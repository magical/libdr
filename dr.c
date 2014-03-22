/*
This file is part of libdr. Copyright 2014 Andrew Ekstedt. See LICENSE for
terms of use. NO WARRANTY.
*/

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <math.h>
#include "dr.h"

enum dr_sample_type {
	int16_type,
	int32_type,
	float_type,
	double_type,
};

static int dr_grow(struct dr_state*);
static int dr_add_frames_internal(struct dr_state*, const void* data, size_t frames, enum dr_sample_type);
static double dr_get_sample(struct dr_state*, const void* data, size_t i, unsigned ch, enum dr_sample_type);

static int debug = 0;

// Terminology note: a frame consisists of one sample for each channel.

struct dr_state {
	unsigned int channels;
	unsigned int rate;

	// partial sum and peak of each channel of current fragment
	double* channel_sum;
	double* channel_peak;

	// rms and peak of each fragment, per channel
	double** frag_rms;
	double** frag_peak;

	uint64_t nfrag; // number of fragments processed
	uint64_t fragalloc; // number of fragments allocated
	uint64_t fragpos; // number of frames processes in the current fragment
	uint64_t fraglen; // number of frames in a fragment
};

// Callback for qsort
static int doublecmp(const void *va, const void *vb)
{
	double a = *(double *)va;
	double b = *(double *)vb;
	return (a<b) - (a>b); // reverse order
}

struct dr_state*
dr_create(unsigned int channels, unsigned int rate)
{
	struct dr_state* m;
	unsigned int ch;

	m = malloc(sizeof *m);
	if(m == NULL) {
		return NULL;
	}

	m->channels = channels;
	m->rate = rate;

	m->frag_rms = malloc(channels * sizeof m->frag_rms[0]);
	m->frag_peak = malloc(channels * sizeof m->frag_peak[0]);
	m->channel_sum = malloc(channels * sizeof m->channel_sum[0]);
	m->channel_peak = malloc(channels * sizeof m->channel_peak[0]);

	if(m->channel_sum == NULL || m->channel_peak == NULL || m->frag_rms == NULL || m->frag_peak == NULL) {
		goto cleanup;
	}

	for(ch = 0; ch < channels; ch++) {
		m->frag_rms[ch] = NULL;
		m->frag_peak[ch] = NULL;
		m->channel_sum[ch] = 0;
		m->channel_peak[ch] = 0;
	}

	m->nfrag = 0;
	m->fragalloc = 0;
	m->fragpos = 0;
	m->fraglen = 3 * rate;
	if(rate == 44100) {
		m->fraglen = 3 * 44160;
	}

	return m;

cleanup:
	free(m->frag_rms);
	free(m->frag_peak);
	free(m->channel_sum);
	free(m->channel_peak);
	free(m);
	return NULL;
}

void
dr_destroy(struct dr_state** mp)
{
	unsigned ch;
	if(mp) {
		struct dr_state* m = *mp;
		for(ch = 0; ch < m->channels; ch++) {
			if(m->frag_rms) free(m->frag_rms[ch]);
			if(m->frag_peak) free(m->frag_peak[ch]);
		}
		free(m->frag_rms);
		free(m->frag_peak);
		free(m->channel_sum);
		free(m->channel_peak);
		free(m);
	}
}

unsigned dr_channels(struct dr_state* m) {
	return m->channels;
}

static double
dr_get_sample(struct dr_state* m, const void* data, size_t i, unsigned ch, enum dr_sample_type type)
{
	i = i * m->channels + ch;
	switch(type) {
	case int16_type: return (double)((int16_t*)data)[i] * (double)INT16_MAX;
	case int32_type: return (double)((int32_t*)data)[i] * (double)INT32_MAX;
	case float_type: return (double)((float*)data)[i];
	case double_type: return ((double*)data)[i];
	}
	abort();
}

int dr_add_frames_int16(struct dr_state* m, const int16_t* data, size_t frames) {
	return dr_add_frames_internal(m, data, frames, int16_type);
}
int dr_add_frames_int32(struct dr_state* m, const int32_t* data, size_t frames) {
	return dr_add_frames_internal(m, data, frames, int32_type);
}
int dr_add_frames_float(struct dr_state* m, const float* data, size_t frames) {
	return dr_add_frames_internal(m, data, frames, float_type);
}
int dr_add_frames_double(struct dr_state* m, const double* data, size_t frames) {
	return dr_add_frames_internal(m, data, frames, double_type);
}

int
dr_add_frames_internal(struct dr_state* m, const void* data, size_t frames, enum dr_sample_type type)
{
	double v;
	unsigned ch;
	size_t i, start, end;
	int err;

	if(debug) {
		fprintf(stderr, "dr_add_frames: Adding %zd frames\n", frames);
	}

	for(i = 0; i < frames;) {
		if(m->fragpos == 0) {
			for(ch = 0; ch < m->channels; ch++) {
				m->channel_sum[ch] = 0;
				m->channel_peak[ch] = 0;
			}
		}

		start = i;
		end = i + (m->fraglen - m->fragpos);
		if(end > frames) {
			end = frames;
		}
		for(; i < end; i++) {
			for(ch = 0; ch < m->channels; ch++) {
				v = dr_get_sample(m, data, i, ch, type);
				m->channel_sum[ch] += v * v;
				v = fabs(v);
				if(m->channel_peak[ch] < v) {
					m->channel_peak[ch] = v;
				}
			}
		}

		m->fragpos += end - start;
		if(debug) {
			fprintf(stderr, "dr_add_frames: pos %ld, len %ld\n", m->fragpos, m->fraglen);
		}
		if(m->fragpos == m->fraglen) {
			err = dr_grow(m);
			if(err) {
				return err;
			}
			for(ch = 0; ch < m->channels; ch++) {
				m->frag_rms[ch][m->nfrag] =
					sqrt(2.0 * m->channel_sum[ch] / (double)m->fraglen);
				m->frag_peak[ch][m->nfrag] = m->channel_peak[ch];
			}
			m->nfrag++;
			m->fragpos = 0;
		}
	}
	return 0;
}

static int
dr_grow(struct dr_state* m)
{
	uint64_t chunk = 10000;
	uint64_t alloc = m->fragalloc;
	unsigned ch;
	void *p;
	if(m->nfrag < alloc) {
		return 0;
	}
	while(alloc <= m->nfrag) {
		alloc += chunk;
	}
	if(debug) {
		fprintf(stderr, "dr_grow: growing to %"PRIu64" fragments\n", alloc);
	}
	for(ch = 0; ch < m->channels; ch++) {
		p = realloc(m->frag_rms[ch], alloc * sizeof m->frag_rms[ch][0]);
		if(p == NULL) {
			return 1;
		}
		m->frag_rms[ch] = p;
		p = realloc(m->frag_peak[ch], alloc * sizeof m->frag_peak[ch][0]);
		if(p == NULL) {
			return 1;
		}
		m->frag_peak[ch] = p;
	}
	m->fragalloc = alloc;
	return 0;
}

int
dr_finish(struct dr_state* m, double* drp, double* peakp)
{
	double dr, ch_dr, peak, ch_peak, v;
	size_t i, n;
	unsigned ch;
	int err;

	if(debug) {
		printf("Total samples: %"PRIu64"\n", m->nfrag*m->fraglen + m->fragpos);
	}

	// Finish incomplete fragment...
	err = dr_grow(m);
	if(err) {
		return err;
	}
	for(ch = 0; ch < m->channels; ch++) {
		m->frag_rms[ch][m->nfrag] =
			sqrt(2.0 * m->channel_sum[ch] / (double)m->fraglen);
		m->frag_peak[ch][m->nfrag] = m->channel_peak[ch];
	}
	m->nfrag++;
	//m->fragpos = 0;

	if(m->nfrag == 0) {
		*drp = 0;
		*peakp = 0;
		return 0;
	}

	dr = 0;
	peak = 0;
	for(ch = 0; ch < m->channels; ch++) {
		qsort(m->frag_rms[ch], m->nfrag, sizeof m->frag_rms[ch][0], doublecmp);
		n = m->nfrag / 5;
		ch_dr = 0;
		ch_peak = 0;
		for(i = 0; i < n; i++) {
			v = m->frag_rms[ch][i];
			ch_dr += v * v;
		}

		qsort(m->frag_peak[ch], m->nfrag, sizeof m->frag_peak[ch][0], doublecmp);
		if(m->nfrag > 1) {
			ch_peak = m->frag_peak[ch][1];
		} else {
			ch_peak = m->frag_peak[ch][0];
		}

		ch_dr = sqrt(ch_dr / (double)n) / ch_peak;
		if(peak < ch_peak) {
			peak = ch_peak;
		}
		if(debug) {
			printf("ch%u DR: %6.2f   Peak: %6.2f\n", ch, 20*log10(ch_dr), 20*log10(ch_peak));
		}
		dr += 20 * log10(ch_dr);
	}

	dr = dr / m->channels;

	// score = (int) round(dr);
	if(drp) *drp = dr;
	if(peakp) *peakp = peak;
	return 0;
}
