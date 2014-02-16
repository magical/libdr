#define DR_VERSION "0.1.0"

struct dr_state;

struct dr_state* dr_create(unsigned channels, unsigned sample_rate);
void dr_destroy(struct dr_state**);
/*int dr_reset(struct dr_state*, unsigned channels, unsigned sample_rate);*/

int dr_add_frames_int16(struct dr_state*, const int16_t* data, size_t frames);
int dr_add_frames_int32(struct dr_state*, const int32_t* data, size_t frames);
int dr_add_frames_float(struct dr_state*, const float* data, size_t frames);
int dr_add_frames_double(struct dr_state*, const double* data, size_t frames);
int dr_finish(struct dr_state*, double* drp, double* peakp);

unsigned dr_channels(struct dr_state*);
