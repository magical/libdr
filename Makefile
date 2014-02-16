CFLAGS=-O3 -Wall -Wextra -Wconversion -Wdeclaration-after-statement -finline-functions

all: libdr.a
libdr.a: libdr.a(dr.o) Makefile
