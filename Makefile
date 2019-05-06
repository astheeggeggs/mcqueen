all: dmodel dsim dgen dscaling dparameters

IDIR_GSL_HEADERS = libs/gsl-1.16
LIB_GSL=libs/gsl-1.16/.libs/libgsl.a

INCLUDES=-I libs/seq_file -I libs/string_buffer -I libs/cJSON -I $(IDIR_GSL_HEADERS)
LIBDIRS=-Llibs/string_buffer 
CFLAGS=-Wall -Wextra
LINK=-lstrbuf -lpthread -lz -lm

SRCS=$(wildcard src/*.c) libs/cJSON/cJSON.c
HDRS=$(wildcard src/*.h)

ifdef VERBOSE
	CFLAGS := $(CFLAGS) -DVERBOSE=1
endif

ifdef DEBUG
	OPT = -O0 -Wstack-protector -fstack-protector
	DEBUG_ARGS = -g -ggdb -DDEBUG=1
else
	OPT = -O3 # -DNDEBUG=1
	DEBUG_ARGS = 
endif

dmodel: src/tools/dmodel.c $(SRCS) $(HDRS)
	$(CC) $(OPT) $(DEBUG_ARGS) $(CFLAGS) $(INCLUDES) $(LIBDIRS) -I src/ -o dmodel src/tools/dmodel.c $(SRCS) $(LIB_GSL) $(LINK)

dsim: src/tools/dsim.c $(SRCS) $(HDRS)
	$(CC) $(OPT) $(DEBUG_ARGS) $(CFLAGS) $(INCLUDES) $(LIBDIRS) -I src/ -o dsim src/tools/dsim.c $(SRCS) $(LIB_GSL) $(LINK)

dgen: src/tools/dgen.c $(SRCS) $(HDRS)
	$(CC) $(OPT) $(DEBUG_ARGS) $(CFLAGS) $(INCLUDES) $(LIBDIRS) -I src/ -o dgen src/tools/dgen.c $(SRCS) $(LIB_GSL) $(LINK)

dscaling: src/tools/dscaling.c $(SRCS) $(HDRS)
	$(CC) $(OPT) $(DEBUG_ARGS) $(CFLAGS) $(INCLUDES) $(LIBDIRS) -I src/ -o dscaling src/tools/dscaling.c $(SRCS) $(LIB_GSL) $(LINK)

dparameters: src/tools/dparameters.c $(SRCS) $(HDRS)
	$(CC) $(OPT) $(DEBUG_ARGS) $(CFLAGS) $(INCLUDES) $(LIBDIRS) -I src/ -o dparameters src/tools/dparameters.c $(SRCS) $(LIB_GSL) $(LINK)

clean:
	rm -rf dmodel dsim dgen dscaling dparameters *.dSYM *.greg

.PHONY: all clean
