compiler := g++
cflags   := -I./include -Wall -Wno-narrowing
ldflags  :=
sseflags :=

ifndef platform
platform := $(shell $(compiler) -dumpmachine)
endif

ifeq (1,$(enable_sse))
sseflags += -DCOMMON_USE_SSE
endif

ifeq (release-lto,$(build))
cflags += -O2 -march=native -flto=4 -fuse-linker-plugin -DNDEBUG $(sseflags)
ldflags += -O2 -march=native -flto=4 -fuse-linker-plugin
endif

ifeq (release,$(build))
cflags += -O2 -march=native -DNDEBUG $(sseflags)
endif

ifeq (debug,$(build))
cflags += -O0 -g
endif

# Enables "debugging" of the local visibility test.
# In fact this intentionally breaks the visibility test for the
# dodecagonal tiling, useful to illustrate the need to properly
# shift the window (see the corresponding Mathematica worksheet).
ifeq (1,$(visibility_debug))
cflags += -DDODECAGONAL_VISIBILITY_DEBUG
endif

# Computes the slope of a point instead of the angle, so omitting
# the application of the arctan().
# Currently we only consider this useful for debugging purpose.
# Also only implemented in some modules.
ifeq (1,$(slope_debug))
cflags += -DSLOPE_ANGLE_DEBUG
endif

objects := cyclotomic_radial chiral_radial corner_cases \
           griddual hexagonal higher_cyclo histogram \
           tuebingen arith_visibility pdf_writer random \
           cyclotomic_random

all: $(objects)

%.o: src/%.cpp
	$(compiler) -c -o $@ $(cflags) $<

cyclotomic_radial: cyclotomic_octagonal.o cyclotomic_decagonal.o \
                   cyclotomic_dodecagonal.o cyclotomic_rhombic_penrose.o \
                   common.o cyclotomic_radial.o; $(compiler) -o $@ $(ldflags) $^
chiral_radial: chiral_radial.o common.o; $(compiler) -o $@ $(ldflags) $^
corner_cases: corner_cases.o common.o; $(compiler) -o $@ $(ldflags) $^
griddual: griddual.o common.o; $(compiler) -o $@ $(ldflags) $^
hexagonal: hexagonal.o common.o; $(compiler) -o $@ $(ldflags) $^
higher_cyclo: higher_cyclo.o common.o; $(compiler) -o $@ $(ldflags) $^
histogram: histogram.o common.o; $(compiler) -o $@ $(ldflags) $^
tuebingen: tuebingen.o common.o; $(compiler) -o $@ $(ldflags) $^
arith_visibility: arith_visibility.o common.o; $(compiler) -o $@ $(ldflags) $^
pdf_writer: pdf_writer.o common.o; $(compiler) -o $@ $(ldflags) $^ -lcairo
random: random.o common.o; $(compiler) -o $@ $(ldflags) $^
cyclotomic_random: cyclotomic_random.o cyclotomic_octagonal.o \
                   cyclotomic_decagonal.o cyclotomic_dodecagonal.o \
                   cyclotomic_rhombic_penrose.o \
                   common.o; $(compiler) -o $@ $(ldflags) $^

clean:
	rm -f *.o
	rm -f $(objects)

strip:
	strip -s $(objects)

rebuild: clean all
