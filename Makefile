compiler := g++
cflags   := -I./include -Wall -Wno-narrowing
ldflags  :=

ifndef platform
platform := $(shell $(compiler) -dumpmachine)
endif

ifeq (release-lto,$(build))
cflags += -O2 -march=native -flto=4 -fuse-linker-plugin -DNDEBUG
ldflags += -O2 -march=native -flto=4 -fuse-linker-plugin
endif

ifeq (release,$(build))
cflags += -O2 -march=native -DNDEBUG
endif

ifeq (debug,$(build))
cflags += -O0 -g
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
