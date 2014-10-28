compiler := g++
flags    := -I./include

ifndef platform
platform := $(shell $(compiler) -dumpmachine)
endif

ifeq (release,$(build))
flags += -O2 -march=native
endif

ifeq (debug,$(build))
flags += -O0 -g
endif

objects := cyclotomic_radial chiral_radial corner_cases \
           griddual hexagonal higher_cyclo histogram \
           tuebingen arith_visibility pdf_writer

all: $(objects)

%.o: src/%.cpp
	$(compiler) -c -o $@ $(flags) $<

cyclotomic_radial: cyclotomic_octagonal.o cyclotomic_decagonal.o \
                   cyclotomic_dodecagonal.o cyclotomic_rhombic_penrose.o \
                   common.o cyclotomic_radial.o; $(compiler) -o $@ $(flags) $^
chiral_radial: chiral_radial.o common.o; $(compiler) -o $@ $(flags) $^
corner_cases: corner_cases.o common.o; $(compiler) -o $@ $(flags) $^
griddual: griddual.o common.o; $(compiler) -o $@ $(flags) $^
hexagonal: hexagonal.o common.o; $(compiler) -o $@ $(flags) $^
higher_cyclo: higher_cyclo.o common.o; $(compiler) -o $@ $(flags) $^
histogram: histogram.o common.o; $(compiler) -o $@ $(flags) $^
tuebingen: tuebingen.o common.o; $(compiler) -o $@ $(flags) $^
arith_visibility: arith_visibility.o common.o; $(compiler) -o $@ $(flags) $^
pdf_writer: pdf_writer.o common.o; $(compiler) -o $@ $(flags) $^ -lcairo

clean:
	rm -f *.o
	rm -f $(objects)

strip:
	strip -s $(objects)
