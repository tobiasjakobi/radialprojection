compiler := g++
flags    :=

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
           tuebingen

all: build;

%: %.cpp; $(compiler) -c -o $@ $(flags) $<

build: $(objects)

clean:
	rm -f *.o
	rm -f $(objects)

strip:
	strip -s $(objects)

