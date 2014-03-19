compiler    := gcc
extra_flags :=

ifndef platform
platform := $(shell $(compiler) -dumpmachine)
endif

ifeq (release,$(build))
extra_flags += -O2 -march=native
endif

ifeq (debug,$(build))
extra_flags += -O0 -g
endif

cc      := $(compiler) -std=gnu99
cpp     := $(subst cc,++,$(compiler)) -std=gnu++0x
flags   := -fPIC $(extra_flags)
ldflags := -shared -Wl,--version-script=link.T
objects :=
           
objects += cyclotomic_radial chiral_radial corner_cases \
           griddual hexagonal higher_cyclo histogram \
           tuebingen

all: build;

%.o: %.cpp
	$(cpp) -c -o $@ $(flags) $<

build: $(objects)

clean:
	rm -f *.o
	rm -f $(objects)

strip:
	strip -s $(objects)

