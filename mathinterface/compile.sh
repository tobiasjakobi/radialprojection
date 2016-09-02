#!/bin/bash

## Adjust the directory depending on which Mathematica
## version you target.
mlinkdir="/var/tmp/Mathematica/9.0.1/SystemFiles/Links/MathLink/DeveloperKit"
#mlinkdir="/var/tmp/Mathematica/10.3/SystemFiles/Links/MathLink/DeveloperKit"
#mlinkdir="/var/tmp/Mathematica/11.0/SystemFiles/Links/MathLink/DeveloperKit"

mlinkfiles="${mlinkdir}/Linux-x86-64/CompilerAdditions"
mprep="${mlinkfiles}/mprep"

includes="${mlinkfiles}"
libdir="${mlinkfiles}"

## Mathematica9 still uses libML64i3, while later
## versions use the i4 version.
libraries="-lML64i3 -lm -lpthread -lrt -lstdc++"
#libraries="-lML64i4 -lm -lpthread -lrt -lstdc++"

# clean
rm -f minterface.c minterface.o feed2mathematica.o feed2mathematica

# compile mathematica glue code
"${mprep}" feed2mathematica.tm -o minterface.c
gcc -c -I"${includes}" minterface.c

g++ -c -I"${includes}" feed2mathematica.cpp

g++ feed2mathematica.o minterface.o -L"${libdir}" $libraries -o feed2mathematica
