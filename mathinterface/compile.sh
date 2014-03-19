#!/bin/bash

mlinkdir="/usr/local/Wolfram/Mathematica/9.0/SystemFiles/Links/MathLink/DeveloperKit"
mlinkfiles="${mlinkdir}/Linux-x86-64/CompilerAdditions"
mprep="${mlinkfiles}/mprep"

includes="${mlinkfiles}"
libdir="${mlinkfiles}"
libraries="-lML64i3 -lm -lpthread -lrt -lstdc++"

# clean
rm -f minterface.c minterface.o feed2mathematica.o feed2mathematica

# compile mathematica glue code
"${mprep}" feed2mathematica.tm -o minterface.c
gcc -c -I"${includes}" minterface.c

g++ -c -I"${includes}" feed2mathematica.cpp

g++ feed2mathematica.o minterface.o -L"${libdir}" $libraries -o feed2mathematica

