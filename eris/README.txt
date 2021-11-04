To use a large number of particles:

Specifically, edit the file include/armadillo_bits/config.hpp and uncomment the line with: // #define ARMA_64BIT_WORD. In version 3.4 this should be near line 59.

Alternatively, you can define ARMA_64BIT_WORD before including the Armadillo header in your program, eg:

#define ARMA_64BIT_WORD
#include <armadillo>
#include <iostream>
...
