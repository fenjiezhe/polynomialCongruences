#ifndef REDUCE_H
#define REDUCE_H

#endif // REDUCE_H

#include "stdinc.h"


extern ZZ lbound;
extern vec_ZZ RL;


long SatNorm(const vec_ZZ& z);


mat_ZZ reduce_lattice(mat_ZZ& B);   //return unimodular matrix

void off_diag_reduce(mat_ZZ& B);

void rounding(mat_ZZ& B, const ZZ X, long delta=3, double cc=1.5);
