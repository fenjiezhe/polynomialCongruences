#ifndef PARAM_H
#define PARAM_H

#endif // PARAM_H

#include "stdinc.h"



long gen_N(long l,ZZ& N); //l=素数bit数，N=pq,返回N的bit数

ZZ gen_X1(const ZZ N,const long h,const long delta=3);
ZZ gen_X2(const ZZ X1, const long h, const long delta=3, double t=1.5);

//可以通过用户输入指定h的值
long set_h(const ZZ N, const long delta=3);


ZZ *gen_param(long &h,const ZZ N,const long delta);


void xX(ZZX& f,ZZ X); //f(x)变为f(xX)

void inverse_xX(vec_ZZ& f,ZZ X); //f(xX)变为f(x)

void gen_poly(ZZ_pX& f,long fd,const ZZ r); //多项式，次数，根

mat_ZZ gen_lattice(const long h, const ZZX ff, const ZZ N, const ZZ X); //返回格基矩阵

