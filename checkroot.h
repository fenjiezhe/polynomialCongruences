#ifndef CHECKROOT_H
#define CHECKROOT_H

#endif // CHECKROOT_H
#include "stdinc.h"

//vec_ZZ b1 为 f(xX)对应的向量
bool check_root(const ZZ X,const ZZ r,const vec_ZZ b1);

ZZ evaluation(ZZX &f, ZZ a);

ZZ norm2(const vec_ZZ& v);

bool GHbound_satisfied(const ZZ N, const long h,const long n,const vec_ZZ g);


