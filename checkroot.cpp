#include "checkroot.h"
#include "gen.h"

bool check_root(const ZZ X,const ZZ r,const vec_ZZ b1)
{
  ZZX v;
  vec_ZZ b(b1); //f(xX)对应的向量
  inverse_xX(b,X);//还原为f(x)对应的向量
  conv(v,b); //向量转换为多项式
  if(IsZero(evaluation(v,r))) return true;
  else return false;
}

//算多项式在a点的值f(a)
ZZ evaluation(ZZX &f, ZZ a)
{
    int n = deg(f);
    int i;
    ZZ sum;
    sum =0;
    for(i=0; i<=n; i++)
    {
        add(sum,sum, coeff(f,i)*power(a,i));
    }
    return sum;
}

ZZ norm2(const vec_ZZ& v)
{
  ZZ sqr_sum(0);
  for(int k=1;k<=v.length();k++){
    add(sqr_sum,sqr_sum,sqr(v(k)));
  }
  return sqr_sum;
}

bool GHbound_satisfied(const ZZ N, const long h, const long n, const vec_ZZ g)
{
  ZZ bnd;
  FloorToZZ(bnd,to_RR(power(N,2*h-2))/to_RR(n));
  ZZ lhs;
  lhs=norm2(g);
  return lhs<bnd;
}
