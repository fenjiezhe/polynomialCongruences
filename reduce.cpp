#include "reduce.h"
#include "stdinc.h"


ZZ lbound;
vec_ZZ RL;



//return unimodular matrix
//append reduced matrix to matrix.txt

mat_ZZ reduce_lattice(mat_ZZ& B)
{
  ofstream fout;
  fout.open("/home/tianyuan/Desktop/qtfile/non-qt-pro/reduced_mat.txt");
  ofstream ofs;
  ofs.open ("/home/tianyuan/Desktop/qtfile/non-qt-pro/param.txt",\
              ofstream::out | ofstream::app);

    double a,b;
    mat_ZZ U;
    a=GetTime();
    G_LLL_RR(B,U);
    b=GetTime();

    fout<<B<<endl;


    //这里考虑子矩阵维度应取多少
    vec_ZZ u1(U(1));
    int count=0;
    for(int k=u1.length();k>0;k--){
      if(u1(k)==0){
       count++;
       cout<<k<<endl;
      }
    }

    cout<<"reduction time: "<<b-a<<"(s)"<<endl;
    cout<<"length of U(1) is "<<u1.length()<<endl;

    cout<<"The first "<<u1.length()-count<<" tuples are non-zero."<<endl;

    ofs<<"\nreduction time: "<<b-a<<"(s)"<<endl;
    ofs<<"The first "<<u1.length()-count<<" tuples are non-zero."<<endl;

    ofs<<endl;
    ofs<<"length of U(1) is "<<u1.length()<<endl;
    //ofs << u1<<endl;

      ofs.close();
      fout.close();

    return U;
}


void off_diag_reduce(mat_ZZ& B) //size-reduce lower-triangular matrix
{
  ZZ q,r;
  ZZ a,b;
  for(long i=B.NumRows()-1;i>0;i--){
    b=B(i,i);
    //run through row(j)
    for(long j=i+1;j<=B.NumRows();j++){
      a=B(j,i);
      DivRem(q,r,a,b);
      vec_ZZ sc_vec;
      mul(sc_vec,q,B(i));
      //sc_vec=q*B(i);
      sub(B(j),B(j),sc_vec);
      //B(j)-=sc_vec;
     }
  }


  ofstream fout;
  fout.open("/home/tianyuan/Desktop/qtfile/non-qt-pro/size_reduced_mat.txt");
  fout<<B<<endl;
  fout.close();
}


void rounding(mat_ZZ& B,const ZZ X,long delta,double cc)
{
  RR c;
  RR t;
  long n=B.NumRows();

  ofstream fout;
  fout.open("/home/tianyuan/Desktop/qtfile/non-qt-pro/reduced_mat.txt");

  //cout<<"X^(n-delta) = "<<power(X,n-delta)<<endl;

  c=power(to_RR(cc),n);   //c=1.5^n
  cout<<"c = "<<c<<endl;

  RR k;
  k=c/to_RR(power(X,n-delta));
  cout<<"c/(X^(n-delta)) = "<<k<<endl;

  for(int i=1;i<=n;i++){
    for(int j=1;j<=i;j++){
      t=k*to_RR(B(i,j));
      RoundToZZ(B(i,j),t);
    }
  }
  fout<<B<<endl;
  fout.close();

}

//提前终止：算各系数绝对值之和
long SatNorm(const vec_ZZ& z)
{
   long n = z.length();
   long j;
   ZZ sum;
   sum = 0;

   for (j = 0; j < n; j++)
   {
      add(sum,sum,abs(z[j]));
   }

   if(sum<lbound)
   {
       RL = z;
       return 1;
   }
   else
       return 0;
}



