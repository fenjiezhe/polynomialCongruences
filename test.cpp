#include "test.h"
#include "gen.h"
#include "reduce.h"
#include "checkroot.h"

//l,N,f_N,h,X,f,B
void param_poly_lattice()
{
  long h;
  ZZ N,X,r;

  cout<<"enter RSA prime length:"<<endl;
  long l;
  cin>>l;

  long delta=3;
  gen_N(l,N);

  ZZ_p::init(N);
  ZZ_pX f_N;

  h=set_h(N);
  X=gen_X1(N,h);
  r = X-1;
  gen_poly(f_N,delta,r);


  mat_ZZ B;
  ZZX f;
  conv(f,f_N);
  B = gen_lattice(h,f,N,X);
  reduce_lattice(B);

  if(check_root(X,r,B(1)))
    cout<<"Coppersmith method succeeds."<<endl;
  else cerr<<"fail!!!!"<<endl;

}


//l,N,f_N,h,X,f,B
void param_poly_lattice2()
{
  ofstream ofs;
  ofs.open ("/home/tianyuan/Desktop/qtfile/non-qt-pro/param.txt",\
              ofstream::out | ofstream::app);
  long h;
  ZZ N,X,X2,r;

  cout<<"enter RSA prime length:"<<endl;
  long l;
  cin>>l;

  long delta=3;
  gen_N(l,N);

  ZZ_p::init(N);
  ZZ_pX f_N;

  h=set_h(N);
  X=gen_X1(N,h);
  X2=gen_X2(X,h);
  r = X2 - 1;
  gen_poly(f_N,delta,r);


  mat_ZZ B;
  ZZX f;
  conv(f,f_N);
  B = gen_lattice(h,f,N,X);
  mat_ZZ BB;
  BB = gen_lattice(h,f,N,X2);


  reduce_lattice(B);

  if(check_root(X,r,B(1)))
   {
    cout<<"\nCoppersmith method succeeds."<<endl;
    ofs<<"\nCoppersmith method succeeds."<<endl;
  }
  else{

    cerr<<"\nfails!!!!"<<endl;
    ofs<<"\nRounding fails!!!!"<<endl;
  }

  off_diag_reduce(BB);
  mat_ZZ BBor(BB);
  rounding(BB,X2);
  mat_ZZ UU;
  vec_ZZ vv;
  UU=reduce_lattice(BB);
  mul(vv,UU(1),BBor);

  if(check_root(X2,r,vv))
   { cout<<"\nCoppersmith method with rounding succeeds."<<endl;
    ofs<<"\nRounding succeeds."<<endl;

  }
  else{

    cerr<<"\nRounding fails!!!!"<<endl;
    ofs<<"\nRounding fails!!!!"<<endl;
    cerr<<"GH bound satisfied ? "<<GHbound_satisfied(N,h,h*delta,vv)<<endl;
    ofs<<"GH bound satisfied ? "<<GHbound_satisfied(N,h,h*delta,vv)<<endl;

  }

  ZZX ff;
  conv(ff,vv);
  //cout<<"short vector "<<ff<<endl;
  Vec< Pair< ZZX, long > > factors;
  ZZ dd;

  factor(dd, factors, ff);
  cout<<"\nterms num: "<<endl;
  cout << factors.length() << "\n";

  ofs<<"\nterms num: "<<endl;
  ofs << factors.length() << "\n";

  ofs.close();

}

/*void test_param()
{
  long h;
  cout<<"enter RSA prime length:"<<endl;
  long l;
  cin>>l;
  ZZ N,X1,X2;
  gen_N(l,N);

  h=set_h(N);
  X1=gen_X1(N,h);
  X2=gen_X2(X1,h);
  r = X1-1;
  gen_poly(f_N,delta,r);


  mat_ZZ B;
  ZZX f;
  conv(f,f_N);
  B = gen_lattice(h,f,N,X);
  reduce_lattice(B);

  if(check_root(X,r,B(1)))
    cout<<"Coppersmith method succeeds."<<endl;
  else cerr<<"fail!!!!"<<endl;


}*/

/*
void poly_param_lattice(){
  long l;
  long delta=3;
  ZZ X,N,r;
  ZZ *pp;
  ZZX f;
  mat_ZZ B;


  cout<<"***********step 0: generate parameters"<<endl;
  gen_N(l,N);
  pp=gen_param(h,N,delta);
  X=pp[0];
  r=pp[1]-1;

  ofstream xr("Xr.txt");
  xr<<X<<"\n"<<r<<endl;
  xr.close();

  ZZ_p::init(N);
  ZZ_pX f_p;

  cout<<"***********step 1: generate polynomial"<<endl;
  gen_poly(f_p,delta,r);
  conv(f,f_p);
  ZZX fc(f);
  cout<<"fc = "<<fc<<endl;
  cout<<"***********step 2: generate lattice basis and reduce it"<<endl;
  B=gen_lattice(h,f,N,pp[0]);
  //cout<<B<<endl;

  mat_ZZ C(B);
cout<<"reduction"<<endl;
  reduce_lattice(B);
}

*/

//测试f(x)=[19 14 1]
//v(x) = [3 8 -24 -8 -1 2]
void test_19141()
{
  ZZX f;
  ZZ X; X=2;
  ZZ N;N=35;
  SetCoeff(f,0,19); SetCoeff(f,1,14); SetCoeff(f,2,1);
  cout<<"f(x) = "<<f<<endl;
  long h=3;

  mat_ZZ p;
  p=gen_lattice(h,f,N,X);

  reduce_lattice(p);




  /*cout<<"before rounding:\n";
  cout<<p<<endl;
  rounding(p,X,2);
  cout<<"after rounding:\n";
  cout<<p<<endl;*/

}

/*
f(x) = [19 14 1]
before reduction:
[[1225 0 0 0 0 0]
[0 2450 0 0 0 0]
[665 980 140 0 0 0]
[0 1330 1960 280 0 0]
[361 1064 936 224 16 0]
[0 722 2128 1872 448 32]
]

after reduction
reduction time: 0
[[3 16 -96 -64 -16 64]
[49 100 0 160 0 64]
[-128 50 -20 80 160 32]
[-201 8 132 -32 -48 32]
[-83 -142 52 8 32 160]
[61 32 148 -128 48 128]
]
v(x) = [3 8 -24 -8 -1 2]
 */
