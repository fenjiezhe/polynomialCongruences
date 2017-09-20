#include "inch.h"
#include "stdinc.h"











mat_ZZ reduce_lattice_earlier(mat_ZZ& B){
    double a,b;
    mat_ZZ U;
    a=GetTime();
    G_LLL_RR(B,U,0.99,0,SatNorm);
    b=GetTime();
    cout<<"early abort reduction time: "<<b-a<<"(s)"<<endl;
    return U;
}






void test_off_diag_reduce(){
  mat_ZZ B;
  ifstream fin("/home/tianyuan/Desktop/qtfile/non-qt-pro/size_reduced_mat.txt");
  fin>>B;
  cout<<"before off_diag reduce: "<<endl;
  cout<<B<<endl;
  fin.close();
  cout<<"after off_diag reduce: "<<endl;
  off_diag_reduce(B);
  cout<<B<<endl;
}



/*

void test1(bool use_round=false){
  long h,l;
  cout<<"please enter the RSA prime bit length: "<<endl;
  cin>>l;
  long delta=3;
  ZZ X,N,r;
  ZZ *pp;
  ZZX f;
  mat_ZZ B,Br;


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
  //cout<<B<<endl;

  reduce_lattice_earlier(C);

  cout<<"***********step 3: check root"<<endl;
  if(check_root(X,r,B(1)))
    cout<<"Coppersmith method succeeds."<<endl;
  else cerr<<"fail!!!!"<<endl;

  if(check_root(X,r,RL))
    cout<<"Early abort succeeds."<<endl;
  else cerr<<"fail!!!!"<<endl;




  /////////////////////////////////////////////////
  if(use_round){

    mat_ZZ U;
    vec_ZZ v;
    cout<<endl;
    cout<<"==============with rounding=============="<<endl;

    Br=gen_lattice(h,fc,N,pp[1]);



    off_diag_reduce(Br);
    mat_ZZ reB(Br);
    rounding(Br,pp[1],delta);
    U=reduce_lattice(Br);
    mul(v,U(1),reB);
    if(v.length()<7){
      cout<<"v(x) = "<<v<<endl;
    }
    if(check_root(pp[1],r,v))
      cout<<"Coppersmith rounding succeeds."<<endl;
    else cerr<<"fail!!!!"<<endl;
  }
}
*/
void  test_fplll_correctness(){
  ifstream fin("Xr.txt");
  ZZ X,r;
  fin>>X>>r;
  cout<<"X = "<<X<<endl;
  cout<<"r = "<<r<<endl;
  fin.close();
  fin.open("time.txt");
  mat_ZZ B;
  fin>>B;
  cout<<"B(1) = "<<B(1)<<endl;
  if(check_root(X,r,B(1)))
       cout<<"Coppersmith method succeeds."<<endl;
  else cerr<<"fail!!!!"<<endl;
}




int main()
{
  //test0();
  //test1(true);
  //test_fplll_correctness();
  //test_19141();

  param_poly_lattice2();
  //cout<<"test X1 and X2\n"<<endl;
  //test_param();
  //test_off_diag_reduce();
   return 0;
}





