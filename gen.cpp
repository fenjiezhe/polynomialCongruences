#include "gen.h"

//l:RSA prime bit length
//delta:polynomial degree
//if asymp = true, then h takes asymptotic value;else h is assigned by user
//return X1 and X2(rounding bound)

//不要在声明和定义都放参数默认值



long gen_N(long l,ZZ& N)
{
  ZZ p,q;
  ofstream ofs;
  ofs.open ("/home/tianyuan/Desktop/qtfile/non-qt-pro/param.txt",\
              ofstream::out | ofstream::trunc);
  //Any contents that existed in the file before it is open are discarded.

  p=GenPrime_ZZ(l);
  q=GenPrime_ZZ(l);
  N=p*q;
  cout<<"N = "<<N<<endl;
  ofs<<"N = "<<N<<endl;

  cout<<"size of N is "<<NumBits(N)<<endl;
  ofs<<"size of N is "<<NumBits(N)<<endl;


  ofs.close();

  return NumBits(N);
}


void gen_poly(ZZ_pX& f,long fd,const ZZ r)
{
  SetCoeff(f,0,to_ZZ_p(-r));
  SetCoeff(f,1,1);


  ofstream ofs;
  ofs.open ("/home/tianyuan/Desktop/qtfile/non-qt-pro/param.txt",\
              ofstream::out | ofstream::app);

  // g is a random polynomial with degree fd-1
  ZZ_pX g;
  while(1)
  {
    random(g,fd);
    if(deg(g)==fd-1)
      break;
  }

  //f=(x-r)*g
  mul(f,f,g);

  //make f monic
  MakeMonic(f);

  cout<<"poly = "<<f<<endl;
  ofs<<"poly = "<<f<<endl;

  cout<<"root = "<<r<<endl;
  ofs<<"root = "<<r<<endl;

  ofs.close();

}


long set_h(const ZZ N,const long delta)
{
  cout<<"set h, if you want h takes asymptotic value, enter 0."<<endl;
  long h;
  cin>>h;


  ofstream ofs;
  ofs.open ("/home/tianyuan/Desktop/qtfile/non-qt-pro/param.txt",\
              ofstream::out | ofstream::app);

  if(!h){
    h=NumBits(N)/delta;
    cout<<"h = O(logN/delta) = "<<h<<endl;
    ofs<<"h = O(logN/delta) = "<<h<<endl;
  }
  else {
    cout<<"assigned by user, h = "<<h<<endl;
    ofs<<"assigned by user, h = "<<h<<endl;

    cout<<"h = O(logN/delta) = "<<NumBits(N)/delta<<endl;
    ofs<<"h = O(logN/delta) = "<<NumBits(N)/delta<<endl;

  }
  cout<<"n = h*delta = "<<h*delta<<endl;
  ofs<<"n = h*delta = "<<h*delta<<endl;

  ofs.close();
  return h;
}

ZZ gen_X1(const ZZ N,const long h,const long delta)
{
  ZZ Y;
  long m=h*delta;


  ofstream ofs;
  ofs.open ("/home/tianyuan/Desktop/qtfile/non-qt-pro/param.txt",\
              ofstream::out | ofstream::app);

  //floor(2^{-0.5}*N^{(h-1)/(m-1)}*m^{-1/(m-1)};

  RoundToZZ(Y,pow(to_RR(N),to_RR(h-1)/to_RR(m-1))/pow(to_RR(m),\
  1/to_RR(m-1))/pow(to_RR(2.0), to_RR(0.5)));

  cout<<"size of X: "<<NumBits(Y)<<endl;
  cout<<"X = "<<Y<<endl;


  ofs<<"size of X: "<<NumBits(Y)<<endl;
  ofs<<"X = "<<Y<<endl;
  ofs.close();

  return Y;

}

ZZ gen_X2(const ZZ X1,const long h,const long delta,double t)
{
  long n=h*delta;
  RR k1,c,ratio;
  RR san(3);
  RR er(2);
  c=pow(to_RR(t),to_RR(n)); //c=1.5^n
  k1=pow(to_RR(n),to_RR(1.5))*\
      pow((san*c-er)/(er*c-er),to_RR(n-1));
  ZZ cc;
  RoundToZZ(cc,c);
  k1/=to_RR(cc);
  k1+=to_RR(1);
  ratio = pow(k1,-er/to_RR(n-1));
  cout<<"ratio = "<<ratio<<endl;
  ZZ X2;
  RoundToZZ(X2,ratio*to_RR(X1));
  cout<<"X2 = "<<X2<<endl;
  return X2;
}

//deprecated
ZZ* gen_param(long &h,const ZZ N,const long delta)
{

  ZZ X,Y;
  RR c;
  ZZ* pp;
  pp = new ZZ[2];

  cout<<"*****************"<<endl;
  //cout<<"delta = "<<delta<<endl;

  if(!h){
    h=NumBits(N)/delta;
    cout<<"h = O(logN/delta) = "<<h<<endl;
  }
  else cout<<"assigned by user, h = "<<h<<endl;

  //dimension of lattice
  long m=h*delta;

  cout<<"dim = h*delta = "<<m<<endl;

    //X = floor(2^{-0.5}*N^{(h-1)/(m-1)}*m^{-1/(m-1)};
RoundToZZ(Y,pow(to_RR(N),to_RR(h-1)/to_RR(m-1))/pow(to_RR(m),1/to_RR(m-1))/pow(to_RR(2.0), to_RR(0.5)));



    pp[0]= Y;


    RR k;
    RR temp1,temp2,temp3,temp4;
    RR RX;
    c=pow(to_RR(1.5),to_RR(m)); //c=1.5^n
    temp1=pow(to_RR(m),to_RR(1.5)); //m^1.5
    temp2=power((to_RR(3)*c-to_RR(2))/(to_RR(2)*c-to_RR(2)),m-1);
    temp3=to_RR(1)/c;

    k=temp1*temp2;
    k*=temp3;
    k+=to_RR(1);
    RX=pow(to_RR(N),to_RR(h-1)/to_RR(m-1))/pow(to_RR(m),1/to_RR(m-1))/pow(to_RR(2.0), to_RR(0.5));
    temp4=pow(k,to_RR(2)/to_RR(m-1));//k^(2/(n-1))
    RX/=temp4;
    RoundToZZ(X,RX);
    cout<<"rounding X = "<<X<<endl;
    cout<<"size of rounding X"<<NumBits(X)<<endl;
    cout<<"X2/X1 = "<<to_RR(X)/to_RR(Y)<<endl;
    pp[1]=X;
    cout<<"*****************\n"<<endl;

  return pp;
}




//由gij(xX)生成矩阵
//将生成的格基矩阵写入matrix_data.txt
//返回格基矩阵
mat_ZZ gen_lattice(const long h,const ZZX ff,const ZZ N,const ZZ X)
{
  ZZX f(ff);
  long delta=deg(f);
  long n=delta*h; //dimension

  ZZX g[h][delta]; //polynomial family

  ZZX fi; set(fi);

  xX(f,X); //f(x)->f(xX)

  ZZ Nh;

  Nh=power(N,h-1);

  lbound=Nh;

  mat_ZZ B; //lattice basis matrix
  B.SetDims(n,n);

  for(long i=0;i<h;i++){
    ZZX xx;
    SetCoeff(xx,0,1);

    for(long j=0;j<delta;j++){

      mul(g[i][j],Nh,fi); //N^(h-1-i)f(xX)^i
      mul(g[i][j],g[i][j],xx); //上式乘以(xX)
      LeftShift(xx,xx,1);
      mul(xx,xx,X);

      vec_ZZ u(g[i][j].rep);

      for(int k=0;k<=deg(g[i][j]);k++){
        B(i*delta+j+1)(k+1)=u[k];
      }
    }
    mul(fi,fi,f);
    Nh/=N;

  }

  ofstream fout("/home/tianyuan/Desktop/qtfile/non-qt-pro/cop_mat.txt");
  fout<<B;
  fout.close();

  return B;

}



void xX(ZZX& f,ZZ X)
{
  long degree = deg(f);
  for(int i=1;i<=degree;i++){
    f.rep[i] *= power(X,i);
  }
}



void inverse_xX(vec_ZZ& f,ZZ X)
{
  for(int i=0;i<f.length();i++){
    f[i] /= power(X,i);
  }
}


