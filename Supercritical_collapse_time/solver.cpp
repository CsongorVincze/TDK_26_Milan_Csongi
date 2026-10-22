// Module 4 finite-difference/RK4 solver; output interpolation does not alter steps.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vector>
using V=std::vector<double>;
int main(int argc,char**argv){
 if(argc!=7 && argc!=8)return 2;
 double A=std::stod(argv[1]),cfl=std::stod(argv[3]),end=std::stod(argv[4]);
 int n=std::stoi(argv[2]); double dx=1./(n-1),dt0=cfl*25*5/std::sinh(5.)*dx;
 V p(n),v(n),a(n),b(n),gam(n),q(n),w(n),kp[4],kv[4];
 for(int k=0;k<4;k++){kp[k].resize(n);kv[k].resize(n);}
 for(int i=0;i<n;i++){
  double x=i*dx,r=25*std::sinh(5*x)/std::sinh(5.),rx=125*std::cosh(5*x)/std::sinh(5.);
  p[i]=A*std::exp(-r*r/8);a[i]=1/(rx*rx*dx*dx);
  b[i]= i ? (2/(r*rx)-25*r/(rx*rx*rx))/(2*dx):0;
  gam[i]=7*std::pow(std::max((r-20)/5,0.),3);
 }
 p.back()=0;
 auto rhs=[&](const V&f,const V&z,V&df,V&dz){
  for(int i=0;i<n-1;i++){
   double lap=i ? a[i]*(f[i+1]-2*f[i]+f[i-1])+b[i]*(f[i+1]-f[i-1]):6*a[0]*(f[1]-f[0]);
   df[i]=z[i];dz[i]=lap-(f[i]-f[i]*f[i]*f[i])-gam[i]*z[i];
  }df.back()=dz.back()=0;
 };
 std::ofstream out;bool record=std::string(argv[6])!="-";
 double spacing=std::stod(argv[5]),next=spacing,t=0;long steps=0;
 if(record){out.open(argv[6]);out<<std::setprecision(17)<<"time,phi_center\n0,"<<A<<"\n";}
 std::ofstream fields;double next_field=.2;
 if(argc==8){fields.open(argv[7]);fields<<std::setprecision(17)<<"time";for(int i=0;i<n;i++)fields<<",phi_"<<i;fields<<"\n0";for(double z:p)fields<<","<<z;fields<<"\n";}
 bool crossed=false;
 while(t<end){
  double peak=1;for(double z:p)peak=std::max(peak,std::abs(z));
  double dt=std::min({dt0,.15/peak,end-t});
  rhs(p,v,kp[0],kv[0]);
  for(int k=1;k<4;k++){
   double f=k==3?1.:.5;
   for(int i=0;i<n;i++){q[i]=p[i]+f*dt*kp[k-1][i];w[i]=v[i]+f*dt*kv[k-1][i];}
   rhs(q,w,kp[k],kv[k]);
  }
  for(int i=0;i<n;i++){
   q[i]=p[i]+dt*(kp[0][i]+2*kp[1][i]+2*kp[2][i]+kp[3][i])/6;
   w[i]=v[i]+dt*(kv[0][i]+2*kv[1][i]+2*kv[2][i]+kv[3][i])/6;
   if(!std::isfinite(q[i])||!std::isfinite(w[i]))return 3;
  }
  crossed=std::abs(q[0])>=12;
  double f=crossed?(12-std::abs(p[0]))/(std::abs(q[0])-std::abs(p[0])):1.;
  double last=t+dt*f;
  if(record){while(next<=last){out<<next<<","<<p[0]+(next-t)/dt*(q[0]-p[0])<<"\n";next+=spacing;}}
  if(fields.is_open()){
   while(next_field<=last){fields<<next_field;for(int i=0;i<n;i++)fields<<","<<p[i]+(next_field-t)/dt*(q[i]-p[i]);fields<<"\n";next_field+=.2;}
   if(crossed){fields<<last;for(int i=0;i<n;i++)fields<<","<<p[i]+f*(q[i]-p[i]);fields<<"\n";}
  }
  steps++;t=last;
  if(crossed){if(record)out<<t<<","<<std::copysign(12.,q[0])<<"\n";break;}
  p.swap(q);v.swap(w);
 }
 if(fields.is_open() && !crossed){fields<<t;for(double z:p)fields<<","<<z;fields<<"\n";}
 std::cout<<std::setprecision(17)<<crossed<<" "<<t<<" "<<steps<<"\n";
}
