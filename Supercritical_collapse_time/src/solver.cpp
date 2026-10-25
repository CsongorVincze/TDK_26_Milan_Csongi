// Radial focusing Klein-Gordon: identical Module 4 stencil and RK4.
// Extend the domain by adding grid nodes; preserve every original inner node.
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#ifdef EXTENDED_PRECISION
using Real=long double;
#else
using Real=double;
#endif
#ifndef SPATIAL_ORDER
#define SPATIAL_ORDER 2
#endif
using V=std::vector<Real>;
int main(int argc,char**argv){
 if(argc!=11){std::cerr<<"A base_N CFL horizon cutoff radius sponge_start strength sample_interval trajectory\n";return 2;}
 Real A=std::stold(argv[1]),cfl=std::stold(argv[3]),end=std::stold(argv[4]),level=std::stold(argv[5]);
 int base=std::stoi(argv[2]);Real radius=std::stold(argv[6]),start=std::stold(argv[7]),strength=std::stold(argv[8]);
 // argv[9]: sample interval; argv[10]: path or '-'.
 if(base<3 || cfl<=0 || end<=0 || level<=A || radius<=start || start<=0)return 2;
 Real dx=Real(1)/(base-1),scale=25/std::sinh(Real(5));
 int n=std::ceil(std::asinh(radius/scale)/5/dx-Real(1e-12))+1;
 if(radius==25)n=base;
 Real actual_radius=scale*std::sinh(5*dx*(n-1));
 Real dt0=cfl*5*scale*dx;
 V p(n),v(n),a(n),b(n),gamma(n),r(n),q(n),w(n),kp[4],kv[4];
 for(int k=0;k<4;k++){kp[k].resize(n);kv[k].resize(n);}
 for(int i=0;i<n;i++){
  Real x=i*dx,rx=5*scale*std::cosh(5*x);r[i]=scale*std::sinh(5*x);
  p[i]=A*std::exp(-r[i]*r[i]/8);a[i]=1/(rx*rx*dx*dx);
  b[i]=i?(2/(r[i]*rx)-25*r[i]/(rx*rx*rx))/(2*dx):0;
  gamma[i]=strength*std::pow(std::max((r[i]-start)/(actual_radius-start),Real(0)),3);
 }
 p.back()=0;
 auto rhs=[&](const V&f,const V&z,V&df,V&dz){
  Real origin=6*a[0]*(f[1]-f[0]);
  if(SPATIAL_ORDER==4)origin=a[0]*(-f[2]+16*f[1]-15*f[0])/2;
  dz[0]=origin-(f[0]-f[0]*f[0]*f[0])-gamma[0]*z[0];df[0]=z[0];
  for(int i=1;i<n-1;i++){
   Real lap=a[i]*(f[i+1]-2*f[i]+f[i-1])+b[i]*(f[i+1]-f[i-1]);
   if(SPATIAL_ORDER==4 && i<n-2){
    Real fm2=f[i==1?1:i-2];
    lap=a[i]*(-f[i+2]+16*f[i+1]-30*f[i]+16*f[i-1]-fm2)/12
        +b[i]*(-f[i+2]+8*f[i+1]-8*f[i-1]+fm2)/6;
   }
   df[i]=z[i];dz[i]=lap-(f[i]-f[i]*f[i]*f[i])-gamma[i]*z[i];
  }df.back()=dz.back()=0;
 };
 std::ofstream out;bool record=std::string(argv[10])!="-";
 Real spacing=std::stold(argv[9]),next=spacing,t=0,peakcenter=std::abs(p[0]);long steps=0;
 auto sample=[&](Real ts,Real f){
  Real center=p[0]+f*(q[0]-p[0]),velocity=v[0]+f*(w[0]-v[0]),half=std::abs(center)/2;
  Real r50=std::numeric_limits<Real>::quiet_NaN(),prev=std::abs(center);
  for(int i=1;i<n;i++){Real value=std::abs(p[i]+f*(q[i]-p[i]));if(value<=half){r50=r[i-1]+(r[i]-r[i-1])*(prev-half)/(prev-value);break;}prev=value;}
  out<<ts<<","<<center<<","<<velocity<<","<<r50<<"\n";
 };
 if(record){out.open(argv[10]);if(!out)return 4;out<<std::setprecision(std::numeric_limits<Real>::max_digits10)<<"time,phi_center,pi_center,r50\n";q=p;w=v;sample(0,0);}
 bool crossed=false;Real t12=-1,t24=-1;int zero_crossings=0;
 while(t<end){
  Real peak=1;for(Real z:p)peak=std::max(peak,std::abs(z));
  Real dt=std::min({dt0,Real(.15)/peak,end-t});
  rhs(p,v,kp[0],kv[0]);
  for(int k=1;k<4;k++){
   Real f=k==3?1:Real(.5);
   for(int i=0;i<n;i++){q[i]=p[i]+f*dt*kp[k-1][i];w[i]=v[i]+f*dt*kv[k-1][i];}
   rhs(q,w,kp[k],kv[k]);
  }
  for(int i=0;i<n;i++){
   q[i]=p[i]+dt*(kp[0][i]+2*kp[1][i]+2*kp[2][i]+kp[3][i])/6;
   w[i]=v[i]+dt*(kv[0][i]+2*kv[1][i]+2*kv[2][i]+kv[3][i])/6;
   if(!std::isfinite(q[i])||!std::isfinite(w[i]))return 3;
  }
  Real oldabs=std::abs(p[0]),newabs=std::abs(q[0]);
  if((p[0]>0 && q[0]<0)||(p[0]<0 && q[0]>0))zero_crossings++;
  if(t12<0 && newabs>=12 && oldabs<12)t12=t+dt*(12-oldabs)/(newabs-oldabs);
  if(t24<0 && newabs>=24 && oldabs<24)t24=t+dt*(24-oldabs)/(newabs-oldabs);
  peakcenter=std::max(peakcenter,newabs);
  crossed=newabs>=level;
  Real f=crossed?(level-oldabs)/(newabs-oldabs):1,last=t+dt*f;
  if(record){while(next<last-Real(1e-12)){sample(next,(next-t)/dt);next+=spacing;}}
  steps++;
  if(crossed){if(record)sample(last,f);t=last;break;}
  t=last;p.swap(q);v.swap(w);
 }
 if(record && !crossed){q=p;w=v;sample(t,0);}
 Real finalcenter=crossed?std::copysign(level,q[0]):p[0];
 std::cout<<std::setprecision(std::numeric_limits<Real>::max_digits10)<<crossed<<" "<<t<<" "<<steps<<" "<<(crossed?(finalcenter>0?1:-1):0)<<" "<<peakcenter<<" "<<finalcenter<<" "<<actual_radius<<" "<<t12<<" "<<t24<<" "<<zero_crossings<<"\n";
}
