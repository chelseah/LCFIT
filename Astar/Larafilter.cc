#include "Larafilter.h"
Larafilter::Larafilter(double *filter, int nf, int nx, int ny, int nz, double x1, double dx, double y1, double dy, double z1, double dz, double fratio,double phi,double Req,  double ggraveq, double groteq):Req_(Req){
  model_ = new LaraModelfilter(filter, nf, nx, ny, nz, x1, dx, y1, dy, z1, dz,fratio,phi,Req,ggraveq,groteq);
  quadpolar_ = new QuadPolar(*model_);
  quad2d_ = new QuadTwoD(*model_);
  pi_=3.1415926535897931;  
}

Larafilter::~Larafilter(){
  delete quadpolar_;
  delete quad2d_;
  delete model_;
}

void Larafilter::Resetphi(double phi){
  model_->Resetphi(phi);
}

double Larafilter::Cal_F0(){
  double xmin[2],xmax[2];
  double F0;
  xmin[0] = 0; xmax[0] = Req_;
  xmin[1] = 0; xmax[1] = pi_*2.; 
  model_->ResetPolarflag(1);
  F0 = quadpolar_->Integrate(xmin,2,xmax,2);
  return F0;
}

double Larafilter::Cal_Lum(){
  double xmin[2],xmax[2];
  double F0;
  //xmin[0] = Req_; xmax[0] = Req_;
  //xmin[1] = 0; xmax[1] = pi_*2.; 
  xmin[0] = -pi_/2.; xmax[0] = pi_/2.; 
  xmin[1] = 0; xmax[1] = pi_; 
  model_->ResetPolarflag(0);
  F0 = quad2d_->Integrate(xmin,2,xmax,2);
  //F0 = quadpolar_->Integrate(xmin,3,xmax,3);
  return F0;
}


