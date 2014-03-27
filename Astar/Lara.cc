#include "Lara.h"
Lara::Lara(double fratio,double phi,double Req,  double ggraveq, double groteq):Req_(Req){
  model_ = new LaraModel(fratio,phi,Req,ggraveq,groteq);
  quadpolar_ = new QuadPolar(*model_);
  quadsphere_ = new QuadSphere(*model_);
  pi_=3.1415926535897931;  
}

Lara::~Lara(){
  delete quadpolar_;
  delete quadsphere_;
  delete model_;
}

void Lara::CalTtheta(double *theta, int nth, double *T, int nT, double absR){
  model_->CalTtheta(theta,nth,T,nT,absR);
}
void Lara::Cal_Txy(double *x, int nx, double *y, int ny, double *Txy, int nT){
  model_-> Cal_Txy(x, nx,y,ny,Txy,nT);
}

void Lara::Cal_F(double *phase, int np, double *F, int nf,double theta, double a, double b){
  double F0;
  model_->Cal_F(phase,np,F,nf,theta,a,b); 
  Cal_F0(&F0,1);
  F0/= (pi_*Req_*Req_);
  //feff=Feff_();
  //printf("F0=%f,feff=%f\n",F0,feff);
  for (int i=0; i<np; i++){
    F[i] /= F0; 
  }
  return;
}


void Lara::Cal_F0(double *F0,int np){
  double xmin[2],xmax[2];
  xmin[0] = 0; xmax[0] = Req_;
  xmin[1] = 0; xmax[1] = pi_*2.; 
  *F0 = quadpolar_->Integrate(xmin,2,xmax,2);
  return;
}

double Lara::Cal_Lum(){
  double xmin[3],xmax[3];
  double F0;
  xmin[0] = 0; xmax[0] = Req_;
  xmin[1] = 0; xmax[1] = pi_*2.; 
  xmin[2] = 0; xmax[2] = pi_; 
  F0 = quadsphere_->Integrate(xmin,3,xmax,3);
  return F0;
}


