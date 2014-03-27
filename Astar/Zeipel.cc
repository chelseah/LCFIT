#include "Zeipel.h"
Zeipel::Zeipel(double fratio,double phi,double Req,  double ggraveq, double groteq, double beta):Req_(Req){
  model_ = new ZeipelModel(fratio,phi,Req,ggraveq,groteq,beta);
  quadpolar_ = new QuadPolar(*model_);
  quadsphere_ = new QuadSphere(*model_);
  pi_=3.1415926535897931;  
}

Zeipel::~Zeipel(){
  delete quadpolar_;
  delete quadsphere_;
  delete model_;
}

void Zeipel::Cal_F(double *phase, int np, double *F, int nf,double theta, double a, double b){
  double F0;
  model_->Cal_F(phase,np,F,nf,theta,a,b); 
  Cal_F0(&F0,1);
  //feff=Feff_();
  //printf("F0=%f,feff=%f\n",F0,feff);
  for (int i=0; i<np; i++){
    F[i] /= F0; 
  }
  return;
}


void Zeipel::Cal_F0(double *F0,int np){
  double xmin[2],xmax[2];
  xmin[0] = 0; xmax[0] = Req_;
  xmin[1] = 0; xmax[1] = pi_*2.; 
  *F0 = quadpolar_->Integrate(xmin,2,xmax,2);
  return;
}

double Zeipel::Cal_Lum(){
  double xmin[3],xmax[3];
  double F0;
  xmin[0] = 0; xmax[0] = Req_;
  xmin[1] = 0; xmax[1] = pi_*2.; 
  xmin[2] = 0; xmax[2] = pi_; 
  F0 = quadsphere_->Integrate(xmin,3,xmax,3);
  return F0;
}


