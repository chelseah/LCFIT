#include "Zeipelmodel.h"
ZeipelModel::ZeipelModel(double fratio,double phi,double Req,  double ggraveq, double groteq, double beta):fratio_(fratio),phi_(phi),Req_(Req),ggraveq_(ggraveq),groteq_(groteq),beta_(beta){
  //quad2d = new QuadTwoD(*this);
  pi_=3.1415926535897931;  
}

ZeipelModel::~ZeipelModel(){
  //delete quad2d;
}

void ZeipelModel::Cal_F(double *phase, int np, double *F, int nf,double theta, double a, double b){
  double *x, *y,*g,*r;
  double F0,feff;
  x = new double [np];
  y = new double [np];
  g = new double [np];
  r = new double [3*np];
  for (int i=0; i<np; i++){
    x[i] = phase[i]*2.*pi_*a*cos(theta)-b*sin(theta);
    y[i] = -phase[i]*2.*pi_*a*sin(theta)+b*cos(theta);
    double d,z;
    //rnew = new double [3]; //x0,y0,z0
    d = Determinant_(x[i],y[i]);
    z = Calzcoord_(x[i],y[i],d);
    r[3*i] = x[i]; r[3*i+1] = y[i]; r[3*i+2] = z;
  }
  Calgeff_(r,np,3,g,np);
  //Cal_F0(&F0);
  feff=Feff_();
  //printf("F0=%f,feff=%f\n",F0,feff);
  for (int i=0; i<np; i++){
    F[i] = pow(g[i],4.*beta_)/((1-feff)/pi_/Req_/Req_); 
  }
  delete [] x;
  delete [] y;
  delete [] g;
  delete [] r;
  return;
}


double ZeipelModel::fx(double *x, int nx) {
  //x[0] = r; x[1] = theta (phi); x[2] = theta; poly coordinate; how do I do the sphere again?
  double g;
  double x0, y0;
  double *r;
  int len = 1;
  r = new double [3]; //x,y,z
  if(nx ==2){
    x0 = x[0]*cos(x[1]);
    y0 = x[0]*sin(x[1]);
    double d,z;
    //rnew = new double [3]; //x0,y0,z0
    d = Determinant_(x0,y0);
    z = Calzcoord_(x0,y0,d);
    r[0] = x0; r[1] = y0; r[2] = z;

    Calgeff_(r, len, 3, &g, len);
  } else {
    if(nx==3){
      double z0;
      x0 = x[0]*cos(x[1])*sin(x[2]);
      y0 = x[0]*sin(x[1])*sin(x[2]);
      z0 = x[0]*cos(x[2]);
      r[0] = x0; r[1] = y0; r[2] = z0; 
      Calgeff_(r, len, 3, &g, len);
    }
  }
  //printf("%f %f\n",x[0],g);
  //return x[0]*pow(g,4.*beta_);
  delete [] r;
  return pow(g,4.*beta_);
}

double ZeipelModel::Determinant_(double x,double y){
  double da, db, dc;
  double f2 = (1-fratio_)*(1-fratio_);
  double sinphi2 = sin(phi_)*sin(phi_);
  double cosphi2 = cos(phi_)*cos(phi_);
  double d;
  da = 4.*y*y * (1-f2)*(1-f2) * sinphi2 * cosphi2;
  db = cosphi2 * f2 + sinphi2;
  dc = (y*y * sinphi2 - Req_*Req_ + x*x) * f2 + y*y * cosphi2;
  d = da - 4*db*dc;
  if(d<0){
    d=0;
  }
  return d;
}

double ZeipelModel::Calzcoord_(double x, double y, double d){
  double za, zb;
  double f2 = (1-fratio_)*(1-fratio_);
  za = -2*y*(1-f2)*sin(phi_)*cos(phi_)+sqrt(d);
  zb = 2*(f2*cos(phi_)*cos(phi_)+sin(phi_)*sin(phi_));
  return za/zb;
}

void ZeipelModel::Rotate_(double *x, double *xnew){
  xnew[0] = x[0];
  xnew[1] = x[1]*cos(phi_)+x[2]*sin(phi_);
  xnew[2] = -x[1]*sin(phi_) + x[2]*cos(phi_);
  return;
}

void ZeipelModel::Calgeff_(double *r, int np, int nd, double *g, int ng){
  for (int i = 0; i<np; i++){
    double absR,absRper,gi,gj,gz,ggrave,grote;
    double *rtemp,*rnew;
    rtemp = new double [3]; //x,y,z
    rnew = new double [3]; //x0,y0,z0
    for (int j=0; j<nd; j++){
      rtemp[j] =r[nd*i+j]; 
    }
    Rotate_(rtemp,rnew);
    absR = sqrt(rnew[0]*rnew[0]+rnew[1]*rnew[1]+rnew[2]*rnew[2]);
    absRper = sqrt(rnew[0]*rnew[0] + rnew[2]*rnew[2]);
    ggrave = -ggraveq_*(Req_/absR)*(Req_/absR);
    grote = groteq_/(Req_/absRper);
    gi = ggrave * rnew[0]/absR + grote*rnew[0]/absRper;
    gj = ggrave * rnew[1]/absR;
    gz = ggrave * rnew[2]/absR + grote*rnew[2]/absRper;
    g[i] = sqrt(gi*gi+gj*gj+gz*gz);
    delete [] rtemp;
    delete [] rnew;
  }
}

double ZeipelModel::Feff_(){
  double feff, f2=(1-fratio_)*(1-fratio_);
  feff = 1- sqrt(f2*cos(phi_)*cos(phi_)+sin(phi_)*sin(phi_));
  return feff;
}
