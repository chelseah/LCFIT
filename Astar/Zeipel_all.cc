#include "Zeipel_all.h"
Zeipel::Zeipel(double fratio,double phi,double Req,  double ggraveq, double groteq, double beta):fratio_(fratio),phi_(phi),Req_(Req),ggraveq_(ggraveq),groteq_(groteq),beta_(beta){
  pi_=3.1415926535897931;  
  absc_[0] = 0.1488743389816312;
  absc_[1] = 0.4333953941292472;
  absc_[2] = 0.6794095682990244;
  absc_[3] = 0.8650633666889845;
  absc_[4] = 0.9739065285171717;
  w_ [0] = 0.2955242247147529;
  w_[1] = 0.2692667193099963;
  w_[2] = 0.2190863625159821;
  w_[3] = 0.1494513491505806;
  w_[4] = 0.0666713443086881;
}

Zeipel::~Zeipel(){
}

void Zeipel::Cal_F(double *phase, int np, double *F, int nf,double theta, double a, double b){
  double *x, *y,*g;
  double F0,feff;
  x = new double [np];
  y = new double [np];
  g = new double [np];
  for (int i=0; i<np; i++){
    x[i] = phase[i]*2.*pi_*a*cos(theta)-b*sin(theta);
    y[i] = -phase[i]*2.*pi_*a*sin(theta)+b*cos(theta);
  }
  Calgeff_(x,np,y,np,g,np);
  Cal_F0(&F0,1);
  feff=Feff_();
  //printf("F0=%f,feff=%f\n",F0,feff);
  for (int i=0; i<np; i++){
    F[i] = pow(g[i],4.*beta_)/(F0/pi_/Req_/Req_*(1-feff)); 
  }
  delete [] x;
  delete [] y;
  delete [] g;
  return;
}

void Zeipel::Cal_F0(double *F0,int np){
  double x1,x2,y1,y2;
  x1 = 0; x2 = Req_;
  y1 = 0; y2 = pi_*2.; 
  *F0 = Integrate(x1,x2,y1,y2);
  return;
}


double Zeipel::fx(double *x,int nx) {
  //x[0] = r; x[1] = theta; poly coordinate
  double g;
  double x0, y0;
  int len = 1;
  x0 = x[0]*cos(x[1]);
  y0 = x[0]*sin(x[1]);
  Calgeff_(&x0, len , &y0, len, &g, len);
  //printf("%f %f\n",x[0],g);
  return x[0]*pow(g,4.*beta_);
}

double Zeipel::Determinant_(double x,double y){
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

double Zeipel::Calzcoord_(double x, double y, double d){
  double za, zb;
  double f2 = (1-fratio_)*(1-fratio_);
  za = -2*y*(1-f2)*sin(phi_)*cos(phi_)+sqrt(d);
  zb = 2*(f2*cos(phi_)*cos(phi_)+sin(phi_)*sin(phi_));
  return za/zb;
}

void Zeipel::Rotate_(double *x, double *xnew){
  xnew[0] = x[0];
  xnew[1] = x[1]*cos(phi_)+x[2]*sin(phi_);
  xnew[2] = -x[1]*sin(phi_) + x[2]*cos(phi_);
  return;
}

void Zeipel::Calgeff_(double *x, int nx, double *y, int ny, double *g, int ng){
  for (int i = 0; i<nx; i++){
    double d,z,absR,absRper,gi,gj,gz,ggrave,grote;
    double *r,*rnew;
    r = new double [3]; //x,y,z
    rnew = new double [3]; //x0,y0,z0
    d = Determinant_(x[i],y[i]);
    z = Calzcoord_(x[i],y[i],d);
    r[0] = x[i]; r[1] = y[i]; r[2] = z;
    Rotate_(r,rnew);
    absR = sqrt(rnew[0]*rnew[0]+rnew[1]*rnew[1]+rnew[2]*rnew[2]);
    absRper = sqrt(rnew[0]*rnew[0] + rnew[2]*rnew[2]);
    ggrave = -ggraveq_*(Req_/absR)*(Req_/absR);
    grote = groteq_/(Req_/absRper);
    gi = ggrave * rnew[0]/absR + grote*rnew[0]/absRper;
    gj = ggrave * rnew[1]/absR;
    gz = ggrave * rnew[2]/absR + grote*rnew[2]/absRper;
    g[i] = sqrt(gi*gi+gj*gj+gz*gz);
    delete [] r;
    delete [] rnew;
  }
}

double Zeipel::Feff_(){
  double feff, f2=(1-fratio_)*(1-fratio_);
  feff = 1- sqrt(f2*cos(phi_)*cos(phi_)+sin(phi_)*sin(phi_));
  return feff;
}

double Zeipel::Integrate(double x1, double x2, double y1,double y2){
  x1_ = x1; x2_=x2;
  y1_ = y1; y2_=y2;
  return Integrate1D_(x1,x2);
}

double Zeipel::Integrate1D_(double x1, double x2){
  //x2>x1
  double xm = 0.5*(x1+x2);
  double xr = 0.5*(x2-x1);
  double s = 0;
  for (int j=0; j<5; j++){
    double dx = xr*absc_[j];
    double xm1 = xm+dx;
    double xm2 = xm-dx;
    s+=w_[j]*(Integrate2D_(y1_,y2_,xm1) + Integrate2D_(y1_,y2_,xm2));
  }
  return s *= xr;
}
double Zeipel::Integrate2D_(double x1, double x2, double y){
  //x2>x1, here x1 is y1, x2 is y2
  double xm = 0.5*(x1+x2);
  double xr = 0.5*(x2-x1);
  double s = 0;
  double *xm1 = new double [2];
  double *xm2 = new double [2];
  xm1[0] = y;
  xm2[0] = y;
  for (int j=0; j<5; j++){
    double dx = xr*absc_[j];
    xm1[1] = xm+dx;
    xm2[1] = xm-dx;
    s+=w_[j]*(fx(xm1,2) + fx(xm2,2));
  }
  delete [] xm1;
  delete [] xm2;
  return s *= xr;
}
