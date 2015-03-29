#include "Zeipelmodel.h"
ZeipelModel::ZeipelModel(double fratio,double phi,double Req,  double ggraveq, double groteq, double beta, double Rp):fratio_(fratio),phi_(phi),Req_(Req),ggraveq_(ggraveq),groteq_(groteq),beta_(beta),Rp_(Rp){
  //quad2d = new QuadTwoD(*this);
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

ZeipelModel::~ZeipelModel(){
  //delete quad2d;
}

void ZeipelModel::Cal_F(double *phase, int np, double *F, int nf,double theta, double a, double b){
  double *x, *y,*g,*r,*area;
  double F0,feff;
  x = new double [np];
  y = new double [np];
  g = new double [np];
  area = new double [np];
  r = new double [3*np];
  feff=Feff_();
  for (int i=0; i<np; i++){
    x[i] = phase[i]*2.*pi_*a*cos(theta)-b*sin(theta);
    y[i] = -phase[i]*2.*pi_*a*sin(theta)+b*cos(theta);
    area[i] = InterArea(x[i],y[i]);
    if (area[i] == 0){
    F[i] = 0;
    } else{
    F[i] = Int_F(x[i],y[i])/(1-feff)*(pi_*Req_*Req_)/area[i];
       // printf("%f %f %f\n",phase[i],area,Int_F(x[i],y[i]));
    }
  //}  
    //double d,z;
    ////rnew = new double [3]; //x0,y0,z0
    //  
    //d = Determinant_(x[i],y[i]);
    //if (d==-1){
    //  d =0;
    //}
    //z = Calzcoord_(x[i],y[i],d);
    //r[3*i] = x[i]; r[3*i+1] = y[i]; r[3*i+2] = z;
  }
  //Calgeff_(r,np,3,g,np);
  ////Cal_F0(&F0);
  //feff=Feff_();
  ////printf("F0=%f,feff=%f\n",F0,feff);
  //for (int i=0; i<np; i++){
  //  printf("%f %f %f %f\n",phase[i],area[i],Int_F(x[i],y[i]),pow(g[i],4.*beta_));
  //  //printf("F[i] = %f %f\n", pow(g[i],4.*beta_),Int_F(x[i],y[i]));
  //  //F[i] = pow(g[i],4.*beta_)/(1-feff)*(pi_*Req_*Req_); 
  //}
  
  delete [] x;
  delete [] y;
  delete [] g;
  delete [] r;
  return;
}


double ZeipelModel::Int_F(double x, double y){
  double xmin[2],xmax[2];
  double F;
  xmin[0] = 0; xmax[0] = Rp_;
  xmin[1] = 0; xmax[1] = pi_*2.;
  F = IntegrateShift_(xmin,2,xmax,2,x,y,false);
  return F;
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
    if(d==-1){
      d=0;
    }
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

double ZeipelModel::InterArea(double x,double y){
  //do simple integration
  double xmin[2],xmax[2];
  double area;
  xmin[0] = 0; xmax[0] = Rp_;
  xmin[1] = 0; xmax[1] = pi_*2.;
  area = IntegrateShift_(xmin,2,xmax,2,x,y,true);
  return area;

}

double ZeipelModel::fxShift(double *x, int nx,double xs,double ys,bool areaflag) {
  //x[0] = r; x[1] = theta (phi); x[2] = theta; poly coordinate; how do I do the sphere again?
  double g;
  double x0, y0;
  double *r;
  int len = 1;
  r = new double [3]; //x,y,z
  x0 = x[0]*cos(x[1])+xs;
  y0 = x[0]*sin(x[1])+ys;
  
  double d,z;
  //rnew = new double [3]; //x0,y0,z0
  d = Determinant_(x0,y0);
  if(d==-1){
    return 0;
  }
  if (!areaflag){
  z = Calzcoord_(x0,y0,d);
  r[0] = x0; r[1] = y0; r[2] = z;

  Calgeff_(r, len, 3, &g, len);
  delete [] r;
  return pow(g,4.*beta_);
  } else{
    return 1.0;
  }
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
    d=-1;
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

double ZeipelModel::IntegrateShift_(double * xmin, int nmin, double *xmax, int nmax, double x, double y, bool areaflag){
  for (int i=0; i<nmin; i++){
    xmin_[i] = xmin[i];
    xmax_[i] = xmax[i];
  }
  return Integrate1DShift_(xmin_[0],xmax_[0], x,  y, areaflag);
}

double ZeipelModel::Integrate1DShift_(double x1, double x2, double x , double y,bool areaflag){
  //x2>x1
  double xm = 0.5*(x1+x2);
  double xr = 0.5*(x2-x1);
  double s = 0;
  for (int j=0; j<5; j++){
    double dx = xr*absc_[j];
    double xm1 = xm+dx;
    double xm2 = xm-dx;
    s+=w_[j]*(Integrate2DShift_(xmin_[1],xmax_[1],xm1, x, y,areaflag) + Integrate2DShift_(xmin_[1],xmax_[1],xm2, x, y, areaflag));
  }
  return s *= xr;
}
double ZeipelModel::Integrate2DShift_(double x1, double x2, double y, double xs, double ys, bool areaflag){
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
    s+=w_[j]*(y*fxShift(xm1,2,xs,ys,areaflag) + y*fxShift(xm2,2,xs,ys,areaflag));
  }
  delete [] xm1;
  delete [] xm2;
  return s *= xr;
}
