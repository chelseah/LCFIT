#include "Laramodelfilter.h"
LaraModelfilter::LaraModelfilter(double *filter, int nf, int nx, int ny, int nz, double x1, double dx, double y1, double dy, double z1, double dz, double fratio,double phi,double Req,  double ggraveq, double groteq):fratio_(fratio),phi_(phi),Req_(Req),ggraveq_(ggraveq),groteq_(groteq){
  //quad2d = new QuadTwoD(*this);
  filter_ = new double [nf];
  nx_ = nx;
  ny_=ny;
  nz_=nz;
  nf_ = nf;
  x1_ = x1;
  y1_=y1;
  z1_=z1;
  dx_ = dx;
  dy_ = dy;
  dz_ = dz;
  for (int k=0;k<nz_;k++){
    for(int j=0;j<ny_; j++){
      for (int i=0; i<nx_; i++){
        filter_[i+j*nx_+k*(nx_*ny_)] = filter[i+j*nx_+k*(nx_*ny_)];
      }
    }
  }
  //printf("in c, %f %f %f", filter_[0],filter_[3*(ny_*nz_)+4*nz_+5],filter_[11*(ny_*nz_)+12*nz_+13]);  
  pi_=3.1415926535897931;
  tol_ = 1.e-5;
  Maxiter_ = 50000;
  polarflag_=0;
  w_ = sqrt(groteq_/ggraveq);
  //printf("w_=%f\n",w_);
  //double tempf = Readfilter_(4.6, log10(1.0), 0.0);
  //printf("tempf=%f\n",tempf);
}

LaraModelfilter::~LaraModelfilter(){
  delete filter_;
}

double LaraModelfilter::Readfilter_(double g, double T, double phi){
  int indx_g, indx_T, indx_phi;
  indx_g = (int)((g-x1_)/dx_+0.5);
  indx_T = (int)((T-y1_)/dy_+0.5);
  indx_phi = (int)((phi-z1_)/dz_+0.5);
  if(indx_g>=nx){
    indx_g=nx-1;
    printf("Warning: g index out of bound!\n")
  }
  if(indx_g<0){
    indx_g=0;
    printf("Warning: g index out of bound!\n")
  }
  //printf("%f %f %f %f %f %d %d %d, %d %d %d\n",g, T, phi,z1_, dz_,indx_g,indx_T,indx_phi, nx_, ny_, nz_);
  //printf("%f %f %f %f %d %d %d, %d %d %d %f\n",g, T, y1_, dy_,indx_g,indx_T,indx_phi, nx_, ny_, nz_,filter_[indx_phi+indx_T*nz_+indx_g*ny_*nz_]);
  return filter_[indx_phi+indx_T*nz_+indx_g*ny_*nz_];
}

void LaraModelfilter::Resetphi(double phi){
  phi_=phi;
}

void LaraModelfilter::ResetPolarflag(int flag){
  polarflag_ = flag;
}

double LaraModelfilter::fx(double *x, int nx) {
  //x[0] = r; x[1] = theta (phi); x[2] = theta; poly coordinate; how do I do the sphere again?
  double T;
  double x0, y0;
  double *r;
  double mu;
  int len = 1;
  r = new double [3]; //x,y,z
  if(polarflag_ ==1){
    //transform to cartitian on skyplane
    x0 = x[0]*cos(x[1]);
    y0 = x[0]*sin(x[1]);
    double d,z;
    //rnew = new double [3]; //x0,y0,z0
    d = Determinant_(x0,y0);
    z = Calzcoord_(x0,y0,d);
    r[0] = x0/Req_; r[1] = y0/Req_; r[2] = z/Req_;
    mu = r[2]/sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); 
    CalTeff_(r, len, 3, &T, len);
    double g = ggraveq_/(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); 
    delete [] r;
    return pow(10,Readfilter_(log10(g),log10(T),mu));

  } else {
    if(polarflag_==0){
      double z0,absR;
      //x[0] is phi, x[1] is theta in spherical coordinate
      //theta
      z0 = Req_*cos(x[0])*sin(x[1]);
      x0 = Req_*sin(x[0])*sin(x[1]);
      y0 = Req_*cos(x[1]);
      r[0] = x0/Req_; r[1] = y0/Req_; r[2] = z0/Req_; 
      absR=sqrt(z0*z0+x0*x0+y0*y0);
      mu = z0/Req_;
      //printf("z0=%f,x0= %f,y0= %f, absR = %f, Req=%f\n",x0,y0,z0,absR,Req_);
      CalTeff_(r, len, 3, &T, len);
      //return x[0]*pow(g,4.*beta_);
    double g = ggraveq_/(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]); 
    delete [] r;
    //printf("%f %f\n",x[0],g);
    //filter is in grid of logg, logT and viewing angel mu,mu is 0 when look 
    //through the surface, 1 when at the edge, the returned value 
    //is in log10(photons/s/cm^2/sr)
    return Req_*Req_*sin(x[1])*pow(10,Readfilter_(log10(g),log10(T),mu));

    }
  }
  }

double LaraModelfilter::Determinant_(double x,double y){
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

double LaraModelfilter::Calzcoord_(double x, double y, double d){
  double za, zb;
  double f2 = (1-fratio_)*(1-fratio_);
  za = -2*y*(1-f2)*sin(phi_)*cos(phi_)+sqrt(d);
  zb = 2*(f2*cos(phi_)*cos(phi_)+sin(phi_)*sin(phi_));
  return za/zb;
}

void LaraModelfilter::Rotate_(double *x, double *xnew){
  xnew[0] = x[0];
  xnew[1] = x[1]*cos(phi_)+x[2]*sin(phi_);
  xnew[2] = -x[1]*sin(phi_) + x[2]*cos(phi_);
  return;
}

double LaraModelfilter::Caltau_(double theta, double r){
  double tau;
  tau = 1./3.*w_*w_*pow(r,3.)*pow(cos(theta),3.)+cos(theta)+log(tan(theta/2.));
  return tau;
}

double LaraModelfilter::Solveangle_(double tau,double theta){
  double x = 0., x0 = tan(theta/2.),h;
  /*newton method*/
  double fx,dfx;
  for (int i=0; i < Maxiter_; i++){
    x = x0;
    //printf("x=%f\n",x);
    fx = AngleFx_(x,tau);
    if(fabs(fx)<1.e-5){
      return x;
    }
    dfx = AngleDfx_(x);
    h = fx/dfx;
    while((x-h)<0){
      h/=2.;
    }
    x-=h;
  }
  printf("in angle,exceed maxiter,fx=%f,x=%f,theta=%f\n",fx,x,theta/pi_);
  return x;
  /*bisec method*/
  //double x1=0.01, x2 = 0.99*pi_/2.,dx;
  //double f = AngleFx_(tan(x1/2.),tau);
  //double fmid = AngleFx_(tan(x2/2.),tau);
  //if (f*fmid >=0.0) {
  //  printf("Root not bracketed,f=%f,fmid=%f\n",f,fmid);
  //  throw("Root must be bracketed fro bisection in rtbis");
  //}
  //double rtb;
  //rtb = f<0.0?(dx=x2-x1,x1):(dx=x1-x2,x2);
  //for (int j=0; j< Maxiter_; j++){
  //  double xmid = rtb+dx*0.5;
  //  dx*=0.5;
  //  fmid = AngleFx_(tan(xmid/2.),tau);
  //  if(fmid<=0.) rtb = xmid;
  //  if(fabs(dx) < tol_ || fmid ==0) return rtb;
  //}
  //throw("Too many bisections in rtbis");
}
double LaraModelfilter::AngleFx_(double x,double tau){
  double fx;
  fx = (1.-x*x)/(1.+x*x)+log(x)-tau;   
  return fx;
}
double LaraModelfilter::AngleDfx_(double x){
  double dfx; 
  dfx = -4.*x/pow((1.+x*x),2.)+1./x;
  return dfx;
}

double LaraModelfilter::Solveradius_(double theta){
  double x = 0., x0 = 1.0,h;
  double fx,dfx;
  /*newton method*/
  for (int i=0; i < Maxiter_; i++){
    x = x0;
    //printf("x=%f\n",x);
    fx = RFx_(x,sin(theta));
    if(fabs(fx)<tol_){
      return x;
    }
    dfx = RDfx_(x,sin(theta));
    h = fx/dfx;
    while((x-h)<0){
      h/=2.;
    }
    x-=h;
  }
  printf("in radius,exceed maxiter,fx=%f,x=%f,theta=%f\n",fx,x,theta/pi_);
  return x;
} 


double LaraModelfilter::RFx_(double x, double sintheta){
  double fx;
  fx = 0.5*pow(x,3.)*sintheta*sintheta - (1./w_/w_+0.5)*x + 1./w_/w_;
  //printf("insolver,%f %f %f\n",fx,x,sintheta);
  return fx;
}

double LaraModelfilter::RDfx_(double x, double sintheta){
  double dfx;
  dfx = 1.5*pow(x,2.)*sintheta*sintheta - (1./w_/w_+0.5);
  return dfx;
}

void LaraModelfilter::CalTtheta(double *theta, int nth, double *T, int nT, double absR){
//double Tpole = pow((1/pow(1-fratio_,4.)+pow(w_,4.)*pow(1-fratio_,2.)-2*w_*w_/(1-fratio_)),1./8.)*pow((1.-w_*w_*pow((1-fratio_),3.)),-1./6.);
  double Tpole = pow(1/(1-fratio_),1./2.)*pow(exp(2./3.*w_*w_*pow((1-fratio_),3)),0.25);
  double tau,Thalf, tanTheta, Ta,Tb;
  double inversR;
  for (int i=0; i<nth; i++){
    
    //absR/=(2./(1-fratio_));
  if(fabs(theta[i])<1.e-3){
      //T[i] = pow(exp(2./3.*w_*w_*absR*absR*absR),0.25)/Tpole;
      T[i] = 1.0;
      //printf("%f %f\n",T[i],w_);
      continue;
    }
    if(fabs(theta[i]-pi_/2.)<1.e-2)
    {
      T[i] = sqrt(2./(2.+w_*w_))*pow((1-w_*w_),1./12.)*exp(-4./3.*w_*w_/pow((2+w_*w_),3));
      continue;
    }
  
    if(cos(theta[i])<1.e-2){
      //T[i] = sqrt(2./(2.+w_*w_))*pow((1-w_*w_),1./12.)*exp(-4./3.*w_*w_/pow((2+w_*w_),3));
      Ta = 1./pow((1-fratio_),4.) + pow(w_,4.)*(1-fratio_)*(1-fratio_)*sin(theta[i])*sin(theta[i]);
      Tb = 2.*w_*w_*sin(theta[i])*sin(theta[i])/(1-fratio_);
      T[i] = pow((Ta-Tb),1./8.) * pow((1-w_*w_*pow((1-fratio_),3)),1./6.)/Tpole;
      continue;
    }
    inversR = pow(cos(theta[i]),2.)/pow((1-fratio_),2)+pow(sin(theta[i]),2.);
    if(w_==0){
      absR=1;
    } else{
      absR = Solveradius_(theta[i]);
    }
    printf("%f %f %f\n", sqrt(1./inversR),absR,theta[i]/pi_);
        //printf("theta=%f\n",theta);
    tau = Caltau_(theta[i],absR);
    //printf("tau=%f\n",tau);
  
    Thalf = Solveangle_(tau,theta[i]);
    //printf("Thalf=%f\n",Thalf);
    tanTheta = (2*Thalf)/(1-Thalf*Thalf);
    //printf("tanTheta=%f,tan(theta)=%f\n",tanTheta,tan(theta));
    Ta = 1./pow(absR,4.) + pow(w_,4.)*absR*absR*sin(theta[i])*sin(theta[i]);
    Tb = 2.*w_*w_*sin(theta[i])*sin(theta[i])/absR;
    T[i] = pow((Ta-Tb),1./8.) * sqrt(fabs(tanTheta/tan(theta[i])))/Tpole;
  
  } 
}
void LaraModelfilter::CalTeff_(double *r, int np, int nd, double *T, int ng){
  
 //double Tpole = pow((1/pow(1-fratio_,4.)+pow(w_,4.)*pow(1-fratio_,2.)-2*w_*w_/(1-fratio_)),1./8.)*pow((1.-w_*w_*pow((1-fratio_),3.)),-1./6.);
  double Tpole = pow(1/pow(1-fratio_,4.),1./8.)*pow(exp(2./3.*w_*w_*pow((1-fratio_),3)),0.25);
  for (int i = 0; i<np; i++){
    double absR,absRper,theta;
    double *rtemp,*rnew;
    if(T[i] == -1){
      T[i]=0.;
      continue;
    }
    rtemp = new double [3]; //x,y,z
    rnew = new double [3]; //x0,y0,z0
    for (int j=0; j<nd; j++){
      rtemp[j] =r[nd*i+j]; 
    }
    Rotate_(rtemp,rnew);
    //double absRtemp = sqrt(rtemp[0]*rtemp[0]+rtemp[1]*rtemp[1]+rtemp[2]*rtemp[2]);
    //printf("rnew=%f %f %f\n",rnew[0],rnew[1],rnew[2]);
    absR = sqrt(rnew[0]*rnew[0]+rnew[1]*rnew[1]+rnew[2]*rnew[2]);
    //absRper = sqrt(rnew[0]*rnew[0] + rnew[2]*rnew[2]);
    //printf("absR=%f,absRper=%f\n",absR,absRtemp);
    theta = acos(rnew[1]/absR);
    //printf("theta=%f\n",theta);
    if (theta<0){
      theta=fabs(theta);
    }
    
    if (theta>pi_/2.){
      theta = pi_-theta;
    }
    if(theta>pi_/2. || theta <0){
      //printf("theta=%f\n",theta);
    }

   double tau,Thalf, tanTheta, Ta,Tb;
  if(fabs(theta)<1.e-3){
      //T[i] = pow(exp(2./3.*w_*w_*absR*absR*absR),0.25)/Tpole;
      T[i] = 1.0;
      //printf("%f %f\n",T[i],w_);
      continue;
    }
    if(fabs(theta-pi_/2.)<1.e-2)
    {
      T[i] = sqrt(2./(2.+w_*w_))*pow((1-w_*w_),1./12.)*exp(-4./3.*w_*w_/pow((2+w_*w_),3));
      continue;
    }
  
    if(cos(theta)<1.e-2){
      //T[i] = sqrt(2./(2.+w_*w_))*pow((1-w_*w_),1./12.)*exp(-4./3.*w_*w_/pow((2+w_*w_),3));
      Ta = 1./pow((1-fratio_),4.) + pow(w_,4.)*(1-fratio_)*(1-fratio_)*sin(theta)*sin(theta);
      Tb = 2.*w_*w_*sin(theta)*sin(theta)/(1-fratio_);
      T[i] = pow((Ta-Tb),1./8.) * pow((1-w_*w_*pow((1-fratio_),3)),1./6.)/Tpole;
      continue;
    }
    //inversR = pow(cos(theta[i]),2.)/pow((1-fratio_),2)+pow(sin(theta[i]),2.);
    //absR = Solveradius_(theta);
    if(w_==0){
      absR=1;
    } else{
      absR = Solveradius_(theta);
    }

   
  
   
   tau = Caltau_(theta,absR);

   Thalf = Solveangle_(tau,theta);
   tanTheta = (2*Thalf)/(1-Thalf*Thalf);
   Ta = 1./pow(absR,4.) + pow(w_,4.)*absR*absR*sin(theta)*sin(theta);
   Tb = 2.*w_*w_*sin(theta)*sin(theta)/absR;
    T[i] = pow((Ta-Tb),1./8.) * sqrt(fabs(tanTheta/tan(theta)))/Tpole;
   //T[i] = pow((Ta-Tb),1./8.) * sqrt(fabs(tanTheta/tan(theta)));
  
   //printf("%f %f\n",theta,T[i]);
    delete [] rtemp;
    delete [] rnew;
  }
}

double LaraModelfilter::Feff_(){
  double feff, f2=(1-fratio_)*(1-fratio_);
  feff = 1- sqrt(f2*cos(phi_)*cos(phi_)+sin(phi_)*sin(phi_));
  return feff;
}
