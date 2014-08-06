#include "quad3d.h"
#include "model.h"
QuadTriD::QuadTriD(Model &model)
  : model_(model){
  if(model_.dimen()!=kdimen_){
    printf("Error: the dimention of model is %d, quad only do 3-D integration!\n ", model_.dimen());
    exit(1);
  }
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

QuadTriD::~QuadTriD(){

}

double QuadTriD::Integrate(double* xmin, int nmin, double* xmax, int nmax){
  for (int i=0; i<kdimen_; i++){
    xmin_[i] = xmin[i];
    xmax_[i] = xmax[i];
  }
  return Integrate1D_(xmin_[0],xmax_[0]);
}

double QuadTriD::Integrate1D_(double x1, double x2){
  //x2>x1
  double xm = 0.5*(x1+x2);
  double xr = 0.5*(x2-x1);
  double s = 0;
  for (int j=0; j<5; j++){
    double dx = xr*absc_[j];
    double xm1 = xm+dx;
    double xm2 = xm-dx;
    s+=w_[j]*(Integrate2D_(xmin_[1],xmax_[1],xm1) + Integrate2D_(xmin_[1],xmax_[1],xm2));
  }
  return s *= xr;
}
double QuadTriD::Integrate2D_(double x1, double x2, double y){
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
    s+=w_[j]*(Integrate3D_(xmin_[2],xmax_[2],xm1) + Integrate3D_(xmin_[2],xmax_[2],xm2));
  }
  delete [] xm1;
  delete [] xm2;
  return s *= xr;
}
double QuadTriD::Integrate3D_(double x1, double x2, double *y){
  //x2>x1, here x1 is y1, x2 is y2
  double xm = 0.5*(x1+x2);
  double xr = 0.5*(x2-x1);
  double s = 0;
  double *xm1 = new double [kdimen_];
  double *xm2 = new double [kdimen_];
  xm1[0] = y[0]; xm1[1] = y[1];
  xm2[0] = y[0]; xm2[1] = y[1];
  for (int j=0; j<5; j++){
    double dx = xr*absc_[j];
    xm1[2] = xm+dx;
    xm2[2] = xm-dx;
    s+=w_[j]*(model_.fx(xm1,kdimen_) + model_.fx(xm2,kdimen_));
  }
  delete [] xm1;
  delete [] xm2;
  return s *= xr;
}