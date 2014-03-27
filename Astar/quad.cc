#include "quad.h"
#include "model.h"
Quad::Quad(Model &model)
  : model_(model){
  if(model_.dimen()!=1){
    printf("Error: the dimention of model is %d, quad only do 1-D integration!\n ", model_.dimen());
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

Quad::~Quad(){

}

double Quad::Integrate(double x1, double x2){
  //x2>x1
  double xm = 0.5*(x1+x2);
  double xr = 0.5*(x2-x1);
  double s = 0;
  for (int j=0; j<5; j++){
    double dx = xr*absc_[j];
    double xm1 = xm+dx;
    double xm2 = xm-dx;
    s+=w_[j]*(model_.fx(&xm1,1) + model_.fx(&xm2,1));
  }
  return s *= xr;
}
