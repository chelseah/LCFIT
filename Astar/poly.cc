#include "poly.h"
//#include <stdio.h>
Polynomial::Polynomial(double* coeff, int lenc )
  :order_(lenc){
  coeff_ = new double [order_];
  for (int i =0; i<order_; i++){
    coeff_[i] = coeff[i];
  }
}

Polynomial::~Polynomial(){
  delete [] coeff_;
}

double Polynomial::fx(double *x, int nx){
  double result=0.,powx=1.;
  for (int i=0; i<order_; i++){
    //printf("%d,%f\n",i,coeff_[i]);
    result+=coeff_[i]*powx;
    powx*=*x;
  }
  //printf("result=%f\n",result);
  return result;
}
