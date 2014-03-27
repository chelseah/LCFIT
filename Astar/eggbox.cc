#include "eggbox.h"
//#include <stdio.h>
#include <math.h>
Eggbox::Eggbox( ){
}

Eggbox::~Eggbox(){
}

double Eggbox::fx(double *x, int nx){
  double result=0.;
//result = exp(pow((2.+cos(x[0]/2.0)*cos(x[1]/2.0)),5.));
  result = pow((2.+cos(x[0]/2.0)*cos(x[1]/2.0)),5.);
  //printf("result=%f\n",result);
  return result;
}
