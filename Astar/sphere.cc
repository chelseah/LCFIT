#include "sphere.h"
#include <stdio.h>
Sphere::Sphere(){
}

Sphere::~Sphere(){
}

double Sphere::fx(double *x, int nx) {
  //x[0] = r; x[1] = phi (0-2pi); x[2] = theta (0-pi); spherical coordinate

  return x[0]*x[0]*sin(x[2]);
}
