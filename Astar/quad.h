#ifndef QUAD_H
#define QUAD_H
#include "cmath"
#include "cstdio"
#include <stdlib.h>
#include <stdio.h>
#include "integrator.h"
class Model;
class Quad : public Integrator {
  public:
    //do integration of \int f(x) dx; with xlim, x1,x2; 
    Quad(Model &model);
    double Integrate(double x1, double x2);
    ~Quad();
  private:
    //const int dimen_;
    Model &model_; //function to evaluate f(x,y)
    //const to be used in functions; abscissas and weights
    double absc_[5],w_[5];
};
#endif
