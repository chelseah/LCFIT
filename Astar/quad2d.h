#ifndef QUAD2D_H
#define QUAD2D_H
#include "integrator.h"
#include <stdio.h>
#include <stdlib.h>
#include "model.h"
class Model;
class QuadTwoD : public Integrator {
  public:
    //do integration of \int f(x,y) dxdy; with xlim, x1,x2; ylim, y1,y2
    QuadTwoD(Model &model);
    double Integrate(double* xmin, int nmin, double* xmax, int nmax);
    ~QuadTwoD();
    Model &model_; //function to evaluate f(x,y)
  private:
    static const int kdimen_ = 2;
    double Integrate1D_(double x1, double x2);
    double Integrate2D_(double x1, double x2,double y);
    double xmin_[kdimen_], xmax_[kdimen_];
    double absc_[5], w_[5];
};
#endif
