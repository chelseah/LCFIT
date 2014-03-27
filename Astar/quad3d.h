#ifndef QUAD3D_H
#define QUAD3D_H
#include "integrator.h"
#include <stdio.h>
#include <stdlib.h>
#include "model.h"
class Model;
class QuadTriD : public Integrator {
  public:
    //do integration of \int f(x,y) dxdy; with xlim, x1,x2; ylim, y1,y2
    QuadTriD(Model &model);
    double Integrate(double *xmin, int nmin, double *xmax, int nmax);
    ~QuadTriD();
    Model &model_; //function to evaluate f(x,y)
  private:
    static const int kdimen_=3;
    double Integrate1D_(double x1, double x2);
    double Integrate2D_(double x1, double x2,double y);
    double Integrate3D_(double x1, double x2, double *y);
    double xmin_[kdimen_],xmax_[kdimen_];
    double absc_[5], w_[5];
};
#endif
