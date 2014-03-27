#ifndef LARA_H
#define LARA_H
#include <stdio.h>
#include <math.h>
#include "Laramodelfilter.h"
#include "quadpolar.h"
#include "quad2d.h"
class Larafilter {
  public: 
    Larafilter(double *filter, int nf, int nx, int ny, int nz,double x1, double dx, double y1, double dy, double z1, double dz, double fratio,double phi,double Req,  double ggraveq, double groteq);
    ~Larafilter();
    //a in cgs/rsun, b is b*Rpole;
    double Cal_F0();
    double Cal_Lum();
    void Resetphi(double phi);
    //void Cal_geff(double x);
    //double fx(double* x) ;
    //int dimen() const {return kDimen_;}
  private:
    //static const int kDimen_ = 2;
    LaraModelfilter* model_;
    QuadPolar* quadpolar_;
    QuadTwoD* quad2d_;
    //double Determinant_(double x,double y);
    //double Calzcoord_(double x,double y,double d);
    //double Feff_();
    //old coordinate is x[0],x[1],x[2]; new coordinate ix xnew[0][1][2]
    //void Rotate_(double *x, double *xnew);
    //void Calgeff_(double *x, int nx, double *y, int ny, double *g, int ng);
    double Req_;
    double pi_;
};
#endif
