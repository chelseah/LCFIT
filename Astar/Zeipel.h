#ifndef ZEIPEL_H
#define ZEIPEL_H
#include <stdio.h>
#include <math.h>
#include "Zeipelmodel.h"
#include "quadpolar.h"
#include "quadsphere.h"
class Zeipel {
  public: 
    Zeipel(double fratio,double phi,double Req,  double ggraveq, double groteq, double beta);
    ~Zeipel();
    //a in cgs/rsun, b is b*Rpole;
    void Cal_F(double *phase, int np, double *F, int nf,double theta, double a, double b);
    void Cal_F0(double *F0, int np);
    double Cal_Lum();
    //void Cal_geff(double x);
    //double fx(double* x) ;
    //int dimen() const {return kDimen_;}
  private:
    //static const int kDimen_ = 2;
    ZeipelModel* model_;
    QuadPolar* quadpolar_;
    QuadSphere* quadsphere_;
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
