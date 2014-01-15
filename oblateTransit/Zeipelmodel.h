#ifndef ZEIPELMODEL_H
#define ZEIPELMODEL_H
#include <stdio.h>
#include <math.h>
#include "model.h"
class ZeipelModel : public Model{
  public: 
    ZeipelModel(double fratio,double phi,double Req,  double ggraveq, double groteq, double beta);
    ~ZeipelModel();
    //a in cgs/rsun, b is b*Rpole;
    void Cal_F(double *phase, int np, double *F, int nf,double theta, double a, double b);
    //void Cal_F0(double *F0);
    double fx(double* x, int nx) ;
    int dimen() const {return kDimen_;}
  private:
    static const int kDimen_ = 2;
    //QuadTwoD* quad2d;
    double Determinant_(double x,double y);
    double Calzcoord_(double x,double y,double d);
    double Feff_();
    //old coordinate is x[0],x[1],x[2]; new coordinate ix xnew[0][1][2]
    void Rotate_(double *x, double *xnew);
    void Calgeff_(double *r, int np, int nd, double *g, int ng);
    double fratio_, phi_, Req_, ggraveq_, groteq_,beta_;
    double pi_;
};
#endif
