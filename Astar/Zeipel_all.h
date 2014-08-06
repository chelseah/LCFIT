#ifndef ZEIPEL_ALL_H
#define ZEIPEL__ALL_H
#include <stdio.h>
#include <math.h>
class Zeipel {
  public: 
    Zeipel(double fratio,double phi,double Req,  double ggraveq, double groteq, double beta);
    ~Zeipel();
    //a in cgs/rsun, b is b*Rpole;
    void Cal_F(double *phase, int np, double *F, int nf,double theta, double a, double b);
    void Cal_F0(double *F0,int np);
    double fx(double* x,int nx) ;
    int dimen() const {return kDimen_;}
  private:
    static const int kDimen_ = 2;
    double Determinant_(double x,double y);
    double Calzcoord_(double x,double y,double d);
    double Integrate(double x1, double x2, double y1, double y2);
    double Integrate1D_(double x1, double x2);
    double Integrate2D_(double x1, double x2,double y);
    double x1_,x2_, y1_,y2_;
    double Feff_();
    //old coordinate is x[0],x[1],x[2]; new coordinate ix xnew[0][1][2]
    void Rotate_(double *x, double *xnew);
    void Calgeff_(double *x, int nx, double *y, int ny, double *g, int ng);
    double fratio_, phi_, Req_, ggraveq_, groteq_,beta_;
    double pi_,absc_[5],w_[5];
};
#endif
