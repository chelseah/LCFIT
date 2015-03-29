#ifndef ZEIPELMODEL_H
#define ZEIPELMODEL_H
#include <stdio.h>
#include <math.h>
#include "model.h"
class ZeipelModel : public Model{
  public: 
    ZeipelModel(double fratio,double phi,double Req,  double ggraveq, double groteq, double beta, double Rp);
    ~ZeipelModel();
    //a in cgs/rsun, b is b*Rpole;
    void Cal_F(double *phase, int np, double *F, int nf,double theta, double a, double b);
    //void Cal_F0(double *F0);
    double fx(double* x, int nx) ;
    double Feff_();
    int dimen() const {return kDimen_;}
  private:
    static const int kDimen_ = 2;
    //QuadTwoD* quad2d;
    double Int_F(double x, double y);
    double fxShift(double* x, int nx,double xs,double ys, bool areaflag) ;
    double Determinant_(double x,double y);
    double Calzcoord_(double x,double y,double d);
    //old coordinate is x[0],x[1],x[2]; new coordinate ix xnew[0][1][2]
    void Rotate_(double *x, double *xnew);
    void Calgeff_(double *r, int np, int nd, double *g, int ng);
    double InterArea(double x, double y); //calculated the fraction intersect area between the planet and the star
    double IntegrateShift_(double * xmin, int nmin, double *xmax, int nmax, double x, double y, bool areaflag);
    double Integrate1DShift_(double x1, double x2, double x , double y, bool areaflag);
    double Integrate2DShift_(double x1, double x2, double y, double xs,double ys, bool areaflag);
    double fratio_, phi_, Req_, ggraveq_, groteq_,beta_, Rp_;
    double pi_;
    double absc_[5], w_[5];
    double xmin_[2], xmax_[2];
};
#endif
