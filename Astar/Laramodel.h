#ifndef LARAMODEL_H
#define LARAMODEL_H
#include <stdio.h>
#include <math.h>
#include "model.h"
class LaraModel : public Model{
  public: 
    LaraModel(double fratio,double phi,double Req,  double ggraveq, double groteq);
    ~LaraModel();
    //a in cgs/rsun, b is b*Rpole;
    void Cal_Txy(double *x, int nx, double *y, int ny, double* Txy, int nT);
    void Cal_F(double *phase, int np, double *F, int nf,double theta, double a, double b);
    //void Cal_F0(double *F0);
    double fx(double* x, int nx) ;

    void CalTtheta(double *theta, int nth, double *T, int nT, double absR);
    int dimen() const {return kDimen_;}
   private:
    static const int kDimen_ = 2;
    //QuadTwoD* quad2d;
    double Determinant_(double x,double y);
    double Calzcoord_(double x,double y,double d);
    double Solveangle_(double tau,double theta);
    double AngleFx_(double x,double tau);
    double AngleDfx_(double x);
    
    double Solveradius_(double theta);
    double RFx_(double x,double sintheta);
    double RDfx_(double x,double sintheta);
    double Feff_();
    //old coordinate is x[0],x[1],x[2]; new coordinate ix xnew[0][1][2]
    void Rotate_(double *x, double *xnew);
    double Caltau_(double theta, double r);
    void CalTeff_(double *r, int np, int nd, double *g, int ng);
    double fratio_, phi_, Req_, ggraveq_, groteq_,w_;
    double pi_,tol_;
    int Maxiter_;
};
#endif
