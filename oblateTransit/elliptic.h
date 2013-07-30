#ifndef ELLIPTIC_H
#define ELLIPTIC_H
#include <math.h>
#include <iostream>
class Elliptic{
  public:
    Elliptic();
    ~Elliptic();
    void Calee(double *k, int lk, double *ee, int lee);
    void Calkk(double *k, int lk, double *kk, int lkk);
    void Ellpic_bulirsch(double *n, int ln, double *k, int lk, double *nk, int lnk);
  private:
     static const double ea1_=0.44325141463;
     static const double ea2_=0.06260601220;
     static const double ea3_=0.04757383546;
     static const double ea4_=0.01736506451;
     static const double eb1_=0.24998368310;
     static const double eb2_=0.09200180037;
     static const double eb3_=0.04069697526;
     static const double eb4_=0.00526449639;

     static const double ka0_=1.38629436112;
     static const double ka1_=0.09666344259;
     static const double ka2_=0.03590092383;
     static const double ka3_=0.03742563713;
     static const double ka4_=0.01451196212;
     static const double kb0_=0.5;
     static const double kb1_=0.12498593597;
     static const double kb2_=0.06880248576;
     static const double kb3_=0.03328355346;
     static const double kb4_=0.00441787012;
     static const double pi_ =3.1415926535898;
};
#endif
