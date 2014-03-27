#ifndef POLY_H
#define POLY_H
#include <stdio.h>
#include "model.h"
class Polynomial : public Model{
  public: 
    Polynomial(double* coeff, int lenc);
    ~Polynomial();
    double fx(double* x, int nx);
    int dimen() const {return kDimen_;}
  private:
    static const int kDimen_ = 1;
    int order_;
    double* coeff_;
};
#endif
