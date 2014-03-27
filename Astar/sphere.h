#ifndef SPHERE_H
#define SPHERE_H
#include <stdio.h>
#include <math.h>
#include "model.h"
class Sphere : public Model{
  public: 
    Sphere();
    ~Sphere();
    double fx(double* x,int nx);
    int dimen() const {return kDimen_;}
  private:
    static const int kDimen_ = 3;
};
#endif
