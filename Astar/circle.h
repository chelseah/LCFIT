#ifndef CIRCLE_H
#define CIRCLE_H
#include <stdio.h>
#include "model.h"
class Circle : public Model{
  public: 
    Circle();
    ~Circle();
    double fx(double* x,int nx);
    int dimen() const {return kDimen_;}
  private:
    static const int kDimen_ = 2;
};
#endif
