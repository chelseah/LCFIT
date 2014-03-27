#ifndef EGGBOX_H
#define EGGBOX_H
#include <stdio.h>
#include "model.h"
class Eggbox : public Model{
  public: 
    Eggbox();
    ~Eggbox();
    double fx(double* x, int nx);
    int dimen() const {return kDimen_;}
  private:
    static const int kDimen_ = 2;
};
#endif
