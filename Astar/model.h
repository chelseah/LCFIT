#ifndef Model_H_
#define Model_H_
class Model{
  public:
    virtual ~Model() {}
    //f(x) to evaluate the function to be integrated at x
    virtual double fx(double *x,int nx) =0;
    virtual int dimen() const = 0;
};
#endif
