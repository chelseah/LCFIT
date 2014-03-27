#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "model.h"
#include "poly.h"
#include "circle.h"
#include "sphere.h"
#include "eggbox.h"
#include "quad.h"
#include "quad2d.h"
#include "quadpolar.h"
#include "quad3d.h"
#include "quadsphere.h"

int main(int argc, char *argv[]){
  //test polynomial bit of the code
  Polynomial* model;
  //first group, order 0
  double *coeff1,*coeff2, *coeff3;
  double x;
  coeff1 = new double [1];
  coeff2 = new double [2];
  coeff3 = new double [3];
  printf("-----------test polynomial--------------------------\n");
  coeff1[0] = 2;
  model = new Polynomial(coeff1,1);
  x = 2.;
  printf("coeff=%f,x=%f,answer=%f\n", coeff1[0],x,model->fx(&x,1));
  coeff2[0] = 2; coeff2[1] = 3;
  delete model;
  model = new Polynomial(coeff2,2);
  x = 2.;
  printf("coeff=%f,%f,x=%f,answer=%f\n", coeff2[0], coeff2[1],x,model->fx(&x,1));
  coeff3[0] = 2; coeff3[1] = 3; coeff3[2] = 4;
  x = 2.;
  delete model;
  model = new Polynomial(coeff3,3);
  printf("coeff=%f,%f,%f,x=%f,answer=%f\n", coeff3[0],coeff3[1],coeff3[2],x,model->fx(&x,1));
  delete model;
  printf("-----------test quad--------------------------\n");
  Quad *quad;
  double fx,x1,x2;
  x1 = 1.; x2 = 2.;
  quad = new Quad(*model);
  fx = quad->Integrate(x1,x2);
  printf("coeff=%f,%f,%f,x1=%f,x2=%f, fx=%f\n", coeff3[0],coeff3[1],coeff3[2],x1,x2,fx);
  printf("-----------test quad2d--------------------------\n");
  QuadTwoD *quad2d;
  double xmin[2],xmax[2];
  Circle* model2D;
  model2D = new Circle();
  quad2d = new QuadTwoD(*model2D);
  xmin[0] = 0; xmax[0] = 2.; xmin[1] = 0; xmax[1] = 3.1415926535897931*2.;
  fx = quad2d->Integrate(xmin,2,xmax,2);
  printf("x1=%f,x2=%f,y1=%f,y2=%f, fx=%f\n", xmin[0],xmax[0],xmin[1],xmax[1],fx);
  delete quad2d;
  //printf("-----------test quad2d eggbox------------------------\n");
  //Eggbox* modeleggbox;
  //modeleggbox = new Eggbox();
  //quad2d = new QuadTwoD(*modeleggbox);
  //xmin[0] = 0; xmax[0] = 3.1415926535897931*4.; xmin[1] = 0; xmax[1] = 3.1415926535897931*4.;
  //fx = quad2d->Integrate(xmin,2,xmax,2);
  //printf("x1=%f,x2=%f,y1=%f,y2=%f, fx=%f\n", xmin[0],xmax[0],xmin[1],xmax[1],log(fx));

  //
  //
  //printf("-----------test quad3d--------------------------\n");
  //QuadTriD *quad3d;
  //double xmin3d[3],xmax3d[3];
  //Sphere* model3D;
  //model3D = new Sphere();
  //quad3d = new QuadTriD(*model3D);
  //xmin3d[0] = 0; xmax3d[0] = 2.; xmin3d[1] = 0; xmax3d[1] = 3.1415926535897931*2.; xmin3d[2] = 0; xmax3d[2] = 3.1415926535897931;
  //fx = quad3d->Integrate(xmin3d,3,xmax3d,3);
  //printf("x1=%f,x2=%f,y1=%f,y2=%f, z1=%f,z2=%f,fx=%f\n", xmin3d[0],xmax3d[0],xmin3d[1],xmax3d[1],xmin3d[2],xmax3d[2],fx);
  printf("-----------test quadpolar--------------------------\n");
  QuadPolar *quadpolar;
  coeff1[0] = 1;
  model = new Polynomial(coeff1,1);
  quadpolar = new QuadPolar(*model);
  xmin[0] = 0; xmax[0] = 2.; xmin[1] = 0; xmax[1] = 3.1415926535897931*2.;
  fx = quadpolar->Integrate(xmin,2,xmax,2);
  printf("x1=%f,x2=%f,y1=%f,y2=%f, fx=%f\n", xmin[0],xmax[0],xmin[1],xmax[1],fx);
  //printf("-----------test quadsphere--------------------------\n");
  //QuadSphere *quadsphere;
  //quadsphere = new QuadSphere(*model);
  //xmin3d[0] = 0; xmax3d[0] = 2.; xmin3d[1] = 0; xmax3d[1] = 3.1415926535897931*2.; xmin3d[2] = 0; xmax3d[2] = 3.1415926535897931;
  //fx = quadsphere->Integrate(xmin3d,3,xmax3d,3);
  //printf("x1=%f,x2=%f,y1=%f,y2=%f, z1=%f,z2=%f,fx=%f\n", xmin[0],xmax[0],xmin[1],xmax[1],xmin[2],xmax[2],fx);



  //delete [] model;
  delete [] coeff1;
  delete [] coeff2;
  delete [] coeff3;
  delete quad;
  delete quad2d;
  //delete quadpolar;
  return 0;
}
