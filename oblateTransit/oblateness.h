#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<time.h>
#include<gsl/gsl_sf_ellint.h>
#define pi 3.1415926535897932384626433832795028841971693993751

void relativeFlux(double *variables, int n, double phi, double *deficitFlux, double *circleAnalogy);
