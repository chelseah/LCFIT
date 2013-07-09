#include"oblateness.h"
#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0 - EPS)

float ran1(long *idum)
{
  int j;
  long k;
  static long iy=0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0 || !iy) {
    if (-(*idum) < 1) *idum=1;
    else *idum = -(*idum);
    for (j=NTAB+7;j>=0;j--){
      k = (*idum)/IQ;
      *idum = IA*(*idum - k*IQ)-IR*k;
      if (*idum<0) *idum += IM;
      if (j<NTAB) iv[j] = *idum;
    }
    iy = iv[0];
  }

  k= (*idum)/IQ;
  *idum = IA*(*idum-k*IQ)-IR*k;
  if (*idum < 0 ) *idum += IM;
  j = iy/NDIV;
  iy = iv[j];
  iv[j] = *idum;
  if ((temp=AM*iy) > RNMX) return RNMX;
  else return temp;
}
static long seed=0;

/* limb darkening function */
double limbDarkening(double u1, double u2, double r2)
{
	double t = 1-r2;
	return((u1+2*u2)*sqrt(t)-u2*t);
}

/* find the intersection point on the line (x1,y1) to (x2,y2) */
void findRoot(double x1, double y1, double x2, double y2, double xc, double yc, double *xval, double *yval, int *flag)
{
	double d1, d2, d;
	d1 = (x1-xc)*(x1-xc)+(y1-yc)*(y1-yc)-1;
	if(fabs(d1)<1e-12) { *xval = x1; *yval = y1; *flag = 1; return;}
	d2 = (x2-xc)*(x2-xc)+(y2-yc)*(y2-yc)-1;
	if(fabs(d2)<1e-12) { *flag = 0; return;} /* to avoid double counting */
	/* if d1*d2>0, no intersection point is there */
	if(d1*d2>0) { *flag = 0; return;}

	double x, y;
	while(fabs(x1-x2)>1e-6 || fabs(y1-y2)>1e-6)
	{
		x = (x1+x2)/2.0; y = (y1+y2)/2.0;
		d = (x-xc)*(x-xc)+(y-yc)*(y-yc)-1;
		if(fabs(d)<1e-12) { 
			*xval = x; *yval = y; 
			*flag = 1; 
			return;
		}
		d1 = (x1-xc)*(x1-xc)+(y1-yc)*(y1-yc)-1.0;
		if(d*d1>0) {x1 = x; y1 = y;}
		else {x2 =x; y2 = y;}
	}
	d1 = (x1-xc)*(x1-xc)+(y1-yc)*(y1-yc)-1.0;
	if(d1>0) {
		*xval = x2; *yval = y2; *flag = 1;
	}
	else {
		*xval = x1; *yval = y1; *flag = 1;
	}
}

/* find the intersection points on the ellipse */
void intersectionPoints(double *x, double *y, int *n, double a, double b, double xc, double yc)
{
	int N=500, i, count=0;
	double phi1, phi2, dphi=2*pi/N;
	double xval, yval;
	double x1, y1, x2, y2;
	int flag;
	for(i=0; i<N; i++)
	{
		phi1 = i*dphi; phi2 = (i+1)*dphi;
		x1 = a*cos(phi1); x2 = a*cos(phi2);
		y1 = b*sin(phi1); y2 = b*sin(phi2);
		findRoot(x1, y1, x2, y2, xc, yc, &xval, &yval, &flag);
		if(flag == 0) continue;
		if(count > 1) {printf("More intersection points than expected!\n"); exit(1);}
		x[count] = xval; y[count] = yval; count++;
		/* there is one intersection point between (x1,y1) and (x2,y2) */
	}
	*n = count;
	return;
}

/* area calculation */
double triangleArea(double *x, double *y)
{
	double dx = x[1]-x[0];
	double dy = y[1]-y[0];
	double d1 = sqrt(dx*dx+dy*dy);
	double k, h;
	if(fabs(dx)<1e-6) h = fabs(x[1]);
	else {
		k = dy/dx;
		h = fabs(y[0]-k*x[0])/sqrt(1+k*k);
	}
	double area = 0.5*d1*h;
	return(area);
}

double circleChordArea(double *x, double *y, double xc, double yc, double rc)
{
	double x1=x[0]-xc, y1=y[0]-yc;
	double x2=x[1]-xc, y2=y[1]-yc;
	double dbeta;
	double l12;
	l12 = (x1-x2)*(x1-x2)+(y1-y2)*(y1-y2);
	double l1, l2;
	l1 = rc*rc; l2 = l1;
	dbeta = acos((l1+l2-l12)/(2*sqrt(l1*l2)));
	double area;
	area = 0.5*rc*rc*(dbeta-sin(dbeta));
	return(area);
}

/* calculate the deficite flux of the star at phase phi */
void relativeFlux(double *variables, int n, double phi, double *deficitFlux, double *circleAnalogy)
{
	/* variables */
/*  variables[0] = Req;
	variables[1] = Rpol;
	varibales[2] = alpha;
	variables[3] = sma;
	variables[4] = period;
	variables[5] = inclination;
	variables[6] = u1; limbdarkening coefficient
	variables[7] = u2; limbdarkening coefficient
	variables[8] = n; observation times
	variables[9] = percentage; phase interval
*/
	double rc=1.0;
	/* distance between centers of the planet and the star */
	double d=variables[3]*sqrt(sin(phi*2*pi)*sin(phi*2*pi)+cos(phi*2*pi)*cos(phi*2*pi)*cos(variables[5])*cos(variables[5]));
	
	/* the embedded sphere's contribution */
	//*circleAnalogy = standardCurve(variables[1], d, variables[6], variables[7]);


	/* if no transit happens */
	if(d >= (1.0+variables[0]))	{ *deficitFlux = 0.0; return; }
	
	/* angle between the major axis and the line connecting two centers */
	double beta;
	double b0=variables[3]*cos(variables[5]);
	if(phi<=0) beta = variables[2]+asin(b0/d); /* ingress */
	if(phi>0)  beta = variables[2]-pi-asin(b0/d); /* egress */

	/* star's position and intersection points */
	double xc, yc;
	xc = d*cos(beta);
	yc = d*sin(beta);
	double xinter[2], yinter[2];
	int npoints;
	intersectionPoints(xinter, yinter, &npoints, variables[0], variables[1], xc, yc);
	
	if(d>1.0 && npoints<2) { /* the planet is outside the stellar plane, and there is no intersectiosn */
		*deficitFlux = 0.0;
		return;
	}

	/* calculate the limitation region */
	double theta1, theta2, t;
	double F00; /* F(x;a,b,alpha,0,0)-F(x;b,0,0) */
	if(d<1.0 && npoints<2) { /* the planet is inside the stellar plane and there is no intersections */
		theta1 = 0.0;
		theta2 = 2*pi;
		F00 = pi*variables[1]*(variables[0]-variables[1]);
	}
	if(npoints==2) { /* two intersections */
		if(yinter[0] >= 0)
	  		theta1 = acos(xinter[0]/variables[0]);
		else
	  		theta1 = 2*pi-acos(xinter[0]/variables[0]);
		if(yinter[1] >= 0)
	  		theta2 = acos(xinter[1]/variables[0]);
		else
	  		theta2 = 2*pi-acos(xinter[1]/variables[0]);
		if(theta1 > theta2){
			t = theta1;
			theta1 = theta2;
			theta2 = t;
		}
	}

	double theta0 = (theta1+theta2)*0.5;
	double xtest, ytest;
	xtest = variables[0]*cos(theta0);
	ytest = variables[1]*sin(theta0);
	double dtest;
	dtest = (xtest-xc)*(xtest-xc)+(ytest-yc)*(ytest-yc)-1.0;
	if(dtest>0) {
		t = theta2-2*pi;
		theta2 = theta1;
		theta1 = t;
	}

	double epsilon = variables[0]/variables[1];
	double epsilon2 = epsilon*epsilon;
	double chordArea1, chordArea2;
	double dtheta=theta2-theta1;
	double ksi, psi, cc, p2;
	if(npoints == 2) {
		if(dtheta<=pi)
		  chordArea1 = dtheta*variables[0]*variables[1]*0.5-triangleArea(xinter, yinter);
		else
		  chordArea1 = dtheta*variables[0]*variables[1]*0.5+triangleArea(xinter, yinter);
		chordArea2 = circleChordArea(xinter, yinter, xc, yc, 1.0);

		if(d>=(variables[1]+1.0)) /* the two circles are separated */
		  F00 = chordArea1+chordArea2;
		else if(d>(1.0-variables[1])) { /* the two circles have intersection points */
			ksi = acos((variables[1]*variables[1]+d*d-1.0)/(2*variables[1]*d));
			psi = acos((1.0+d*d-variables[1]*variables[1])/(2.0*d));
			p2 = variables[1]*variables[1];
			cc = ksi*p2+psi-sqrt(4*d*d-(1+d*d-p2)*(1+d*d-p2))*0.5;
			F00 = chordArea1+chordArea2-cc;
		}
		else /* the small circle is included in the larger one */
		  F00 = chordArea1+chordArea2-pi*variables[1]*variables[1];
	}

	int Nrays, i;
	Nrays = 3000;
	double eta1, eta2, theta, R, x, y;
	double contribution=0.0, separation1, separation2;
	double total=0.0;
	
	for(i=1; i<=Nrays; i++)
	{
		eta1 = ran1(&seed);
		eta2 = ran1(&seed);
		theta = (1.0-eta2)*theta1+eta2*theta2;
		R = sqrt(eta1+(1-eta1)/epsilon2);
		x = variables[0]*R*cos(theta);
		y = variables[1]*R*sin(theta);
		/* distance to the center of the star */
		separation1 = (x-xc)*(x-xc)+(y-yc)*(y-yc);
		/* if the ray locates outside the stellar plane, or inside the circular plane of planet, ignore it */
		if(separation1 > 1.0) continue;
		/* distance to the center of the planet */
		separation2 = x*x+y*y;
		if(separation2 < variables[1]*variables[1]) continue;
		/* otherwise, add its contribution on */
		contribution += limbDarkening(variables[6], variables[7], separation1);
	}
	
	double totalArea, Istar;
	totalArea = dtheta*(variables[0]*variables[0]-variables[1]*variables[1])*0.5/epsilon;
	Istar = contribution*totalArea/Nrays;
	*deficitFlux = (1.0-variables[6]-variables[7])*F00+Istar;
}


