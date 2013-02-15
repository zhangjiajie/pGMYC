#include <math.h>
#include <stdlib.h>
#include <stdio.h>

double triadLike(double t1, double t2, double t3, 
					double *xt, double L1, double L2, double L3);
int feasible(double p[]);
double objective(double p[]);
double triadObs(double t1, double t2, double t3, double tint, 
	double L1, double L2, double L3);
double penalty(double p[]);
double addconstr(double x);
double BranchLike(double rate, double timeLength, double charLength);

