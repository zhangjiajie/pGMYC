#include "TreeUtils.h"
double Min1D(double (*func)(double x),double *xmin,
						double ax,double bx, double cx,double ftol );
void plot2d(double low1, double high1, double low2, double high2, int gridSize);
double MinND(TREE t,int method, int algorithm, double (*func)(double p[]),void (*grad)(double [], double []),
		double p[], int numvar,int *iter, double ftol,
		double linMinDelta,int *success);
