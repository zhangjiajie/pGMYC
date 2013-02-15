int dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	double(*func)(double []), void (*dfunc)(double [], double []));

double f1dim1(double x);
void mnbrak1(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
	double (*func)(double));
void linmin1(double p[], double xi[], int n, double *fret, 
	double (*func)(double []));
double brent1(double ax, double bx, double cx, double (*f)(double), double tol,
	double *xmin);
int powell1(double p[], double **xi, int n, double ftol, int *iter, double *fret,
	double (*func)(double []));
void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []));

