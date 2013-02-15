#include "TimeAlgorithms.h"
int BFGSwrapper
	(
	int numvar,
	double x[],
	double 		(*objective)(double []),
	void		(*gradient)(double [], double []),
	double		*max_obj
	);
void sfun_(int *N,double X[],double *F, double G[]);
int TNwrapper
	(
	int numvar,
	double x[],
	double 		(*objective)(double []),
	void		(*gradient)(double [], double []),
	double		*max_obj
	);
