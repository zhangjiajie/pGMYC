#include "TimeAlgorithms.h"
int BFGSwrapper
	(
	int numvar,
	double x[],
	double 		(*objective)(double []),
	void		(*gradient)(double [], double []),
	double		*max_obj
	);
