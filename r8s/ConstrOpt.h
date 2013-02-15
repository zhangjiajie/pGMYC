#include "nexus.h"
int ConstrOpt(
	TREE		t,
	struct	
	  NexDataType	*Nex,
	int		isConstrained,	
	int		numVar,
	double		p[],
	double 		(*objective)(double p[]),
	void		(*gradient)(double [], double []),
	int		method,
	int		algorithm,
	double 		(*penalty)(double p[]),
	
	double		*maxLike,
	int			*numIter,
	int			*numRestartIter,
	int			*numBarrierIter
	
	)
;

