#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "TreeUtils.h"
#include "Maximize.h"
#include "powell.h"
#include "NRCvectorUtils.h"
#include "DistrFuncs.h"
#include "memory.h"
#include "TimeAlgorithms.h"
#include "ObjFunc.h"
#include "TNwrapper.h"

/*****************************************************************************************
/*****************************************************************************************
/* N DIMENSIONAL MINIMIZATION BY POWELL'S METHOD AND Quasi-NEWTON METHOD or TN method

	PARAMETERS:
		 'func'			-	name of function to be minimized, of the form 
		 				double func(double p[])
		 'p[numvar]'	- 	a 1-offset vector containing the estimates on return
		 					and containing a guess on entry	
		 'ftol'			-	fractional tolerance required
		 'fmin'			-	function value at minimum
		 
		
	COMMENTS:
		make sure to negate the function if you want to maximize it;
		no end of trouble if parameters don't match because of careless 
			defining of ints to longs in some of these files
		
	Q-NEWTON:
	Imposes the ancillary termination criterion of Gill et al. p. 306:U3, that is
	checks to see if the gradient is near 0.

*/



double MinND(TREE t,int method, int algorithm,double (*func)(double p[]),void (*gradient)(double [], double []),
	double p[], int numvar,int *numIter, double ftol,double linMinDelta,int *success )

// 'method' seems irrelevant at this point, not referred to below

{

int i,j,k,ierror;
double **xdir,fmin,*y, *pSave, *D,nm,crit;
NODETYPE** nodeArray;

switch (algorithm)

{
case POWELL:

	xdir=matrix(1,numvar,1,numvar);
	for (i=1;i<=numvar;i++) 
		for (j=1;j<=numvar;j++) 
		{
		if (i==j) xdir[i][i]=1.0; 
		else xdir[i][j]=0.0;
		}
	*success=powell1(p,xdir,numvar,ftol,numIter, &fmin, func);

	free_matrix(xdir, 1,numvar,1,numvar);
	break;
 
case QNEWT:
	*success=dfpmin(p,numvar,ftol,numIter,&fmin,func,gradient);
	break;

case TN:
	fmin=func(p); /* TN wants an initial guess at the function value */ 
	ierror=TNwrapper
		(
		numvar,
		p,
		func,	
		gradient,
		&fmin
		);
	if (ierror==0 || ierror==3 ) /* NB! Sometimes an error return code of 3 gives an OK result...see tn.f docs */
		*success = 1;
	else
		*success = 0;
	break;
}
return fmin;

}


