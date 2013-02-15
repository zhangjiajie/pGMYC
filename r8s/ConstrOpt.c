/*   General nonlinear minimization algorithm with or without constraints.
	See further details under 'penalty' module.
	This should work with any objective function! 

Returns the number of iterations used in the first optimization;
and the number of iterations used in the LAST of the restarts.


*/

#define	CONSOLE 0	/* set to 0 for inclusion in MacRate */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Maximize.h"
#include "NRCvectorUtils.h"
#include "ConstrOpt.h"
#include "MyUtilities.h"
#include "TreeUtils.h"
#include "memory.h"
#include "ObjFunc.h"
#include "penalty.h" // don't delete

/***************************************************************/


double 		rk;		/* barrier factor */

/***************************************************************/

int ConstrOpt(
	TREE		t,
	struct	
	  NexDataType	*Nex,
	int		isConstrained,		/* 0 = unconstrained internal times */
	int		numVar,



	double		p[],
	double 		(*objective)(double p[]),
	void		(*gradient)(double [], double []),
	int		method,	
	int		algorithm,
	double 		(*penalty)(double p[]),
	
	double		*max_obj,
	int		*numPowellIter,
	int		*numRestartIter,
	int		*numBarrierIter
	
	)
	
{


	int		verbose;
	double		PowellTol;
	double		barrierTol;
	int		maxPowellIter;
	int		maxBarrierIter;
	double		initBarrierFactor;
	double		barrierMultiplier;
	double		linminOffset;
	double		contractFactor;
	int		maxContractIter;
	int		restarts;
	double		perturb_factor;
	double		local_factor; 

extern int powellMode;

extern	double		gContractFactor;	
extern	int		gMaxContractIter;	/* declared in NRCoptimize module */
extern 	int		isFeasible, gNVar;
extern NODETYPE*	gRoot;
extern int		gPowellTrace;
double	Ftest,fmin,Fsave,test;
double *pSave;
int i,j, success, k, jj, fail_flag, fl;
double minfmin;

verbose=Nex->RateBlockParms.verbose;
PowellTol=Nex->RateBlockParms.ftol;
barrierTol=Nex->RateBlockParms.barrierTol;
maxPowellIter=Nex->RateBlockParms.maxIter;
maxBarrierIter=Nex->RateBlockParms.maxBarrierIter;
initBarrierFactor=Nex->RateBlockParms.initBarrierFactor;
barrierMultiplier=Nex->RateBlockParms.barrierMultiplier;
linminOffset=Nex->RateBlockParms.linminOffset;
contractFactor=Nex->RateBlockParms.contractFactor;
maxContractIter=Nex->RateBlockParms.maxContractIter;
restarts=Nex->RateBlockParms.num_restarts;
perturb_factor=Nex->RateBlockParms.perturb_factor;
local_factor=Nex->RateBlockParms.local_factor; 

/* some globals needed in opt routines */
gContractFactor=contractFactor;
gMaxContractIter=maxContractIter;

pSave=vector(1, gNVar);

rk=initBarrierFactor;		

if (!isConstrained)
	maxBarrierIter=1;	/* for unconstrained iteration, just do the following once */

if (isConstrained && maxBarrierIter<2)
	{
	doGenericAlert("Too few iterations for constrained optimization");
	return 0;
	}
	
for (i=1;i<=maxBarrierIter;i++)
	{
	    if (verbose > 0)
	      {
	      if (isConstrained)
		printf("\n[Barrier iteration %i]\n",i);
	      printf("...Checking the starting point...(for some barrier constant)...\n");
	      }
	    *numPowellIter=maxPowellIter; /* don't remove */
	    if (check_initial_point(objective, p))
		    {
		    if (verbose > 0)
		    	printf("...Passed...now optimizing...\n");
		    fmin=MinND(t,method,algorithm, objective,gradient,p,numVar,numPowellIter, PowellTol, 
			linminOffset, &success);
		    if (gPowellTrace && powellMode==1)
			{
			printf("\nTRACE:The proposed soln:\n");
			for (j=1;j<=gNVar;j++)
			    printf("p[%i] %f\n", j, p[j]);
			printf("TRACE: Objective function value = %f\n\n", fmin);
			}
		    Fsave=fmin;
		    for (j=1;j<=gNVar;j++)
			pSave[j]=p[j];
		    if (!success)
			{
			printf("MinND returned failure in ConstrOpt on initial search (may restart!)\n");
		/* 	return 0; *//* FATAL ERROR */
			}
		    else
			{
			if (verbose > 0)
			  printf("...Initial solution=%f\n", Fsave);
			}
		    }
	    else
		    {
		    printf("...WARNING: Trial point not feasible (barrier iter %i restart %i)\n", i, k);
		    return 0; /* FATAL ERROR */
		    }
		
	for (k=1;k<=restarts;k++)  /* Will do from 0 to r restarts from solution just found, each time
				    checking to see if it matches previous solution.  As soon as it matches
				    the routine terminates,  on the assumption that this is a true local soln;
				    otherwise after r restarts it will report a fatal error (r=# restarts)*/
	    {
	    fail_flag=0;
	    if (gPowellTrace && powellMode==1)
		printf("\nTRACE: Perturbing the trial solution and retrying...\n");
	    *numRestartIter=maxPowellIter; /* don't remove */
	/**    for (j=1;j<=gNVar;j++)
		pSave[j]=p[j]; **/
	    if (verbose > 0)
	      printf("......starting perturbation %i\n",k);
	    if (!perturb_p(p,numVar,perturb_factor)) /* perturb the soln */
		printf("......perturbation %i failed!\n",k);
	    if (gPowellTrace && powellMode==1)
		{
		printf("\nTRACE:The point perturbed from previous soln:\n");
		for (j=1;j<=gNVar;j++)
		printf("p[%i] %f\n", j, p[j]);
		}
	    if (check_initial_point(objective, p)) /* perturbed point feasible?*/
		    {				    /* ...yes, optimize */
		    fmin=MinND(t,method,algorithm,objective,gradient,p,numVar,numRestartIter, PowellTol, 
			linminOffset, &success);
		    if (!success)
			{
			printf("MinND returned failure in ConstrOpt (while retrying)\n");
			return 0; /* FATAL ERROR */
			}
		    else /* soln. found, check if its the same as saved point */
			{
			if (fmin < Fsave)  /* this is a better soln;, save it */
			    {
			    Fsave = fmin;
			    for (j=1;j<=gNVar;j++)
				pSave[j]=p[j];
						
			    }
			if (verbose > 0)
			  printf("......solution for perturbation %i=%f...best=%f\n", k,fmin, Fsave);


#if 0
			if (!same_points(pSave, p, gNVar, local_factor))
			    {
			    fail_flag=1;
			    }
			 else
			    continue; /* if are the same points, don't do any more retries */
#endif
			}
		    }
	    else
		    {
		    printf("WARNING: ConstrOpt: Initial RETRY point not feasible...keeping initial soln. (barrier iter %i restart %i)\n", i, k);
		    *max_obj = Fsave;
		    }
	    }
	    for (j=1;j<=gNVar;j++)
		p[j]=pSave[j];
	    *max_obj = Fsave;    
	    
#if 0
	if(fail_flag)    
	    {
	    printf("FATAL ERROR--Inflection point or non-optimum still found after %i retries\n", k);
	    for (j=1;j<=gNVar;j++)
				printf("%i %f (%f)\n", j, p[j], pSave[j]);
	    return 0; /* FATAL ERROR (can only get to here if last retry was different than previous soln. */
	    }	
#endif			
	*numBarrierIter=i;
	if (isConstrained)  /** ...some of this may need midification ***/
			{
			if  (i==1)
				Fsave=2*fmin;	    /* this will force the routine to examine at least
							two values of the barrier contract factor */
			Ftest= fmin-penalty(p);	    /* this is the true value of the function */
			test=(Ftest-Fsave)/Ftest;   /* calculates a fractional tolerance */
			if (fabs(test) < barrierTol)
			    break;
			Fsave=Ftest;		
			rk*=barrierMultiplier;		/* Adjust factor in penalty function*/
			*max_obj = Ftest;
			}
		/* Note that each iteration uses the p[] estimate of the previous iteration
			as its initial guess */

/* NB.! THIS WHOLE ROUTINE HAS BEEN MODIFED WITHOUT CHECKING THE
 * CONSTRAINED PART OF THE ALGORITHM!
 */

	}



free_vector(pSave, 1, gNVar);
return 1;
}

