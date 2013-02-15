#include "NRCvectorUtils.h"
#include <stdio.h>
#include <stdlib.h>
#include "ObjFunc.h"
#include "TreeUtils.h"
#include "TimeAlgorithms.h"
#include "ConstrOpt.h"
#include "objective.h"
#include "MyUtilities.h"
#include "penalty.h"
#include "DistrFuncs.h"
#include "nexus.h"
// #include "malloc.h"
#include "structures.h"
#include "TNwrapper.h"
#include "memory.h"
#include "DrawTree.h"
#include "TreeSim.h"

/* private functions */

static void tree2pTimeArray_helper(NODETYPE *node,double pTime[]);
static int children_tips_are_zeroes(NODETYPE * parent);
static void setUpLowHigh_helper(NODETYPE * node, double LOW[],double HIGH[],double minDur);
static void setUpLowHighTN(TREE t,int nvar,int tvar, double initRate,double LOW[],double HIGH[]);
static void newAllNodeArray_helper(NODETYPE *node,NODETYPE ** nodeArray);
static void newNodeArray_helper(NODETYPE *node,NODETYPE ** nodeArray);

static void assignArrayTimesToLL_helper(NODETYPE * node, double lp[]);


/* GLOBALS */

double		*gLOW,*gHIGH;
StackPtr	gFStack,gPStack;
NODETYPE * 	gRoot,*gRootDesc;
int		gRatesAreGamma;
double		gAlpha;
long		gNumSites;
int		gClampRoot;
int		gisConstrained;		/* 0 = unconstrained internal times */
double		gPowellTol;
double		gbarrierTol;
int		gmaxPowellIter;
int		gmaxBarrierIter;
double		ginitBarrierFactor;
double		gbarrierMultiplier;
double		glinminOffset;
double		gcontractFactor;
int		gmaxContractIter;

#define		MAX_FEASIBLE_TRIES 25	/* If a perturbation is not not feasible, retry this many times */
#define		PRINT_TREES	0
#define		MAX_NODES	200	/***** KLUDGE ******/
#define		INIT_RATE	50
#define		INIT_RATE2	100 
#define		INIT_RATIO	2.5
#define		INIT_RATE_FRACTION 0.75
#define		MAX_MULT_TRIES  25	/* kludge :: This must be larger than #rate* # time guesses */
double		gInitTimeFudge;		/* These are used to perturb initial
					conditions */

extern int N;
extern double S[],a[][2];
double chiSq;
double gGamma_b;
double gGamma_c;
int	gIndex;

static void doObjFuncGuts(TREE Tree, int method, int nRates, int algorithm, int * success, double initRateArg);

// ***************** Below is the new wrapper function around the older 'guts' function **********

void doObjFunc(TREE Tree, int method, int nRates, int algorithm, int * success)
{
extern struct NexDataType *gNexDataPtr;
int verbose;
if (method==PENLIKE)
	{
	
	// For PL only, use as an estimate of initRate, the value obtained by Langley Fitch clock analysis

	double initRateArg;

	traverseMultiplyLength(Tree->root, 1,1); // this just forces rounding of input branch lengths ANY TIME divtime functions are called! To make sure the gradients don't go wonky: i.e. likelihood functions convert lengths to integers but gradients don't; this way we make sure the data are ALWAYS COMPLIANT

	verbose = gNexDataPtr->RateBlockParms.verbose;
	if (verbose>0)
		printf("Using DIVTIME with clock (Langley-Fitch) options to obtain an initial guess at rates for PL search\n");
       	gNexDataPtr->RateBlockParms.verbose=0;
       	doObjFuncGuts(Tree,LaF,nRates,TN,success,0.0); // can pass dumb 0 initRate for LaF
       	gNexDataPtr->RateBlockParms.verbose= verbose; /* restore output verbositude to current value */

	if (*success)
		initRateArg=Tree->estRate;	// this was stored during the optimization
	else
		{
		printf("Could not obtain a clock-based initial rate estimate in PL setup\n");
		return;
		}

	doObjFuncGuts(Tree,method,nRates,algorithm,success,initRateArg);  
	}
else
	doObjFuncGuts(Tree,method,nRates,algorithm,success,0.0);  // for other methods, pass dumb initRate, which will be ignored
return;
}
/*****************************************************************************/

static void doObjFuncGuts(TREE Tree, int method, int nRates, int algorithm, int * success, double initRateArg)

{

  extern int 	gNVar,gFloatRoot,gIndex; /* initialized in convertLLtoArray */
  extern double gSmoothing,gFit,gLike; /* defined in TimeAlgorithms.c */
  extern int 	gVarMinFlag,gEstRoot;
  extern double /* pTime,*/  gnpexp;	/* 'gnpexp' defined in ReadNexusFile2.c */
  extern NODETYPE *gRootDescPenalty;
  extern NODETYPE * gRoot;               /* field of node structure */
  extern struct NexDataType *gNexDataPtr;
  struct NexDataType *Nex;
  struct RBP * rbp;
	
  double (*obj_func_array[10])(double[]); /* array of pointers to the various
		objective functions...indexed by 'method' which is set up below */
  void (*gradient_array[10])(double[],double []); /* array of pointers to the various
		gradient functions...indexed by 'method' which is set up below */

  NODETYPE 	*r1,*r2,*r,*child1,*child2;
  NODE		root;
  char		*TreeName;

  double	*pTime=NULL,*D=NULL;
  double	sum,max,b,c,fObj;
  double	*p;	/* pointer to the array of times, etc. */
  double	maxLike,minRate;
  double 	local_factor, perturb_factor; /* used to adjust initial guess for times */
  double	ftol,	initRate;
  
  /* ...a bunch of arrays to be...*/
  double	*Mult_Like;
  double	*Z_Dist;
  double	*Rate_hat, *Rate_hat2;
  double 	*p1,*p2;
  double 	*lp,**storeParm,**storeParm_init,**storeGrad,**storeGradInit, *approxGrad;

  int 		i,j,k,m,jj, ii,nvar,tvar,count=0,maxi,chiSqDF,was_a_failure=0,verbose,KFR,k1,k2, storek1,success_init;
  int 		rateType,penaltyType,neighborPenalty;
  int		numIter,numRestartIter;
  int		numBarrierIter;
  int	    	NUM_TIME_GUESSES;
  int	    	NUM_RATE_GUESSES;
  int	    	NUM_RESTARTS;
  int	    	NUM_MULT_TRIES;
  int		*Fail_Flag;
  int		anyZeroTerminal,anyZeroInternal;


/********************************************************************************/
/*******************************   Main code ************************************/

obj_func_array[LaF]=objLangFitch;
obj_func_array[LFLOCAL]=objLangFitchLocal;
obj_func_array[NP]=objNP;
obj_func_array[PENLIKE]=objPenLike;

gradient_array[LaF]=GradientLF;
gradient_array[PENLIKE]=GradientPL;

root=Tree->root;
TreeName=Tree->name;

  gRoot=root;	/* Global is needed for tree-based objective functions */
  gRootDescPenalty=root; /* Used in the penalty function too */
  Nex=gNexDataPtr;
  rbp=&(Nex->RateBlockParms);
  penaltyType=rbp->PenaltyType;
  neighborPenalty=rbp->NeighborPenalty;
  verbose=rbp->verbose;
  gnpexp=rbp->npexp; /***KLUDGE***/
  gClampRoot=rbp->clampRoot;
  gSmoothing=rbp->smoothing; /* Global is needed in tree-based objective functions */
  gNumSites=rbp->numSites;
  ObjFuncDefaults();

  NUM_TIME_GUESSES=rbp->num_time_guesses;
  NUM_RESTARTS=rbp->num_restarts;

  NUM_MULT_TRIES=NUM_TIME_GUESSES;
  if (NUM_MULT_TRIES>MAX_MULT_TRIES)
    fatal("Too many initial starts for static arrays");


/** check if time constraints have been set, and set global appropriately **/


/** NOTE! If using TN, we don't want to set gisConstrained, because this will calc a modified Obj func! **/

if (constraintsPresent(root) && algorithm != TN)
	gisConstrained=1;
else
	gisConstrained=0;

/** check if rates are variable across sites, and set globals for use in obj funcs**/

if(rbp->RatesAreGamma)
	{
	gRatesAreGamma=1;
	gAlpha=rbp->alpha;
	}
else
	gRatesAreGamma=0;


/****************** ...SOME WARNINGS....   *************************/
if (rbp->roundFlag==0 && (method != NP)) /* This gets set when the input trees are ultrametric */
	{
	doGenericAlert("Model-based DIVTIME routines require rounding of input branch lengths");
	*success=0;
	return;
	}
if (method==USER) /* This gets set when the input trees are ultrametric */
	{
	doGenericAlert("Tree is already ultrametric! No need for DIVTIME");
	*success=0;
	return;
	}
	

if (algorithm==QNEWT && gRatesAreGamma)
	{
	doGenericAlert("QNEWT cannot be used if rates are gamma distributed; use POWELL");
	*success=0;
	return;
	}
if ((algorithm==QNEWT) && gisConstrained)
	{
	doGenericAlert("QNEWT cannot be used if time constraints are present; use POWELL");
	*success=0;
	return;
	}

if ((algorithm==TN) && (method == LFLOCAL))
	{
	doGenericAlert("Can only use algorithm=POWELL with local clock method");
	*success=0;
	return;
	}
if ((algorithm==QNEWT) && (method == NP))
	{
	doGenericAlert("QNEWT only available for LF and PL");
	*success=0;
	return;
	}
if ((algorithm==TN) && (method == NP))
	{
	doGenericAlert("TN only available for LF and PL");
	*success=0;
	return;
	}


/* 
QNEWT method fails often with 0-length terminals. (maybe also internals). This is because
	the rate wants to go to zero and beyond, but the derivative at that point is still non-
	zero. This will have to be ultimately fixed by proper invocation of constraints and
	boundaries. Makes me worry about using fossils or terminals > age 0.

	In these cases POWELL often performs better. At least it can find a point where the 
	rest of the gradient is near 0. QNEWT doesn't even get close!
*/

anyZeroTerminal=any_zero_terminal_branches(root);
anyZeroInternal=any_zero_internal_branches(root);
if (algorithm==QNEWT && anyZeroTerminal)
	{
	doGenericAlert("ZERO-LENGTH TERMINAL BRANCHES IN TREE (QNEWT will fail for PL and low smoothing values)\nTry algorithm=TN or powell!");
	*success=0;
	return;
	}
if (algorithm==TN && anyZeroTerminal)
	{
	doGenericAlert("ZERO-LENGTH TERMINAL BRANCHES IN TREE: TN will impose small nonzero bounds on parameters to overcome");
	}
if (anyZeroInternal)
	{
	doGenericAlert("ZERO-LENGTH BRANCHES IN TREE (you should run COLLAPSE first)");
	*success=0;
	return;
	}
if (penaltyType==1 && (anyZeroInternal || anyZeroTerminal))
	{
	doGenericAlert("Log penalty does not permit zero-length branches (try log(0) yourself)");
	*success=0;
	return;
	}
#if 0
if (penaltyType==1 && neighborPenalty==0 && method==PENLIKE && algorithm != POWELL)
	{
	doGenericAlert("At the moment can only do POWELL on log (anc/desc) penalty for PL!");
	*success=0;
	return;
	}
if (penaltyType==0 && neighborPenalty==1 && method==PENLIKE && algorithm != POWELL)
	{
	doGenericAlert("At the moment can only do POWELL on additive (neighbor) penalty for PL!");
	*success=0;
	return;
	}
#endif
if (node_tomy(root) > 2)
	doGenericAlert("ROOT IS A BASAL POLYCHOTOMY (is the tree UNROOTED?)");

switch (warnEstRoot(root))
	{
	case 1:
		doGenericAlert("You are trying to estimate the age of the root\nbut there is probably insufficient information\n(Try using FIXAGE or enforcing time constraints)\n...bailing on search!");
		*success=0;
		return;		// this is the only case where we bail
	case 2:
		if (verbose >0)
			doGenericAlert("You are trying to estimate the age of the root\nbut with the given constraints it is possible that a range of solutions exist");
	case 0:			// everything OK
		;
	}


/***** ..... Initial header output....see this file for details */




#include "ObjFuncHeader.h"




/****************** Allocate some arrays ************************/

if (!pTime  && !D)  /* only do this the first time...otherwise reuse these arrays...*/
           pTime=allocateTimeArray(root,method,nRates,&tvar,&nvar,&D); 
					/* This sets up gNVav, allocates pTIME, and calcs number of parameters */
zeroEstRate(root);	/* just zero out this value at each node */

lp = vector(1,nvar);
Fail_Flag = (int *)myMalloc((NUM_MULT_TRIES+1)*sizeof (int)); /* make this a 1-off array */
Mult_Like = vector(1,NUM_MULT_TRIES);
Z_Dist = vector(1,NUM_MULT_TRIES);
Rate_hat = vector(1,NUM_MULT_TRIES); 
Rate_hat2 = vector(1,NUM_MULT_TRIES);
p1 = vector(1,NUM_MULT_TRIES);
p2 = vector(1,NUM_MULT_TRIES);
storeParm = matrix(1,NUM_MULT_TRIES,1,nvar);
storeParm_init = matrix(1,NUM_MULT_TRIES,1,nvar);
storeGrad = matrix(1,NUM_MULT_TRIES,1,nvar);
approxGrad=vector(1,nvar);
storeGradInit = matrix(1,NUM_MULT_TRIES,1,nvar);
if (algorithm==TN)
	{
	gLOW=vector(0,nvar-1);
	gHIGH=vector(0,nvar-1);
	}


/****************** ...loop over multiple initial guesses ...   ******/

for (m=1;m<=NUM_TIME_GUESSES;m++) 
  {
  ++count;
  ii=(m-1); /* index into a 1-d array */

  if (verbose>0)
    printf("Starting optimization (random starting point %i)\n", m);
  maxorder(root);
  if (!setupFeasibleTimes(root))   /* Calculate an initial FEASIBLE guess at times and put on tree */
	{
	*success=0;			/* the constraints provided were probably invalid...bail*/
	printf ("...bailing...\n");
	return;
	}

  tree2pTimeArray(root,pTime);		/* Copies tree times to pTime array */
  initRate=treeLength(root)/treeDurLength(root); /* I use this instead of mean_rate(), because the latter is very sensitive to outliers that pop up on some branches that setUpFeasibleTimes initialized to be very short */
  minRate=rbp->minRateFactor*initRate;


  switch (method) /* Do miscellaneous set up stuff 
		    ** tvar = number of times
		    ** nvar = tvar + number of rate variables **/
	{
	case PENLIKE:// notice that we use the passed argument initRateArg here ONLY!	
			initRate=initRateArg;
			/*initTreeRates(root,gEstRoot,initRate);*/
			for (i=tvar+1;i<=nvar;i++)
				pTime[i]=initRate*(1+INIT_RATE_FRACTION*(0.5-myRand()));
  			minRate=rbp->minRateFactor*initRate;
			assignArrayRatesToLL2(root,pTime);
			check_if_exceed_factorial(root); /* make sure branch lengths
				aren't toooo long on this tree. AND make sure they are not between
				0 and 1, which would suggest that they are not in units of numbers
				of substitutions.*/
			break;

	case LaF:	
			pTime[gNVar]=initRate;
			check_if_exceed_factorial(root); 
			break;
	case LFLOCAL:
			for (i=gNVar-nRates+1;i<=gNVar;i++)	
				pTime[i]=initRate;
			check_if_exceed_factorial(root);
			break;
	case NP:	break;
	}

 /*save the initial point and gradient*/
 	for (i=1;i<=nvar;i++)
		storeParm_init[count][i]=pTime [i];
	if (method==PENLIKE || method==LaF)
		{
		gradient_array[method](pTime,D);/* get gradient at solution */
 		for (i=1;i<=nvar;i++)
			{
			storeGradInit[count][i]=D[i];
			}
#if 0 // useful when checking on gradient calculations 
		gradient_array[method](pTime,D);/* get gradient at solution */
		Dapprox(pTime,approxGrad, nvar, obj_func_array[method],0.00001);
		printf("Numerical gradient calculation prior to search:\n");
		for (i=1;i<=nvar;i++)
			printf("[%i]\t%e\t%e\n",i,D[i],approxGrad[i]);
#endif
		}


if (verbose>=2)
	printf(" Some initial conditions:\n  Root age = %f\n  Init rate per site = %e\n  MinRate per site = %e\n",root->time,initRate/gNumSites,minRate/gNumSites);

  if (algorithm==TN) /* set up the vectors containing lower and upper bounds for node times, used only by TN, also note
				0-length branches and fix these to have minimum non-zero durations */
	{
	setUpLowHighTN(Tree,nvar,tvar,minRate,gLOW,gHIGH);
	}

/*  } */

/*
 * Here is the call to the optimization routine
 */

	if (!ConstrOpt(
		    Tree,
		    Nex,
		    gisConstrained,		
		    gNVar /*nvar*/,	/* set above */
		    pTime,
		    obj_func_array[method],
		    gradient_array[method],
		    method,
		    algorithm,
		    penalty,
		    &maxLike,
		    &numIter,
		    &numRestartIter,
		    &numBarrierIter
		    ))
		{
		was_a_failure=1;
		Fail_Flag[m]=1; /* Set the failure code to 1: Failure in ConstrOpt or lower level routine */
		Tree->timesAvailable=0; /* confim that times have not been estimated */
		}

	else /* optimization returned OK */
		{
//		if (rbp->checkGradient && (method==PENLIKE || method==LaF)) /* ...but check the gradient if requested...*/
		if (rbp->checkGradient) /* ...but check the gradient if requested...*/
			{
// note that I *could* use the exact gradients here for some objective funcs; should allow user to choose
// ...this is to permit me to develop and debug the log penalty function, which as of yet does not have a gradient func 
// (except for neighbor penalty, which seems flaky)
#if 0 // useful when checking gradient calcs
		gradient_array[method](pTime,D);/* get gradient at solution */
		Dapprox(pTime,approxGrad, nvar, obj_func_array[method],0.00001);
		printf("Analytic and Numerical gradient calculation at solution:\n");
		for (i=1;i<=nvar;i++)
			printf("[%i]\t%e\t%e\n",i,D[i],approxGrad[i]);
#endif
			if ( method==NP || method == LFLOCAL/* ||  (method==PENLIKE && penaltyType==1) */) 
				{
				if (verbose>0) 
					printf("*** Analytical gradient not available: using numerical approximation to gradient in check gradient step ***\n    (may be inaccurate around 0, which will give wrong active[] sign)\n");
				Dapprox(pTime,D, nvar, obj_func_array[method],0.00001);
				}
			else
				{
				if (verbose>0)
					printf("*** Using analytical formula for gradient in check gradient step ***\n");
				gradient_array[method](pTime,D);/* get exact gradient at solution */
				}
			if (method != NP) // the true gradient is negative of the value we calculated for these methods
				{
				for (i=1;i<=nvar;i++)
					D[i]=-D[i];	
				}
			if(checkGradient(Tree,pTime,D,maxLike,rbp->ftol,verbose))
				{
				Fail_Flag[m]=0; 	/* keep track of whether this rep succeeded */ 
				Tree->timesAvailable=1; /* note that times have been constructed */
				Tree->method=method;	/* ..and how..*/
				}
			else
				Fail_Flag[m]=2; /* Failure code=2 means gradient was not 0 at proposed soln */
			}
		else
			{
			Fail_Flag[m]=0; 	/* keep track of whether this rep succeeded */ 
			Tree->timesAvailable=1; /* note that times have been constructed */
			Tree->method=method;	/* ..and how..*/
			}
		if (algorithm==TN)	/*...trap for estimated rates that run into the lower bound we impose */
			for (i=tvar+1;i<=nvar;i++)
				{
				if (pTime[i]==minRate)
					{
					doGenericAlert("Warning: An estimated rate crashed into the imposed lower bound on rates (see MINRATEFACTOR)\n\
You may be extrapolating too deep in tree for too low a smoothing value\n");
					break;
					}
				}
		}
	if ( method==NP || method==LFLOCAL /* || (method==PENLIKE && penaltyType==1) */) 
		Dapprox(pTime,D, nvar, obj_func_array[method],0.00001);
	else
		gradient_array[method](pTime,D);/* get gradient at solution */
	for (i=1;i<=gNVar /*nvar*/;i++)
		{
		storeParm[count][i]=pTime [i]; /* Save solutions and gradients*/
		storeGrad[count][i]=D[i];
		}
		
	/**** Peak diagnostic does brute force search around peak  ****
	 *
	 * peak_peek(objective,pTime,gNVar,0.01, 2); 
	 *
	 *************************/
 
/*
 * Save the objective function and some other stuff
 */


  switch (method) 
	{

	case PENLIKE:
		Mult_Like[m]=-maxLike;
		break;
	case LaF:
   		Mult_Like[m]=-maxLike;
   		Rate_hat[m]=pTime[gNVar];
		break;
	case LFLOCAL:
   		Mult_Like[m]=-maxLike;
   /** FIX 		Rate_hat[m]=pTime[gNVar]; **/
		break;
	case NP: 
   		Mult_Like[m]=maxLike;
//		Mult_Like[m]-=1.0;	/* Corrects for the amount added in obj func */
		Mult_Like[m]*=(-1);	 /* make this temporarily negative so we 
			can look for the maximum value across runs for any objective function*/
	}

  if (verbose > 0 && algorithm !=TN)
	{
  	printf("Optimization replicates used in first pass...%i\n",numIter);
  	printf("Optimization replicates used in LAST restart...%i\n",numRestartIter);
	}

   } /* end m (loop over multiple initial starting points)*/
  

  for (i=1,max = -1e100,maxi=1;i<=NUM_MULT_TRIES;i++)
	if (Mult_Like[i]>max)
		max=Mult_Like[i],maxi=i; /* finds max likelihood in array */

   Tree->obj=Mult_Like[maxi];

/* Put the times corresponding to the best estimate 
    back onto the tree data structure,  and do some other stuff with it. */
 
    for(j=1;j<=gNVar;j++) lp[j]=storeParm[maxi][j];


    if (verbose>0)
       printf("\nUsing optimization from starting point %i as best estimate\n", maxi);
    pTimeArray2tree(root,lp); 

/*
 * Print out results specific to chosen method
 */

    switch (method)
    {
    case PENLIKE:
	{
   	assignArrayRatesToLL2(root,lp); 
	fObj=objPenLike(lp);
	chiSq = LFuncons(root);    /* Check the Chi-sq test on this best tree */
	chiSqDF=numBranches(root)-(numIntNodes(root)-1)-1;
			/* df are number of branches (the data) minus the number of 
			estimated parameters: there are number of interior node times
			(-1 for the root), and one rate parameter */
	if (verbose>=1)
	  {
  	  printf("\nGoodness of fit test for soln. %i (best): Chi squared = %6f (df=%i)\n",maxi, chiSq,chiSqDF);
	  printf("Objective function value:%f\n",fObj);
	  printf("Likelihood portion of objective function=%f\n",gLike);
	  printf("Penalty portion (divided by smoothing factor):%f\n",(fObj-gLike)/gSmoothing);
	  }

	pTimeArray2tree(root,lp);
    	assignArrayRatesToLL2(root,lp); 
	
	break;
	}

    case LaF:
	b=0.0;
	c=0.0;
	chiSq = LFchiSq(root, pTime[gNVar]);    /* Check the Chi-sq test on this best tree */
	chiSqDF=numBranches(root)-(numIntNodes(root)-1)-1;
			/* df are number of branches (the data) minus the number of 
			estimated parameters: there are number of interior node times
			(-1 for the root), and one rate parameter */
	if (verbose>=1)
	  {
  	  printf("Test of molecular clock on soln. %i (best): Chi squared = %6f (df=%i)\n\n",maxi, chiSq,chiSqDF);
	  }
	rateType=1;
	Tree->estRate=lp[gNVar];	/* save the rate */
/**/ 	assignArrayRatesToLL_LF(Tree->root,Tree->estRate); /* MODIFY THIS IN LOCAL CLOCK MODEL !! */
	break;
    case LFLOCAL:
	rateType=1;
/**/ 	assignArrayRatesToLL_LFLOCAL(Tree->root,lp); /* MODIFY THIS IN LOCAL CLOCK MODEL !! */
	break;
    case NP:
	{
	set_est_rates(root,0.0,0.0,1); /* sets these up in case needed for ratograms */
	break;
	}
    }



/*
 * Output the parameter estimates
 */


if (verbose >= 2)
  {

  printf("\nStarting points for searches:\n\n\t");
  for (i=1;i<=NUM_MULT_TRIES;i++) printf("%8i  ",i);printf("\n\n");
  for (i=1;i<=nvar;i++)
	{
	if (i==1)printf("Times\n");
	if (i==tvar+1)printf("\nRates (substitutions per site per unit time)\n");
	printf("p[%2i]\t",i);
	for (j=1;j<=NUM_MULT_TRIES;j++)
		printf("% 8.3f  ",storeParm_init[j][i]/gNumSites);
	printf("\n");
	}

  }

if (verbose > 0 )
  {
  printf("\nParameter estimates:\n\n\t");
  for (i=1;i<=NUM_MULT_TRIES;i++) printf("%8i  ",i);printf("\n\n");
  for (i=1;i<= gNVar /* nvar*/;i++)
	{
	if (i==1)printf("Times\n");
	if (i==tvar+1)printf("\nRates (substitutions per site per unit time)\n");
	printf("p[%2i]\t",i);
	for (j=1;j<=NUM_MULT_TRIES;j++)
		printf("% 8.3g  ",storeParm[j][i]/gNumSites);
	printf("\n");
	}
  printf("\n\t");
  for (i=1;i<=NUM_MULT_TRIES;i++)
	{
	switch (Fail_Flag[i])
		{
		case 0: printf("  PASSED  ");break;
		case 1: printf("FAILED(1) ");break;
		case 2: printf("FAILED(2) ");break;
		} 
	}	
  printf("\n");
  printf("__________________\nResult codes:\nPASSED = OK\nFAILED(1) = Optimization routine returned error\nFAILED(2) = Solution's gradient was not 0\n");
  printf("__________________\n\n");
  printf("\nObj->\t");
  for (i=1;i<=NUM_MULT_TRIES;i++)
	{ 
	if (method==NP)
		printf("%+8.3g  ",-Mult_Like[i]);
	else
		printf("%+8.3f  ",Mult_Like[i]);
	}	
  printf("\n");
  }


/*
 *  Show the gradient if requested
 */


if (verbose > 0 && gNexDataPtr->RateBlockParms.showGradient)
	if (method==PENLIKE || method == LaF)
		{
		printf("\nGradient at starting points:\n\n\t");
  		for (i=1;i<=NUM_MULT_TRIES;i++) printf("%8i  ",i);printf("\n\n");
  		for (i=1;i<=nvar;i++)
			{
			if (i==1)printf("Times\n");
			if (i==tvar+1)printf("\nRates\n");
			printf("p[%2i]\t",i);
			for (j=1;j<=NUM_MULT_TRIES;j++)
				printf("% 8.3g  ",storeGradInit[j][i]);
			printf("\n");
			}
		printf("\n\n2-Norm = ");
		for (j=1;j<=NUM_MULT_TRIES;j++)
			printf("%8.3g  ",norm(storeGradInit[j],1,nvar));
		printf("\n");

  		printf("\nGradient at solutions:\n\n\t");
  		for (i=1;i<=NUM_MULT_TRIES;i++) printf("%8i  ",i);printf("\n\n");
  		for (i=1;i<=nvar;i++)
			{
			if (i==1)printf("Times\n");
			if (i==tvar+1)printf("\nRates\n");
			printf("p[%2i]\t",i);
			for (j=1;j<=NUM_MULT_TRIES;j++)
				printf("% 8.3g  ",storeGrad[j][i]);
			printf("\n");
			}
		printf("\n\n2-Norm = ");
		for (j=1;j<=NUM_MULT_TRIES;j++)
			printf("%8.3g  ",norm(storeGrad[j],1,nvar));
		printf("\n");
		}
/*
 * Show the convergence rate report if requested
 */


if (gNexDataPtr->RateBlockParms.showConvergence )
	{
	printf("Powell Convergence Diagnostics\nObjective:\n");
	while (hasElements(gFStack))
		printf("F=%e\n",popD(gFStack));
	printf("Convergence Diagnostics\nParameters:\n");
	while (hasElements(gPStack))
		printf("Norm=%e\n",popD(gPStack));
	}

if (NUM_MULT_TRIES>1 && verbose>=2)
	{
	printf("CHRONOGRAMS FOR  MULTIPLE SOLUTIONS:\n");
	printf("\nbegin trees;\n");
	for (i=1;i<=NUM_MULT_TRIES;i++)
		{
		if (i==maxi)
			printf("Tree  BEST%i=",i);
		else
			printf("Tree  8s%i=",i);
		for(j=1;j<=tvar;j++) lp[j]=storeParm[i][j];
		gIndex=1;
		gRoot=root;
		pTimeArray2tree(root,lp);
		make_parens(root,1);
		printf(";\n");
		}
	printf("end;\n");
	for (i=1;i<=NUM_MULT_TRIES;i++)
		{
		if (i==maxi)
			printf("Tree  BEST%i\n",i);
		else
			printf("Tree  r8s%i\n",i);
		for(j=1;j<=tvar;j++) lp[j]=storeParm[i][j];
		gIndex=1;
		gRoot=root;
		pTimeArray2tree(root,lp);
		DrawTree(root,2,80);
		}
/* MODIFIED 01/21/01 ** CAREFUL HERE ! */ 
/* NOW PUT THE OPTIMAL TIMES BACK ON THE TREE !! */
   	 for(j=1;j<=tvar;j++) lp[j]=storeParm[maxi][j];
   	 gIndex=1;
   	 gRoot=root;
 /*  	 pTimeArray2tree(root,lp);
    	assignArrayRatesToLL2T(root,lp); */
	}

/*
 * Free up everything
 */

free_vector(D,1,nvar);
free_vector(pTime,1,nvar);
free_vector(lp,1,nvar);
free_vector(Z_Dist,1,NUM_MULT_TRIES);
free_vector(Mult_Like,1,NUM_MULT_TRIES);
free_vector(Rate_hat,1,NUM_MULT_TRIES); 
free_vector(Rate_hat2,1,NUM_MULT_TRIES);
free_vector(p1,1,NUM_MULT_TRIES);
free_vector(p2,1,NUM_MULT_TRIES);

free_matrix(storeParm,1,NUM_MULT_TRIES,1,nvar);
free_matrix(storeParm_init,1,NUM_MULT_TRIES,1,nvar);
free_matrix(storeGrad,1,NUM_MULT_TRIES,1,nvar);
myFree(Fail_Flag);
if (algorithm==TN)
	{
	free_vector(gLOW,0,nvar-1);
	free_vector(gHIGH,0,nvar-1);
	}
*success=!was_a_failure;

return;
}
/************************************************************/
void ObjFuncDefaults(void) /* hopefully deprecated now */
{
gPowellTol=0.0000001;
gbarrierTol=0.0001;
gmaxPowellIter=500;
gmaxBarrierIter=10;
ginitBarrierFactor=.25;
gbarrierMultiplier=0.10;
glinminOffset=0.05;
gcontractFactor=0.1;
gmaxContractIter=10;

return;
}


/************************************************************/

int perturb_p(double p[], int n, double perturb_factor)

/* PERTURBATION OF p[] VECTOR 

For each component of a n-dimensional point,  perturbs it by an amount of
+ or -perturb_factor; checks to see if this new point is feasible,  and does
the same for all other components.  If any component change causes infeasibility, 
then the original component is restored.  IF NO CHANGE IN ANY COMPONENT OCCURS
then an error return is passed. 'errcount' records the number of components that
could not be changed.


*/


{
int i,j,binary, errcount=0;
extern int isFeasible,gIndex,gPowellTrace;
double ps, r,x;
for (i=1;i<=n;i++)
	{
	ps=p[i];
	r=2*(myRand()-0.5);
/* printf("r=%f\n",r);*/
	p[i]*=(1+r*perturb_factor);
	gIndex=1;
	pTimeArray2tree(gRoot,p);
	isFeasible=1;
	check_feasible(gRoot);
	if (!isFeasible)
	    {
	    p[i]=ps;
	    ++errcount;
	    if (gPowellTrace)
		debug_check_feasible(gRoot);
	    }
	}


isFeasible=1;
if (errcount==n)
    return 0; /* error return if all components stayed the same !*/
else
    return 1;		
}

/************************************************************/

int same_points(double p1[], double p2[], int n, double tolerance)

/* Tests whether any coordinates in two vectors differ by more than fractional
    tolerance;
 * if so,  this returns 0; otherwise the points are the "same" and returns 1;
 */

{
int ix;
for (ix=1;ix<=n;ix++)
    {
    if (fabs(p1[n]+p2[n])>0.01) /* ignore if ages are too close to zero (roundoff) */
	if (2*fabs(p1[n]-p2[n])/(p1[n]+p2[n]) > tolerance)
	    return 0;
    }
return 1;
}

void peak_peek(objfunc objective,double p[],int nvar,double sizeFactor, int grid)

/* Evaluates the objective function on lattice neighborhood around the point p.
   Tests whether any of the points on the lattice have a better score than the original point
   The lattice has a dimension for each parameter in p, given by nvar. It has a width
   given by p[k]*sizefactor in the direction of parameter k.  It has a number of lattice
   points in each dimension given by 'grid'.
*/

{
int next_ix(int ix[],int n,int top);
double *pCopy,min=+1e100,gg,obj;
int i;
int *ix;
/*printf("Center of lattice:\n");
for (i=1;i<=nvar;i++)
	  {
	  printf("p[%2i]\t%8.3f\n",i,p[i]);
	  }
*/
if (grid==1)
	{
	printf("Grid must be >1!\n");
	return;
	}
pCopy=vector(1,nvar);
ix=(int *)myMalloc((nvar+1)*sizeof(int));
gg=(grid+1)/2.0;
for (i=1;i<=nvar;i++)
	{
	ix[i]=1;
	pCopy[i]=p[i];
	}
printf("\nPeeking at the Peak\n\n");
do
	{
/*	for (i=1;i<=nvar;i++)
		printf("%1i",ix[i]);
	printf("\n");
*/
	for (i=1;i<=nvar;i++)
		{
		pCopy[i]=p[i]*( 1+2*sizeFactor*((float)(ix[i]-1)/(grid-1)-0.5));
		}
	
	obj=objective(pCopy);
	if (obj<min)min=obj;

/*  	for (i=1;i<=nvar;i++)
	  {
	  printf("p[%2i]\t%8.3f\n",i,pCopy[i]);
	  }
*/
	printf("Peek/Objective = %f\tbest=%f\n",obj,min);
	}
	while(next_ix(ix,nvar,grid));
printf("Minimum objective value on lattice = %f\n",min);
return;
}

int next_ix(int ix[],int n,int top)
{
int k=1;
while (k<=n)
  {
  if (ix[k]<=top-1)
	{
	++ix[k];
	return 1;
	}
  else
	{
	ix[k]=1;
	++k;
	}
  }
return 0;
}
/*******************************************************************************/
static double * allocateTimeArray(NODETYPE * root, int method,int nRates, int *tvar, int *nvar, double **D)

/* ALLOCATION OF PARAMETER  ARRAY (...and derivative array, D)!   */

{
extern int gNVar;
NODETYPE *child;
double *p;

if (root)
	{
	*tvar=numFreeNodes(root);		/* number of time parameters */
	if (method==PENLIKE)
		*nvar= *tvar+numBranches(root);	/* add a rate parameter for every branch in the tree */
	if (method==PENLIKET)
		*nvar= *tvar+numBranches(root)+1; /* every node gets a parameter, incl. root */
	if (method==LaF)
		*nvar= *tvar+1;
	if (method==LFLOCAL)
		*nvar= *tvar+nRates;
	if (method==GAMMA)
		*nvar= *tvar+2;		/* add two rate parameters */
	if (method==NP)
		*nvar= *tvar;
	gNVar= *nvar;
	}
else
	return NULL; /*error*/

p=vector(1,*nvar);
*D=vector(1,*nvar);
if (p)
	return p;
else
	fatal("Couldn't allocate solution array p[]\n");
}
/***********************************************************************************/

static void setUpLowHighTN(TREE t,int nvar, int tvar, double minRate, double LOW[],double HIGH[])

/* Sets up gLOW and gHIGH arrays for TN bound constraint algorithm; recurses through tree and scans
	for min/max constraints.
   Also performs special handling on any 0-length terminal branches: sets up a min duration, so that TN
	does not get mucked up
   Also puts a minimum value on rate parameters. Under PL and low smoothing these go to 0 otherwise and fail to converge
*/

{
#define LARGE_VAL	1e20
	extern struct NexDataType *gNexDataPtr;
	double minDur;
	int i;
	minDur=t->root->time * gNexDataPtr->RateBlockParms.minDurFactor;	
	for (i=0;i<tvar;i++) /* set up default values, including those for both times and rates */
		{
		LOW[i]=0.0;
		HIGH[i]=LARGE_VAL;
		}
	for (i=tvar;i<nvar;i++) /* set up default values, including those for both times and rates */
		{
		LOW[i]=minRate;
		HIGH[i]=LARGE_VAL;
		}
	gIndex=0; /* unlike elsewhere these are 0-offset arrays for FORTRAN calls in TN */
	setUpLowHigh_helper(t->root,LOW,HIGH,minDur);
	return;	
}
static void setUpLowHigh_helper(NODETYPE *node,double LOW[],double HIGH[],double minDur)

/* NB! DOESN'T YET WORK WITH NON-EXTANT TERMINALS */

{
	NODETYPE *child;
	if (isFree(node))
		{
		if (node->nodeIsConstrainedMin) 
			{
			if (node->minAge == minDur) /* resets to zero in case it has been set on a previous CV run perhaps (because of pruning taxa, which ought not to be persistent) */
				node->minAge=0.0;
			else
				LOW[gIndex]=node->minAge;
			}
		if (node->nodeIsConstrainedMax) 
			{
			HIGH[gIndex]=node->maxAge;
			}
		if (children_tips_are_zeroes(node)) /* only important case is when all child branches are terminal and 0; then the node needs to have a minimum slightly above 0 */
			{
			/* if the node is already constrained to a minimum, assume that that min is larger than minDur and do nothing*/
			if (!node->nodeIsConstrainedMin) 
					{
					node->nodeIsConstrainedMin=1;
					node->minAge=minDur;
					LOW[gIndex]=minDur;
					}
			}			
		++gIndex;
		}
	if (isTip(node))
			return;
	child=node->firstdesc;
	SIBLOOP(child)
			{
			setUpLowHigh_helper(child,LOW,HIGH,minDur);	
			}
	return;	
}
static int children_tips_are_zeroes(NODETYPE * parent)
{
/* return 1 if all children of node are tips AND all of child branches have 0-length */
NODETYPE *  n;
if (isTip(parent))return 0;
n=parent->firstdesc;
do
	{
	if (!isTip(n) ||  n->length != 0.0) return 0;
	n=n->sib;
	} while (n);
return 1;

}
/***********************************************************************************/

static void tree2pTimeArray(NODETYPE *node,double pTime[])

/* modern code 10.29.00 */

{
	gIndex=1;
	tree2pTimeArray_helper(node,pTime);
	return;	
}
static void tree2pTimeArray_helper(NODETYPE *node,double pTime[])

/* modern code 10.29.00 */

{
	NODETYPE *child;
	if (isFree(node)) 
		{
		pTime[gIndex++]=node->time;
		}
	if (isTip(node))
		return;
	child=node->firstdesc;
	SIBLOOP(child)
			{
			tree2pTimeArray_helper(child,pTime);
			}
	return;	
}
/* modern code 10.7.99 */

/***********************************************************************************/
int checkGradientSimple(double solution[],double gradient[],double Obj,double ftol,int nParm)

/* Examines optimality conditions for linear inequality bound constraints,
	and returns 1 if conditions are satisfied. See Gill et al., for conditions.

   solution, gradient, and Obj are all evaluated at the parameter estimate.

*/

{
extern struct NexDataType *gNexDataPtr;
extern int gNVar;
int N,i,j,numConstr=0,success=1, *active;
double nm,crit,active_eps;

active=(int*)myMalloc((gNVar+1)*sizeof(int));

for (i=1;i<=nParm;i++) active[i]=0;



/*
   Check if gradient for all non-active (free) parameters is 0 (condition J2) 
   Termination criterion of Gill et al. p. 306:U3 
   This includes all rate and time parameters.
*/

nm=norm_not_active(gradient,active,1,nParm);
crit= pow(ftol,0.333)*(1+fabs(Obj)); /* need to consider rate variables too? */
printf("Checking gradient norm: norm=%g critical value not to exceed = %g ... ",nm, crit);
if (nm > crit)
	{ printf ("FAILED\n"); return 0; }
else
	{ printf ("PASSED\n"); return 1; }
}
/***********************************************************************************/
int checkGradient(TREE t,double solution[],double gradient[],double Obj,double ftol,int verbose)

/* Examines optimality conditions for linear inequality bound constraints,
	and returns 1 if conditions are satisfied. See Gill et al., for conditions.

   solution, gradient, and Obj are all evaluated at the parameter estimate.

*/

{
extern struct NexDataType *gNexDataPtr;
extern int gNVar;
NODETYPE ** nodeArray,*node,*n;
int N,i,j,numConstr=0,success=1, *active;
double nm,crit,active_eps;
active_eps=gNexDataPtr->RateBlockParms.activeEpsilon * t->root->time;
N=numFreeNodes(t->root);
nodeArray=newNodeArray(t);
active=(int*)myMalloc((gNVar+1)*sizeof(int));

set_active(solution, nodeArray,N,active,active_eps);
/*
   Check to see if all *time* parameters satisfy bounds (Gill et al. condition J1, p. 77). 
   Not necessary. All the feasible checking makes sure this condition holds!

*/


/*
   Check if gradient for all non-active (free) parameters is 0 (condition J2) 
   Termination criterion of Gill et al. p. 306:U3 
   This includes all rate and time parameters.
*/

nm=norm_not_active(gradient,active,1,gNVar);
crit= pow(ftol,0.333)*(1+fabs(Obj)); /* need to consider rate variables too? */
if (nm > crit)
	{
	success=0;
	}

/* check if gradient for active constraints are either positive or negative
(Gill et al., condition J3). Notice that within this routine we assume the sign of the gradients is correct
with respect to the real objective function. In other words, since we usually negate everything for LaF and PL
this routine requires the input gradient to be corrected back to its true sign. 

I ignore parameters whose gradients are nearly zero in this test. It seems that roundoff error in these
cases may often produce the wrong sign. Below I set a fairly arbitrary criterion to determine whether a
gradient is big enough at an active parameter to worry about.

   At the moment we ignore possible bounds on the rate parameters. These only come into play in the
   pathological case of zero-length terminals. Consider treating this later...
*/

if (verbose>0)
	{
	printf("...checking gradient of the solution...\n\n");
	printf("\tNorm for free parameters in gradient (%f) should be less than cutoff (%f)\n",nm,crit);
	printf("\tChecking active constraints for node times only (activeEpsilon=%f and window [actEps*rootage]=%f)\n",gNexDataPtr->RateBlockParms.activeEpsilon,active_eps);
	printf("\tParam\tEstimate\tGradient\tActive*\tTaxon\n");
	for (i=1;i<=N;i++)
			{
			n=nodeArray[i];
			printf("\t[%i]\t%f\t%f\t%i\t%s\n",i,solution[i],gradient[i],active[i],n->taxon_name);
			}
	printf("\n\tGradients for rate parameters (if any)\n");
	printf("\tParam\tEstimate\tGradient\n");
	for (i=N+1;i<=gNVar;i++)
			printf("\t[%i]\t%f\t%f\n",i,solution[i],gradient[i]);
	// we aren't treating the rate parameters as active or not, so don't print out active[] for them
	
	printf("------------------------------------------------------------------------------\n");
	printf("*Key\n");
	printf("\t+1 = maximum age constraint is reached (gradient may not be 0)\n");
	printf("\t-1 = minimum age constraint is reached (gradient may not be 0)\n");
	printf("\t 0 = no constraint present or constraint not reached (gradient should be 0)\n\n");
	printf("\t (note that small gradient values (< |%f|) are not examined with respect to active constraints to avoid spurious roundoff issues)\n\n",crit*0.0001);
	}
for (i=1;i<=gNVar;i++)
	{
	if (active[i] != 0)
	    if (fabs(gradient[i]) > 0.0001*crit) // we ignore this test when the gradient is approx 0 anyway 
		{
		if ((active[i]== +1 && gradient[i]<0) ||
		    (active[i]== -1 && gradient[i]>0))
			{
			if (verbose>0)
				printf("Active parameter [%i] gradient has wrong sign at solution\n",i);
			success=0;
			}
		}
	}
if (verbose>0)
	{
	if (success)
		printf("*** Gradient check passed ***\n");
	else
		printf("*** Gradient check FAILED ***\n");
	}
myFree(nodeArray);
myFree(active);
return success;
}
/***********************************************************************************/


int set_active(double * solution, NODETYPE **nodeArray,int nNodes,int active[],double active_eps)

/* Checks all node parameters to see which are on boundaries and sets these as active constraints. Leaves others alone.
   If parameter is close to upper boundary, active[]=+1
   If parameter is close to lower boundary, active[]=-1
   If parameter is interior to boundaries,  active[]= 0
   Closeness is decided by 'active_eps'

   Note that if the user sets a min and max constraint to be very close to each other, then it might be possible
   for both constraints to be "active", which causes a problematic indeterminacy. In such a case, we issue a warning
   telling the user to either make activeEpsilon smaller or the constraints farther appart, and we set active[]->0. This
   may cause the gradient norm check to fail, but so be it...the user ought to make changes. If it doesn't fail, no harm
   done.
*/

{
extern int gNVar;
int i,numActive=0;
NODETYPE *node;
for (i=1;i<=gNVar;i++) 
	active[i]=0; // the active array covers all the parameters for future work
for (i=1;i<=nNodes;i++)
	{
	node=nodeArray[i];
	if (node->nodeIsConstrainedMax)
		if (fabs(solution[i]-node->maxAge)<active_eps)
			{
			active[i]=+1;
			++numActive;
			}
	if (node->nodeIsConstrainedMin)
		if (fabs(solution[i]-node->minAge)<active_eps)
			{
			active[i]=-1;
			++numActive;
			}
	if (node->nodeIsConstrainedMin && node->nodeIsConstrainedMax)
		if (fabs(node->maxAge-node->minAge)<=2*active_eps)
			{
			doGenericAlert("Check gradient results are problematic. See below...");	
			printf("Check gradient problem at node %s\n",node->taxon_name);
			printf("Node's min and max constraints are too close together for the current value of activeEpsilon\n");
			printf("The direct solution is to either make the constraints further apart or make activeEpsilon smaller\n");
			printf("An even better approach might be to use FIXAGE instead of constraints if they are so close!\n");
			active[i]=0; // remains this way by default
			}
	}

return numActive;
}

/***********************************************************************************/
NODETYPE** newAllNodeArray(TREE t)

{
long  N;
NODETYPE ** nodeArray;
N=numNodes(t->root);
nodeArray=(NODETYPE**)myMalloc((N+1)*sizeof(NODETYPE*)); /* 1-offset array */
gIndex=1;
newAllNodeArray_helper(t->root,nodeArray);
return nodeArray;
}
NODETYPE** newNodeArray(TREE t)
/* this one just makes a node array for free nodes */
{
long  N;
NODETYPE ** nodeArray;
N=numFreeNodes(t->root);
nodeArray=(NODETYPE**)myMalloc((N+1)*sizeof(NODETYPE*)); /* 1-offset array */
gIndex=1;
newNodeArray_helper(t->root,nodeArray);
return nodeArray;
}

static void newAllNodeArray_helper(NODETYPE *node,NODETYPE ** nodeArray)
{
	NODETYPE *child;
	nodeArray[gIndex++]=node;
	child=node->firstdesc;
	SIBLOOP(child)
			{
			newAllNodeArray_helper(child,nodeArray);
			}
	return;	
}
static void newNodeArray_helper(NODETYPE *node,NODETYPE ** nodeArray)
{
	NODETYPE *child;
	if (isFree(node)) 
		{
		nodeArray[gIndex++]=node;
		}
	if (isTip(node))
		return;
	child=node->firstdesc;
	SIBLOOP(child)
			{
			newNodeArray_helper(child,nodeArray);
			}
	return;	
}


/***********************************************************************************/


void pTimeArray2tree(NODETYPE *root,double lp[])
{
NODETYPE * child;
gIndex=1;
assignArrayTimesToLL_helper(root,lp);
return;
}

static void assignArrayTimesToLL_helper(NODETYPE * node, double lp[])
{
NODETYPE *child;
if (isFree(node))
    node->time=lp[gIndex++];
if (isTip(node))
	return;
child=node->firstdesc;
SIBLOOP(child)
    assignArrayTimesToLL_helper(child,lp);
return;
}

void Dapprox(double p[],double grad[],int n, double (* obj)(double p[]),double h)

/*
Finite difference approximation to the gradient using central difference approximation,

	f(x+h)-f(x-h) / 2h

Technical issue here. If the step size is too large, the estimated derivative is wrong because the
central difference approx is off. If it is too small, then especially near the optimum where the gradient
should be zero, we run into roundoff errors, dividing small difference by small differences.

At the moment this routine overrides the value of h passed to it! The current value is based on experiments with some sample data sets, but it may not be perfect for every data set or method/algorithm. The step length should be chosen according to
some more clever scheme...

NB. Some comments from Gill et al. on this stuff: don't be tempted to use this routine as a plug in for an optimization
that requires gradients. It is very difficult to get high precision near the gradient of zero.

*/

{
int i;
double f,f1,f2,dif,psave,p1,p2;
//h=0.00001;
for (i=1;i<=n;i++)
	{
	psave=p[i];
//	p1=psave-psave*h;   // scaling this way turns out to be a bad idea! Causes flip-flop of sign sometimes
	p1=psave-h;
	p[i]=p1;
	f1=obj(p);
//	p2=psave+psave*h;
	p2=psave+h;
	p[i]=p2;
	f2=obj(p);
	p[i]=psave;
//	dif=(f2-f1)/(2*psave*h);
	dif=(f2-f1)/(2*h);
	grad[i]=dif;
//	printf("***[%i] %16.12f %16.12f %16.12f %16.12f %16.12f %16.12f %f\n",i,p[1],p[2],p1,p2,f1,f2,dif);
	}
return;
}
