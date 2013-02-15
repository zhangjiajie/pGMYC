/**

Idea for another smoothing function: calculate a term for each node on the unrooted tree, the term
corresponding to the variance of the rates of all incident branches. This makes every node equivalent 
in some sense, but also may add too much weight to sister group rates being similar (why should they be
as similar as ancestor/descendant branches? they shouldn't)

Could also expand the "window" size to look at the variance among 2-off, etc., neighbors. Before when I
tried this approach, I don't think I explicitly calculated variances...

Started to work on this with function NeighborVariance...see below

**/


/* This module has a bunch of nasty stuff to actually implement the objective
function on a tree, among other things.  It also has routines to find an initial
feasible point. 

*/
#include "DistrFuncs.h"
#include "TreeUtils.h"
#include "TreeSim.h"
#include "TimeAlgorithms.h"
#include "memory.h"
#include "math.h"
#include "penalty.h"
#include "stdio.h"
#include "stdlib.h"
#include "ObjFunc.h"
#include "nexus.h"

#define SQR(x) 		((x)*(x))
#define ZERO(x)   	(fabs(x) < 0.0001)
#define MAX_FACTORIAL 	500			/* precompute values up to this point */
#define DEBUG 		0			/* 0-2 to display different levels of debugging info */
#define LARGE_VAL 	1e20
/* define LARGE_VAL carefully.  It must be larger than any likely value of the objective
function at the solution, but it must not be too close to HUGE_VAL, which throws
an exception on some machines */

static double NeighborSum(NODETYPE *n, int * numBranches);
static double NeighborVariance(NODETYPE *n);
static double penalizedRatesNeighbor(NODETYPE *n);
void derivRateNeighbor(NODETYPE * n, double p[], double grad[],int *ixPtr);

static double BranchLikeSumNegBinomial(double rate, long nSites,double alpha, double T, double k);
static void tree2pTimeArray(NODETYPE *node,double pTime[]);
static double recurseSldWin(NODETYPE* node);
static void assignArrayRatesToLL2_helper(NODETYPE *root,double lp[], int *ix);
static void initTreeRates_helper(NODETYPE * node, int *index,double rate);
static double recurseLangFitchLocal(NODETYPE *node, NODETYPE * itsAncestor, double p[]);
static double recursePenLike(NODETYPE *node, NODETYPE * itsAncestor);
static double recursePenLikeT(NODETYPE *node);
static double penalizedRates(NODETYPE *root);
static double recursePenalizeRates(NODETYPE *node, NODETYPE * itsAncestor);
static double penalizedRatesT(NODETYPE *root);
static double recursePenalizeRatesT(NODETYPE *node);
static double recursePenalizeRates2(NODETYPE *node, NODETYPE * itsAncestor);
static void derivTimeLF(NODETYPE * n, double p[], double grad[],int *ixPtr);

int		powellMode;
int		gNVar;
double 		gSmoothing,gFit,gLike;
double 		FactLookup[MAX_FACTORIAL+1];
double 		logFactLookup[MAX_FACTORIAL+1];
int 		gRootFlag, 
    		gFloatRoot;
int 		gVarMinFlag;


/************************************************************/
/************************************************************/

			/* Gradients */

/************************************************************/
/************************************************************/

/* Langley fitch */

/************************************************************/

void GradientLF(double p[], double grad[])

// To generalize this to the LFLOCAL model just write a function that operates on the subtree defined by the different rate models

{
extern int gNVar;
extern NODETYPE * gRoot;
int index=1;
double rate;
rate=p[gNVar];


/**!!! Following two calls may be too expensive during a search ***/
pTimeArray2tree(gRoot,p); 
assignArrayRatesToLL2(gRoot,p); 
/* do the partial Derivs wrt time parameters */
derivTimeLF(gRoot,p,grad,&index);
/* now do the partial Deriv wrt rate parameter */
grad[index]=- (treeLength(gRoot)/rate - treeDurLength(gRoot)); /* minimize it, stupid */

//printf("---index=%i grad=%f rate=%f treeLength=%f TreeDur=%f\n",index,grad[index],rate,treeLength(gRoot),treeDurLength(gRoot));

return;


}
static void derivTimeLF(NODETYPE * n, double p[], double grad[],int *ixPtr)

/* Calculates the derivatives of log likelihood with respect to the time variables
   for the LF  method */



{
extern int gNVar;
NODETYPE * child;
double g=0.0,rate,gp;
rate=p[gNVar];
if (isFree(n))
	{
	if (!isRoot(n))
		{
		if (n->length==0.0)
			g=rate;
		else
			g = -n->length/(n->anc->time-n->time)+rate;

		}
	if (!isTip(n))
		{
		child=n->firstdesc;
		SIBLOOP(child)
			{
			if (child->length==0.0)
				gp=-rate;
			else
				gp=(child->length)/(n->time-child->time)-rate;
#if 0
printf("$$:%f %f %f %f %f %f\n",child->length,n->time,child->time,rate,g,gp);
#endif
			g+=gp;
			}
		}
	grad[*ixPtr]=-g;
	++(*ixPtr);
	}
child=n->firstdesc;
SIBLOOP(child)
			{
			derivTimeLF(child,p,grad,ixPtr);
			}
return;
}

/************************************************************/

/* Penalized Likelihood */

/************************************************************/

void GradientPL(double p[], double grad[])

{
extern struct NexDataType *gNexDataPtr;
extern NODETYPE * gRoot;
int index=1;



/**!!! Following two calls may be too expensive during a search ***/
pTimeArray2tree(gRoot,p); 
assignArrayRatesToLL2(gRoot,p); 



derivTime(gRoot,p,grad,&index);

if (gNexDataPtr->RateBlockParms.PenaltyType==0)
	derivRate(gRoot,p,grad,&index);
else
	derivRateLog(gRoot,p,grad,&index);
return;
}

void derivTime(NODETYPE * n, double p[], double grad[],int *ixPtr)

/* Calculates the derivatives of log likelihood with respect to the time variables
   for the PL method */

/* I think I can cut the time in ca. 1/2 for these gradient routines by precomputing the 
differences between node times and ancestor times (and rates) and storing them on trees.
Looks like we often calculate these things twice as part of the child loops in these two 
routines.*/


{
NODETYPE * child;
double g=0.0;
if (isFree(n))
	{
	if (!isRoot(n))
		{
		if (n->length ==0.0)
			g=n->estRate;
		else
			g = -n->length/(n->anc->time-n->time)+n->estRate;
		}
	if (!isTip(n))
		{
		child=n->firstdesc;
		SIBLOOP(child)
			{
			if (child->length ==0.0)
				g-=child->estRate;
			else
				g+=(child->length)/(n->time-child->time)-child->estRate;
			}
		}
	grad[*ixPtr]=-g;
	++(*ixPtr);
	}
child=n->firstdesc;
SIBLOOP(child)
			{
			derivTime(child,p,grad,ixPtr);
			}
return;
}


void derivRate(NODETYPE * n, double p[], double grad[],int *ixPtr)

/* Calculates the derivatives of log likelihood with respect to the rate variables
   for the PL method*/


{
NODETYPE * child;
double g=0.0,meanr=0.0;
int tomy=0;
extern double gSmoothing;

if (!isRoot(n))
	{
	g=n->length/n->estRate-(n->anc->time-n->time); /* part due to likelihood */

	if (isRoot(n->anc)) /* node is immediate desc of root: special case */
		{
		child=n->anc->firstdesc;
		SIBLOOP(child)
			{
			++tomy;
			meanr+=child->estRate;
			}
		meanr/=tomy;
 		g+= (-2*gSmoothing)*(n->estRate-meanr)/tomy;

		child=n->firstdesc;
		SIBLOOP(child)
				g+= 2*gSmoothing*(child->estRate-n->estRate);
		}
	else
		{
		g+=(-2*gSmoothing)*(n->estRate-n->anc->estRate);
		if (!isTip(n))
			{
			child=n->firstdesc;
			SIBLOOP(child)
				g+=2*gSmoothing*(child->estRate-n->estRate);
			}
		}
	grad[*ixPtr]=-g;
	++(*ixPtr);
	}

child=n->firstdesc;
SIBLOOP(child)
	{
	derivRate(child,p,grad,ixPtr);
	}
return;
}
void derivRateLog(NODETYPE * n, double p[], double grad[],int *ixPtr)

/* Calculates the derivatives of log likelihood with respect to the rate variables
   for the PL method USING A LOG PENALTY ON THE RATES*/


{
NODETYPE * child;
double g=0.0,meanr=0.0,lognrate;
int tomy=0;
extern double gSmoothing;

if (!isRoot(n))
	{
	g=n->length/n->estRate-(n->anc->time-n->time); /* part due to likelihood */

	lognrate=log(n->estRate);
	if (isRoot(n->anc)) /* node is immediate desc of root: special case */
		{
		child=n->anc->firstdesc;
		SIBLOOP(child)
			{
			++tomy;
			meanr+=log(child->estRate);
			}
		meanr/=tomy;
 		g+= (-2*gSmoothing/n->estRate)*(lognrate-meanr)/tomy;

		child=n->firstdesc;
		SIBLOOP(child)
				g+= 2*gSmoothing*(log(child->estRate)-lognrate)/n->estRate;
		}
	else
		{
		g+=(-2*gSmoothing)*(lognrate-log(n->anc->estRate))/n->estRate;
		if (!isTip(n))
			{
			child=n->firstdesc;
			SIBLOOP(child)
				g+= 2*gSmoothing*(log(child->estRate)-lognrate)/n->estRate;
			}
		}
	grad[*ixPtr]=-g;
	++(*ixPtr);
	}

child=n->firstdesc;
SIBLOOP(child)
	{
	derivRateLog(child,p,grad,ixPtr);
	}
return;
}
/************************************************************/
void derivRateNeighbor(NODETYPE * n, double p[], double grad[],int *ixPtr)

/* Calculates the derivatives of log likelihood with respect to the rate variables
   for the PL method using log penalty and neighbor variance*/
   
/* Numerous experiments show problems at high smoothing values for TN routine here. The TN
	routine requires that the function to be minimized be bounded below. This is not true
	for this function, since log(r) goes to negative infinity as r goes to zero. QNEWT seems
	to work better at high smoothing values--but then it croaks at low smoothing values! */


{
NODETYPE * child, *anc;
double gradLike,g=0.0,meanr=0.0,logsum1,logsum2,nRate,ancRate;
int nbranch1,nbranch2;
extern double gSmoothing;

if (!isRoot(n))
	{
	anc=n->anc;
	nRate=n->estRate;
	logsum1=NeighborSum(anc,&nbranch1);
	g+= (2*log(nRate)/nRate-2*logsum1/(nbranch1*nRate))/nbranch1;
	
	if (!isTip(n))  // this is the case of an internal rate (which has two terms instead of just one for terminal rates)
		{
		logsum2=NeighborSum(n,&nbranch2);
		g+= (2*log(nRate)/nRate-2*logsum2/(nbranch2*nRate))/nbranch2;
		}
	gradLike=n->length/nRate-(n->anc->time-n->time); /* part due to likelihood */
	grad[*ixPtr]=-(gradLike-gSmoothing*g);  /* it's a minimization */
//printf ("GRAD[%i]: %e %e %i %i %f %f %f %f\n",*ixPtr,grad[*ixPtr],nRate,nbranch1,nbranch2,gradLike,g,logsum1,logsum2);
	++(*ixPtr);
	}

child=n->firstdesc;
SIBLOOP(child)
	{
	derivRateNeighbor(child,p,grad,ixPtr);
	}
return;
}
/************************************************************/

void assignArrayRatesToLL_LF(NODETYPE *node,double rate)

/* Assigns all nodes a single rate; ignores root rate */

{
NODETYPE *child;
node->estRate=rate;
child=node->firstdesc;
SIBLOOP(child)
    assignArrayRatesToLL_LF(child,rate);
return;
}
void assignArrayRatesToLL_LFLOCAL(NODETYPE *node,double p[])

/* Assigns all nodes rates according to local model; ignores root rate */

{
NODETYPE *child;
extern int gNVar;
int rateIndex;
rateIndex=gNVar-node->modelID;
node->estRate=p[rateIndex];
child=node->firstdesc;
SIBLOOP(child)
    assignArrayRatesToLL_LFLOCAL(child,p);
return;
}

    
void assignArrayRatesToLL2T(NODETYPE *root,double lp[])

/* includes root rate */

{
NODETYPE *child;
int index=numFreeNodes(root)+1;  /* set index to one after last time in array */
assignArrayRatesToLL2_helper(root,lp, &index);
return;
}

void assignArrayRatesToLL2(NODETYPE *root,double lp[])

/* ignores root rate */

{
NODETYPE *child;
int index=numFreeNodes(root)+1;
child=root->firstdesc;
SIBLOOP(child)
    assignArrayRatesToLL2_helper(child,lp, &index);
return;
}

static void assignArrayRatesToLL2_helper(NODETYPE * node, double lp[], int *index)
{
NODETYPE *child;
node->estRate=lp[(*index)++];
child=node->firstdesc;
SIBLOOP(child)
    assignArrayRatesToLL2_helper(child,lp,index);
return;
    
}

void initTreeRates(NODETYPE *root, int includeRootFlag,double rate)

/* Initialize all the branch's rates for the penalized likelihood method */

{
NODETYPE *child;
int index=numIntNodes(root);
if (!includeRootFlag)
	--index;	/* points to last time in pTime array */
++index;		 /* set index to one after last time in array */
child=root->firstdesc;
SIBLOOP(child)
    initTreeRates_helper(child,&index,rate);
return;
}

static void initTreeRates_helper(NODETYPE * node, int *index,double rate)
{
NODETYPE *child;
node->estRate=rate;
++(*index);
child=node->firstdesc;
SIBLOOP(child)
    initTreeRates_helper(child,index,rate);
return;
    
}

/*******************************************************/

int warnEstRoot(NODETYPE * root)
/*
	1.  By default the program tries to estimate all internal nodes including the root.
	2.  However, this is only possible under the following conditions:
		A.  At least one internal node time is fixed (t>0?) with setage command
		B.  Not all the tips are extant (or the same age)
		C.  Some node (other than root) has a maximum age constraint AND the search is constrained
			(NB.  This often constrains internal times to a range of values only)
	3.  If none of these conditions are met, the search will bail with a warning that
		further age information must be supplied, forcing the user, e.g., to set the
		root age to 1.0.
*/
{
if (isFree(root))	/* default set by Tree_Init, unless changed by setage command */
	{
	if 	(
		(numFreeNodes(root)<numIntNodes(root)) ||  /* if some internal nodes have been fixed... */
	 	(tipsDifferentAges(root)) ||
		(maxAgePresent(root))
		)
			{
//			doGenericAlert("You are trying to estimate the age of the root\nbut with the given constraints it is possible that a range of solutions exist");
			return 2;
			}
	else
		{
//		doGenericAlert("You are trying to estimate the age of the root\nbut there is probably insufficient information\n(Try using FIXAGE or enforcing time constraints)\n...bailing on search!");
		return 1;			/* none of conditions hold */
		}
	}
else
	return 0;	/* root is already set; no warning necessary */

}


/*******************************************************/
/******* ROUTINES FOR FINDING A FEASIBLE POINT *********/
/*******************************************************/


int setupFeasibleTimes(NODETYPE * root)

/* Set up some feasible times and stores them in tree. Present-day is time 0.*/

{
extern int isFeasible;
double Length, minTime=0.0,maxTime=0.0;
NODETYPE *child;

descMinAge(root,&minTime,&maxTime); /* get the LARGEST minimum and max age of all descendants */

if (isFree(root))
	{
	if (root->nodeIsConstrainedMax)   /* if there is a root max age constraint */
        	  root->time=minTime+(root->maxAge-minTime)*(0.02+myRand()*0.96); 
	else
		{
		if (minTime !=0 )
		  root->time=minTime*1.25;
		else	
		  {
		  if (maxAgePresent(root))
		      root->time=maxTime*1.25; /* but what if maxTime = 0? */
		  else
		      root->time=1.0; /* if no minages or maxages...shoultn't happen,
			currently precluded by 'warnEstRoot' */
		  }
		}
	}
child=root->firstdesc;
SIBLOOP(child)
    {
    if (!aFeasibleTime(child,root->time))
	return 0;		
    }	
isFeasible=1;	/* global */
return 1;
}
/*******************/
int aFeasibleTime(NODETYPE *node,double timeAnc)

/* Moves through the tree, checking the constraints, and sets times of each node so
that they are also feasible according to constraints.  Constraints include explicit
minage and maxage statements for internal nodes and times of leaf nodes.  */

{
double minTime=0.0,maxTime=0.0;		/* important that this be set to 0.0; see 'descMinAge' */
NODETYPE *child;
descMinAge(node,&minTime,&maxTime);/* this is the largest minimum age of this node and ALL descendants */
if (isFree(node))
  {
	
  if (node->nodeIsConstrainedMax)
     if (node->maxAge < timeAnc)
	timeAnc=node->maxAge;   /* the age of this node must be <= to its maxAge */


  node->time = timeAnc - (timeAnc-minTime)*(0.02+myRand()*0.96)/log(node->order+3);

/* ...the idea here is to assign a random time to node that is between timeAnc and timeAnc-minTime, but not equal to either 
   ...however, when minTime=0 this often makes the MRCA of a big basal clade way too recent, if the random number happens to call
      for it. Therefore I divide by a monotonic function of the node->order to try to correct.  This is all clunky
      but there is no clear solution.  We just want a series of trial points anyway, although the more diffuse they are, the better.
      If there is an immediate desc node with a mintime, then in that case, we should throw a uniform random number down,
      but tricky to identify, because such a constraint may stem from a shallower node. Note that the function log(node->order+3)
      is just to keep it from being less than 1 but not as large as node->order.*/

/*printf("feasible time=%f timeAnc=%f minTime=%f order=%i\n",node->time,timeAnc,minTime,node->order);*/
  }

if (!isTip(node))
 {
 child=node->firstdesc;
 SIBLOOP(child)
    {
    if (!aFeasibleTime(child,node->time))
	return 0;
    }
 }
return 1;
}
/*******************/
double nodeLowerBound(NODETYPE *node)

/* Finds the lower bound on a node's age. This is the LARGER of
	(1) the oldest FIXED node age among descendants, and
	(2) the oldest 'minimum age' constraint among descendants

   If the node itself is fixed, we return its age.
*/

{
double minTime=0.0,maxTime=0.0;
descMinAge(node,&minTime,&maxTime);
return minTime;
}
/*******************/
double nodeUpperBound(NODETYPE *node)

/* Finds the upper bound on a node's age. This is the SMALLER of
	(1) the youngest FIXED node age among ancestors back to the root, and
	(2) the youngest 'maximum age' constraint among ancestors back to root

   NB! If there is no upper bound, return a value of 1e20 (i.e., big)
*/

{
double maxTime=1e20;

for (;node;node=node->anc) /* go through all the ancestors incl. root */
		{
		if (!isFree(node)) /* FIXED */
			{	
			if (node->time<maxTime)
				maxTime=node->time;
			}
		else
			if (node->nodeIsConstrainedMax && node->maxAge < maxTime)
				maxTime = node->maxAge;

		} 
return maxTime;
}


void descMinAge(NODETYPE *node, double *curMin,double *curMax)

	/* finds the largest minimum age and max age constraint among all descendants of node
	INCLUDING THIS NODE!.
	NB.  Second parm must point to 0.0 on first call (and 3rd to a large val) */

{
	int copyindex;
	NODETYPE *child;

	if (!node) return;

	if (!isFree(node)/* || isTip(node) */)
		{
		if (node->time > *curMin)
			*curMin=node->time;  /* This does not prevent a desc age from being
						OLDER than this node */
		if (node->time > *curMax)
			*curMax=node->time;
		}

	else
		{
		if ((node->nodeIsConstrainedMin) && (node->minAge > *curMin))
			*curMin = node->minAge;
		if ((node->nodeIsConstrainedMax) && (node->maxAge > *curMax)) /* ??? */
			*curMax = node->maxAge;
		}
	if (!isTip(node))
		{
		child=node->firstdesc;
		SIBLOOP(child)
				{
				descMinAge(child,curMin,curMax);
				}
		}	
	return;
}
/********************************************************************/
/**************** OBJECTIVE FUNCTIONS *******************************/
/********************************************************************/
double objPenLike(double p[])

/* A penalized likelihood objective function, consisting of two terms:
	(1) a likelihood calculated via LF but with each branch having a different
		rate parameter
	(2) a penalty deducted from the likelihood term, comprising a smoothing
		factor mutliplied by squared differences between neighboring
		branch's rates

NOTES: Large values of the smoothing parameter lead to frequent problems of
	nonconvergence of Powell.  Often these can be proven by using the
	peak_peek function to show that there are neighboring points that
	are more optimal than the proposed solution.  Restarts and perturbations
	often do not help!  At this time I don't have a solution.  Make sure to
	do lots of time guesses (which now begin with randomly perturbed rate
	guesses too). On the other hand, reasonable values of smoothing seem to work.

*/

{
  extern struct NexDataType *gNexDataPtr;
  extern NODETYPE * gRoot;	    /* This global is declared when the whole 
					    algorithm is called */
  extern int gisConstrained;	       /* are we doing a constrained optimization? */
  extern int isFeasible;
  extern int gEstRoot,gRatesAreGamma;
  extern long gNumSites;
  double obj=0.0,Like=0.0,rt, K, k1, k2,PY;
  int i;
  NODETYPE *child;
  extern int gFirstDesc;
  static int f=1;

  if (DEBUG > 1)
    printf("--DEBUG-- Entry to ObjPenLike\n");
  if (gisConstrained)
    Like+=penalty(p); 
  
  if (f) {setupLogFactLookup(); f=0;}       /* temporary kludge to call this init once only*/
					    /**  !! I DONT NEED FACTORIALS--THEY ARE JUST CONSTANTS IN THE ML */

  pTimeArray2tree(gRoot, p); /* put times from p[] onto tree */
  assignArrayRatesToLL2(gRoot, p); /* put all the rates from p[] onto tree */

  child=gRoot->firstdesc;
  SIBLOOP(child)
	{
	  rt=recursePenLike(child, gRoot);
	  if (rt == LARGE_VAL)
	    return rt;
	  else
	    Like+=rt; /* do it on one clade descended from root*/
	  /*i.e., make a tree in which root has only  one descendent node*/
	}

if (gNexDataPtr->RateBlockParms.NeighborPenalty==1)
	PY = penalizedRatesNeighbor(gRoot);
else
	PY = penalizedRates(gRoot);


gFit = -Like + PY;	/* ...this is currently bogus */
gLike=Like;
obj = Like + gSmoothing*PY;  /* remember we are minimizing everything */
if (DEBUG>1)
    printLikes(gRoot);
if (DEBUG>0)
    printf("AT RETURN in objPenLike (obj,penalty): %e\t%e\n",obj,PY);
 return obj;
}
/**********************/

static double recursePenLike(NODETYPE *node, NODETYPE * itsAncestor)

/* Recursively calculates the likelihood part of the penalized like */

{
  extern int gRatesAreGamma;
  extern double gAlpha;
  extern long gNumSites;
  NODETYPE *child;
  double obj,d,rt;
  
  if (!node) 
	return(0.0);
  d=itsAncestor->time-node->time;

  if (gRatesAreGamma)
    	obj=BranchLikeSumNegBinomial(node->estRate,gNumSites,gAlpha,d,node->length);
  else
  	obj=BranchLike(node->estRate,d,node->length); /* first arg is the branch's rate, stored in estRate location */
  node->like = obj;                            /* store this for later display */
  if (obj == LARGE_VAL)
	return obj;
  if (isTip(node)) 
    return obj;
  child=node->firstdesc;
  SIBLOOP(child)
    {
      obj+=recursePenLike(child,node);
    }
  return obj;
}
/**********************/

static double penalizedRates(NODETYPE *root)

/* Calculates the penalty for rate variation across branches.
   Penalty consists of squared deviations of rates between ancestor and descendant branchs
   AND squared deviations among the immediate descendants of the root node.
*/

{
  extern struct NexDataType *gNexDataPtr;
  extern NODETYPE * gRoot;    /* This global is declared when the whole 
					    algorithm is called */

  int tomy=0;
  double obj=0.0,thisTime,ancTime,thisLength,basal_rate,lr,rt,s=0.0,ss=0.0,r;
  NODETYPE *child;
  child=root->firstdesc;
  SIBLOOP(child)
    {
      rt=recursePenalizeRates(child, root);
      if (rt==LARGE_VAL)
		return rt;
      else
      	obj+=rt;
      ++tomy;
/*** ...NB! if one of children of root is a tip, and we PRUNE it during CV, possible trouble here 
	(currently not a problem, because we don't predict the length along that branch) */
      r=child->estRate;
      if (gNexDataPtr->RateBlockParms.PenaltyType==1)
	r=log(r);
      s+=r;
      ss+=r*r; 
    }
#if 0
  obj+=2*( ss-s*s/tomy) ;	/* this is basically the variance of the rates descended immediately from the root node */

#else
/** NB. the factor of two is needed to get the gradient correct (the gradient assumes we are
minimizing the simple squared deviation; the code above minimizes the variance, off by a factor of two here */

  obj += (ss-s*s/tomy)/tomy;	/* exactly the variance AS OF 5/26/01*/
#endif



/* printf("**%e\n",obj);*/
  return obj; 
}
/**********************/

static double recursePenalizeRates(NODETYPE *node, NODETYPE * itsAncestor)
{
  extern struct NexDataType *gNexDataPtr;
  int copyindex;
  NODETYPE *child;
  double obj=0.0,d,ranc,rdesc,o, oVar;
  if (!node) return(0.0);
  ranc=node->estRate;
  if (ranc < 0.0)
	return LARGE_VAL;  /* clamp down on time violations */
  else
      {
      child=node->firstdesc;
      SIBLOOP(child)
	{
	  rdesc=child->estRate;
	  if (rdesc < 0.0)
		return LARGE_VAL; /* signal from local_rate:clamp down on time violations */
	  else
            switch (gNexDataPtr->RateBlockParms.PenaltyType)
		{
		case 0: o= ranc-rdesc;break; // normal "additive"
		case 1: o= log(ranc)-log(rdesc);break; // logarithmic
    		}
	   obj+=o*o; // These are all x^2 terms
    /*     printf("--%f %f %f %f %f %f\n",
	    itsAncestor->time,node->time,child->time,ranc,rdesc,o);*/
	}
     }  

  if (isTip(node)) 
    return obj;
  child=node->firstdesc;
  SIBLOOP(child)
    {
      obj+=recursePenalizeRates(child,node);
    }
  return obj;	
}

/**********************/
static double penalizedRatesNeighbor(NODETYPE *n)

/* Calculates the penalty for rate variation across branches.
   Penalty is the variance of rates around a node. 
*/

{
  double obj=0.0,rt;
  NODETYPE *child;
  if (isTip(n))
  	return 0.0;
  rt=NeighborVariance(n);
  if (rt==LARGE_VAL)
	return rt;
  else
	obj=rt;
  child=n->firstdesc;
  SIBLOOP(child)
    {
	obj+=penalizedRatesNeighbor(child);
    }
  return obj; 
}
static double NeighborVariance(NODETYPE *n)
{
/* Calculate the variance in rates around node n. If n is the root, just do the descendants. If n is a tip, ignore.*/

extern struct NexDataType * gNexDataPtr;
double var,r,s=0.0,ss=0.0,rt;
int numNeighbor=0;
NODETYPE *child;
if (isTip(n))
	return 0.0;  // is this OK for a log penalty!?
child=n->firstdesc;
SIBLOOP(child) // get the neighbors who are descendants
    {
	++numNeighbor;
	r=child->estRate;
	if (r<0) 
		return LARGE_VAL;
	if (gNexDataPtr->RateBlockParms.PenaltyType == 1) //...for log penalties
		{
		if (r==0.0) return LARGE_VAL;
		r = log(r);
		}
	s+=r;
	ss+=r*r;
    }

if (!isRoot(n)) // Unless it's the root, also count the ancestral branch
    {
	++numNeighbor;
	r=n->estRate;
	if (r<0)
		return LARGE_VAL;
	if (gNexDataPtr->RateBlockParms.PenaltyType == 1) //...for log penalties
		{
		if (r==0.0) return LARGE_VAL;
		r = log(r);
		}
	s+=r;
	ss+=r*r;
    }
var = (ss-s*s/numNeighbor)/numNeighbor;	
return var;
}
static double NeighborSum(NODETYPE *n, int * numBranches)
{
/* Calculate the sum of log rates around node n. If n is the root, just do the descendants. If n is a tip, ignore.
   Also returns the number of incident branches on that node (useful for later calcs of variance, etc.)*/

double s=0.0;
int numNeighbor=0;
NODETYPE *child;
*numBranches=0;
if (isTip(n))
	return 0.0;
child=n->firstdesc;
SIBLOOP(child) // get the neighbors who are descendants
    {
	++(*numBranches);
//printf("1###%f\n",child->estRate);
	s+=log(child->estRate);
    }

if (!isRoot(n)) // Unless it's the root, also count the ancestral branch
    {
	++(*numBranches);
//printf("2###%f\n",n->estRate);
	s+=log(n->estRate);
    }
return s;
}

/**********************/

static double recursePenalizeRates2(NODETYPE *node, NODETYPE * itsAncestor)

/* Attempt to do curvature (second derivative) minimization with discrete estimate of curvature*/

{
  int copyindex;
  NODETYPE *child;
  double obj=0.0,d,ranc,rdesc,o, oVar,rancanc,oo;
  if (!node) return(0.0);
  ranc=node->estRate;
  if (!isRoot(node->anc))
	{
  rancanc=node->anc->estRate;
  if (ranc < 0.0  || rancanc < 0)
	return LARGE_VAL;  /* clamp down on time violations */
  else
      {
      child=node->firstdesc;
      SIBLOOP(child)
	{
    
	  rdesc=child->estRate;
	  if (rdesc < 0.0)
		return LARGE_VAL; /* signal from local_rate:clamp down on time violations */
	  else
		{
          	oo=2*ranc-rancanc-rdesc;
		o= oo*oo;
    		}
    
	   obj+=o;
     /*    printf("--%f %f %f %f %f %f\n",
	    itsAncestor->time,node->time,child->time,ranc,rdesc,o);*/
	}
     }  
	}
  if (isTip(node))
/** former code 
    return obj;
  ***/
	{
	oo=2*node->estRate-ranc-rancanc;
	obj=oo*oo;
	} /** new code **/
  child=node->firstdesc;
  SIBLOOP(child)
    {
      obj+=recursePenalizeRates2(child,node);
    }
  return obj;	
}


/********************************************************************/
/********************************************************************/

double objLangFitch(double p[])

/* This is the objective function for the Langley and Fitch clock method.
 *
 *	-If gisConstrained is set, the objective function has a penalty added to it.
 *	-If gEstRoot is set, the root node is additionally estimated.  However, this can
 *	only be done accurately if terminal tips have times > 0, or if there are maximum
 *	age constraints

 */

{
  extern NODETYPE * gRoot,*gRootDesc;    /* This global is declared when the whole 
					    algorithm is called */
  extern int gisConstrained;	       /* are we doing a constrained optimization? */
  extern int isFeasible;
  extern int gEstRoot;
  double obj=0.0,rootObj, d1, d2, thisTime,ancTime,thisLength,rt, K, k1, k2;
  int nnodes,i;
  NODETYPE *child, *child1, *child2;
  extern int gFirstDesc;
  static int f=1;

  if (DEBUG > 1)
    printf("--DEBUG-- Entry to ObjLangFitch\nRate=%g\n", p[gNVar]);
  if (gisConstrained)
    obj+=penalty(p); 
  
  if (f) {setupLogFactLookup(); f=0;}       /* temporary kludge to call this init once only*/
					    /**  !! I DONT NEED FACTORIALS--THEY ARE JUST CONSTANTS IN THE ML */

  pTimeArray2tree(gRoot, p); /* put times from p[] onto tree */

  child=gRoot->firstdesc;
  SIBLOOP(child)
	{
	  rt=recurseLangFitch(child, gRoot, p);
	  if (rt == LARGE_VAL)
	    return rt;
	  else
	    obj+=rt; /* do it on one clade descended from root*/
	  /*i.e., make a tree in which root has only  one descendent node*/
	}
if (DEBUG>1)
    printLikes(gRoot);
if (DEBUG>0)
    printf("RETURN FROM ObjLangFitch: %e\n",obj);
 return obj;
}
/**********************/

double recurseLangFitch(NODETYPE *node, NODETYPE * itsAncestor, double p[])
{
  extern int gRatesAreGamma;		/* are rates gamma distributed across sites? */
  extern double gAlpha;
  extern long gNumSites;
  int copyindex;
  NODETYPE *child;
  double obj,d,rt;
  
  if (!node) 
	return(0.0);
  d=itsAncestor->time-node->time;
  if (gRatesAreGamma)
    obj=BranchLikeSumNegBinomial(p[gNVar],gNumSites,gAlpha,d,node->length);
  else
    obj=BranchLike(p[gNVar],d,node->length);     /* first arg is the rate, stored in last element of pTime array */
  node->like = obj;                            /* store this for later display */
  if (obj == LARGE_VAL)
	return obj;
  if (isTip(node)) 
    return obj;
  child=node->firstdesc;
  SIBLOOP(child)
    {
      obj+=recurseLangFitch(child,node,p);
    }
  return obj;
}
/**********************/
double objLangFitchLocal(double p[])

/* This is the objective function for the Langley and Fitch LOCAL clock method.
 *
 *	-If gisConstrained is set, the objective function has a penalty added to it.
 *	-If gEstRoot is set, the root node is additionally estimated.  However, this can
 *	only be done accurately if terminal tips have times > 0, or if there are maximum
 *	age constraints

 
  The rates are stored in the p[] array starting BACKWARDS from the last array element

*/

{
  extern NODETYPE * gRoot,*gRootDesc;    /* This global is declared when the whole 
					    algorithm is called */
  extern int gisConstrained;	       /* are we doing a constrained optimization? */
  extern int isFeasible;
  extern int gEstRoot;
  double obj=0.0,rootObj, d1, d2, thisTime,ancTime,thisLength,rt, K, k1, k2;
  int nnodes,i;
  NODETYPE *child, *child1, *child2;
  extern int gFirstDesc;
  static int f=1;

  if (DEBUG > 1)
    printf("--DEBUG-- Entry to ObjLangFitch\nRate=%g\n", p[gNVar]);
  if (gisConstrained)
    obj+=penalty(p); 
  
  if (f) {setupLogFactLookup(); f=0;}       /* temporary kludge to call this init once only*/
					    /**  !! I DONT NEED FACTORIALS--THEY ARE JUST CONSTANTS IN THE ML */

  pTimeArray2tree(gRoot, p); /* put times from p[] onto tree */
  assignArrayRatesToLL_LFLOCAL(gRoot,p);

  child=gRoot->firstdesc;
  SIBLOOP(child)
	{
	  rt=recurseLangFitchLocal(child, gRoot, p);
	  if (rt == LARGE_VAL)
	    return rt;
	  else
	    obj+=rt; /* do it on one clade descended from root*/
	  /*i.e., make a tree in which root has only  one descendent node*/
	}
if (DEBUG>1)
    printLikes(gRoot);
if (DEBUG>0)
    printf("RETURN FROM ObjLangFitch: %e\n",obj);
 return obj;
}

static double recurseLangFitchLocal(NODETYPE *node, NODETYPE * itsAncestor, double p[])
{
  extern int gRatesAreGamma;		/* are rates gamma distributed across sites? */
  extern double gAlpha;
  extern long gNumSites;
  int copyindex,rateIndex;
  NODETYPE *child;
  double obj,d,rt;
  
  if (!node) 
	return(0.0);
  d=itsAncestor->time-node->time;
  rateIndex=gNVar-node->modelID; /* gets the proper rateindex  for this node */
  if (gRatesAreGamma)
    obj=BranchLikeSumNegBinomial(p[rateIndex],gNumSites,gAlpha,d,node->length);
  else
    obj=BranchLike(p[rateIndex],d,node->length);     /* first arg is the rate, stored in last element of pTime array */
  node->like = obj;                            /* store this for later display */
  if (obj == LARGE_VAL)
	return obj;
  if (isTip(node)) 
    return obj;
  child=node->firstdesc;
  SIBLOOP(child)
    {
      obj+=recurseLangFitchLocal(child,node,p);
    }
  return obj;
}
/***********************************************************************************/
/***********************************************************************************/

double objNP(double p[])

/* This is the objective function for the NPRS method.
'npexp' is the exponent in the smoothing function. Make it global below */

/* NB. DOES NOT YET IMPLEMENT 'estroot' OPTION */

/* At one time I added 1.0 to the objective function always to avoid the case where obj=0.0 for clocklike data.
   This required subtracting 1.0  at the end of the day. I've changed the latter back, but is this right? */
{
  extern struct NexDataType *gNexDataPtr;	
  extern NODETYPE * gRoot,*gRootDesc;    /* This global is declared when the whole 
					    algorithm is called */
  extern int	gisConstrained;	       /* are we doing a constrained optimization? */
  extern double gnpexp; /* global from ReadNexusFile */ 
  extern int	gClampRoot; /* global from ReadNexusFile */

  extern int gEstRoot;
  extern int isFeasible;
  static int firstTime=1,num_branches;
  double obj=0.0,thisTime,ancTime,thisLength,basal_rate,lr,rt,s=0.0,ss=0.0,r;
  int nnodes,i,tomy=0;
  NODETYPE *child;
  
  if (firstTime)
	{
	num_branches=numBranches(gRoot);
	firstTime=0;
	}

  if (gisConstrained)
    obj+=penalty(p); 
  

  pTimeArray2tree(gRoot,p);


/*** Now find objective function over rest of tree ***/

 child=gRoot->firstdesc;
  SIBLOOP(child)
    {
	++tomy;
	r=local_rate(child);
	if (gNexDataPtr->RateBlockParms.PenaltyType == 1) //...for log penalties
		r = log(r);
	s+=r;
	ss+=r*r;
      rt=recurseNP(child, gRoot, p);
      if (rt==LARGE_VAL)
		return rt;
      else
      	obj+=rt; /* do it on one clade descended from root*/
      /*i.e., make a tree in which root has only  one descendent node*/
    }

  obj += (ss-s*s/tomy)/tomy;	/* exactly the variance AS OF 5/26/01*/


  return obj; 
}
/**********************/

double recurseNP(NODETYPE *node, NODETYPE * itsAncestor, double p[])
{
  extern struct NexDataType *gNexDataPtr;	
  int copyindex;
  NODETYPE *child;
  double obj=0.0,d,ranc,rdesc,o, oVar;
  extern double gnpexp; /* global from ReadNexusFile */ 
  if (!node) return(0.0);
  ranc=local_rate(node);
  if (ranc < 0.0)
	return LARGE_VAL;  /* clamp down on time violations */

  if (gVarMinFlag)  /* if we are minimizing the variance of rates! */
  	obj+= /* SQR */ (ranc); /*** IS THIS A MISTAKE?  ***/
  else
      {
      child=node->firstdesc;
      SIBLOOP(child)
	{
    
	  rdesc=local_rate(child);

	  if (rdesc < 0.0)
		    return LARGE_VAL; /* signal from local_rate:clamp down on time violations */
	  if (gNexDataPtr->RateBlockParms.PenaltyType == 0)
	  	o= pow(fabs(ranc-rdesc),gnpexp);
	  else
	  	o= pow(fabs(log(ranc)-log(rdesc)),gnpexp);
	  obj+=o;
	}
     }  

  if (isTip(node)) 
    return obj;
  child=node->firstdesc;
  SIBLOOP(child)
    {
      obj+=recurseNP(child,node,p);
    }
  return obj;	
}

double local_rate(NODETYPE *node)

/* Estimates the local rate of evolution for the branch subtending this node. */

{
double rlocal,cumul=0.0,rlocal_desc,rl,rt,rL,rT;
int numrates=0;
NODETYPE *nodes_anc, *child,*this;
if (isRoot(node))
	{
	 fatal ("attempted to estimate local rate at root");
	}
rT = node->anc->time-node->time;
if (rT <= 0.0) 
			return -1.0; 	/* force a LARGE_VAL in calling routine */
rL = node->length;
return rL/rT;
}



/***********************************************************************************/
double mean_rate(NODETYPE *node)

/* Estimates the summed rate of evolution over all branches descended from node.
 Have to divide by number of branches after exiting
*/

{
double rlocal,cumul=0.0,rlocal_desc,rl,rt,rL,rT;
int numrates=0;
NODETYPE *nodes_anc, *child,*this;

if (!isRoot(node))
	{ 
	rL = node->length;
	rT = node->anc->time-node->time;
	cumul+= rL/rT;
	}

if (!isTip(node))
{
child=node->firstdesc;
SIBLOOP(child)
    {
    cumul+=mean_rate(child);
    }
}
return cumul;

}

/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
/***********************************************************************************/
double LFchiSq(NODETYPE *node,  double rate)

/* NB!  PAUP phylograms will have branch lengths on a scale from 0 to something less than
 * 1.0 (usually).  These are frequencies,  rather than numbers of changes.   The chi-sq test
 * needs numbers of changes!  You have to multiply the reported chi-sq value by the sequence length!
 */

{
  double cs=0.0;
  NODETYPE *child;
  child=node->firstdesc;
  SIBLOOP (child)
    {
      cs+=LFcs1(child,node,rate);
    }
  return cs;
}

double LFcs1(NODETYPE *node,  NODETYPE *itsAncestor, double rate)
{
  int copyindex;
  double chiSq=0.0,expected;
  NODETYPE *child;
  
  if (!node) return 0.0;
  expected = rate * (itsAncestor->time - node->time);
  if (fabs(expected) > 0.0001)  /* if duration is zero we might have a problem; for now this would
	mean that the expected change is zero.  Generally the observed change is zero too.  I think
	in a chi-squared test we would just not count this cell! */
	{
	chiSq= SQR(node->length - expected)/expected;
	node->chiSq = chiSq; /* if its zero, this doesn't get recorded...problem elsewhere? */
	}
  if (isTip(node)) 
    {
      return chiSq;
    }
  child=node->firstdesc;
  SIBLOOP(child)
    {
      chiSq+=LFcs1(child,node,rate);
    }
  return chiSq;	
}
/***********************************************************************************/
void printnodeLike(NODETYPE *node)
{
    double duration;
    NODETYPE *anc;
    if (!isRoot(node))
	{
	anc=node->anc;
	duration=node->anc->time-node->time;
	printf("node %3i (%s) age=%4.2f | anc %3i (%s) age=%4.2f | dur=%4.2f len=%4.2f like=%e\n",
	    node->id, node->taxon_name,node->time, 
	    anc->id, anc->taxon_name, anc->time, 
	    duration,node->length, node->like);
	}
    return;
}
/***********************************************************************************/
double BranchLike(double rate, double timeLength, double charLength)
{
  /* calculates NEGATIVE log likelihood of a branch whose characters are evolving according
     to a Poisson process (negative only because we are minimizing!)*/
/* NOTE I CAN IGNORE THE FACTORIAL STUFF AS LONG AS I JUST NEED ML ESTIMATES AND LR TESTS;
   HOWEVER, FOR THE TIME BEING, I'M LEAVING IT IN, AS IT SHOULD BE PRETTY FAST ANYWAY. ON
   8.9.00 I LOOKED INTO DELETING IT IN THE CALCULATION OF Z, BELOW, AND INDEED IT STILL WORKS. */

// Important: This routine converts charLength to an integral value! 
  
  double x,z, z1;
//*******  TEMPORARY
// if (rate == 0.0)
//  	return LARGE_VAL;
//*******
	
  if (timeLength<0.0 || rate<0.0) 
	return  LARGE_VAL;  /* must negate -infinity to minimize */	
  x= rate*timeLength;
  if (x > 0.0)
  	z =  -((int)charLength*log(x)-x   - logFact((long)charLength)  );
  if (x == 0.0 && charLength > 0.0)
	z =  LARGE_VAL;
  if (x == 0.0 && charLength == 0.0)
	z = 0.0;	
/* printf("BL:Poisson:(r,T,k,L)%e %e %e %e\n",rate,timeLength,charLength,-z);*/ 
  return z;
}

static double BranchLikeSumNegBinomial(double rate, long nSites,double alpha, double T, double k)

/* 
Log likelihood of a sum of nSites negative binomial distributions, each of which is governed by
a NB distribution with parameters alpha, and rate*T/alpha. Each of these results from a compounding
of a Poisson distribution with parameter rate*T and a gamma distribution of rates with parameter 
alpha. The gamma is normalized so that it has mean of rate*T (by using rate*T/alpha as the first
paramater and alpha as the second parameter of a gamma distribution.

NOTE! This function divides the input arg rate by nSites so that it works in terms of substitutions
per site internally, but the rest of the program sees it in units of total substitutions. This keeps it
consistent with Langley Fitch w/o gamma and PL w/o gamma.

*/

{
  
  double lz,zzz,p,q,x,na;

/***/

rate/=nSites;

/***/

  if (T <=0.0) return LARGE_VAL;	/*to accomodate log(0) or log (-x). */
                                        /*Must negate this because we are minimizing!*/
  if (T==0.0  && k==0.0 && (rate > 0.0 && alpha > 0.0) )
/*  if (ZERO(T) && ZERO(k) && (rate > 0.0 && alpha > 0.0) ) */
    return (0.0);
  if (rate <=0.0 || alpha<= 0.0) return LARGE_VAL;

  x=rate*T/alpha;
  p=1/(1+x);
  q=1-p;
  na=nSites*alpha;
  lz=gammln(k+na)+k*log(q)+na*log(p)   -gammln(na)-gammln(k+1);

/* printf("**NegBin:%e %e %e %e %e %e %e\n",rate,p,gammln(k+na),k*log(q),-gammln(na),-gammln(k+1),na*log(p));*/

/*zzz=BranchLike(b*c,T,k);
printf("BLgamma:(b,c,T,k,L) %e %e %e %e %e %e\n",b,c,T,k,lz,-zzz);*/

  return -lz; /*(its a minimization function, stupid)*/
}

double BranchLikeGamma(double b, double c, double T, double k)
{
  /* calculates log likelihood of a branch whose characters are evolving according
     to a Poisson compounded with a gamma distribution (params b,c), which is negative binomial.
	k = number of substitutions
	T = branch duration */
  
  double lz,zzz;
  if (T <=0.0) return LARGE_VAL;	/*to accomodate log(0) or log (-x). */
                                        /*Must negate this because we are minimizing!*/
  if (T==0.0 && k==0.0 && (b > 0.0 && c > 0.0) )
    return (0.0);
  if (b <=0.0 || c<= 0.0) return LARGE_VAL;

/*k=4;b=0.001;c=1000.0;T=1.0; ...test values...*/


  lz=gammln(k+c)+k*log(T*b)-gammln(c)-gammln(k+1)-(c+k)*log(1+T*b);

/*printf("**BLgamma:%e %e %e %e %e\n",gammln(k+c),k*log(T*b),-gammln(c),-gammln(k+1),-(c+k)*log(1+T*b));*/

/*zzz=BranchLike(b*c,T,k);
printf("BLgamma:(b,c,T,k,L) %e %e %e %e %e %e\n",b,c,T,k,lz,-zzz);*/
  return -lz; /*(its a minimization function, stupid)*/
}

double BL(double rate, double timeLength, double charLength)
{
/* calculates likelihood of a branch whose characters are evolving according
to a Poisson process */

  double x,z;
  if (timeLength<=0.0) return 0.0;
  if (timeLength==0.0 && charLength==0.0  && (rate > 0.0) )
    return (1.0);
  if (rate <= 0.0 ) return 0.0;  
  x= rate*timeLength;
  /* z=exp(-x)*pow(x,charLength)/FactLookup[(int)charLength]; */

  z= -x + charLength*log(x) - logFact((long)charLength);
  z= exp(z);

  return z;
}

double Factorial(double x)  /* never gets called now ! */
{
  extern double FactLookup[];
  if (x <= MAX_FACTORIAL)
    return FactLookup[(int)x];
  else
    doGenericAlert ("Factorial size exceeded");
  return 0.0;
}

void setupFactLookup(void)
{
#if 0   /* Some hardware throws an exception for overflows below
  extern double FactLookup[];
  int i;
  FactLookup[0]=1.0;
  for(i=1;i<=50;i++)
    FactLookup[i]=i*FactLookup[i-1];
  for (i=51;i<=MAX_FACTORIAL;i++)
    FactLookup[i]=pow(i/2.7182818,(double)i)*sqrt(2*3.14159*i);
/* use stirling formula for 51 < n < Max */
#endif
  return;
}

void setupLogFactLookup(void)
{
  extern double FactLookup[];
  extern double logFactLookup[];
  double z;
  int i;
  FactLookup[0]=1.0;
  for(i=1;i<=50;i++)
    FactLookup[i]=i*FactLookup[i-1];
/********** 6.10.02 I had forgot to include the following statement! Important for 0-length branches ******/ 
  logFactLookup[0]=0.0;
  for(i=1;i<=50;i++)
    logFactLookup[i]=log(FactLookup[i]); /* precompute the log of these */
  for (i=51;i<=MAX_FACTORIAL;i++)
    {
    z=i*(log((double)i)-1)+log(sqrt(2*3.14159*i));
    logFactLookup[i]=z;
    }
/* use stirling formula for 51 < n < Max */

/*for (i=1;i<=MAX_FACTORIAL;i++)
	printf("%e\n",logFactLookup[i]);*/

  return;
}

double logFact(long k) /* calculate factorials */
{
  if (k<=MAX_FACTORIAL)
    return logFactLookup[k]; /* use precomputed values if arg < MAX, otherwise use Stirling's formula */
  else
    return k*(log((double)k)-1)+log(sqrt(2*3.14159*k));
}

void plotOpt(double p[],int grid,double p1low,double p1high,
	double p2low,double p2high, char *p1label,  char *p2label)
{
int i,j;
double p1interval,p2interval,p1,p2,obj;
p1interval=(p1high-p1low)/grid;
p2interval=(p2high-p2low)/grid;

printf("\t\t\t%s\n", p2label);

for (i=1;i<=grid;i++)
	{
	p1=p1low+p1interval*(i-1);
	if (i==grid/2)
		printf("%6s\t\t", p1label);
	else
		printf("\t\t");
	for (j=1;j<=grid;j++)
		{
		p2=p2low+p2interval*(j-1);
		p[1]=p1;
		p[2]=p2;
		obj=BD_Like(p);
		printf("%6f  ",obj);
		}
	printf("\n");
	}
}

void check_if_exceed_factorial(NODETYPE *node)
{
    
  NODETYPE *child;
  
  if (!isRoot(node))
/*...no longer needed 
    if (node->length > MAX_FACTORIAL) 
	{
	printf("%e ---->",node->length);
	fatal("Factorial size exceeded\n");
	}
*/
    if (node->length > 0.0 && node->length < 1.0) 
	{
	printf("%e ---->",node->length);
	printf("WARNING! A branch length was a real number between 0.0 and 1.0;\nlengths must be integers representing absolute numbers of substitutions\n"); 
	printf("Attempting to fix by rounding\n");
	if (node->length <0.5)node->length=0.0;
	else node->length=1.0;
	}
  child=node->firstdesc;
  SIBLOOP(child)
    {
    check_if_exceed_factorial(child);
    }
  return;	
    
}


