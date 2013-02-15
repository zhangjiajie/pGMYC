#include "TreeUtils.h"
#include "penalty.h"
#include "math.h"
#include "nexus.h"
/* Module for penalty/barrier function */

/* When do we need to worry about constraints on the tree?  It turns out, at least
for Langley and Fitch algorithm, that the tree-constraints (that is, a descendant
must be younger than an ancestor, and all tips are at time 0), do NOT have to be
included explicitly via a penalty function.  I THINK this is because in the
absence of fossils, the optimum
of the function is always on the interior of the feasible space--it is never on a
boundary (i.e., we never reconstruct any branch duration to be zero under L-F (although
I haven't expressly looked at cases where the number of substitutions on a branch is
zero!).

However, when fossils are included, a time constraint may cram the reconstructed time
of a node right up against a constraint, and the standard POWELL seems to go goofy. Note
that the calculation of the likelihood currently blows up (returns HUGE_VAL) when any
duration becomes <= 0.  This is OK when no constraints are in force, but seems to be
inadequate to get POWELL to work properly when a constraint is enforced on an internal
node.  I think that is because the likelihood is going to increase smoothely to infinity
as the branch shrinks to zero under non-constraint circumstances, but when there is a constraint,
the likelihood blows up somewhere short of a branch length of zero, as soon as the fossil
constraint is reached, so there is a jump.  Point is that we can't just, say, check to 
see if we're violating the constraint, and then add some big number to the likelihood
as seen as we are violating it.

*/
/* In this kind of constrained optimization we maximize F(x) subject to some constraints,
g(x).  We do this by maximizing a different function R(x) = F(x) + k G(x), where G(x) is
a penalty or barrier function based on g(x).   For
example, here we will use a reciprocal function G(x) = 1/g(x), where we write all
constraints as g(x) >= 0.  See the function 'addconstr' below. We start with some reasonable
value of k and repeat the optimization each time reducing k appropriately.  See
the constants in 'constrOpt'
*/



#define LARGE_VAL 1e+30

int isFeasible,gPenaltyIx;
NODETYPE * gRootDescPenalty;	/* initialized in 'ObjFunc.c'*/

/*********************************************************************/

double penalty(double p[])


{
extern int gEstRoot;
extern double rk;	/* this is the constant k as described above */
extern int isFeasible;
isFeasible=1;

/* 7.24.00: 

	If we are estimating the root node, we assume that there are no
	constraints on its ages.  However, there is one more node, the root,
	in the p[] array.  To ignore it, we start indexing one past it at [2].
*/

gPenaltyIx=1;
return rk*traversePenalty(gRootDescPenalty,p);
}
/*********************************************************************/

double penalty_rate(double p[])

/* This clamps the rate ratio to less than or equal some value MAX_RATIO */
#define MAX_RATIO 10.0
{
extern double rk;	/* this is the constant k as described above */
extern int isFeasible;
extern int gNVar;

double penalty=0.0;
isFeasible=1;

/*penalty = addconstr(MAX_RATIO-p[gNVar+1]) + addconstr(p[gNVar+1]);*/
/*printf("##%e\t%e\t%e\n",p[gNVar+1],penalty,rk*penalty);*/
return rk*penalty;
}
/*********************************************************************/
double traversePenalty(NODETYPE *node, double p[])
{
	double penalty=0.0,temp;
	NODETYPE *child;
	if (!node) return(-1);
	if (isFree(node))
	  {
	  if (node->nodeIsConstrainedMin)
			penalty += addconstr( p[gPenaltyIx] - node->minAge);
	  if (node->nodeIsConstrainedMax)
			penalty += addconstr( node->maxAge - p[gPenaltyIx]);
	  ++gPenaltyIx;
	  }
	if (!isTip(node) ) 
	  {
	  child=node->firstdesc;
	  SIBLOOP(child)
			penalty+=traversePenalty(child,p);
	  }
	return penalty;
}

/*********************************************************************/
double addconstr(double x)

/* Adds a reciprocal constraint to the penalty function; 
parameter 'x' means that x > 0 is a constraint */

{
extern int isFeasible;
if (x>0.0) return 1/x;
else 
	{
	isFeasible=0;
	return LARGE_VAL;	/* shouldn't ever need this value */
	}

}


/**********************************************************************/

void check_feasible(NODETYPE *node)

/* Checks if the times currently set on tree are feasible. This is done in the context
	of whether the search is calling for time constraints or not. If it is, then 
	every point must satisfy relevant min and max age constrains.  If not, then it
	must merely obey the tree constraints (age can't be older than its ancestor).
	The next routine is identical, except it moves through the WHOLE tree and prints
	error messages.  Current routine bails at first violation.

 NB!!! Have to set isFeasible to 1 prior to call, then check it */

{
	extern int isFeasible,gisConstrained;
	NODETYPE *child;
	if (!isFeasible)
		return;
	if (!isRoot(node))
				{
				if(node->time > node->anc->time)
					{
					isFeasible=0;
					return;
					}
				}
	if ((node->nodeIsConstrainedMin) && (node->time < node->minAge))
				{
				isFeasible=0;
				return;
				}
	if ((node->nodeIsConstrainedMax)&&(node->time > node->maxAge))
				{
				isFeasible=0;
				return;
				}
	if (!isTip(node))
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			check_feasible(child);
		}
	return;
}
/**********************************************************************/
/**********************************************************************/

void debug_check_feasible(NODETYPE *node)
/* prints out useful stuff if desired when a point is not feasible */

{
	extern int isFeasible,gisConstrained;
	NODETYPE *child;
	if (!isRoot(node))
				{
				if(node->time >= node->anc->time)
					{
					printf("FEASIBLE VIOLATION: node %s:%i (%f) is older than ancestor %s:%i (%f)\n", 
					    node->taxon_name,node->id,node->time, node->anc->taxon_name,node->anc->id,node->anc->time);
					isFeasible=0;
					}
				}
	if ((node->nodeIsConstrainedMin) && (node->time <= node->minAge))
				{
				printf("FEASIBLE VIOLATION: node %s:%i (%f) is younger than its min age (%f)\n", 
					    node->taxon_name,node->id,node->time, node->minAge);
				isFeasible=0;
				}
	if ((node->nodeIsConstrainedMax)&&(node->time >= node->maxAge))
				{
				printf("FEASIBLE VIOLATION: node %s:%i (%f) older than its max age (%f)\n", 
					    node->taxon_name,node->id,node->time, node->maxAge);
				isFeasible=0;
				}
	if (!isTip(node))
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			debug_check_feasible(child);
		}
	return;
}
/**********************************************************************/

int check_initial_point(double (*objective)(double p[]),  double p[])

/**** only works if the tree structure has the p[] times on it! ***/

{
    extern struct NexDataType *gNexDataPtr;
    extern int gNVar;
    extern NODETYPE * gRoot;	
    extern int isFeasible;
    int i;
    double f_init;
    f_init=(objective)(p);
/*
    if (gNexDataPtr->RateBlockParms.verbose)
	printf("Objective function at initial feasible point=%f\n", f_init);
*/ 
    isFeasible=1;
    check_feasible(gRoot);
    if (!isFeasible)
			{
			doGenericAlert("A point was NOT feasible");
			printf("The point:\n");
			for (i=1;i<=gNVar-1;i++)
			    printf("p[%2i] %6f\n",i,p [i]);
debug_check_feasible(gRoot);
			return 0;
			}
    return 1;
}
