/*

Module to implement a simple vaguely covarion-like switch model with binary switch,S, and binary trait, T.
Modified massively from ancestral.c

7.24.2012. Worked on getting a better initial rate estimate
	- and the zero-length branch problem (see below, grep for 'cherry')

Contains several routines from NRC and some minor modifications to same...

*/
#include "ObjFunc.h"
#include "Maximize.h"
#include "structures.h"
#include "nexus.h"
#include "covarion.h"
#include "TreeUtils.h"
#include "DistrFuncs.h"

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SQR(x)        ((x)*(x))
#define MAX(a,b) ((a) >(b) ? (a):(b))
#define LARGE_VAL 1e100


#define MIN_BRANCH_LENGTH 1e-6  // we do not allow calculations on 0-length branches
				// imagine a cherry with a 0 in one leaf and 1 in the other
				// and 0-length branches to the mrca. That's a 0 likelihood

double rootPrior[3];

double gR, gS, gQ01;  //switch rate = s; trait rate = r; fixed s rate = qQ01 for that special case

double ** covProb;

void writeNexusChar(NODETYPE *node);
void writeNexusCharRecurse(NODETYPE *node);

void simulatePrecursorCharRecurse(NODETYPE *node, int nodeState, double s, double r);
void simulatePrecursorChar(TREE t, double s, double r );
void surface_peek(double (*obj)(double p[]), int dim, double x1low, double x1high, double x2low,double x2high);
void peak_peek_2(double p[], double (*obj)(double p[]), int dim, double scaleFactor, double ftol);
void simulateBinaryCharRecurse(NODETYPE *node, int nodeState, double q01, double q10);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
static double ** transitionProb2stateSymmetric(double mu, double t);
static double ** transitionProb2state(double pi1, double pi2, double t);
static double ** transitionProb(double s, double r, double t);
static double ** transitionProb4(double s,  double t);
static double ** mat_transpose(double **A);
static double ** mat_mult(double **A,double **B);
static void mat_print(double **m, int maxi, int maxj);
double **dmatrix(long nrl, long nrh, long ncl, long nch);	// defined in continuousML.c
static void uppassCovarion(NODETYPE *node);
static void setupCLmaxopt(NODETYPE *node);
static void setupCL(NODETYPE *node); // conditional likelihoods



double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);

int nParm, gNT, p2tindex; 

int NSTATES;
int model;
int verbose;


void covarionOptimize(TREE t,int *numIter, double ftol,double linMinDelta,int *success, int nmodel, int doMarginals, int estimateFlag, int doRecon )
{

extern struct NexDataType *gNexDataPtr;
extern NODETYPE * gRoot;
StrListPtr DM, TL;
PtrList nodeList;
int NT,i,j, ii,ixTL, num_retry=5,count1,count0;
double *p, *pretry,*D, obj, logobj, logobjretry,**PT, AIC, rate_init_estimate;
char * tip_name,*tip_value_str, *dummy,*found_tip;
NODE a;
float meanTip=0.0;
int tip_value,verbose;
double minobj,minS,minR;
double **m1, **m2;
struct RBP * rbp;

int impute=0;

//**** Turn on to print most output ****//

verbose=gNexDataPtr->RateBlockParms.verbose;

//**** ******************************* ****//

double (*obj_func_array[10])(double[]);
obj_func_array[0] = objBinaryTraitSymmetric;
obj_func_array[1] = objBinaryTrait;
obj_func_array[2] = objCovarion;
obj_func_array[3] = objCovarion; 	// this is the precursor_1 model which uses same objective with a tweak to set two parms equal
obj_func_array[4] = objCovarion4;	// covarion 4 state 1 parm model


// ***!!!! FIGURE OUT WHY ITS OK TO USE NSTATES IN ARRAY DECLARATIONS THIS SEEMS BOGUS TO ME  ***

rbp=&(gNexDataPtr->RateBlockParms);


model=nmodel; //make it global

switch (model)  
			{
			case 0: // binary trait, symmetric, one parm
				NSTATES=2;
				nParm=1;
				rootPrior[0]=1/2.; rootPrior[1]=1/2.; 
				break;
			case 1: // binary trait, asymmetric, two parms
				NSTATES=2;
				nParm=2;
				// root Prior is set to estimated stationary frequency in objBinaryTrait
				break;
			case 2: // three-state model, two parms
				NSTATES=3;
				nParm=2;
				rootPrior[0]=1/3.; rootPrior[1]=1/3.; rootPrior[2]=1/3.;
				break;
			case 3: // three-state model, one parm
				NSTATES=3;
				nParm=1;
				rootPrior[0]=1/3.; rootPrior[1]=1/3.; rootPrior[2]=1/3.;
				break;
			case 4: // four state covarionish, one parm
				NSTATES=4;
				nParm=1;
				rootPrior[0]=1/4.; rootPrior[1]=1/4.; rootPrior[2]=1/4.; rootPrior[3]=1/4.;
				break;
			}
					


DM=gNexDataPtr->DMList;
TL=gNexDataPtr->TaxaList;
gRoot = t->root;
gRoot -> length = 0.0;
if (verbose) 
	{
	printf ("Warning! Setting the root's subtending branch length to ZERO whether it is previously set or not\n");
	if (gNexDataPtr->RateBlockParms.cov_brlens==1)
		printf("All branch lengths are set to = 1.0\n");
	if (gNexDataPtr->RateBlockParms.cov_brlens==0)
		printf("Using user supplied branch lengths\n");
	printf("Using model %i with %i states and %i parameters\n",model, NSTATES,nParm);
	}
nodeList=pNewList();
TreeToTaxaPtrList(t->root,nodeList);




NT=t->numTaxa;
for (i=1;i<=NT;i++)
	{
	a=(NODE)(pListgetkthNode(nodeList, i)->item);
	tip_name=a->taxon_name;
	ixTL=findMatchStr(TL,tip_name); // lookup the taxon name from the matrix ordering and get the relevant data matrix ordering corresponding name
	if (ixTL == 0)
		{
		printf("Could not find tree taxon name %s in matrix\n",tip_name);
		exit(1);
		}
	tip_value_str=getkthStr(DM,(long)(ixTL));
	found_tip=getkthStr(TL,(long)(ixTL));
	tip_value=strtod(tip_value_str,NULL);
// printf("%i\t%s\t%s\t%s\t%f\n",i,tip_name,found_tip,tip_value_str,tip_value);

// set up conditional likelihoods of leaves based on data

	a->state=*tip_value_str; // store this for later display conveniently in DrawTree
	if (model==4)
		{
		switch (*tip_value_str)  // this is a single char
			{
			case '0': 
				(a->CL)[0] = 1.0;
				(a->CL)[1] = 1.0;
				(a->CL)[2] = 0.0;
				(a->CL)[3] = 0.0;
				++count0;	
				break;
			case '1': 
				(a->CL)[0] = 0.0;
				(a->CL)[1] = 0.0;  
				(a->CL)[2] = 1.0;
				(a->CL)[3] = 1.0;
				++count1;	
				break;
			case '?': 
				(a->CL)[0] = 1.0;
				(a->CL)[1] = 1.0;  
				(a->CL)[2] = 1.0;
				(a->CL)[3] = 1.0;
				break;
			}
		}
	if (model==2 || model==3)
		{
		switch (*tip_value_str)  // this is a single char
			{
			case '0': 
				(a->CL)[0] = 1.0;
				(a->CL)[1] = 1.0;
				(a->CL)[2] = 0.0;
				++count0;	
	
				break;
			case '1': 
				(a->CL)[0] = 0.0;
				(a->CL)[1] = 0.0;  
				(a->CL)[2] = 1.0;
				++count1;	
				break;
			case '?': 
			
				if (impute)
					{
					if (a->opt == 0)  // leftover from prev run
						{
						(a->CL)[0] = 1.0;
						(a->CL)[1] = 1.0;  
						(a->CL)[2] = 0.0;
						}
					if (a->opt == 1)
						{
						(a->CL)[0] = 0.0;
						(a->CL)[1] = 0.0;  
						(a->CL)[2] = 1.0;
						}
					}
				else
					{
					(a->CL)[0] = 1.0;
					(a->CL)[1] = 1.0;  
					(a->CL)[2] = 1.0;
					}
				break;
			}
		}
	if (model==0 || model ==1)
		{
		switch (*tip_value_str)  
			{
			case '0': 
				(a->CL)[0] = 1.0;
				(a->CL)[1] = 0.0;
				++count0;	
	
				break;
			case '1': 
				(a->CL)[0] = 0.0;
				(a->CL)[1] = 1.0;  
				++count1;	
				break;
			case '?':
				if (impute && model==1)
					{
					(a->CL)[0] = 0.0;
					(a->CL)[1] = 0.0;  
					(a->CL)[a->opt]=1.0;  // this will be a leftover from any previous run
					
					printf("Imputing taxon %s, setting a->CL[%i] to 1.0\n",a->taxon_name,a->opt);
					
					}
				else
					{
					(a->CL)[0] = 1.0;
					(a->CL)[1] = 1.0;
					}
				break;
			}

		}
	} // end NT loop

p = dvector(1,nParm);
pretry = dvector(1,nParm);
D = dvector(1,nParm); //gradient


if (estimateFlag==1)
	{
	// Set up initial search start

if (gNexDataPtr->RateBlockParms.cov_brlens==1)
	rate_init_estimate = MIN(count0,count1)/(float)numBranches(t->root);	
else	
	rate_init_estimate = MIN(count0,count1)/treeLength(t->root); // roughly an upper bound...
	gS=gR=rate_init_estimate;
	if (verbose)
		printf("Using initial estimate of rates = %g\n",gS);
	
	p[1] = gS; 
	if (nParm==2) p[2] = gR; 
	
	
	// do the search
	
	logobj= -MinND(t,0, POWELL,obj_func_array[model],NULL,p, nParm,numIter, ftol,linMinDelta,success);
	
	// check the solution if requested
	
	if (model==2 || model==1) // retry to look for boundary estimate for s_rate >> 1
	  {
	  num_retry=2;
	  for (ii=1;ii<=num_retry;ii++)
		{
		*numIter=100;
		//pretry[1]=10;pretry[2]=0.1;
		pretry[1]=0.01*2*ii;pretry[2]=0.01*2*ii;
		logobjretry= -MinND(t,0, POWELL,obj_func_array[model],NULL,pretry, nParm,numIter, ftol,linMinDelta,success);
		if ( (logobjretry >  logobj) && fabs ((logobjretry - logobj)/logobjretry) > ftol)
			{
			if (verbose) printf("Found a better solution on retry %i of model: initial estimate p[1]=%10.8e p[2]=%10.8e obj=%10.8e\n",ii,p[1],p[2],logobj);
			p[1]=pretry[1];
			p[2]=pretry[2];
			logobj=logobjretry;
			}
		else
			if (verbose) printf("Retry solution: p[1]=%10.8e p[2]=%10.8e obj=%10.8e\n",pretry[1],pretry[2],logobjretry);
		}
	  } // end if model 1,2
	if (rbp->checkGradient && verbose==1) // silly to check this if we are not going to print it...
			{
			printf("NumIter done = %i\n",*numIter);
			Dapprox(p,D,nParm, obj_func_array[model], 1e-11); // approximate the gradient
			printf("Approximate gradient: %g %g\n",D[1],D[2]);
			if (!checkGradientSimple(p,D, -logobj, ftol, nParm) ) // do a formal check by the books
				peak_peek_2(p, obj_func_array[model], 6, 0.01,ftol); // and if we fail, just look at the local surface
			}
	

	
	// Report model results	
	}
if (estimateFlag==0) // just calculate the likelihood
	{
	p[1] = rbp->s_rate;
	if (nParm==2) p[2] = rbp->r_rate;
	logobj=-(obj_func_array[model])(p);
	if (verbose) printf("Using user supplied rate parameters!\n");
	}



AIC = 2*nParm - 2*logobj;
if (verbose)
	{
	switch (model)
		{
		case 0: 
			printf("Estimated parameters: %f ; Maximum likelihood=%g\n",p[1],logobj); break;
		case 1:
			printf("Estimated parameters: %f %f; Maximum likelihood=%g\n",p[1],p[2],logobj);
			printf("Stationary frequencies: %f %f\n",p[2]/(p[1]+p[2]),p[1]/(p[1]+p[2]) );
			break;
		case 2:
			printf("Estimated parameters: %f %f; Maximum likelihood=%g\n",p[1],p[2],logobj); break; 
		case 3: 
			printf("Estimated parameters: %f ; Maximum likelihood=%g\n",p[1],logobj); break;
		case 4: 
			printf("Estimated parameters: %f ; Maximum likelihood=%g\n",p[1],logobj); break;
		}
	printf("Model log Like = %f;  AIC value for %i parameters = %f\n",logobj,nParm,AIC);
	}




// Possibly do the reconstruction using these rate(s)
	
gS = p[1]; if (nParm==2) gR = p[2];

if (doRecon)
	{
	 setupCLmaxopt(t->root);
	 uppassCovarion(t->root);
 	}




if (doMarginals)
	{
	if (verbose) 
		printf("Calculating marginals...\n");
	findMarginals(t->root);
	t->root=gRoot; //kludge : see marginals code below
	}

if (verbose)
	{
	if (model==0)
		printf("%i\t%i\t%6.3g\t   --   \t%6.3f\t%6.3f\n",NSTATES,nParm,p[1],logobj,AIC);
	else
		printf("%i\t%i\t%6.3g\t%6.3g\t%6.3f\t%6.3f\n",NSTATES,nParm,p[1],p[2],logobj,AIC);
	}
if (model==1 || model==2)
	{
//	if (verbose) printf("\n\nLikelihood surface\n");
//	surface_peek(obj_func_array[model], 50, 0.001, 0.051, 0.001, 0.021);
	}


free_dvector(p,1,nParm);
free_dvector(pretry,1,nParm);


//printf ("begin r8s; covarion estimate=no simulate=yes model=binary_2 brlens=user q01=%f q10=%f; end;\n",p[1],p[2]);

//writeNexusChar(gRoot);

}

// **********************************************************************************

//// Following does a one-char simulation using the 2-state assym model estimates and writes a nexus file
//// Root prior is assumed to be the stationary distribution based on the two rates
//// unless FIXROOT is set to 1, in which case we always fix the root to state 0.

#define FIXROOT 1

void simulateBinaryChar(TREE t, double q01, double q10 )
{
double rootPriorState1;

NSTATES=2;
if (FIXROOT)
	rootPriorState1=0;
else
	rootPriorState1 = q10/(q01 + q10);
simulateBinaryCharRecurse(t->root,rootPriorState1,q01,q10);
//printCovarion(t->root,0);
writeNexusChar(t->root);

//printf ("begin r8s;\ncovarion estimate=yes model=binary_2 brlens=user;\ncovarion estimate=yes model=switch brlens=user;\ncovarion estimate=yes model=switch_1 brlens=user;\ncovarion estimate=yes model=binary_1 brlens=user;\nend;\n");

}
// **********************************************************************************
void simulatePrecursorChar(TREE t, double s, double r )
{
double rootPriorState;
NSTATES=3;

rootPriorState = randTrinary(1/3.,1/3.); // uniform on 1/3,1/3,1/3
simulatePrecursorCharRecurse(t->root,rootPriorState,s,r);
writeNexusChar(t->root);
printf ("begin r8s;\ncovarion estimate=yes model=binary_2 brlens=user;\ncovarion estimate=yes model=switch brlens=user;\ncovarion estimate=yes model=switch_1 brlens=user;\ncovarion estimate=yes model=binary_1 brlens=user;\nend;\n");

}

void simulatePrecursorCharRecurse(NODETYPE *node, int nodeState, double s, double r)

// Simulate a character under the 2-state assym model, starting with some initial state, and storing
// trait in node->opt field.

{
    NODETYPE *child;
    double **PT;
	int child_state;
	node->opt = nodeState; // store the character state in this field
	child=node->firstdesc;
	SIBLOOP(child)
			{
			PT = transitionProb(s,r,child->length);
			child_state = randTrinary(PT[nodeState][0],PT[nodeState][1]); // returns a 0,1 or 2 according to rate matrix
//printf("**%f %f %i %i %f %f\n",s,r,nodeState,child_state,PT[nodeState][0],PT[nodeState][1]);
			simulatePrecursorCharRecurse(child,child_state,s,r);
  			free_dmatrix(PT,0,NSTATES-1,0,NSTATES-1);
			}
	return;
    
}




/***********************************************************************************/

// Find marginal ancestral state probabilities at every internal node.
// To do this, we have to reroot at each internal node, redo the calculations and report

void findMarginals(NODETYPE *root)
{
PtrList nodeList;
long NT;
int i,j,saveFlag;
NODETYPE * a, *rroot,*rootfirstdesc;
extern NODETYPE * gRoot;
double saveCL[NSTATES];

rootfirstdesc = root->firstdesc;
nodeList=pNewList();
TreeToNodePtrList(root,nodeList);
NT=numNodes(root);
for (i=1;i<=NT;i++)
	{
	saveFlag=0;
	a=(NODE)(pListgetkthNode(nodeList, i)->item);
	
//	if (isTip(a)) continue;  // skip rerooting at leaves. Careful here--if we were to reroot on a leaf (which is ok in principle),
							// then the setupCL function rewrites the CL values for this leaf, because it
							// is now being treated as the root. We'd have to reinitialize all the leaf
							// values every time to do something like that...
							
							// OK, really not sure whether to infer states at leaves. So I'll ignore for now
#if 1
	if (isTip(a))
		{
		saveFlag=1;
		for (j=0;j<NSTATES;j++)
			saveCL[j] = (a->CL)[j];
		}
	//printf("Rooting on node %li\n",a->id);
#endif
	ReRoot2(a); // use this kind of rerooting (type 2) which roots at a node.

//printf("Rerooting at node %li\n",a->id);
//DrawTree(a,1,100);

	gRoot = a;

	/*** Now do conditional likelihoods on tree ***/

  	setupCL(a);
  
	// Find the best weighted by priors; this will correspond to the max like across the tree (maybe slow
	// compared to calculating conditional likelihoods???

  	for (j=0;j<NSTATES;j++)
		{
		(a->CLmarg)[j] = rootPrior[j]*(a->CL)[j];
		}

#if 1
	if (saveFlag==1)  // all because a is no longer a tip but we want to restore it if it WAS a tip above (before it got rerooted)
		{
		for (j=0;j<NSTATES;j++)
			(a->CL)[j] = saveCL[j] ;
		saveFlag=0;
		}
#endif

	}
	

ReRoot2(root); // reroot on original root node and pray that everything is still the same
gRoot=root;

return;
}


/***********************************************************************************/

double objCovarion(double p[])


{
  extern struct NexDataType *gNexDataPtr;	
  extern NODETYPE * gRoot;    /* This global is declared when the whole algorithm is called */

  static int firstTime=1,num_branches;
  int i;
  double obj=-1e100,max=0.0,L;
  NODETYPE *child;
  
 if (p[1] <0 ) {return LARGE_VAL;};
 gS = p[1];
 
 if (nParm == 2)
 		{
 		if (p[2] <0 ) {return LARGE_VAL;};
 		gR = p[2];
 		}
 		
  //gS = gR = p[1]; 
  
/*** Now do conditional likelihoods on tree ***/

  setupCL(gRoot);
  
// Find the best weighted by priors; this will correspond to the max like across the tree (maybe slow
// compared to calculating conditional likelihoods???

  obj=0.0;
  for (i=0;i<NSTATES;i++)
		{
		obj += rootPrior[i]*(gRoot->CL)[i];
		}
//printf ("--------------------------------->%f %f %f\n",p[1],p[2],-obj);
//printCovarion(gRoot,0);
  return -log(obj); // it's a minimization
}


double objCovarionFixed(double p[]) // deprecated
{
  extern struct NexDataType *gNexDataPtr;	
  extern NODETYPE * gRoot;    /* This global is declared when the whole algorithm is called */

  static int firstTime=1,num_branches;
  int i;
  double obj=-1e100,max=0.0,L;
  NODETYPE *child;
  
 if (p[1] <0 ) {return LARGE_VAL;};
 
  gS = gQ01 ; gR = p[1];

//printf ("%f %f\n",gS,gR);
  
/*** Now do conditional likelihoods on tree ***/

  setupCL(gRoot);
  
// Find the best weighted by priors; this will correspond to the max like across the tree (maybe slow
// compared to calculating conditional likelihoods???

  obj=0.0;
  for (i=0;i<NSTATES;i++)
		{
		obj += rootPrior[i]*(gRoot->CL)[i];
		}
//printf ("--------------------------------->%f %f %f\n",p[1],p[2],-obj);
//printCovarion(gRoot);
  return -log(obj); // it's a minimization
}
/**********************/

static void setupCLmaxopt(NODETYPE *node) // Pupko's 2002 algorithm
{
  NODETYPE *child;
  double lsum,cl,max, brlen;
  int i,j,k,opt;
  double **PT;
  if (!node) return;
  child=node->firstdesc;
  SIBLOOP(child)
    {
      setupCLmaxopt(child);
    }
//  PT = transitionProb(gS,gR,node->length);
  
brlen = MAX(node->length, MIN_BRANCH_LENGTH);
if (model==4) 
 	PT = transitionProb4(gS,brlen);
if (model==3)
	PT = transitionProb(gS,gS,brlen);
if (model==2)
	PT =transitionProb(gS,gR,brlen);
if (model==1)
	PT = transitionProb2state(gS,gR,brlen);
if (model==0)
	PT = transitionProb2stateSymmetric(gS,brlen);

  
//printf("Transition matrix for node %li\n",node->id);
//mat_print(PT, NSTATES,NSTATES);
  for (i=0;i<NSTATES;i++) // i is the state of this node's parent
	{
  	max=-1.0;
  	for (j=0;j<NSTATES;j++) // j is the state of this node; find max over j
		{
		cl = PT[i][j];

		if (isTip(node)) // This agrees with Cecile's formulation of how to deal with tip states here
			{
			cl *= (node->CL)[j];  // oddly, this special case needs this other value of CL
			}
		else
			{
  			child=node->firstdesc;
			SIBLOOP(child)
  				{
				cl *= (child->CLmax)[j];
				}
			}
		if (cl > max)  // strict inequality means we favor the lowest state if tied
			{
			max=cl;
			opt=j;
			}
		}
	(node->CLmax)[i]=max;
	(node->CLopt)[i]=opt;
	}
  free_dmatrix(PT,0,NSTATES-1,0,NSTATES-1);
  return ;	
}
/**********************/

static void uppassCovarion(NODETYPE *node)
{

// Take the conditional likelihood scores calculated in the Pupko algorithm, add in the root prior
// and recurse from root to tip choosing best state.

// if all root likelihoods are the same, this chooses state 0 by default

  NODETYPE *child;
  int i,opt;
  double max,L;
  child=node->firstdesc;
  if (isRoot(node))
	{
	max=0;
	opt=0;
	for (i=0;i<NSTATES;i++)
		{
		L = rootPrior[i]*(node->CLmax)[i];
		(node->CLmax)[i]=L;
		if (L > max) {opt=i ; max = L; };
		}
	node->opt=opt;
	}
  else
	{
	node->opt=(node->CLopt)[ node->anc->opt ];
	
	
	}
  SIBLOOP(child)
    {
      uppassCovarion(child);
    }
  return ;	
}

/***********************************************************************************/
void printCovarion(NODETYPE *node, int doMarginals)
{
    float diff;
    double sum;
    int i;
    NODETYPE *child;

      if(*(node->taxon_name))
	    printf("%.12s ", node->taxon_name);
      else
	    printf("    Node %li   ", node->id);

	if (!isRoot(node))
		{
		
		printf("(anc:%li)\t",node->anc->id);
		
		}


	sum=0;
	for (i=0;i<NSTATES;i++)
		sum += (node->CLmarg)[i] ;
	// Order of output: CL, CLmax, CLmarg, CLnorm, CLopt
	if (model==2 || model==3)
		{
		printf(": %4.2e %4.2e %4.2e",(node->CL)[0],(node->CL)[1],(node->CL)[2]);
		printf(": %4.2e %4.2e %4.2e",(node->CLmax)[0],(node->CLmax)[1],(node->CLmax)[2]);
		if (doMarginals)
			{
			printf(": %6.2f %6.2f %6.2f",log((node->CLmarg)[0]),log((node->CLmarg)[1]),log((node->CLmarg)[2]));
			if (sum > 0.0)
				printf(": %4.2f %4.2f %4.2f",(node->CLmarg)[0]/sum,(node->CLmarg)[1]/sum,(node->CLmarg)[2]/sum);
			else
				printf(":               ");
			}
		printf(": %2i %2i %2i => %i\n",(node->CLopt)[0],(node->CLopt)[1],(node->CLopt)[2],node->opt);
		}
	if (model==0 || model==1)
		{
		printf(": %4.2e %4.2e",(node->CL)[0],(node->CL)[1]);
		printf(": %4.2e %4.2e",(node->CLmax)[0],(node->CLmax)[1]);
		if (doMarginals)
			{
			printf(": %6.2f %6.2f",log((node->CLmarg)[0]),log((node->CLmarg)[1]));
			if (sum > 0.0)
				printf(": %4.2f %4.2f",(node->CLmarg)[0]/sum,(node->CLmarg)[1]/sum);
			else
				printf(":               ");
			}
		printf(": %2i %2i => %i\n",(node->CLopt)[0],(node->CLopt)[1],node->opt);
		}

    if(!isTip(node))
	{
    	child=node->firstdesc;
    	SIBLOOP(child) 
		printCovarion(child,doMarginals);
	}
    return;
    
}

void printChanges(NODETYPE *node)
{
    float diff;
    NODETYPE *child;

	if (!isRoot(node))
		{
		if (node->opt != node->anc->opt)
			{
			printf ("%i -> %i\t::\t", node->anc->opt,node->opt);
			printf("%li ==> ",node->anc->id);
			if(*(node->taxon_name))
				printf("%s (%li)\n", node->taxon_name, node->id);
			else
				printf("%li\n", node->id);
		
			}
		}
    if(!isTip(node))
		{
			child=node->firstdesc;
			SIBLOOP(child) 
			printChanges(child);
		}
    return;
    
}

static double ** mat_transpose(double **A)
{
double sum;
int i,j,k;
double **p;
p=dmatrix(0,NSTATES-1,0,NSTATES-1);
for (i=0;i<NSTATES;i++)
	{
	for (j=0;j<NSTATES;j++)
		{
		p[i][j]=A[j][i];
		}
	}	
return p;
}


static double ** mat_mult(double **A,double **B)
{
double sum;
int i,j,k;
double **p;
p=dmatrix(0,NSTATES-1,0,NSTATES-1);
for (i=0;i<NSTATES;i++)
	{
	for (j=0;j<NSTATES;j++)
		{
		sum=0;
		for (k=0;k<NSTATES;k++)
			sum += A[i][k]*B[k][j];
		p[i][j]=sum;
		}
	}	
return p;
}

static void mat_print(double **m, int maxi, int maxj)
{
int i,j;
for (i=0;i<maxi;i++)
	{
	printf("%i:\t",i);
	for (j=0;j<maxj;j++)
		printf("%f\t",m[i][j]);
	printf("\n");	
	}
printf("\n");
}

static double ** transitionProb(double s, double r, double t)

{
extern struct NexDataType *gNexDataPtr;
double L1,L2,L3; // eigenvalues
double **V, **VT, **P, **D, cc, cd, c1, c2, L[NSTATES], C[NSTATES];
int i,j;
V=dmatrix(0,NSTATES-1,0,NSTATES-1);
VT=dmatrix(0,NSTATES-1,0,NSTATES-1);
//P=dmatrix(0,NSTATES-1,0,NSTATES-1);
D=dmatrix(0,NSTATES-1,0,NSTATES-1);

// Hijack the branch lengths if called for by the r8s block.
if (gNexDataPtr->RateBlockParms.cov_brlens==1) t=1.0;


cc = 1/sqrt(3);

//if (r != s)
if (fabs (r-s) > (r+s)*0.00001)  // wow, if I check for absolute inequality, I get crazy likelihoods in the case when r and s are just a little off...
	{

	L[0] = 0;
	L[1] = -r - s + sqrt(SQR(r) + SQR(s) - r*s);
	L[2] = -r - s - sqrt(SQR(r) + SQR(s) - r*s);
	
	//printf("eigenvalues: %f %f %f\n", L[0],L[1],L[2]);
	
	// make the eigenmatrix
	for (j=0;j<NSTATES;j++)
		{
		C[j] = sqrt (SQR(s)*SQR(r+L[j]) + SQR(s+L[j])*SQR(r+L[j])+SQR(r)*SQR(s+L[j]));
		for (i=0;i<NSTATES; i++)
			{
			switch (i)
				{
				case 0: V[i][j] = s*(r+L[j])/C[j]; break;
				case 1: V[i][j] = (s+L[j])*(r+L[j])/C[j]; break;
				case 2: V[i][j] = (s+L[j])*r/C[j]; break;
				}
			}
		}
	
	
	
	}
	
else  // r==s

	{
//printf ("Assuming r==s\n");
	L[0] = 0;
	L[1] = -r ;
	L[2] = -3*r ;
	
	V[0][0] = 1/sqrt(3);   V[0][1] = +1/sqrt(2); 	V[0][2] = -1/sqrt(6);
	V[1][0] = 1/sqrt(3);   V[1][1] = 0;  			V[1][2] = +2/sqrt(6);
	V[2][0] = 1/sqrt(3);   V[2][1] = -1/sqrt(2); 	V[2][2] = -1/sqrt(6);
	
	//printf("eigenvalues: %f %f %f\n", L[0],L[1],L[2]);
	
	
	}

//printf("normalizing constants: %f %f %f\n", C[0],C[1],C[2]);

//printf("Eigenmatrix\n");
//mat_print(V, NSTATES,NSTATES);

// diagonal matrix
D[0][0] = 1; D[0][1]=0; D[0][2]=0;
D[1][0] = 0; D[1][1]=exp(t*L[1]); D[1][2]=0;
D[2][0] = 0; D[2][1]=0; D[2][2]=exp(t*L[2]);

//printf("Diagonal matrix\n");
//mat_print(D, NSTATES,NSTATES);

VT= mat_transpose(V);
P = mat_mult(mat_mult(V,D),VT);
free_dmatrix(D,0,NSTATES-1,0,NSTATES-1);
free_dmatrix(V,0,NSTATES-1,0,NSTATES-1);
free_dmatrix(VT,0,NSTATES-1,0,NSTATES-1);
return P;
}

/**********************/

static void setupCL(NODETYPE *node) // conditional likelihoods
{
  extern NODETYPE * gRoot;    /* This global is declared when the whole algorithm is called */
  NODETYPE *child;
  double lsum,cl, **PT, brlen;
  int i,k;
  if (!node) return;
  if (isTip(node)) return;
  child=node->firstdesc;
  SIBLOOP(child)
    {
      setupCL(child);
    }
  for (k=0;k<NSTATES;k++)  // node's state
        	{
        	cl = 1.0;
  			child=node->firstdesc;
  			SIBLOOP(child)
        		{
			brlen = MAX(child->length, MIN_BRANCH_LENGTH);
         		if (model==4) // Just set the two rates equal and let this function do the rest...
 					PT = transitionProb4(gS,brlen);
       			if (model==3) // Just set the two rates equal and let this function do the rest...
 					PT = transitionProb(gS,gS,brlen);
        		if (model==2)
 					PT = transitionProb(gS,gR,brlen);
         		if (model==1)
 					PT = transitionProb2state(gS,gR,brlen);
			if (model==0)
 					PT = transitionProb2stateSymmetric(gS,brlen);
//mat_print(PT,NSTATES,NSTATES);
            	lsum = 0.0;
            	for (i=0;i<NSTATES;i++)  // child's state
                        {
                        //lsum += PT[i][k]*((child->CL)[i]); THAT WAS A BAD BUG
                        lsum += PT[k][i]*((child->CL)[i]);
//printf("Partial Like sums at node %i, child %i : i=%i k=%i PT=%g childCL=%g lsum=%g\n",node->id,child->id,i,k,PT[i][k],(child->CL)[i],lsum);

                        }
            	cl *= lsum;
  				free_dmatrix(PT,0,NSTATES-1,0,NSTATES-1);
            	}
        	(node->CL)[k] = cl;
//printf ("Node %li k=%i CL[k]=%g\n",node->id,k,(node->CL)[k]);
        	}
 //printf ("===>%i  %g %g\n",node->id,(node->CL)[0],(node->CL)[1]);
 //printf("==--==%li %li\n",node,gRoot);
  return ;
}
/**********************/
//static double ** transitionProb2state(double pi0, double mu, double t)
static double ** transitionProb2state(double mu1, double mu2, double t)

{
extern struct NexDataType *gNexDataPtr;
double **P, z, pi1,pi2, scale;

//mu2=mu1;

pi1 = mu2/(mu1+mu2);
pi2 = 1 - pi1;

//scale = 2*mu1*mu2/(mu1 + mu2); // normalizes this so mu will read out in terms of substitutions / site


P=dmatrix(0,1,0,1);

// Hijack the branch lengths if called for by the r8s block.
if (gNexDataPtr->RateBlockParms.cov_brlens==1) t=1.0;

z = exp(-(mu1+mu2) * t );   

P[0][0] = pi1 + (1 - pi1) * z;		P[0][1] = pi2 * (1-z);
P[1][0] = pi1 * (1 - z);			P[1][1] = pi2 + (1 - pi2) * z;

//printf ("...... %f %f %f %f %f %f\n",mu1,mu2,pi1,pi2,z,t);
return P;
}
/***********************************************************************************/
static double ** transitionProb2stateSymmetric(double mu, double t)

{
extern struct NexDataType *gNexDataPtr;
double **P, z, pi0,pi1, scale;

//mu2=mu1;

pi0 = 0.5;
pi1 = 0.5;


P=dmatrix(0,1,0,1);

// Hijack the branch lengths if called for by the r8s block.
if (gNexDataPtr->RateBlockParms.cov_brlens==1) t=1.0;

//scale = 0.5; // 2*1/2*1/2 normalizes this so mu will read out in terms of substitutions / site
z = exp(- 2 * mu * t );   

P[0][0] = pi0 + (1 - pi0) * z;		P[0][1] = pi1 * (1-z);
P[1][0] = pi0 * (1 - z);			P[1][1] = pi1 + (1 - pi1) * z;


return P;
}
/***********************************************************************************/

double objBinaryTrait(double p[])


{
  extern struct NexDataType *gNexDataPtr;	
  extern NODETYPE * gRoot;    /* This global is declared when the whole algorithm is called */

  static int firstTime=1,num_branches;
  int i;
  double obj=-1e100,max=0.0,L,pi;
  NODETYPE *child;
  
 if (p[1] <0 ) {return LARGE_VAL;};
 if (p[2] <0 ) {return LARGE_VAL;};
 
  gS = p[1]; gR = p[2];
  
  
  
/*** Now do conditional likelihoods on tree ***/

  setupCL(gRoot);
 
// Find the best weighted by priors; this will correspond to the max like across the tree (maybe slow
// compared to calculating conditional likelihoods???

	pi = (p[2]/(p[1]+p[2]));
	rootPrior[0]=pi; rootPrior[1]=1-pi;
	obj = rootPrior[0]*(gRoot->CL)[0] + rootPrior[1]*(gRoot->CL)[1];  // Using the stationary freqs as priors
//	obj = .5*(gRoot->CL)[0] + .5*(gRoot->CL)[1];  
	
//printf ("--------------------------------->%i  %g %g %g  %g  %g\n",gRoot->id,(gRoot->CL)[0],(gRoot->CL)[1],p[1],p[2],-obj);
//printCovarion(gRoot);
  return -log(obj); // it's a minimization
}
/**********************/

/***********************************************************************************/

double objBinaryTraitSymmetric(double p[])


{
  extern struct NexDataType *gNexDataPtr;	
  extern NODETYPE * gRoot;    /* This global is declared when the whole algorithm is called */

  static int firstTime=1,num_branches;
  int i;
  double obj=-1e100,max=0.0,L;
  NODETYPE *child;
  
 if (p[1] <0 ) {return LARGE_VAL;};
 
  gS = p[1]; 
  
  
  
/*** Now do conditional likelihoods on tree ***/

  setupCL(gRoot);
 
// Find the best weighted by priors; this will correspond to the max like across the tree (maybe slow
// compared to calculating conditional likelihoods???

	obj = rootPrior[0]*(gRoot->CL)[0] + rootPrior[1]*(gRoot->CL)[1];  
	
//printf ("--------------------------------->%i  %g %g %g  %g  %g\n",gRoot->id,(gRoot->CL)[0],(gRoot->CL)[1],p[1],p[2],-obj);
//printCovarion(gRoot,0); exit(1);
  return -log(obj); // it's a minimization
}
/**********************/

void simulateBinaryCharRecurse(NODETYPE *node, int nodeState, double q01, double q10)

// Simulate a character under the 2-state assym model, starting with some initial state, and storing
// trait in node->opt field.

{
    NODETYPE *child;
    double **PT;
	int child_state;
	node->opt = nodeState; // store the character state in this field
	child=node->firstdesc;
	SIBLOOP(child)
			{
			//if (child->length==FLT_MAX) // this is a legacy ... and stupid indicatory of no length
			// OK interim fix works, but now if you READ a phylogram from file, you need to execute
			// a 'convert_branchlength_to_time' command to get the times for the statements below
			if (0) // this is a legacy ... and stupid indicatory of no length
				PT = transitionProb2state(q01,q10,child->length);
			else
				PT = transitionProb2state(q01,q10,node->time-child->time);
			// Previous is useful for trees generated by r8s tree sim code, which have times
			// but no lengths
			child_state = randBinary(PT[nodeState][1]); // returns a 1 with this probability
//printf("**%f %f %i %i %f %f %f\n",q01,q10,nodeState,child_state,node->time,child->time,PT[nodeState][1]);
			simulateBinaryCharRecurse(child,child_state,q01,q10);
  			free_dmatrix(PT,0,NSTATES-1,0,NSTATES-1);
			}
	return;
    
}

void writeNexusChar(NODETYPE *node)
{
printf("#NEXUS\nBegin data;\ndimensions ntax=%i nchar=1;\nformat symbols=\"01\";\nmatrix\n",numdesc(node));
writeNexusCharRecurse(node);
printf("\n;\nend;\n");

/*
printf("begin trees;\ntree sim = ");
make_parens(node, 0);
printf(";\nend;\n");
*/
}
void writeNexusCharRecurse(NODETYPE *node)
{
	int state,width;
	if (NSTATES>2) // for these models the genotype - phyenotype mapping is more complicated...
		{
		if (node->opt == 2)
			state = 1;
		else
			state = 0;
		}
	else
		state=node->opt;
    NODETYPE *child;
    if (isTip(node))
	 if (*(node->taxon_name)=='\0')
		{
		width = log10(node->id)+1; 
		printf("tx%-*i\t%i\n", width, node->id,state);
		}
	else
    		printf("%s\t%i\n",node->taxon_name,state);
	child=node->firstdesc;
	SIBLOOP(child)
		writeNexusCharRecurse(child);

}

void peak_peek_2(double p[], double (*obj)(double p[]), int dim, double scaleFactor, double ftol)
{
int i,j;
double p1,p2,scale,f;
double pp[3], min=1e20,supposed_optimum;
supposed_optimum = obj(p);

printf("Examining neighborhood of solution by brute force...\nNegative log likelihood; Columns are p[1]; rows are p[2]\n");
for (i=1;i<=dim+1;i++)
	{
	for (j=1;j<=dim+1;j++)
		{
		pp[1] = p[1] * ( 1  +  (j-1-dim/2) * scaleFactor / dim );
		pp[2] = p[2] * ( 1  +  (i-1-dim/2) * scaleFactor / dim );
		f = obj(pp);
		if (f < min ) min = f; 
		printf ("%f ",f); 
		}
	printf("\n");
	}

if (min < supposed_optimum && fabs((min - supposed_optimum)/min) >= ftol)
	printf ("WARNING: BETTER solution found with optimum = %10.8e somewhere in the local search grid (supposed optimum = %10.8e)\n",min, supposed_optimum);
else
	printf ("No better solution found in the local search grid\n");

}


void surface_peek(double (*obj)(double p[]), int dim, double x1low, double x1high, double x2low,double x2high)
{

// calc the obj function on an even grid of dim x dim points
// assuming this is a 2-parm objective

int i,j;
double pp[3],p1,p2,scale1,scale2,f,min=1e20,p1save,p2save;

for (i=1;i<=dim;i++)
	{
	for (j=1;j<=dim;j++)
		{
		if (dim > 1)
			{
			pp[1] = x1low + (j-1) * (x1high - x1low) / (dim-1) ;
			pp[2] = x2low + (i-1) * (x2high - x2low) / (dim-1) ;
			}
		else
			{
			pp[1] = x1low;
			pp[2] = x2low;
			}
		f = obj(pp);
		if (f < min ) { min = f; p1save=pp[1]; p2save=pp[2] ;}

//printf("%f %f : %f\n",pp[1],pp[2],f);

		printf ("%f ",f); 
		}
	printf("\n");
	}
printf ("Best solution found on this grid is at (%f , %f) with obj=%f\n",p1save,p2save,min);

}
static double ** transitionProb4(double s,  double t)

{
extern struct NexDataType *gNexDataPtr;
double L1,L2,L3; // eigenvalues
double **V, **VT, **P, **D, *N, cc, cd, c1, c2, L[NSTATES], C[NSTATES];
int i,j;
V=dmatrix(0,NSTATES-1,0,NSTATES-1);
VT=dmatrix(0,NSTATES-1,0,NSTATES-1);
//P=dmatrix(0,NSTATES-1,0,NSTATES-1);
D=dmatrix(0,NSTATES-1,0,NSTATES-1);

N=dvector(0,NSTATES-1);

// Hijack the branch lengths if called for by the r8s block.
if (gNexDataPtr->RateBlockParms.cov_brlens==1) t=1.0;

#if 1

	L[0] = 0.0;
	L[1] = -2*s ;
	L[2] = -2*s - sqrt(2)*s ;
	L[3] = -2*s + sqrt(2)*s ;
	
	V[0][0] = 1;   V[0][1] = +1; 	V[0][2] = -1;			V[0][3] = -1;
	V[1][0] = 1;   V[1][1] = -1;  	V[1][2] = +1+sqrt(2);	V[1][3] = +1-sqrt(2);
	V[2][0] = 1;   V[2][1] = -1; 	V[2][2] = -1-sqrt(2);	V[2][3] = -1+sqrt(2);
	V[3][0] = 1;   V[3][1] = +1; 	V[3][2] = +1;			V[3][3] = +1;

	N[0]=2;	N[1]=2; N[2] = sqrt(2 + 2*SQR(1+sqrt(2)) ); N[3] = sqrt(2 + 2*SQR(-1+sqrt(2)) );


	for (i=0;i<NSTATES;i++)
		for (j=0;j<NSTATES;j++)
			V[i][j] /= N[j];
	
	//printf("eigenvalues: %f %f %f\n", L[0],L[1],L[2]);
	
	

//printf("normalizing constants: %f %f %f\n", C[0],C[1],C[2]);

//printf("Eigenmatrix\n");
//mat_print(V, NSTATES,NSTATES);

// diagonal matrix

for (i=0;i<NSTATES;i++)
	for (j=0;j<NSTATES;j++)
		{
		if (i==j) 
			D[i][j] = exp(t*L[i]);
		else
			D[i][j] = 0;
		}
//printf("Diagonal matrix\n");
//mat_print(D, NSTATES,NSTATES);






VT= mat_transpose(V);
P = mat_mult(mat_mult(V,D),VT);
#endif

#if 0   // Jukes - Cantor like model
//  HACK FOR NOW ....

P=dmatrix(0,NSTATES-1,0,NSTATES-1); // DELETE WHEN YOU DELETE THIS CODE!
for (i=0;i<NSTATES;i++)
	for (j=0;j<NSTATES;j++)
		{
		if (i==j) 
			P[i][j] = 0.25+0.75*exp(-t*s);
		else
			P[i][j] = 0.25*(1-exp(-t*s));
		}
#endif


free_dmatrix(D,0,NSTATES-1,0,NSTATES-1);
free_dmatrix(V,0,NSTATES-1,0,NSTATES-1);
free_dmatrix(VT,0,NSTATES-1,0,NSTATES-1);
free_dvector(N,0,NSTATES-1);
return P;
}

double objCovarion4(double p[])


{
  extern struct NexDataType *gNexDataPtr;	
  extern NODETYPE * gRoot;    /* This global is declared when the whole algorithm is called */

  static int firstTime=1,num_branches;
  int i;
  double obj=-1e100,max=0.0,L;
  NODETYPE *child;
  
 if (p[1] <0 ) {return LARGE_VAL;};
 
  gS = p[1]; 
  
  
  
/*** Now do conditional likelihoods on tree ***/

  setupCL(gRoot);
 
// Find the best weighted by priors; this will correspond to the max like across the tree (maybe slow
// compared to calculating conditional likelihoods???

  obj=0.0;
  for (i=0;i<NSTATES;i++)
		{
		obj += rootPrior[i]*(gRoot->CL)[i];
		}
	
//printf ("--------------------------------->%i  %g %g %g  %g  %g\n",gRoot->id,(gRoot->CL)[0],(gRoot->CL)[1],p[1],p[2],-obj);
//printCovarion(gRoot);
  return -log(obj); // it's a minimization
}
