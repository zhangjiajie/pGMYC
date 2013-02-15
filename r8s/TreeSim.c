#define BIG_VAL 1e20 
#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#include "TreeSim.h"
#include "MyUtilities.h"
#include "TreeUtils.h"
#include "DistrFuncs.h"
#include "ObjFunc.h"
#include "memory.h"

#define SQR(x)	((x)*(x))
/*#define         drand48()       ((double)rand()/RAND_MAX)*/

/* Note the seed for drand48 is set in the ReadNexusFile caller, doSim */

long name_index;


/* private functions */

static int isOrphan(NODETYPE *node);
static long findArrayElem(float *A, long N, float X);
void connect_nodes(NODETYPE * desc1, NODETYPE * desc2, NODETYPE * anc);
int event(double par,double n);
void normalize_ages(NODETYPE * node, double age);
char * next_name(void);
void setup_auto_vectors(NODETYPE * node, int *index,
			double array1[], double array2[]);

/****************************************************************************/
/****************************************************************************/

/*			RANDOM BRANCH SELECTION ROUTINES		    */

/****************************************************************************/
/****************************************************************************/

void markRandomNodes(TREE Tree, long nMark, NODETYPE ** markedNodes)

/* Mark a random sample of size nMark of the nodes of a tree (without duplication ...
	...and without creating any orphaned grades with no descendants!?)*/
/* Array 'markedNodes' must be allocated to size nMark+1 in calling routine; this maintains
	a list of marked nodes */

{
NODETYPE * node, **nodeArray;
long i,nNodes;
if (markedNodes == NULL)
	fatal("markedNodes array not allocated");
if (nMark > Tree->numTaxa/2)
	doGenericAlert("You may be trying to sample too many nodes on these trees for algorithm");
nodeArray=newAllNodeArray(Tree);
nNodes=numNodes(Tree->root);
unMarkTree(Tree);
for (i=1;i<=nMark;i++)
	{
	node=nextRndNode(nNodes, nodeArray);
/*	do
		{
		node=nextRndNode(nNodes, nodeArray);
		}
		while (isOrphan(node) ); 
printf("NUMDN:%i\n",numUnMarkedDesc(node));*/		
	markNode(node);
	markedNodes[i]=node;
	}
myFree(nodeArray);
return;
}

void RandomBranches(TREE Tree, long nNodes, NODETYPE ** nodeArray, long nMark, NODETYPE ** markedNodes, int withReplace)
{

/* makes a selection of randomly chosen nodes where branches of longer durations are preferentially sampled */
/* Sampling is without replacement--no node is sampled more than once */


float *cumul, T, P,TipDurSum=0.0,TipFract;
long i,j, *count,rawTerms=0,finalTerms=0;
NODETYPE * node;
cumul = (float *)myMalloc(nNodes*sizeof(float)); /* 0-offset array */
count = (long *)myMalloc((nNodes+1)*sizeof(long)); /* 1-offset array */
i=1;
node = nodeArray[i];
if (isRoot(node))
	cumul[i-1] = 0;
else
	cumul[i-1] = node->anc->time - node->time;  
if (isTip(node))
		TipDurSum=(node->anc->time - node->time);
for  (i=2;i<=nNodes;i++)
	{
	node = nodeArray[i];
	if (!isRoot(node))
		cumul[i-1] = cumul[i-2] + (node->anc->time - node->time);  
	if (isTip(node))
		TipDurSum+=(node->anc->time - node->time);
	}
/*
for  (i=1;i<=nNodes;i++)
	printf("Cumul:%li\t%f\n",i,cumul[i-1]);
*/
T = cumul[nNodes-1]; /* Sum of all durations is just in the last bin */
TipFract=TipDurSum/T;
// printf("Fraction of total length in terminal branches, T: %f\t%f\n",TipFract,T);
for (i=1;i<=nNodes;i++) 
	count[i]=0;
unMarkTree(Tree);
for (i=1;i<=nMark;i++)  
	{
	if (withReplace)
		{
		P = myRand() * T; /* get a random variate between 0..T */
		j = findArrayElem(cumul,nNodes,P);
		node=nodeArray[j+1]; /* add one because nodeArray is 1-offset array */
		}
	else
	do 
		{
		P = myRand() * T; /* get a random variate between 0..T */
	/* 	for (j=0;j<=nNodes-1;j++)
			if (P < cumul[j]) break; legacy code known to work; next line is quicker though */
		j = findArrayElem(cumul,nNodes,P);
		node=nodeArray[j+1]; /* add one because nodeArray is 1-offset array */

		} while  (isNodeMarked(node)) ; /* takes care of w/o replacement */
	++count[j+1];
	markNode(node);
	markedNodes[i]=node;
	}


/*
printf("Total count:\n");
for (i=1;i<=nNodes;i++) 
	printf("NODE:%li\t%s\t%i\t%li\n",i,nodeArray[i]->taxon_name,nodeArray[i]->id,count[i]);
*/
return;
}




static int isOrphan(NODETYPE *node)
{
int flag=0;
NODETYPE *n;
unMarkNode(node);
if (numUnMarkedDesc(node) == 0 ) 
	return 1;  // this node is an orphaned node
markNode(node); // temporaily mark this node 
n=node->anc;
unMarkNode(n);
if (numUnMarkedDesc(n) == 0 ) 
		{
		unMarkNode(node);
		markNode(n);
		return 1;  // the ancestor is NOW an orphaned node, so our marked candidate is bad; unmark it and return
		}
else
	{
	unMarkNode(node);
	markNode(n);
	return 0; // the ancestor is not an orphan either, so return its mark and return good
	}
}


/****************************************************************************/
static long findArrayElem(float *A, long N, float X)

/* in array A[0..N-1] which has a sorted, increasing series of values, 
	find the index, J, such that A[J-1] < X < A[J] */ 

{
long j, jlow, jhigh, jtry,jsize;
jlow=0;
jhigh=N-1;
jsize=jhigh-jlow;
do 
	{
	jtry=jlow+jsize/2;
	if (X<A[jtry])
		jhigh=jtry;
	else
		jlow=jtry;
	jsize=jhigh-jlow;
	} while (jsize > 1);
return jhigh;
}


/****************************************************************************/
/****************************************************************************/

/*			CHARACTER EVOLUTION ROUTINES			    */

/****************************************************************************/
/****************************************************************************/

void set_branch_rates(NODETYPE *node, double curRate, double rateChangeAmt,
		double minRate, double maxRate,double transitionProb, int gradual,
		int model)

/* Recursively moves through tree assigning rates to nodes (stored in node->nodeReal).

There are two MODELS:

MODEL=NORMAL	

rates are randomly chosen around a mean of curRate, and with a standard
deviation of rateChangeAmt.  Cannot exceed bounds given by minRate and
maxRate.

MODEL=AUTOCORR

 **NOPE THE FOLLOWING HAS BEEN CHANGED TO DISABLE THE USE OF transitionProb
Begins at root with 'curRate', and then this rate evolves on tree.  With
some fixed probability, 'transitionProb' it changes to another rate. If 
'gradual'==1 this new rate is picked so as to be distributed uniformly on 
[curRate+rateChangeAmt, curRate-rateChangeAmt].  If gradual==0, then a stepwise model 
is used so that EITHER curRate+rateChangeAmt or curRate-rateChangeAmt is chosen
with equal probability.  

Note that rates are not allowed to exceed the bounds specified by minRate 
and maxRate.

 */

{
    NODETYPE *child;
    double newRate,r;
    if (model==1) 	/* normal model */
	{
	newRate=curRate+normdev()*rateChangeAmt;
	if (newRate < minRate)
		newRate=minRate;
	if (newRate > maxRate)
		newRate=maxRate;
	child=node->firstdesc;
	SIBLOOP(child)
		set_branch_rates(child,curRate,rateChangeAmt,
				minRate,maxRate,transitionProb,gradual, model);
	node->estRate=node->nodeReal=newRate;	/* general storage place */
	}

#if 0
	r=myRand();
	newRate=curRate;
	if (r < transitionProb) /* do a transition; i.e., change the rate */
	   {
	   r=myRand();
	
	   newRate+=2*rateChangeAmt*(2*r-1);
	   if(gradual)
	
		newRate+=2*rateChangeAmt*(2*r-1);
	
	
	   else	/*stepwise changes in rate*/
		{
		if (r<0.5)
			newRate+=rateChangeAmt;
		else
			newRate-=rateChangeAmt;
		}
	    if (newRate < minRate)
		    newRate=minRate;
	    if (newRate > maxRate)
		    newRate=maxRate;
	   }
#endif	
	
    if (model==2)
	{
	newRate=curRate;
	node->estRate=node->nodeReal=curRate+normdev()*rateChangeAmt;
	if(node->nodeReal < minRate) node->estRate=node->nodeReal=minRate;	
	if(node->nodeReal > maxRate) node->estRate=node->nodeReal=maxRate;	
	child=node->firstdesc;
	SIBLOOP(child)
		set_branch_rates(child,newRate,rateChangeAmt,
				minRate,maxRate,transitionProb,gradual, model);
	}

    
    
    /*printf("RATES:%f %f %f\n",curRate,newRate,node->nodeReal);*/
    return;


}
void set_branch_lengths(NODETYPE *node, int infinite)

/* Recursively moves through tree assigning branch lengths by looking at the
branch rate,r, stored in nodeReal, getting the duration,d, and generating
a poisson variate with mean r*d.  HOWEVER, if infinite is true, we generate
the EXPECTED branch lengths, rather than a poisson deviate*/
{
NODETYPE *child;
double mu;

if (!isRoot(node))
	{
	if (infinite)
		{  
		if (!isRoot(node))
			node->length = node->nodeReal  *(node->anc->time-node->time);  
			/* assumes rate has been stored! */
		}
	else
		{
		mu=node->nodeReal*(node->anc->time-node->time);
		node->length = poidev(mu);  
/*printf("rate=%f rd=%f length=%f\n",node->nodeReal,mu,node->length);*/
		}
	}
child=node->firstdesc;
SIBLOOP(child)
	set_branch_lengths(child,infinite);
return;


}
/***********************************************************************************/
void set_est_rates(NODETYPE *node,double b,double c,int rateType)

/* b and c are shape parameters of a gamma distribution.  

rateType=1;  LF or NP
rateType=2;  GAMMA

*/

{
	NODETYPE *child;
	double T,k;
	if(!isRoot(node))
	  {
	  T=node->anc->time-node->time;
	  k=node->length;
	  if (rateType==1)
	    node->estRate=k/T;
	  else
		{  
	    node->estRate=b*(k+c)/(1+T*b);
		}
/*printf("%f %f %f %f %f %f\n",b,c,T,k,node->estRate,k/T);*/
	  }
	child=node->firstdesc;
	SIBLOOP(child) 
		set_est_rates(child,b,c,rateType);
	return;
}

/****************************************************************************/
/****************************************************************************/

/*			BIRTH AND DEATH ROUTINES			    */

/****************************************************************************/
/****************************************************************************/
int BDDiversity(long n_taxa, double spec_rate, double extinct_rate,
		double char_rate, int interval)
{

  long MAX_ALLOWABLE_SIZE;
  long i,reps,ix,iy,iz;
  long *lineage;
  long ntaxa;
  long s1,s2;
  double time,xinc;
  long synaps,sz1;
  double dt;
  long num_dts=0;

  MAX_ALLOWABLE_SIZE=3*n_taxa;
  lineage = (long *)myMalloc(MAX_ALLOWABLE_SIZE*sizeof(long));
  dt = 0.1/n_taxa;
/*  NOTES ON SPECIFYING AN APPROPRIATE dt VALUE: ** IMPORTANT **

Evidently, in order for the discrete approximation to the Poisson processes to work
in this simulation, especially the 1-(1-x)^n term in 'event', we have to make dt smaller
as number of taxa, NTAXA, gets bigger.  That's because the probability of a coalescent event
in dt goes up with NTAXA.  However, 'event' only allows one such event in the interval dt.

Similar considerations may hold for the rate parameters, lambda, mu, and chi.  These should
probably not exceed 1.0 or so.
*/

      synaps=0;
      sz1=0;
      time=0.0;
      ntaxa=n_taxa;  /* Note this is the number of taxa at any time
		--NOT the number of taxa with surviving descendants;
		It can therefore hover at 2 toward the root as extinct
		lineages come and go, leaving a long root branch that eventually
		terminates at the real root without leaving any surviving clades
		except one that is nested up further*/

      for (i=0;i<ntaxa;i++)
        lineage[i]=1;

      while (ntaxa>1)
        {
          if (event(dt*spec_rate,ntaxa))
            {
              /*printf("speciation\n");*/
              s1 = (long)(ntaxa*myRand());
              do {s2=(long)(ntaxa*myRand());} while (s1 == s2);
              lineage[s1] += lineage[s2];
              lineage[s2] = lineage[ntaxa-1];
              ntaxa--;
		/*printf("%li %li\n",lineage[s1],ntaxa);*/
            }

          if (event(dt*extinct_rate,ntaxa+1))
            {
              /*          printf("extinction\n");*/
              lineage[ntaxa]=0;
              ntaxa++;
            }

          if (event(dt*char_rate,ntaxa) )
            {
              /*        printf("**synapomorphy\n");*/
              do {s1 = (int)(ntaxa*myRand());} while (0==lineage[s1]);
		/* printf("%li %li\n",lineage[s1],ntaxa);*/
	      if (lineage[s1]<n_taxa) /* don't count full clade size--they're phony
				due to extinct sister group at root */
		      {
	              synaps++;
	              if (1==lineage[s1]) sz1++;
		    /*  ++CladeSizeHisto[lineage[s1]];*/ /* lineage[si] is the clade size 
				at this point; increment the appropriate histo bin;
			N.B.!  We exclude any innovations subtending the root node by
			the if statement.  You might think that stopping at ntaxa =1
			would suffice, but this doesn't prevent innovations from
			accumulating along the branch subtending the root in the case where
			the root's left (or right) descendant is extinct (and hence the 
			reconstructed root is the right descendant of that)*/
		      }
            }

          time += dt;
	++num_dts;
	
	if ( (num_dts/interval) == (num_dts/(double)interval) )
		printf("%f %li\n",time, ntaxa);
        }

printf("tree %d, total time = %f\n",reps,time);

fflush(NULL);

}

/************************************************************************************
 * 
 * DISCRETE TIME Routine for constructing a tree structure according to a 
 * birth-death process backward in time;  i.e., conditional on ntaxa.
 * 
 */

NODETYPE* BDTree(long n_taxa, double spec_rate, double extinct_rate,
		double char_rate)
{
	char *cn;
  long MAX_ALLOWABLE_SIZE;
  long i,reps,ix,iy,iz;
  NODETYPE *(*nodeArray), *root;
  long ntaxa;
  long s1,s2;
  double time,xinc;
  long synaps,sz1;
  double dt;
  long num_dts=0;

  MAX_ALLOWABLE_SIZE=3*n_taxa;
  nodeArray = (NODETYPE **)myMalloc(MAX_ALLOWABLE_SIZE*sizeof(NODETYPE *));
  dt = 0.1/n_taxa;

/*  NOTES ON SPECIFYING AN APPROPRIATE dt VALUE: ** IMPORTANT ** SEE ABOVE */

      name_index=0;
      time=dt;		/* initialize so that time refers 
			to end of the dt increment */
      ntaxa=n_taxa;  /* Note this is the number of taxa at any time
		--NOT the number of taxa with surviving descendants;
		It can therefore hover at 2 toward the root as extinct
		lineages come and go, leaving a long root branch that eventually
		terminates at the real root without leaving any surviving clades
		except one that is nested up further*/

      for (i=0;i<ntaxa;i++)
	{
        nodeArray[i]=newnode(); /* store an array of n_taxa new nodes */
	cn=next_name();
	setNodeName(nodeArray[i],cn);  	
	myFree(cn);
	}

      while (ntaxa>1)
        {
          if (event(dt*spec_rate,ntaxa))
            {
              /*printf("speciation\n");*/
              s1 = (long)(ntaxa*myRand());
              do {s2=(long)(ntaxa*myRand());} while (s1 == s2);
	      nodeArray[s1]=coalesce_nodes(nodeArray[s1],nodeArray[s2]);
	      nodeArray[s1]->time = time;
	      nodeArray[s2]=nodeArray[ntaxa-1]; /*contract array by one */
              ntaxa--;
		/*printf("%li %li\n",lineage[s1],ntaxa);*/
            }

          if (event(dt*extinct_rate,ntaxa+1))
            {
              /*          printf("extinction\n");*/
	      nodeArray[ntaxa]=newnode();
	      nodeArray[ntaxa]->time=time;
		cn=next_name();
		setNodeName(nodeArray[ntaxa],cn);	
		myFree(cn);
              ntaxa++;
            }


          time += dt;
        } /* end while */

root=nodeArray[0];
#if 1
normalize_ages(root,time-dt);/* divide all ages by the age of the root 
	(subtract dt because of time+= statement above at end of loop)*/
#endif
myFree(nodeArray);
return root;	/* This is the root node...the last coalescence */
}
/*************************************************/
NODETYPE* BDback(long n_taxa, double rate, int normalFlag)

/* Simulation of backward Yule (coalescent process) using calls
 * to exponential distribution of waiting times.  If normalFlag ==1 then Root node age is normalized
 * to 1. Also, in that case, the speciation rate does not need to be specified
 */

{
  char *cn;
  long MAX_ALLOWABLE_SIZE;
  long i,reps,ix,iy,iz;
  NODETYPE *(*nodeArray), *root;
  long ntaxa;
  long s1,s2;
  double time,xinc;
  double dt; 

  MAX_ALLOWABLE_SIZE=2*n_taxa;
  nodeArray = (NODETYPE **)myMalloc(MAX_ALLOWABLE_SIZE*sizeof(NODETYPE *));

      name_index=0;
      time=0;		/* initialize so that time refers 
			to end of the dt increment */
      ntaxa=n_taxa;  /* Note this is the number of taxa at any time
		--NOT the number of taxa with surviving descendants;
		It can therefore hover at 2 toward the root as extinct
		lineages come and go, leaving a long root branch that eventually
		terminates at the real root without leaving any surviving clades
		except one that is nested up further*/

      for (i=0;i<ntaxa;i++)
	{
        nodeArray[i]=newnode(); /* store an array of n_taxa new nodes */
	cn=next_name();
	setNodeName(nodeArray[i],cn);  	
	myFree(cn);
	}

      while (ntaxa>1)
        {
	  s1 = (long)(ntaxa*myRand());
	  do {s2=(long)(ntaxa*myRand());} while (s1 == s2);
	  nodeArray[s1]=coalesce_nodes(nodeArray[s1],nodeArray[s2]);
	  dt = hexp(ntaxa*rate);
          time += dt;
	  nodeArray[s1]->time = time;
	  nodeArray[s2]=nodeArray[ntaxa-1]; /*contract array by one */
	  ntaxa--;
	    /*printf("%li %li\n",lineage[s1],ntaxa);*/
        } /* end while */

root=nodeArray[0];

if (normalFlag)
	normalize_ages(root,time);/* divide all ages by the age of the root */
myFree(nodeArray);
return root;	/* This is the root node...the last coalescence */
}


/*************************************************/
NODETYPE* Yule_C(long n_taxa, double speciation)

/* Simulation of backward Yule (coalescent process) using calls
 * to birthDist distribution, which provides a stable age structure to nodes 
 */

{
  char *cn;
  long MAX_ALLOWABLE_SIZE;
  long i,reps,ix,iy,iz;
  NODETYPE *(*nodeArray), *root;
  long ntaxa;
  long s1,s2;
  double time,xinc;
  double dt, rate=1.0; 
  double *times;

  ntaxa=n_taxa;
  MAX_ALLOWABLE_SIZE=2*n_taxa;
  nodeArray = (NODETYPE **)myMalloc(MAX_ALLOWABLE_SIZE*sizeof(NODETYPE *));
  times = (double*)myMalloc((ntaxa-2)*sizeof(double));

/*  NOTES ON SPECIFYING AN APPROPRIATE dt VALUE: ** IMPORTANT ** SEE ABOVE */

      name_index=0;
      time=0;		/* initialize so that time refers 
			to end of the dt increment */

      for (i=0;i<ntaxa;i++)
	{
        nodeArray[i]=newnode(); /* create n_taxa new leaf nodes */
	cn=next_name();
	setNodeName(nodeArray[i],cn);  	
	myFree(cn);
	}
      for (i=0;i<ntaxa-2;i++)
	{
        times[i]=birthDist(speciation); /* draw a random internal time*/
// printf("Time=%f\n",times[i]);
	}

      qsort((void *)times,ntaxa-2,sizeof (double),compar); 
		/* sorts the time in increasing order from present back */

      ix=0;
      while (ntaxa>1)
        {
	  s1 = (long)(ntaxa*myRand());
	  do {s2=(long)(ntaxa*myRand());} while (s1 == s2);
	  nodeArray[s1]=coalesce_nodes(nodeArray[s1],nodeArray[s2]);
// Nex two lines seem like bogus legacies
	  dt = hexp(ntaxa*rate);
          time += dt;
	  nodeArray[s1]->time = times[ix];
	  nodeArray[s2]=nodeArray[ntaxa-1]; /*contract array by one */
	  ntaxa--;
	  ++ix;
	    /*printf("%li %li\n",lineage[s1],ntaxa);*/
        } /* end while */

root=nodeArray[0];
root->time=1.0;
myFree(nodeArray);
myFree(times);
return root;	/* This is the root node...the last coalescence */
}

double PH_gamma (long n, double *times, double T)
/* gamma of Pybus-Harvey:
 times[] is 0-offset and is ordered from recent to past on the n-1 internal nodes, opposite their sense...
 d = is the difference between ordered node ages
 T = age of root

NB! If n = 2, we return 0, otherwise we get a NAN. THIS HAS IMPLICATIONS, and the original paper of PH does not
consider this boundary condition.

 */
{
long i,j,k,ll;
double A,B,C,gamma,d;
if (n == 2) return 0.0;
B=sqrt(1.0/(12*(n-2)));
A=0;
for (j=2;j<=n;j++)
	{
	if (j==2)
		d=T-times[n-3];
	else if (j==n)
		d=times[0];
	else
		d=times[n-j]-times[n-j-1];
	A+=j*d;
	}
C=0;
for (i=2;i<=n-1;i++)
	for (k=2;k<=i;k++)
		{
		if (k==2)
			d=T-times[n-3];
		else if (k==n)
			d=times[0];
		else
			d=times[n-k]-times[n-k-1];
		C+=k*d;
		}

//printf ("%f %f %f\n",A,B,C);

C/=(n-2);
gamma = (C-A/2)/(A*B);
return gamma;
}

/*************************************************/
NODETYPE* RY_1997(long n_taxa, double T, double speciation, double extinction, double sampling)
/*

Simulation of Rannala Yang (1997) kernel function to generate trees based on birth-death-sampling
rates.

 */

{
  char *cn;
  long MAX_ALLOWABLE_SIZE;
  long i,reps,ix,iy,iz;
  NODETYPE *(*nodeArray), *root;
  long ntaxa;
  long s1,s2;
  double time,xinc;
  double dt, rate=1.0, gamma; 
  double *times;

  ntaxa=n_taxa;
  MAX_ALLOWABLE_SIZE=2*n_taxa;
  nodeArray = (NODETYPE **)myMalloc(MAX_ALLOWABLE_SIZE*sizeof(NODETYPE *));
  times = (double*)myMalloc((ntaxa-2)*sizeof(double));

/*  NOTES ON SPECIFYING AN APPROPRIATE dt VALUE: ** IMPORTANT ** SEE ABOVE */

      name_index=0;
      time=0;		/* initialize so that time refers 
			to end of the dt increment */

      for (i=0;i<ntaxa;i++)
	{
        nodeArray[i]=newnode(); /* create n_taxa new leaf nodes */
	cn=next_name();
	setNodeName(nodeArray[i],cn);  	
	myFree(cn);
	}
      for (i=0;i<ntaxa-2;i++)
	{
        times[i]=T*RY_1997_Dist(speciation,extinction,sampling); /* draw a random internal time scaled on 0..1
		and scale by T, the age of the clade */
// printf("Time=%f\n",times[i]);
	}

      qsort((void *)times,ntaxa-2,sizeof (double),compar); 
		/* sorts the time in increasing order from present back */
	
	  gamma = PH_gamma(ntaxa,times,T);
	  printf("[Pybus-Harvey gamma statistic for this clade is %f]\n",gamma);

      ix=0;
      while (ntaxa>1)
        {
	  s1 = (long)(ntaxa*myRand());
	  do {s2=(long)(ntaxa*myRand());} while (s1 == s2);
	  nodeArray[s1]=coalesce_nodes(nodeArray[s1],nodeArray[s2]);
	  nodeArray[s1]->time = times[ix];
	  nodeArray[s2]=nodeArray[ntaxa-1]; /*contract array by one */
	  ntaxa--;
	  ++ix;
	    /*printf("%li %li\n",lineage[s1],ntaxa);*/
        } /* end while */

root=nodeArray[0];
//root->time=1.0;
root->time=T;
myFree(nodeArray);
myFree(times);
return root;	/* This is the root node...the last coalescence */
}

/************************************************************************************
 * 
 * DISCRETE TIME Routine for constructing a tree structure according to a birth-death process
 * forward in time.  Interval in time is T.  Time is measured from 0 = present to T = root time.
 * 
 */

NODETYPE* BDTreeForward(double T, double spec_rate, double extinct_rate,
		double char_rate)
{
	char *cn;
  long MAX_ALLOWABLE_SIZE;
  long i,j, reps,ix,iy,iz, num_intervals;
  NODETYPE *(*nodeArray), *root,  *temp_internal,*n;
  long ntaxa;
  long s1,s2, nlineage;
  double time=0.0,xinc;
  long synaps,sz1;
  double dt;
  long num_dts=0;

  MAX_ALLOWABLE_SIZE=50*exp((spec_rate+extinct_rate)*T);  /* just a guess */
  nodeArray = (NODETYPE **)myMalloc(MAX_ALLOWABLE_SIZE*sizeof(NODETYPE *));
  dt = 0.001/exp((spec_rate-extinct_rate)*T); /*  NOTES ON SPECIFYING AN APPROPRIATE dt VALUE: ** IMPORTANT ** SEE ABOVE */
  num_intervals=T/dt;
/*   printf("T=%f dt=%f num_intervals=%li\n", T, dt, num_intervals);*/
  time=T;
  root=newnode();   /* make the root node and its first two descendants*/
  root->time=time;
  nodeArray[0]=newnode();

if (1)
  {
  nodeArray[1]=newnode();
  connect_nodes(nodeArray[0], nodeArray[1], root);  
  nlineage=2;
printf ("blech\n");
  }
else
  {
n=nodeArray[0];
n->anc=root;
root->firstdesc=n;
nlineage=1;
  }

  for (i=0;i<num_intervals;i++)
    {
    time-=dt;
          if (event(dt*extinct_rate,nlineage))
            {
             /* printf("extinction: time=%f\n", time);*/
              s1 = (long)(nlineage*myRand()); /* returns a rand on 0..nlineage-1; i.e. a proper index into nodeArray?*/
	      nodeArray[s1]->time = time;
	      for (j=s1;j<nlineage-1;j++)
		nodeArray[j]=nodeArray[j+1]; /*contract array by one */
              nlineage--;
            }

          if (event(dt*spec_rate,nlineage))
            {
             /* printf("speciation: time=%f\n", time);*/
              s1 = (long)(nlineage*myRand()); 
	      nodeArray[s1]->time = time;
	      temp_internal=nodeArray[s1];
	      for (j=s1;j<nlineage-1;j++)
		nodeArray[j]=nodeArray[j+1]; /*contract array by one */
              nlineage--;

	      if (nlineage+2 > MAX_ALLOWABLE_SIZE)
		fatal("Too many taxa generated in BDforward\n");
	      nodeArray[nlineage]=newnode();
	      nodeArray[nlineage+1]=newnode();
	      connect_nodes(nodeArray[nlineage],nodeArray[nlineage+1],temp_internal );
	      nlineage+=2;
            }
    
    
    }
  for (i=0;i<nlineage;i++)
	      nodeArray[i]->time = time;
    

myFree(nodeArray);
return root;	/* This is the root node...the last coalescence */
}



int event(double par,double n)
{
  double r;
  double P = 1.0 - pow(1.0-par,n);
  r=myRand();
/*  printf("param=%f n=%f test=%f rand=%f\n",par, n,  P, r);*/
  if (r<P) return 1;
  return 0;
}
void connect_nodes(NODETYPE * desc1, NODETYPE * desc2, NODETYPE * anc)
{
NODETYPE * ancestor;
desc1->anc=anc;
desc2->anc=anc;
anc->firstdesc=desc1;
desc1->sib=desc2;
return;
}


NODETYPE * coalesce_nodes(NODETYPE * node1, NODETYPE * node2)
{
NODETYPE * ancestor;
ancestor = newnode();
node1->anc=ancestor;
node2->anc=ancestor;
ancestor->firstdesc=node1;
node1->sib=node2;
return ancestor;
}

char * next_name(void)
{
char  taxon_name[10], *c;
++name_index;
sprintf(taxon_name,"t%li",name_index);
c=DupStr(taxon_name);
return c;
}

void normalize_ages(NODETYPE * node, double age)
{
NODETYPE *child;
node->time/=age;
child=node->firstdesc;
SIBLOOP(child)
	normalize_ages(child,age);
return;
}

/****************************************************************************/
/****************************************************************************/
/*			MISCELLANY					    */
/****************************************************************************/
/****************************************************************************/

double angle(double *vec1, double *vec2, int arraySize)
{
double ang, size1=0.0, size2=0.0, cross=0.0;
int i;
for (i=0;i<arraySize;i++)
	{
	size1+=vec1[i]*vec1[i];
	size2+=vec2[i]*vec2[i];
	cross+=vec1[i]*vec2[i];
	}
size1=sqrt(size1);
size2=sqrt(size2);
ang = acos(cross/(size1*size2));
return ang;
}
double euclid_distance(double *vec1, double *vec2, int arraySize)
{
#define SMALL 0.01
double ang, size=0.0;
int i;
for (i=0;i<arraySize;i++)
	{
	if (fabs(vec1[i]+vec2[i]) < SMALL)
		--arraySize; /* ignore info for nodes that are very nearly zero
			as their percent errors will tend to be bogus */
	else
		size+=2.0*fabs(vec1[i]-vec2[i])/(vec1[i]+vec2[i]);
	}

return size/arraySize;  /* average percent time difference */
}

double correlation(double x[], double y[], unsigned long n)
{
        unsigned long j;
        double yt,xt,r;
        double syy=0.0,sxy=0.0,sxx=0.0,ay=0.0,ax=0.0;

        for (j=1;j<=n;j++) {
                ax += x[j];
                ay += y[j];
        }
        ax /= n;
        ay /= n;
        for (j=1;j<=n;j++) {
                xt=x[j]-ax;
                yt=y[j]-ay;
                sxx += xt*xt;
                syy += yt*yt;
                sxy += xt*yt;
        }
        r=sxy/sqrt(sxx*syy);
return r;
}


void setup_auto_vectors(NODETYPE * node, int *index,
			double array1[], double array2[])
{
NODETYPE *child;

if(isRoot(node))
	*index=0; /* arrays are 0-offset so start indexing at 1 */
else
	if (!isTip(node))
		{
		child=node->firstdesc;
		SIBLOOP(child)
			{
			array1[*index]=node->nodeReal;
			array2[*index]=child->nodeReal;
			++(*index);
			}
		}
child=node->firstdesc;
SIBLOOP(child)
	setup_auto_vectors(child,index,array1,array2);
return;
}
double tree_auto_correlation(NODETYPE * root)
{
int k,nb, index;
double *array1, *array2,r;
k=numdesc(root);
nb=2*k-2;	/* maximum num of branches in tree, and array size for autocor*/
array1=(double *)myMalloc(nb*sizeof(double));
array2=(double *)myMalloc(nb*sizeof(double));
setup_auto_vectors(root,&index,array1,array2);
r=correlation(array1-1,array2-1,(unsigned long)index); /* sub 1 for 1-offsets*/
return r;
}

/***********************************************************************************/
double BD_Like(double params[])

// IGNORE following comments. This code now merely does a pure birth likelihood!

/* Calculates the **log**        likelihood under a birth death model according to Nee et al.'s
 * equation 21 in their 1994 "Reconstructed evolutionary process" paper. 
 * NO!! Now we use the kernel function below to calculate based on waiting times.  Either way should
 * be the same.
 * 
 * params[2] = a = mu/lambda
 * params[1] = r = lambda-mu
 * 
 * Remember this is a one-offset array.  
 *
 * 
 */

{
    double BDkernel_func(double a, double r, double xn, double xnwait, int n);
    extern int gNtips, gBDnumvar,gStemFlag; // gStemFlag=1 means we will include the branch subtending the root node in calcs.
    extern double *gTimes; // array of times for each node in a binary tree (n-1 of these) including root, increasing order with
							// root the last element (see more comments under function defn
    extern double gLogNm1, gRootDur;
    extern double gB;
    extern NODETYPE* gRoot;
	
    double sum=0.0, a, r, s=0.0, prod=1.0, xnwait, xn,p, like, D;
    int i, N, n ;
    r=params[1];
    if (gBDnumvar>1)
        a=params[2];
    else
	a=0.0;	/*pure birth model*/
    
   
#if 0

    if (a<=0.0 || r<=0.0 || a>=1.0)
	return BIG_VAL; /* corresponding to log(0) case (and accounting for need
			    to negate for MinND */
     N=gNtips-2;

    for (i=0;i<=gNtips-3;i++) 
	{
	s+= gTimes[i];
	printf("[1] time=%f, sum=%f\n", gTimes[i], s);
	}
    for (i=0;i<=gNtips-2;i++) /* just the internal node times */
	{
	sum-= 2*log(exp(r*gTimes[i])-a);
	printf("[2] time=%f, sum=%f\n", gTimes[i], sum);
	}
	
    sum+=N*log(r)+r*s+gNtips*log(1-a) + log (24.);
    return -sum; /* negate to send to minimize function */
#endif

#if 0
     N=gNtips-2;

    s=0.0;
    for (i=0;i<=gNtips-3;i++) /* all internal ages except root's (should be outside of function) */
	{
	s+= gTimes[i];
	printf("[1] gNtips=%i time=%f, sum=%e\n", gNtips, gTimes[i], s);
	}
	
    prod=1.0;
    for (i=0;i<=gNtips-2;i++) /* just the internal node times */
	{
	p= SQR(exp(r*gTimes[i])-a);
	prod*=p;
	printf("[2] time=%f, p=%e prod=%e\n", gTimes[i], p, prod);
	}

    like=pow(r, N) * exp(r*s)*pow(1-a, gNtips) /prod;

like=N*log(r)+r*s+gNtips*log(1-a)-log(prod);

    printf("a=%f r=%f r^n=%e exp=%e 1-a^n=%e prod=%e like=%e\n",a, r, pow(r, N), exp(r*s), pow(1-a, gNtips), prod, like );
    return -like; /* negate to send to minimize function */
#endif

#if 0
    if (a<0.0 || r<=0.0 || a>=1.0)
	return BIG_VAL;
    N=gNtips-2;
    sum=0.0;
    for (i=1;i<=gNtips-2;i++) 
	{
	xn=gTimes[i];
	xnwait=xn-gTimes[i-1];
	n=gNtips-i;
//printf ("In BDLike: %i %f %f %f %f %i\n",gNtips,a,r,xn,xnwait,n);
	p=BDkernel_func(a, r, xn, xnwait, n);
	sum+=p;

	if (gStemFlag==1)
		{
		sum+=BDkernel_func(a,r,gTimes[gNtips-1],gRootDur,1); // add the term for branch subtending the root of this tree
		}

/*	printf("a=%f r=%f xn=%f xnwait=%f n=%i p=%f prod=%f YuleProd=%f \n",
	    a, r, xn, xnwait, n,  p, sum,n*r*exp(-n*r*xnwait) );*/
	}
#endif	

#if 1
// whole thing is a stupid brute force optimization, when we could do it analytically, but useful 
// as check...

    if (a<0.0 || r<=0.0 || a>=1.0)
	return BIG_VAL;
	D = get_sum_durations(gRoot); // stupid to recalc: should cache somewhere
	sum = (gNtips-2)*log(r)-r*D;

	if (gStemFlag==1)
		{
		sum += log(r)-r*(gRoot->anc->time-gRoot->time);
			//printf ("D=%f sum=%f r=%f rootancT=%f rootT=%f stemFlag=%i\n",D,sum,r,gRoot->anc->time,gRoot->time,gStemFlag);	
		}
	//else
			//printf ("D=%f sum=%f r=%f rootT=%f stemFlag=%i\n",D,sum,r,gRoot->time,gStemFlag);	
#endif	
	
    return -sum;
}

double YuleLike(NODETYPE* root, int stemFlag, double *speciation)

// analytic results for likelihood and ML estimate of rate in Yule
// Handles the case of a single branch descended from the root

/* A technical issue on the LR test for sister group diversities. The likelihood for a stem or crown clade is 

		L(n) = (n-1)! * rate ^ (n-n0) * exp(-rate*sum_dur).

The LR test is  

		L(n1)*L(n2) /L(n1+n2).

It seems like the factorial coefficients might be something like (n1-1)!(n2-1)!/(n1+n2-1)!, which would be very different from 1.
However, in the denominator, we are ALSO looking only at diversification processes that have produced two sister clades with 
exactly n1 and n2 species, so the correct coefficient for the denominator is also (n1-1)!(n2-1)!. Actually, it seems like there 
should be an additional factor of 2 in both numerator and denominator, for the two events of n1,n2 or n2,n1 (if n1 ne n2).
In any case, the coefficients drop out, and we can actually calculate the YuleLike without them, at least for this sole purpose of
doing an LR test. Note that the function does NOT calculate the coefficients.
*/

{
double D, r, logLike;
long n,n0;
n=numdesc(root);
D=get_sum_durations(root);
if (stemFlag) //stem
	{
	n0=1;
	D += root->anc->time - root->time;
	}
else  //crown
	n0=2;
r = (n-n0)/D;
*speciation = r;
if (  (stemFlag==1 && n==1 ) || (stemFlag==0 && n==2) )// Model: single terminal branch or crown clade with n=2!
	logLike = 0; // because likelihood is 1 that a rate of 0 will generate observed data!
else
	logLike = (n-n0)*log(r)-r*D;
return logLike;	
}

double BDkernel_func(double a, double r, double xn, double xnwait, int n)
{

double w;
 
/* a=0.0;*//************!!!!!!!!!!!!!!!!!!!*****************/

//printf ("In kernel: %f %f %f %f %i\n",a,r,xn,xnwait,n);

if (xnwait > xn)
	{
    printf("ERROR in BDkernel_func: waiting time too large\n");
	printf ("In kernel: %f %f %f %f %i\n",a,r,xn,xnwait,n);
	}

else
    {
  /*  w=n*r*exp(-n*r*xnwait)*pow(1-a*exp(-r*(xn-xnwait)), n-1)/pow(1-a*exp(-r*xn), n);*/
    w=log(n)+log(r)-n*r*xnwait+(n-1)*log(1-a*exp(-r*(xn-xnwait)))-n*log(1-a*exp(-r*xn));
#if 0
   printf("***[p =  %e X %e / %e = %e]\n",n*r*exp(-n*r*xnwait), 
	    pow(1-a*exp(-r*(xn-xnwait)), n-1), pow(1-a*exp(-r*xn), n) , w);
#endif 
    return w;  
    }
    
}
/***************************************************/
// Just generate the clade size, not the tree
/** starting from a split ***/

double Yule_forward(double rate,  double T, double *sum_durations, int stemFlag)
{
double curTime=0.0, wait_time;
double N;
if (stemFlag==1) N=1;
else N=2;
*sum_durations=0.0;
while(1)
    {
    wait_time=hexp(N*rate);
    curTime+=wait_time;
    if (curTime>T)
	{
	(*sum_durations)+=N*(T-(curTime-wait_time));
	break;
	}
    (*sum_durations)+=N*wait_time;
    N+=1.0;
    }  
return N;
    
}

/************************************************************************************/

// Build a random tree from a split such that its basal two sister groups are generate with
// rates spec1 and spec2 with a forward Yule process over time T.

NODETYPE* SisterGroupYule(double T, double spec1, double spec2, double *Ntips1, double *Ntips2, double *sum_durations)
{
NODETYPE *r1, *r2, *first1, *first2;
double sum_durations1, sum_durations2;

r1 = YuleTreeForward(T, spec1, Ntips1, &sum_durations1,1);
r2 = YuleTreeForward(T, spec2, Ntips2, &sum_durations2,1);

// make r1 the root and get rid of r2, updating node pointers...

first1 = r1->firstdesc;
first2 = r2->firstdesc;
first1->sib = first2;
first2->anc=r1;

*sum_durations = sum_durations1 + sum_durations2;
return r1;
}


/************************************************************************************
 * 
 * Routine for constructing a tree structure according to a Yule/pure birth process
 * forward in time.  Interval in time is T.  Time is measured from 0 = present to T = root time.
 * INITIAL TREE HAS A ROOT NODE AND ITS TWO IMMEDIATE DESCENDANTS! THUS, THIS IS A SIMULATION
 * OF A CROWN CLADE 
 */

NODETYPE* YuleTreeForward(double T, double spec_rate, double *Ntips, double *sum_durations,int stemFlag)
{
	char *cn;
  long MAX_ALLOWABLE_SIZE;
  long i,j, reps,ix,iy,iz, num_intervals;
  NODETYPE *(*nodeArray), *root,  *temp_internal, *n;
  long ntaxa;
  long s1,s2, nlineage;
  double time=0.0,xinc, wait_time;
int stem=0;

  MAX_ALLOWABLE_SIZE=50*exp((spec_rate)*T);  /* just a guess */
  nodeArray = (NODETYPE **)myMalloc(MAX_ALLOWABLE_SIZE*sizeof(NODETYPE *));
  if (!nodeArray) fatal ("Couldn't allocate enough memory in BDforward\n");
/*   printf("T=%f dt=%f num_intervals=%li\n", T, dt, num_intervals);*/
  time=T;
  root=newnode();   /* make the root node and its first two descendants*/
  root->time=time;
  nodeArray[0]=newnode();
if (!stemFlag)  // for crown group simulation
  {
  nodeArray[1]=newnode();
  connect_nodes(nodeArray[0], nodeArray[1], root);  
  nlineage=2;
  }
else // for stem group simulation
  {
n=nodeArray[0];
n->anc=root;
root->firstdesc=n;
nlineage=1;
  }

*sum_durations=0.0;
  while(1)
    {
             /* printf("speciation: time=%f\n", time);*/
    wait_time=hexp(nlineage*spec_rate);
    time-=wait_time;
    if (time<0.0)
	{
	(*sum_durations)+=nlineage*(time+wait_time);
	time=0.0;
	break;
	}
    (*sum_durations)+=nlineage*wait_time;
    s1 = (long)(nlineage*myRand());  // pick one of the candidate nodes at random
    nodeArray[s1]->time = time;	// it will have this split time
    temp_internal=nodeArray[s1]; // save it
    for (j=s1;j<nlineage-1;j++) // delete it from candidate nodes
        nodeArray[j]=nodeArray[j+1]; /*contract array by one */
    nlineage--;
    
    if (nlineage+2 > MAX_ALLOWABLE_SIZE)
    fatal("Too many taxa generated in BDforward\n");
    nodeArray[nlineage]=newnode(); // make two new nodes and add to tree
    nodeArray[nlineage+1]=newnode();
    connect_nodes(nodeArray[nlineage],nodeArray[nlineage+1],temp_internal );
    nlineage+=2;
    }
  for (i=0;i<nlineage;i++)
	      nodeArray[i]->time = time;
*Ntips=nlineage;    

myFree(nodeArray);
return root;	/* This is the root node...the last coalescence */
}
