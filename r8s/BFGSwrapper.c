r8s/BFGSwrapper.h                                                                                   0000644 0000766 0000120 00000000254 07565204055 013547  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "TimeAlgorithms.h"
int BFGSwrapper
	(
	int numvar,
	double x[],
	double 		(*objective)(double []),
	void		(*gradient)(double [], double []),
	double		*max_obj
	);
                                                                                                                                                                                                                                                                                                                                                    r8s/ConstrOpt.c                                                                                     0000644 0000766 0000120 00000016061 10357561057 013357  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /*   General nonlinear minimization algorithm with or without constraints.
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

                                                                                                                                                                                                                                                                                                                                                                                                                                                                               r8s/ConstrOpt.h                                                                                     0000644 0000766 0000120 00000000550 07565204055 013357  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "nexus.h"
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

                                                                                                                                                        r8s/DistrFuncs.c                                                                                    0000644 0000766 0000024 00000036635 11210541233 013516  0                                                                                                    ustar   sandermj                        staff                                                                                                                                                                                                                  #include "DistrFuncs.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "memory.h"
#include "NRCvectorUtils.h"
#include "MyUtilities.h"

#define AUTO_SEED_PROMPT 0	/* set to 1 if myRand() will always prompt for seed */

long iix;	/* seed */


/*****************************************************************************************************/
double RY_1997_Dist(double speciation, double extinction, double sampling) 

{
double U,phi,y,d;
U=myRand();
d=exp (extinction - speciation);

phi = (sampling*speciation*(d-1) + (extinction-speciation)*d)/(d-1);

y=(log(phi-U*sampling*speciation)-log(phi-U*sampling*speciation+U*(speciation-extinction)))/(extinction-speciation);


return y;
}
/*****************************************************************************************************/
double birthDist(double lambda) 

/* returns a random number from a birth process 
waiting time with fixed interval 1.0. Times are measured from 0 = present
to 1=root; cf Ross, Stochastic Processes, p. 145, for the density.  The
equation below is then the inverse function of the cdf, obtained by integration.
THen it is subtracted from 1.0 to get the time sense right! 

It's easy to derive that X = (lambda + log(y-y*exp(-lambda)+exp(-lambda))/lambda;
Now just right the first lambda as log(exp(lambda)) and combine terms in the 
numerator.  This reduces to what's below.

NB. Watch out for the overflow case of lambda very large

*/
{
double r;
r=myRand();
if (lambda > 20) // this is a quick first test
		{ // if r*exp(lambda) >> 1 then r*(exp(lambda)-1)+1 reduces to just r*exp(lambda)
		  // Therefore below we test if log(r)+lambda is big, which is less prone to overflow
		if (log(r)+lambda > 20) // this is a "proper" test if we'll overflow
			return -log(r)/lambda; // easier calculation for overflow case
		}
return 1-(log(r*(exp(lambda)-1)+1)/lambda);
}
/*****************************************************************************************************/

double hgeom(double param)

	/*..............returns a random geometric variate distributed with parameter
	param, according to c.d.f. 1-(1-param)^n  ......(returns a double inc ase of BIG numbers)*/

{
double z, den;


if (param<1.0e-8)  den = (-param) - (param*param)/2;  /*taylor series to avoid roundoff error when
					subtracting a very small param from 1 and then taking logarithm:
					log (1 + x) = x - x^2/2 +x^3/3 + ....        */
					
else den=log(1-param);

z= log (1-(double)myRand())  / den;





return (floor(z)  +1.0 );	
	/* is this the right truncation of the real-valued expression above? YES
	Checked by reference to the expected mean of the distribution in 100,000 replicates
	EX = 1/param   Var = (1-param)/param^2  See Olkin, Gleser, and Derman, p. 193ff. Probability
	Models and Applications, 1980.*/
}

/*****************************************************************************************************/
double hexp(double param)

	/*..............returns a random exponential variate distributed with parameter
	param, according to c.d.f. 1-exp(-param*x) 
	Function has been validated by checking that the median deviate is equal to log 2/param, and
	the mean deviate is 1/param */

{
double z, den;


z= -log (1-myRand())  / param;


return (z);	

}
#define M 714025
#define IA 1366
#define IC 150889      		 
/***********************************************************************/
double myRand(void)
{
/* procedure returns a real random value on [0,1]. The variable iix must be
declared globally as a long integer by the main program and must be
initialized as a seed integer   */

extern long iix;
static long iy,ir[98];
static int iff=0;
int j;



/********/

return rand()/(double)RAND_MAX;


/*******/


#if AUTO_SEED_PROMPT
static int flag=1;
if (flag)
	{
	getseed();
	flag=0;
	}	/* Initialize for random number generator*/
#endif
/*return((double)rand()/RAND_MAX);*/

if (iix<0 || iff==0) {
	iff=1;
	if((iix=(IC-(iix)) % M) < 0) iix= -iix;
	for (j=1;j<=97;j++) {
		iix=(IA*(iix)+IC) % M;
		ir[j]=iix;
	}
	iix=(IA*iix+IC) % M;
	iy=iix;
   }
j= 1 + 97.0*iy/M;
if (j>97 || j < 1) printf("error in rand\n");
iy=ir[j];
iix=(IA*iix+IC) % M;
ir[j]=iix;
return (  (double)iy/M);


}
/*****************************************************************************************************/
void getseed(void)
{
	printf("Please type a seed for the random number generator\n");
	printf("XXXXX\n");
	scanf("%5li", &iix);
	return;

}
/*****************************************************************************************************/
long rndLong(long maxLong)  /* rand int on [1..maxLong] */
{
long i;
i=myRand()*maxLong+1;
if (i>maxLong)
	return maxLong;
else
	return i;


}

/*****************************************************************************************************/
void bshuf(int *targetArray, int *excludeArray, long numChars, long includedChars)

/* Generates a bootstrap weight array, 'targetArray', which can be used by various routines.
THis array has 'numChars' sites in it.  Some of these sites may be excluded.  This information
is provided in 'excludeArray', which MUST be present.  By default this array (also of length 
'numChars') has all 1's in it.  Any site can be excluded from bootstrapping by setting to zero.
The actual number of non-zero sites is 'includedChars' <= 'numChars'.
 */
{
long ix, choice,validCount=0;
for (ix=0;ix<numChars;ix++)
	{
	targetArray[ix]=0;
	}
for (ix=0;ix<includedChars;ix++) /* we only want this many replicates--which is possibly less than numChars
					when some chars have been excluded */
	{
	for (;;)
		{
		choice = rndLong(numChars);
		if (excludeArray[choice-1]>0) break;
		}
	++targetArray[choice-1];
	}
/*for (ix=0;ix<N;ix++)
	printf("%i %i\n",targetArray[ix],excludeArray[ix]);*/
return;
}
/*****************************************************************************************************/
void bshuf2(int *targetArray, long numChars)

/* fills a supplied bootstrap weight array, 'targetArray', which can be used by various routines.
THis array has 'numChars' sites in it. 
 */
{
long ix, choice,validCount=0;
for (ix=0;ix<numChars;ix++)
	{
	targetArray[ix]=0;
	}
for (ix=0;ix<numChars;ix++) 
	{
	choice = rndLong(numChars);
	++targetArray[choice-1];
	}
return;
}
/*****************************************************************************************************/
void bshuf3(float *weightArray, long numChars,int nreps,  char * buffer1,char * buffer2)

/* Sets up a series of bootstrap arrays, targetArray (addree supplied by user), by multinomial sampling
from weights in real-valued weightArray.  We do this setting up a cumulative distribution function
for the multinomial weights and then throwing random numbers at the ordinate.  Note that weightArray
should have weights in the range 0..numChars, which below we normalize to probabilities by dividing by numChars
 */
 
 /* buffer1 is printed for the first replicate, buffer2 for all the remaining replicates
  * 
  */
 
{
int *targetArray,k;
long ix,j,choice,validCount=0;
float *pCumul,*pMean, p;
targetArray=(int *)myMalloc((numChars)*sizeof(int));
pCumul=(float *)myMalloc((numChars)*sizeof(float));
pMean=(float *)myMalloc((numChars)*sizeof(float));/* Just to check if algorithm OK */
pCumul[0]= weightArray[0]/numChars;
for (ix=1;ix<numChars;ix++) /* set up cumul dist function */
	{
	pCumul[ix]=pCumul[ix-1]+weightArray[ix]/numChars;
	pMean[ix]=0.0;
	/*printf("%i %f\n", ix, pCumul[ix]);*/
	}
printf("[Weighted Bootstrap Resamples:]\n");
for (k=1;k<=nreps;k++)
    {
    printf("begin paup;\n");
    for (ix=0;ix<numChars;ix++) 
	    {
	    targetArray[ix]=0;
	    }
    for (ix=0;ix<numChars;ix++) 
	    {
	    p = myRand(); /* rand on [0,1] */
	    for (j=0;j<numChars;j++)
		    if (p<=pCumul[j])
			    break;
	    ++targetArray[j];
	    }
    printf("weights ");
    for (j=0;j<numChars-1;j++)
	    {
	    if ((j>0)&& ((j/20)==(j/20.0)))
		printf("\n");
	    printf("%i:%i, ", targetArray[j],j+1);
	    }
    printf("%i:%i;\n",  targetArray[j],numChars);
    printf("end;\n");
    if (k==1)
	{
	if (buffer1)
	    printf("%s\n", buffer1);
/* printf("hs;contree all/majrule=yes strict=no file=the_middle_trees replace=yes append=no;");*/
	}
    else
        if (buffer2)
	    printf("%s\n", buffer2);
    printf("\n");
   for (j=0;j<numChars;j++)
	    {
	    pMean[j]+=targetArray[j]/(float)nreps;
	    }
    }
printf("[Mean vector\n");
for (j=0;j<numChars;j++)
	{
	if ((j>0)&& ((j/10)==(j/10.0)))
	    printf("\n");
	printf("%f:%i, ", pMean[j],j+1);
	}
printf("]\n");
myFree(pMean);
myFree(pCumul);
myFree(targetArray);
return;
}
/*****************************************************************************************************/

/* returns a standard normal deviate */

float normdev(void)
{
        static int iset=0;
        static float gset;
        float fac,rsq,v1,v2;

        if  (iset == 0) {
                do {
                        v1=2.0*myRand()-1.0;
                        v2=2.0*myRand()-1.0;
                        rsq=v1*v1+v2*v2;
                } while (rsq >= 1.0 || rsq == 0.0);
                fac=sqrt(-2.0*log(rsq)/rsq);
                gset=v1*fac;
                iset=1;
                return v2*fac;
        } else {
                iset=0;
                return gset;
        }
}
/*****************************************************************************************************/

#define PI 3.141592654

/* returns a poisson deviate with mean xm */

double poidev(double xm)
{
        double gammln(double xx);
        static double sq,alxm,g,oldm=(-1.0);
        double em,t,y;

        if (xm < 12.0) {
                if (xm != oldm) {
                        oldm=xm;
                        g=exp(-xm);
                }
                em = -1;
                t=1.0;
                do {
                        ++em;
                        t *= myRand();
                } while (t > g);
        } else {
                if (xm != oldm) {
                        oldm=xm;
                        sq=sqrt(2.0*xm);
                        alxm=log(xm);
                        g=xm*alxm-gammln(xm+1.0);
                }
                do {
                        do {
                                y=tan(PI*myRand());
                                em=sq*y+xm;
                        } while (em < 0.0);
                        em=floor(em);
                        t=0.9*(1.0+y*y)*exp(em*alxm-gammln(em+1.0)-g);
                } while (myRand() > t);
        }
        return em;
}
#undef PI
double gammln(double xx)
{
        double x,y,tmp,ser;
        static double cof[6]={76.18009172947146,-86.50532032941677,
                24.01409824083091,-1.231739572450155,
                0.1208650973866179e-2,-0.5395239384953e-5};
        int j;

        y=x=xx;
        tmp=x+5.5;
        tmp -= (x+0.5)*log(tmp);
        ser=1.000000000190015;
        for (j=0;j<=5;j++) ser += cof[j]/++y;
        return -tmp+log(2.5066282746310005*ser/x);
}
float factrl(int n)
{
        static int ntop=4;
        static float a[33]={1.0,1.0,2.0,6.0,24.0};
        int j;

        if (n < 0) nrerror("Negative factorial in routine factrl");
        if (n > 32) return exp(gammln(n+1.0));
        while (ntop<n) {
                j=ntop++;
                a[ntop]=a[j]*ntop;
        }
        return a[n];
}

/*********************************************************************/
int * taxon_sample_simple(int numtaxa, int numrandom)
{
int n,i,j, *temp,*sample,rand_ix ;
temp=(int *)myMalloc(numtaxa*sizeof(int));
sample=(int *)myMalloc(numrandom*sizeof(int));
for (j=0;j<numtaxa;j++)
	temp[j]=j+1; /* this is a 0-offset vector of all the valid ids to be sampled here. These
			    are basically every id from 1..numtaxa */
n=numtaxa;
for (i=0;i<numrandom;i++) 
	{
	rand_ix=rndLong((long)(n))-1;
	sample[i]=temp[rand_ix];
	remove_array_item(temp, n, sample[i]);
	--n;
	}

return sample;
}

/*********************************************************************/
void taxon_sample(int numtaxa, int numfixed, int numrandom, 
		int fixed[], int sample[],
		int nstart,int nstop,int nstart2,int nstop2,int numrandom2)
/* numtaxa = total number of taxa in analysis (PAUP ntaxa)
   numfixed = number of taxa that will be kept constant across ALL replicates
   numrandom = number of randomly chosen taxa in each replicate
   fixed[]= integer array of ids (on 1..numtaxa) of the fixed taxa
   sample[]= integer array of ids (on 0..numtaxa-1) of the randomly selected taxa PLUS
	the fixed taxa!  i.e.,  the final chosen sample
   
    DAMMIT,  let's keep ALL ids on [1..n]!! But all arrays are 0-offset

   Both the arrays are passed to the function and must therefore be initialized
   correctly
   
   ** This function now operates in three modes
	(1) If nstart==0, the routine samples from ALL of the ids from 1..numtaxa, excluding
		the fixed taxa.    The number of taxa sampled is numrandom
	(2) if nstart>0, it samples only from the set of id's ranging from [nstart..nstop], and
		from these it samples exactly numtaxa (numtaxa<=nstop-nstart+1)
	(3) if nstart2>0, it samples from both ranges [nstart..nstop] and [nstart2..nstop2], 
		sampling numtaxa from the first range and numtaxa2 from the second.
	NOTE! These n* parameters are all on conventional [1..n] range
	NOTE! Should initialize fixed[] to all zeros in the case where no fixed list present; then
	    nothing will happen!
	NOTE! Specifying a fixed taxon within the range of [nstart..nstop] will lead to unpredictable events
*/		    

{

    int i,j,k,n,   nsample, temp[MAX_TAXON_ARRAY], temp2[MAX_TAXON_ARRAY], rand_ix;
    nsample=numrandom+numrandom2+numfixed; /* final number of taxa in the sample */
    if ( nsample > MAX_TAXON_ARRAY || numtaxa > MAX_TAXON_ARRAY)
	fatal("Array bounds exceeded in taxon_sample");
    for (i=0;i<numfixed;i++)
	sample[i]=fixed[i]; /* add the fixed taxa to the final sample list */


/* set up the temp list for case (1) above */

    if (nstart==0)
	{
    	for (j=0;j<numtaxa;j++)
		temp[j]=j+1; /* this is a 0-offset vector of all the valid ids to be sampled here. These
			    are basically every id from 1..numtaxa */
	}

/* BUT now have to remove the fixed taxa id's from the temp list, because
we don't want to consider sampling from them (they're already there!) */

    n=numtaxa;
    for (j=0;j<numfixed;j++)
	{
	remove_array_item(temp, n, fixed[j]); 
	--n;
	}
/*printf("temp vector for case(1)\n");
  for (k=0;k<n;k++)
	    printf("%i %i\n", k, temp[k]);*/
	

/* set up the temp list for case (2) above */

    if (nstart>0)
	for (j=0;j<nstop-nstart+1;j++)
		temp[j]=nstart+j;

    if (nstart2>0)
	{
	for (j=0;j<nstop2-nstart2+1;j++)
		temp2[j]=nstart2+j;
	/*printf("temp vector for case(3)\n");
	for (k=0;k<nstop2-nstart2+1;k++)
		    printf("%i %i\n", k, temp2[k]);*/
	}

	    
    if (nstart==0)
	n=numtaxa-numfixed;
    else
	n=nstop-nstart+1;
	
    for (i=0;i<numrandom;i++) 
	{
	rand_ix=rndLong((long)(n))-1;
	/*printf("rand=%i\n", rand_ix);*/
	sample[i+numfixed]=temp[rand_ix];
	remove_array_item(temp, n, sample[i+numfixed]);
	--n;
	}

if (nstart2>0)
    {
    n=nstop2-nstart2+1;
    for (i=0;i<numrandom2;i++) 
	{
	rand_ix=rndLong((long)(n))-1;
	/*printf("rand=%i\n", rand_ix);*/
	sample[i+numfixed+numrandom]=temp2[rand_ix];
	remove_array_item(temp2, n, sample[i+numfixed+numrandom]);
	--n;
	}
    }
	
	
 /*   printf("[The taxon sample:\n");
    for (j=0;j<nsample;j++)
	printf("%i %i\n", j, sample[j]);*/
    return;
    
    
}
void remove_array_item(int array[], int num_elements, int item)

/* works on 0-offset arrays ! If item is not found, nothing happens*/

{
    int i,ix,  last;
    last=num_elements-1;
    for (ix=0;ix<=last;ix++)
	if (array[ix]==item)
	    {
	    for (i=ix;i<last;i++)
		array[i]=array[i+1];
	    }
    return;
}

                                                                                                   r8s/DistrFuncs.h                                                                                    0000644 0000766 0000024 00000001504 11210535760 013516  0                                                                                                    ustar   sandermj                        staff                                                                                                                                                                                                                  #define MAX_TAXON_ARRAY 500

double RY_1997_Dist(double speciation, double extinction, double sampling) ;
double birthDist(double lambda);
float normdev(void);
float factrl(int n);
double gammln(double xx);
double poidev(double xm);
long rndLong(long maxLong);
double hgeom(double param);
double hexp(double param);
double myRand(void);
void getseed(void);
void bshuf3(float *weightArray, long numChars,int nreps,  char * buffer1, 
	char * buffer2);
void bshuf2(int *targetArray, long numChars);
void bshuf(int *targetArray, int *excludeArray, 
	long numChars, long includedChars);
void taxon_sample(int numtaxa, int numfixed, int numrandom, 
		int fixed[], int sample[],
		int nstart,int nstop,int nstart2,int stop2,int numrandom2);
int * taxon_sample_simple(int, int);
void remove_array_item(int array[], int length, int item_index);
                                                                                                                                                                                            r8s/DrawTree.c                                                                                      0000644 0000766 0000120 00000026573 11643634325 013150  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "TreeUtils.h"
#include "DrawTree.h"
#include "memory.h"
#include <math.h>
#include "MyUtilities.h"

/* private functions */

void assignX(NODETYPE *node,  int Xleft,  int Xright, int Xwidth, int treeMode);
void Assign_XY_Tree(NODETYPE *root,  int xUpLeft,int yUpLeft, int xLowRight,
			int yLowRight,int treemode);
void  MakeTree(NODETYPE* root,	 int xUpLeft, int yUpLeft, int xLowRight, 
		int yLowRight, int treeWidth, int nameWidth, int treemode);
int 		assignY2(NODETYPE *node,  double *YcurPtr,  double yinc);
void		Tprint(NODETYPE *r),
		HorizLine( int x1,  int x2,  int y),
		VertLine( int y1,  int y2,  int x),
		DrawIt(NODETYPE *node),
		DrawNames(NODETYPE *node, int treemode),
		swap( int *x1,  int *x2);
void		MaxTaxLength(NODETYPE *node);
char*		TreeToString(void);

/* GLOBALS */

static double gAge;
int gHyp,gMax=0;
double gMaxToTipLength, gMaxToTipLengthRate;
char matrx[MAXHEIGHT][MAXWIDTH];

int		WINWIDTH,WINHEIGHT;

NODETYPE*	internalList[MAXINTERNALNODENAMES];	/* fixed linear array of pointers to nodes */
int		internalListix=0;
/***********************************************************************************/

void DrawTree(NODETYPE *root, int treemode, int userMaxWidth)


/*

	treemode	type of tree
	   0		   cladogram
	   1		   phylogram
	   2 		   chronogram
	   4		   ratogram
	   9		   cladogram with taxon names replaced with character state in column 1!
	   10		   phylogram with taxon names replaced with character state in column 1!
*/

{


	int numtax,treeWidth,nameWidth;

	numtax=numdesc(root);
	gAge=root->time; /* used in chronogram draws */
	if ((gAge==0.0) && (treemode==2))
	    fatal("Times do not appear to have been set on trees; try converting lengths to time");
	gMax=0;		/* width of longest taxon name */
	MaxTaxLength(root);	  /* put the width of the longest taxon name in gMax */
	gMax+=2; 	/* Add two characters for possible special characters '%' */
	gMax+=6; 	/* Add 6 characters for possible terminal character states */
	nameWidth=gMax;
	gMaxToTipLength=calcMaxToTip(root);
	gMaxToTipLengthRate=calcMaxToTipRate(root);
if (userMaxWidth !=0)
	treeWidth=userMaxWidth;
else
	treeWidth=max(max(2*numtax,MINWIDTH),userMaxWidth);
	WINWIDTH=treeWidth+nameWidth; /* window is big enough to a name area of size gMax and a tree area that
						might be as small as MINWIDTH but is larger for big trees */ 


	if (WINWIDTH > MAXWIDTH)
		{
		WINWIDTH = MAXWIDTH - nameWidth -1;
		treeWidth = WINWIDTH - nameWidth;
		printf("Tree is so large that it has been compressed horizontally and resolution may be lost\n");
		}
	WINHEIGHT=min(2*numtax-1,MAXHEIGHT);

	if (treemode == 9 || treemode == 10) printf("Currently printing marginal probs that switch is ON at each node\n");
	MakeTree(root,0,0,WINWIDTH-1,WINHEIGHT-1,treeWidth,nameWidth,treemode);

}
/***********************************************************************************/

double calcMaxToTip(NODETYPE* node)

/* Calculates maximum distance from root to tip when lengths are available */

{
	double max=0.0,temp,thisLength;
	NODETYPE *child;
	if (!node) return(0.0);

	if (isRoot(node))
		{
		thisLength=0.0;
		}
	else
		{
		thisLength=node->length;	/* don't add length under the root */
		}
	if (isTip(node)) 
			{
			return (thisLength);  /* treat a tip and a compact node the same way */
			}
	child=node->firstdesc;
	SIBLOOP(child) {
			temp=calcMaxToTip(child);
			if (temp > max) max = temp;
			}
	return thisLength+max;
}
/***********************************************************************************/

double calcMaxToTipRate(NODETYPE* node)

/* Calculates maximum distance from root to tip when lengths are available */

{
	double max=0.0,temp,thisLength;
	NODETYPE *child;
	if (!node) return(0.0);

	if (isRoot(node))
		{
		thisLength=0.0;
		}
	else
		{
		thisLength=node->estRate;	/* don't add length under the root */
		}
	if (isTip(node)) 
			{
			return (thisLength);  /* treat a tip and a compact node the same way */
			}
	child=node->firstdesc;
	SIBLOOP(child) {
			temp=calcMaxToTipRate(child);
			if (temp > max) max = temp;
			}
	return thisLength+max;
}
/************************************************************************************/


void  MakeTree(NODETYPE* root,	 int xUpLeft, int yUpLeft, int xLowRight, 
		int yLowRight, int treeWidth, int nameWidth, int treemode)
										 
	/*  Takes a string tree description and plots the branches
	of it in a box-- NOW ONLY USED FOR DUMB TERMINAL DRAWING */										 
										 
{
  
	NODETYPE *ax,*bx;
	int x,y, windowWidth,treeAreaWidth;
	char* s;	


	windowWidth=xLowRight-xUpLeft+1; /* character or pixel width of window */
	Assign_XY_Tree(root,xUpLeft,yUpLeft,xUpLeft+treeWidth-1,yLowRight,treemode);	
	for (y=0;y<MAXHEIGHT;y++)
		for (x=0;x<MAXWIDTH;x++)  matrx[y][x]=SPACE;

	DrawIt(root);
	DrawNames(root,treemode);
	s=TreeToString();
	printf("%s\n",s);
	myFree(s);
	return;


}  
/***********************************************************************************/
char* TreeToString(void)
{
	char *s;
	int width[MAXHEIGHT],
		TotChars,
		NChars,
		row,
		col;
	char *ptr;
	
	TotChars=0;
	for (row=0;row<=WINHEIGHT-1;row++)
		{
		NChars=WINWIDTH;
		for (col=WINWIDTH-1;col>=0;col--,NChars--)
			{
			if (matrx[row][col] != SPACE) break;
			}
		width[row]=NChars;
		TotChars+=(NChars+1); /* the extra '1' is for the CR added on each line
								4D WANTS TO SEE CR ONLY, NOT CR-LF */
		}

	s=(char*)myMalloc(sizeof(char)*(TotChars+1));
	if (!s) return NULL;


	ptr=s;
	
	for (row=0;row<=WINHEIGHT-1;row++)
		{
		for (col=0;col<=width[row]-1;col++)
			{
			*ptr=matrx[row][col];
		/*	if ((*ptr == 10) || (*ptr == 13)) don't store extraneous
						LF's or CR's that may have
						crept in via makegroup  
				*ptr=SPACE;  instead replace with space --FIXED */
			++ptr;
			}
	/*	*ptr=RETURN;
		++ptr;*/
		*ptr=LF;
		++ptr;	/* UNIX and WWW recognizes LF as CR-LF */
		}

	*ptr=0;

return (s);
}


/***********************************************************************************/
void Assign_XY_Tree(NODETYPE *root,  int xUpLeft,int yUpLeft, int xLowRight,
						 int yLowRight, int treemode	)
						  
/*	Assigns X and Y display coordinates to the nodes in a tree structure.  Uses
	the  integer coordinates of the upper left and lower right of the
	box in which the tree should be displayed.  The X,Y values are stored in the
	tree structure */


{
		 int N,yinc;
		 double yupleft;
		 yupleft=yUpLeft;
		if(root==NULL) 
						return;
			
		(void)maxorder(root);
		(void)numdesc(root);
		assignX(root,xUpLeft,xLowRight,xLowRight-xUpLeft+1,treemode);
		N=root->numdesc;
		if (N==1) yinc=0;
		else yinc = (yLowRight-yUpLeft)/(float)(N-1);
		assignY2(root,&yupleft,yinc);						
		return;
}
/***********************************************************************************/

void MaxTaxLength(NODETYPE *node)

	/* Finds the length of the longest string contained in the tree structure */
{
	NODETYPE *child;
	extern int gMax;
	int temp, length, max=0;
	if (!node) return;
	child=node->firstdesc;
	SIBLOOP(child) 
		{
		if (*(child->taxon_name) == '\0')
			length=4; /* use max likely number of digits in the id
					instead of the absent string */
		else
			length=strlen(child->taxon_name);
		if (length>gMax) 
			gMax=length;
		MaxTaxLength(child);
		}
	return;


}
/***********************************************************************************/
int assignY2(NODETYPE *node,  double *YcurPtr,  double yinc)
{
	NODETYPE *child;
	int sum=0,count=0;
	
	if (!node) return 0;
	if(isTip(node)  || (node->isCompactNode) )
		{
		node->Y= *YcurPtr;
		(*YcurPtr)+=yinc;
		return(node->Y);
		}

	child=node->firstdesc;
	
	SIBLOOP(child) {
		sum+=assignY2(child,YcurPtr,yinc);
		++count;
		}
	sum/=count;
	node->Y=sum;
	return(sum);
}
/***********************************************************************************/
void assignX(NODETYPE *node,  int Xleft,  int Xright, int Xwidth, int treeMode)
{
	if (isRoot(node)) {
			node->X = Xleft;
			}
	else
			switch (treeMode)
			{
			case 9:
			case 0:
				node->X = Xleft + (Xright - Xleft-1)/(float)(node->order + 1);
				break;
			case 10:
			case 1:
				node->X = Xleft + (Xwidth-1)*node->length/gMaxToTipLength;
				break;
			case 2:
				node->X = Xleft + (Xwidth-1)*(gAge-node->time)/gAge;
				break;
			
			case 4:
				node->X = Xleft + (Xwidth-1)*node->estRate/gMaxToTipLengthRate;
				break;
			
			}


	if (node->sib) assignX(node->sib,Xleft,Xright,Xwidth,treeMode);
	
	if (node->isCompactNode) return;  

	switch (treeMode)
		{
		case 2:
		if (node->firstdesc) 
			assignX(node->firstdesc,Xleft,Xright,Xwidth,treeMode);
		break;
		default:
		if (node->firstdesc) 
			assignX(node->firstdesc,node->X,Xright,Xwidth,treeMode);
			
		}
	return;
}
/***********************************************************************************/

void DrawIt(NODETYPE *node)
{
	NODETYPE *child;
	 int x,count=0;
	 char *name;
	if (!node) return;

/* 'draw' the taxon name for the node */

	child=node->firstdesc;
	SIBLOOP(child) {  /* if node is already a tip, this will just fall through */
		HorizLine(node->X,child->X,child->Y);
			/* horiz line from the node's X to its child */
		VertLine(node->Y,child->Y,node->X);
			/* vert line from the node's Y to meet the previous horiz  line */
		DrawIt(child);
		}
	return;
}

/*********************************************************************/
#define ccheck(c) (isalpha(c) || ((c)=='%'))
void DrawNames(NODETYPE *node, int treemode)
{
	extern int gLabel;
	NODETYPE *child;
	 int x,room,slength,offset,count=0,code=0,width;
	 char *name, *nameFromShell, s[15],ss[20],c;
	 double sum, switch_on;
	if (!node) return;

/* 'draw' the taxon name for the node */
	if (treemode==9 || treemode==10)  // this is the oddball case of printing the character reconstruction instead of name
		{
		if (isTip(node))
			{
			if (  (node->CL)[2]	== 1.0) c='+';
			else						c='-';
//			sprintf(s, "%*i (%c)", 1,node->opt,c);
//			name=s;
			snprintf(ss, 20,"%*i (%c) %.10s", 1,node->opt,c,node->taxon_name);
// had no end of trouble with abort traps until I switched over snprintf!
			name=ss;
			
			}
		else
			{
			sum = (node->CLmarg)[0] + (node->CLmarg)[1] + (node->CLmarg)[2];
			switch_on =     1 - (node->CLmarg)[0]/sum ;
			snprintf(s, 15, "%1i %4.2f (%li)", node->opt, switch_on, node->id);
			name=s;
			}
		}
	else
		name=node->taxon_name;
	if (gLabel)
	 if (!(*name)) /* if an internal node name doesn't exist ...*/
	    {
	    name=s;  /*...then set up a temp string using fixed array */
	    if(node->id == 0)
		width=1;
	    else
	        width=log10(node->id)+1;/* get the number of digits in the id */
	    sprintf(name, "%*li", width,node->id); /* just use its id number */
	    }
	slength=strlen(name)+2; /* this is how much room we need */
	x=(node->X)+1;  /* start writing  one character to the right of x */
	while(*name)   /* loop through the string */
		  {
			if (count>gMax) break;
		  	matrx[node->Y][x]=*name;
		  	++count;
		  	++x;
		  	++name;
		  }
	child=node->firstdesc;
	SIBLOOP(child) {  /* if node is already a tip, this will just fall through */
		DrawNames(child,treemode);
		}
	return;
}

/*********************************************************************/


void HorizLine( int x1,  int x2,  int y)
{
 int i;
if (x1>x2) swap(&x1,&x2);
for (i=x1;i<=x2;i++) matrx[y][i]=DASH;
return;
}
void VertLine( int y1,  int y2,  int x)
{
 int i;
if (y1>y2) swap(&y1,&y2);
for (i=y1+1;i<=y2-1;i++) matrx[i][x]=BAR;
matrx[y1][x]=PLUS;
matrx[y2][x]=PLUS;
return;
}
void swap( int *x1,  int *x2)
{
 int temp;
temp=*x1;
*x1=*x2;
*x2=temp;
return;
}

                                                                                                                                     r8s/DrawTree.h                                                                                      0000644 0000766 0000120 00000001773 07565204055 013151  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #define MAXINSTRING 5000	/*max size of tree description input string */
#define MAXfSTRING 100
#define MAXINTERNALNODENAMES 100	/* needed for a global fixed length array */
#define LF 10
#define RETURN 13
#define COLON ':'
#define BAR	'|'
#define PLUS '+'
#define DASH '-'
#define	SPACE	' '
#define COMMA	','
#define	RIGHTPARENS	')'
#define	LEFTPARENS	'('
#define	isvalidtaxchar(c)	isalnum(c) || (ispunct(c) && ((c)!=COMMA) && ((c)!=RIGHTPARENS) && ((c)!=LEFTPARENS) && ((c)!=SPACE))
#define min(a,b)			( (a)<=(b) )  ? (a):(b)
#define max(a,b)			( (a)>=(b) )  ? (a):(b)
#define SIBLOOP(c)			for (; (c); (c)=(c)->sib)
#define isTip(c)			 ( (c)->firstdesc == NULL )
#define MINWIDTH	20	/* this is the minimum width in the window 
						allowed for the tree itself (i.e, minus the taxon names) */ 
#define MAXWIDTH 150
#define MAXHEIGHT 500



/* STRUCTURES AND PROTOTYPES */




double calcMaxToTip(NODETYPE* node);
double calcMaxToTipRate(NODETYPE* node);
void DrawTree(NODETYPE *root, int treemode,int userMaxWidth);


     r8s/Maximize.h                                                                                      0000644 0000766 0000120 00000000623 07565204055 013210  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "TreeUtils.h"
double Min1D(double (*func)(double x),double *xmin,
						double ax,double bx, double cx,double ftol );
void plot2d(double low1, double high1, double low2, double high2, int gridSize);
double MinND(TREE t,int method, int algorithm, double (*func)(double p[]),void (*grad)(double [], double []),
		double p[], int numvar,int *iter, double ftol,
		double linMinDelta,int *success);
                                                                                                             r8s/MinND.c                                                                                         0000644 0000766 0000120 00000004270 11217270432 012356  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include <stdio.h>
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


                                                                                                                                                                                                                                                                                                                                        r8s/MyUtilities.c                                                                                   0000644 0000766 0000120 00000012251 07565204055 013701  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "MyUtilities.h"
#include "memory.h"

/****  Miscellaneous utility commands ****/




char * slurpFile (FILE * inFileStream, long maxSize)

{
        char    *BigBuffer;
        int             c;
        long    count=0,i=0;

        
        
        BigBuffer=(char*)malloc(maxSize*sizeof(char));
        if (!BigBuffer) 
                {
                doGenericAlert("Could not allocate file buffer");
                return NULL;            
                }

        while ((c=fgetc(inFileStream)) != EOF)  /* have to define c as int so that EOF can be detected*/
                {
                if (i >= maxSize-1) /* have to save room for terminating null */
                        {
                        doGenericAlert("Slurped file exceeds allotted maximum");
                        return NULL;
                        }
                BigBuffer[i]=c;
                ++i;
                
                }
                BigBuffer[i]='\0';

return BigBuffer;       
                
}



void fatal(char *s)
{
	long i;
	printf("!FATAL ERROR!");
	doGenericAlert(s);
	for (i=1;i<100000;i++);
	exit(1);
	return;
}

void doGenericAlert(char* errorMsg)
{
char* s;
printf("********************* WARNING **********************\n\n");
printf("%s\n" ,errorMsg);
printf("\n****************************************************\n");
return;
}

void strtoupper(char *s)
{
	char *temps;
	temps=s; 
	while(  *temps ) 
		{
		*temps=toupper(*temps); 
		++temps;
		}
	/* converts string to upper case */
	return;
}
/***************************************************/

void concat(char **destHndl, char *source)

/* concatenates source into destination; have to be careful, because 'realloc' might
create a new destination pointer and we have to make sure to copy that to the old
'destHndl' */

{
char *tempPtr, *destPtr;
long lengthDest,lengthSource, length;

destPtr=*destHndl;

lengthDest=strlen(destPtr);
lengthSource=strlen(source);
length=lengthDest+lengthSource+1;

tempPtr=(char *)myRealloc(destPtr,(length*sizeof(char))); /* "myRealloc" */
if (tempPtr==NULL) 
	fatal("myReallocation error in concat");
if (tempPtr != destPtr) /* myRealloc had to move pointer */
	{
	destPtr=tempPtr;
	*destHndl=destPtr;	/* make sure to save the new pointer */
	}
strcat(destPtr,source);
return;
}

/***************************************************/
		
char*	DupStr(char* s)
{

/* makes a dynamic memory copy of a string -- returns NULL on error*/

char* sNew;
sNew=(char*)myMalloc((strlen(s)+1)*sizeof(char));
if (sNew != NULL)
	strcpy(sNew,s);
return sNew;

}
/***************************************************/

char*	makeEmptyString(void)
{
char* s;
s=(char*)myMalloc(sizeof(char));
if (s != NULL)
	*s = '\0';
return  s;
}
/***************************************************/
FILE* PromptFileName(char* promptMsg, char* mode)

{
FILE* fpntr;
char fnIn[FILENAME_MAX];	/* defined in stdio.h */
printf("%s",promptMsg);
scanf("%s",fnIn);
if (  (fpntr=fopen(fnIn,mode)) )
	return fpntr;
else
	fatal("Error in file handling");

}
/***************************************************/
int isStrInteger(char* s)

/* Checks to see if a string represents an arbitrary length integer number */

{
char * p;
p=s;
while (*p)
	{
	if (!isdigit(*p))
		return 0;
	++p;
	}
return 1;
}

/***************************************************/
#define numX	100
#define numY	60

void dumbPlot(double X[],  double Y[], int N)

/* X and Y are 0-offset arrays */


{
    double Xmax, Ymax, Xmin, Ymin, Xdif, Ydif, Xintv, Yintv;
    char m[numX+1][numY+1];
    int ix, iy, Xa, Ya;
    for (ix=0;ix<numX+1;ix++)
	for (iy=0;iy<numY+1;iy++)
	    m[ix][iy]=' ';
    array_minmax(X, N, &Xmin, &Xmax);
    array_minmax(Y, N, &Ymin, &Ymax);
    Xdif=Xmax-Xmin;
    Ydif=Ymax-Ymin;
    Xintv=Xdif/numX;
    Yintv=Ydif/numY;

    printf("Ascii Plot of %i Points\n\n", N);

#if 0
    for (ix=0;ix<N;ix++)
	{
	printf("%f\t%f\n", X[ix], Y[ix]);
	}
#endif
    for (ix=0;ix<N;ix++)
	{
	    Xa=(X[ix]-Xmin)/Xintv;
	    Ya=(Y[ix]-Ymin)/Yintv;
	    m[Xa][Ya]='*';
	}
    
    for (iy=numY;iy>=0;iy--)
	{
	if (iy==numY)
	    printf("%6.2f", Ymax);
	if (iy==0)
	    printf("%6.2f", Ymin);
	printf("\t|");
	for (ix=0;ix<=numX;ix++)
	   printf("%c",  m[ix][iy]);
	printf ("\n");
	}
	printf("\t");
	for (ix=0;ix<=numX;ix++)
	   printf("-");
	printf ("\n");
    printf("%6.2f", Xmin);
    for (ix=0;ix<=numX-12;ix++)
	   printf(" ");
    printf("%6f\n", Xmax);
	      
	
    
}

void array_minmax(double X[], int N,  double *min,  double *max)
{
    int i;
    *min=+1e100;
    *max=-1e100;
    for (i=0;i<N;i++)
	{
	    if (X[i]>*max) *max=X[i];
	    if (X[i]<*min) *min=X[i];
	}
    return;
}

void binHisto(long histo[], long N, long binSize)

/* expects a histogram array (1-off) in which element histo[k] is the count for a variate value of k.
   Puts sum of counts into bins of size binSize. N is the number of elements (array length -1)  */

{
long *binHisto, i,j,base,numBins;
numBins = N/binSize+1;
binHisto=(long *)myMalloc((numBins+1)*sizeof(long));
base=1;
for (i=1;i<=numBins;i++)
	{
	base=1+(i-1)*binSize;
	binHisto[i]=0;
	for (j=base;j<base+binSize;j++)
		binHisto[i]+=histo[j];
	}
printf("Binned histogram\n");
for (i=1;i<=numBins;i++)
	if (binHisto[i]>0)
		printf("%li\t%li\n",i*binSize,binHisto[i]);

}
                                                                                                                                                                                                                                                                                                                                                       r8s/MyUtilities.h                                                                                   0000644 0000766 0000120 00000001123 07565204055 013702  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdlib.h>

void			 doGenericAlert(char* errorMsg);
void 			fatal(char *s);
void			strtoupper(char *s);
void			concat(char **destHndl, char *source);
char*			DupStr(char* s);
char*			makeEmptyString(void);
FILE* 			PromptFileName(char* promptMsg, char* mode);
int				isStrInteger(char* s);
void array_minmax(double X[], int N,  double *min,  double *max);
void dumbPlot(double X[],  double Y[], int N);
char * slurpFile (FILE * inFileStream, long maxSize);
void binHisto(long * histo, long N, long binSize);
                                                                                                                                                                                                                                                                                                                                                                                                                                             r8s/NRCvectorUtils.c                                                                                0000644 0000766 0000120 00000003460 10113434106 014271  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "NRCvectorUtils.h"
#include "stdio.h"
#include "stdlib.h"
#include "memory.h"
#include "math.h"

#define SQR(a) ((a)*(a))
#define PNORM 2

double norm(double v[],int nl, int nh)

/* finds the p-norm of a 1-off vector using components from nl..nh */

{
int k;
double z=0.0;
for (k=nl;k<=nh;k++)
	{
	if (PNORM==2)
#if PNORM==2
		z+=SQR(v[k]);
#else
		z+=pow(v[k],PNORM);
#endif
	}
return pow(z,1.0/PNORM);

}
double norm_not_active(double v[],int active[],int nl, int nh)

/* finds the p-norm of a 1-off vector using components from nl..nh. 
   Only include elements not in the active set: these have active[] == 0*/

{
int k;
double z=0.0;
for (k=nl;k<=nh;k++)
	{
	if (PNORM==2)
#if PNORM==2
	    if (active[k] ==0)
		z+=SQR(v[k]);
#else
		z+=pow(v[k],PNORM);
#endif
	}
return pow(z,1.0/PNORM);

}

double *vector(int nl, int nh)
{
	double *v;
	v=(double *)myMalloc((unsigned) (nh-nl+1)*sizeof(double));
	if (!v) nrerror("allocation failure in vector()");
	return (v-nl);
}
void free_vector(double *v, int nl, int nh)
{
	myFree((char*)(v+nl));

	return;

}
void nrerror(char error_text[])
{
	void exit();long j;
	printf("Numerical recipes run-time error...\n");
	printf("%s\n",error_text);
	printf("...now exiting to system...\n");
     exit(1);
}
double **matrix(int nrl, int nrh, int ncl, int nch)
{
	int i;
	double **m;
	m=(double **)myMalloc((unsigned)(nrh-nrl+1)*sizeof(double*));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m -= nrl;
	for (i=nrl;i<=nrh;i++)  {
		m[i]=(double *) myMalloc((unsigned)(nch-ncl+1)*sizeof(double));
		if (!m[i]) nrerror("allocation failure 2 in matrix()");
		m[i] -= ncl;
	}
	return (m);
}

void free_matrix(double **m, int nrl, int nrh, int ncl, int nch)
/* free a float matrix allocated by matrix() */
{
int i;
for (i=nrl;i<=nrh;i++)
    {
    myFree(m[i]+ncl);
    }
myFree(m+nrl);
}
                                                                                                                                                                                                                r8s/NRCvectorUtils.h                                                                                0000644 0000766 0000120 00000000512 07565204055 014310  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  double norm(double v[],int nl, int nh);
double norm_not_active(double v[],int ac[],int nl, int nh);

double *vector(int nl, int nh);
void free_vector(double *v, int nl, int nh);
void nrerror(char error_text[]);
double **matrix(int nrl, int nrh, int ncl, int nch);
void free_matrix(double **m, int nrl, int nrh, int ncl, int nch);
                                                                                                                                                                                      r8s/ObjFunc.c                                                                                       0000644 0000766 0000120 00000122273 10637044177 012755  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "NRCvectorUtils.h"
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
h=0.00001;
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
//	printf("***[%i] %12.10f %12.10f %12.7f %12.7f %f\n",i,p1,p2,f1,f2,dif);
	}
return;
}
                                                                                                                                                                                                                                                                                                                                     r8s/ObjFunc.h                                                                                       0000644 0000766 0000120 00000002466 10155467117 012761  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #ifndef _LANG_FITCH
#define _LANG_FITCH
#include "TreeUtils.h"

typedef double (*objfunc)(double[]); 

#define USER -1
#define LaF 0
#define HMM 1
#define NP  2
#define GAMMA 3
#define SLDWIN 4
#define PENLIKE 5
#define PENLIKET 6
#define LFLOCAL 7

#define POWELL 	0
#define QNEWT 	1
#define TN	2

NODETYPE** newAllNodeArray(TREE t);
NODETYPE** newNodeArray(TREE t);
void Dapprox(double p[],double grad[],int n, double (* obj)(double p[]),double h);
int set_active(double * solution, NODETYPE **nodeArray,int nNodes,int active[],double active_eps);
int checkGradient(TREE t,double solution[],double gradient[],double Obj,double ftol,int verbose);
void pTimeArray2tree(NODETYPE *node,double pTime[]);
static void tree2pTimeArray(NODETYPE *node,double pTime[]);

static double * allocateTimeArray(NODETYPE * root, int method,int nRates,int *tvar, int *nvar,double **D);

void peak_peek(objfunc objective,double p[],int nvar,double sizeFactor, int grid);
int perturb(
	double p[], 
	int nvar,
	int numperts, 
	double perturb_factor,
	double local_factor, 
	double unpertOpt,
	objfunc objective
	);
int perturb_p(double p[], int n, double perturb_factor);
int same_points(double p1[], double p2[], int n, double tolerance);
void doObjFunc(TREE t, int method, int nRates,int algorithm, int * success);

void ObjFuncDefaults(void);
#endif
                                                                                                                                                                                                          r8s/ObjFuncHeader.h                                                                                 0000644 0000766 0000120 00000007427 10157130415 014062  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /**** This is just the header output...got tired of looking at it ****/


if (verbose > 0)
  {
  printf("\n\n----------------------------------------\n\n\n");
  printf("LINEAGE RATE/TIME ANALYSIS FOR TREE %s\n\n",TreeName);
  switch (method)
	{
	case PENLIKET:printf("Method = Penalized Likelihood(T)\n");printf("Smoothing factor = %f\n",rbp->smoothing);break;
	case PENLIKE:printf("Method = Penalized Likelihood\n");printf("Smoothing factor = %f\n",rbp->smoothing);
		switch (rbp->NeighborPenalty)
			{
			case 0:printf("Penalty function = Ancestor-Descendant\n");break;
			case 1:printf("Penalty function = Neighbor variance\n");break;
			}
		switch (rbp->PenaltyType)
			{
			case 0:printf("Scale for rate penalty = ADDITIVE\n");break;
			case 1:printf("Scale for rate penalty = LOG\n");break;
			}
		printf("Minimum allowed rate = %f of initial average rate estimate\n",rbp->minRateFactor);
		printf("Minimum allowed duration on 0-length terminal branches = %f of root's age\n",rbp->minDurFactor);
		break;
	case LaF:printf("Method = Langley and Fitch\n");break;
	case LFLOCAL:printf("Method = Langley and Fitch (with %i local rates)\n",nRates);break;
	case GAMMA:printf("Method = Gamma-Negative-Binomial\n");break;
	case HMM:printf("Method = Hidden Markov\n");break;
	case NP:
		{
		printf("Method = Non-parametric (exp=%f)\n",rbp->npexp);
		switch (rbp->PenaltyType)
			{
			case 0:printf("Scale for rate penalty = ADDITIVE\n");break;
			case 1:printf("Scale for rate penalty = LOG\n");break;
			}
		if (gVarMinFlag)
		    printf("(Minimizing variance of local rates)\n");
		else
		    printf("(Minimizing local transformations (NPRS))\n");
		}
	}
  switch (algorithm)
	{
	case POWELL: printf("Optimization via Powell's method\n");break;
	case QNEWT:  printf("Optimization via quasi-Newton method with analytical gradients\n");break;
	case TN:     printf("Optimization via Truncated-Newton (TN) method with bound constraints\n");break;
	}

  printf("\n----------------------------------------\n");

    printf("Substitution Model\n");
    if (rbp->RatesAreGamma)
    	{
	printf("\tRates are gamma distributed across sites\n");
	printf("\tShape parameter of gamma distribution = %6.2f\n",rbp->alpha);
	}
    else
	printf("\tRates are equal across sites\n");
	

    printf("Global/Local Search Parameters\n");
    printf("\tNumber of searches from random starts = %i\n",NUM_TIME_GUESSES);
    printf("\tNumber of restarts after each search = %i\n",rbp->num_restarts);
    printf("\tLocal perturbation on restarts = %4.3g\n", rbp->perturb_factor);
    printf("\tLocal fractional tolerance after restarts = %4.3g\n", rbp->local_factor);
	

	if (algorithm==TN)
		printf("Optimization parameters set automatically by TN routine\n");
	else
		{
	  	printf("Optimization parameters\n");
		printf("\tMaximum number of iterations  = %i\n", rbp->maxIter);
		printf("\tFunction Tolerance = %4.3g\n", rbp->ftol);
		printf("\tlinminOffset = %4.3g\n",rbp->linminOffset);
		printf("\tContract Factor = %4.3g\n", rbp->contractFactor);
		printf("\tMax number of contract iterations  = %i\n",rbp->maxContractIter);
    	}
  if (method==NP)
	{
  	if (gClampRoot)
		printf("Root rates are CLAMPED\n");
  	else
        	printf("Root rates are FREE\n");
	}
  if (gisConstrained && algorithm != TN)
    {
    printf("Time constraints are enforced using barrier optimization\n");
    printf("\tBarrier Tol = %4.3g\n", rbp->barrierTol);
    printf("\tMax Barrier Iter = %6i\n", rbp->maxBarrierIter);
    printf("\tInit Barrier Factor = %4.3g\n", rbp->initBarrierFactor);
    printf("\tBarrier Multiplier = %4.3g\n", rbp->barrierMultiplier);
    }
  printf("Length of tree input =  %g\n",treeLength(root));
  printf("Number of taxa  =  %i\n",numdesc(root));
  printf("Number of sites in sequences =  %li\n",rbp->numSites);
  printf("\n----------------------------------------\n");
  }
                                                                                                                                                                                                                                         r8s/ReadNexusFile2.c                                                                                0000644 0000766 0000024 00000561644 11644142002 014217  0                                                                                                    ustar   sandermj                        staff                                                                                                                                                                                                                  /* REVISION HISTORY
 * 
 * 8.14.99.   Radical.  Changed everything to case insensitive by redefining isEqual macro
 *	    in structures.h, and changing the strtoupper() statements in nexttoken2.c
 *	    This is all more or less togled by the STU macro defn in structures.h
 * 9.3.99  Radical. Changed parse-assignment function to use a static
	    char buffer for LocalToken. New function is parse_assignment2

   9.13.99 Fixed 'compar()' function used in qsort routines TWICE in this program
    It was bogus and worked on SGI (coincidentally) but not on LINUX

   10.13.99 Changed 'BranchLike' in the LF algorithm so it converts real branch lengths
	to integer branh lengths before calculating likelihoods...formerly this was allowing
	real values for branch lengths...although this should never have come up
    1.16.00 Fixed nasty bug in setupFeasibleTimes. Had failed to initialize
	minTime=0.0 prior to descMinAge call.  Only a problem on some hardware.
    1.16.00 Fixed bug in ABCSuperTree routines that assigned a large double (1e100) as a
	branch length to the wtSet[] array, which is a Float, causing an exception on some
	hardware.  In the newNode function, I now use FLT_MAX instead of 1e100 (KLUDGE!) as
	the nulll branch length, and in TreeDrawUtils, which is the only place I need the
	stupid thing, I recognize FLT_MAX too.
    4.19.00 Changed the multiply_branchlength_by command so that final branch lengths are rounded
	to the nearest integer.  This allows me to run the collapse command afterwards, and take
	care of those nasty zero-length branches. (otherwise, I could collapse, but this would be prior
	to rounding down to the integer of zero, and then I was stuck with zero-length branches again
    4.21.00  Added a feature allowing 'setage taxon=root age=xx' to permit fixing the root node age of the tree.
	This corrects a bug that allowed the default root age to be greater than 1.0 whenever internal constraints
	were set with the constrain_time commands.  Now you have to explicitly set the root's age, or watch out!
    9.19.00  Changed command syntax for several commands so that options are remembered between runs (by
	changing to static values of parameters in functions. Incl.: divtime, describe.

   11.1.00 Series of significant modifications to code. Working on release 1.0.

    5.26.01 Changed PL and  NP routines to use the variance of rates across the root's children's branches
		as the quantity to be minimized. This may have slight effects on previous runs. Updated gradients
		for PL and objfuncs for PL and NP to do this.
    5.26.01 Added gamma distributed rates functionality for PL (w/POWELL only)--gradients not done yet...Note there
		is an issue in that BranchLike uses units of raw substitutions, whereas BranchLikeNegBinom uses subst/site.
    7.9.01 Added estimates of confidence intervals on a single node time using curvature of
	likelihood surface (see doDivTime:confidence)
	12.11.01 Fixed factorial calculations to allow any argument
	12.18.01 Fixed nasty bug in gradient calculations. Had not initialized g=0, which
			was a problem sometimes in estimating root node (long story)
	2.8.02 Fixed collision between CALIBRATE command and FIXAGE command. Used to permit fixage on an
		ultrametric tree, but calibrate command expects that ages were all scaled on 0 to 1 and not
		messed with after that initialization. Added warnings.
	2.16.02 Dumb bug: NPRS wanted to do a gradient check in MinND, but I don't have the gradient! Added test
		(should add gradient and make termination criteria a clean separate routine
	2.17.02 Added ROUND=YES|NO option to BRLEN command. This permits user to NOT round branch lengths
		to nearest integer. For input ultrametric trees, we'd really like to keep exact real lengths
		for use in the CALIBRATE command
	3.28.02 (Matt Lavin) Gradient check was too strict when using constrained nodes (should ignore non-zero
		gradient in the direction of the constraint if the constraint is active). Deciding that the 
		constraint is active relies on ACTIVE_EPS tolerance factor, which I've increased to 0.01
	4.15.02 (Ben Warren) Roundoff error bugs:
			1. In Langley-Fitch, use of very large ages (e.g., 500000) caused inconsistent results.
			The ZERO() macro in the BranchLike function was stupidly rounding arguments to 0.0 when
			they were merely small. I've just removed this macro.
			2. In NPRS, the same inputs were causing problems. Here the culprit was adding +1.0 to the
			objective function in ObjNP(). Originally I'd added this because Powell's termination criteria
			get fooled whenever there is truly a clock for NPRS (then the objective function = 0). Temp
			solution is to drop the addition of 1.0 but this will probably re-introduce the other bug. Need
			to improve Powell or use other optimization engines. 
	6.4.02 (Torsten Eriksson) Fixed very stupid memory leak in Powell. I had left some debugging code in there which
			allocated but did not free meory in linmin
		Also disabled the convergence diagnostics for Powell; probably slow down program a lot; can easily uncomment
		the changes in 'powell1'
	6.10.02 Several changes to BranchLike and LogFactLookup. Bug fixes concerning 0-length branches. See code for details 

Some insights on zero-length branches under PL: The likelihood component for a branch with k=0 substitutions is exp(-rT), 
which is a maximum along the directions of both r=0 and T=0. This could easily cause optimization methods to have problems,
although it seems that gradient methods perform worst. Once r =~ 0 or T =~ 0, the gradient will be 0 for the other variable
and the function will be maximized. For high smoothing values, the neighboring branches with k>0 substitutions keep these problems
in check, but as smoothing gets low, there comes a point at which r can be 0 without too  much of a penalty...then we have 
problems optimizing. The solution apparently is to set a small but finite minimum lower bound on rates.
	6.11.02. Rewrote and cleaned up tests in BranchLike. Better check to make sure this doesn't mess everything up downstream 
	6.22.02. Modified the check_feasible routine! Now a point is considered feasible even if a node time is EQUAL to a min or max
		 age constraint. Prior to this, I required strictly greater or less than. This is necessary because under TN algorithm
		 we frequently find points exactly satisfying constraint. Then when 'perturb_p' tries to work, it checks every node in 
		 tree for strict inequality, failed and bombed. [can also solve just by not restarting...]
	6.23.02. 'Bug' fixed in TNBC line search. For 0-length branches, the linesearch kept iterating. Now we force termination in the
			hopes that a restart perturbation will fix things up!
	6.24.02. Bug fixed in cross-validation. Zero-length branches were sometimes causing spurious calculations of chi-squared values.
		Needed to check for zero-length branches in cvSquareErrorBranch routine and prevent division by zero.																			
	7.25.02  See fix on 6.22.02. Didn't quite fix this right. Made sure all tests were for strict inequality
	7.26.02  Removed "include <malloc.h>" in ReadNex... and memory.c, the only two places it is invoked. Couldn't find this include
		 file under Mac OS X.
	7.26.02  Fixed bug in pruneTaxon command, which did not correctly update all pointers in tree structure when a taxon was deleted,
		 causing erroneous print outs under "describe plot=xxx_description". Other option values in describe worked OK.
	11.6.02  Eliminated all drand48 and srand48 references, and replaced with myRand() function which uses stdlib rand(). Also call srand() rather than srand48 everywhere. Hopefully this will make cross-platform development easier.
	11.6.02  (Torsten Eriksson) Eliminated, hopefully, a memory leak in TN (TNwrapper), where I failed to deallocate arrays.
	11.?.? [sometime while in Germany] Added a log transformation to the penalty function in NP and PL. Invoked by SET PENALTY=ADD|LOG;
		Note that this only works with POWELL-haven't calculated derivatives yet.
		Also began to implement a general 'neighbor variance' penalty, not done, see TimeAlgorithms...
	4.3.03  Upon migrating to MAC OS 10.2 found a bug in the COLLAPSE command which is only fixed by REMOVING the node 
		destructor (was screwing up the recursion). Now I've got a bunch of dangling nodes still! Figure out how to
		deallocate them (perhaps maintain a global junk list)
	4.31.03 Added feature. A likelihood ratio based relative rate test at a user supplied node in the tree: RRLIKE TAXON=XXX;
	5.13.03 Modified warnEstRoot so that it does not issue warnings itself but merely returns error code.
		Continued work on RRLIKE command. It now always takes time constraints or fixed ages into account if they are available for the clades
		in question. If no constraints or fixed ages are available, then it sets the root of the focal subtree to an arbitrary
		value.
	6.3.03 Fixed the 'birthDist' function to permit large values of lambda. Currently these overflowed. Can replace the
		density with a simpler form when lambda is large
	6.8.03 In TreeSim, I now require a new seed for each run; otherwise it uses a default seed and issues a warning. I didn't
		like having a seed held over from previous runs.
	7.16.03 Lots of stuff:
			working on a fossil cross validation scheme
			added log penalty scaling with gradient this time to penalized likelihood
			added neighbor variance penalty to PL in conjunction with log scaling
	9.19.03 Added the VCV function which calculates the variance covariance matrix of a tree based on the lengths
			of the branches subtending the MRCA of each pair of taxa.
		Syntax is VCV taxon = name; which works on the subtree descended from name
	3.3.04  Added another fossil cross validation scheme, fossilcrossvfixed, which uses only a set of fixed ages
	6.14.04 Made two submodels of BDback, depending on whether we normalize the root age to 1. 
			diversemodel=bdback command does not; 'bdbacknormal' does.
	8.26.04 Fixed bug in profile branch command which did not correct for trees in which branch was missing (A. Antonelli)
	8.26.04 Fixed bug in doObjFunc which mistakenly set algorithm to 'TN' in multiple time guesses by an expression 'if (algorithm=TN)' blunder (changed to '==' : also A. Antonelli).
	8.26.04 Important changes to PL routines: now estimating the initial starting guess on the rate by doing a LaF clock
		analysis first and taking that estimator as the guess. This required putting a wrapper around the doObjFunc
		routine. Other methods unaffected, still using crude guess.
	8.26.04 Begun some error reporting. Will now signal error if the basic command is wrong, but still overlooks merely wrong option setting
		syntax
	...?	Implemented CO command (continuous character rate estimation)
	12.4.04 Beginning work on checkGradient routine. Problems identified in correct setting for ACTIVE_EPS. Warning now
		issued if a min and max constraint might BOTH be treated as generating an active constraint (and thus the program
		arbitrarily picks one, leading to a false conclusion that the gradient's sign is wrong at the constraint)
	12.6.04 Problem noted in optimization when user specifies ROUND=NO. The calculation of the objective function rounds
		the character lengths on a branch always, but the calculation of the gradient (at least in LF), uses real 
		arithmetic. This causes routines to come to different conclusions depending on whether it uses gradient or 
		non gradient methods. Effect is slight on param estimates or obj func estimate, but in case examined when 
		checking the gradient using numerical approx it is quite important. Added warning message in the BLFORMAT 
		command, and now DIVTIME will bail with a warning if you try to use any method other than NPRS without
		rounding input (NPRS doesn't rely on calculation of a likelihood, so not an issue
	12.8.04 Fixed bug: ObjFunc.c was subtracting 1.0 from obj() under NPRS because at one point I had added +1.0 to it
		in TimeAlgorithms.c (to fix what I thought was a problem with clocklike data sets returning an objective
		function value of 0.0). However, I had stopped adding the 1.0 in the obj() in TimeAlgorithms.c.
		This was sometimes generating NPRS obj() values less than zero. Now the two modules consistently do not
		modify the objective.
	12.8.05 Changed 'execute' command so that we can read multiple blocks from different files, for example, to read
		a trees block from one file and a r8s block from another
	12.10.05 Implemented an analytical gradient for log penalty function under PL.
	1.8.05 Added charset ouput to MRP command.
	4.29.05. Added warning about only using POWELL with LFLOCAL in response to 'bug report'
	12.05. Added ancestral state reconstruction using squared change parsimony
	1.06.06. Removed all traces of 4PL and jk4PL code. Mix of C and C++ causing some problems on some compilers
	1.06.06. Using -pedantic gcc option to find misc bugs. Removed HUGE_VAL from two functions (this is a reserved word)
	7.25.06. False bug in user supplied ultrametric tree stuff...make sure user sets ROUND=NO and lengths=persite if latter
			is appropriate. Added routing rootToTips to print out those distances. Found slight roundoff error in PAUP
			tree descriptions in this respect.
Find Hilmar Lapp's email, where he found some bugs in interactive mode that need fixing....
	6/22/07. The rounding issue in branch lengths came up again. To reiterate from above, *all* branch lengths input to doObjFunc...must
			be integers. This is because the likelihood function converts them to integers but the gradient functions do not, leading
			to frequent catastrophic convergence failures...Major change to doObjFunc is to call 'traverseMultiplyLength', which is used
			to force rounding of all branch lengths before any optimization is done. Now there is no choice. Previously, I warned the user
			about the issue but didn't actually fix it.	
    9/2011.  Adding a new feature on a 3-state markov character model (covarion.c).
    		 Added a feature for options on rerooting, and fixed some rerooting bugs. Now we can
    		 reroot at a NODE, instead of just former behavior where we always had a binary root.
    		 */

/****  Module for Nexus File functions  *******/

#define EFRON1996	1  /* only needed for this flavor of bootstrapping */


/*******************************************************************************/

#include "continuousML.h"
#include "NRCvectorUtils.h"
#include "storeNexusFile.h"
#include <sys/types.h>
/* #include <malloc.h> */
#include "Maximize.h"
#include "WuLi.h"
#include "nexus.h"
#include "MyUtilities.h"
#include "memory.h"
#include "ObjFunc.h"
#include "TimeAlgorithms.h"
#include "myOutput.h"
#include "TreeSim.h"
#include "ObjFunc.h"
#include "DrawTree.h"
#include "moment.h"
#include "DistrFuncs.h"
#include "distance.h"
#include "structures.h"
#include "TreeUtils.h"
#include "ancestral.h"
#include "covarion.h"
#include "relativeRates.h"
#include <math.h>
#include <ctype.h>


/*****  private functions ******/

int parse_assignment2(char * target);
static void doVCVCommand(void);

static void doCovarionCommand(void);
static void doAncestralCommand(void);
static void doContOptCommand(void);
static void histoStat(long h[], long N, long nTaxa,long *count, double *mean, double *freq1class, long *maxS, double *dominance);
static void doRRLikeTestCommand(void);
static void doConfidence(TREE thisTree,char * nodeName,int method,int nRates,int algorithm,double cutoff,int JMAX);
static void doLocalModelCommand(void);
static void doBLFormatCommand(void);
static void doFossilCrossVfixed(TREE tree, int method,int nRates,int algorithm,double cvStart,double cvInc,double cvNum, int cvSample);
static void doFossilCrossV(TREE tree, int method,int nRates,int algorithm,double cvStart,double cvInc,double cvNum, int cvSample);
static float doCrossV(TREE tree, int method, int nRates,int algorithm, double c1, double c2, double c3, int);
static void doUnSetAgeCommand(void);
static void doShowAgeCommand(void);
static void doReRootCommand(void);
static void doPruneTaxonCommand(void);
static void doSetAgeCommand(void);
static void doClusterHistogramCommand(void);
static void doCollapseCommand(void);
static void doClearTreesCommand(void);
static void doExecuteCommand(void);
static void printHelp(void);
static void doDivTimeCommand(void);
static void doSimpleCladeCheckCommand(void);
static void doEFRON_Weights_Command(char *buffer);
static void doB_Weights_Command(char * buffer, char *buffer1, char* buffer2);
static void efron1996(int *weightArray,int nchars,int num_points,char *buffer, 
    long *index);
static void doCladeCheckCommand(void);
static void doBootCharCommand(char* buffer);
static void doBranchProfileCommand(void);
static void doClade_Set_Command(void);
static void printNexus(int ntaxa, int nchars, StrListPtr taxaList,  char **matrix);
	
static void doSuperCommand(void);  
static void doCalibrateCommand(void);
static void doPrintCommand(void);
static void doSimBlock(void);
static void doSimCommand(void);
static void doTaxaSetCommand(void);
static void doSetCommand(void);
static void doDistanceCommand(void);
static void doBSCommand(void);
static void doError(char* p[], int which);
static void doDataBlock(void);
static void doBootBlock(void);
static void doTranslateCommand(void);
static void doWuLiCommand(NODETYPE * root);
static void doCharsBlock(void);
static void doFormat(void);
static void doMatrix(void);
static void doIndel(void);
static void doTaxaBlock(void);
static void doCharDimensions(void);
static void doTaxDimensions(void);
static void doTaxLabels(void);
static void doTreeBlock(void);
static void doTreeCommand(void);
static void doExSets(void);
static void checkMatrix(void);
static void doUnrecognizedBlock(void);
static void doSitesCommand(int); /* exclude third positions */
static void doConstrain_TimeCommand();
static void doBootCommand(StrListPtr b, char* bu);
static char * doIncludeCommand(void);
static StrListPtr doFixedTaxaListCommand(void);
static void doMRCACommand(void);
static void doLengthMultiplyCommand();

static void doSaveTree(NODETYPE *root);
static int parse_assignment(char * target,char ** token);
static void doBD(void);
void doDimensions(void);
void doMatrixGeneral(void);

/******   globals   ********/

StackPtr gFStack,gPStack;
char LocalToken[MAX_LOCAL_TOKEN_SIZE];

int gEstRoot;
int gInteractive,gLabel;
int gSeedisSet=0;

StrListPtr	gTaxaList;
int	gNewLine;
int	gColumn;
int 	gFirstDesc;
char 	*bufPtr;
double 	gnpexp;
StrListPtr gTaxaSet;	
char  	*aTokenPtr;
int 	curTree=0;		/* index for the current tree description being parsed */
	
char	*nexError[2]= {
						"Error: Not a NEXUS file",		/* 0 */
						"Error opening NEXUS file"		/* 1 */
						};

struct NexDataType gNexData;	/* This is THE data structure for the NEXUS data */
struct NexDataType *gNexDataPtr;	/* This is THE data structure for the NEXUS data */

/**************************/
/**************************/

/* 
	Read a NEXUS file buffer and set up a global data structure containing everything. 
	See nexus.h for that data structure.
	Returns NULL on error. 
*/

void readNexusFile(char * theNexusFileBuffer)
{
	char *stemp;
	int c;
	int flag;
	int ix;
	long bufLength;

	struct NexDataType *nptr;
	
/*mallopt (M_DEBUG, 1);*/

	bufPtr=theNexusFileBuffer;	/* Initialize this global to beginning of the 
					buffer and will sweep through it until end of buffer */
	
	bufLength=strlen(theNexusFileBuffer);	

	if ( bufPtr != NULL )
		{
		aTokenPtr=nextToken();
		if(!isEqual(aTokenPtr,"#NEXUS"))
			{
			doError(nexError,0);
			return;		/* not a NEXUS file */
			}
		while (aTokenPtr=nextToken(), *aTokenPtr)
			{
			if (isEqual(aTokenPtr,"BEGIN"))
				{
				nextToken();if (!*aTokenPtr) return;
				stemp=DupStr(aTokenPtr);	/* get the block name and store in 'stemp'*/
				if (!stemp) 
					fatal ("Error reading block name");
				
				if (isEqual(aTokenPtr=nextToken(),";")) /* pop the terminating semicolon */
					{
					if (isEqual(stemp,"TAXA"))
							doTaxaBlock();
					else				
					if (isEqual(stemp,"CHARACTERS"))
							doCharsBlock();
					else				
					if (isEqual(stemp,"TREES"))
							doTreeBlock();
					else				
					if (isEqual(stemp,"RATES"))
							doRateBlock();
					else				
					if (isEqual(stemp,"R8S"))
							doRateBlock();
					else				
					if (isEqual(stemp,"BOOTSTRAP"))
							doBootBlock();
					else				
					if (isEqual(stemp,"SIMULATION"))
							doSimBlock();
					else				
					if (isEqual(stemp,"DATA"))
							doDataBlock();
					else 
						{  /* token is not a recognized block */
							doUnrecognizedBlock();
						}
					/*if (!*aTokenPtr)break;*/
					}
				free(stemp);				
				}
			}
		return;	/* normal return */
		}
	else
		{
		doError(nexError,1);
		return;
		}
	
}
/****************************************************************/
void doInteractive(void)
{
#define BUFSIZE 500
char Token[MAX_TOKEN_SIZE];
char inputBuffer[BUFSIZE],*pLast;
aTokenPtr=Token;
    for (;;)
	{
	printf("\nr8s>");
	fgets(inputBuffer,BUFSIZE,stdin);
#if STU
	strtoupper(inputBuffer);
#endif
	if (strlen(inputBuffer)>0)  /* if something in the buffer */
		{
		pLast=inputBuffer+strlen(inputBuffer)-1;
		if (*pLast != ';')
			{
			*(pLast+1)=';';
			*(pLast+2)='\0';
			}
		bufPtr=inputBuffer;
			doRateBlock();
		}
	    
	}
return;
}
void doCommandLineControl(char *inputBuffer)
{
char Token[MAX_TOKEN_SIZE];
char *pLast;
aTokenPtr=Token;
#if STU
strtoupper(inputBuffer);
#endif
if (strlen(inputBuffer)>0)  /* if something in the buffer */
		{
		pLast=inputBuffer+strlen(inputBuffer)-1;
		if (*pLast != ';')
			{
			*(pLast+1)=';';
			*(pLast+2)='\0';
			}
		bufPtr=inputBuffer;
		doRateBlock();
		}
  
    
}
static void doExecuteCommand(void)
{

char *theNexusFileBuffer, fnInput[FILENAME_MAX];
FILE * inStream =NULL;
strcpy(fnInput, aTokenPtr=nextToken());   /* set file name */
#if 0
if (gNexDataPtr)
	freeNexusStructure(gNexDataPtr);
gNexDataPtr=initialize_nexus();
#endif

// *** try this to allow multiple reads from different files... 
if (!gNexDataPtr)
	gNexDataPtr=initialize_nexus();
// ***
if (!gNexDataPtr)
	fatal("Failure to allocate nexus data structure in main.c");
if (!(inStream=fopen(fnInput,"r")) )
		    {
		    printf("file=%s\n", fnInput);
		    fatal("Error in file handling\n");
		    }
if (inStream)
	    {
	    theNexusFileBuffer=storeNexusFile(inStream);
	    readNexusFile(theNexusFileBuffer);
	    };
return;
}


static void printHelp(void)
{

	char * filename="r8s.helpfile";
	FILE* fpntr;
	char * buffer;
	buffer=NULL;
	if (  (fpntr=fopen(filename,"r")) )
		{
		buffer=slurpFile (fpntr, 10000);
		printf("%s\n", buffer);
		}
	else
	    printf("Failed to open file '%s'\n", filename);
	
	
	return;

}
	
/****************************************************************/
struct NexDataType * initialize_nexus(void)
{

struct NexDataType *p;

p=(struct NexDataType *)myMalloc(sizeof(struct NexDataType ));
if (!p)
    return NULL;
aTokenPtr=(char *)myMalloc(MAX_TOKEN_SIZE*sizeof(char));
if (!aTokenPtr)
    fatal ("Couldn't allocate aTokenPtr");
    
gTaxaSet=NULL;
gLabel=1;

p->TDList = newStrList(); /* initialize the list of tree descriptions */
p->TDLabelList = newStrList(); /* initialize the list of tree labels */
p->TaxaList = newStrList();
p->TransList=newStrList();
p->inTrees=NULL;
p->TaxSetNameList=NULL;
p->TaxSetLists=NULL;
p->excArray=NULL; /* can't be initialized further until we know num of chars */

p->isChars=0;
p->isTrees=0;
p->isTaxa=0;
p->isTranslate=0;	/*...flags for when these elements are read */
p->NTaxa=0;
p->NChars=0;
p->NumTrees=0;
p->matchchar='.';
p->gapchar='-';
p->missingchar='?';

p->RateBlockParms.NeighborPenalty=0;
p->RateBlockParms.PenaltyType=0;
p->RateBlockParms.checkGradient=0;
p->RateBlockParms.clampRoot=1;
p->RateBlockParms.isBS=0;
p->RateBlockParms.NReps=1;
p->RateBlockParms.seed=1;
p->RateBlockParms.RRtype=WULI;
p->RateBlockParms.npexp=2.0;
p->RateBlockParms.verbose=1;
p->RateBlockParms.num_restarts=1;
p->RateBlockParms.num_time_guesses=1;
p->RateBlockParms.num_rate_guesses=1;
p->RateBlockParms.smoothing=1.0;
p->RateBlockParms.showGradient=0;
p->RateBlockParms.showConvergence=0;
p->RateBlockParms.ftol=1e-6;
p->RateBlockParms.barrierTol=0.0001;
p->RateBlockParms.activeEpsilon=0.001;
p->RateBlockParms.maxIter=500;
p->RateBlockParms.maxBarrierIter=10;
p->RateBlockParms.initBarrierFactor=.25;
p->RateBlockParms.barrierMultiplier=0.10;
p->RateBlockParms.linminOffset=0.05;
p->RateBlockParms.contractFactor=0.1;
p->RateBlockParms.maxContractIter=10;
p->RateBlockParms.local_factor=0.01;
p->RateBlockParms.perturb_factor=0.01;
p->RateBlockParms.RatesAreGamma=0;
p->RateBlockParms.alpha=1.0;
p->RateBlockParms.numSites=1;
p->RateBlockParms.lengthFmt=0;  //'total'
p->RateBlockParms.roundFlag=1;  // round input branch lengths by default
p->RateBlockParms.clockFmt=0;
p->RateBlockParms.minRateFactor=0.05;
p->RateBlockParms.minDurFactor=0.001;

return p;
}
/****************************************************************/
/****************   BLOCK PROCESSING FUNCTIONS ******************/
/****************************************************************/


void doDataBlock(void)
{
	do 				/* need to put in error checking in case no DIMENSIONS statement */
		{
		aTokenPtr=nextToken();
		if (isEqual(aTokenPtr,"DIMENSIONS"))  
			doDimensions();
		if (isEqual(aTokenPtr,"FORMAT"))  
			doFormat();
		if (isEqual(aTokenPtr,"MATRIX"))
			doMatrixGeneral();

		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));
	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
return;
}
/****************************************************************/
void doUnrecognizedBlock(void)
{
	do 				
		{
		aTokenPtr=nextToken();
		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));
	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
return;
}
/****************************************************************/
void doCharsBlock(void)
{
	do 				/* need to put in error checking in case no DIMENSIONS statement */
		{
		aTokenPtr=nextToken();
		if (isEqual(aTokenPtr,"DIMENSIONS"))  
			doCharDimensions();
		if (isEqual(aTokenPtr,"FORMAT"))  
			doFormat();
		if (isEqual(aTokenPtr,"MATRIX"))
			doMatrix();

		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));
	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
return;
}
/****************************************************************/

static void doBootBlock(void)
{
StrListPtr fixedList=NULL;
char *buffer, *buffer1, *buffer2;
buffer=NULL;
buffer1=NULL;
buffer2=NULL;
	do 				/* need to put in error checking in case no DIMENSIONS statement */
		{
		aTokenPtr=nextToken();
		if (isEqual(aTokenPtr,"CLADE_SETS"))  
			doClade_Set_Command();
		if (isEqual(aTokenPtr,"MRCA"))  
			doMRCACommand();
		if (isEqual(aTokenPtr,"FIXED"))  
			fixedList=doFixedTaxaListCommand();
		if (isEqual(aTokenPtr,"BOOT"))  
			doBootCommand(fixedList, buffer);
		if (isEqual(aTokenPtr,"BOOTCHARS"))  
			doBootCharCommand(buffer);
		if (isEqual(aTokenPtr,"INCLUDE"))  
			buffer=doIncludeCommand();
		if (isEqual(aTokenPtr,"INCLUDE1"))  
			buffer1=doIncludeCommand();
		if (isEqual(aTokenPtr,"INCLUDE2"))  
			buffer2=doIncludeCommand();
		if (isEqual(aTokenPtr,"WEIGHTS"))  
			doB_Weights_Command(buffer,buffer1,buffer2);
		if (isEqual(aTokenPtr,"EFRON_WEIGHTS"))  
			doEFRON_Weights_Command(buffer);
	     	if (isEqual(aTokenPtr,"CLADE_CHECK"))  
			doCladeCheckCommand();
	     	if (isEqual(aTokenPtr,"SIMPLE_CLADE_CHECK"))  
			doSimpleCladeCheckCommand();

		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));
	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
return;
}
/****************************************************************/
static void  doClade_Set_Command()
{
/* for each tree, creates a clade list and stores a pointer to that list
   in the tree structure */
       
	int matchCount;
	StrListPtr allTaxaList;
	double ovl, ovlmax, ovlmin;
	PtrList SetList;
	PtrList lnode, pnode, focal_cladeSet, cur_cladeSet, fnode, clnode;
	TREE thisTree, focal_tree, curTree;
	Set fclade, cclade;
	int ntaxa, ntrees;
	while (!isEqual(aTokenPtr=nextToken(),";")); /* just the one word command */
	if (gNexDataPtr->isTrees)
		{
		lnode=gNexDataPtr->inTrees;
		allTaxaList=newStrList();
		thisTree=lnode->item;
		TreeToTaxaList(thisTree->root,allTaxaList); /* important to do this ONCE only! */
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			thisTree->cladeSet = Tree2CladeSet(thisTree, allTaxaList);
			}
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			SetList=thisTree->cladeSet;
			printCladeSets(SetList);
			}
		}
	freeStrList(allTaxaList);

#if 1 
/*      do the fuzzy bootstrap */

	if (gNexDataPtr->isTrees)
		{
		lnode=gNexDataPtr->inTrees;
		ntrees=pLengthList(lnode);
		if (ntrees >= 2) /* must have a focal tree and one other to continue! */
		    {
		    focal_tree=lnode->item;
		    focal_cladeSet=focal_tree->cladeSet;
		    printf("\nFuzzy Clade Analysis\nMinimax 	 BP  \t\t   Clade\n-----------------------------------\n");
		    LISTLOOP(focal_cladeSet) /* loop through all the clades in the focal tree */
			{
			fclade=focal_cladeSet->item;
			lnode=(gNexDataPtr->inTrees)->next;
			ovlmin=+9999.999;
			matchCount=0;
			 LISTLOOP(lnode)  /* loop through the other trees ...*/
			    {
			    curTree=lnode->item;
			    cur_cladeSet=curTree->cladeSet;
			    ovlmax=-9999.999;
			    LISTLOOP(cur_cladeSet) /*...checking each of their clades */
				{
				cclade=cur_cladeSet->item;
				/* test_set(fclade, cclade); */
				ovl=set_overlap(fclade, cclade);
				/*printf("overlap... %f\n", ovl);*/
				if (ovl>ovlmax) ovlmax=ovl;
				if (ovl==1.0)++matchCount; /* there was a clade matching the focal clade */   
				}
			    /*printf("max overlap for this clade = %f\n", ovlmax);*/
			    if (ovlmax<ovlmin) ovlmin=ovlmax;
			    }
			printf("%f\t%f\t", ovlmin, matchCount/(float)(ntrees-1));   
			print_set(fclade);
			}
			
		    }
		}
#endif	
    
}
static void doCladeCheckCommand()
/** Reads a list of taxa and checks to see if that group is a clade on all trees.
Reports the proportion of trees in which this group is a clade.
ADDED.  Lots of junk to do Efron 1996 bootstrap.  A real hack.  We generate sets of NUMBER_IN_EFRON_SET trees. The
first in every set is the P(j),  

To avoid a bias in which we always take the tree-weights over the boundary in the direction toward pj,  we
filp a coin to choose between that value of w and the next lower one,  which has pj on the other side of the
boundary toward pcent. **/
{

#define NUMBER_IN_EFRON_SET 26

	int isClade[NUMBER_IN_EFRON_SET];
	StrListPtr aTaxaList, txPtr, nLptr;
	PtrList nodeList, mrcaPtr;
	PtrList lnode;
	TREE thisTree;
	NODETYPE *mrca, *node;
	int i, ix=0, kix=-1,nList,counter=0, first_of_set, jix,flipCount,
	    last_of_set, watch_this_block=0, counter2=0, coin, cindex;
	aTaxaList=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	
			appendStrList(aTaxaList,aTokenPtr); /* store the label */
	nList=lengthList(aTaxaList);
	if (nList < 2)
	    fatal("Must have at least two names in CLADE_CHECK command");
	printf("cat phase2_header >> paup_phase2\n");
	if (gNexDataPtr->isTrees)
		{
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			++ix;
			++kix;
			if ( (ix-1)/NUMBER_IN_EFRON_SET == (float)(ix-1)/NUMBER_IN_EFRON_SET )
			    {
			    kix=0;
			    printf("# .................................\n");
			    }
			printf("# Tree %i:Specified group IS ",ix);
			thisTree=lnode->item;
			if (group_a_clade(thisTree->root, aTaxaList))
				{
				++counter;
				printf("a clade\n");
				isClade[kix]=1;

				}
			else
				{
				printf("NOT a clade\n");
				isClade[kix]=0;
				}
			if ( (ix)/NUMBER_IN_EFRON_SET == (float)(ix)/NUMBER_IN_EFRON_SET ) /* last of set */
			    	{
					flipCount=0;
					for (jix=0;jix<NUMBER_IN_EFRON_SET-1;jix++)
						if (isClade[jix] != isClade[jix+1])
							++flipCount;
					if (flipCount==1)
					  {
					  if (isClade[0]==0)
					    {
					    for (jix=0;jix<NUMBER_IN_EFRON_SET-1;jix++)
						if (isClade[jix]) /* here is the transition */
							{
							if (myRand()>0.5)coin=1;else coin=0;
							if (coin) 
								cindex=ix-NUMBER_IN_EFRON_SET+jix;
							else 
								cindex=ix-NUMBER_IN_EFRON_SET+jix+1;
							printf("agrep -d \';\' \'Weight set %li:\' theWsearchWeights >> paup_phase2\n", cindex);
				/* printf("cat PHASE2_INCLUDE >> paup_phase2\n");*/
							break;
							}
					    }
					  else
					    printf("# ONE CHANGE BUT WRONG DIRECTION\n");
					  }
					if (flipCount==0)
						printf("# NO changes\n");
					if (flipCount>1)
						printf("# WARNING! Uneven boundary conditions\n");

				}
			}
		}
	    printf("# Proportion (%i) of %i trees with group monophyletic (EFRON):%f\n",counter2, ix, 
			(float)NUMBER_IN_EFRON_SET*counter2/ix);
	    printf("cat phase2_tailer >> paup_phase2\n"); /*wrapper for nexus syntax */

	return;
    
    
}
static void doSimpleCladeCheckCommand()
/** Reads a list of taxa and checks to see if that group is a clade on all trees.
Reports the proportion of trees in which this group is a clade.

Writes a script in shell language to allow agreping from the weights file **/

{


	StrListPtr aTaxaList, txPtr, nLptr;
	PtrList nodeList, mrcaPtr;
	PtrList lnode;
	TREE thisTree;
	NODETYPE *mrca, *node;
	int i, ix=0, nList,counter=0, first_of_set, 
	    last_of_set, watch_this_block=0, counter2=0, coin, cindex;
	aTaxaList=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	
			appendStrList(aTaxaList,aTokenPtr); /* store the label */
	nList=lengthList(aTaxaList);
	if (nList < 2)
	    fatal("Must have at least two names in SIMPLE_CLADE_CHECK command");

	printf("echo \"#nexus\" > theNULLHweights\n");
	printf("echo \"begin bootstrap;\" >> theNULLHweights\n");
	printf("echo \"include file=PHASE1A_INCLUDE;\" >> theNULLHweights\n");
	printf("#Check for group monophyly...\n");
	if (gNexDataPtr->isTrees)
		{
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			++ix;
			printf("# Tree %i:Specified group IS ",ix);
			thisTree=lnode->item;
			if (group_a_clade(thisTree->root, aTaxaList))
				{
				++counter;
				printf("a clade\n");
				}
			else
				{
				printf("NOT a clade\n");
				printf("agrep -d \';\' \'Weight set %li:\' the_phase1_weights >> theNULLHweights\n", ix);
				}
			}
		}
	printf("# Proportion (%i) of %i trees with group monophyletic:%f\n",counter, ix, 
			(float)counter/ix);
	printf("echo \";end;\" >> theNULLHweights\n");

	return;
    
    
}
static void doB_Weights_Command(char *buffer, char *buffer1,char *buffer2)

/* For every weight statement in the input NEXUS file, we generate NREPS
new weight statements to allow further bootstrapping.  The weights can be real,
for example in EFRON resampling, which requires multinomial resampling*/

#define NREPS 100
#define MAX_WEIGHT_ARRAY 2500
{
	float XweightArray[MAX_WEIGHT_ARRAY];
	float weight;
	int character, j;
	char * dummy;
	for (j=0;j<MAX_WEIGHT_ARRAY;j++)
	    {
	    XweightArray[j]=0.0;
	    }
	while (!isEqual(aTokenPtr=nextToken(), ";"))
	    {
	    if (isEqual(aTokenPtr, ","))
		aTokenPtr=nextToken(); /* skip commas */
	    weight=strtod(aTokenPtr, &dummy);
	    aTokenPtr=nextToken();
	    if (!isEqual(aTokenPtr,":")) fatal("Improperly formatted b_weights statement");
	    aTokenPtr=nextToken();
	    character=strtod(aTokenPtr, &dummy);
	    if (character >=MAX_WEIGHT_ARRAY)
		fatal("Too many characters in B_WEIGHT_COMMAND: recompile with larger array");
	    XweightArray[character-1]=weight;
	    }	
	/* printf("TEST OF WEIGHTS COMMAND\n");
	for (j=0;j<MAX_WEIGHT_ARRAY;j++)
	    {
	    if ((j>0)&& ((j/10)==(j/10.0)))
		printf("\n");
	    
	    printf("%4.2f:%i, ", XweightArray[j],j+1);
	    }*/
	bshuf3(XweightArray, character, NREPS, buffer1,buffer2); /* assumes the last character read is the highest; ie. weights
		    read sequentially--violates NEXUS format,  but this command is only in my BOOT block anyway */
	if (buffer)
	    printf("%s\n", buffer);
	return;    
    
}
static void doEFRON_Weights_Command(char *buffer)

/* For every INTEGER weight statement in the input NEXUS file, we generate a new set of REAL weights on the
boundary between R1 and R2; see EFRON 1996 */

{
	int weightArray[MAX_WEIGHT_ARRAY];
	float weight;
	long index;
	int character, j;
	char * dummy;
	for (j=0;j<MAX_WEIGHT_ARRAY;j++)
	    {
	    weightArray[j]=0.0;
	    }
	while (!isEqual(aTokenPtr=nextToken(), ";"))
	    {
	    if (isEqual(aTokenPtr, ","))
		aTokenPtr=nextToken(); /* skip commas */
	    weight=strtod(aTokenPtr, &dummy);
	    aTokenPtr=nextToken();
	    if (!isEqual(aTokenPtr,":")) fatal("Improperly formatted b_weights statement");
	    aTokenPtr=nextToken();
	    character=strtod(aTokenPtr, &dummy);
	    if (character >=MAX_WEIGHT_ARRAY)
		fatal("Too many characters in EFRON_WEIGHT_COMMAND: recompile with larger array");
	    weightArray[character-1]=weight;
	    }
/* NB!  careful to make sure that efron1996 actually generates the right number(NUMBER_IN_EFRON_SET)  of increments ! */	
	efron1996(weightArray,character,NUMBER_IN_EFRON_SET,buffer, &index);

	return;    
    
}

static void doMRCACommand(void)

/** assigns an internal name to the MRCA of a set of taxa: 
Usage MRCA new_internal_name taxon1 ...taxonN ; if only one taxonname is given,  then 
assume it is an internal name and replace it with the newname
 **/

{
	StrListPtr aTaxaList, txPtr, nLptr;
StrListPtr mrcaTaxa;
	PtrList nodeList, mrcaPtr;
	PtrList lnode;
	TREE thisTree;
	NODETYPE *mrca, *node;
	char *new_internal_name, *old_internal_name;
	int i, ix=0, nList;
	gTaxaList=txPtr=aTaxaList=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	
			{
			appendStrList(aTaxaList,aTokenPtr); /* store the label */
			}
	nList=lengthList(aTaxaList);
	if (nList < 2)
		{
	    	doGenericAlert("Must have at least two names in MRCA command");
		return;
		}
	new_internal_name=aTaxaList->s; /* this is the first name in the list */
	nLptr=aTaxaList->next; /*points to taxon names*/
	if (gNexDataPtr->isTrees)
		{
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			if (nList == 2)
			    {
			    node=find_taxon_name(thisTree->root,nLptr->s);
			    if (node)
				{
				printf("Redefining node name: %s to %s\n",node->taxon_name,new_internal_name);
				setNodeName(node, new_internal_name);
				}
			    else
				{
				doGenericAlert("BAD MRCA COMMAND: Taxon name misspelled or not on tree");
				return;
				}
			    }
			else
			    {
			    mrca=MRCA(thisTree->root, nLptr);
			    if (mrca)
				{
				if (mrca->taxon_name)
				    if(!isEqual(mrca->taxon_name,"")) /* careful, name initialized to "" */
					{
					doGenericAlert("MRCA is overwriting an existing node name");
					printf("[** The overwritten node is %s **]\n",mrca->taxon_name);
					}
				setNodeName(mrca, new_internal_name);
				printf("Defining clade name: %s\n",new_internal_name);
				}
			    else
				{
				doGenericAlert("BAD MRCA COMMAND: Taxon name misspelled or not on tree");
				return;
				}
			    }
			}
		}


	return;
    
    
}

static void doBootCommand(StrListPtr fixedList, char* buffer)

/* Process the taxon bootstrap command, write the appropriate
NEXUS syntax to delete taxa and, also write stuff from a character
buffer that might have been included with the 'include' command */

{
	char  *dummy;
	int sample[MAX_TAXON_ARRAY], fixed[MAX_TAXON_ARRAY], included[MAX_TAXON_ARRAY], 
			    ntaxa=0, nrandom=0, nfixed=0,nsample=0,
			nstart=0,nstop=0,nstart2=0,nstop2=0,nrandom2=0;
	long aSeed=1, nreps, i, j, k;
	for (k=0;k<MAX_TAXON_ARRAY;k++)
		fixed[k]=0; /* necessary default for later 
				calls to 'taxon_sample' */
	if (fixedList) /* if there are fixed taxa */
	    {
	    nfixed=lengthList(fixedList);
	    for (k=0;k<nfixed;k++)
		fixed[k]=strtod(getkthStr(fixedList, k+1), &dummy);
		    /* internal representation of taxon ids is on 1..n */
	    }
		
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NSTART"))
			nstart=strtod(LocalToken,&dummy);
		if (parse_assignment2("NSTOP"))
			nstop=strtod(LocalToken,&dummy);
		if (parse_assignment2("NSTART2"))
			nstart2=strtod(LocalToken,&dummy);
		if (parse_assignment2("NSTOP2"))
			nstop2=strtod(LocalToken,&dummy);
		if (parse_assignment2("NREPS"))
			nreps=strtod(LocalToken,&dummy);
		if (parse_assignment2("NRANDOM"))
			nrandom=strtod(LocalToken,&dummy);
		if (parse_assignment2("NRANDOM2"))
			nrandom2=strtod(LocalToken,&dummy);
		if (parse_assignment2("NTAXA"))
			ntaxa=strtod(LocalToken,&dummy);
		if (parse_assignment2("SEED"))
			aSeed=strtod(LocalToken,&dummy);
		}
	nsample=nfixed+nrandom+nrandom2; /* total number to be in sample */
	srand(aSeed);
	printf("begin paup;\n");
	for (k=1;k<=nreps;k++)
	    {
	    for (i=0;i<ntaxa;i++)
		included[i]=0;
	    taxon_sample(ntaxa,nfixed, nrandom, fixed, sample,
			nstart,nstop,nstart2,nstop2,nrandom2);
	    for (j=0;j<nsample;j++)
		    included[sample[j]-1]=1;
	    printf("[The taxon sample is:");
	    for (j=0;j<ntaxa;j++)
		if (included[j]) 
		    printf("%i ", j+1); 
	    printf("]\n");
	    printf("delete ");
	    for (j=0;j<ntaxa;j++)
		if (!included[j]) 
		    printf("%i ", j+1); /* +1 is to reconvert back to external
					representation of taxon ids on 1..n */
	    printf("/prune=yes;\n");
	    if (buffer)
		    printf("%s\n", buffer);

	    }
	printf("end;\n");
	return;
}
static void doBootCharCommand(char* buffer)

/* Process the character bootstrap command, write the appropriate
NEXUS syntax to weight characters and, also write stuff from a character
buffer that might have been included with the 'include' command */


{
	char  *dummy;
	double u, u1=0.0, u2=0.0, a;
	int *weightArray;
	float * pMean,pSum=0.0;
	long aSeed=1, nreps=0, i, j, k,nchars=0, index=0,ix;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NREPS"))
			nreps=strtod(LocalToken,&dummy);
		if (parse_assignment2("NCHARS"))
			nchars=strtod(LocalToken,&dummy);
		if (parse_assignment2("SEED"))
			aSeed=strtod(LocalToken,&dummy);
		}
	srand(aSeed);
	if (nchars>0)
		{
		pMean=(float *)myMalloc((nchars)*sizeof(float));
		for (ix=1;ix<nchars;ix++)
			pMean[ix]=0.0;
		weightArray=(int *)myMalloc(nchars*sizeof(int));
		printf("begin paup;\n");
		for (k=1;k<=nreps;k++)
		    {
		    u1=0.0;
		    u2=0.0;
		    ++index;
		    bshuf2(weightArray,nchars);
		    printf("\n\n[******************************]\n");
		    printf("[*** Bootstrap replicate %i ***]\n\n", k);
		/*    for (j=0;j<nchars;j++)
			    {
			    u=weightArray[j]-1.0; 
			    u1+=u*u*u;
			    u2+=u*u;
			    } *//* for EFRON 96 algorithm */
		    a= (1/6.0)*u1/pow(u2, 1.5);
		    printf("[Weight set %li:][w=1.0][a=%f]weights ", index, a);
		    for (j=0;j<nchars-1;j++)
			    {
			    if ((j>0)&& ((j/10)==(j/10.0)))
				printf("\n");
			    printf("%i:%i, ", weightArray[j],j+1);
			    }
		    printf("%i:%i;\n", weightArray[nchars-1],nchars);
		    if (buffer)
			    printf("%s\n", buffer);
		   for (j=0;j<nchars;j++)
			    {
			    pMean[j]+=weightArray[j]/(float)nreps;
			    }
	
	/*	    if (EFRON1996)
			efron1996(weightArray,nchars,NUMBER_IN_EFRON_SET,buffer, &index);*/

		    }
		printf("[Mean vector\n");
		for (j=0;j<nchars;j++)
			{
			if ((j>0)&& ((j/10)==(j/10.0)))
			    printf("\n");
			printf("%f:%i, ", pMean[j],j+1);
			pSum+=pMean[j];
			}
		printf(" sum=%f]\n",pSum);
		printf("end;\n");
		myFree(weightArray);
		}
	return;
}
static void efron1996(int *weightArray,int nchars,int num_points, char *buffer, 
    long *index)

/* writes PAUP code to implement part of Efron et al., 1996 boot algorithm.
 * This component writes commands that generate a set of weight statements
 * corresponding to the search for the boundary vectors.  Calculates proper
 * weights for all points,  w,  such that w is an element of [0,xinc, 2*xinc, 
 * ..., 1].  For each point the p vectors wP(j)+(1-w)P(cent) are calculated 
 * (see bottom left of p. 7089 of paper).  A row of weights is printed along
 * with any commands stored in the buffer from a previous include.
 *
 * 'weightArray' contains the bootstrap weight vector, P(j) on entry.
 * Spits out NUMBER_IN_EFRON_SET rows of weights.
 */

{
double w, u, u1=0.0, u2=0.0, a, pj,xinc;
int j,k;
if (num_points<2) fatal("Too few points in efron1996");
xinc=1/(num_points-1.0);
printf("begin paup;\n");
for (k=1;k<=num_points;k++) 
    {
    w=1-(k-1)*xinc;
    u1=0.0;
    u2=0.0;
    ++(*index);
    printf("\n[Weight set %li:][w=%f] ", *index, w);
    for (j=0;j<nchars;j++)
	    {
	    u= (1-w)+w*weightArray[j]-1.0;
	    u1+=u*u*u;
	    u2+=u*u;
	    }
    a= (1/6.0)*u1/pow(u2, 1.5);
    printf(" [a=%f] ", a);
    printf("weights ");
    for (j=0;j<nchars-1;j++)
	    {
	    if ((j>0)&& ((j/30)==(j/30.0)))
		printf("\n");
	    pj=(1-w)+w*weightArray[j];
	    
	    printf("%4.2f:%i, ", pj,j+1);
	    }
    pj=(1-w)+w*weightArray[nchars-1];
    printf("%f:%i;\n", pj,nchars);
    if (buffer)
	    printf("%s\n", buffer);
    printf("\n");
    } 
printf("end;\n");
return;
}


static char * doIncludeCommand(void)

/* Reads the specified file and passes it to the 'boot' command so 
that its text can be printed for each boot replicate.
NB!  The include must be executed BEFORE THE BOOT command.
AARRPHHH: demands that filenames be in caps!*/

{
	FILE* fpntr;
	char  *filename, *buffer;
	buffer=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("FILE"))
			filename=DupStr(LocalToken);
		}
	if (  (fpntr=fopen(filename,"r")) )
		{
		buffer=slurpFile (fpntr, 10000);
		}
	else
	    printf("Failed to open file '%s'\n", filename);
	return buffer;
}

/*
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		appendStrList(gNexDataPtr->TaxaList,aTokenPtr);
		}

*/
static StrListPtr doFixedTaxaListCommand(void)
{
	StrListPtr aTaxaList; /* probably just numbers, but treat them
		as strings */

	aTaxaList=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	
			appendStrList(aTaxaList,aTokenPtr); /* store the label */
	return aTaxaList;

}



/****************************************************************/

void doTaxaBlock(void)
{
	do 				/* need to put in error checking in case no DIMENSIONS statement */
		{
		aTokenPtr=nextToken();
		if (isEqual(aTokenPtr,"DIMENSIONS"))  
			doTaxDimensions();
		if (isEqual(aTokenPtr,"TAXLABELS"))
			doTaxLabels();

		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));
	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
return;
}
/**************************************************************/
void doTreeBlock(void)
{
	do 
		{
		aTokenPtr=nextToken();
		if (isEqual(aTokenPtr,"TREE") || isEqual(aTokenPtr,"UTREE") )  
									/* process TREE command */
			doTreeCommand();
		if (isEqual(aTokenPtr,"TRANSLATE") )  
									/* process TREE command */
			doTranslateCommand();
			
		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));

	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
if (gNexDataPtr->isTrees)
	{
	/*printf("[Number of trees read = %i]\n",gNexDataPtr->NumTrees);		
	print_tree_list(gNexDataPtr->inTrees);*/
	}
return;
}
/**************************************************************/
void doSimBlock(void)
{
	do 
		{
		aTokenPtr=nextToken();
		if (isEqual(aTokenPtr,"SIMULATE") )  
									/* process TREE command */
			doSimCommand();
		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));

	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
return;
}
/****************************************************************/

void doDivTimeCommand(void)
{
extern int powellMode;
extern double gGamma_c;
extern double gGamma_b;
extern int gisConstrained, gVarMinFlag,gEstRoot,gFloatRoot;
	double (*obj_func_array[10])(double[]); /* array of pointers to the various
		objective functions...indexed by 'method' which is set up below */
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD ,*method_string;
double  EstMult=1,PrdMult=1,cutoff=2.0;
static	int method=LaF,
	    algorithm=POWELL;
static	long iTree=0;	
	long numTrees;
	int j,success,crossv=0,crossv2=0,cvSample=0,nRates=1,maxBisect=20,fossilFlag=0,fossilFixedFlag=0;
	NODETYPE *root, *found_node;
static	double cvStart=0.0,cvInc=1.0,cvNum=1;
	int confidence=0;
StackPtr S;

	powellMode=1;
	method_string=NULL; /*default*/
/*
	obj_func_array[LaF]=objLangFitch;
	obj_func_array[NP]=objNP;
	obj_func_array[GAMMA]=objGamma;
	obj_func_array[PENLIKE]=objPenLike;
	obj_func_array[PENLIKET]=objPenLikeT;
*/

#define POWELL_STACK_SIZE 35
gFStack=newStack(POWELL_STACK_SIZE);
gPStack=newStack(POWELL_STACK_SIZE);


	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("ALGORITHM"))
			{
			if (isEqual(LocalToken,"POWELL"))
				algorithm=POWELL;
			if (isEqual(LocalToken,"QNEWT"))
				algorithm=QNEWT;
			if (isEqual(LocalToken,"TN"))
				algorithm=TN;
			}
		if (parse_assignment2("METHOD"))
			{
			method_string=DupStr(LocalToken);
			}
		if (parse_assignment2("TAXON"))
			{
			taxon=DupStr(LocalToken);
			}
		if (parse_assignment2("TREE"))
			{
			iTree=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("NRATES"))
			nRates=strtod(LocalToken,&dummy);
		if (parse_assignment2("MAXBISECT"))
			maxBisect=strtod(LocalToken,&dummy);
		if (parse_assignment2("CONFIDENCE"))
			{
			if (isEqual(LocalToken,"YES"))
				confidence=1;
			else
				confidence=0;
			}
		if (parse_assignment2("CROSSV"))
			{
			if (isEqual(LocalToken,"YES"))
				crossv=1;
			else
				crossv=0;
			}
		if (parse_assignment2("FOSSILFIXED"))
			{
			if (isEqual(LocalToken,"YES"))
				fossilFixedFlag=1;
			else
				fossilFixedFlag=0;
			}
		if (parse_assignment2("FOSSILCONSTRAINED"))
			{
			if (isEqual(LocalToken,"YES"))
				fossilFlag=1;
			else
				fossilFlag=0;
			}
		if (parse_assignment2("CUTOFF"))
			cutoff=strtod(LocalToken,&dummy);
		if (parse_assignment2("CVSTART"))
			cvStart=strtod(LocalToken,&dummy);
		if (parse_assignment2("CVINC"))
			cvInc=strtod(LocalToken,&dummy);
		if (parse_assignment2("CVNUM"))
			cvNum=strtod(LocalToken,&dummy);
		if (parse_assignment2("CVSAMPLE"))
			cvSample=strtod(LocalToken,&dummy);
		
		}
	if (method_string!=NULL)
		{
		if (isEqual(method_string,"PL"))method=PENLIKE; 
		if (isEqual(method_string,"LF"))method=LaF; 
		if (isEqual(method_string,"NP") || isEqual(method_string,"NPRS"))
		    {
		    method=NP;
		    gVarMinFlag=0; /* this is a kludgy way to distinguish between two flavors
					of NPRS alogorithm */
		    }
		if (isEqual(method_string,"NPVAR"))
		    {
		    method=NP;
		    gVarMinFlag=1;
		    }

		}

		if ((method==LaF)  && (nRates > 1) )method=LFLOCAL; /* temp..I can just use this always for LF */



	if (gNexDataPtr->isTrees)
	    {
	      lnode=gNexDataPtr->inTrees;
	      if (iTree>0) /* a specific tree was indicated */
		{
		if (iTree > pLengthList(lnode))
			{
			doGenericAlert("Invalid tree specified");
			return;
			}
		else
			{
			thisTree=(pListgetkthNode(lnode,iTree))->item;
			doObjFunc(thisTree,method,nRates,algorithm,&success);
			}
		}
	      else
		{
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			if (crossv)
				{
				if (fossilFlag)
					doFossilCrossV(thisTree,method,nRates,algorithm,cvStart,cvInc,cvNum,cvSample);				
				else if (fossilFixedFlag)
					doFossilCrossVfixed(thisTree,method,nRates,algorithm,cvStart,cvInc,cvNum,cvSample);				
				else
					doCrossV(thisTree,method,nRates,algorithm,cvStart,cvInc,cvNum,cvSample);
				}
			else
				{
				doObjFunc(thisTree,method,nRates,algorithm,&success);
				if (confidence)
					{
					doConfidence(thisTree,taxon,method,nRates,algorithm,cutoff,maxBisect);
					}
				}
			}
		}
	      if (method==GAMMA)
				{
				thisTree->est_b=gGamma_b;
				thisTree->est_c=gGamma_c;
				}
	    }
	if (method_string)
		myFree(method_string);
	freeStack(gFStack);
	freeStack(gPStack);
	return;						
}

/***************/
static void
doConfidence(TREE T,char * nodeName,int method,int nRates,int algorithm,double cutoff,int JMAX)

/*
Find the confidence interval on an estimated node time. Construct the interval by
finding the values of that node time at which the likelihood drops by an amount 'cutoff'.
This is done by examining a range of possible times roughly between an upper and lower bound determined
from the age constraints and fixed node times. The focal node is fixed at various times across this
range, and the search is restarted (cf Cutler, MBE, 2000) estimating all other parameters as before. Rather
than search the whole range, a bisection strategy is used with an NRC function.


nodeName -	Determine confidence interval for this node
cutoff	--	Target of the search is (Max Like - cutoff)
JMAX	--	Maximum number of bisections allowed
*/


{
int i,maxPts=10,success,j;
NODETYPE *n;
double upper,lower,R,t,tSoln,factor,solnObj,targetObj,bump,low,high,tLow,objLow,tHigh,objHigh;
double x1,x2,dx,f,fmid,xmid,rtb,xacc;
solnObj=T->obj;	/* value of the obj function at the solution */
targetObj=solnObj-cutoff;
if (!(n=find_taxon_name(T->root,nodeName)))
	{
	doGenericAlert("Failed to find node name in 'confidence'");
	return;
	}
if (!isFree(n))
	{
	doGenericAlert("Cannot estimate confidence limit on FIXED node");
	return;
	}
tSoln=n->time;			/* save the estimated age of node */	
upper=nodeUpperBound(n);
if (upper>=1e20)
	upper=2*n->time;	/* if no upper bound, arbitrarily put it at 2X the node's age, but
					we'll check to see if this accomodates results */
lower=nodeLowerBound(n);
factor=0.9;			/* let's squeeze the search interval by this amount to prevent bumping
					against bounds! */
R=upper-lower;
bump=R*(1-factor)/2;
lower+=bump;
upper-=bump;
low=lower;
high=tSoln;
xacc=(upper-lower)*0.01;
x1=low;
x2=tSoln;
	{/* modifed from NRC 'rtbis' */
	n->free=0;
	n->time=x1;
	doObjFunc(T,method,nRates,algorithm,&success);
	f=T->obj - targetObj;
	fmid=solnObj - targetObj;
	if (f*fmid >= 0.0)
		{
		doGenericAlert ("Confidence search failed: no crossover point");return;
		}
	rtb = f < 0.0 ? (dx=x2-x1,x1):(dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++)
		{
		xmid=rtb+(dx*=0.5);
		n->time=xmid;
		doObjFunc(T,method,nRates,algorithm,&success);
		fmid=T->obj - targetObj;
		if (fmid <= 0.0) rtb=xmid;
		if (/*  fabs(dx) < xacc || fmid == 0.0 */ fabs(fmid) < 0.1) /* my termination criterion! */
			{
			/* printf ("**** Lower t = %f objective function - targetObj= %f [iters=%i]\n",rtb,fmid,j);*/
			tLow=rtb;
			objLow=fmid+targetObj;
			break;
			}
		}
	if (j>=JMAX)
		doGenericAlert("Confidence search failed");
	}
x2=tSoln;	/* lazy-ass copy of the previous code! */
x1=upper;
	{/* modifed from NRC 'rtbis' */
	n->free=0;
	n->time=x1;
	doObjFunc(T,method,nRates,algorithm,&success);
	f=T->obj - targetObj;
	fmid=solnObj - targetObj;
	if (f*fmid >= 0.0)
		{
		doGenericAlert ("Confidence search failed: no crossover point");return;
		}
	rtb = f < 0.0 ? (dx=x2-x1,x1):(dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++)
		{
		xmid=rtb+(dx*=0.5);
		n->time=xmid;
		doObjFunc(T,method,nRates,algorithm,&success);
		fmid=T->obj - targetObj;
		if (fmid <= 0.0) rtb=xmid;
		if (/*  fabs(dx) < xacc || fmid == 0.0 */ fabs(fmid) < 0.1) /* my termination criterion! */
			{
			/* printf ("**** Higher t = %f objective function - targetObj= %f [iters=%i]\n",rtb,fmid,j);*/
			tHigh=rtb;
			objHigh=fmid+targetObj;
			break;
			}
		}
	if (j>=JMAX)
		doGenericAlert("Confidence search failed");
	}
printf("\nConfidence interval for node %s using cutoff value of %f\n",n->taxon_name,cutoff);
printf("Point\t\tAge\t\tObj\n");
printf("Lower\t\t%6.2f\t\t%6.2f\n",tLow,objLow);
printf("Higher\t\t%6.2f\t\t%6.2f\n",tHigh,objHigh);
printf("Soln\t\t%6.2f\t\t%6.2f\n",tSoln,solnObj);



n->free=1;
return;
}
/***************/

static float 
doCrossV(TREE tree, int method,int nRates,int algorithm,double cvStart,double cvInc,double cvNum, int cvSample)

/*  
	Does a cross validation analysis in which we 
	(1)sequentially remove each tip (leaving the tip's ancestor in place),
	(2)do a full estimation on remaining subtree, 
	(3)then calculate a prediction error for that removed terminal branch. 
	(4) Then puts the terminal back on the tree.

	If the method is LaF or NP then one round of CV is invoked.
	If the method is PENLIKE, then analysis is repeated with the smoothing parameter chosen from a range
	from [cvStart, cvStart+cvInc, ...,cvStart+cvInc*(cvNum-1)].
	If cvSample==0 then we cross validate on all the taxa; if cvRep>0, we randomly sample that many taxa
		and use only those taxa. We use the same random sample for ALL smoothing levels, however!


*/
{
char *Result, *Good="Good", *Failed="Failed";
int i,j,k,success,collFlag=0,ntips,verbose,overallGood=1;
double * cvScore,*chiSqScore, cvSum,chiSq,chiSqSum,*cvTotalScore,*cvTotalScoreChiSq,bestChiSq,bestSmooth;
int * sample, *cvResult, *cvResultFinal;
int numSuccess, numFail,bestJ;
double smooth;
PtrList tipNodeList;
NODETYPE *CVNode;
ntips=numdesc(tree->root);
cvResult = (int *)myMalloc((ntips+1)*sizeof(int));
cvResultFinal = (int *)myMalloc((cvNum+1)*sizeof(int));
cvScore = vector(1,ntips);
chiSqScore = vector(1,ntips);
cvTotalScore=vector(1,cvNum);
cvTotalScoreChiSq=vector(1,cvNum);
tipNodeList = pNewList();
TreeToTaxaPtrList(tree->root,tipNodeList); /* get a list of all the tip nodes */

/* NOTES: This will not work when tips have ages > 0. In that case, we must enforce constraints on that tips
	ancestor node; otherwise, once it is pruned, that ancestor might be inferred to be younger than the 
	pruned tip, causing all hell to break loose when doing the predicted value */

if (cvSample>0)
	{
	srand(gNexDataPtr->RateBlockParms.seed);	/* sets up random seed */
	sample=taxon_sample_simple(ntips,cvSample);
	ntips=cvSample;
	}

if (method==LaF || method==NP)
	cvNum=1;	/* just overrides the default value of cvNum for these methods, forcing them to only go once */
verbose=gNexDataPtr->RateBlockParms.verbose;
if (verbose > 0)
	printf("Begin cross-validation analyses...\n");
for (j=0;j<cvNum;j++)
	{
	numSuccess=0;
	numFail=0;
        cvSum=0.0;
	chiSqSum=0.0;
	smooth=pow(10.0,j*cvInc+cvStart);
	gNexDataPtr->RateBlockParms.smoothing=smooth;
	for (k=1;k<=ntips;k++)
		{
		if (cvSample>0)
			i=sample[k-1];
		else
			i=k;
		CVNode=(NODE)(pListgetkthNode(tipNodeList,(long)i)->item);
		RemoveTaxonLvAnc(CVNode);
		gNexDataPtr->RateBlockParms.verbose=0; 
	/* suppress all output from the actual optimization run...may want to allow it for debugging though! */
		doObjFunc(tree,method,nRates,algorithm,&success);
		gNexDataPtr->RateBlockParms.verbose= verbose; /* restore output verbositude to current value */
		AddChild(CVNode->anc,CVNode); /* important to reattach before next call */
		if (success)
			{
			cvResult[k]=1; /* good */
			++numSuccess;				
			cvScore[k]=cvSquareErrorBranch(tree,CVNode,method,&chiSq); 
			cvSum+=cvScore[k];
			chiSqScore[k]=chiSq;
			chiSqSum+=chiSq;
			printf("+\n");
			}
		else
			{
			cvResult[k]=0; /*failed */
			cvScore[k]=0.0;
			chiSqScore[k]=0.0; /* just put default values in */
			++numFail;
			printf("-\n");
			}
		}
	cvTotalScore[j+1]=cvSum;
	cvTotalScoreChiSq[j+1]=chiSqSum;
	if (numFail==0)
		cvResultFinal[j+1]=1; /* all prunings led to successful optimizations */
	else
		cvResultFinal[j+1]=0; /* some prunings had failed optimizations */
	printf("\n");
	if (verbose>0)
		{	
		printf(".....................................................................\n");
		printf("\nCV Results for smoothing = %f:\nPruned taxon\tSq\t\tChiSq\t\tResult\n",smooth);
		for (i=1;i<=ntips;i++)
			{
			CVNode=(NODE)(pListgetkthNode(tipNodeList,(long)i)->item);
			if (cvResult[i])
				Result=Good;
			else
				Result=Failed;
			printf("%8.8s\t%8.2f\t%8.2f\t%s\n",CVNode->taxon_name,cvScore[i],chiSqScore[i],Result);
/*			printf ("Cross Validation Score [%2i] = %f\t[%f]\n",i,cvScore[i],chiSqScore[i]);
		printf("Cross Validation Score Total (%i pruned terminals):smoothing = %f CV=%f chiSq=%f\n",numSuccess,smooth,cvSum/numSuccess,chiSqSum/numSuccess);
*/
			}
		if (numFail>0)
			printf("** Note that %i failed prunings occurred **\n",numFail);
		}

/* REMEMBER WE ARE OFTEN NOT COUNTING A TIP DESCENDED FROM THE ROOT */

	}
printf("********************************************************************************\n\n");
printf("Results of cross validation analysis for tree %s\n",tree->name);
  switch (method)
	{
	case PENLIKE:printf("Method = Penalized Likelihood\n");break;
	case LaF:printf("Method = Langley and Fitch\n");break;
	case LFLOCAL:printf("Method = Langley and Fitch (with %i local rates)\n",nRates);break;
	case NP:printf("Method = Non-parametric\n");
	}
  switch (algorithm)
	{
	case POWELL: printf("Optimization via Powell's method\n");break;
	case QNEWT:  printf("Optimization via quasi-Newton method with analytical gradients\n");
	}
printf("\nlog10\n");
printf("smooth\tsmooth\t\tSq Error\tChi Square Error\n");
printf("--------------------------------------------------------------------------------\n");
overallGood=1;
bestChiSq=1e20;
bestJ=0;
for (j=0;j<cvNum;j++)
	{
	if (cvResultFinal[j+1])
		{
		if (cvTotalScoreChiSq[j+1]<bestChiSq)
			{
			bestChiSq=cvTotalScoreChiSq[j+1];
			bestJ=j;
			}
		Result=Good;
		}
	else
		{
		overallGood=0;
		Result=Failed;
		}
	smooth=pow(10.0,j*cvInc+cvStart);
	printf("%6.2f\t%6.2g\t\t%6.2f\t%6.2f\t(%s)\n",j*cvInc+cvStart,smooth,cvTotalScore[j+1],cvTotalScoreChiSq[j+1],Result);
	}
printf("********************************************************************************\n\n");

bestSmooth=pow(10.0,bestJ*cvInc+cvStart);
printf("Optimum: %6.2f\t%6.2g\t\t%6.2f\t%6.2f\n",bestJ*cvInc+cvStart,bestSmooth,cvTotalScore[bestJ+1],cvTotalScoreChiSq[bestJ+1]);
if (!overallGood)
	printf("WARNING: Cross validation procedure had errors: optimum may be incorrect\n");

printf("********************************************************************************\n\n");
myFree(cvResult);
myFree(cvResultFinal);
free_vector(cvScore,1,ntips);
free_vector(chiSqScore,1,ntips);
free_vector(cvTotalScore,1,cvNum);
free_vector(cvTotalScoreChiSq,1,cvNum);
freepList(tipNodeList);
return bestSmooth;
}

/***************/




static void 
doFossilCrossV(TREE tree, int method,int nRates,int algorithm,double cvStart,double cvInc,double cvNum, int cvSample)

/*  
	Does a cross validation analysis in which we 
	(1)sequentially unconstrain each node that has a constraint (fixed nodes are not affected)
	(2)do a full estimation on the tree, 
	(3)then calculate the deviation of the estimate for that node versus the constraint (if the constraint is now violated) 
	(4) sums these errors across all constrained nodes

	If the method is LaF or NP then one round of CV is invoked.
	If the method is PENLIKE, then analysis is repeated with the smoothing parameter chosen from a range
	from [cvStart, cvStart+cvInc, ...,cvStart+cvInc*(cvNum-1)].

	Reports two kinds of error, a fractional value per constrained node, and a raw value per constrained node in units of time.

*/
{
char *Result, *Good="Good", *Failed="Failed";
int i,j,k,success,collFlag=0,ntips,verbose,numFixed,numConstrained,numNodes,curIndex;
double * cvScore,*cvScoreRaw, cvSum,chiSq,cvRawSum,*cvTotalScore,*cvTotalScoreRaw,*fixedTime, *estTime;
int * sample, *cvResult, *cvResultFinal;
int numSuccess, numFail,wasConstrainedMin,wasConstrainedMax,wasFixed,wasConstrained;
double smooth,saveTime;
PtrList tipNodeList;
NODETYPE *CVNode;
NODE *nodeArray;
ntips=numdesc(tree->root);
cvResult = (int *)myMalloc((ntips+1)*sizeof(int));
cvResultFinal = (int *)myMalloc((cvNum+1)*sizeof(int));
cvScore = vector(1,ntips);
cvScoreRaw = vector(1,ntips);
cvTotalScore=vector(1,cvNum);
cvTotalScoreRaw=vector(1,cvNum);
tipNodeList = pNewList();
TreeToTaxaPtrList(tree->root,tipNodeList); /* get a list of all the tip nodes */


if (method==LaF || method==NP)
	cvNum=1;	/* just overrides the default value of cvNum for these methods, forcing them to only go once */



numFixed=numFixedNodes(tree->root);
numConstrained=numConstrainedNodes(tree->root);
numNodes=numFixed+numConstrained;
if (numFixed <1 || numConstrained <2)
	{
	doGenericAlert("Must have at least one fixed and two constrained nodes for fossil cross validation");
	return;
	}
fixedTime=(double *)myMalloc(numNodes*sizeof(double));
estTime=(double *)myMalloc(numNodes*sizeof(double));

	// makes a single array beginning with fixed nodes and ending with constrained ones (if any)
nodeArray=(NODE *)myMalloc(numConstrained*sizeof(NODE));
curIndex=0;
//setupFixedNodeArray(tree->root, nodeArray, &curIndex);

setupConstrainedNodeArray(tree->root, nodeArray, &curIndex); 


verbose=gNexDataPtr->RateBlockParms.verbose;
if (verbose > 0)
	printf("Begin fossil cross-validation analyses...\n");
for (j=0;j<cvNum;j++)
	{
	numSuccess=0;
	numFail=0;
    cvSum=0.0;
	cvRawSum=0.0;
	smooth=pow(10.0,j*cvInc+cvStart);
	gNexDataPtr->RateBlockParms.smoothing=smooth;
	for (k=1;k<=numConstrained;k++)
		{
		i=k-1;
		CVNode=nodeArray[i];
		
		if (isConstrainedMax(CVNode))
			{
			wasConstrainedMax=1;
			CVNode->nodeIsConstrainedMax=0;
			}
		else
			wasConstrainedMax=0;
		if (isConstrainedMin(CVNode))
			{
			wasConstrainedMin=1;
			CVNode->nodeIsConstrainedMin=0;
			}
		else
			wasConstrainedMin=0;
			
		gNexDataPtr->RateBlockParms.verbose=0; 
	/* suppress all output from the actual optimization run...may want to allow it for debugging though! */
		doObjFunc(tree,method,nRates,algorithm,&success);
		gNexDataPtr->RateBlockParms.verbose= verbose; /* restore output verbositude to current value */
		if (success)
			{
			estTime[i]=CVNode->time;
			cvResult[k]=1; /* good */
			++numSuccess;
	// Below I deal with all three possible cases: the two simple constraints, or both together;
	// The latter is handled by setting up the next condition and then going through both tests below it
	// Notice that any given time cannot violate both constraints simultaneously

			cvScore[k]=0.0;
			cvScoreRaw[k]=0.0;

			if (wasConstrainedMin)
				{
				if (estTime[i]<CVNode->minAge) // constraint violated, so calculate departure
					{
					cvScore[k]+=2*fabs(CVNode->minAge-estTime[i])/(CVNode->minAge+estTime[i]); 
					cvScoreRaw[k]+=fabs(CVNode->minAge-estTime[i]);
					}
				}
			if (wasConstrainedMax)
				{
				if (estTime[i]>CVNode->maxAge) // constraint violated, so calculate departure
					{
					cvScore[k]+=2*fabs(CVNode->maxAge-estTime[i])/(CVNode->maxAge+estTime[i]); 
					cvScoreRaw[k]+=fabs(CVNode->maxAge-estTime[i]);
					}
				}
			cvRawSum+=cvScoreRaw[k];
			cvSum+=cvScore[k];
			}
		else
			{
			cvResult[k]=0; /*failed */
			cvScore[k]=-99.9;
			cvScoreRaw[k]=0.0; /* should probably set this to some bogus value! */
			++numFail;
			}

	// Restore constraints for this node
		if (wasConstrainedMax)
			CVNode->nodeIsConstrainedMax=1; // Notice, unconstraining does not delete constraint times in struct
		if (wasConstrainedMin)
			CVNode->nodeIsConstrainedMin=1;

		}
		
	printf("\nFossil-constrained cross validation analysis\n\n");
	printf("\tNode\t\tEst Age\t\tMin\tMax\t\tFract Score\tRaw Score\n");
	printf("----------------------------------------------------------------------------------------\n");
	for (i=0;i<numConstrained;i++)
		{
		CVNode=nodeArray[i];
		printf("%i\t%.6s\t\t%6.2f\t\t",i+1,CVNode->taxon_name,estTime[i]);
		if (isConstrainedMin(CVNode))
			printf("%6.2f\t",CVNode->minAge);
		else
			printf("   --   ");
		if (isConstrainedMax(CVNode))
			printf("%6.2f\t",CVNode->maxAge);
		else
			printf("   --   ");
		printf("\t%f\t%6.2f\n",cvScore[i+1],cvScoreRaw[i+1]);
		}

	cvTotalScore[j+1]=cvSum/(numNodes-numFail); /* on the off chance that some reps failed, don't count their contributions */
	cvTotalScoreRaw[j+1]=cvRawSum/(numNodes-numFail);
	
	if (numFail==0)
		cvResultFinal[j+1]=1; /* all prunings led to successful optimizations */
	else
		cvResultFinal[j+1]=0; /* some prunings had failed optimizations */

/* REMEMBER WE ARE OFTEN NOT COUNTING A TIP DESCENDED FROM THE ROOT */

	}
printf("********************************************************************************\n\n");
printf("Summary results of fossil-constrained cross validation analysis for tree %s\n",tree->name);
  switch (method)
	{
	case PENLIKE:printf("Method = Penalized Likelihood\n");break;
	case LaF:printf("Method = Langley and Fitch\n");break;
	case LFLOCAL:printf("Method = Langley and Fitch (with %i local rates)\n",nRates);break;
	case NP:printf("Method = Non-parametric\n");
	}
  switch (algorithm)
	{
	case POWELL: printf("Optimization via Powell's method\n");break;
	case QNEWT:  printf("Optimization via quasi-Newton method with analytical gradients\n");
	}
printf("\nFixed nodes:%i\nConstrained nodes:%i\n",numFixed,numConstrained);
printf("\nlog10\n");
printf("smooth\tsmooth\t\tFract Error\tRaw Error\n");
printf("--------------------------------------------------------------------------------\n");
for (j=0;j<cvNum;j++)
	{
	if (cvResultFinal[j+1])
		Result=Good;
	else
		Result=Failed;
	smooth=pow(10.0,j*cvInc+cvStart);
	printf("%6.2f\t%6.2g\t\t%6.4f\t\t%6.4f\t(%s)\n",j*cvInc+cvStart,smooth,cvTotalScore[j+1],cvTotalScoreRaw[j+1],Result);
	}
printf("********************************************************************************\n\n");
myFree(cvResult);
myFree(cvResultFinal);
myFree(fixedTime);
myFree(estTime);
myFree(nodeArray);
free_vector(cvScore,1,ntips);
free_vector(cvScoreRaw,1,ntips);
free_vector(cvTotalScore,1,cvNum);
free_vector(cvTotalScoreRaw,1,cvNum);
freepList(tipNodeList);


return;
}


/****************************************************************/
static void 
doFossilCrossVfixed(TREE tree, int method,int nRates,int algorithm,double cvStart,double cvInc,double cvNum, int cvSample)

/*  
	Does a cross validation analysis in which we 
	(1)sequentially fix each node that is fixed 
	(2)do a full estimation on the tree, 
	(3)then calculate the deviation of the estimate for that node versus the original fixed value 
	(4) sums these errors across all fixed nodes 

	If the method is LaF or NP then one round of CV is invoked.
	If the method is PENLIKE, then analysis is repeated with the smoothing parameter chosen from a range
	from [cvStart, cvStart+cvInc, ...,cvStart+cvInc*(cvNum-1)].

	Reports two kinds of error, a fractional value per constrained node, and a raw value per constrained node in units of time.

*/
{
char *Result, *Good="Good", *Failed="Failed";
int i,j,k,success,collFlag=0,ntips,verbose,numFixed,numConstrained,numNodes,curIndex;
double * cvScore,*cvScoreRaw, cvSum,chiSq,cvRawSum,*cvTotalScore,*cvTotalScoreRaw,*fixedTime, *estTime;
int * sample, *cvResult, *cvResultFinal;
int numSuccess, numFail,wasConstrainedMin,wasConstrainedMax,wasFixed,wasConstrained;
double smooth,saveTime;
PtrList tipNodeList;
NODETYPE *CVNode;
NODE *nodeArray;
ntips=numdesc(tree->root);
cvResult = (int *)myMalloc((ntips+1)*sizeof(int));
cvResultFinal = (int *)myMalloc((cvNum+1)*sizeof(int));
cvScore = vector(1,ntips);
cvScoreRaw = vector(1,ntips);
cvTotalScore=vector(1,cvNum);
cvTotalScoreRaw=vector(1,cvNum);
tipNodeList = pNewList();
TreeToTaxaPtrList(tree->root,tipNodeList); /* get a list of all the tip nodes */


if (method==LaF || method==NP)
	cvNum=1;	/* just overrides the default value of cvNum for these methods, forcing them to only go once */



numFixed=numFixedNodes(tree->root);
numNodes=numFixed;
if (numFixed <2 )
	{
	doGenericAlert("Must have at least two fixed nodes for fossil cross validation");
	return;
	}
fixedTime=(double *)myMalloc(numNodes*sizeof(double));
estTime=(double *)myMalloc(numNodes*sizeof(double));

nodeArray=(NODE *)myMalloc(numNodes*sizeof(NODE));
curIndex=0;
setupFixedNodeArray(tree->root, nodeArray, &curIndex);

verbose=gNexDataPtr->RateBlockParms.verbose;
if (verbose > 0)
	printf("Begin fossil cross-validation analyses...\n");
for (j=0;j<cvNum;j++)
	{
	numSuccess=0;
	numFail=0;
    cvSum=0.0;
	cvRawSum=0.0;
	smooth=pow(10.0,j*cvInc+cvStart);
	gNexDataPtr->RateBlockParms.smoothing=smooth;
	for (k=1;k<=numNodes;k++)
		{
		i=k-1;
		CVNode=nodeArray[i];
		
		CVNode->free=1;
		saveTime=CVNode->time;	
		gNexDataPtr->RateBlockParms.verbose=0; 
	/* suppress all output from the actual optimization run...may want to allow it for debugging though! */
		doObjFunc(tree,method,nRates,algorithm,&success);
		gNexDataPtr->RateBlockParms.verbose= verbose; /* restore output verbositude to current value */
		if (success)
			{
			estTime[i]=CVNode->time;
			cvResult[k]=1; /* good */
			++numSuccess;

			cvScore[k]=2*fabs(saveTime-estTime[i])/(saveTime+estTime[i]); 
			cvScoreRaw[k]=fabs(saveTime-estTime[i]);
			cvRawSum+=cvScoreRaw[k];
			cvSum+=cvScore[k];
			}
		else
			{
			cvResult[k]=0; /*failed */
			cvScore[k]=-99.9;
			cvScoreRaw[k]=0.0; /* should probably set this to some bogus value! */
			++numFail;
			}
		CVNode->free=0; // restore fixity and age
		CVNode->time=saveTime;
		}
		
	printf("\nFossil-fixed cross validation analysis\n\n");
	printf("Node\tTaxon\t\tAge\t\tEst.Age\t\tFract.Score\tRaw.Score\n");
	printf("---------------------------------------------------------------------------\n");
	for (i=0;i<numNodes;i++)
		{
		CVNode=nodeArray[i];
		printf("%i\t%.6s\t\t%6.2f\t\t%6.2f\t",i+1,CVNode->taxon_name,CVNode->time,estTime[i]);
		printf("\t%f\t%6.2f\n",cvScore[i+1],cvScoreRaw[i+1]);
		}

	cvTotalScore[j+1]=cvSum/(numNodes-numFail); /* on the off chance that some reps failed, don't count their contributions */
	cvTotalScoreRaw[j+1]=cvRawSum/(numNodes-numFail);
	
	if (numFail==0)
		cvResultFinal[j+1]=1; /* all prunings led to successful optimizations */
	else
		cvResultFinal[j+1]=0; /* some prunings had failed optimizations */

/* REMEMBER WE ARE OFTEN NOT COUNTING A TIP DESCENDED FROM THE ROOT */

	}
printf("********************************************************************************\n\n");
printf("Summary results of fossil-fixed cross validation analysis for tree %s\n",tree->name);
  switch (method)
	{
	case PENLIKE:printf("Method = Penalized Likelihood\n");break;
	case LaF:printf("Method = Langley and Fitch\n");break;
	case LFLOCAL:printf("Method = Langley and Fitch (with %i local rates)\n",nRates);break;
	case NP:printf("Method = Non-parametric\n");
	}
  switch (algorithm)
	{
	case POWELL: printf("Optimization via Powell's method\n");break;
	case QNEWT:  printf("Optimization via quasi-Newton method with analytical gradients\n");
	}
printf("\nFixed nodes:%i\n",numFixed);
printf("\nlog10\n");
printf("smooth\tsmooth\t\tFract Error\tRaw Error\n");
printf("--------------------------------------------------------------------------------\n");
for (j=0;j<cvNum;j++)
	{
	if (cvResultFinal[j+1])
		Result=Good;
	else
		Result=Failed;
	smooth=pow(10.0,j*cvInc+cvStart);
	printf("%6.2f\t%6.2g\t\t%6.4f\t\t%6.4f\t(%s)\n",j*cvInc+cvStart,smooth,cvTotalScore[j+1],cvTotalScoreRaw[j+1],Result);
	}
printf("********************************************************************************\n\n");
myFree(cvResult);
myFree(cvResultFinal);
myFree(fixedTime);
myFree(estTime);
myFree(nodeArray);
free_vector(cvScore,1,ntips);
free_vector(cvScoreRaw,1,ntips);
free_vector(cvTotalScore,1,cvNum);
free_vector(cvTotalScoreRaw,1,cvNum);
freepList(tipNodeList);


return;
}


/****************************************************************/

/* Following routines handle the custom 'Rate' block */

void doRateBlock(void)
{
  NODETYPE *root;
  char* TD, *TreeName;
  int ix;
  long numTrees,j;
extern int gisConstrained;
extern double gGamma_c;
extern double gGamma_b;
PtrList lnode;
TREE thisTree;
 

gTaxaList=NULL;	   /* initialize this somewhere else ? */
  

/* 
 * Following code sets up an exclusion array the FIRST time this block is called.
 * After that it is assumed to be there.  Therefore,  if you want to re-execute
 * a new matrix with different num of chars,  we'll have to set the pointer to 
 * the exclusion array to NULL!
 */

if ((gNexDataPtr->NChars > 0 ) && (gNexDataPtr->excArray == NULL))
 {
 gNexDataPtr->excArray=(int*)myMalloc(gNexDataPtr->NChars*sizeof(int)); 
 if(gNexDataPtr->excArray ==NULL)
    fatal("Allocation error in excArray");
 }
if (gNexDataPtr->excArray)
  for (ix=0;ix<gNexDataPtr->NChars;ix++)
    gNexDataPtr->excArray[ix]=1;



  do 				/* need to put in error checking in case no DIMENSIONS statement */
	    {
	     aTokenPtr=nextToken();
	     if (*aTokenPtr=='\0')
		{
		return;
		}
	     
	     else if   (isEqual(aTokenPtr, "QUIT") || 
			isEqual(aTokenPtr, "Q") ||
			isEqual(aTokenPtr, "BYE"))
			exit(1);
	     else if(isEqual(aTokenPtr, "COVARION"))
			doCovarionCommand();
	     else if(isEqual(aTokenPtr, "ANC"))
			doAncestralCommand();
	     else if(isEqual(aTokenPtr, "CO"))
			doContOptCommand();
	     else if(isEqual(aTokenPtr, "VCV"))
			doVCVCommand();
	     else if(isEqual(aTokenPtr, "RRLIKE"))
			doRRLikeTestCommand();
	     else if(isEqual(aTokenPtr, "LOCALMODEL"))
			doLocalModelCommand();
	     else if(isEqual(aTokenPtr, "BLFORMAT"))
			doBLFormatCommand();
	     else if(isEqual(aTokenPtr, "SHOWAGE") || isEqual(aTokenPtr, "SHOW") )
			doShowAgeCommand();
	     else if(isEqual(aTokenPtr, "UNFIXAGE"))
			doUnSetAgeCommand();
	     else if(isEqual(aTokenPtr, "FIXAGE"))
			doSetAgeCommand();
	     else if(isEqual(aTokenPtr, "PRUNE"))
			doPruneTaxonCommand();
	     else if(isEqual(aTokenPtr, "REROOT"))
			doReRootCommand();
	     else if(isEqual(aTokenPtr, "CLEARTREES"))
			doClearTreesCommand();
	     else if(isEqual(aTokenPtr, "CLUSTER_HISTOGRAM"))
			doClusterHistogramCommand();
	     else if(isEqual(aTokenPtr, "COLLAPSE"))
			doCollapseCommand();
	     else if(isEqual(aTokenPtr, "EXECUTE"))
			doExecuteCommand();
	     else if (isEqual(aTokenPtr,"PROFILE"))  
			doBranchProfileCommand();
	     else if (isEqual(aTokenPtr,"MRCA"))  
			doMRCACommand();
	     else if (isEqual(aTokenPtr,"CLADE_CHECK"))  
			doCladeCheckCommand();
	     else if (isEqual(aTokenPtr,"CONSTRAIN"))  
		doConstrain_TimeCommand();
	     else if (isEqual(aTokenPtr,"CALIBRATE"))  
		doCalibrateCommand();
	     else if (isEqual(aTokenPtr,"DESCRIBE") || isEqual(aTokenPtr,"DESC"))  
		doPrintCommand();
	     else if (isEqual(aTokenPtr,"TAXASET"))  
		doTaxaSetCommand();
	     else if (isEqual(aTokenPtr,"DISTANCE"))  
		doDistanceCommand();
	     else if (isEqual(aTokenPtr,"BS"))  
		doBSCommand();
	     else if (isEqual(aTokenPtr,"SET"))  
		doSetCommand();
	     else if (isEqual(aTokenPtr,"MRP"))
		doSuperCommand();  
	     else if (isEqual(aTokenPtr,"DIVTIME") || isEqual(aTokenPtr,"DIV"))
		doDivTimeCommand();  
	     else if (isEqual(aTokenPtr,"BD"))
		doBD();
	     else if (isEqual(aTokenPtr,"SIMULATE") )  
			doSimCommand();
	     else if (isEqual(aTokenPtr,"MULTIPLY_BRANCHLENGTH_BY"))
		doLengthMultiplyCommand();
	     else if (isEqual(aTokenPtr,"MULT"))
		doLengthMultiplyCommand();
	     else if (isEqual(aTokenPtr,"CONVERT_BRANCHLENGTH_TO_TIME"))
		{
	        if (gNexDataPtr->isTrees)
			{
			lnode=gNexDataPtr->inTrees;
			LISTLOOP (lnode)
				{
				thisTree=lnode->item;
				convert_branchlength_to_time(thisTree->root);
				}
			}
		}
	     else if (isEqual(aTokenPtr,"RR"))
		{  
	        if (gNexDataPtr->isTrees)
			{
			lnode=gNexDataPtr->inTrees;
			LISTLOOP (lnode)
				{
				thisTree=lnode->item;
				doWuLiCommand(thisTree->root);
				}
			}
		}
    }  while (!isEqual(aTokenPtr,"END") && !isEqual(aTokenPtr,"ENDBLOCK") && aTokenPtr);
    aTokenPtr=nextToken();
    if (!isEqual(aTokenPtr,";"))
	    doGenericAlert("Block not terminated with semicolon");
return;
}


/****************************************************************/
/**************** COMMAND PROCESSING FUNCTIONS ******************/
/****************************************************************/
int gNtips, gBDnumvar,gStemFlag;
double *gTimes,gRootDur;
#define POLYTOMY 1
#define SQR(x) (x)*(x)

void doBD(void) /* do birth-death stats on a specified clade */
    {
    extern NODETYPE* gRoot;
    long minN=1000000000,maxN=-1000000000;
    PtrList lnode;
    TREE thisTree;
    int N0,success, numIter, i, n, nint, ntrees, ix=1, profile_flag=0, dp_flag=0, Nee_flag=0, sumerr=0, ixm=0, numvar;
    double K1,K1var,K2,K2var,pure_birth_estimate, B, Optimum, a, r, AgeRoot, AgeFoundNode, CalAgeRoot, CalAgeFoundNode=-1.0, 
	    Kendall_var, Moran_var, mean1,adev1,sdev1,var1,skew1,curt1, mean2,adev2,sdev2,var2,skew2,curt2, 
	    Kendall_estimate, Raup_estimate, Like1, Like2, Like2Param,Like1Param,LR,SG,L1P,L2first,L2second,rate_ml;
    double params[3], *Y, *X,  *data1, *data2, *data3, *data4, *data5, *data6, *data7, *data8;
	double r1,r2,r0,LR_other;
    int root_id=1; /* as default use the root node */
    int stemFlag=0; /* default assume crown group */
	int SisterLRflag=0;
	int countLR=0,countSG=0,n1,n2,smaller;
    char *dummy,  *taxon;
    double rootAge;
    NODETYPE * found_node, *root, *first, *second;
    taxon=NULL;
    while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			{ /* must accept either a name(string) or taxon id (int) */
			if (isdigit(*LocalToken))
			    root_id=strtod(LocalToken,&dummy);
			else
			    taxon=DupStr(LocalToken);
			}
		if (parse_assignment2("AGE"))
			CalAgeFoundNode=strtod(LocalToken,&dummy);
		if (parse_assignment2("PROFILE"))
			if (isEqual(LocalToken,"YES"))
				profile_flag=1;
			else
				profile_flag=0;
		if (parse_assignment2("DIVPLOT"))
			if (isEqual(LocalToken,"YES"))
				dp_flag=1;
			else
				dp_flag=0;
		if (parse_assignment2("NEE"))
			if (isEqual(LocalToken,"YES"))
				Nee_flag=1;
			else
				Nee_flag=0;
		if (parse_assignment2("SISTER_LR"))
			if (isEqual(LocalToken,"YES"))
				SisterLRflag=1;
			else
				SisterLRflag=0;
		if (parse_assignment2("STEM"))
			if (isEqual(LocalToken,"YES"))
				stemFlag=1;
			else
				stemFlag=0;
		
		}

if (gNexDataPtr->isTrees)
		{
		lnode=gNexDataPtr->inTrees;
		ntrees=pLengthList(lnode);
		/*if(profile_flag)*/ /* always need some of these now */
		    {
		    data1=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data2=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data3=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data4=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data5=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data6=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data7=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data8=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    }
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;



    if (taxon)
	    found_node=find_taxon_name(root,taxon);
    else
	{ 
	if (root_id>1)
	    found_node=find_id(root, root_id);
	else if (root_id==1)
	    found_node=root;	/* the default setting; look at whole tree */
	}
    if(found_node)
	{
	root_id=found_node->id;
	gTimes=sort_node_time_array(found_node); /* zero-offset array gTimes */
	
	
n=numdesc(found_node)-1; /* See comments under sort_node_time_array() */
if (n>=2)
    if(gTimes[n-1]==gTimes[n-2])
	++sumerr;

	if(dp_flag)  /* Plot diversity over time.  WORKS ONLY FOR RECONSTRUCTED PROCESS (ASSUMES NO EXTINCTION)*/
	    {
	    Y=(double*)myMalloc(n*sizeof(double));
	    X=(double*)myMalloc(n*sizeof(double));
	    for (i=0;i<n;i++)
		{
		Y[i]=log(i+2); 
		X[n-1-i]=found_node->time-gTimes[i];
		}
	    dumbPlot(X,  Y, n);
	    myFree(X);
	    myFree(Y);
	    }

// ** Seems like the calibration stuff is a bit mucked up, check Kendal_est..

	AgeRoot=root->time;
	AgeFoundNode=found_node->time;
	if (CalAgeFoundNode == -1.0)
	    CalAgeFoundNode = found_node->time; /* if we didn't read a calibrated age for this, just set to internal value */
	CalAgeRoot=AgeRoot*CalAgeFoundNode/AgeFoundNode;
	B=get_sum_durations(found_node);
	gNtips=numdesc(found_node);
	if (stemFlag)
		N0=1;
	else
		N0=2;
	Kendall_estimate=((gNtips-N0 )/B)        /*/  CalAgeRoot   */;
	Kendall_var=SQR(Kendall_estimate)/(N0*(exp(Kendall_estimate *CalAgeFoundNode -1)));
	Moran_var=SQR(Kendall_estimate)/(gNtips-N0);
	Raup_estimate=(log(gNtips)-log(N0))/CalAgeFoundNode;
	if(profile_flag)
	    {
	    data1[ix]=Kendall_estimate;
	    data2[ix]=Raup_estimate;
	    data3[ix]=gNtips;
		if (gNtips<minN)minN=gNtips;
		if (gNtips>maxN)maxN=gNtips;
	    data4[ix]=B;
	    }

	if (!profile_flag)
	    {
	    printf("Age of root = %f\n", AgeRoot);
	    printf("ML estimate of lineage birth rate under Yule model using durations (node %i numtips=%i B=%f ):\n", 
			    root_id, gNtips, B);
	    printf("Whole Tree Root internal age = %f....calibrated age = %f\n", AgeRoot, CalAgeRoot);
	    printf("Subtree for BD root internal age = %f....calibrated age = %f\n", AgeFoundNode, CalAgeFoundNode);
	    printf("Kendall estimate of lineage birth rate under Yule model = %f (%f)\n",Kendall_estimate );
	    printf("Kendall's 1949 estimate of variance and std dev of rate estimate: %f\t%f\n",
		    Kendall_var, sqrt(Kendall_var) );
	    printf("Moran's 1951 estimate of variance and std dev of rate estimate: %f\t%f\n",
		    Moran_var, sqrt(Moran_var) );
	    printf("'Raup' estimate of lineage birth rate under Yule model (log(N)/t))= %f\n", Raup_estimate);
	    }
	if (Nee_flag) 
	    {
	    gBDnumvar=2; /* 1 = pure birth model (seems to agree with Kendall estimator in practice!), 2= birth-death model */
	    params[1]=2.0;
	    params[2]=0.0;
	    numIter=100;
	    Like2 = -MinND(thisTree,-1,POWELL,BD_Like,NULL,params,gBDnumvar,&numIter, 0.00001,0.0001, &success );
		// the -1 is a dummy argument, MinND doesn't need it
	    r=params[1]/CalAgeRoot;
	    a=params[2];
	    data6[ix]=r;
	    data7[ix]=a;
	    printf("ML estimate in the BD model:a=%f \tr=%f\tspec=%f\textinct=%f\n", a, r, r/(1-a), a*r/(1-a));
	    gBDnumvar=1;
	    Like1 = -MinND(thisTree,-1,POWELL,BD_Like,NULL,params,gBDnumvar,&numIter, 0.00001,0.0001, &success );
	    r=params[1]/CalAgeRoot;
	    data8[ix]=r;
	    printf("ML estimate in the B (Yule) model:r=%f\n", r);
	    LR=2*(Like2-Like1);
	    data5[ix]=LR;
	    printf("Like1=%f Like2=%f LR=%f\n", Like1, Like2, LR);

	/*    plotOpt(params,10,0.0,1.0,0.0,5.0, "a", "r");*/
	    }
	if (SisterLRflag) 
	    {
		
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
	    //gBDnumvar=1; /* 1 = pure birth model (seems to agree with Kendall estimator in practice!), 2= birth-death model */
	    //params[1]=1.0;
	    //params[2]=0.0;
	    //numIter=100;
		
		root=thisTree->root;
		first=root->firstdesc;
		second=first->sib;
	// two param model
		n1=gNtips=numdesc(first);
		//gStemFlag = 1;
		//gRoot=first;
	    //numIter=100;
	    //Like1 = -MinND(thisTree,-1,POWELL,BD_Like,NULL,params,gBDnumvar,&numIter, 0.00001,0.0001, &success );
	    //r=params[1]/CalAgeRoot;
		//data6[ix]=r;
	    //a=params[2];
	//    printf("ML estimate in the first sister:a=%f \tr=%f\tspec=%f\textinct=%f\n", a, r, r/(1-a), a*r/(1-a));
			// the -1 is a dummy argument, MinND doesn't need it
		n2=gNtips=numdesc(second);
		//gStemFlag = 1;
		//gRoot=second;
	    //numIter=100;
	    //Like2 = -MinND(thisTree,-1,POWELL,BD_Like,NULL,params,gBDnumvar,&numIter, 0.00001,0.0001, &success );
	    //r=params[1]/CalAgeRoot;
		//data7[ix]=r;
	    //a=params[2];
	    //data8[ix]=r;
	 //   printf("ML estimate in the second sister:a=%f \tr=%f\tspec=%f\textinct=%f\n", a, r, r/(1-a), a*r/(1-a));
		//Like2Param = Like1+Like2;
	// one param model
		//gNtips=numdesc(root);
		//gStemFlag = 0;
		//gRoot=root;
	    //numIter=100;
		//Like1Param = -MinND(thisTree,-1,POWELL,BD_Like,NULL,params,gBDnumvar,&numIter, 0.00001,0.0001, &success );
	    //r=params[1]/CalAgeRoot;
		//data8[ix]=r;
	    //a=params[2];
		L2first=YuleLike(first, 1, &rate_ml);
		data6[ix]=rate_ml;
		r1=rate_ml;
		printf("L2first=%f mlrate2param=%f\n",L2first,rate_ml);
		L2second=YuleLike(second, 1, &rate_ml);
		data7[ix]=rate_ml;
		r2=rate_ml;
		printf("L2second=%f mlrate2param=%f\n",L2second,rate_ml);
		L1P=YuleLike(root, 0, &rate_ml);
		data8[ix]=rate_ml;
		r0=rate_ml;
	//    printf("ML estimate in the whole tree:a=%f \tr=%f\tspec=%f\textinct=%f\n", a, r, r/(1-a), a*r/(1-a));

	   // LR=2*(Like2Param-Like1Param);
		LR = 2*(L2first+L2second-L1P);
		//LR_other =  2*(log(r1/r0)*(n1-1)+log(r2/r0)*(n2-1)) ; alternative sweet and correct formula for the LR test (but undefined if n = 2 for crown)
		//printf("*************** %f %f\n",LR,LR_other);
// LR /= 1.29; chi-square correction factor...
	    data5[ix]=LR;
		if (LR>3.841) ++countLR;
		if (n1<n2) smaller=n1; else smaller=n2;
		SG = 2.0*smaller/(n1+n2-1);
		if (n1==n2) SG=1; //boundary case
		if (SG<0.05) ++countSG;
		//printf("n1=%i n2=%i LR=%f SG=%f\n",n1,n2,LR,SG);
	    //printf("Like1=%f Like2=%f Like1Param=%f Like2Param=%f LR=%f\n", Like1,Like2, Like1Param, Like2Param, LR);

		printf("L1P=%f mlrate1param=%f\n",L1P,rate_ml);



	/*    plotOpt(params,10,0.0,1.0,0.0,5.0, "a", "r");*/
	    }
	    
	} /* end found_node */
	++ix;
			} /* end LISTLOOP */
		}
	 if (stemFlag)
		printf("(Stem group simulation:N0=1)\n");
	 else
		printf("(Crown group simulation:N0=2)\n");
    if(profile_flag)
	{
	printf("\n****************\nProfile Analysis of %i trees\n", ntrees);
	moment(data3,ntrees,&mean1,&adev1,&sdev1,&var1,&skew1,&curt1);
	printf("Mean clade size = %f (range=[%li,%li])\n", mean1,minN,maxN);
	moment(data4,ntrees,&mean1,&adev1,&sdev1,&var1,&skew1,&curt1);
	printf("Mean B = %f\n", mean1);
	moment(data1,ntrees,&mean1,&adev1,&sdev1,&var1,&skew1,&curt1); K1=mean1;K1var=var1;
	printf("Summary statistics on Yule(Kendall perfect information) estimator for %i trees\n");
	printf("Mean diversification rate = %f\n", mean1);
	printf("Variance and standard deviation of diversification rate = %f \t%f\n", var1, sdev1);
	moment(data2,ntrees,&mean2,&adev2,&sdev2,&var2,&skew2,&curt2);K2=mean2;K2var=var2;
	printf("Summary statistics on Yule(Raup minimal information) estimator for %i trees\n");
	printf("Mean diversification rate = %f\n", mean2);
	printf("Variance and standard deviation of diversification rate = %f \t%f\n", var2, sdev2);
	moment(data5,ntrees,&mean2,&adev2,&sdev2,&var2,&skew2,&curt2);
	printf("Summary statistics on LR Test %i trees\n");
	printf("Mean and variance = %f \t%f\n",mean2, var2);
	moment(data6,ntrees,&mean2,&adev2,&sdev2,&var2,&skew2,&curt2);
	printf("Mean and sdev of r in Yule model for clade 1 = %f \t%f\n",mean2, sdev2);
	moment(data7,ntrees,&mean2,&adev2,&sdev2,&var2,&skew2,&curt2);
	printf("Mean and sdev of r in Yule model for clade 2 %f \t%f\n",mean2, var2);
	moment(data8,ntrees,&mean2,&adev2,&sdev2,&var2,&skew2,&curt2);
	printf("Mean and sdev of r  in Yule model for whole clade %f \t%f\n",mean2, var2);

	if (SisterLRflag) 
		printf ("reject(LR|SG)?\t%i\t%i\n",countLR,countSG);
//	printf("stats\t%f\t%f\t%f\t%f\n",K1,K1var,K2,K2var);

	myFree(data1);
	myFree(data2);
	}
    if (taxon)
	    myFree(taxon);
    return;
    }
/****************************************************************/
void doPrintCommand(void)  /* Draws ASCII version of trees */
{
int 	likeFlag=0;
static int plotwidth=0;
int 	whichTrees=0;	/*0=input; 1=output*/
static int	treemode=0;	/*0=cladogram, etc.*/
int	numTrees,j, treeix=1;
static  int tree=0;
char *	dummy;
NODETYPE *root;
PtrList lnode;
TREE thisTree;
char*	TD, *TreeName;
char 	*clado="CLADOGRAM", *phylo="PHYLOGRAM",*chrono="CHRONOGRAM", *rato="RATOGRAM", 
	*td="TREE DESCRIPTION", *pd="PHYLO DESCRIPTION",  *rd=
	"RATO DESCRIPTION", *ni="NODE INFORMATION", *id="ID INFORMATION", *trace="TRACE", 
	*tracephy="TRACEPHY", *marg="MARGINALS DESCRIPTION", *modeString;

	while (!isEqual(aTokenPtr=nextToken(),";")) 
		{
		if (parse_assignment2("TREE"))
			tree=strtod(LocalToken,&dummy);
		if (parse_assignment2("PLOTWIDTH"))
			plotwidth=strtod(LocalToken,&dummy);
		if (parse_assignment2("DISPLAY_ANCESTOR"))
			if (isEqual(LocalToken,"YES"))
				gLabel=1;
			else
				gLabel=0;
		if (parse_assignment2("SOURCE"))
			{
			if (isEqual(LocalToken,"INPUT"))
				whichTrees=0;
			else
				whichTrees=1;
			}
		if (parse_assignment2("PLOT"))
			{
			if (isEqual(LocalToken,"CLADOGRAM"))
				treemode=0;
			if (isEqual(LocalToken,"PHYLOGRAM"))
				treemode=1;
			if (isEqual(LocalToken,"CHRONOGRAM"))
				treemode=2;
			if (isEqual(LocalToken,"NODE_INFO"))
				treemode=3;
			if (isEqual(LocalToken,"RATOGRAM"))
				treemode=4;
			if (isEqual(LocalToken,"CHRONO_DESCRIPTION"))
				treemode=5; /* also prints durations */
			if (isEqual(LocalToken,"TREE_DESCRIPTION"))
				treemode=5; /* also prints durations */
			if (isEqual(LocalToken,"PHYLO_DESCRIPTION"))
				treemode=6;
			if (isEqual(LocalToken,"RATO_DESCRIPTION"))
				treemode=7;
			if (isEqual(LocalToken,"ID_DESCRIPTION"))
				treemode=8;
			if (isEqual(LocalToken,"TRACE"))
				treemode=9;
			if (isEqual(LocalToken,"TRACEPHY"))
				treemode=10;
			if (isEqual(LocalToken,"MARG_DESCRIPTION"))
				treemode=11;
			}
		if (parse_assignment2("LIKE"))
			{
			if (isEqual(LocalToken,"YES"))
				likeFlag=1;
			else
				likeFlag=0;
			}
		}
/* now do it */

	switch (treemode) 
		{
		case 0: modeString=clado;break;
		case 1: modeString=phylo;break;
		case 2: modeString=chrono;break;
		case 3: modeString=ni;break;
		case 4: modeString=rato;break;
		case 5: modeString=td;break;
		case 6: modeString=pd;break;
		case 7: modeString=rd;break;
		case 8: modeString=id;break;
		case 9: modeString=trace;break;
		case 10: modeString=tracephy;break;
		case 11: modeString=marg;break;
		}
	if (whichTrees==0)
		{
		if (gNexDataPtr->isTrees)
			{
			lnode=gNexDataPtr->inTrees;
			LISTLOOP (lnode)
				{
				if ((tree==0) || (tree==treeix)) 
				    /* if no tree specified OR a specific tree specified */
				  {
				  thisTree=lnode->item;
				  printf("[Printing tree %i]\n", treeix);
				  printf("\n[%s of tree %s]\n",modeString, thisTree->name);
				  if (treemode==8)
				    {
				    printf("\ntree %s = ", thisTree->name);
				    make_parens(thisTree->root, 3); /* TD with ids for node names*/
				    printf(";\n");
				    }
				  if (treemode==5)
				    {
				    printf("\ntree %s = ", thisTree->name);
				    make_parens(thisTree->root, 1); /* TD of chronogram*/
				    printf(";\n");
				    }
				  if (treemode==6)
				    {
				    printf("\ntree %s = ", thisTree->name);
				    make_parens(thisTree->root, 0); /* TD of phylogram*/
				    printf(";\n");
				    }
				  if (treemode==11)
				    {
				    printf("\ntree %s = ", thisTree->name);
				    make_parens(thisTree->root, 4); /* marginals description string*/
				    printf(";\n");
				    }
				  if (treemode==7)
				    {
				 /*   set_est_rates(thisTree->root,thisTree->est_b,thisTree->est_c);*/
				    printf("\ntree %s = ", thisTree->name);
				    make_parens(thisTree->root, 2); /* TD of phylogram*/
				    printf(";\n");
				    }
				  if (treemode==3)
				    {
				    printtree(thisTree->root);
				    }
				  if (  (treemode==0) || 
					(treemode==1) || 
					(treemode==9) || 
					(treemode==10) || 
					(treemode==2) || 
					(treemode==4))
				    DrawTree(thisTree->root,treemode, plotwidth);
				  if (likeFlag)
					printLikes(thisTree->root);
				  }
				++treeix;
				}
			}
		else
			doGenericAlert("No input trees available\n");
		}
	/* else.....*/
	return;
}
/****************************************************************/
static void doBLFormatCommand(void)
{
PtrList lnode;
TREE thisTree;
char * dummy;
extern long gNumSites;
struct RBP * rbp;
int roundflag=1;

rbp=&(gNexDataPtr->RateBlockParms);
printf("Executing blformat command...\n");
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("ROUND"))
			{
			if (isEqual(LocalToken,"YES"))
				rbp->roundFlag=1;
			else
				rbp->roundFlag=0;
			}
		if (parse_assignment2("NSITES"))
			{
			rbp->numSites=strtod(LocalToken,&dummy);
			gNumSites=rbp->numSites;
			printf ("Number of sites in sequences set to %li\n",rbp->numSites);
			}
		if (parse_assignment2("ULTRAMETRIC"))
			{
			if (isEqual(LocalToken,"YES"))
				{
				rbp->clockFmt=1;
				}
			else
				rbp->clockFmt=0;
			}
		if (parse_assignment2("LENGTHS"))
			{
			if (isEqual(LocalToken,"TOTAL"))
				{
				rbp->lengthFmt=0;
				printf("Branch lengths assumed to be in units of raw numbers of substitutions\n");
				}
			else
			  if (isEqual(LocalToken,"PERSITE"))
				{
				rbp->lengthFmt=1;
				printf("Branch lengths assumed to be in units of numbers of substitutions per site\n");
				}			
			}

		}



	if (gNexDataPtr->isTrees )
		{
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			if (rbp->lengthFmt  == 1)  /* if lengths are per site, convert to total numbers of substs. */
				traverseMultiplyLength(thisTree->root, (double)rbp->numSites,rbp->roundFlag);
			if (rbp->lengthFmt  == 0 && rbp->roundFlag)  
				traverseMultiplyLength(thisTree->root, 1,rbp->roundFlag); // this just forces a rounding by default in case stupid user inputs such
			if (rbp->clockFmt == 1)	   /* if tree is ultrametric, calculate times on scale of [0,1] */
				{
//				rootToTips(thisTree->root,0.0); ...checks for ultrametricity... for debugging mostly
				convert_branchlength_to_time(thisTree->root);
				thisTree->timesAvailable=1;
				thisTree->method=USER;	/* save the fact that USER supplied times */
				}
				
			}
		}
	if (rbp->lengthFmt  == 0 && rbp->roundFlag)  
		printf("All branch lengths were rounded to the nearest integer\n");
	if (rbp->lengthFmt == 1)
		{
		printf("All branch lengths multipled by the %li sites in the sequence\n",rbp->numSites);
		if (rbp->roundFlag)
			printf("Branch lengths rounded to nearest integer\n(This may not be a good idea when using CALIBRATE on ultrametric user supplied input trees. Use ROUND=NO then)\n");
		else
			printf("Branch lengths not rounded on input\n(If doing DIVTIME, you should SET ROUND=YES, the default)\n");
		}
	if (rbp->clockFmt == 1)
		printf("Tree is assumed to be ultrametric (Use CALIBRATE command to scale it)\n");


	return;
}
/****************************************************************/
void doFormat(void)
{

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (isEqual(aTokenPtr,"INTERLEAVE"))
			gNexDataPtr->Intlvflag=1;
		if (parse_assignment2("MISSING"))
			gNexDataPtr->missingchar = *LocalToken;
		if (parse_assignment2("GAP"))
			gNexDataPtr->gapchar = *LocalToken;
		if (parse_assignment2("MATCHCHAR"))
			gNexDataPtr->matchchar = *LocalToken;
		}
/* Report */
	printf("Missing character = %c\n",gNexDataPtr->missingchar);
	printf("Gap character = %c\n",gNexDataPtr->gapchar);
	printf("Match character = %c\n",gNexDataPtr->matchchar);


	return;
}
/*----------------------------------------------------------------*/
void doSetCommand(void)
{
extern int gPowellTrace;
extern long gNumSites;
char  *dummy;
long aSeed;
struct RBP * rbp;
rbp=&(gNexDataPtr->RateBlockParms);
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("MINDURFACTOR"))
			rbp->minDurFactor=strtod(LocalToken,&dummy);
		if (parse_assignment2("MINRATEFACTOR"))
			rbp->minRateFactor=strtod(LocalToken,&dummy);
		if (parse_assignment2("NEIGHBORPENALTY"))
		    {
		    if (isEqual(LocalToken, "YES"))
			rbp->NeighborPenalty=1;
		    if (isEqual(LocalToken, "NO"))
			rbp->NeighborPenalty=0;
		    }
		if (parse_assignment2("PENALTY"))
		    {
		    if (isEqual(LocalToken, "ADD"))
			rbp->PenaltyType=0;
		    if (isEqual(LocalToken, "LOG"))
			rbp->PenaltyType=1;
		    }
		if (parse_assignment2("RATES"))
		    {
		    if (isEqual(LocalToken, "EQUAL"))
			rbp->RatesAreGamma=0;
		    else
		      if (isEqual(LocalToken, "GAMMA"))
			rbp->RatesAreGamma=1;
		    }
		if (parse_assignment2("ACTIVEEPSILON"))
			rbp->activeEpsilon=strtod(LocalToken,&dummy);
		if (parse_assignment2("SHAPE"))
			rbp->alpha=strtod(LocalToken,&dummy);
		if (parse_assignment2("BARRIERTOL"))
			rbp->barrierTol=strtod(LocalToken,&dummy);
		if (parse_assignment2("MAXITER"))
			rbp->maxIter=strtod(LocalToken,&dummy);
		if (parse_assignment2("MAXBARRIERITER"))
			rbp->maxBarrierIter=strtod(LocalToken,&dummy);
		if (parse_assignment2("INITBARRIERFACTOR"))
			rbp->initBarrierFactor=strtod(LocalToken,&dummy);
		if (parse_assignment2("BARRIERMULTIPLIER"))
			rbp->barrierMultiplier=strtod(LocalToken,&dummy);
		if (parse_assignment2("LINMINOFFSET"))
			rbp->linminOffset=strtod(LocalToken,&dummy);
		if (parse_assignment2("CONTRACTFACTOR"))
			rbp->contractFactor=strtod(LocalToken,&dummy);
		if (parse_assignment2("MAXCONTRACTITER"))
			rbp->maxContractIter=strtod(LocalToken,&dummy);
		if (parse_assignment2("FTOL"))
			rbp->ftol=strtod(LocalToken,&dummy);
		if (parse_assignment2("SHOWCONVERGENCE"))
		    {
		    if (isEqual(LocalToken, "YES"))
			gNexDataPtr->RateBlockParms.showConvergence=1;
		    else
			gNexDataPtr->RateBlockParms.showConvergence=0;
		    }
		if (parse_assignment2("SHOWGRADIENT"))
		    {
		    if (isEqual(LocalToken, "YES"))
			gNexDataPtr->RateBlockParms.showGradient=1;
		    else
			gNexDataPtr->RateBlockParms.showGradient=0;
		    }
		if (parse_assignment2("CHECKGRADIENT"))
		    {
		    if (isEqual(LocalToken, "YES"))
			gNexDataPtr->RateBlockParms.checkGradient=1;
		    else
			gNexDataPtr->RateBlockParms.checkGradient=0;
		    }
		if (parse_assignment2("CLAMPROOT"))
		    {
		    if (isEqual(LocalToken, "YES"))
			gNexDataPtr->RateBlockParms.clampRoot=1;
		    else
			gNexDataPtr->RateBlockParms.clampRoot=0;
		    }
		if (parse_assignment2("TRACE"))
		    {
		    if (isEqual(LocalToken, "YES") || isEqual(LocalToken, "ON"))
			gPowellTrace=1;
		    else
			gPowellTrace=0;
		    }
		if (parse_assignment2("VERBOSE"))
			gNexDataPtr->RateBlockParms.verbose=strtod(LocalToken,&dummy);
		if (parse_assignment2("LOCAL_FACTOR"))
			gNexDataPtr->RateBlockParms.local_factor=strtod(LocalToken,&dummy);
		if (parse_assignment2("PERTURB_FACTOR"))
			gNexDataPtr->RateBlockParms.perturb_factor=strtod(LocalToken,&dummy);
		if (parse_assignment2("NPEXP"))
			{
			gNexDataPtr->RateBlockParms.npexp=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("SMOOTHING"))
			{
			gNexDataPtr->RateBlockParms.smoothing=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("SITES"))
			{
	     		if (isEqual(LocalToken,"ALL"))  
				doSitesCommand(0);
	     		if (isEqual(LocalToken,"EXCLUDE12"))  
				doSitesCommand(1);
	     		if (isEqual(LocalToken,"EXCLUDE3"))  
				doSitesCommand(3);

			}
		if (parse_assignment2("SEED"))
			{
			aSeed=strtod(LocalToken,&dummy);
			srand(aSeed);
			}
		if (parse_assignment2("NUM_RESTARTS"))
			gNexDataPtr->RateBlockParms.num_restarts=strtod(LocalToken,&dummy);
		if (parse_assignment2("NUM_RATE_GUESSES"))
			gNexDataPtr->RateBlockParms.num_rate_guesses=strtod(LocalToken,&dummy);
		if (parse_assignment2("NUM_TIME_GUESSES"))
			gNexDataPtr->RateBlockParms.num_time_guesses=strtod(LocalToken,&dummy);
		}





return;
}

/*----------------------------------------------------------------*/
void doMatrixGeneral(void)

/* reads a NEXUS data matrix and stores in global dmatrix as a set of strings corresponding 
to the rows of the matrix; interleaved and whitespace is allowed */

#define isNewLine(c) (((c)=='\n') || ((c)=='\r'))		
{
	int 	
			taxon, 
			saveix=0, 
			firstTaxon=1, 
			lastTaxon, 
			nextTaxon;
						
	double 	z;
	StrListPtr DM, TL;	
	
    printf("\n\nWARNING: reading DATA block is under construction!\n");

	if (gNexDataPtr->NTaxa == 0 || gNexDataPtr->NChars ==0)
	    fatal("Must specify matrix dimensions prior to reading matrix");

	lastTaxon=gNexDataPtr->NTaxa;	/* get from global ntaxa */
	firstTaxon=1;
	taxon=1;
	lastTaxon=gNexDataPtr->NTaxa;	/*from global data */
	
/* initialize the global lists */

	DM=gNexDataPtr->DMList = newStrListN(gNexDataPtr->NTaxa);
	TL=gNexDataPtr->TaxaList=newStrList();



/* Loop through all tokens in the command */
	gNewLine=1;
	for (;;) 
	    {
	    aTokenPtr=nextToken();
	    if (isNewLine(*aTokenPtr))
		continue;
	    if (isEqual(aTokenPtr, ";"))
		break;
	    appendStrList(TL,aTokenPtr); /*Store the taxon label */
	    aTokenPtr=nextToken();
	    while(!isNewLine(*aTokenPtr))
		{
		catkthStr(DM, aTokenPtr, (long)taxon);
		aTokenPtr=nextToken() ; /* skip over possible new lines */
		}
	    printf("%s\t%s\n", getkthStr(TL,(long)(taxon)), getkthStr(DM,(long)(taxon)));
	    ++taxon;
	    }



	gNewLine=0;
	checkMatrix();
	gNexDataPtr->isTaxa=1;	/* Set flag showing that labels read */
	gNexDataPtr->isChars=1;
	return;
}

/*----------------------------------------------------------------*/
void doMatrix(void)

/* reads a NEXUS data matrix and stores in global dmatrix as a set of strings corresponding 
to the rows of the matrix; interleaved and whitespace is allowed, but taxa must
be in same order as in TAXA block; polymorphisms as sets are simply stored as in the originial 
matrix--not parsed in any way--watch out! */

{
	int 	
			taxon, 
			saveix=0, 
			firstTaxon=1, 
			lastTaxon, 
			nextTaxon;
						
	double 	z;
	StrListPtr DM;	
	

	lastTaxon=gNexDataPtr->NTaxa;	/* get from global ntaxa */
	firstTaxon=1;
	taxon=1;
	lastTaxon=gNexDataPtr->NTaxa;	/*from global data */
	
/* initialize the global list for DM */

	gNexDataPtr->DMList = newStrListN(gNexDataPtr->NTaxa);
	DM=gNexDataPtr->DMList;

/* Loop through all tokens in the command */


	for (;;)
		{

/* Get a token ...and test it for a variety of conditions */

		aTokenPtr=nextToken();
		
/* bail at end of command  */

		if (isEqual(aTokenPtr,";"))
			break;	
									
/* reset index when the first taxon name is encountered (e.g. in interleaved) and 
start at top of loop again with new token */ 

		if (isEqual(aTokenPtr,getkthStr(gNexDataPtr->TaxaList,(long)firstTaxon)))
			{
			taxon=1;	
			continue;
			}
			
/* if token is the NEXT taxon name then on next loop begin storing in next row */
/* ... otherwise store data on present row */
		
		if (taxon<lastTaxon)		
			{
			if (isEqual(aTokenPtr,getkthStr(gNexDataPtr->TaxaList,(long)(taxon+1)))) /* taxon+1 is the next taxon */

				++taxon;

			else

				catkthStr(DM, aTokenPtr, taxon);

			}

/*... but note that we have to be careful not to check past the end of the labels array */ 

		else			

			catkthStr(DM, aTokenPtr, taxon);

		
		}	/* end for */


	gNexDataPtr->isChars=1;	/* This flag is set if a matrix command is read; might cause
							problems if the matrix command is empty of a mstrix! */
	checkMatrix();
	return;
}

/*----------------------------------------------------------------*/
void doDimensions(void)
{
	char  *dummy;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NCHAR"))
			{
			gNexDataPtr->NChars=(int)strtod(LocalToken,&dummy);
			printf("Number of characters in matrix = %i\n",gNexDataPtr->NChars);		
			}
		if (parse_assignment2("NTAX"))
			{
			gNexDataPtr->NTaxa=(int)strtod(LocalToken,&dummy);
			printf("Number of taxa in matrix = %i\n",gNexDataPtr->NTaxa);		
			}
		}
	return;
}
/*----------------------------------------------------------------*/
void doCharDimensions(void)
{
	char  *dummy;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NCHAR"))
			{
			gNexDataPtr->NChars=(int)strtod(LocalToken,&dummy);
			}
		}
	printf("Number of characters in matrix = %i\n",gNexDataPtr->NChars);		
	return;
}
/*----------------------------------------------------------------*/
void doTaxDimensions(void)
{
	char  *dummy;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NTAX"))
			{
			gNexDataPtr->NTaxa=(int)strtod(LocalToken,&dummy);
			}
		}
	printf("Number of taxa in matrix = %i\n",gNexDataPtr->NTaxa);		
	return;
}
/*----------------------------------------------------------------*/

void doTaxLabels(void)
{
	int ix=0,len;

	gNexDataPtr->TaxaList=newStrList();

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		appendStrList(gNexDataPtr->TaxaList,aTokenPtr);
		}
	if (lengthList(gNexDataPtr->TaxaList)<gNexDataPtr->NTaxa) 
		fatal ("Too few taxon labels");
	else
		gNexDataPtr->isTaxa=1;	/* Set flag showing that labels read */
	return;
}


/**************************************************************/
void doTranslateCommand(void)

/* Currently the program WILL read trees with numbers as taxon names and just store those numbers
 * as strings.  All we have to do below is recurse and replace tip numbers with stuff from the list
 * below...so do it someday
 */

{
	gNexDataPtr->TransList=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	/* if its not a ';' it should be the number*/
		{
	
		if (  (!isdigit(*aTokenPtr)) && (!isEqual(aTokenPtr,","))  )
			appendStrList(gNexDataPtr->TransList,aTokenPtr); /* store the label */
		}
	gNexDataPtr->isTranslate=1;	/* set flag */
	printf("Trees stored WITH translation table\n");		
	return;
}

/****************************************************************/

void doTreeCommand(void)
{
		long size;
		extern int curTree;
		int flag=0;
		char *stemp, *tree_name, *TD, *p;
		PtrList aTreeList;

	if (gNexDataPtr->inTrees == NULL)  /* if this is the first tree */
		{
		gNexDataPtr->inTrees=pNewListAlt(sizeof(struct treetype));
		aTreeList=gNexDataPtr->inTrees;
		}
	else
		{
						/* if a later tree */
		aTreeList=pListAddNode(gNexDataPtr->inTrees,sizeof(struct treetype));
		if (aTreeList==NULL)fatal("Couldn't allocate tree list properly");
		}
/*printf("%li\n", (long)aTreeList);*/
	TD=makeEmptyString();
	tree_name=makeEmptyString();
	if ( isEqual(aTokenPtr=nextToken(),";")) return;

	if (isEqual(aTokenPtr,"*")) aTokenPtr=nextToken(); /* first token might be an
			asterisk; if so ignore and get the next token */

	concat(&tree_name,aTokenPtr);	/* this token should be the tree label */
	(void)appendStrList(gNexDataPtr->TDLabelList,tree_name);
	if (  !isEqual(aTokenPtr=nextToken(),"=")) 
			return; /* error; '=' should be next */

p=strchr(bufPtr,';');
if (p)
	{
	size=(long)(p-bufPtr);
/*printf(":%10.10s:%li:%c:\n",bufPtr,size,*p);*/
	*p='\0';
	TD=DupStr(bufPtr);
	(void)appendStrList(gNexDataPtr->TDList,TD);
	bufPtr=++p;
	}
#if 0
	while (!isEqual(aTokenPtr=nextToken(),";"))
		concat(&TD,aTokenPtr);	/* get the tree string */
	(void)appendStrList(gNexDataPtr->TDList,TD);
#endif			
	Tree_Initialize(aTreeList->item, TD, tree_name);
	printf("Reading tree %s\n",tree_name);
	++curTree;
	gNexDataPtr->NumTrees=curTree;
	gNexDataPtr->isTrees=1;
/*print_tree_list(gNexDataPtr->inTrees);*/
	myFree(tree_name);
	myFree(TD);

	return;						
}
/****************************************************************/

void doDistanceCommand(void)
{
	StrListPtr aTaxaList;


	if ((gNexDataPtr->isChars==0) || ( gNexDataPtr->isTaxa==0))
		return;	/* don't have the right data from NEXUS file, so bail */


	aTaxaList=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	
			appendStrList(aTaxaList,aTokenPtr); /* store the label */
	doDistance(aTaxaList);
	freeStrList(aTaxaList);
	return;

}
/****************************************************************/
void doLengthMultiplyCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD;
	int roundflag=1; /* force it to round to nearest integer for now */	
	double multiplier=0.0;

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		multiplier=strtod(aTokenPtr, &dummy);
		}
		printf("All branch lengths multipled by AND ROUNDED TO NEAREST INTEGER %f\n", multiplier);
	if ((gNexDataPtr->isTrees) && (multiplier >  0.0) )
		{

		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			traverseMultiplyLength(thisTree->root, multiplier,roundflag);
			}
		}
	return;						
}
/****************************************************************/
static void doClusterHistogramCommand(void)

/* Default mode is to to print a cumulative histogram of unrooted cluster sizes across
 * a set of trees.  Also bins this histogram into two larger bins: "shallow" clades
 * of size 2-3 and "deep" clusters of size >3.
 * 
 * A cluster size is the size of the smaller partition in any bipartition of the taxa
 * 
 * If option NORMALIZE=YES is invoked,  then the program expects a set of N model
 * trees interleaved with a set of N strict consenus trees made from consensing the model
 * trees and the estimated trees.  It divides the number of clades in the latter by the nuber
 * of clades in the former and spits this out.  In other words,  the list goes 
 *	model tree1
 *	consensus tree1
 *	model tree2
 *	consensus tree 2,  etc.
 * 
 * NB.  All trees are assumed to be of the same size! Unprdictable results otherwise.
 * NB.  The clade of ALL taxa is ignored.
 */

{
	float b0=0.0, b1=0.0;
	PtrList lnode;
	TREE thisTree, modelTree, strictTree;
	char * dummy, *taxon, *TD;
	
	double multiplier=0.0;
	long * histo = NULL,*histo1, *histo2, TSize, bin1[2], bin2[2];
	double * cumHisto=NULL;
	int arraySize,i,nTrees,normFlag=0, ix;

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NORMALIZE"))
			{
			if (isEqual(LocalToken,"YES"))
				normFlag=1;
			else
				normFlag=0;
			}
		}
	if (normFlag==0)
	    {
	    if (gNexDataPtr->isTrees )
		    {
		    nTrees=pLengthList(gNexDataPtr->inTrees);
		    lnode=gNexDataPtr->inTrees;
		    LISTLOOP (lnode)
			    {
			    thisTree=lnode->item;
			    if (histo == NULL)
				    {
				    TSize=numdesc(thisTree->root);
				    arraySize=floor(LOG2(TSize))-1;
				    if (arraySize>0)
					{
					histo=(long *)myMalloc(arraySize*sizeof(long));
					cumHisto=(double *)myMalloc(arraySize*sizeof(double));
					for (i=0;i<arraySize;i++)
					    cumHisto[i]=0.0;
					}
				    }
			    for (i=0;i<arraySize;i++)
				histo[i]=0; 
			    ClusterHistogram(thisTree->root,histo,TSize);
			    for (i=0;i<arraySize;i++)
				cumHisto[i]+=histo[i];
			    }
		    printf("[\nMean histogram of cluster sizes\nSize of tree=%li\n\n", (long)TSize);
		    printf("Number of trees=%i\n",nTrees);
		    for (i=0;i<arraySize-1;i++)
		      if (i==arraySize -2 ) 
			printf("%li-%li:\t\t%f\n", 
			    (long)pow(2, i+1), (long)pow(2, i+2), cumHisto[i]/nTrees);
		      else
			printf("%li-%li:\t\t%f\n", 
			    (long)pow(2, i+1), (long)pow(2, i+2)-1, cumHisto[i]/nTrees);
		    printf("\n]\n");
    
		    }
	    else
		    printf("No trees currently in memory\n");
	    }
	else if (normFlag==1)
	    {
	    if (gNexDataPtr->isTrees )
		    {
		    nTrees=pLengthList(gNexDataPtr->inTrees)/2; /* each list is half as long */
		    lnode=gNexDataPtr->inTrees;
		    thisTree=lnode->item;
		    TSize=numdesc(thisTree->root);
		    arraySize=floor(LOG2(TSize))-1;
		    if (arraySize>0)
			{
			histo1=(long *)myMalloc(arraySize*sizeof(long));
			histo2=(long *)myMalloc(arraySize*sizeof(long));
			cumHisto=(double *)myMalloc(arraySize*sizeof(double));
			for (i=0;i<arraySize;i++)
			    cumHisto[i]=0.0;
			}

		    for (ix=1;ix<=nTrees;ix++) 
			    {
			    modelTree=lnode->item;
			    strictTree=lnode->next->item;
			    bin1[0]=0;bin1[1]=0;
			    bin2[0]=0;bin2[1]=0;
			    for (i=0;i<arraySize;i++)
				{
				histo1[i]=0;
				histo2[i]=0;
				}
			    ClusterHistogram(modelTree->root,histo1,TSize);
			    ClusterHistogram(strictTree->root,histo2,TSize);
			    for (i=0;i<arraySize-1;i++)
				{
				if (!((histo1[i] == 0) && (histo2[i] == 0)))
					{ 
					cumHisto[i]+=(float)histo2[i]/histo1[i];
	    /* IMPORTANT --> */		if (i==0) /*  "shallow" bins */
					    {
					    bin1[0]+=histo1[i];
					    bin2[0]+=histo2[i];
					    }
					else  /*..remaining "not shallow" bins */
					    {
					    bin1[1]+=histo1[i];
					    bin2[1]+=histo2[i];
					    }
					}
				}
			    b0+=(float)bin2[0]/(bin1[0]);
			    b1+=(float)bin2[1]/(bin1[1]);
			    lnode=lnode->next->next;  /* skip two */
			    }
		    printf("[\nNormalized histogram of cluster sizes\nSize of tree=%li\n\n", (long)TSize);
		    for (i=0;i<arraySize-1;i++)
		      if (i==arraySize -2 ) 
			printf("%li-%li:\t\t%f\n", 
			    (long)pow(2, i+1), (long)pow(2, i+2), cumHisto[i]/nTrees);
		      else
			printf("%li-%li:\t\t%f\n", 
			    (long)pow(2, i+1), (long)pow(2, i+2)-1, cumHisto[i]/nTrees);
		    
		    printf("[\nSummary cluster sizes\n");
		    printf("Shallow (2-3):%f\nNot shallow (4..):%f\n",
				 b0/nTrees, b1/nTrees);
		    printf("Deep(%li-%li):%f\n",
			(long)pow(2, arraySize-1), (long)pow(2,  arraySize),cumHisto[arraySize-2]/nTrees);
		    printf("\n]\n");
    
		    }
	    else
		    printf("No trees currently in memory\n");
	    }
	return;						
}
/****************************************************************/
static void doClearTreesCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD;
	
	double multiplier=0.0;

	if (gNexDataPtr->isTrees )
		{
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			Tree_Destructor(thisTree);
			}
		freepList(gNexDataPtr->inTrees);
		gNexDataPtr->inTrees=NULL;
		gNexDataPtr->isTrees=0;
		printf("All trees cleared from memory\n");
		}
	else
		printf("No trees currently in memory\n");
	return;						
}
/****************************************************************/
static void doCollapseCommand(void)

/* Collapses any zero-length branches to polytomies.

	** NB.! A BUG EXISTS: If an internal node has a name, but COLLAPSE
	is run, the name may get overwritten!
*/

{
	PtrList lnode;
	TREE thisTree;
	int num;

	if (gNexDataPtr->isTrees )
		{
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			num=0;
			thisTree=lnode->item;
			while(any_zero_internal_branches(thisTree->root))
				{
				collapse_zero(thisTree->root);
			/* have to call this repeatedly because collapse
				only works on first zero branch; then the tree
				is different and a node is missing, so...*/
				++num;
				}
			printf("%i zero-length branches collapsed\n",num);
			thisTree->numBranches=numBranches(thisTree->root);
			}
		}
	else
		printf("No trees currently in memory\n");
	return;						
}
/****************************************************************/
void doBSCommand(void)
{
	char * dummy, *taxon, *TD;
	

	gNexDataPtr->RateBlockParms.isBS=1;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NREPS"))
			gNexDataPtr->RateBlockParms.NReps=strtod(LocalToken,&dummy);
		if (parse_assignment2("SEED"))
			gNexDataPtr->RateBlockParms.seed=strtod(LocalToken,&dummy);
		}

		return;						
}
/****************************************************************/
static void doSuperCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD;
	
	long numTrees, numTaxa,treeNum,charNum,taxNum;
	int j, ii=1, ll, mm, icur, numInt,method, nn, wtFlag=0;
	int maxClades=0, maxTaxa;
	NODETYPE *root;
/* Fix these fixed length arrays!! */

	char **matrix /* [MAXTAX][MAXCLADES] */;
	float *wtset;
	StrListPtr aTaxaList,firstTaxaList;
	


	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("METHOD"))
			{
			if (isEqual(LocalToken,"BAUM"))
				method=0;
			else
				method=1; /* PURVIS */
			}
		if (parse_assignment2("WEIGHTS"))
			{
			if (isEqual(LocalToken,"YES"))
				wtFlag=1;
			else
				wtFlag=0;
			}
		}
	if (gNexDataPtr->isTrees)
		{

		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode) /* sets up a list of all unique taxa names
				    in taxaList[1] */
			{
			thisTree=lnode->item;
			root=thisTree->root;
			aTaxaList=newStrList();
			TreeToTaxaList(root, aTaxaList);
			if (ii==1)
			    firstTaxaList=aTaxaList;
			maxClades+=numIntNodes(root)-1; /* allow one char (clade)
			    for each interior node in tree (less the root)*/
			if (ii>1)
				{
				glomStrLists(firstTaxaList, aTaxaList);
				freeStrList(aTaxaList);
				}
			++ii;
			}
		maxTaxa=lengthList(firstTaxaList);

		matrix=(char **)myMalloc(maxTaxa*sizeof(char *));
		wtset=(float *)myMalloc(maxClades*sizeof(float));
		for (ll=0;ll<maxTaxa;ll++)  /* initialize the ABC matrix */
		    {
		    matrix[ll]=(char *)myMalloc(maxClades*sizeof(char));
		    for (mm=0;mm<maxClades;mm++)
			matrix[ll][mm]='?';
		    }
		    
   
    
		if (method==0)  /* Do the Baum and Ragan supertree */
		    {
		    icur=0;   
		    gColumn=0;/* global must be set prior to following */
		    lnode=gNexDataPtr->inTrees;
		    LISTLOOP (lnode)
			    {
			    thisTree=lnode->item;
			    root=thisTree->root;
			    numInt=numIntNodes(root)-1;
			    aTaxaList=newStrList();
			    TreeToTaxaList(root, aTaxaList);
			    numTaxa=lengthList(aTaxaList);
			    for(ll=1;ll<=numTaxa;ll++)
				{
				taxon=getkthStr(aTaxaList, ll);
				mm=findMatchStr(firstTaxaList, taxon);
				if (mm)
				   for (j=0;j<numInt;j++)
				    matrix[mm-1][icur+j]='0';
				}
			    j++;
			    icur+=numInt;
			    freeStrList(aTaxaList);
			    ABCSuperTree(root, firstTaxaList, matrix, wtset);			    
			    }
			    
	
		    printf("[Baum,  Ragan Supertree]\n\n");
		    printNexus(maxTaxa, maxClades,firstTaxaList, matrix );
		    if (wtFlag)
			{
			printf("weights ");
			for (nn=0;nn<maxClades-1;nn++)
				{
				if ((nn>0)&& ((nn/10)==(nn/10.0)))
				    printf("\n");
				printf("%5.2f:%i, ", wtset[nn],nn+1);
				}
			printf("%5.2f:%i;\n", wtset[maxClades-1],maxClades);
			}
		    }
		
  
		if (method==1) /* Do the Purvis kind of supertree */
		    {
		    for (ll=0;ll<maxTaxa;ll++)  /* initialize the ABC matrix */
			for (mm=0;mm<maxClades;mm++)
			    matrix[ll][mm]='?';
		    gColumn=0;/* global must be set prior to following */
		    lnode=gNexDataPtr->inTrees;
		    LISTLOOP (lnode)
			    {
			    thisTree=lnode->item;
			    root=thisTree->root;
			    ABCSuperTreePurvis(root, firstTaxaList, matrix, wtset);			    
			    }
			    
	
		    printf("[Purvis Supertree]\n\n");
		    printNexus(maxTaxa, maxClades,firstTaxaList, matrix );
		    }
	    freeStrList(firstTaxaList);

	    treeNum=1;   
	    charNum=1;
	    lnode=gNexDataPtr->inTrees;
	    printf ("begin paup;\n");
	    LISTLOOP (lnode)
			    {
			    thisTree=lnode->item;
			    root=thisTree->root;
			    numInt=numIntNodes(root)-1;
			    numTaxa=numdesc(root);
			    if (numInt>0)
			   	printf("charset input%i = %i-%i;\n",treeNum++,charNum,charNum+numInt-1);
			    else
				printf("[skipping input%i (no clades)]\n",treeNum++);
			    charNum+=numInt;
			    }
	    printf("end;\n");
	    }

	return;						

}  
/****************************************************************/
void printNexus(int ntaxa, int nchars, 	
	StrListPtr taxaList,  char **matrix)
{
int ll, mm;
printf("#Nexus\n");
printf("Begin data;\ndimensions ntax=%i nchar=%i;\n", ntaxa, nchars);
printf("format symbols=\"01\" missing=?;\nmatrix\n");

		for (ll=0;ll<ntaxa;ll++)
		    {
		    printf("%s\t", getkthStr(taxaList, ll+1));
		    for (mm=0;mm<nchars;mm++)
			printf("%c", matrix[ll][mm]);
		    printf("\n");
		    }
printf(";\nend;\n");
return;    
}

/****************************************************************/

void doConstrain_TimeCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD;
	
	long numTrees;
	int j, removeAll=0,fixAll=0;
	int root_id=0;
	int flagMax=-1,flagMin=-1;
	double maxAge=1.0e20, minAge=0.0; /* standard defaults; same as newnode*/
	NODETYPE *root, *found_node;

	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("REMOVE"))
			if (isEqual(LocalToken,"ALL"))
				removeAll=1;
		
		if (parse_assignment2("TAXON"))
			{
			if (isdigit(*LocalToken))
			    root_id=strtod(LocalToken,&dummy);
			else
			    taxon=DupStr(LocalToken);
			}
		if (parse_assignment2("MAXAGE"))
			{
			if (isEqual(LocalToken,"NONE"))
				flagMax=0;
			else
				{
				maxAge=strtod(LocalToken,&dummy);
				flagMax=1;
				}
			}
		if (parse_assignment2("MINAGE"))
			{
			if (isEqual(LocalToken,"NONE"))
				flagMin=0;
			else
				{
				minAge=strtod(LocalToken,&dummy);
				flagMin=1;
				}
			}
		if (parse_assignment2("MAX_AGE"))
			{
			if (isEqual(LocalToken,"NONE"))
				flagMax=0;
			else
				{
				maxAge=strtod(LocalToken,&dummy);
				flagMax=1;
				}
			}
		if (parse_assignment2("MIN_AGE"))
			{
			if (isEqual(LocalToken,"NONE"))
				flagMin=0;
			else
				{
				minAge=strtod(LocalToken,&dummy);
				flagMin=1;
				}
			}
		}
	if (gNexDataPtr->isTrees)
		{

		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (removeAll)
				unSetConstraints(root); /* remove all constraints from tree */
			else
			 if (root_id>0 || taxon != NULL) 
			  {
			  if (root_id>0)
			    found_node=find_id(root, root_id);
			  else
			    found_node=find_taxon_name(root,taxon);
			  if(found_node)
			     {
			     if (isFree(found_node))
				{
				if (flagMin==1)	/* if flags remain at -1, then no action taken */
					{
					printf("Setting minimum age constraint for taxon %s to %f\n",found_node->taxon_name,minAge);
					found_node->nodeIsConstrainedMin=1;
					found_node->minAge=minAge;
					}
				if (flagMin==0)
					{
					printf("Removing any minimum age constraint for taxon %s\n",found_node->taxon_name);
					found_node->nodeIsConstrainedMin=0;
					}
				if (flagMax==1)
					{
					printf("Setting maximum age constraint for taxon %s to %f\n",found_node->taxon_name,maxAge);
					found_node->nodeIsConstrainedMax=1;
					found_node->maxAge=maxAge;
					}
				if (flagMax==0)
					{
					printf("Removing any maximum age constraint for taxon %s\n",found_node->taxon_name);
					found_node->nodeIsConstrainedMax=0;
					}
				}
			     }
			  else
			     doGenericAlert("Cannot assign a constraint to a node that is fixed\nUse UNFIXAGE on this node");
			  }
			}
		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
void doCalibrateCommand(void)

/* Spits out the ages of all nodes relative to a given age of a specific node. On exit,
	all nodes will have been rescaled according to the one given node's age */

{
	PtrList lnode;
	TREE thisTree;
	char * dummy,  *TD, *profile_taxon;
	

static  double calAge=1.0;
	char * taxon=NULL;
	long numTrees;
	int i, j, profileFlag=0, id=1, ix, profile_node_id=0, taxon_id=0, calflag=0;
	double time, scalefactor=1.0; /*default calibration makes it equivalent to no correction */
	double ave, adev, sdev, var, skew, curt;
	double *data, *data_one_node;
	NODETYPE *root, *found_node, *node;
	StrListPtr aTaxaList, txPtr;
	PtrList nodeList, nLptr, mrcaPtr;
	NODETYPE *mrca, *s;
	extern long gNumSites;

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			{ 
			calflag=1;
			taxon=DupStr(LocalToken);
			}
		if (parse_assignment2("AGE"))
			calAge=strtod(LocalToken,&dummy);
		if (parse_assignment2("PROFILE_NODE"))
			{
			profile_taxon=DupStr(LocalToken);
			profileFlag=1;
			}
		
		}
	
	/*..............do the work...........*/	
		
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees and just spits out 
		 *calibrated ages.....
		 */
		ix=1;
		numTrees=pLengthList(gNexDataPtr->inTrees);
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			printf("Tree %i\n", ix);
			thisTree=lnode->item;
			if (thisTree->method != USER)
				doGenericAlert("WARNING: Calibrate command is designed for user-supplied ultrametric trees only!");
			root=thisTree->root;
			if (thisTree->timesAvailable)
			  {
			  if (taxon)
			    {
			    found_node=find_taxon_name(root,taxon);
			    printf("\nCalibrated Ages based on taxon %s @ age %f\n", 
				taxon, calAge);
			    }
			  else
			    {
			    found_node=root;
			    printf("\nCalibrated Ages based on ROOT age @ age %f\n", 
				calAge);
			    }
			  time=found_node->time;
			  scaleTree(root,calAge,found_node);
/*			  print_ages(root, time, calAge,thisTree->method);*/
			  }
			else
			  doGenericAlert("Times unavailable");
			++ix;
			}
		/*....works on a profile of identical topology trees and does
		/* summary stats by node...
		*/	
			
		if(profileFlag)
		    {
			data=(double*)myMalloc((numTrees+1)*sizeof(double));
			ix=1; /*init the index of trees */
			lnode=gNexDataPtr->inTrees;
			LISTLOOP (lnode)
				{
				thisTree=lnode->item;
				if (thisTree->timesAvailable)
				  {
				root=thisTree->root;
				/*...first get the scalefactor if calibration is internal*/
				found_node=find_taxon_name(root,taxon); 
					/* taxon is the calibrated taxon */
				if(found_node)
				    {
				    scalefactor=calAge/found_node->time;
				    }
				/*...now get node corresponding to id...*/
				else
				    printf("Couldn't find a calibration:Using root=1.0");
				
				node=find_taxon_name(root,profile_taxon); 
				if(node)
					{
					data[ix]=node->time*scalefactor;
					++ix;
// printf("Time of node for tree %i = %f\n",ix,data[ix]);
					}
				  }
				else
				  doGenericAlert("Times unavailable");
				}
			moment(data, numTrees, &ave, &adev,&sdev,
				&var, &skew, &curt); /* remember a 1-offset array */
			printf("Profile information for node across stored trees:\n"); 
			printf("Node=%s : Num trees=%i Mean time=%f  Standard deviation=%f Skewness=%f\n", 
					profile_taxon, numTrees,ave, sdev, skew);
			myFree(data);
		    }


		    


		}
	if (taxon)
		myFree(taxon);
	return;						
}
static void doShowAgeCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy,  *TD, *profile_taxon;
	

static  double calAge=1.0;
static  long iTree;
	char * taxon=NULL;
	long numTrees;
	int i, j, profileFlag=0, id=1, ix, profile_node_id=0, taxon_id=0, calflag=0,showNamed=0;
	double time, scalefactor=1.0; /*default calibration makes it equivalent to no correction */
	double ave, adev, sdev, var, skew, curt;
	double *data, *data_one_node;
	NODETYPE *root, *found_node, *node;
	StrListPtr aTaxaList, txPtr;
	PtrList nodeList, nLptr, mrcaPtr;
	NODETYPE *mrca, *s;

	
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TREE"))
			{
			iTree=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("SHOWNAMED")) // prints out a brief list of ages of named internal nodes only
			{
			if (isEqual(LocalToken,"YES"))
				showNamed=1;
			else
				showNamed=0;
			}
		}

	/*..............do the work...........*/	
		
	if (!gNexDataPtr->isTrees)
		printf("No input trees available\n");
	else
		{

		/*....works on any collection of trees and just spits out 
		 *calibrated ages.....
		 */
		ix=1;
		numTrees=pLengthList(gNexDataPtr->inTrees);
		lnode=gNexDataPtr->inTrees;
	        if (iTree>0) /* a specific tree was indicated */
			{
			if (iTree > pLengthList(lnode))
				{
				doGenericAlert("Invalid tree specified");
				return;
				}
			else
				{
				thisTree=(pListgetkthNode(lnode,iTree))->item;
				if (thisTree->timesAvailable)
					print_ages(thisTree->root,1.0,1.0,thisTree->method);
				}
			}
	        else
			LISTLOOP (lnode)
				{
				thisTree=lnode->item;
				root=thisTree->root;
				if (thisTree->timesAvailable)
				    if (showNamed)
					{
					printf("-----------------------------------------------------------------------------------------\n");
					printf("\nAges of internal named nodes only:\n");
					print_named_ages(root);
					}
				    else // the usual setting...
					{
					printf("-----------------------------------------------------------------------------------------\n");
					printf("Estimated ages and substitution rates for tree %s\n\n",thisTree->name);
					switch(thisTree->method)
						{
						case USER:printf("Reconstruction method: User-supplied ultrametric tree\n");break;
						case LaF:printf("Reconstruction method: Langley-Fitch (clock)\n");break;
						case NP:printf("Reconstruction method: NPRS\n");break;
						case PENLIKE:printf("Reconstruction method: Penalized likelihood\n");break;
						}
					printf("Named internal nodes indicated by [*]\n");
					printf("Rates are for branches subtending indicated node\n");
					printf("Rates are in units of substitutions per site per unit time\n\n");
					printf("\t\t\t\t  Constraints\t\t\t\tRates\n");
					printf("\tNode\t   Fix [Mod]\t  Min\t  Max\t  Age\t\tEstimated\tLocal\n");
					printf("-----------------------------------------------------------------------------------------\n");

					print_ages(root, 1.0,1.0,thisTree->method); /* use this more complex function to do something simple here! */
					printf("-----------------------------------------------------------------------------------------\n");

					summarize_rates(thisTree);

					}
				else
					doGenericAlert("Times and rates unavailable");
				++ix;
				}





		/*....works on a profile of identical topology trees and does
		/* summary stats by node...
		*/	
			
		if(profileFlag) /* currently never invokes...leave for use later perhaps */
		    {
			data=(double*)myMalloc((numTrees+1)*sizeof(double));
			ix=1; /*init the index of trees */
			lnode=gNexDataPtr->inTrees;
			LISTLOOP (lnode)
				{
				thisTree=lnode->item;
				root=thisTree->root;
				/*...first get the scalefactor if calibration is internal*/
				found_node=find_taxon_name(root,taxon); 
					/* taxon is the calibrated taxon */
				if(found_node)
				    {
				    scalefactor=calAge/found_node->time;
				    }
				/*...now get node corresponding to id...*/
				else
				    printf("Couldn't find a calibration:Using root=1.0");
				
				node=find_taxon_name(root,profile_taxon); 
				if(node)
					{
					data[ix]=node->time*scalefactor;
					++ix;
					}
				}
			moment(data, numTrees, &ave, &adev,&sdev,
				&var, &skew, &curt); /* remember a 1-offset array */
			printf("Profile information for node across stored trees:\n"); 
			printf("Node=%s : Num trees=%i Mean time=%f  Standard deviation=%f Skewness=%f\n", 
					profile_taxon, numTrees,ave, sdev, skew);
			myFree(data);
		    }


		    


		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
static void doRRLikeTestCommand(void)

/* 

Do a likelihood ratio test on the N clades descended from one node. Use the localmodel feature.

Comments: Interesting issues here, because we might want to use time constraint information and we might not.
When the focal node is the root of a clade without constraints at all, the DIVTIME routine will bail because
it will think the root cannot be estimated. In that case, it seems reasonable to set the age of the (local) root
node temporarily to 1.0. Alternatively for each node we might make a time-informed and a time-uninformed relative
rate test. In the latter, we ignore all time information of descendant nodes. In the former we take info into account,
and if the problem is inestimable, let the chips fall where they may! 
*/

{
	PtrList lnode;
	TREE thisTree,subtree;
	char * dummy, *taxon, *TD;
 	double Like0,Like1,LR;	
	long numTrees;
	int constrain,warn;
	int i, j, id=1, ix, stemFlag=1, model=0; //init. stemFlag to include stem lineage of each child
	int success=0,method,algorithm,nRates,allFlag=0;
	NODETYPE *root, *child, *node,*found_node,*saveNode;

	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees  */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root; // eventually do this for an arbitrary node!
			if (taxon)
			    {
				if (isEqual(taxon,"ALL"))
					{
					allFlag=1; // implement later
					}
				else
					{
			   		if (isEqual(taxon,"ROOT"))
						found_node=root;
			    		else
			    			found_node=find_taxon_name(root,taxon);
					}
			    	if (found_node)
					{
					warn=warnEstRoot(found_node);
					if (warn == 1)  // only when the divtime would normally bail do we consider this an unconstr search
						{
						constrain=0;	
						found_node->free=0;
						found_node->time=100.0; // unfixed so fix it at some arbitray age
						printf("...Insufficient time constraints present in RRLike...fixing age of %s at 100.0\n",found_node->taxon_name);
						}
					else
						constrain=1;
					saveNode=found_node->anc;
					subtree=Subtree_Initialize(thisTree,found_node); // pull off this subtree and work on it only
					algorithm=TN;
					method=LaF;
					nRates=1;
					doObjFunc(subtree,method,nRates,algorithm,&success);
					Like0=subtree->obj;

					child=found_node->firstdesc;
					SIBLOOP(child)
						{
						setLocalModel(child,model,stemFlag);
						++model;
						}
					nRates=model;
					algorithm=POWELL; // note we have to switch to POWELL for the Local clock methods!
					method=LFLOCAL;
					doObjFunc(subtree,method,nRates,algorithm,&success);
					Like1=subtree->obj;
					LR = -2 * (Like0-Like1);
					printf("\nRELATIVE RATE TEST USING LIKELIHOOD RATIO\n\n");
					printf("    Node       Clock   Non-clock     LR Stat    d.f.    Constrained?\n");
					printf("--------------------------------------------------------------------\n");
					printf("%8s    %8.2f    %8.2f    %8.2f    %4i        %4i\n",found_node->taxon_name,Like0,Like1,LR,nRates-1,constrain);
					found_node->anc=saveNode; // restore the subtree's link to original tree
					if (warn == 1)  
						{
						found_node->free=1;
						found_node->time=0.0; // reset values 
						}
					}
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
static void doLocalModelCommand(void)

/* 

Takes instructions on setting up a local clock model with finite number of rate parameters distributed across tree.
Parameters are indexed from 0..N-1, if there are N rates. Assigns all branches in clade 'Taxon' some index value. If
Taxon is a tip, then its subtending branch is used. If 'stem' is set to yes, then the stem branch is assigned as well
as all members of the clade.

 */

{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD;
	
	long numTrees;
	int i, j, id=1, ix, stemFlag=0, model;
	NODETYPE *root, *found_node, *node;

	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		if (parse_assignment2("STEM"))
			{
			if (isEqual(LocalToken,"YES"))
				stemFlag=1;
			else
				stemFlag=0;
			}
		if (parse_assignment2("RATEINDEX"))
			model=strtod(LocalToken,&dummy);

			if (isEqual(LocalToken,"ALL"))
				(void)preOrder(root,fixNodeAge); 
					/* force all nodes to have their ages fixed (hopefully ages are set somehow!) */
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees  */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (taxon)
			    {
			   	if (isEqual(taxon,"ROOT"))
					found_node=root;
			    	else
			    		found_node=find_taxon_name(root,taxon);
			    	if (found_node)
					{
					setLocalModel(found_node,model,stemFlag);
					}
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
void doSetAgeCommand(void)

/* Sets the age of any node in the tree, but this is transient if node is internal */
/* DO NOT PERMIT THE USE OF this command (now 'FixAge') on ultrametric trees...collides
   with 'Calibrate' command ...probably should just have one command, and prevent user
   from fixage on more than one node for ultrametric trees LATER...*/

{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD, *profile_taxon;
	
	long numTrees;
	int i, j, profileFlag=0, id=1, ix, profile_node_id=0, taxon_id=0, calflag=0;
	double age, time, scalefactor=1.0; /*default calibration makes it equivalent to no correction */
	double *data, *data_one_node;
	NODETYPE *root, *found_node, *node;
	StrListPtr aTaxaList, txPtr;
	PtrList nodeList, nLptr, mrcaPtr;
	NODETYPE *mrca, *s;
	struct RBP * rbp;
	rbp=&(gNexDataPtr->RateBlockParms);

	taxon=NULL;

	if (rbp->clockFmt)
		{
		doGenericAlert("Can't SETAGE on ultrametric trees; use CALIBRATE command instead!");
		return;	
		}

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		if (parse_assignment2("AGE"))
			age=strtod(LocalToken,&dummy);
		if (parse_assignment2("FIX"))
			if (isEqual(LocalToken,"ALL"))
				(void)preOrder(root,fixNodeAge); 
					/* force all nodes to have their ages fixed (hopefully ages are set somehow!) */
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees  */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (taxon)
			    {
			    if (isEqual(taxon,"ALL"))
				(void)preOrder(root,fixNodeAge); 
				/* If only say 'setAge taxon=all' then fix all times to whatever they are. 
				   Doesn't permit us to set all nodes to some one age--that'd be dumb*/
			    else
				{
			   	if (isEqual(taxon,"ROOT"))
					found_node=root;
			    	else
			    		found_node=find_taxon_name(root,taxon);
			    	if (found_node)
					{
					found_node->free=0;
			    		found_node->time=age;
					found_node->nodeIsConstrainedMax=0;
					found_node->nodeIsConstrainedMin=0; /* overrides preexisting min or max constraints ! */
					printf("Fixing age of %s at %f\n",taxon,age);
					printf(" (The age of this node will no longer be estimated.)\n");
					printf(" (This command overides any previous age constraints for this node.)\n");
					printf(" (The total number of fixed ages is now %i)\n",numFixedNodes(root));
					}
				}
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
static void doUnSetAgeCommand(void)

/* Frees the age of a node, which will subsequently be estimated */

{
	PtrList lnode;
	TREE thisTree;
	char *taxon;
	NODETYPE *root, *found_node, *node;

	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees  */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (taxon)
			    {
			    if (isEqual(taxon,"ALL"))
				(void)preOrder(root,unFixNodeAge);
			    else
				{
			   	 if (isEqual(taxon,"ROOT"))
					found_node=root;
			   	 else
			    		found_node=find_taxon_name(root,taxon);
			   	 if (found_node)
					{
					found_node->free=1;
					printf("Unfixing age of %s\n",taxon);
					printf(" (the age of this node WILL be estimated in future searches)\n");
					}
				}
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
static void doCovarionCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	int success, maxIterations,nParm, showProbs=0,showChanges=0;
	double ftol=0.0001, linMinDelta=0.05;
	char *dummy;
	struct RBP * rbp;
	gNexDataPtr->RateBlockParms.estCov=0;
	gNexDataPtr->RateBlockParms.cov_brlens=0;
	rbp=&(gNexDataPtr->RateBlockParms);
	maxIterations=rbp->maxIter;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("ESTIMATE"))
		    {
		    if (isEqual(LocalToken,"YES"))
				gNexDataPtr->RateBlockParms.estCov=1;
		    else
				gNexDataPtr->RateBlockParms.estCov=0;

		    }
		if (parse_assignment2("SHOWPROBS"))
		    {
		    if (isEqual(LocalToken,"YES"))
				showProbs=1;
		    else
				showProbs=0;
		    }
		if (parse_assignment2("SHOWCHANGES"))
		    {
		    if (isEqual(LocalToken,"YES"))
				showChanges=1;
		    else
				showChanges=0;
		    }
		if (parse_assignment2("BRLENS"))
		    {
		    if (isEqual(LocalToken,"ONE"))
				gNexDataPtr->RateBlockParms.cov_brlens=1;
		    if (isEqual(LocalToken,"USER"))
				gNexDataPtr->RateBlockParms.cov_brlens=0;

		    }
		if (parse_assignment2("S_RATE"))
			gNexDataPtr->RateBlockParms.s_rate=strtod(LocalToken,&dummy);
		if (parse_assignment2("R_RATE"))
			{
			gNexDataPtr->RateBlockParms.r_rate=strtod(LocalToken,&dummy);
			}
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			printf("...doing 3-state optimization...\n");
			printf("Optimization parameters:\n  ftol...%f\n");
			covarionOptimize(thisTree,&maxIterations, rbp->ftol,rbp->linminOffset,&success );
			if (showProbs) printCovarion(thisTree->root);
			if (showChanges) printChanges(thisTree->root);
			}

		}
	return;		
}				
/****************************************************************/
static void doAncestralCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	int success, maxIterations,nParm;
	double ftol=0.0001, linMinDelta=0.05;
	char *dummy;
	struct RBP * rbp;
	rbp=&(gNexDataPtr->RateBlockParms);
	maxIterations=rbp->maxIter;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			printf("...doing ancestral state squared change parsimony optimization...\n");
			printf("Optimization parameters:\n  ftol...%f\n");
			ancestralOptimize(thisTree,&maxIterations, rbp->ftol,rbp->linminOffset,&success );
			printf("Node\t\tState\t\tAnc State\tDiff\tSign of Difference\n");
			printAncestral(thisTree->root);
			}

		}
	return;		
}				
/****************************************************************/
static void doContOptCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	int nRates=1, success, maxIterations,nParm;
	double ftol=0.0001, linMinDelta=0.05;
	char *dummy;
	maxIterations=1000;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NRATES"))
			nRates=strtod(LocalToken,&dummy);
		}
	nParm=nRates+1;
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			printf("...doing continuous character rate optimization...\n");
			contOptimize(thisTree,nParm,&maxIterations, ftol,linMinDelta,&success );
			}

		}
	return;						
}

/****************************************************************/
void doPruneTaxonCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD, *profile_taxon;
	
	long numTrees;
	int i, j, id=1, ix,  taxon_id=0;
	NODETYPE *root, *found_node, *node;
	StrListPtr aTaxaList, txPtr;
	PtrList nodeList, nLptr, mrcaPtr;
	NODETYPE *mrca, *s;

	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (taxon)
			    {
			    found_node=find_taxon_name(root,taxon);
			    if (found_node)
				{
			    	thisTree->root=RemoveTaxon(thisTree,found_node);
				thisTree->numBranches=numBranches(thisTree->root);
				printf("Pruning taxon %s\n",taxon);
				}
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
static void doVCVCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD, *profile_taxon;
	double T;
	long i, j;
	NODETYPE *root, *found_node, *node;
	NODE a,b,c;
	PtrList nodeList;
	long lengthList;
	
	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (taxon)
			    {
			    found_node=find_taxon_name(root,taxon);
			    if (found_node)
					{
					nodeList=pNewList();
					TreeToTaxaPtrList(root,nodeList);
					lengthList=pLengthList(nodeList);
					for (i=1;i<=lengthList;i++)
						{
						a=(NODE)(pListgetkthNode(nodeList, i)->item);
						printf("%s\t",a->taxon_name);
						for (j=1;j<=lengthList;j++)
							{
							b=(NODE)(pListgetkthNode(nodeList, j)->item);
							c=mrca(a,b);
							T=pathLengthTime(c,found_node);
							printf("%f\t",T);
							}
						printf("\n");
						}
					freepList(nodeList);
					}
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
void doReRootCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD, *profile_taxon;
	
	long numTrees;
	int i, j, id=1, ix,  taxon_id=0, atNode=0;
	NODETYPE *root, *found_node, *node;
	StrListPtr aTaxaList, txPtr;
	PtrList nodeList, nLptr, mrcaPtr;
	NODETYPE *mrca, *s;

	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		if (parse_assignment2("ATNODE"))  // reroots at a node rather than maintaining a binary root
			if (isEqual(LocalToken,"YES"))
				atNode=1;
			else
				atNode=0;
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (taxon)
			    {
			    found_node=find_taxon_name(root,taxon);
			    if (atNode==0)
			    	thisTree->root=ReRoot(found_node);
			    if (atNode==1)
			    	thisTree->root=ReRoot2(found_node);
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
static void doBranchProfileCommand(void)
{
	extern long gNumSites;
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD, *profile_taxon;
	
	long numTrees;
	int i, j, profileFlag=0, id=1, ix, profile_node_id=0, taxon_id=0, calflag=0,parmFlag=2;
	double calAge=1.0, time, scalefactor=1.0; /*default calibration makes it equivalent to no correction */
	double ave, adev, sdev, var, skew, curt,min=1e20,max=-1e20;
	double *data, *data_one_node;
	NODETYPE *root, *found_node, *node;
	StrListPtr aTaxaList, txPtr;
	PtrList nodeList, nLptr, mrcaPtr;
	NODETYPE *mrca, *s;

	taxon=NULL;
	profile_taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			{
			profile_taxon=DupStr(LocalToken);
			}
		if (parse_assignment2("PARAMETER"))
			{
			if (isEqual(LocalToken,"LENGTH"))
				parmFlag=1;
			else if (isEqual(LocalToken,"AGE"))
				parmFlag=2;
			else if (isEqual(LocalToken,"RATE"))
				parmFlag=3;

			}
		
		}
	
	/*..............do the work...........*/	
		
	if (gNexDataPtr->isTrees)
		{
		numTrees=pLengthList(gNexDataPtr->inTrees);
		data=(double*)myMalloc((numTrees+1)*sizeof(double));
		ix=0; /*init the index of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			node=find_taxon_name(root,profile_taxon); 
			if(node)
				{
				++ix;
				switch (parmFlag)
					{
					case 1: data[ix]=node->length;break;
					case 2: data[ix]=node->time;break;
					case 3: data[ix]=node->estRate/gNumSites;break;
					}
			/*	printf("%f\n",data[ix]);*/
				}
			else
				{
				printf("WARNING! Profiled node not found on tree (may have been collapsed?)\n");
				}
			}
		for (i=1;i<=ix;i++) // remember ix may be less than numTrees if some nodes are collapsed and not found
			{
			if (data[i]>max)max=data[i];
			if (data[i]<min)min=data[i];
			}
		printf("Profile information for node %s across %i tree(s) out of %i trees:\n",profile_taxon,ix,numTrees); 
		if (ix >=2)
			{
			moment(data, ix, &ave, &adev,&sdev,
				&var, &skew, &curt); /* remember a 1-offset array */
			printf("Mean = %f  Std dev = %f Min = %f Max = %f\n", ave, sdev,min,max);
			}
		else
			printf("Profile cannot be obtained when number of trees with given node < 2\n");
		myFree(data);
		}
	if(profile_taxon)
		myFree(profile_taxon);
	return;						
}
/****************************************************************/
void doSimCommand(void)
{
float T_LF[25][25];
float T_PL[25][25];
NODETYPE *nodea, *nodeb, *nodec, *nodeInt;
float bx,cx,lb,lc,trueAge,rangeFactor;
int ksteps,ib,ic;  // ksteps an even number please
	PtrList lnode;
	TREE thisTree;
	int verbose = 0, resetSeed=0;
	float bestSmooth;
	extern NODETYPE * gRoot;

	char * 	dummy;
	
	double			*RepCount,*RepMean,*RepDominance,*RepFreq1Class,*RepFractMonophyletic,*RepMonophyleticSpecies;
	double  *X, *Y, * time1, *time2,*time3,*data1,*data2,*chiSqArray, *data3, 
		ang,mean,adev,LFsdev,NPsdev, sdev, var,skew,curt,av, 
		Kendall_var1, Kendall_var2, kappa, B, NN1,NN2,freq1class,dominance;
	double 	NPmean,LFmean,chiSqmean;
	int	whichBetter;
	double	NN;
	NODETYPE * root, *node, **markedNodes, **nodeArray, *saveTree;
	extern 	int gIndex;
	extern 	double chiSq; /* declared in ObjFunc */
	int	i,j,k,success1,success2,kk,jj,irepcount,TotalReps;

	long    MaxGroupSize,nMark=0,count,s1,s2,maxS,
		countTaxa,countExc,binSize=10,nTaxa,nNodes,size,size2,countMonotypes, countMonophyletic, 
		countMonophyleticSpecies,*h,*histo,*histo2,*histoMonophyletic,
		*histoTotal,*histo2Total,theSeed=1,*histoB,sizeB;
	int    	
		N0=1,
		stemFlag=0,
		exclusive=1,
		rndBranchDur=0,
		withReplace=1,
		diversemodel=1,
		Yule_flag=0,
		CharEvol_flag=0, 
		save_flag=1,
		nreps=1,
		nrepsPerTree=1,
		nrepsPerBrRate=1,
		ntaxa=10,
		interval=10,
		infinite_flag=0,
		ratemodel=1, 
		gradual_flag=1, 
		silentFlag=0;
	double 	
		speciation=1.0,
		speciation2=1.0,
		extinction=0.0,
		sampling_fraction=1.0,
		start_rate=1.0,
		change_rate=0.0,
		min_rate=0.1,
		max_rate=2.0,
		rate_transition=0.0, 
		T=1.0;

/*print_mem_dbg(__FILE__,__LINE__);*/
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXONOMY"))
			{
			if (isEqual(LocalToken,"EXCLUSIVE"))
				exclusive=1;
			if (isEqual(LocalToken,"NESTED"))
				exclusive=0;
			}
		if (parse_assignment2("RNDBRANCHDUR"))
			if (isEqual(LocalToken,"YES"))
				rndBranchDur=1;
			else
				rndBranchDur=0;
		if (parse_assignment2("WITHREPLACE"))
			if (isEqual(LocalToken,"YES"))
				withReplace=1;
			else
				withReplace=0;
		if (parse_assignment2("CHAREVOL"))
			if (isEqual(LocalToken,"YES"))
				CharEvol_flag=1;
			else
				CharEvol_flag=0;
		if (parse_assignment2("NREPSPERBRRATE"))
			nrepsPerBrRate=strtod(LocalToken,&dummy);
		if (parse_assignment2("NREPSPERTREE"))
			nrepsPerTree=strtod(LocalToken,&dummy);
		if (parse_assignment2("BINSIZE"))
			binSize=strtod(LocalToken,&dummy);
		if (parse_assignment2("NMARK"))
			nMark=strtod(LocalToken,&dummy);
		if (parse_assignment2("NREPS"))
			nreps=strtod(LocalToken,&dummy);
		if (parse_assignment2("NTAXA"))
			ntaxa=strtod(LocalToken,&dummy);
		if (parse_assignment2("SPECIATION"))
			speciation=strtod(LocalToken,&dummy);
		if (parse_assignment2("SPECIATION2"))
			speciation2=strtod(LocalToken,&dummy);
		if (parse_assignment2("EXTINCTION"))
			extinction=strtod(LocalToken,&dummy);
		if (parse_assignment2("SAMPLING"))
			sampling_fraction=strtod(LocalToken,&dummy);
		if (parse_assignment2("INTERVAL"))
			interval=strtod(LocalToken,&dummy);
		if (parse_assignment2("SEED"))
			{
			theSeed=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("STARTRATE"))
			start_rate=strtod(LocalToken,&dummy);
		if (parse_assignment2("CHANGERATE"))
			change_rate=strtod(LocalToken,&dummy);
		if (parse_assignment2("MINRATE"))
			min_rate=strtod(LocalToken,&dummy);
		if (parse_assignment2("MAXRATE"))
			max_rate=strtod(LocalToken,&dummy);
		if (parse_assignment2("RATETRANSITION"))
			rate_transition=strtod(LocalToken,&dummy);
		if (parse_assignment2("T"))
			{
			T=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("INFINITE"))
			if (isEqual(LocalToken,"NO"))
				infinite_flag=0;
			else
				infinite_flag=1;
		if (parse_assignment2("GRADUAL"))
			if (isEqual(LocalToken,"NO"))
				gradual_flag=0;
		if (parse_assignment2("SAVE_FLAG"))
			if (isEqual(LocalToken,"YES"))
				save_flag=1;
			else
				save_flag=0;
		if (parse_assignment2("SILENT"))
			if (isEqual(LocalToken,"YES"))
				silentFlag=1;
			else
				silentFlag=0;
		if (parse_assignment2("STEM"))
			if (isEqual(LocalToken,"YES"))
				stemFlag=1;
			else
				stemFlag=0;
		if (parse_assignment2("RATEMODEL"))
			{
			if (isEqual(LocalToken,"NORMAL"))
				ratemodel=1;
			if (isEqual(LocalToken,"AUTOCORR"))
				ratemodel=2;
			}
		if (parse_assignment2("DIVERSEMODEL"))
			{
			Yule_flag=1;
			if (isEqual(LocalToken,"YULE"))
				diversemodel=1;
			if (isEqual(LocalToken,"BDBACKNORMAL"))
				diversemodel=5;
			if (isEqual(LocalToken,"BDBACK"))
				diversemodel=2;
			if (isEqual(LocalToken,"YULE_C"))
				diversemodel=3;
			if (isEqual(LocalToken,"BDFORWARD")) 
				diversemodel=4;
			if (isEqual(LocalToken,"RY1997")) 
				diversemodel=6;
			if (isEqual(LocalToken,"YULE_SISTERS")) 
				diversemodel=7;

			}
		}

	if (diversemodel==7)
				printf("[Expected ratio of sister group diversities=%f]\n",exp(T*(speciation2-speciation)));

	srand(theSeed);
	if (theSeed==1)
		{
		if (!silentFlag)doGenericAlert("WARNING: YOU ARE USING A DEFAULT SEED FOR RANDOM NUMBERS");
		}
	if (!silentFlag)printf("\n\n** r8s simulation run **\n\n");
	verbose=gNexDataPtr->RateBlockParms.verbose;

	kappa=speciation*T;
	Kendall_var1=SQR(speciation)/(2*(exp(kappa)-1.));
	Kendall_var2=Kendall_var1*SQR(sinh(0.5*kappa)/(0.5*kappa));
	time1=(double*)myMalloc((ntaxa-2)*sizeof(double));
	time2=(double*)myMalloc((ntaxa-2)*sizeof(double));
	time3=(double*)myMalloc((ntaxa-2)*sizeof(double));
	data1=(double*)myMalloc((nreps+1)*sizeof(double)); /* 1-offset array */
	data2=(double*)myMalloc((nreps+1)*sizeof(double)); /* 1-offset array */
	data3=(double*)myMalloc((nreps+1)*sizeof(double)); /* 1-offset array */
	chiSqArray=(double*)myMalloc((nreps+1)*sizeof(double)); /* 1-offset array */

#define SIM_LOOP 0	/* for doing lots of simulations */
#if SIM_LOOP

printf("ChangeRate,Transition,ChiSq,LF,NP,Which\n");

/* for (jj=0;jj<=10;jj++) */ 
  for (kk=0;kk<=10;kk++)
	{
	/* change_rate=jj*max_rate/20.; */
	rate_transition=kk/10.0;
#endif




/* start simulating */

#if 0  /* TEST  */
	printf("Simulation of Yule_forward\n");
	for (i=1;i<=nreps;i++)
		{
		NN1=Yule_forward(speciation, T, &B,stemFlag);
		data1[i]=NN1;
		if (stemFlag) N0=1; else N0=2;	
		data2[i]=(NN1-N0)/B;
		data3[i]=(log(NN1)-log(N0))/T;
		}
	moment(data1,nreps,&mean,&adev,&sdev,&var,&skew,&curt);
	printf("Test of Yule Forward routine: mean=%fvar=%f\n",mean,var);
	moment(data2,nreps,&mean,&adev,&sdev,&var,&skew,&curt);
	printf("Test of K-infinite estimator: mean=%fvar=%f\n",mean,var);
	moment(data3,nreps,&mean,&adev,&sdev,&var,&skew,&curt);
	printf("Test of K-1 routine: mean=%fvar=%f\n",mean,var);
	/*return;*/
#endif

	count=1;
	if (Yule_flag)
	    {
/*print_mem_dbg(__FILE__,__LINE__);*/
	    if (!silentFlag)printf("Diversification simulation:\nseed = %li\nnreps=%i\nntaxa=%i\nspec rate = %f\nextinct rate=%f\n",
		    theSeed,nreps,ntaxa,speciation,extinction);
	    switch (diversemodel)
			{
			case 1:
		    	  if (!silentFlag)printf("MODEL=Forward Yule model\n");
			  if (stemFlag)
				printf("(Stem group simulation:N0=1)\n");
			  else
				printf("(Crown group simulation:N0=2)\n");
			  break;
			case 2:
		    	  if (!silentFlag)printf("MODEL=Backward birth-death model\n");
		    	  if (!silentFlag)printf("(Root node time normalized to one)\n");
			  break;
			}
	  /*printf("Predicted estimator variance:K1(infinite)=%f K2 (k=1) = %f\n",
	     Kendall_var1, Kendall_var2);*/
	    for (i=1;i<=nreps;i++)
		{
		if (!silentFlag)
			printf ("...generating replicate tree number %i\n",i);
		/*root=BDTree(ntaxa,speciation, extinction,0.1);*/
		/*root = BDTreeForward(T,speciation, extinction,0.1);*/
		    {
		    switch (diversemodel)
			{
			case 1:
		    	  root = YuleTreeForward(T, speciation, &NN1, &B,stemFlag);
			  break;
			case 2: /* this is bdback without normalizing root to 1 */
			 /* root=BDTree(ntaxa,speciation, extinction,0.1);*/
			  root=BDback(ntaxa,speciation,0);
			 /* data1[i]=treeDurLength(root);*/
			  data1[i]=treeAgeSum(root)/numIntNodes(root);
			  printf("Duration of tree = %f age=%f\n", data1[i], root->time);
			  break;
			case 5: /* this is bdback normalizing root to 1 */
			 /* root=BDTree(ntaxa,speciation, extinction,0.1);*/
			  root=BDback(ntaxa,1.0,1); /* set speciation to 1; doesn't really matter anyway given the renormalization */
			 /* data1[i]=treeDurLength(root);*/
			  data1[i]=treeAgeSum(root)/numIntNodes(root);
			  printf("Duration of tree = %f age=%f\n", data1[i], root->time);
			  break;
			case 3:
			    root=Yule_C(ntaxa, speciation);
			    break;
			case 4:	
			    root=BDTreeForward(T, speciation, extinction,0.0);
			    break;
			case 6:	
			    root=RY_1997(ntaxa, T, speciation,extinction,sampling_fraction); // Rannala Yang 1997 model
			    break;
			case 7:	
			    root=SisterGroupYule(T, speciation, speciation2, &NN1, &NN2, &B);
			    break;


			}
		 /*   data1[i]=NN1;
		    data2[i]=(NN1-2)/B;
		    data3[i]=B;*/
		    if(save_flag) /* now do this by default! */
			doSaveTree(root);
		    else
			DisposeTree(root);
		    }
/*****/
#if SIMLOOP

/*********** !!!!!!!!!! note that the following doObjFunc calls have invalid arg lists !!!!!!! **/

		gIndex=0;
		tree2aTimeArray(root,time1);

		gnpexp=gNexDataPtr->RateBlockParms.npexp;  /* KLUDGE */
		(void)doObjFunc(objLangFitch,root,"Simulated",LaF,&success1);
		gIndex=0;
		tree2aTimeArray(root,time2);


		(void)doObjFunc(objNP,root,"Simulated",NP,&success2);
		gIndex=0;
		tree2aTimeArray(root,time3);



		/*for (k=0;k<ntaxa-2;k++)
		    printf("%f %f %f\n",time1[k],time2[k],time3[k]);*/
		if (success1 && success2) /* store rep results only if both opts work*/
			{
			data1[count]=euclid_distance(time1,time2,ntaxa-2);
			data2[count]=euclid_distance(time1,time3,ntaxa-2);
			chiSqArray[count]=chiSq; /* a global */
			if (verbose)
			    printf("%f %f\n",data1[count],data2[count]);
		/*if (data2[count]>0.01)
		  {
		  printf("$$$%f\n", data2[count]);
		  for (k=0;k<ntaxa-2;k++)
		    printf("%f %f %f\n",time1[k],time2[k],time3[k]);
		  gNexDataPtr->RateBlockParms.verbose=2;
		  (void)doObjFunc(objNP,root,"Simulated",NP,&success2);
		  DrawTree(root, 1, 0);
		  DrawTree(root, 2, 0);
		  DrawTree(root, 4, 0);
		  printtree(root);
		  make_parens(root, 0);
		  exit(1);
		  }*/
			++count;
			}
		else
			doGenericAlert("WARNING: LF or NP failed\n");

		DisposeTree(root); /*******!!  !!********/
#endif
		} /* end nreps */
	       } /* endif */
		/* BDDiversity(ntaxa,speciation, extinction,0.1,interval); */

	--count;
/***************************

	Do the character evolution simulation


	For the normally-distributed model, the parameters start_rate and change_rate correspond to
	the mean and standard deviation of the normal respectively. The min and max values are still respected.


****************************/

		if(CharEvol_flag)
		  {
		  if (gNexDataPtr->isTrees)
		    {
			if (!silentFlag)
				{
				printf("\nBranch evolution simulation:\nseed=%li\n\nrate transition=%f\n",
					theSeed, start_rate,change_rate,min_rate,max_rate,rate_transition);
				printf("Gradual rate change=%i\n",gradual_flag);
				printf("Infinite=%i\n",infinite_flag);
				if(ratemodel==1)
					{
					printf("RATE MODEL:Normally distributed\n");
					printf("with parameters: mean=%f, sdev=%f, minrate=%f, maxrate=%f\n", 
						start_rate, change_rate, min_rate,max_rate);
					}
				else if (ratemodel==2)
					{
					printf("RATE MODEL:Autocorrelated\n");
					printf("with parameters: startrate=%f change rate=%f minrate=%f maxrate=%f\n", 
					start_rate,change_rate,min_rate,max_rate);
					printf("transition probability=%f change amount=%f\n", rate_transition, change_rate);
					}
				}
			lnode=gNexDataPtr->inTrees;
			LISTLOOP (lnode)
				{
				i=1;
				thisTree=lnode->item;
				for (j=1;j<=nrepsPerTree;j++)
					{
//					set_branch_rates(thisTree->root,start_rate,change_rate, min_rate,max_rate,rate_transition,gradual_flag, ratemodel);
					saveTree=copyTree(thisTree->root); // makes a deep copy because stuff will be overwr
					for (jj=1;jj<=nrepsPerBrRate;jj++)
					  {
//					  thisTree->root=copyTree(saveTree);
//					  set_branch_lengths(thisTree->root,infinite_flag);
//				    DrawTree(thisTree->root,1, 0);
//					printtree(thisTree->root);	
					  if (!silentFlag)
						{
						printf("\n ** Tree %i (Rate Replicate %i, Branch Length Rep %i)\n", i,j,jj);
						printf("tree SIMTREE = ");
							make_parens(thisTree->root, 0); /* TD of phylogram*/
						printf(";\n");
						}

#if 1
					// note if you estimate times, this will bollocks up the set rates, branches above..
trueAge=50.0;
rangeFactor=0.05; // increase to make the range of min and max branch lengths larger
ksteps=20;
gNexDataPtr->RateBlockParms.num_time_guesses=2;
for (ib=-ksteps/2;ib<=+ksteps/2;ib++)
	for (ic=-ksteps/2;ic<=+ksteps/2;ic++)
		{
lb=floor(pow(10, log10(5000)+ib*rangeFactor)); // make sure these are ints because of that lousy rounding problem with gradients
lc=floor(pow(10, log10(5000)+ic*rangeFactor));
printf("\n--%6.1f %6.1f\n\n",lb,lc);
thisTree->root=copyTree(saveTree);
(thisTree->root)->free=0;
(thisTree->root)->time=T;
(thisTree->root)->nodeIsConstrainedMax=0;
(thisTree->root)->nodeIsConstrainedMin=0; 
root=thisTree->root;
nodea=find_taxon_name(root,"a");
nodeb=find_taxon_name(root,"b");
nodec=find_taxon_name(root,"c");
nodeInt=nodec->anc;
nodea->length=10000;
nodeInt->length=5000;
		nodeb->length=lb;
		nodec->length=lc;
		printtree(root);
		DrawTree(root,1, 0);

					  doObjFunc(thisTree,LaF,1,TN,&success1);
					  print_ages(thisTree->root, 1.0,1.0,thisTree->method); 
		T_LF[ib+ksteps/2][ic+ksteps/2]=nodeInt->time;
DisposeTree(thisTree->root);

/*
thisTree->root=copyTree(saveTree);
(thisTree->root)->free=0;
(thisTree->root)->time=T;
(thisTree->root)->nodeIsConstrainedMax=0;
(thisTree->root)->nodeIsConstrainedMin=0; 
root=thisTree->root;
nodea=find_taxon_name(root,"a");
nodeb=find_taxon_name(root,"b");
nodec=find_taxon_name(root,"c");
nodeInt=nodec->anc;
nodea->length=10000;
nodeInt->length=5000;
		nodeb->length=lb;
		nodec->length=lc;
//					  bestSmooth=doCrossV(thisTree,LaF,1,TN,1.0,0.5,1,0);
					  bestSmooth=doCrossV(thisTree,PENLIKE,1,TN,0,0.5,8,0);
//					  gNexDataPtr->RateBlockParms.smoothing=bestSmooth;
DisposeTree(thisTree->root);
*/

thisTree->root=copyTree(saveTree);
(thisTree->root)->free=0;
(thisTree->root)->time=T;
(thisTree->root)->nodeIsConstrainedMax=0;
(thisTree->root)->nodeIsConstrainedMin=0; 
root=thisTree->root;
nodea=find_taxon_name(root,"a");
nodeb=find_taxon_name(root,"b");
nodec=find_taxon_name(root,"c");
nodeInt=nodec->anc;
nodea->length=10000;
nodeInt->length=5000;
		nodeb->length=lb;
		nodec->length=lc;
					  gNexDataPtr->RateBlockParms.smoothing=0.0001;
					  doObjFunc(thisTree,PENLIKE,1,TN,&success1);
		T_PL[ib+ksteps/2][ic+ksteps/2]=nodeInt->time;
					  print_ages(thisTree->root, 1.0,1.0,thisTree->method); 
DisposeTree(thisTree->root);

		}
//					  DisposeTree(thisTree->root);

for (ib=0;ib<=ksteps;ib++)
	{
	for (ic=0;ic<=+ksteps;ic++)
		printf("%6.1f\t",fabs(trueAge-T_LF[ib][ic]));
	printf("\n");
	}
printf("\n");
for (ib=0;ib<=ksteps;ib++)
	{
	for (ic=0;ic<=+ksteps;ic++)
		printf("%6.1f\t",fabs(trueAge-T_PL[ib][ic]));
	printf("\n");
	}
printf("\n");

printf("\t");
for (ic=0;ic<=ksteps;ic++)
	{
	lc=pow(10, log10(5000)+(ic-ksteps/2)*rangeFactor);
	printf("%6.1f\t",lc);
	}
printf("\n");

for (ib=0;ib<=ksteps;ib++)
	{
	lb=pow(10, log10(5000)+(ib-ksteps/2)*rangeFactor);
	printf("%6.1f\t",lb);
	for (ic=0;ic<=+ksteps;ic++)
		{
		printf("%6.1f\t",fabs(trueAge-T_LF[ib][ic])-fabs(trueAge-T_PL[ib][ic]));
		}
	printf("\n");
	}
printf("\n");
for (ib=0;ib<=ksteps;ib++)
	{
	for (ic=0;ic<=+ksteps;ic++)
		{
		printf("%6.1f\t",fabs(trueAge-T_LF[ib][ic])-fabs(trueAge-T_PL[ib][ic]));
		}
	printf("\n");
	}
printf("\n");
for (ib=0;ib<=ksteps;ib++)
	{
	for (ic=0;ic<=ib;ic++)
		if (fabs(trueAge-T_PL[ib][ic]) < fabs(trueAge-T_LF[ib][ic]))
			printf("1");
		else
			printf("0");
	printf("\n");
	}
printf("\n");
#endif
					  }
					thisTree->root=saveTree;
					}
				++i;
				}
		    }
		  else
		    printf("No trees currently in memory\n");
		  }





#if 0
	moment(data1,count,&LFmean,&adev,&LFsdev,&var,&skew,&curt);
	if (nreps-count)
		printf("There were %i failed replicates\n",nreps-count);
	/*printf("distance sim to LF: mean=%f stdev=%f\n",mean,sdev);*/
	moment(data2,count,&NPmean,&adev,&NPsdev,&var,&skew,&curt);
	/*printf("distance sim to NP: mean=%f stdev=%f\n",mean,sdev);*/
	moment(chiSqArray,count,&chiSqmean,&adev,&sdev,&var,&skew,&curt);
	if (NPmean<LFmean)
		whichBetter=1; /* NP has lower mean error */
	else
		whichBetter=0;

	printf("%f\t%f\t%f\t%f (+-%f)\t%f (+-%f)\t%i\n",
		change_rate, rate_transition,chiSqmean,LFmean,LFsdev/sqrt(nreps), 
		NPmean,NPsdev/sqrt(nreps), whichBetter);
#endif
#if SIM_LOOP
	} 
#endif


/***** RANDOM BRANCH SAMPLING SIMULATION ****/

//  NB. The nested taxon model always picks nodes without replacement...see nextRndNode()

/*****/

#if 0
		// Allocation and initialization

		    markedNodes=(NODETYPE**)myMalloc((nMark+1)*sizeof(NODETYPE*)); /* 1-offset array */
		    lnode=gNexDataPtr->inTrees;
			/* following lines assume that all trees in list have some number of nodes! */
			/* And we're going to add all the histogram entries together across trees */
		    thisTree=lnode->item;
		    nNodes=numNodes(thisTree->root);
		    nTaxa=numdesc(thisTree->root);
			 /* Following are all (1-off) histos: # one bin for each possible taxon size */
		    MaxGroupSize=nTaxa; /* Following arrays need to range from a group size of 0 to nTaxa*/
		    histo=(long *)myMalloc((1+MaxGroupSize)*sizeof(long));
		    histoB=(long *)myMalloc((1+MaxGroupSize)*sizeof(long));
		    histoMonophyletic=(long *)myMalloc((1+MaxGroupSize)*sizeof(long)); 
		    histo2=(long *)myMalloc((1+MaxGroupSize)*sizeof(long)); 
		    histoTotal=(long *)myMalloc((1+MaxGroupSize)*sizeof(long));
		    histo2Total=(long *)myMalloc((1+MaxGroupSize)*sizeof(long)); 

			TotalReps=nreps*nrepsPerTree;
		    RepCount=(double *)myMalloc((TotalReps+1)*sizeof(double)); // 1-offset for moment function
		    RepMean=(double *)myMalloc((TotalReps+1)*sizeof(double)); // 1-offset for moment function
		    RepDominance=(double *)myMalloc((TotalReps+1)*sizeof(double));
		    RepFreq1Class=(double *)myMalloc((TotalReps+1)*sizeof(double));
		    RepFractMonophyletic=(double *)myMalloc((TotalReps+1)*sizeof(double));
		    RepMonophyleticSpecies=(double *)myMalloc((TotalReps+1)*sizeof(double));

		    irepcount=0;
			for (i=0;i<=nTaxa;i++)
							{
							histoTotal[i]=0;
							histo2Total[i]=0;
							histoB[i]=0;
							}

			if (exclusive)	// for exclusive models we never want to sample with replacement!
				withReplace=0;
				
		    LISTLOOP (lnode)
			    	{
			    	thisTree=lnode->item;
					nodeArray = newAllNodeArray(thisTree);
					for (j=1;j<=nrepsPerTree;j++)
					   {
			    		++irepcount;
			    		printf("...working on replicate %i\n",irepcount);
						for (i=0;i<=MaxGroupSize;i++)
							{
							histo[i]=0;
							histoMonophyletic[i]=0;
							histo2[i]=0;
							}
					   if (exclusive) // under exclusive model, the root node is always part of a taxon; mark it in 
					   				  // 0-th element of this array
					   		markedNodes[0]=thisTree->root;
					   if (rndBranchDur) // the other nodes marked are in elements 1..N	
							RandomBranches(thisTree,nNodes,nodeArray,nMark,markedNodes,withReplace); /* NB!  If we have a very long branch, the rand chars will hit it often, and the naive w/o replacement algorithm will thrash */
					   else
							markRandomNodes(thisTree,nMark,markedNodes);  /* NB! Change the dynamic allocation in this routine to static...*/
								
					   countTaxa=0;
					   countMonotypes=0;
					   countMonophyletic=0;
					   countExc=0;
					   countMonophyleticSpecies=0;
/*
	Under the exclusive model there will be N+1 groups recognized, where N is the number of marks;
	Under the nested model there will be N groups recognized
*/
					   for (i=0;i<=nMark;i++)
								{
								if (i==0 && !exclusive) // skip the root node (in ..[0]) if nested model
									continue;
								++countTaxa;
								node=markedNodes[i];
								unMarkNode(node);/* have to unmark this node for next function to work right */
								size = numUnMarkedDesc(node); // NB! sometimes you can get a size of zero! Then this is an "orphan"	
								size2= numdesc(node); // size of this taxon assuming it's monophyletic
								++histo[size]; // these are the group sizes possibly paraphyletic (for exclusive sampling)
								if (size<=9) sizeB=0;
								if (size>9 && size<=99) sizeB=1;
								if (size>99 && size<=999) sizeB=2;
								if (size>999) sizeB=3;
								/*sizeB=(long)floor((size-1)/10.0)+1;*/
								++histoB[sizeB];
								++histoTotal[size]; // keep track across reps
								countExc+=size;
								++histo2[size2]; // these are the monophyletic group sizes (for clade sampling)
								++histo2Total[size2]; // keep track across reps
								if (size2==1)
									++countMonotypes;
								if (size>1 && size==size2) /* don't count the "monophyletic" monotypes*/
									{
									++countMonophyletic;
									++histoMonophyletic[size];
									countMonophyleticSpecies+=size; // num of species in monophyletic higher taxa
									}
								markNode(node);
//if (size2>10) printf("Mono:%i Exc:%i\n",size2,size);
								}


// At this point store the stats for each replicate (regardless of whether replicate across one tree or many)							
						if (exclusive)
							h=histo;
						else
							h=histo2;							
						histoStat(h, MaxGroupSize,nTaxa, &count, &mean, &freq1class, &maxS, &dominance);
						RepCount[irepcount]=count; // number of taxa > 0 found 
						RepMean[irepcount]=mean;
						RepDominance[irepcount]=dominance;
						RepFreq1Class[irepcount]=freq1class;
						RepFractMonophyletic[irepcount] = (double)countMonophyletic/(countTaxa-countMonotypes);
						RepMonophyleticSpecies[irepcount]= (double)countMonophyleticSpecies/(nTaxa-countMonotypes);
printf("Monotypes and max size:%i %i %i\n",irepcount,histo[1],maxS);
// printf("%i %f %f %f\n",irepcount,RepMean[irepcount],RepDominance[irepcount],RepFreq1Class[irepcount]);
					   } // end nrepspertree
				myFree(nodeArray);
				} // end listloop

		// Print results

		    if (!silentFlag)
			{
		    	printf("******************************************\n\n");
		    	printf("Results of simulation of random taxonomies\n\n");
		    	if (exclusive) printf ("Taxonomy consists of exclusive groups\n"); 
		    	else printf ("Taxonomy consists of nested groups\n");
		    	if (rndBranchDur) printf ("Branches are sampled according to durations\n"); 
		    	else printf ("Branches are sampled equally\n");
		    	if (withReplace) printf("Branches are sampled with replacement\n");
		    	else printf ("Branches are sampled without replacement\n");

		    	printf("Average histogram for model across %li replicate tree simulations\n",nreps);
		    	for (i=1;i<=MaxGroupSize;i++)
			    if (exclusive)
				{
				if (histoTotal[i]>0)
					printf("%li\t%f\n",i,histoTotal[i]/(float)TotalReps);
				}
			    else
				if (histo2Total[i]>0)
					printf("%li\t%f\n",i,histo2Total[i]/(float)TotalReps);

			}


		// Make a plot of the sum total histogram across all replicates

		    if (exclusive)
			h=histoTotal;
		    else
			h=histo2Total;

		    histoStat(h, MaxGroupSize,nTaxa, &count, &mean, &freq1class, &maxS, &dominance); // just to get count
		    X = (double *)myMalloc(count*sizeof(double));
		    Y = (double *)myMalloc(count*sizeof(double));
		    count=0;
		    for (i=1;i<=MaxGroupSize;i++)
				if (h[i]>0) 
					{
					X[count]=log10(i);
					Y[count]=log10(h[i]/(float)TotalReps);
					++count;
					}
		    for (i=0;i<count;i++)
			printf ("%f\t%f\n",X[i],Y[i]);
		    if (!silentFlag)			
		    	dumbPlot(X,Y,count);

			moment(RepCount,TotalReps,&mean,&adev,&sdev,&var,&skew,&curt);
			printf("Number of higher units generated across reps: mean=%f var=%f\n",mean,var);
			moment(RepMean,TotalReps,&mean,&adev,&sdev,&var,&skew,&curt);
			printf("Taxon size across reps: mean=%f var=%f\n",mean,var);
			moment(RepDominance,TotalReps,&mean,&adev,&sdev,&var,&skew,&curt);
			printf("Dominance across reps: mean=%f var=%f\n",mean,var);
			moment(RepFreq1Class,TotalReps,&mean,&adev,&sdev,&var,&skew,&curt);
			printf("Size 1 class (monotypes) across reps: mean=%f var=%f\n",mean,var);
			moment(RepFractMonophyletic,TotalReps,&mean,&adev,&sdev,&var,&skew,&curt);
			printf("Fraction of monophyletic nonmonotypes: mean=%f var=%f\n",mean,var);
			moment(RepMonophyleticSpecies,TotalReps,&mean,&adev,&sdev,&var,&skew,&curt);
			printf("Fraction of all species found in nonmonotypic monophyletic higher taxa: mean=%f var=%f\n",mean,var);
		    	for (i=0;i<=3;i++)
				printf("%li\t%f\n",i,histoB[i]/(float)TotalReps);

			myFree(X);
			myFree(Y);
			myFree(histo);
			myFree(histo2);
			myFree(histoTotal);
			myFree(histo2Total);
			myFree(histoMonophyletic);
			myFree(markedNodes);
			myFree(RepMean);
			myFree(RepDominance);
			myFree(RepFreq1Class);
#endif

/*****/
/*****/



/*  put this as a command such as 'statistics'
		moment(data1,nreps,&mean,&adev,&sdev,&var,&skew,&curt);
		printf("Test of Mean Duration in BDBack: mean=%f var=%f\n",mean,var);
		moment(data2,nreps,&mean,&adev,&sdev,&var,&skew,&curt);
		printf("Test of K-infinite estimator in YuleTreeForward: mean=%fvar=%f\n",mean,var);
		moment(data3,nreps,&mean,&adev,&sdev,&var,&skew,&curt);
		printf("Test of mean-B in YuleTreeForward: mean=%fvar=%f\n",mean,var);
*/
	myFree(data1);
	myFree(data2);
	myFree(chiSqArray);
	myFree(time1);
	myFree(time2);
	myFree(time3);
	return;						
}


/****************************************************************/
static void histoStat(long h[], long N, long nTaxa, long *count, double *mean, double *freq1class, long *maxS, double *dominance)
{
// NB! The histo array, h, is of size N+1 and is zero-offset, because the min size class is zero for our stuff
// Thus, note that it should be called by ...(h,N,...), rather than (...,N+1).

// Note that h[0] contains the count of those taxa orphaned with no terminals in them at all. 
// The number of actual taxa observed in the classification is therefore slightly less than the number of randomly selected nodes for marking


long i,s2=0;
*maxS = -1;
*count=0;
for (i=1;i<=N;i++)
	{
	
	if (h[i]>0) 
		{
		*maxS = i;  // this is the maximum size observed in the histogram
		(*count)+=h[i];	// total number of points (higher taxa)
		s2+=h[i]*i;
		}
	}
	
*mean=(double)s2 /(*count);
*freq1class = (double) h[1]/(*count);
*dominance = (double)(*maxS)/nTaxa;

return;
}

/****************************************************************/

static void doSaveTree(NODETYPE *root)  /* adds a simulated tree to the tree list */
{
extern int curTree;
int flag=0, numtips, numinternals, roottomy;
char *stemp, *tree_name="SIMTREE", *TD="";
PtrList aTreeList;
TREE aTree;


	if (gNexDataPtr->inTrees == NULL)  /* if this is the first tree */
	    {
	    gNexDataPtr->inTrees=pNewListAlt(sizeof(struct treetype));
	    aTreeList=gNexDataPtr->inTrees;
	    }
	else				/* if a later tree */
	    aTreeList=pListAddNode(gNexDataPtr->inTrees,sizeof(struct treetype));
		
	(void)appendStrList(gNexDataPtr->TDLabelList,tree_name);
	(void)appendStrList(gNexDataPtr->TDList,TD);			
	++curTree;
	gNexDataPtr->NumTrees=curTree;
	gNexDataPtr->isTrees=1;
	aTree=aTreeList->item;
	if (root)
	    {
	    init_node_ids(root, 0); 
	    init_free(root); /* sets default to estimate all internal nodes but root */ 
	    aTree->root=root;
	    aTree->name=DupStr(tree_name);
	    aTree->TD=DupStr(TD);
	    TreeStats(root, &numtips, &numinternals, &roottomy);
	    aTree->numTaxa=numtips;
	    aTree->numBranches=numBranches(root);
	    aTree->basalTomy=roottomy;
	    }

return;						
}
/****************************************************************/
void doWuLiCommand(NODETYPE * root)
{
	int id[3], ix,jx;
	char c, *dummy;
	

	if ((gNexDataPtr->isChars==0) || ( gNexDataPtr->isTaxa==0)
		|| gTaxaSet == NULL)
		return;	/* don't have the right data from NEXUS file, so bail */

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("RRTYPE"))
			{
			if (isEqual(LocalToken,"WULI"))
				gNexDataPtr->RateBlockParms.RRtype = WULI;
			if (isEqual(LocalToken,"STEEL"))
				gNexDataPtr->RateBlockParms.RRtype = STEEL;
			if (isEqual(LocalToken,"TAJIMA"))
				gNexDataPtr->RateBlockParms.RRtype = TAJIMA;
			if (isEqual(LocalToken,"MIKE"))
				gNexDataPtr->RateBlockParms.RRtype = MIKE;
			}
		if (parse_assignment2("BS"))
		    {
		    if (isEqual(LocalToken,"YES"))
			gNexDataPtr->RateBlockParms.isBS=1;
		    else
			gNexDataPtr->RateBlockParms.isBS=0;

		    }
		if (parse_assignment2("NREPS"))
			gNexDataPtr->RateBlockParms.NReps=strtod(LocalToken,&dummy);
		if (parse_assignment2("SEED"))
			{
			gNexDataPtr->RateBlockParms.seed=strtod(LocalToken,&dummy);
			srand(gNexDataPtr->RateBlockParms.seed);
			}
		}

	doRelativeRates(gTaxaSet,root);
	
	return;
}
/****************************************************************/
void doTaxaSetCommand(void)
{
	int id[3], ix,jx;
	char c, *dummy;

	if ( gNexDataPtr->isTaxa==0)
		return;	/* don't have the right data from NEXUS file, so bail */

	if(gTaxaSet)  /* get rid of old list if it's present */
		    {
		    freeStrList(gTaxaSet);
		    }

	gTaxaSet=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	
			appendStrList(gTaxaSet,aTokenPtr); /* store the label */
	
	printf("Using the following taxa:\n");xprintStrList(gTaxaSet);
	return;
}

/****************************************************************/
void doExSets(void)

/* Sets up an exclusion set array  in which a zero means excluded and 1 means included.
Format is 'exsets n1 n2 n3 n4 - n5 n6;'  NOTE THAT THERE MUST BE A SPACE BEFORE AND AFTER
THE DASH--THIS IS A NON NEXUS COMPLIANT WORKAROUND, but the NEXUS format does not recognize
dashes as punctuation since they can also represent gaps (?) FIX Later
Each time an exclusion set is invoked, the array is reset to match command */
{
	char* dummy;
	long icur,ilast,ix;
	int *excArray;
	excArray=gNexDataPtr->excArray;
	for (ix=0;ix<gNexDataPtr->NChars;ix++)
		excArray[ix]=1;	/* initializes exclusion set array for use in block*/
/*	fprintf(fpOut2,"[!NOTE: Some sites excluded in following analyses]\n");*/
	while (!isEqual(aTokenPtr=nextToken(),";"))	/* if its not a ';' it should be a number*/
		{
	
		if (  isdigit(*aTokenPtr)  )
			{
				icur=strtod(aTokenPtr,&dummy);
				ilast=icur;
				excArray[icur-1]=0;	/* this is a zero offset array */
			}
		else
			if (isEqual(aTokenPtr,"-"))
				{
				aTokenPtr=nextToken();
				icur=strtod(aTokenPtr,&dummy);
				for (ix=ilast;ix<=icur;ix++)
					excArray[ix-1]=0;	/* this is a zero offset array */
					
				}
			
		

		}
			
					
return;
}
/****************************************************************/
void doSitesCommand(int what) /* include or exclude positions */
{
long ix;

switch (what) 
	{
	case 0:
		for (ix=0;ix<gNexDataPtr->NChars;ix++)
			gNexDataPtr->excArray[ix]=1;
		printf("\n\n*** All sites included from now on ***\n\n\n");
		break;
	case 1:
		for (ix=0;ix<gNexDataPtr->NChars;ix++)
			if ( (ix+1)/3 != (ix+1)/3.0)
				gNexDataPtr->excArray[ix]=0;
		printf("\n\n*** First and second positions excluded from now on ***\n\n\n");
		break;
	case 3:
		for (ix=0;ix<gNexDataPtr->NChars;ix++)
			if ( (ix+1)/3 == (ix+1)/3.0)
				gNexDataPtr->excArray[ix]=0;
		printf("\n\n*** Third positions excluded from now on ***\n\n\n");
		break;

	}
return;
}

/****************************************************************/
/****************  MISCELLANEOUS FUNCTIONS **********************/
/****************************************************************/

void freeNexusStructure(struct NexDataType *nex)
{
freeStrList(nex->TaxaList);
freeStrList(nex->TDList);
freeStrList(nex->TDLabelList);
if (nex->isChars)
	freeStrList(nex->DMList);	/* this won't be allocated if no characters */
freeStrList(nex->TransList);
myFree(nex);


return;
}
/****************************************************************/

void doError(char* p[], int which)
{
doGenericAlert(p[which]);
}
/****************************************************************/
void TreeSummary(int whichTree)
{
	NODETYPE *root;
	char * TreeName, *TD;
	int numTips,numInternals, roottomy;
	TreeName=getkthStr(gNexDataPtr->TDLabelList,whichTree);
	TD=getkthStr(gNexDataPtr->TDList,whichTree);
	root=string_to_tree(TD);
	if (root != NULL)
			{
			TreeStats(root,&numTips,&numInternals, &roottomy);
			DisposeTree(root);
			}
	printf("Processing tree %i (%s) (taxa=%i; No. internal nodes = %i; Basal tomy=%i)\n",
				whichTree, TreeName,numTips,numInternals,roottomy); 

	return;

}
/****************************************************************/
int parse_assignment(char * target,char ** token)

/* on entry 'aTokenPtr' contains the putative first word of a three token
assignment statement of the form 'word1=word2'.  This function checks to see
if word1 is the same as 'target' and if so, it returns the address of a string
containing 'word2' or NULL if an error occurs.  aTokenPtr is
set to the last token in the assignment statement
If no match, aTokenPtr is left unchanged!! */

/*** BAD CODE *** causes memory leaks, probably failing to 
		free LocalTokens */

{
		if (isEqual(aTokenPtr,target))
			{
			aTokenPtr=nextToken();
			/*if (aTokenPtr==NULL) return 0;*/
			if (!isEqual(aTokenPtr,"="))
				{
				printf("Bad assignment statement=(%s)\n",aTokenPtr);
				fatal("exiting...");
				}
			aTokenPtr=nextToken();
			*token = DupStr(aTokenPtr);
			return 1;
			}
	return 0;
}

/****************************************************************/
int parse_assignment2(char * target)

/* on entry 'aTokenPtr' contains the putative first word of a three token
assignment statement of the form 'word1=word2'.  This function checks to see
if word1 is the same as 'target' and if so, it returns the address of a string
containing 'word2' or NULL if an error occurs.  aTokenPtr is
set to the last token in the assignment statement
If no match, aTokenPtr is left unchanged!! */

{
		if (isEqual(aTokenPtr,target))
			{
			aTokenPtr=nextToken();
			if (!isEqual(aTokenPtr,"="))
				{
				printf("Bad assignment statement=(%s)\n",aTokenPtr);
				fatal("exiting...");
				}
			aTokenPtr=nextToken();
			if (strlen(aTokenPtr)< MAX_LOCAL_TOKEN_SIZE -1)
				strcpy(LocalToken,aTokenPtr);
			else
				fatal("local token size exceeded\n");
			return 1;
			}
	return 0;
}

/****************************************************************/


void checkMatrix(void)
{
int itaxa;
char* c;
for (itaxa=0;itaxa<gNexDataPtr->NTaxa;itaxa++)
	{
	c=getkthStr(gNexDataPtr->DMList,(long)(itaxa+1));
	while(*c++)
		if ( (*c=='{') || (*c=='}'))
			{
			doGenericAlert("Polymorphism not allowed: Do not invoke rate tests");
			return;
			}
	}
	
return;	
}
/***************/
#if 0
int gNComp;
static void doCrossV(PtrList TreeList, int method,double EstMult,double PrdMult,double cvStart,double cvInc,double cvNum)

/*  does a cross validation with the range of the tuning parameter set to run from
	[cvStart, cvStart+cvInc, ...,cvStart+cvInc*(cvNum-1)]

NB! Haven't added variable 'algorithm' here yet

*/
{
#define isEven(k) ((k)/2 == ((k)/2.0))
TREE tree1,tree2;
PtrList p;
int i,j,success,nTrees,collFlag=0;
double * cvScore,cvSum;

double smooth;

nTrees=pLengthList(TreeList);
if (!isEven(nTrees))
	{
	doGenericAlert ("Must have even number of trees to do cross validation");
	return;
	}
else
     {
     gNComp=nTrees/2;
     cvScore = vector(1,gNComp);
     for (j=0;j<cvNum;j++)
	{
        cvSum=0.0;
	smooth=pow(10.0,j*cvInc+cvStart);
	gNexDataPtr->RateBlockParms.smoothing=smooth;
	p=TreeList;
	for (i=1;i<=gNComp;i++)
		{
		tree1=p->item;
		tree2=p->next->item;
		if (j==0)		/* only do this on the first pass; otherwise we keep multiplying the lengths */
			{
			if (EstMult != 1.0)traverseMultiplyLength(tree1->root, EstMult);
			if (PrdMult != 1.0)traverseMultiplyLength(tree2->root, PrdMult);
			}
		if (collapseLengthsTree2Tree(tree1,tree2))
			collFlag=1;				/* set this if some tree had to collapse a branch */
		copyLengthsTree2Tree(tree1->root,tree2->root); /* put lengths from tree2 onto field nodeReal of tree1 */
		doObjFunc(tree1,method,0,&success);
		cvScore[i]=cvSquareError(tree1,method);
		cvSum+=cvScore[i];
		p=p->next->next;
		} 	
	printf("\nCross Validation Analysis\n\n");
	for (i=1;i<=gNComp;i++)
		printf ("Cross Validation Score [%2i] = %f\n",i,cvScore[i]);
	printf("Cross Validation Score Total:smoothing = %f CV=%f\n",smooth,cvSum/gNComp);
	if (collFlag)
		printf("NOTE: Some partitions had 0-length branches that had to be collapsed to estimate cv scores\n");
	}
    }
}

#endif
                                                                                            r8s/TNwrapper.c                                                                                     0000644 0000766 0000120 00000005067 10357562456 013356  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "Maximize.h"
#include "NRCvectorUtils.h"
#include "ConstrOpt.h"
#include "MyUtilities.h"
#include "TreeUtils.h"
#include "memory.h"
#include "ObjFunc.h"
#include "TNwrapper.h"
double 		(*gObj)(double []);
void		(*gGrad)(double [], double []);

int TNwrapper
	(
	int numvar,
	double x[],
	double 		(*objective)(double []),
	void		(*gradient)(double [], double []),
	double		*max_obj
	)

{
int IERROR,LW,*IPIVOT,i;
double f;
double *W;
extern double *gLOW,*gHIGH;
double *g;


void tnbc_(int *, int *,double [],double *,double [],double [],int *k, 
	void (*)(int *,double [],double *, double []),
	double [], double [], int []);

LW=14*numvar;
W=(double*)myMalloc(LW*sizeof(double));
IPIVOT=(int*)myMalloc(numvar*sizeof(int));
g=vector(1,numvar);

f=(*objective)(x); /* f at starting point */
gradient(x,g); /* g at starting point */  
gObj=objective; /* globals used by sfun */
gGrad=gradient;

/*
for (i=0;i<numvar;i++)
	printf("[%2i]:%f\t%f\n",i,gLOW[i],gHIGH[i]);
*/
tnbc_(&IERROR, &numvar, x+1, &f, g+1, W, &LW, sfun_,gLOW,gHIGH,IPIVOT);
/*
printf("BOUNDS REACHED:\n");	
for (i=0;i<numvar;i++)
	printf("[%2i]:%i\n",i,IPIVOT[i]);
*/
*max_obj=f;
free_vector(g,1,numvar);
myFree(W);
myFree(IPIVOT);
return IERROR;

}

/* whew NASTY to get the 1-offset/0-offset shit between NRC, C, and FORTRAN */
/* Integers and doubles get passed as pointers, as do arrays. */


void sfun_(int *N,double X[],double *F, double G[])

{
*F=(*gObj)(X-1);
(*gGrad)(X-1,G-1);
return;
}




#if 0
C***********************************************************************
C EASY TO USE, NO BOUNDS
C***********************************************************************
C MAIN PROGRAM TO MINIMIZE A FUNCTION (REPRESENTED BY THE ROUTINE SFUN)
C OF N VARIABLES X
C
      DOUBLE PRECISION  X(50), F, G(50), W(700)
      EXTERNAL          SFUN
C
C DEFINE SUBROUTINE PARAMETERS
C N  - NUMBER OF VARIABLES
C X  - INITIAL ESTIMATE OF THE SOLUTION
C F  - ROUGH ESTIMATE OF FUNCTION VALUE AT SOLUTION
C LW - DECLARED LENGTH OF THE ARRAY W
C
      N  = 10
      DO 10 I = 1,N
         X(I) = I / FLOAT(N+1)
10      CONTINUE
      F  = 1.D0
      LW = 700
      CALL TN (IERROR, N, X, F, G, W, LW, SFUN)
      STOP
      END
C
C
C     SUBROUTINE SFUN (N, X, F, G)
C     DOUBLE PRECISION  X(N), G(N), F, T
C
C ROUTINE TO EVALUATE FUNCTION (F) AND GRADIENT (G) OF THE OBJECTIVE
C FUNCTION AT THE POINT X
C
C     F = 0.D0
C     DO 10 I = 1,N
C        T    = X(I) - I
C        F    = F + T*T
C        G(I) = 2.D0 * T
C10    CONTINUE
C      RETURN
C      END


#endif
                                                                                                                                                                                                                                                                                                                                                                                                                                                                         r8s/TNwrapper.h                                                                                     0000644 0000766 0000120 00000000557 10357562531 013354  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "TimeAlgorithms.h"
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
                                                                                                                                                 r8s/TimeAlgorithms.c                                                                                0000644 0000766 0000120 00000131336 10637045116 014351  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /**

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


                                                                                                                                                                                                                                                                                                  r8s/TimeAlgorithms.h                                                                                0000644 0000766 0000120 00000004421 10357545737 014365  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  int warnEstRoot(NODETYPE * root);
double nodeLowerBound(NODETYPE *node);
double nodeUpperBound(NODETYPE *node);
void assignArrayRatesToLL_LFLOCAL(NODETYPE *node,double p[]);

void GradientLF(double p[], double grad[]);
void GradientPL(double p[], double grad[]);
void derivTime(NODETYPE * n, double p[], double grad[],int *ixPtr);
void derivRate(NODETYPE * n, double p[], double grad[],int *ixPtr);
void derivRateLog(NODETYPE * n, double p[], double grad[],int *ixPtr);

void printnodeLike(NODETYPE *node);
double mean_rate(NODETYPE *node);
void tree2aTimeArray(NODETYPE *node, double *array);
void plotOpt(double p[],int grid,double p1low,double p1high,
	double p2low,double p2high, char *p1label,  char *p2label);
double LFcs1(NODETYPE *node,  NODETYPE *itsAncestor, double rate);
int setupParrays(NODETYPE *node, int itsAncestor);
int setupFeasibleTimes(NODETYPE * root);
void traverseSetUpFeasible(NODETYPE * node,double maxLength);
void spewtree(NODETYPE *node);

void assignArrayRatesToLL_LF(NODETYPE *root,double rate);

void assignArrayTimesToLL(NODETYPE *node,double lp[]);
/*
void assignArrayTimesToLL2(NODETYPE *root,double lp[], int includeRootFlag);
void assignArrayTimesToLL2_helper(NODETYPE * node, double lp[], int *index);
*/
void assignArrayRatesToLL2(NODETYPE *root,double lp[]);
void assignArrayRatesToLL2T(NODETYPE *root,double lp[]);
void initTreeRates(NODETYPE *root, int includeRootFlag,double rate);

double recurseLangFitch(NODETYPE *node, NODETYPE * itsAncestor, double p[]);
double LFchiSq(NODETYPE *node, double rate);
double BranchLike(double rate, double timeLength, double charLength);
double BranchLikeGamma(double rate1, double rate2,double timeLength, double charLength);
double BL(double rate, double timeLength, double charLength);
void descMinAge(NODETYPE *node, double *curMin,double *curMax);
int aFeasibleTime(NODETYPE *node,double timeAnc);
double Factorial(double x);
void setupFactLookup(void);
void setupLogFactLookup(void);
double logFact(long);
double recurseNP(NODETYPE *node, NODETYPE * itsAncestor, double p[]);
double local_rate(NODETYPE *node);
void check_if_exceed_factorial(NODETYPE *node);

/* These all have to have the same arg list */
double objNP(double p[]);
double objLangFitch(double p[]);
double objLangFitchLocal(double p[]);
double objPenLike(double p[]);
                                                                                                                                                                                                                                               r8s/TreeDrawUtils.c                                                                                 0000644 0000766 0000120 00000046413 07565204055 014165  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include <limits.h>
#include "TreeUtils.h"
#include "MacApp.h"
#include "TreeDrawUtils.h"

double 	gMaxToTipLength;	/* used in phylogram draws */
int		gRootFlag;

/***********************************************************************************/

void  MakeTreeStruct(char *TreeDescriptionPtr,struct MyRecType *theDataPtr)
										 
	/*  Takes a string tree description and creates the linked list tree structure,
	a pointer to which is then stored in the global data structure 'theData' */										 
										 
{
	NODETYPE *root;
	Str255 temp;
	long i;


	if (TreeDescriptionPtr==NULL) {
		doGenericAlert("BAD tree description");
		theDataPtr->treeStructisMade=0; /* note that it has NOT been made */
		return;
		};

	

	root=string_to_tree(TreeDescriptionPtr);

		

	if(root ==NULL) {  	/* error: tree description bad */
		doGenericAlert("BAD tree description");
		theDataPtr->treeStructisMade=0; /* note that it has NOT been made */
		return;
		}
		
	theDataPtr->theRoot=root;	/* store the pointer to the tree structure */
	theDataPtr->displayRoot=theDataPtr->theRoot; /* initially drawn tree will draw from root */
	
	theDataPtr->treeStructisMade=1; /* note that it has been made */
	SetInternalNodeCompact(theDataPtr->theRoot,theDataPtr->collapseInternalNode);
	if (theDataPtr->allowCompact)
		SetCompactNodes(theDataPtr->theRoot);
	theDataPtr->QListHasChanged=1;
	
	gRootFlag=1; 	/* So we don't count the root */
	gMaxToTipLength=calcMaxToTip(root);
	if ((gMaxToTipLength< FLT_MAX) && (gMaxToTipLength>0.0) )
		theDataPtr->lengthsAvailable=1;
		/* previous code checks to see if ALL of the initial lengths in the tree have 
		been replaced by smaller and more reasonable numbers (they were all initialized
		to a huge value).  If something went wrong on reading a tree description list,
		particularly if soome branch did not have a length, the flag will stay zero*/
	return;
}  
/**********************************************************************************/
void MyDrawString(char *s, int flag)
{
int length;
length=strlen(s);
if (length>0)
	{
	if (!flag) 			/* draw normal */
		TextMode(srcCopy);
	else				/* invert */
		TextMode(notSrcCopy);
	DrawText(s,0,length);  
	TextMode(srcOr);
	}
return;
}
/***********************************************************************************/
void 	DrawIntName(NODETYPE *node,Rect *contRectPtr)
{
	int width,length;
	length=strlen(node->taxon_name);
	width=TextWidth(node->taxon_name,0,length);
	MoveTo(contRectPtr->right-kScrollbarAdjust-width-3,contRectPtr->bottom-3);
	DrawText(node->taxon_name,0,length);
/*	MyDrawString(node->taxon_name,0);
		 For now (!) always write this as normal uninverted text */
	TextFace(FONTSTYLE);
	return;
}
/***********************************************************************************/

void Assign_XY_Tree(NODETYPE *root,  Rect *TreeRectPtr, int treeMode)
						  
/*	Assigns X and Y display coordinates to the nodes in a tree structure.  Uses
	the  integer coordinates of the upper left and lower right of the
	box in which the tree should be displayed.  The X,Y values are stored in the
	tree structure.
	treeMode is the type of tree: 0 = cladogram, 1 = phylogram , 2= chronogram*/


{
		extern double 	gMaxToTipLength;
		int N;
		double yinc,yUpLeft;  /*have to be double for accurate placement of lines */
		if(root==NULL) 
						return;
			
		(void)maxorder(root);
		(void)numdesc(root);
		gRootFlag=1;		/* needed in recursive 'calcMaxToTip'  */
		if (treeMode == 1)
			gMaxToTipLength=calcMaxToTip(root); /*for phylogram find longest path to tip*/

		gRootFlag=1;		/* needed in recursive 'assignX'  */
		assignX(root,TreeRectPtr->left,TreeRectPtr->right,
					TreeRectPtr->right-TreeRectPtr->left+1,treeMode);
		N=root->numdesc;
		if (N==1) yinc=0;
		else yinc = (TreeRectPtr->bottom-TreeRectPtr->top)/(double)(N-1);
		yUpLeft=TreeRectPtr->top;
		assignY2(root,&yUpLeft,yinc);						
		return;
}

/***********************************************************************************/

double calcMaxToTip(NODETYPE* node)

/* Calculates maximum distance from root to tip when lengths are available */

{
	double max=0.0,temp,thisLength;
	NODETYPE *child;
	if (!node) return(0.0);

	if (gRootFlag == 1)
		{
		gRootFlag = 0;
		thisLength=0.0;
		}
	else
		{
		thisLength=node->length;	/* don't add length under the root */
		}
	if (isTip(node)  ||  (node->isCompactNode)  ) 
			{
			return (thisLength);  /* treat a tip and a compact node the same way */
			}
	child=node->firstdesc;
	SIBLOOP(child) {
			temp=calcMaxToTip(child);
			if (temp > max) max = temp;
			}
	return thisLength+max;
}
/***********************************************************************************/

int assignY2(NODETYPE *node,  double *YcurPtr,  double yinc)
{
	NODETYPE *child;
	int sum=0,count=0;
	
	if (!node) return;
	if(isTip(node)  || (node->isCompactNode) )
		{
		node->Y= *YcurPtr;
		(*YcurPtr)+=yinc;
		return(node->Y);
		}

	child=node->firstdesc;
	
	SIBLOOP(child) {
		sum+=assignY2(child,YcurPtr,yinc);
		++count;
		}
	sum/=count;
	node->Y=sum;
	return(sum);
}
/***********************************************************************************/

void assignX(NODETYPE *node,  int Xleft,  int Xright, int Xwidth, int treeMode)
{
	extern double 	gMaxToTipLength;
	if (gRootFlag) {
			node->X = Xleft;
			gRootFlag = 0;
			}
	else
			switch (treeMode)
			{
			case 0:
				node->X = Xleft + (Xright - Xleft-1)/(float)(node->order + 1);
				break;
			case 1:
				node->X = Xleft + (Xwidth-1)*node->length/gMaxToTipLength;
				break;
			case 2:
				node->X = Xleft + (Xwidth-1)*(1-node->time);
				break;
			
			}


	if (node->sib) assignX(node->sib,Xleft,Xright,Xwidth,treeMode);
	
	if (node->isCompactNode) return;  

	switch (treeMode)
		{
		case 2:
		if (node->firstdesc) 
			assignX(node->firstdesc,Xleft,Xright,Xwidth,treeMode);
		break;
		default:
		if (node->firstdesc) 
			assignX(node->firstdesc,node->X,Xright,Xwidth,treeMode);
			
		}
	return;
}
/**********************************************************************************/
void DrawHigherName(struct MyRecType *DataPtr, Rect *locContentRectPtr)
{
char *r="Root";
int length;
NODETYPE *node;
node=DataPtr->displayRoot;
if (node==DataPtr->theRoot)
	{
		MoveTo(locContentRectPtr->left+2,locContentRectPtr->bottom-3);
		MyDrawString(r,0);	
	}
else
	if (*(node->taxon_name))
		{
		MoveTo(locContentRectPtr->left+2,locContentRectPtr->bottom-3);
		MyDrawString(node->taxon_name,0);
		}
return;
}

/**********************************************************************************/
void DrawTree(WindowPtr theWindow)

										 
	/*  Plots the branches
	of it in a box in a Mac Window,the current GrafPort.  Does not yet free
	up the allocated space for the tree structure.  Nor does it write the taxon
	names!*/										 
										 
{
 	int Top,Left,Bottom,Right; 
	int x,y, windowWidth,treeAreaWidth;
	Str255 temp;
	Rect TreeRect;
	Rect *DrawRectPtr;
	int treeMode,offsetRuler,j;
	NODETYPE *root;
	extern int gMax;
	double maxLength;
	struct MyRecType * dptr;
	
	dptr=(struct MyRecType *)GetWRefCon(theWindow); /* window's data */
	root=dptr->displayRoot;
	treeMode=dptr->treeMode;
/*	ContToDrawRect(&theWindow->portRect,DrawRectPtr);*/


	switch (treeMode)	/* Checks to see if necessary info is availabe for this tree view type*/
		{
		case 0:			/* Can always show cladogram */
			break;		
		case 1:			/* phylogram*/
			if (!dptr->lengthsAvailable)
				{
				doGenericAlert("Branch lengths are not currently available");
				dptr->treeMode=0;	/* Restore default */
				return;
				}
			break;
		case 2:			/* chronogram */
			if (!dptr->timesAvailable)
				{
				doGenericAlert("Branch times are not currently available");
				dptr->treeMode=0;
				return;
				}
			break;
		
		}



	gMax=0;
	
	if (treeMode == 0) 
		offsetRuler=0;
	else
		offsetRuler=dOFFSETRULER;



	(dptr->TreeRectPtr)->top=(theWindow->portRect).top + BORDER + iTextHalfHt;
	(dptr->TreeRectPtr)->left=(theWindow->portRect).left +BORDER;
	(dptr->TreeRectPtr)->bottom=(theWindow->portRect).bottom -BORDER - iTextHalfHt - offsetRuler-kScrollbarAdjust;
	(dptr->TreeRectPtr)->right=(theWindow->portRect).right -BORDER-kScrollbarAdjust;
	
	windowWidth=dptr->TreeRectPtr->right-dptr->TreeRectPtr->left+1; /* character or pixel width of window */
				
	if (root->isCompactNode) {
		root->isCompactNode=0;
			MaxTaxLength(root);
						 /* put the width of the longest taxon name in gMax */
			treeAreaWidth=max(MINWIDTH,windowWidth-gMax-6) ;
						/* tree area is the larger of MINWIDTH and the specified
						window minus the taxon areas;
						the '-6' gives space for the space between tip and name */
			gMax=min(gMax+6,windowWidth-MINWIDTH);
					/* reset gMax to be the display width ALLOWED by the difference between
					the actual size of the window and the minimum size of the tree;  this gets
					used to determine how much of the possibly length taxon name gets displayed*/		
		dptr->TreeRectPtr->right=dptr->TreeRectPtr->left+treeAreaWidth-1;
		Assign_XY_Tree(root,dptr->TreeRectPtr,treeMode);	
		MacDrawTree(root); 
		root->isCompactNode=1;
		}
	else
		{
		MaxTaxLength(root);
		treeAreaWidth=max(MINWIDTH,windowWidth-gMax-6) ;
		gMax=min(gMax+6,windowWidth-MINWIDTH);
		dptr->TreeRectPtr->right=dptr->TreeRectPtr->left+treeAreaWidth-1;
		Assign_XY_Tree(root,dptr->TreeRectPtr,treeMode);	
		MacDrawTree(root);
		} 
		/* all this is to allow drawing of the inside of a compact node, which would
		be precluded by the fact that the root of that subtree has its flag set to compact,
		and this would force AssignXY to bail out immediately */

	if (treeMode > 0)
		{
		gRootFlag=1; 	/* So we don't count the root */
		maxLength=calcMaxToTip(root);
		drawRuler(*(dptr->TreeRectPtr),treeAreaWidth,treeMode,maxLength);
		}

	return;


}  
/***********************************************************************************/
#define dINSET_RULER 15
#define dTICKMARK_LENGTH 3
void drawRuler(Rect TreeRect,int width, int treeMode, double maxLength)
{
double length=1.0,interval;	/* factor controls roughly how long the phylogram bar is */
Str255 numStr;
int	leftRuler, rightRuler,yRuler,i,x;
switch (treeMode)
	{
	case 1:
		if (maxLength > 1.0)
			{
			while ( length< 0.1* maxLength)
				length *= 10.0;
			}
		x=length;
		length=length*width/maxLength;
		leftRuler=TreeRect.left;
		rightRuler=TreeRect.left+length;
		yRuler=TreeRect.bottom+dINSET_RULER;
		NumToString((long)x,numStr);
		MoveTo(leftRuler,yRuler);
		LineTo(rightRuler,yRuler);
		MoveTo(rightRuler + 2, yRuler);
		DrawString(numStr);
		interval=length/10.;
		for (i=0;i<=10;i++)
			{
			x=leftRuler+i*interval;
			MoveTo(x,yRuler); 
			LineTo(x, yRuler-dTICKMARK_LENGTH);
			}
		
		
		break;
	case 2:
		leftRuler=TreeRect.left;
		rightRuler=TreeRect.right;
		yRuler=TreeRect.bottom+dINSET_RULER;
		MoveTo(leftRuler,yRuler);
		LineTo(rightRuler,yRuler);
		interval=width/10.;
		for (i=0;i<=10;i++)
			{
			x=leftRuler+i*interval;
			MoveTo(x,yRuler); 
			LineTo(x, yRuler-dTICKMARK_LENGTH);
			}
		MoveTo(leftRuler-2, yRuler+FONTSIZE);
		DrawString("\p1.0");
		MoveTo(rightRuler-2, yRuler+FONTSIZE);
		DrawString("\p0.0");


		break;
	}

}
/***********************************************************************************/
void MacDrawTree(NODETYPE *node)
{
	NODETYPE *child;
	Rect theRect;
	int length,invertFlag;
	char *asterisk="*";
	char s[20];
	if (!node) return;
	if (!isTip(node)  ) 
		{
			SetRect(&theRect,node->X-NODEBOX,node->Y-NODEBOX,node->X+NODEBOX,node->Y+NODEBOX);
			FrameRect(&theRect);
			PaintRect(&theRect);
			if (node->nodeIsConstrained)
				{
				MoveTo(node->X+4,node->Y+3);
				LineTo(node->X+4,node->Y-3);
				MoveTo(node->X-5,node->Y+3);
				LineTo(node->X-5,node->Y-3);
				}
			MoveTo(node->X+3,node->Y+5);
			if  (*(node->taxon_name))  
					{
			/*		MyDrawString(asterisk,node->isQueryNode);*/ /* asterisk */
					MyDrawString(node->taxon_name,node->isQueryNode);
					
					TextFace(FONTSTYLE);
					}
		}
	if (!(node->isCompactNode) )  /* only draw branches to descendants if node is not compact */
		{
			child=node->firstdesc;
			SIBLOOP(child) 
				{
				MoveTo(node->X,child->Y);
				LineTo(child->X,child->Y);
				
				Move(TNOFFSETX,TNOFFSETY);  /*  offset from tip of branch*/
				if (isTip(child))
					{
					length=strlen(child->taxon_name);
					if (TextWidth(child->taxon_name,0,length) <= gMax)
						MyDrawString(child->taxon_name,child->isQueryNode);  /* as long as name isnt too 
																	long...*/
					TextFace(FONTSTYLE);
					}	
				MoveTo(node->X,node->Y);
				LineTo(node->X,child->Y);
				MacDrawTree(child);
				}
		}
	return;
}
/***********************************************************************************/
int maxorder(NODETYPE *node)
{
	 int max,temp;
	NODETYPE *child;
	if (!node) return(-1);
	if (isTip(node)  ||  (node->isCompactNode)  ) {node->order=0; return (0);}
			/* treat a tip and a compact node the same way */
	max=0;
	child=node->firstdesc;
	SIBLOOP(child) {
			temp=maxorder(child);
			if (temp > max) max = temp;
			}
	node->order=max+1;
	return (max+1);
}
/****************************************************/
void MaxTaxLength(NODETYPE *node)

	/* Finds the length of the longest string contained in the tree structure */
{
	NODETYPE *child;
	extern int gMax;
	int temp, length, max=0;
	if (!node) return;

	if (isTip(node))
		{
		length=strlen(node->taxon_name);
		temp=(short)TextWidth(node->taxon_name,0,(short)length);
		if (temp>gMax) gMax=temp;
		}


	if (node->isCompactNode) return;
	child=node->firstdesc;
	SIBLOOP(child) 
		{
		MaxTaxLength(child);
		}
	return;


}
/****************************************************/
NODETYPE *SearchTreeNodes(Point localPt, int radius, NODETYPE *root)

	/* Checks to see if local mouse-coordinates in 'P' are within some distance 'radius'
	from a node in the window's tree structure.  If so the function returns a pointer to
	the node; if not it returns NULL.  */
{
	NODETYPE	*foundNode;

	foundNode=TraverseScanPt(localPt, radius, root);
	return(foundNode);
}
/****************************************************/
void SetCompactNodes(NODETYPE *rootPtr)
{
extern int kMAXSCRTAXA;
int ix=0;
do	{
	numdesc(rootPtr);
	SetOneNode(rootPtr);
	++ix;
	}
	while (rootPtr->numdesc > kMAXSCRTAXA && ix< kMAXSCRTAXA);
return;
}

int SetOneNode(NODETYPE *node)
{
if (node->isCompactNode) return (0);
if ( (node->sib) && SetOneNode(node->sib) ) return (1);
if ( (node->firstdesc) && SetOneNode(node->firstdesc) ) return (1);
if (node->numdesc >kMAXSCRTAXA)
	{
	node->isCompactNode=1;
	return (1);
	}
else
	return (0);
}

void ClearCompactNodes(NODETYPE *node)
{
if (!node) return;
node->isCompactNode=0;
if (node->sib) ClearCompactNodes(node->sib);
if (node->firstdesc) ClearCompactNodes(node->firstdesc);
return;
}

void SetInternalNodeCompact(NODETYPE *node,int mode)
{
/* if the collapse internal node mode flag is set, traverse tree and set to compact
	node any node that is (a) internal, and (b) has a taxon name associated */

if (!node) return;
if 	( 
	(*(node->taxon_name))  
	&& (!isTip(node)) 
	&& mode
	)
	node->isCompactNode=1;
if (node->sib) SetInternalNodeCompact(node->sib,mode);
if (node->firstdesc) SetInternalNodeCompact(node->firstdesc,mode);
return;
}



/***********************************************************************************/

NODETYPE *TraverseScanPt(Point localPt, int radius, NODETYPE *node)

{
	NODETYPE	*foundNode;
	
	if (!node) return (NULL);
	if (  (node->X - radius) < localPt.h && (node->X + radius) > localPt.h &&
	    	  (node->Y - radius) < localPt.v && (node->Y + radius) > localPt.v ) return(node);
	else 
		{
			if (node->sib) 
				{
				foundNode=TraverseScanPt(localPt,radius, node->sib);
				if(foundNode) return (foundNode);
				}
			if (node->firstdesc) 
				{
				foundNode=TraverseScanPt(localPt,radius,node->firstdesc);
				if(foundNode) return (foundNode);
				}
			return (NULL);
		}
}
/***********************************************************************************/

NODETYPE *QTreeNodes(Point localPt, int radius, NODETYPE *node)

/* scans through the tree and returns the node corresponding to either a node that is clicked
or the node next to the taxon name that is clicked */

{
	NODETYPE *child;
	NODETYPE	*foundNode;
	Rect		taxRect;
	int			length,width,height,starth,startv;
	char		*asterisk="*";

/* check for click in node box */
	if (!node  || node->isCompactNode) return (NULL);
	if (  (node->X - radius) < localPt.h && (node->X + radius) > localPt.h &&
	    	  (node->Y - radius) < localPt.v && (node->Y + radius) > localPt.v ) return(node);

/* check for click in taxon name area, only if there is a taxon name for that node
	and ONLY for terminal taxa */
	
	if (*(node->taxon_name)  && (isTip(node)))
		{
		width=TextWidth(node->taxon_name,0,strlen(node->taxon_name));
		height=FONTSIZE/2;
		SetRect(&taxRect,node->X +4, node->Y-height,node->X +4+width, node->Y+height);
		if (PtInRect(localPt,&taxRect)) 
			{
			return (node);
		/*	doQToggleTaxon(node);
			return((NODETYPE *)-1);	*/	/* indicates a taxon was toggled */
			}
		}

/* check for click in asterisk that corresponds to higher node name, but only
for internal node that has a name */

	starth=node->X+3;
	startv=node->Y+5;
	if (*(node->taxon_name)  && (!isTip(node)))
		if (  (node->X +3) < localPt.h && (node->X +9) > localPt.h &&
	    	  (node->Y - 3) < localPt.v && (node->Y + 3) > localPt.v ) 
	    	 {
			 FLAGFLIP(node->isQueryNode);
	    	 MoveTo(starth,startv);
			 MyDrawString(asterisk,node->isQueryNode); /* asterisk */
			 return((NODETYPE *)-1);	/* just toggle the status and let 'DrawIntName handle it */
	    	 }

 
	child=node->firstdesc;
	SIBLOOP(child) 
			{
			foundNode=QTreeNodes(localPt,radius, child);
			if(foundNode) return (foundNode);
			}
	return (NULL);	/* Haven't found either a node or a taxon name so far */

}
void doQToggleTaxon(NODETYPE *node)
{
FontInfo	fi;
int 		height,width,starth,startv;
Rect		taxRect;

/*GetFontInfo( &fi );*/

starth=node->X+TNOFFSETX;
startv=node->Y+TNOFFSETY;
/*width=TextWidth(node->taxon_name,0,length);
SetRect(&taxRect,starth, startv-fi.ascent,starth+width, startv+fi.descent);
EraseRect(&taxRect);*/

MoveTo( starth,startv );

FLAGFLIP(node->isQueryNode);
MyDrawString(node->taxon_name,node->isQueryNode);  
TextFace(FONTSTYLE);	

return;

}
void QToggleClade(NODETYPE *node)

/* sets or resets the query status of all of the tips descended from node.
DOES NOT alter the internal nodes.  Decision on whether to set or reset
is based on node->toggleDesc flag. Thus, clicking on an internal node will
either set all descendants or reset all descendants depending on the history of clicks
on that internal node*/
{
	int flag;
	if (node->toggleDesc) 
		flag=0; 
	else 
		flag=1;
	QSetAllNodes(node,flag);
	node->toggleDesc=flag;
	return;
}
void QSetAllNodes(NODETYPE *node,int flag)
{

	NODETYPE *child;
	if (!node) return;
	if (isTip(node))
		node->isQueryNode=flag;		
	child=node->firstdesc;
	SIBLOOP(child) {
			QSetAllNodes(child,flag);
			}
	return;
}

/*********************************************************/
NODETYPE *ScanInternNodes(Point localPt, int radius, NODETYPE *node)


{
	NODETYPE	*foundNode;
	int			length,width,height;

/* check for click in asterisk */
	if (!node) return (NULL);
	if (!isTip(node))
		if (  (node->X +6- radius) < localPt.h && (node->X +6+ radius) > localPt.h &&
	    	  (node->Y - radius) < localPt.v && (node->Y + radius) > localPt.v ) 
	    	 {
	    	 return(node);
	    	 }



 
	if (node->sib) 
		{
		foundNode=ScanInternNodes(localPt,radius, node->sib);
		if(foundNode) return (foundNode);
		}
	if (node->firstdesc) 
		{
		foundNode=ScanInternNodes(localPt,radius,node->firstdesc);
		if(foundNode) return (foundNode);
		}
	return (NULL);

}

                                                                                                                                                                                                                                                     r8s/TreeDrawUtils.h                                                                                 0000644 0000766 0000120 00000002614 07565204055 014165  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  double 		calcMaxToTip(NODETYPE* root);

int 		maxorder(NODETYPE *),
			assignY2(NODETYPE *node,  double *YcurPtr,  double yinc),
			SetOneNode(NODETYPE *),
			doQWriteFile(struct MyRecType *DataPtr)
			;

NODETYPE	
		 	*TraverseScanPt(Point localPt, int radius, NODETYPE *node),
		 	*QTreeNodes(Point localPt, int radius, NODETYPE *node),
		 	*ScanInternNodes(Point localPt, int radius, NODETYPE *node),
		 	*SearchTreeNodes(Point globalPt, int radius, NODETYPE *root);

void		Tprint(NODETYPE *),
			MyDrawString(char *s, int flag),
			DrawIntName(NODETYPE *node,Rect *contRectPtr),
			DrawHigherName(struct MyRecType *DataPtr, Rect *locContentRectPtr),
			SetInternalNodeCompact(NODETYPE *node,int mode),
			ClearCompactNodes(NODETYPE *node),
			QSetAllNodes(NODETYPE *node,int flag),
			QToggleClade(NODETYPE *node),
			doQToggleTaxon(NODETYPE *node),
			assignX(NODETYPE *,  int,  int,int width,int treeMode),
			assignY(NODETYPE *,  int,  int),
			Assign_XY_Tree(NODETYPE *root,  Rect *TreeRectPtr, int treeMode),
			MagnifyTree(NODETYPE *node, int xLeft, int yLeft, double xFactor, double yFactor),
			MaxTaxLength(NODETYPE *node),
	 		MacDrawTree(NODETYPE *node),
			DrawTree(WindowPtr w),

			SetCompactNodes(NODETYPE *),
			MakeTreeStruct(char *TreeDescriptionPtr,struct MyRecType *theDataPtr),

			Dispose(struct MyRecType *theDataPtr);

void drawRuler(Rect TreeRect,int width, int treeMode, double maxLength);
                                                                                                                    r8s/TreeSim.c                                                                                       0000644 0000766 0000024 00000111072 11354151762 013004  0                                                                                                    ustar   sandermj                        staff                                                                                                                                                                                                                  #define BIG_VAL 1e20 
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
 times[] is 0-offset and is ordered from recent to past on the n-2 internal nodes, opposite their sense...
 d = is the difference between ordered node ages
 T = age of root
 */
{
long i,j,k,ll;
double A,B,C,gamma,d;
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
                                                                                                                                                                                                                                                                                                                                                                                                                                                                      r8s/TreeSim.h                                                                                       0000644 0000766 0000024 00000003314 11220055016 012773  0                                                                                                    ustar   sandermj                        staff                                                                                                                                                                                                                  #include "TreeUtils.h"
double YuleLike(NODETYPE* root, int stemFlag, double *speciation);
NODETYPE* SisterGroupYule(double T, double spec1, double spec2, double *Ntips1, double *Ntips2, double *sum_durations);
NODETYPE* RY_1997(long n_taxa, double T, double speciation, double extinction, double sampling);
double PH_gamma (long n, double *times, double T);

void RandomBranches(TREE Tree, long nNodes, NODETYPE ** nodeArray, long nMark, NODETYPE ** markedNodes, int  withReplace);
void markRandomNodes(TREE Tree, long nMark, NODETYPE ** markedNodes);
int BDDiversity(long n_taxa, double spec_rate, double extinct_rate,
		double char_rate, int interval);
NODETYPE* BDTree(long n_taxa, double spec_rate, double extinct_rate,
		double char_rate);
NODETYPE * coalesce_nodes(NODETYPE * node1, NODETYPE * node2);
void set_branch_rates(NODETYPE *node, double curRate, double rateChangeAmt,
	double minRate, double maxRate,double transitionProb, int gradual,
	int model);
double angle(double *vec1, double *vec2, int arraySize);
double euclid_distance(double *vec1, double *vec2, int arraySize);
double correlation(double x[], double y[], unsigned long n);
double tree_auto_correlation(NODETYPE * root);
void set_branch_lengths(NODETYPE *node, int infinite);
void set_est_rates(NODETYPE *node,double b, double c, int d);
double BD_Like(double params[]);
NODETYPE* BDTreeForward(double T, double spec_rate, double extinct_rate,
		double char_rate);
double Yule_forward(double rate,  double T, double *sum_durations,int stemFlag);
NODETYPE* YuleTreeForward(double T, double spec_rate, double *Ntips, double *sum_durations,int stemFlag);
NODETYPE* BDback(long n_taxa,double specRate,int normalFlag);
NODETYPE* Yule_C(long n_taxa, double speciation);
                                                                                                                                                                                                                                                                                                                    r8s/TreeUtils.c                                                                                     0000644 0000766 0000120 00000215735 11650316467 013355  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /******************        Tree Utility module    ********************

					
	NOTE: The tree/node data structure is defined in TreeUtils.h
					
**********************************************************************/
#if 0
#define isEqual(a,b)		(!strcmp((a),(b)))
#endif

#include <limits.h>
#include "ObjFunc.h"
#include "TreeUtils.h"
#include "MyUtilities.h"
#include "NRCvectorUtils.h"
#include "memory.h"
#include "DrawTree.h"
#include "structures.h"
#include <stdlib.h>
#include <math.h>
#include "TimeAlgorithms.h"
#include "moment.h"
#include "DistrFuncs.h"

/* ...Global variables throughout module */

char 	*gStringptr;
int	gId, gFlag;
int	gCount;

double sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)
void Flip2(NODETYPE *a);
void printnode(NODETYPE *node);
static void setModel(NODE n,int model);
static void preOrderIntArg(NODETYPE *node,void (*func)(NODE node, int iarg),int iarg);
static double SThelper(NODETYPE * node,double factor);
static void collapse_zero_2trees(NODE node1, NODE node2);
static double cvCS(NODETYPE * node);
static double cvCST(NODETYPE * node);
static double cvSQNP(NODE node);
static double cvSELF(NODE n,double rate);
static double moveNER(NODE node);
static double setNER(NODE node);
static double cvSQET(NODETYPE * node);
static double zeroER(NODETYPE *node);
static void insertNode(NODETYPE *node,  NODETYPE* anc);
static void updateOneSubtree(NODETYPE *subRoot);
NODETYPE * prevSib(NODETYPE* node);

/****************************************************************/

// returns a 1-off variance-covariance matrix based on the ultrametric distances on tree t including only branches with model ID matching model
// note that the order of indices and terminals is determined by the call to 'TreeToTaxaPtrList'

double ** tree2VCV(TREE t, int model)
{
	PtrList lnode;
	long i, j;
	NODETYPE *root, *node;
	NODE a,b,c;
	PtrList nodeList;
	long lengthList,n;
	double ** vcv,T;
	root=t->root;
	nodeList=pNewList();
	TreeToTaxaPtrList(root,nodeList);
	n=pLengthList(nodeList); // the number of taxa!
	vcv = matrix(1,n,1,n);
printf("\nVariance-covariance matrix for model %i\n",model);
	for (i=1;i<=n;i++)
			{
			a=(NODE)(pListgetkthNode(nodeList, i)->item);
			printf("%s\t",a->taxon_name);
			for (j=1;j<=n;j++)
					{
					b=(NODE)(pListgetkthNode(nodeList, j)->item);
					c=mrca(a,b);
					T=pathLengthTimeModel(c,root,model);
					printf("%f\t",T);
					vcv[i][j]=T;
					}
			printf("\n");
			}
	freepList(nodeList);
	return vcv;						
}
/*****************************************************************************************************/

NODETYPE * nextRndNode(long nNodes, NODETYPE ** nodeArray)

/* return a randomly selected node which has not yet been marked. Dumbly keeps looking for available nodes,
	and will do so forever if they've all been sampled. Very brute force routine--stupid even! Better check in calling routine ! */

{
NODETYPE * node;
long rn;
do 
	{
	rn=rndLong(nNodes);
	node = nodeArray[rn];
	} while (isNodeMarked(node));
// markNode(node); TEMPORARY, hope this works
return node;
}

/*****************************************************************************************************/
void markNode( NODETYPE  * n)
{
n->isQueryNode=1;
return;
}
void unMarkNode( NODETYPE  * n)
{
n->isQueryNode=0;
return;
}
int isNodeMarked( NODETYPE  *n)
{
return n->isQueryNode;
}
void unMarkTree(TREE T)
{
preOrderVoid(T->root,unMarkNode);
return;
}

void setLocalModel(NODETYPE *n,int model,int stemFlag)
{
int saveModel;
saveModel=n->modelID;
	
preOrderIntArg(n,setModel,model);
if (/* !isTip(n) &&*/ stemFlag==0)
	n->modelID=saveModel; /* for interior nodes, if we DONT want the stem rate assigned we better have saved it,
					because the recursion will assign it */
return;
}

static void setModel(NODE n,int model)
{
n->modelID=model;
return;
}


NODETYPE * sister(NODETYPE * n)
{
if (n->anc == NULL)
	return NULL;
if (n->anc->firstdesc==n)
	return n->sib;
else
	return prevSib(n);


}

void copyLengthsTree2Tree(NODETYPE * node1,NODETYPE * node2)
{
NODETYPE * child1, *child2;
node1->nodeReal=node2->length;
if (!isTip(node1))
	{
	child1=node1->firstdesc;
	child2=node2->firstdesc;
	for (;child1;child1=child1->sib,child2=child2->sib)
		copyLengthsTree2Tree(child1,child2);
	}
return;
}

/***********************************************************************************/
void preOrderVoid(NODETYPE *node,void (*f)(NODETYPE *))
{
	NODETYPE *child;
	(*f)(node);
	if (!isTip(node))
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			preOrderVoid(child,f);
		}
	return ;	
}
double preOrder(NODETYPE *node,double (*f)(NODETYPE *))
{
	double sum=0;
	NODETYPE *child;
	sum+=(*f)(node);
	if (!isTip(node))
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			sum += preOrder(child,f);
		}
	return (sum);	
}
double preOrderArg(NODETYPE *node,double (*func)(NODE node, double farg),double farg)
{
/* */

	double sum=0;
	NODETYPE *child;
	sum+=(*func)(node,farg);
	if (!isTip(node))
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			sum += preOrderArg(child,func,farg);
		}
	return (sum);	
}
static void preOrderIntArg(NODETYPE *node,void (*func)(NODE node, int iarg),int iarg)
{
/* */

	NODETYPE *child;
	(*func)(node,iarg);
	if (!isTip(node))
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			preOrderIntArg(child,func,iarg);
		}
	return;	
}
/***********************************************************************************/
void setNodeEstRate(NODE node)

/* for PENLIKET method, this is a kludge which ends up calculating the mean branch rates
	and storing them in the usual place, so ratogram draws are possible. Also stores
	node rate estimates in node->nodeEstRate field */

{
(void)preOrder(node,moveNER);
(void)preOrder(node,setNER);
return;
}
static double moveNER(NODE node)
{
node->nodeEstRate=node->estRate;
return 0.0;
}
static double setNER(NODE node)
{
if (!isRoot(node))
	node->estRate=(node->nodeEstRate+node->anc->nodeEstRate)/2.0;
return 0.0;
}

/***********************************************************************************/


void zeroEstRate(NODETYPE *node)
{
(void)preOrder(node,zeroER);
return;
}
static double zeroER(NODETYPE *node)
{
node->estRate=0.0;
node->nodeEstRate=0.0;
return 0.0;
} 

/***********************************************************************************/

double LFuncons(NODETYPE *node)
{

return preOrder(node,LFuncons1);

}

double LFuncons1(NODETYPE *node)
{
double expected, chiSq=0.0;
if (!isRoot(node))
	{
	expected=node->estRate*(node->anc->time-node->time);
  	if (fabs(expected) > 0.0001)  
		{
		chiSq= SQR(node->length - expected)/expected;
		}
/*printf("SAME expected=%f nodelength=%f node->nodeReal=%f chiSq=%f\n",expected,node->length,node->nodeReal,chiSq);*/
	}
return chiSq;
}
double LFunconsT(NODETYPE *node)
{

return preOrder(node,LFuncons1T);

}

double LFuncons1T(NODETYPE *node)
{
double expected, chiSq=0.0,r;
if (!isRoot(node))
	{
	r=(node->estRate+node->anc->estRate)/2;
	expected=r*(node->anc->time-node->time);
  	if (fabs(expected) > 0.0001)  
		{
		chiSq= SQR(node->length - expected)/expected;
		}
/*printf("SAME expected=%f nodelength=%f node->nodeReal=%f chiSq=%f\n",expected,node->length,node->nodeReal,chiSq);*/
	}
return chiSq;
}

/***********************************************************************************/
double cvSquareErrorBranch(TREE t, NODE node,int method,double *chiSq)

/* 

For all methods, uses estimated rates and times to calc the prediction error of branch subtending node n.

Prediction error calculated as follows: 

LF: uses the overall tree-wide estimated rate

NP and PL: uses the estimated rate of the branch BELOW the removed branch, OR
		if the pruned branch's ancestor is the root, uses the estimated rate of
		the pruned branch's sister branch.

For 0-length branches, just returns 0 for ChiSq if expected is 0.

 */

{
NODETYPE *n;
double sq=0.0;
double d,expected;
if (isRoot(node))return 0.0;
*chiSq=0.0;
d=node->anc->time-node->time;
switch (method)
	{
	case LaF:
		expected=t->estRate*d;
		sq=SQR(expected - node->length);
		if (expected == 0.0)
			*chiSq=0.0;
		else
			*chiSq=sq/expected;
		break;
	case NP:
		if (!isRoot(node->anc)) 
			{
			expected=node->anc->length/(node->anc->anc->time - node->anc->time)*d;
			sq=SQR(expected - node->length);
			if (expected == 0.0)
				*chiSq=0.0;
			else
				*chiSq=sq/expected;
			}
		else
			{
			n=sister(node);
			expected=n->length/(n->anc->time - n->time)*d;
			sq=SQR(expected - node->length);
			if (expected == 0.0)
				*chiSq=0.0;
			else
				*chiSq=sq/expected;

			/*!!!!!!!!  Put in code to handle case where rate is at the root */
			}
		break;


		break;
	case PENLIKE:
		if (!isRoot(node->anc)) 
			{
			expected=node->anc->estRate*d;
			sq=SQR(expected - node->length);
			if (expected == 0.0)
				*chiSq=0.0;
			else
				*chiSq=sq/expected;
/*
printf("%s rate=%f expected=%f dur=%f length=%f SqEr=%f chiSq=%f\n",node->taxon_name,node->anc->estRate,expected,d,node->length,sq,*chiSq);
*/
			}
		else
			{
			expected=sister(node)->estRate*d;
			sq=SQR(expected - node->length);
			if (expected == 0.0)
				*chiSq=0.0;
			else
				*chiSq=sq/expected;
			}
		break;
	case PENLIKET: /* NOT UPDATED YET TO HANDLE ROOT ISSUES */
		expected=node->anc->estRate*d;
		sq=SQR(expected - node->length)/expected;
		/* IMPORTANT: For this method, we use this node's ancestor for the rate. The node
			itself was a tip, and was deleted from the analysis ! */
		break;
	}
return sq;

}
#if 0
/***********************************************************************************/
double cvSquareError(TREE t, int method)

/* For all methods, uses estimated rates and times to calc the prediction error of branch
	lengths that are stored in the variable node->nodeReal */

{
double sq=0.0;
extern int gNComp;
switch (method)
	{
	case LaF:
		(void)preOrderArg(t->root,cvSELF,t->estRate/(gNComp-1));
		sq=preOrder(t->root,cvCS);
		break;
	case NP:
		sq=preOrder(t->root,cvSQNP);
		break;
	case PENLIKE:
		sq=preOrder(t->root,cvCS);
		break;
	case PENLIKET:
		sq=preOrder(t->root,cvSQET);
		break;
	}
printf("FINAL CALC ON SQ ERR= %f\n",sq);
return sq;

}
/***********************************************************************************/


static double cvSQNP(NODE node)
{
double expected, Sq=0.0,d,estRate;
extern int gNComp;
if (!isRoot(node))
	{
	expected=node->length/(gNComp-1); /* seems trivial, but correct for the NP model */
  	if (fabs(expected) > 0.0001)  
		{
		Sq= SQR(node->nodeReal - expected);
		}
	}
return Sq;
}

/***********************************************************************************/

static double cvSELF(NODE n,double rate)
{
n->estRate=rate;
return 0.0;
}

/***********************************************************************************/

static double cvCS(NODETYPE * node)
{
double expected, chiSq=0.0,d;
extern int gNComp;
if (!isRoot(node))
	{
	d=node->anc->time-node->time;
	expected=node->estRate*d/(gNComp-1);
  	if (fabs(expected) > 0.0001)  
		{
		chiSq= SQR(node->nodeReal - expected);
		}

printf("rate=%f expected=%f dur=%f length=%f node->nodeReal=%f SqEr=%f\n",node->estRate,expected,d,node->length,node->nodeReal,chiSq);


	}
return chiSq;
}
/***********************************************************************************/


double cvCST(NODETYPE * node)
{
double expected, chiSq=0.0,d,r;
if (!isRoot(node))
	{
	d=node->anc->time-node->time;
	r=(node->estRate+node->anc->estRate)/2;
	expected=r*d;
  	if (fabs(expected) > 0.0001)  
		{
		chiSq= SQR(node->nodeReal - expected) /expected ;
		}

/*printf("rate=%f expected=%f dur=%f length=%f node->nodeReal=%f chiSq=%f\n",node->estRate,expected,d,node->length,node->nodeReal,chiSq);
*/

	}
return chiSq;
}
/***********************************************************************************/

static double cvSQET(NODETYPE * node)
{
double expected, sq=0.0,d,r;
extern int gNComp;
if (!isRoot(node))
	{
	d=node->anc->time-node->time;
	r=(node->estRate+node->anc->estRate)/2;
	expected=r*d/(gNComp-1);
  	if (fabs(expected) > 0.0001)  
		{
		sq= SQR(node->nodeReal - expected);
		}
	}
return sq;
}
#endif

/********************************************/

void unSetConstraints(NODETYPE * node)


{
	NODETYPE *child;

	node->nodeIsConstrainedMax=0;
	node->nodeIsConstrainedMin=0;
	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			unSetConstraints(child);
		}
	return;
}


double unFixNodeAge(NODETYPE *node) 

	/* only allowed on internals for convenience*/
{
if (!isTip(node))
	node->free=1;
return 0.0;
} 


double fixNodeAge(NODETYPE *node)
{
node->free=0;
return 0.0;
} 

int numFixedNodes(NODETYPE * node)

/* Returns number of (internal) descendents of the clade at node with fixed ages */

{
	NODETYPE *child;
	int numFix=0;
	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			numFix+=numFixedNodes(child);
		if (isFixed(node))
			++numFix;
		}
	return numFix;
}
int numConstrainedNodes(NODETYPE * node)

/* Returns number of (internal) descendents of the clade at node with constrained (not fixed) ages */

{
	NODETYPE *child;
	int numCon=0;
	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			numCon+=numConstrainedNodes(child);
		if (isConstrained(node))
			++numCon;
		}
	return numCon;
}

void setupFixedNodeArray(NODETYPE * node, NODE nodeArray[], int *curIndex)

/* Populates a node array with the nodes that are currently constrained. */

{
	NODETYPE *child;
	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			setupFixedNodeArray(child,nodeArray,curIndex);
		if (!isFree(node))
			nodeArray[(*curIndex)++]=node;
		}
	return;
}

void setupConstrainedNodeArray(NODETYPE * node, NODE nodeArray[], int *curIndex)

/* Populates a node array with the nodes that are currently constrained. */

{
	NODETYPE *child;
	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			setupConstrainedNodeArray(child,nodeArray,curIndex);
		if (isConstrained(node))
			nodeArray[(*curIndex)++]=node;
		}
	return;
}

int maxAgePresent(NODETYPE * node)

/* Returns a 1 if any descendents of the clade at node have a max age constraint set */

{
	NODETYPE *child;

	if (node->nodeIsConstrainedMax)
		return 1;

	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			if(maxAgePresent(child))
				return 1;
		}
	return 0;
}
int constraintsPresent(NODETYPE * node)

/* returns a 1 if any descendents of the clade at node have time constraints set */

{
	NODETYPE *child;

	if (node->nodeIsConstrainedMax || node->nodeIsConstrainedMin)
		return 1;

	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			if(constraintsPresent(child))
				return 1;
		}
	return 0;
}
int tipsDifferentAges(NODETYPE *node)

/* determines whether all the tips are the same age: 1=different, 0=same */

{
	static int first=1;
	static double save;
	NODETYPE *child;
	if (isTip(node)) 
		{
		if (first)
			{
			save=node->time;
			first=0;
			}
		else
			if (save != node->time)
				return 1;
		}
	child=node->firstdesc;
	SIBLOOP(child) 
		if(tipsDifferentAges(child))
			return 1;
	return 0;
}

/***********************************/
void updateSubtrees(NODETYPE *srcNode)

/* Copy all the current times on the source tree to the various subtree data structures 
	defined from it */

{
	NODETYPE *child,*subRoot;
	if (!isTip(srcNode) && !isRoot(srcNode))
			{
			subRoot=srcNode->nodePtr;
			if (subRoot)
				updateOneSubtree(subRoot);	/* update the subtree for this source tree node */
			}
	child=srcNode->firstdesc;
	SIBLOOP(child)
			updateSubtrees(child);
	return;
}

static void updateOneSubtree(NODETYPE *subNode)
{
	NODETYPE *child,*srcNode;
	subNode->time=(subNode->nodePtr)->time;		/* subNode->nodePtr maintains a pointer to the corresponding source tree node */
	child=subNode->firstdesc;
	SIBLOOP(child)
		updateOneSubtree(child);
	return;


}


/***********************************/

void setupSubtrees(NODETYPE * srcNode)

/* for each internal node of tree rooted at srcNode, setup a subtree and store a pointer to this
subtree in the internal node's nodePtr location */

{
	NODETYPE *child;
		{
		if (!isTip(srcNode) && !isRoot(srcNode))
			srcNode->nodePtr=createSubtree(srcNode,0);
		child=srcNode->firstdesc;
		SIBLOOP(child)
			setupSubtrees(child);
		}
	return;
}

/***********************************/

NODETYPE *createSubtree(NODETYPE *srcNode, int SubtreeSize)

/* Returns a pointer to a newly allocated tree, which is created by copying a 
subtree from tree srcRoot. Copies pertinent time information from source nodes.
Each node also stores a pointer to the node on the source tree from whence it came.
This permits rapid updating of information about time.

At the moment this routine ignores SubtreeSize, and makes a subtree from three branches
surrounding srcNode.

*/
{
NODETYPE *root, *cnode,*node,*child;
if (!isTip(srcNode) && !isRoot(srcNode)) /* only allow subtrees from internal nodes! */
	{
	root=newnode();
	copyNodeInfo(srcNode->anc,root);
	root->nodePtr=srcNode->anc;
	cnode=newnode();
	AddChild(root,cnode);
	copyNodeInfo(srcNode,cnode);
	cnode->nodePtr=srcNode;
	child=srcNode->firstdesc;
	SIBLOOP(child)  /* for each child of the source node, create a child on the copied tree */
		{
		node=newnode();
		AddChild(cnode,node);
		copyNodeInfo(child,node);
		node->nodePtr=child;		
		}
	return root;
	}
else
	return NULL;
}

/**********************************/
NODETYPE*  copyTree(NODETYPE* a)
// returns a node that is either a tip, or the root of a properly formatted tree, but its ancestor and sibs are undefined
{
NODETYPE* child,*first,*newfirst,*newn,*n,*prev;
newn = newnode();
copyNodeInfo(a,newn);
if(!isTip(a))
	{
	first=a->firstdesc;
	newfirst=copyTree(first);
	newn->firstdesc=newfirst;
	newfirst->anc = newn;
	prev=newfirst;
	child=first->sib; // start loop with the second sib in the sib list...
	SIBLOOP(child)
		{
		n = copyTree(child);
		prev->sib=n;
		prev=prev->sib;
		n->anc = newn;
		}
	}
return  newn;
}

/**********************************/



void AddChild(NODETYPE * parent, NODETYPE * theChild)
        {
	NODETYPE *aChild;
	if (parent->firstdesc)
	    {
	    aChild=parent->firstdesc;
	    if (aChild)
		    {
		    while(aChild->sib)
			    aChild=aChild->sib;
		    aChild->sib=theChild;
		    }
	    }
	else
	    parent->firstdesc=theChild;
	theChild->anc=parent;
	theChild->sib=NULL;
        return;
        }
void RemoveTaxonLvAnc(NODETYPE * rmTaxon)

/* remove a tip or clade, but leave its ancestor node in place */

{
NODETYPE * prev;
prev=prevSib(rmTaxon);
if(prev)	/* either rmTaxon is the firstdesc or its got a prev sib */
	prev->sib=rmTaxon->sib;
else
	rmTaxon->anc->firstdesc=rmTaxon->sib;
rmTaxon->anc=NULL;
rmTaxon->sib=NULL;
}


NODETYPE * RemoveTaxon(TREE T,NODETYPE * rmTaxon)

/* Removes a taxon, or clade, including the stem lineage
 * Does not remove the node below the stem lineage if that becomes of degree two
 * Does not deallocate memory for the subtree that is deleted, or change links on that subtree.
 * Won't allow removal of the root node
 * If the node is one of only two children of the root, the root is removed as well.
 * RETURNS A POINTER TO THE PRUNED TREE!
 */
     {
     NODETYPE *n, *prev, *parent,*grandparent,*sis,*root;
     if (T)
     	root=T->root;
     else
	root=NULL;	/* used for cases in which don't know the tree and don't care */
     if (rmTaxon==NULL)
	return root;
     if (!isRoot(rmTaxon))
        {
	parent=rmTaxon->anc;
	grandparent=parent->anc; /* might be NULL if parent is the root */
	if (node_tomy(parent)==2)
		{
		sis=sister(rmTaxon);
		if (isRoot(parent)) 
			{
			sis->anc=NULL;	/* make sure this node acquires 'root' status */
			sis->sib=NULL;
			return sis; /* new root of tree is this sister node */
			}
		else
			{
			sis->anc=grandparent;
			prev=prevSib(parent);
			if (prev)
				{
				prev->sib=sis;
				}
			else
				{
				grandparent->firstdesc=sis;
				}
			sis->sib=parent->sib;
			sis->length+=parent->length;
			}
		}
	if (node_tomy(parent)>2)
		{
		prev=prevSib(rmTaxon);
		if(prev)	/* either rmTaxon is the firstdesc or its got a prev sib */
			prev->sib=rmTaxon->sib;
		else
			parent->firstdesc=rmTaxon->sib;
		}

        return root;
        }
     }

NODETYPE * prevSib(NODETYPE* node)

/* returns the sib that points to this sib, or null if this sib is the first desc or if this sib is root */

	{
	NODETYPE *prev, *n;
	prev=NULL;
	if(!isRoot(node))
	    {
	    n=node->anc->firstdesc;
	    while(n != node)  
			    {
			    if (n->sib == NULL)
				     return NULL;   
			    prev=n;
			    n=n->sib;
			    }
	    return prev;
	    }
	else
	    return NULL;
	}

NODETYPE * ReRoot(NODETYPE * atNode)  // Sept 2011: fixed some bugs in here that I fixed in the
									// parallel code in mysmalltreelib.c

/* Reroots a tree in place, returning the node pointer to the new root. The old root node is deleted
	and a new root node is instantiated in its place. New root has id=1 and length=0.0. Nothing else
	is changed. Any time a node becomes a 1-tomy in this process, however, the branch lengths are
	combined.
	*/

    {
	NODETYPE *n, *r;
	if(isRoot(atNode))
	    return atNode; /* don't change the root */
	n=atNode->anc;
	if (!isRoot(n))
	    {
	    r=newnode();
	    r->id = 1; // this will be the new root. By convention its id=1; old root is deleted
	    r->length = 0.0; // also convention
//	    RemoveTaxon(NULL,atNode);
		RemoveTaxonLvAnc(atNode);
	    AddChild(r, atNode);
	    Flip(n);
	    AddChild(r, n);
	    n->length=0; /* leave all the length on the left root's branch */
//	    init_node_ids(r, 0);
	    return r;
	    }
	else
	    return n; /* don't change the root here either */
	   
    }
NODETYPE * ReRoot2(NODETYPE * atNode)  // Sept 2011: fixed some bugs in here that I fixed in the
									// parallel code in mysmalltreelib.c

// RoDo of ReRoot! 
// Make the node the root of the tree (as opposed to making it the sister node of the rest of the tree!).

// Note if we start with a binary root node, we keep that node, even as we reroot on other internal
// nodes. That means the rerootings often have a degree one node on one branch, but should be fine 
// for calculations

    {
	NODETYPE *n, *r;
	if(isRoot(atNode))
	    return atNode; /* don't change the root */
	Flip2(atNode);
	return atNode;
    }
void Flip2(NODETYPE *a)

// the subtree "below" node a, becomes one of the children of 'a' now.

    {
	NODETYPE * b,  *saveAnc, *parent;
	float saveLength;
	b=a->anc;	
	if (!isRoot(b))
	    {
	    Flip(b);  /* recurse until the root, then back up */		
	    }
//	RemoveTaxon(NULL,a);
	RemoveTaxonLvAnc(a);
	AddChild(a, b);
	b->length=a->length; /* flip the branch lengths too */
	return;
    }
void Flip(NODETYPE *a)

// the subtree "below" node a, becomes one of the children of 'a' now.

    {
	NODETYPE * b,  *saveAnc, *parent;
	float saveLength;
	b=a->anc;	
	if (!isRoot(b))
	    {
	    Flip(b);  /* recurse until the root, then back up */		
	    }
//	RemoveTaxon(NULL,a);
	RemoveTaxonLvAnc(a);
	AddChild(a, b);
	b->length=a->length; /* flip the branch lengths too */
#if 1
	if (node_tomy(b)==1)  /* then delete this node and combine branch lengths*/
	    {
	    saveLength=b->length;
//	    RemoveTaxon(NULL,b);
		parent = b->anc;
	    RemoveTaxonLvAnc(b);
	    AddChild(parent, b->firstdesc);
	    b->firstdesc->length+=saveLength;
	    Node_Destructor(b);
	    /* deallocate node b HERE */
	    };
#endif
	return;
    }
/***********************************************************************************/
void traverseMultiplyLength(NODETYPE * node, double multiplier,int roundflag)

/* multiply all branch lengths by a constant and round to nearest integer */

{
	NODETYPE *child;
	double value=0;
//printf("node %s:%f %f %f %i\n",node->taxon_name, node->length,value,multiplier,roundflag);
	value=node->length*multiplier;
	if (roundflag)
		{
		if (value-floor(value)<0.5)
			node->length=floor(value);
		else
			node->length=ceil(value); /* rounding to nearest integer */
		}
	else
		node->length=value;
//printf("node %s:%f %f %f %i\n",node->taxon_name, node->length,value,multiplier,roundflag);
	child=node->firstdesc;
	SIBLOOP(child) 
		traverseMultiplyLength(child, multiplier,roundflag);

	return;
}
/***********************************************************************************/
double treeDurLength(NODETYPE * node)

/* sums the branch durations over tree */

{
	NODETYPE *child;
	double dur;
	if (isRoot(node))
		dur=0.0;
	else
		dur=node->anc->time - node->time;
	child=node->firstdesc;
	SIBLOOP(child)
		dur+=treeDurLength(child);
	return dur;
}
/***********************************************************************************/
double treeLength(NODETYPE * node)

/* sums the branch lengths over tree */

{
	NODETYPE *child;
	double dur;
	if (isRoot(node))
		dur=0.0;
	else
		dur=node->length;
	child=node->firstdesc;
	SIBLOOP(child)
		dur+=treeLength(child);
	return dur;
}
/***********************************************************************************/
double treeAgeSum(NODETYPE * node)

/* sums the node ages over tree */

{
	NODETYPE *child;
	double dur;
	dur=node->time;
	child=node->firstdesc;
	SIBLOOP(child)
		dur+=treeAgeSum(child);
	return dur;
}
/**********************************************************************/

int isNodeDescendant(NODETYPE *nodeA, NODETYPE *nodeB)
/*
 * Is nodeA the strict descendant of nodeB or identical to node B?
 * Returns 1 if it is, 0 if it is not 
 * 
 */

{
NODETYPE *node;
for(node=nodeA;node;node=node->anc) 
	/* worst case, terminates when node = NULL at ancestor of root */
	{
	if (node==nodeB) return 1;
	}  
return 0;    
}
/**********************************************************************/
int group_a_clade(NODETYPE *root, StrListPtr taxaList)

/* is the specified list of taxa a clade on tree 'root'? 
	Note that the list might contain MORE taxa than are found on the tree (i.e.,
	it might be a pruned tree.  We allow this.  First, we make a new
	taxaList that contains the intersection of the taxaList and the MRCA taxa
	on the tree.  This we check to see if it's identical to the MRCA list,
	which only happens when group is consistent with that clade
*/

{
NODETYPE *mrca;   
StrListPtr mrcaTaxa,intersectTaxaList;
mrcaTaxa=newStrList();
mrca=MRCA(root, taxaList);	/* mrca of the taxa list */
if (mrca)
    {
    TreeToTaxaList(mrca,  mrcaTaxa); /* set up list of taxa actually descended from mrca node */
/*    printf("group:\n");xprintStrList(taxaList);
    printf("clade:\n");xprintStrList(mrcaTaxa);*/

#if 0  /* enable this for pruning as described above */
    intersectTaxaList=string_list_intersect(mrcaTaxa,taxaList);
    if (string_lists_same(intersectTaxaList, mrcaTaxa))
#else
    if (string_lists_same(taxaList, mrcaTaxa))
#endif
	return 1;/* are they the same?; if so, is monophyletic */
    else 
	return 0;
    }  
}


/**********************************************************************/

NODETYPE * MRCA(NODETYPE *root, StrListPtr taxaList)
/*
 * On tree with 'root',  returns node of the MRCA of taxa in taxa List (a list of strings) 
 * NOTE: some taxa in list may not be on tree!  In that case we find the MRCA of those that
 * ARE on the tree.  If none of the taxa are on the tree return NULL.(BOMBS)
 * YIKES, does this really work? */

{
NODETYPE *node, *firstTaxonNode=NULL, *otherTaxonNode;
PtrList pOther, p, nodeList=NULL, nLptr;
int nList, k, i;
StrListPtr txPtr;
NODETYPE *s;

nList=lengthList(taxaList);
if (nList<2)
	{
    	doGenericAlert("taxa list has fewer than two taxa!");
	return NULL;
	}


/** convert the taxa list to a (possibly smaller) list of corresponding nodes **/ 

 
for(txPtr=taxaList;txPtr;txPtr=txPtr->next)
    {
    s=find_taxon_name(root,txPtr->s);
    if (s) /* don't include taxa that aren't on tree */
	{
	if(!nodeList) /* create first node, or...*/
	    {
		nodeList=pNewListAlt(sizeof(NODETYPE*));
		nLptr=nodeList;
	    }
	else	    /* add a new node, if list is already there */
	    {
		pListAddNode(nodeList, sizeof(NODETYPE*));
		nLptr=nLptr->next;    
	    }
	nLptr->item=s;
	}
    }

if (!nodeList)
	return NULL;	/* BAIL IF THERE WERE NO TAXA and hence NO NODES */

if (pLengthList(nodeList)<lengthList(taxaList))
	doGenericAlert("MRCA COMMAND: Num nodes less than num labels: You probably have misspelled a taxon name!");

p=nodeList;
for (firstTaxonNode=(NODETYPE *)(p->item);!isRoot(firstTaxonNode);firstTaxonNode=firstTaxonNode->anc)
	/* traverse the ancestry path starting from taxon 1... */
	{
	 /*...and check at each node whether the other taxa are descendants of that node...*/	
	for (pOther=p->next;pOther;pOther=pOther->next)
	    {
	    otherTaxonNode=(NODETYPE *)(pOther->item);
	    if (!isNodeDescendant(otherTaxonNode, firstTaxonNode)) 
		goto a1; /* ...at least one taxon was not a descendant, so...---> */
	    }
	freepList(nodeList);
	return firstTaxonNode; /* all members of list were descendants of this node, so return */
a1:	;   /* ----> ....so traverse to its ancestor and repeat...*/
	}  
freepList(nodeList);
return firstTaxonNode;   /* now the root: just return that by default */ 
}

/**********************************************************************/
NODE  mrca(NODE a, NODE b )

// the mrca of two nodes

{
NODE p,psave;
for (p=a;p;p=p->anc)
	p->nodeFlag=1;
for (p=b;p;p=p->anc)
	if (p->nodeFlag==1)
		goto a1;
a1:psave=p;
for (p=a;p;p=p->anc)
	p->nodeFlag=0;

return psave; 
}
double pathLengthTimeModel(NODE a, NODE b, int model)

// The sum of durations between two nodes: only include branches that match the modelID
{
NODE  p,anc,c;
double T=0.0;
c = mrca (a,b);
for (p=a;p!=c;p=p->anc)
	if (p->modelID == model)
		T += p->anc->time - p->time;
for (p=b;p!=c;p=p->anc)
	if (p->modelID == model)
		T += p->anc->time - p->time;
return T;
}
double pathLengthTime(NODE a, NODE b)

// The sum of durations between two nodes
{
NODE  p,anc,c;
double T=0.0;
c = mrca (a,b);
for (p=a;p!=c;p=p->anc)
	T += p->anc->time - p->time;
for (p=b;p!=c;p=p->anc)
	T += p->anc->time - p->time;
return T;
}
/**********************************************************************/
void setNodeName(NODETYPE *node, char *name)
{
char *copy;
copy=DupStr(name);
myFree(node->taxon_name);
node->taxon_name=copy;
return;

}

/**********************************************************************/
void make_parens(NODETYPE *node, int flag)

/* writes a parens formatted tree description with labels and durations or
lengths.  flag=0: print lengths; flag =1: print durations as lengths,  
flag=2: print rates as lengths, flag=3: print node id's as lengths,
flag=4: print normalized marginal of CLmarg[0] as length, and anc state as second number*/

{
  extern long gNumSites;
  double value, duration;
  int width;
  
  if (flag==4)
  	value = (node->CLmarg)[0]/((node->CLmarg)[0]+(node->CLmarg)[1]+(node->CLmarg)[2]);
  if (!isRoot(node))
	{
	if (flag==0)
		value = node->length;
	else if (flag == 1)
		value = node->anc->time - node->time; /* duration */
	else if (flag == 2)
		value = node->estRate/gNumSites;		/*rate*/
	}

  if (isTip(node)) 
    {
    if (*(node->taxon_name)=='\0')
		{
		width = log10(node->id)+1; 
		printf("tx%-*i", width, node->id);
		}
    else
      	printf("%s",node->taxon_name);
    if (flag == 3)
		printf(":%i",node->id);
    if (flag < 3)
    	printf(":%-8.6f",value);
    if (flag == 4)
		printf(":%-8.6f:%i",value,node->opt);
    }
  else printf("(");

  if (node->firstdesc) make_parens(node->firstdesc,flag);

  if (!isTip(node))
    {
      printf(")");
      if (*(node->taxon_name)!='\0') 
	    printf("%s",node->taxon_name);
      if (!isRoot(node)) 
	 	{
         if (flag == 3)
			printf(":%i",node->id);
         if (flag < 3)
	    	printf(":%-8.6f",value);
		}
	  if (flag==4)
			printf(":%-8.6f:%i",value,node->opt);
    }

  if (node->sib) printf(","),make_parens(node->sib,flag);

}
/***********************************************************************************/
void TreeStats(NODETYPE *root, int * numtips, int * numinternals, int * roottomy)

/* gets some info on a tree, including the number of tips, internal nodes (incl. root),
and the number of immediate descendants of the root node, the roottomy level */

{
	NODETYPE *child;
	*roottomy=0;
	child=root->firstdesc;
	SIBLOOP(child) 
		++(*roottomy);
	*numtips = numdesc(root);
	*numinternals = numIntNodes(root);

	return;
}
/***********************************************************************************/
int node_tomy(NODETYPE *node)

/* number of immediate descendants of this node (including this one!) */

{
	NODETYPE *child;
	int tomy=0;
	child=node->firstdesc;
	SIBLOOP(child) 
		++tomy;

	return tomy;
}
/***********************************************************************************/
int maxorder(NODETYPE *node)
{
	int max,temp;
	NODETYPE *child;
	if (!node) return(-1);
	if (isTip(node)  ) {node->order=0; return (0);}
	max=0;
	child=node->firstdesc;
	SIBLOOP(child) {
			temp=maxorder(child);
			if (temp > max) max = temp;
			}
	node->order=max+1;
	return (max+1);
}
/***********************************************************************************/
void init_free(NODETYPE *node)

/* Initializes the free flag for each node.  By default free is set to 0 for tips, 1 for internal nodes and root*/

{
	NODETYPE *child;
	if (isTip(node)) 
		{
		node->free=0; 
		return;
		}
	else 
		node->free=1;		
	child=node->firstdesc;
	SIBLOOP(child) {
			init_free(child);
			}
	return;
}
/***********************************************************************************/
int numFreeNodes(NODETYPE *node)
{
/* returns number of  nodes in the tree that have their free flags set, meaning
that we are estimating their ages.  Counts the root and tips too! */

	int sum=0;
	NODETYPE *child;
	if (isFree(node)) 
		++sum;
	if (isTip(node))
		return sum;
	child=node->firstdesc;
	SIBLOOP(child) 
		sum += numFreeNodes(child);
	return (sum);	
}

/***********************************************************************************/
void Tree_Initialize(TREE aTree, char *TD, char *name)
{
NODETYPE * root;
int numtips, numinternals, roottomy;
root=string_to_tree(TD);
if (root)
	{
	init_node_ids(root, 0);
	init_free(root); /* sets default to estimate all internal nodes but root */ 
	aTree->root=root;
	aTree->name=DupStr(name);
	aTree->TD=DupStr(TD);
	TreeStats(root, &numtips, &numinternals, &roottomy);
	aTree->numTaxa=numtips;
	aTree->numBranches=numBranches(root);
	aTree->basalTomy=roottomy;
	aTree->cladeSet=NULL;
	aTree->est_b=0.0;
	aTree->est_c=0.0;
	root->anc=NULL;
	aTree->timesAvailable=0;
	aTree->method=USER;
	}
return;
}
/***********************************************************************************/
TREE Subtree_Initialize(TREE T,NODETYPE *node)

// Creates a "tree" by using some subclade of an existing tree. 'node' is the node on existing tree
// that will become the root. Does not change any of the information on the existing tree, but does
// NOT allocate a copy of this subtree--merely uses existing data structure for tree and allocates extra info 
// Careful allocating the atree. Still getting used to not using TREE as object init, as one can in C++
{
TREE aTree;
aTree=(struct treetype *)myMalloc(sizeof(struct treetype));
if (aTree)
	{
	aTree->name=T->name;
	aTree->root=node;
	aTree->numTaxa=T->numTaxa;
	aTree->numBranches=T->numBranches;
	aTree->basalTomy=T->basalTomy;
	aTree->est_b=0.0;
	aTree->est_c=0.0;
	aTree->root->anc=NULL;
	aTree->timesAvailable=0;
	aTree->method=USER;
	}
return aTree;
}
/***********************************************************************************/
void Tree_Destructor(TREE aTree)
{
DisposeTree(aTree->root);
myFree(aTree->name);
myFree(aTree->TD);
/* should free the cladeSet array if present!!! */
myFree(aTree);
return;
}
/***********************************************************************************/
void Node_Destructor(NODETYPE *node)
{
    
	myFree(node->taxon_name);	
	myFree(node);
	return;    
}




/***********************************************************************************/
void DisposeTree(NODETYPE *node)

// This used to be broken! I think this works now but haven't tested it...

	/* Frees up the tree memory and its taxon names */
{
	NODETYPE *child;
	if (!node) return;

	DisposeTree(node->firstdesc);
	DisposeTree(node->sib);
	Node_Destructor(node);
	return;
}
/***********************************************************************************/
NODETYPE *makegroup(void)  
{

	/* Returns a pointer to a tree structure corresponding to everything within a
	parentheses formatted string whose first character is at address 'gStringptr'.
	This function is called recursively each time a left parenthesis is encountered.
	NOTE: the 'gStringptr' MUST point to the first left parens on entry to this function.
	On exit 'gStringptr' should point to the rightmost right parens in the group--
	or to the last character in the name or number after colon; this
	makes it ready to skip to the next character and continue parsing 
	
	NOTE THAT THIS ROUTINE CANNOT HANDLE AN INTERNAL NODE WITH BOTH A NAME AND A LENGTH

	STUPIDLY, THIS ROUTINE DOES NOT USE TOKENS, SO IT DOESN'T PROPERLY TAKE CARE OF IMBEDDED
	SINGLE QUOTES OR BRACKETS, ETC.


	EVEN WORSE, the storage of numbers is corrupted if there is a space between colon and number
	or if the number starts with .xxx rather than 0.xxx.

*/

	NODETYPE *root, *currnode, *prevnode;
	extern char *gStringptr;	/* points to current position in string tree description
								and must be global for recursive calls to work right */
	char *character, *name, *delim=" ,):"; /* taxon name delimiters include space, comma, parens,colon */
	char *singleQuote = "'";
	char* dummy;
	extern int gCount;
	size_t length;
	int first;
	root=newnode();if (root==NULL) return(NULL);
	currnode=root;
	first=1;
	while (*gStringptr != '\0') {
		++gStringptr;
		switch (*gStringptr)  {
			case(LEFTPARENS):{    /* recursively go down into next clade */
				prevnode=currnode;
				currnode=makegroup();
				if(currnode==NULL) return(NULL);
				if (first) {
							prevnode->firstdesc=currnode;
							first=0;
							}
				else prevnode->sib=currnode;
				currnode->anc=root;
				break;
				}
			case(RIGHTPARENS): /* check to see if there is a taxon name after the parens
								OR a number after a colon, and store. 
									First letter must follow colon */
				{
					++gStringptr;  /* look ahead */
					if (*gStringptr == COLON)
						{
							++gStringptr;
							if (isdigit(*gStringptr) || (*gStringptr == '-') || (*gStringptr == '+') ) {  /* only checks first char !!!*/
								root->length=strtod(gStringptr,&dummy);
								length=strcspn(gStringptr,delim);
								gStringptr+=length-1;	/* increment but 
									leave at last character rather than
									one past, to fulfill the definition of function
									above */
								if (root->length <0.0)
								    {
								    printf("** WARNING: A negative branch length was set to ZERO\n");
								    root->length=0.0;
								    }
							}
				
						}
					else
						if (isalnum(*gStringptr))
						/* RECENTLY CHANGED THIS TO 'isalnum'from isalpha*/
							{
							length=strcspn(gStringptr,delim);
							name=(char *)myMalloc((length+1)*sizeof(char));
		
							myFree(root->taxon_name); 
		
							root->taxon_name=name;	
							if (name==NULL) return(NULL);
							strncpy(root->taxon_name,gStringptr,length);
							*((root->taxon_name)+length)='\0';  
							
							gStringptr+=length-1;	/* see comment above */

							}
						else
							--gStringptr;	/* it was neither a name or number */
					return(root);
				}
			default:{  /* check for valid taxon name or number after name and store */
					if (*gStringptr == COLON)
						{
							++gStringptr;
							if (isdigit(*gStringptr) || (*gStringptr == '-') || (*gStringptr == '+')) {  /* only checks first char !!!*/
								currnode->length=strtod(gStringptr,&dummy);
								if (currnode->length <0.0)
								    {
								    printf("** WARNING: A negative branch length was set to ZERO\n");
								    currnode->length=0.0;
								    }
								length=strcspn(gStringptr,delim);
								gStringptr+=length-1;	/* increment but
									see comment above for explanation of -1 */
							}
				
						}






			else
				if (IsItAlphaNum(*gStringptr)) {  /* only checks first letter !!!*/
					prevnode=currnode;
					currnode=newnode();
					if (currnode==NULL) return(NULL);
					if (first) {
								prevnode->firstdesc=currnode;
								first=0;
								}
					else prevnode->sib=currnode;
					currnode->anc=root;
					length=strcspn(gStringptr,delim);
					name=(char *)myMalloc((length+1)*sizeof(char));

					myFree(currnode->taxon_name); 

					currnode->taxon_name=name;	
					if (name==NULL) return(NULL);
					strncpy(currnode->taxon_name,gStringptr,length);
					*((currnode->taxon_name)+length)='\0';  
					
					gStringptr+=length-1;	/* increment but only if two or more characters */
					}
				}
		
		}
	}
}

/***********************************************************************************/
NODETYPE *makegroup2(void)  
{

	/* Returns a pointer to a tree structure corresponding to everything within a
	parentheses formatted string whose first character is at address 'gStringptr'.
	This function is called recursively each time a left parenthesis is encountered.
	NOTE: the 'gStringptr' MUST point to the first left parens on entry to this function.
	On exit 'gStringptr' should point to the rightmost right parens in the group--
	or to the last character in the name or number after colon; this
	makes it ready to skip to the next character and continue parsing 
	
	NOTE THAT THIS ROUTINE CANNOT HANDLE AN INTERNAL NODE WITH BOTH A NAME AND A LENGTH

	STUPIDLY, THIS ROUTINE DOES NOT USE TOKENS, SO IT DOESN'T PROPERLY TAKE CARE OF IMBEDDED
	SINGLE QUOTES OR BRACKETS, ETC.

	BUG! If there is a blank in a taxon name within single quotes, the name will not be parse right.
*/

	NODETYPE *root, *currnode, *prevnode;
	extern char *gStringptr;	/* points to current position in string tree description
								and must be global for recursive calls to work right */
	char *character, *name, *delim=" ,):"; /* taxon name delimiters include space, comma, parens,colon */
	char* dummy;
	extern int gCount;
	size_t length;
	int first;
	root=newnode();if (root==NULL) return(NULL);
	currnode=root;
	first=1;
	while (*gStringptr != '\0') {
		++gStringptr;
		switch (*gStringptr)  {
			case(LEFTPARENS):{    /* recursively go down into next clade */
				prevnode=currnode;
				currnode=makegroup2();
				if(currnode==NULL) return(NULL);
				insertNode(currnode, prevnode);	/* add the new node (actually,
				    the whole subtree!) to prevnode,  its ancestor */
				break;
				}
			case(RIGHTPARENS): /* check to see if there is a taxon name after the parens
								OR a number after a colon, and store. 
									First letter must follow colon */
				{
					++gStringptr;  /* look ahead */
					if (*gStringptr == COLON)
						{
							++gStringptr;
							if (isdigit(*gStringptr) || (*gStringptr == '-') || (*gStringptr == '+') ) {  /* only checks first char !!!*/
								root->length=strtod(gStringptr,&dummy);
								length=strcspn(gStringptr,delim);
								gStringptr+=length-1;	/* increment but 
									leave at last character rather than
									one past, to fulfill the definition of function
									above */
								if (root->length <0.0)
								    {
								    printf("** WARNING: A negative branch length was set to ZERO\n");
								    root->length=0.0;
								    }
							}
				
						}
					else
						if (isalnum(*gStringptr))
						/* RECENTLY CHANGED THIS TO 'isalnum'from isalpha*/
							{
							length=strcspn(gStringptr,delim);
							name=(char *)myMalloc((length+1)*sizeof(char));
		
							myFree(root->taxon_name); 
		
							root->taxon_name=name;	
							if (name==NULL) return(NULL);
							strncpy(root->taxon_name,gStringptr,length);
							*((root->taxon_name)+length)='\0';  
							
							gStringptr+=length-1;	/* see comment above */

							}
						else
							--gStringptr;	/* it was neither a name or number */
					return(root); /*... of the current subtree...*/
				}
			default:
			    {  /* check for valid taxon name or number after name and store */
				    if (*gStringptr == COLON)
						{
							++gStringptr;
							if (isdigit(*gStringptr) || (*gStringptr == '-') || (*gStringptr == '+')) {  /* only checks first char !!!*/
								currnode->length=strtod(gStringptr,&dummy);
								if (currnode->length <0.0)
								    {
								    printf("** WARNING: A negative branch length was set to ZERO\n");
								    currnode->length=0.0;
								    }
								length=strcspn(gStringptr,delim);
								gStringptr+=length-1;	/* increment but
									see comment above for explanation of -1 */
							}
				
						}
				    else
					if (IsItAlphaNum(*gStringptr))  /*  its a terminal, add it */
						{  /* only checks first letter !!!*/
						prevnode=currnode;
						currnode=newnode();
						if (currnode==NULL) return(NULL);
						insertNode(currnode, prevnode);
						length=strcspn(gStringptr,delim);
						name=(char *)myMalloc((length+1)*sizeof(char));
	
						myFree(currnode->taxon_name); 
	
						currnode->taxon_name=name;	
						if (name==NULL) return(NULL);
						strncpy(currnode->taxon_name,gStringptr,length);
						*((currnode->taxon_name)+length)='\0';  
						
						gStringptr+=length-1;	/* increment but only if two or more characters */
						}
			}
		
		}
	}
}
static void insertNode(NODETYPE *node,  NODETYPE* anc)

/* ....looks dangerous, if the node is a polytomy, this seems to delete some children! */

{
                node->anc=anc;
                if (anc->firstdesc == NULL)  /* this is first child of anc */
                        anc->firstdesc=node;
                else                    /* this is nth child and has a sib */
                        anc->firstdesc->sib=node;
}
/***********************************************************************************/
NODETYPE *newnode(void)
{

	NODETYPE *node;
	node=(NODETYPE *)myMalloc(sizeof(NODETYPE));		

	if (node==NULL) fatal("Toast");
	node->anc=NULL;
	node->firstdesc=NULL;
	node->sib=NULL;
	node->nodePtr=NULL;
	node->isCompactNode=0;
	node->isQueryNode=0;
	node->toggleDesc=0;
	node->taxon_name=(char *)myMalloc(sizeof(char));
	node->length=FLT_MAX;		/* Big number lets us check later */ 
	node->time=0.0; /* NOTE: 'drawtree' checks this value at the root node
		to determine if times have been set */
	node->minAge=0.0;
	node->maxAge=/*  1.0  */  1.0e20;
	node->nodeIsConstrainedMax=0;
	node->nodeIsConstrainedMin=0;
	node->free=0;
	node->like = 0.0;
	node->id=0;
	node->modelID=0;
	node->nodeFlag=0;
	

	if (node->taxon_name ==NULL) fatal("Couldn't allocate name in node");;
	*(node->taxon_name)='\0'; 	/* store a null string for now */
	return (node);
}


/***********************************************************************************/
void copyNodeInfo(NODETYPE *source,NODETYPE *dest)

/* Copies SOME information about one node to another node */

{

	dest->taxon_name=DupStr(source->taxon_name);
	dest->length=source->length;
	dest->time=source->time;
	dest->minAge=source->minAge;
	dest->maxAge=source->maxAge;
	dest->id=source->id;
	dest->free=source->free;
	dest->numdesc=source->numdesc;
	dest->estRate=source->estRate;
	dest->nodeReal=source->nodeReal;
	dest->nodeIsConstrainedMax=source->nodeIsConstrainedMax;
	dest->nodeIsConstrainedMin=source->nodeIsConstrainedMin;
	return;
}


/***********************************************************************************/
void ClusterHistogram(NODETYPE * node, long *histo,long TSize)

/* moves through a tree setting up a histogram of cluster sizes in which bins are 
 * on a log 2 scale: thus for a tree of size 32.  Note that the last bin is slightly larger...
 * 
 *	2-3
 *	4-7
 *	8-16
 */

{
	NODETYPE *child;
	int ix;
	long N;
	if (!node || isTip(node))
		return;
	if (!isRoot(node))
	    {
	    N=node->numdesc;
	    ix=floor(LOG2(MIN(N,TSize-N)))-1;
	    if (N==TSize/2)
		--ix;		/* put this special cluster size in next lowest bin...otherwise
				    it sits there almost alone in the last bin */
	    ++histo[ix];
	    }
	child=node->firstdesc;
	SIBLOOP(child) 
		ClusterHistogram(child, histo,TSize);
	return;
}


/***********************************************************************************/

int numNodes(NODETYPE *node)
{

return numdesc(node)+numIntNodes(node);


}

/***********************************************************************************/

long numUnMarkedDesc(NODETYPE *node)

/* determines the number of leaves descended from a node EXCEPT for descendant clades
	that are 'marked' */ 
/* Careful! Must unmark the root node of this subtree in the caller and then remark when done */

{
	long sum=0;
	NODETYPE *child;
	if (!node) return(-1);
	if (isNodeMarked(node))
		return 0;
	if (isTip(node)) 
		{
		return 1;
		}
	child=node->firstdesc;
	SIBLOOP(child) 
		sum+=numUnMarkedDesc(child);
	return (sum);
}
/***********************************************************************************/

int numdesc(NODETYPE *node)

/* determines the number of leaves descended from every node and stores them at node */

{
	long sum=0;
	NODETYPE *child;
	if (!node) return(-1);
	if (isTip(node)) 
		{
		node->numdesc=1; 
		return (1);
		}
	child=node->firstdesc;
	SIBLOOP(child) 
		sum+=numdesc(child);
	node->numdesc=sum;
	return (sum);
}
/***********************************************************************************/
void printnode(NODETYPE *node)
{
    double duration, estR;
    NODETYPE *anc;
    if (isNodeMarked(node))
	printf("*");
    if (!isRoot(node))
	{
	anc=node->anc;
	duration=node->anc->time-node->time;
	printf("node %3i (%s) age=%4.2g | anc %3i (%s) age=%4.2g | dur=%4.2g len=%4.2g rate=%6.3g nodeReal=%6.3g age bounds=[%g..%g]\n",
	    node->id, node->taxon_name,node->time, 
	    anc->id, anc->taxon_name, anc->time, 
	    duration,node->length,node->estRate,node->nodeReal,nodeUpperBound(node),nodeLowerBound(node));
	}
    else
	printf("node %3i (%s) age=%4.2f len=%4.2g\n",
	    node->id, node->taxon_name,node->time, node->length);
    return;
}
/***********************************************************************************/
void printtree(NODETYPE *node)
{
	NODETYPE *child;
	if (!node) 
		return;
	printnode(node);
#if 0
	if (node->nodePtr)
		{
		printf("\n\n --->Subtree info for this node:\n\n");
		DrawTree(node->nodePtr,0,0); /*printtree(node->nodePtr);*/
		printf("\n\n --->End of subtree info for this node:\n\n");
		}
#else
	if (node->nodePtr)
		printf("-->Points to node with label %s\n",node->nodePtr->taxon_name);
#endif
	child=node->firstdesc;
	SIBLOOP(child) 
		printtree(child);
	return;
}
/***********************************************************************************/
void printLikes(NODETYPE *node)
{
	NODETYPE *child;
	if (!node) 
		return;
	printnodeLike(node);
	child=node->firstdesc;
	SIBLOOP(child) 
		printLikes(child);
	return;
}
/***********************************************************************************/

int numIntNodes(NODETYPE *node)
{
/* returns number of internal nodes in the tree.  Counts the root too, so this must
be subtracted later! */

	int sum=0;
	NODETYPE *child;
	if (!node) return(-1);
	if (isTip(node)) 
		return (0);
	child=node->firstdesc;
	SIBLOOP(child) 
		sum += numIntNodes(child); /* add one for each child and all that children's*/
	return (1+sum);	/* the 1 is to count this node, which must be internal */
}
/***********************************************************************************/
int numBranches(NODETYPE *node)
{
/* returns number of branches in the tree.  Does NOT count a branch subtending
the root.  Note that this number is <= 2*ntaxa-2 because of polytomies */

	return numdesc(node)+numIntNodes(node)-1; /* for the root */
}
/***********************************************************************************/
/***********************************************************************************/

NODETYPE *string_to_tree(char *tree_description)

/*  Takes a string tree description in NEXUS parens format and returns the root
	node of a linked-list tree structure.  The names of the taxa are stored
	in a string pointed to by node->taxon_name.  The routine first checks to
	see if the string description is valid--hopefully it's good at this.
	Stores lengths of branches and internal taxon names.
	
	Returns NULL if fails. 
*/ 


{
	NODETYPE *root;
	extern char *gStringptr;
	int ix;
	gStringptr=tree_description;
	if (stringcheck(tree_description)) {;
		gStringptr=strchr(tree_description,LEFTPARENS); /* move to first occurrence of left paren*/

		root=makegroup();		
		}
	else 
		return(NULL);
	return(root);
}
/***********************************************************************************/
int stringcheck(char *td)
	/*  Checks a tree description statement and returns error 
	if it has unbalanced parentheses, or if a right parens precedes a left out of turn.
	It does NOT catch a premature termination caused by an early right parens that matches
	up with a left parens */

{
	long count=0,parenscount=0;
	if (*td == '\0') return (0);  /* string is empty */
	while(*td)  {  /* while character is not null */
		if (*td == LEFTPARENS) ++parenscount;
		if (*td == RIGHTPARENS) --parenscount;				
		if (parenscount < 0) return (0); /* right parens preceded left */
		td++;  /* used to be *td++  */
		}
	if (parenscount != 0) return (0); /* unbalanced parens */
	return (1);	/* string OK  */
}
/***********************************************************************************/
NODETYPE * find_taxon_name(NODETYPE *node,char *name)
/* returns the node that has a taxon name that matches 'name' or NULL
if it fails to find */


{
	NODETYPE *child, *found_node;

	if (node->taxon_name)
		if (isEqual(name,node->taxon_name))
			return node;
	child=node->firstdesc;
	SIBLOOP(child) 
		if (found_node=find_taxon_name(child,name) )
			return found_node;
	return NULL;
}
/***********************************************************************************/

int collapseLengthsTree2Tree(TREE t1,TREE t2)

/* Moves through two identical trees simultaneously, and whenever it finds a zero-length
	branch on tree1, it collapses that branch AND the corrsponding one on tree2.
	Information about the length of the branch on tree2 is discarded! If the trees are not
	isomorphic, the results are unpredictable.
*/
{
int retFlag=0;
if (any_zero_internal_branches(t1->root))
	retFlag=1;
while(any_zero_internal_branches(t1->root))
				collapse_zero_2trees(t1->root,t2->root);

return retFlag;
}
static void collapse_zero_2trees(NODETYPE * node1, NODETYPE * node2)
{
NODETYPE *  child1,* child2;

if(!isRoot(node1) && !isTip(node1))
	{
	if (node1->length==0.0)
			{
			collapse_node(node1);
			collapse_node(node2);
			return;
			}
	}
child1=node1->firstdesc;
child2=node2->firstdesc;
for (;child1;child1=child1->sib,child2=child2->sib)
	collapse_zero_2trees(child1,child2);

return;
}
/***********************************************************************************/
int any_zero_internal_branches(NODETYPE *node)
/* returns 1 if any INTERNAL nodes, other than the root node, have zero-length branches */


{
	NODETYPE *child;

	if(!isRoot(node) && !isTip(node))
	   if (node->length==0.0)
			return 1;
	child=node->firstdesc;
	SIBLOOP(child) 
		if (any_zero_internal_branches(child) )
			return 1;
	return 0;
}
/***********************************************************************************/
int any_zero_terminal_branches(NODETYPE *node)
/* returns 1 if any TERMINAL nodes,  have zero-length branches */


{
	NODETYPE *child;

	if(isTip(node))
	   if (node->length==0.0)
			return 1;
	child=node->firstdesc;
	SIBLOOP(child) 
		if (any_zero_terminal_branches(child) )
			return 1;
	return 0;
}
/***********************************************************************************/
void collapse_zero(NODETYPE *node)
/* collapses THE FIRST FOUND zero-length branch to polytomies*/


{
	NODETYPE *child;

	if(!isRoot(node) && !isTip(node))
	   if (node->length==0.0)
			{
			collapse_node(node);
			return;
			}
	child=node->firstdesc;
	SIBLOOP(child) 
		collapse_zero(child);
	return;
}
/***********************************************************************************/
void collapse_node(NODETYPE *node)
{

// Remove the branch subtending node, collapsing to a polytomy
// Ignore if root or tip

NODETYPE *anc, *right, *left, *first_desc, *last_desc, *nd;

if (isTip(node) && isRoot(node))return;
anc=node->anc;
first_desc=node->firstdesc;
right=node->sib;  

for (nd=first_desc;nd->sib;nd=nd->sib);
last_desc=nd;/* this fragment finds node's last immediate desc. */
last_desc->sib=right;

if (anc->firstdesc==node)
    anc->firstdesc=first_desc; /* this node is the leftmost */
else			/*find the node to the left of it */
    {
    for (nd=anc->firstdesc;nd->sib!=node;nd=nd->sib);
    left=nd;
    left->sib=first_desc;
    }
for (nd=first_desc;nd!=last_desc->sib;nd=nd->sib)
    nd->anc=anc;
// Node_Destructor(node); /* WATCH OUT! THAT's for damn sure; this screws up the recursion in outer routine*/
return;    
}

/***********************************************************************************/
NODETYPE * find_id(NODETYPE *node,int id)
/* returns the node that has an id that matches 'id' or NULL
if it fails to find */


{
	NODETYPE *child, *found_node=NULL;

	if (node->id == id)
			return node;
	child=node->firstdesc;
	SIBLOOP(child) 
		if (found_node=find_id(child,id) )
			return found_node;
	return NULL;
}

/***********************************************************************************/
void init_node_ids(NODETYPE *node, int id)
{
	NODETYPE *child;
	if (!node) 
		return;
	gId=id+1;
	node->id=gId;
	child=node->firstdesc;
	SIBLOOP(child) 
		init_node_ids(child, gId);
	return;
}
/***********************************************************************************/
void print_ages(NODETYPE *node, double time, double calAge,int method)
{
    NODETYPE *child;
    extern long gNumSites;

/*...ages and misc...*/

    if (time != 0.0)
    	node->time*= (calAge/time); /* trap this occasional error */
    if (!isTip(node))
      {
      if(*(node->taxon_name))
	    printf(" [*]  %.8s\t", node->taxon_name);
      else
	    printf("     (%i)\t", node->id);
      }
    else
	printf("      %.8s\t",node->taxon_name);
    if (node->free)
	    printf("      ");
    else
	    printf("     *");
    printf(" [%1i] ",node->modelID);
    if (node->nodeIsConstrainedMin)
    	printf("\t%7.2f\t",node->minAge);
    else
	printf("\t   --   ");
    if (node->nodeIsConstrainedMax)
    	printf("%7.2f\t",node->maxAge);
    else
	printf("   --\t");

    printf("%7.2f\t\t",node->time);
    if (!isRoot(node))
      switch (method)
		{
		case USER:   /* user supplied chronogram */  	
			printf("   --   \t%.4e\n",node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case LaF:    	
			printf("%.4e\t%.4e\n",node->estRate/gNumSites,node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case LFLOCAL:    	
			printf("%.4e\t%.4e\n",node->estRate/gNumSites,node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case NP:
			printf("   --   \t%.4e\n",node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case PENLIKE:
			printf("%.4e\t%.4e\n",node->estRate/gNumSites,node->length/(node->anc->time-node->time)/gNumSites);
			break;
		default:;
		}
     else
	printf("\n");



    if(!isTip(node))
	{
    	child=node->firstdesc;
    	SIBLOOP(child) 
		print_ages(child, time, calAge,method);
	}
    return;
    
}
/***********************************************************************************/
void print_named_ages(NODETYPE *node)
{
// prints out the ages of internal named nodes only...
    NODETYPE *child;
    if (!isTip(node))
      {
      if(*(node->taxon_name))
	    printf(" [**]  %.8s\t%7.2f\n", node->taxon_name,node->time);
      }

    if(!isTip(node))
	{
    	child=node->firstdesc;
    	SIBLOOP(child) 
		print_named_ages(child);
	}
    return;
    
}
/***********************************************************************************/
void summarize_rates(TREE t)
{
extern long gNumSites;
double *r,min=1e20,max=-1e20,mean,sdev,adev,var,skew,curt;
long ix=0,i;
NODETYPE * root;
root=t->root;
r=(double *)myMalloc((t->numBranches)*sizeof(double));
recurse_summarize_rates(t->root,&ix,r);

moment(&r[0]-1,t->numBranches,&mean,&adev,&sdev,&var,&skew,&curt);
for (i=0;i<t->numBranches;i++)
		{
		if (r[i]>max)max=r[i];
		if (r[i]<min)min=r[i];
		}
printf("\nSummary of rate variation (substitutions per site per unit time)\n  Mean    = %.4g\n  Std Dev = %.4g\n  Min     = %.4g\n  Max     = %.4g\n  Range   = %.4g\n  Ratio   =  %.4g\n",mean/gNumSites,sdev/gNumSites,min/gNumSites,max/gNumSites,(max-min)/gNumSites,max/min);

return;

}

static void recurse_summarize_rates(NODETYPE * n, long * ix, double r[])
{
NODETYPE *child;
if (!isRoot(n))
	{
	r[*ix]=n->estRate;
	++(*ix);
	}
if (!isTip(n))
	{
	child=n->firstdesc;
	SIBLOOP(child)
		recurse_summarize_rates(child,ix,r);
	}
return;
}


void print_rates(NODETYPE *node,int method)
{
    NODETYPE *child;
    extern long gNumSites;
   

/* ...now rates...*/


    if (!isRoot(node))
	{
    	if(*(node->taxon_name))
	    printf("\t%.7s\t", node->taxon_name);
   	 else
	    printf("\t%i\t", node->id);

	switch (method)
		{
		case USER:   /* user supplied chronogram */  	
			printf("\t\t%.4g\n",node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case LaF:    	
			printf("\t\t%.4g\t%.4g\n",node->estRate/gNumSites,node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case LFLOCAL:    	
			printf("\t\t%.4g\t%.4g\n",node->estRate/gNumSites,node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case NP:
			printf("\t\t%.4g\n",node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case PENLIKE:
			printf("\t\t%.4g\t%.4g\n",node->estRate/gNumSites,node->length/(node->anc->time-node->time)/gNumSites);
			break;
		default:;
		}
	}
    if(!isTip(node))
	{
    	child=node->firstdesc;
    	SIBLOOP(child) 
		print_rates(child,method);
	}
    return;
    
}
/***********************************************************************************/
void convert_branchlength_to_time(NODETYPE *root)

/* takes a set of ultrametric branch lengths and just scales them to [0,1] over tree and plops
them into the time field of all nodes. 
DO THIS WHEN THE BLFORMAT COMMAND HAS OPTION CLOCK=YES
*/

{
void recurseBTT(NODETYPE *node,  double scalefactor, double tcurr);
double scalefactor;
scalefactor=calcMaxToTip (root);
//printf("max to tip:%f\n",scalefactor);
recurseBTT(root,  scalefactor,  1.0);   
}
void recurseBTT(NODETYPE *node,  double scalefactor, double tcurr)
{
	NODETYPE *child;
	if (isRoot(node))
			node->time = tcurr;
	else
			node->time = tcurr - node->length/scalefactor;
	child=node->firstdesc;
	SIBLOOP(child) 
		recurseBTT(child, scalefactor, node->time);
	return;
    
}

void scaleTree(NODETYPE * root, double calAge, NODETYPE * calNode)

/* Scale all times on a tree according to given calibration time for one node */
{
double scaleFactor;

scaleFactor=calAge/calNode->time;
(void)preOrderArg(root,SThelper,scaleFactor);
return;
}
static double SThelper(NODETYPE * node,double factor)
{
node->time *= factor;
return 0.0;  // required by the prototype...
}

/***********************************************************************************/
static int gC;

double * sort_node_time_array(NODETYPE *root)

/* returns a pointer to an array that has the sorted INTERNAL node times in increasing order
from the present backward, including the root node time.
N.B.  There will be one element of this array for each speciation event, even if there
is a polytomy, and the diversity increases several steps at one instant in time.  Therefore,
the number of elements is equal to ntaxa-1, which happens to be the number of internal nodes
in a fully bifurcating tree.  It does not matter if the cladogram is actually bifurcating
or not! (This was a headache to resolve)*/

{
void recurse_time_get(NODETYPE *node,  double times[]);
int compar(const void *v1,  const void *v2);
int n, i;
double *times;
n=numdesc(root)-1; 
gC=0;
times=(double*)myMalloc(n*sizeof(double));
recurse_time_get(root, times);
qsort((void *)times, n, sizeof(double), compar);

/*
printf("Sorted times of internal nodes\n");
printf ("%i internal times\n", n);
for (i=0;i<n;i++)
    printf("%f\n", times[i]);
*/
return times;
}

void recurse_time_get(NODETYPE *node,  double times[])
{
	NODETYPE *child;
	int i;
	if (!isTip(node))
		{
		for (i=1;i<=node_tomy(node)-1;i++) 
		    {
		    times[gC]=node->time;
		    ++gC;
		    }
		}
	child=node->firstdesc;
	SIBLOOP(child) 
		recurse_time_get(child,times);
	return;
    
}
int compar(const void *v1,  const void *v2)
{
double V;
V= *(double *)v1 - *(double *)v2;
if (V>0.0)
	return 1;
else if (V <0.0)
	return -1;
else if (V==0.0)
	return 0;

}
double get_sum_durations(NODETYPE *node)
{
	NODETYPE *child;
	double dur=0;
	if (isTip(node))
	    return 0.0;
	child=node->firstdesc;
	SIBLOOP(child)
		{
		dur+=node->time-child->time;
		/*printf("%f %f %f\n", node->time, child->time, dur);*/
		dur+=get_sum_durations(child);
		}
	return dur;
    
}
/***********************************************************************************/
void print_tree_list(PtrList treeList)
{
PtrList lnode;
TREE thisTree;
lnode=treeList;
LISTLOOP (lnode)
	{
	thisTree=lnode->item;
	printf("Tree %s\nNum taxa = %i\nNum branches = %i\nBasal tomy = %i\n\n",
		thisTree->name,thisTree->numTaxa,thisTree->numBranches,
		thisTree->basalTomy);
	DrawTree(thisTree->root,0, 0);
	}
return;
}
/***********************************************************************************/
void TreeToTaxaList(NODETYPE *node,  StrListPtr taxaList)	
{

/* Moves through clade from node, compiling list of descendants;
on entry taxaList must be a valid pointer to a possibly empty list */

	NODETYPE *child;
	if (isTip(node)) 
		{
		appendStrList(taxaList, node->taxon_name);
		return;
		}

	else	/* interior node */
		{
		child=node->firstdesc;
		SIBLOOP(child)
			{
			TreeToTaxaList(child, taxaList);
			}
		return;
		}
}
/***********************************************************************************/
void TreeToTaxaPtrList(NODETYPE *node,  PtrList NodeList)	
{

/* Moves through clade from node, compiling list of terminals NODES (!);
on entry taxaList must be a valid pointer to a possibly empty list */

	NODETYPE *child;
	if (isTip(node)) 
		{
		pListAddItem(NodeList, node);
		return;
		}

	else	/* interior node */
		{
		child=node->firstdesc;
		SIBLOOP(child)
			{
			TreeToTaxaPtrList(child, NodeList);
			}
		return;
		}
}
/***********************************************************************************/
void TreeToNodePtrList(NODETYPE *node,  PtrList NodeList)	
{

/* Moves through clade from node, compiling list of all nodes;
on entry nodeList must be a valid pointer to a possibly empty list */

	NODETYPE *child;
	pListAddItem(NodeList, node);
	child=node->firstdesc;
	SIBLOOP(child)
		{
		TreeToNodePtrList(child, NodeList);
		}
	return;
}
/***********************************************************************************/
void ABCSuperTree(NODETYPE *node, StrListPtr UniqueList, 
		    char **dataMatrix,float *wtset)	
{

/* on entry external gColumn MUST point to a valid column in dataMatrix, zero-offset;
i.e.,  start the thang at ZERO
'UniqueList' is the list of all the taxa in the analysis
 */

	extern int gColumn; /* declared in ReadNexusFile2 */
	NODETYPE *child;
	int ll, mm, j, numTaxa;
	char * taxon;
	StrListPtr aTaxaList;
	if (isTip(node)) 
		return;
	else	/* interior node */
		{
		if (!isRoot(node)) /* Don't add a char for whole tree */
		    {
		    wtset[gColumn]=node->length;
		    aTaxaList=newStrList();
		    TreeToTaxaList(node, aTaxaList);
		    numTaxa=lengthList(aTaxaList);
		    for(ll=1;ll<=numTaxa;ll++)
			    {
			    taxon=getkthStr(aTaxaList, ll);
			    mm=findMatchStr(UniqueList, taxon);
			    if (mm)
				dataMatrix[mm-1][gColumn]='1';
			    }
		    freeStrList(aTaxaList);
		    ++gColumn;
		    }
		child=node->firstdesc;
		SIBLOOP(child)
			{
			ABCSuperTree(child, UniqueList, dataMatrix,wtset);
			}
		return;
		}
}
/***********************************************************************************/
void ABCSuperTreePurvis(NODETYPE *node, StrListPtr UniqueList, 
		    char **dataMatrix,float *wtset)	
{
/* Sets up a datamatrix coded according to Purvis' supertree method.
 * Visits each internal node.  Constructs a set of characters for that node consisting of
 * 1's for a subclade and 0's for the other taxa in the clade,  not in that subclade.
 * Does a character for all the immediate subclades of that node,  then recurses in.
 * 
 * 
 */


/* on entry external gColumn MUST point to a valid column in dataMatrix, zero-offset;
i.e.,  start the thang at ZERO */

	extern int gColumn; /* declared in ReadNexusFile2 */
	NODETYPE *child;
	int ll, mm, j, numTaxa, numTaxa2;
	char * taxon;
	StrListPtr aTaxaList, aTaxaList2;
	if (isTip(node)) 
		return;
	else	/* interior node */
		{
		aTaxaList=newStrList();
		TreeToTaxaList(node, aTaxaList);
		numTaxa=lengthList(aTaxaList);
		child=node->firstdesc;
		SIBLOOP(child)
			{
			if (!isTip(child)) 
			    {
			    aTaxaList2=newStrList();
			    TreeToTaxaList(child,  aTaxaList2);
			    numTaxa2=lengthList(aTaxaList2);
			    for(ll=1;ll<=numTaxa;ll++)
				{
				 taxon=getkthStr(aTaxaList, ll);
				 mm=findMatchStr(UniqueList, taxon);
				 if (mm)
				    dataMatrix[mm-1][gColumn]='0';
				 }
			    for(ll=1;ll<=numTaxa2;ll++)
				{
				 taxon=getkthStr(aTaxaList2, ll);
				 mm=findMatchStr(UniqueList, taxon);
				 if (mm)
				    dataMatrix[mm-1][gColumn]='1';
				 }
				 
			
			    ++gColumn;
			    }
			
			ABCSuperTreePurvis(child, UniqueList, dataMatrix,wtset);
			}
		freeStrList(aTaxaList);
		return;
		}
}
/***********************************************************************************/
PtrList Tree2CladeSet(TREE thisTree, StrListPtr allTaxaList)
    {
    PtrList CladeSetList;
    CladeSetList=pNewList(); 
    Tree2CladeSets(thisTree->root, allTaxaList, thisTree->numTaxa, CladeSetList);
    return CladeSetList;
    }
void Tree2CladeSets(NODETYPE *node, StrListPtr allTaxaList, int nTaxa, 
		    PtrList SetList)	
{

/*  Recurses through a tree, obtaining a list of all the clades in the tree.
 *  A generic list,  'SetList',  of pointers to set vectors is repeatedly
 *  added to as we traverse the tree. 'allTaxaList' is a string list containing
 *  all the taxon names.  Sets are represented as integer (binary) vectors of
 *  size 'nTaxa'.  Membership in a clade for some taxon is signified by a 1 in 
 *  that position.  Note that 'allTaxaList' HAS TO BE CREATED ONCE and then used
 *  for all trees,  otherwise each set ordering will be unique and sets won't
 *  make sense.
 * 
 */


	NODETYPE *child;
	int ll, mm, j, numTaxa, i;
	char * taxon;
	Set cladeSet;
	StrListPtr aTaxaList;
	if (isTip(node)) 
		return;
	else	/* interior node */
		{
		if (!isRoot(node)) /* Don't add a char for whole tree */
		    {
		    cladeSet=newSet(nTaxa);
		    aTaxaList=newStrList();
		    TreeToTaxaList(node, aTaxaList); /* get taxa in this clade */
		    numTaxa=lengthList(aTaxaList); /* how many? */
		    for(ll=1;ll<=numTaxa;ll++)
			    {
			    taxon=getkthStr(aTaxaList, ll);
			    mm=findMatchStr(allTaxaList, taxon); /* find the position
					in the vector for this taxon */
			    if (mm)
				add_to_set(cladeSet, mm);
			    }
		    freeStrList(aTaxaList);
		    pListAddItem(SetList, cladeSet);
		    }
		child=node->firstdesc;
		SIBLOOP(child)
			{
			Tree2CladeSets(child, allTaxaList, nTaxa, SetList);
			}
		return;
		}
}
void printCladeSets(PtrList SetList)
    {
    int i, nTaxa;
    Set cladeSet2, cladeSet1;
    PtrList curP2, curP1;
    curP2=SetList;
    while(curP2)
	{
	    cladeSet2=(Set)(curP2->item);
	    print_set(cladeSet2);
	    curP2=curP2->next;
	}
#if 0	
    curP1=SetList;
    while(curP1)
      {
	cladeSet1=(Set)(curP1->item);
	curP2=curP1->next;
	while(curP2)
	    {
		cladeSet2=(Set)(curP2->item);
		test_set(cladeSet1, cladeSet2);
		curP2=curP2->next;
	    }
	curP1=curP1->next;
      }
#endif
    return;
    }	


void rootToTips(NODETYPE* node,double curLen)

/* Calculates distances from root to each tip and prints...*/

{
	NODETYPE *child;

	if (!isRoot(node))
		{
		curLen+=node->length;	/* don't add length under the root */
		}
	if (isTip(node)) 
			{
			printf("Root to tip:%s: %f\n",node->taxon_name,curLen);
			}
	child=node->firstdesc;
	SIBLOOP(child) {
			rootToTips(child,curLen);
			}
	return ;
}
/***********************************************************************************/


                                   r8s/TreeUtils.h                                                                                     0000644 0000766 0000120 00000021724 11637721346 013354  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #ifndef _TREEUTILS
#define _TREEUTILS

#define LARGE_NODES	0		/* set to 1 iff we want to use all fields in every node,
					if this is 0, we cannot do HMM now! */
#define FLT_MAX 1e35	/* no longer defined in limits.h--do it here temporarily */
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "structures.h"
#define LF 10
#define RETURN 13
#define COLON ':'
#define BAR	'|'
#define PLUS '+'
#define DASH '-'
#define	SPACE	' '
#define COMMA	','
#define	RIGHTPARENS	')'
#define	LEFTPARENS	'('
#define MAXSTRING	5000		/* maximum length of string (INCREASE LATER) */
#define IsItAlphaNum(c)	  ( (c) >= 48 && (c)<=57 ) || ( (c) >= 65 && (c)<=90 ) || ( (c) >= 97 && (c)<=122 ) 
#define min(a,b)			( (a)<=(b) )  ? (a):(b)
#define max(a,b)			( (a)>=(b) )  ? (a):(b)
#define SIBLOOP(c)			for (; (c); (c)=(c)->sib)
#define isTip(c)			 ( (c)->firstdesc == NULL )
#define isConstrained(c)		( (c)->nodeIsConstrainedMax || (c)->nodeIsConstrainedMin)
#define isConstrainedMax(c)		( (c)->nodeIsConstrainedMax)
#define isConstrainedMin(c)		( (c)->nodeIsConstrainedMin)
#define isFree(c)			( (c)->free == 1)
#define isFixed(c)			( (c)->free == 0)
#define isRoot(c)			( (c)->anc == NULL )
#define FLAGFLIP(c)			( (c) = (c) ^ 1)  /* use XOR operator */
#define MIN(a,b) (((a)<(b))?(a):(b))
#define LN2 0.69314718
#define LOG2(x) (log((double)(x))/LN2)  /*base 2 logarithm */

#define MAXCLADES 200 /* also defined in ReadNexusFile2 */
#define MAXTAX	  100

/* STRUCTURES AND PROTOTYPES */

/**********************************************/

/*  Node structure */

struct nodetype {
				struct nodetype 		*anc;
				struct nodetype 		*firstdesc;
				struct nodetype 		*sib;
				struct nodetype			*nodePtr;	/* generic pointer to some other node */
				char 				*taxon_name;
				double				length;		/* length of subtending branch */
				int 				order;
				long 				numdesc;
				int				numSelectedDesc;/* number of selected nodes 
									    below this one (including this one) */
				long				id;	
				int 				X;		/* positions on screen */
				int 				Y;
				double				time;		/* current time of node... */
				double				nodeReal;	/* Let's use this for various real numbers */
				short				nodeFlag;	/* for various flags */
				double				estRate;	/* estimated rate, usually for branch */
				double 				nodeEstRate;	/* estimated rate,special for node method */

				int				isQueryNode;	/* 1 if node to be used in query; USED IN NODE MARKING ROUTINES */
				
				int				isCompactNode;	/* 1 if this node is displayed
										as a clade of all its descendants*/
				int				toggleDesc;	/* 1 if all descendants should
										be queried */
				int				nodeIsConstrainedMax;
				int 				nodeIsConstrainedMin;
  				int				modelID;	/* takes integer values for different rate parms under 
										local clock model */
				int				free;	/* 1 if we estimate this node's time */
				double				minAge;		/* present = 0; root = 1 */
				double				maxAge;		/* ...These are constraints on ages*/
					
				double				cumulProb;	/* Used in RandomTree modules */
				double				like;		/* likelihood of subtending branch */
				double				chiSq;


				double				CL[4];		/* conditional likelihood for four states */
				double				CLmax[4];	/* Pupko max of 4 conditional likelihoods*/
				double				CLmarg[4];	/* Marginal likelihoods computed by rerooting */
				int				CLopt[4];	/* Pupko optimal state choice */
				int				opt;		/* optimal state at this node */

#if LARGE_NODES
                                double beta[2];
                                double beta_sum;
                                double delta[2];
                                int psi[2];
#endif
				};
typedef 	struct nodetype NODETYPE;
typedef		struct nodetype * NODE;

/**********************************************/

/*  Tree structure */

struct treetype 	{
			char *		name;
			char *		TD;
			long		numTaxa;
			long		numBranches;
			long		basalTomy;
			NODETYPE 	*root;
			PtrList		cladeSet;
			double 		est_b;
			double		est_c;	/*estimates of gamma rate parms */
			double		estRate;
			int		timesAvailable;	/* 1 if times have been estimated */
			int		method;		/* method used to estimate times */
			double		obj;		/* value of obj func at soln */
			};
typedef 	struct treetype  * TREE;




/*************************************************/

void TreeToNodePtrList(NODETYPE *node,  PtrList NodeList)	;
double ** tree2VCV(TREE t, int i);
void print_named_ages(NODETYPE *node);
double pathLengthTime(NODE a, NODE b);
double pathLengthTimeModel(NODE a, NODE b, int model);
NODE  mrca(NODE a, NODE b );
TREE Subtree_Initialize(TREE T, NODETYPE *node);
long numUnMarkedDesc(NODETYPE *node);
NODETYPE * nextRndNode(long nNodes, NODETYPE ** nodeArray);
void markNode( NODETYPE  * n);
void unMarkNode( NODETYPE  * n);
int isNodeMarked( NODETYPE  *n);
void unMarkTree(TREE T);

void setLocalModel(NODETYPE *n,int model,int stemFlag);
void summarize_rates(TREE t);
static void recurse_summarize_rates(NODETYPE * n, long * ix, double r[]);

void scaleTree(NODETYPE * root, double calAge, NODETYPE * calNode);

void setupConstrainedNodeArray(NODETYPE * node, NODE nodeArray[], int *curIndex);
void setupFixedNodeArray(NODETYPE * node, NODE nodeArray[], int *curIndex);
int numConstrainedNodes(NODETYPE * node);
int numFixedNodes(NODETYPE * node);
double unFixNodeAge(NODETYPE *node);
double fixNodeAge(NODETYPE *node);
void print_rates(NODETYPE *n,int method);
void RemoveTaxonLvAnc(NODETYPE * n);
int collapseLengthsTree2Tree(TREE t1,TREE t2);
void setNodeEstRate(NODE node);
void zeroEstRate(NODETYPE *node);
double cvSquareError(TREE t, int method);
double cvSquareErrorBranch(TREE t, NODE n,int method,double *chiSq);


void copyLengthsTree2Tree(NODETYPE * node1,NODETYPE * node2);
double LFunconsT(NODETYPE *node);
double LFuncons1T(NODETYPE *node);
double LFuncons(NODETYPE *node);
double LFuncons1(NODETYPE *node);
void preOrderVoid(NODETYPE *node,void (*f)(NODETYPE *));
double preOrderArg(NODETYPE *node,double (*func)(NODE node, double farg),double farg);
double preOrder(NODETYPE *node,double (*f)(NODETYPE *));
void unSetConstraints(NODETYPE * node);
int maxAgePresent(NODETYPE * node);
int constraintsPresent(NODETYPE * node);

int tipsDifferentAges(NODETYPE *node);

void		ABCSuperTreePurvis(NODETYPE *node, StrListPtr UniqueList, 
		    char **dataMatrix,float *wtset)	;
void		ABCSuperTree(NODETYPE *node, StrListPtr UniqueList, 
		    char **dataMatrix,float *wtset);	
void		TreeToTaxaList(NODETYPE *node,  StrListPtr taxaList);	
void		TreeToTaxaPtrList(NODETYPE *node,  PtrList NodeList);	
int			 maxorder(NODETYPE *node);
int 			numdesc(NODETYPE *),
			stringcheck(char *td);
			

int numFreeNodes(NODETYPE *node);
NODETYPE	* sister(NODETYPE * n);
NODETYPE 	* ReRoot(NODETYPE * atNode);
NODETYPE 	* ReRoot2(NODETYPE * atNode);
void 		Flip(NODETYPE *a);
NODETYPE * 		RemoveTaxon(TREE t,NODETYPE * theChild);
void 		AddChild(NODETYPE * parent, NODETYPE * theChild);

void		 updateSubtrees(NODETYPE *srcNode);
NODETYPE	 *createSubtree(NODETYPE *srcNode, int SubtreeSize);
void 		copyNodeInfo(NODETYPE *source,NODETYPE *dest);
double 		*sort_node_time_array(NODETYPE *root);
double 		get_sum_durations(NODETYPE *node);
NODETYPE	*newnode(void),
		*makegroup(void),
		*string_to_tree(char *tree_description);
void 		collapse_node(NODETYPE *node);			
void 		collapse_zero(NODETYPE *node);
void 		Node_Destructor(NODETYPE *node);
void 		make_parens(NODETYPE *root, int flag);
void 		Tree_Destructor(TREE aTree);
void		DisposeTree(NODETYPE *node);
int 		numIntNodes(NODETYPE *node);
int 		numBranches(NODETYPE *node);
void 		TreeStats(NODETYPE *root, int * numtips, 
				int * numinternals, int * roottomy);
void		setNodeName(NODETYPE *node, char *name);
void 		printtree(NODETYPE *node);
NODETYPE * 	find_taxon_name(NODETYPE *node,char *name);
void 		Tree_Initialize(TREE aTree, char *TD, char *name);
void 		print_tree_list(PtrList treeList);
void		print_ages(NODETYPE *node, double time, double calAge,int method);
void		init_node_ids(NODETYPE *node, int gId);
void		convert_branchlength_to_time(NODETYPE *root);
NODETYPE	*find_id(NODETYPE *node,int id);
int		node_tomy(NODETYPE *node);
int		isNodeDescendant(NODETYPE *nodeA, 
		    NODETYPE *nodeB);
NODETYPE *	MRCA(NODETYPE *, StrListPtr taxaList);
void		traverseMultiplyLength(NODETYPE *, double x,int round);
void 		Tree2CladeSets(NODETYPE *node, StrListPtr allTaxaList, int nTaxa, 
		    PtrList SetList);	
void 		printCladeSets(PtrList SetList);
PtrList		Tree2CladeSet(TREE thisTree, StrListPtr allTaxaList);
int 		group_a_clade(NODETYPE *root, StrListPtr taxaList);
int 		any_zero_internal_branches(NODETYPE *node);
int             any_zero_terminal_branches(NODETYPE *node);
void 		printLikes(NODETYPE *node);
void 		printnodeLike(NODETYPE *node);
void 		ClusterHistogram(NODETYPE * node, long * array,long TSize);
double 		treeDurLength(NODETYPE * node);
double 		treeLength(NODETYPE * node);
double 		treeAgeSum(NODETYPE * node);
int 		compar(const void *v1,  const void *v2);
int 		numNodes(NODETYPE *node);
void 		init_free(NODETYPE *node);
void rootToTips(NODETYPE* node,double curLen);
NODETYPE*  copyTree(NODETYPE* a);
#endif
                                            r8s/WuLi.c                                                                                          0000644 0000766 0000120 00000023750 07565204055 012306  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "moment.h"
#include "nexus.h"
#include "WuLi.h"
#include "MyUtilities.h"
#include "memory.h"
#include "DistrFuncs.h"
#include "distance.h"
#define SMALL	0.00001
double Sqr(double);


/*********************************************************************/
void WuLiStub(int inGroup1,int inGroup2, int outGroup)
{
	extern FILE * outstream;
	int* excArray; /* a local exclusion array.  NOTE THAT routines in distance, and WuLi use the
		local and globabl exclusion array inconsistently.  CHeck this when redoing bootstrap*/
	int bs;
	extern struct NexDataType *gNexDataPtr;
	extern long iix;
	double zscore[3]={1.96,2.57,3.30}; /* cutoffs for P=0.05,0.01, and 0.001 in
					a two-tailed z-test */
	int 	i,j,id[4],ih,error,errorbs=0,errorsdevbs=0,errorsdev=0;
	double	z,zbs,delta,stddev,*data,dif,difbs,Poa,Pob,doa,dob;
	double	mean,adev,sdev,svar,skew,curt;
	long	bufPos=0,ix;
	long	N;
	long	NChars,actualNChars=0;
	int	*saveExcArray;	/* saves the current exclusion array while bootstrapping */


	NChars = gNexDataPtr->NChars;
	excArray=gNexDataPtr->excArray;
	
/**** Following is code to implement the bootstrap estimate of variance ****/	

	bs=gNexDataPtr->RateBlockParms.isBS;	/* just a flag */
	N=gNexDataPtr->RateBlockParms.NReps;
	if (bs && (NChars >0) && (N >0) )
		{
		saveExcArray=(int *)myMalloc(NChars*sizeof(int));
		for (ix=0;ix<NChars;ix++)
			{
			saveExcArray[ix]=excArray[ix];
			if(excArray[ix]>0)
				++actualNChars;/* count the number of included chars */
			}
		iix = gNexDataPtr->RateBlockParms.seed;
		data = (double *)myMalloc(N*sizeof(double));
		
		for(i=0;i<N;i++)  /* Do the resampling and store results of each rep */
			{
			bshuf(excArray,saveExcArray,NChars,actualNChars);
/************>>>>>>*/
			error=WuLiTest(inGroup1,inGroup2,outGroup,&difbs,&stddev,
				       &Poa,&Pob,&doa,&dob);
			if (error)
				errorbs=1;	/* set this if any error among bs reps */
			data[i]=difbs;
			}
		moment(data-1,(int)N,&mean,&adev,&sdev,&svar,&skew,&curt);
		difbs=mean;	/* these are "bias corrected" bootstrap differences. 
					Am I sure I want this ? */

	/*** Error handling ***/
		if (sdev==0.0)
			errorsdevbs=1;
		if (!errorsdevbs && !errorbs)
			zbs=difbs/sdev;
	/********/	
		for (ix=0;ix<NChars;ix++)
			excArray[ix]=saveExcArray[ix];	
				/* returns exclusion set array to original value*/
		myFree(data);
		myFree(saveExcArray);
		}
/***************************** Non bootstrap code ******************/		

	error=WuLiTest(inGroup1,inGroup2,outGroup,&dif,&stddev,&Poa,&Pob,&doa,&dob);

	/*** Error handling ***/
		if (stddev==0.0)
			errorsdev=1;
		if (!errorsdev && !error)
			z=dif/stddev;
	/**********/


	printf("(%15.15s (%15.15s %15.15s))\t\t",
		getkthStr(gNexDataPtr->TaxaList,outGroup),
		getkthStr(gNexDataPtr->TaxaList,inGroup1),
		getkthStr(gNexDataPtr->TaxaList,inGroup2)
		);



	if (error || errorsdev)
		printf("Error%1i%1i\t",error,errorsdev);
	else
		{
	/*	printf("%+6.4f",z);
		for (i=0;i<3;i++) if (fabs(z) < zscore[i]) break; 
		for (j=0;j<i;j++) printf("*");
		for (j=i;j<3;j++) printf(" ");
	*/	}
	if (bs)
		{
		if (errorbs || errorsdevbs)
			printf("Error%1i%1i\t",errorbs,errorsdevbs);
		else
			{
		/*	printf("\t%+6.4f",zbs);
			for (i=0;i<3;i++) if (fabs(zbs) < zscore[i]) break;
			for (j=0;j<i;j++) printf("*");
			for (j=i;j<3;j++) printf(" ");
			printf("[% +6.4f % +6.4f] (% 6.4f % 6.4f) ",dif,difbs,stddev,sdev);*/
                        
			zbs=dif/sdev; /* this is the z value using observed difference and bs
					estimate of std dev ! */
			printf("\nRR: P(oa)=%f P(ob)= %f d(oa)=%f d(ob)=%f z=%6.4f",Poa,Pob,doa,dob,zbs);
                        for (i=0;i<3;i++) if (fabs(zbs) < zscore[i]) break; 
                        for (j=0;j<i;j++) printf("*");
			printf("\n");
			}
		}
		
	printf("\n");



/* free(excArray);*/  /* can't be constantly creating and destroying this array!! */
return;
}


/*********************************************************************/
int WuLiTest(int inGroup1,int inGroup2, int outGroup,
		 double *dif, double *stddev,double *Poa, double *Pob, double *doa, double *dob)

/* Performs the relative rate test of Wu and Li (1985) as described in Muse and Weir 
(1992), Genetics 132:269-276.  

Requires a NEXUS file stored in Character Block/Tree Block format (this can be enforced
in the NEXUS options in MacClade.  Also, the data matrix must NOT contain explicit
polymorphisms, e.g., sites stored as "{a/c}".  Use the IUPAC codes instead.  Polymorphic
sites, question marks, and sites with gaps are ignored in the pairwise distance calculations
in PQCalc().
You MAY use the dot format for sequences that are the same as the first line of data.

Usage: Following the Data block in the NEXUS file add the following block:

	begin rates;
		wuli  ingroup1 ingroup2 outgroup;
		......
	end;
	
The ingroup and outgroups may be NEXUS taxon numbers or 
names (exactly as they appear in the taxa block).  The outgroup MUST BE LAST.
Any number of 'wuli' commands may appear in the block.

*/



/* all id numbers passed as parameters are 1..ntaxa, ie. unit offset */

{
	
	extern struct NexDataType *gNexDataPtr;
	int 	i,j,id[4],ih,error=0,RRtype;
	long	seqLength,nP,nQ,mA,mB;
	double	z,num,denom,dij,dijk,ssqr,delta,SQdif;
	double p0a,p0b,d1,d2;
	char *pi, *pj, *pRow1;
	const float K = 20.0/19;
	double	P[4][4],
			Q[4][4],
			a[4][4],
			b[4][4],
			c[4][4],
			d[4][4],
			Ahat[4][4],
			Bhat[4][4],
			VarK[4][4];
			


	id[1]=inGroup1;
	id[2]=inGroup2;
	id[3]=outGroup;
	
	RRtype=gNexDataPtr->RateBlockParms.RRtype;	/* just a flag */

	pRow1=getkthStr(gNexDataPtr->DMList,1);	
	for (i=1;i<=2;i++)
		for (j=i+1;j<=3;j++)
			{
			pi=getkthStr(gNexDataPtr->DMList,id[i]);
			pj=getkthStr(gNexDataPtr->DMList,id[j]);
			if (RRtype==MIKE)
				seqLength=aaCalc1(pi,pj,pRow1,&P[i][j],&nP);	
			else
				seqLength=PQCalc1(pi,pj,pRow1,&P[i][j],&Q[i][j],&nP,&nQ);
			
			if (seqLength <= 1) 
				{
				return 1;   /* Need two or more valid sites in all sequence comparisons */
				}

			if (RRtype == WULI)
				{
				if ((2*Q[i][j] >= 1.0) || (2*P[i][j]+Q[i][j] >= 1.0) )
					{
					return 2;  /* Divergence too large for Wu Li to handle */
					}
				a[i][j]=1.0/(1-2*P[i][j]-Q[i][j]);
				b[i][j]=1.0/(1-2*Q[i][j]);
				d[i][j]=0.5*(a[i][j]+b[i][j]);
				Ahat[i][j]=0.5*log(a[i][j])-0.25*log(b[i][j]);
				Bhat[i][j]=0.5*log(b[i][j]);
				VarK[i][j]= (Sqr(a[i][j])*P[i][j]+Sqr(d[i][j])*Q[i][j]
					-Sqr(a[i][j]*P[i][j]+d[i][j]*Q[i][j]))/seqLength;
				}
			}

	if (RRtype == WULI)
		{
		if (P[1][2]+Q[1][2] < SMALL)
			return 9;  /* Divergence too small -- can cause negative arguments to variance
							calculations if missing data are a problem */



		Bhat[0][3]=0.5*(Bhat[1][3]+Bhat[2][3]-Bhat[1][2]);
		Ahat[0][3]=0.5*(Ahat[1][3]+Ahat[2][3]-Ahat[1][2]);
		Q[0][3]=0.5*(1-exp(-2*Bhat[0][3]));
		P[0][3]=0.5*(1-Q[0][3]-exp(-2*Ahat[0][3]-Bhat[0][3]));
		a[0][3]=1.0/(1-2*P[0][3]-Q[0][3]);
		b[0][3]=1.0/(1-2*Q[0][3]);
		d[0][3]=0.5*(a[0][3]+b[0][3]);
		VarK[0][3]=(Sqr(a[0][3])*P[0][3]+Sqr(d[0][3])*Q[0][3]
				-Sqr(a[0][3]*P[0][3]+d[0][3]*Q[0][3]))/seqLength;
		*dif=Ahat[1][3]+Bhat[1][3]-(Ahat[2][3]+Bhat[2][3]);
		SQdif=VarK[1][3]+VarK[2][3]-2*VarK[0][3];
		if (SQdif < 0.0)
			{
			error = 3;	/* Apparently this can happen when d12 is zero (or small)
				but d13 and d23 are not equal because of missing data */
			*stddev=0.0;	/* just set it to zero */
			}
		else
			*stddev=sqrt(SQdif);
		}

	  if (RRtype == STEEL)	/* Next is the Steel et al. nonparametric method */
		{
		(void)tripletSites(inGroup1,inGroup2,outGroup,&dijk,&mA,&mB);
		*dif=P[1][3]+Q[1][3]-(P[2][3]+Q[2][3]);
		dij=P[1][2]+Q[1][2];
		SQdif=(dij-Sqr(*dif)-dijk)/(seqLength-1);
		if (SQdif < 0.0)
			{
			error = 4+error; /* See above... */
			*stddev=0.0;	/* just set it to zero */
			}
		else
			*stddev=sqrt(SQdif);
		}

	  if (RRtype == MIKE)	/* just doodling...didn't work */ 
		{
		if (P[1][3]>=0.95 || P[2][3]>=0.95)
			doGenericAlert("Argument out of bounds in AA distance");
		else
			{
			d1=-0.95*log(1-K*P[1][3]);
			d2=-0.95*log(1-K*P[2][3]);
			/*
			d1=1.6*(pow(1-P[1][3],-1/1.6)-1);
			d2=1.6*(pow(1-P[2][3],-1/1.6)-1);
			*/
			*dif = d1-d2;
			*doa=d1;
			*dob=d2;
			*Poa=P[1][3];
			*Pob=P[2][3];
/*
printf("P=%f\td=%f\n",P[1][3],d1);
printf("P=%f\td=%f\n",P[2][3],d2);
*/
			}	
		}
	  if (RRtype == TAJIMA)	
		{
		(void)tripletSites(inGroup1,inGroup2,outGroup,&dijk,&mA,&mB);
		*dif=Sqr(mA-mB)/(mA+mB);
		*stddev=1.0;
		}
return error;	/* error free return */
	
}






long tripletSites(int i, int j, int k, double *P, long *MA, long *MB)

/*  Returns the number of valid sites (nongap, nonmissing) in sequences
that are different in each of the three taxa, where i and j and k are NEXUS 
taxon numbers.  Used by the Steel test. Stores P which is the proportion  */

{

	extern struct NexDataType *gNexDataPtr;
	int* excArray;		/* array for exclusion set */
	long isite, validcount=0,difcount=0,mA=0,mB=0;
	char *pi, *pj, *pk,*pRow1, ci,cj,ck,missing,gap,match;
	
	excArray=gNexDataPtr->excArray;
	gap=gNexDataPtr->gapchar;
	missing=gNexDataPtr->missingchar;
	match=gNexDataPtr->matchchar;

	pi=getkthStr(gNexDataPtr->DMList,i);
	pj=getkthStr(gNexDataPtr->DMList,j);
	pk=getkthStr(gNexDataPtr->DMList,k);
	pRow1=getkthStr(gNexDataPtr->DMList,1);	
	for (isite=0;isite<gNexDataPtr->NChars;isite++)
	  if (excArray[isite])  /* process site only if not in exclusion set */
		{
		ci=pi[isite];
		cj=pj[isite];
		ck=pk[isite];
		if (ci==match) ci = pRow1[isite];  /* check for 'period' format in sequences */
		if (cj==match) cj = pRow1[isite];  /* if present, set to data for first row */
		if (ck==match) ck = pRow1[isite];  /* if present, set to data for first row */
		
		if (  strchr("ACGT",ci) &&  strchr("ACGT",cj) &&  strchr("ACGT",ck) ) 
				/* only consider when the three sites are in ACGT */
			{
			++validcount;
			if ((ci != cj )	&& (ci != ck) && (ck != cj)  )
				++difcount;

			if( ( cj == ck ) && (cj != ci) )
				++mA;
			if( ( ci == ck ) && (cj != ci) )
				++mB;	/* these are for TAJIMA's method */

			}	
		}
	if (validcount > 0) 
		*P=(double)difcount/validcount;
	*MA=mA;*MB=mB;
	return validcount;

}

double Sqr(double x)
{
return x*x;
}
                        r8s/WuLi.h                                                                                          0000644 0000766 0000120 00000000510 07565204055 012300  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #define WULI	0
#define STEEL	1
#define TAJIMA	2
#define MIKE	3
int WuLiTest(int inGroup1,int inGroup2, int outGroup,
		 double *dif, double *stddev,double *Poa, double *Pob, double *doa, double *dob);
long tripletSites(int i, int j, int k, double *P, long *MA, long *MB);
void WuLiStub(int inGroup1,int inGroup2, int outGroup);
                                                                                                                                                                                        r8s/ancestral.c                                                                                     0000644 0000766 0000120 00000007207 10414517145 013373  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /*

Module to reconstruct ancestral states of continuous traits via squared change parsimony.
Traits are stored in node->nodeReal

Contains several routines from NRC and some minor modifications to same...

*/
#include "ObjFunc.h"
#include "Maximize.h"
#include "structures.h"
#include "nexus.h"
#include "ancestral.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SQR(x)        ((x)*(x))

static void pTimeArray2treeAncestral(NODETYPE * node, double lp[]);
double recurseAncestral(NODETYPE *node);



double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);

int gNParm, gNT, p2tindex; 

void ancestralOptimize(TREE t,int *numIter, double ftol,double linMinDelta,int *success )
{
extern struct NexDataType *gNexDataPtr;
extern NODETYPE * gRoot;
StrListPtr DM, TL;
PtrList nodeList;
int nParm,NT,i,j, ixTL;
double *p, obj,tip_value;
char * tip_name,*tip_value_str, *dummy,*found_tip;
NODE a;
float meanTip=0.0;

DM=gNexDataPtr->DMList;
TL=gNexDataPtr->TaxaList;
gRoot = t->root;
nodeList=pNewList();
TreeToTaxaPtrList(t->root,nodeList);

nParm=numIntNodes(t->root);
NT=t->numTaxa;
p = dvector(1,nParm);
for (i=1;i<=NT;i++)
	{
	a=(NODE)(pListgetkthNode(nodeList, i)->item);
	tip_name=a->taxon_name;
	ixTL=findMatchStr(TL,tip_name); // lookup the taxon name from the matrix ordering and get the relevant data matrix ordering corresponding name
	tip_value_str=getkthStr(DM,(long)(ixTL));
	found_tip=getkthStr(TL,(long)(ixTL));
	tip_value=strtod(tip_value_str,NULL);
// printf("%i\t%s\t%s\t%s\t%f\n",i,tip_name,found_tip,tip_value_str,tip_value);
	a->nodeReal=tip_value;
	meanTip+=tip_value;
	}
meanTip/=NT;
for (i=1;i<=nParm;i++)
	p[i]=meanTip; // use mean tip value as guess for all nodes...


obj=MinND(t,0, POWELL,objAncestral,NULL,p, nParm,numIter, ftol,linMinDelta,success );



free_dvector(p,1,nParm);
}




/***********************************************************************************/

double objAncestral(double p[])


{
  extern struct NexDataType *gNexDataPtr;	
  extern NODETYPE * gRoot;    /* This global is declared when the whole algorithm is called */

  static int firstTime=1,num_branches;
  double obj;
  NODETYPE *child;
  
  p2tindex=1;
  pTimeArray2treeAncestral(gRoot,p);  // puts the array values onto the tree 


/*** Now find objective function over rest of tree ***/

  obj=recurseAncestral(gRoot);

  return obj; 
}
/**********************/

double recurseAncestral(NODETYPE *node)
{
  NODETYPE *child;
  double obj=0.0;
  if (!node) return(0.0);
  child=node->firstdesc;
  SIBLOOP(child)
    {
      obj += SQR(child->nodeReal-node->nodeReal);
      obj += recurseAncestral(child);
    }
  return obj;	
}

static void pTimeArray2treeAncestral(NODETYPE * node, double lp[])
{
NODETYPE *child;
if (!isTip(node))
	node->nodeReal=lp[p2tindex++];
child=node->firstdesc;
SIBLOOP(child)
	pTimeArray2treeAncestral(child,lp);
return;
}
/***********************************************************************************/
void printAncestral(NODETYPE *node)
{
    float diff;
    NODETYPE *child;

    if (!isTip(node))
      {
      if(*(node->taxon_name))
	    printf(" [*]  %.8s\t", node->taxon_name);
      else
	    printf("     (%i)\t", node->id);
      }
    else
	printf("      %.8s\t",node->taxon_name);

//    printf("%7.2f\t\t",node->time);
    if (!isRoot(node))
	{
	diff = node->nodeReal - node->anc->nodeReal;
	printf("%.4f\t\t%.4f\t\t%.4f",node->nodeReal,node->anc->nodeReal,diff);
	if (diff >= 0.0) printf("\t\t+\n");
	else printf("\t\t-\n");
	}
     else
	printf("%.4f\n",node->nodeReal);



    if(!isTip(node))
	{
    	child=node->firstdesc;
    	SIBLOOP(child) 
		printAncestral(child);
	}
    return;
    
}
                                                                                                                                                                                                                                                                                                                                                                                         r8s/ancestral.h                                                                                     0000644 0000766 0000120 00000000242 10343376052 013371  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  void printAncestral(NODETYPE *node);
void ancestralOptimize(TREE t,int *numIter, double ftol,double linMinDelta,int *success );
double objAncestral(double p[]);

                                                                                                                                                                                                                                                                                                                                                              r8s/blas.f                                                                                          0000644 0000766 0000120 00000021606 07565204055 012350  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  C%% TRUNCATED-NEWTON METHOD: BLAS
C   NOTE: ALL ROUTINES HERE ARE FROM LINPACK WITH THE EXCEPTION
C         OF DXPY (A VERSION OF DAXPY WITH A=1.0)
C   WRITTEN BY:  STEPHEN G. NASH
C                OPERATIONS RESEARCH AND APPLIED STATISTICS DEPT.
C                GEORGE MASON UNIVERSITY
C                FAIRFAX, VA 22030
C****************************************************************
      DOUBLE PRECISION FUNCTION DDOT(N,DX,INCX,DY,INCY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     FORMS THE DOT PRODUCT OF TWO VECTORS.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DTEMP
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      DDOT = 0.0D0
      DTEMP = 0.0D0
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DTEMP = DTEMP + DX(IX)*DY(IY)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      DDOT = DTEMP
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,5)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DTEMP = DTEMP + DX(I)*DY(I)
   30 CONTINUE
      IF( N .LT. 5 ) GO TO 60
   40 MP1 = M + 1
      DO 50 I = MP1,N,5
        DTEMP = DTEMP + DX(I)*DY(I) + DX(I + 1)*DY(I + 1) +
     *   DX(I + 2)*DY(I + 2) + DX(I + 3)*DY(I + 3) + DX(I + 4)*DY(I + 4)
   50 CONTINUE
   60 DDOT = DTEMP
      RETURN
      END
      SUBROUTINE DAXPY(N,DA,DX,INCX,DY,INCY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     CONSTANT TIMES A VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1),DA
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF (DA .EQ. 0.0D0) RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DA*DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DA*DX(I)
        DY(I + 1) = DY(I + 1) + DA*DX(I + 1)
        DY(I + 2) = DY(I + 2) + DA*DX(I + 2)
        DY(I + 3) = DY(I + 3) + DA*DX(I + 3)
   50 CONTINUE
      RETURN
      END
      DOUBLE PRECISION FUNCTION DNRM2 ( N, DX, INCX)
      IMPLICIT         DOUBLE PRECISION (A-H,O-Z)
      INTEGER          NEXT
      DOUBLE PRECISION DX(1),CUTLO,CUTHI,HITEST,SUM,XMAX,ZERO,ONE
      DATA   ZERO, ONE /0.0D0, 1.0D0/
C
C     EUCLIDEAN NORM OF THE N-VECTOR STORED IN DX() WITH STORAGE
C     INCREMENT INCX .
C     IF    N .LE. 0 RETURN WITH RESULT = 0.
C     IF N .GE. 1 THEN INCX MUST BE .GE. 1
C
C           C.L.LAWSON, 1978 JAN 08
C
C     FOUR PHASE METHOD     USING TWO BUILT-IN CONSTANTS THAT ARE
C     HOPEFULLY APPLICABLE TO ALL MACHINES.
C         CUTLO = MAXIMUM OF  DSQRT(U/EPS)  OVER ALL KNOWN MACHINES.
C         CUTHI = MINIMUM OF  DSQRT(V)      OVER ALL KNOWN MACHINES.
C     WHERE
C         EPS = SMALLEST NO. SUCH THAT EPS + 1. .GT. 1.
C         U   = SMALLEST POSITIVE NO.   (UNDERFLOW LIMIT)
C         V   = LARGEST  NO.            (OVERFLOW  LIMIT)
C
C     BRIEF OUTLINE OF ALGORITHM..
C
C     PHASE 1    SCANS ZERO COMPONENTS.
C     MOVE TO PHASE 2 WHEN A COMPONENT IS NONZERO AND .LE. CUTLO
C     MOVE TO PHASE 3 WHEN A COMPONENT IS .GT. CUTLO
C     MOVE TO PHASE 4 WHEN A COMPONENT IS .GE. CUTHI/M
C     WHERE M = N FOR X() REAL AND M = 2*N FOR COMPLEX.
C
C     VALUES FOR CUTLO AND CUTHI..
C     FROM THE ENVIRONMENTAL PARAMETERS LISTED IN THE IMSL CONVERTER
C     DOCUMENT THE LIMITING VALUES ARE AS FOLLOWS..
C     CUTLO, S.P.   U/EPS = 2**(-102) FOR  HONEYWELL.  CLOSE SECONDS ARE
C                   UNIVAC AND DEC AT 2**(-103)
C                   THUS CUTLO = 2**(-51) = 4.44089E-16
C     CUTHI, S.P.   V = 2**127 FOR UNIVAC, HONEYWELL, AND DEC.
C                   THUS CUTHI = 2**(63.5) = 1.30438E19
C     CUTLO, D.P.   U/EPS = 2**(-67) FOR HONEYWELL AND DEC.
C                   THUS CUTLO = 2**(-33.5) = 8.23181D-11
C     CUTHI, D.P.   SAME AS S.P.  CUTHI = 1.30438D19
C     DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C     DATA CUTLO, CUTHI / 4.441E-16,  1.304E19 /
      DATA CUTLO, CUTHI / 8.232D-11,  1.304D19 /
C
      IF(N .GT. 0) GO TO 10
         DNRM2  = ZERO
         GO TO 300
C
   10 ASSIGN 30 TO NEXT
      SUM = ZERO
      NN = N * INCX
C                                                 BEGIN MAIN LOOP
      I = 1
   20    GO TO NEXT,(30, 50, 70, 110)
   30 IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
      ASSIGN 50 TO NEXT
      XMAX = ZERO
C
C                        PHASE 1.  SUM IS ZERO
C
   50 IF( DX(I) .EQ. ZERO) GO TO 200
      IF( DABS(DX(I)) .GT. CUTLO) GO TO 85
C
C                                PREPARE FOR PHASE 2.
      ASSIGN 70 TO NEXT
      GO TO 105
C
C                                PREPARE FOR PHASE 4.
C
  100 I = J
      ASSIGN 110 TO NEXT
      SUM = (SUM / DX(I)) / DX(I)
  105 XMAX = DABS(DX(I))
      GO TO 115
C
C                   PHASE 2.  SUM IS SMALL.
C                             SCALE TO AVOID DESTRUCTIVE UNDERFLOW.
C
   70 IF( DABS(DX(I)) .GT. CUTLO ) GO TO 75
C
C                     COMMON CODE FOR PHASES 2 AND 4.
C                     IN PHASE 4 SUM IS LARGE.  SCALE TO AVOID OVERFLOW.
C
  110 IF( DABS(DX(I)) .LE. XMAX ) GO TO 115
         SUM = ONE + SUM * (XMAX / DX(I))**2
         XMAX = DABS(DX(I))
         GO TO 200
C
  115 SUM = SUM + (DX(I)/XMAX)**2
      GO TO 200
C
C
C                  PREPARE FOR PHASE 3.
C
   75 SUM = (SUM * XMAX) * XMAX
C
C
C     FOR REAL OR D.P. SET HITEST = CUTHI/N
C     FOR COMPLEX      SET HITEST = CUTHI/(2*N)
C
   85 HITEST = CUTHI/FLOAT( N )
C
C                   PHASE 3.  SUM IS MID-RANGE.  NO SCALING.
C
      DO 95 J =I,NN,INCX
      IF(DABS(DX(J)) .GE. HITEST) GO TO 100
   95    SUM = SUM + DX(J)**2
      DNRM2 = DSQRT( SUM )
      GO TO 300
C
  200 CONTINUE
      I = I + INCX
      IF ( I .LE. NN ) GO TO 20
C
C              END OF MAIN LOOP.
C
C              COMPUTE SQUARE ROOT AND ADJUST FOR SCALING.
C
      DNRM2 = XMAX * DSQRT(SUM)
  300 CONTINUE
      RETURN
      END
      SUBROUTINE DCOPY(N,DX,INCX,DY,INCY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     COPIES A VECTOR, X, TO A VECTOR, Y.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     JACK DONGARRA, LINPACK, 3/11/78.
C
      DOUBLE PRECISION DX(1),DY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,7)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DX(I)
   30 CONTINUE
      IF( N .LT. 7 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,7
        DY(I) = DX(I)
        DY(I + 1) = DX(I + 1)
        DY(I + 2) = DX(I + 2)
        DY(I + 3) = DX(I + 3)
        DY(I + 4) = DX(I + 4)
        DY(I + 5) = DX(I + 5)
        DY(I + 6) = DX(I + 6)
   50 CONTINUE
      RETURN
      END
C******************************************************************
C SPECIAL BLAS FOR Y = X+Y
C******************************************************************
      SUBROUTINE DXPY(N,DX,INCX,DY,INCY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
C
C     VECTOR PLUS A VECTOR.
C     USES UNROLLED LOOPS FOR INCREMENTS EQUAL TO ONE.
C     STEPHEN G. NASH 5/30/89.
C
      DOUBLE PRECISION DX(1),DY(1)
      INTEGER I,INCX,INCY,IX,IY,M,MP1,N
C
      IF(N.LE.0)RETURN
      IF(INCX.EQ.1.AND.INCY.EQ.1)GO TO 20
C
C        CODE FOR UNEQUAL INCREMENTS OR EQUAL INCREMENTS
C          NOT EQUAL TO 1
C
      IX = 1
      IY = 1
      IF(INCX.LT.0)IX = (-N+1)*INCX + 1
      IF(INCY.LT.0)IY = (-N+1)*INCY + 1
      DO 10 I = 1,N
        DY(IY) = DY(IY) + DX(IX)
        IX = IX + INCX
        IY = IY + INCY
   10 CONTINUE
      RETURN
C
C        CODE FOR BOTH INCREMENTS EQUAL TO 1
C
C
C        CLEAN-UP LOOP
C
   20 M = MOD(N,4)
      IF( M .EQ. 0 ) GO TO 40
      DO 30 I = 1,M
        DY(I) = DY(I) + DX(I)
   30 CONTINUE
      IF( N .LT. 4 ) RETURN
   40 MP1 = M + 1
      DO 50 I = MP1,N,4
        DY(I) = DY(I) + DX(I)
        DY(I + 1) = DY(I + 1) + DX(I + 1)
        DY(I + 2) = DY(I + 2) + DX(I + 2)
        DY(I + 3) = DY(I + 3) + DX(I + 3)
   50 CONTINUE
      RETURN
      END
                                                                                                                          r8s/continuousML.c                                                                                  0000644 0000766 0000120 00000023761 10357556225 014071  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /*

Module to implement ML estimation of rate parameter(s) for a multivariate normal
model of continuous trait evolution.


Contains several routines from NRC and some minor modifications to same...

*/

#define BIG_VAL 1e20

#include "continuousML.h"
#include "ObjFunc.h"
#include "Maximize.h"
#include "structures.h"
#include "nexus.h"
#include "NRCvectorUtils.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>


double invertmatrix(double **a,int N, double **inv);
void ludcmpdouble(double **r, int n, int *indx, double *f);
double **dmatrix(long nrl, long nrh, long ncl, long nch);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);
int *ivector(long nl, long nh);
void free_ivector(int *v, long nl, long nh);

double **gVCV, **gV, **gInv, *gX, *gXm, * gW; // globals
int gNParm, gNT; 
typedef double ** DoubleMatrix;
#define MAX_MODELS 32
DoubleMatrix * gDMVCV;

/**********************************************************************************

Returns the negative log of the likelihood function of the MV normal model 

L = (2 Pi)^(-0.5k) {det V}^(-0.5) exp ( - 0.5 * (X-m)'vInv (X-m) )

log L = constant - 0.5 logDet(V) - 0.5 * {(X-m)'vInv (X-m)} 

	where the constant can be ignored for our purposes (estimation and LR tests)

Have to use a bunch of global variables and arrays, matrices to do this calculation since the prototype is fixed above.
All such structures follow NRC protocols as below. 

nParm = number of free parameters to be estimated (minimum of two)
NT    = number of terminal taxa

p[] is the vector of current parameter values: p[1]...p[nParm-1] contain rates, p[nParm] contains the mean value of the trait

gVCV[1...NT][1...NT] is the variance covariance matrix

gVCM[1..NT][1...NT] is a matrix of integers which indicates which rate parameter is associated with which cell in the VCV matrix.
	These integers range over [0...Nparm-2]. They are determined by the user setting up the model. Clunky, but otherwise I have
	to traverse the tree and recalculate distances all the time. VCM will be set up once, as is VCV. The matrix to be inverted
	is then the direct product of VCV and VCM...see below.
gV[1..NT][1..NT] should be allocated globally in advance; this is the final variance covariance matrix
gInv[1..NT][1..NT] should be allocated globally in advance; this is the inverse of the final variance covariance matrix
gX[1..NT] should be allocated globally in advance; this is the X vector of observations
gXm[1..NT] should be allocated globally in advance; this is the (X - mean) vector
gW[1..NT] should be allocated in advance, a working temp vector
*/
double contObj(double *p)  // the standard objective function prototype for use throughout r8s
{
double logDet,logL,mean,M=0.0;
int nParm,NT,numModels;
int r,c,i,modelIx;
NT=gNT;
nParm=gNParm;
numModels=nParm-1;

for (i=1;i<=nParm-1;i++)
	if (p[i] <0.0) return (BIG_VAL); // if any of the rates go negative, return a crazy neg log likelihood

mean = p[gNParm]; 
for (r=1;r<=NT;r++)
	{
	for (c=1;c<=NT;c++) 
		{
		if (r > c) // since the matrix is symmetric, save this step for lower triangle
			gV[r][c]=gV[c][r];
		else	   // or, actually do the calculation...
			{
			gV[r][c]=0.0;
			for (modelIx=0;modelIx<numModels;modelIx++)
				gV[r][c] += gDMVCV[modelIx][r][c]*p[modelIx+1]; // modelIx+1 is the correct index into the p[]
			}
		}
	gXm[r]=gX[r]-mean;
	}
logDet=invertmatrix(gV,NT,gInv);

// *** do the (X-m)'gInv (X-m) multiplications

for (c=1;c<=NT;c++)
	{
	gW[c]=0.0;
	for (r=1;r<=NT;r++) 
		gW[c] += gInv[r][c] * gXm[r];
	}
for (c=1;c<=NT;c++)
	M += gXm[c]*gW[c];

// ***

logL = -0.5 * logDet - 0.5* M;

return -logL;
}

void contOptimize(TREE t,int nParm,int *numIter, double ftol,double linMinDelta,int *success )
{
extern struct NexDataType *gNexDataPtr;
StrListPtr DM, TL;
PtrList nodeList;
int NT,i,j, model,numModels,modelIx,ixTL;
double *p, obj,tip_value;
char * tip_name,*tip_value_str, *dummy,*found_tip;
NODE a;
DoubleMatrix dmVCV[MAX_MODELS]; 

DM=gNexDataPtr->DMList;
TL=gNexDataPtr->TaxaList;

nodeList=pNewList();
TreeToTaxaPtrList(t->root,nodeList);


gDMVCV = dmVCV;
NT=t->numTaxa;
gNParm=nParm;
gNT=NT;
gV = dmatrix(1,NT,1,NT);
gInv = dmatrix(1,NT,1,NT);
gX = dvector(1,NT); 
gXm = dvector(1,NT); 
gW = dvector(1,NT); 
p = dvector(1,nParm);
for (i=1;i<=NT;i++)
	{
	a=(NODE)(pListgetkthNode(nodeList, i)->item);
	tip_name=a->taxon_name;
	ixTL=findMatchStr(TL,tip_name); // lookup the taxon name from the matrix ordering and get the relevant data matrix ordering corresponding name
	tip_value_str=getkthStr(DM,(long)(ixTL));
	found_tip=getkthStr(TL,(long)(ixTL));
	tip_value=strtod(tip_value_str,NULL);
// printf("%i\t%s\t%s\t%s\t%f\n",i,tip_name,found_tip,tip_value_str,tip_value);
	gX[i]=tip_value;
	}
for (i=1;i<=nParm;i++)
	p[i]=1.0; //bad first guess!

numModels = nParm - 1; // this is the number of rates
if (numModels > MAX_MODELS)
	fatal("Number of rate categories exceeds limits in continuousML.c\n");
for (modelIx=0;modelIx<numModels;modelIx++)
	dmVCV[modelIx] = tree2VCV(t, modelIx);  // allocate the global matrix and set it up with values based on tree

obj=MinND(t,0, POWELL,contObj,NULL,p, nParm,numIter, ftol,linMinDelta,success );

printf("\n\nParameter estimates:\n");
for (i=1;i<=nParm-1;i++)
	printf("Model %2i rate = \t%f\n",i-1,p[i]);
printf("Mean trait = \t%f\n",p[nParm]);

printf ("Log Likelihood = %f\n",-obj);


}




/*
 *  invertmatrix.c
 *  Brian O'Meara 26viii04
 */

#define NR_END 1  


/*  Function to invert a matrix
 *  Input is a two dimensional double array of size N x N in standard C format with indices 0 to N-1, as well as an integer (N) describing the number of rows. Void output, as it changes the pre-existing matrix. 
 *
 *  Usage: "invertmatrix(&inmatrix[0][0],N);"
 *
 *  Note: To get N, use "int N=sizeof(inmatrix)/sizeof(inmatrix[0]);"
 *  Note: "invertmatrix" overwrites the matrix you pass with with the matrix's inverse, 
 *         so make sure you have stored a copy of the original matrix if you want to keep it.
 *
 * Dependencies: The function uses numerical recipes in c code. Since we want double precision, not float precision, we had to convert the convert_matrix.c function to dconvert_matrix, ludcmp.c to ludcmpdouble, and lubksb.c to lubksbdouble. dconvert_matrix and lubksdouble are included in this file, but ludcmpdouble.c must be kept separate, as it gives an error otherwise ["warning: passing arg 4 of `ludcmpdouble' from incompatible pointer type"]. Also include nrutil.c to compile.
  */

void lubksbdouble(double **q, int n, int *indx, double b[])
{
	int i,ii=0,ip,j;
	double sum;

	for (i=1;i<=n;i++) {
		ip=indx[i];
		sum=b[ip];
		b[ip]=b[i];
		if (ii)
			for (j=ii;j<=i-1;j++) sum -= q[i][j]*b[j];
		else if (sum) ii=i;
		b[i]=sum;
	}
	for (i=n;i>=1;i--) {
		sum=b[i];
		for (j=i+1;j<=n;j++) sum -= q[i][j]*b[j];
		b[i]=sum/q[i][i];
	}
}
 



double invertmatrix(double **a,int N, double **y) 
{
    double *col,*d,logDet=0.0;
    int r,c,step;
    int i,j,*indx;
    d=dvector(1,N);
    col=dvector(1,N);
    indx=ivector(1,N);
    ludcmpdouble(a,N,indx,d);
    for(j=1;j<=N;j++) {
        for(i=1;i<=N;i++) col[i]=0.0;
        col[j]=1.0;
        lubksbdouble(a,N,indx,col);
        for(i=1;i<=N;i++) y[i][j]=col[i];
    }
    for (j=1;j<=N;j++)
	{
	logDet += log (fabs(a[j][j])); // we could keep track of the sign in 'd', but taking abs val has same effect...
	}
    return logDet; 
    
}

#define NRANSI
#define TINY 1.0e-20;

void ludcmpdouble(double **r, int n, int *indx, double *f)
{
	int i,imax,j,k;
	double big,dum,sum,temp;
	double *vv;

	vv=dvector(1,n);
	*f=1.0;
	for (i=1;i<=n;i++) {
		big=0.0;
		for (j=1;j<=n;j++)
			if ((temp=fabs(r[i][j])) > big) big=temp;
		if (big == 0.0) nrerror("Singular matrix in routine ludcmpdouble");
		vv[i]=1.0/big;
	}
	for (j=1;j<=n;j++) {
		for (i=1;i<j;i++) {
			sum=r[i][j];
			for (k=1;k<i;k++) sum -= r[i][k]*r[k][j];
			r[i][j]=sum;
		}
		big=0.0;
		for (i=j;i<=n;i++) {
			sum=r[i][j];
			for (k=1;k<j;k++)
				sum -= r[i][k]*r[k][j];
			r[i][j]=sum;
			if ( (dum=vv[i]*fabs(sum)) >= big) {
				big=dum;
				imax=i;
			}
		}
		if (j != imax) {
			for (k=1;k<=n;k++) {
				dum=r[imax][k];
				r[imax][k]=r[j][k];
				r[j][k]=dum;
			}
			*f = -(*f);
			vv[imax]=vv[j];
		}
		indx[j]=imax;
		if (r[j][j] == 0.0) r[j][j]=TINY;
		if (j != n) {
			dum=1.0/(r[j][j]);
			for (i=j+1;i<=n;i++) r[i][j] *= dum;
		}
	}
	free_dvector(vv,1,n);
}
#undef TINY
#undef NRANSI
#define FREE_ARG void*

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}
double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}
void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}
void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}
#undef NR_END

               r8s/continuousML.h                                                                                  0000644 0000766 0000120 00000000223 10121427437 014051  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "TreeUtils.h"
void contOptimize(TREE t,int nParm,int *numIter, double ftol,double linMinDelta,int *success );
double contObj(double p[]);
                                                                                                                                                                                                                                                                                                                                                                             r8s/covarion.c                                                                                      0000644 0000766 0000120 00000032700 11645740567 013250  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /*

Module to implement a simple covarion model with binary switch,S, and binary trait, T.
Modified massively from ancestral.c

We transform the 2 2-state traits to a single four state trait:

	(S,T) => Z

	(0,0) => 0
	(0,1) => 1
	(1,0) => 2
	(1,1) => 3
	


Contains several routines from NRC and some minor modifications to same...

*/
#include "ObjFunc.h"
#include "Maximize.h"
#include "structures.h"
#include "nexus.h"
#include "covarion.h"
#include "TreeUtils.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SQR(x)        ((x)*(x))
#define MAX(a,b) ((a) >(b) ? (a):(b))
#define LARGE_VAL 1e100

#define NSTATES 3

double rootPrior[NSTATES]={1/3.,1/3.,1/3.};

double gR, gS;  //switch rate = s; trait rate = r

double ** covProb;

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
static double ** transitionProb(double s, double r, double t);
static double ** mat_transpose(double **A);
static double ** mat_mult(double **A,double **B);
static void mat_print(double **m, int maxi, int maxj);
double **dmatrix(long nrl, long nrh, long ncl, long nch);	// defined in continuousML.c
static void uppassCovarion(NODETYPE *node);
static void setupCLmaxopt(NODETYPE *node);
static void setupCL(NODETYPE *node); // conditional likelihoods



double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);

int gNParm, gNT, p2tindex; 

void covarionOptimize(TREE t,int *numIter, double ftol,double linMinDelta,int *success )
{
extern struct NexDataType *gNexDataPtr;
extern NODETYPE * gRoot;
StrListPtr DM, TL;
PtrList nodeList;
int nParm,NT,i,j, ixTL;
double *p, obj, **PT;
char * tip_name,*tip_value_str, *dummy,*found_tip;
NODE a;
float meanTip=0.0;
int tip_value;
double minobj,minS,minR;

double **m1, **m2;



DM=gNexDataPtr->DMList;
TL=gNexDataPtr->TaxaList;
gRoot = t->root;
gRoot -> length = 0.0;
printf ("Warning! Setting the root's subtending branch length to ZERO whether it is previously set or not\n");
if (gNexDataPtr->RateBlockParms.cov_brlens==1)
	printf("All branch lengths are set to = 1.0\n");
if (gNexDataPtr->RateBlockParms.cov_brlens==0)
	printf("Using user supplied branch lengths\n");
nodeList=pNewList();
TreeToTaxaPtrList(t->root,nodeList);




NT=t->numTaxa;
for (i=1;i<=NT;i++)
	{
	a=(NODE)(pListgetkthNode(nodeList, i)->item);
	tip_name=a->taxon_name;
	ixTL=findMatchStr(TL,tip_name); // lookup the taxon name from the matrix ordering and get the relevant data matrix ordering corresponding name
	tip_value_str=getkthStr(DM,(long)(ixTL));
	found_tip=getkthStr(TL,(long)(ixTL));
	tip_value=strtod(tip_value_str,NULL);
// printf("%i\t%s\t%s\t%s\t%f\n",i,tip_name,found_tip,tip_value_str,tip_value);

	// set up conditional likelihoods of leaves based on data
	switch (*tip_value_str)  // this is a single char
		{
		case '0': 
			(a->CL)[0] = 1.0;
			(a->CL)[1] = 1.0;
			(a->CL)[2] = 0.0;

			break;
		case '1': 
			(a->CL)[0] = 0.0;
			(a->CL)[1] = 0.0;  
			(a->CL)[2] = 1.0;
			break;
		case '?': 
			(a->CL)[0] = 1.0;
			(a->CL)[1] = 1.0;  
			(a->CL)[2] = 1.0;
			break;


		}


	}
nParm=2;
p = dvector(1,nParm);




if (gNexDataPtr->RateBlockParms.estCov==1)
	{
	gR=0.2; gS=0.3;
	
	nParm = 2;
	p[1] = gS; 
	p[2] = gR; 
	
	obj=MinND(t,0, POWELL,objCovarion,NULL,p, nParm,numIter, ftol,linMinDelta,success );
	

	}
else
	{
	p[1] = gNexDataPtr->RateBlockParms.s_rate;
	p[2] = gNexDataPtr->RateBlockParms.r_rate;
	obj=objCovarion(p); // just get this likelihood for these parameters
	}
gS = p[1]; gR = p[2];

if (gS > 0 && gR > 0)
	{
	setupCLmaxopt(t->root);
	uppassCovarion(t->root);
	
	if (gNexDataPtr->RateBlockParms.estCov==1)
		printf("Estimated parameters: %f %f; Maximum likelihood=%g\n",p[1],p[2],-obj);
	else
		printf("User supplied rate parameters: %f %f; Likelihood at optimal solution=%g\n",p[1],p[2],-obj);
	}
	
else
	printf ("Cov rate parameters were <= 0\n");

//printf("Currently NOT calculating marginals\n");
findMarginals(t->root);

t->root=gRoot; //kludge : see marginals code below

free_dvector(p,1,nParm);


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
	
#if 0
	if (isTip(a))  // don't reroot but do this hack for the moment: FIX THIS LATER...
		{
		if ((a->CL)[2] == 1.0) 
			{(a->CLmarg)[0]=0; (a->CLmarg)[1]=0; (a->CLmarg)[2]=1;}
		else
			{(a->CLmarg)[0]=0.5; (a->CLmarg)[1]=0.5; (a->CLmarg)[2]=0;}
		continue;
		}
#endif
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
 if (p[2] <0 ) {return LARGE_VAL;};
 
  gS = p[1]; gR = p[2];
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
//printCovarion(gRoot);
  return -obj; // it's a minimization
}
/**********************/

static void setupCLmaxopt(NODETYPE *node) // Pupko's 2002 algorithm
{
  NODETYPE *child;
  double lsum,cl,max, brLen;
  int i,j,k,opt;
  double **PT;
  if (!node) return;
  child=node->firstdesc;
  SIBLOOP(child)
    {
      setupCLmaxopt(child);
    }
  PT = transitionProb(gS,gR,node->length);
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
		if (cl > max)
			{
			max=cl;
			opt=j;
			}
		}
	(node->CLmax)[i]=max;
	(node->CLopt)[i]=opt;
	}
  free_dmatrix(PT,0,NSTATES,0,NSTATES);
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
void printCovarion(NODETYPE *node)
{
    float diff;
    double sum;
    NODETYPE *child;

      if(*(node->taxon_name))
	    printf("%.12s\t", node->taxon_name);
      else
	    printf("    Node %li  \t", node->id);

	sum = (node->CLmarg)[0] + (node->CLmarg)[1] + (node->CLmarg)[2];
	// Order of output: CL, CLmax, CLmarg, CLnorm, CLopt
	printf(": %4.2e %4.2e %4.2e",(node->CL)[0],(node->CL)[1],(node->CL)[2]);
	printf(": %4.2e %4.2e %4.2e",(node->CLmax)[0],(node->CLmax)[1],(node->CLmax)[2]);
	printf(": %4.2e %4.2e %4.2e",(node->CLmarg)[0],(node->CLmarg)[1],(node->CLmarg)[2]);
	if (sum > 0.0)
		printf(": %4.2f %4.2f %4.2f",(node->CLmarg)[0]/sum,(node->CLmarg)[1]/sum,(node->CLmarg)[2]/sum);
	else
		printf(":               ");
	printf(": %2i %2i %2i => %i\n",(node->CLopt)[0],(node->CLopt)[1],(node->CLopt)[2],node->opt);


    if(!isTip(node))
	{
    	child=node->firstdesc;
    	SIBLOOP(child) 
		printCovarion(child);
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
V=dmatrix(0,NSTATES,0,NSTATES);
VT=dmatrix(0,NSTATES,0,NSTATES);
P=dmatrix(0,NSTATES,0,NSTATES);
D=dmatrix(0,NSTATES,0,NSTATES);

// Hijack the branch lengths if called for by the r8s block.
if (gNexDataPtr->RateBlockParms.cov_brlens==1) t=1.0;


cc = 1/sqrt(3);

if (r != s)
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
return P;
}

/**********************/

static void setupCL(NODETYPE *node) // conditional likelihoods
{
  NODETYPE *child;
  double lsum,cl, **PT;
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
 				PT = transitionProb(gS,gR,child->length);
  				//mat_print(PT,NSTATES,NSTATES);
            	lsum = 0.0;
            	for (i=0;i<NSTATES;i++)  // child's state
                        {
                        lsum += PT[i][k]*((child->CL)[i]);
                        }
            	cl *= lsum;
  				free_dmatrix(PT,0,NSTATES,0,NSTATES);
            	}
        	(node->CL)[k] = cl;
        	}
  return ;
}
/**********************/

                                                                r8s/covarion.h                                                                                      0000644 0000766 0000120 00000000346 11637722625 013252  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  void findMarginals(NODETYPE *root);
void printChanges(NODETYPE *node);
void printCovarion(NODETYPE *node);
void covarionOptimize(TREE t,int *numIter, double ftol,double linMinDelta,int *success );
double objCovarion(double p[]);

                                                                                                                                                                                                                                                                                          r8s/distance.c                                                                                      0000644 0000766 0000120 00000016223 10357547124 013215  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "distance.h"
#include "nexus.h"
#include "string.h"
#include "myOutput.h"
#include "MyUtilities.h"
#include "memory.h"
#include <math.h>


/* Calculates hamming and K2P distances between two sequences;
also returns the absolute number of transitions and transversions in nP and nQ respectively*/


/*
Hamming K2P distances validated by comparison of two sample data sets to PAUP* d47.
*/

double distance(char * pi, char * pj, char * pRow1, int kind, long * nP, long * nQ)
{
double a,b,Ahat,Bhat,P,Q;
const float A = 20/19;
long seqLength;

seqLength = PQCalc1(pi,pj,pRow1,&P,&Q, nP,nQ);

switch (kind)
	{
	case 0: /* Hamming distance */
		return P+Q;
		
	case 1:	/* K2P distance */
		a=1.0/(1-2*P-Q);
		b=1.0/(1-2*Q);
		Ahat=0.5*log(a)-0.25*log(b);
		Bhat=0.5*log(b);
		return Ahat+Bhat;

	}

}




long PQCalc1(char * pi, char * pj, char * pRow1, double *P, double *Q, long *nP, long * nQ)

/*  Returns the number of valid sites (nongap, nonmissing) in sequences
between taxa i and j, where i and j are NEXUS taxon numbers.  Stores P and Q,
which are the proportion of transitions and transversions, respectively.
On input takes three pointers: two to the relevant char strings and the third to the
first row of the matrix (in the case where '.' is used) or NULL if '.' is not used.
 */

{

	extern struct NexDataType *gNexDataPtr;
	int* excArray;		/* array for exclusion set */
	long isite, validcount=0,transitcount=0,transvcount=0, slength;
	long	missingcount1=0,missingcount2=0,missingcountsame=0,
			gapcount1=0,gapcount2=0,gapcountsame=0;
	char ci,cj,gap,match,missing;
	excArray=gNexDataPtr->excArray;
	gap=gNexDataPtr->gapchar;
	missing=gNexDataPtr->missingchar;
	match=gNexDataPtr->matchchar;
	
	slength = strlen(pi);
	for (isite=0; isite < slength; isite++)
	  if (excArray[isite] > 0)  /* process site only if weight positive */
		{
		ci=pi[isite];
		cj=pj[isite];
		
		if (ci==missing) ++missingcount1;
		if (cj==missing) ++missingcount2;
		if ((ci==missing) && (cj==missing)) ++missingcountsame;
		if (ci==gap) ++gapcount1;
		if (cj==gap) ++gapcount2;
		if ((ci==gap) && (cj==gap)) ++gapcountsame;
		
		
		if (pRow1)  /* if a first row in matrix was passed to us (rather than NULL) */
			{ 
			if (ci==match) ci = pRow1[isite];  /* check for 'period' format in sequences */
			if (cj==match) cj = pRow1[isite];  /* if present, set to data for first row */
			}
		if (  strchr("ACGT",ci) &&  strchr("ACGT",cj)  ) /* only consider when both sites
														are in ACGT */
			{
			validcount+=excArray[isite]; /* weight site */
			if (ci != cj )	/* only consider variable sites in the following */
				if (
					((ci=='C') && (cj=='T'))  ||
					((ci=='T') && (cj=='C'))  ||
					((ci=='A') && (cj=='G'))  ||
					((ci=='G') && (cj=='A'))
					)
					transitcount+=excArray[isite];/* weight site */
				else
					transvcount+=excArray[isite];/* weight site */
			}
		}
	if (validcount > 0) 
		{
		*P=(double)transitcount/validcount;
		*Q=(double)transvcount/validcount;
		*nP=transitcount;
		*nQ=transvcount;
		}
	return validcount;

}

long aaCalc1(char * pi, char * pj, char * pRow1, double *P,long *n)

/*  Returns the number of valid sites (nongap, nonmissing) in sequences
between taxa i and j, where i and j are NEXUS taxon numbers.  Stores P and Q,
which are the proportion of transitions and transversions, respectively.
On input takes three pointers: two to the relevant char strings and the third to the
first row of the matrix (in the case where '.' is used) or NULL if '.' is not used.
 */

{

	extern struct NexDataType *gNexDataPtr;
	int* excArray;		/* array for exclusion set */
	long isite, validcount=0,transitcount=0,transvcount=0, slength;
	long	missingcount1=0,missingcount2=0,missingcountsame=0,
			gapcount1=0,gapcount2=0,gapcountsame=0;
	char ci,cj,gap,match,missing;
	*n=0;
	excArray=gNexDataPtr->excArray;
	gap=gNexDataPtr->gapchar;
	missing=gNexDataPtr->missingchar;
	match=gNexDataPtr->matchchar;
	
	slength = strlen(pi);
	for (isite=0; isite < slength; isite++)
	  if (excArray[isite] > 0)  /* process site only if weight positive */
		{
		ci=pi[isite];
		cj=pj[isite];
		
		if (ci==missing) ++missingcount1;
		if (cj==missing) ++missingcount2;
		if ((ci==missing) && (cj==missing)) ++missingcountsame;
		if (ci==gap) ++gapcount1;
		if (cj==gap) ++gapcount2;
		if ((ci==gap) && (cj==gap)) ++gapcountsame;
		
		
		if (pRow1)  /* if a first row in matrix was passed to us (rather than NULL) */
			{ 
			if (ci==match) ci = pRow1[isite];  /* check for 'period' format in sequences */
			if (cj==match) cj = pRow1[isite];  /* if present, set to data for first row */
			}
		if (!(ci==missing || cj==missing || ci==gap || cj==gap)) /* also check for  valid AA codes in here! */
			{
			validcount+=excArray[isite]; 
			if (ci != cj )	
				(*n)+=excArray[isite];
			}
		}
	if (validcount > 0) 
		*P=(double)(*n)/validcount;
	else 
		*P=0;
	return validcount;

}


void doDistance(StrListPtr aTaxaList)
{
/* do distance matrix */

int kk,kind,j,i,ix,jx,taxonID;
double *X, *Y;
int numDistances, idist=0;
char *taxon1,*taxon2, *taxon, *dummy,*pi,*pj,*pRow1;
long NList;
double d;
long nP,nQ;
extern struct NexDataType *gNexDataPtr;	

NList=lengthList(aTaxaList);
numDistances=NList*NList/2.-NList/2.; /* upper triangular nondiagonal entries */
X=(double *)myMalloc(numDistances*sizeof(double));
Y=(double *)myMalloc(numDistances*sizeof(double));
for (ix=1;ix<=NList;ix++) /* convert any taxon ids to taxon names 
				unless already stored that way*/
		{
		taxon=getkthStr(aTaxaList,ix);
		if(isStrInteger(taxon))
			{
			taxonID=strtod(taxon,&dummy);
			setkthNode(aTaxaList, ix, getkthStr(gNexDataPtr->TaxaList,taxonID));
			}
		}


PRINT_LINE;
printf("\nDistances for selected taxa\n\n");
PRINT_LINE;

for (kk=0;kk<=2;kk++)
	{
	switch (kk) 
		{
		case 0:
			kind=0;
			printf("\n\nAbsolute Distance Matrix:\n\n          ");
			break;
		case 1:
			kind=1;
			printf("\n\nKimura 2-parameter Distance Matrix:\n\n          ");
			break;
		case 2:
			kind=1;
			printf("\n\nTransition/Transversion Matrix (transitions above diagonal):\n\n          ");
			break;
		}
	
	for (j=1;j<=NList;j++)
		{
		printf("%8.8s  ",getkthStr(aTaxaList,(long)j ));
		}
	printf("\n");
	for (i=1;i<=NList;i++)
		{
		taxon1 = getkthStr(aTaxaList,(long)i );
		printf("%8.8s  ",taxon1);
		for (j=1;j<=NList;j++)
				{ /* set up the unsorted 3-list */
				taxon2 = getkthStr(aTaxaList,(long)j );
				ix = findMatchStr(gNexDataPtr->TaxaList, taxon1);
				if (ix ==0)
					doGenericAlert ("Matching taxon label not found in WuLi");
				jx = findMatchStr(gNexDataPtr->TaxaList, taxon2);
				if (jx ==0)
					doGenericAlert ("Matching taxon label not found in WuLi");
				pRow1=getkthStr(gNexDataPtr->DMList,1);	
				pi=getkthStr(gNexDataPtr->DMList,ix);
				pj=getkthStr(gNexDataPtr->DMList,jx);
				d=distance(pi, pj, pRow1, kind, &nP, &nQ);
				if (i<j)
					{
					if (kk==2)
						{
						printf("%8li  ",nP); /* transitions */
						Y[idist]=nP;
						X[idist]=nQ;
						idist++;
						}
					else
						printf("%8f  ",d);
					}
				else
					if (i>j)
						{
						if (kk==2)
							printf("%8li  ",nQ); /* transversion*/
						else
							printf("%8li  ",nP+nQ);
						}
					else
						printf("      --  ");
				}
		printf("\n");
		}
	}
PRINT_LINE;

dumbPlot(X, Y, numDistances);

return;
}
                                                                                                                                                                                                                                                                                                                                                                             r8s/distance.h                                                                                      0000644 0000766 0000120 00000000547 07565204055 013224  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "structures.h"
void doDistance(StrListPtr aTaxaList);
long aaCalc1(char * pi, char * pj, char * pRow1, double *P,long *n);
long PQCalc(int i, int j, double *P, double *Q);
long PQCalc1(char * pi, char * pj, char * pRow1, double *P, double *Q, long *nP, long * nQ);
double distance(char * pi, char * pj, char * pRow1, int kind, long * nP, long * nQ);
                                                                                                                                                         r8s/main.c                                                                                          0000644 0000766 0000120 00000004155 11153632653 012345  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /********************************************************

		r8s 
*/

#define VERSION 1.72

/*

********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "storeNexusFile.h"
#include "nexus.h"



int main(int argc,char * argv[])
    {
    char *p;
    extern int gInteractive;
    extern struct NexDataType *gNexDataPtr;	/* This is THE data structure for the NEXUS data */
    char *theNexusFileBuffer, fnInput[FILENAME_MAX],theArg,CLBuf[256];
    FILE * inStream =NULL;
    int cFlag=0,c=0;
    long l;
    
    gInteractive=1; /* default is interactive mode */
    gNexDataPtr=initialize_nexus();

    fprintf(stderr,"r8s version %4.2f (compiled %s)\n",VERSION,__DATE__);
    if (argc == 1)
	{
	    doInteractive();
	    return 1;
	}
    
    else
    for (++argv, c=1;c<argc;c++)
	{
	if (**argv=='-')
		{
		p=*argv;
		++p;
		switch(tolower(*p))
				{
				case 'b':
					gInteractive=0;break; /* set to batch mode */
				case 'c':
					++argv;
					strcpy(CLBuf, *argv);
					gInteractive=0;
					cFlag=1;
					break;
			        case 'f':
					++argv;
					strcpy(fnInput, *argv);   /* set file name */
					if (  !(inStream=fopen(fnInput,"r")) )
						{
						printf("Error opening %s\n", fnInput);
						exit(1);
						}
					else
						fprintf(stderr, "[...reading file %s]\n", fnInput);
					break;
			        case 'v':
					printf("r8s version %4.2f (%s)\n",VERSION,__DATE__);
					break;
			        case 'h':
					printf("Usage: r8s [-b] [-h] [-v] [-f datafile] [-c commandstring]\n");
					printf("\t-b\tBatch process the datafile\n");
					printf("\t-h\tThis information...\n");
					printf("\t-v\tPrint version and compilation date\n");
					printf("\t-c\tOpen and execute commandstring immediately\n");
				}
		}
	++argv; if (!*argv) break;
	}
    
    if (!gNexDataPtr)
	fatal("Failure to allocate nexus data structure in main.c");
    if (inStream)
	    {
	    theNexusFileBuffer=storeNexusFile(inStream);
	    readNexusFile(theNexusFileBuffer);
	    };
    
    if(gInteractive)
	doInteractive();
    if(cFlag)
	doCommandLineControl(CLBuf);    
    
    
    return 1;
    }
                                                                                                                                                                                                                                                                                                                                                                                                                   r8s/memory.c                                                                                        0000644 0000766 0000120 00000003564 07565204055 012737  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include <stdlib.h>
#include <stdio.h>
#include "memory.h"
#include "errno.h"
/* #include <malloc.h> */

/* Modified July 1996 back to standard C without Mac toolbox calls */



void * myMalloc(size_t theSize)
{
void * p;

errno = 0;
p = (char *)malloc(theSize);
		/* print_mem_dbg(); */
if (errno)
	{
	perror("Low Level allocation error in myMalloc");
	exit(1);
	}
return p;
}

void myFree(void * p)
{
errno=0;

free(p);
if (errno)
	{
	perror("Low Level free error in myFree");
	exit(1);
	}
return;

}

void * myRealloc(void * p, size_t theSize)
{
long i,j;
void * pp;

errno=0;
pp = (char *)realloc(p,theSize);
if (errno)
	{
	perror("Low Level reallocation error in myRealloc");
	exit(1);
	}


return pp;
}

#if MEM_DBG
/* if you want to do some serious memory debugging, set this to 1 in header 
but you may have to link to library with -lmalloc in the Makefile.
PROBABLY SGI specific */


void print_mem_dbg(char *file_name,int line)
{
struct mallinfo mi;
long netSpace;
mi=mallinfo();
printf("%s:%d\n[%i  %i  %i  %i  %i  %i  %i]\n", file_name,line, 
	mi.uordblks,mi.usmblks, 
	mi.arena, mi.ordblks, mi.smblks, mi.fsmblks,mi.fordblks);

return;
}



#endif
#if 0 
     struct mallinfo  {
             int arena;         /* total space in arena */
             int ordblks;       /* number of ordinary blocks */
             int smblks;        /* number of small blocks */
             int hblkhd;        /* space in holding block headers */
             int hblks;         /* number of holding blocks */
             int usmblks;       /* space in small blocks in use */
             int fsmblks;       /* space in free small blocks */
             int uordblks;      /* space in ordinary blocks in use */

             int fordblks;      /* space in free ordinary blocks */
             int keepcost;      /* space penalty if keep option */
                                /* is used */
     }
 
#endif
                                                                                                                                            r8s/memory.h                                                                                        0000644 0000766 0000120 00000000324 07565204055 012733  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include <stdlib.h>
#define MEM_DBG 0	/* for memory debugging */

void * myMalloc(size_t theSize);
void myFree(void * p);
void * myRealloc(void * p, size_t theSize);
void print_mem_dbg(char *file_name,int line);
                                                                                                                                                                                                                                                                                                            r8s/moment.c                                                                                        0000644 0000766 0000120 00000001313 07565204055 012714  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include <math.h>
#include "MyUtilities.h"
#include "moment.h"

void moment(double data[], int n, double *ave, double *adev, double *sdev,
	double *var, double *skew, double *curt)
{
	int j;
	float ep=0.0,s,p;

	if (n <= 1) doGenericAlert("n must be at least 2 in moment");
	s=0.0;
	for (j=1;j<=n;j++) s += data[j];
	*ave=s/n;
	*adev=(*var)=(*skew)=(*curt)=0.0;
	for (j=1;j<=n;j++) {
		*adev += fabs(s=data[j]-(*ave));
		*var += (p=s*s);
		*skew += (p *= s);
		*curt += (p *= s);
	}
	*adev /= n;
	*var=(*var-ep*ep/n)/(n-1);
	*sdev=sqrt(*var);
	if (*var) {
		*skew /= (n*(*var)*(*sdev));
		*curt=(*curt)/(n*(*var)*(*var))-3.0;
	} /* 
	  else doGenericAlert("No skew/kurtosis when variance = 0 (in moment)");
	  */
}
                                                                                                                                                                                                                                                                                                                     r8s/moment.h                                                                                        0000644 0000766 0000120 00000000166 07565204055 012726  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  void moment(double data[], int n, double *ave, double *adev, double *sdev,
	double *var, double *skew, double *curt);
                                                                                                                                                                                                                                                                                                                                                                                                          r8s/myOutput.h                                                                                      0000644 0000766 0000120 00000000111 07565204055 013263  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #define PRINT_LINE printf("-----------------------------------------\n")
                                                                                                                                                                                                                                                                                                                                                                                                                                                       r8s/nextToken2.c                                                                                    0000644 0000766 0000120 00000013573 07565204055 013471  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "nexus.h"
#include "MyUtilities.h"
#include "memory.h"
#include "structures.h"



#define isNL_NEXUSwhiteSpace(c)  ( strchr(" \t\v\f", (c)) || (((c) <= 6) && ((c) >= 0)))
#define isNL_NEXUSpunct(c) ( strchr(NL_punct,(c)) )
#define NULL_RETURN {*aTokenPtr='\0';return aTokenPtr;}

#define CHECK_OVERFLOW  if (cix>=MAX_TOKEN_SIZE-1) doGenericAlert("Token Size Exceeded in nextToken")
/********************************************************/


/*	
	Gets the next token from input stream 'fpointer', and copies it onto the global
	buffer pointed to by 'aTokenPtr'.  If there is NO next token, we copy a null
	string onto that buffer.  That's a signal for the main caller routine...

	If the global variable gNewLine=1 then the newline characters, '\n' and '\r'
	ARE returned as individual tokens,  when encountered.  The normal state is
	gNewLine=0,  which treats these as white space delimiters too.  The only time
	NEXUS file needs to think about newlines is when reading interleaved matrices!

*/

char *nextToken(void)


	{
	extern char *aTokenPtr; 

	extern char * bufPtr;	/*declared and initialized in readNexusFile.c */
	extern int gNewLine;	/*declared and set in readNexusFile.c */

	char *punct="()[]{}/\\,;:=*\'\"`+";	/* these are NEXUS definitions */
	char *NL_punct="()[]{}/\\,;:=*\'\"`+\r\n";	/* NEXUS definitions plus stuff for newlines*/
	char c;
	int cix=0;	/* counter to monitor token size */
	
	*aTokenPtr='\0';
	
	if  ((c=*bufPtr++) == '\0') NULL_RETURN
	
	/* First block below handles the case where newline characters must be reported*/
	
	if (gNewLine)
	    {
	    while (( isNL_NEXUSwhiteSpace(c) ) || (c=='['))  
			    /* this whole loop is in case multiple comments separated by whitespace */ 
		    {
		    while ( isNL_NEXUSwhiteSpace(c) )  /* skip white space and land on next 'c'*/
				    {
				    c=*bufPtr++;
				    if (c=='\0')  NULL_RETURN;/* check for embedded EOF */
				    }
			    
		    if (c=='[')		/* skip the comment and land on next 'c' after comment */
			    {
			    while (c !=']')
				    {
				    c=*bufPtr++;
				    if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
				    }
			    c=*bufPtr++;	/* get next char after ']' */
			    if (c=='\0') NULL_RETURN;
			    }
		    }
	    
    
	    if (c=='\'')		/* deal with single-quoted tokens */
		    {
    
		    aTokenPtr[cix++]=c;
		    CHECK_OVERFLOW;
		    while (  (c=*bufPtr++) != '\'')
			    {
			    if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
			    aTokenPtr[cix++]=c;	/* this is a valid character in the word, add to token */
			    CHECK_OVERFLOW;
			    }
		    aTokenPtr[cix++]=c;	/* add the terminating quote too */
		    CHECK_OVERFLOW;
		    aTokenPtr[cix]='\0'; /* null terminate the string */
#if STU
		    strtoupper(aTokenPtr);
#endif
		    return(aTokenPtr);
		    }	/* return everything between single quotes, including the quotes, as a token*/		
	    aTokenPtr[cix++]=c;		/* char is either punctuation or part of word, so add it to token */
	    CHECK_OVERFLOW;
	    
	    if (!isNL_NEXUSpunct(c))	/* next char is part of word, so add all of word until white,punct,eof,
									    or Token size exceeded */
		    {
		    for (;;)
			    {
			    c=*bufPtr++;  
			    if (  isNL_NEXUSpunct(c) || isNL_NEXUSwhiteSpace(c) 
						    ||  (c == '\0') )
				    {
				    --bufPtr; /* word is terminated by some c that is not part of word;
								     push c back into stream and deal with it on
								    next call to this function; meantime, break out, 
								    and return this token*/
				    break;
				    };
			    aTokenPtr[cix++]=c;	/* this is a valid character in the word, add to token */
			    CHECK_OVERFLOW;
			    }
		    }
	    }
	else  /* identical to block above except for character test definitions! */
	    {
	    while (( isNEXUSwhiteSpace(c) ) || (c=='['))  
			    /* this whole loop is in case multiple comments separated by whitespace */ 
		    {
		    while ( isNEXUSwhiteSpace(c) )  /* skip white space and land on next 'c'*/
				    {
				    c=*bufPtr++;
				    if (c=='\0')  NULL_RETURN;/* check for embedded EOF */
				    }
			    
		    if (c=='[')		/* skip the comment and land on next 'c' after comment */
			    {
			    while (c !=']')
				    {
				    c=*bufPtr++;
				    if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
				    }
			    c=*bufPtr++;	/* get next char after ']' */
			    if (c=='\0') NULL_RETURN;
			    }
		    }
	    
    
	    if (c=='\'')		/* deal with single-quoted tokens */
		    {
		    aTokenPtr[cix++]=c;
		    CHECK_OVERFLOW;
		    while (  (c=*bufPtr++) != '\'')
			    {
			    if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
			    aTokenPtr[cix++]=c;	/* this is a valid character in the word, add to token */
			    CHECK_OVERFLOW;
			    }
		    aTokenPtr[cix++]=c;	/* add the terminating quote too */
		    CHECK_OVERFLOW;
		    aTokenPtr[cix]='\0'; /* null terminate the string */
#if STU
		    strtoupper(aTokenPtr);
#endif
		    return(aTokenPtr);
		    }	/* return everything between single quotes, including the quotes, as a token*/		
	    aTokenPtr[cix++]=c;		/* char is either punctuation or part of word, so add it to token */
	    CHECK_OVERFLOW;
	    
	    if (!isNEXUSpunct(c))	/* next char is part of word, so add all of word until white,punct,eof,
									    or Token size exceeded */
		    {
		    for (;;)
			    {
			    c=*bufPtr++; 
			    if (  isNEXUSpunct(c) || isNEXUSwhiteSpace(c) 
						    ||  (c == '\0') )
				    {
				    --bufPtr; /* word is terminated by some c that is not part of word;
								     push c back into stream and deal with it on
								    next call to this function; meantime, break out, 
								    and return this token*/
				    break;
				    };
			    aTokenPtr[cix++]=c;	/* this is a valid character in the word, add to token */
			    CHECK_OVERFLOW;
			    }
		    }
	    }




	
	
		
	aTokenPtr[cix]='\0'; /* null terminate the string */
#if STU
	strtoupper(aTokenPtr);
#endif
	return(aTokenPtr);
	}
                                                                                                                                     r8s/nexus.h                                                                                         0000644 0000766 0000120 00000010223 11637204150 012553  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #ifndef _NEXUS
#define _NEXUS

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdlib.h>
#include "TreeUtils.h"
#include "structures.h"
#include "MyUtilities.h"

#define MAX_TOKEN_SIZE 10000		/* we've got room */
#define MAX_LOCAL_TOKEN_SIZE 128
#define MAXTREES	100		/* maximum number of tree descriptions stored from file */

#if 0
#define isEqual(a,b)		(!strcmp((a),(b)))
#define isEqualUL(a,b)		(!strcmp((a),(b)))
#endif

/* (!strcmp((strtoupper(a), a),(strtoupper(b), b)))*/

#define isNEXUSpunct(c) 		( strchr(punct,(c)) )
#define isNEXUSwhiteSpace(c)	( isspace((c)) || (((c) <= 6) && ((c) >= 0)))
	/* current NEXUS format also excludes ASCII 0-6 */



struct	RBP 	{
				int	clampRoot; /* 0 = separate subtree
						optimizations */
				int	isBS;	/* toggle bootstrap */
				long 	NReps;	/* bootstrap replicates */
				long 	seed;	/* random number seed */
				int	RRtype; /* 0=WuLi; 1=Steel et al.*/
				double	npexp;	/* exponent in the NP optimization */
				int	verbose;/* verbosity for rate block */
				int	num_restarts;
				int	num_rate_guesses;
				int	num_time_guesses;
				double	local_factor; /* fractional tolerance for
					two points to be considered the same */
				double perturb_factor; /*fractional displacement to look
					for another optimum */
				double	smoothing;	/* smoothing factor in penalized like*/
				double  ftol;		/* fraction func tolerance */
				double	barrierTol;
				int	maxIter;
				int	maxBarrierIter;
				double	initBarrierFactor;
				double  barrierMultiplier;
				double	linminOffset;
				double	contractFactor;
				int	maxContractIter;
				int	showConvergence;
				int	checkGradient;
				int	showGradient;
				int	RatesAreGamma;	/* across sites */
				double  alpha;		/* shape param */
				double activeEpsilon;	/* fractional distance from a constraint; if closer than this distance, a solution is said to be "active" */
				long	numSites;
				int	clockFmt;	/* 1 = trees assumed to be ultrametric on input */
				int	lengthFmt;	/* 0 = branch lengths are in numbers of subst.; 1 = subst/site */
				int	roundFlag;	/* 0 = branch lengths are not rounded on input; 1 = rounded */
				int	PenaltyType;
				int	NeighborPenalty; // 1=penalize with neighbor variance; 0=old style ancestor/desc squared
				float minRateFactor; // a fraction of the average rate to impose a min on all rates under PL
				float minDurFactor; // a fraction of the root's age to impose a min duration for 0-length 
							// terminal branches 
				int estCov;	// should we try to estimate the covarion matrix rate parameters?
				double s_rate;
				double r_rate; // s and r rates in covarion matrix
				int	cov_brlens; // set to 1 if we will set all branch lengths to 1; otherwise use supplied values
				};	
				

/* This is the data structure containing all the information for a NEXUS file */

struct 	NexDataType {
			int		isChars;
			int		isTrees;
			int		isTaxa;
			int		isTranslate;	/*...flags for when these elements are read */
			int 		NTaxa;			/* number of taxa */
			int 		NChars;			/* number of characters */
			int 		Intlvflag;		/* flag is set if data matrix is interleaved */
			char		matchchar;
			char		gapchar;
			char		missingchar;
			int			NumTrees;		/* number of trees in data structure */
			StrListPtr 	TaxaList;		/* list of taxon names */
			StrListPtr 	TDList;			/* list of tree descriptions */
			StrListPtr 	TDLabelList;		/* list of tree description labels */
			StrListPtr 	DMList;			/* The data matrix as a list of row strings*/
			StrListPtr	TransList;		/* Translation table list */
			StrListPtr 	TaxSetNameList;		/* list of names of taxsets*/
			int*		excArray;		/* array of flags for excluding characters */
			struct RBP	RateBlockParms;	/* the rate block parameters */
			PtrList		inTrees;	/* list of trees */
			PtrList		TaxSetLists;	/* list of pointers to the taxsets, each of which
								is a StrList */
				};
				
/* All lists use my 'list' data structure */
void readNexusFile(char * buffer);

void TreeSummary(int whichTree);
char *nextToken(void);
void freeNexusStructure(struct NexDataType *nex);
struct NexDataType * initialize_nexus(void);
void doRateBlock(void);
void doCommandLineControl(char *buffer);
void doInteractive(void);
#endif
                                                                                                                                                                                                                                                                                                                                                                             r8s/nr.h                                                                                            0000755 0000766 0000120 00000075357 07565204055 012067  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #ifndef _NR_H_
#define _NR_H_

#ifndef _FCOMPLEX_DECLARE_T_
typedef struct FCOMPLEX {float r,i;} fcomplex;
#define _FCOMPLEX_DECLARE_T_
#endif /* _FCOMPLEX_DECLARE_T_ */

#ifndef _ARITHCODE_DECLARE_T_
typedef struct {
	unsigned long *ilob,*iupb,*ncumfq,jdif,nc,minint,nch,ncum,nrad;
} arithcode;
#define _ARITHCODE_DECLARE_T_
#endif /* _ARITHCODE_DECLARE_T_ */

#ifndef _HUFFCODE_DECLARE_T_
typedef struct {
	unsigned long *icod,*ncod,*left,*right,nch,nodemax;
} huffcode;
#define _HUFFCODE_DECLARE_T_
#endif /* _HUFFCODE_DECLARE_T_ */

#include <stdio.h>

#if defined(__STDC__) || defined(ANSI) || defined(NRANSI) /* ANSI */

void addint(double **uf, double **uc, double **res, int nf);
void airy(float x, float *ai, float *bi, float *aip, float *bip);
void amebsa(float **p, float y[], int ndim, float pb[],	float *yb,
	float ftol, float (*funk)(float []), int *iter, float temptr);
void amoeba(float **p, float y[], int ndim, float ftol,
	float (*funk)(float []), int *iter);
float amotry(float **p, float y[], float psum[], int ndim,
	float (*funk)(float []), int ihi, float fac);
float amotsa(float **p, float y[], float psum[], int ndim, float pb[],
	float *yb, float (*funk)(float []), int ihi, float *yhi, float fac);
void anneal(float x[], float y[], int iorder[], int ncity);
double anorm2(double **a, int n);
void arcmak(unsigned long nfreq[], unsigned long nchh, unsigned long nradd,
	arithcode *acode);
void arcode(unsigned long *ich, unsigned char **codep, unsigned long *lcode,
	unsigned long *lcd, int isign, arithcode *acode);
void arcsum(unsigned long iin[], unsigned long iout[], unsigned long ja,
	int nwk, unsigned long nrad, unsigned long nc);
void asolve(unsigned long n, double b[], double x[], int itrnsp);
void atimes(unsigned long n, double x[], double r[], int itrnsp);
void avevar(float data[], unsigned long n, float *ave, float *var);
void balanc(float **a, int n);
void banbks(float **a, unsigned long n, int m1, int m2, float **al,
	unsigned long indx[], float b[]);
void bandec(float **a, unsigned long n, int m1, int m2, float **al,
	unsigned long indx[], float *d);
void banmul(float **a, unsigned long n, int m1, int m2, float x[], float b[]);
void bcucof(float y[], float y1[], float y2[], float y12[], float d1,
	float d2, float **c);
void bcuint(float y[], float y1[], float y2[], float y12[],
	float x1l, float x1u, float x2l, float x2u, float x1,
	float x2, float *ansy, float *ansy1, float *ansy2);
void beschb(double x, double *gam1, double *gam2, double *gampl,
	double *gammi);
float bessi(int n, float x);
float bessi0(float x);
float bessi1(float x);
void bessik(float x, float xnu, float *ri, float *rk, float *rip,
	float *rkp);
float bessj(int n, float x);
float bessj0(float x);
float bessj1(float x);
void bessjy(float x, float xnu, float *rj, float *ry, float *rjp,
	float *ryp);
float bessk(int n, float x);
float bessk0(float x);
float bessk1(float x);
float bessy(int n, float x);
float bessy0(float x);
float bessy1(float x);
float beta(float z, float w);
float betacf(float a, float b, float x);
float betai(float a, float b, float x);
float bico(int n, int k);
void bksub(int ne, int nb, int jf, int k1, int k2, float ***c);
float bnldev(float pp, int n, long *idum);
float brent(float ax, float bx, float cx,
	float (*f)(float), float tol, float *xmin);
void broydn(float x[], int n, int *check,
	void (*vecfunc)(int, float [], float []));
void bsstep(float y[], float dydx[], int nv, float *xx, float htry,
	float eps, float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void caldat(long julian, int *mm, int *id, int *iyyy);
void chder(float a, float b, float c[], float cder[], int n);
float chebev(float a, float b, float c[], int m, float x);
void chebft(float a, float b, float c[], int n, float (*func)(float));
void chebpc(float c[], float d[], int n);
void chint(float a, float b, float c[], float cint[], int n);
float chixy(float bang);
void choldc(float **a, int n, float p[]);
void cholsl(float **a, int n, float p[], float b[], float x[]);
void chsone(float bins[], float ebins[], int nbins, int knstrn,
	float *df, float *chsq, float *prob);
void chstwo(float bins1[], float bins2[], int nbins, int knstrn,
	float *df, float *chsq, float *prob);
void cisi(float x, float *ci, float *si);
void cntab1(int **nn, int ni, int nj, float *chisq,
	float *df, float *prob, float *cramrv, float *ccc);
void cntab2(int **nn, int ni, int nj, float *h, float *hx, float *hy,
	float *hygx, float *hxgy, float *uygx, float *uxgy, float *uxy);
void convlv(float data[], unsigned long n, float respns[], unsigned long m,
	int isign, float ans[]);
void copy(double **aout, double **ain, int n);
void correl(float data1[], float data2[], unsigned long n, float ans[]);
void cosft(float y[], int n, int isign);
void cosft1(float y[], int n);
void cosft2(float y[], int n, int isign);
void covsrt(float **covar, int ma, int ia[], int mfit);
void crank(unsigned long n, float w[], float *s);
void cyclic(float a[], float b[], float c[], float alpha, float beta,
	float r[], float x[], unsigned long n);
void daub4(float a[], unsigned long n, int isign);
float dawson(float x);
float dbrent(float ax, float bx, float cx,
	float (*f)(float), float (*df)(float), float tol, float *xmin);
void ddpoly(float c[], int nc, float x, float pd[], int nd);
int decchk(char string[], int n, char *ch);
void derivs(float x, float y[], float dydx[]);
float df1dim(float x);
void dfour1(double data[], unsigned long nn, int isign);
void dfpmin(float p[], int n, float gtol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []));
float dfridr(float (*func)(float), float x, float h, float *err);
void dftcor(float w, float delta, float a, float b, float endpts[],
	float *corre, float *corim, float *corfac);
void dftint(float (*func)(float), float a, float b, float w,
	float *cosint, float *sinint);
void difeq(int k, int k1, int k2, int jsf, int is1, int isf,
	int indexv[], int ne, float **s, float **y);
void dlinmin(float p[], float xi[], int n, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float[]));
double dpythag(double a, double b);
void drealft(double data[], unsigned long n, int isign);
void dsprsax(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n);
void dsprstx(double sa[], unsigned long ija[], double x[], double b[],
	unsigned long n);
void dsvbksb(double **u, double w[], double **v, int m, int n, double b[],
	double x[]);
void dsvdcmp(double **a, int m, int n, double w[], double **v);
void eclass(int nf[], int n, int lista[], int listb[], int m);
void eclazz(int nf[], int n, int (*equiv)(int, int));
float ei(float x);
void eigsrt(float d[], float **v, int n);
float elle(float phi, float ak);
float ellf(float phi, float ak);
float ellpi(float phi, float en, float ak);
void elmhes(float **a, int n);
float erfcc(float x);
float erff(float x);
float erffc(float x);
void eulsum(float *sum, float term, int jterm, float wksp[]);
float evlmem(float fdt, float d[], int m, float xms);
float expdev(long *idum);
float expint(int n, float x);
float f1(float x);
float f1dim(float x);
float f2(float y);
float f3(float z);
float factln(int n);
float factrl(int n);
void fasper(float x[], float y[], unsigned long n, float ofac, float hifac,
	float wk1[], float wk2[], unsigned long nwk, unsigned long *nout,
	unsigned long *jmax, float *prob);
void fdjac(int n, float x[], float fvec[], float **df,
	void (*vecfunc)(int, float [], float []));
void fgauss(float x, float a[], float *y, float dyda[], int na);
void fill0(double **u, int n);
void fit(float x[], float y[], int ndata, float sig[], int mwt,
	float *a, float *b, float *siga, float *sigb, float *chi2, float *q);
void fitexy(float x[], float y[], int ndat, float sigx[], float sigy[],
	float *a, float *b, float *siga, float *sigb, float *chi2, float *q);
void fixrts(float d[], int m);
void fleg(float x, float pl[], int nl);
void flmoon(int n, int nph, long *jd, float *frac);
float fmin(float x[]);
void four1(float data[], unsigned long nn, int isign);
void fourew(FILE *file[5], int *na, int *nb, int *nc, int *nd);
void fourfs(FILE *file[5], unsigned long nn[], int ndim, int isign);
void fourn(float data[], unsigned long nn[], int ndim, int isign);
void fpoly(float x, float p[], int np);
void fred2(int n, float a, float b, float t[], float f[], float w[],
	float (*g)(float), float (*ak)(float, float));
float fredin(float x, int n, float a, float b, float t[], float f[], float w[],
	float (*g)(float), float (*ak)(float, float));
void frenel(float x, float *s, float *c);
void frprmn(float p[], int n, float ftol, int *iter, float *fret,
	float (*func)(float []), void (*dfunc)(float [], float []));
void ftest(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *f, float *prob);
float gamdev(int ia, long *idum);
float gammln(float xx);
float gammp(float a, float x);
float gammq(float a, float x);
float gasdev(long *idum);
void gaucof(int n, float a[], float b[], float amu0, float x[], float w[]);
void gauher(float x[], float w[], int n);
void gaujac(float x[], float w[], int n, float alf, float bet);
void gaulag(float x[], float w[], int n, float alf);
void gauleg(float x1, float x2, float x[], float w[], int n);
void gaussj(float **a, int n, float **b, int m);
void gcf(float *gammcf, float a, float x, float *gln);
float golden(float ax, float bx, float cx, float (*f)(float), float tol,
	float *xmin);
void gser(float *gamser, float a, float x, float *gln);
void hpsel(unsigned long m, unsigned long n, float arr[], float heap[]);
void hpsort(unsigned long n, float ra[]);
void hqr(float **a, int n, float wr[], float wi[]);
void hufapp(unsigned long index[], unsigned long nprob[], unsigned long n,
	unsigned long i);
void hufdec(unsigned long *ich, unsigned char *code, unsigned long lcode,
	unsigned long *nb, huffcode *hcode);
void hufenc(unsigned long ich, unsigned char **codep, unsigned long *lcode,
	unsigned long *nb, huffcode *hcode);
void hufmak(unsigned long nfreq[], unsigned long nchin, unsigned long *ilong,
	unsigned long *nlong, huffcode *hcode);
void hunt(float xx[], unsigned long n, float x, unsigned long *jlo);
void hypdrv(float s, float yy[], float dyyds[]);
fcomplex hypgeo(fcomplex a, fcomplex b, fcomplex c, fcomplex z);
void hypser(fcomplex a, fcomplex b, fcomplex c, fcomplex z,
	fcomplex *series, fcomplex *deriv);
unsigned short icrc(unsigned short crc, unsigned char *bufptr,
	unsigned long len, short jinit, int jrev);
unsigned short icrc1(unsigned short crc, unsigned char onech);
unsigned long igray(unsigned long n, int is);
void iindexx(unsigned long n, long arr[], unsigned long indx[]);
void indexx(unsigned long n, float arr[], unsigned long indx[]);
void interp(double **uf, double **uc, int nf);
int irbit1(unsigned long *iseed);
int irbit2(unsigned long *iseed);
void jacobi(float **a, int n, float d[], float **v, int *nrot);
void jacobn(float x, float y[], float dfdx[], float **dfdy, int n);
long julday(int mm, int id, int iyyy);
void kendl1(float data1[], float data2[], unsigned long n, float *tau, float *z,
	float *prob);
void kendl2(float **tab, int i, int j, float *tau, float *z, float *prob);
void kermom(double w[], double y, int m);
void ks2d1s(float x1[], float y1[], unsigned long n1,
	void (*quadvl)(float, float, float *, float *, float *, float *),
	float *d1, float *prob);
void ks2d2s(float x1[], float y1[], unsigned long n1, float x2[], float y2[],
	unsigned long n2, float *d, float *prob);
void ksone(float data[], unsigned long n, float (*func)(float), float *d,
	float *prob);
void kstwo(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *d, float *prob);
void laguer(fcomplex a[], int m, fcomplex *x, int *its);
void lfit(float x[], float y[], float sig[], int ndat, float a[], int ia[],
	int ma, float **covar, float *chisq, void (*funcs)(float, float [], int));
void linbcg(unsigned long n, double b[], double x[], int itol, double tol,
	 int itmax, int *iter, double *err);
void linmin(float p[], float xi[], int n, float *fret,
	float (*func)(float []));
void lnsrch(int n, float xold[], float fold, float g[], float p[], float x[],
	 float *f, float stpmax, int *check, float (*func)(float []));
void load(float x1, float v[], float y[]);
void load1(float x1, float v1[], float y[]);
void load2(float x2, float v2[], float y[]);
void locate(float xx[], unsigned long n, float x, unsigned long *j);
void lop(double **out, double **u, int n);
void lubksb(float **a, int n, int *indx, float b[]);
void ludcmp(float **a, int n, int *indx, float *d);
void machar(int *ibeta, int *it, int *irnd, int *ngrd,
	int *machep, int *negep, int *iexp, int *minexp, int *maxexp,
	float *eps, float *epsneg, float *xmin, float *xmax);
void matadd(double **a, double **b, double **c, int n);
void matsub(double **a, double **b, double **c, int n);
void medfit(float x[], float y[], int ndata, float *a, float *b, float *abdev);
void memcof(float data[], int n, int m, float *xms, float d[]);
int metrop(float de, float t);
void mgfas(double **u, int n, int maxcyc);
void mglin(double **u, int n, int ncycle);
float midexp(float (*funk)(float), float aa, float bb, int n);
float midinf(float (*funk)(float), float aa, float bb, int n);
float midpnt(float (*func)(float), float a, float b, int n);
float midsql(float (*funk)(float), float aa, float bb, int n);
float midsqu(float (*funk)(float), float aa, float bb, int n);
void miser(float (*func)(float []), float regn[], int ndim, unsigned long npts,
	float dith, float *ave, float *var);
void mmid(float y[], float dydx[], int nvar, float xs, float htot,
	int nstep, float yout[], void (*derivs)(float, float[], float[]));
void mnbrak(float *ax, float *bx, float *cx, float *fa, float *fb,
	float *fc, float (*func)(float));
void mnewt(int ntrial, float x[], int n, float tolx, float tolf);
void moment(float data[], int n, float *ave, float *adev, float *sdev,
	float *var, float *skew, float *curt);
void mp2dfr(unsigned char a[], unsigned char s[], int n, int *m);
void mpadd(unsigned char w[], unsigned char u[], unsigned char v[], int n);
void mpdiv(unsigned char q[], unsigned char r[], unsigned char u[],
	unsigned char v[], int n, int m);
void mpinv(unsigned char u[], unsigned char v[], int n, int m);
void mplsh(unsigned char u[], int n);
void mpmov(unsigned char u[], unsigned char v[], int n);
void mpmul(unsigned char w[], unsigned char u[], unsigned char v[], int n,
	int m);
void mpneg(unsigned char u[], int n);
void mppi(int n);
void mprove(float **a, float **alud, int n, int indx[], float b[],
	float x[]);
void mpsad(unsigned char w[], unsigned char u[], int n, int iv);
void mpsdv(unsigned char w[], unsigned char u[], int n, int iv, int *ir);
void mpsmu(unsigned char w[], unsigned char u[], int n, int iv);
void mpsqrt(unsigned char w[], unsigned char u[], unsigned char v[], int n,
	int m);
void mpsub(int *is, unsigned char w[], unsigned char u[], unsigned char v[],
	int n);
void mrqcof(float x[], float y[], float sig[], int ndata, float a[],
	int ia[], int ma, float **alpha, float beta[], float *chisq,
	void (*funcs)(float, float [], float *, float [], int));
void mrqmin(float x[], float y[], float sig[], int ndata, float a[],
	int ia[], int ma, float **covar, float **alpha, float *chisq,
	void (*funcs)(float, float [], float *, float [], int), float *alamda);
void newt(float x[], int n, int *check,
	void (*vecfunc)(int, float [], float []));
void odeint(float ystart[], int nvar, float x1, float x2,
	float eps, float h1, float hmin, int *nok, int *nbad,
	void (*derivs)(float, float [], float []),
	void (*rkqs)(float [], float [], int, float *, float, float,
	float [], float *, float *, void (*)(float, float [], float [])));
void orthog(int n, float anu[], float alpha[], float beta[], float a[],
	float b[]);
void pade(double cof[], int n, float *resid);
void pccheb(float d[], float c[], int n);
void pcshft(float a, float b, float d[], int n);
void pearsn(float x[], float y[], unsigned long n, float *r, float *prob,
	float *z);
void period(float x[], float y[], int n, float ofac, float hifac,
	float px[], float py[], int np, int *nout, int *jmax, float *prob);
void piksr2(int n, float arr[], float brr[]);
void piksrt(int n, float arr[]);
void pinvs(int ie1, int ie2, int je1, int jsf, int jc1, int k,
	float ***c, float **s);
float plgndr(int l, int m, float x);
float poidev(float xm, long *idum);
void polcoe(float x[], float y[], int n, float cof[]);
void polcof(float xa[], float ya[], int n, float cof[]);
void poldiv(float u[], int n, float v[], int nv, float q[], float r[]);
void polin2(float x1a[], float x2a[], float **ya, int m, int n,
	float x1, float x2, float *y, float *dy);
void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
void powell(float p[], float **xi, int n, float ftol, int *iter, float *fret,
	float (*func)(float []));
void predic(float data[], int ndata, float d[], int m, float future[], int nfut);
float probks(float alam);
void psdes(unsigned long *lword, unsigned long *irword);
void pwt(float a[], unsigned long n, int isign);
void pwtset(int n);
float pythag(float a, float b);
void pzextr(int iest, float xest, float yest[], float yz[], float dy[],
	int nv);
float qgaus(float (*func)(float), float a, float b);
void qrdcmp(float **a, int n, float *c, float *d, int *sing);
float qromb(float (*func)(float), float a, float b);
float qromo(float (*func)(float), float a, float b,
	float (*choose)(float (*)(float), float, float, int));
void qroot(float p[], int n, float *b, float *c, float eps);
void qrsolv(float **a, int n, float c[], float d[], float b[]);
void qrupdt(float **r, float **qt, int n, float u[], float v[]);
float qsimp(float (*func)(float), float a, float b);
float qtrap(float (*func)(float), float a, float b);
float quad3d(float (*func)(float, float, float), float x1, float x2);
void quadct(float x, float y, float xx[], float yy[], unsigned long nn,
	float *fa, float *fb, float *fc, float *fd);
void quadmx(float **a, int n);
void quadvl(float x, float y, float *fa, float *fb, float *fc, float *fd);
float ran0(long *idum);
float ran1(long *idum);
float ran2(long *idum);
float ran3(long *idum);
float ran4(long *idum);
void rank(unsigned long n, unsigned long indx[], unsigned long irank[]);
void ranpt(float pt[], float regn[], int n);
void ratint(float xa[], float ya[], int n, float x, float *y, float *dy);
void ratlsq(double (*fn)(double), double a, double b, int mm, int kk,
	double cof[], double *dev);
double ratval(double x, double cof[], int mm, int kk);
float rc(float x, float y);
float rd(float x, float y, float z);
void realft(float data[], unsigned long n, int isign);
void rebin(float rc, int nd, float r[], float xin[], float xi[]);
void red(int iz1, int iz2, int jz1, int jz2, int jm1, int jm2, int jmf,
	int ic1, int jc1, int jcf, int kc, float ***c, float **s);
void relax(double **u, double **rhs, int n);
void relax2(double **u, double **rhs, int n);
void resid(double **res, double **u, double **rhs, int n);
float revcst(float x[], float y[], int iorder[], int ncity, int n[]);
void reverse(int iorder[], int ncity, int n[]);
float rf(float x, float y, float z);
float rj(float x, float y, float z, float p);
void rk4(float y[], float dydx[], int n, float x, float h, float yout[],
	void (*derivs)(float, float [], float []));
void rkck(float y[], float dydx[], int n, float x, float h,
	float yout[], float yerr[], void (*derivs)(float, float [], float []));
void rkdumb(float vstart[], int nvar, float x1, float x2, int nstep,
	void (*derivs)(float, float [], float []));
void rkqs(float y[], float dydx[], int n, float *x,
	float htry, float eps, float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void rlft3(float ***data, float **speq, unsigned long nn1,
	unsigned long nn2, unsigned long nn3, int isign);
float rofunc(float b);
void rotate(float **r, float **qt, int n, int i, float a, float b);
void rsolv(float **a, int n, float d[], float b[]);
void rstrct(double **uc, double **uf, int nc);
float rtbis(float (*func)(float), float x1, float x2, float xacc);
float rtflsp(float (*func)(float), float x1, float x2, float xacc);
float rtnewt(void (*funcd)(float, float *, float *), float x1, float x2,
	float xacc);
float rtsafe(void (*funcd)(float, float *, float *), float x1, float x2,
	float xacc);
float rtsec(float (*func)(float), float x1, float x2, float xacc);
void rzextr(int iest, float xest, float yest[], float yz[], float dy[], int nv);
void savgol(float c[], int np, int nl, int nr, int ld, int m);
void score(float xf, float y[], float f[]);
void scrsho(float (*fx)(float));
float select(unsigned long k, unsigned long n, float arr[]);
float selip(unsigned long k, unsigned long n, float arr[]);
void shell(unsigned long n, float a[]);
void shoot(int n, float v[], float f[]);
void shootf(int n, float v[], float f[]);
void simp1(float **a, int mm, int ll[], int nll, int iabf, int *kp,
	float *bmax);
void simp2(float **a, int n, int l2[], int nl2, int *ip, int kp, float *q1);
void simp3(float **a, int i1, int k1, int ip, int kp);
void simplx(float **a, int m, int n, int m1, int m2, int m3, int *icase,
	int izrov[], int iposv[]);
void simpr(float y[], float dydx[], float dfdx[], float **dfdy,
	int n, float xs, float htot, int nstep, float yout[],
	void (*derivs)(float, float [], float []));
void sinft(float y[], int n);
void slvsm2(double **u, double **rhs);
void slvsml(double **u, double **rhs);
void sncndn(float uu, float emmc, float *sn, float *cn, float *dn);
double snrm(unsigned long n, double sx[], int itol);
void sobseq(int *n, float x[]);
void solvde(int itmax, float conv, float slowc, float scalv[],
	int indexv[], int ne, int nb, int m, float **y, float ***c, float **s);
void sor(double **a, double **b, double **c, double **d, double **e,
	double **f, double **u, int jmax, double rjac);
void sort(unsigned long n, float arr[]);
void sort2(unsigned long n, float arr[], float brr[]);
void sort3(unsigned long n, float ra[], float rb[], float rc[]);
void spctrm(FILE *fp, float p[], int m, int k, int ovrlap);
void spear(float data1[], float data2[], unsigned long n, float *d, float *zd,
	float *probd, float *rs, float *probrs);
void sphbes(int n, float x, float *sj, float *sy, float *sjp, float *syp);
void splie2(float x1a[], float x2a[], float **ya, int m, int n, float **y2a);
void splin2(float x1a[], float x2a[], float **ya, float **y2a, int m, int n,
	float x1, float x2, float *y);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);
void spread(float y, float yy[], unsigned long n, float x, int m);
void sprsax(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n);
void sprsin(float **a, int n, float thresh, unsigned long nmax, float sa[],
	unsigned long ija[]);
void sprspm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[],
	float sc[], unsigned long ijc[]);
void sprstm(float sa[], unsigned long ija[], float sb[], unsigned long ijb[],
	float thresh, unsigned long nmax, float sc[], unsigned long ijc[]);
void sprstp(float sa[], unsigned long ija[], float sb[], unsigned long ijb[]);
void sprstx(float sa[], unsigned long ija[], float x[], float b[],
	unsigned long n);
void stifbs(float y[], float dydx[], int nv, float *xx,
	float htry, float eps, float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void stiff(float y[], float dydx[], int n, float *x,
	float htry, float eps, float yscal[], float *hdid, float *hnext,
	void (*derivs)(float, float [], float []));
void stoerm(float y[], float d2y[], int nv, float xs,
	float htot, int nstep, float yout[],
	void (*derivs)(float, float [], float []));
void svbksb(float **u, float w[], float **v, int m, int n, float b[],
	float x[]);
void svdcmp(float **a, int m, int n, float w[], float **v);
void svdfit(float x[], float y[], float sig[], int ndata, float a[],
	int ma, float **u, float **v, float w[], float *chisq,
	void (*funcs)(float, float [], int));
void svdvar(float **v, int ma, float w[], float **cvm);
void toeplz(float r[], float x[], float y[], int n);
void tptest(float data1[], float data2[], unsigned long n, float *t, float *prob);
void tqli(float d[], float e[], int n, float **z);
float trapzd(float (*func)(float), float a, float b, int n);
void tred2(float **a, int n, float d[], float e[]);
void tridag(float a[], float b[], float c[], float r[], float u[],
	unsigned long n);
float trncst(float x[], float y[], int iorder[], int ncity, int n[]);
void trnspt(int iorder[], int ncity, int n[]);
void ttest(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *t, float *prob);
void tutest(float data1[], unsigned long n1, float data2[], unsigned long n2,
	float *t, float *prob);
void twofft(float data1[], float data2[], float fft1[], float fft2[],
	unsigned long n);
void vander(double x[], double w[], double q[], int n);
void vegas(float regn[], int ndim, float (*fxn)(float [], float), int init,
	unsigned long ncall, int itmx, int nprn, float *tgral, float *sd,
	float *chi2a);
void voltra(int n, int m, float t0, float h, float *t, float **f,
	float (*g)(int, float), float (*ak)(int, int, float, float));
void wt1(float a[], unsigned long n, int isign,
	void (*wtstep)(float [], unsigned long, int));
void wtn(float a[], unsigned long nn[], int ndim, int isign,
	void (*wtstep)(float [], unsigned long, int));
void wwghts(float wghts[], int n, float h,
	void (*kermom)(double [], double ,int));
int zbrac(float (*func)(float), float *x1, float *x2);
void zbrak(float (*fx)(float), float x1, float x2, int n, float xb1[],
	float xb2[], int *nb);
float zbrent(float (*func)(float), float x1, float x2, float tol);
void zrhqr(float a[], int m, float rtr[], float rti[]);
float zriddr(float (*func)(float), float x1, float x2, float xacc);
void zroots(fcomplex a[], int m, fcomplex roots[], int polish);

#else /* ANSI */
/* traditional - K&R */

void addint();
void airy();
void amebsa();
void amoeba();
float amotry();
float amotsa();
void anneal();
double anorm2();
void arcmak();
void arcode();
void arcsum();
void asolve();
void atimes();
void avevar();
void balanc();
void banbks();
void bandec();
void banmul();
void bcucof();
void bcuint();
void beschb();
float bessi();
float bessi0();
float bessi1();
void bessik();
float bessj();
float bessj0();
float bessj1();
void bessjy();
float bessk();
float bessk0();
float bessk1();
float bessy();
float bessy0();
float bessy1();
float beta();
float betacf();
float betai();
float bico();
void bksub();
float bnldev();
float brent();
void broydn();
void bsstep();
void caldat();
void chder();
float chebev();
void chebft();
void chebpc();
void chint();
float chixy();
void choldc();
void cholsl();
void chsone();
void chstwo();
void cisi();
void cntab1();
void cntab2();
void convlv();
void copy();
void correl();
void cosft();
void cosft1();
void cosft2();
void covsrt();
void crank();
void cyclic();
void daub4();
float dawson();
float dbrent();
void ddpoly();
int decchk();
void derivs();
float df1dim();
void dfour1();
void dfpmin();
float dfridr();
void dftcor();
void dftint();
void difeq();
void dlinmin();
double dpythag();
void drealft();
void dsprsax();
void dsprstx();
void dsvbksb();
void dsvdcmp();
void eclass();
void eclazz();
float ei();
void eigsrt();
float elle();
float ellf();
float ellpi();
void elmhes();
float erfcc();
float erff();
float erffc();
void eulsum();
float evlmem();
float expdev();
float expint();
float f1();
float f1dim();
float f2();
float f3();
float factln();
float factrl();
void fasper();
void fdjac();
void fgauss();
void fill0();
void fit();
void fitexy();
void fixrts();
void fleg();
void flmoon();
float fmin();
void four1();
void fourew();
void fourfs();
void fourn();
void fpoly();
void fred2();
float fredin();
void frenel();
void frprmn();
void ftest();
float gamdev();
float gammln();
float gammp();
float gammq();
float gasdev();
void gaucof();
void gauher();
void gaujac();
void gaulag();
void gauleg();
void gaussj();
void gcf();
float golden();
void gser();
void hpsel();
void hpsort();
void hqr();
void hufapp();
void hufdec();
void hufenc();
void hufmak();
void hunt();
void hypdrv();
fcomplex hypgeo();
void hypser();
unsigned short icrc();
unsigned short icrc1();
int igray();
void iindexx();
void indexx();
void interp();
int irbit1();
int irbit2();
void jacobi();
void jacobn();
long julday();
void kendl1();
void kendl2();
void kermom();
void ks2d1s();
void ks2d2s();
void ksone();
void kstwo();
void laguer();
void lfit();
void linbcg();
void linmin();
void lnsrch();
void load();
void load1();
void load2();
void locate();
void lop();
void lubksb();
void ludcmp();
void machar();
void matadd();
void matsub();
void medfit();
void memcof();
int metrop();
void mgfas();
void mglin();
float midexp();
float midinf();
float midpnt();
float midsql();
float midsqu();
void miser();
void mmid();
void mnbrak();
void mnewt();
void moment();
void mp2dfr();
void mpadd();
void mpdiv();
void mpinv();
void mplsh();
void mpmov();
void mpmul();
void mpneg();
void mppi();
void mprove();
void mpsad();
void mpsdv();
void mpsmu();
void mpsqrt();
void mpsub();
void mrqcof();
void mrqmin();
void newt();
void odeint();
void orthog();
void pade();
void pccheb();
void pcshft();
void pearsn();
void period();
void piksr2();
void piksrt();
void pinvs();
float plgndr();
float poidev();
void polcoe();
void polcof();
void poldiv();
void polin2();
void polint();
void powell();
void predic();
float probks();
void psdes();
void pwt();
void pwtset();
float pythag();
void pzextr();
float qgaus();
void qrdcmp();
float qromb();
float qromo();
void qroot();
void qrsolv();
void qrupdt();
float qsimp();
float qtrap();
float quad3d();
void quadct();
void quadmx();
void quadvl();
float ran0();
float ran1();
float ran2();
float ran3();
float ran4();
void rank();
void ranpt();
void ratint();
void ratlsq();
double ratval();
float rc();
float rd();
void realft();
void rebin();
void red();
void relax();
void relax2();
void resid();
float revcst();
void reverse();
float rf();
float rj();
void rk4();
void rkck();
void rkdumb();
void rkqs();
void rlft3();
float rofunc();
void rotate();
void rsolv();
void rstrct();
float rtbis();
float rtflsp();
float rtnewt();
float rtsafe();
float rtsec();
void rzextr();
void savgol();
void score();
void scrsho();
float select();
float selip();
void shell();
void shoot();
void shootf();
void simp1();
void simp2();
void simp3();
void simplx();
void simpr();
void sinft();
void slvsm2();
void slvsml();
void sncndn();
double snrm();
void sobseq();
void solvde();
void sor();
void sort();
void sort2();
void sort3();
void spctrm();
void spear();
void sphbes();
void splie2();
void splin2();
void spline();
void splint();
void spread();
void sprsax();
void sprsin();
void sprspm();
void sprstm();
void sprstp();
void sprstx();
void stifbs();
void stiff();
void stoerm();
void svbksb();
void svdcmp();
void svdfit();
void svdvar();
void toeplz();
void tptest();
void tqli();
float trapzd();
void tred2();
void tridag();
float trncst();
void trnspt();
void ttest();
void tutest();
void twofft();
void vander();
void vegas();
void voltra();
void wt1();
void wtn();
void wwghts();
int zbrac();
void zbrak();
float zbrent();
void zrhqr();
float zriddr();
void zroots();

#endif /* ANSI */

#endif /* _NR_H_ */
                                                                                                                                                                                                                                                                                 r8s/nrutil.c                                                                                        0000755 0000766 0000120 00000020326 07565204055 012742  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /* CAUTION: This is the ANSI C (only) version of the Numerical Recipes
   utility file nrutil.c.  Do not confuse this file with the same-named
   file nrutil.c that is supplied in the same subdirectory or archive
   as the header file nrutil.h.  *That* file contains both ANSI and
   traditional K&R versions, along with #ifdef macros to select the
   correct version.  *This* file contains only ANSI C.               */

#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#define NR_END 1
#define FREE_ARG char*

void nrerror(char error_text[])
/* Numerical Recipes standard error handler */
{
	fprintf(stderr,"Numerical Recipes run-time error...\n");
	fprintf(stderr,"%s\n",error_text);
	fprintf(stderr,"...now exiting to system...\n");
	exit(1);
}

float *vector(long nl, long nh)
/* allocate a float vector with subscript range v[nl..nh] */
{
	float *v;

	v=(float *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(float)));
	if (!v) nrerror("allocation failure in vector()");
	return v-nl+NR_END;
}

int *ivector(long nl, long nh)
/* allocate an int vector with subscript range v[nl..nh] */
{
	int *v;

	v=(int *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(int)));
	if (!v) nrerror("allocation failure in ivector()");
	return v-nl+NR_END;
}

unsigned char *cvector(long nl, long nh)
/* allocate an unsigned char vector with subscript range v[nl..nh] */
{
	unsigned char *v;

	v=(unsigned char *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(unsigned char)));
	if (!v) nrerror("allocation failure in cvector()");
	return v-nl+NR_END;
}

unsigned long *lvector(long nl, long nh)
/* allocate an unsigned long vector with subscript range v[nl..nh] */
{
	unsigned long *v;

	v=(unsigned long *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(long)));
	if (!v) nrerror("allocation failure in lvector()");
	return v-nl+NR_END;
}

double *dvector(long nl, long nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
	double *v;

	v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
	if (!v) nrerror("allocation failure in dvector()");
	return v-nl+NR_END;
}

float **matrix(long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t)((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(float *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

double **dmatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a double matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	double **m;

	/* allocate pointers to rows */
	m=(double **) malloc((size_t)((nrow+NR_END)*sizeof(double*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;

	/* allocate rows and set pointers to them */
	m[nrl]=(double *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(double)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

int **imatrix(long nrl, long nrh, long ncl, long nch)
/* allocate a int matrix with subscript range m[nrl..nrh][ncl..nch] */
{
	long i, nrow=nrh-nrl+1,ncol=nch-ncl+1;
	int **m;

	/* allocate pointers to rows */
	m=(int **) malloc((size_t)((nrow+NR_END)*sizeof(int*)));
	if (!m) nrerror("allocation failure 1 in matrix()");
	m += NR_END;
	m -= nrl;


	/* allocate rows and set pointers to them */
	m[nrl]=(int *) malloc((size_t)((nrow*ncol+NR_END)*sizeof(int)));
	if (!m[nrl]) nrerror("allocation failure 2 in matrix()");
	m[nrl] += NR_END;
	m[nrl] -= ncl;

	for(i=nrl+1;i<=nrh;i++) m[i]=m[i-1]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **submatrix(float **a, long oldrl, long oldrh, long oldcl, long oldch,
	long newrl, long newcl)
/* point a submatrix [newrl..][newcl..] to a[oldrl..oldrh][oldcl..oldch] */
{
	long i,j,nrow=oldrh-oldrl+1,ncol=oldcl-newcl;
	float **m;

	/* allocate array of pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in submatrix()");
	m += NR_END;
	m -= newrl;

	/* set pointers to rows */
	for(i=oldrl,j=newrl;i<=oldrh;i++,j++) m[j]=a[i]+ncol;

	/* return pointer to array of pointers to rows */
	return m;
}

float **convert_matrix(float *a, long nrl, long nrh, long ncl, long nch)
/* allocate a float matrix m[nrl..nrh][ncl..nch] that points to the matrix
declared in the standard C manner as a[nrow][ncol], where nrow=nrh-nrl+1
and ncol=nch-ncl+1. The routine should be called with the address
&a[0][0] as the first argument. */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1;
	float **m;

	/* allocate pointers to rows */
	m=(float **) malloc((size_t) ((nrow+NR_END)*sizeof(float*)));
	if (!m) nrerror("allocation failure in convert_matrix()");
	m += NR_END;
	m -= nrl;

	/* set pointers to rows */
	m[nrl]=a-ncl;
	for(i=1,j=nrl+1;i<nrow;i++,j++) m[j]=m[j-1]+ncol;
	/* return pointer to array of pointers to rows */
	return m;
}

float ***f3tensor(long nrl, long nrh, long ncl, long nch, long ndl, long ndh)
/* allocate a float 3tensor with range t[nrl..nrh][ncl..nch][ndl..ndh] */
{
	long i,j,nrow=nrh-nrl+1,ncol=nch-ncl+1,ndep=ndh-ndl+1;
	float ***t;

	/* allocate pointers to pointers to rows */
	t=(float ***) malloc((size_t)((nrow+NR_END)*sizeof(float**)));
	if (!t) nrerror("allocation failure 1 in f3tensor()");
	t += NR_END;
	t -= nrl;

	/* allocate pointers to rows and set pointers to them */
	t[nrl]=(float **) malloc((size_t)((nrow*ncol+NR_END)*sizeof(float*)));
	if (!t[nrl]) nrerror("allocation failure 2 in f3tensor()");
	t[nrl] += NR_END;
	t[nrl] -= ncl;

	/* allocate rows and set pointers to them */
	t[nrl][ncl]=(float *) malloc((size_t)((nrow*ncol*ndep+NR_END)*sizeof(float)));
	if (!t[nrl][ncl]) nrerror("allocation failure 3 in f3tensor()");
	t[nrl][ncl] += NR_END;
	t[nrl][ncl] -= ndl;

	for(j=ncl+1;j<=nch;j++) t[nrl][j]=t[nrl][j-1]+ndep;
	for(i=nrl+1;i<=nrh;i++) {
		t[i]=t[i-1]+ncol;
		t[i][ncl]=t[i-1][ncl]+ncol*ndep;
		for(j=ncl+1;j<=nch;j++) t[i][j]=t[i][j-1]+ndep;
	}

	/* return pointer to array of pointers to rows */
	return t;
}

void free_vector(float *v, long nl, long nh)
/* free a float vector allocated with vector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_ivector(int *v, long nl, long nh)
/* free an int vector allocated with ivector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_cvector(unsigned char *v, long nl, long nh)
/* free an unsigned char vector allocated with cvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_lvector(unsigned long *v, long nl, long nh)
/* free an unsigned long vector allocated with lvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_dvector(double *v, long nl, long nh)
/* free a double vector allocated with dvector() */
{
	free((FREE_ARG) (v+nl-NR_END));
}

void free_matrix(float **m, long nrl, long nrh, long ncl, long nch)
/* free a float matrix allocated by matrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch)
/* free a double matrix allocated by dmatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_imatrix(int **m, long nrl, long nrh, long ncl, long nch)
/* free an int matrix allocated by imatrix() */
{
	free((FREE_ARG) (m[nrl]+ncl-NR_END));
	free((FREE_ARG) (m+nrl-NR_END));
}

void free_submatrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a submatrix allocated by submatrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_convert_matrix(float **b, long nrl, long nrh, long ncl, long nch)
/* free a matrix allocated by convert_matrix() */
{
	free((FREE_ARG) (b+nrl-NR_END));
}

void free_f3tensor(float ***t, long nrl, long nrh, long ncl, long nch,
	long ndl, long ndh)
/* free a float f3tensor allocated by f3tensor() */
{
	free((FREE_ARG) (t[nrl][ncl]+ndl-NR_END));
	free((FREE_ARG) (t[nrl]+ncl-NR_END));
	free((FREE_ARG) (t+nrl-NR_END));
}
                                                                                                                                                                                                                                                                                                          r8s/objective.h                                                                                     0000644 0000766 0000120 00000000671 07565204055 013402  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include <math.h>
#include <stdlib.h>
#include <stdio.h>

double triadLike(double t1, double t2, double t3, 
					double *xt, double L1, double L2, double L3);
int feasible(double p[]);
double objective(double p[]);
double triadObs(double t1, double t2, double t3, double tint, 
	double L1, double L2, double L3);
double penalty(double p[]);
double addconstr(double x);
double BranchLike(double rate, double timeLength, double charLength);

                                                                       r8s/parse_command.c                                                                                 0000644 0000766 0000120 00000004515 07565204055 014234  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #define INTEGER	    0
#define	REAL	    1
#define FLAG	    2
#define STRING	    3
#define CHARACTER   4
#define MAX_COMMANDS 25
struct cList {
		int variable_type;
		char * option_name;
		void * variable_address;
		};



/*
 * Processes a command consisting of 'option=value' strings separated by spaces and ending in ';'.
 * The array of structures 'comList' contains the syntax.  Each element of the array has the three
 * structure members corresponding to the option name string that is looked for,  the type of variable
 * it is,  and the address of a variable that the 'value' will be stored in.  For integers and reals this
 * member is just a pointer to the stored values.  For type FLAG,  the parser will look for 'option=YES' or
 * 'option=NO' and store an integer 1 or 0 respectively.  For type STRING,  the pointer points to a location
 * where a pointer to the string will be stored.  Hence to retrieve that string we have to dereference twice.
  */
  
void dummy(void)
{
char *tHndl;
double min_age, max_age; 

 
  
struct cList b[MAX_COMMANDS] =
    {
    {STRING, "TAXON", &tHndl}, 
    {REAL, "MIN_AGE", &min_age}, 
    {REAL, "MAX_AGE", &max_age} 
    };
}

void parse_command(struct cList comList[], int ncommands)
{
char * localTokenPtr,  *dummy;
extern char * aTokenPtr;
int i;
if (ncommands > MAX_COMMANDS)
    fatal("Too many commands in parse_command\n");
while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		for (i=0;i<ncommands;i++)
		    if (isEqual(aTokenPtr,   (comList[i].option_name)))
			{
			if (parse_assignment(comList[i].option_name,&localTokenPtr))
			    {
			    switch (comList[i].variable_type)
				{
				 case INTEGER: *(int *)(comList[i].variable_address)=strtod(localTokenPtr,&dummy);
					myFree(localTokenPtr);
					break;
				 case REAL: *(double *)(comList[i].variable_address)=strtod(localTokenPtr,&dummy); 
					myFree(localTokenPtr);
					break;
				 case STRING: *(char**)(comList[i].variable_address)=localTokenPtr; /* don't free! */
					break;
				 case CHARACTER: *(char *)(comList[i].variable_address)=strtod(localTokenPtr,&dummy); 
					myFree(localTokenPtr);
					break;
				 case FLAG: 
					if (isEqual(localTokenPtr, "YES"))
						*(int *)(comList[i].variable_address)=1;
					else
						*(int *)(comList[i].variable_address)=0;
					myFree(localTokenPtr);
					break;
				}
			    }
			}
			

		
		}

 
 return;   
    
}
                                                                                                                                                                                   r8s/penalty.c                                                                                       0000644 0000766 0000120 00000016331 07565204055 013077  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "TreeUtils.h"
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
                                                                                                                                                                                                                                                                                                       r8s/penalty.h                                                                                       0000644 0000766 0000120 00000000444 07565204055 013102  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  double penalty_rate(double p[]);
void debug_check_feasible(NODETYPE *node);
void check_feasible(NODETYPE *node);
double traversePenalty(NODETYPE *node, double p[]);
double addconstr(double x);
double penalty(double p[]);
int check_initial_point(double (*objective)(double p[]),  double p[]);
                                                                                                                                                                                                                            r8s/powell.c                                                                                        0000644 0000766 0000120 00000025264 10357560066 012732  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /* next two are  TEMPORARY to debug LNSRCH */
#include "TreeUtils.h"
#include "DrawTree.h"
#include "MyUtilities.h"

#include "structures.h"
#include <math.h>
#define NRANSI
#include "NRCvectorUtils.h"
#include "powell.h"
#define ITMAX 3000
static float sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)    /*...from Powell */
#define SIGN(a,b) ((b)>0.0 ? fabs(a):-fabs(a))
#define MAX(a,b) ((a) >(b) ? (a):(b))


double          gContractFactor=0.10; /* this is old crap */
int             gMaxContractIter=10;
int		gPowellTrace=0;


double *pcom,*xicom,(*nrfunc)(double []);
int ncom;

int powell1(double p[], double **xi, int n, double ftol, int *iter, double *fret,
	double (*func)(double []))
{
extern	int		gmaxPowellIter,powellMode;
extern StackPtr gFStack,gPStack,gTestStack;
	int i,itmax,ibig,j;
	double del,fp,fptt,t,*pt,*ptt,*xit,*pdif;
	itmax=*iter;
	pt=vector(1,n);
	ptt=vector(1,n);
	xit=vector(1,n);
/*	pdif=vector(1,n);*/
	*fret=(*func)(p);
	for (j=1;j<=n;j++) pt[j]=p[j];
	for (*iter=1;;++(*iter)) {
		if(gPowellTrace)
		    {
		     printf("TRACE (MODE=%i)(Powell iteration %i)(start)\n", powellMode,*iter);
		     for (i=1;i<=n;i++)
			 printf("p[%i] %f\n",i, p[i]);
		     printf("	Objective function value = %f\n", (*func)(p));
		    }
	
		fp=(*fret);
		ibig=0;
		del=0.0;
		for (i=1;i<=n;i++) {
			for (j=1;j<=n;j++) xit[j]=xi[j][i];
			fptt=(*fret);
			linmin1(p,xit,n,fret,func);
			if (fabs(fptt-(*fret)) > del) {
				del=fabs(fptt-(*fret));
				ibig=i;
			}
		}
/* NEW ADDITIONS */

/*
for(j=1;j<=n;j++)pdif[j]=pt[j]-p[j];
pushD(gFStack,fp-(*fret));
pushD(gPStack,norm(pdif,1,n));
*/
/********/
		if (2.0*fabs(fp-(*fret)) <= ftol*(fabs(fp)+fabs(*fret))) {
			free_vector(xit,1,n);
			free_vector(ptt,1,n);
			free_vector(pt,1,n);

			return 1;
		}
		if (*iter >= itmax) 
			{
			printf("powell1 exceeding maximum iterations.\n");
			return 0;
			}
		for (j=1;j<=n;j++) {
			ptt[j]=2.0*p[j]-pt[j];
			xit[j]=p[j]-pt[j];
			pt[j]=p[j];
		}
		fptt=(*func)(ptt);
		if (fptt < fp) {
			t=2.0*(fp-2.0*(*fret)+fptt)*SQR(fp-(*fret)-del)-del*SQR(fp-fptt);
			if (t < 0.0) {
				linmin1(p,xit,n,fret,func);
				for (j=1;j<=n;j++) {
					xi[j][ibig]=xi[j][n];
					xi[j][n]=xit[j];
				}
			}
		}
	}
}
#undef NRANSI
#define NRANSI
#define CGOLD 0.3819660
#define ZEPS 1.0e-10
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

double brent1(double ax, double bx, double cx, double (*f)(double), double tol,
	double *xmin)
{
	int iter;
	double a,b,d,etemp,fu,fv,fw,fx,p,q,r,tol1,tol2,u,v,w,x,xm;
	double e=0.0;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=(*f)(x);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol2=2.0*(tol1=tol*fabs(x)+ZEPS);
		if (fabs(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (fabs(e) > tol1) {
			r=(x-w)*(fx-fv);
			q=(x-v)*(fx-fw);
			p=(x-v)*q-(x-w)*r;
			q=2.0*(q-r);
			if (q > 0.0) p = -p;
			q=fabs(q);
			etemp=e;
			e=d;
			if (fabs(p) >= fabs(0.5*q*etemp) || p <= q*(a-x) || p >= q*(b-x))
				d=CGOLD*(e=(x >= xm ? a-x : b-x));
			else {
				d=p/q;
				u=x+d;
				if (u-a < tol2 || b-u < tol2)
					d=SIGN(tol1,xm-x);
			}
		} else {
			d=CGOLD*(e=(x >= xm ? a-x : b-x));
		}
		u=(fabs(d) >= tol1 ? x+d : x+SIGN(tol1,d));
		fu=(*f)(u);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			SHFT(v,w,x,u)
			SHFT(fv,fw,fx,fu)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				v=w;
				w=u;
				fv=fw;
				fw=fu;
			} else if (fu <= fv || v == x || v == w) {
				v=u;
				fv=fu;
			}
		}
	}
	nrerror("Too many iterations in brent1");
	*xmin=x;
	return fx;
}
#undef CGOLD
#undef ZEPS
#undef SHFT
#undef NRANSI
#define NRANSI
#define TOL 2.0e-4


void linmin1(double p[], double xi[], int n, double *fret, double (*func)(double []))
{
	int j;
	double xx,xmin,fx,fb,fa,bx,ax;

	ncom=n;
	pcom=vector(1,n);
	xicom=vector(1,n);
	nrfunc=func;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	mnbrak1(&ax,&xx,&bx,&fa,&fx,&fb,f1dim1);
	*fret=brent1(ax,xx,bx,f1dim1,TOL,&xmin);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free_vector(xicom,1,n);
	free_vector(pcom,1,n);
}
#undef TOL
#undef NRANSI
#define NRANSI
#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak1(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
	double (*func)(double))
{
	double ulim,u,r,q,fu,dum;

	*fa=(*func)(*ax);
	*fb=(*func)(*bx);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=(*func)(*cx);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(MAX(fabs(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=(*func)(u);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,(*func)(u))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=(*func)(u);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=(*func)(u);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}
#undef GOLD
#undef GLIMIT
#undef TINY
#undef SHFT
#undef NRANSI
#define NRANSI

double f1dim1(double x)
{
	int j;
	double f,*xt;

	xt=vector(1,ncom);
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=(*nrfunc)(xt);
	free_vector(xt,1,ncom);
	return f;
}
#undef NRANSI


#define NRANSI
#define EPS 1.0e-10
#define FREEALL free_vector(xi,1,n);free_vector(h,1,n);free_vector(g,1,n);

void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []))
{
	int j,its,i;
	double gg,gam,fp,dgg;
	double *g,*h,*xi;

	g=vector(1,n);
	h=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,xi);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		linmin1(p,xi,n,fret,func);
		if (2.0*fabs(*fret-fp) <= ftol*(fabs(*fret)+fabs(fp)+EPS)) {
			FREEALL
			return;
		}
		fp=(*func)(p);
		(*dfunc)(p,xi);
		if(gPowellTrace)
		    {
		     printf("TRACE (Powell iteration %i)(start)\n", *iter);
		     for (i=1;i<=n;i++)
			 printf("p[%i] %f %f\n",i, p[i],xi[i]);
		     printf("	Objective function value = %f\n", (*func)(p));
		    }
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	nrerror("Too many iterations in frprmn");
}
#undef EPS
#undef FREEALL
#undef NRANSI

#define NRANSI
#define EPS 3.0e-8
#define TOLX (4*EPS)
#define STPMX 100.0
#define FMAX(a,b) ((a) >(b) ? (a):(b))

#define FREEALL free_vector(xi,1,n);free_vector(pnew,1,n); \
free_matrix(hessin,1,n,1,n);free_vector(hdg,1,n);free_vector(g,1,n); \
free_vector(dg,1,n);

int dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	double(*func)(double []), void (*dfunc)(double [], double []))
{
	int lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
		 double *f, double stpmax, int *check, double (*func)(double []));
	int check,i,its,j;
	int itmax;
	double den,fac,fad,fae,fp,stpmax,sum=0.0,sumdg,sumxi,temp,test;
	double *dg,*g,*hdg,**hessin,*pnew,*xi;

	itmax=*iter; /* my code */
	dg=vector(1,n);
	g=vector(1,n);
	hdg=vector(1,n);
	hessin=matrix(1,n,1,n);
	pnew=vector(1,n);
	xi=vector(1,n);
	fp=(*func)(p);
	(*dfunc)(p,g);
	for (i=1;i<=n;i++) {
		for (j=1;j<=n;j++) hessin[i][j]=0.0;
		hessin[i][i]=1.0;
		xi[i] = -g[i];
		sum += p[i]*p[i];
	}
	stpmax=STPMX*FMAX(sqrt(sum),(double)n);
	for (its=1;its<=itmax;its++) {
		*iter=its;
		if (lnsrch(n,p,fp,g,xi,pnew,fret,stpmax,&check,func)==0)
			return 0;
		fp = *fret;
		for (i=1;i<=n;i++) {
			xi[i]=pnew[i]-p[i];
			p[i]=pnew[i];
		}
		test=0.0;
		for (i=1;i<=n;i++) {
			temp=fabs(xi[i])/FMAX(fabs(p[i]),1.0);
			if (temp > test) test=temp;
		}
		if (test < TOLX) {
			FREEALL
			return 1;
		}
		for (i=1;i<=n;i++) dg[i]=g[i];
		(*dfunc)(p,g);
		if(gPowellTrace)
		    {
		     printf("TRACE (dfpmin iteration %i)(start)\n", *iter);
		     for (i=1;i<=n;i++)
			 printf("p[%i] %f %f\n",i, p[i],g[i]);
		     printf("	Objective function value = %f\n", (*func)(p));
		    }
		test=0.0;
		den=FMAX(*fret,1.0);
		for (i=1;i<=n;i++) {
			temp=fabs(g[i])*FMAX(fabs(p[i]),1.0)/den;
			if (temp > test) test=temp;
		}
		if (test < gtol) {
			FREEALL
			return 1;
		}
		for (i=1;i<=n;i++) dg[i]=g[i]-dg[i];
		for (i=1;i<=n;i++) {
			hdg[i]=0.0;
			for (j=1;j<=n;j++) hdg[i] += hessin[i][j]*dg[j];
		}
		fac=fae=sumdg=sumxi=0.0;
		for (i=1;i<=n;i++) {
			fac += dg[i]*xi[i];
			fae += dg[i]*hdg[i];
			sumdg += SQR(dg[i]);
			sumxi += SQR(xi[i]);
		}
		if (fac*fac > EPS*sumdg*sumxi) {
			fac=1.0/fac;
			fad=1.0/fae;
			for (i=1;i<=n;i++) dg[i]=fac*xi[i]-fad*hdg[i];
			for (i=1;i<=n;i++) {
				for (j=1;j<=n;j++) {
					hessin[i][j] += fac*xi[i]*xi[j]
					-fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
				}
			}
		}
		for (i=1;i<=n;i++) {
			xi[i]=0.0;
			for (j=1;j<=n;j++) xi[i] -= hessin[i][j]*g[j];
		}
	}
/*	nrerror("too many iterations in dfpmin");*/
	FREEALL
	return 0;
}
#undef EPS
#undef TOLX
#undef STPMX
#undef FREEALL
#undef NRANSI
#define NRANSI
#define ALF 1.0e-4
#define TOLX 1.0e-7

int lnsrch(int n, double xold[], double fold, double g[], double p[], double x[],
	double *f, double stpmax, int *check, double (*func)(double []))
{
extern NODETYPE *gRoot;
	int i;
	double a,alam,alam2,alamin,b,disc,f2,fold2,rhs1,rhs2,slope,sum,temp,
		test,tmplam;

	*check=0;
	for (sum=0.0,i=1;i<=n;i++) sum += p[i]*p[i];
	sum=sqrt(sum);
	if (sum > stpmax)
		for (i=1;i<=n;i++) p[i] *= stpmax/sum;
	for (slope=0.0,i=1;i<=n;i++)
		slope += g[i]*p[i];
	test=0.0;
	for (i=1;i<=n;i++) {
		temp=fabs(p[i])/FMAX(fabs(xold[i]),1.0);
		if (temp > test) test=temp;
	}
	alamin=TOLX/test;
	alam=1.0;
	for (;;) {
		for (i=1;i<=n;i++) x[i]=xold[i]+alam*p[i];
		*f=(*func)(x);
		if (alam < alamin) {
			for (i=1;i<=n;i++) x[i]=xold[i];
			*check=1;
			return 1;
		} else if (*f <= fold+ALF*alam*slope) return 1;
		else {
			if (alam == 1.0)
				tmplam = -slope/(2.0*(*f-fold-slope));
			else {
				rhs1 = *f-fold-alam*slope;
				rhs2=f2-fold2-alam2*slope;
				a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
				b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
				if (a == 0.0) tmplam = -slope/(2.0*b);
				else {
					disc=b*b-3.0*a*slope;
					if (disc<0.0)
						{
						 doGenericAlert("Roundoff problem in lnsrch.");
						 return 0;
						}
					else tmplam=(-b+sqrt(disc))/(3.0*a);
				}
				if (tmplam>0.5*alam)
					tmplam=0.5*alam;
			}
		}
		alam2=alam;
		f2 = *f;
		fold2=fold;
		alam=FMAX(tmplam,0.1*alam);
	}
}
#undef ALF
#undef TOLX
#undef NRANSI
                                                                                                                                                                                                                                                                                                                                            r8s/powell.h                                                                                        0000644 0000766 0000120 00000001276 07565204055 012734  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  int dfpmin(double p[], int n, double gtol, int *iter, double *fret,
	double(*func)(double []), void (*dfunc)(double [], double []));

double f1dim1(double x);
void mnbrak1(double *ax, double *bx, double *cx, double *fa, double *fb, double *fc,
	double (*func)(double));
void linmin1(double p[], double xi[], int n, double *fret, 
	double (*func)(double []));
double brent1(double ax, double bx, double cx, double (*f)(double), double tol,
	double *xmin);
int powell1(double p[], double **xi, int n, double ftol, int *iter, double *fret,
	double (*func)(double []));
void frprmn(double p[], int n, double ftol, int *iter, double *fret,
	double (*func)(double []), void (*dfunc)(double [], double []));

                                                                                                                                                                                                                                                                                                                                  r8s/relativeRates.c                                                                                 0000644 0000766 0000120 00000011003 07565204055 014224  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include "nexus.h"
#include "TreeUtils.h"
#include "relativeRates.h"
#include "WuLi.h"
#include "MyUtilities.h"
#include "root3taxa.h"
#include "distance.h"

/* Global */


/***** COMMENTS *******

Note that missing data can cause problems in the following way:  Currently distances are calculated
pairwise.  This means that the relative rates tests can use sites that do not have data in all
three taxa.  Sometimes the distance between the two ingroups can therefore be zero but the two 
distances to the outgroups might not be equal!  When this occurs the variance calculations used
in WuLi and Steel et al. get mucked up and can try to take the square root of a negative number.

Long term solution is to give user a choice between (a) including only sites where data is present
in all taxa or (b) explaining the source of these errors.


*/


/**************************************************************************/
int doRelativeRates(StrListPtr aTaxaList, NODETYPE * root)
{
	extern FILE * outstream;
	extern struct NexDataType *gNexDataPtr;	
	int id[3], ix,jx,taxonID,kk,kind;
	char c, *dummy,*taxon,*tmp, *pi,*pj,*pRow1, *taxon1,*taxon2;
	StrListPtr WholeSelList,rootedList, unrootedList;
	struct MyRecType * dptr;
	struct NexDataType * nexPtr;	/* This is THE data structure for the NEXUS data */
	int i,j,k,NList;


	unrootedList=newStrListN(3);
	NList=lengthList(aTaxaList);
	WholeSelList=newStrListN(NList); /* for some reason I need to copy
			to a new list before writing to it; I suspect a problem
			in the 'setkthnode' line below but OK for now*/
	for (ix=1;ix<=NList;ix++) /* convert any taxon ids to taxon names 
								unless already stored that way*/
		{
		taxon=getkthStr(aTaxaList,ix);
		if(isStrInteger(taxon))
			{
			taxonID=strtod(taxon,&dummy);
			setkthNode(WholeSelList, ix, getkthStr(gNexDataPtr->TaxaList,taxonID));
			}
		else
			setkthNode(WholeSelList, ix, taxon);
		}




		
	nexPtr=gNexDataPtr;	

/* Preliminaries: do some checking to see if we can proceed and open an output file */

	if (nexPtr->isChars==0)
		{
		doGenericAlert("Characters not available in NEXUS file");
		return 0;	
		}
		
	if (nexPtr->isTaxa==0) 
		{
		doGenericAlert("Taxa not available in NEXUS file");
		return 0;	
		}
#if 0		
	printf("RELATIVE RATE TESTS: Method = ");
	switch (gNexDataPtr->RateBlockParms.RRtype)
		{
		case WULI:printf("Wu & Li\n"); break;
		case MIKE:printf("Mike's method\n"); break;
		case STEEL:printf("Steel et al.\n"); break;
		case TAJIMA:printf("Tajima\n"); break;
		}
	printf("(* = P<0.05; ** = P<0.01; *** = P<0.001)\n");
	printf("(Positive z-score means higher rate in first ingroup)\n\n");
	if (gNexDataPtr->RateBlockParms.isBS)
		printf("Bootstrap estimates of variance: N replicates = %li; Seed = %li\n",
		gNexDataPtr->RateBlockParms.NReps,
		gNexDataPtr->RateBlockParms.seed);
	printf("(Outgroup        (Ingroup1          Ingroup2     ))\t\tz (\"exact\")");

	if (gNexDataPtr->RateBlockParms.isBS)
		printf("\tz (bootstrap)\t[mean:an bs] (sdev: an bs)\n");	/* shift the column heading over if showing bs results */
	else printf("\n");
	printf("---------------------------------------------------------------------------------------------------\n");
#endif
	i=1;j=2;k=3;
	
	for (i=1;i<j;i++)
		for (j=i+1;j<k;j++)
			for (k=j+1;k<=NList;k++)
				{ /* set up the unsorted 3-list */
				setkthNode(unrootedList,1,getkthStr(WholeSelList,(long)i) );
				setkthNode(unrootedList,2,getkthStr(WholeSelList,(long)j) );
				setkthNode(unrootedList,3,getkthStr(WholeSelList,(long)k) );
				rootedList=root3taxa(unrootedList,root); /* returns the properly sorted 3-list */
				/* Now convert the three taxon names in List to ID numbers and call WuLi */
				if (rootedList)  /* if null it means a polytomy or other error */
					{
					for (ix=0;ix<3;ix++)
						{
							taxon=getkthStr(rootedList,ix+1);
							jx = findMatchStr(nexPtr->TaxaList, taxon);
							if (jx ==0)
								doGenericAlert ("Matching taxon label not found in WuLi");
							id[ix]=jx;	/* make sure ids are on [1..ntaxa] */
						}
					(void)WuLiStub(id[0],id[1],id[2]);
					freeStrList(rootedList);
					}
				}




				
freeStrList(unrootedList);			
freeStrList(WholeSelList);	
return 1;
}
/**************************************************************************/
int doGroupRelRates(StrListPtr ig1List,StrListPtr ig2List,StrListPtr ogList)
{
long ig1Size,ig2Size, ogSize, i,j;

ig1Size=lengthList(ig1List);
ig2Size=lengthList(ig2List);
ogSize=lengthList(ogList);
if ((ig1Size==0) | (ig2Size==0) | (ogSize==0)) 
	{
	doGenericAlert("At least one of taxa lists is empty");
	return (0);
	}

} 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             r8s/relativeRates.h                                                                                 0000644 0000766 0000120 00000000074 07565204055 014237  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  int doRelativeRates(StrListPtr aTaxaList, NODETYPE * root);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                    r8s/root3taxa.c                                                                                     0000644 0000766 0000120 00000011641 07565204055 013346  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #if 0
#define isEqual(a,b)		(!strcmp((a),(b)))
#endif

#include "TreeUtils.h"
#include "memory.h"
#include "root3taxa.h"
#include "nexus.h"
#include "structures.h"
#include <stdlib.h>

int		gbSLcount;
StrListPtr 	gThreeList,g3SelList;
/***********************************************************************************/
/***********************************************************************************/

StrListPtr root3taxa(StrListPtr unsortedList, NODETYPE *root)	

/* 	
	If the three taxa form a polytomy, bails. Otherwise...
	Builds a list of three taxon names using my list structure.
	First and second elements are the ingroup, third is the outgroup.
	Returns the address of the list.  NULL is the error return.
	
	The next THREE functions are all necessary for this routine.
	
	Note that the algorithm for determining which taxon of the three is an outgroup is
	not intuitive.  Each taxon in order that it is encountered on a traversal is assigned 
	either 1,2,4.  When these are added up as we proceed higher in the tree, they make for
	unique pairs of sister group numbers.  Thus, if (4,(1,2)), we have 4 versus 3 and can tell
	that the single number 4 is the outgroup.  i.e., 1,2,4 are never obtained by addition of
	themselves--this is equivalent to binary coding.
	
	NOT A GOOD X-TREE ROUTINE.  IT TRAVERSES THE ENTIRE TREE EVEN IF THREE TAXA ARE CLOSELY RELATED	

	*/



{
	extern StrListPtr g3SelList;
	int total,OG,IG1,IG2;
	char* dummy;
	StrListPtr localThreeList;
	
	g3SelList=unsortedList;
/*	if (numSelected(root)==3)  */ /* no longer confine to 3 taxa in selection */
		{
		gbSLcount=1;
		gThreeList=newStrListN(3);
		localThreeList=newStrListN(3);
		if (gThreeList && localThreeList)
			total=bSLrecurse(root);  /* should work now */
		else
			{
		/*	doGenericAlert("Allocation error in String List(build list)");*/
			return NULL;			/* Avoids allocation errors */
			}
		OG=bSLrecurse2(root);
		if ((-OG >=1) && (-OG <=3))
			{
			switch (-OG)
				{
				case 1:
					IG1=2;IG2=3;OG=1; break;
				case 2:
					IG1=1;IG2=3;OG=2; break;
				case 3:
					IG1=1;IG2=2;OG=3; break;
				}
				
			setkthNode(localThreeList,1, getkthStr(gThreeList,IG1) );
			setkthNode(localThreeList,2, getkthStr(gThreeList,IG2) );
			setkthNode(localThreeList,3, getkthStr(gThreeList,OG) );			
			freeStrList(gThreeList);
			return localThreeList;
			}
		else if (OG == -2000) /* polytomy */
			return NULL;
		}
/*	else
		return NULL;*/

}

int bSLrecurse(NODETYPE *node)

/* first pass through the tree doing some stuff when it finds any of the three taxa whose
names are contained in the global list 'g3SelList'.*/
	
{
	extern struct NexDataType *gNexDataPtr;	
	char *dummy, *taxon;
	int sum=0,k,ix;
	NODETYPE *child;
	if (!node) return(-1);
	if (isTip(node) ) 
		{
		if (gNexDataPtr->isTranslate)	/* trees stored with translation table */
			{
			ix=strtod(node->taxon_name,&dummy);/* this is the taxon code */
			taxon=getkthStr(gNexDataPtr->TransList,ix);
			}
		else
			taxon=node->taxon_name;
		if 		(
		  		(isEqual(taxon,getkthStr(g3SelList,1))) ||
		  		(isEqual(taxon,getkthStr(g3SelList,2))) ||
		  		(isEqual(taxon,getkthStr(g3SelList,3)))
		  	 	)  /**** new ****/
					{
					switch (gbSLcount)
						{
						case 1: k=1;break;
						case 2:	k=2;break;
						case 3:	k=4;break;	/* these codes get assigned to the three taxa */
						}
					setkthNode(gThreeList, gbSLcount, taxon);	/*put in list at appropriate place*/
					node->numSelectedDesc=k; 
					++gbSLcount;
					return (k);
					}
	    else
			{
			node->numSelectedDesc=0; 
			return (0);
			}
		}
	else	/* interior node */
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			sum+=bSLrecurse(child);
		node->numSelectedDesc=sum;
		return (sum);
		}
}

int bSLrecurse2(NODETYPE *node)	  /* this routine should not need changing */

/* second pass through the tree making decisions about which of three taxon is
the outgroup.  Returns a code to 'root3taxa' */

{
	int sum=0,lastsum,count=0;
	NODETYPE *child;
	if (!node) return(-1);
	if (isTip(node) ) 
		return (node->numSelectedDesc);

	else	/* interior node */
		{
		child=node->firstdesc;
		SIBLOOP(child)
			{
			if (child->numSelectedDesc > 0)  /* can ignore children with no selected grandchildren*/
				{
				++count;
				lastsum=sum;		/* we need to track up to two numbers */
				sum=bSLrecurse2(child);
				if (sum<0) return (sum); /* This is how we bail out of the recursion when OG IS FOUND*/
				}
			}
		if (count >= 3)
			return (-2000); /* this a a polytomy */

		switch (count) 
			{
			case 0: return (-1000); /* this means none selected--shouldn't happen */
			case 1: break;			/* do nothing; just continue looking */
			case 2: 
				if (sum+lastsum==7)
					{
					if ((sum==1) || (lastsum==1)) return (-1);
					if ((sum==2) || (lastsum==2)) return (-2);
					if ((sum==4) || (lastsum==4)) return (-3);	/* decode the codes */
					/* abs of return value will tell us which element of list is the OG! MAGIC, huh? */
					
					}
			}
		
		return (node->numSelectedDesc);
		}
}

                                                                                               r8s/root3taxa.h                                                                                     0000644 0000766 0000120 00000000172 07565204055 013350  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  int			bSLrecurse(NODETYPE *node),
			bSLrecurse2(NODETYPE *node);

StrListPtr root3taxa(StrListPtr list, NODETYPE *root);
                                                                                                                                                                                                                                                                                                                                                                                                      r8s/sfun_.c                                                                                         0000644 0000766 0000120 00000000303 07565204055 012525  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  void sfun_(int *N,double X[],double *F, double G[])
{
int i;
double T;
*F=0;
for (i=1;i<=*N;i++)
	{

         T    = X[i-1] - i;
         *F    = *F + T*T;
         G[i-1] = 2 * T;
	}
return;
}

                                                                                                                                                                                                                                                                                                                             r8s/storeNexusFile.c                                                                                0000644 0000766 0000120 00000001760 07565204055 014402  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include <stdio.h>
#include <stdlib.h>
#include "storeNexusFile.h"
#include "MyUtilities.h"

#define	MAX_BUFFER_SIZE	20000000  /* THIS ERROR ISN"T ALWAYS CAUGHT APPARENTLY!!*/


/***************

Prompts for a file name, reads the file into a large, fixed length buffer, and
returns a pointer to that buffer.

NEED TO SLURP ANY SIZE FILE; REWRITE WITH A SMALL BUFFER THAT KEEPS ADDING
DYNAMIC SPACE TO THE BIG BUFFER

*/ 


char * storeNexusFile (FILE * inFileStream)

{
	char	*BigBuffer;
	int		c;
	long	count=0,i=0;

	
	
	BigBuffer=(char*)malloc(MAX_BUFFER_SIZE*sizeof(char));
	if (!BigBuffer)	
		{
		doGenericAlert("Could not allocate file buffer");
		return NULL;		
		}

	while ((c=fgetc(inFileStream)) != EOF)	/* have to define c as int so that EOF can be detected*/
		{
		if (i >= MAX_BUFFER_SIZE-1) /* have to save room for terminating null */
			{
			doGenericAlert("Nexus file exceeds 500k maximum");
			return NULL;
			}
		BigBuffer[i]=c;
		++i;
		
		}
		BigBuffer[i]='\0';

return BigBuffer;	
		
}
                r8s/storeNexusFile.h                                                                                0000644 0000766 0000120 00000000063 07565204055 014402  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #include <stdio.h>
char * storeNexusFile (FILE *);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                             r8s/structures.c                                                                                    0000644 0000766 0000120 00000037256 07565204055 013657  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /*  MODULE FOR STRING LIST STRUCTURES




StrListPtr 		newStrList(void);
StrListPtr 		lastStrNode(StrListPtr node);
StrListPtr 		kthStrNode(StrListPtr node, long k);
int			setkthNode(StrListPtr node, long k, char* s);
int			appendStrList(StrListPtr firstNode, char *s);
void			xprintStrList(StrListPtr aList);
long 			lengthList(StrListPtr node);
char*			getkthStr(StrListPtr node, long k);
void 			catkthStr(StrListPtr list, char* s, long i);
StrListPtr 		newStrListN(long numElements);
void 			freeStrList(StrListPtr node);
long 			findMatchStr(StrListPtr List, char * target)


*/



#include "structures.h"
#include "MyUtilities.h"
#include "memory.h"

StackPtr newStack(int maxElements)
{
StackPtr s;
s=(StackPtr)myMalloc(sizeof(struct Stack));
s->maxElements=maxElements;
s->numElements=0;
s->data=(double*)myMalloc(maxElements*sizeof(double));
return s;
}
void freeStack(StackPtr S)
{
myFree(S->data);
myFree (S);
return;
}
int hasElements(StackPtr S)
{
if (S->numElements>0)
	return 1;
else
	return 0;
}

#define DMIN(a,b)			( (a)<(b) )  ? (a):(b)
void pushD(StackPtr S, double x)
{
int i,top;
double *s;
s=S->data;
top=DMIN(S->numElements,S->maxElements-1);
for (i=top;i>=1;i--)
	s[i]=s[i-1];
s[0]=x;
if (S->numElements < S->maxElements)
	++(S->numElements);
return;
}
double popD(StackPtr S)
{
int i,top;
double *s,x;
s=S->data;
x=s[0];
top=S->numElements;
for (i=0;i<=top-2;i++)
	s[i]=s[i+1];
if (S->numElements > 0)
	--(S->numElements);
return x;
}

/***************************************************/
StrListPtr 
string_list_intersect(StrListPtr s1, StrListPtr s2)


{
long L1,L2;
StrListPtr snew;
snew = NULL;
L1=lengthList(s1);
L2=lengthList(s2);
while(s1)
	{
	if (findMatchStr(s2,s1->s))
		{
		if (!snew)
			snew = newStrList();
		appendStrList(snew, s1->s);
		}
	s1=s1->next;
	}
return snew;
} 

/***************************************************/
int 
string_lists_same(StrListPtr s1, StrListPtr s2)

/* returns 1 if the string lists are the same (i.e., have the same
elements in any order), otherwise returns 0 */

{
long L1,L2;
L1=lengthList(s1);
L2=lengthList(s2);
if (L1 != L2)
	return 0;
while(s1)
	{
	if (!findMatchStr(s2,s1->s)) return 0;
	s1=s1->next;
	}
return 1;
} 

/***************************************************/

/* finds the first occurrence of string target in List, and returns the index for that 
string (on 1..n).  Returns 0 if not found. */

long 
findMatchStr(StrListPtr List, char * target)
{
long ix=1;
while (List)
	{
	if (isEqual(target,List->s)) return (ix);
	++ix;
	List = List ->next;
	}
	
return 0;
}

/***************************************************/
void glomStrLists(StrListPtr A,  StrListPtr B)
{
/* (destructively merges two string lists.  Finds the union of the list and
MAKES THE FIRST LIST (A) THIS UNION.  If you want a NEW list write something*/

long lengthB, ixB;
char * Belement;

lengthB = lengthList(B); 
for (ixB=1;ixB<=lengthB; ixB++)
    {
    Belement=getkthStr(B, ixB);
    if (!findMatchStr(A, Belement))
	appendStrList(A, Belement);
    }  
    
}
/***************************************************/
		
/* Prints out the contents of a string list, one element per row */

void 			
xprintStrList(StrListPtr node)
	{
	long count=0;
	
	
	while ((node != NULL)  )
		{
		if(node->s != NULL)
			printf("%s\n",node->s);
		node=node->next;
		}
		
	return;
	}

/***************************************************/
		
/* makes head node of a linked list and initializes 
string pointer and next pointer to NULL

--returns NULL on error*/

StrListPtr 		
newStrList(void)

	{
	
	StrListPtr node;
	node = (StrListPtr)myMalloc(sizeof(struct StrList));
	if (node != NULL)
		{
		node->next=NULL;
		node->s=NULL;
		}
	return node;
	}
/***************************************************/
		
/* Makes linked list of numElements elements and initializes them all to NULL

--returns NULL on error*/

StrListPtr 		
newStrListN(long numElements)

	{
	
	StrListPtr node;
	long k;
	
	if (numElements<=0 )
		return NULL;
	node = newStrList();
	for (k=1;k<=numElements-1;k++)
		appendStrList(node,NULL);
	
	return node;
	}
/***************************************************/		

/* returns last node in linked list or NULL if error */

StrListPtr 		
lastStrNode(StrListPtr node)

	{
	long count=0;
	if (node != NULL)
		while (node->next != NULL)
			{
			++count;
			node=node->next;
			}
	return node;
	}
/***************************************************/		

/* returns the kth node, or NULL if past end or bad k; k is on [1..size of list] */

StrListPtr 		
kthStrNode(StrListPtr node, long k)

{
	if ( k<0) return NULL; /* Error */
	while (  (node->next != NULL)  && (--k > 0)  )
		{
		node=node->next;	
		}
	return node;
}
/***************************************************/
		
/* Sets the kth element in the list to the value of string s;
	--returns 1 if OK, 0 if error*/

int				
setkthNode(StrListPtr list, long k, char* s)
{
	StrListPtr node;
	node = kthStrNode(list,k);
	if (node != NULL)
		{
		/** major modification in next line to allow clean overwrite of string**/
		if (node->s)
			myFree(node->s);
		node->s=DupStr(s);	/* make a persistent version of this string */
		if (node->s != NULL)
			return 1;
		else
			return 0; /* error in DupStr */
		}
	else
		return 0; /* Error */

}
/***************************************************/
		
/* Returns a ptr to the kth string element in the list, or NULL on error 

NB!  If you are going to stash this string somewhere, use DupStr to make a persistent
copy of the string itself.  Otherwise, some other routine may free the original locatiions
that this pointer points to.  */



char*			
getkthStr(StrListPtr node, long k)
{
	char* s;
	if ( k<0) return NULL; /* Error */
	while (  (node->next != NULL)  && (--k > 0)  )
		{
		node=node->next;	
		};
	if (node != NULL)
		{
		s=node->s;
		return /* DupStr */(s);	/* Return a persistent version of this or else all hell breaks
								loose when we dispose of the string list */
		}
	else
		return NULL; /* Error */

}
/***************************************************/		

/* Adds a new item to the list, but...
	If the last node in the list has a null string, it puts the item
	in that last item rather than making a new node.  This allows
	one to append to a newly created empty list.  
	A duplicate of the string is made before storing.
	If the parameter 'str' is NULL, just adds an empty node with a NULL string;
	Returns 1 if OK, 0 if error */



int		
appendStrList(StrListPtr firstNode, char *str)
{

	StrListPtr	lastNode, newNode;
	char* 		cpStr;
	
	if (firstNode == NULL ) return 0; 	/* error  */
	
	lastNode=lastStrNode(firstNode);
	if (lastNode != NULL)
		{
		if (str == NULL)	/* add a null new node */
			{
			newNode=newStrList();
			lastNode->next=newNode;
			return 1;		
			}
		else				/* add a new node that has a string */
			{
			cpStr=DupStr(str);/* make a persistent copy of the string and store */
			if (cpStr == NULL)
				return 0;	/* error */
			if (lastNode->s == NULL)  /* there's no item stored at this node, so put it here */
				{
				lastNode->s=cpStr;
				return 1;
				}
			else				/* there IS an item at this node, make a new node and store*/
				{
				newNode=newStrList();
				lastNode->next=newNode;
				newNode->s=cpStr;
				return 1;
				}
			}
		}
	else
		return 0; 	/* error */
	}
/***************************************************/		
/* returns number of elements of list */

long 		lengthList(StrListPtr node)

{
	long length=1;
	if (node == NULL) return 0; /* Error */
	while (  node->next != NULL )
		{
		node=node->next;
		++length;	
		}
	return length;
}
/***************************************************/		
/* Finds the kth string in the list and concatenates string 'ss' to it.
If that string is NULL, it just puts 'ss' there. */

void catkthStr(StrListPtr list, char* ss, long i)
	{
	StrListPtr 		ithnode;
	ithnode=kthStrNode(list,i);
	if (ithnode->s == NULL) 
		ithnode->s=DupStr(ss);		/* destination string is null, so just put string here */
	else
		concat(&(ithnode->s),ss);	/* destination exists so concatenate */
	return;
	}
/***************************************************/		
void freeStrList(StrListPtr node)
{
if (node == NULL) return;
if (node->next)
	freeStrList(node->next);
myFree(node->s);
myFree (node);
return;


}

/***************************************************/
/***************************************************/
/* GENERIC POINTER LISTS */

/* NB!  Sizeof refers to the size of the object being pointed to
    by the pointers in the list. */

/* NB!  A list is defined if it has been created by pNewList
 *	A list is defined but empty if it has no items,  which
	is indicated by having a NULL pointer as its .item
 */


/***************************************************/
/***************************************************/
		
/* makes head node of a linked list of void pointers 

--returns NULL on error*/

PtrList pNewList(void)

/*  USE THIS ONE FROM NOW ON! 
    Returns a pointer to the headnode of a generic pointer list.
 * Sets the first item to NULL,  making it an "empty" list
 */

{
	
	PtrList node;
	node = (PtrList)myMalloc(sizeof(struct PtrListStruct));
	node->next=NULL;
	node->item=NULL;
	return node;
}

PtrList pNewListAlt(size_t size) /* OLD! probably should NOT initialize item here! 
	(it makes adding items later tricky. FIX some day.  Used to initialize
	tree lists and the like,  and currently sets up their first nodes) */
{
	
	PtrList node;
	node = (PtrList)myMalloc(sizeof(struct PtrListStruct));
	node->next=NULL;
	node->item=(void *)myMalloc(size);
	
return node;
}
/***************************************************/
long pLengthList(PtrList node)

/** ! seems like this will return a length of 1 even if the first node is empty*/
{
	long length=1;
	if (node == NULL) fatal("list is empty"); /* Error */
	while (  node->next)
		{
		node=node->next;
		++length;	
		}
	return length;
}
PtrList pListLastNode(PtrList node)

	{
	long count=0;
	if (node != NULL)
		while (node->next != NULL)
			{
			++count;
			node=node->next;
			}
	return node;
	}
PtrList pListgetkthNode(PtrList node, long k)

{
	if ( k<0 ) return NULL; /* Error */
	while (  (node->next != NULL)  && (--k > 0)  )
		{
		node=node->next;	
		}
	return node;
}


PtrList pListAddNode(PtrList firstNode, size_t size)

/*** Appends a node of the given size on the list and returns a pointer
	to that node. Item is set to NULL. If firstNode is NULL, 
	create a new EMPTY list ***/

{

	PtrList	lastNode, newNode;
	
	if (firstNode == NULL ) 
		return pNewListAlt(size);
	
	lastNode=pListLastNode(firstNode);
	if (lastNode != NULL)
		{
			newNode=pNewListAlt(size);
			if (newNode)
				{
				lastNode->next=newNode;
				return newNode;	
				}
			else
				fatal("Couldn't allocate newnode in pListAddNode");	
		}
	else
		{
		fatal("Couldn't get lastNode right in pListAddNode");
		return NULL; 	/* error */
		    
		}
	}
/***************************************************/		
void pListAddItem2(PtrList firstNode, size_t size, void * ptrItem)

/** CLUNKY CODE DO NOT USE **/
/*** given a pointer to an item of size 'size', this adds this item
    to the list.  Does this either at the current last node if that has no item
    or at a new next node ***/

	{
	PtrList	newNode,  lastNode;
	if (firstNode==NULL)
	    return; /* error, no list */
	lastNode=pListLastNode(firstNode);
	if (lastNode->item == NULL)
		{
			lastNode->item=ptrItem;
		}
	else
		{
		newNode=pListAddNode(firstNode, size);
		newNode->item=ptrItem;
		}
	return;
	}
/***************************************************/		
void pListAddItem(PtrList firstNode, void * ptrItem)

/*** given a pointer to an item of size 'size', this adds this item
    to the list.  Does this either at the current last node if that has no item
    or at a new next node ***/

	{
	PtrList	newNode,  lastNode;
	if (firstNode==NULL)
	    return; /* error, no list */
	lastNode=pListLastNode(firstNode);
	if (lastNode->item == NULL)
		{
			lastNode->item=ptrItem;
		}
	else
		{
		newNode=pNewList();
		lastNode->next=newNode;
		newNode->item=ptrItem;
		}
	return;
	}
/***************************************************/		
/***************************************************/		
void DfreepList(PtrList node)

/* EXTREME free of generic list */

{
if (node == NULL) return;
if (node->next)
	freepList(node->next);
myFree(node->item);  /** NB! This may be inadequate! Because item may itself
			have pointers to dynamically created objects. OR
			IT MAY BE DESTRUCTIVE,  say,  if we just want to maintain
			a list of pointers to things.  This now ALSO DELETES
			the things pointed TOO! */
myFree (node);
return;


}
void freepList(PtrList node)

/* more GENTLE free of generic list (doesn't destroy the elements of the list,
 * just the pointers to the elements 
 */

{
if (node == NULL) return;
if (node->next)
	freepList(node->next);
myFree (node);
return;
}
/**********************************/
/**********************************/
/* Kludgy utils for set-vectors */
/* items in a set are numbered from 1..size, where size is the max possible
 * number of items.
 * NB! It is assumed that any operations on two sets forces both to have
 * the same 'size'.  Results are undefined otherwise.
 */


void test_set(Set a,  Set b)
{
    printf("\n");
    printf("set 1:       ");print_set(a);
    printf("set 2:       ");print_set(b);
    printf("intersection:");print_set(intersect_set(a, b));
    printf("union:       ");print_set(union_set(a, b));
    if (sets_Equal(a, b))
	printf("Sets EQUAL\n");
    if (is_subset(a, b))
	printf("1 is a subset of 2\n");
    if (is_superset(a, b))
	printf("1 is a superset of 2\n");
    printf("Percent overlap=%f\n", set_overlap(a, b));
    return;
	

}
double set_overlap(Set a, Set b)
    {
	int csize, dsize;
	Set c, d;
	c=intersect_set(a, b);
	d=union_set(a, b);
		
	csize=sizeof_set(c);
	dsize=sizeof_set(d);
	if (dsize==0)
	    return 0.0;
	else
	    return (double)csize/dsize;
	
    }

int is_subset(Set a,  Set b)
    { /* is a a subset of b? can also be equal to b! */
    if (sets_Equal(a, intersect_set(a, b)))
	return 1;
    else 
	return 0;
    }
int is_superset(Set a,  Set b)
    { /* is a superset of b? can also be equal to b! */
    if (sets_Equal(b, intersect_set(a, b)))
	return 1;
    else 
	return 0;
    }

Set newSet(int size)
    {
	Set theSet;
	int j;
	theSet = (Set)myMalloc(sizeof(struct SetStruct));
	theSet->size=size;
	theSet->element=(int *)myMalloc(size*sizeof(int));
	for (j=1;j<=size;j++)
			 remove_from_set(theSet,  j);
	return theSet;
    }
void add_to_set(Set theSet,  int item_id)
    {
    (theSet->element)[item_id-1]=1;
    return;
   }
void remove_from_set(Set theSet,  int item_id)
    {
    (theSet->element)[item_id-1]=0;
    return;
   }

int is_set_member(Set theSet,  int item_id)
    {
	if ((theSet->element)[item_id-1]==1)
	    return 1;
	else 
	    return 0;
    }

void print_set(Set theSet)
    {
    int i;
    for (i=1;i<=theSet->size;i++)
    if (is_set_member(theSet, i))
	printf("1");
    else
	printf("0");
    printf("\n");
    return;
    }

Set union_set(Set a,  Set b)
{
    Set c;
    int item, csize;
    c=newSet(a->size);
    for (item=1;item<=a->size;item++)
	if ((is_set_member(a, item)) || (is_set_member(b, item)))
	    add_to_set(c, item);
    return c;
}
Set intersect_set(Set a,  Set b)
{
    Set c;
    int item;
    c=newSet(a->size);
    for (item=1;item<=a->size;item++)
	if ((is_set_member(a, item)) && (is_set_member(b, item)))
	    add_to_set(c, item);
    return c;
}

int sizeof_set(Set a)
{
int n=0, item;
for (item=1;item<=a->size;item++)
    if (is_set_member(a, item))
	++n;
return n;    
}
int is_empty_set(Set a)
{
int n=0, item;
for (item=1;item<=a->size;item++)
    if (is_set_member(a, item))
	return 0;
return 1;    
}

int sets_Equal( Set a,   Set b)
/* note that two empty sets are considered equal */
/* returns 1 if two sets are equal, zero otherwise */

{
int i;
for (i=1;i<=a->size;i++)
    if (is_set_member(a, i) != is_set_member(b, i))
	return 0;   
return 1;   
}

/**********************************/
                                                                                                                                                                                                                                                                                                                                                  r8s/structures.h                                                                                    0000644 0000766 0000120 00000005112 10357560240 013640  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  #ifndef _STRUCTURES
#define _STRUCTURES
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define MAX_LIST_SIZE	1000

#define STU	0	/* STRING TO UPPER? */
#define isEqual(a,b)		(!strcasecmp((a),(b))) 
#define LISTLOOP(c)		for (; (c); (c)=(c)->next)


/* My list structure */

struct Stack {
				int maxElements;
				int numElements;
				double *data;
		};

struct StrList {
				char* s;
				struct StrList* next;
				};
struct PtrListStruct {
				void * item;
				struct PtrListStruct* next;
				};

struct SetStruct {
				int size;	/* max possible elements */
				int *element;	/* pointer to an int array of
						0's and 1's */
		   };

typedef struct Stack		* StackPtr;
typedef struct PtrListStruct	* PtrList;			
typedef struct SetStruct	* Set;
typedef struct StrList		* StrListPtr;

/*************************************/

StackPtr		newStack(int maxElements);
void 			freeStack(StackPtr S);
void 			pushD(StackPtr S, double x);
double			popD(StackPtr S);
int 			hasElements(StackPtr S);

StrListPtr 		string_list_intersect(StrListPtr s1, StrListPtr s2);
int			is_subset(Set a,  Set b);
int			is_superset(Set a,  Set b);
void			test_set(Set a,  Set b);
int			is_empty_set(Set a);
int			sizeof_set(Set a);
Set			intersect_set(Set a,  Set b);
Set			union_set(Set a,  Set b);
void			print_set(Set theSet);
int			is_set_member(Set theSet,  int item_id);
void			add_to_set(Set theSet,  int item_id);
void			remove_from_set(Set theSet,  int item_id);
Set			newSet(int size);
void			glomStrLists(StrListPtr A,  StrListPtr B);
StrListPtr 		newStrList(void);
StrListPtr 		lastStrNode(StrListPtr node);
StrListPtr 		kthStrNode(StrListPtr node, long k);
int			setkthNode(StrListPtr node, long k, char* s);
int			appendStrList(StrListPtr firstNode, char *s);
void			xprintStrList(StrListPtr aList);
long 			lengthList(StrListPtr node);
char*			getkthStr(StrListPtr node, long k);
void 			catkthStr(StrListPtr list, char* s, long i);
StrListPtr 		newStrListN(long numElements);
void 			freeStrList(StrListPtr node);
long			findMatchStr(StrListPtr List, char * target);
int			string_lists_same(StrListPtr s1, StrListPtr s2);

void			pListAddItem(PtrList firstNode, void * ptrItem);
void			pListAddItem2(PtrList firstNode, size_t size, void * ptrItem);
void			DfreepList(PtrList node);
void			freepList(PtrList node);
PtrList			pNewListAlt(size_t size);
PtrList			pNewList(void);
PtrList			pListLastNode(PtrList node);
PtrList			pListgetkthNode(PtrList node, long k);
PtrList			pListAddNode(PtrList firstNode, size_t size);
long 			pLengthList(PtrList node);
double			set_overlap(Set a, Set b);
int sets_Equal( Set a,   Set b);
#endif

                                                                                                                                                                                                                                                                                                                                                                                                                                                      r8s/tn.f                                                                                            0000644 0000766 0000120 00000154525 10170322670 012045  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  C%% TRUNCATED-NEWTON METHOD:  SUBROUTINES
C   FOR OTHER MACHINES, MODIFY ROUTINE MCHPR1 (MACHINE EPSILON)
c   WRITTEN BY:  STEPHEN G. NASH
C                OPERATIONS RESEARCH AND APPLIED STATISTICS DEPT.
C                GEORGE MASON UNIVERSITY
C                FAIRFAX, VA 22030
C******************************************************************
      SUBROUTINE TN (IERROR, N, X, F, G, W, LW, SFUN)
      IMPLICIT          DOUBLE PRECISION (A-H,O-Z)
      INTEGER           IERROR, N, LW
      DOUBLE PRECISION  X(N), G(N), F, W(LW)
C
C THIS ROUTINE SOLVES THE OPTIMIZATION PROBLEM
C
C            MINIMIZE F(X)
C               X
C
C WHERE X IS A VECTOR OF N REAL VARIABLES.  THE METHOD USED IS
C A TRUNCATED-NEWTON ALGORITHM (SEE "NEWTON-TYPE MINIMIZATION VIA
C THE LANCZOS METHOD" BY S.G. NASH (SIAM J. NUMER. ANAL. 21 (1984),
C PP. 770-778).  THIS ALGORITHM FINDS A LOCAL MINIMUM OF F(X).  IT DOES
C NOT ASSUME THAT THE FUNCTION F IS CONVEX (AND SO CANNOT GUARANTEE A
C GLOBAL SOLUTION), BUT DOES ASSUME THAT THE FUNCTION IS BOUNDED BELOW.
C IT CAN SOLVE PROBLEMS HAVING ANY NUMBER OF VARIABLES, BUT IT IS
C ESPECIALLY USEFUL WHEN THE NUMBER OF VARIABLES (N) IS LARGE.
C
C SUBROUTINE PARAMETERS:
C
C IERROR - (INTEGER) ERROR CODE
C          ( 0 => NORMAL RETURN)
C          ( 2 => MORE THAN MAXFUN EVALUATIONS)
C          ( 3 => LINE SEARCH FAILED TO FIND
C          (          LOWER POINT (MAY NOT BE SERIOUS)
C          (-1 => ERROR IN INPUT PARAMETERS)
C N      - (INTEGER) NUMBER OF VARIABLES
C X      - (REAL*8) VECTOR OF LENGTH AT LEAST N; ON INPUT, AN INITIAL
C          ESTIMATE OF THE SOLUTION; ON OUTPUT, THE COMPUTED SOLUTION.
C G      - (REAL*8) VECTOR OF LENGTH AT LEAST N; ON OUTPUT, THE FINAL
C          VALUE OF THE GRADIENT
C F      - (REAL*8) ON INPUT, A ROUGH ESTIMATE OF THE VALUE OF THE
C          OBJECTIVE FUNCTION AT THE SOLUTION; ON OUTPUT, THE VALUE
C          OF THE OBJECTIVE FUNCTION AT THE SOLUTION
C W      - (REAL*8) WORK VECTOR OF LENGTH AT LEAST 14*N
C LW     - (INTEGER) THE DECLARED DIMENSION OF W
C SFUN   - A USER-SPECIFIED SUBROUTINE THAT COMPUTES THE FUNCTION
C          AND GRADIENT OF THE OBJECTIVE FUNCTION.  IT MUST HAVE
C          THE CALLING SEQUENCE
C             SUBROUTINE SFUN (N, X, F, G)
C             INTEGER           N
C             DOUBLE PRECISION  X(N), G(N), F
C
C THIS IS AN EASY-TO-USE DRIVER FOR THE MAIN OPTIMIZATION ROUTINE
C LMQN.  MORE EXPERIENCED USERS WHO WISH TO CUSTOMIZE PERFORMANCE
C OF THIS ALGORITHM SHOULD CALL LMQN DIRECTLY.
C
C----------------------------------------------------------------------
C THIS ROUTINE SETS UP ALL THE PARAMETERS FOR THE TRUNCATED-NEWTON
C ALGORITHM.  THE PARAMETERS ARE:
C
C ETA    - SEVERITY OF THE LINESEARCH
C MAXFUN - MAXIMUM ALLOWABLE NUMBER OF FUNCTION EVALUATIONS
C XTOL   - DESIRED ACCURACY FOR THE SOLUTION X*
C STEPMX - MAXIMUM ALLOWABLE STEP IN THE LINESEARCH
C ACCRCY - ACCURACY OF COMPUTED FUNCTION VALUES
C MSGLVL - DETERMINES QUANTITY OF PRINTED OUTPUT
C          0 = NONE, 1 = ONE LINE PER MAJOR ITERATION.
C MAXIT  - MAXIMUM NUMBER OF INNER ITERATIONS PER STEP
C
      DOUBLE PRECISION ETA, ACCRCY, XTOL, STEPMX, DSQRT, MCHPR1
      EXTERNAL         SFUN
C
C SET UP PARAMETERS FOR THE OPTIMIZATION ROUTINE
C
      MAXIT = N/2
      IF (MAXIT .GT. 50) MAXIT = 50
      IF (MAXIT .LE. 0) MAXIT = 1
      MSGLVL = 0
      MAXFUN = 150*N
      ETA = .25D0
      STEPMX = 1.D1
      ACCRCY = 1.D2*MCHPR1()
      XTOL = DSQRT(ACCRCY)
C
C MINIMIZE THE FUNCTION
C
      CALL LMQN (IERROR, N, X, F, G, W, LW, SFUN,
     *     MSGLVL, MAXIT, MAXFUN, ETA, STEPMX, ACCRCY, XTOL)
C
C PRINT THE RESULTS
C
      IF (IERROR .NE. 0) WRITE(*,800) IERROR
      WRITE(*,810) F
      IF (MSGLVL .LT. 1) RETURN
      WRITE(*,820)
      NMAX = 10
      IF (N .LT. NMAX) NMAX = N
      WRITE(*,830) (I,X(I),I=1,NMAX)
      RETURN
800   FORMAT(//,' ERROR CODE =', I3)
810   FORMAT(//,' OPTIMAL FUNCTION VALUE = ', 1PD22.15)
820   FORMAT(10X, 'CURRENT SOLUTION IS (AT MOST 10 COMPONENTS)', /,
     *       14X, 'I', 11X, 'X(I)')
830   FORMAT(10X, I5, 2X, 1PD22.15)
      END
C
C
      SUBROUTINE TNBC (IERROR, N, X, F, G, W, LW, SFUN, LOW, UP, IPIVOT)
      IMPLICIT          DOUBLE PRECISION (A-H,O-Z)
      INTEGER           IERROR, N, LW, IPIVOT(N)
      DOUBLE PRECISION  X(N), G(N), F, W(LW), LOW(N), UP(N)
C
C THIS ROUTINE SOLVES THE OPTIMIZATION PROBLEM
C
C   MINIMIZE     F(X)
C      X
C   SUBJECT TO   LOW <= X <= UP
C
C WHERE X IS A VECTOR OF N REAL VARIABLES.  THE METHOD USED IS
C A TRUNCATED-NEWTON ALGORITHM (SEE "NEWTON-TYPE MINIMIZATION VIA
C THE LANCZOS ALGORITHM" BY S.G. NASH (TECHNICAL REPORT 378, MATH.
C THE LANCZOS METHOD" BY S.G. NASH (SIAM J. NUMER. ANAL. 21 (1984),
C PP. 770-778).  THIS ALGORITHM FINDS A LOCAL MINIMUM OF F(X).  IT DOES
C NOT ASSUME THAT THE FUNCTION F IS CONVEX (AND SO CANNOT GUARANTEE A
C GLOBAL SOLUTION), BUT DOES ASSUME THAT THE FUNCTION IS BOUNDED BELOW.
C IT CAN SOLVE PROBLEMS HAVING ANY NUMBER OF VARIABLES, BUT IT IS
C ESPECIALLY USEFUL WHEN THE NUMBER OF VARIABLES (N) IS LARGE.
C
C SUBROUTINE PARAMETERS:
C
C IERROR  - (INTEGER) ERROR CODE
C           ( 0 => NORMAL RETURN
C           ( 2 => MORE THAN MAXFUN EVALUATIONS
C           ( 3 => LINE SEARCH FAILED TO FIND LOWER
C           (          POINT (MAY NOT BE SERIOUS)
C           (-1 => ERROR IN INPUT PARAMETERS
C N       - (INTEGER) NUMBER OF VARIABLES
C X       - (REAL*8) VECTOR OF LENGTH AT LEAST N; ON INPUT, AN INITIAL
C           ESTIMATE OF THE SOLUTION; ON OUTPUT, THE COMPUTED SOLUTION.
C G       - (REAL*8) VECTOR OF LENGTH AT LEAST N; ON OUTPUT, THE FINAL
C           VALUE OF THE GRADIENT
C F       - (REAL*8) ON INPUT, A ROUGH ESTIMATE OF THE VALUE OF THE
C           OBJECTIVE FUNCTION AT THE SOLUTION; ON OUTPUT, THE VALUE
C           OF THE OBJECTIVE FUNCTION AT THE SOLUTION
C W       - (REAL*8) WORK VECTOR OF LENGTH AT LEAST 14*N
C LW      - (INTEGER) THE DECLARED DIMENSION OF W
C SFUN    - A USER-SPECIFIED SUBROUTINE THAT COMPUTES THE FUNCTION
C           AND GRADIENT OF THE OBJECTIVE FUNCTION.  IT MUST HAVE
C           THE CALLING SEQUENCE
C             SUBROUTINE SFUN (N, X, F, G)
C             INTEGER           N
C             DOUBLE PRECISION  X(N), G(N), F
C LOW, UP - (REAL*8) VECTORS OF LENGTH AT LEAST N CONTAINING
C           THE LOWER AND UPPER BOUNDS ON THE VARIABLES.  IF
C           THERE ARE NO BOUNDS ON A PARTICULAR VARIABLE, SET
C           THE BOUNDS TO -1.D38 AND 1.D38, RESPECTIVELY.
C IPIVOT  - (INTEGER) WORK VECTOR OF LENGTH AT LEAST N, USED
C           TO RECORD WHICH VARIABLES ARE AT THEIR BOUNDS.
C
C THIS IS AN EASY-TO-USE DRIVER FOR THE MAIN OPTIMIZATION ROUTINE
C LMQNBC.  MORE EXPERIENCED USERS WHO WISH TO CUSTOMIZE PERFORMANCE
C OF THIS ALGORITHM SHOULD CALL LMQBC DIRECTLY.
C
C----------------------------------------------------------------------
C THIS ROUTINE SETS UP ALL THE PARAMETERS FOR THE TRUNCATED-NEWTON
C ALGORITHM.  THE PARAMETERS ARE:
C
C ETA    - SEVERITY OF THE LINESEARCH
C MAXFUN - MAXIMUM ALLOWABLE NUMBER OF FUNCTION EVALUATIONS
C XTOL   - DESIRED ACCURACY FOR THE SOLUTION X*
C STEPMX - MAXIMUM ALLOWABLE STEP IN THE LINESEARCH
C ACCRCY - ACCURACY OF COMPUTED FUNCTION VALUES
C MSGLVL - CONTROLS QUANTITY OF PRINTED OUTPUT
C          0 = NONE, 1 = ONE LINE PER MAJOR ITERATION.
C MAXIT  - MAXIMUM NUMBER OF INNER ITERATIONS PER STEP
C
      DOUBLE PRECISION  ETA, ACCRCY, XTOL, STEPMX, DSQRT, MCHPR1
      EXTERNAL          SFUN
C
C SET PARAMETERS FOR THE OPTIMIZATION ROUTINE
C
      MAXIT = N/2
      IF (MAXIT .GT. 50) MAXIT = 50
      IF (MAXIT .LE. 0) MAXIT = 1
      MSGLVL = 0 
      MAXFUN = 150*N
      ETA = .25D0
      STEPMX = 1.D1
      ACCRCY = 1.D2*MCHPR1()



      XTOL = DSQRT(ACCRCY)
C
C MINIMIZE FUNCTION
C
      CALL LMQNBC (IERROR, N, X, F, G, W, LW, SFUN, LOW, UP, IPIVOT,
     *            MSGLVL, MAXIT, MAXFUN, ETA, STEPMX, ACCRCY, XTOL)
C
C PRINT RESULTS
C
      IF (MSGLVL .LT. 1) RETURN
C **MJS** I moved the next 2 lines of code from just before the previous test to where it is now
      IF (IERROR .NE. 0) WRITE(*,800) IERROR
      WRITE(*,810) F
      WRITE(*,820)
      NMAX = 10
      IF (N .LT. NMAX) NMAX = N
      WRITE(*,830) (I,X(I),I=1,NMAX)
      RETURN
800   FORMAT(//,' ERROR CODE =', I3)
810   FORMAT(//,' OPTIMAL FUNCTION VALUE = ', 1PD22.15)
820   FORMAT(10X, 'CURRENT SOLUTION IS (AT MOST 10 COMPONENTS)', /,
     *       14X, 'I', 11X, 'X(I)')
830   FORMAT(10X, I5, 2X, 1PD22.15)
      END
C
C
      SUBROUTINE LMQN (IFAIL, N, X, F, G, W, LW, SFUN,
     *            MSGLVL, MAXIT, MAXFUN, ETA, STEPMX, ACCRCY, XTOL)
      IMPLICIT          DOUBLE PRECISION (A-H,O-Z)
      INTEGER           MSGLVL, N, MAXFUN, IFAIL, LW
      DOUBLE PRECISION  X(N), G(N), W(LW), ETA, XTOL, STEPMX, F, ACCRCY
C
C THIS ROUTINE IS A TRUNCATED-NEWTON METHOD.
C THE TRUNCATED-NEWTON METHOD IS PRECONDITIONED BY A LIMITED-MEMORY
C QUASI-NEWTON METHOD (THIS PRECONDITIONING STRATEGY IS DEVELOPED
C IN THIS ROUTINE) WITH A FURTHER DIAGONAL SCALING (SEE ROUTINE NDIA3).
C FOR FURTHER DETAILS ON THE PARAMETERS, SEE ROUTINE TN.
C
      INTEGER I, ICYCLE, IOLDG, IPK, IYK, LOLDG, LPK, LSR,
     *     LWTEST, LYK, LYR, NFTOTL, NITER, NM1, NUMF, NWHY
      DOUBLE PRECISION ABSTOL, ALPHA, DIFNEW, DIFOLD, EPSMCH,
     *     EPSRED, FKEEP, FM, FNEW, FOLD, FSTOP, FTEST, GNORM, GSK,
     *     GTG, GTPNEW, OLDF, OLDGTP, ONE, PE, PEPS, PNORM, RELTOL,
     *     RTEPS, RTLEPS, RTOL, RTOLSQ, SMALL, SPE, TINY,
     *     TNYTOL, TOLEPS, XNORM, YKSK, YRSR, ZERO
      LOGICAL LRESET, UPD1
C
C THE FOLLOWING IMSL AND STANDARD FUNCTIONS ARE USED
C
      DOUBLE PRECISION DABS, DDOT, DSQRT, STEP1, DNRM2
      EXTERNAL SFUN
      COMMON /SUBSCR/ LGV,LZ1,LZK,LV,LSK,LYK,LDIAGB,LSR,LYR,
     *     LOLDG,LHG,LHYK,LPK,LEMAT,LWTEST
C
C INITIALIZE PARAMETERS AND CONSTANTS
C
      IF (MSGLVL .GE. -2) WRITE(*,800)
      CALL SETPAR(N)
      UPD1 = .TRUE.
      IRESET = 0
      NFEVAL = 0
      NMODIF = 0
      NLINCG = 0
      FSTOP = F
      ZERO = 0.D0
      ONE = 1.D0
      NM1 = N - 1
C
C WITHIN THIS ROUTINE THE ARRAY W(LOLDG) IS SHARED BY W(LHYR)
C
      LHYR = LOLDG
C
C CHECK PARAMETERS AND SET CONSTANTS
C
      CALL CHKUCP(LWTEST,MAXFUN,NWHY,N,ALPHA,EPSMCH,
     *     ETA,PEPS,RTEPS,RTOL,RTOLSQ,STEPMX,FTEST,
     *     XTOL,XNORM,X,LW,SMALL,TINY,ACCRCY)
      IF (NWHY .LT. 0) GO TO 120
      CALL SETUCR(SMALL,NFTOTL,NITER,N,F,FNEW,
     *     FM,GTG,OLDF,SFUN,G,X)
      FOLD = FNEW
      IF (MSGLVL .GE. 1) WRITE(*,810) NITER,NFTOTL,NLINCG,FNEW,GTG
C
C CHECK FOR SMALL GRADIENT AT THE STARTING POINT.
C
      FTEST = ONE + DABS(FNEW)
      IF (GTG .LT. 1.D-4*EPSMCH*FTEST*FTEST) GO TO 90
C
C SET INITIAL VALUES TO OTHER PARAMETERS
C
      ICYCLE = NM1
      TOLEPS = RTOL + RTEPS
      RTLEPS = RTOLSQ + EPSMCH
      GNORM  = DSQRT(GTG)
      DIFNEW = ZERO
      EPSRED = 5.0D-2
      FKEEP  = FNEW
C
C SET THE DIAGONAL OF THE APPROXIMATE HESSIAN TO UNITY.
C
      IDIAGB = LDIAGB
      DO 10 I = 1,N
         W(IDIAGB) = ONE
         IDIAGB = IDIAGB + 1
10    CONTINUE
C
C ..................START OF MAIN ITERATIVE LOOP..........
C
C COMPUTE THE NEW SEARCH DIRECTION
C
      MODET = MSGLVL - 3
      CALL MODLNP(MODET,W(LPK),W(LGV),W(LZ1),W(LV),
     *     W(LDIAGB),W(LEMAT),X,G,W(LZK),
     *     N,W,LW,NITER,MAXIT,NFEVAL,NMODIF,
     *     NLINCG,UPD1,YKSK,GSK,YRSR,LRESET,SFUN,.FALSE.,IPIVOT,
     *     ACCRCY,GTPNEW,GNORM,XNORM)
20    CONTINUE
      CALL DCOPY(N,G,1,W(LOLDG),1)
      PNORM = DNRM2(N,W(LPK),1)
      OLDF = FNEW
      OLDGTP = GTPNEW
C
C PREPARE TO COMPUTE THE STEP LENGTH
C
      PE = PNORM + EPSMCH
C
C COMPUTE THE ABSOLUTE AND RELATIVE TOLERANCES FOR THE LINEAR SEARCH
C
      RELTOL = RTEPS*(XNORM + ONE)/PE
      ABSTOL = - EPSMCH*FTEST/(OLDGTP - EPSMCH)
C
C COMPUTE THE SMALLEST ALLOWABLE SPACING BETWEEN POINTS IN
C THE LINEAR SEARCH
C
      TNYTOL = EPSMCH*(XNORM + ONE)/PE
      SPE = STEPMX/PE
C
C SET THE INITIAL STEP LENGTH.
C
      ALPHA = STEP1(FNEW,FM,OLDGTP,SPE)
C
C PERFORM THE LINEAR SEARCH
C
      CALL LINDER(N,SFUN,SMALL,EPSMCH,RELTOL,ABSTOL,TNYTOL,
     *     ETA,ZERO,SPE,W(LPK),OLDGTP,X,FNEW,ALPHA,G,NUMF,
     *     NWHY,W,LW)
C
      FOLD = FNEW
      NITER = NITER + 1
      NFTOTL = NFTOTL + NUMF
      GTG = DDOT(N,G,1,G,1)
      IF (MSGLVL .GE. 1) WRITE(*,810) MAXFUN,NITER,NFTOTL,NLINCG,FNEW,GTG
      IF (NWHY .LT. 0) GO TO 120
      IF (NWHY .EQ. 0 .OR. NWHY .EQ. 2) GO TO 30
C
C THE LINEAR SEARCH HAS FAILED TO FIND A LOWER POINT
C
      NWHY = 3
      GO TO 100
30    IF (NWHY .LE. 1) GO TO 40
      CALL SFUN(N,X,FNEW,G)
      NFTOTL = NFTOTL + 1
C
C TERMINATE IF MORE THAN MAXFUN EVALUTATIONS HAVE BEEN MADE
C
40    NWHY = 2
      IF (NFTOTL .GT. MAXFUN) GO TO 110
      NWHY = 0
C
C SET UP PARAMETERS USED IN CONVERGENCE AND RESETTING TESTS
C
      DIFOLD = DIFNEW
      DIFNEW = OLDF - FNEW
C
C IF THIS IS THE FIRST ITERATION OF A NEW CYCLE, COMPUTE THE
C PERCENTAGE REDUCTION FACTOR FOR THE RESETTING TEST.
C
      IF (ICYCLE .NE. 1) GO TO 50
      IF (DIFNEW .GT. 2.0D0 *DIFOLD) EPSRED = EPSRED + EPSRED
      IF (DIFNEW .LT. 5.0D-1*DIFOLD) EPSRED = 5.0D-1*EPSRED
50    CONTINUE
      GNORM = DSQRT(GTG)
      FTEST = ONE + DABS(FNEW)
      XNORM = DNRM2(N,X,1)
C
C TEST FOR CONVERGENCE
C
      IF ((ALPHA*PNORM .LT. TOLEPS*(ONE + XNORM)
     *     .AND. DABS(DIFNEW) .LT. RTLEPS*FTEST
     *     .AND. GTG .LT. PEPS*FTEST*FTEST)
     *     .OR. GTG .LT. 1.D-4*ACCRCY*FTEST*FTEST) GO TO 90
C
C COMPUTE THE CHANGE IN THE ITERATES AND THE CORRESPONDING CHANGE
C IN THE GRADIENTS
C
      ISK = LSK
      IPK = LPK
      IYK = LYK
      IOLDG = LOLDG
      DO 60 I = 1,N
         W(IYK) = G(I) - W(IOLDG)
         W(ISK) = ALPHA*W(IPK)
         IPK = IPK + 1
         ISK = ISK + 1
         IYK = IYK + 1
         IOLDG = IOLDG + 1
60    CONTINUE
C
C SET UP PARAMETERS USED IN UPDATING THE DIRECTION OF SEARCH.
C
      YKSK = DDOT(N,W(LYK),1,W(LSK),1)
      LRESET = .FALSE.
      IF (ICYCLE .EQ. NM1 .OR. DIFNEW .LT.
     *     EPSRED*(FKEEP-FNEW)) LRESET = .TRUE.
      IF (LRESET) GO TO 70
      YRSR = DDOT(N,W(LYR),1,W(LSR),1)
      IF (YRSR .LE. ZERO) LRESET = .TRUE.
70    CONTINUE
      UPD1 = .FALSE.
C
C      COMPUTE THE NEW SEARCH DIRECTION
C
      MODET = MSGLVL - 3
      CALL MODLNP(MODET,W(LPK),W(LGV),W(LZ1),W(LV),
     *     W(LDIAGB),W(LEMAT),X,G,W(LZK),
     *     N,W,LW,NITER,MAXIT,NFEVAL,NMODIF,
     *     NLINCG,UPD1,YKSK,GSK,YRSR,LRESET,SFUN,.FALSE.,IPIVOT,
     *     ACCRCY,GTPNEW,GNORM,XNORM)
      IF (LRESET) GO TO 80
C
C      STORE THE ACCUMULATED CHANGE IN THE POINT AND GRADIENT AS AN
C      "AVERAGE" DIRECTION FOR PRECONDITIONING.
C
      CALL DXPY(N,W(LSK),1,W(LSR),1)
      CALL DXPY(N,W(LYK),1,W(LYR),1)
      ICYCLE = ICYCLE + 1
      GOTO 20
C
C RESET
C
80    IRESET = IRESET + 1
C
C INITIALIZE THE SUM OF ALL THE CHANGES IN X.
C
      CALL DCOPY(N,W(LSK),1,W(LSR),1)
      CALL DCOPY(N,W(LYK),1,W(LYR),1)
      FKEEP = FNEW
      ICYCLE = 1
      GO TO 20
C
C ...............END OF MAIN ITERATION.......................
C
90    IFAIL = 0
      F = FNEW
      RETURN
100   OLDF = FNEW
C
C LOCAL SEARCH HERE COULD BE INSTALLED HERE
C
110    F = OLDF
C
C SET IFAIL
C
120   IFAIL = NWHY
      RETURN
800   FORMAT(//' NIT   NF   CG', 9X, 'F', 21X, 'GTG',//)
810   FORMAT(' ',I4,1X,I3,1X,I4,1X,I4,1X,1PD22.15,2X,1PD15.8)
      END
C
C
      SUBROUTINE LMQNBC (IFAIL, N, X, F, G, W, LW, SFUN, LOW, UP,
     *   IPIVOT, MSGLVL, MAXIT, MAXFUN, ETA, STEPMX, ACCRCY, XTOL)
      IMPLICIT         DOUBLE PRECISION (A-H,O-Z)
      INTEGER          MSGLVL,N,MAXFUN,IFAIL,LW
      INTEGER          IPIVOT(N)
      DOUBLE PRECISION ETA,XTOL,STEPMX,F,ACCRCY
      DOUBLE PRECISION X(N),G(N),W(LW),LOW(N),UP(N)
C
C THIS ROUTINE IS A BOUNDS-CONSTRAINED TRUNCATED-NEWTON METHOD.
C THE TRUNCATED-NEWTON METHOD IS PRECONDITIONED BY A LIMITED-MEMORY
C QUASI-NEWTON METHOD (THIS PRECONDITIONING STRATEGY IS DEVELOPED
C IN THIS ROUTINE) WITH A FURTHER DIAGONAL SCALING (SEE ROUTINE NDIA3).
C FOR FURTHER DETAILS ON THE PARAMETERS, SEE ROUTINE TNBC.
C
      INTEGER I, ICYCLE, IOLDG, IPK, IYK, LOLDG, LPK, LSR,
     *     LWTEST, LYK, LYR, NFTOTL, NITER, NM1, NUMF, NWHY
      DOUBLE PRECISION ABSTOL, ALPHA, DIFNEW, DIFOLD, EPSMCH, EPSRED,
     *     FKEEP, FLAST, FM, FNEW, FOLD, FSTOP, FTEST, GNORM, GSK,
     *     GTG, GTPNEW, OLDF, OLDGTP, ONE, PE, PEPS, PNORM, RELTOL,
     *     RTEPS, RTLEPS, RTOL, RTOLSQ, SMALL, SPE, TINY,
     *     TNYTOL, TOLEPS, XNORM, YKSK, YRSR, ZERO
      LOGICAL CONV, LRESET, UPD1, NEWCON
C
C THE FOLLOWING STANDARD FUNCTIONS AND SYSTEM FUNCTIONS ARE USED
C
      DOUBLE PRECISION DABS, DDOT, DNRM2, DSQRT, STEP1
      EXTERNAL SFUN
      COMMON/SUBSCR/ LGV, LZ1, LZK, LV, LSK, LYK, LDIAGB, LSR, LYR,
     *     LOLDG, LHG, LHYK, LPK, LEMAT, LWTEST
C
C CHECK THAT INITIAL X IS FEASIBLE AND THAT THE BOUNDS ARE CONSISTENT
C
      CALL CRASH(N,X,IPIVOT,LOW,UP,IER)
      IF (IER .NE. 0) WRITE(*,800)
      IF (IER .NE. 0) RETURN
      IF (MSGLVL .GE. 1) WRITE(*,810)
C
C INITIALIZE VARIABLES
C
      CALL SETPAR(N)
      UPD1 = .TRUE.
      IRESET = 0
      NFEVAL = 0
      NMODIF = 0
      NLINCG = 0
      FSTOP = F
      CONV = .FALSE.
      ZERO = 0.D0
      ONE = 1.D0
      NM1 = N - 1
C
C WITHIN THIS ROUTINE THE ARRAY W(LOLDG) IS SHARED BY W(LHYR)
C
      LHYR = LOLDG
C
C CHECK PARAMETERS AND SET CONSTANTS
C
      CALL CHKUCP(LWTEST,MAXFUN,NWHY,N,ALPHA,EPSMCH,
     *     ETA,PEPS,RTEPS,RTOL,RTOLSQ,STEPMX,FTEST,
     *     XTOL,XNORM,X,LW,SMALL,TINY,ACCRCY)
      IF (NWHY .LT. 0) GO TO 160
      CALL SETUCR(SMALL,NFTOTL,NITER,N,F,FNEW,
     *     FM,GTG,OLDF,SFUN,G,X)
      FOLD = FNEW
      FLAST = FNEW
C
C TEST THE LAGRANGE MULTIPLIERS TO SEE IF THEY ARE NON-NEGATIVE.
C BECAUSE THE CONSTRAINTS ARE ONLY LOWER BOUNDS, THE COMPONENTS
C OF THE GRADIENT CORRESPONDING TO THE ACTIVE CONSTRAINTS ARE THE
C LAGRANGE MULTIPLIERS.  AFTERWORDS, THE PROJECTED GRADIENT IS FORMED.
C
      DO 10 I = 1,N
         IF (IPIVOT(I) .EQ. 2) GO TO 10
         IF (-IPIVOT(I)*G(I) .GE. 0.D0) GO TO 10
         IPIVOT(I) = 0
10    CONTINUE
      CALL ZTIME(N,G,IPIVOT)
      GTG = DDOT(N,G,1,G,1)
      IF (MSGLVL .GE. 1)
     *    CALL MONIT(N,X,FNEW,G,NITER,NFTOTL,NFEVAL,LRESET,IPIVOT)
C
C CHECK IF THE INITIAL POINT IS A LOCAL MINIMUM.
C
      FTEST = ONE + DABS(FNEW)
      IF (GTG .LT. 1.D-4*EPSMCH*FTEST*FTEST) GO TO 130
C
C SET INITIAL VALUES TO OTHER PARAMETERS
C
      ICYCLE = NM1
      TOLEPS = RTOL + RTEPS
      RTLEPS = RTOLSQ + EPSMCH
      GNORM  = DSQRT(GTG)
      DIFNEW = ZERO
      EPSRED = 5.0D-2
      FKEEP  = FNEW
C
C SET THE DIAGONAL OF THE APPROXIMATE HESSIAN TO UNITY.
C
      IDIAGB = LDIAGB
      DO 15 I = 1,N
         W(IDIAGB) = ONE
         IDIAGB = IDIAGB + 1
15    CONTINUE
C
C ..................START OF MAIN ITERATIVE LOOP..........
C
C COMPUTE THE NEW SEARCH DIRECTION
C
      MODET = MSGLVL - 3
      CALL MODLNP(MODET,W(LPK),W(LGV),W(LZ1),W(LV),
     *     W(LDIAGB),W(LEMAT),X,G,W(LZK),
     *     N,W,LW,NITER,MAXIT,NFEVAL,NMODIF,
     *     NLINCG,UPD1,YKSK,GSK,YRSR,LRESET,SFUN,.TRUE.,IPIVOT,
     *     ACCRCY,GTPNEW,GNORM,XNORM)
20    CONTINUE
      CALL DCOPY(N,G,1,W(LOLDG),1)
      PNORM = DNRM2(N,W(LPK),1)
      OLDF = FNEW
      OLDGTP = GTPNEW
C
C PREPARE TO COMPUTE THE STEP LENGTH
C
      PE = PNORM + EPSMCH
C
C COMPUTE THE ABSOLUTE AND RELATIVE TOLERANCES FOR THE LINEAR SEARCH
C
      RELTOL = RTEPS*(XNORM + ONE)/PE
      ABSTOL = - EPSMCH*FTEST/(OLDGTP - EPSMCH)
C
C COMPUTE THE SMALLEST ALLOWABLE SPACING BETWEEN POINTS IN
C THE LINEAR SEARCH
C
      TNYTOL = EPSMCH*(XNORM + ONE)/PE
      CALL STPMAX(STEPMX,PE,SPE,N,X,W(LPK),IPIVOT,LOW,UP)
C
C SET THE INITIAL STEP LENGTH.
C
      ALPHA = STEP1(FNEW,FM,OLDGTP,SPE)
C
C PERFORM THE LINEAR SEARCH
C
      CALL LINDER(N,SFUN,SMALL,EPSMCH,RELTOL,ABSTOL,TNYTOL,
     *     ETA,ZERO,SPE,W(LPK),OLDGTP,X,FNEW,ALPHA,G,NUMF,
     *     NWHY,W,LW)
      NEWCON = .FALSE.
      IF (DABS(ALPHA-SPE) .GT. 1.D1*EPSMCH) GO TO 30
      NEWCON = .TRUE.
      NWHY   = 0
      CALL MODZ(N,X,W(LPK),IPIVOT,EPSMCH,LOW,UP,FLAST,FNEW)
      FLAST = FNEW
C
30    IF (MSGLVL .GE. 3) WRITE(*,820) ALPHA,PNORM
      FOLD = FNEW
      NITER = NITER + 1
      NFTOTL = NFTOTL + NUMF
C
C IF REQUIRED, PRINT THE DETAILS OF THIS ITERATION
C
      IF (MSGLVL .GE. 1)
     *    CALL MONIT(N,X,FNEW,G,NITER,NFTOTL,NFEVAL,LRESET,IPIVOT)
      IF (NWHY .LT. 0) GO TO 160
      IF (NWHY .EQ. 0 .OR. NWHY .EQ. 2) GO TO 40
C
C THE LINEAR SEARCH HAS FAILED TO FIND A LOWER POINT
C
      NWHY = 3
      GO TO 140
40    IF (NWHY .LE. 1) GO TO 50
      CALL SFUN(N,X,FNEW,G)
      NFTOTL = NFTOTL + 1
C
C TERMINATE IF MORE THAN MAXFUN EVALUATIONS HAVE BEEN MADE
C
50    NWHY = 2

C *** MJS ***
C     WRITE(*,57)NFTOTL,MAXFUN,NUMF
57    FORMAT(1X,I6,1X,I6,1X,I6)

      IF (NFTOTL .GT. MAXFUN) GO TO 150
C NEXT LINE ***MJS*** THIS WAS ADDED TO PREVENT PROBLEMS WITH 0-LENGTH BRANCHES. WITH TWO 0-LENGTH SISTERS UNDER LF/PL OR 
C ONE TERMINAL UNDER PL, THE LINESEARCH KEPT ITERATING WITHOUT MAKING PROGRESS FOR MANY STARTING POINTS. I THINK WE CAN BEAT
C THIS BY TERMINATING AND THEN DOING A PERTURBATION AND RESTART
      IF (NITER .GT. MAXFUN) GO TO 150
      NWHY = 0
C
C SET UP PARAMETERS USED IN CONVERGENCE AND RESETTING TESTS
C
      DIFOLD = DIFNEW
      DIFNEW = OLDF - FNEW
C
C IF THIS IS THE FIRST ITERATION OF A NEW CYCLE, COMPUTE THE
C PERCENTAGE REDUCTION FACTOR FOR THE RESETTING TEST.
C
      IF (ICYCLE .NE. 1) GO TO 60
      IF (DIFNEW .GT. 2.D0*DIFOLD) EPSRED = EPSRED + EPSRED
      IF (DIFNEW .LT. 5.0D-1*DIFOLD) EPSRED = 5.0D-1*EPSRED
60    CALL DCOPY(N,G,1,W(LGV),1)
      CALL ZTIME(N,W(LGV),IPIVOT)
      GTG = DDOT(N,W(LGV),1,W(LGV),1)
      GNORM = DSQRT(GTG)
      FTEST = ONE + DABS(FNEW)
      XNORM = DNRM2(N,X,1)
C
C TEST FOR CONVERGENCE
C
      CALL CNVTST(CONV,ALPHA,PNORM,TOLEPS,XNORM,DIFNEW,RTLEPS,
     *     FTEST,GTG,PEPS,EPSMCH,GTPNEW,FNEW,FLAST,G,IPIVOT,N,ACCRCY)
      IF (CONV) GO TO 130
      CALL ZTIME(N,G,IPIVOT)
C
C COMPUTE THE CHANGE IN THE ITERATES AND THE CORRESPONDING CHANGE
C IN THE GRADIENTS
C
      IF (NEWCON) GO TO 90
      ISK = LSK
      IPK = LPK
      IYK = LYK
      IOLDG = LOLDG
      DO 70 I = 1,N
         W(IYK) = G(I) - W(IOLDG)
         W(ISK) = ALPHA*W(IPK)
         IPK = IPK + 1
         ISK = ISK + 1
         IYK = IYK + 1
         IOLDG = IOLDG + 1
70    CONTINUE
C
C SET UP PARAMETERS USED IN UPDATING THE PRECONDITIONING STRATEGY.
C
      YKSK = DDOT(N,W(LYK),1,W(LSK),1)
      LRESET = .FALSE.
      IF (ICYCLE .EQ. NM1 .OR. DIFNEW .LT.
     *     EPSRED*(FKEEP-FNEW)) LRESET = .TRUE.
      IF (LRESET) GO TO 80
      YRSR = DDOT(N,W(LYR),1,W(LSR),1)
      IF (YRSR .LE. ZERO) LRESET = .TRUE.
80    CONTINUE
      UPD1 = .FALSE.
C
C      COMPUTE THE NEW SEARCH DIRECTION
C
90    IF (UPD1 .AND. MSGLVL .GE. 3) WRITE(*,830)
      IF (NEWCON .AND. MSGLVL .GE. 3) WRITE(*,840)
      MODET = MSGLVL - 3
      CALL MODLNP(MODET,W(LPK),W(LGV),W(LZ1),W(LV),
     *     W(LDIAGB),W(LEMAT),X,G,W(LZK),
     *     N,W,LW,NITER,MAXIT,NFEVAL,NMODIF,
     *     NLINCG,UPD1,YKSK,GSK,YRSR,LRESET,SFUN,.TRUE.,IPIVOT,
     *     ACCRCY,GTPNEW,GNORM,XNORM)
      IF (NEWCON) GO TO 20
      IF (LRESET) GO TO 110
C
C COMPUTE THE ACCUMULATED STEP AND ITS CORRESPONDING
C GRADIENT DIFFERENCE.
C
      CALL DXPY(N,W(LSK),1,W(LSR),1)
      CALL DXPY(N,W(LYK),1,W(LYR),1)
      ICYCLE = ICYCLE + 1
      GOTO 20
C
C RESET
C
110   IRESET = IRESET + 1
C
C INITIALIZE THE SUM OF ALL THE CHANGES IN X.
C
      CALL DCOPY(N,W(LSK),1,W(LSR),1)
      CALL DCOPY(N,W(LYK),1,W(LYR),1)
      FKEEP = FNEW
      ICYCLE = 1
      GO TO 20
C
C ...............END OF MAIN ITERATION.......................
C
130   IFAIL = 0
      F = FNEW
      RETURN
140   OLDF = FNEW
C
C LOCAL SEARCH COULD BE INSTALLED HERE
C
150   F = OLDF
      IF (MSGLVL .GE. 1) CALL MONIT(N,X,
     *     F,G,NITER,NFTOTL,NFEVAL,IRESET,IPIVOT)
C
C SET IFAIL
C
160   IFAIL = NWHY
      RETURN
800   FORMAT(' THERE IS NO FEASIBLE POINT; TERMINATING ALGORITHM')
810   FORMAT(//'  NIT   NF   CG', 9X, 'F', 21X, 'GTG',//)
820   FORMAT('        LINESEARCH RESULTS:  ALPHA,PNORM',2(1PD12.4))
830   FORMAT(' UPD1 IS TRUE - TRIVIAL PRECONDITIONING')
840   FORMAT(' NEWCON IS TRUE - CONSTRAINT ADDED IN LINESEARCH')
      END
C
C
      SUBROUTINE MONIT(N,X,F,G,NITER,NFTOTL,NFEVAL,IRESET,IPIVOT)
C
C PRINT RESULTS OF CURRENT ITERATION
C
      IMPLICIT         DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N),F,G(N),GTG
      INTEGER          IPIVOT(N)
C     LOGICAL		IRESET  *** MJS added this line. BUG HERE. Most of the callers pass a logical IRESET, but some pass an integer

C
      GTG = 0.D0
      DO 10 I = 1,N
         IF (IPIVOT(I) .NE. 0) GO TO 10
         GTG = GTG + G(I)*G(I)
10    CONTINUE
      WRITE(*,800) NITER,NFTOTL,NFEVAL,F,GTG
      RETURN
800   FORMAT(' ',I4,1X,I4,1X,I4,1X,1PD22.15,2X,1PD15.8)
      END
C
C
      SUBROUTINE ZTIME(N,X,IPIVOT)
      IMPLICIT         DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N)
      INTEGER          IPIVOT(N)
C
C THIS ROUTINE MULTIPLIES THE VECTOR X BY THE CONSTRAINT MATRIX Z
C
      DO 10 I = 1,N
         IF (IPIVOT(I) .NE. 0) X(I) = 0.D0
10    CONTINUE
      RETURN
      END
C
C
      SUBROUTINE STPMAX(STEPMX,PE,SPE,N,X,P,IPIVOT,LOW,UP)
      IMPLICIT         DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION LOW(N),UP(N),X(N),P(N),STEPMX,PE,SPE,T
      INTEGER          IPIVOT(N)
C
C COMPUTE THE MAXIMUM ALLOWABLE STEP LENGTH
C
      SPE = STEPMX / PE
C SPE IS THE STANDARD (UNCONSTRAINED) MAX STEP
      DO 10 I = 1,N
         IF (IPIVOT(I) .NE. 0) GO TO 10
         IF (P(I) .EQ. 0.D0) GO TO 10
         IF (P(I) .GT. 0.D0) GO TO 5
         T = LOW(I) - X(I)
         IF (T .GT. SPE*P(I)) SPE = T / P(I)
         GO TO 10
5        T = UP(I) - X(I)
         IF (T .LT. SPE*P(I)) SPE = T / P(I)
10    CONTINUE
      RETURN
      END
C
C
      SUBROUTINE MODZ(N,X,P,IPIVOT,EPSMCH,LOW,UP,FLAST,FNEW)
      IMPLICIT         DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N), P(N), EPSMCH, DABS, TOL, LOW(N), UP(N),
     *                 FLAST, FNEW
      INTEGER          IPIVOT(N)
C
C UPDATE THE CONSTRAINT MATRIX IF A NEW CONSTRAINT IS ENCOUNTERED
C
      DO 10 I = 1,N
         IF (IPIVOT(I) .NE. 0) GO TO 10
         IF (P(I) .EQ. 0.D0) GO TO 10
         IF (P(I) .GT. 0.D0) GO TO 5
         TOL = 1.D1 * EPSMCH * (DABS(LOW(I)) + 1.D0)
         IF (X(I)-LOW(I) .GT. TOL) GO TO 10
         FLAST = FNEW
         IPIVOT(I) = -1
         X(I) = LOW(I)
         GO TO 10
5        TOL = 1.D1 * EPSMCH * (DABS(UP(I)) + 1.D0)
         IF (UP(I)-X(I) .GT. TOL) GO TO 10
         FLAST = FNEW
         IPIVOT(I) = 1
         X(I) = UP(I)
10    CONTINUE
      RETURN
      END
C
C
      SUBROUTINE CNVTST(CONV,ALPHA,PNORM,TOLEPS,XNORM,DIFNEW,RTLEPS,
     *     FTEST,GTG,PEPS,EPSMCH,GTPNEW,FNEW,FLAST,G,IPIVOT,N,ACCRCY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL CONV,LTEST
      INTEGER IPIVOT(N)
      DOUBLE PRECISION G(N), ALPHA, PNORM, TOLEPS, XNORM, DIFNEW,
     *     RTLEPS, FTEST, GTG, PEPS, EPSMCH, GTPNEW, FNEW, FLAST, ONE,
     *     CMAX, T, ACCRCY
C
C TEST FOR CONVERGENCE
C
      IMAX = 0
      CMAX = 0.D0
      LTEST = FLAST - FNEW .LE. -5.D-1*GTPNEW
      DO 10 I = 1,N
         IF (IPIVOT(I) .EQ. 0 .OR. IPIVOT(I) .EQ. 2) GO TO 10
         T = -IPIVOT(I)*G(I)
         IF (T .GE. 0.D0) GO TO 10
         CONV = .FALSE.
         IF (LTEST) GO TO 10
         IF (CMAX .LE. T) GO TO 10
         CMAX = T
         IMAX = I
10    CONTINUE
      IF (IMAX .EQ. 0) GO TO 15
      IPIVOT(IMAX) = 0
      FLAST = FNEW
      RETURN
15    CONTINUE
      CONV = .FALSE.
      ONE = 1.D0
      IF ((ALPHA*PNORM .GE. TOLEPS*(ONE + XNORM)
     *     .OR. DABS(DIFNEW) .GE. RTLEPS*FTEST
     *     .OR. GTG .GE. PEPS*FTEST*FTEST)
     *     .AND. GTG .GE. 1.D-4*ACCRCY*FTEST*FTEST) RETURN
      CONV = .TRUE.
C
C FOR DETAILS, SEE GILL, MURRAY, AND WRIGHT (1981, P. 308) AND
C FLETCHER (1981, P. 116).  THE MULTIPLIER TESTS (HERE, TESTING
C THE SIGN OF THE COMPONENTS OF THE GRADIENT) MAY STILL NEED TO
C MODIFIED TO INCORPORATE TOLERANCES FOR ZERO.
C
      RETURN
      END
C
C
      SUBROUTINE CRASH(N,X,IPIVOT,LOW,UP,IER)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X(N),LOW(N),UP(N)
      INTEGER IPIVOT(N)
C
C THIS INITIALIZES THE CONSTRAINT INFORMATION, AND ENSURES THAT THE
C INITIAL POINT SATISFIES  LOW <= X <= UP.
C THE CONSTRAINTS ARE CHECKED FOR CONSISTENCY.
C
      IER = 0
      DO 30 I = 1,N
         IF (X(I) .LT. LOW(I)) X(I) = LOW(I)
         IF (X(I) .GT. UP(I)) X(I) = UP(I)
         IPIVOT(I) = 0
         IF (X(I) .EQ. LOW(I)) IPIVOT(I) = -1
         IF (X(I) .EQ. UP(I)) IPIVOT(I) = 1
         IF (UP(I) .EQ. LOW(I)) IPIVOT(I) = 2
         IF (LOW(I) .GT. UP(I)) IER = -I
30    CONTINUE
      RETURN
      END
C
C THE VECTORS SK AND YK, ALTHOUGH NOT IN THE CALL,
C ARE USED (VIA THEIR POSITION IN W) BY THE ROUTINE MSOLVE.
C
      SUBROUTINE MODLNP(MODET,ZSOL,GV,R,V,DIAGB,EMAT,
     *     X,G,ZK,N,W,LW,NITER,MAXIT,NFEVAL,NMODIF,NLINCG,
     *     UPD1,YKSK,GSK,YRSR,LRESET,SFUN,BOUNDS,IPIVOT,ACCRCY,
     *     GTP,GNORM,XNORM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER MODET,N,NITER,IPIVOT(1)
      DOUBLE PRECISION ZSOL(N),G(N),GV(N),R(N),V(N),DIAGB(N),W(LW)
      DOUBLE PRECISION EMAT(N),ZK(N),X(N),ACCRCY
      DOUBLE PRECISION ALPHA,BETA,DELTA,GSK,GTP,PR,
     *     QOLD,QNEW,QTEST,RHSNRM,RNORM,RZ,RZOLD,TOL,VGV,YKSK,YRSR
      DOUBLE PRECISION GNORM,XNORM
      DOUBLE PRECISION DDOT,DNRM2
      LOGICAL FIRST,UPD1,LRESET,BOUNDS
      EXTERNAL SFUN
C
C THIS ROUTINE PERFORMS A PRECONDITIONED CONJUGATE-GRADIENT
C ITERATION IN ORDER TO SOLVE THE NEWTON EQUATIONS FOR A SEARCH
C DIRECTION FOR A TRUNCATED-NEWTON ALGORITHM.  WHEN THE VALUE OF THE
C QUADRATIC MODEL IS SUFFICIENTLY REDUCED,
C THE ITERATION IS TERMINATED.
C
C PARAMETERS
C
C MODET       - INTEGER WHICH CONTROLS AMOUNT OF OUTPUT
C ZSOL        - COMPUTED SEARCH DIRECTION
C G           - CURRENT GRADIENT
C GV,GZ1,V    - SCRATCH VECTORS
C R           - RESIDUAL
C DIAGB,EMAT  - DIAGONAL PRECONDITONING MATRIX
C NITER       - NONLINEAR ITERATION #
C FEVAL       - VALUE OF QUADRATIC FUNCTION
C
C *************************************************************
C INITIALIZATION
C *************************************************************
C
C GENERAL INITIALIZATION
C
      IF (MODET .GT. 0) WRITE(*,800)
      IF (MAXIT .EQ. 0) RETURN
      FIRST = .TRUE.
      RHSNRM = GNORM
      TOL = 1.D-12
      QOLD = 0.D0
C
C INITIALIZATION FOR PRECONDITIONED CONJUGATE-GRADIENT ALGORITHM
C
      CALL INITPC(DIAGB,EMAT,N,W,LW,MODET,
     *            UPD1,YKSK,GSK,YRSR,LRESET)
      DO 10 I = 1,N
         R(I) = -G(I)
         V(I) = 0.D0
         ZSOL(I) = 0.D0
10    CONTINUE
C
C ************************************************************
C MAIN ITERATION
C ************************************************************
C
      DO 30 K = 1,MAXIT
         NLINCG = NLINCG + 1
         IF (MODET .GT. 1) WRITE(*,810) K
C
C CG ITERATION TO SOLVE SYSTEM OF EQUATIONS
C
         IF (BOUNDS) CALL ZTIME(N,R,IPIVOT)
         CALL MSOLVE(R,ZK,N,W,LW,UPD1,YKSK,GSK,
     *                 YRSR,LRESET,FIRST)
         IF (BOUNDS) CALL ZTIME(N,ZK,IPIVOT)
         RZ = DDOT(N,R,1,ZK,1)
         IF (RZ/RHSNRM .LT. TOL) GO TO 80
         IF (K .EQ. 1) BETA = 0.D0
         IF (K .GT. 1) BETA = RZ/RZOLD
         DO 20 I = 1,N
            V(I) = ZK(I) + BETA*V(I)
20       CONTINUE
         IF (BOUNDS) CALL ZTIME(N,V,IPIVOT)
         CALL GTIMS(V,GV,N,X,G,W,LW,SFUN,FIRST,DELTA,ACCRCY,XNORM)
         IF (BOUNDS) CALL ZTIME(N,GV,IPIVOT)
         NFEVAL = NFEVAL + 1
         VGV = DDOT(N,V,1,GV,1)
         IF (VGV/RHSNRM .LT. TOL) GO TO 50
         CALL NDIA3(N,EMAT,V,GV,R,VGV,MODET)
C
C COMPUTE LINEAR STEP LENGTH
C
         ALPHA = RZ / VGV
         IF (MODET .GE. 1) WRITE(*,820) ALPHA
C
C COMPUTE CURRENT SOLUTION AND RELATED VECTORS
C
         CALL DAXPY(N,ALPHA,V,1,ZSOL,1)
         CALL DAXPY(N,-ALPHA,GV,1,R,1)
C
C TEST FOR CONVERGENCE
C
         GTP = DDOT(N,ZSOL,1,G,1)
         PR = DDOT(N,R,1,ZSOL,1)
         QNEW = 5.D-1 * (GTP + PR)
         QTEST = K * (1.D0 - QOLD/QNEW)
         IF (QTEST .LT. 0.D0) GO TO 70
         QOLD = QNEW
         IF (QTEST .LE. 5.D-1) GO TO 70
C
C PERFORM CAUTIONARY TEST
C
         IF (GTP .GT. 0) GO TO 40
         RZOLD = RZ
30    CONTINUE
C
C TERMINATE ALGORITHM
C
      K = K-1
      GO TO 70
C
C TRUNCATE ALGORITHM IN CASE OF AN EMERGENCY
C
40    IF (MODET .GE. -1) WRITE(*,830) K
      CALL DAXPY(N,-ALPHA,V,1,ZSOL,1)
      GTP = DDOT(N,ZSOL,1,G,1)
      GO TO 90
50    CONTINUE
      IF (MODET .GT. -2) WRITE(*,840)
60    IF (K .GT. 1) GO TO 70
      CALL MSOLVE(G,ZSOL,N,W,LW,UPD1,YKSK,GSK,YRSR,LRESET,FIRST)
      CALL NEGVEC(N,ZSOL)
      IF (BOUNDS) CALL ZTIME(N,ZSOL,IPIVOT)
      GTP = DDOT(N,ZSOL,1,G,1)
70    CONTINUE
      IF (MODET .GE. -1) WRITE(*,850) K,RNORM
      GO TO 90
80    CONTINUE
      IF (MODET .GE. -1) WRITE(*,860)
      IF (K .GT. 1) GO TO 70
      CALL DCOPY(N,G,1,ZSOL,1)
      CALL NEGVEC(N,ZSOL)
      IF (BOUNDS) CALL ZTIME(N,ZSOL,IPIVOT)
      GTP = DDOT(N,ZSOL,1,G,1)
      GO TO 70
C
C STORE (OR RESTORE) DIAGONAL PRECONDITIONING
C
90    CONTINUE
      CALL DCOPY(N,EMAT,1,DIAGB,1)
      RETURN
800   FORMAT(' ',//,' ENTERING MODLNP')
810   FORMAT(' ',//,' ### ITERATION ',I2,' ###')
820   FORMAT(' ALPHA',1PD16.8)
830   FORMAT(' G(T)Z POSITIVE AT ITERATION ',I2,
     *     ' - TRUNCATING METHOD',/)
840   FORMAT(' ',10X,'HESSIAN NOT POSITIVE-DEFINITE')
850   FORMAT(' ',/,8X,'MODLAN TRUNCATED AFTER ',I3,' ITERATIONS',
     *     '  RNORM = ',1PD14.6)
860   FORMAT(' PRECONDITIONING NOT POSITIVE-DEFINITE')
      END
C
C
      SUBROUTINE NDIA3(N,E,V,GV,R,VGV,MODET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION E(N),V(N),GV(N),R(N),VGV,VR,DDOT
C
C UPDATE THE PRECONDITIOING MATRIX BASED ON A DIAGONAL VERSION
C OF THE BFGS QUASI-NEWTON UPDATE.
C
      VR = DDOT(N,V,1,R,1)
      DO 10 I = 1,N
         E(I) = E(I) - R(I)*R(I)/VR + GV(I)*GV(I)/VGV
         IF (E(I) .GT. 1.D-6) GO TO 10
         IF (MODET .GT. -2) WRITE(*,800) E(I)
         E(I) = 1.D0
10    CONTINUE
      RETURN
800   FORMAT(' *** EMAT NEGATIVE:  ',1PD16.8)
      END
C
C      SERVICE ROUTINES FOR OPTIMIZATION
C
      SUBROUTINE NEGVEC(N,V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N
      DOUBLE PRECISION V(N)
C
C NEGATIVE OF THE VECTOR V
C
      INTEGER I
      DO 10 I = 1,N
         V(I) = -V(I)
10    CONTINUE
      RETURN
      END
C
C
      SUBROUTINE LSOUT(ILOC,ITEST,XMIN,FMIN,GMIN,XW,FW,GW,U,A,
     *     B,TOL,EPS,SCXBD,XLAMDA)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION XMIN,FMIN,GMIN,XW,FW,GW,U,A,B,
     *     TOL,EPS,SCXBD,XLAMDA
C
C ERROR PRINTOUTS FOR GETPTC
C
      DOUBLE PRECISION YA,YB,YBND,YW,YU
      YU = XMIN + U
      YA = A + XMIN
      YB = B + XMIN
      YW = XW + XMIN
      YBND = SCXBD + XMIN
      WRITE(*,800)
      WRITE(*,810) TOL,EPS
      WRITE(*,820) YA,YB
      WRITE(*,830) YBND
      WRITE(*,840) YW,FW,GW
      WRITE(*,850) XMIN,FMIN,GMIN
      WRITE(*,860) YU
      WRITE(*,870) ILOC,ITEST
      RETURN
800   FORMAT(///' OUTPUT FROM LINEAR SEARCH')
810   FORMAT('  TOL AND EPS'/2D25.14)
820   FORMAT('  CURRENT UPPER AND LOWER BOUNDS'/2D25.14)
830   FORMAT('  STRICT UPPER BOUND'/D25.14)
840   FORMAT('  XW, FW, GW'/3D25.14)
850   FORMAT('  XMIN, FMIN, GMIN'/3D25.14)
860   FORMAT('  NEW ESTIMATE'/2D25.14)
870   FORMAT('  ILOC AND ITEST'/2I3)
      END
C
C
      DOUBLE PRECISION FUNCTION STEP1(FNEW,FM,GTP,SMAX)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION FNEW,FM,GTP,SMAX
C
C ********************************************************
C STEP1 RETURNS THE LENGTH OF THE INITIAL STEP TO BE TAKEN ALONG THE
C VECTOR P IN THE NEXT LINEAR SEARCH.
C ********************************************************
C
      DOUBLE PRECISION ALPHA,D,EPSMCH
      DOUBLE PRECISION DABS,MCHPR1
      EPSMCH = MCHPR1()
      D = DABS(FNEW-FM)
      ALPHA = 1.D0
      IF (2.D0*D .LE. (-GTP) .AND. D .GE. EPSMCH)
     *     ALPHA = -2.D0*D/GTP
      IF (ALPHA .GE. SMAX) ALPHA = SMAX
      STEP1 = ALPHA
      RETURN
      END
C
C
      DOUBLE PRECISION FUNCTION MCHPR1()
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION X
C
C RETURNS THE VALUE OF EPSMCH, WHERE EPSMCH IS THE SMALLEST POSSIBLE
C REAL NUMBER SUCH THAT 1.0 + EPSMCH .GT. 1.0
C
C FOR VAX
C
      MCHPR1 = 1.D-17
C
C FOR SUN
C
C     MCHPR1 = 1.0842021724855D-19
      RETURN
      END
C
C
      SUBROUTINE CHKUCP(LWTEST,MAXFUN,NWHY,N,ALPHA,EPSMCH,
     *     ETA,PEPS,RTEPS,RTOL,RTOLSQ,STEPMX,TEST,
     *     XTOL,XNORM,X,LW,SMALL,TINY,ACCRCY)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LW,LWTEST,MAXFUN,NWHY,N
      DOUBLE PRECISION ACCRCY,ALPHA,EPSMCH,ETA,PEPS,RTEPS,RTOL,
     *     RTOLSQ,STEPMX,TEST,XTOL,XNORM,SMALL,TINY
      DOUBLE PRECISION X(N)
C
C CHECKS PARAMETERS AND SETS CONSTANTS WHICH ARE COMMON TO BOTH
C DERIVATIVE AND NON-DERIVATIVE ALGORITHMS
C
      DOUBLE PRECISION DABS,DSQRT,MCHPR1
      EPSMCH = MCHPR1()
      SMALL = EPSMCH*EPSMCH
      TINY = SMALL
      NWHY = -1
      RTEPS = DSQRT(EPSMCH)
      RTOL = XTOL
      IF (DABS(RTOL) .LT. ACCRCY) RTOL = 1.D1*RTEPS
C
C CHECK FOR ERRORS IN THE INPUT PARAMETERS
C
      IF (LW .LT. LWTEST
     *      .OR. N .LT. 1 .OR. RTOL .LT. 0.D0 .OR. ETA .GE. 1.D0 .OR.
     *      ETA .LT. 0.D0 .OR. STEPMX .LT. RTOL .OR.
     *      MAXFUN .LT. 1) RETURN
      NWHY = 0
C
C SET CONSTANTS FOR LATER
C
      RTOLSQ = RTOL*RTOL
      PEPS = ACCRCY**0.6666D0
      XNORM = DNRM2(N,X,1)
      ALPHA = 0.D0
      TEST = 0.D0
      RETURN
      END
C
C
      SUBROUTINE SETUCR(SMALL,NFTOTL,NITER,N,F,FNEW,
     *            FM,GTG,OLDF,SFUN,G,X)
      IMPLICIT         DOUBLE PRECISION (A-H,O-Z)
      INTEGER          NFTOTL,NITER,N
      DOUBLE PRECISION F,FNEW,FM,GTG,OLDF,SMALL
      DOUBLE PRECISION G(N),X(N)
      EXTERNAL         SFUN
C
C CHECK INPUT PARAMETERS, COMPUTE THE INITIAL FUNCTION VALUE, SET
C CONSTANTS FOR THE SUBSEQUENT MINIMIZATION
C
      FM = F
C
C COMPUTE THE INITIAL FUNCTION VALUE
C
      CALL SFUN(N,X,FNEW,G)
      NFTOTL = 1
C
C SET CONSTANTS FOR LATER
C
      NITER = 0
      OLDF = FNEW
      GTG = DDOT(N,G,1,G,1)
      RETURN
      END
C
C
      SUBROUTINE GTIMS(V,GV,N,X,G,W,LW,SFUN,FIRST,DELTA,ACCRCY,XNORM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION V(N),GV(N),DINV,DELTA,G(N)
      DOUBLE PRECISION F,X(N),W(LW),ACCRCY,DSQRT,XNORM
      LOGICAL FIRST
      EXTERNAL SFUN
      COMMON/SUBSCR/ LGV,LZ1,LZK,LV,LSK,LYK,LDIAGB,LSR,LYR,
     *     LHYR,LHG,LHYK,LPK,LEMAT,LWTEST
C
C THIS ROUTINE COMPUTES THE PRODUCT OF THE MATRIX G TIMES THE VECTOR
C V AND STORES THE RESULT IN THE VECTOR GV (FINITE-DIFFERENCE VERSION)
C
      IF (.NOT. FIRST) GO TO 20
      DELTA = DSQRT(ACCRCY)*(1.D0+XNORM)
      FIRST = .FALSE.
20    CONTINUE
      DINV = 1.D0/DELTA
      IHG = LHG
      DO 30 I = 1,N
         W(IHG) = X(I) + DELTA*V(I)
         IHG = IHG + 1
30    CONTINUE
      CALL SFUN(N,W(LHG),F,GV)
      DO 40 I = 1,N
         GV(I) = (GV(I) - G(I))*DINV
40    CONTINUE
      RETURN
      END
C
C
      SUBROUTINE MSOLVE(G,Y,N,W,LW,UPD1,YKSK,GSK,
     *     YRSR,LRESET,FIRST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION G(N),Y(N),W(LW),YKSK,GSK,YRSR
      LOGICAL UPD1,LRESET,FIRST
C
C THIS ROUTINE SETS UPT THE ARRAYS FOR MSLV
C
      COMMON/SUBSCR/ LGV,LZ1,LZK,LV,LSK,LYK,LDIAGB,LSR,LYR,
     *     LHYR,LHG,LHYK,LPK,LEMAT,LWTEST
      CALL MSLV(G,Y,N,W(LSK),W(LYK),W(LDIAGB),W(LSR),W(LYR),W(LHYR),
     *     W(LHG),W(LHYK),UPD1,YKSK,GSK,YRSR,LRESET,FIRST)
      RETURN
      END
      SUBROUTINE MSLV(G,Y,N,SK,YK,DIAGB,SR,YR,HYR,HG,HYK,
     *     UPD1,YKSK,GSK,YRSR,LRESET,FIRST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION G(N),Y(N)
C
C THIS ROUTINE ACTS AS A PRECONDITIONING STEP FOR THE
C LINEAR CONJUGATE-GRADIENT ROUTINE.  IT IS ALSO THE
C METHOD OF COMPUTING THE SEARCH DIRECTION FROM THE
C GRADIENT FOR THE NON-LINEAR CONJUGATE-GRADIENT CODE.
C IT REPRESENTS A TWO-STEP SELF-SCALED BFGS FORMULA.
C
      DOUBLE PRECISION DDOT,YKSK,GSK,YRSR,RDIAGB,YKHYK,GHYK,
     *     YKSR,YKHYR,YRHYR,GSR,GHYR
      DOUBLE PRECISION SK(N),YK(N),DIAGB(N),SR(N),YR(N),HYR(N),HG(N),
     *     HYK(N),ONE
      LOGICAL LRESET,UPD1,FIRST
      IF (UPD1) GO TO 100
      ONE = 1.D0
      GSK = DDOT(N,G,1,SK,1)
      IF (LRESET) GO TO 60
C
C COMPUTE HG AND HY WHERE H IS THE INVERSE OF THE DIAGONALS
C
      DO 57 I = 1,N
         RDIAGB = 1.0D0/DIAGB(I)
         HG(I) = G(I)*RDIAGB
         IF (FIRST) HYK(I) = YK(I)*RDIAGB
         IF (FIRST) HYR(I) = YR(I)*RDIAGB
57    CONTINUE
      IF (FIRST) YKSR = DDOT(N,YK,1,SR,1)
      IF (FIRST) YKHYR = DDOT(N,YK,1,HYR,1)
      GSR = DDOT(N,G,1,SR,1)
      GHYR = DDOT(N,G,1,HYR,1)
      IF (FIRST) YRHYR = DDOT(N,YR,1,HYR,1)
      CALL SSBFGS(N,ONE,SR,YR,HG,HYR,YRSR,
     *     YRHYR,GSR,GHYR,HG)
      IF (FIRST) CALL SSBFGS(N,ONE,SR,YR,HYK,HYR,YRSR,
     *     YRHYR,YKSR,YKHYR,HYK)
      YKHYK = DDOT(N,HYK,1,YK,1)
      GHYK = DDOT(N,HYK,1,G,1)
      CALL SSBFGS(N,ONE,SK,YK,HG,HYK,YKSK,
     *     YKHYK,GSK,GHYK,Y)
      RETURN
60    CONTINUE
C
C COMPUTE GH AND HY WHERE H IS THE INVERSE OF THE DIAGONALS
C
      DO 65 I = 1,N
         RDIAGB = 1.D0/DIAGB(I)
         HG(I) = G(I)*RDIAGB
         IF (FIRST) HYK(I) = YK(I)*RDIAGB
65    CONTINUE
      IF (FIRST) YKHYK = DDOT(N,YK,1,HYK,1)
      GHYK = DDOT(N,G,1,HYK,1)
      CALL SSBFGS(N,ONE,SK,YK,HG,HYK,YKSK,
     *     YKHYK,GSK,GHYK,Y)
      RETURN
100   CONTINUE
      DO 110 I = 1,N
110      Y(I) = G(I) / DIAGB(I)
      RETURN
      END
C
C
      SUBROUTINE SSBFGS(N,GAMMA,SJ,YJ,HJV,HJYJ,YJSJ,YJHYJ,
     *     VSJ,VHYJ,HJP1V)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N
      DOUBLE PRECISION GAMMA,YJSJ,YJHYJ,VSJ,VHYJ
      DOUBLE PRECISION SJ(N),YJ(N),HJV(N),HJYJ(N),HJP1V(N)
C
C SELF-SCALED BFGS
C
      INTEGER I
      DOUBLE PRECISION BETA,DELTA
      DELTA = (1.D0 + GAMMA*YJHYJ/YJSJ)*VSJ/YJSJ
     *     - GAMMA*VHYJ/YJSJ
      BETA = -GAMMA*VSJ/YJSJ
      DO 10 I = 1,N
         HJP1V(I) = GAMMA*HJV(I) + DELTA*SJ(I) + BETA*HJYJ(I)
10    CONTINUE
      RETURN
      END
C
C ROUTINES TO INITIALIZE PRECONDITIONER
C
      SUBROUTINE INITPC(DIAGB,EMAT,N,W,LW,MODET,
     *     UPD1,YKSK,GSK,YRSR,LRESET)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DIAGB(N),EMAT(N),W(LW)
      DOUBLE PRECISION YKSK,GSK,YRSR
      LOGICAL LRESET,UPD1
      COMMON/SUBSCR/ LGV,LZ1,LZK,LV,LSK,LYK,LDIAGB,LSR,LYR,
     *     LHYR,LHG,LHYK,LPK,LEMAT,LWTEST
      CALL INITP3(DIAGB,EMAT,N,LRESET,YKSK,YRSR,W(LHYK),
     *     W(LSK),W(LYK),W(LSR),W(LYR),MODET,UPD1)
      RETURN
      END
      SUBROUTINE INITP3(DIAGB,EMAT,N,LRESET,YKSK,YRSR,BSK,
     *     SK,YK,SR,YR,MODET,UPD1)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DOUBLE PRECISION DIAGB(N),EMAT(N),YKSK,YRSR,BSK(N),SK(N),
     *     YK(N),COND,SR(N),YR(N),DDOT,SDS,SRDS,YRSK,TD,D1,DN
      LOGICAL LRESET,UPD1
      IF (UPD1) GO TO 90
      IF (LRESET) GO TO 60
      DO 10 I = 1,N
         BSK(I) = DIAGB(I)*SR(I)
10    CONTINUE
      SDS = DDOT(N,SR,1,BSK,1)
      SRDS = DDOT(N,SK,1,BSK,1)
      YRSK = DDOT(N,YR,1,SK,1)
      DO 20 I = 1,N
         TD = DIAGB(I)
         BSK(I) = TD*SK(I) - BSK(I)*SRDS/SDS+YR(I)*YRSK/YRSR
         EMAT(I) = TD-TD*TD*SR(I)*SR(I)/SDS+YR(I)*YR(I)/YRSR
20    CONTINUE
      SDS = DDOT(N,SK,1,BSK,1)
      DO 30 I = 1,N
         EMAT(I) = EMAT(I) - BSK(I)*BSK(I)/SDS+YK(I)*YK(I)/YKSK
30    CONTINUE
      GO TO 110
60    CONTINUE
      DO 70 I = 1,N
         BSK(I) = DIAGB(I)*SK(I)
70    CONTINUE
      SDS = DDOT(N,SK,1,BSK,1)
      DO 80 I = 1,N
         TD = DIAGB(I)
         EMAT(I) = TD - TD*TD*SK(I)*SK(I)/SDS + YK(I)*YK(I)/YKSK
80    CONTINUE
      GO TO 110
90    CONTINUE
      CALL DCOPY(N,DIAGB,1,EMAT,1)
110   CONTINUE
      IF (MODET .LT. 1) RETURN
      D1 = EMAT(1)
      DN = EMAT(1)
      DO 120 I = 1,N
         IF (EMAT(I) .LT. D1) D1 = EMAT(I)
         IF (EMAT(I) .GT. DN) DN = EMAT(I)
120   CONTINUE
      COND = DN/D1
      WRITE(*,800) D1,DN,COND
800   FORMAT(' ',//8X,'DMIN =',1PD12.4,'  DMAX =',1PD12.4,
     *     ' COND =',1PD12.4,/)
      RETURN
      END
C
C
      SUBROUTINE SETPAR(N)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER LSUB(14)
      COMMON/SUBSCR/ LSUB,LWTEST
C
C SET UP PARAMETERS FOR THE OPTIMIZATION ROUTINE
C
      DO 10 I = 1,14
          LSUB(I) = (I-1)*N + 1
10    CONTINUE
      LWTEST = LSUB(14) + N - 1
      RETURN
      END
C
C      LINE SEARCH ALGORITHMS OF GILL AND MURRAY
C
      SUBROUTINE LINDER(N,SFUN,SMALL,EPSMCH,RELTOL,ABSTOL,
     *     TNYTOL,ETA,SFTBND,XBND,P,GTP,X,F,ALPHA,G,NFTOTL,
     *     IFLAG,W,LW)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      INTEGER N,NFTOTL,IFLAG,LW
      DOUBLE PRECISION SMALL,EPSMCH,RELTOL,ABSTOL,TNYTOL,ETA,
     *     SFTBND,XBND,GTP,F,ALPHA
      DOUBLE PRECISION P(N),X(N),G(N),W(LW)
C
C
      INTEGER I,IENTRY,ITEST,L,LG,LX,NUMF,ITCNT
      DOUBLE PRECISION A,B,B1,BIG,E,FACTOR,FMIN,FPRESN,FU,
     *     FW,GMIN,GTEST1,GTEST2,GU,GW,OLDF,SCXBND,STEP,
     *     TOL,U,XMIN,XW,RMU,RTSMLL,UALPHA
      LOGICAL BRAKTD
C
C      THE FOLLOWING STANDARD FUNCTIONS AND SYSTEM FUNCTIONS ARE
C      CALLED WITHIN LINDER
C
      DOUBLE PRECISION DDOT,DSQRT
      EXTERNAL SFUN
C
C      ALLOCATE THE ADDRESSES FOR LOCAL WORKSPACE
C
      LX = 1
      LG = LX + N
      LSPRNT = 0
      NPRNT  = 10000
      RTSMLL = DSQRT(SMALL)
      BIG = 1.D0/SMALL
      ITCNT = 0
C
C      SET THE ESTIMATED RELATIVE PRECISION IN F(X).
C
      FPRESN = 10.D0*EPSMCH
      NUMF = 0
      U = ALPHA
      FU = F
      FMIN = F
      GU = GTP
      RMU = 1.0D-4
C
C      FIRST ENTRY SETS UP THE INITIAL INTERVAL OF UNCERTAINTY.
C
      IENTRY = 1
10    CONTINUE
C
C TEST FOR TOO MANY ITERATIONS
C
      ITCNT = ITCNT + 1
      IFLAG = 1
      IF (ITCNT .GT. 20) GO TO 50
      IFLAG = 0
      CALL GETPTC(BIG,SMALL,RTSMLL,RELTOL,ABSTOL,TNYTOL,
     *     FPRESN,ETA,RMU,XBND,U,FU,GU,XMIN,FMIN,GMIN,
     *     XW,FW,GW,A,B,OLDF,B1,SCXBND,E,STEP,FACTOR,
     *     BRAKTD,GTEST1,GTEST2,TOL,IENTRY,ITEST)
CLSOUT
      IF (LSPRNT .GE. NPRNT) CALL LSOUT(IENTRY,ITEST,XMIN,FMIN,GMIN,
     *     XW,FW,GW,U,A,B,TOL,RELTOL,SCXBND,XBND)
C
C      IF ITEST=1, THE ALGORITHM REQUIRES THE FUNCTION VALUE TO BE
C      CALCULATED.
C



      IF (ITEST .NE. 1) GO TO 30
      UALPHA = XMIN + U
      L = LX
      DO 20 I = 1,N
         W(L) = X(I) + UALPHA*P(I)
         L = L + 1
20    CONTINUE
      CALL SFUN(N,W(LX),FU,W(LG))
      NUMF = NUMF + 1
      GU = DDOT(N,W(LG),1,P,1)
C
C      THE GRADIENT VECTOR CORRESPONDING TO THE BEST POINT IS
C      OVERWRITTEN IF FU IS LESS THAN FMIN AND FU IS SUFFICIENTLY
C      LOWER THAN F AT THE ORIGIN.
C
      IF (FU .LE. FMIN .AND. FU .LE. OLDF-UALPHA*GTEST1)
     *     CALL DCOPY(N,W(LG),1,G,1)
      GOTO 10
C
C      IF ITEST=2 OR 3 A LOWER POINT COULD NOT BE FOUND
C
30    CONTINUE
      NFTOTL = NUMF
      IFLAG = 1
      IF (ITEST .NE. 0) GO TO 50
C
C      IF ITEST=0 A SUCCESSFUL SEARCH HAS BEEN MADE
C
      IFLAG = 0
      F = FMIN
      ALPHA = XMIN
      DO 40 I = 1,N
         X(I) = X(I) + ALPHA*P(I)
40    CONTINUE
50    RETURN
      END
C
C
      SUBROUTINE GETPTC(BIG,SMALL,RTSMLL,RELTOL,ABSTOL,TNYTOL,
     *     FPRESN,ETA,RMU,XBND,U,FU,GU,XMIN,FMIN,GMIN,
     *     XW,FW,GW,A,B,OLDF,B1,SCXBND,E,STEP,FACTOR,
     *     BRAKTD,GTEST1,GTEST2,TOL,IENTRY,ITEST)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      LOGICAL BRAKTD
      INTEGER IENTRY,ITEST
      DOUBLE PRECISION BIG,SMALL,RTSMLL,RELTOL,ABSTOL,TNYTOL,
     *     FPRESN,ETA,RMU,XBND,U,FU,GU,XMIN,FMIN,GMIN,
     *     XW,FW,GW,A,B,OLDF,B1,SCXBND,E,STEP,FACTOR,
     *     GTEST1,GTEST2,TOL,DENOM
C
C ************************************************************
C GETPTC, AN ALGORITHM FOR FINDING A STEPLENGTH, CALLED REPEATEDLY BY
C ROUTINES WHICH REQUIRE A STEP LENGTH TO BE COMPUTED USING CUBIC
C INTERPOLATION. THE PARAMETERS CONTAIN INFORMATION ABOUT THE INTERVAL
C IN WHICH A LOWER POINT IS TO BE FOUND AND FROM THIS GETPTC COMPUTES A
C POINT AT WHICH THE FUNCTION CAN BE EVALUATED BY THE CALLING PROGRAM.
C THE VALUE OF THE INTEGER PARAMETERS IENTRY DETERMINES THE PATH TAKEN
C THROUGH THE CODE.
C ************************************************************
C
      LOGICAL CONVRG
      DOUBLE PRECISION ABGMIN,ABGW,ABSR,A1,CHORDM,CHORDU,
     *     D1,D2,P,Q,R,S,SCALE,SUMSQ,TWOTOL,XMIDPT
      DOUBLE PRECISION ZERO, POINT1,HALF,ONE,THREE,FIVE,ELEVEN
C
C THE FOLLOWING STANDARD FUNCTIONS AND SYSTEM FUNCTIONS ARE CALLED
C WITHIN GETPTC
C
      DOUBLE PRECISION DABS, DSQRT
C
      ZERO = 0.D0
      POINT1 = 1.D-1
      HALF = 5.D-1
      ONE = 1.D0
      THREE = 3.D0
      FIVE = 5.D0
      ELEVEN = 11.D0
C
C      BRANCH TO APPROPRIATE SECTION OF CODE DEPENDING ON THE
C      VALUE OF IENTRY.
C
      GOTO (10,20), IENTRY
C
C      IENTRY=1
C      CHECK INPUT PARAMETERS
C
10      ITEST = 2
      IF (U .LE. ZERO .OR. XBND .LE. TNYTOL .OR. GU .GT. ZERO)
     *     RETURN
      ITEST = 1
      IF (XBND .LT. ABSTOL) ABSTOL = XBND
      TOL = ABSTOL
      TWOTOL = TOL + TOL
C
C A AND B DEFINE THE INTERVAL OF UNCERTAINTY, X AND XW ARE POINTS
C WITH LOWEST AND SECOND LOWEST FUNCTION VALUES SO FAR OBTAINED.
C INITIALIZE A,SMIN,XW AT ORIGIN AND CORRESPONDING VALUES OF
C FUNCTION AND PROJECTION OF THE GRADIENT ALONG DIRECTION OF SEARCH
C AT VALUES FOR LATEST ESTIMATE AT MINIMUM.
C
      A = ZERO
      XW = ZERO
      XMIN = ZERO
      OLDF = FU
      FMIN = FU
      FW = FU
      GW = GU
      GMIN = GU
      STEP = U
      FACTOR = FIVE
C
C      THE MINIMUM HAS NOT YET BEEN BRACKETED.
C
      BRAKTD = .FALSE.
C
C SET UP XBND AS A BOUND ON THE STEP TO BE TAKEN. (XBND IS NOT COMPUTED
C EXPLICITLY BUT SCXBND IS ITS SCALED VALUE.)  SET THE UPPER BOUND
C ON THE INTERVAL OF UNCERTAINTY INITIALLY TO XBND + TOL(XBND).
C
      SCXBND = XBND
      B = SCXBND + RELTOL*DABS(SCXBND) + ABSTOL
      E = B + B
      B1 = B
C
C COMPUTE THE CONSTANTS REQUIRED FOR THE TWO CONVERGENCE CRITERIA.
C
      GTEST1 = -RMU*GU
      GTEST2 = -ETA*GU
C
C SET IENTRY TO INDICATE THAT THIS IS THE FIRST ITERATION
C
      IENTRY = 2
      GO TO 210
C
C IENTRY = 2
C
C UPDATE A,B,XW, AND XMIN
C
20      IF (FU .GT. FMIN) GO TO 60
C
C IF FUNCTION VALUE NOT INCREASED, NEW POINT BECOMES NEXT
C ORIGIN AND OTHER POINTS ARE SCALED ACCORDINGLY.
C
      CHORDU = OLDF - (XMIN + U)*GTEST1
      IF (FU .LE. CHORDU) GO TO 30
C
C THE NEW FUNCTION VALUE DOES NOT SATISFY THE SUFFICIENT DECREASE
C CRITERION. PREPARE TO MOVE THE UPPER BOUND TO THIS POINT AND
C FORCE THE INTERPOLATION SCHEME TO EITHER BISECT THE INTERVAL OF
C UNCERTAINTY OR TAKE THE LINEAR INTERPOLATION STEP WHICH ESTIMATES
C THE ROOT OF F(ALPHA)=CHORD(ALPHA).
C
      CHORDM = OLDF - XMIN*GTEST1
      GU = -GMIN
      DENOM = CHORDM-FMIN
      IF (DABS(DENOM) .GE. 1.D-15) GO TO 25
          DENOM = 1.D-15
          IF (CHORDM-FMIN .LT. 0.D0)  DENOM = -DENOM
25    CONTINUE
      IF (XMIN .NE. ZERO) GU = GMIN*(CHORDU-FU)/DENOM
      FU = HALF*U*(GMIN+GU) + FMIN
      IF (FU .LT. FMIN) FU = FMIN
      GO TO 60
30      FW = FMIN
      FMIN = FU
      GW = GMIN
      GMIN = GU
      XMIN = XMIN + U
      A = A-U
      B = B-U
      XW = -U
      SCXBND = SCXBND - U
      IF (GU .LE. ZERO) GO TO 40
      B = ZERO
      BRAKTD = .TRUE.
      GO TO 50
40    A = ZERO
50    TOL = DABS(XMIN)*RELTOL + ABSTOL
      GO TO 90
C
C IF FUNCTION VALUE INCREASED, ORIGIN REMAINS UNCHANGED
C BUT NEW POINT MAY NOW QUALIFY AS W.
C
60    IF (U .LT. ZERO) GO TO 70
      B = U
      BRAKTD = .TRUE.
      GO TO 80
70    A = U
80    XW = U
      FW = FU
      GW = GU
90    TWOTOL = TOL + TOL
      XMIDPT = HALF*(A + B)
C
C CHECK TERMINATION CRITERIA
C
      CONVRG = DABS(XMIDPT) .LE. TWOTOL - HALF*(B-A) .OR.
     *     DABS(GMIN) .LE. GTEST2 .AND. FMIN .LT. OLDF .AND.
     *     (DABS(XMIN - XBND) .GT. TOL .OR. .NOT. BRAKTD)
      IF (.NOT. CONVRG) GO TO 100
      ITEST = 0
      IF (XMIN .NE. ZERO) RETURN
C
C IF THE FUNCTION HAS NOT BEEN REDUCED, CHECK TO SEE THAT THE RELATIVE
C CHANGE IN F(X) IS CONSISTENT WITH THE ESTIMATE OF THE DELTA-
C UNIMODALITY CONSTANT, TOL.  IF THE CHANGE IN F(X) IS LARGER THAN
C EXPECTED, REDUCE THE VALUE OF TOL.
C
      ITEST = 3
      IF (DABS(OLDF-FW) .LE. FPRESN*(ONE + DABS(OLDF))) RETURN
      TOL = POINT1*TOL
      IF (TOL .LT. TNYTOL) RETURN
      RELTOL = POINT1*RELTOL
      ABSTOL = POINT1*ABSTOL
      TWOTOL = POINT1*TWOTOL
C
C CONTINUE WITH THE COMPUTATION OF A TRIAL STEP LENGTH
C
100   R = ZERO
      Q = ZERO
      S = ZERO
      IF (DABS(E) .LE. TOL) GO TO 150
C
C FIT CUBIC THROUGH XMIN AND XW
C
      R = THREE*(FMIN-FW)/XW + GMIN + GW
      ABSR = DABS(R)
      Q = ABSR
      IF (GW .EQ. ZERO .OR. GMIN .EQ. ZERO) GO TO 140
C
C COMPUTE THE SQUARE ROOT OF (R*R - GMIN*GW) IN A WAY
C WHICH AVOIDS UNDERFLOW AND OVERFLOW.
C
      ABGW = DABS(GW)
      ABGMIN = DABS(GMIN)
      S = DSQRT(ABGMIN)*DSQRT(ABGW)
      IF ((GW/ABGW)*GMIN .GT. ZERO) GO TO 130
C
C COMPUTE THE SQUARE ROOT OF R*R + S*S.
C
      SUMSQ = ONE
      P = ZERO
      IF (ABSR .GE. S) GO TO 110
C
C THERE IS A POSSIBILITY OF OVERFLOW.
C
      IF (S .GT. RTSMLL) P = S*RTSMLL
      IF (ABSR .GE. P) SUMSQ = ONE +(ABSR/S)**2
      SCALE = S
      GO TO 120
C
C THERE IS A POSSIBILITY OF UNDERFLOW.
C
110   IF (ABSR .GT. RTSMLL) P = ABSR*RTSMLL
      IF (S .GE. P) SUMSQ = ONE + (S/ABSR)**2
      SCALE = ABSR
120   SUMSQ = DSQRT(SUMSQ)
      Q = BIG
      IF (SCALE .LT. BIG/SUMSQ) Q = SCALE*SUMSQ
      GO TO 140
C
C COMPUTE THE SQUARE ROOT OF R*R - S*S
C
130   Q = DSQRT(DABS(R+S))*DSQRT(DABS(R-S))
      IF (R .GE. S .OR. R .LE. (-S)) GO TO 140
      R = ZERO
      Q = ZERO
      GO TO 150
C
C COMPUTE THE MINIMUM OF FITTED CUBIC
C
140   IF (XW .LT. ZERO) Q = -Q
      S = XW*(GMIN - R - Q)
      Q = GW - GMIN + Q + Q
      IF (Q .GT. ZERO) S = -S
      IF (Q .LE. ZERO) Q = -Q
      R = E
      IF (B1 .NE. STEP .OR. BRAKTD) E = STEP
C
C CONSTRUCT AN ARTIFICIAL BOUND ON THE ESTIMATED STEPLENGTH
C
150   A1 = A
      B1 = B
      STEP = XMIDPT
      IF (BRAKTD) GO TO 160
      STEP = -FACTOR*XW
      IF (STEP .GT. SCXBND) STEP = SCXBND
      IF (STEP .NE. SCXBND) FACTOR = FIVE*FACTOR
      GO TO 170
C
C IF THE MINIMUM IS BRACKETED BY 0 AND XW THE STEP MUST LIE
C WITHIN (A,B).
C
160   IF ((A .NE. ZERO .OR. XW .GE. ZERO) .AND. (B .NE. ZERO .OR.
     *     XW .LE. ZERO)) GO TO 180
C
C IF THE MINIMUM IS NOT BRACKETED BY 0 AND XW THE STEP MUST LIE
C WITHIN (A1,B1).
C
      D1 = XW
      D2 = A
      IF (A .EQ. ZERO) D2 = B
C THIS LINE MIGHT BE
C     IF (A .EQ. ZERO) D2 = E
      U = - D1/D2
      STEP = FIVE*D2*(POINT1 + ONE/U)/ELEVEN
      IF (U .LT. ONE) STEP = HALF*D2*DSQRT(U)
170   IF (STEP .LE. ZERO) A1 = STEP
      IF (STEP .GT. ZERO) B1 = STEP
C
C REJECT THE STEP OBTAINED BY INTERPOLATION IF IT LIES OUTSIDE THE
C REQUIRED INTERVAL OR IT IS GREATER THAN HALF THE STEP OBTAINED
C DURING THE LAST-BUT-ONE ITERATION.
C
180   IF (DABS(S) .LE. DABS(HALF*Q*R) .OR.
     *     S .LE. Q*A1 .OR. S .GE. Q*B1) GO TO 200
C
C A CUBIC INTERPOLATION STEP
C
      STEP = S/Q
C
C THE FUNCTION MUST NOT BE EVALUTATED TOO CLOSE TO A OR B.
C
      IF (STEP - A .GE. TWOTOL .AND. B - STEP .GE. TWOTOL) GO TO 210
      IF (XMIDPT .GT. ZERO) GO TO 190
      STEP = -TOL
      GO TO 210
190   STEP = TOL
      GO TO 210
200   E = B-A
C
C IF THE STEP IS TOO LARGE, REPLACE BY THE SCALED BOUND (SO AS TO
C COMPUTE THE NEW POINT ON THE BOUNDARY).
C
210   IF (STEP .LT. SCXBND) GO TO 220
      STEP = SCXBND
C
C MOVE SXBD TO THE LEFT SO THAT SBND + TOL(XBND) = XBND.
C
      SCXBND = SCXBND - (RELTOL*DABS(XBND)+ABSTOL)/(ONE + RELTOL)
220   U = STEP
      IF (DABS(STEP) .LT. TOL .AND. STEP .LT. ZERO) U = -TOL
      IF (DABS(STEP) .LT. TOL .AND. STEP .GE. ZERO) U = TOL
      ITEST = 1
      RETURN
      END
                                                                                                                                                                           r8s/tnc.c                                                                                           0000644 0000766 0000120 00000137446 10357567246 012231  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /* tnc : truncated newton bound contrained minimization
         using gradient information, in C */

/*
 * Copyright (c) 2002-2004, Jean-Sebastien Roy (js@jeannot.org)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * This software is a C implementation of TNBC, a truncated newton minimization
 * package originally developed by Stephen G. Nash in Fortran.
 * 
 * The original source code can be found at :
 * http://iris.gmu.edu/~snash/nash/software/software.html
 * 
 * Copyright for the original TNBC fortran routines:
 * 
 *   TRUNCATED-NEWTON METHOD:  SUBROUTINES
 *     WRITTEN BY:  STEPHEN G. NASH
 *           SCHOOL OF INFORMATION TECHNOLOGY & ENGINEERING
 *           GEORGE MASON UNIVERSITY
 *           FAIRFAX, VA 22030
 */

/*
 * Conversion into C by Elisabeth Nguyen & Jean-Sebastien Roy
 * Modifications by Jean-Sebastien Roy, 2001-2002
 */

static char const rcsid[] =
  "@(#) $Jeannot: tnc.c,v 1.202 2004/04/18 10:32:30 js Exp $";

static char const copyright[] =
  "(c) 2002-2003, Jean-Sebastien Roy (js@jeannot.org)";

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "tnc.h"

typedef enum
{
  TNC_FALSE = 0,
  TNC_TRUE
} logical;

/*
 * Return code strings
 */

char *tnc_rc_string[10] =
{
  "Memory allocation failed",
  "Invalid parameters (n<0)",
  "Infeasible (low bound > up bound)",
  "Local minima reach (|pg| ~= 0)",
  "Converged (|f_n-f_(n-1)| ~= 0)",
  "Maximum number of function evaluations reached",
  "Linear search failed",
  "All lower bounds are equal to the upper bounds",
  "Unable to progress",
  "User requested end of minimization"
};

/*
 * getptc return codes
 */
typedef enum
{
  GETPTC_OK     = 0, /* Suitable point found */
  GETPTC_EVAL   = 1, /* Function evaluation required */
  GETPTC_EINVAL = 2, /* Bad input values */
  GETPTC_FAIL   = 3  /* No suitable point found */
} getptc_rc;

/*
 * linearSearch return codes
 */
typedef enum
{
  LS_OK        = 0, /* Suitable point found */
  LS_MAXFUN    = 1, /* Max. number of function evaluations reach */
  LS_FAIL      = 2, /* No suitable point found */
  LS_USERABORT = 3, /* User requested end of minimization */
  LS_ENOMEM    = 4  /* Memory allocation failed */
} ls_rc;

/*
 * Prototypes
 */
static tnc_rc tnc_minimize(int n, double x[], double *f, double g[],
  tnc_function *function, void *state, 
  double xscale[], double *fscale,
  double low[], double up[], tnc_message messages,
  int maxCGit, int maxnfeval, int *nfeval,
  double eta, double stepmx, double accuracy,
  double fmin, double ftol, double rescale);

static getptc_rc getptcInit(double *reltol, double *abstol, double tnytol, 
  double eta, double rmu, double xbnd, 
  double *u, double *fu, double *gu, double *xmin,
  double *fmin, double *gmin, double *xw, double *fw, 
  double *gw, double *a, double *b, double *oldf, 
  double *b1, double *scxbnd, double *e, double *step, 
  double *factor, logical *braktd, double *gtest1, 
  double *gtest2, double *tol);

static getptc_rc getptcIter(double big, double 
  rtsmll, double *reltol, double *abstol, double tnytol, 
  double fpresn, double xbnd, 
  double *u, double *fu, double *gu, double *xmin,
  double *fmin, double *gmin, double *xw, double *fw, 
  double *gw, double *a, double *b, double *oldf, 
  double *b1, double *scxbnd, double *e, double *step, 
  double *factor, logical *braktd, double *gtest1, 
  double *gtest2, double *tol);

static void printCurrentIteration(int n, double f, double g[], int niter,
  int nfeval, int pivot[]);

static double initialStep(double fnew, double fmin, double gtp, double smax);

static ls_rc linearSearch(int n, tnc_function *function, void *state,
  double low[], double up[],
  double xscale[], double fscale, int pivot[],
  double eta, double ftol, double xbnd, 
  double p[], double x[], double *f, 
  double *alpha, double gfull[], int maxnfeval, int *nfeval);

static int tnc_direction(double *zsol, double *diagb,
  double *x, double *g, int n,
  int maxCGit, int maxnfeval, int *nfeval, 
  logical upd1, double yksk, double yrsr,
  double *sk, double *yk, double *sr, double *yr,
  logical lreset, tnc_function *function, void *state, 
  double xscale[], double fscale,
  int *pivot, double accuracy,
  double gnorm, double xnorm, double *low, double *up);

static double stepMax(double step, int n, double x[], double p[], int pivot[], 
  double low[], double up[], double xscale[]);

/* Active set of constraints */
static void setContraints(int n, double x[], int pivot[], double xscale[],
  double low[], double up[]);

static logical addConstraint(int n, double x[], double p[], int pivot[],
  double low[], double up[], double xscale[]);

static logical removeConstraint(double gtpnew, double f, 
  double *fLastConstraint, double g[], int pivot[], int n);

static void project(int n, double x[], int pivot[]);

static int hessianTimesVector(double v[], double gv[], int n, 
  double x[], double g[], tnc_function *function, void *state, 
  double xscale[], double fscale,
  double accuracy, double xnorm, double low[], double up[]);

static int msolve(double g[], double *y, int n, 
  double sk[], double yk[], double diagb[], double sr[], 
  double yr[], logical upd1, double yksk, double yrsr, 
  logical lreset);

static void diagonalScaling(int n, double e[], double v[], double gv[],
  double r[]);

static void ssbfgs(int n, double gamma, double sj[], double *hjv,
  double hjyj[], double yjsj, 
  double yjhyj, double vsj, double vhyj, double hjp1v[]);

static int initPreconditioner(double diagb[], double emat[], int n, 
  logical lreset, double yksk, double yrsr,
  double sk[], double yk[], double sr[], double yr[], 
  logical upd1);

/* Scaling */
static void coercex(int n, double x[], double low[], double up[]);
static void unscalex(int n, double x[], double xscale[]);
static void scaleg(int n, double g[], double xscale[], double fscale);
static void scalex(int n, double x[], double xscale[]);
static void projectConstants(int n, double x[], double xscale[]);

/* Machine precision */
static double mchpr1(void);

/* Special blas for incx=incy=1 */
static double ddot1(int n, double dx[], double dy[]);
static void dxpy1(int n, double dx[], double dy[]);
static void daxpy1(int n, double da, double dx[], double dy[]);
static void dcopy1(int n, double dx[], double dy[]);
static double dnrm21(int n, double dx[]);

/* additionnal blas-like functions */
static void dneg1(int n, double v[]);
static double dnrmi1(int n, double v[]);

/*
 * This routine solves the optimization problem
 * 
 *   minimize   f(x)
 *     x
 *   subject to   low <= x <= up
 * 
 * where x is a vector of n real variables. The method used is
 * a truncated-newton algorithm (see "newton-type minimization via
 * the lanczos algorithm" by s.g. nash (technical report 378, math.
 * the lanczos method" by s.g. nash (siam j. numer. anal. 21 (1984),
 * pp. 770-778).  this algorithm finds a local minimum of f(x). It does
 * not assume that the function f is convex (and so cannot guarantee a
 * global solution), but does assume that the function is bounded below.
 * it can solve problems having any number of variables, but it is
 * especially useful when the number of variables (n) is large.
 * 
 */
extern int tnc(int n, double x[], double *f, double g[], tnc_function *function,
  void *state, double low[], double up[], double scale[], int messages, 
  int maxCGit, int maxnfeval, double eta, double stepmx, 
  double accuracy, double fmin, double ftol, double rescale, int *nfeval)
{
  int rc, frc, i, nc, nfeval_local,
    free_low = TNC_FALSE, free_up = TNC_FALSE,
    free_g = TNC_FALSE;
  double *xscale = NULL, fscale, epsmch, rteps;

  if(nfeval==NULL)
  {
    /* Ignore nfeval */
    nfeval = &nfeval_local;
  }
  *nfeval = 0;
  
  /* Version info */
  if (messages & TNC_MSG_VERS)
  {
    fprintf(stderr, "tnc: Version %s, %s\n",TNC_VERSION,copyright);
    fprintf(stderr, "tnc: RCS ID: %s\n",rcsid);
  }

  /* Check for errors in the input parameters */
  if (n == 0)
  {
    rc = TNC_CONSTANT;
    goto cleanup;
  }
  
  if (n < 0)
  {
    rc = TNC_EINVAL;
    goto cleanup;
  }
  
  /* Check bounds arrays */
  if (low == NULL)
  {
    low = malloc(n*sizeof(*low));
    if (low == NULL)
    {
      rc = TNC_ENOMEM;
      goto cleanup;
    }
    free_low = TNC_TRUE;
    for (i = 0 ; i < n ; i++) low[i] = -HUGE_VAL;
  }
  if (up == NULL)
  {
    up = malloc(n*sizeof(*up));
    if (up == NULL)
    {
      rc = TNC_ENOMEM;
      goto cleanup;
    }
    free_up = TNC_TRUE;
    for (i = 0 ; i < n ; i++) up[i] = HUGE_VAL;
  }

  /* Coherency check */
  for (i = 0 ; i < n ; i++)
  {
    if (low[i] > up [i])
    {
      rc = TNC_INFEASIBLE;
      goto cleanup;
    }
  }

  /* Coerce x into bounds */
  coercex(n, x, low, up);

  if (maxnfeval < 1)
  {
    rc = TNC_MAXFUN;
    goto cleanup;
  }
  
  /* Allocate g if necessary */
  if(g == NULL)
  {
    g = malloc(n*sizeof(*g));
    if (g == NULL)
    {
      rc = TNC_ENOMEM;
      goto cleanup;
    }
    free_g = TNC_TRUE;
  }

  /* Initial function evaluation */  
  frc = function(x, f, g, state);
  (*nfeval) ++;
  if (frc)
  {
    rc = TNC_USERABORT;
    goto cleanup;
  }
    
  /* Constant problem ? */
  for (nc = 0, i = 0 ; i < n ; i++)
    if ((low[i] == up[i]) || (scale != NULL && scale[i] == 0.0))
      nc ++;

  if (nc == n)
  {
    rc = TNC_CONSTANT;
    goto cleanup;
  }

  /* Scaling parameters */  
  xscale = malloc(sizeof(*xscale)*n);
  if (xscale == NULL)
  {
    rc = TNC_ENOMEM;
    goto cleanup;
  }
  fscale = 1.0;

  for (i = 0 ; i < n ; i++)
  {
    if (scale != NULL)
    {
      xscale[i] = fabs(scale[i]);
      if (xscale[i] == 0.0)
        low[i] = up[i] = x[i];
    }
    else if (low[i] != -HUGE_VAL && up[i] != HUGE_VAL)
      xscale[i] = up[i] - low[i];
    else
      xscale[i] = 1.0+fabs(x[i]);
  }

  /* Default values for parameters */
  epsmch = mchpr1();
  rteps = sqrt(epsmch);

  if (stepmx < rteps * 10.0) stepmx = 1.0e1;
  if (eta < 0.0 || eta >= 1.0) eta = 0.25;
  if (rescale < 0) rescale = 1.3;
  if (maxCGit < 0) /* maxCGit == 0 is valid */
  {
    maxCGit = n / 2;
    if (maxCGit < 1) maxCGit = 1;
    else if (maxCGit > 50) maxCGit = 50;
  }
  if (maxCGit > n) maxCGit = n;
  if (ftol < 0.0) ftol = 0.0;
  if (accuracy <= epsmch) accuracy = rteps;

  /* Optimisation */
  rc = tnc_minimize(n, x, f, g, function, state,
    xscale, &fscale, low, up, messages, 
    maxCGit, maxnfeval, nfeval, eta, stepmx, accuracy, fmin, ftol, rescale);

cleanup:
  if (messages & TNC_MSG_EXIT)
    fprintf(stderr, "tnc: %s\n", tnc_rc_string[rc - TNC_MINRC]);

  if (xscale) free(xscale);
  if (free_low) free(low);
  if (free_up) free(up);
  if (free_g) free(g);
  
  return rc;
}

/* Coerce x into bounds */
static void coercex(int n, double x[], double low[], double up[])
{
  int i;
  
  for (i = 0 ; i < n ; i++)
  {
    if (x[i]<low[i]) x[i] = low[i];
    else if (x[i]>up[i]) x[i] = up[i];
  }
}

/* Unscale x */
static void unscalex(int n, double x[], double xscale[])
{
  int i;
  for (i = 0 ; i < n ; i++)
    x[i] *= xscale[i];
}

/* Scale x */
static void scalex(int n, double x[], double xscale[])
{
  int i;
  for (i = 0 ; i < n ; i++)
    if (xscale[i]>0.0)
      x[i] /= xscale[i];
}

/* Scale g */
static void scaleg(int n, double g[], double xscale[], double fscale)
{
  int i;
  for (i = 0 ; i < n ; i++)
    g[i] *= xscale[i]*fscale;
}

/* Caculate the pivot vector */
static void setContraints(int n, double x[], int pivot[], double xscale[],
  double low[], double up[])
{
  int i;
  double epsmch;
  
  epsmch = mchpr1();

  for (i = 0; i < n; i++)
  {
    double tol;

    /* tolerances should be better ajusted */
    if (xscale[i] == 0.0)
    {
      pivot[i] = 2;
    }
    else
    {
       tol = epsmch * 10.0 * (fabs(low[i]) + 1.0);
      if ((x[i]*xscale[i] - low[i] <= tol) && low[i] != - HUGE_VAL)
        pivot[i] = -1;
      else
      {
         tol = epsmch * 10.0 * (fabs(up[i]) + 1.0);
        if ((x[i]*xscale[i] - up[i] >= tol) && up[i] != HUGE_VAL)
          pivot[i] = 1;
        else
          pivot[i] = 0;
      }
    }
  }
}

/*
 * This routine is a bounds-constrained truncated-newton method.
 * the truncated-newton method is preconditioned by a limited-memory
 * quasi-newton method (this preconditioning strategy is developed
 * in this routine) with a further diagonal scaling
 * (see routine diagonalscaling).
 */
static tnc_rc tnc_minimize(int n, double x[], 
  double *f, double gfull[], tnc_function *function, void *state, 
  double xscale[], double *fscale,
  double low[], double up[], tnc_message messages, 
  int maxCGit, int maxnfeval, int *nfeval, double eta, double stepmx, 
  double accuracy, double fmin, double ftol, double rescale)
{
  double fLastReset, difnew, epsmch, epsred, oldgtp,
    difold, oldf, rteps, xnorm, newscale,
    gnorm, ustpmax, fLastConstraint, spe, yrsr, yksk,
    *temp = NULL, *sk = NULL, *yk = NULL, *diagb = NULL, *sr = NULL,
    *yr = NULL, *oldg = NULL, *pk = NULL, *g = NULL;
  double alpha = 0.0; /* Default unused value */
  int i, icycle, niter = 0, oldnfeval, *pivot = NULL, frc;
  logical lreset, newcon, upd1, remcon;
  tnc_rc rc = TNC_ENOMEM; /* Default error */

  /* Allocate temporary vectors */  
  oldg = malloc(sizeof(*oldg)*n);
   if (oldg == NULL) goto cleanup;
  g = malloc(sizeof(*g)*n);
   if (g == NULL) goto cleanup;
  temp = malloc(sizeof(*temp)*n);
   if (temp == NULL) goto cleanup;
  diagb = malloc(sizeof(*diagb)*n);
   if (diagb == NULL) goto cleanup;
  pk = malloc(sizeof(*pk)*n);
   if (pk == NULL) goto cleanup;

  sk = malloc(sizeof(*sk)*n);
   if (sk == NULL) goto cleanup;
  yk = malloc(sizeof(*yk)*n);
   if (yk == NULL) goto cleanup;
  sr = malloc(sizeof(*sr)*n);
   if (sr == NULL) goto cleanup;
  yr = malloc(sizeof(*yr)*n);
   if (yr == NULL) goto cleanup;

  pivot = malloc(sizeof(*pivot)*n);
   if (pivot == NULL) goto cleanup;

  /* Initialize variables */
  epsmch = mchpr1();
  rteps = sqrt(epsmch);

  difnew = 0.0;
  epsred = 0.05;
  upd1 = TNC_TRUE;
  icycle = n - 1;
  newcon = TNC_TRUE;
  
  /* Uneeded initialisations */
  lreset = TNC_FALSE;
  yrsr = 0.0;
  yksk = 0.0;
  
  /* Initial scaling */
  scalex(n, x, xscale);
  (*f) *= *fscale;

  /* initial pivot calculation */
  setContraints(n, x, pivot, xscale, low, up);

  dcopy1(n, gfull, g);
  scaleg(n, g, xscale, *fscale);

  /* Test the lagrange multipliers to see if they are non-negative. */
  for (i = 0; i < n; i++)
    if (-pivot[i] * g[i] < 0.0)
      pivot[i] = 0;

  project(n, g, pivot);

  /* Set initial values to other parameters */
  gnorm = dnrm21(n, g);

  fLastConstraint = *f; /* Value at last constraint */
  fLastReset = *f; /* Value at last reset */

  if (messages & TNC_MSG_ITER) fprintf(stderr,
    "  NIT   NF   F                       GTG\n");
  if (messages & TNC_MSG_ITER) printCurrentIteration(n, *f / *fscale, gfull,
    niter, *nfeval, pivot);

  /* Set the diagonal of the approximate hessian to unity. */
  for (i = 0; i < n; i++) diagb[i] = 1.0;

  /* Start of main iterative loop */
  while(TNC_TRUE)
  {
    /* Tolerance should be user modifiable */
    if (dnrmi1(n, g) <= 1.0e-2*rteps*fabs(*f))
    {
      /* |PG| == 0.0 => local minimum */
      dcopy1(n, gfull, g);
      project(n, g, pivot);
      if (messages & TNC_MSG_INFO) fprintf(stderr,
        "tnc: |pg| = %g -> local minimum\n",dnrmi1(n, g));
      rc = TNC_LOCALMINIMUM;
      break;
    }

    /* Terminate if more than maxnfeval evaluations have been made */
    if (*nfeval >= maxnfeval)
    {
      rc = TNC_MAXFUN;
      break;
    }

    /* Rescale function if necessary */
    newscale = dnrmi1(n, g);
    if ((newscale > epsmch) && (fabs(log10(newscale)) > rescale))
    {
      newscale = 1.0/newscale;
      
      *f *= newscale;
      *fscale *= newscale;
      gnorm *= newscale;
      fLastConstraint *= newscale;
      fLastReset *= newscale;
      difnew *= newscale;

      for (i = 0; i < n; i++) g[i] *= newscale;
      for (i = 0; i < n; i++) diagb[i] = 1.0;

      upd1 = TNC_TRUE;
      icycle = n - 1;
      newcon = TNC_TRUE;

      if (messages & TNC_MSG_INFO) fprintf(stderr, 
        "tnc: fscale = %g\n", *fscale);
    }

    dcopy1(n, x, temp);
    project(n, temp, pivot);
    xnorm = dnrm21(n, temp);
    oldnfeval = *nfeval;

    /* Compute the new search direction */
    frc = tnc_direction(pk, diagb, x, g, n, maxCGit, maxnfeval, nfeval,
      upd1, yksk, yrsr, sk, yk, sr, yr,
      lreset, function, state, xscale, *fscale,
      pivot, accuracy, gnorm, xnorm, low, up);

    if (frc == -1)
    {
      rc = TNC_ENOMEM;
      break;
    }

    if (frc)
    {
      rc = TNC_USERABORT;
      break;
    }

    if (!newcon)
    {
      if (!lreset)
      {
        /* Compute the accumulated step and its corresponding gradient
          difference. */
        dxpy1(n, sk, sr);
        dxpy1(n, yk, yr);
        icycle++;
      }
      else
      {
        /* Initialize the sum of all the changes */
        dcopy1(n, sk, sr);
        dcopy1(n, yk, yr);
        fLastReset = *f;
        icycle = 1;
      }
    }

    dcopy1(n, g, oldg);
    oldf = *f;
    oldgtp = ddot1(n, pk, g);

    /* Maximum unconstrained step length */
    ustpmax = stepmx / (dnrm21(n, pk) + epsmch);

    /* Maximum constrained step length */
    spe = stepMax(ustpmax, n, x, pk, pivot, low, up, xscale);
    
    if (spe > 0.0)
    {
      ls_rc lsrc;
      /* Set the initial step length */
      alpha = initialStep(*f, fmin / (*fscale), oldgtp, spe);

      /* Perform the linear search */
      lsrc = linearSearch(n, function, state, low, up,
        xscale, *fscale, pivot,
        eta, ftol, spe, pk, x, f, &alpha, gfull, maxnfeval, nfeval);

      if (lsrc == LS_ENOMEM)
      {
        rc = TNC_ENOMEM;
        break;
      }

      if (lsrc == LS_USERABORT)
      {
        rc = TNC_USERABORT;
        break;
      }

      /* If we went up to the maximum unconstrained step, increase it */
      if (alpha >= 0.9 * ustpmax)
      {
        stepmx *= 1e2;
        if (messages & TNC_MSG_INFO) fprintf(stderr,
          "tnc: stepmx = %g\n", stepmx);
      }

      /* If we went up to the maximum constrained step,
         a new constraint was encountered */
      if (alpha - spe >= -epsmch * 10.0)
      {
        newcon = TNC_TRUE;
      }
      else
      {
        /* Break if the linear search has failed to find a lower point */
        if (lsrc != LS_OK)
        {
          if (lsrc == LS_MAXFUN) rc = TNC_MAXFUN;
          else rc = TNC_LSFAIL;
          break;
        }
        newcon = TNC_FALSE;
      }
    }
    else
    {
      /* Maximum constrained step == 0.0 => new constraint */
      newcon = TNC_TRUE;
    }

    if (newcon)
    {
      if(!addConstraint(n, x, pk, pivot, low, up, xscale))
      {
        if(*nfeval == oldnfeval)
        {
          rc = TNC_NOPROGRESS;
          break;
        }
      }
      fLastConstraint = *f;
    }

    niter++;

    /* Set up parameters used in convergence and resetting tests */
    difold = difnew;
    difnew = oldf - *f;

    /* If this is the first iteration of a new cycle, compute the
       percentage reduction factor for the resetting test */
    if (icycle == 1)
    {
      if (difnew > difold * 2.0) epsred += epsred;
      if (difnew < difold * 0.5) epsred *= 0.5;
    }

    dcopy1(n, gfull, g);
    scaleg(n, g, xscale, *fscale);

    dcopy1(n, g, temp);
    project(n, temp, pivot);
    gnorm = dnrm21(n, temp);

    /* Reset pivot */
    remcon = removeConstraint(oldgtp, *f, &fLastConstraint, g, pivot, n);

    if (!remcon && !newcon)
    {
      /* No constraint removed & no new constraint : test for convergence */
      if (fabs(difnew) <= ftol*epsmch*0.5*(fabs(oldf)+fabs(*f)))
      {
        if (messages & TNC_MSG_INFO) fprintf(stderr, 
          "tnc: |fn-fn-1] = %g -> convergence\n",fabs(difnew));
        rc = TNC_CONVERGED;
        break;
      }
    }

    project(n, g, pivot);

    if (messages & TNC_MSG_ITER) printCurrentIteration(n, *f / *fscale, gfull,
      niter, *nfeval, pivot);

    /* Compute the change in the iterates and the corresponding change in the
      gradients */
    if (!newcon)
    {
      for (i = 0; i < n; i++)
      {
        yk[i] = g[i] - oldg[i];
        sk[i] = alpha * pk[i];
      }

      /* Set up parameters used in updating the preconditioning strategy */
      yksk = ddot1(n, yk, sk);
      
      if (icycle == (n - 1) || difnew < epsred * (fLastReset - *f))
        lreset = TNC_TRUE;
      else
      {
        yrsr = ddot1(n, yr, sr);
        if (yrsr <= 0.0) lreset = TNC_TRUE;
        else lreset = TNC_FALSE;
      }
      upd1 = TNC_FALSE;
    }
  }

  if (messages & TNC_MSG_ITER) printCurrentIteration(n, *f / *fscale, gfull,
    niter, *nfeval, pivot);

  /* Unscaling */
  unscalex(n, x, xscale);
  coercex(n, x, low, up);
  (*f) /= *fscale;

cleanup: 
  if (oldg) free(oldg);
  if (g) free(g);
  if (temp) free(temp);
  if (diagb) free(diagb);
  if (pk) free(pk);

  if (sk) free(sk);
  if (yk) free(yk);
  if (sr) free(sr);
  if (yr) free(yr);
  
  if (pivot) free(pivot);

  return rc;
}

/* Print the results of the current iteration */
static void printCurrentIteration(int n, double f, double g[], int niter,
  int nfeval, int pivot[])
{
  int i;
  double gtg;

  gtg = 0.0;
  for (i = 0; i < n; i++)
    if (pivot[i] == 0)
      gtg += g[i] * g[i];

  fprintf(stderr, " %4d %4d %22.15E  %15.8E\n", niter, nfeval, f, gtg);
}

/*
 * Set x[i] = 0.0 if direction i is currently constrained 
 */
static void project(int n, double x[], int pivot[])
{
  int i;
  for (i = 0; i < n; i++)
    if (pivot[i] != 0)
      x[i] = 0.0;
}

/*
 * Set x[i] = 0.0 if direction i is constant
 */
static void projectConstants(int n, double x[], double xscale[])
{
  int i;
  for (i = 0; i < n; i++)
    if (xscale[i] == 0.0)
      x[i] = 0.0;
}

/*
 * Compute the maximum allowable step length
 */
static double stepMax(double step, int n, double x[], double dir[],
  int pivot[], double low[], double up[], double xscale[])
{
  int i;
  double t;

  /* Constrained maximum step */
  for (i = 0; i < n; i++)
  {
    if ((pivot[i] == 0) && (dir[i] != 0.0))
    {
      if (dir[i] < 0.0)
      {
        t = low[i]/xscale[i] - x[i];
        if (t > step * dir[i]) step = t / dir[i];
      }
      else
      {
        t = up[i]/xscale[i] - x[i];
        if (t < step * dir[i]) step = t / dir[i];
      }
    }
  }
  
  return step;
}

/*
 * Update the constraint vector pivot if a new constraint is encountered
 */
static logical addConstraint(int n, double x[], double p[], int pivot[],
  double low[], double up[], double xscale[])
{
  int i, newcon = TNC_FALSE;
  double tol, epsmch;

  epsmch = mchpr1();

  for (i = 0; i < n; i++)
  {
    if ((pivot[i] == 0) && (p[i] != 0.0))
    {
       if (p[i] < 0.0 && low[i] != - HUGE_VAL)
      {
         tol = epsmch * 10.0 * (fabs(low[i]) + 1.0);
        if (x[i]*xscale[i] - low[i] <= tol)
        {
          pivot[i] = -1;
          x[i] = low[i]/xscale[i];
          newcon = TNC_TRUE;
        }
      }
      else if (up[i] != HUGE_VAL)
      {
        tol = epsmch * 10.0 * (fabs(up[i]) + 1.0);
        if (up[i] - x[i]*xscale[i] <= tol)
        {
          pivot[i] = 1;
          x[i] = up[i]/xscale[i];
          newcon = TNC_TRUE;
        }
      }
    }
  }
  return newcon;
}

/*
 * Check if a constraint is no more active
 */
static logical removeConstraint(double gtpnew, double f, 
  double *fLastConstraint, double g[], int pivot[], int n)
{
  double cmax, t;
  int imax, i;
  logical ltest;

  imax = -1;
  cmax = 0.0;
  ltest = (*fLastConstraint - f) <= (gtpnew * -0.5);
  for (i = 0; i < n; i++)
  {
    if (pivot[i] != 2)
    {
      t = -pivot[i] * g[i];
      if (t < 0.0)
      {
        if ((!ltest) && (cmax > t))
        {
          cmax = t;
          imax = i;
        }
      }
    }
  }

  if (imax != -1)
  {
    pivot[imax] = 0;
    *fLastConstraint = f;
    return TNC_TRUE;
  }
  else
    return TNC_FALSE;

/*
 * For details, see gill, murray, and wright (1981, p. 308) and
 * fletcher (1981, p. 116). The multiplier tests (here, testing
 * the sign of the components of the gradient) may still need to
 * modified to incorporate tolerances for zero.
 */
}

/*
 * This routine performs a preconditioned conjugate-gradient
 * iteration in order to solve the newton equations for a search
 * direction for a truncated-newton algorithm.
 * When the value of the quadratic model is sufficiently reduced,
 * the iteration is terminated.
 */
static int tnc_direction(double *zsol, double *diagb,
  double *x, double g[], int n,
  int maxCGit, int maxnfeval, int *nfeval, 
  logical upd1, double yksk, double yrsr,
  double *sk, double *yk, double *sr, double *yr,
  logical lreset, tnc_function *function, void *state, 
  double xscale[], double fscale,
  int *pivot, double accuracy,
  double gnorm, double xnorm, double low[], double up[])
{
  double alpha, beta, qold, qnew, rhsnrm, tol, vgv, rz, rzold, qtest, pr, gtp;
  int i, k, frc;
  /* Temporary vectors */
  double *r = NULL, *zk = NULL, *v = NULL, *emat = NULL, *gv = NULL;

  /* No CG it. => dir = -grad */
  if (maxCGit == 0)
  {
    dcopy1(n, g, zsol);
    dneg1(n, zsol);
    project(n, zsol, pivot);
    return 0;
  }

  /* General initialization */
  rhsnrm = gnorm;
  tol = 1e-12;
  qold = 0.0;
  rzold = 0.0; /* Uneeded */

  frc = -1; /* ENOMEM here */
  r = malloc(sizeof(*r)*n); /* Residual */
  if (r == NULL) goto cleanup;
  v = malloc(sizeof(*v)*n);
  if (v == NULL) goto cleanup;
  zk = malloc(sizeof(*zk)*n);
  if (zk == NULL) goto cleanup;
  emat = malloc(sizeof(*emat)*n); /* Diagonal preconditoning matrix */
  if (emat == NULL) goto cleanup;
  gv = malloc(sizeof(*gv)*n); /* hessian times v */
  if (gv == NULL) goto cleanup;

  /* Initialization for preconditioned conjugate-gradient algorithm */
  frc = initPreconditioner(diagb, emat, n, lreset, yksk, yrsr, sk, yk, sr, yr,
    upd1);
  if (frc) goto cleanup;

  for (i = 0; i < n; i++)
  {
    r[i] = -g[i];
    v[i] = 0.0;
    zsol[i] = 0.0; /* Computed search direction */
  }

  /* Main iteration */
  for (k = 0; k < maxCGit; k++)
  {
    /* CG iteration to solve system of equations */
    project(n, r, pivot);
    frc = msolve(r, zk, n, sk, yk, diagb, sr, yr, upd1, yksk, yrsr, lreset);
    if (frc) goto cleanup;
    project(n, zk, pivot);
    rz = ddot1(n, r, zk);

    if ((rz / rhsnrm < tol) || ((*nfeval) >= (maxnfeval-1)))
    {
      /* Truncate algorithm in case of an emergency
         or too many function evaluations */
      if (k == 0)
      {
        dcopy1(n, g, zsol);
        dneg1(n, zsol);
        project(n, zsol, pivot);
      }
      break;
    }
    if (k == 0) beta = 0.0;
    else beta = rz / rzold;

    for (i = 0; i < n; i++)
      v[i] = zk[i] + beta * v[i];

    project(n, v, pivot);
    frc = hessianTimesVector(v, gv, n, x, g, function, state,
      xscale, fscale, accuracy, xnorm, low, up);
    ++(*nfeval);
    if (frc) goto cleanup;
    project(n, gv, pivot);

    vgv = ddot1(n, v, gv);
    if (vgv / rhsnrm < tol)
    {
      /* Truncate algorithm in case of an emergency */
      if (k == 0)
      {
        frc = msolve(g, zsol, n, sk, yk, diagb, sr, yr, upd1, yksk, yrsr,
          lreset);
        if (frc) goto cleanup;
        dneg1(n, zsol);
        project(n, zsol, pivot);
      }
      break;
    }
    diagonalScaling(n, emat, v, gv, r);

    /* Compute linear step length */
    alpha = rz / vgv;

    /* Compute current solution and related vectors */
    daxpy1(n, alpha, v, zsol);
    daxpy1(n, -alpha, gv, r);

    /* Test for convergence */
    gtp = ddot1(n, zsol, g);
    pr = ddot1(n, r, zsol);
    qnew = (gtp + pr) * 0.5;
    qtest = (k + 1) * (1.0 - qold / qnew);
    if (qtest <= 0.5) break;

    /* Perform cautionary test */
    if (gtp > 0.0)
    {
      /* Truncate algorithm in case of an emergency */
      daxpy1(n, -alpha, v, zsol);
      break;
    }

    qold = qnew;
    rzold = rz;
  }

  /* Terminate algorithm */
  /* Store (or restore) diagonal preconditioning */
  dcopy1(n, emat, diagb);

cleanup:
  if (r) free(r);
  if (v) free(v);
  if (zk) free(zk);
  if (emat) free(emat);
  if (gv) free(gv);
  return frc;
}

/* 
 * Update the preconditioning matrix based on a diagonal version
 * of the bfgs quasi-newton update.
 */
static void diagonalScaling(int n, double e[], double v[], double gv[],
  double r[])
{
  int i;
  double vr, vgv;

  vr = 1.0/ddot1(n, v, r);
  vgv = 1.0/ddot1(n, v, gv);
  for (i = 0; i < n; i++)
  {
    e[i] += - r[i]*r[i]*vr + gv[i]*gv[i]*vgv;
    if (e[i] <= 1e-6) e[i] = 1.0;
  }
}

/*
 * Returns the length of the initial step to be taken along the
 * vector p in the next linear search.
 */
static double initialStep(double fnew, double fmin, double gtp, double smax)
{
  double d, alpha;

  d = fabs(fnew - fmin);
  alpha = 1.0;
  if (d * 2.0 <= -(gtp) && d >= mchpr1()) alpha = d * -2.0 / gtp;
  if (alpha >= smax) alpha = smax;

  return alpha;
}

/*
 * Hessian vector product through finite differences
 */
static int hessianTimesVector(double v[], double gv[], int n, 
  double x[], double g[], tnc_function *function, void *state, 
  double xscale[], double fscale,
  double accuracy, double xnorm, double low[], double up[])
{
  double dinv, f, delta, *xv;
  int i, frc;
  
  xv = malloc(sizeof(*xv)*n);
  if (xv == NULL) return -1;

  delta = accuracy * (xnorm + 1.0);
  for (i = 0; i < n; i++)
    xv[i] = x[i] + delta * v[i];

  unscalex(n, xv, xscale);
  coercex(n, xv, low, up);
  frc = function(xv, &f, gv, state);
  free(xv);
  if (frc) return 1;
  scaleg(n, gv, xscale, fscale);

  dinv = 1.0 / delta;
  for (i = 0; i < n; i++)
    gv[i] = (gv[i] - g[i]) * dinv;
    
  projectConstants(n, gv, xscale);

  return 0;
}

/*
 * This routine acts as a preconditioning step for the 
 * linear conjugate-gradient routine. It is also the 
 * method of computing the search direction from the 
 * gradient for the non-linear conjugate-gradient code. 
 * It represents a two-step self-scaled bfgs formula. 
 */
static int msolve(double g[], double y[], int n, 
  double sk[], double yk[], double diagb[], double sr[], 
  double yr[], logical upd1, double yksk, double yrsr, 
  logical lreset)
{
  double ghyk, ghyr, yksr, ykhyk, ykhyr, yrhyr, rdiagb, gsr, gsk;
  int i, frc;
  double *hg = NULL, *hyk = NULL, *hyr = NULL;

  if (upd1)
  {
    for (i = 0; i < n; i++) y[i] = g[i] / diagb[i];
    return 0;
  }

  frc = -1;
  gsk = ddot1(n, g, sk);
  hg = malloc(sizeof(*hg)*n);
  if (hg == NULL) goto cleanup;
  hyr = malloc(sizeof(*hyr)*n);
  if (hyr == NULL) goto cleanup;
  hyk = malloc(sizeof(*hyk)*n);
  if (hyk == NULL) goto cleanup;
  frc = 0;

  /* Compute gh and hy where h is the inverse of the diagonals */
  if (lreset)
  {
    for (i = 0; i < n; i++)
    {
      rdiagb = 1.0 / diagb[i];
      hg[i] = g[i] * rdiagb;
      hyk[i] = yk[i] * rdiagb;
    }
    ykhyk = ddot1(n, yk, hyk);
    ghyk = ddot1(n, g, hyk);
    ssbfgs(n, 1.0, sk, hg, hyk, yksk, ykhyk, gsk, ghyk, y);
  }
  else
  {
    for (i = 0; i < n; i++)
    {
      rdiagb = 1.0 / diagb[i];
      hg[i] = g[i] * rdiagb;
      hyk[i] = yk[i] * rdiagb;
      hyr[i] = yr[i] * rdiagb;
    }
    gsr = ddot1(n, g, sr);
    ghyr = ddot1(n, g, hyr);
    yrhyr = ddot1(n, yr, hyr);
    ssbfgs(n, 1.0, sr, hg, hyr, yrsr, yrhyr, gsr, ghyr, hg);
    yksr = ddot1(n, yk, sr);
    ykhyr = ddot1(n, yk, hyr);
    ssbfgs(n, 1.0, sr, hyk, hyr, yrsr, yrhyr, yksr, ykhyr, hyk);
    ykhyk = ddot1(n, hyk, yk);
    ghyk = ddot1(n, hyk, g);
    ssbfgs(n, 1.0, sk, hg, hyk, yksk, ykhyk, gsk, ghyk, y);
  }
  
cleanup:
  if (hg) free(hg);
  if (hyk) free(hyk);
  if (hyr) free(hyr);
  
  return frc;
}

/*
 * Self-scaled BFGS
 */
static void ssbfgs(int n, double gamma, double sj[], double hjv[],
  double hjyj[], double yjsj, 
  double yjhyj, double vsj, double vhyj, double hjp1v[])
{
  double beta, delta;
  int i;

  if (yjsj == 0.0)
  {
    delta = 0.0;
    beta = 0.0;
  }
  else
  {
    delta = (gamma * yjhyj / yjsj + 1.0) * vsj / yjsj - gamma * vhyj / yjsj;
    beta = -gamma * vsj / yjsj;
  }

  for (i = 0; i < n; i++)
    hjp1v[i] = gamma * hjv[i] + delta * sj[i] + beta * hjyj[i];
}

/*
 * Initialize the preconditioner
 */
static int initPreconditioner(double diagb[], double emat[], int n, 
  logical lreset, double yksk, double yrsr,
  double sk[], double yk[], double sr[], double yr[], 
  logical upd1)
{
  double srds, yrsk, td, sds;
  int i;
  double *bsk;

  if (upd1)
  {
    dcopy1(n, diagb, emat);
    return 0;
  }
  
  bsk = malloc(sizeof(*bsk)*n);
  if (bsk == NULL) return -1;
  
  if (lreset) 
  {
    for (i = 0; i < n; i++) bsk[i] = diagb[i] * sk[i];
    sds = ddot1(n, sk, bsk);
    if (yksk == 0.0) yksk = 1.0;
    if (sds == 0.0) sds = 1.0;
    for (i = 0; i < n; i++)
    {
      td = diagb[i];
      emat[i] = td - td * td * sk[i] * sk[i] / sds + yk[i] * yk[i] / yksk;
    }
  }
  else
  {
    for (i = 0; i < n; i++) bsk[i] = diagb[i] * sr[i];
    sds = ddot1(n, sr, bsk);
    srds = ddot1(n, sk, bsk);
    yrsk = ddot1(n, yr, sk);
    if (yrsr == 0.0) yrsr = 1.0;
    if (sds == 0.0) sds = 1.0;
    for (i = 0; i < n; i++)
    {
      td = diagb[i];
      bsk[i] = td * sk[i] - bsk[i] * srds / sds + yr[i] * yrsk / yrsr;
      emat[i] = td - td * td * sr[i] * sr[i] / sds + yr[i] * yr[i] / yrsr;
    }
    sds = ddot1(n, sk, bsk);
    if (yksk == 0.0) yksk = 1.0;
    if (sds == 0.0) sds = 1.0;
    for (i = 0; i < n; i++)
      emat[i] = emat[i] - bsk[i] * bsk[i] / sds + yk[i] * yk[i] / yksk;
  }
  
  free(bsk);
  return 0;
}


/*
 * Line search algorithm of gill and murray
 */
static ls_rc linearSearch(int n, tnc_function *function, void *state,
  double low[], double up[],
  double xscale[], double fscale, int pivot[],
  double eta, double ftol, double xbnd,
  double p[], double x[], double *f,
  double *alpha, double gfull[], int maxnfeval, int *nfeval)
{
  double b1, big, tol, rmu, fpresn, fu, gu, fw, gw, gtest1, gtest2,
    oldf, fmin, gmin, rtsmll, step, a, b, e, u, ualpha, factor, scxbnd, xw,
    epsmch, reltol, abstol, tnytol, pe, xnorm, rteps;
  double *temp = NULL, *tempgfull = NULL, *newgfull = NULL;
  int maxlsit = 64, i, itcnt, frc;
  ls_rc rc;
  getptc_rc itest;
  logical braktd;
  
  rc = LS_ENOMEM;
  temp = malloc(sizeof(*temp)*n);
  if (temp == NULL) goto cleanup;
  tempgfull = malloc(sizeof(*tempgfull)*n);
  if (tempgfull == NULL) goto cleanup;
  newgfull = malloc(sizeof(*newgfull)*n);
  if (newgfull == NULL) goto cleanup;

  dcopy1(n, gfull, temp);
  scaleg(n, temp, xscale, fscale);
  gu = ddot1(n, temp, p);

  dcopy1(n, x, temp);
  project(n, temp, pivot);
  xnorm = dnrm21(n, temp);

  /* Compute the absolute and relative tolerances for the linear search */
  epsmch = mchpr1();
  rteps = sqrt(epsmch);
  pe = dnrm21(n, p) + epsmch;
  reltol = rteps * (xnorm + 1.0) / pe;
  abstol = -epsmch * (1.0 + fabs(*f)) / (gu - epsmch);

  /* Compute the smallest allowable spacing between points in the linear
    search */
  tnytol = epsmch * (xnorm + 1.0) / pe;

  rtsmll = epsmch;
  big = 1.0 / (epsmch * epsmch);
  itcnt = 0;

  /* Set the estimated relative precision in f(x). */
  fpresn = epsmch * ftol;

  u = *alpha;
  fu = *f;
  fmin = *f;
  rmu = 1e-4;

  /* Setup */
  itest = getptcInit(&reltol, &abstol, tnytol, eta, rmu, 
    xbnd, &u, &fu, &gu, alpha, &fmin, &gmin, &xw, &fw, &gw, &a, &b,
    &oldf, &b1, &scxbnd, &e, &step, &factor, &braktd, &gtest1, &gtest2, &tol);

  /* If itest == GETPTC_EVAL, the algorithm requires the function value to be 
    calculated */
  while(itest == GETPTC_EVAL)
  {
    /* Test for too many iterations or too many function evals */
    if ((++itcnt > maxlsit) || ((*nfeval) >= maxnfeval)) break;

    ualpha = *alpha + u;
    for (i = 0; i < n; i++)
      temp[i] = x[i] + ualpha * p[i];

    /* Function evaluation */
    unscalex(n, temp, xscale);
    coercex(n, temp, low, up);

    frc = function(temp, &fu, tempgfull, state);
    ++(*nfeval);
    if (frc)
    {
      rc = LS_USERABORT;
      goto cleanup;
    }

    fu *= fscale;

    dcopy1(n, tempgfull, temp);
    scaleg(n, temp, xscale, fscale);
    gu = ddot1(n, temp, p);

    itest = getptcIter(big, rtsmll, &reltol, &abstol, tnytol, fpresn,
      xbnd, &u, &fu, &gu, alpha, &fmin, &gmin, &xw, &fw, &gw, &a, &b,
      &oldf, &b1, &scxbnd, &e, &step, &factor, &braktd, &gtest1, &gtest2, &tol);
    
    /* New best point ? */
    if (*alpha == ualpha)
      dcopy1(n, tempgfull, newgfull);
  }

  if (itest == GETPTC_OK)
  {
    /* A successful search has been made */
    *f = fmin;
    daxpy1(n, *alpha, p, x);
    dcopy1(n, newgfull, gfull);
    rc = LS_OK;
  }
  /* Too many iterations ? */
  else if (itcnt > maxlsit) rc = LS_FAIL;
  /* If itest=GETPTC_FAIL or GETPTC_EINVAL a lower point could not be found */
  else if (itest != GETPTC_EVAL) rc = LS_FAIL;
  /* Too many function evaluations */
  else rc = LS_MAXFUN;

cleanup:
  if (temp) free(temp);
  if (tempgfull) free(tempgfull);
  if (newgfull) free(newgfull);

  return rc;
}

/*
 * getptc, an algorithm for finding a steplength, called repeatedly by
 * routines which require a step length to be computed using cubic
 * interpolation. The parameters contain information about the interval
 * in which a lower point is to be found and from this getptc computes a
 * point at which the function can be evaluated by the calling program.
 */
static getptc_rc getptcInit(double *reltol, double *abstol, double tnytol, 
  double eta, double rmu, double xbnd, 
  double *u, double *fu, double *gu, double *xmin,
  double *fmin, double *gmin, double *xw, double *fw, 
  double *gw, double *a, double *b, double *oldf, 
  double *b1, double *scxbnd, double *e, double *step, 
  double *factor, logical *braktd, double *gtest1, 
  double *gtest2, double *tol)
{
  /* Check input parameters */
  if (*u <= 0.0 || xbnd <= tnytol || *gu > 0.0) return GETPTC_EINVAL;
  if (xbnd < *abstol) *abstol = xbnd;
  *tol = *abstol;

  /* a and b define the interval of uncertainty, x and xw are points */
  /* with lowest and second lowest function values so far obtained. */
  /* initialize a,smin,xw at origin and corresponding values of */
  /* function and projection of the gradient along direction of search */
  /* at values for latest estimate at minimum. */

  *a = 0.0;
  *xw = 0.0;
  *xmin = 0.0;
  *oldf = *fu;
  *fmin = *fu;
  *fw = *fu;
  *gw = *gu;
  *gmin = *gu;
  *step = *u;
  *factor = 5.0;

  /* The minimum has not yet been bracketed. */
  *braktd = TNC_FALSE;

  /* Set up xbnd as a bound on the step to be taken. (xbnd is not computed */
  /* explicitly but scxbnd is its scaled value.) Set the upper bound */
  /* on the interval of uncertainty initially to xbnd + tol(xbnd). */
  *scxbnd = xbnd;
  *b = *scxbnd + *reltol * fabs(*scxbnd) + *abstol;
  *e = *b + *b;
  *b1 = *b;

  /* Compute the constants required for the two convergence criteria. */
  *gtest1 = -rmu * *gu;
  *gtest2 = -eta * *gu;

  /* If the step is too large, replace by the scaled bound (so as to */
  /* compute the new point on the boundary). */
  if (*step >= *scxbnd)
  {  
    *step = *scxbnd;
    /* Move sxbd to the left so that sbnd + tol(xbnd) = xbnd. */
    *scxbnd -= (*reltol * fabs(xbnd) + *abstol) / (1.0 + *reltol);
  }
  *u = *step;
  if (fabs(*step) < *tol && *step < 0.0) *u = -(*tol);
  if (fabs(*step) < *tol && *step >= 0.0) *u = *tol;
  return GETPTC_EVAL;
}

static getptc_rc getptcIter(double big, double 
  rtsmll, double *reltol, double *abstol, double tnytol, 
  double fpresn, double xbnd,
  double *u, double *fu, double *gu, double *xmin,
  double *fmin, double *gmin, double *xw, double *fw, 
  double *gw, double *a, double *b, double *oldf, 
  double *b1, double *scxbnd, double *e, double *step, 
  double *factor, logical *braktd, double *gtest1, 
  double *gtest2, double *tol)
{
  double abgw, absr, p, q, r, s, scale, denom, 
    a1, d1, d2, sumsq, abgmin, chordm, chordu,
    xmidpt, twotol;
  logical convrg;

  /* Update a,b,xw, and xmin */
  if (*fu <= *fmin)
  {
    /* If function value not increased, new point becomes next */
    /* origin and other points are scaled accordingly. */
    chordu = *oldf - (*xmin + *u) * *gtest1;
    if (*fu > chordu)
    {
      /* The new function value does not satisfy the sufficient decrease */
      /* criterion. prepare to move the upper bound to this point and */
      /* force the interpolation scheme to either bisect the interval of */
      /* uncertainty or take the linear interpolation step which estimates */
      /* the root of f(alpha)=chord(alpha). */

      chordm = *oldf - *xmin * *gtest1;
      *gu = -(*gmin);
      denom = chordm - *fmin;
      if (fabs(denom) < 1e-15)
      {
        denom = 1e-15;
        if (chordm - *fmin < 0.0) denom = -denom;
      }
      if (*xmin != 0.0) *gu = *gmin * (chordu - *fu) / denom;
      *fu = 0.5 * *u * (*gmin + *gu) + *fmin;
      if (*fu < *fmin) *fu = *fmin;
    }
    else
    {
      *fw = *fmin;
      *fmin = *fu;
      *gw = *gmin;
      *gmin = *gu;
      *xmin += *u;
      *a -= *u;
      *b -= *u;
      *xw = -(*u);
      *scxbnd -= *u;
      if (*gu <= 0.0)
      {
        *a = 0.0;
      }
      else
      {
        *b = 0.0;
        *braktd = TNC_TRUE;
      }
      *tol = fabs(*xmin) * *reltol + *abstol;
      goto ConvergenceCheck;
    }
  }

  /* If function value increased, origin remains unchanged */
  /* but new point may now qualify as w. */
  if (*u < 0.0)
    *a = *u;
  else
  {
    *b = *u;
    *braktd = TNC_TRUE;
  }
  *xw = *u;
  *fw = *fu;
  *gw = *gu;

ConvergenceCheck:
  twotol = *tol + *tol;
  xmidpt = 0.5 * (*a + *b);

  /* Check termination criteria */
  convrg = (fabs(xmidpt) <= twotol - 0.5 * (*b - *a)) || 
    (fabs(*gmin) <= *gtest2 && *fmin < *oldf && ((fabs(*xmin - xbnd) > *tol) ||
    (! (*braktd))));
  if (convrg)
  {
    if (*xmin != 0.0) return GETPTC_OK;

    /*
     * If the function has not been reduced, check to see that the relative
     * change in f(x) is consistent with the estimate of the delta-
     * unimodality constant, tol. If the change in f(x) is larger than
     * expected, reduce the value of tol.
     */
    if (fabs(*oldf - *fw) <= fpresn * 0.5 * (fabs(*fw) + fabs(*oldf)))
      return GETPTC_FAIL;
    *tol = 0.1 * *tol;
    if (*tol < tnytol) return GETPTC_FAIL;
    *reltol = 0.1 * *reltol;
    *abstol = 0.1 * *abstol;
    twotol = 0.1 * twotol;
  }

  /* Continue with the computation of a trial step length */
  r = 0.0;
  q = 0.0;
  s = 0.0;
  if (fabs(*e) > *tol)
  {
    /* Fit cubic through xmin and xw */
    r = 3.0 * (*fmin - *fw) / *xw + *gmin + *gw;
    absr = fabs(r);
    q = absr;
    if (*gw != 0.0 && *gmin != 0.0)
    {
      /* Compute the square root of (r*r - gmin*gw) in a way
         which avoids underflow and overflow. */
      abgw = fabs(*gw);
      abgmin = fabs(*gmin);
      s = sqrt(abgmin) * sqrt(abgw);
      if (*gw / abgw * *gmin > 0.0) 
      {
        if (r >= s || r <= -s)
        {
          /* Compute the square root of r*r - s*s */
          q = sqrt(fabs(r + s)) * sqrt(fabs(r - s));
        }
        else
        {
          r = 0.0;
          q = 0.0;
          goto MinimumFound;
        }
      }
      else
      {
        /* Compute the square root of r*r + s*s. */
        sumsq = 1.0;
        p = 0.0;
        if (absr >= s)
        {
          /* There is a possibility of underflow. */
          if (absr > rtsmll) p = absr * rtsmll;
          if (s >= p)
          {
            double value = s / absr;
            sumsq = 1.0 + value * value;
          }
          scale = absr;
        }
        else
        {
          /* There is a possibility of overflow. */
          if (s > rtsmll) p = s * rtsmll;
          if (absr >= p)
          {
            double value = absr / s;
            sumsq = 1.0 + value * value;
          }
          scale = s;
        }
        sumsq = sqrt(sumsq);
        q = big;
        if (scale < big / sumsq) q = scale * sumsq;
      }
    }

    /* Compute the minimum of fitted cubic */
    if (*xw < 0.0) q = -q;
    s = *xw * (*gmin - r - q);
    q = *gw - *gmin + q + q;
    if (q > 0.0) s = -s;
    if (q <= 0.0) q = -q;
    r = *e;
    if (*b1 != *step || *braktd) *e = *step;
  }

MinimumFound:
  /* Construct an artificial bound on the estimated steplength */
  a1 = *a;
  *b1 = *b;
  *step = xmidpt;
  if ( (! *braktd) || ((*a == 0.0 && *xw < 0.0) || (*b == 0.0 && *xw > 0.0)) )
  {
    if (*braktd)
    {
      /* If the minimum is not bracketed by 0 and xw the step must lie
         within (a1,b1). */
      d1 = *xw;
      d2 = *a;
      if (*a == 0.0) d2 = *b;
      /* This line might be : */
      /* if (*a == 0.0) d2 = *e */
      *u = -d1 / d2;
      *step = 5.0 * d2 * (0.1 + 1.0 / *u) / 11.0;
      if (*u < 1.0) *step = 0.5 * d2 * sqrt(*u);
    }
    else
    {
      *step = -(*factor) * *xw;
      if (*step > *scxbnd) *step = *scxbnd;
      if (*step != *scxbnd) *factor = 5.0 * *factor;
    }
    /* If the minimum is bracketed by 0 and xw the step must lie within (a,b) */
    if (*step <= 0.0) a1 = *step;
    if (*step > 0.0) *b1 = *step;
  }

/*
 *   Reject the step obtained by interpolation if it lies outside the
 *   required interval or it is greater than half the step obtained
 *   during the last-but-one iteration.
 */
  if (fabs(s) <= fabs(0.5 * q * r) || s <= q * a1 || s >= q * *b1)
    *e = *b - *a;
  else
  {
    /* A cubic interpolation step */
    *step = s / q;

    /* The function must not be evaluated too close to a or b. */
    if (*step - *a < twotol || *b - *step < twotol)
    {
      if (xmidpt <= 0.0)
        *step = -(*tol);
      else
        *step = *tol;
    }
  }

  /* If the step is too large, replace by the scaled bound (so as to */
  /* compute the new point on the boundary). */
  if (*step >= *scxbnd)
  {  
    *step = *scxbnd;
    /* Move sxbd to the left so that sbnd + tol(xbnd) = xbnd. */
    *scxbnd -= (*reltol * fabs(xbnd) + *abstol) / (1.0 + *reltol);
  }
  *u = *step;
  if (fabs(*step) < *tol && *step < 0.0) *u = -(*tol);
  if (fabs(*step) < *tol && *step >= 0.0) *u = *tol;
  return GETPTC_EVAL;
}

/*
 * Return epsmch, where epsmch is the smallest possible
 * power of 2 such that 1.0 + epsmch > 1.0
 */
static double mchpr1(void)
{
  static double epsmch = 0.0;
  
  if (epsmch == 0.0)
  {
    double eps = 1.0;
    while((1.0 + (eps*0.5)) > 1.0)
      eps *= 0.5;
    epsmch = eps;
  }

  return epsmch;
}

/* Blas like routines */

/* dy+=dx */
static void dxpy1(int n, double dx[], double dy[])
{
  int i;
  for (i = 0; i < n; i++)
    dy[i] += dx[i];
}

/* dy+=da*dx */
static void daxpy1(int n, double da, double dx[], double dy[])
{
  int i;
  for (i = 0; i < n; i++)
    dy[i] += da*dx[i];
}

/* Copy dx -> dy */
/* Could use memcpy */
static void dcopy1(int n, double dx[], double dy[])
{
  int i;
  for (i = 0; i < n; i++)
    dy[i] = dx[i];
}

/* Negate */
static void dneg1(int n, double v[])
{
  int i;
  for (i = 0; i < n; i++)
    v[i] = -v[i];
}

/* Dot product */
static double ddot1(int n, double dx[], double dy[])
{
  int i;
  double dtemp = 0.0;
  for (i = 0; i < n; i++)
    dtemp += dy[i]*dx[i];
  return dtemp;
}

/* Infinity norm */
static double dnrmi1(int n, double v[])
{
  int i;
  double dtemp, dmax;
  for (dmax = fabs(v[0]), i = 1; i < n; i++)
    if ((dtemp = fabs(v[i])) > dmax) dmax = dtemp;
  return dmax;
}

/* Euclidian norm */
static double dnrm21(int n, double dx[])
{
  int i;
  double dssq = 1.0, dscale = 0.0;

  for (i = 0; i < n; i++)
  {
    if (dx[i] != 0.0)
    {
      double dabsxi = fabs(dx[i]);
      if (dscale<dabsxi)
      {
        /* Normalization to prevent overflow */
        double ratio = dscale/dabsxi;
        dssq = 1.0 + dssq*ratio*ratio;
        dscale = dabsxi;
      }
      else
      {
        double ratio = dabsxi/dscale;
        dssq += ratio*ratio;
      }
    }
  }

  return dscale*sqrt(dssq);
}
                                                                                                                                                                                                                          r8s/tnc.h                                                                                           0000644 0000766 0000120 00000014435 10357567246 012226  0                                                                                                    ustar   sandermj                        admin                                                                                                                                                                                                                  /* tnc : truncated newton bound contrained minimization
         using gradient information, in C */

/*
 * Copyright (c) 2002-2004, Jean-Sebastien Roy (js@jeannot.org)
 * 
 * Permission is hereby granted, free of charge, to any person obtaining a
 * copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/*
 * This software is a C implementation of TNBC, a truncated newton minimization
 * package originally developed by Stephen G. Nash in Fortran.
 * 
 * The original source code can be found at :
 * http://iris.gmu.edu/~snash/nash/software/software.html
 * 
 * Copyright for the original TNBC fortran routines:
 * 
 *   TRUNCATED-NEWTON METHOD:  SUBROUTINES
 *     WRITTEN BY:  STEPHEN G. NASH
 *           SCHOOL OF INFORMATION TECHNOLOGY & ENGINEERING
 *           GEORGE MASON UNIVERSITY
 *           FAIRFAX, VA 22030
 */

/* $Jeannot: tnc.h,v 1.53 2004/04/18 10:32:30 js Exp $ */

#ifndef _TNC_
#define _TNC_

#define TNC_VERSION "1.2.5"

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Verbosity level
 */
typedef enum {
  TNC_MSG_NONE = 0, /* No messages */
  TNC_MSG_ITER = 1, /* One line per iteration */
  TNC_MSG_INFO = 2, /* Informational messages */
  TNC_MSG_VERS = 4, /* Version info */
  TNC_MSG_EXIT = 8, /* Exit reasons */

  TNC_MSG_ALL = TNC_MSG_ITER | TNC_MSG_INFO
    | TNC_MSG_VERS | TNC_MSG_EXIT /* All messages */
} tnc_message;

/*
 * Possible return values for tnc
 */
typedef enum
{
  TNC_MINRC        = -3, /* Constant to add to get the rc_string */
  TNC_ENOMEM       = -3, /* Memory allocation failed */
  TNC_EINVAL       = -2, /* Invalid parameters (n<0) */
  TNC_INFEASIBLE   = -1, /* Infeasible (low bound > up bound) */
  TNC_LOCALMINIMUM =  0, /* Local minima reach (|pg| ~= 0) */
  TNC_CONVERGED    =  1, /* Converged (|f_n-f_(n-1)| ~= 0) */
  TNC_MAXFUN       =  2, /* Max. number of function evaluations reach */
  TNC_LSFAIL       =  3, /* Linear search failed */
  TNC_CONSTANT     =  4, /* All lower bounds are equal to the upper bounds */
  TNC_NOPROGRESS   =  5, /* Unable to progress */
  TNC_USERABORT    =  6  /* User requested end of minization */
} tnc_rc;

/*
 * Return code strings
 * use tnc_rc_string[rc - TNC_MINRC] to get the message associated with
 * return code rc.
 */

extern char *tnc_rc_string[10];

/*
 * A function as required by tnc
 * state is a void pointer provided to the function at each call
 *
 * x     : on input, then vector of variables (should not be modified)
 * f     : on output, the value of the function
 * g     : on output, the value of the gradient
 * state : on input, the value of the state variable as provided to tnc
 *
 * must returns 0 if no error occurs or 1 to immediately end the minimization.
 *
 */
typedef int tnc_function(double x[], double *f, double g[], void *state);

/*
 * tnc : minimize a function with variables subject to bounds, using
 *       gradient information.
 *
 * n         : number of variables (must be >= 0)
 * x         : on input, initial estimate ; on output, the solution
 * f         : on output, the function value at the solution
 * g         : on output, the gradient value at the solution
 *             g should be an allocated vector of size n or NULL,
 *             in which case the gradient value is not returned.
 * function  : the function to minimize (see tnc_function)
 * state     : used by function (see tnc_function)
 * low, up   : the bounds
 *             set low[i] to -HUGE_VAL to remove the lower bound
 *             set up[i] to HUGE_VAL to remove the upper bound
 *             if low == NULL, the lower bounds are removed.
 *             if up == NULL, the upper bounds are removed.
 * scale     : scaling factors to apply to each variable
 *             if NULL, the factors are up-low for interval bounded variables
 *             and 1+|x] fo the others.
 * messages  : see the tnc_message enum
 * maxCGit   : max. number of hessian*vector evaluation per main iteration
 *             if maxCGit == 0, the direction chosen is -gradient
 *             if maxCGit < 0, maxCGit is set to max(1,min(50,n/2))
 * maxnfeval : max. number of function evaluation
 * eta       : severity of the line search. if < 0 or > 1, set to 0.25
 * stepmx    : maximum step for the line search. may be increased during call
 *             if too small, will be set to 10.0
 * accuracy  : relative precision for finite difference calculations
 *             if <= machine_precision, set to sqrt(machine_precision)
 * fmin      : minimum function value estimate
 * ftol      : precision goal for the value of f in the stoping criterion
 *             relative to the machine precision and the value of f.
 *             if ftol < 0.0, ftol is set to 0.0
 * rescale   : Scaling factor (in log10) used to trigger rescaling
 *             if 0, rescale at each iteration
 *             if a big value, never rescale
 *             if < 0, rescale is set to 1.3
 * nfeval    : on output, the number of function evaluations.
 *             ignored if nfeval==NULL.
 *
 * The tnc function returns a code defined in the tnc_rc enum.
 * On output, x, f and g may be very slightly out of sync because of scaling.
 *
 */
extern int tnc(int n, double x[], double *f, double g[],
  tnc_function *function, void *state,
  double low[], double up[], double scale[],
  int messages, int maxCGit, int maxnfeval, double eta, double stepmx,
  double accuracy, double fmin, double ftol, double rescale, int *nfeval);

#ifdef __cplusplus
}
#endif

#endif /* _TNC_ */
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   