#include "DistrFuncs.h"
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

int randBinary (double p) // Return a 1 with probability p; a 0 with prob 1-p.
{
if (myRand() < p)
	return 1;
else
	return 0;

}
int randTrinary (double p0,double p1) // Return a 0,1, or 2 with probability of p0,p1, or 1-p0-p1.
{
double r;
r=myRand();
if (r < p0) return 0;
if (r < p0 + p1) return 1;
return 2;
}
