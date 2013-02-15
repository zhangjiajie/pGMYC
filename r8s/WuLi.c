#include <math.h>
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
		printf("%+6.4f",z);
		for (i=0;i<3;i++) if (fabs(z) < zscore[i]) break; 
		for (j=0;j<i;j++) printf("*");
		for (j=i;j<3;j++) printf(" ");
		}
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
	long	seqLength,nP[4][4],nQ[4][4],mA,mB;
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
				seqLength=PQCalc1(pi,pj,pRow1,&P[i][j],&Q[i][j],&nP[i][j],&nQ[i][j]);
			

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

// ** new stuff 

// ratio = ()/();

//


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
