#include "distance.h"
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
