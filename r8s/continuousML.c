/*

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

