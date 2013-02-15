#include "NRCvectorUtils.h"
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
