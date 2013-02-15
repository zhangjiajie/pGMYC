
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
