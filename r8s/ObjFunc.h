#ifndef _LANG_FITCH
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
int checkGradientSimple(double solution[],double gradient[],double Obj,double ftol,int nParm);
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
