int warnEstRoot(NODETYPE * root);
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
