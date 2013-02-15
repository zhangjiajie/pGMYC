void simulateBinaryChar(TREE t, double q01, double q10 );
void findMarginals(NODETYPE *root);
void printChanges(NODETYPE *node);
void printCovarion(NODETYPE *node, int doMarginals);
void covarionOptimize(TREE t,int *numIter, double ftol,double linMinDelta,int *success, int nstates, int doMarginals, int estimate, int doRecon );
double objCovarion(double p[]);
double objBinaryTrait(double p[]);
double objBinaryTraitSymmetric(double p[]);
double objCovarionFixed(double p[]);
double objCovarion4(double p[]);

