#include "TreeUtils.h"
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
