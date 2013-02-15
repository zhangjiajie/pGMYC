#include "structures.h"
void doDistance(StrListPtr aTaxaList);
long aaCalc1(char * pi, char * pj, char * pRow1, double *P,long *n);
long PQCalc(int i, int j, double *P, double *Q);
long PQCalc1(char * pi, char * pj, char * pRow1, double *P, double *Q, long *nP, long * nQ);
double distance(char * pi, char * pj, char * pRow1, int kind, long * nP, long * nQ);
