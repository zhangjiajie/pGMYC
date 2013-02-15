#define WULI	0
#define STEEL	1
#define TAJIMA	2
#define MIKE	3
int WuLiTest(int inGroup1,int inGroup2, int outGroup,
		 double *dif, double *stddev,double *Poa, double *Pob, double *doa, double *dob);
long tripletSites(int i, int j, int k, double *P, long *MA, long *MB);
void WuLiStub(int inGroup1,int inGroup2, int outGroup);
