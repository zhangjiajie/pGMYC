#define MAX_TAXON_ARRAY 500

int randTrinary (double p0,double p1); // Return a 0,1, or 2 with probability of p0,p1, or 1-p0-p1.
int randBinary (double p); // Return a 1 with probability p; a 0 with prob 1-p.
double RY_1997_Dist(double speciation, double extinction, double sampling) ;
double birthDist(double lambda);
float normdev(void);
float factrl(int n);
double gammln(double xx);
double poidev(double xm);
long rndLong(long maxLong);
double hgeom(double param);
double hexp(double param);
double myRand(void);
void getseed(void);
void bshuf3(float *weightArray, long numChars,int nreps,  char * buffer1, 
	char * buffer2);
void bshuf2(int *targetArray, long numChars);
void bshuf(int *targetArray, int *excludeArray, 
	long numChars, long includedChars);
void taxon_sample(int numtaxa, int numfixed, int numrandom, 
		int fixed[], int sample[],
		int nstart,int nstop,int nstart2,int stop2,int numrandom2);
int * taxon_sample_simple(int, int);
void remove_array_item(int array[], int length, int item_index);
