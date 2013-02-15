#ifndef _NEXUS
#define _NEXUS

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdlib.h>
#include "TreeUtils.h"
#include "structures.h"
#include "MyUtilities.h"

#define MAX_TOKEN_SIZE 10000		/* we've got room */
#define MAX_LOCAL_TOKEN_SIZE 128
#define MAXTREES	100		/* maximum number of tree descriptions stored from file */

#if 0
#define isEqual(a,b)		(!strcmp((a),(b)))
#define isEqualUL(a,b)		(!strcmp((a),(b)))
#endif

/* (!strcmp((strtoupper(a), a),(strtoupper(b), b)))*/

#define isNEXUSpunct(c) 		( strchr(punct,(c)) )
#define isNEXUSwhiteSpace(c)	( isspace((c)) || (((c) <= 6) && ((c) >= 0)))
	/* current NEXUS format also excludes ASCII 0-6 */



struct	RBP 	{
				int	clampRoot; /* 0 = separate subtree
						optimizations */
				int	isBS;	/* toggle bootstrap */
				long 	NReps;	/* bootstrap replicates */
				long 	seed;	/* random number seed */
				int	RRtype; /* 0=WuLi; 1=Steel et al.*/
				double	npexp;	/* exponent in the NP optimization */
				int	verbose;/* verbosity for rate block */
				int	num_restarts;
				int	num_rate_guesses;
				int	num_time_guesses;
				double	local_factor; /* fractional tolerance for
					two points to be considered the same */
				double perturb_factor; /*fractional displacement to look
					for another optimum */
				double	smoothing;	/* smoothing factor in penalized like*/
				double  ftol;		/* fraction func tolerance */
				double	barrierTol;
				int	maxIter;
				int	maxBarrierIter;
				double	initBarrierFactor;
				double  barrierMultiplier;
				double	linminOffset;
				double	contractFactor;
				int	maxContractIter;
				int	showConvergence;
				int	checkGradient;
				int	showGradient;
				int	RatesAreGamma;	/* across sites */
				double  alpha;		/* shape param */
				double activeEpsilon;	/* fractional distance from a constraint; if closer than this distance, a solution is said to be "active" */
				long	numSites;
				int	clockFmt;	/* 1 = trees assumed to be ultrametric on input */
				int	lengthFmt;	/* 0 = branch lengths are in numbers of subst.; 1 = subst/site */
				int	roundFlag;	/* 0 = branch lengths are not rounded on input; 1 = rounded */
				int	PenaltyType;
				int	NeighborPenalty; // 1=penalize with neighbor variance; 0=old style ancestor/desc squared
				float minRateFactor; // a fraction of the average rate to impose a min on all rates under PL
				float minDurFactor; // a fraction of the root's age to impose a min duration for 0-length 
							// terminal branches 
				int estCov;	// should we try to estimate the covarion matrix rate parameters?
				double s_rate;
				double r_rate; // s and r rates in covarion matrix
				int	cov_brlens; // set to 1 if we will set all branch lengths to 1; otherwise use supplied values
				};	
				

/* This is the data structure containing all the information for a NEXUS file */

struct 	NexDataType {
			int		isChars;
			int		isTrees;
			int		isTaxa;
			int		isTranslate;	/*...flags for when these elements are read */
			int 		NTaxa;			/* number of taxa */
			int 		NChars;			/* number of characters */
			int 		Intlvflag;		/* flag is set if data matrix is interleaved */
			char		matchchar;
			char		gapchar;
			char		missingchar;
			int			NumTrees;		/* number of trees in data structure */
			StrListPtr 	TaxaList;		/* list of taxon names */
			StrListPtr 	TDList;			/* list of tree descriptions */
			StrListPtr 	TDLabelList;		/* list of tree description labels */
			StrListPtr 	DMList;			/* The data matrix as a list of row strings*/
			StrListPtr	TransList;		/* Translation table list */
			StrListPtr 	TaxSetNameList;		/* list of names of taxsets*/
			int*		excArray;		/* array of flags for excluding characters */
			struct RBP	RateBlockParms;	/* the rate block parameters */
			PtrList		inTrees;	/* list of trees */
			PtrList		TaxSetLists;	/* list of pointers to the taxsets, each of which
								is a StrList */
				};
				
/* All lists use my 'list' data structure */
void readNexusFile(char * buffer);

void TreeSummary(int whichTree);
char *nextToken(void);
void freeNexusStructure(struct NexDataType *nex);
struct NexDataType * initialize_nexus(void);
void doRateBlock(void);
void doCommandLineControl(char *buffer);
void doInteractive(void);
#endif
