/* REVISION HISTORY
 * 
 * 8.14.99.   Radical.  Changed everything to case insensitive by redefining isEqual macro
 *	    in structures.h, and changing the strtoupper() statements in nexttoken2.c
 *	    This is all more or less togled by the STU macro defn in structures.h
 * 9.3.99  Radical. Changed parse-assignment function to use a static
	    char buffer for LocalToken. New function is parse_assignment2

   9.13.99 Fixed 'compar()' function used in qsort routines TWICE in this program
    It was bogus and worked on SGI (coincidentally) but not on LINUX

   10.13.99 Changed 'BranchLike' in the LF algorithm so it converts real branch lengths
	to integer branh lengths before calculating likelihoods...formerly this was allowing
	real values for branch lengths...although this should never have come up
    1.16.00 Fixed nasty bug in setupFeasibleTimes. Had failed to initialize
	minTime=0.0 prior to descMinAge call.  Only a problem on some hardware.
    1.16.00 Fixed bug in ABCSuperTree routines that assigned a large double (1e100) as a
	branch length to the wtSet[] array, which is a Float, causing an exception on some
	hardware.  In the newNode function, I now use FLT_MAX instead of 1e100 (KLUDGE!) as
	the nulll branch length, and in TreeDrawUtils, which is the only place I need the
	stupid thing, I recognize FLT_MAX too.
    4.19.00 Changed the multiply_branchlength_by command so that final branch lengths are rounded
	to the nearest integer.  This allows me to run the collapse command afterwards, and take
	care of those nasty zero-length branches. (otherwise, I could collapse, but this would be prior
	to rounding down to the integer of zero, and then I was stuck with zero-length branches again
    4.21.00  Added a feature allowing 'setage taxon=root age=xx' to permit fixing the root node age of the tree.
	This corrects a bug that allowed the default root age to be greater than 1.0 whenever internal constraints
	were set with the constrain_time commands.  Now you have to explicitly set the root's age, or watch out!
    9.19.00  Changed command syntax for several commands so that options are remembered between runs (by
	changing to static values of parameters in functions. Incl.: divtime, describe.

   11.1.00 Series of significant modifications to code. Working on release 1.0.

    5.26.01 Changed PL and  NP routines to use the variance of rates across the root's children's branches
		as the quantity to be minimized. This may have slight effects on previous runs. Updated gradients
		for PL and objfuncs for PL and NP to do this.
    5.26.01 Added gamma distributed rates functionality for PL (w/POWELL only)--gradients not done yet...Note there
		is an issue in that BranchLike uses units of raw substitutions, whereas BranchLikeNegBinom uses subst/site.
    7.9.01 Added estimates of confidence intervals on a single node time using curvature of
	likelihood surface (see doDivTime:confidence)
	12.11.01 Fixed factorial calculations to allow any argument
	12.18.01 Fixed nasty bug in gradient calculations. Had not initialized g=0, which
			was a problem sometimes in estimating root node (long story)
	2.8.02 Fixed collision between CALIBRATE command and FIXAGE command. Used to permit fixage on an
		ultrametric tree, but calibrate command expects that ages were all scaled on 0 to 1 and not
		messed with after that initialization. Added warnings.
	2.16.02 Dumb bug: NPRS wanted to do a gradient check in MinND, but I don't have the gradient! Added test
		(should add gradient and make termination criteria a clean separate routine
	2.17.02 Added ROUND=YES|NO option to BRLEN command. This permits user to NOT round branch lengths
		to nearest integer. For input ultrametric trees, we'd really like to keep exact real lengths
		for use in the CALIBRATE command
	3.28.02 (Matt Lavin) Gradient check was too strict when using constrained nodes (should ignore non-zero
		gradient in the direction of the constraint if the constraint is active). Deciding that the 
		constraint is active relies on ACTIVE_EPS tolerance factor, which I've increased to 0.01
	4.15.02 (Ben Warren) Roundoff error bugs:
			1. In Langley-Fitch, use of very large ages (e.g., 500000) caused inconsistent results.
			The ZERO() macro in the BranchLike function was stupidly rounding arguments to 0.0 when
			they were merely small. I've just removed this macro.
			2. In NPRS, the same inputs were causing problems. Here the culprit was adding +1.0 to the
			objective function in ObjNP(). Originally I'd added this because Powell's termination criteria
			get fooled whenever there is truly a clock for NPRS (then the objective function = 0). Temp
			solution is to drop the addition of 1.0 but this will probably re-introduce the other bug. Need
			to improve Powell or use other optimization engines. 
	6.4.02 (Torsten Eriksson) Fixed very stupid memory leak in Powell. I had left some debugging code in there which
			allocated but did not free meory in linmin
		Also disabled the convergence diagnostics for Powell; probably slow down program a lot; can easily uncomment
		the changes in 'powell1'
	6.10.02 Several changes to BranchLike and LogFactLookup. Bug fixes concerning 0-length branches. See code for details 

Some insights on zero-length branches under PL: The likelihood component for a branch with k=0 substitutions is exp(-rT), 
which is a maximum along the directions of both r=0 and T=0. This could easily cause optimization methods to have problems,
although it seems that gradient methods perform worst. Once r =~ 0 or T =~ 0, the gradient will be 0 for the other variable
and the function will be maximized. For high smoothing values, the neighboring branches with k>0 substitutions keep these problems
in check, but as smoothing gets low, there comes a point at which r can be 0 without too  much of a penalty...then we have 
problems optimizing. The solution apparently is to set a small but finite minimum lower bound on rates.
	6.11.02. Rewrote and cleaned up tests in BranchLike. Better check to make sure this doesn't mess everything up downstream 
	6.22.02. Modified the check_feasible routine! Now a point is considered feasible even if a node time is EQUAL to a min or max
		 age constraint. Prior to this, I required strictly greater or less than. This is necessary because under TN algorithm
		 we frequently find points exactly satisfying constraint. Then when 'perturb_p' tries to work, it checks every node in 
		 tree for strict inequality, failed and bombed. [can also solve just by not restarting...]
	6.23.02. 'Bug' fixed in TNBC line search. For 0-length branches, the linesearch kept iterating. Now we force termination in the
			hopes that a restart perturbation will fix things up!
	6.24.02. Bug fixed in cross-validation. Zero-length branches were sometimes causing spurious calculations of chi-squared values.
		Needed to check for zero-length branches in cvSquareErrorBranch routine and prevent division by zero.																			
	7.25.02  See fix on 6.22.02. Didn't quite fix this right. Made sure all tests were for strict inequality
	7.26.02  Removed "include <malloc.h>" in ReadNex... and memory.c, the only two places it is invoked. Couldn't find this include
		 file under Mac OS X.
	7.26.02  Fixed bug in pruneTaxon command, which did not correctly update all pointers in tree structure when a taxon was deleted,
		 causing erroneous print outs under "describe plot=xxx_description". Other option values in describe worked OK.
	11.6.02  Eliminated all drand48 and srand48 references, and replaced with myRand() function which uses stdlib rand(). Also call srand() rather than srand48 everywhere. Hopefully this will make cross-platform development easier.
	11.6.02  (Torsten Eriksson) Eliminated, hopefully, a memory leak in TN (TNwrapper), where I failed to deallocate arrays.
	11.?.? [sometime while in Germany] Added a log transformation to the penalty function in NP and PL. Invoked by SET PENALTY=ADD|LOG;
		Note that this only works with POWELL-haven't calculated derivatives yet.
		Also began to implement a general 'neighbor variance' penalty, not done, see TimeAlgorithms...
	4.3.03  Upon migrating to MAC OS 10.2 found a bug in the COLLAPSE command which is only fixed by REMOVING the node 
		destructor (was screwing up the recursion). Now I've got a bunch of dangling nodes still! Figure out how to
		deallocate them (perhaps maintain a global junk list)
	4.31.03 Added feature. A likelihood ratio based relative rate test at a user supplied node in the tree: RRLIKE TAXON=XXX;
	5.13.03 Modified warnEstRoot so that it does not issue warnings itself but merely returns error code.
		Continued work on RRLIKE command. It now always takes time constraints or fixed ages into account if they are available for the clades
		in question. If no constraints or fixed ages are available, then it sets the root of the focal subtree to an arbitrary
		value.
	6.3.03 Fixed the 'birthDist' function to permit large values of lambda. Currently these overflowed. Can replace the
		density with a simpler form when lambda is large
	6.8.03 In TreeSim, I now require a new seed for each run; otherwise it uses a default seed and issues a warning. I didn't
		like having a seed held over from previous runs.
	7.16.03 Lots of stuff:
			working on a fossil cross validation scheme
			added log penalty scaling with gradient this time to penalized likelihood
			added neighbor variance penalty to PL in conjunction with log scaling
	9.19.03 Added the VCV function which calculates the variance covariance matrix of a tree based on the lengths
			of the branches subtending the MRCA of each pair of taxa.
		Syntax is VCV taxon = name; which works on the subtree descended from name
	3.3.04  Added another fossil cross validation scheme, fossilcrossvfixed, which uses only a set of fixed ages
	6.14.04 Made two submodels of BDback, depending on whether we normalize the root age to 1. 
			diversemodel=bdback command does not; 'bdbacknormal' does.
	8.26.04 Fixed bug in profile branch command which did not correct for trees in which branch was missing (A. Antonelli)
	8.26.04 Fixed bug in doObjFunc which mistakenly set algorithm to 'TN' in multiple time guesses by an expression 'if (algorithm=TN)' blunder (changed to '==' : also A. Antonelli).
	8.26.04 Important changes to PL routines: now estimating the initial starting guess on the rate by doing a LaF clock
		analysis first and taking that estimator as the guess. This required putting a wrapper around the doObjFunc
		routine. Other methods unaffected, still using crude guess.
	8.26.04 Begun some error reporting. Will now signal error if the basic command is wrong, but still overlooks merely wrong option setting
		syntax
	...?	Implemented CO command (continuous character rate estimation)
	12.4.04 Beginning work on checkGradient routine. Problems identified in correct setting for ACTIVE_EPS. Warning now
		issued if a min and max constraint might BOTH be treated as generating an active constraint (and thus the program
		arbitrarily picks one, leading to a false conclusion that the gradient's sign is wrong at the constraint)
	12.6.04 Problem noted in optimization when user specifies ROUND=NO. The calculation of the objective function rounds
		the character lengths on a branch always, but the calculation of the gradient (at least in LF), uses real 
		arithmetic. This causes routines to come to different conclusions depending on whether it uses gradient or 
		non gradient methods. Effect is slight on param estimates or obj func estimate, but in case examined when 
		checking the gradient using numerical approx it is quite important. Added warning message in the BLFORMAT 
		command, and now DIVTIME will bail with a warning if you try to use any method other than NPRS without
		rounding input (NPRS doesn't rely on calculation of a likelihood, so not an issue
	12.8.04 Fixed bug: ObjFunc.c was subtracting 1.0 from obj() under NPRS because at one point I had added +1.0 to it
		in TimeAlgorithms.c (to fix what I thought was a problem with clocklike data sets returning an objective
		function value of 0.0). However, I had stopped adding the 1.0 in the obj() in TimeAlgorithms.c.
		This was sometimes generating NPRS obj() values less than zero. Now the two modules consistently do not
		modify the objective.
	12.8.05 Changed 'execute' command so that we can read multiple blocks from different files, for example, to read
		a trees block from one file and a r8s block from another
	12.10.05 Implemented an analytical gradient for log penalty function under PL.
	1.8.05 Added charset ouput to MRP command.
	4.29.05. Added warning about only using POWELL with LFLOCAL in response to 'bug report'
	12.05. Added ancestral state reconstruction using squared change parsimony
	1.06.06. Removed all traces of 4PL and jk4PL code. Mix of C and C++ causing some problems on some compilers
	1.06.06. Using -pedantic gcc option to find misc bugs. Removed HUGE_VAL from two functions (this is a reserved word)
	7.25.06. False bug in user supplied ultrametric tree stuff...make sure user sets ROUND=NO and lengths=persite if latter
			is appropriate. Added routing rootToTips to print out those distances. Found slight roundoff error in PAUP
			tree descriptions in this respect.
Find Hilmar Lapp's email, where he found some bugs in interactive mode that need fixing....
	6/22/07. The rounding issue in branch lengths came up again. To reiterate from above, *all* branch lengths input to doObjFunc...must
			be integers. This is because the likelihood function converts them to integers but the gradient functions do not, leading
			to frequent catastrophic convergence failures...Major change to doObjFunc is to call 'traverseMultiplyLength', which is used
			to force rounding of all branch lengths before any optimization is done. Now there is no choice. Previously, I warned the user
			about the issue but didn't actually fix it.	
    9/2011.  Adding a new feature on a 3-state markov character model (covarion.c).
    		 Added a feature for options on rerooting, and fixed some rerooting bugs. Now we can
    		 reroot at a NODE, instead of just former behavior where we always had a binary root.
    		 
    12/6/2010		 NEW FEATURE SKIPS ANYTHING BEFORE '#NEXUS' in input file. WATCH OUT!
    11/14/2012	Well, after a bunch of additions for precursor code, I added a bug to the RemoveTaxonLvAnc function. It was ok for the new
		code, but it broke the cross validation code. A simple fix seems to work. [In the function, I had set the deleted subtree's		   ancestor to NULL (reasonable, eh?), but forgot that the cross val code, needed that anc to be save, so now I save it in 
		another variable. Thanks to Tomas Flouri.]

    		 */



/****  Module for Nexus File functions  *******/

#define EFRON1996	1  /* only needed for this flavor of bootstrapping */


/*******************************************************************************/

#include "continuousML.h"
#include "NRCvectorUtils.h"
#include "storeNexusFile.h"
#include <sys/types.h>
/* #include <malloc.h> */
#include "Maximize.h"
#include "WuLi.h"
#include "nexus.h"
#include "MyUtilities.h"
#include "memory.h"
#include "ObjFunc.h"
#include "TimeAlgorithms.h"
#include "myOutput.h"
#include "TreeSim.h"
#include "ObjFunc.h"
#include "DrawTree.h"
#include "moment.h"
#include "DistrFuncs.h"
#include "distance.h"
#include "structures.h"
#include "TreeUtils.h"
#include "ancestral.h"
#include "covarion.h"
#include "relativeRates.h"
#include <math.h>
#include <ctype.h>


/*****  private functions ******/

int parse_assignment2(char * target);
static void doVCVCommand(void);

static void doCovarionCommand(void);
static void doAncestralCommand(void);
static void doContOptCommand(void);
static void histoStat(long h[], long N, long nTaxa,long *count, double *mean, double *freq1class, long *maxS, double *dominance);
static void doRRLikeTestCommand(void);
static void doConfidence(TREE thisTree,char * nodeName,int method,int nRates,int algorithm,double cutoff,int JMAX);
static void doLocalModelCommand(void);
static void doBLFormatCommand(void);
static void doFossilCrossVfixed(TREE tree, int method,int nRates,int algorithm,double cvStart,double cvInc,double cvNum, int cvSample);
static void doFossilCrossV(TREE tree, int method,int nRates,int algorithm,double cvStart,double cvInc,double cvNum, int cvSample);
static float doCrossV(TREE tree, int method, int nRates,int algorithm, double c1, double c2, double c3, int);
static void doUnSetAgeCommand(void);
static void doShowAgeCommand(void);
static void doReRootCommand(void);
static void doPruneTaxonCommand(void);
static void doSetAgeCommand(void);
static void doClusterHistogramCommand(void);
static void doCollapseCommand(void);
static void doClearTreesCommand(void);
static void doExecuteCommand(void);
static void printHelp(void);
static void doDivTimeCommand(void);
static void doSimpleCladeCheckCommand(void);
static void doEFRON_Weights_Command(char *buffer);
static void doB_Weights_Command(char * buffer, char *buffer1, char* buffer2);
static void efron1996(int *weightArray,int nchars,int num_points,char *buffer, 
    long *index);
static void doCladeCheckCommand(void);
static void doBootCharCommand(char* buffer);
static void doBranchProfileCommand(void);
static void doClade_Set_Command(void);
static void printNexus(int ntaxa, int nchars, StrListPtr taxaList,  char **matrix);
	
static void doSuperCommand(void);  
static void doCalibrateCommand(void);
static void doPrintCommand(void);
static void doSimBlock(void);
static void doSimCommand(void);
static void doTaxaSetCommand(void);
static void doSetCommand(void);
static void doDistanceCommand(void);
static void doBSCommand(void);
static void doError(char* p[], int which);
static void doDataBlock(void);
static void doBootBlock(void);
static void doTranslateCommand(void);
static void doWuLiCommand(NODETYPE * root);
static void doCharsBlock(void);
static void doFormat(void);
static void doMatrix(void);
static void doIndel(void);
static void doTaxaBlock(void);
static void doCharDimensions(void);
static void doTaxDimensions(void);
static void doTaxLabels(void);
static void doTreeBlock(void);
static void doTreeCommand(void);
static void doExSets(void);
static void checkMatrix(void);
static void doUnrecognizedBlock(void);
static void doSitesCommand(int); /* exclude third positions */
static void doConstrain_TimeCommand();
static void doBootCommand(StrListPtr b, char* bu);
static char * doIncludeCommand(void);
static StrListPtr doFixedTaxaListCommand(void);
static void doMRCACommand(void);
static void doLengthMultiplyCommand();

static void doSaveTree(NODETYPE *root);
static int parse_assignment(char * target,char ** token);
static void doBD(void);
void doDimensions(void);
void doMatrixGeneral(void);

/******   globals   ********/

StackPtr gFStack,gPStack;
char LocalToken[MAX_LOCAL_TOKEN_SIZE];

int gEstRoot;
int gInteractive,gLabel;
int gSeedisSet=0;

StrListPtr	gTaxaList;
int	gNewLine;
int	gColumn;
int 	gFirstDesc;
char 	*bufPtr;
double 	gnpexp;
StrListPtr gTaxaSet;	
char  	*aTokenPtr;
int 	curTree=0;		/* index for the current tree description being parsed */
	
char	*nexError[2]= {
						"Error: Not a NEXUS file",		/* 0 */
						"Error opening NEXUS file"		/* 1 */
						};

struct NexDataType gNexData;	/* This is THE data structure for the NEXUS data */
struct NexDataType *gNexDataPtr;	/* This is THE data structure for the NEXUS data */

/**************************/
/**************************/

/* 
	Read a NEXUS file buffer and set up a global data structure containing everything. 
	See nexus.h for that data structure.
	Returns NULL on error. 
*/

void readNexusFile(char * theNexusFileBuffer)
{
	char *stemp;
	int c;
	int flag;
	int ix;
	long bufLength;

	struct NexDataType *nptr;
	
/*mallopt (M_DEBUG, 1);*/

	bufPtr=theNexusFileBuffer;	/* Initialize this global to beginning of the 
					buffer and will sweep through it until end of buffer */
	
	bufLength=strlen(theNexusFileBuffer);	

	if ( bufPtr != NULL )
		{
		aTokenPtr=nextToken();
		if(!isEqual(aTokenPtr,"#NEXUS"))
			{ 
			if (!aTokenPtr)   // NEW FEATURE SKIPS ANYTHING BEFORE '#NEXUS' in file
				{			 // NEW
			
			doError(nexError,0);
			return;		/* not a NEXUS file */
			
				}		//NEW
			else		// NEW
				aTokenPtr=nextToken();  // NEW
			}
		while (aTokenPtr=nextToken(), *aTokenPtr)
			{
			if (isEqual(aTokenPtr,"BEGIN"))
				{
				nextToken();if (!*aTokenPtr) return;
				stemp=DupStr(aTokenPtr);	/* get the block name and store in 'stemp'*/
				if (!stemp) 
					fatal ("Error reading block name");
				
				if (isEqual(aTokenPtr=nextToken(),";")) /* pop the terminating semicolon */
					{
					if (isEqual(stemp,"TAXA"))
							doTaxaBlock();
					else				
					if (isEqual(stemp,"CHARACTERS"))
							doCharsBlock();
					else				
					if (isEqual(stemp,"TREES"))
							doTreeBlock();
					else				
					if (isEqual(stemp,"RATES"))
							doRateBlock();
					else				
					if (isEqual(stemp,"R8S"))
							doRateBlock();
					else				
					if (isEqual(stemp,"BOOTSTRAP"))
							doBootBlock();
					else				
					if (isEqual(stemp,"SIMULATION"))
							doSimBlock();
					else				
					if (isEqual(stemp,"DATA"))
							doDataBlock();
					else 
						{  /* token is not a recognized block */
							doUnrecognizedBlock();
						}
					/*if (!*aTokenPtr)break;*/
					}
				free(stemp);				
				}
			}
		return;	/* normal return */
		}
	else
		{
		doError(nexError,1);
		return;
		}
	
}
/****************************************************************/
void doInteractive(void)
{
#define BUFSIZE 500
char Token[MAX_TOKEN_SIZE];
char inputBuffer[BUFSIZE],*pLast;
aTokenPtr=Token;
    for (;;)
	{
	printf("\nr8s>");
	fgets(inputBuffer,BUFSIZE,stdin);
#if STU
	strtoupper(inputBuffer);
#endif
	if (strlen(inputBuffer)>0)  /* if something in the buffer */
		{
		pLast=inputBuffer+strlen(inputBuffer)-1;
		if (*pLast != ';')
			{
			*(pLast+1)=';';
			*(pLast+2)='\0';
			}
		bufPtr=inputBuffer;
			doRateBlock();
		}
	    
	}
return;
}
void doCommandLineControl(char *inputBuffer)
{
char Token[MAX_TOKEN_SIZE];
char *pLast;
aTokenPtr=Token;
#if STU
strtoupper(inputBuffer);
#endif
if (strlen(inputBuffer)>0)  /* if something in the buffer */
		{
		pLast=inputBuffer+strlen(inputBuffer)-1;
		if (*pLast != ';')
			{
			*(pLast+1)=';';
			*(pLast+2)='\0';
			}
		bufPtr=inputBuffer;
		doRateBlock();
		}
  
    
}
static void doExecuteCommand(void)
{

char *theNexusFileBuffer, fnInput[FILENAME_MAX];
FILE * inStream =NULL;
strcpy(fnInput, aTokenPtr=nextToken());   /* set file name */
#if 0
if (gNexDataPtr)
	freeNexusStructure(gNexDataPtr);
gNexDataPtr=initialize_nexus();
#endif

// *** try this to allow multiple reads from different files... 
if (!gNexDataPtr)
	gNexDataPtr=initialize_nexus();
// ***
if (!gNexDataPtr)
	fatal("Failure to allocate nexus data structure in main.c");
if (!(inStream=fopen(fnInput,"r")) )
		    {
		    printf("file=%s\n", fnInput);
		    fatal("Error in file handling\n");
		    }
if (inStream)
	    {
	    theNexusFileBuffer=storeNexusFile(inStream);
	    readNexusFile(theNexusFileBuffer);
	    };
return;
}


static void printHelp(void)
{

	char * filename="r8s.helpfile";
	FILE* fpntr;
	char * buffer;
	buffer=NULL;
	if (  (fpntr=fopen(filename,"r")) )
		{
		buffer=slurpFile (fpntr, 10000);
		printf("%s\n", buffer);
		}
	else
	    printf("Failed to open file '%s'\n", filename);
	
	
	return;

}
	
/****************************************************************/
struct NexDataType * initialize_nexus(void)
{

struct NexDataType *p;

p=(struct NexDataType *)myMalloc(sizeof(struct NexDataType ));
if (!p)
    return NULL;
aTokenPtr=(char *)myMalloc(MAX_TOKEN_SIZE*sizeof(char));
if (!aTokenPtr)
    fatal ("Couldn't allocate aTokenPtr");
    
gTaxaSet=NULL;
gLabel=1;

p->TDList = newStrList(); /* initialize the list of tree descriptions */
p->TDLabelList = newStrList(); /* initialize the list of tree labels */
p->TaxaList = newStrList();
p->TransList=newStrList();
p->inTrees=NULL;
p->TaxSetNameList=NULL;
p->TaxSetLists=NULL;
p->excArray=NULL; /* can't be initialized further until we know num of chars */

p->isChars=0;
p->isTrees=0;
p->isTaxa=0;
p->isTranslate=0;	/*...flags for when these elements are read */
p->NTaxa=0;
p->NChars=0;
p->NumTrees=0;
p->matchchar='.';
p->gapchar='-';
p->missingchar='?';

p->RateBlockParms.NeighborPenalty=0;
p->RateBlockParms.PenaltyType=0;
p->RateBlockParms.checkGradient=1;
p->RateBlockParms.clampRoot=1;
p->RateBlockParms.isBS=0;
p->RateBlockParms.NReps=1;
p->RateBlockParms.seed=1;
p->RateBlockParms.RRtype=WULI;
p->RateBlockParms.npexp=2.0;
p->RateBlockParms.verbose=1;
p->RateBlockParms.num_restarts=1;
p->RateBlockParms.num_time_guesses=1;
p->RateBlockParms.num_rate_guesses=1;
p->RateBlockParms.smoothing=1.0;
p->RateBlockParms.showGradient=0;
p->RateBlockParms.showConvergence=0;
p->RateBlockParms.ftol=1e-8;
p->RateBlockParms.barrierTol=0.0001;
p->RateBlockParms.activeEpsilon=0.001;
p->RateBlockParms.maxIter=500;
p->RateBlockParms.maxBarrierIter=10;
p->RateBlockParms.initBarrierFactor=.25;
p->RateBlockParms.barrierMultiplier=0.10;
p->RateBlockParms.linminOffset=0.05;
p->RateBlockParms.contractFactor=0.1;
p->RateBlockParms.maxContractIter=10;
p->RateBlockParms.local_factor=0.01;
p->RateBlockParms.perturb_factor=0.01;
p->RateBlockParms.RatesAreGamma=0;
p->RateBlockParms.alpha=1.0;
p->RateBlockParms.numSites=1;
p->RateBlockParms.lengthFmt=0;  //'total'
p->RateBlockParms.roundFlag=1;  // round input branch lengths by default
p->RateBlockParms.clockFmt=0;
p->RateBlockParms.minRateFactor=0.05;
p->RateBlockParms.minDurFactor=0.001;

return p;
}
/****************************************************************/
/****************   BLOCK PROCESSING FUNCTIONS ******************/
/****************************************************************/


void doDataBlock(void)
{
	do 				/* need to put in error checking in case no DIMENSIONS statement */
		{
		aTokenPtr=nextToken();
		if (isEqual(aTokenPtr,"DIMENSIONS"))  
			doDimensions();
		if (isEqual(aTokenPtr,"FORMAT"))  
			doFormat();
		if (isEqual(aTokenPtr,"MATRIX"))
			doMatrixGeneral();

		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));
	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
return;
}
/****************************************************************/
void doUnrecognizedBlock(void)
{
	do 				
		{
		aTokenPtr=nextToken();
		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));
	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
return;
}
/****************************************************************/
void doCharsBlock(void)
{
	do 				/* need to put in error checking in case no DIMENSIONS statement */
		{
		aTokenPtr=nextToken();
		if (isEqual(aTokenPtr,"DIMENSIONS"))  
			doCharDimensions();
		if (isEqual(aTokenPtr,"FORMAT"))  
			doFormat();
		if (isEqual(aTokenPtr,"MATRIX"))
			doMatrix();

		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));
	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
return;
}
/****************************************************************/

static void doBootBlock(void)
{
StrListPtr fixedList=NULL;
char *buffer, *buffer1, *buffer2;
buffer=NULL;
buffer1=NULL;
buffer2=NULL;
	do 				/* need to put in error checking in case no DIMENSIONS statement */
		{
		aTokenPtr=nextToken();
		if (isEqual(aTokenPtr,"CLADE_SETS"))  
			doClade_Set_Command();
		if (isEqual(aTokenPtr,"MRCA"))  
			doMRCACommand();
		if (isEqual(aTokenPtr,"FIXED"))  
			fixedList=doFixedTaxaListCommand();
		if (isEqual(aTokenPtr,"BOOT"))  
			doBootCommand(fixedList, buffer);
		if (isEqual(aTokenPtr,"BOOTCHARS"))  
			doBootCharCommand(buffer);
		if (isEqual(aTokenPtr,"INCLUDE"))  
			buffer=doIncludeCommand();
		if (isEqual(aTokenPtr,"INCLUDE1"))  
			buffer1=doIncludeCommand();
		if (isEqual(aTokenPtr,"INCLUDE2"))  
			buffer2=doIncludeCommand();
		if (isEqual(aTokenPtr,"WEIGHTS"))  
			doB_Weights_Command(buffer,buffer1,buffer2);
		if (isEqual(aTokenPtr,"EFRON_WEIGHTS"))  
			doEFRON_Weights_Command(buffer);
	     	if (isEqual(aTokenPtr,"CLADE_CHECK"))  
			doCladeCheckCommand();
	     	if (isEqual(aTokenPtr,"SIMPLE_CLADE_CHECK"))  
			doSimpleCladeCheckCommand();

		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));
	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
return;
}
/****************************************************************/
static void  doClade_Set_Command()
{
/* for each tree, creates a clade list and stores a pointer to that list
   in the tree structure */
       
	int matchCount;
	StrListPtr allTaxaList;
	double ovl, ovlmax, ovlmin;
	PtrList SetList;
	PtrList lnode, pnode, focal_cladeSet, cur_cladeSet, fnode, clnode;
	TREE thisTree, focal_tree, curTree;
	Set fclade, cclade;
	int ntaxa, ntrees;
	while (!isEqual(aTokenPtr=nextToken(),";")); /* just the one word command */
	if (gNexDataPtr->isTrees)
		{
		lnode=gNexDataPtr->inTrees;
		allTaxaList=newStrList();
		thisTree=lnode->item;
		TreeToTaxaList(thisTree->root,allTaxaList); /* important to do this ONCE only! */
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			thisTree->cladeSet = Tree2CladeSet(thisTree, allTaxaList);
			}
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			SetList=thisTree->cladeSet;
			printCladeSets(SetList);
			}
		}
	freeStrList(allTaxaList);

#if 1 
/*      do the fuzzy bootstrap */

	if (gNexDataPtr->isTrees)
		{
		lnode=gNexDataPtr->inTrees;
		ntrees=pLengthList(lnode);
		if (ntrees >= 2) /* must have a focal tree and one other to continue! */
		    {
		    focal_tree=lnode->item;
		    focal_cladeSet=focal_tree->cladeSet;
		    printf("\nFuzzy Clade Analysis\nMinimax 	 BP  \t\t   Clade\n-----------------------------------\n");
		    LISTLOOP(focal_cladeSet) /* loop through all the clades in the focal tree */
			{
			fclade=focal_cladeSet->item;
			lnode=(gNexDataPtr->inTrees)->next;
			ovlmin=+9999.999;
			matchCount=0;
			 LISTLOOP(lnode)  /* loop through the other trees ...*/
			    {
			    curTree=lnode->item;
			    cur_cladeSet=curTree->cladeSet;
			    ovlmax=-9999.999;
			    LISTLOOP(cur_cladeSet) /*...checking each of their clades */
				{
				cclade=cur_cladeSet->item;
				/* test_set(fclade, cclade); */
				ovl=set_overlap(fclade, cclade);
				/*printf("overlap... %f\n", ovl);*/
				if (ovl>ovlmax) ovlmax=ovl;
				if (ovl==1.0)++matchCount; /* there was a clade matching the focal clade */   
				}
			    /*printf("max overlap for this clade = %f\n", ovlmax);*/
			    if (ovlmax<ovlmin) ovlmin=ovlmax;
			    }
			printf("%f\t%f\t", ovlmin, matchCount/(float)(ntrees-1));   
			print_set(fclade);
			}
			
		    }
		}
#endif	
    
}
static void doCladeCheckCommand()
/** Reads a list of taxa and checks to see if that group is a clade on all trees.
Reports the proportion of trees in which this group is a clade.
ADDED.  Lots of junk to do Efron 1996 bootstrap.  A real hack.  We generate sets of NUMBER_IN_EFRON_SET trees. The
first in every set is the P(j),  

To avoid a bias in which we always take the tree-weights over the boundary in the direction toward pj,  we
filp a coin to choose between that value of w and the next lower one,  which has pj on the other side of the
boundary toward pcent. **/
{

#define NUMBER_IN_EFRON_SET 26

	int isClade[NUMBER_IN_EFRON_SET];
	StrListPtr aTaxaList, txPtr, nLptr;
	PtrList nodeList, mrcaPtr;
	PtrList lnode;
	TREE thisTree;
	NODETYPE *mrca, *node;
	int i, ix=0, kix=-1,nList,counter=0, first_of_set, jix,flipCount,
	    last_of_set, watch_this_block=0, counter2=0, coin, cindex;
	aTaxaList=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	
			appendStrList(aTaxaList,aTokenPtr); /* store the label */
	nList=lengthList(aTaxaList);
	if (nList < 2)
	    fatal("Must have at least two names in CLADE_CHECK command");
	printf("cat phase2_header >> paup_phase2\n");
	if (gNexDataPtr->isTrees)
		{
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			++ix;
			++kix;
			if ( (ix-1)/NUMBER_IN_EFRON_SET == (float)(ix-1)/NUMBER_IN_EFRON_SET )
			    {
			    kix=0;
			    printf("# .................................\n");
			    }
			printf("# Tree %i:Specified group IS ",ix);
			thisTree=lnode->item;
			if (group_a_clade(thisTree->root, aTaxaList))
				{
				++counter;
				printf("a clade\n");
				isClade[kix]=1;

				}
			else
				{
				printf("NOT a clade\n");
				isClade[kix]=0;
				}
			if ( (ix)/NUMBER_IN_EFRON_SET == (float)(ix)/NUMBER_IN_EFRON_SET ) /* last of set */
			    	{
					flipCount=0;
					for (jix=0;jix<NUMBER_IN_EFRON_SET-1;jix++)
						if (isClade[jix] != isClade[jix+1])
							++flipCount;
					if (flipCount==1)
					  {
					  if (isClade[0]==0)
					    {
					    for (jix=0;jix<NUMBER_IN_EFRON_SET-1;jix++)
						if (isClade[jix]) /* here is the transition */
							{
							if (myRand()>0.5)coin=1;else coin=0;
							if (coin) 
								cindex=ix-NUMBER_IN_EFRON_SET+jix;
							else 
								cindex=ix-NUMBER_IN_EFRON_SET+jix+1;
							printf("agrep -d \';\' \'Weight set %li:\' theWsearchWeights >> paup_phase2\n", cindex);
				/* printf("cat PHASE2_INCLUDE >> paup_phase2\n");*/
							break;
							}
					    }
					  else
					    printf("# ONE CHANGE BUT WRONG DIRECTION\n");
					  }
					if (flipCount==0)
						printf("# NO changes\n");
					if (flipCount>1)
						printf("# WARNING! Uneven boundary conditions\n");

				}
			}
		}
	    printf("# Proportion (%i) of %i trees with group monophyletic (EFRON):%f\n",counter2, ix, 
			(float)NUMBER_IN_EFRON_SET*counter2/ix);
	    printf("cat phase2_tailer >> paup_phase2\n"); /*wrapper for nexus syntax */

	return;
    
    
}
static void doSimpleCladeCheckCommand()
/** Reads a list of taxa and checks to see if that group is a clade on all trees.
Reports the proportion of trees in which this group is a clade.

Writes a script in shell language to allow agreping from the weights file **/

{
	StrListPtr aTaxaList, txPtr, nLptr;
	PtrList nodeList, mrcaPtr;
	PtrList lnode;
	TREE thisTree;
	NODETYPE *mrca, *node;
	int i, ix=0, nList,counter=0, first_of_set, 
	    last_of_set, watch_this_block=0, counter2=0, coin, cindex;
	aTaxaList=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	
			appendStrList(aTaxaList,aTokenPtr); /* store the label */
	nList=lengthList(aTaxaList);
	if (nList < 2)
	    fatal("Must have at least two names in SIMPLE_CLADE_CHECK command");

	printf("echo \"#nexus\" > theNULLHweights\n");
	printf("echo \"begin bootstrap;\" >> theNULLHweights\n");
	printf("echo \"include file=PHASE1A_INCLUDE;\" >> theNULLHweights\n");
	printf("#Check for group monophyly...\n");
	if (gNexDataPtr->isTrees)
		{
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			++ix;
			printf("# Tree %i:Specified group IS ",ix);
			thisTree=lnode->item;
			if (group_a_clade(thisTree->root, aTaxaList))
				{
				++counter;
				printf("a clade\n");
				}
			else
				{
				printf("NOT a clade\n");
				printf("agrep -d \';\' \'Weight set %li:\' the_phase1_weights >> theNULLHweights\n", ix);
				}
			}
		}
	printf("# Proportion (%i) of %i trees with group monophyletic:%f\n",counter, ix, 
			(float)counter/ix);
	printf("echo \";end;\" >> theNULLHweights\n");

	return;
    
    
}
static void doB_Weights_Command(char *buffer, char *buffer1,char *buffer2)

/* For every weight statement in the input NEXUS file, we generate NREPS
new weight statements to allow further bootstrapping.  The weights can be real,
for example in EFRON resampling, which requires multinomial resampling*/

#define NREPS 100
#define MAX_WEIGHT_ARRAY 2500
{
	float XweightArray[MAX_WEIGHT_ARRAY];
	float weight;
	int character, j;
	char * dummy;
	for (j=0;j<MAX_WEIGHT_ARRAY;j++)
	    {
	    XweightArray[j]=0.0;
	    }
	while (!isEqual(aTokenPtr=nextToken(), ";"))
	    {
	    if (isEqual(aTokenPtr, ","))
		aTokenPtr=nextToken(); /* skip commas */
	    weight=strtod(aTokenPtr, &dummy);
	    aTokenPtr=nextToken();
	    if (!isEqual(aTokenPtr,":")) fatal("Improperly formatted b_weights statement");
	    aTokenPtr=nextToken();
	    character=strtod(aTokenPtr, &dummy);
	    if (character >=MAX_WEIGHT_ARRAY)
		fatal("Too many characters in B_WEIGHT_COMMAND: recompile with larger array");
	    XweightArray[character-1]=weight;
	    }	
	/* printf("TEST OF WEIGHTS COMMAND\n");
	for (j=0;j<MAX_WEIGHT_ARRAY;j++)
	    {
	    if ((j>0)&& ((j/10)==(j/10.0)))
		printf("\n");
	    
	    printf("%4.2f:%i, ", XweightArray[j],j+1);
	    }*/
	bshuf3(XweightArray, character, NREPS, buffer1,buffer2); /* assumes the last character read is the highest; ie. weights
		    read sequentially--violates NEXUS format,  but this command is only in my BOOT block anyway */
	if (buffer)
	    printf("%s\n", buffer);
	return;    
    
}
static void doEFRON_Weights_Command(char *buffer)

/* For every INTEGER weight statement in the input NEXUS file, we generate a new set of REAL weights on the
boundary between R1 and R2; see EFRON 1996 */

{
	int weightArray[MAX_WEIGHT_ARRAY];
	float weight;
	long index;
	int character, j;
	char * dummy;
	for (j=0;j<MAX_WEIGHT_ARRAY;j++)
	    {
	    weightArray[j]=0.0;
	    }
	while (!isEqual(aTokenPtr=nextToken(), ";"))
	    {
	    if (isEqual(aTokenPtr, ","))
		aTokenPtr=nextToken(); /* skip commas */
	    weight=strtod(aTokenPtr, &dummy);
	    aTokenPtr=nextToken();
	    if (!isEqual(aTokenPtr,":")) fatal("Improperly formatted b_weights statement");
	    aTokenPtr=nextToken();
	    character=strtod(aTokenPtr, &dummy);
	    if (character >=MAX_WEIGHT_ARRAY)
		fatal("Too many characters in EFRON_WEIGHT_COMMAND: recompile with larger array");
	    weightArray[character-1]=weight;
	    }
/* NB!  careful to make sure that efron1996 actually generates the right number(NUMBER_IN_EFRON_SET)  of increments ! */	
	efron1996(weightArray,character,NUMBER_IN_EFRON_SET,buffer, &index);

	return;    
    
}

static void doMRCACommand(void)

/** assigns an internal name to the MRCA of a set of taxa: 
Usage MRCA new_internal_name taxon1 ...taxonN ; if only one taxonname is given,  then 
assume it is an internal name and replace it with the newname
 **/

{
	StrListPtr aTaxaList, txPtr, nLptr;
StrListPtr mrcaTaxa;
	PtrList nodeList, mrcaPtr;
	PtrList lnode;
	TREE thisTree;
	NODETYPE *mrca, *node;
	char *new_internal_name, *old_internal_name;
	int i, ix=0, nList;
	gTaxaList=txPtr=aTaxaList=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	
			{
			appendStrList(aTaxaList,aTokenPtr); /* store the label */
			}
	nList=lengthList(aTaxaList);
	if (nList < 2)
		{
	    	doGenericAlert("Must have at least two names in MRCA command");
		return;
		}
	new_internal_name=aTaxaList->s; /* this is the first name in the list */
	nLptr=aTaxaList->next; /*points to taxon names*/
	if (gNexDataPtr->isTrees)
		{
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			if (nList == 2)
			    {
			    node=find_taxon_name(thisTree->root,nLptr->s);
			    if (node)
				{
				printf("Redefining node name: %s to %s\n",node->taxon_name,new_internal_name);
				setNodeName(node, new_internal_name);
				}
			    else
				{
				doGenericAlert("BAD MRCA COMMAND: Taxon name misspelled or not on tree");
				return;
				}
			    }
			else
			    {
			    mrca=MRCA(thisTree->root, nLptr);
			    if (mrca)
				{
				if (mrca->taxon_name)
				    if(!isEqual(mrca->taxon_name,"")) /* careful, name initialized to "" */
					{
					doGenericAlert("MRCA is overwriting an existing node name");
					printf("[** The overwritten node is %s **]\n",mrca->taxon_name);
					}
				setNodeName(mrca, new_internal_name);
				printf("Defining clade name: %s\n",new_internal_name);
				}
			    else
				{
				doGenericAlert("BAD MRCA COMMAND: Taxon name misspelled or not on tree");
				return;
				}
			    }
			}
		}


	return;
    
    
}

static void doBootCommand(StrListPtr fixedList, char* buffer)

/* Process the taxon bootstrap command, write the appropriate
NEXUS syntax to delete taxa and, also write stuff from a character
buffer that might have been included with the 'include' command */

{
	char  *dummy;
	int sample[MAX_TAXON_ARRAY], fixed[MAX_TAXON_ARRAY], included[MAX_TAXON_ARRAY], 
			    ntaxa=0, nrandom=0, nfixed=0,nsample=0,
			nstart=0,nstop=0,nstart2=0,nstop2=0,nrandom2=0;
	long aSeed=1, nreps, i, j, k;
	for (k=0;k<MAX_TAXON_ARRAY;k++)
		fixed[k]=0; /* necessary default for later 
				calls to 'taxon_sample' */
	if (fixedList) /* if there are fixed taxa */
	    {
	    nfixed=lengthList(fixedList);
	    for (k=0;k<nfixed;k++)
		fixed[k]=strtod(getkthStr(fixedList, k+1), &dummy);
		    /* internal representation of taxon ids is on 1..n */
	    }
		
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NSTART"))
			nstart=strtod(LocalToken,&dummy);
		if (parse_assignment2("NSTOP"))
			nstop=strtod(LocalToken,&dummy);
		if (parse_assignment2("NSTART2"))
			nstart2=strtod(LocalToken,&dummy);
		if (parse_assignment2("NSTOP2"))
			nstop2=strtod(LocalToken,&dummy);
		if (parse_assignment2("NREPS"))
			nreps=strtod(LocalToken,&dummy);
		if (parse_assignment2("NRANDOM"))
			nrandom=strtod(LocalToken,&dummy);
		if (parse_assignment2("NRANDOM2"))
			nrandom2=strtod(LocalToken,&dummy);
		if (parse_assignment2("NTAXA"))
			ntaxa=strtod(LocalToken,&dummy);
		if (parse_assignment2("SEED"))
			aSeed=strtod(LocalToken,&dummy);
		}
	nsample=nfixed+nrandom+nrandom2; /* total number to be in sample */
	srand(aSeed);
	printf("begin paup;\n");
	for (k=1;k<=nreps;k++)
	    {
	    for (i=0;i<ntaxa;i++)
		included[i]=0;
	    taxon_sample(ntaxa,nfixed, nrandom, fixed, sample,
			nstart,nstop,nstart2,nstop2,nrandom2);
	    for (j=0;j<nsample;j++)
		    included[sample[j]-1]=1;
	    printf("[The taxon sample is:");
	    for (j=0;j<ntaxa;j++)
		if (included[j]) 
		    printf("%i ", j+1); 
	    printf("]\n");
	    printf("delete ");
	    for (j=0;j<ntaxa;j++)
		if (!included[j]) 
		    printf("%i ", j+1); /* +1 is to reconvert back to external
					representation of taxon ids on 1..n */
	    printf("/prune=yes;\n");
	    if (buffer)
		    printf("%s\n", buffer);

	    }
	printf("end;\n");
	return;
}
static void doBootCharCommand(char* buffer)

/* Process the character bootstrap command, write the appropriate
NEXUS syntax to weight characters and, also write stuff from a character
buffer that might have been included with the 'include' command */


{
	char  *dummy;
	double u, u1=0.0, u2=0.0, a;
	int *weightArray;
	float * pMean,pSum=0.0;
	long aSeed=1, nreps=0, i, j, k,nchars=0, index=0,ix;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NREPS"))
			nreps=strtod(LocalToken,&dummy);
		if (parse_assignment2("NCHARS"))
			nchars=strtod(LocalToken,&dummy);
		if (parse_assignment2("SEED"))
			aSeed=strtod(LocalToken,&dummy);
		}
	srand(aSeed);
	if (nchars>0)
		{
		pMean=(float *)myMalloc((nchars)*sizeof(float));
		for (ix=1;ix<nchars;ix++)
			pMean[ix]=0.0;
		weightArray=(int *)myMalloc(nchars*sizeof(int));
		printf("begin paup;\n");
		for (k=1;k<=nreps;k++)
		    {
		    u1=0.0;
		    u2=0.0;
		    ++index;
		    bshuf2(weightArray,nchars);
		    printf("\n\n[******************************]\n");
		    printf("[*** Bootstrap replicate %i ***]\n\n", k);
		/*    for (j=0;j<nchars;j++)
			    {
			    u=weightArray[j]-1.0; 
			    u1+=u*u*u;
			    u2+=u*u;
			    } *//* for EFRON 96 algorithm */
		    a= (1/6.0)*u1/pow(u2, 1.5);
		    printf("[Weight set %li:][w=1.0][a=%f]weights ", index, a);
		    for (j=0;j<nchars-1;j++)
			    {
			    if ((j>0)&& ((j/10)==(j/10.0)))
				printf("\n");
			    printf("%i:%i, ", weightArray[j],j+1);
			    }
		    printf("%i:%i;\n", weightArray[nchars-1],nchars);
		    if (buffer)
			    printf("%s\n", buffer);
		   for (j=0;j<nchars;j++)
			    {
			    pMean[j]+=weightArray[j]/(float)nreps;
			    }
	
	/*	    if (EFRON1996)
			efron1996(weightArray,nchars,NUMBER_IN_EFRON_SET,buffer, &index);*/

		    }
		printf("[Mean vector\n");
		for (j=0;j<nchars;j++)
			{
			if ((j>0)&& ((j/10)==(j/10.0)))
			    printf("\n");
			printf("%f:%i, ", pMean[j],j+1);
			pSum+=pMean[j];
			}
		printf(" sum=%f]\n",pSum);
		printf("end;\n");
		myFree(weightArray);
		}
	return;
}
static void efron1996(int *weightArray,int nchars,int num_points, char *buffer, 
    long *index)

/* writes PAUP code to implement part of Efron et al., 1996 boot algorithm.
 * This component writes commands that generate a set of weight statements
 * corresponding to the search for the boundary vectors.  Calculates proper
 * weights for all points,  w,  such that w is an element of [0,xinc, 2*xinc, 
 * ..., 1].  For each point the p vectors wP(j)+(1-w)P(cent) are calculated 
 * (see bottom left of p. 7089 of paper).  A row of weights is printed along
 * with any commands stored in the buffer from a previous include.
 *
 * 'weightArray' contains the bootstrap weight vector, P(j) on entry.
 * Spits out NUMBER_IN_EFRON_SET rows of weights.
 */

{
double w, u, u1=0.0, u2=0.0, a, pj,xinc;
int j,k;
if (num_points<2) fatal("Too few points in efron1996");
xinc=1/(num_points-1.0);
printf("begin paup;\n");
for (k=1;k<=num_points;k++) 
    {
    w=1-(k-1)*xinc;
    u1=0.0;
    u2=0.0;
    ++(*index);
    printf("\n[Weight set %li:][w=%f] ", *index, w);
    for (j=0;j<nchars;j++)
	    {
	    u= (1-w)+w*weightArray[j]-1.0;
	    u1+=u*u*u;
	    u2+=u*u;
	    }
    a= (1/6.0)*u1/pow(u2, 1.5);
    printf(" [a=%f] ", a);
    printf("weights ");
    for (j=0;j<nchars-1;j++)
	    {
	    if ((j>0)&& ((j/30)==(j/30.0)))
		printf("\n");
	    pj=(1-w)+w*weightArray[j];
	    
	    printf("%4.2f:%i, ", pj,j+1);
	    }
    pj=(1-w)+w*weightArray[nchars-1];
    printf("%f:%i;\n", pj,nchars);
    if (buffer)
	    printf("%s\n", buffer);
    printf("\n");
    } 
printf("end;\n");
return;
}


static char * doIncludeCommand(void)

/* Reads the specified file and passes it to the 'boot' command so 
that its text can be printed for each boot replicate.
NB!  The include must be executed BEFORE THE BOOT command.
AARRPHHH: demands that filenames be in caps!*/

{
	FILE* fpntr;
	char  *filename, *buffer;
	buffer=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("FILE"))
			filename=DupStr(LocalToken);
		}
	if (  (fpntr=fopen(filename,"r")) )
		{
		buffer=slurpFile (fpntr, 10000);
		}
	else
	    printf("Failed to open file '%s'\n", filename);
	return buffer;
}

/*
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		appendStrList(gNexDataPtr->TaxaList,aTokenPtr);
		}

*/
static StrListPtr doFixedTaxaListCommand(void)
{
	StrListPtr aTaxaList; /* probably just numbers, but treat them
		as strings */

	aTaxaList=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	
			appendStrList(aTaxaList,aTokenPtr); /* store the label */
	return aTaxaList;

}



/****************************************************************/

void doTaxaBlock(void)
{
	do 				/* need to put in error checking in case no DIMENSIONS statement */
		{
		aTokenPtr=nextToken();
		if (isEqual(aTokenPtr,"DIMENSIONS"))  
			doTaxDimensions();
		if (isEqual(aTokenPtr,"TAXLABELS"))
			doTaxLabels();

		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));
	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
return;
}
/**************************************************************/
void doTreeBlock(void)
{
	do 
		{
		aTokenPtr=nextToken();
		if (isEqual(aTokenPtr,"TREE") || isEqual(aTokenPtr,"UTREE") )  
									/* process TREE command */
			doTreeCommand();
		if (isEqual(aTokenPtr,"TRANSLATE") )  
									/* process TREE command */
			doTranslateCommand();
			
		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));

	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
if (gNexDataPtr->isTrees)
	{
	/*printf("[Number of trees read = %i]\n",gNexDataPtr->NumTrees);		
	print_tree_list(gNexDataPtr->inTrees);*/
	}
return;
}
/**************************************************************/
void doSimBlock(void)
{
	do 
		{
		aTokenPtr=nextToken();
		if (isEqual(aTokenPtr,"SIMULATE") )  
									/* process TREE command */
			doSimCommand();
		}  while (!isEqual(aTokenPtr,"END")  &&
						(!isEqual(aTokenPtr,"ENDBLOCK") ));

	aTokenPtr=nextToken();
	if (!isEqual(aTokenPtr,";"))
		doGenericAlert("Block not terminated with semicolon");
return;
}
/****************************************************************/

void doDivTimeCommand(void)
{
extern int powellMode;
extern double gGamma_c;
extern double gGamma_b;
extern int gisConstrained, gVarMinFlag,gEstRoot,gFloatRoot;
	double (*obj_func_array[10])(double[]); /* array of pointers to the various
		objective functions...indexed by 'method' which is set up below */
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD ,*method_string;
double  EstMult=1,PrdMult=1,cutoff=2.0;
static	int method=LaF,
	    algorithm=POWELL;
static	long iTree=0;	
	long numTrees;
	int j,success,crossv=0,crossv2=0,cvSample=0,nRates=1,maxBisect=20,fossilFlag=0,fossilFixedFlag=0;
	NODETYPE *root, *found_node;
static	double cvStart=0.0,cvInc=1.0,cvNum=1;
	int confidence=0;
StackPtr S;

	powellMode=1;
	method_string=NULL; /*default*/
/*
	obj_func_array[LaF]=objLangFitch;
	obj_func_array[NP]=objNP;
	obj_func_array[GAMMA]=objGamma;
	obj_func_array[PENLIKE]=objPenLike;
	obj_func_array[PENLIKET]=objPenLikeT;
*/

#define POWELL_STACK_SIZE 35
gFStack=newStack(POWELL_STACK_SIZE);
gPStack=newStack(POWELL_STACK_SIZE);


	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("ALGORITHM"))
			{
			if (isEqual(LocalToken,"POWELL"))
				algorithm=POWELL;
			if (isEqual(LocalToken,"QNEWT"))
				algorithm=QNEWT;
			if (isEqual(LocalToken,"TN"))
				algorithm=TN;
			}
		if (parse_assignment2("METHOD"))
			{
			method_string=DupStr(LocalToken);
			}
		if (parse_assignment2("TAXON"))
			{
			taxon=DupStr(LocalToken);
			}
		if (parse_assignment2("TREE"))
			{
			iTree=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("NRATES"))
			nRates=strtod(LocalToken,&dummy);
		if (parse_assignment2("MAXBISECT"))
			maxBisect=strtod(LocalToken,&dummy);
		if (parse_assignment2("CONFIDENCE"))
			{
			if (isEqual(LocalToken,"YES"))
				confidence=1;
			else
				confidence=0;
			}
		if (parse_assignment2("CROSSV"))
			{
			if (isEqual(LocalToken,"YES"))
				crossv=1;
			else
				crossv=0;
			}
		if (parse_assignment2("FOSSILFIXED"))
			{
			if (isEqual(LocalToken,"YES"))
				fossilFixedFlag=1;
			else
				fossilFixedFlag=0;
			}
		if (parse_assignment2("FOSSILCONSTRAINED"))
			{
			if (isEqual(LocalToken,"YES"))
				fossilFlag=1;
			else
				fossilFlag=0;
			}
		if (parse_assignment2("CUTOFF"))
			cutoff=strtod(LocalToken,&dummy);
		if (parse_assignment2("CVSTART"))
			cvStart=strtod(LocalToken,&dummy);
		if (parse_assignment2("CVINC"))
			cvInc=strtod(LocalToken,&dummy);
		if (parse_assignment2("CVNUM"))
			cvNum=strtod(LocalToken,&dummy);
		if (parse_assignment2("CVSAMPLE"))
			cvSample=strtod(LocalToken,&dummy);
		
		}
	if (method_string!=NULL)
		{
		if (isEqual(method_string,"PL"))method=PENLIKE; 
		if (isEqual(method_string,"LF"))method=LaF; 
		if (isEqual(method_string,"NP") || isEqual(method_string,"NPRS"))
		    {
		    method=NP;
		    gVarMinFlag=0; /* this is a kludgy way to distinguish between two flavors
					of NPRS alogorithm */
		    }
		if (isEqual(method_string,"NPVAR"))
		    {
		    method=NP;
		    gVarMinFlag=1;
		    }

		}

		if ((method==LaF)  && (nRates > 1) )method=LFLOCAL; /* temp..I can just use this always for LF */



	if (gNexDataPtr->isTrees)
	    {
	      lnode=gNexDataPtr->inTrees;
	      if (iTree>0) /* a specific tree was indicated */
		{
		if (iTree > pLengthList(lnode))
			{
			doGenericAlert("Invalid tree specified");
			return;
			}
		else
			{
			thisTree=(pListgetkthNode(lnode,iTree))->item;
			doObjFunc(thisTree,method,nRates,algorithm,&success);
			}
		}
	      else
		{
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			if (crossv)
				{
				if (fossilFlag)
					doFossilCrossV(thisTree,method,nRates,algorithm,cvStart,cvInc,cvNum,cvSample);				
				else if (fossilFixedFlag)
					doFossilCrossVfixed(thisTree,method,nRates,algorithm,cvStart,cvInc,cvNum,cvSample);				
				else
					doCrossV(thisTree,method,nRates,algorithm,cvStart,cvInc,cvNum,cvSample);
				}
			else
				{
				doObjFunc(thisTree,method,nRates,algorithm,&success);
				if (confidence)
					{
					doConfidence(thisTree,taxon,method,nRates,algorithm,cutoff,maxBisect);
					}
				}
			}
		}
	      if (method==GAMMA)
				{
				thisTree->est_b=gGamma_b;
				thisTree->est_c=gGamma_c;
				}
	    }
	if (method_string)
		myFree(method_string);
	freeStack(gFStack);
	freeStack(gPStack);
	return;						
}

/***************/
static void
doConfidence(TREE T,char * nodeName,int method,int nRates,int algorithm,double cutoff,int JMAX)

/*
Find the confidence interval on an estimated node time. Construct the interval by
finding the values of that node time at which the likelihood drops by an amount 'cutoff'.
This is done by examining a range of possible times roughly between an upper and lower bound determined
from the age constraints and fixed node times. The focal node is fixed at various times across this
range, and the search is restarted (cf Cutler, MBE, 2000) estimating all other parameters as before. Rather
than search the whole range, a bisection strategy is used with an NRC function.


nodeName -	Determine confidence interval for this node
cutoff	--	Target of the search is (Max Like - cutoff)
JMAX	--	Maximum number of bisections allowed
*/


{
int i,maxPts=10,success,j;
NODETYPE *n;
double upper,lower,R,t,tSoln,factor,solnObj,targetObj,bump,low,high,tLow,objLow,tHigh,objHigh;
double x1,x2,dx,f,fmid,xmid,rtb,xacc;
solnObj=T->obj;	/* value of the obj function at the solution */
targetObj=solnObj-cutoff;
if (!(n=find_taxon_name(T->root,nodeName)))
	{
	doGenericAlert("Failed to find node name in 'confidence'");
	return;
	}
if (!isFree(n))
	{
	doGenericAlert("Cannot estimate confidence limit on FIXED node");
	return;
	}
tSoln=n->time;			/* save the estimated age of node */	
upper=nodeUpperBound(n);
if (upper>=1e20)
	upper=2*n->time;	/* if no upper bound, arbitrarily put it at 2X the node's age, but
					we'll check to see if this accomodates results */
lower=nodeLowerBound(n);
factor=0.9;			/* let's squeeze the search interval by this amount to prevent bumping
					against bounds! */
R=upper-lower;
bump=R*(1-factor)/2;
lower+=bump;
upper-=bump;
low=lower;
high=tSoln;
xacc=(upper-lower)*0.01;
x1=low;
x2=tSoln;
	{/* modifed from NRC 'rtbis' */
	n->free=0;
	n->time=x1;
	doObjFunc(T,method,nRates,algorithm,&success);
	f=T->obj - targetObj;
	fmid=solnObj - targetObj;
	if (f*fmid >= 0.0)
		{
		doGenericAlert ("Confidence search failed: no crossover point");return;
		}
	rtb = f < 0.0 ? (dx=x2-x1,x1):(dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++)
		{
		xmid=rtb+(dx*=0.5);
		n->time=xmid;
		doObjFunc(T,method,nRates,algorithm,&success);
		fmid=T->obj - targetObj;
		if (fmid <= 0.0) rtb=xmid;
		if (/*  fabs(dx) < xacc || fmid == 0.0 */ fabs(fmid) < 0.1) /* my termination criterion! */
			{
			/* printf ("**** Lower t = %f objective function - targetObj= %f [iters=%i]\n",rtb,fmid,j);*/
			tLow=rtb;
			objLow=fmid+targetObj;
			break;
			}
		}
	if (j>=JMAX)
		doGenericAlert("Confidence search failed");
	}
x2=tSoln;	/* lazy-ass copy of the previous code! */
x1=upper;
	{/* modifed from NRC 'rtbis' */
	n->free=0;
	n->time=x1;
	doObjFunc(T,method,nRates,algorithm,&success);
	f=T->obj - targetObj;
	fmid=solnObj - targetObj;
	if (f*fmid >= 0.0)
		{
		doGenericAlert ("Confidence search failed: no crossover point");return;
		}
	rtb = f < 0.0 ? (dx=x2-x1,x1):(dx=x1-x2,x2);
	for (j=1;j<=JMAX;j++)
		{
		xmid=rtb+(dx*=0.5);
		n->time=xmid;
		doObjFunc(T,method,nRates,algorithm,&success);
		fmid=T->obj - targetObj;
		if (fmid <= 0.0) rtb=xmid;
		if (/*  fabs(dx) < xacc || fmid == 0.0 */ fabs(fmid) < 0.1) /* my termination criterion! */
			{
			/* printf ("**** Higher t = %f objective function - targetObj= %f [iters=%i]\n",rtb,fmid,j);*/
			tHigh=rtb;
			objHigh=fmid+targetObj;
			break;
			}
		}
	if (j>=JMAX)
		doGenericAlert("Confidence search failed");
	}
printf("\nConfidence interval for node %s using cutoff value of %f\n",n->taxon_name,cutoff);
printf("Point\t\tAge\t\tObj\n");
printf("Lower\t\t%6.2f\t\t%6.2f\n",tLow,objLow);
printf("Higher\t\t%6.2f\t\t%6.2f\n",tHigh,objHigh);
printf("Soln\t\t%6.2f\t\t%6.2f\n",tSoln,solnObj);



n->free=1;
return;
}
/***************/

static float 
doCrossV(TREE tree, int method,int nRates,int algorithm,double cvStart,double cvInc,double cvNum, int cvSample)

/*  
	Does a cross validation analysis in which we 
	(1)sequentially remove each tip (leaving the tip's ancestor in place),
	(2)do a full estimation on remaining subtree, 
	(3)then calculate a prediction error for that removed terminal branch. 
	(4) Then puts the terminal back on the tree.

	If the method is LaF or NP then one round of CV is invoked.
	If the method is PENLIKE, then analysis is repeated with the smoothing parameter chosen from a range
	from [cvStart, cvStart+cvInc, ...,cvStart+cvInc*(cvNum-1)].
	If cvSample==0 then we cross validate on all the taxa; if cvRep>0, we randomly sample that many taxa
		and use only those taxa. We use the same random sample for ALL smoothing levels, however!


*/
{
char *Result, *Good="Good", *Failed="Failed";
int i,j,k,success,collFlag=0,ntips,verbose,overallGood=1;
double * cvScore,*chiSqScore, cvSum,chiSq,chiSqSum,*cvTotalScore,*cvTotalScoreChiSq,bestChiSq,bestSmooth;
int * sample, *cvResult, *cvResultFinal;
int numSuccess, numFail,bestJ;
double smooth;
PtrList tipNodeList;
NODETYPE *CVNode, *saveAncNode;
ntips=numdesc(tree->root);
cvResult = (int *)myMalloc((ntips+1)*sizeof(int));
cvResultFinal = (int *)myMalloc((cvNum+1)*sizeof(int));
cvScore = vector(1,ntips);
chiSqScore = vector(1,ntips);
cvTotalScore=vector(1,cvNum);
cvTotalScoreChiSq=vector(1,cvNum);
tipNodeList = pNewList();
TreeToTaxaPtrList(tree->root,tipNodeList); /* get a list of all the tip nodes */

/* NOTES: This will not work when tips have ages > 0. In that case, we must enforce constraints on that tips
	ancestor node; otherwise, once it is pruned, that ancestor might be inferred to be younger than the 
	pruned tip, causing all hell to break loose when doing the predicted value */

if (cvSample>0)
	{
	srand(gNexDataPtr->RateBlockParms.seed);	/* sets up random seed */
	sample=taxon_sample_simple(ntips,cvSample);
	ntips=cvSample;
	}

if (method==LaF || method==NP)
	cvNum=1;	/* just overrides the default value of cvNum for these methods, forcing them to only go once */
verbose=gNexDataPtr->RateBlockParms.verbose;
if (verbose > 0)
	printf("Begin cross-validation analyses...\n");
for (j=0;j<cvNum;j++)
	{
	numSuccess=0;
	numFail=0;
        cvSum=0.0;
	chiSqSum=0.0;
	smooth=pow(10.0,j*cvInc+cvStart);
	gNexDataPtr->RateBlockParms.smoothing=smooth;
	for (k=1;k<=ntips;k++)
		{
		if (cvSample>0)
			i=sample[k-1];
		else
			i=k;
		CVNode=(NODE)(pListgetkthNode(tipNodeList,(long)i)->item);
		saveAncNode = CVNode->anc; // this should be ok as long as CVNode is never the root, which it should node be from prev line
		RemoveTaxonLvAnc(CVNode);
		gNexDataPtr->RateBlockParms.verbose=0; 
	/* suppress all output from the actual optimization run...may want to allow it for debugging though! */
		doObjFunc(tree,method,nRates,algorithm,&success);
		gNexDataPtr->RateBlockParms.verbose= verbose; /* restore output verbositude to current value */
		//AddChild(CVNode->anc,CVNode); /* important to reattach before next call. [!this was buggy after I mucked with the 'RemoveTaxonLvAnc' function and made the ancestor NULL.] */
		AddChild(saveAncNode,CVNode); /* important to reattach before next call */
		if (success)
			{
			cvResult[k]=1; /* good */
			++numSuccess;				
			cvScore[k]=cvSquareErrorBranch(tree,CVNode,method,&chiSq); 
			cvSum+=cvScore[k];
			chiSqScore[k]=chiSq;
			chiSqSum+=chiSq;
			printf("+\n");
			}
		else
			{
			cvResult[k]=0; /*failed */
			cvScore[k]=0.0;
			chiSqScore[k]=0.0; /* just put default values in */
			++numFail;
			printf("-\n");
			}
		}
	cvTotalScore[j+1]=cvSum;
	cvTotalScoreChiSq[j+1]=chiSqSum;
	if (numFail==0)
		cvResultFinal[j+1]=1; /* all prunings led to successful optimizations */
	else
		cvResultFinal[j+1]=0; /* some prunings had failed optimizations */
	printf("\n");
	if (verbose>0)
		{	
		printf(".....................................................................\n");
		printf("\nCV Results for smoothing = %f:\nPruned taxon\tSq\t\tChiSq\t\tResult\n",smooth);
		for (i=1;i<=ntips;i++)
			{
			CVNode=(NODE)(pListgetkthNode(tipNodeList,(long)i)->item);
			if (cvResult[i])
				Result=Good;
			else
				Result=Failed;
			printf("%8.8s\t%8.2f\t%8.2f\t%s\n",CVNode->taxon_name,cvScore[i],chiSqScore[i],Result);
/*			printf ("Cross Validation Score [%2i] = %f\t[%f]\n",i,cvScore[i],chiSqScore[i]);
		printf("Cross Validation Score Total (%i pruned terminals):smoothing = %f CV=%f chiSq=%f\n",numSuccess,smooth,cvSum/numSuccess,chiSqSum/numSuccess);
*/
			}
		if (numFail>0)
			printf("** Note that %i failed prunings occurred **\n",numFail);
		}

/* REMEMBER WE ARE OFTEN NOT COUNTING A TIP DESCENDED FROM THE ROOT */

	}
printf("********************************************************************************\n\n");
printf("Results of cross validation analysis for tree %s\n",tree->name);
  switch (method)
	{
	case PENLIKE:printf("Method = Penalized Likelihood\n");break;
	case LaF:printf("Method = Langley and Fitch\n");break;
	case LFLOCAL:printf("Method = Langley and Fitch (with %i local rates)\n",nRates);break;
	case NP:printf("Method = Non-parametric\n");
	}
  switch (algorithm)
	{
	case POWELL: printf("Optimization via Powell's method\n");break;
	case QNEWT:  printf("Optimization via quasi-Newton method with analytical gradients\n");
	}
printf("\nlog10\n");
printf("smooth\tsmooth\t\tSq Error\tChi Square Error\n");
printf("--------------------------------------------------------------------------------\n");
overallGood=1;
bestChiSq=1e20;
bestJ=0;
for (j=0;j<cvNum;j++)
	{
	if (cvResultFinal[j+1])
		{
		if (cvTotalScoreChiSq[j+1]<bestChiSq)
			{
			bestChiSq=cvTotalScoreChiSq[j+1];
			bestJ=j;
			}
		Result=Good;
		}
	else
		{
		overallGood=0;
		Result=Failed;
		}
	smooth=pow(10.0,j*cvInc+cvStart);
	printf("%6.2f\t%6.2g\t\t%6.2f\t%6.2f\t(%s)\n",j*cvInc+cvStart,smooth,cvTotalScore[j+1],cvTotalScoreChiSq[j+1],Result);
	}
printf("********************************************************************************\n\n");

bestSmooth=pow(10.0,bestJ*cvInc+cvStart);
printf("Optimum: %6.2f\t%6.2g\t\t%6.2f\t%6.2f\n",bestJ*cvInc+cvStart,bestSmooth,cvTotalScore[bestJ+1],cvTotalScoreChiSq[bestJ+1]);
if (!overallGood)
	printf("WARNING: Cross validation procedure had errors: optimum may be incorrect\n");

printf("********************************************************************************\n\n");
myFree(cvResult);
myFree(cvResultFinal);
free_vector(cvScore,1,ntips);
free_vector(chiSqScore,1,ntips);
free_vector(cvTotalScore,1,cvNum);
free_vector(cvTotalScoreChiSq,1,cvNum);
freepList(tipNodeList);
return bestSmooth;
}

/***************/




static void 
doFossilCrossV(TREE tree, int method,int nRates,int algorithm,double cvStart,double cvInc,double cvNum, int cvSample)

/*  
	Does a cross validation analysis in which we 
	(1)sequentially unconstrain each node that has a constraint (fixed nodes are not affected)
	(2)do a full estimation on the tree, 
	(3)then calculate the deviation of the estimate for that node versus the constraint (if the constraint is now violated) 
	(4) sums these errors across all constrained nodes

	If the method is LaF or NP then one round of CV is invoked.
	If the method is PENLIKE, then analysis is repeated with the smoothing parameter chosen from a range
	from [cvStart, cvStart+cvInc, ...,cvStart+cvInc*(cvNum-1)].

	Reports two kinds of error, a fractional value per constrained node, and a raw value per constrained node in units of time.

*/
{
char *Result, *Good="Good", *Failed="Failed";
int i,j,k,success,collFlag=0,ntips,verbose,numFixed,numConstrained,numNodes,curIndex;
double * cvScore,*cvScoreRaw, cvSum,chiSq,cvRawSum,*cvTotalScore,*cvTotalScoreRaw,*fixedTime, *estTime;
int * sample, *cvResult, *cvResultFinal;
int numSuccess, numFail,wasConstrainedMin,wasConstrainedMax,wasFixed,wasConstrained;
double smooth,saveTime;
PtrList tipNodeList;
NODETYPE *CVNode;
NODE *nodeArray;
ntips=numdesc(tree->root);
cvResult = (int *)myMalloc((ntips+1)*sizeof(int));
cvResultFinal = (int *)myMalloc((cvNum+1)*sizeof(int));
cvScore = vector(1,ntips);
cvScoreRaw = vector(1,ntips);
cvTotalScore=vector(1,cvNum);
cvTotalScoreRaw=vector(1,cvNum);
tipNodeList = pNewList();
TreeToTaxaPtrList(tree->root,tipNodeList); /* get a list of all the tip nodes */


if (method==LaF || method==NP)
	cvNum=1;	/* just overrides the default value of cvNum for these methods, forcing them to only go once */



numFixed=numFixedNodes(tree->root);
numConstrained=numConstrainedNodes(tree->root);
numNodes=numFixed+numConstrained;
if (numFixed <1 || numConstrained <2)
	{
	doGenericAlert("Must have at least one fixed and two constrained nodes for fossil cross validation");
	return;
	}
fixedTime=(double *)myMalloc(numNodes*sizeof(double));
estTime=(double *)myMalloc(numNodes*sizeof(double));

	// makes a single array beginning with fixed nodes and ending with constrained ones (if any)
nodeArray=(NODE *)myMalloc(numConstrained*sizeof(NODE));
curIndex=0;
//setupFixedNodeArray(tree->root, nodeArray, &curIndex);

setupConstrainedNodeArray(tree->root, nodeArray, &curIndex); 


verbose=gNexDataPtr->RateBlockParms.verbose;
if (verbose > 0)
	printf("Begin fossil cross-validation analyses...\n");
for (j=0;j<cvNum;j++)
	{
	numSuccess=0;
	numFail=0;
    cvSum=0.0;
	cvRawSum=0.0;
	smooth=pow(10.0,j*cvInc+cvStart);
	gNexDataPtr->RateBlockParms.smoothing=smooth;
	for (k=1;k<=numConstrained;k++)
		{
		i=k-1;
		CVNode=nodeArray[i];
		
		if (isConstrainedMax(CVNode))
			{
			wasConstrainedMax=1;
			CVNode->nodeIsConstrainedMax=0;
			}
		else
			wasConstrainedMax=0;
		if (isConstrainedMin(CVNode))
			{
			wasConstrainedMin=1;
			CVNode->nodeIsConstrainedMin=0;
			}
		else
			wasConstrainedMin=0;
			
		gNexDataPtr->RateBlockParms.verbose=0; 
	/* suppress all output from the actual optimization run...may want to allow it for debugging though! */
		doObjFunc(tree,method,nRates,algorithm,&success);
		gNexDataPtr->RateBlockParms.verbose= verbose; /* restore output verbositude to current value */
		if (success)
			{
			estTime[i]=CVNode->time;
			cvResult[k]=1; /* good */
			++numSuccess;
	// Below I deal with all three possible cases: the two simple constraints, or both together;
	// The latter is handled by setting up the next condition and then going through both tests below it
	// Notice that any given time cannot violate both constraints simultaneously

			cvScore[k]=0.0;
			cvScoreRaw[k]=0.0;

			if (wasConstrainedMin)
				{
				if (estTime[i]<CVNode->minAge) // constraint violated, so calculate departure
					{
					cvScore[k]+=2*fabs(CVNode->minAge-estTime[i])/(CVNode->minAge+estTime[i]); 
					cvScoreRaw[k]+=fabs(CVNode->minAge-estTime[i]);
					}
				}
			if (wasConstrainedMax)
				{
				if (estTime[i]>CVNode->maxAge) // constraint violated, so calculate departure
					{
					cvScore[k]+=2*fabs(CVNode->maxAge-estTime[i])/(CVNode->maxAge+estTime[i]); 
					cvScoreRaw[k]+=fabs(CVNode->maxAge-estTime[i]);
					}
				}
			cvRawSum+=cvScoreRaw[k];
			cvSum+=cvScore[k];
			}
		else
			{
			cvResult[k]=0; /*failed */
			cvScore[k]=-99.9;
			cvScoreRaw[k]=0.0; /* should probably set this to some bogus value! */
			++numFail;
			}

	// Restore constraints for this node
		if (wasConstrainedMax)
			CVNode->nodeIsConstrainedMax=1; // Notice, unconstraining does not delete constraint times in struct
		if (wasConstrainedMin)
			CVNode->nodeIsConstrainedMin=1;

		}
		
	printf("\nFossil-constrained cross validation analysis\n\n");
	printf("\tNode\t\tEst Age\t\tMin\tMax\t\tFract Score\tRaw Score\n");
	printf("----------------------------------------------------------------------------------------\n");
	for (i=0;i<numConstrained;i++)
		{
		CVNode=nodeArray[i];
		printf("%i\t%.6s\t\t%6.2f\t\t",i+1,CVNode->taxon_name,estTime[i]);
		if (isConstrainedMin(CVNode))
			printf("%6.2f\t",CVNode->minAge);
		else
			printf("   --   ");
		if (isConstrainedMax(CVNode))
			printf("%6.2f\t",CVNode->maxAge);
		else
			printf("   --   ");
		printf("\t%f\t%6.2f\n",cvScore[i+1],cvScoreRaw[i+1]);
		}

	cvTotalScore[j+1]=cvSum/(numNodes-numFail); /* on the off chance that some reps failed, don't count their contributions */
	cvTotalScoreRaw[j+1]=cvRawSum/(numNodes-numFail);
	
	if (numFail==0)
		cvResultFinal[j+1]=1; /* all prunings led to successful optimizations */
	else
		cvResultFinal[j+1]=0; /* some prunings had failed optimizations */

/* REMEMBER WE ARE OFTEN NOT COUNTING A TIP DESCENDED FROM THE ROOT */

	}
printf("********************************************************************************\n\n");
printf("Summary results of fossil-constrained cross validation analysis for tree %s\n",tree->name);
  switch (method)
	{
	case PENLIKE:printf("Method = Penalized Likelihood\n");break;
	case LaF:printf("Method = Langley and Fitch\n");break;
	case LFLOCAL:printf("Method = Langley and Fitch (with %i local rates)\n",nRates);break;
	case NP:printf("Method = Non-parametric\n");
	}
  switch (algorithm)
	{
	case POWELL: printf("Optimization via Powell's method\n");break;
	case QNEWT:  printf("Optimization via quasi-Newton method with analytical gradients\n");
	}
printf("\nFixed nodes:%i\nConstrained nodes:%i\n",numFixed,numConstrained);
printf("\nlog10\n");
printf("smooth\tsmooth\t\tFract Error\tRaw Error\n");
printf("--------------------------------------------------------------------------------\n");
for (j=0;j<cvNum;j++)
	{
	if (cvResultFinal[j+1])
		Result=Good;
	else
		Result=Failed;
	smooth=pow(10.0,j*cvInc+cvStart);
	printf("%6.2f\t%6.2g\t\t%6.4f\t\t%6.4f\t(%s)\n",j*cvInc+cvStart,smooth,cvTotalScore[j+1],cvTotalScoreRaw[j+1],Result);
	}
printf("********************************************************************************\n\n");
myFree(cvResult);
myFree(cvResultFinal);
myFree(fixedTime);
myFree(estTime);
myFree(nodeArray);
free_vector(cvScore,1,ntips);
free_vector(cvScoreRaw,1,ntips);
free_vector(cvTotalScore,1,cvNum);
free_vector(cvTotalScoreRaw,1,cvNum);
freepList(tipNodeList);


return;
}


/****************************************************************/
static void 
doFossilCrossVfixed(TREE tree, int method,int nRates,int algorithm,double cvStart,double cvInc,double cvNum, int cvSample)

/*  
	Does a cross validation analysis in which we 
	(1)sequentially fix each node that is fixed 
	(2)do a full estimation on the tree, 
	(3)then calculate the deviation of the estimate for that node versus the original fixed value 
	(4) sums these errors across all fixed nodes 

	If the method is LaF or NP then one round of CV is invoked.
	If the method is PENLIKE, then analysis is repeated with the smoothing parameter chosen from a range
	from [cvStart, cvStart+cvInc, ...,cvStart+cvInc*(cvNum-1)].

	Reports two kinds of error, a fractional value per constrained node, and a raw value per constrained node in units of time.

*/
{
char *Result, *Good="Good", *Failed="Failed";
int i,j,k,success,collFlag=0,ntips,verbose,numFixed,numConstrained,numNodes,curIndex;
double * cvScore,*cvScoreRaw, cvSum,chiSq,cvRawSum,*cvTotalScore,*cvTotalScoreRaw,*fixedTime, *estTime;
int * sample, *cvResult, *cvResultFinal;
int numSuccess, numFail,wasConstrainedMin,wasConstrainedMax,wasFixed,wasConstrained;
double smooth,saveTime;
PtrList tipNodeList;
NODETYPE *CVNode;
NODE *nodeArray;
ntips=numdesc(tree->root);
cvResult = (int *)myMalloc((ntips+1)*sizeof(int));
cvResultFinal = (int *)myMalloc((cvNum+1)*sizeof(int));
cvScore = vector(1,ntips);
cvScoreRaw = vector(1,ntips);
cvTotalScore=vector(1,cvNum);
cvTotalScoreRaw=vector(1,cvNum);
tipNodeList = pNewList();
TreeToTaxaPtrList(tree->root,tipNodeList); /* get a list of all the tip nodes */


if (method==LaF || method==NP)
	cvNum=1;	/* just overrides the default value of cvNum for these methods, forcing them to only go once */



numFixed=numFixedNodes(tree->root);
numNodes=numFixed;
if (numFixed <2 )
	{
	doGenericAlert("Must have at least two fixed nodes for fossil cross validation");
	return;
	}
fixedTime=(double *)myMalloc(numNodes*sizeof(double));
estTime=(double *)myMalloc(numNodes*sizeof(double));

nodeArray=(NODE *)myMalloc(numNodes*sizeof(NODE));
curIndex=0;
setupFixedNodeArray(tree->root, nodeArray, &curIndex);

verbose=gNexDataPtr->RateBlockParms.verbose;
if (verbose > 0)
	printf("Begin fossil cross-validation analyses...\n");
for (j=0;j<cvNum;j++)
	{
	numSuccess=0;
	numFail=0;
    cvSum=0.0;
	cvRawSum=0.0;
	smooth=pow(10.0,j*cvInc+cvStart);
	gNexDataPtr->RateBlockParms.smoothing=smooth;
	for (k=1;k<=numNodes;k++)
		{
		i=k-1;
		CVNode=nodeArray[i];
		
		CVNode->free=1;
		saveTime=CVNode->time;	
		gNexDataPtr->RateBlockParms.verbose=0; 
	/* suppress all output from the actual optimization run...may want to allow it for debugging though! */
		doObjFunc(tree,method,nRates,algorithm,&success);
		gNexDataPtr->RateBlockParms.verbose= verbose; /* restore output verbositude to current value */
		if (success)
			{
			estTime[i]=CVNode->time;
			cvResult[k]=1; /* good */
			++numSuccess;

			cvScore[k]=2*fabs(saveTime-estTime[i])/(saveTime+estTime[i]); 
			cvScoreRaw[k]=fabs(saveTime-estTime[i]);
			cvRawSum+=cvScoreRaw[k];
			cvSum+=cvScore[k];
			}
		else
			{
			cvResult[k]=0; /*failed */
			cvScore[k]=-99.9;
			cvScoreRaw[k]=0.0; /* should probably set this to some bogus value! */
			++numFail;
			}
		CVNode->free=0; // restore fixity and age
		CVNode->time=saveTime;
		}
		
	printf("\nFossil-fixed cross validation analysis\n\n");
	printf("Node\tTaxon\t\tAge\t\tEst.Age\t\tFract.Score\tRaw.Score\n");
	printf("---------------------------------------------------------------------------\n");
	for (i=0;i<numNodes;i++)
		{
		CVNode=nodeArray[i];
		printf("%i\t%.6s\t\t%6.2f\t\t%6.2f\t",i+1,CVNode->taxon_name,CVNode->time,estTime[i]);
		printf("\t%f\t%6.2f\n",cvScore[i+1],cvScoreRaw[i+1]);
		}

	cvTotalScore[j+1]=cvSum/(numNodes-numFail); /* on the off chance that some reps failed, don't count their contributions */
	cvTotalScoreRaw[j+1]=cvRawSum/(numNodes-numFail);
	
	if (numFail==0)
		cvResultFinal[j+1]=1; /* all prunings led to successful optimizations */
	else
		cvResultFinal[j+1]=0; /* some prunings had failed optimizations */

/* REMEMBER WE ARE OFTEN NOT COUNTING A TIP DESCENDED FROM THE ROOT */

	}
printf("********************************************************************************\n\n");
printf("Summary results of fossil-fixed cross validation analysis for tree %s\n",tree->name);
  switch (method)
	{
	case PENLIKE:printf("Method = Penalized Likelihood\n");break;
	case LaF:printf("Method = Langley and Fitch\n");break;
	case LFLOCAL:printf("Method = Langley and Fitch (with %i local rates)\n",nRates);break;
	case NP:printf("Method = Non-parametric\n");
	}
  switch (algorithm)
	{
	case POWELL: printf("Optimization via Powell's method\n");break;
	case QNEWT:  printf("Optimization via quasi-Newton method with analytical gradients\n");
	}
printf("\nFixed nodes:%i\n",numFixed);
printf("\nlog10\n");
printf("smooth\tsmooth\t\tFract Error\tRaw Error\n");
printf("--------------------------------------------------------------------------------\n");
for (j=0;j<cvNum;j++)
	{
	if (cvResultFinal[j+1])
		Result=Good;
	else
		Result=Failed;
	smooth=pow(10.0,j*cvInc+cvStart);
	printf("%6.2f\t%6.2g\t\t%6.4f\t\t%6.4f\t(%s)\n",j*cvInc+cvStart,smooth,cvTotalScore[j+1],cvTotalScoreRaw[j+1],Result);
	}
printf("********************************************************************************\n\n");
myFree(cvResult);
myFree(cvResultFinal);
myFree(fixedTime);
myFree(estTime);
myFree(nodeArray);
free_vector(cvScore,1,ntips);
free_vector(cvScoreRaw,1,ntips);
free_vector(cvTotalScore,1,cvNum);
free_vector(cvTotalScoreRaw,1,cvNum);
freepList(tipNodeList);


return;
}


/****************************************************************/

/* Following routines handle the custom 'Rate' block */

void doRateBlock(void)
{
  NODETYPE *root;
  char* TD, *TreeName;
  int ix;
  long numTrees,j;
extern int gisConstrained;
extern double gGamma_c;
extern double gGamma_b;
PtrList lnode;
TREE thisTree;
 

gTaxaList=NULL;	   /* initialize this somewhere else ? */
  

/* 
 * Following code sets up an exclusion array the FIRST time this block is called.
 * After that it is assumed to be there.  Therefore,  if you want to re-execute
 * a new matrix with different num of chars,  we'll have to set the pointer to 
 * the exclusion array to NULL!
 */

if ((gNexDataPtr->NChars > 0 ) && (gNexDataPtr->excArray == NULL))
 {
 gNexDataPtr->excArray=(int*)myMalloc(gNexDataPtr->NChars*sizeof(int)); 
 if(gNexDataPtr->excArray ==NULL)
    fatal("Allocation error in excArray");
 }
if (gNexDataPtr->excArray)
  for (ix=0;ix<gNexDataPtr->NChars;ix++)
    gNexDataPtr->excArray[ix]=1;



  do 				/* need to put in error checking in case no DIMENSIONS statement */
	    {
	     aTokenPtr=nextToken();
	     if (*aTokenPtr=='\0')
		{
		return;
		}
	     
	     else if   (isEqual(aTokenPtr, "QUIT") || 
			isEqual(aTokenPtr, "Q") ||
			isEqual(aTokenPtr, "BYE"))
			exit(1);
	     else if(isEqual(aTokenPtr, "COVARION"))
			doCovarionCommand();
	     else if(isEqual(aTokenPtr, "TRAITEVOL"))
			doCovarionCommand();
	     else if(isEqual(aTokenPtr, "ANC"))
			doAncestralCommand();
	     else if(isEqual(aTokenPtr, "CO"))
			doContOptCommand();
	     else if(isEqual(aTokenPtr, "VCV"))
			doVCVCommand();
	     else if(isEqual(aTokenPtr, "RRLIKE"))
			doRRLikeTestCommand();
	     else if(isEqual(aTokenPtr, "LOCALMODEL"))
			doLocalModelCommand();
	     else if(isEqual(aTokenPtr, "BLFORMAT"))
			doBLFormatCommand();
	     else if(isEqual(aTokenPtr, "SHOWAGE") || isEqual(aTokenPtr, "SHOW") )
			doShowAgeCommand();
	     else if(isEqual(aTokenPtr, "UNFIXAGE"))
			doUnSetAgeCommand();
	     else if(isEqual(aTokenPtr, "FIXAGE"))
			doSetAgeCommand();
	     else if(isEqual(aTokenPtr, "PRUNE"))
			doPruneTaxonCommand();
	     else if(isEqual(aTokenPtr, "REROOT"))
			doReRootCommand();
	     else if(isEqual(aTokenPtr, "CLEARTREES"))
			doClearTreesCommand();
	     else if(isEqual(aTokenPtr, "CLUSTER_HISTOGRAM"))
			doClusterHistogramCommand();
	     else if(isEqual(aTokenPtr, "COLLAPSE"))
			doCollapseCommand();
	     else if(isEqual(aTokenPtr, "EXECUTE"))
			doExecuteCommand();
	     else if (isEqual(aTokenPtr,"PROFILE"))  
			doBranchProfileCommand();
	     else if (isEqual(aTokenPtr,"MRCA"))  
			doMRCACommand();
	     else if (isEqual(aTokenPtr,"CLADE_CHECK"))  
			doCladeCheckCommand();
	     else if (isEqual(aTokenPtr,"CONSTRAIN"))  
		doConstrain_TimeCommand();
	     else if (isEqual(aTokenPtr,"CALIBRATE"))  
		doCalibrateCommand();
	     else if (isEqual(aTokenPtr,"DESCRIBE") || isEqual(aTokenPtr,"DESC"))  
		doPrintCommand();
	     else if (isEqual(aTokenPtr,"TAXASET"))  
		doTaxaSetCommand();
	     else if (isEqual(aTokenPtr,"DISTANCE"))  
		doDistanceCommand();
	     else if (isEqual(aTokenPtr,"BS"))  
		doBSCommand();
	     else if (isEqual(aTokenPtr,"SET"))  
		doSetCommand();
	     else if (isEqual(aTokenPtr,"MRP"))
		doSuperCommand();  
	     else if (isEqual(aTokenPtr,"DIVTIME") || isEqual(aTokenPtr,"DIV"))
		doDivTimeCommand();  
	     else if (isEqual(aTokenPtr,"BD"))
		doBD();
	     else if (isEqual(aTokenPtr,"SIMULATE") )  
			doSimCommand();
	     else if (isEqual(aTokenPtr,"MULTIPLY_BRANCHLENGTH_BY"))
		doLengthMultiplyCommand();
	     else if (isEqual(aTokenPtr,"MULT"))
		doLengthMultiplyCommand();
	     else if (isEqual(aTokenPtr,"CONVERT_BRANCHLENGTH_TO_TIME"))
		{
	        if (gNexDataPtr->isTrees)
			{
			lnode=gNexDataPtr->inTrees;
			LISTLOOP (lnode)
				{
				thisTree=lnode->item;
				convert_branchlength_to_time(thisTree->root);
				}
			}
		}
	     else if (isEqual(aTokenPtr,"RR"))
		{  
	        if (gNexDataPtr->isTrees)
			{
			lnode=gNexDataPtr->inTrees;
			LISTLOOP (lnode)
				{
				thisTree=lnode->item;
				doWuLiCommand(thisTree->root);
				}
			}
		}
    }  while (!isEqual(aTokenPtr,"END") && !isEqual(aTokenPtr,"ENDBLOCK") && aTokenPtr);
    aTokenPtr=nextToken();
    if (!isEqual(aTokenPtr,";"))
	    doGenericAlert("Block not terminated with semicolon");
return;
}


/****************************************************************/
/**************** COMMAND PROCESSING FUNCTIONS ******************/
/****************************************************************/
int gNtips, gBDnumvar,gStemFlag;
double *gTimes,gRootDur;
#define POLYTOMY 1
#define SQR(x) (x)*(x)

void doBD(void) /* do birth-death stats on a specified clade */
    {
    extern NODETYPE* gRoot;
    long minN=1000000000,maxN=-1000000000;
    PtrList lnode;
    TREE thisTree;
    int N0,success, numIter, i, n, nint, ntrees, ix=1, profile_flag=0, dp_flag=0, Nee_flag=0, sumerr=0, ixm=0, numvar;
    double K1,K1var,K2,K2var,pure_birth_estimate, B, Optimum, a, r, AgeRoot, AgeFoundNode, CalAgeRoot, CalAgeFoundNode=-1.0, 
	    Kendall_var, Moran_var, mean1,adev1,sdev1,var1,skew1,curt1, mean2,adev2,sdev2,var2,skew2,curt2, 
	    Kendall_estimate, Raup_estimate, Like1, Like2, Like2Param,Like1Param,LR,SG,L1P,L2first,L2second,rate_ml;
    double params[3], *Y, *X,  *data1, *data2, *data3, *data4, *data5, *data6, *data7, *data8,*data9;
	double r1,r2,r0,LR_other,gamma;
    int root_id=1; /* as default use the root node */
    int stemFlag=0, /* default assume crown group */
	gammaFlag=0;
	int SisterLRflag=0;
	int countLR=0,countSG=0,n1,n2,smaller;
    char *dummy,  *taxon;
    double rootAge;
    NODETYPE * found_node, *root, *first, *second;
    taxon=NULL;
    while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			{ /* must accept either a name(string) or taxon id (int) */
			if (isdigit(*LocalToken))
			    root_id=strtod(LocalToken,&dummy);
			else
			    taxon=DupStr(LocalToken);
			}
		if (parse_assignment2("GAMMA"))
			if (isEqual(LocalToken,"YES"))
				gammaFlag=1;
			else
				gammaFlag=0;
		if (parse_assignment2("AGE"))
			CalAgeFoundNode=strtod(LocalToken,&dummy);
		if (parse_assignment2("PROFILE"))
			if (isEqual(LocalToken,"YES"))
				profile_flag=1;
			else
				profile_flag=0;
		if (parse_assignment2("DIVPLOT"))
			if (isEqual(LocalToken,"YES"))
				dp_flag=1;
			else
				dp_flag=0;
		if (parse_assignment2("NEE"))
			if (isEqual(LocalToken,"YES"))
				Nee_flag=1;
			else
				Nee_flag=0;
		if (parse_assignment2("SISTER_LR"))
			if (isEqual(LocalToken,"YES"))
				SisterLRflag=1;
			else
				SisterLRflag=0;
		if (parse_assignment2("STEM"))
			if (isEqual(LocalToken,"YES"))
				stemFlag=1;
			else
				stemFlag=0;
		
		}

if (gNexDataPtr->isTrees)
		{
		lnode=gNexDataPtr->inTrees;
		ntrees=pLengthList(lnode);
		/*if(profile_flag)*/ /* always need some of these now */
		    {
		    data1=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data2=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data3=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data4=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data5=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data6=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data7=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data8=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    data9=(double*)myMalloc((ntrees+1)*sizeof(double)); /* 1-offset array */
		    }
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;



    if (taxon)
	    found_node=find_taxon_name(root,taxon);
    else
	{ 
	if (root_id>1)
	    found_node=find_id(root, root_id);
	else if (root_id==1)
	    found_node=root;	/* the default setting; look at whole tree */
	}
    if(found_node)
	{
	root_id=found_node->id;
	gTimes=sort_node_time_array(found_node); /* zero-offset array gTimes */
	
	
	n=numdesc(found_node)-1; /* See comments under sort_node_time_array() */
	if (n>=2)
    		if(gTimes[n-1]==gTimes[n-2])
			++sumerr;

	if (gammaFlag)
		{
	  	gamma = PH_gamma(n+1,gTimes,root->time);
	    	data9[ix]=gamma;
		printf("Gamma statistic = %f\n",gamma);
		}

	if(dp_flag)  /* Plot diversity over time.  WORKS ONLY FOR RECONSTRUCTED PROCESS (ASSUMES NO EXTINCTION)*/
	    {
	    Y=(double*)myMalloc(n*sizeof(double));
	    X=(double*)myMalloc(n*sizeof(double));
	    for (i=0;i<n;i++)
		{
		Y[i]=log(i+2); 
		X[n-1-i]=found_node->time-gTimes[i];
		}
	    dumbPlot(X,  Y, n);
	    myFree(X);
	    myFree(Y);
	    }

// ** Seems like the calibration stuff is a bit mucked up, check Kendal_est..

	AgeRoot=root->time;
	AgeFoundNode=found_node->time;
	if (CalAgeFoundNode == -1.0)
	    CalAgeFoundNode = found_node->time; /* if we didn't read a calibrated age for this, just set to internal value */
	CalAgeRoot=AgeRoot*CalAgeFoundNode/AgeFoundNode;
	B=get_sum_durations(found_node);
	gNtips=numdesc(found_node);
	if (stemFlag)
		N0=1;
	else
		N0=2;
	Kendall_estimate=((gNtips-N0 )/B)        /*/  CalAgeRoot   */;
	Kendall_var=SQR(Kendall_estimate)/(N0*(exp(Kendall_estimate *CalAgeFoundNode -1)));
	Moran_var=SQR(Kendall_estimate)/(gNtips-N0);
	Raup_estimate=(log(gNtips)-log(N0))/CalAgeFoundNode;
	if(profile_flag)
	    {
	    data1[ix]=Kendall_estimate;
	    data2[ix]=Raup_estimate;
	    data3[ix]=gNtips;
		if (gNtips<minN)minN=gNtips;
		if (gNtips>maxN)maxN=gNtips;
	    data4[ix]=B;
	    }

	if (!profile_flag)
	    {
	    printf("Age of root = %f\n", AgeRoot);
	    printf("ML estimate of lineage birth rate under Yule model using durations (node %i numtips=%i B=%f ):\n", 
			    root_id, gNtips, B);
	    printf("Whole Tree Root internal age = %f....calibrated age = %f\n", AgeRoot, CalAgeRoot);
	    printf("Subtree for BD root internal age = %f....calibrated age = %f\n", AgeFoundNode, CalAgeFoundNode);
	    printf("Kendall estimate of lineage birth rate under Yule model = %f (%f)\n",Kendall_estimate );
	    printf("Kendall's 1949 estimate of variance and std dev of rate estimate: %f\t%f\n",
		    Kendall_var, sqrt(Kendall_var) );
	    printf("Moran's 1951 estimate of variance and std dev of rate estimate: %f\t%f\n",
		    Moran_var, sqrt(Moran_var) );
	    printf("'Raup' estimate of lineage birth rate under Yule model (log(N)/t))= %f\n", Raup_estimate);
	    }
	if (Nee_flag) 
	    {
	    gBDnumvar=2; /* 1 = pure birth model (seems to agree with Kendall estimator in practice!), 2= birth-death model */
	    params[1]=2.0;
	    params[2]=0.0;
	    numIter=100;
	    Like2 = -MinND(thisTree,-1,POWELL,BD_Like,NULL,params,gBDnumvar,&numIter, 0.00001,0.0001, &success );
		// the -1 is a dummy argument, MinND doesn't need it
	    r=params[1]/CalAgeRoot;
	    a=params[2];
	    data6[ix]=r;
	    data7[ix]=a;
	    printf("ML estimate in the BD model:a=%f \tr=%f\tspec=%f\textinct=%f\n", a, r, r/(1-a), a*r/(1-a));
	    gBDnumvar=1;
	    Like1 = -MinND(thisTree,-1,POWELL,BD_Like,NULL,params,gBDnumvar,&numIter, 0.00001,0.0001, &success );
	    r=params[1]/CalAgeRoot;
	    data8[ix]=r;
	    printf("ML estimate in the B (Yule) model:r=%f\n", r);
	    LR=2*(Like2-Like1);
	    data5[ix]=LR;
	    printf("Like1=%f Like2=%f LR=%f\n", Like1, Like2, LR);

	/*    plotOpt(params,10,0.0,1.0,0.0,5.0, "a", "r");*/
	    }
	if (SisterLRflag) 
	    {
		
/* A technical issue on the LR test for sister group diversities. The likelihood for a stem or crown clade is 

		L(n) = (n-1)! * rate ^ (n-n0) * exp(-rate*sum_dur).

The LR test is  

		L(n1)*L(n2) /L(n1+n2).

It seems like the factorial coefficients might be something like (n1-1)!(n2-1)!/(n1+n2-1)!, which would be very different from 1.
However, in the denominator, we are ALSO looking only at diversification processes that have produced two sister clades with 
exactly n1 and n2 species, so the correct coefficient for the denominator is also (n1-1)!(n2-1)!. Actually, it seems like there 
should be an additional factor of 2 in both numerator and denominator, for the two events of n1,n2 or n2,n1 (if n1 ne n2).
In any case, the coefficients drop out, and we can actually calculate the YuleLike without them, at least for this sole purpose of
doing an LR test. Note that the function does NOT calculate the coefficients.
*/
	    //gBDnumvar=1; /* 1 = pure birth model (seems to agree with Kendall estimator in practice!), 2= birth-death model */
	    //params[1]=1.0;
	    //params[2]=0.0;
	    //numIter=100;
		
		root=thisTree->root;
		first=root->firstdesc;
		second=first->sib;
	// two param model
		n1=gNtips=numdesc(first);
		//gStemFlag = 1;
		//gRoot=first;
	    //numIter=100;
	    //Like1 = -MinND(thisTree,-1,POWELL,BD_Like,NULL,params,gBDnumvar,&numIter, 0.00001,0.0001, &success );
	    //r=params[1]/CalAgeRoot;
		//data6[ix]=r;
	    //a=params[2];
	//    printf("ML estimate in the first sister:a=%f \tr=%f\tspec=%f\textinct=%f\n", a, r, r/(1-a), a*r/(1-a));
			// the -1 is a dummy argument, MinND doesn't need it
		n2=gNtips=numdesc(second);
		//gStemFlag = 1;
		//gRoot=second;
	    //numIter=100;
	    //Like2 = -MinND(thisTree,-1,POWELL,BD_Like,NULL,params,gBDnumvar,&numIter, 0.00001,0.0001, &success );
	    //r=params[1]/CalAgeRoot;
		//data7[ix]=r;
	    //a=params[2];
	    //data8[ix]=r;
	 //   printf("ML estimate in the second sister:a=%f \tr=%f\tspec=%f\textinct=%f\n", a, r, r/(1-a), a*r/(1-a));
		//Like2Param = Like1+Like2;
	// one param model
		//gNtips=numdesc(root);
		//gStemFlag = 0;
		//gRoot=root;
	    //numIter=100;
		//Like1Param = -MinND(thisTree,-1,POWELL,BD_Like,NULL,params,gBDnumvar,&numIter, 0.00001,0.0001, &success );
	    //r=params[1]/CalAgeRoot;
		//data8[ix]=r;
	    //a=params[2];
		L2first=YuleLike(first, 1, &rate_ml);
		data6[ix]=rate_ml;
		r1=rate_ml;
		printf("L2first=%f mlrate2param=%f\n",L2first,rate_ml);
		L2second=YuleLike(second, 1, &rate_ml);
		data7[ix]=rate_ml;
		r2=rate_ml;
		printf("L2second=%f mlrate2param=%f\n",L2second,rate_ml);
		L1P=YuleLike(root, 0, &rate_ml);
		data8[ix]=rate_ml;
		r0=rate_ml;
	//    printf("ML estimate in the whole tree:a=%f \tr=%f\tspec=%f\textinct=%f\n", a, r, r/(1-a), a*r/(1-a));

	   // LR=2*(Like2Param-Like1Param);
		LR = 2*(L2first+L2second-L1P);
		//LR_other =  2*(log(r1/r0)*(n1-1)+log(r2/r0)*(n2-1)) ; alternative sweet and correct formula for the LR test (but undefined if n = 2 for crown)
		//printf("*************** %f %f\n",LR,LR_other);
// LR /= 1.29; chi-square correction factor...
	    data5[ix]=LR;
		if (LR>3.841) ++countLR;
		if (n1<n2) smaller=n1; else smaller=n2;
		SG = 2.0*smaller/(n1+n2-1);
		if (n1==n2) SG=1; //boundary case
		if (SG<0.05) ++countSG;
		//printf("n1=%i n2=%i LR=%f SG=%f\n",n1,n2,LR,SG);
	    //printf("Like1=%f Like2=%f Like1Param=%f Like2Param=%f LR=%f\n", Like1,Like2, Like1Param, Like2Param, LR);

		printf("L1P=%f mlrate1param=%f\n",L1P,rate_ml);



	/*    plotOpt(params,10,0.0,1.0,0.0,5.0, "a", "r");*/
	    }
	    
	} /* end found_node */
	++ix;
			} /* end LISTLOOP */
		}
	 if (stemFlag)
		printf("(Stem group simulation:N0=1)\n");
	 else
		printf("(Crown group simulation:N0=2)\n");
    if(profile_flag)
	{
	printf("\n****************\nProfile Analysis of %i trees\n", ntrees);
	moment(data3,ntrees,&mean1,&adev1,&sdev1,&var1,&skew1,&curt1);
	printf("Mean clade size = %f (range=[%li,%li])\n", mean1,minN,maxN);
	moment(data4,ntrees,&mean1,&adev1,&sdev1,&var1,&skew1,&curt1);
	printf("Mean B = %f\n", mean1);
	moment(data1,ntrees,&mean1,&adev1,&sdev1,&var1,&skew1,&curt1); K1=mean1;K1var=var1;
	printf("Summary statistics on Yule(Kendall perfect information) estimator for %i trees\n");
	printf("Mean diversification rate = %f\n", mean1);
	printf("Variance and standard deviation of diversification rate = %f \t%f\n", var1, sdev1);
	moment(data2,ntrees,&mean2,&adev2,&sdev2,&var2,&skew2,&curt2);K2=mean2;K2var=var2;
	printf("Summary statistics on Yule(Raup minimal information) estimator for %i trees\n");
	printf("Mean diversification rate = %f\n", mean2);
	printf("Variance and standard deviation of diversification rate = %f \t%f\n", var2, sdev2);
	moment(data5,ntrees,&mean2,&adev2,&sdev2,&var2,&skew2,&curt2);
	printf("Summary statistics on LR Test %i trees\n");
	printf("Mean and variance = %f \t%f\n",mean2, var2);
	moment(data6,ntrees,&mean2,&adev2,&sdev2,&var2,&skew2,&curt2);
	printf("Mean and sdev of r in Yule model for clade 1 = %f \t%f\n",mean2, sdev2);
	moment(data7,ntrees,&mean2,&adev2,&sdev2,&var2,&skew2,&curt2);
	printf("Mean and sdev of r in Yule model for clade 2 %f \t%f\n",mean2, var2);
	moment(data8,ntrees,&mean2,&adev2,&sdev2,&var2,&skew2,&curt2);
	printf("Mean and sdev of r  in Yule model for whole clade %f \t%f\n",mean2, var2);
	moment(data9,ntrees,&mean2,&adev2,&sdev2,&var2,&skew2,&curt2);
	printf("Mean and sdev of gamma  %f \t%f\n",mean2, var2);

	if (SisterLRflag) 
		printf ("reject(LR|SG)?\t%i\t%i\n",countLR,countSG);
//	printf("stats\t%f\t%f\t%f\t%f\n",K1,K1var,K2,K2var);

	myFree(data1);
	myFree(data2);
	}
    if (taxon)
	    myFree(taxon);
    return;
    }
/****************************************************************/
void doPrintCommand(void)  /* Draws ASCII version of trees */
{
int 	likeFlag=0;
static int plotwidth=0;
int 	whichTrees=0;	/*0=input; 1=output*/
static int	treemode=0;	/*0=cladogram, etc.*/
int	numTrees,j, treeix=1;
static  int tree=0;
char *	dummy;
NODETYPE *root;
PtrList lnode;
TREE thisTree;
char*	TD, *TreeName;
char 	*clado="CLADOGRAM", *phylo="PHYLOGRAM",*chrono="CHRONOGRAM", *rato="RATOGRAM", 
	*td="TREE DESCRIPTION", *pd="PHYLO DESCRIPTION",  *rd=
	"RATO DESCRIPTION", *ni="NODE INFORMATION", *id="ID INFORMATION", *trace="TRACE", 
	*tracephy="TRACEPHY", *marg="MARGINALS DESCRIPTION", *modeString;

	while (!isEqual(aTokenPtr=nextToken(),";")) 
		{
		if (parse_assignment2("TREE"))
			tree=strtod(LocalToken,&dummy);
		if (parse_assignment2("PLOTWIDTH"))
			plotwidth=strtod(LocalToken,&dummy);
		if (parse_assignment2("DISPLAY_ANCESTOR"))
			if (isEqual(LocalToken,"YES"))
				gLabel=1;
			else
				gLabel=0;
		if (parse_assignment2("SOURCE"))
			{
			if (isEqual(LocalToken,"INPUT"))
				whichTrees=0;
			else
				whichTrees=1;
			}
		if (parse_assignment2("PLOT"))
			{
			if (isEqual(LocalToken,"CLADOGRAM"))
				treemode=0;
			if (isEqual(LocalToken,"PHYLOGRAM"))
				treemode=1;
			if (isEqual(LocalToken,"CHRONOGRAM"))
				treemode=2;
			if (isEqual(LocalToken,"NODE_INFO"))
				treemode=3;
			if (isEqual(LocalToken,"RATOGRAM"))
				treemode=4;
			if (isEqual(LocalToken,"CHRONO_DESCRIPTION"))
				treemode=5; /* also prints durations */
			if (isEqual(LocalToken,"TREE_DESCRIPTION"))
				treemode=5; /* also prints durations */
			if (isEqual(LocalToken,"PHYLO_DESCRIPTION"))
				treemode=6;
			if (isEqual(LocalToken,"RATO_DESCRIPTION"))
				treemode=7;
			if (isEqual(LocalToken,"ID_DESCRIPTION"))
				treemode=8;
			if (isEqual(LocalToken,"TRACE"))
				treemode=9;
			if (isEqual(LocalToken,"TRACEPHY"))
				treemode=10;
			if (isEqual(LocalToken,"MARG_DESCRIPTION"))
				treemode=11;
			}
		if (parse_assignment2("LIKE"))
			{
			if (isEqual(LocalToken,"YES"))
				likeFlag=1;
			else
				likeFlag=0;
			}
		}
/* now do it */

	switch (treemode) 
		{
		case 0: modeString=clado;break;
		case 1: modeString=phylo;break;
		case 2: modeString=chrono;break;
		case 3: modeString=ni;break;
		case 4: modeString=rato;break;
		case 5: modeString=td;break;
		case 6: modeString=pd;break;
		case 7: modeString=rd;break;
		case 8: modeString=id;break;
		case 9: modeString=trace;break;
		case 10: modeString=tracephy;break;
		case 11: modeString=marg;break;
		}
	if (whichTrees==0)
		{
		if (gNexDataPtr->isTrees)
			{
			lnode=gNexDataPtr->inTrees;
			LISTLOOP (lnode)
				{
				if ((tree==0) || (tree==treeix)) 
				    /* if no tree specified OR a specific tree specified */
				  {
				  thisTree=lnode->item;
				  printf("[Printing tree %i]\n", treeix);
				  printf("\n[%s of tree %s]\n",modeString, thisTree->name);
				  if (treemode==8)
				    {
				    printf("\ntree %s = ", thisTree->name);
				    make_parens(thisTree->root, 3); /* TD with ids for node names*/
				    printf(";\n");
				    }
				  if (treemode==5)
				    {
				    printf("\ntree %s = ", thisTree->name);
				    make_parens(thisTree->root, 1); /* TD of chronogram*/
				    printf(";\n");
				    }
				  if (treemode==6)
				    {
				    printf("\ntree %s = ", thisTree->name);
				    make_parens(thisTree->root, 0); /* TD of phylogram*/
				    printf(";\n");
				    }
				  if (treemode==11)
				    {
				    printf("#nexus\nbegin trees;\n\ntree %s = ", thisTree->name);
				    make_parens(thisTree->root, 4); /* marginals description string*/
				    printf(";\nend;\n");
				    }
				  if (treemode==7)
				    {
				 /*   set_est_rates(thisTree->root,thisTree->est_b,thisTree->est_c);*/
				    printf("\ntree %s = ", thisTree->name);
				    make_parens(thisTree->root, 2); /* TD of phylogram*/
				    printf(";\n");
				    }
				  if (treemode==3)
				    {
				    printtree(thisTree->root);
				    }
				  if (  (treemode==0) || 
					(treemode==1) || 
					(treemode==9) || 
					(treemode==10) || 
					(treemode==2) || 
					(treemode==4))
				    DrawTree(thisTree->root,treemode, plotwidth);
				  if (likeFlag)
					printLikes(thisTree->root);
				  }
				++treeix;
				}
			}
		else
			doGenericAlert("No input trees available\n");
		}
	/* else.....*/
	return;
}
/****************************************************************/
static void doBLFormatCommand(void)
{
PtrList lnode;
TREE thisTree;
char * dummy;
extern long gNumSites;
struct RBP * rbp;
int roundflag=1;

rbp=&(gNexDataPtr->RateBlockParms);
printf("Executing blformat command...\n");
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("ROUND"))
			{
			if (isEqual(LocalToken,"YES"))
				rbp->roundFlag=1;
			else
				rbp->roundFlag=0;
			}
		if (parse_assignment2("NSITES"))
			{
			rbp->numSites=strtod(LocalToken,&dummy);
			gNumSites=rbp->numSites;
			printf ("Number of sites in sequences set to %li\n",rbp->numSites);
			}
		if (parse_assignment2("ULTRAMETRIC"))
			{
			if (isEqual(LocalToken,"YES"))
				{
				rbp->clockFmt=1;
				}
			else
				rbp->clockFmt=0;
			}
		if (parse_assignment2("LENGTHS"))
			{
			if (isEqual(LocalToken,"TOTAL"))
				{
				rbp->lengthFmt=0;
				printf("Branch lengths assumed to be in units of raw numbers of substitutions\n");
				}
			else
			  if (isEqual(LocalToken,"PERSITE"))
				{
				rbp->lengthFmt=1;
				printf("Branch lengths assumed to be in units of numbers of substitutions per site\n");
				}			
			}

		}



	if (gNexDataPtr->isTrees )
		{
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			if (rbp->lengthFmt  == 1)  /* if lengths are per site, convert to total numbers of substs. */
				traverseMultiplyLength(thisTree->root, (double)rbp->numSites,rbp->roundFlag);
			if (rbp->lengthFmt  == 0 && rbp->roundFlag)  
				traverseMultiplyLength(thisTree->root, 1,rbp->roundFlag); // this just forces a rounding by default in case stupid user inputs such
			if (rbp->clockFmt == 1)	   /* if tree is ultrametric, calculate times on scale of [0,1] */
				{
//				rootToTips(thisTree->root,0.0); ...checks for ultrametricity... for debugging mostly
				convert_branchlength_to_time(thisTree->root);
				thisTree->timesAvailable=1;
				thisTree->method=USER;	/* save the fact that USER supplied times */
				}
				
			}
		}
	if (rbp->lengthFmt  == 0 && rbp->roundFlag)  
		printf("All branch lengths were rounded to the nearest integer\n");
	if (rbp->lengthFmt == 1)
		{
		printf("All branch lengths multipled by the %li sites in the sequence\n",rbp->numSites);
		if (rbp->roundFlag)
			printf("Branch lengths rounded to nearest integer\n(This may not be a good idea when using CALIBRATE on ultrametric user supplied input trees. Use ROUND=NO then)\n");
		else
			printf("Branch lengths not rounded on input\n(If doing DIVTIME, you should SET ROUND=YES, the default)\n");
		}
	if (rbp->clockFmt == 1)
		printf("Tree is assumed to be ultrametric (Use CALIBRATE command to scale it)\n");


	return;
}
/****************************************************************/
void doFormat(void)
{

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (isEqual(aTokenPtr,"INTERLEAVE"))
			gNexDataPtr->Intlvflag=1;
		if (parse_assignment2("MISSING"))
			gNexDataPtr->missingchar = *LocalToken;
		if (parse_assignment2("GAP"))
			gNexDataPtr->gapchar = *LocalToken;
		if (parse_assignment2("MATCHCHAR"))
			gNexDataPtr->matchchar = *LocalToken;
		}
/* Report */
	if (gNexDataPtr->RateBlockParms.verbose)
		{
		printf("Missing character = %c\n",gNexDataPtr->missingchar);
		printf("Gap character = %c\n",gNexDataPtr->gapchar);
		printf("Match character = %c\n",gNexDataPtr->matchchar);
		}

	return;
}
/*----------------------------------------------------------------*/
void doSetCommand(void)
{
extern int gPowellTrace;
extern long gNumSites;
char  *dummy;
long aSeed;
struct RBP * rbp;
rbp=&(gNexDataPtr->RateBlockParms);
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("MINDURFACTOR"))
			rbp->minDurFactor=strtod(LocalToken,&dummy);
		if (parse_assignment2("MINRATEFACTOR"))
			rbp->minRateFactor=strtod(LocalToken,&dummy);
		if (parse_assignment2("NEIGHBORPENALTY"))
		    {
		    if (isEqual(LocalToken, "YES"))
			rbp->NeighborPenalty=1;
		    if (isEqual(LocalToken, "NO"))
			rbp->NeighborPenalty=0;
		    }
		if (parse_assignment2("PENALTY"))
		    {
		    if (isEqual(LocalToken, "ADD"))
			rbp->PenaltyType=0;
		    if (isEqual(LocalToken, "LOG"))
			rbp->PenaltyType=1;
		    }
		if (parse_assignment2("RATES"))
		    {
		    if (isEqual(LocalToken, "EQUAL"))
			rbp->RatesAreGamma=0;
		    else
		      if (isEqual(LocalToken, "GAMMA"))
			rbp->RatesAreGamma=1;
		    }
		if (parse_assignment2("ACTIVEEPSILON"))
			rbp->activeEpsilon=strtod(LocalToken,&dummy);
		if (parse_assignment2("SHAPE"))
			rbp->alpha=strtod(LocalToken,&dummy);
		if (parse_assignment2("BARRIERTOL"))
			rbp->barrierTol=strtod(LocalToken,&dummy);
		if (parse_assignment2("MAXITER"))
			rbp->maxIter=strtod(LocalToken,&dummy);
		if (parse_assignment2("MAXBARRIERITER"))
			rbp->maxBarrierIter=strtod(LocalToken,&dummy);
		if (parse_assignment2("INITBARRIERFACTOR"))
			rbp->initBarrierFactor=strtod(LocalToken,&dummy);
		if (parse_assignment2("BARRIERMULTIPLIER"))
			rbp->barrierMultiplier=strtod(LocalToken,&dummy);
		if (parse_assignment2("LINMINOFFSET"))
			rbp->linminOffset=strtod(LocalToken,&dummy);
		if (parse_assignment2("CONTRACTFACTOR"))
			rbp->contractFactor=strtod(LocalToken,&dummy);
		if (parse_assignment2("MAXCONTRACTITER"))
			rbp->maxContractIter=strtod(LocalToken,&dummy);
		if (parse_assignment2("FTOL"))
			{
			rbp->ftol=strtod(LocalToken,&dummy);
			printf("Parameter ftol set to %g\n",rbp->ftol);
			}
		if (parse_assignment2("SHOWCONVERGENCE"))
		    {
		    if (isEqual(LocalToken, "YES"))
			gNexDataPtr->RateBlockParms.showConvergence=1;
		    else
			gNexDataPtr->RateBlockParms.showConvergence=0;
		    }
		if (parse_assignment2("SHOWGRADIENT"))
		    {
		    if (isEqual(LocalToken, "YES"))
			gNexDataPtr->RateBlockParms.showGradient=1;
		    else
			gNexDataPtr->RateBlockParms.showGradient=0;
		    }
		if (parse_assignment2("CHECKGRADIENT"))
		    {
		    if (isEqual(LocalToken, "YES"))
				gNexDataPtr->RateBlockParms.checkGradient=1;
		    else
				gNexDataPtr->RateBlockParms.checkGradient=0;
		    }
		if (parse_assignment2("CLAMPROOT"))
		    {
		    if (isEqual(LocalToken, "YES"))
			gNexDataPtr->RateBlockParms.clampRoot=1;
		    else
			gNexDataPtr->RateBlockParms.clampRoot=0;
		    }
		if (parse_assignment2("TRACE"))
		    {
		    if (isEqual(LocalToken, "YES") || isEqual(LocalToken, "ON"))
			gPowellTrace=1;
		    else
			gPowellTrace=0;
		    }
		if (parse_assignment2("VERBOSE"))
			{
			gNexDataPtr->RateBlockParms.verbose=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("LOCAL_FACTOR"))
			gNexDataPtr->RateBlockParms.local_factor=strtod(LocalToken,&dummy);
		if (parse_assignment2("PERTURB_FACTOR"))
			gNexDataPtr->RateBlockParms.perturb_factor=strtod(LocalToken,&dummy);
		if (parse_assignment2("NPEXP"))
			{
			gNexDataPtr->RateBlockParms.npexp=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("SMOOTHING"))
			{
			gNexDataPtr->RateBlockParms.smoothing=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("SITES"))
			{
	     		if (isEqual(LocalToken,"ALL"))  
				doSitesCommand(0);
	     		if (isEqual(LocalToken,"EXCLUDE12"))  
				doSitesCommand(1);
	     		if (isEqual(LocalToken,"EXCLUDE3"))  
				doSitesCommand(3);

			}
		if (parse_assignment2("SEED"))
			{
			aSeed=strtod(LocalToken,&dummy);
			srand(aSeed);
			printf("Setting seed for random numbers to %i\n",aSeed);
			}
		if (parse_assignment2("NUM_RESTARTS"))
			gNexDataPtr->RateBlockParms.num_restarts=strtod(LocalToken,&dummy);
		if (parse_assignment2("NUM_RATE_GUESSES"))
			gNexDataPtr->RateBlockParms.num_rate_guesses=strtod(LocalToken,&dummy);
		if (parse_assignment2("NUM_TIME_GUESSES"))
			gNexDataPtr->RateBlockParms.num_time_guesses=strtod(LocalToken,&dummy);
		}





return;
}

/*----------------------------------------------------------------*/
void doMatrixGeneral(void)

/* reads a NEXUS data matrix and stores in global dmatrix as a set of strings corresponding 
to the rows of the matrix; interleaved and whitespace is allowed */

#define isNewLine(c) (((c)=='\n') || ((c)=='\r'))		
{
	int 	
			taxon, 
			saveix=0, 
			firstTaxon=1, 
			lastTaxon, 
			nextTaxon;
						
	double 	z;
	StrListPtr DM, TL;	
	
    //printf("\n\nWARNING: reading DATA block is under construction!\n");

	if (gNexDataPtr->NTaxa == 0 || gNexDataPtr->NChars ==0)
	    fatal("Must specify matrix dimensions prior to reading matrix");

	lastTaxon=gNexDataPtr->NTaxa;	/* get from global ntaxa */
	firstTaxon=1;
	taxon=1;
	lastTaxon=gNexDataPtr->NTaxa;	/*from global data */
	
/* initialize the global lists */

	DM=gNexDataPtr->DMList = newStrListN(gNexDataPtr->NTaxa);
	TL=gNexDataPtr->TaxaList=newStrList();



/* Loop through all tokens in the command */
	gNewLine=1;
	for (;;) 
	    {
	    aTokenPtr=nextToken();
	    if (isNewLine(*aTokenPtr))
		continue;
	    if (isEqual(aTokenPtr, ";"))
		break;
	    appendStrList(TL,aTokenPtr); /*Store the taxon label */
	    aTokenPtr=nextToken();
	    while(!isNewLine(*aTokenPtr))
		{
		catkthStr(DM, aTokenPtr, (long)taxon);
		aTokenPtr=nextToken() ; /* skip over possible new lines */
		}
	    //printf("%s\t%s\n", getkthStr(TL,(long)(taxon)), getkthStr(DM,(long)(taxon)));
	    ++taxon;
	    }



	gNewLine=0;
	checkMatrix();
	gNexDataPtr->isTaxa=1;	/* Set flag showing that labels read */
	gNexDataPtr->isChars=1;
	return;
}

/*----------------------------------------------------------------*/
void doMatrix(void)

/* reads a NEXUS data matrix and stores in global dmatrix as a set of strings corresponding 
to the rows of the matrix; interleaved and whitespace is allowed, but taxa must
be in same order as in TAXA block; polymorphisms as sets are simply stored as in the originial 
matrix--not parsed in any way--watch out! */

{
	int 	
			taxon, 
			saveix=0, 
			firstTaxon=1, 
			lastTaxon, 
			nextTaxon;
						
	double 	z;
	StrListPtr DM;	
	

	lastTaxon=gNexDataPtr->NTaxa;	/* get from global ntaxa */
	firstTaxon=1;
	taxon=1;
	lastTaxon=gNexDataPtr->NTaxa;	/*from global data */
	
/* initialize the global list for DM */

	gNexDataPtr->DMList = newStrListN(gNexDataPtr->NTaxa);
	DM=gNexDataPtr->DMList;

/* Loop through all tokens in the command */


	for (;;)
		{

/* Get a token ...and test it for a variety of conditions */

		aTokenPtr=nextToken();
		
/* bail at end of command  */

		if (isEqual(aTokenPtr,";"))
			break;	
									
/* reset index when the first taxon name is encountered (e.g. in interleaved) and 
start at top of loop again with new token */ 

		if (isEqual(aTokenPtr,getkthStr(gNexDataPtr->TaxaList,(long)firstTaxon)))
			{
			taxon=1;	
			continue;
			}
			
/* if token is the NEXT taxon name then on next loop begin storing in next row */
/* ... otherwise store data on present row */
		
		if (taxon<lastTaxon)		
			{
			if (isEqual(aTokenPtr,getkthStr(gNexDataPtr->TaxaList,(long)(taxon+1)))) /* taxon+1 is the next taxon */

				++taxon;

			else

				catkthStr(DM, aTokenPtr, taxon);

			}

/*... but note that we have to be careful not to check past the end of the labels array */ 

		else			

			catkthStr(DM, aTokenPtr, taxon);

		
		}	/* end for */


	gNexDataPtr->isChars=1;	/* This flag is set if a matrix command is read; might cause
							problems if the matrix command is empty of a mstrix! */
	checkMatrix();
	return;
}

/*----------------------------------------------------------------*/
void doDimensions(void)
{
	char  *dummy;
	int verbose;
	verbose=gNexDataPtr->RateBlockParms.verbose;

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NCHAR"))
			{
			gNexDataPtr->NChars=(int)strtod(LocalToken,&dummy);
			if (verbose) printf("Number of characters in matrix = %i\n",gNexDataPtr->NChars);		
			}
		if (parse_assignment2("NTAX"))
			{
			gNexDataPtr->NTaxa=(int)strtod(LocalToken,&dummy);
			if (verbose) printf("Number of taxa in matrix = %i\n",gNexDataPtr->NTaxa);		
			}
		}
	return;
}
/*----------------------------------------------------------------*/
void doCharDimensions(void)
{
	char  *dummy;
	int verbose;
	verbose=gNexDataPtr->RateBlockParms.verbose;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NCHAR"))
			{
			gNexDataPtr->NChars=(int)strtod(LocalToken,&dummy);
			}
		}
	if (verbose) printf("Number of characters in matrix = %i\n",gNexDataPtr->NChars);		
	return;
}
/*----------------------------------------------------------------*/
void doTaxDimensions(void)
{
	char  *dummy;
	int verbose;
	verbose=gNexDataPtr->RateBlockParms.verbose;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NTAX"))
			{
			gNexDataPtr->NTaxa=(int)strtod(LocalToken,&dummy);
			}
		}
	if (verbose) printf("Number of taxa in matrix = %i\n",gNexDataPtr->NTaxa);		
	return;
}
/*----------------------------------------------------------------*/

void doTaxLabels(void)
{
	int ix=0,len;

	gNexDataPtr->TaxaList=newStrList();

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		appendStrList(gNexDataPtr->TaxaList,aTokenPtr);
		}
	if (lengthList(gNexDataPtr->TaxaList)<gNexDataPtr->NTaxa) 
		fatal ("Too few taxon labels");
	else
		gNexDataPtr->isTaxa=1;	/* Set flag showing that labels read */
	return;
}


/**************************************************************/
void doTranslateCommand(void)

/* Currently the program WILL read trees with numbers as taxon names and just store those numbers
 * as strings.  All we have to do below is recurse and replace tip numbers with stuff from the list
 * below...so do it someday
 */

{
	gNexDataPtr->TransList=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	/* if its not a ';' it should be the number*/
		{
	
		if (  (!isdigit(*aTokenPtr)) && (!isEqual(aTokenPtr,","))  )
			appendStrList(gNexDataPtr->TransList,aTokenPtr); /* store the label */
		}
	gNexDataPtr->isTranslate=1;	/* set flag */
	printf("Trees stored WITH translation table\n");		
	return;
}

/****************************************************************/

void doTreeCommand(void)
{
		long size;
		extern int curTree;
		int flag=0;
		char *stemp, *tree_name, *TD, *p;
		PtrList aTreeList;

	if (gNexDataPtr->inTrees == NULL)  /* if this is the first tree */
		{
		gNexDataPtr->inTrees=pNewListAlt(sizeof(struct treetype));
		aTreeList=gNexDataPtr->inTrees;
		}
	else
		{
						/* if a later tree */
		aTreeList=pListAddNode(gNexDataPtr->inTrees,sizeof(struct treetype));
		if (aTreeList==NULL)fatal("Couldn't allocate tree list properly");
		}
/*printf("%li\n", (long)aTreeList);*/
	TD=makeEmptyString();
	tree_name=makeEmptyString();
	if ( isEqual(aTokenPtr=nextToken(),";")) return;

	if (isEqual(aTokenPtr,"*")) aTokenPtr=nextToken(); /* first token might be an
			asterisk; if so ignore and get the next token */

	concat(&tree_name,aTokenPtr);	/* this token should be the tree label */
	(void)appendStrList(gNexDataPtr->TDLabelList,tree_name);
	if (  !isEqual(aTokenPtr=nextToken(),"=")) 
			return; /* error; '=' should be next */

p=strchr(bufPtr,';');
if (p)
	{
	size=(long)(p-bufPtr);
/*printf(":%10.10s:%li:%c:\n",bufPtr,size,*p);*/
	*p='\0';
	TD=DupStr(bufPtr);
	(void)appendStrList(gNexDataPtr->TDList,TD);
	bufPtr=++p;
	}
#if 0
	while (!isEqual(aTokenPtr=nextToken(),";"))
		concat(&TD,aTokenPtr);	/* get the tree string */
	(void)appendStrList(gNexDataPtr->TDList,TD);
#endif			
	Tree_Initialize(aTreeList->item, TD, tree_name);
	if (gNexDataPtr->RateBlockParms.verbose) printf("Reading tree %s\n",tree_name);
	++curTree;
	gNexDataPtr->NumTrees=curTree;
	gNexDataPtr->isTrees=1;
/*print_tree_list(gNexDataPtr->inTrees);*/
	myFree(tree_name);
	myFree(TD);

	return;						
}
/****************************************************************/

void doDistanceCommand(void)
{
	StrListPtr aTaxaList;


	if ((gNexDataPtr->isChars==0) || ( gNexDataPtr->isTaxa==0))
		return;	/* don't have the right data from NEXUS file, so bail */


	aTaxaList=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	
			appendStrList(aTaxaList,aTokenPtr); /* store the label */
	doDistance(aTaxaList);
	freeStrList(aTaxaList);
	return;

}
/****************************************************************/
void doLengthMultiplyCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD;
	int roundflag=1; /* force it to round to nearest integer for now */	
	double multiplier=0.0;

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		multiplier=strtod(aTokenPtr, &dummy);
		}
		printf("All branch lengths multipled by AND ROUNDED TO NEAREST INTEGER %f\n", multiplier);
	if ((gNexDataPtr->isTrees) && (multiplier >  0.0) )
		{

		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			traverseMultiplyLength(thisTree->root, multiplier,roundflag);
			}
		}
	return;						
}
/****************************************************************/
static void doClusterHistogramCommand(void)

/* Default mode is to to print a cumulative histogram of unrooted cluster sizes across
 * a set of trees.  Also bins this histogram into two larger bins: "shallow" clades
 * of size 2-3 and "deep" clusters of size >3.
 * 
 * A cluster size is the size of the smaller partition in any bipartition of the taxa
 * 
 * If option NORMALIZE=YES is invoked,  then the program expects a set of N model
 * trees interleaved with a set of N strict consenus trees made from consensing the model
 * trees and the estimated trees.  It divides the number of clades in the latter by the nuber
 * of clades in the former and spits this out.  In other words,  the list goes 
 *	model tree1
 *	consensus tree1
 *	model tree2
 *	consensus tree 2,  etc.
 * 
 * NB.  All trees are assumed to be of the same size! Unprdictable results otherwise.
 * NB.  The clade of ALL taxa is ignored.
 */

{
	float b0=0.0, b1=0.0;
	PtrList lnode;
	TREE thisTree, modelTree, strictTree;
	char * dummy, *taxon, *TD;
	
	double multiplier=0.0;
	long * histo = NULL,*histo1, *histo2, TSize, bin1[2], bin2[2];
	double * cumHisto=NULL;
	int arraySize,i,nTrees,normFlag=0, ix;

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NORMALIZE"))
			{
			if (isEqual(LocalToken,"YES"))
				normFlag=1;
			else
				normFlag=0;
			}
		}
	if (normFlag==0)
	    {
	    if (gNexDataPtr->isTrees )
		    {
		    nTrees=pLengthList(gNexDataPtr->inTrees);
		    lnode=gNexDataPtr->inTrees;
		    LISTLOOP (lnode)
			    {
			    thisTree=lnode->item;
			    if (histo == NULL)
				    {
				    TSize=numdesc(thisTree->root);
				    arraySize=floor(LOG2(TSize))-1;
				    if (arraySize>0)
					{
					histo=(long *)myMalloc(arraySize*sizeof(long));
					cumHisto=(double *)myMalloc(arraySize*sizeof(double));
					for (i=0;i<arraySize;i++)
					    cumHisto[i]=0.0;
					}
				    }
			    for (i=0;i<arraySize;i++)
				histo[i]=0; 
			    ClusterHistogram(thisTree->root,histo,TSize);
			    for (i=0;i<arraySize;i++)
				cumHisto[i]+=histo[i];
			    }
		    printf("[\nMean histogram of cluster sizes\nSize of tree=%li\n\n", (long)TSize);
		    printf("Number of trees=%i\n",nTrees);
		    for (i=0;i<arraySize-1;i++)
		      if (i==arraySize -2 ) 
			printf("%li-%li:\t\t%f\n", 
			    (long)pow(2, i+1), (long)pow(2, i+2), cumHisto[i]/nTrees);
		      else
			printf("%li-%li:\t\t%f\n", 
			    (long)pow(2, i+1), (long)pow(2, i+2)-1, cumHisto[i]/nTrees);
		    printf("\n]\n");
    
		    }
	    else
		    printf("No trees currently in memory\n");
	    }
	else if (normFlag==1)
	    {
	    if (gNexDataPtr->isTrees )
		    {
		    nTrees=pLengthList(gNexDataPtr->inTrees)/2; /* each list is half as long */
		    lnode=gNexDataPtr->inTrees;
		    thisTree=lnode->item;
		    TSize=numdesc(thisTree->root);
		    arraySize=floor(LOG2(TSize))-1;
		    if (arraySize>0)
			{
			histo1=(long *)myMalloc(arraySize*sizeof(long));
			histo2=(long *)myMalloc(arraySize*sizeof(long));
			cumHisto=(double *)myMalloc(arraySize*sizeof(double));
			for (i=0;i<arraySize;i++)
			    cumHisto[i]=0.0;
			}

		    for (ix=1;ix<=nTrees;ix++) 
			    {
			    modelTree=lnode->item;
			    strictTree=lnode->next->item;
			    bin1[0]=0;bin1[1]=0;
			    bin2[0]=0;bin2[1]=0;
			    for (i=0;i<arraySize;i++)
				{
				histo1[i]=0;
				histo2[i]=0;
				}
			    ClusterHistogram(modelTree->root,histo1,TSize);
			    ClusterHistogram(strictTree->root,histo2,TSize);
			    for (i=0;i<arraySize-1;i++)
				{
				if (!((histo1[i] == 0) && (histo2[i] == 0)))
					{ 
					cumHisto[i]+=(float)histo2[i]/histo1[i];
	    /* IMPORTANT --> */		if (i==0) /*  "shallow" bins */
					    {
					    bin1[0]+=histo1[i];
					    bin2[0]+=histo2[i];
					    }
					else  /*..remaining "not shallow" bins */
					    {
					    bin1[1]+=histo1[i];
					    bin2[1]+=histo2[i];
					    }
					}
				}
			    b0+=(float)bin2[0]/(bin1[0]);
			    b1+=(float)bin2[1]/(bin1[1]);
			    lnode=lnode->next->next;  /* skip two */
			    }
		    printf("[\nNormalized histogram of cluster sizes\nSize of tree=%li\n\n", (long)TSize);
		    for (i=0;i<arraySize-1;i++)
		      if (i==arraySize -2 ) 
			printf("%li-%li:\t\t%f\n", 
			    (long)pow(2, i+1), (long)pow(2, i+2), cumHisto[i]/nTrees);
		      else
			printf("%li-%li:\t\t%f\n", 
			    (long)pow(2, i+1), (long)pow(2, i+2)-1, cumHisto[i]/nTrees);
		    
		    printf("[\nSummary cluster sizes\n");
		    printf("Shallow (2-3):%f\nNot shallow (4..):%f\n",
				 b0/nTrees, b1/nTrees);
		    printf("Deep(%li-%li):%f\n",
			(long)pow(2, arraySize-1), (long)pow(2,  arraySize),cumHisto[arraySize-2]/nTrees);
		    printf("\n]\n");
    
		    }
	    else
		    printf("No trees currently in memory\n");
	    }
	return;						
}
/****************************************************************/
static void doClearTreesCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD;
	
	double multiplier=0.0;

	if (gNexDataPtr->isTrees )
		{
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			Tree_Destructor(thisTree);
			}
		freepList(gNexDataPtr->inTrees);
		gNexDataPtr->inTrees=NULL;
		gNexDataPtr->isTrees=0;
		printf("All trees cleared from memory\n");
		}
	else
		printf("No trees currently in memory\n");
	return;						
}
/****************************************************************/
static void doCollapseCommand(void)

/* Collapses any zero-length branches to polytomies.

	** NB.! A BUG EXISTS: If an internal node has a name, but COLLAPSE
	is run, the name may get overwritten!
*/

{
	PtrList lnode;
	TREE thisTree;
	int num;

	if (gNexDataPtr->isTrees )
		{
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			num=0;
			thisTree=lnode->item;
			while(any_zero_internal_branches(thisTree->root))
				{
				collapse_zero(thisTree->root);
			/* have to call this repeatedly because collapse
				only works on first zero branch; then the tree
				is different and a node is missing, so...*/
				++num;
				}
			printf("%i zero-length branches collapsed\n",num);
			thisTree->numBranches=numBranches(thisTree->root);
			}
		}
	else
		printf("No trees currently in memory\n");
	return;						
}
/****************************************************************/
void doBSCommand(void)
{
	char * dummy, *taxon, *TD;
	

	gNexDataPtr->RateBlockParms.isBS=1;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NREPS"))
			gNexDataPtr->RateBlockParms.NReps=strtod(LocalToken,&dummy);
		if (parse_assignment2("SEED"))
			gNexDataPtr->RateBlockParms.seed=strtod(LocalToken,&dummy);
		}

		return;						
}
/****************************************************************/
static void doSuperCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD;
	
	long numTrees, numTaxa,treeNum,charNum,taxNum;
	int j, ii=1, ll, mm, icur, numInt,method, nn, wtFlag=0;
	int maxClades=0, maxTaxa;
	NODETYPE *root;
/* Fix these fixed length arrays!! */

	char **matrix /* [MAXTAX][MAXCLADES] */;
	float *wtset;
	StrListPtr aTaxaList,firstTaxaList;
	


	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("METHOD"))
			{
			if (isEqual(LocalToken,"BAUM"))
				method=0;
			else
				method=1; /* PURVIS */
			}
		if (parse_assignment2("WEIGHTS"))
			{
			if (isEqual(LocalToken,"YES"))
				wtFlag=1;
			else
				wtFlag=0;
			}
		}
	if (gNexDataPtr->isTrees)
		{

		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode) /* sets up a list of all unique taxa names
				    in taxaList[1] */
			{
			thisTree=lnode->item;
			root=thisTree->root;
			aTaxaList=newStrList();
			TreeToTaxaList(root, aTaxaList);
			if (ii==1)
			    firstTaxaList=aTaxaList;
			maxClades+=numIntNodes(root)-1; /* allow one char (clade)
			    for each interior node in tree (less the root)*/
			if (ii>1)
				{
				glomStrLists(firstTaxaList, aTaxaList);
				freeStrList(aTaxaList);
				}
			++ii;
			}
		maxTaxa=lengthList(firstTaxaList);

		matrix=(char **)myMalloc(maxTaxa*sizeof(char *));
		wtset=(float *)myMalloc(maxClades*sizeof(float));
		for (ll=0;ll<maxTaxa;ll++)  /* initialize the ABC matrix */
		    {
		    matrix[ll]=(char *)myMalloc(maxClades*sizeof(char));
		    for (mm=0;mm<maxClades;mm++)
			matrix[ll][mm]='?';
		    }
		    
   
    
		if (method==0)  /* Do the Baum and Ragan supertree */
		    {
		    icur=0;   
		    gColumn=0;/* global must be set prior to following */
		    lnode=gNexDataPtr->inTrees;
		    LISTLOOP (lnode)
			    {
			    thisTree=lnode->item;
			    root=thisTree->root;
			    numInt=numIntNodes(root)-1;
			    aTaxaList=newStrList();
			    TreeToTaxaList(root, aTaxaList);
			    numTaxa=lengthList(aTaxaList);
			    for(ll=1;ll<=numTaxa;ll++)
				{
				taxon=getkthStr(aTaxaList, ll);
				mm=findMatchStr(firstTaxaList, taxon);
				if (mm)
				   for (j=0;j<numInt;j++)
				    matrix[mm-1][icur+j]='0';
				}
			    j++;
			    icur+=numInt;
			    freeStrList(aTaxaList);
			    ABCSuperTree(root, firstTaxaList, matrix, wtset);			    
			    }
			    
	
		    printf("[Baum,  Ragan Supertree]\n\n");
		    printNexus(maxTaxa, maxClades,firstTaxaList, matrix );
		    if (wtFlag)
			{
			printf("weights ");
			for (nn=0;nn<maxClades-1;nn++)
				{
				if ((nn>0)&& ((nn/10)==(nn/10.0)))
				    printf("\n");
				printf("%5.2f:%i, ", wtset[nn],nn+1);
				}
			printf("%5.2f:%i;\n", wtset[maxClades-1],maxClades);
			}
		    }
		
  
		if (method==1) /* Do the Purvis kind of supertree */
		    {
		    for (ll=0;ll<maxTaxa;ll++)  /* initialize the ABC matrix */
			for (mm=0;mm<maxClades;mm++)
			    matrix[ll][mm]='?';
		    gColumn=0;/* global must be set prior to following */
		    lnode=gNexDataPtr->inTrees;
		    LISTLOOP (lnode)
			    {
			    thisTree=lnode->item;
			    root=thisTree->root;
			    ABCSuperTreePurvis(root, firstTaxaList, matrix, wtset);			    
			    }
			    
	
		    printf("[Purvis Supertree]\n\n");
		    printNexus(maxTaxa, maxClades,firstTaxaList, matrix );
		    }
	    freeStrList(firstTaxaList);

	    treeNum=1;   
	    charNum=1;
	    lnode=gNexDataPtr->inTrees;
	    printf ("begin paup;\n");
	    LISTLOOP (lnode)
			    {
			    thisTree=lnode->item;
			    root=thisTree->root;
			    numInt=numIntNodes(root)-1;
			    numTaxa=numdesc(root);
			    if (numInt>0)
			   	printf("charset input%i = %i-%i;\n",treeNum++,charNum,charNum+numInt-1);
			    else
				printf("[skipping input%i (no clades)]\n",treeNum++);
			    charNum+=numInt;
			    }
	    printf("end;\n");
	    }

	return;						

}  
/****************************************************************/
void printNexus(int ntaxa, int nchars, 	
	StrListPtr taxaList,  char **matrix)
{
int ll, mm;
printf("#Nexus\n");
printf("Begin data;\ndimensions ntax=%i nchar=%i;\n", ntaxa, nchars);
printf("format symbols=\"01\" missing=?;\nmatrix\n");

		for (ll=0;ll<ntaxa;ll++)
		    {
		    printf("%s\t", getkthStr(taxaList, ll+1));
		    for (mm=0;mm<nchars;mm++)
			printf("%c", matrix[ll][mm]);
		    printf("\n");
		    }
printf(";\nend;\n");
return;    
}

/****************************************************************/

void doConstrain_TimeCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD;
	
	long numTrees;
	int j, removeAll=0,fixAll=0;
	int root_id=0;
	int flagMax=-1,flagMin=-1;
	double maxAge=1.0e20, minAge=0.0; /* standard defaults; same as newnode*/
	NODETYPE *root, *found_node;

	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("REMOVE"))
			if (isEqual(LocalToken,"ALL"))
				removeAll=1;
		
		if (parse_assignment2("TAXON"))
			{
			if (isdigit(*LocalToken))
			    root_id=strtod(LocalToken,&dummy);
			else
			    taxon=DupStr(LocalToken);
			}
		if (parse_assignment2("MAXAGE"))
			{
			if (isEqual(LocalToken,"NONE"))
				flagMax=0;
			else
				{
				maxAge=strtod(LocalToken,&dummy);
				flagMax=1;
				}
			}
		if (parse_assignment2("MINAGE"))
			{
			if (isEqual(LocalToken,"NONE"))
				flagMin=0;
			else
				{
				minAge=strtod(LocalToken,&dummy);
				flagMin=1;
				}
			}
		if (parse_assignment2("MAX_AGE"))
			{
			if (isEqual(LocalToken,"NONE"))
				flagMax=0;
			else
				{
				maxAge=strtod(LocalToken,&dummy);
				flagMax=1;
				}
			}
		if (parse_assignment2("MIN_AGE"))
			{
			if (isEqual(LocalToken,"NONE"))
				flagMin=0;
			else
				{
				minAge=strtod(LocalToken,&dummy);
				flagMin=1;
				}
			}
		}
	if (gNexDataPtr->isTrees)
		{

		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (removeAll)
				unSetConstraints(root); /* remove all constraints from tree */
			else
			 if (root_id>0 || taxon != NULL) 
			  {
			  if (root_id>0)
			    found_node=find_id(root, root_id);
			  else
			    found_node=find_taxon_name(root,taxon);
			  if(found_node)
			     {
			     if (isFree(found_node))
				{
				if (flagMin==1)	/* if flags remain at -1, then no action taken */
					{
					printf("Setting minimum age constraint for taxon %s to %f\n",found_node->taxon_name,minAge);
					found_node->nodeIsConstrainedMin=1;
					found_node->minAge=minAge;
					}
				if (flagMin==0)
					{
					printf("Removing any minimum age constraint for taxon %s\n",found_node->taxon_name);
					found_node->nodeIsConstrainedMin=0;
					}
				if (flagMax==1)
					{
					printf("Setting maximum age constraint for taxon %s to %f\n",found_node->taxon_name,maxAge);
					found_node->nodeIsConstrainedMax=1;
					found_node->maxAge=maxAge;
					}
				if (flagMax==0)
					{
					printf("Removing any maximum age constraint for taxon %s\n",found_node->taxon_name);
					found_node->nodeIsConstrainedMax=0;
					}
				}
			     }
			  else
			     doGenericAlert("Cannot assign a constraint to a node that is fixed\nUse UNFIXAGE on this node");
			  }
			}
		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
void doCalibrateCommand(void)

/* Spits out the ages of all nodes relative to a given age of a specific node. On exit,
	all nodes will have been rescaled according to the one given node's age */

{
	PtrList lnode;
	TREE thisTree;
	char * dummy,  *TD, *profile_taxon;
	

static  double calAge=1.0;
	char * taxon=NULL;
	long numTrees;
	int i, j, profileFlag=0, id=1, ix, profile_node_id=0, taxon_id=0, calflag=0;
	double time, scalefactor=1.0; /*default calibration makes it equivalent to no correction */
	double ave, adev, sdev, var, skew, curt;
	double *data, *data_one_node;
	NODETYPE *root, *found_node, *node;
	StrListPtr aTaxaList, txPtr;
	PtrList nodeList, nLptr, mrcaPtr;
	NODETYPE *mrca, *s;
	extern long gNumSites;

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			{ 
			calflag=1;
			taxon=DupStr(LocalToken);
			}
		if (parse_assignment2("AGE"))
			calAge=strtod(LocalToken,&dummy);
		if (parse_assignment2("PROFILE_NODE"))
			{
			profile_taxon=DupStr(LocalToken);
			profileFlag=1;
			}
		
		}
	
	/*..............do the work...........*/	
		
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees and just spits out 
		 *calibrated ages.....
		 */
		ix=1;
		numTrees=pLengthList(gNexDataPtr->inTrees);
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			printf("Tree %i\n", ix);
			thisTree=lnode->item;
			if (thisTree->method != USER)
				doGenericAlert("WARNING: Calibrate command is designed for user-supplied ultrametric trees only!");
			root=thisTree->root;
			if (thisTree->timesAvailable)
			  {
			  if (taxon)
			    {
			    found_node=find_taxon_name(root,taxon);
			    printf("\nCalibrated Ages based on taxon %s @ age %f\n", 
				taxon, calAge);
			    }
			  else
			    {
			    found_node=root;
			    printf("\nCalibrated Ages based on ROOT age @ age %f\n", 
				calAge);
			    }
			  time=found_node->time;
			  scaleTree(root,calAge,found_node);
/*			  print_ages(root, time, calAge,thisTree->method);*/
			  }
			else
			  doGenericAlert("Times unavailable");
			++ix;
			}
		/*....works on a profile of identical topology trees and does
		/* summary stats by node...
		*/	
			
		if(profileFlag)
		    {
			data=(double*)myMalloc((numTrees+1)*sizeof(double));
			ix=1; /*init the index of trees */
			lnode=gNexDataPtr->inTrees;
			LISTLOOP (lnode)
				{
				thisTree=lnode->item;
				if (thisTree->timesAvailable)
				  {
				root=thisTree->root;
				/*...first get the scalefactor if calibration is internal*/
				found_node=find_taxon_name(root,taxon); 
					/* taxon is the calibrated taxon */
				if(found_node)
				    {
				    scalefactor=calAge/found_node->time;
				    }
				/*...now get node corresponding to id...*/
				else
				    printf("Couldn't find a calibration:Using root=1.0");
				
				node=find_taxon_name(root,profile_taxon); 
				if(node)
					{
					data[ix]=node->time*scalefactor;
					++ix;
// printf("Time of node for tree %i = %f\n",ix,data[ix]);
					}
				  }
				else
				  doGenericAlert("Times unavailable");
				}
			moment(data, numTrees, &ave, &adev,&sdev,
				&var, &skew, &curt); /* remember a 1-offset array */
			printf("Profile information for node across stored trees:\n"); 
			printf("Node=%s : Num trees=%i Mean time=%f  Standard deviation=%f Skewness=%f\n", 
					profile_taxon, numTrees,ave, sdev, skew);
			myFree(data);
		    }


		    


		}
	if (taxon)
		myFree(taxon);
	return;						
}
static void doShowAgeCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy,  *TD, *profile_taxon;
	

static  double calAge=1.0;
static  long iTree;
	char * taxon=NULL;
	long numTrees;
	int i, j, profileFlag=0, id=1, ix, profile_node_id=0, taxon_id=0, calflag=0,showNamed=0;
	double time, scalefactor=1.0; /*default calibration makes it equivalent to no correction */
	double ave, adev, sdev, var, skew, curt;
	double *data, *data_one_node;
	NODETYPE *root, *found_node, *node;
	StrListPtr aTaxaList, txPtr;
	PtrList nodeList, nLptr, mrcaPtr;
	NODETYPE *mrca, *s;

	
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TREE"))
			{
			iTree=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("SHOWNAMED")) // prints out a brief list of ages of named internal nodes only
			{
			if (isEqual(LocalToken,"YES"))
				showNamed=1;
			else
				showNamed=0;
			}
		}

	/*..............do the work...........*/	
		
	if (!gNexDataPtr->isTrees)
		printf("No input trees available\n");
	else
		{

		/*....works on any collection of trees and just spits out 
		 *calibrated ages.....
		 */
		ix=1;
		numTrees=pLengthList(gNexDataPtr->inTrees);
		lnode=gNexDataPtr->inTrees;
	        if (iTree>0) /* a specific tree was indicated */
			{
			if (iTree > pLengthList(lnode))
				{
				doGenericAlert("Invalid tree specified");
				return;
				}
			else
				{
				thisTree=(pListgetkthNode(lnode,iTree))->item;
				if (thisTree->timesAvailable)
					print_ages(thisTree->root,1.0,1.0,thisTree->method);
				}
			}
	        else
			LISTLOOP (lnode)
				{
				thisTree=lnode->item;
				root=thisTree->root;
				if (thisTree->timesAvailable)
				    if (showNamed)
					{
					printf("-----------------------------------------------------------------------------------------\n");
					printf("\nAges of internal named nodes only:\n");
					print_named_ages(root);
					}
				    else // the usual setting...
					{
					printf("-----------------------------------------------------------------------------------------\n");
					printf("Estimated ages and substitution rates for tree %s\n\n",thisTree->name);
					switch(thisTree->method)
						{
						case USER:printf("Reconstruction method: User-supplied ultrametric tree\n");break;
						case LaF:printf("Reconstruction method: Langley-Fitch (clock)\n");break;
						case NP:printf("Reconstruction method: NPRS\n");break;
						case PENLIKE:printf("Reconstruction method: Penalized likelihood\n");break;
						}
					printf("Named internal nodes indicated by [*]\n");
					printf("Rates are for branches subtending indicated node\n");
					printf("Rates are in units of substitutions per site per unit time\n\n");
					printf("\t\t\t\t  Constraints\t\t\t\tRates\n");
					printf("\tNode\t   Fix [Mod]\t  Min\t  Max\t  Age\t\tEstimated\tLocal\n");
					printf("-----------------------------------------------------------------------------------------\n");

					print_ages(root, 1.0,1.0,thisTree->method); /* use this more complex function to do something simple here! */
					printf("-----------------------------------------------------------------------------------------\n");

					summarize_rates(thisTree);

					}
				else
					doGenericAlert("Times and rates unavailable");
				++ix;
				}





		/*....works on a profile of identical topology trees and does
		/* summary stats by node...
		*/	
			
		if(profileFlag) /* currently never invokes...leave for use later perhaps */
		    {
			data=(double*)myMalloc((numTrees+1)*sizeof(double));
			ix=1; /*init the index of trees */
			lnode=gNexDataPtr->inTrees;
			LISTLOOP (lnode)
				{
				thisTree=lnode->item;
				root=thisTree->root;
				/*...first get the scalefactor if calibration is internal*/
				found_node=find_taxon_name(root,taxon); 
					/* taxon is the calibrated taxon */
				if(found_node)
				    {
				    scalefactor=calAge/found_node->time;
				    }
				/*...now get node corresponding to id...*/
				else
				    printf("Couldn't find a calibration:Using root=1.0");
				
				node=find_taxon_name(root,profile_taxon); 
				if(node)
					{
					data[ix]=node->time*scalefactor;
					++ix;
					}
				}
			moment(data, numTrees, &ave, &adev,&sdev,
				&var, &skew, &curt); /* remember a 1-offset array */
			printf("Profile information for node across stored trees:\n"); 
			printf("Node=%s : Num trees=%i Mean time=%f  Standard deviation=%f Skewness=%f\n", 
					profile_taxon, numTrees,ave, sdev, skew);
			myFree(data);
		    }


		    


		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
static void doRRLikeTestCommand(void)

/* 

Do a likelihood ratio test on the N clades descended from one node. Use the localmodel feature.

Comments: Interesting issues here, because we might want to use time constraint information and we might not.
When the focal node is the root of a clade without constraints at all, the DIVTIME routine will bail because
it will think the root cannot be estimated. In that case, it seems reasonable to set the age of the (local) root
node temporarily to 1.0. Alternatively for each node we might make a time-informed and a time-uninformed relative
rate test. In the latter, we ignore all time information of descendant nodes. In the former we take info into account,
and if the problem is inestimable, let the chips fall where they may! 
*/

{
	PtrList lnode;
	TREE thisTree,subtree;
	char * dummy, *taxon, *TD;
 	double Like0,Like1,LR;	
	long numTrees;
	int constrain,warn;
	int i, j, id=1, ix, stemFlag=1, model=0; //init. stemFlag to include stem lineage of each child
	int success=0,method,algorithm,nRates,allFlag=0;
	NODETYPE *root, *child, *node,*found_node,*saveNode;

	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees  */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root; // eventually do this for an arbitrary node!
			if (taxon)
			    {
				if (isEqual(taxon,"ALL"))
					{
					allFlag=1; // implement later
					}
				else
					{
			   		if (isEqual(taxon,"ROOT"))
						found_node=root;
			    		else
			    			found_node=find_taxon_name(root,taxon);
					}
			    	if (found_node)
					{
					warn=warnEstRoot(found_node);
					if (warn == 1)  // only when the divtime would normally bail do we consider this an unconstr search
						{
						constrain=0;	
						found_node->free=0;
						found_node->time=100.0; // unfixed so fix it at some arbitray age
						printf("...Insufficient time constraints present in RRLike...fixing age of %s at 100.0\n",found_node->taxon_name);
						}
					else
						constrain=1;
					saveNode=found_node->anc;
					subtree=Subtree_Initialize(thisTree,found_node); // pull off this subtree and work on it only
					algorithm=TN;
					method=LaF;
					nRates=1;
					doObjFunc(subtree,method,nRates,algorithm,&success);
					Like0=subtree->obj;

					child=found_node->firstdesc;
					SIBLOOP(child)
						{
						setLocalModel(child,model,stemFlag);
						++model;
						}
					nRates=model;
					algorithm=POWELL; // note we have to switch to POWELL for the Local clock methods!
					method=LFLOCAL;
					doObjFunc(subtree,method,nRates,algorithm,&success);
					Like1=subtree->obj;
					LR = -2 * (Like0-Like1);
					printf("\nRELATIVE RATE TEST USING LIKELIHOOD RATIO\n\n");
					printf("    Node       Clock   Non-clock     LR Stat    d.f.    Constrained?\n");
					printf("--------------------------------------------------------------------\n");
					printf("%8s    %8.2f    %8.2f    %8.2f    %4i        %4i\n",found_node->taxon_name,Like0,Like1,LR,nRates-1,constrain);
					found_node->anc=saveNode; // restore the subtree's link to original tree
					if (warn == 1)  
						{
						found_node->free=1;
						found_node->time=0.0; // reset values 
						}
					}
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
static void doLocalModelCommand(void)

/* 

Takes instructions on setting up a local clock model with finite number of rate parameters distributed across tree.
Parameters are indexed from 0..N-1, if there are N rates. Assigns all branches in clade 'Taxon' some index value. If
Taxon is a tip, then its subtending branch is used. If 'stem' is set to yes, then the stem branch is assigned as well
as all members of the clade.

 */

{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD;
	
	long numTrees;
	int i, j, id=1, ix, stemFlag=0, model;
	NODETYPE *root, *found_node, *node;

	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		if (parse_assignment2("STEM"))
			{
			if (isEqual(LocalToken,"YES"))
				stemFlag=1;
			else
				stemFlag=0;
			}
		if (parse_assignment2("RATEINDEX"))
			model=strtod(LocalToken,&dummy);

			if (isEqual(LocalToken,"ALL"))
				(void)preOrder(root,fixNodeAge); 
					/* force all nodes to have their ages fixed (hopefully ages are set somehow!) */
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees  */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (taxon)
			    {
			   	if (isEqual(taxon,"ROOT"))
					found_node=root;
			    	else
			    		found_node=find_taxon_name(root,taxon);
			    	if (found_node)
					{
					setLocalModel(found_node,model,stemFlag);
					}
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
void doSetAgeCommand(void)

/* Sets the age of any node in the tree, but this is transient if node is internal */
/* DO NOT PERMIT THE USE OF this command (now 'FixAge') on ultrametric trees...collides
   with 'Calibrate' command ...probably should just have one command, and prevent user
   from fixage on more than one node for ultrametric trees LATER...*/

{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD, *profile_taxon;
	
	long numTrees;
	int i, j, profileFlag=0, id=1, ix, profile_node_id=0, taxon_id=0, calflag=0;
	double age, time, scalefactor=1.0; /*default calibration makes it equivalent to no correction */
	double *data, *data_one_node;
	NODETYPE *root, *found_node, *node;
	StrListPtr aTaxaList, txPtr;
	PtrList nodeList, nLptr, mrcaPtr;
	NODETYPE *mrca, *s;
	struct RBP * rbp;
	rbp=&(gNexDataPtr->RateBlockParms);

	taxon=NULL;

	if (rbp->clockFmt)
		{
		doGenericAlert("Can't SETAGE on ultrametric trees; use CALIBRATE command instead!");
		return;	
		}

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		if (parse_assignment2("AGE"))
			age=strtod(LocalToken,&dummy);
		if (parse_assignment2("FIX"))
			if (isEqual(LocalToken,"ALL"))
				(void)preOrder(root,fixNodeAge); 
					/* force all nodes to have their ages fixed (hopefully ages are set somehow!) */
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees  */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (taxon)
			    {
			    if (isEqual(taxon,"ALL"))
				(void)preOrder(root,fixNodeAge); 
				/* If only say 'setAge taxon=all' then fix all times to whatever they are. 
				   Doesn't permit us to set all nodes to some one age--that'd be dumb*/
			    else
				{
			   	if (isEqual(taxon,"ROOT"))
					found_node=root;
			    	else
			    		found_node=find_taxon_name(root,taxon);
			    	if (found_node)
					{
					found_node->free=0;
			    		found_node->time=age;
					found_node->nodeIsConstrainedMax=0;
					found_node->nodeIsConstrainedMin=0; /* overrides preexisting min or max constraints ! */
					printf("Fixing age of %s at %f\n",taxon,age);
					printf(" (The age of this node will no longer be estimated.)\n");
					printf(" (This command overides any previous age constraints for this node.)\n");
					printf(" (The total number of fixed ages is now %i)\n",numFixedNodes(root));
					}
				}
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
static void doUnSetAgeCommand(void)

/* Frees the age of a node, which will subsequently be estimated */

{
	PtrList lnode;
	TREE thisTree;
	char *taxon;
	NODETYPE *root, *found_node, *node;

	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees  */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (taxon)
			    {
			    if (isEqual(taxon,"ALL"))
				(void)preOrder(root,unFixNodeAge);
			    else
				{
			   	 if (isEqual(taxon,"ROOT"))
					found_node=root;
			   	 else
			    		found_node=find_taxon_name(root,taxon);
			   	 if (found_node)
					{
					found_node->free=1;
					printf("Unfixing age of %s\n",taxon);
					printf(" (the age of this node WILL be estimated in future searches)\n");
					}
				}
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
static void doCovarionCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	int success, maxIterations,nParm, showProbs=0,showChanges=0,model=0,doRecon=0,doMarginals=0,estimate=0,simulate=0;
	double ftol=0.0001, linMinDelta=0.05;
	char *dummy;
	struct RBP * rbp;
double p1low,p1high,p2low,p2high;
int dim=8,likeSurface=0;
double (*obj_func_array[4])(double[]);
obj_func_array[0] = objBinaryTraitSymmetric;
obj_func_array[1] = objBinaryTrait;
obj_func_array[2] = objCovarion;
obj_func_array[3] = objCovarion; // this is the switch_1 model which uses same objective with a tweak to set two parms equal

	gNexDataPtr->RateBlockParms.estCov=0;
	gNexDataPtr->RateBlockParms.cov_brlens=0;
	rbp=&(gNexDataPtr->RateBlockParms);
	maxIterations=rbp->maxIter;
	int nstates=3;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("MODEL"))
		    {
		    if (isEqual(LocalToken,"BINARY_1"))
				model=0;
		    if (isEqual(LocalToken,"BINARY_2"))
				model=1;
		    if (isEqual(LocalToken,"SWITCH")) // deprecated
				model=2;
		    if (isEqual(LocalToken,"SWITCH_2"))
				model=2;
		    if (isEqual(LocalToken,"SWITCH_1"))
				model=3;
		    if (isEqual(LocalToken,"PRECURSOR_2"))
				model=2;
		    if (isEqual(LocalToken,"PRECURSOR_1"))
				model=3;
		    if (isEqual(LocalToken,"MULT_1"))
				model=4;

		    }
		if (parse_assignment2("MARGINALS"))
		    {
		    if (isEqual(LocalToken,"YES"))
				doMarginals=1;
		    else
				doMarginals=0;

		    }
		if (parse_assignment2("ANCSTATES"))
		    {
		    if (isEqual(LocalToken,"YES"))
				doRecon=1;
		    else
				doRecon=0;

		    }
		if (parse_assignment2("ESTIMATE"))
		    {
		    if (isEqual(LocalToken,"YES"))
				estimate=1;
		    else
				estimate=0;
		    }
		if (parse_assignment2("SIMULATE"))
		    {
		    if (isEqual(LocalToken,"YES"))
				simulate=1;
		    else
				simulate=0;
		    }
		if (parse_assignment2("SHOWPROBS"))
		    {
		    if (isEqual(LocalToken,"YES"))
				showProbs=1;
		    else
				showProbs=0;
		    }
		if (parse_assignment2("SHOWCHANGES"))
		    {
		    if (isEqual(LocalToken,"YES"))
				showChanges=1;
		    else
				showChanges=0;
		    }
		if (parse_assignment2("BRLENS"))
		    {
		    if (isEqual(LocalToken,"ONE"))
				gNexDataPtr->RateBlockParms.cov_brlens=1;
		    if (isEqual(LocalToken,"USER"))
				gNexDataPtr->RateBlockParms.cov_brlens=0;

		    }
		if (parse_assignment2("Q01"))
			gNexDataPtr->RateBlockParms.s_rate=strtod(LocalToken,&dummy);
		if (parse_assignment2("Q10"))
			gNexDataPtr->RateBlockParms.r_rate=strtod(LocalToken,&dummy);
		if (parse_assignment2("S_RATE"))
			gNexDataPtr->RateBlockParms.s_rate=strtod(LocalToken,&dummy);
		if (parse_assignment2("R_RATE"))
			gNexDataPtr->RateBlockParms.r_rate=strtod(LocalToken,&dummy);
			
			
			
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			if (simulate)
				simulateBinaryChar(thisTree, gNexDataPtr->RateBlockParms.s_rate, gNexDataPtr->RateBlockParms.r_rate );
//				simulatePrecursorChar(thisTree, gNexDataPtr->RateBlockParms.s_rate, gNexDataPtr->RateBlockParms.r_rate );
			else
				covarionOptimize(thisTree,&maxIterations, rbp->ftol,rbp->linminOffset,&success,model,doMarginals, estimate, doRecon );
			if (showProbs) printCovarion(thisTree->root,doMarginals);
			if (showChanges) printChanges(thisTree->root);
			}

		}
	return;		
}				
/****************************************************************/
static void doAncestralCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	int success, maxIterations,nParm;
	double ftol=0.0001, linMinDelta=0.05;
	char *dummy;
	struct RBP * rbp;
	rbp=&(gNexDataPtr->RateBlockParms);
	maxIterations=rbp->maxIter;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			printf("...doing ancestral state squared change parsimony optimization...\n");
			printf("Optimization parameters:\n  ftol...%f\n");
			ancestralOptimize(thisTree,&maxIterations, rbp->ftol,rbp->linminOffset,&success );
			printf("Node\t\tState\t\tAnc State\tDiff\tSign of Difference\n");
			printAncestral(thisTree->root);
			}

		}
	return;		
}				
/****************************************************************/
static void doContOptCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	int nRates=1, success, maxIterations,nParm;
	double ftol=0.0001, linMinDelta=0.05;
	char *dummy;
	maxIterations=1000;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("NRATES"))
			nRates=strtod(LocalToken,&dummy);
		}
	nParm=nRates+1;
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			printf("...doing continuous character rate optimization...\n");
			contOptimize(thisTree,nParm,&maxIterations, ftol,linMinDelta,&success );
			}

		}
	return;						
}

/****************************************************************/
void doPruneTaxonCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD, *profile_taxon;
	
	long numTrees;
	int i, j, id=1, ix,  taxon_id=0;
	NODETYPE *root, *found_node, *node;
	StrListPtr aTaxaList, txPtr;
	PtrList nodeList, nLptr, mrcaPtr;
	NODETYPE *mrca, *s;

	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (taxon)
			    {
			    found_node=find_taxon_name(root,taxon);
			    if (found_node)
				{
			    	thisTree->root=RemoveTaxon(thisTree,found_node);
				thisTree->numBranches=numBranches(thisTree->root);
				printf("Pruning taxon %s\n",taxon);
				}
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
static void doVCVCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD, *profile_taxon;
	double T;
	long i, j;
	NODETYPE *root, *found_node, *node;
	NODE a,b,c;
	PtrList nodeList;
	long lengthList;
	
	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (taxon)
			    {
			    found_node=find_taxon_name(root,taxon);
			    if (found_node)
					{
					nodeList=pNewList();
					TreeToTaxaPtrList(root,nodeList);
					lengthList=pLengthList(nodeList);
					for (i=1;i<=lengthList;i++)
						{
						a=(NODE)(pListgetkthNode(nodeList, i)->item);
						printf("%s\t",a->taxon_name);
						for (j=1;j<=lengthList;j++)
							{
							b=(NODE)(pListgetkthNode(nodeList, j)->item);
							c=mrca(a,b);
							T=pathLengthTime(c,found_node);
							printf("%f\t",T);
							}
						printf("\n");
						}
					freepList(nodeList);
					}
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
void doReRootCommand(void)
{
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD, *profile_taxon;
	
	long numTrees;
	int i, j, id=1, ix,  taxon_id=0, atNode=0;
	NODETYPE *root, *found_node, *node;
	StrListPtr aTaxaList, txPtr;
	PtrList nodeList, nLptr, mrcaPtr;
	NODETYPE *mrca, *s;

	taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			taxon=DupStr(LocalToken);
		if (parse_assignment2("ATNODE"))  // reroots at a node rather than maintaining a binary root
			if (isEqual(LocalToken,"YES"))
				atNode=1;
			else
				atNode=0;
		}
	
	if (gNexDataPtr->isTrees)
		{

		/*....works on any collection of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			if (taxon)
			    {
			    found_node=find_taxon_name(root,taxon);
			    if (atNode==0)
			    	thisTree->root=ReRoot(found_node);
			    if (atNode==1)
			    	thisTree->root=ReRoot2(found_node);
			    }
			}

		}
	if (taxon)
		myFree(taxon);
	return;						
}
/****************************************************************/
static void doBranchProfileCommand(void)
{
	extern long gNumSites;
	PtrList lnode;
	TREE thisTree;
	char * dummy, *taxon, *TD, *profile_taxon;
	
	long numTrees;
	int i, j, profileFlag=0, id=1, ix, profile_node_id=0, taxon_id=0, calflag=0,parmFlag=2;
	double calAge=1.0, time, scalefactor=1.0; /*default calibration makes it equivalent to no correction */
	double ave, adev, sdev, var, skew, curt,min=1e20,max=-1e20;
	double *data, *data_one_node;
	NODETYPE *root, *found_node, *node;
	StrListPtr aTaxaList, txPtr;
	PtrList nodeList, nLptr, mrcaPtr;
	NODETYPE *mrca, *s;

	taxon=NULL;
	profile_taxon=NULL;
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXON"))
			{
			profile_taxon=DupStr(LocalToken);
			}
		if (parse_assignment2("PARAMETER"))
			{
			if (isEqual(LocalToken,"LENGTH"))
				parmFlag=1;
			else if (isEqual(LocalToken,"AGE"))
				parmFlag=2;
			else if (isEqual(LocalToken,"RATE"))
				parmFlag=3;

			}
		
		}
	
	/*..............do the work...........*/	
		
	if (gNexDataPtr->isTrees)
		{
		numTrees=pLengthList(gNexDataPtr->inTrees);
		data=(double*)myMalloc((numTrees+1)*sizeof(double));
		ix=0; /*init the index of trees */
		lnode=gNexDataPtr->inTrees;
		LISTLOOP (lnode)
			{
			thisTree=lnode->item;
			root=thisTree->root;
			node=find_taxon_name(root,profile_taxon); 
			if(node)
				{
				++ix;
				switch (parmFlag)
					{
					case 1: data[ix]=node->length;break;
					case 2: data[ix]=node->time;break;
					case 3: data[ix]=node->estRate/gNumSites;break;
					}
			/*	printf("%f\n",data[ix]);*/
				}
			else
				{
				printf("WARNING! Profiled node not found on tree (may have been collapsed?)\n");
				}
			}
		for (i=1;i<=ix;i++) // remember ix may be less than numTrees if some nodes are collapsed and not found
			{
			if (data[i]>max)max=data[i];
			if (data[i]<min)min=data[i];
			}
		printf("Profile information for node %s across %i tree(s) out of %i trees:\n",profile_taxon,ix,numTrees); 
		if (ix >=2)
			{
			moment(data, ix, &ave, &adev,&sdev,
				&var, &skew, &curt); /* remember a 1-offset array */
			printf("Mean = %f  Std dev = %f Min = %f Max = %f\n", ave, sdev,min,max);
			}
		else
			printf("Profile cannot be obtained when number of trees with given node < 2\n");
		myFree(data);
		}
	if(profile_taxon)
		myFree(profile_taxon);
	return;						
}
/****************************************************************/
void doSimCommand(void)
{
float T_LF[25][25];
float T_PL[25][25];
NODETYPE *nodea, *nodeb, *nodec, *nodeInt;
float bx,cx,lb,lc,trueAge,rangeFactor;
int ksteps,ib,ic;  // ksteps an even number please
	PtrList lnode;
	TREE thisTree;
	int verbose = 0, resetSeed=0;
	float bestSmooth;
	extern NODETYPE * gRoot;

	char * 	dummy;
	
	double			*RepCount,*RepMean,*RepDominance,*RepFreq1Class,*RepFractMonophyletic,*RepMonophyleticSpecies;
	double  *X, *Y, * time1, *time2,*time3,*data1,*data2,*chiSqArray, *data3, 
		ang,mean,adev,LFsdev,NPsdev, sdev, var,skew,curt,av, 
		Kendall_var1, Kendall_var2, kappa, B, NN1,NN2,freq1class,dominance;
	double 	NPmean,LFmean,chiSqmean;
	int	whichBetter;
	double	NN;
	NODETYPE * root, *node, **markedNodes, **nodeArray, *saveTree;
	extern 	int gIndex;
	extern 	double chiSq; /* declared in ObjFunc */
	int	i,j,k,success1,success2,kk,jj,irepcount,TotalReps;

	long    MaxGroupSize,nMark=0,count,s1,s2,maxS,
		countTaxa,countExc,binSize=10,nTaxa,nNodes,size,size2,countMonotypes, countMonophyletic, 
		countMonophyleticSpecies,*h,*histo,*histo2,*histoMonophyletic,
		*histoTotal,*histo2Total,theSeed=1,*histoB,sizeB;
	int    	
		N0=1,
		stemFlag=0,
		exclusive=1,
		rndBranchDur=0,
		withReplace=1,
		diversemodel=1,
		Yule_flag=0,
		CharEvol_flag=0, 
		save_flag=1,
		nreps=1,
		nrepsPerTree=1,
		nrepsPerBrRate=1,
		ntaxa=10,
		interval=10,
		infinite_flag=0,
		ratemodel=1, 
		gradual_flag=1, 
		silentFlag=0;
	double 	
		speciation=1.0,
		speciation2=1.0,
		extinction=0.0,
		sampling_fraction=1.0,
		start_rate=1.0,
		change_rate=0.0,
		min_rate=0.1,
		max_rate=2.0,
		rate_transition=0.0, 
		T=1.0;

/*print_mem_dbg(__FILE__,__LINE__);*/
	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("TAXONOMY"))
			{
			if (isEqual(LocalToken,"EXCLUSIVE"))
				exclusive=1;
			if (isEqual(LocalToken,"NESTED"))
				exclusive=0;
			}
		if (parse_assignment2("RNDBRANCHDUR"))
			if (isEqual(LocalToken,"YES"))
				rndBranchDur=1;
			else
				rndBranchDur=0;
		if (parse_assignment2("WITHREPLACE"))
			if (isEqual(LocalToken,"YES"))
				withReplace=1;
			else
				withReplace=0;
		if (parse_assignment2("CHAREVOL"))
			if (isEqual(LocalToken,"YES"))
				CharEvol_flag=1;
			else
				CharEvol_flag=0;
		if (parse_assignment2("NREPSPERBRRATE"))
			nrepsPerBrRate=strtod(LocalToken,&dummy);
		if (parse_assignment2("NREPSPERTREE"))
			nrepsPerTree=strtod(LocalToken,&dummy);
		if (parse_assignment2("BINSIZE"))
			binSize=strtod(LocalToken,&dummy);
		if (parse_assignment2("NMARK"))
			nMark=strtod(LocalToken,&dummy);
		if (parse_assignment2("NREPS"))
			nreps=strtod(LocalToken,&dummy);
		if (parse_assignment2("NTAXA"))
			ntaxa=strtod(LocalToken,&dummy);
		if (parse_assignment2("SPECIATION"))
			speciation=strtod(LocalToken,&dummy);
		if (parse_assignment2("SPECIATION2"))
			speciation2=strtod(LocalToken,&dummy);
		if (parse_assignment2("EXTINCTION"))
			extinction=strtod(LocalToken,&dummy);
		if (parse_assignment2("SAMPLING"))
			sampling_fraction=strtod(LocalToken,&dummy);
		if (parse_assignment2("INTERVAL"))
			interval=strtod(LocalToken,&dummy);
		if (parse_assignment2("SEED"))
			{
			theSeed=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("STARTRATE"))
			start_rate=strtod(LocalToken,&dummy);
		if (parse_assignment2("CHANGERATE"))
			change_rate=strtod(LocalToken,&dummy);
		if (parse_assignment2("MINRATE"))
			min_rate=strtod(LocalToken,&dummy);
		if (parse_assignment2("MAXRATE"))
			max_rate=strtod(LocalToken,&dummy);
		if (parse_assignment2("RATETRANSITION"))
			rate_transition=strtod(LocalToken,&dummy);
		if (parse_assignment2("T"))
			{
			T=strtod(LocalToken,&dummy);
			}
		if (parse_assignment2("INFINITE"))
			if (isEqual(LocalToken,"NO"))
				infinite_flag=0;
			else
				infinite_flag=1;
		if (parse_assignment2("GRADUAL"))
			if (isEqual(LocalToken,"NO"))
				gradual_flag=0;
		if (parse_assignment2("SAVE_FLAG"))
			if (isEqual(LocalToken,"YES"))
				save_flag=1;
			else
				save_flag=0;
		if (parse_assignment2("SILENT"))
			if (isEqual(LocalToken,"YES"))
				silentFlag=1;
			else
				silentFlag=0;
		if (parse_assignment2("STEM"))
			if (isEqual(LocalToken,"YES"))
				stemFlag=1;
			else
				stemFlag=0;
		if (parse_assignment2("RATEMODEL"))
			{
			if (isEqual(LocalToken,"NORMAL"))
				ratemodel=1;
			if (isEqual(LocalToken,"AUTOCORR"))
				ratemodel=2;
			}
		if (parse_assignment2("DIVERSEMODEL"))
			{
			Yule_flag=1;
			if (isEqual(LocalToken,"YULE"))
				diversemodel=1;
			if (isEqual(LocalToken,"BDBACKNORMAL"))
				diversemodel=5;
			if (isEqual(LocalToken,"BDBACK"))
				diversemodel=2;
			if (isEqual(LocalToken,"YULE_C"))
				diversemodel=3;
			if (isEqual(LocalToken,"BDFORWARD")) 
				diversemodel=4;
			if (isEqual(LocalToken,"RY1997")) 
				diversemodel=6;
			if (isEqual(LocalToken,"YULE_SISTERS")) 
				diversemodel=7;

			}
		}

	if (diversemodel==7)
				printf("[Expected ratio of sister group diversities=%f]\n",exp(T*(speciation2-speciation)));

	srand(theSeed);
	if (theSeed==1)
		{
		if (!silentFlag)doGenericAlert("WARNING: YOU ARE USING A DEFAULT SEED FOR RANDOM NUMBERS");
		}
	if (!silentFlag)printf("\n\n** r8s simulation run **\n\n");
	verbose=gNexDataPtr->RateBlockParms.verbose;

	kappa=speciation*T;
	Kendall_var1=SQR(speciation)/(2*(exp(kappa)-1.));
	Kendall_var2=Kendall_var1*SQR(sinh(0.5*kappa)/(0.5*kappa));
	time1=(double*)myMalloc((ntaxa-2)*sizeof(double));
	time2=(double*)myMalloc((ntaxa-2)*sizeof(double));
	time3=(double*)myMalloc((ntaxa-2)*sizeof(double));
	data1=(double*)myMalloc((nreps+1)*sizeof(double)); /* 1-offset array */
	data2=(double*)myMalloc((nreps+1)*sizeof(double)); /* 1-offset array */
	data3=(double*)myMalloc((nreps+1)*sizeof(double)); /* 1-offset array */
	chiSqArray=(double*)myMalloc((nreps+1)*sizeof(double)); /* 1-offset array */

#define SIM_LOOP 0	/* for doing lots of simulations */
#if SIM_LOOP

printf("ChangeRate,Transition,ChiSq,LF,NP,Which\n");

/* for (jj=0;jj<=10;jj++) */ 
  for (kk=0;kk<=10;kk++)
	{
	/* change_rate=jj*max_rate/20.; */
	rate_transition=kk/10.0;
#endif




/* start simulating */

#if 0  /* TEST  */
	printf("Simulation of Yule_forward\n");
	for (i=1;i<=nreps;i++)
		{
		NN1=Yule_forward(speciation, T, &B,stemFlag);
		data1[i]=NN1;
		if (stemFlag) N0=1; else N0=2;	
		data2[i]=(NN1-N0)/B;
		data3[i]=(log(NN1)-log(N0))/T;
		}
	moment(data1,nreps,&mean,&adev,&sdev,&var,&skew,&curt);
	printf("Test of Yule Forward routine: mean=%fvar=%f\n",mean,var);
	moment(data2,nreps,&mean,&adev,&sdev,&var,&skew,&curt);
	printf("Test of K-infinite estimator: mean=%fvar=%f\n",mean,var);
	moment(data3,nreps,&mean,&adev,&sdev,&var,&skew,&curt);
	printf("Test of K-1 routine: mean=%fvar=%f\n",mean,var);
	/*return;*/
#endif

	count=1;
	if (Yule_flag)
	    {
/*print_mem_dbg(__FILE__,__LINE__);*/
	    if (!silentFlag)printf("Diversification simulation:\nseed = %li\nnreps=%i\nntaxa=%i\nspec rate = %f\nextinct rate=%f\n",
		    theSeed,nreps,ntaxa,speciation,extinction);
	    switch (diversemodel)
			{
			case 1:
		    	  if (!silentFlag)printf("MODEL=Forward Yule model\n");
			  if (stemFlag)
				printf("(Stem group simulation:N0=1)\n");
			  else
				printf("(Crown group simulation:N0=2)\n");
			  break;
			case 2:
		    	  if (!silentFlag)printf("MODEL=Backward birth-death model\n");
		    	  if (!silentFlag)printf("(Root node time normalized to one)\n");
			  break;
			}
	  /*printf("Predicted estimator variance:K1(infinite)=%f K2 (k=1) = %f\n",
	     Kendall_var1, Kendall_var2);*/
	    for (i=1;i<=nreps;i++)
		{
		if (!silentFlag)
			printf ("...generating replicate tree number %i\n",i);
		/*root=BDTree(ntaxa,speciation, extinction,0.1);*/
		/*root = BDTreeForward(T,speciation, extinction,0.1);*/
		    {
		    switch (diversemodel)
			{
			case 1:
		    	  root = YuleTreeForward(T, speciation, &NN1, &B,stemFlag);
			  break;
			case 2: /* this is bdback without normalizing root to 1 */
			 /* root=BDTree(ntaxa,speciation, extinction,0.1);*/
			  root=BDback(ntaxa,speciation,0);
			 /* data1[i]=treeDurLength(root);*/
			  data1[i]=treeAgeSum(root)/numIntNodes(root);
			  printf("Duration of tree = %f age=%f\n", data1[i], root->time);
			  break;
			case 5: /* this is bdback normalizing root to 1 */
			 /* root=BDTree(ntaxa,speciation, extinction,0.1);*/
			  root=BDback(ntaxa,1.0,1); /* set speciation to 1; doesn't really matter anyway given the renormalization */
			 /* data1[i]=treeDurLength(root);*/
			  data1[i]=treeAgeSum(root)/numIntNodes(root);
			  printf("Duration of tree = %f age=%f\n", data1[i], root->time);
			  break;
			case 3:
			    root=Yule_C(ntaxa, speciation);
			    break;
			case 4:	
			    root=BDTreeForward(T, speciation, extinction,0.0);
			    break;
			case 6:	
			    root=RY_1997(ntaxa, T, speciation,extinction,sampling_fraction); // Rannala Yang 1997 model
			    break;
			case 7:	
			    root=SisterGroupYule(T, speciation, speciation2, &NN1, &NN2, &B);
			    break;


			}
		 /*   data1[i]=NN1;
		    data2[i]=(NN1-2)/B;
		    data3[i]=B;*/
		    if(save_flag) /* now do this by default! */
			doSaveTree(root);
		    else
			DisposeTree(root);
		    }
/*****/
#if SIMLOOP

/*********** !!!!!!!!!! note that the following doObjFunc calls have invalid arg lists !!!!!!! **/

		gIndex=0;
		tree2aTimeArray(root,time1);

		gnpexp=gNexDataPtr->RateBlockParms.npexp;  /* KLUDGE */
		(void)doObjFunc(objLangFitch,root,"Simulated",LaF,&success1);
		gIndex=0;
		tree2aTimeArray(root,time2);


		(void)doObjFunc(objNP,root,"Simulated",NP,&success2);
		gIndex=0;
		tree2aTimeArray(root,time3);



		/*for (k=0;k<ntaxa-2;k++)
		    printf("%f %f %f\n",time1[k],time2[k],time3[k]);*/
		if (success1 && success2) /* store rep results only if both opts work*/
			{
			data1[count]=euclid_distance(time1,time2,ntaxa-2);
			data2[count]=euclid_distance(time1,time3,ntaxa-2);
			chiSqArray[count]=chiSq; /* a global */
			if (verbose)
			    printf("%f %f\n",data1[count],data2[count]);
		/*if (data2[count]>0.01)
		  {
		  printf("$$$%f\n", data2[count]);
		  for (k=0;k<ntaxa-2;k++)
		    printf("%f %f %f\n",time1[k],time2[k],time3[k]);
		  gNexDataPtr->RateBlockParms.verbose=2;
		  (void)doObjFunc(objNP,root,"Simulated",NP,&success2);
		  DrawTree(root, 1, 0);
		  DrawTree(root, 2, 0);
		  DrawTree(root, 4, 0);
		  printtree(root);
		  make_parens(root, 0);
		  exit(1);
		  }*/
			++count;
			}
		else
			doGenericAlert("WARNING: LF or NP failed\n");

		DisposeTree(root); /*******!!  !!********/
#endif
		} /* end nreps */
	       } /* endif */
		/* BDDiversity(ntaxa,speciation, extinction,0.1,interval); */

	--count;
/***************************

	Do the character evolution simulation


	For the normally-distributed model, the parameters start_rate and change_rate correspond to
	the mean and standard deviation of the normal respectively. The min and max values are still respected.


****************************/

		if(CharEvol_flag)
		  {
		  if (gNexDataPtr->isTrees)
		    {
			if (!silentFlag)
				{
				printf("\nBranch evolution simulation:\nseed=%li\n\nrate transition=%f\n",
					theSeed, start_rate,change_rate,min_rate,max_rate,rate_transition);
				printf("Gradual rate change=%i\n",gradual_flag);
				printf("Infinite=%i\n",infinite_flag);
				if(ratemodel==1)
					{
					printf("RATE MODEL:Normally distributed\n");
					printf("with parameters: mean=%f, sdev=%f, minrate=%f, maxrate=%f\n", 
						start_rate, change_rate, min_rate,max_rate);
					}
				else if (ratemodel==2)
					{
					printf("RATE MODEL:Autocorrelated\n");
					printf("with parameters: startrate=%f change rate=%f minrate=%f maxrate=%f\n", 
					start_rate,change_rate,min_rate,max_rate);
					printf("transition probability=%f change amount=%f\n", rate_transition, change_rate);
					}
				}
			lnode=gNexDataPtr->inTrees;
			LISTLOOP (lnode)
				{
				i=1;
				thisTree=lnode->item;
				for (j=1;j<=nrepsPerTree;j++)
					{
//					set_branch_rates(thisTree->root,start_rate,change_rate, min_rate,max_rate,rate_transition,gradual_flag, ratemodel);
					saveTree=copyTree(thisTree->root); // makes a deep copy because stuff will be overwr
					for (jj=1;jj<=nrepsPerBrRate;jj++)
					  {
//					  thisTree->root=copyTree(saveTree);
//					  set_branch_lengths(thisTree->root,infinite_flag);
//				    DrawTree(thisTree->root,1, 0);
//					printtree(thisTree->root);	
					  if (!silentFlag)
						{
						printf("\n ** Tree %i (Rate Replicate %i, Branch Length Rep %i)\n", i,j,jj);
						printf("tree SIMTREE = ");
							make_parens(thisTree->root, 0); /* TD of phylogram*/
						printf(";\n");
						}

#if 1
					// note if you estimate times, this will bollocks up the set rates, branches above..
trueAge=50.0;
rangeFactor=0.05; // increase to make the range of min and max branch lengths larger
ksteps=20;
gNexDataPtr->RateBlockParms.num_time_guesses=2;
for (ib=-ksteps/2;ib<=+ksteps/2;ib++)
	for (ic=-ksteps/2;ic<=+ksteps/2;ic++)
		{
lb=floor(pow(10, log10(5000)+ib*rangeFactor)); // make sure these are ints because of that lousy rounding problem with gradients
lc=floor(pow(10, log10(5000)+ic*rangeFactor));
printf("\n--%6.1f %6.1f\n\n",lb,lc);
thisTree->root=copyTree(saveTree);
(thisTree->root)->free=0;
(thisTree->root)->time=T;
(thisTree->root)->nodeIsConstrainedMax=0;
(thisTree->root)->nodeIsConstrainedMin=0; 
root=thisTree->root;
nodea=find_taxon_name(root,"a");
nodeb=find_taxon_name(root,"b");
nodec=find_taxon_name(root,"c");
nodeInt=nodec->anc;
nodea->length=10000;
nodeInt->length=5000;
		nodeb->length=lb;
		nodec->length=lc;
		printtree(root);
		DrawTree(root,1, 0);

					  doObjFunc(thisTree,LaF,1,TN,&success1);
					  print_ages(thisTree->root, 1.0,1.0,thisTree->method); 
		T_LF[ib+ksteps/2][ic+ksteps/2]=nodeInt->time;
DisposeTree(thisTree->root);

/*
thisTree->root=copyTree(saveTree);
(thisTree->root)->free=0;
(thisTree->root)->time=T;
(thisTree->root)->nodeIsConstrainedMax=0;
(thisTree->root)->nodeIsConstrainedMin=0; 
root=thisTree->root;
nodea=find_taxon_name(root,"a");
nodeb=find_taxon_name(root,"b");
nodec=find_taxon_name(root,"c");
nodeInt=nodec->anc;
nodea->length=10000;
nodeInt->length=5000;
		nodeb->length=lb;
		nodec->length=lc;
//					  bestSmooth=doCrossV(thisTree,LaF,1,TN,1.0,0.5,1,0);
					  bestSmooth=doCrossV(thisTree,PENLIKE,1,TN,0,0.5,8,0);
//					  gNexDataPtr->RateBlockParms.smoothing=bestSmooth;
DisposeTree(thisTree->root);
*/

thisTree->root=copyTree(saveTree);
(thisTree->root)->free=0;
(thisTree->root)->time=T;
(thisTree->root)->nodeIsConstrainedMax=0;
(thisTree->root)->nodeIsConstrainedMin=0; 
root=thisTree->root;
nodea=find_taxon_name(root,"a");
nodeb=find_taxon_name(root,"b");
nodec=find_taxon_name(root,"c");
nodeInt=nodec->anc;
nodea->length=10000;
nodeInt->length=5000;
		nodeb->length=lb;
		nodec->length=lc;
					  gNexDataPtr->RateBlockParms.smoothing=0.0001;
					  doObjFunc(thisTree,PENLIKE,1,TN,&success1);
		T_PL[ib+ksteps/2][ic+ksteps/2]=nodeInt->time;
					  print_ages(thisTree->root, 1.0,1.0,thisTree->method); 
DisposeTree(thisTree->root);

		}
//					  DisposeTree(thisTree->root);

for (ib=0;ib<=ksteps;ib++)
	{
	for (ic=0;ic<=+ksteps;ic++)
		printf("%6.1f\t",fabs(trueAge-T_LF[ib][ic]));
	printf("\n");
	}
printf("\n");
for (ib=0;ib<=ksteps;ib++)
	{
	for (ic=0;ic<=+ksteps;ic++)
		printf("%6.1f\t",fabs(trueAge-T_PL[ib][ic]));
	printf("\n");
	}
printf("\n");

printf("\t");
for (ic=0;ic<=ksteps;ic++)
	{
	lc=pow(10, log10(5000)+(ic-ksteps/2)*rangeFactor);
	printf("%6.1f\t",lc);
	}
printf("\n");

for (ib=0;ib<=ksteps;ib++)
	{
	lb=pow(10, log10(5000)+(ib-ksteps/2)*rangeFactor);
	printf("%6.1f\t",lb);
	for (ic=0;ic<=+ksteps;ic++)
		{
		printf("%6.1f\t",fabs(trueAge-T_LF[ib][ic])-fabs(trueAge-T_PL[ib][ic]));
		}
	printf("\n");
	}
printf("\n");
for (ib=0;ib<=ksteps;ib++)
	{
	for (ic=0;ic<=+ksteps;ic++)
		{
		printf("%6.1f\t",fabs(trueAge-T_LF[ib][ic])-fabs(trueAge-T_PL[ib][ic]));
		}
	printf("\n");
	}
printf("\n");
for (ib=0;ib<=ksteps;ib++)
	{
	for (ic=0;ic<=ib;ic++)
		if (fabs(trueAge-T_PL[ib][ic]) < fabs(trueAge-T_LF[ib][ic]))
			printf("1");
		else
			printf("0");
	printf("\n");
	}
printf("\n");
#endif
					  }
					thisTree->root=saveTree;
					}
				++i;
				}
		    }
		  else
		    printf("No trees currently in memory\n");
		  }





#if 0
	moment(data1,count,&LFmean,&adev,&LFsdev,&var,&skew,&curt);
	if (nreps-count)
		printf("There were %i failed replicates\n",nreps-count);
	/*printf("distance sim to LF: mean=%f stdev=%f\n",mean,sdev);*/
	moment(data2,count,&NPmean,&adev,&NPsdev,&var,&skew,&curt);
	/*printf("distance sim to NP: mean=%f stdev=%f\n",mean,sdev);*/
	moment(chiSqArray,count,&chiSqmean,&adev,&sdev,&var,&skew,&curt);
	if (NPmean<LFmean)
		whichBetter=1; /* NP has lower mean error */
	else
		whichBetter=0;

	printf("%f\t%f\t%f\t%f (+-%f)\t%f (+-%f)\t%i\n",
		change_rate, rate_transition,chiSqmean,LFmean,LFsdev/sqrt(nreps), 
		NPmean,NPsdev/sqrt(nreps), whichBetter);
#endif
#if SIM_LOOP
	} 
#endif


/***** RANDOM BRANCH SAMPLING SIMULATION ****/

//  NB. The nested taxon model always picks nodes without replacement...see nextRndNode()

/*****/

#if 0
		// Allocation and initialization

		    markedNodes=(NODETYPE**)myMalloc((nMark+1)*sizeof(NODETYPE*)); /* 1-offset array */
		    lnode=gNexDataPtr->inTrees;
			/* following lines assume that all trees in list have some number of nodes! */
			/* And we're going to add all the histogram entries together across trees */
		    thisTree=lnode->item;
		    nNodes=numNodes(thisTree->root);
		    nTaxa=numdesc(thisTree->root);
			 /* Following are all (1-off) histos: # one bin for each possible taxon size */
		    MaxGroupSize=nTaxa; /* Following arrays need to range from a group size of 0 to nTaxa*/
		    histo=(long *)myMalloc((1+MaxGroupSize)*sizeof(long));
		    histoB=(long *)myMalloc((1+MaxGroupSize)*sizeof(long));
		    histoMonophyletic=(long *)myMalloc((1+MaxGroupSize)*sizeof(long)); 
		    histo2=(long *)myMalloc((1+MaxGroupSize)*sizeof(long)); 
		    histoTotal=(long *)myMalloc((1+MaxGroupSize)*sizeof(long));
		    histo2Total=(long *)myMalloc((1+MaxGroupSize)*sizeof(long)); 

			TotalReps=nreps*nrepsPerTree;
		    RepCount=(double *)myMalloc((TotalReps+1)*sizeof(double)); // 1-offset for moment function
		    RepMean=(double *)myMalloc((TotalReps+1)*sizeof(double)); // 1-offset for moment function
		    RepDominance=(double *)myMalloc((TotalReps+1)*sizeof(double));
		    RepFreq1Class=(double *)myMalloc((TotalReps+1)*sizeof(double));
		    RepFractMonophyletic=(double *)myMalloc((TotalReps+1)*sizeof(double));
		    RepMonophyleticSpecies=(double *)myMalloc((TotalReps+1)*sizeof(double));

		    irepcount=0;
			for (i=0;i<=nTaxa;i++)
							{
							histoTotal[i]=0;
							histo2Total[i]=0;
							histoB[i]=0;
							}

			if (exclusive)	// for exclusive models we never want to sample with replacement!
				withReplace=0;
				
		    LISTLOOP (lnode)
			    	{
			    	thisTree=lnode->item;
					nodeArray = newAllNodeArray(thisTree);
					for (j=1;j<=nrepsPerTree;j++)
					   {
			    		++irepcount;
			    		printf("...working on replicate %i\n",irepcount);
						for (i=0;i<=MaxGroupSize;i++)
							{
							histo[i]=0;
							histoMonophyletic[i]=0;
							histo2[i]=0;
							}
					   if (exclusive) // under exclusive model, the root node is always part of a taxon; mark it in 
					   				  // 0-th element of this array
					   		markedNodes[0]=thisTree->root;
					   if (rndBranchDur) // the other nodes marked are in elements 1..N	
							RandomBranches(thisTree,nNodes,nodeArray,nMark,markedNodes,withReplace); /* NB!  If we have a very long branch, the rand chars will hit it often, and the naive w/o replacement algorithm will thrash */
					   else
							markRandomNodes(thisTree,nMark,markedNodes);  /* NB! Change the dynamic allocation in this routine to static...*/
								
					   countTaxa=0;
					   countMonotypes=0;
					   countMonophyletic=0;
					   countExc=0;
					   countMonophyleticSpecies=0;
/*
	Under the exclusive model there will be N+1 groups recognized, where N is the number of marks;
	Under the nested model there will be N groups recognized
*/
					   for (i=0;i<=nMark;i++)
								{
								if (i==0 && !exclusive) // skip the root node (in ..[0]) if nested model
									continue;
								++countTaxa;
								node=markedNodes[i];
								unMarkNode(node);/* have to unmark this node for next function to work right */
								size = numUnMarkedDesc(node); // NB! sometimes you can get a size of zero! Then this is an "orphan"	
								size2= numdesc(node); // size of this taxon assuming it's monophyletic
								++histo[size]; // these are the group sizes possibly paraphyletic (for exclusive sampling)
								if (size<=9) sizeB=0;
								if (size>9 && size<=99) sizeB=1;
								if (size>99 && size<=999) sizeB=2;
								if (size>999) sizeB=3;
								/*sizeB=(long)floor((size-1)/10.0)+1;*/
								++histoB[sizeB];
								++histoTotal[size]; // keep track across reps
								countExc+=size;
								++histo2[size2]; // these are the monophyletic group sizes (for clade sampling)
								++histo2Total[size2]; // keep track across reps
								if (size2==1)
									++countMonotypes;
								if (size>1 && size==size2) /* don't count the "monophyletic" monotypes*/
									{
									++countMonophyletic;
									++histoMonophyletic[size];
									countMonophyleticSpecies+=size; // num of species in monophyletic higher taxa
									}
								markNode(node);
//if (size2>10) printf("Mono:%i Exc:%i\n",size2,size);
								}


// At this point store the stats for each replicate (regardless of whether replicate across one tree or many)							
						if (exclusive)
							h=histo;
						else
							h=histo2;							
						histoStat(h, MaxGroupSize,nTaxa, &count, &mean, &freq1class, &maxS, &dominance);
						RepCount[irepcount]=count; // number of taxa > 0 found 
						RepMean[irepcount]=mean;
						RepDominance[irepcount]=dominance;
						RepFreq1Class[irepcount]=freq1class;
						RepFractMonophyletic[irepcount] = (double)countMonophyletic/(countTaxa-countMonotypes);
						RepMonophyleticSpecies[irepcount]= (double)countMonophyleticSpecies/(nTaxa-countMonotypes);
printf("Monotypes and max size:%i %i %i\n",irepcount,histo[1],maxS);
// printf("%i %f %f %f\n",irepcount,RepMean[irepcount],RepDominance[irepcount],RepFreq1Class[irepcount]);
					   } // end nrepspertree
				myFree(nodeArray);
				} // end listloop

		// Print results

		    if (!silentFlag)
			{
		    	printf("******************************************\n\n");
		    	printf("Results of simulation of random taxonomies\n\n");
		    	if (exclusive) printf ("Taxonomy consists of exclusive groups\n"); 
		    	else printf ("Taxonomy consists of nested groups\n");
		    	if (rndBranchDur) printf ("Branches are sampled according to durations\n"); 
		    	else printf ("Branches are sampled equally\n");
		    	if (withReplace) printf("Branches are sampled with replacement\n");
		    	else printf ("Branches are sampled without replacement\n");

		    	printf("Average histogram for model across %li replicate tree simulations\n",nreps);
		    	for (i=1;i<=MaxGroupSize;i++)
			    if (exclusive)
				{
				if (histoTotal[i]>0)
					printf("%li\t%f\n",i,histoTotal[i]/(float)TotalReps);
				}
			    else
				if (histo2Total[i]>0)
					printf("%li\t%f\n",i,histo2Total[i]/(float)TotalReps);

			}


		// Make a plot of the sum total histogram across all replicates

		    if (exclusive)
			h=histoTotal;
		    else
			h=histo2Total;

		    histoStat(h, MaxGroupSize,nTaxa, &count, &mean, &freq1class, &maxS, &dominance); // just to get count
		    X = (double *)myMalloc(count*sizeof(double));
		    Y = (double *)myMalloc(count*sizeof(double));
		    count=0;
		    for (i=1;i<=MaxGroupSize;i++)
				if (h[i]>0) 
					{
					X[count]=log10(i);
					Y[count]=log10(h[i]/(float)TotalReps);
					++count;
					}
		    for (i=0;i<count;i++)
			printf ("%f\t%f\n",X[i],Y[i]);
		    if (!silentFlag)			
		    	dumbPlot(X,Y,count);

			moment(RepCount,TotalReps,&mean,&adev,&sdev,&var,&skew,&curt);
			printf("Number of higher units generated across reps: mean=%f var=%f\n",mean,var);
			moment(RepMean,TotalReps,&mean,&adev,&sdev,&var,&skew,&curt);
			printf("Taxon size across reps: mean=%f var=%f\n",mean,var);
			moment(RepDominance,TotalReps,&mean,&adev,&sdev,&var,&skew,&curt);
			printf("Dominance across reps: mean=%f var=%f\n",mean,var);
			moment(RepFreq1Class,TotalReps,&mean,&adev,&sdev,&var,&skew,&curt);
			printf("Size 1 class (monotypes) across reps: mean=%f var=%f\n",mean,var);
			moment(RepFractMonophyletic,TotalReps,&mean,&adev,&sdev,&var,&skew,&curt);
			printf("Fraction of monophyletic nonmonotypes: mean=%f var=%f\n",mean,var);
			moment(RepMonophyleticSpecies,TotalReps,&mean,&adev,&sdev,&var,&skew,&curt);
			printf("Fraction of all species found in nonmonotypic monophyletic higher taxa: mean=%f var=%f\n",mean,var);
		    	for (i=0;i<=3;i++)
				printf("%li\t%f\n",i,histoB[i]/(float)TotalReps);

			myFree(X);
			myFree(Y);
			myFree(histo);
			myFree(histo2);
			myFree(histoTotal);
			myFree(histo2Total);
			myFree(histoMonophyletic);
			myFree(markedNodes);
			myFree(RepMean);
			myFree(RepDominance);
			myFree(RepFreq1Class);
#endif

/*****/
/*****/



/*  put this as a command such as 'statistics'
		moment(data1,nreps,&mean,&adev,&sdev,&var,&skew,&curt);
		printf("Test of Mean Duration in BDBack: mean=%f var=%f\n",mean,var);
		moment(data2,nreps,&mean,&adev,&sdev,&var,&skew,&curt);
		printf("Test of K-infinite estimator in YuleTreeForward: mean=%fvar=%f\n",mean,var);
		moment(data3,nreps,&mean,&adev,&sdev,&var,&skew,&curt);
		printf("Test of mean-B in YuleTreeForward: mean=%fvar=%f\n",mean,var);
*/
	myFree(data1);
	myFree(data2);
	myFree(chiSqArray);
	myFree(time1);
	myFree(time2);
	myFree(time3);
	return;						
}


/****************************************************************/
static void histoStat(long h[], long N, long nTaxa, long *count, double *mean, double *freq1class, long *maxS, double *dominance)
{
// NB! The histo array, h, is of size N+1 and is zero-offset, because the min size class is zero for our stuff
// Thus, note that it should be called by ...(h,N,...), rather than (...,N+1).

// Note that h[0] contains the count of those taxa orphaned with no terminals in them at all. 
// The number of actual taxa observed in the classification is therefore slightly less than the number of randomly selected nodes for marking


long i,s2=0;
*maxS = -1;
*count=0;
for (i=1;i<=N;i++)
	{
	
	if (h[i]>0) 
		{
		*maxS = i;  // this is the maximum size observed in the histogram
		(*count)+=h[i];	// total number of points (higher taxa)
		s2+=h[i]*i;
		}
	}
	
*mean=(double)s2 /(*count);
*freq1class = (double) h[1]/(*count);
*dominance = (double)(*maxS)/nTaxa;

return;
}

/****************************************************************/

static void doSaveTree(NODETYPE *root)  /* adds a simulated tree to the tree list */
{
extern int curTree;
int flag=0, numtips, numinternals, roottomy;
char *stemp, *tree_name="SIMTREE", *TD="";
PtrList aTreeList;
TREE aTree;


	if (gNexDataPtr->inTrees == NULL)  /* if this is the first tree */
	    {
	    gNexDataPtr->inTrees=pNewListAlt(sizeof(struct treetype));
	    aTreeList=gNexDataPtr->inTrees;
	    }
	else				/* if a later tree */
	    aTreeList=pListAddNode(gNexDataPtr->inTrees,sizeof(struct treetype));
		
	(void)appendStrList(gNexDataPtr->TDLabelList,tree_name);
	(void)appendStrList(gNexDataPtr->TDList,TD);			
	++curTree;
	gNexDataPtr->NumTrees=curTree;
	gNexDataPtr->isTrees=1;
	aTree=aTreeList->item;
	if (root)
	    {
	    init_node_ids(root, 0); 
	    init_free(root); /* sets default to estimate all internal nodes but root */ 
	    aTree->root=root;
	    aTree->name=DupStr(tree_name);
	    aTree->TD=DupStr(TD);
	    TreeStats(root, &numtips, &numinternals, &roottomy);
	    aTree->numTaxa=numtips;
	    aTree->numBranches=numBranches(root);
	    aTree->basalTomy=roottomy;
	    }

return;						
}
/****************************************************************/
void doWuLiCommand(NODETYPE * root)
{
	int id[3], ix,jx;
	char c, *dummy;
	

	if ((gNexDataPtr->isChars==0) || ( gNexDataPtr->isTaxa==0)
		|| gTaxaSet == NULL)
		return;	/* don't have the right data from NEXUS file, so bail */

	while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		if (parse_assignment2("RRTYPE"))
			{
			if (isEqual(LocalToken,"WULI"))
				gNexDataPtr->RateBlockParms.RRtype = WULI;
			if (isEqual(LocalToken,"STEEL"))
				gNexDataPtr->RateBlockParms.RRtype = STEEL;
			if (isEqual(LocalToken,"TAJIMA"))
				gNexDataPtr->RateBlockParms.RRtype = TAJIMA;
			if (isEqual(LocalToken,"MIKE"))
				gNexDataPtr->RateBlockParms.RRtype = MIKE;
			}
		if (parse_assignment2("BS"))
		    {
		    if (isEqual(LocalToken,"YES"))
			gNexDataPtr->RateBlockParms.isBS=1;
		    else
			gNexDataPtr->RateBlockParms.isBS=0;

		    }
		if (parse_assignment2("NREPS"))
			gNexDataPtr->RateBlockParms.NReps=strtod(LocalToken,&dummy);
		if (parse_assignment2("SEED"))
			{
			gNexDataPtr->RateBlockParms.seed=strtod(LocalToken,&dummy);
			srand(gNexDataPtr->RateBlockParms.seed);
			}
		}

	doRelativeRates(gTaxaSet,root);
	
	return;
}
/****************************************************************/
void doTaxaSetCommand(void)
{
	int id[3], ix,jx;
	char c, *dummy;

	if ( gNexDataPtr->isTaxa==0)
		return;	/* don't have the right data from NEXUS file, so bail */

	if(gTaxaSet)  /* get rid of old list if it's present */
		    {
		    freeStrList(gTaxaSet);
		    }

	gTaxaSet=newStrList();
	while (!isEqual(aTokenPtr=nextToken(),";"))	
			appendStrList(gTaxaSet,aTokenPtr); /* store the label */
	
	printf("Using the following taxa:\n");xprintStrList(gTaxaSet);
	return;
}

/****************************************************************/
void doExSets(void)

/* Sets up an exclusion set array  in which a zero means excluded and 1 means included.
Format is 'exsets n1 n2 n3 n4 - n5 n6;'  NOTE THAT THERE MUST BE A SPACE BEFORE AND AFTER
THE DASH--THIS IS A NON NEXUS COMPLIANT WORKAROUND, but the NEXUS format does not recognize
dashes as punctuation since they can also represent gaps (?) FIX Later
Each time an exclusion set is invoked, the array is reset to match command */
{
	char* dummy;
	long icur,ilast,ix;
	int *excArray;
	excArray=gNexDataPtr->excArray;
	for (ix=0;ix<gNexDataPtr->NChars;ix++)
		excArray[ix]=1;	/* initializes exclusion set array for use in block*/
/*	fprintf(fpOut2,"[!NOTE: Some sites excluded in following analyses]\n");*/
	while (!isEqual(aTokenPtr=nextToken(),";"))	/* if its not a ';' it should be a number*/
		{
	
		if (  isdigit(*aTokenPtr)  )
			{
				icur=strtod(aTokenPtr,&dummy);
				ilast=icur;
				excArray[icur-1]=0;	/* this is a zero offset array */
			}
		else
			if (isEqual(aTokenPtr,"-"))
				{
				aTokenPtr=nextToken();
				icur=strtod(aTokenPtr,&dummy);
				for (ix=ilast;ix<=icur;ix++)
					excArray[ix-1]=0;	/* this is a zero offset array */
					
				}
			
		

		}
			
					
return;
}
/****************************************************************/
void doSitesCommand(int what) /* include or exclude positions */
{
long ix;

switch (what) 
	{
	case 0:
		for (ix=0;ix<gNexDataPtr->NChars;ix++)
			gNexDataPtr->excArray[ix]=1;
		printf("\n\n*** All sites included from now on ***\n\n\n");
		break;
	case 1:
		for (ix=0;ix<gNexDataPtr->NChars;ix++)
			if ( (ix+1)/3 != (ix+1)/3.0)
				gNexDataPtr->excArray[ix]=0;
		printf("\n\n*** First and second positions excluded from now on ***\n\n\n");
		break;
	case 3:
		for (ix=0;ix<gNexDataPtr->NChars;ix++)
			if ( (ix+1)/3 == (ix+1)/3.0)
				gNexDataPtr->excArray[ix]=0;
		printf("\n\n*** Third positions excluded from now on ***\n\n\n");
		break;

	}
return;
}

/****************************************************************/
/****************  MISCELLANEOUS FUNCTIONS **********************/
/****************************************************************/

void freeNexusStructure(struct NexDataType *nex)
{
freeStrList(nex->TaxaList);
freeStrList(nex->TDList);
freeStrList(nex->TDLabelList);
if (nex->isChars)
	freeStrList(nex->DMList);	/* this won't be allocated if no characters */
freeStrList(nex->TransList);
myFree(nex);


return;
}
/****************************************************************/

void doError(char* p[], int which)
{
doGenericAlert(p[which]);
}
/****************************************************************/
void TreeSummary(int whichTree)
{
	NODETYPE *root;
	char * TreeName, *TD;
	int numTips,numInternals, roottomy;
	TreeName=getkthStr(gNexDataPtr->TDLabelList,whichTree);
	TD=getkthStr(gNexDataPtr->TDList,whichTree);
	root=string_to_tree(TD);
	if (root != NULL)
			{
			TreeStats(root,&numTips,&numInternals, &roottomy);
			DisposeTree(root);
			}
	printf("Processing tree %i (%s) (taxa=%i; No. internal nodes = %i; Basal tomy=%i)\n",
				whichTree, TreeName,numTips,numInternals,roottomy); 

	return;

}
/****************************************************************/
int parse_assignment(char * target,char ** token)

/* on entry 'aTokenPtr' contains the putative first word of a three token
assignment statement of the form 'word1=word2'.  This function checks to see
if word1 is the same as 'target' and if so, it returns the address of a string
containing 'word2' or NULL if an error occurs.  aTokenPtr is
set to the last token in the assignment statement
If no match, aTokenPtr is left unchanged!! */

/*** BAD CODE *** causes memory leaks, probably failing to 
		free LocalTokens */

{
		if (isEqual(aTokenPtr,target))
			{
			aTokenPtr=nextToken();
			/*if (aTokenPtr==NULL) return 0;*/
			if (!isEqual(aTokenPtr,"="))
				{
				printf("Bad assignment statement=(%s)\n",aTokenPtr);
				fatal("exiting...");
				}
			aTokenPtr=nextToken();
			*token = DupStr(aTokenPtr);
			return 1;
			}
	return 0;
}

/****************************************************************/
int parse_assignment2(char * target)

/* on entry 'aTokenPtr' contains the putative first word of a three token
assignment statement of the form 'word1=word2'.  This function checks to see
if word1 is the same as 'target' and if so, it returns the address of a string
containing 'word2' or NULL if an error occurs.  aTokenPtr is
set to the last token in the assignment statement
If no match, aTokenPtr is left unchanged!! */

{
		if (isEqual(aTokenPtr,target))
			{
			aTokenPtr=nextToken();
			if (!isEqual(aTokenPtr,"="))
				{
				printf("Bad assignment statement=(%s)\n",aTokenPtr);
				fatal("exiting...");
				}
			aTokenPtr=nextToken();
			if (strlen(aTokenPtr)< MAX_LOCAL_TOKEN_SIZE -1)
				strcpy(LocalToken,aTokenPtr);
			else
				fatal("local token size exceeded\n");
			return 1;
			}
	return 0;
}

/****************************************************************/


void checkMatrix(void)
{
int itaxa;
char* c;
for (itaxa=0;itaxa<gNexDataPtr->NTaxa;itaxa++)
	{
	c=getkthStr(gNexDataPtr->DMList,(long)(itaxa+1));
	while(*c++)
		if ( (*c=='{') || (*c=='}'))
			{
			doGenericAlert("Polymorphism not allowed: Do not invoke rate tests");
			return;
			}
	}
	
return;	
}
/***************/
#if 0
int gNComp;
static void doCrossV(PtrList TreeList, int method,double EstMult,double PrdMult,double cvStart,double cvInc,double cvNum)

/*  does a cross validation with the range of the tuning parameter set to run from
	[cvStart, cvStart+cvInc, ...,cvStart+cvInc*(cvNum-1)]

NB! Haven't added variable 'algorithm' here yet

*/
{
#define isEven(k) ((k)/2 == ((k)/2.0))
TREE tree1,tree2;
PtrList p;
int i,j,success,nTrees,collFlag=0;
double * cvScore,cvSum;

double smooth;

nTrees=pLengthList(TreeList);
if (!isEven(nTrees))
	{
	doGenericAlert ("Must have even number of trees to do cross validation");
	return;
	}
else
     {
     gNComp=nTrees/2;
     cvScore = vector(1,gNComp);
     for (j=0;j<cvNum;j++)
	{
        cvSum=0.0;
	smooth=pow(10.0,j*cvInc+cvStart);
	gNexDataPtr->RateBlockParms.smoothing=smooth;
	p=TreeList;
	for (i=1;i<=gNComp;i++)
		{
		tree1=p->item;
		tree2=p->next->item;
		if (j==0)		/* only do this on the first pass; otherwise we keep multiplying the lengths */
			{
			if (EstMult != 1.0)traverseMultiplyLength(tree1->root, EstMult);
			if (PrdMult != 1.0)traverseMultiplyLength(tree2->root, PrdMult);
			}
		if (collapseLengthsTree2Tree(tree1,tree2))
			collFlag=1;				/* set this if some tree had to collapse a branch */
		copyLengthsTree2Tree(tree1->root,tree2->root); /* put lengths from tree2 onto field nodeReal of tree1 */
		doObjFunc(tree1,method,0,&success);
		cvScore[i]=cvSquareError(tree1,method);
		cvSum+=cvScore[i];
		p=p->next->next;
		} 	
	printf("\nCross Validation Analysis\n\n");
	for (i=1;i<=gNComp;i++)
		printf ("Cross Validation Score [%2i] = %f\n",i,cvScore[i]);
	printf("Cross Validation Score Total:smoothing = %f CV=%f\n",smooth,cvSum/gNComp);
	if (collFlag)
		printf("NOTE: Some partitions had 0-length branches that had to be collapsed to estimate cv scores\n");
	}
    }
}

#endif
