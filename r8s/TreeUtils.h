#ifndef _TREEUTILS
#define _TREEUTILS

#define LARGE_NODES	0		/* set to 1 iff we want to use all fields in every node,
					if this is 0, we cannot do HMM now! */
#define FLT_MAX 1e35	/* no longer defined in limits.h--do it here temporarily */
#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "structures.h"
#define LF 10
#define RETURN 13
#define COLON ':'
#define BAR	'|'
#define PLUS '+'
#define DASH '-'
#define	SPACE	' '
#define COMMA	','
#define	RIGHTPARENS	')'
#define	LEFTPARENS	'('
#define MAXSTRING	5000		/* maximum length of string (INCREASE LATER) */
#define IsItAlphaNum(c)	  ( (c) >= 48 && (c)<=57 ) || ( (c) >= 65 && (c)<=90 ) || ( (c) >= 97 && (c)<=122 ) 
#define min(a,b)			( (a)<=(b) )  ? (a):(b)
#define max(a,b)			( (a)>=(b) )  ? (a):(b)
#define SIBLOOP(c)			for (; (c); (c)=(c)->sib)
#define isTip(c)			 ( (c)->firstdesc == NULL )
#define isConstrained(c)		( (c)->nodeIsConstrainedMax || (c)->nodeIsConstrainedMin)
#define isConstrainedMax(c)		( (c)->nodeIsConstrainedMax)
#define isConstrainedMin(c)		( (c)->nodeIsConstrainedMin)
#define isFree(c)			( (c)->free == 1)
#define isFixed(c)			( (c)->free == 0)
#define isRoot(c)			( (c)->anc == NULL )
#define FLAGFLIP(c)			( (c) = (c) ^ 1)  /* use XOR operator */
#define MIN(a,b) (((a)<(b))?(a):(b))
#define LN2 0.69314718
#define LOG2(x) (log((double)(x))/LN2)  /*base 2 logarithm */

#define MAXCLADES 200 /* also defined in ReadNexusFile2 */
#define MAXTAX	  100

/* STRUCTURES AND PROTOTYPES */

/**********************************************/

/*  Node structure */

struct nodetype {
				struct nodetype 		*anc;
				struct nodetype 		*firstdesc;
				struct nodetype 		*sib;
				struct nodetype			*nodePtr;	/* generic pointer to some other node */
				char 				*taxon_name;
				double				length;		/* length of subtending branch */
				int 				order;
				long 				numdesc;
				int				numSelectedDesc;/* number of selected nodes 
									    below this one (including this one) */
				long				id;	
				int 				X;		/* positions on screen */
				int 				Y;
				double				time;		/* current time of node... */
				double				nodeReal;	/* Let's use this for various real numbers */
				short				nodeFlag;	/* for various flags */
				double				estRate;	/* estimated rate, usually for branch */
				double 				nodeEstRate;	/* estimated rate,special for node method */

				int				isQueryNode;	/* 1 if node to be used in query; USED IN NODE MARKING ROUTINES */
				
				int				isCompactNode;	/* 1 if this node is displayed
										as a clade of all its descendants*/
				int				toggleDesc;	/* 1 if all descendants should
										be queried */
				int				nodeIsConstrainedMax;
				int 				nodeIsConstrainedMin;
  				int				modelID;	/* takes integer values for different rate parms under 
										local clock model */
				int				free;	/* 1 if we estimate this node's time */
				double				minAge;		/* present = 0; root = 1 */
				double				maxAge;		/* ...These are constraints on ages*/
					
				double				cumulProb;	/* Used in RandomTree modules */
				double				like;		/* likelihood of subtending branch */
				double				chiSq;


				double				CL[5];		/* conditional likelihood for four states */
				double				CLmax[5];	/* Pupko max of 4 conditional likelihoods*/
				double				CLmarg[5];	/* Marginal likelihoods computed by rerooting */
				int				CLopt[5];	/* Pupko optimal state choice */
				int				opt;		/* optimal state at this node */
				char			state;		/* character state for discrete char algorithms */

#if LARGE_NODES
                                double beta[2];
                                double beta_sum;
                                double delta[2];
                                int psi[2];
#endif
				};
typedef 	struct nodetype NODETYPE;
typedef		struct nodetype * NODE;

/**********************************************/

/*  Tree structure */

struct treetype 	{
			char *		name;
			char *		TD;
			long		numTaxa;
			long		numBranches;
			long		basalTomy;
			NODETYPE 	*root;
			PtrList		cladeSet;
			double 		est_b;
			double		est_c;	/*estimates of gamma rate parms */
			double		estRate;
			int		timesAvailable;	/* 1 if times have been estimated */
			int		method;		/* method used to estimate times */
			double		obj;		/* value of obj func at soln */
			};
typedef 	struct treetype  * TREE;




/*************************************************/

void TreeToNodePtrList(NODETYPE *node,  PtrList NodeList)	;
double ** tree2VCV(TREE t, int i);
void print_named_ages(NODETYPE *node);
double pathLengthTime(NODE a, NODE b);
double pathLengthTimeModel(NODE a, NODE b, int model);
NODE  mrca(NODE a, NODE b );
TREE Subtree_Initialize(TREE T, NODETYPE *node);
long numUnMarkedDesc(NODETYPE *node);
NODETYPE * nextRndNode(long nNodes, NODETYPE ** nodeArray);
void markNode( NODETYPE  * n);
void unMarkNode( NODETYPE  * n);
int isNodeMarked( NODETYPE  *n);
void unMarkTree(TREE T);

void setLocalModel(NODETYPE *n,int model,int stemFlag);
void summarize_rates(TREE t);
static void recurse_summarize_rates(NODETYPE * n, long * ix, double r[]);

void scaleTree(NODETYPE * root, double calAge, NODETYPE * calNode);

void setupConstrainedNodeArray(NODETYPE * node, NODE nodeArray[], int *curIndex);
void setupFixedNodeArray(NODETYPE * node, NODE nodeArray[], int *curIndex);
int numConstrainedNodes(NODETYPE * node);
int numFixedNodes(NODETYPE * node);
double unFixNodeAge(NODETYPE *node);
double fixNodeAge(NODETYPE *node);
void print_rates(NODETYPE *n,int method);
void RemoveTaxonLvAnc(NODETYPE * n);
int collapseLengthsTree2Tree(TREE t1,TREE t2);
void setNodeEstRate(NODE node);
void zeroEstRate(NODETYPE *node);
double cvSquareError(TREE t, int method);
double cvSquareErrorBranch(TREE t, NODE n,int method,double *chiSq);


void copyLengthsTree2Tree(NODETYPE * node1,NODETYPE * node2);
double LFunconsT(NODETYPE *node);
double LFuncons1T(NODETYPE *node);
double LFuncons(NODETYPE *node);
double LFuncons1(NODETYPE *node);
void preOrderVoid(NODETYPE *node,void (*f)(NODETYPE *));
double preOrderArg(NODETYPE *node,double (*func)(NODE node, double farg),double farg);
double preOrder(NODETYPE *node,double (*f)(NODETYPE *));
void unSetConstraints(NODETYPE * node);
int maxAgePresent(NODETYPE * node);
int constraintsPresent(NODETYPE * node);

int tipsDifferentAges(NODETYPE *node);

void		ABCSuperTreePurvis(NODETYPE *node, StrListPtr UniqueList, 
		    char **dataMatrix,float *wtset)	;
void		ABCSuperTree(NODETYPE *node, StrListPtr UniqueList, 
		    char **dataMatrix,float *wtset);	
void		TreeToTaxaList(NODETYPE *node,  StrListPtr taxaList);	
void		TreeToTaxaPtrList(NODETYPE *node,  PtrList NodeList);	
int			 maxorder(NODETYPE *node);
int 			numdesc(NODETYPE *),
			stringcheck(char *td);
			

int numFreeNodes(NODETYPE *node);
NODETYPE	* sister(NODETYPE * n);
NODETYPE 	* ReRoot(NODETYPE * atNode);
NODETYPE 	* ReRoot2(NODETYPE * atNode);
void 		Flip(NODETYPE *a);
NODETYPE * 		RemoveTaxon(TREE t,NODETYPE * theChild);
void 		AddChild(NODETYPE * parent, NODETYPE * theChild);

void		 updateSubtrees(NODETYPE *srcNode);
NODETYPE	 *createSubtree(NODETYPE *srcNode, int SubtreeSize);
void 		copyNodeInfo(NODETYPE *source,NODETYPE *dest);
double 		*sort_node_time_array(NODETYPE *root);
double 		get_sum_durations(NODETYPE *node);
NODETYPE	*newnode(void),
		*makegroup(void),
		*string_to_tree(char *tree_description);
void 		collapse_node(NODETYPE *node);			
void 		collapse_zero(NODETYPE *node);
void 		Node_Destructor(NODETYPE *node);
void 		make_parens(NODETYPE *root, int flag);
void 		Tree_Destructor(TREE aTree);
void		DisposeTree(NODETYPE *node);
int 		numIntNodes(NODETYPE *node);
int 		numBranches(NODETYPE *node);
void 		TreeStats(NODETYPE *root, int * numtips, 
				int * numinternals, int * roottomy);
void		setNodeName(NODETYPE *node, char *name);
void 		printtree(NODETYPE *node);
NODETYPE * 	find_taxon_name(NODETYPE *node,char *name);
void 		Tree_Initialize(TREE aTree, char *TD, char *name);
void 		print_tree_list(PtrList treeList);
void		print_ages(NODETYPE *node, double time, double calAge,int method);
void		init_node_ids(NODETYPE *node, int gId);
void		convert_branchlength_to_time(NODETYPE *root);
NODETYPE	*find_id(NODETYPE *node,int id);
int		node_tomy(NODETYPE *node);
int		isNodeDescendant(NODETYPE *nodeA, 
		    NODETYPE *nodeB);
NODETYPE *	MRCA(NODETYPE *, StrListPtr taxaList);
void		traverseMultiplyLength(NODETYPE *, double x,int round);
void 		Tree2CladeSets(NODETYPE *node, StrListPtr allTaxaList, int nTaxa, 
		    PtrList SetList);	
void 		printCladeSets(PtrList SetList);
PtrList		Tree2CladeSet(TREE thisTree, StrListPtr allTaxaList);
int 		group_a_clade(NODETYPE *root, StrListPtr taxaList);
int 		any_zero_internal_branches(NODETYPE *node);
int             any_zero_terminal_branches(NODETYPE *node);
void 		printLikes(NODETYPE *node);
void 		printnodeLike(NODETYPE *node);
void 		ClusterHistogram(NODETYPE * node, long * array,long TSize);
double 		treeDurLength(NODETYPE * node);
double 		treeLength(NODETYPE * node);
double 		treeAgeSum(NODETYPE * node);
int 		compar(const void *v1,  const void *v2);
int 		numNodes(NODETYPE *node);
void 		init_free(NODETYPE *node);
void rootToTips(NODETYPE* node,double curLen);
NODETYPE*  copyTree(NODETYPE* a);
#endif
