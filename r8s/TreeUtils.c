/******************        Tree Utility module    ********************

					
	NOTE: The tree/node data structure is defined in TreeUtils.h
					
**********************************************************************/
#if 0
#define isEqual(a,b)		(!strcmp((a),(b)))
#endif

#include <limits.h>
#include "ObjFunc.h"
#include "TreeUtils.h"
#include "MyUtilities.h"
#include "NRCvectorUtils.h"
#include "memory.h"
#include "DrawTree.h"
#include "structures.h"
#include <stdlib.h>
#include <math.h>
#include "TimeAlgorithms.h"
#include "moment.h"
#include "DistrFuncs.h"

/* ...Global variables throughout module */

char 	*gStringptr;
int	gId, gFlag;
int	gCount;

double sqrarg;
#define SQR(a) (sqrarg=(a),sqrarg*sqrarg)
void Flip2(NODETYPE *a);
void printnode(NODETYPE *node);
static void setModel(NODE n,int model);
static void preOrderIntArg(NODETYPE *node,void (*func)(NODE node, int iarg),int iarg);
static double SThelper(NODETYPE * node,double factor);
static void collapse_zero_2trees(NODE node1, NODE node2);
static double cvCS(NODETYPE * node);
static double cvCST(NODETYPE * node);
static double cvSQNP(NODE node);
static double cvSELF(NODE n,double rate);
static double moveNER(NODE node);
static double setNER(NODE node);
static double cvSQET(NODETYPE * node);
static double zeroER(NODETYPE *node);
static void insertNode(NODETYPE *node,  NODETYPE* anc);
static void updateOneSubtree(NODETYPE *subRoot);
NODETYPE * prevSib(NODETYPE* node);

/****************************************************************/

// returns a 1-off variance-covariance matrix based on the ultrametric distances on tree t including only branches with model ID matching model
// note that the order of indices and terminals is determined by the call to 'TreeToTaxaPtrList'

double ** tree2VCV(TREE t, int model)
{
	PtrList lnode;
	long i, j;
	NODETYPE *root, *node;
	NODE a,b,c;
	PtrList nodeList;
	long lengthList,n;
	double ** vcv,T;
	root=t->root;
	nodeList=pNewList();
	TreeToTaxaPtrList(root,nodeList);
	n=pLengthList(nodeList); // the number of taxa!
	vcv = matrix(1,n,1,n);
printf("\nVariance-covariance matrix for model %i\n",model);
	for (i=1;i<=n;i++)
			{
			a=(NODE)(pListgetkthNode(nodeList, i)->item);
			printf("%s\t",a->taxon_name);
			for (j=1;j<=n;j++)
					{
					b=(NODE)(pListgetkthNode(nodeList, j)->item);
					c=mrca(a,b);
					T=pathLengthTimeModel(c,root,model);
					printf("%f\t",T);
					vcv[i][j]=T;
					}
			printf("\n");
			}
	freepList(nodeList);
	return vcv;						
}
/*****************************************************************************************************/

NODETYPE * nextRndNode(long nNodes, NODETYPE ** nodeArray)

/* return a randomly selected node which has not yet been marked. Dumbly keeps looking for available nodes,
	and will do so forever if they've all been sampled. Very brute force routine--stupid even! Better check in calling routine ! */

{
NODETYPE * node;
long rn;
do 
	{
	rn=rndLong(nNodes);
	node = nodeArray[rn];
	} while (isNodeMarked(node));
// markNode(node); TEMPORARY, hope this works
return node;
}

/*****************************************************************************************************/
void markNode( NODETYPE  * n)
{
n->isQueryNode=1;
return;
}
void unMarkNode( NODETYPE  * n)
{
n->isQueryNode=0;
return;
}
int isNodeMarked( NODETYPE  *n)
{
return n->isQueryNode;
}
void unMarkTree(TREE T)
{
preOrderVoid(T->root,unMarkNode);
return;
}

void setLocalModel(NODETYPE *n,int model,int stemFlag)
{
int saveModel;
saveModel=n->modelID;
	
preOrderIntArg(n,setModel,model);
if (/* !isTip(n) &&*/ stemFlag==0)
	n->modelID=saveModel; /* for interior nodes, if we DONT want the stem rate assigned we better have saved it,
					because the recursion will assign it */
return;
}

static void setModel(NODE n,int model)
{
n->modelID=model;
return;
}


NODETYPE * sister(NODETYPE * n)
{
if (n->anc == NULL)
	return NULL;
if (n->anc->firstdesc==n)
	return n->sib;
else
	return prevSib(n);


}

void copyLengthsTree2Tree(NODETYPE * node1,NODETYPE * node2)
{
NODETYPE * child1, *child2;
node1->nodeReal=node2->length;
if (!isTip(node1))
	{
	child1=node1->firstdesc;
	child2=node2->firstdesc;
	for (;child1;child1=child1->sib,child2=child2->sib)
		copyLengthsTree2Tree(child1,child2);
	}
return;
}

/***********************************************************************************/
void preOrderVoid(NODETYPE *node,void (*f)(NODETYPE *))
{
	NODETYPE *child;
	(*f)(node);
	if (!isTip(node))
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			preOrderVoid(child,f);
		}
	return ;	
}
double preOrder(NODETYPE *node,double (*f)(NODETYPE *))
{
	double sum=0;
	NODETYPE *child;
	sum+=(*f)(node);
	if (!isTip(node))
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			sum += preOrder(child,f);
		}
	return (sum);	
}
double preOrderArg(NODETYPE *node,double (*func)(NODE node, double farg),double farg)
{
/* */

	double sum=0;
	NODETYPE *child;
	sum+=(*func)(node,farg);
	if (!isTip(node))
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			sum += preOrderArg(child,func,farg);
		}
	return (sum);	
}
static void preOrderIntArg(NODETYPE *node,void (*func)(NODE node, int iarg),int iarg)
{
/* */

	NODETYPE *child;
	(*func)(node,iarg);
	if (!isTip(node))
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			preOrderIntArg(child,func,iarg);
		}
	return;	
}
/***********************************************************************************/
void setNodeEstRate(NODE node)

/* for PENLIKET method, this is a kludge which ends up calculating the mean branch rates
	and storing them in the usual place, so ratogram draws are possible. Also stores
	node rate estimates in node->nodeEstRate field */

{
(void)preOrder(node,moveNER);
(void)preOrder(node,setNER);
return;
}
static double moveNER(NODE node)
{
node->nodeEstRate=node->estRate;
return 0.0;
}
static double setNER(NODE node)
{
if (!isRoot(node))
	node->estRate=(node->nodeEstRate+node->anc->nodeEstRate)/2.0;
return 0.0;
}

/***********************************************************************************/


void zeroEstRate(NODETYPE *node)
{
(void)preOrder(node,zeroER);
return;
}
static double zeroER(NODETYPE *node)
{
node->estRate=0.0;
node->nodeEstRate=0.0;
return 0.0;
} 

/***********************************************************************************/

double LFuncons(NODETYPE *node)
{

return preOrder(node,LFuncons1);

}

double LFuncons1(NODETYPE *node)
{
double expected, chiSq=0.0;
if (!isRoot(node))
	{
	expected=node->estRate*(node->anc->time-node->time);
  	if (fabs(expected) > 0.0001)  
		{
		chiSq= SQR(node->length - expected)/expected;
		}
/*printf("SAME expected=%f nodelength=%f node->nodeReal=%f chiSq=%f\n",expected,node->length,node->nodeReal,chiSq);*/
	}
return chiSq;
}
double LFunconsT(NODETYPE *node)
{

return preOrder(node,LFuncons1T);

}

double LFuncons1T(NODETYPE *node)
{
double expected, chiSq=0.0,r;
if (!isRoot(node))
	{
	r=(node->estRate+node->anc->estRate)/2;
	expected=r*(node->anc->time-node->time);
  	if (fabs(expected) > 0.0001)  
		{
		chiSq= SQR(node->length - expected)/expected;
		}
/*printf("SAME expected=%f nodelength=%f node->nodeReal=%f chiSq=%f\n",expected,node->length,node->nodeReal,chiSq);*/
	}
return chiSq;
}

/***********************************************************************************/
double cvSquareErrorBranch(TREE t, NODE node,int method,double *chiSq)

/* 

For all methods, uses estimated rates and times to calc the prediction error of branch subtending node n.

Prediction error calculated as follows: 

LF: uses the overall tree-wide estimated rate

NP and PL: uses the estimated rate of the branch BELOW the removed branch, OR
		if the pruned branch's ancestor is the root, uses the estimated rate of
		the pruned branch's sister branch.

For 0-length branches, just returns 0 for ChiSq if expected is 0.

 */

{
NODETYPE *n;
double sq=0.0;
double d,expected;
if (isRoot(node))return 0.0;
*chiSq=0.0;
d=node->anc->time-node->time;
switch (method)
	{
	case LaF:
		expected=t->estRate*d;
		sq=SQR(expected - node->length);
		if (expected == 0.0)
			*chiSq=0.0;
		else
			*chiSq=sq/expected;
		break;
	case NP:
		if (!isRoot(node->anc)) 
			{
			expected=node->anc->length/(node->anc->anc->time - node->anc->time)*d;
			sq=SQR(expected - node->length);
			if (expected == 0.0)
				*chiSq=0.0;
			else
				*chiSq=sq/expected;
			}
		else
			{
			n=sister(node);
			expected=n->length/(n->anc->time - n->time)*d;
			sq=SQR(expected - node->length);
			if (expected == 0.0)
				*chiSq=0.0;
			else
				*chiSq=sq/expected;

			/*!!!!!!!!  Put in code to handle case where rate is at the root */
			}
		break;


		break;
	case PENLIKE:
		if (!isRoot(node->anc)) 
			{
			expected=node->anc->estRate*d;
			sq=SQR(expected - node->length);
			if (expected == 0.0)
				*chiSq=0.0;
			else
				*chiSq=sq/expected;
/*
printf("%s rate=%f expected=%f dur=%f length=%f SqEr=%f chiSq=%f\n",node->taxon_name,node->anc->estRate,expected,d,node->length,sq,*chiSq);
*/
			}
		else
			{
			expected=sister(node)->estRate*d;
			sq=SQR(expected - node->length);
			if (expected == 0.0)
				*chiSq=0.0;
			else
				*chiSq=sq/expected;
			}
		break;
	case PENLIKET: /* NOT UPDATED YET TO HANDLE ROOT ISSUES */
		expected=node->anc->estRate*d;
		sq=SQR(expected - node->length)/expected;
		/* IMPORTANT: For this method, we use this node's ancestor for the rate. The node
			itself was a tip, and was deleted from the analysis ! */
		break;
	}
return sq;

}
#if 0
/***********************************************************************************/
double cvSquareError(TREE t, int method)

/* For all methods, uses estimated rates and times to calc the prediction error of branch
	lengths that are stored in the variable node->nodeReal */

{
double sq=0.0;
extern int gNComp;
switch (method)
	{
	case LaF:
		(void)preOrderArg(t->root,cvSELF,t->estRate/(gNComp-1));
		sq=preOrder(t->root,cvCS);
		break;
	case NP:
		sq=preOrder(t->root,cvSQNP);
		break;
	case PENLIKE:
		sq=preOrder(t->root,cvCS);
		break;
	case PENLIKET:
		sq=preOrder(t->root,cvSQET);
		break;
	}
printf("FINAL CALC ON SQ ERR= %f\n",sq);
return sq;

}
/***********************************************************************************/


static double cvSQNP(NODE node)
{
double expected, Sq=0.0,d,estRate;
extern int gNComp;
if (!isRoot(node))
	{
	expected=node->length/(gNComp-1); /* seems trivial, but correct for the NP model */
  	if (fabs(expected) > 0.0001)  
		{
		Sq= SQR(node->nodeReal - expected);
		}
	}
return Sq;
}

/***********************************************************************************/

static double cvSELF(NODE n,double rate)
{
n->estRate=rate;
return 0.0;
}

/***********************************************************************************/

static double cvCS(NODETYPE * node)
{
double expected, chiSq=0.0,d;
extern int gNComp;
if (!isRoot(node))
	{
	d=node->anc->time-node->time;
	expected=node->estRate*d/(gNComp-1);
  	if (fabs(expected) > 0.0001)  
		{
		chiSq= SQR(node->nodeReal - expected);
		}

printf("rate=%f expected=%f dur=%f length=%f node->nodeReal=%f SqEr=%f\n",node->estRate,expected,d,node->length,node->nodeReal,chiSq);


	}
return chiSq;
}
/***********************************************************************************/


double cvCST(NODETYPE * node)
{
double expected, chiSq=0.0,d,r;
if (!isRoot(node))
	{
	d=node->anc->time-node->time;
	r=(node->estRate+node->anc->estRate)/2;
	expected=r*d;
  	if (fabs(expected) > 0.0001)  
		{
		chiSq= SQR(node->nodeReal - expected) /expected ;
		}

/*printf("rate=%f expected=%f dur=%f length=%f node->nodeReal=%f chiSq=%f\n",node->estRate,expected,d,node->length,node->nodeReal,chiSq);
*/

	}
return chiSq;
}
/***********************************************************************************/

static double cvSQET(NODETYPE * node)
{
double expected, sq=0.0,d,r;
extern int gNComp;
if (!isRoot(node))
	{
	d=node->anc->time-node->time;
	r=(node->estRate+node->anc->estRate)/2;
	expected=r*d/(gNComp-1);
  	if (fabs(expected) > 0.0001)  
		{
		sq= SQR(node->nodeReal - expected);
		}
	}
return sq;
}
#endif

/********************************************/

void unSetConstraints(NODETYPE * node)


{
	NODETYPE *child;

	node->nodeIsConstrainedMax=0;
	node->nodeIsConstrainedMin=0;
	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			unSetConstraints(child);
		}
	return;
}


double unFixNodeAge(NODETYPE *node) 

	/* only allowed on internals for convenience*/
{
if (!isTip(node))
	node->free=1;
return 0.0;
} 


double fixNodeAge(NODETYPE *node)
{
node->free=0;
return 0.0;
} 

int numFixedNodes(NODETYPE * node)

/* Returns number of (internal) descendents of the clade at node with fixed ages */

{
	NODETYPE *child;
	int numFix=0;
	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			numFix+=numFixedNodes(child);
		if (isFixed(node))
			++numFix;
		}
	return numFix;
}
int numConstrainedNodes(NODETYPE * node)

/* Returns number of (internal) descendents of the clade at node with constrained (not fixed) ages */

{
	NODETYPE *child;
	int numCon=0;
	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			numCon+=numConstrainedNodes(child);
		if (isConstrained(node))
			++numCon;
		}
	return numCon;
}

void setupFixedNodeArray(NODETYPE * node, NODE nodeArray[], int *curIndex)

/* Populates a node array with the nodes that are currently constrained. */

{
	NODETYPE *child;
	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			setupFixedNodeArray(child,nodeArray,curIndex);
		if (!isFree(node))
			nodeArray[(*curIndex)++]=node;
		}
	return;
}

void setupConstrainedNodeArray(NODETYPE * node, NODE nodeArray[], int *curIndex)

/* Populates a node array with the nodes that are currently constrained. */

{
	NODETYPE *child;
	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			setupConstrainedNodeArray(child,nodeArray,curIndex);
		if (isConstrained(node))
			nodeArray[(*curIndex)++]=node;
		}
	return;
}

int maxAgePresent(NODETYPE * node)

/* Returns a 1 if any descendents of the clade at node have a max age constraint set */

{
	NODETYPE *child;

	if (node->nodeIsConstrainedMax)
		return 1;

	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			if(maxAgePresent(child))
				return 1;
		}
	return 0;
}
int constraintsPresent(NODETYPE * node)

/* returns a 1 if any descendents of the clade at node have time constraints set */

{
	NODETYPE *child;

	if (node->nodeIsConstrainedMax || node->nodeIsConstrainedMin)
		return 1;

	if (!isTip(node)) 
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			if(constraintsPresent(child))
				return 1;
		}
	return 0;
}
int tipsDifferentAges(NODETYPE *node)

/* determines whether all the tips are the same age: 1=different, 0=same */

{
	static int first=1;
	static double save;
	NODETYPE *child;
	if (isTip(node)) 
		{
		if (first)
			{
			save=node->time;
			first=0;
			}
		else
			if (save != node->time)
				return 1;
		}
	child=node->firstdesc;
	SIBLOOP(child) 
		if(tipsDifferentAges(child))
			return 1;
	return 0;
}

/***********************************/
void updateSubtrees(NODETYPE *srcNode)

/* Copy all the current times on the source tree to the various subtree data structures 
	defined from it */

{
	NODETYPE *child,*subRoot;
	if (!isTip(srcNode) && !isRoot(srcNode))
			{
			subRoot=srcNode->nodePtr;
			if (subRoot)
				updateOneSubtree(subRoot);	/* update the subtree for this source tree node */
			}
	child=srcNode->firstdesc;
	SIBLOOP(child)
			updateSubtrees(child);
	return;
}

static void updateOneSubtree(NODETYPE *subNode)
{
	NODETYPE *child,*srcNode;
	subNode->time=(subNode->nodePtr)->time;		/* subNode->nodePtr maintains a pointer to the corresponding source tree node */
	child=subNode->firstdesc;
	SIBLOOP(child)
		updateOneSubtree(child);
	return;


}


/***********************************/

void setupSubtrees(NODETYPE * srcNode)

/* for each internal node of tree rooted at srcNode, setup a subtree and store a pointer to this
subtree in the internal node's nodePtr location */

{
	NODETYPE *child;
		{
		if (!isTip(srcNode) && !isRoot(srcNode))
			srcNode->nodePtr=createSubtree(srcNode,0);
		child=srcNode->firstdesc;
		SIBLOOP(child)
			setupSubtrees(child);
		}
	return;
}

/***********************************/

NODETYPE *createSubtree(NODETYPE *srcNode, int SubtreeSize)

/* Returns a pointer to a newly allocated tree, which is created by copying a 
subtree from tree srcRoot. Copies pertinent time information from source nodes.
Each node also stores a pointer to the node on the source tree from whence it came.
This permits rapid updating of information about time.

At the moment this routine ignores SubtreeSize, and makes a subtree from three branches
surrounding srcNode.

*/
{
NODETYPE *root, *cnode,*node,*child;
if (!isTip(srcNode) && !isRoot(srcNode)) /* only allow subtrees from internal nodes! */
	{
	root=newnode();
	copyNodeInfo(srcNode->anc,root);
	root->nodePtr=srcNode->anc;
	cnode=newnode();
	AddChild(root,cnode);
	copyNodeInfo(srcNode,cnode);
	cnode->nodePtr=srcNode;
	child=srcNode->firstdesc;
	SIBLOOP(child)  /* for each child of the source node, create a child on the copied tree */
		{
		node=newnode();
		AddChild(cnode,node);
		copyNodeInfo(child,node);
		node->nodePtr=child;		
		}
	return root;
	}
else
	return NULL;
}

/**********************************/
NODETYPE*  copyTree(NODETYPE* a)
// returns a node that is either a tip, or the root of a properly formatted tree, but its ancestor and sibs are undefined
{
NODETYPE* child,*first,*newfirst,*newn,*n,*prev;
newn = newnode();
copyNodeInfo(a,newn);
if(!isTip(a))
	{
	first=a->firstdesc;
	newfirst=copyTree(first);
	newn->firstdesc=newfirst;
	newfirst->anc = newn;
	prev=newfirst;
	child=first->sib; // start loop with the second sib in the sib list...
	SIBLOOP(child)
		{
		n = copyTree(child);
		prev->sib=n;
		prev=prev->sib;
		n->anc = newn;
		}
	}
return  newn;
}

/**********************************/



void AddChild(NODETYPE * parent, NODETYPE * theChild)
        {
	NODETYPE *aChild;
	if (parent->firstdesc)
	    {
	    aChild=parent->firstdesc;
	    if (aChild)
		    {
		    while(aChild->sib)
			    aChild=aChild->sib;
		    aChild->sib=theChild;
		    }
	    }
	else
	    parent->firstdesc=theChild;
	theChild->anc=parent;
	theChild->sib=NULL;
        return;
        }
void RemoveTaxonLvAnc(NODETYPE * rmTaxon)

/* remove a tip or clade, but leave its ancestor node in place */

{
NODETYPE * prev;
prev=prevSib(rmTaxon);
if(prev)	/* either rmTaxon is the firstdesc or its got a prev sib */
	prev->sib=rmTaxon->sib;
else
	rmTaxon->anc->firstdesc=rmTaxon->sib;
rmTaxon->anc=NULL;
rmTaxon->sib=NULL;
}


NODETYPE * RemoveTaxon(TREE T,NODETYPE * rmTaxon)

/* Removes a taxon, or clade, including the stem lineage
 * Does not remove the node below the stem lineage if that becomes of degree two
 * Does not deallocate memory for the subtree that is deleted, or change links on that subtree.
 * Won't allow removal of the root node
 * If the node is one of only two children of the root, the root is removed as well.
 * RETURNS A POINTER TO THE PRUNED TREE!
 */
     {
     NODETYPE *n, *prev, *parent,*grandparent,*sis,*root;
     if (T)
     	root=T->root;
     else
	root=NULL;	/* used for cases in which don't know the tree and don't care */
     if (rmTaxon==NULL)
	return root;
     if (!isRoot(rmTaxon))
        {
	parent=rmTaxon->anc;
	grandparent=parent->anc; /* might be NULL if parent is the root */
	if (node_tomy(parent)==2)
		{
		sis=sister(rmTaxon);
		if (isRoot(parent)) 
			{
			sis->anc=NULL;	/* make sure this node acquires 'root' status */
			sis->sib=NULL;
			return sis; /* new root of tree is this sister node */
			}
		else
			{
			sis->anc=grandparent;
			prev=prevSib(parent);
			if (prev)
				{
				prev->sib=sis;
				}
			else
				{
				grandparent->firstdesc=sis;
				}
			sis->sib=parent->sib;
			sis->length+=parent->length;
			}
		}
	if (node_tomy(parent)>2)
		{
		prev=prevSib(rmTaxon);
		if(prev)	/* either rmTaxon is the firstdesc or its got a prev sib */
			prev->sib=rmTaxon->sib;
		else
			parent->firstdesc=rmTaxon->sib;
		}

        return root;
        }
     }

NODETYPE * prevSib(NODETYPE* node)

/* returns the sib that points to this sib, or null if this sib is the first desc or if this sib is root */

	{
	NODETYPE *prev, *n;
	prev=NULL;
	if(!isRoot(node))
	    {
	    n=node->anc->firstdesc;
	    while(n != node)  
			    {
			    if (n->sib == NULL)
				     return NULL;   
			    prev=n;
			    n=n->sib;
			    }
	    return prev;
	    }
	else
	    return NULL;
	}

NODETYPE * ReRoot(NODETYPE * atNode)  // Sept 2011: fixed some bugs in here that I fixed in the
									// parallel code in mysmalltreelib.c

/* Reroots a tree in place, returning the node pointer to the new root. The old root node is deleted
	and a new root node is instantiated in its place. New root has id=1 and length=0.0. Nothing else
	is changed. Any time a node becomes a 1-tomy in this process, however, the branch lengths are
	combined.
	*/

    {
	NODETYPE *n, *r;
	if(isRoot(atNode))
	    return atNode; /* don't change the root */
	n=atNode->anc;
	if (!isRoot(n))
	    {
	    r=newnode();
	    r->id = 1; // this will be the new root. By convention its id=1; old root is deleted
	    r->length = 0.0; // also convention
//	    RemoveTaxon(NULL,atNode);
		RemoveTaxonLvAnc(atNode);
	    AddChild(r, atNode);
	    Flip(n);
	    AddChild(r, n);
	    n->length=0; /* leave all the length on the left root's branch */
//	    init_node_ids(r, 0);
	    return r;
	    }
	else
	    return n; /* don't change the root here either */
	   
    }
NODETYPE * ReRoot2(NODETYPE * atNode)  // Sept 2011: fixed some bugs in here that I fixed in the
									// parallel code in mysmalltreelib.c

// RoDo of ReRoot! 
// Make the node the root of the tree (as opposed to making it the sister node of the rest of the tree!).

// Note if we start with a binary root node, we keep that node, even as we reroot on other internal
// nodes. That means the rerootings often have a degree one node on one branch, but should be fine 
// for calculations

    {
	NODETYPE *n, *r;
	if(isRoot(atNode))
	    return atNode; /* don't change the root */
	Flip2(atNode);
	return atNode;
    }
void Flip2(NODETYPE *a)

// the subtree "below" node a, becomes one of the children of 'a' now.

    {
	NODETYPE * b,  *saveAnc, *parent;
	float saveLength;
	b=a->anc;	
	if (!isRoot(b))
	    {
	    Flip(b);  /* recurse until the root, then back up */		
	    }
//	RemoveTaxon(NULL,a);
	RemoveTaxonLvAnc(a);
	AddChild(a, b);
	b->length=a->length; /* flip the branch lengths too */
	return;
    }
void Flip(NODETYPE *a)

// the subtree "below" node a, becomes one of the children of 'a' now.

    {
	NODETYPE * b,  *saveAnc, *parent;
	float saveLength;
	b=a->anc;	
	if (!isRoot(b))
	    {
	    Flip(b);  /* recurse until the root, then back up */		
	    }
//	RemoveTaxon(NULL,a);
	RemoveTaxonLvAnc(a);
	AddChild(a, b);
	b->length=a->length; /* flip the branch lengths too */
#if 1
	if (node_tomy(b)==1)  /* then delete this node and combine branch lengths*/
	    {
	    saveLength=b->length;
//	    RemoveTaxon(NULL,b);
		parent = b->anc;
	    RemoveTaxonLvAnc(b);
	    AddChild(parent, b->firstdesc);
	    b->firstdesc->length+=saveLength;
	    Node_Destructor(b);
	    /* deallocate node b HERE */
	    };
#endif
	return;
    }
/***********************************************************************************/
void traverseMultiplyLength(NODETYPE * node, double multiplier,int roundflag)

/* multiply all branch lengths by a constant and round to nearest integer */

{
	NODETYPE *child;
	double value=0;
//printf("node %s:%f %f %f %i\n",node->taxon_name, node->length,value,multiplier,roundflag);
	value=node->length*multiplier;
	if (roundflag)
		{
		if (value-floor(value)<0.5)
			node->length=floor(value);
		else
			node->length=ceil(value); /* rounding to nearest integer */
		}
	else
		node->length=value;
//printf("node %s:%f %f %f %i\n",node->taxon_name, node->length,value,multiplier,roundflag);
	child=node->firstdesc;
	SIBLOOP(child) 
		traverseMultiplyLength(child, multiplier,roundflag);

	return;
}
/***********************************************************************************/
double treeDurLength(NODETYPE * node)

/* sums the branch durations over tree */

{
	NODETYPE *child;
	double dur;
	if (isRoot(node))
		dur=0.0;
	else
		dur=node->anc->time - node->time;
	child=node->firstdesc;
	SIBLOOP(child)
		dur+=treeDurLength(child);
	return dur;
}
/***********************************************************************************/
double treeLength(NODETYPE * node)

/* sums the branch lengths over tree */

{
	NODETYPE *child;
	double dur;
	if (isRoot(node))
		dur=0.0;
	else
		dur=node->length;
	child=node->firstdesc;
	SIBLOOP(child)
		dur+=treeLength(child);
	return dur;
}
/***********************************************************************************/
double treeAgeSum(NODETYPE * node)

/* sums the node ages over tree */

{
	NODETYPE *child;
	double dur;
	dur=node->time;
	child=node->firstdesc;
	SIBLOOP(child)
		dur+=treeAgeSum(child);
	return dur;
}
/**********************************************************************/

int isNodeDescendant(NODETYPE *nodeA, NODETYPE *nodeB)
/*
 * Is nodeA the strict descendant of nodeB or identical to node B?
 * Returns 1 if it is, 0 if it is not 
 * 
 */

{
NODETYPE *node;
for(node=nodeA;node;node=node->anc) 
	/* worst case, terminates when node = NULL at ancestor of root */
	{
	if (node==nodeB) return 1;
	}  
return 0;    
}
/**********************************************************************/
int group_a_clade(NODETYPE *root, StrListPtr taxaList)

/* is the specified list of taxa a clade on tree 'root'? 
	Note that the list might contain MORE taxa than are found on the tree (i.e.,
	it might be a pruned tree.  We allow this.  First, we make a new
	taxaList that contains the intersection of the taxaList and the MRCA taxa
	on the tree.  This we check to see if it's identical to the MRCA list,
	which only happens when group is consistent with that clade
*/

{
NODETYPE *mrca;   
StrListPtr mrcaTaxa,intersectTaxaList;
mrcaTaxa=newStrList();
mrca=MRCA(root, taxaList);	/* mrca of the taxa list */
if (mrca)
    {
    TreeToTaxaList(mrca,  mrcaTaxa); /* set up list of taxa actually descended from mrca node */
/*    printf("group:\n");xprintStrList(taxaList);
    printf("clade:\n");xprintStrList(mrcaTaxa);*/

#if 0  /* enable this for pruning as described above */
    intersectTaxaList=string_list_intersect(mrcaTaxa,taxaList);
    if (string_lists_same(intersectTaxaList, mrcaTaxa))
#else
    if (string_lists_same(taxaList, mrcaTaxa))
#endif
	return 1;/* are they the same?; if so, is monophyletic */
    else 
	return 0;
    }  
}


/**********************************************************************/

NODETYPE * MRCA(NODETYPE *root, StrListPtr taxaList)
/*
 * On tree with 'root',  returns node of the MRCA of taxa in taxa List (a list of strings) 
 * NOTE: some taxa in list may not be on tree!  In that case we find the MRCA of those that
 * ARE on the tree.  If none of the taxa are on the tree return NULL.(BOMBS)
 * YIKES, does this really work? */

{
NODETYPE *node, *firstTaxonNode=NULL, *otherTaxonNode;
PtrList pOther, p, nodeList=NULL, nLptr;
int nList, k, i;
StrListPtr txPtr;
NODETYPE *s;

nList=lengthList(taxaList);
if (nList<2)
	{
    	doGenericAlert("taxa list has fewer than two taxa!");
	return NULL;
	}


/** convert the taxa list to a (possibly smaller) list of corresponding nodes **/ 

 
for(txPtr=taxaList;txPtr;txPtr=txPtr->next)
    {
    s=find_taxon_name(root,txPtr->s);
    if (s) /* don't include taxa that aren't on tree */
	{
	if(!nodeList) /* create first node, or...*/
	    {
		nodeList=pNewListAlt(sizeof(NODETYPE*));
		nLptr=nodeList;
	    }
	else	    /* add a new node, if list is already there */
	    {
		pListAddNode(nodeList, sizeof(NODETYPE*));
		nLptr=nLptr->next;    
	    }
	nLptr->item=s;
	}
    }

if (!nodeList)
	return NULL;	/* BAIL IF THERE WERE NO TAXA and hence NO NODES */

if (pLengthList(nodeList)<lengthList(taxaList))
	doGenericAlert("MRCA COMMAND: Num nodes less than num labels: You probably have misspelled a taxon name!");

p=nodeList;
for (firstTaxonNode=(NODETYPE *)(p->item);!isRoot(firstTaxonNode);firstTaxonNode=firstTaxonNode->anc)
	/* traverse the ancestry path starting from taxon 1... */
	{
	 /*...and check at each node whether the other taxa are descendants of that node...*/	
	for (pOther=p->next;pOther;pOther=pOther->next)
	    {
	    otherTaxonNode=(NODETYPE *)(pOther->item);
	    if (!isNodeDescendant(otherTaxonNode, firstTaxonNode)) 
		goto a1; /* ...at least one taxon was not a descendant, so...---> */
	    }
	freepList(nodeList);
	return firstTaxonNode; /* all members of list were descendants of this node, so return */
a1:	;   /* ----> ....so traverse to its ancestor and repeat...*/
	}  
freepList(nodeList);
return firstTaxonNode;   /* now the root: just return that by default */ 
}

/**********************************************************************/
NODE  mrca(NODE a, NODE b )

// the mrca of two nodes

{
NODE p,psave;
for (p=a;p;p=p->anc)
	p->nodeFlag=1;
for (p=b;p;p=p->anc)
	if (p->nodeFlag==1)
		goto a1;
a1:psave=p;
for (p=a;p;p=p->anc)
	p->nodeFlag=0;

return psave; 
}
double pathLengthTimeModel(NODE a, NODE b, int model)

// The sum of durations between two nodes: only include branches that match the modelID
{
NODE  p,anc,c;
double T=0.0;
c = mrca (a,b);
for (p=a;p!=c;p=p->anc)
	if (p->modelID == model)
		T += p->anc->time - p->time;
for (p=b;p!=c;p=p->anc)
	if (p->modelID == model)
		T += p->anc->time - p->time;
return T;
}
double pathLengthTime(NODE a, NODE b)

// The sum of durations between two nodes
{
NODE  p,anc,c;
double T=0.0;
c = mrca (a,b);
for (p=a;p!=c;p=p->anc)
	T += p->anc->time - p->time;
for (p=b;p!=c;p=p->anc)
	T += p->anc->time - p->time;
return T;
}
/**********************************************************************/
void setNodeName(NODETYPE *node, char *name)
{
char *copy;
copy=DupStr(name);
myFree(node->taxon_name);
node->taxon_name=copy;
return;

}

/**********************************************************************/
void make_parens(NODETYPE *node, int flag)

/* writes a parens formatted tree description with labels and durations or
lengths.  flag=0: print lengths; flag =1: print durations as lengths,  
flag=2: print rates as lengths, flag=3: print node id's as lengths,
flag=4: print normalized marginal of CLmarg[0] as length, and anc state as second number*/

{
  extern long gNumSites;
  double value, duration;
  int width;
  
  if (flag==4)
  	value = (node->CLmarg)[0]/((node->CLmarg)[0]+(node->CLmarg)[1]+(node->CLmarg)[2]);
  if (!isRoot(node))
	{
	if (flag==0)
		value = node->length;
	else if (flag == 1)
		value = node->anc->time - node->time; /* duration */
	else if (flag == 2)
		value = node->estRate/gNumSites;		/*rate*/
	}

  if (isTip(node)) 
    {
    if (*(node->taxon_name)=='\0')
		{
		width = log10(node->id)+1; 
		printf("tx%-*i", width, node->id);
		}
    else
      	printf("%s",node->taxon_name);
    if (flag == 3)
		printf(":%i",node->id);
    if (flag < 3)
    	printf(":%-8.6f",value);
    if (flag == 4)
		printf(":%-8.6f:%i",value,node->opt);
    }
  else printf("(");

  if (node->firstdesc) make_parens(node->firstdesc,flag);

  if (!isTip(node))
    {
      printf(")");
      if (*(node->taxon_name)!='\0') 
	    printf("%s",node->taxon_name);
      if (!isRoot(node)) 
	 	{
         if (flag == 3)
			printf(":%i",node->id);
         if (flag < 3)
	    	printf(":%-8.6f",value);
		}
	  if (flag==4)
			printf(":%-8.6f:%i",value,node->opt);
    }

  if (node->sib) printf(","),make_parens(node->sib,flag);

}
/***********************************************************************************/
void TreeStats(NODETYPE *root, int * numtips, int * numinternals, int * roottomy)

/* gets some info on a tree, including the number of tips, internal nodes (incl. root),
and the number of immediate descendants of the root node, the roottomy level */

{
	NODETYPE *child;
	*roottomy=0;
	child=root->firstdesc;
	SIBLOOP(child) 
		++(*roottomy);
	*numtips = numdesc(root);
	*numinternals = numIntNodes(root);

	return;
}
/***********************************************************************************/
int node_tomy(NODETYPE *node)

/* number of immediate descendants of this node (including this one!) */

{
	NODETYPE *child;
	int tomy=0;
	child=node->firstdesc;
	SIBLOOP(child) 
		++tomy;

	return tomy;
}
/***********************************************************************************/
int maxorder(NODETYPE *node)
{
	int max,temp;
	NODETYPE *child;
	if (!node) return(-1);
	if (isTip(node)  ) {node->order=0; return (0);}
	max=0;
	child=node->firstdesc;
	SIBLOOP(child) {
			temp=maxorder(child);
			if (temp > max) max = temp;
			}
	node->order=max+1;
	return (max+1);
}
/***********************************************************************************/
void init_free(NODETYPE *node)

/* Initializes the free flag for each node.  By default free is set to 0 for tips, 1 for internal nodes and root*/

{
	NODETYPE *child;
	if (isTip(node)) 
		{
		node->free=0; 
		return;
		}
	else 
		node->free=1;		
	child=node->firstdesc;
	SIBLOOP(child) {
			init_free(child);
			}
	return;
}
/***********************************************************************************/
int numFreeNodes(NODETYPE *node)
{
/* returns number of  nodes in the tree that have their free flags set, meaning
that we are estimating their ages.  Counts the root and tips too! */

	int sum=0;
	NODETYPE *child;
	if (isFree(node)) 
		++sum;
	if (isTip(node))
		return sum;
	child=node->firstdesc;
	SIBLOOP(child) 
		sum += numFreeNodes(child);
	return (sum);	
}

/***********************************************************************************/
void Tree_Initialize(TREE aTree, char *TD, char *name)
{
NODETYPE * root;
int numtips, numinternals, roottomy;
root=string_to_tree(TD);
if (root)
	{
	init_node_ids(root, 0);
	init_free(root); /* sets default to estimate all internal nodes but root */ 
	aTree->root=root;
	aTree->name=DupStr(name);
	aTree->TD=DupStr(TD);
	TreeStats(root, &numtips, &numinternals, &roottomy);
	aTree->numTaxa=numtips;
	aTree->numBranches=numBranches(root);
	aTree->basalTomy=roottomy;
	aTree->cladeSet=NULL;
	aTree->est_b=0.0;
	aTree->est_c=0.0;
	root->anc=NULL;
	aTree->timesAvailable=0;
	aTree->method=USER;
	}
return;
}
/***********************************************************************************/
TREE Subtree_Initialize(TREE T,NODETYPE *node)

// Creates a "tree" by using some subclade of an existing tree. 'node' is the node on existing tree
// that will become the root. Does not change any of the information on the existing tree, but does
// NOT allocate a copy of this subtree--merely uses existing data structure for tree and allocates extra info 
// Careful allocating the atree. Still getting used to not using TREE as object init, as one can in C++
{
TREE aTree;
aTree=(struct treetype *)myMalloc(sizeof(struct treetype));
if (aTree)
	{
	aTree->name=T->name;
	aTree->root=node;
	aTree->numTaxa=T->numTaxa;
	aTree->numBranches=T->numBranches;
	aTree->basalTomy=T->basalTomy;
	aTree->est_b=0.0;
	aTree->est_c=0.0;
	aTree->root->anc=NULL;
	aTree->timesAvailable=0;
	aTree->method=USER;
	}
return aTree;
}
/***********************************************************************************/
void Tree_Destructor(TREE aTree)
{
DisposeTree(aTree->root);
myFree(aTree->name);
myFree(aTree->TD);
/* should free the cladeSet array if present!!! */
myFree(aTree);
return;
}
/***********************************************************************************/
void Node_Destructor(NODETYPE *node)
{
    
	myFree(node->taxon_name);	
	myFree(node);
	return;    
}




/***********************************************************************************/
void DisposeTree(NODETYPE *node)

// This used to be broken! I think this works now but haven't tested it...

	/* Frees up the tree memory and its taxon names */
{
	NODETYPE *child;
	if (!node) return;

	DisposeTree(node->firstdesc);
	DisposeTree(node->sib);
	Node_Destructor(node);
	return;
}
/***********************************************************************************/
NODETYPE *makegroup(void)  
{

	/* Returns a pointer to a tree structure corresponding to everything within a
	parentheses formatted string whose first character is at address 'gStringptr'.
	This function is called recursively each time a left parenthesis is encountered.
	NOTE: the 'gStringptr' MUST point to the first left parens on entry to this function.
	On exit 'gStringptr' should point to the rightmost right parens in the group--
	or to the last character in the name or number after colon; this
	makes it ready to skip to the next character and continue parsing 
	
	NOTE THAT THIS ROUTINE CANNOT HANDLE AN INTERNAL NODE WITH BOTH A NAME AND A LENGTH

	STUPIDLY, THIS ROUTINE DOES NOT USE TOKENS, SO IT DOESN'T PROPERLY TAKE CARE OF IMBEDDED
	SINGLE QUOTES OR BRACKETS, ETC.


	EVEN WORSE, the storage of numbers is corrupted if there is a space between colon and number
	or if the number starts with .xxx rather than 0.xxx.

*/

	NODETYPE *root, *currnode, *prevnode;
	extern char *gStringptr;	/* points to current position in string tree description
								and must be global for recursive calls to work right */
	char *character, *name, *delim=" ,):"; /* taxon name delimiters include space, comma, parens,colon */
	char *singleQuote = "'";
	char* dummy;
	extern int gCount;
	size_t length;
	int first;
	root=newnode();if (root==NULL) return(NULL);
	currnode=root;
	first=1;
	while (*gStringptr != '\0') {
		++gStringptr;
		switch (*gStringptr)  {
			case(LEFTPARENS):{    /* recursively go down into next clade */
				prevnode=currnode;
				currnode=makegroup();
				if(currnode==NULL) return(NULL);
				if (first) {
							prevnode->firstdesc=currnode;
							first=0;
							}
				else prevnode->sib=currnode;
				currnode->anc=root;
				break;
				}
			case(RIGHTPARENS): /* check to see if there is a taxon name after the parens
								OR a number after a colon, and store. 
									First letter must follow colon */
				{
					++gStringptr;  /* look ahead */
					if (*gStringptr == COLON)
						{
							++gStringptr;
							if (isdigit(*gStringptr) || (*gStringptr == '-') || (*gStringptr == '+') ) {  /* only checks first char !!!*/
								root->length=strtod(gStringptr,&dummy);
								length=strcspn(gStringptr,delim);
								gStringptr+=length-1;	/* increment but 
									leave at last character rather than
									one past, to fulfill the definition of function
									above */
								if (root->length <0.0)
								    {
								    printf("** WARNING: A negative branch length was set to ZERO\n");
								    root->length=0.0;
								    }
							}
				
						}
					else
						if (isalnum(*gStringptr))
						/* RECENTLY CHANGED THIS TO 'isalnum'from isalpha*/
							{
							length=strcspn(gStringptr,delim);
							name=(char *)myMalloc((length+1)*sizeof(char));
		
							myFree(root->taxon_name); 
		
							root->taxon_name=name;	
							if (name==NULL) return(NULL);
							strncpy(root->taxon_name,gStringptr,length);
							*((root->taxon_name)+length)='\0';  
							
							gStringptr+=length-1;	/* see comment above */

							}
						else
							--gStringptr;	/* it was neither a name or number */
					return(root);
				}
			default:{  /* check for valid taxon name or number after name and store */
					if (*gStringptr == COLON)
						{
							++gStringptr;
							if (isdigit(*gStringptr) || (*gStringptr == '-') || (*gStringptr == '+')) {  /* only checks first char !!!*/
								currnode->length=strtod(gStringptr,&dummy);
								if (currnode->length <0.0)
								    {
								    printf("** WARNING: A negative branch length was set to ZERO\n");
								    currnode->length=0.0;
								    }
								length=strcspn(gStringptr,delim);
								gStringptr+=length-1;	/* increment but
									see comment above for explanation of -1 */
							}
				
						}






			else
				if (IsItAlphaNum(*gStringptr)) {  /* only checks first letter !!!*/
					prevnode=currnode;
					currnode=newnode();
					if (currnode==NULL) return(NULL);
					if (first) {
								prevnode->firstdesc=currnode;
								first=0;
								}
					else prevnode->sib=currnode;
					currnode->anc=root;
					length=strcspn(gStringptr,delim);
					name=(char *)myMalloc((length+1)*sizeof(char));

					myFree(currnode->taxon_name); 

					currnode->taxon_name=name;	
					if (name==NULL) return(NULL);
					strncpy(currnode->taxon_name,gStringptr,length);
					*((currnode->taxon_name)+length)='\0';  
					
					gStringptr+=length-1;	/* increment but only if two or more characters */
					}
				}
		
		}
	}
}

/***********************************************************************************/
NODETYPE *makegroup2(void)  
{

	/* Returns a pointer to a tree structure corresponding to everything within a
	parentheses formatted string whose first character is at address 'gStringptr'.
	This function is called recursively each time a left parenthesis is encountered.
	NOTE: the 'gStringptr' MUST point to the first left parens on entry to this function.
	On exit 'gStringptr' should point to the rightmost right parens in the group--
	or to the last character in the name or number after colon; this
	makes it ready to skip to the next character and continue parsing 
	
	NOTE THAT THIS ROUTINE CANNOT HANDLE AN INTERNAL NODE WITH BOTH A NAME AND A LENGTH

	STUPIDLY, THIS ROUTINE DOES NOT USE TOKENS, SO IT DOESN'T PROPERLY TAKE CARE OF IMBEDDED
	SINGLE QUOTES OR BRACKETS, ETC.

	BUG! If there is a blank in a taxon name within single quotes, the name will not be parse right.
*/

	NODETYPE *root, *currnode, *prevnode;
	extern char *gStringptr;	/* points to current position in string tree description
								and must be global for recursive calls to work right */
	char *character, *name, *delim=" ,):"; /* taxon name delimiters include space, comma, parens,colon */
	char* dummy;
	extern int gCount;
	size_t length;
	int first;
	root=newnode();if (root==NULL) return(NULL);
	currnode=root;
	first=1;
	while (*gStringptr != '\0') {
		++gStringptr;
		switch (*gStringptr)  {
			case(LEFTPARENS):{    /* recursively go down into next clade */
				prevnode=currnode;
				currnode=makegroup2();
				if(currnode==NULL) return(NULL);
				insertNode(currnode, prevnode);	/* add the new node (actually,
				    the whole subtree!) to prevnode,  its ancestor */
				break;
				}
			case(RIGHTPARENS): /* check to see if there is a taxon name after the parens
								OR a number after a colon, and store. 
									First letter must follow colon */
				{
					++gStringptr;  /* look ahead */
					if (*gStringptr == COLON)
						{
							++gStringptr;
							if (isdigit(*gStringptr) || (*gStringptr == '-') || (*gStringptr == '+') ) {  /* only checks first char !!!*/
								root->length=strtod(gStringptr,&dummy);
								length=strcspn(gStringptr,delim);
								gStringptr+=length-1;	/* increment but 
									leave at last character rather than
									one past, to fulfill the definition of function
									above */
								if (root->length <0.0)
								    {
								    printf("** WARNING: A negative branch length was set to ZERO\n");
								    root->length=0.0;
								    }
							}
				
						}
					else
						if (isalnum(*gStringptr))
						/* RECENTLY CHANGED THIS TO 'isalnum'from isalpha*/
							{
							length=strcspn(gStringptr,delim);
							name=(char *)myMalloc((length+1)*sizeof(char));
		
							myFree(root->taxon_name); 
		
							root->taxon_name=name;	
							if (name==NULL) return(NULL);
							strncpy(root->taxon_name,gStringptr,length);
							*((root->taxon_name)+length)='\0';  
							
							gStringptr+=length-1;	/* see comment above */

							}
						else
							--gStringptr;	/* it was neither a name or number */
					return(root); /*... of the current subtree...*/
				}
			default:
			    {  /* check for valid taxon name or number after name and store */
				    if (*gStringptr == COLON)
						{
							++gStringptr;
							if (isdigit(*gStringptr) || (*gStringptr == '-') || (*gStringptr == '+')) {  /* only checks first char !!!*/
								currnode->length=strtod(gStringptr,&dummy);
								if (currnode->length <0.0)
								    {
								    printf("** WARNING: A negative branch length was set to ZERO\n");
								    currnode->length=0.0;
								    }
								length=strcspn(gStringptr,delim);
								gStringptr+=length-1;	/* increment but
									see comment above for explanation of -1 */
							}
				
						}
				    else
					if (IsItAlphaNum(*gStringptr))  /*  its a terminal, add it */
						{  /* only checks first letter !!!*/
						prevnode=currnode;
						currnode=newnode();
						if (currnode==NULL) return(NULL);
						insertNode(currnode, prevnode);
						length=strcspn(gStringptr,delim);
						name=(char *)myMalloc((length+1)*sizeof(char));
	
						myFree(currnode->taxon_name); 
	
						currnode->taxon_name=name;	
						if (name==NULL) return(NULL);
						strncpy(currnode->taxon_name,gStringptr,length);
						*((currnode->taxon_name)+length)='\0';  
						
						gStringptr+=length-1;	/* increment but only if two or more characters */
						}
			}
		
		}
	}
}
static void insertNode(NODETYPE *node,  NODETYPE* anc)

/* ....looks dangerous, if the node is a polytomy, this seems to delete some children! */

{
                node->anc=anc;
                if (anc->firstdesc == NULL)  /* this is first child of anc */
                        anc->firstdesc=node;
                else                    /* this is nth child and has a sib */
                        anc->firstdesc->sib=node;
}
/***********************************************************************************/
NODETYPE *newnode(void)
{

	NODETYPE *node;
	node=(NODETYPE *)myMalloc(sizeof(NODETYPE));		

	if (node==NULL) fatal("Toast");
	node->anc=NULL;
	node->firstdesc=NULL;
	node->sib=NULL;
	node->nodePtr=NULL;
	node->isCompactNode=0;
	node->isQueryNode=0;
	node->toggleDesc=0;
	node->taxon_name=(char *)myMalloc(sizeof(char));
	node->length=FLT_MAX;		/* Big number lets us check later */ 
	node->time=0.0; /* NOTE: 'drawtree' checks this value at the root node
		to determine if times have been set */
	node->minAge=0.0;
	node->maxAge=/*  1.0  */  1.0e20;
	node->nodeIsConstrainedMax=0;
	node->nodeIsConstrainedMin=0;
	node->free=0;
	node->like = 0.0;
	node->id=0;
	node->modelID=0;
	node->nodeFlag=0;
	

	if (node->taxon_name ==NULL) fatal("Couldn't allocate name in node");;
	*(node->taxon_name)='\0'; 	/* store a null string for now */
	return (node);
}


/***********************************************************************************/
void copyNodeInfo(NODETYPE *source,NODETYPE *dest)

/* Copies SOME information about one node to another node */

{

	dest->taxon_name=DupStr(source->taxon_name);
	dest->length=source->length;
	dest->time=source->time;
	dest->minAge=source->minAge;
	dest->maxAge=source->maxAge;
	dest->id=source->id;
	dest->free=source->free;
	dest->numdesc=source->numdesc;
	dest->estRate=source->estRate;
	dest->nodeReal=source->nodeReal;
	dest->nodeIsConstrainedMax=source->nodeIsConstrainedMax;
	dest->nodeIsConstrainedMin=source->nodeIsConstrainedMin;
	return;
}


/***********************************************************************************/
void ClusterHistogram(NODETYPE * node, long *histo,long TSize)

/* moves through a tree setting up a histogram of cluster sizes in which bins are 
 * on a log 2 scale: thus for a tree of size 32.  Note that the last bin is slightly larger...
 * 
 *	2-3
 *	4-7
 *	8-16
 */

{
	NODETYPE *child;
	int ix;
	long N;
	if (!node || isTip(node))
		return;
	if (!isRoot(node))
	    {
	    N=node->numdesc;
	    ix=floor(LOG2(MIN(N,TSize-N)))-1;
	    if (N==TSize/2)
		--ix;		/* put this special cluster size in next lowest bin...otherwise
				    it sits there almost alone in the last bin */
	    ++histo[ix];
	    }
	child=node->firstdesc;
	SIBLOOP(child) 
		ClusterHistogram(child, histo,TSize);
	return;
}


/***********************************************************************************/

int numNodes(NODETYPE *node)
{

return numdesc(node)+numIntNodes(node);


}

/***********************************************************************************/

long numUnMarkedDesc(NODETYPE *node)

/* determines the number of leaves descended from a node EXCEPT for descendant clades
	that are 'marked' */ 
/* Careful! Must unmark the root node of this subtree in the caller and then remark when done */

{
	long sum=0;
	NODETYPE *child;
	if (!node) return(-1);
	if (isNodeMarked(node))
		return 0;
	if (isTip(node)) 
		{
		return 1;
		}
	child=node->firstdesc;
	SIBLOOP(child) 
		sum+=numUnMarkedDesc(child);
	return (sum);
}
/***********************************************************************************/

int numdesc(NODETYPE *node)

/* determines the number of leaves descended from every node and stores them at node */

{
	long sum=0;
	NODETYPE *child;
	if (!node) return(-1);
	if (isTip(node)) 
		{
		node->numdesc=1; 
		return (1);
		}
	child=node->firstdesc;
	SIBLOOP(child) 
		sum+=numdesc(child);
	node->numdesc=sum;
	return (sum);
}
/***********************************************************************************/
void printnode(NODETYPE *node)
{
    double duration, estR;
    NODETYPE *anc;
    if (isNodeMarked(node))
	printf("*");
    if (!isRoot(node))
	{
	anc=node->anc;
	duration=node->anc->time-node->time;
	printf("node %3i (%s) age=%4.2g | anc %3i (%s) age=%4.2g | dur=%4.2g len=%4.2g rate=%6.3g nodeReal=%6.3g age bounds=[%g..%g]\n",
	    node->id, node->taxon_name,node->time, 
	    anc->id, anc->taxon_name, anc->time, 
	    duration,node->length,node->estRate,node->nodeReal,nodeUpperBound(node),nodeLowerBound(node));
	}
    else
	printf("node %3i (%s) age=%4.2f len=%4.2g\n",
	    node->id, node->taxon_name,node->time, node->length);
    return;
}
/***********************************************************************************/
void printtree(NODETYPE *node)
{
	NODETYPE *child;
	if (!node) 
		return;
	printnode(node);
#if 0
	if (node->nodePtr)
		{
		printf("\n\n --->Subtree info for this node:\n\n");
		DrawTree(node->nodePtr,0,0); /*printtree(node->nodePtr);*/
		printf("\n\n --->End of subtree info for this node:\n\n");
		}
#else
	if (node->nodePtr)
		printf("-->Points to node with label %s\n",node->nodePtr->taxon_name);
#endif
	child=node->firstdesc;
	SIBLOOP(child) 
		printtree(child);
	return;
}
/***********************************************************************************/
void printLikes(NODETYPE *node)
{
	NODETYPE *child;
	if (!node) 
		return;
	printnodeLike(node);
	child=node->firstdesc;
	SIBLOOP(child) 
		printLikes(child);
	return;
}
/***********************************************************************************/

int numIntNodes(NODETYPE *node)
{
/* returns number of internal nodes in the tree.  Counts the root too, so this must
be subtracted later! */

	int sum=0;
	NODETYPE *child;
	if (!node) return(-1);
	if (isTip(node)) 
		return (0);
	child=node->firstdesc;
	SIBLOOP(child) 
		sum += numIntNodes(child); /* add one for each child and all that children's*/
	return (1+sum);	/* the 1 is to count this node, which must be internal */
}
/***********************************************************************************/
int numBranches(NODETYPE *node)
{
/* returns number of branches in the tree.  Does NOT count a branch subtending
the root.  Note that this number is <= 2*ntaxa-2 because of polytomies */

	return numdesc(node)+numIntNodes(node)-1; /* for the root */
}
/***********************************************************************************/
/***********************************************************************************/

NODETYPE *string_to_tree(char *tree_description)

/*  Takes a string tree description in NEXUS parens format and returns the root
	node of a linked-list tree structure.  The names of the taxa are stored
	in a string pointed to by node->taxon_name.  The routine first checks to
	see if the string description is valid--hopefully it's good at this.
	Stores lengths of branches and internal taxon names.
	
	Returns NULL if fails. 
*/ 


{
	NODETYPE *root;
	extern char *gStringptr;
	int ix;
	gStringptr=tree_description;
	if (stringcheck(tree_description)) {;
		gStringptr=strchr(tree_description,LEFTPARENS); /* move to first occurrence of left paren*/

		root=makegroup();		
		}
	else 
		return(NULL);
	return(root);
}
/***********************************************************************************/
int stringcheck(char *td)
	/*  Checks a tree description statement and returns error 
	if it has unbalanced parentheses, or if a right parens precedes a left out of turn.
	It does NOT catch a premature termination caused by an early right parens that matches
	up with a left parens */

{
	long count=0,parenscount=0;
	if (*td == '\0') return (0);  /* string is empty */
	while(*td)  {  /* while character is not null */
		if (*td == LEFTPARENS) ++parenscount;
		if (*td == RIGHTPARENS) --parenscount;				
		if (parenscount < 0) return (0); /* right parens preceded left */
		td++;  /* used to be *td++  */
		}
	if (parenscount != 0) return (0); /* unbalanced parens */
	return (1);	/* string OK  */
}
/***********************************************************************************/
NODETYPE * find_taxon_name(NODETYPE *node,char *name)
/* returns the node that has a taxon name that matches 'name' or NULL
if it fails to find */


{
	NODETYPE *child, *found_node;

	if (node->taxon_name)
		if (isEqual(name,node->taxon_name))
			return node;
	child=node->firstdesc;
	SIBLOOP(child) 
		if (found_node=find_taxon_name(child,name) )
			return found_node;
	return NULL;
}
/***********************************************************************************/

int collapseLengthsTree2Tree(TREE t1,TREE t2)

/* Moves through two identical trees simultaneously, and whenever it finds a zero-length
	branch on tree1, it collapses that branch AND the corrsponding one on tree2.
	Information about the length of the branch on tree2 is discarded! If the trees are not
	isomorphic, the results are unpredictable.
*/
{
int retFlag=0;
if (any_zero_internal_branches(t1->root))
	retFlag=1;
while(any_zero_internal_branches(t1->root))
				collapse_zero_2trees(t1->root,t2->root);

return retFlag;
}
static void collapse_zero_2trees(NODETYPE * node1, NODETYPE * node2)
{
NODETYPE *  child1,* child2;

if(!isRoot(node1) && !isTip(node1))
	{
	if (node1->length==0.0)
			{
			collapse_node(node1);
			collapse_node(node2);
			return;
			}
	}
child1=node1->firstdesc;
child2=node2->firstdesc;
for (;child1;child1=child1->sib,child2=child2->sib)
	collapse_zero_2trees(child1,child2);

return;
}
/***********************************************************************************/
int any_zero_internal_branches(NODETYPE *node)
/* returns 1 if any INTERNAL nodes, other than the root node, have zero-length branches */


{
	NODETYPE *child;

	if(!isRoot(node) && !isTip(node))
	   if (node->length==0.0)
			return 1;
	child=node->firstdesc;
	SIBLOOP(child) 
		if (any_zero_internal_branches(child) )
			return 1;
	return 0;
}
/***********************************************************************************/
int any_zero_terminal_branches(NODETYPE *node)
/* returns 1 if any TERMINAL nodes,  have zero-length branches */


{
	NODETYPE *child;

	if(isTip(node))
	   if (node->length==0.0)
			return 1;
	child=node->firstdesc;
	SIBLOOP(child) 
		if (any_zero_terminal_branches(child) )
			return 1;
	return 0;
}
/***********************************************************************************/
void collapse_zero(NODETYPE *node)
/* collapses THE FIRST FOUND zero-length branch to polytomies*/


{
	NODETYPE *child;

	if(!isRoot(node) && !isTip(node))
	   if (node->length==0.0)
			{
			collapse_node(node);
			return;
			}
	child=node->firstdesc;
	SIBLOOP(child) 
		collapse_zero(child);
	return;
}
/***********************************************************************************/
void collapse_node(NODETYPE *node)
{

// Remove the branch subtending node, collapsing to a polytomy
// Ignore if root or tip

NODETYPE *anc, *right, *left, *first_desc, *last_desc, *nd;

if (isTip(node) && isRoot(node))return;
anc=node->anc;
first_desc=node->firstdesc;
right=node->sib;  

for (nd=first_desc;nd->sib;nd=nd->sib);
last_desc=nd;/* this fragment finds node's last immediate desc. */
last_desc->sib=right;

if (anc->firstdesc==node)
    anc->firstdesc=first_desc; /* this node is the leftmost */
else			/*find the node to the left of it */
    {
    for (nd=anc->firstdesc;nd->sib!=node;nd=nd->sib);
    left=nd;
    left->sib=first_desc;
    }
for (nd=first_desc;nd!=last_desc->sib;nd=nd->sib)
    nd->anc=anc;
// Node_Destructor(node); /* WATCH OUT! THAT's for damn sure; this screws up the recursion in outer routine*/
return;    
}

/***********************************************************************************/
NODETYPE * find_id(NODETYPE *node,int id)
/* returns the node that has an id that matches 'id' or NULL
if it fails to find */


{
	NODETYPE *child, *found_node=NULL;

	if (node->id == id)
			return node;
	child=node->firstdesc;
	SIBLOOP(child) 
		if (found_node=find_id(child,id) )
			return found_node;
	return NULL;
}

/***********************************************************************************/
void init_node_ids(NODETYPE *node, int id)
{
	NODETYPE *child;
	if (!node) 
		return;
	gId=id+1;
	node->id=gId;
	child=node->firstdesc;
	SIBLOOP(child) 
		init_node_ids(child, gId);
	return;
}
/***********************************************************************************/
void print_ages(NODETYPE *node, double time, double calAge,int method)
{
    NODETYPE *child;
    extern long gNumSites;

/*...ages and misc...*/

    if (time != 0.0)
    	node->time*= (calAge/time); /* trap this occasional error */
    if (!isTip(node))
      {
      if(*(node->taxon_name))
	    printf(" [*]  %.8s\t", node->taxon_name);
      else
	    printf("     (%i)\t", node->id);
      }
    else
	printf("      %.8s\t",node->taxon_name);
    if (node->free)
	    printf("      ");
    else
	    printf("     *");
    printf(" [%1i] ",node->modelID);
    if (node->nodeIsConstrainedMin)
    	printf("\t%7.2f\t",node->minAge);
    else
	printf("\t   --   ");
    if (node->nodeIsConstrainedMax)
    	printf("%7.2f\t",node->maxAge);
    else
	printf("   --\t");

    printf("%7.2f\t\t",node->time);
    if (!isRoot(node))
      switch (method)
		{
		case USER:   /* user supplied chronogram */  	
			printf("   --   \t%.4e\n",node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case LaF:    	
			printf("%.4e\t%.4e\n",node->estRate/gNumSites,node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case LFLOCAL:    	
			printf("%.4e\t%.4e\n",node->estRate/gNumSites,node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case NP:
			printf("   --   \t%.4e\n",node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case PENLIKE:
			printf("%.4e\t%.4e\n",node->estRate/gNumSites,node->length/(node->anc->time-node->time)/gNumSites);
			break;
		default:;
		}
     else
	printf("\n");



    if(!isTip(node))
	{
    	child=node->firstdesc;
    	SIBLOOP(child) 
		print_ages(child, time, calAge,method);
	}
    return;
    
}
/***********************************************************************************/
void print_named_ages(NODETYPE *node)
{
// prints out the ages of internal named nodes only...
    NODETYPE *child;
    if (!isTip(node))
      {
      if(*(node->taxon_name))
	    printf(" [**]  %.8s\t%7.2f\n", node->taxon_name,node->time);
      }

    if(!isTip(node))
	{
    	child=node->firstdesc;
    	SIBLOOP(child) 
		print_named_ages(child);
	}
    return;
    
}
/***********************************************************************************/
void summarize_rates(TREE t)
{
extern long gNumSites;
double *r,min=1e20,max=-1e20,mean,sdev,adev,var,skew,curt;
long ix=0,i;
NODETYPE * root;
root=t->root;
r=(double *)myMalloc((t->numBranches)*sizeof(double));
recurse_summarize_rates(t->root,&ix,r);

moment(&r[0]-1,t->numBranches,&mean,&adev,&sdev,&var,&skew,&curt);
for (i=0;i<t->numBranches;i++)
		{
		if (r[i]>max)max=r[i];
		if (r[i]<min)min=r[i];
		}
printf("\nSummary of rate variation (substitutions per site per unit time)\n  Mean    = %.4g\n  Std Dev = %.4g\n  Min     = %.4g\n  Max     = %.4g\n  Range   = %.4g\n  Ratio   =  %.4g\n",mean/gNumSites,sdev/gNumSites,min/gNumSites,max/gNumSites,(max-min)/gNumSites,max/min);

return;

}

static void recurse_summarize_rates(NODETYPE * n, long * ix, double r[])
{
NODETYPE *child;
if (!isRoot(n))
	{
	r[*ix]=n->estRate;
	++(*ix);
	}
if (!isTip(n))
	{
	child=n->firstdesc;
	SIBLOOP(child)
		recurse_summarize_rates(child,ix,r);
	}
return;
}


void print_rates(NODETYPE *node,int method)
{
    NODETYPE *child;
    extern long gNumSites;
   

/* ...now rates...*/


    if (!isRoot(node))
	{
    	if(*(node->taxon_name))
	    printf("\t%.7s\t", node->taxon_name);
   	 else
	    printf("\t%i\t", node->id);

	switch (method)
		{
		case USER:   /* user supplied chronogram */  	
			printf("\t\t%.4g\n",node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case LaF:    	
			printf("\t\t%.4g\t%.4g\n",node->estRate/gNumSites,node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case LFLOCAL:    	
			printf("\t\t%.4g\t%.4g\n",node->estRate/gNumSites,node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case NP:
			printf("\t\t%.4g\n",node->length/(node->anc->time-node->time)/gNumSites);
			break;
		case PENLIKE:
			printf("\t\t%.4g\t%.4g\n",node->estRate/gNumSites,node->length/(node->anc->time-node->time)/gNumSites);
			break;
		default:;
		}
	}
    if(!isTip(node))
	{
    	child=node->firstdesc;
    	SIBLOOP(child) 
		print_rates(child,method);
	}
    return;
    
}
/***********************************************************************************/
void convert_branchlength_to_time(NODETYPE *root)

/* takes a set of ultrametric branch lengths and just scales them to [0,1] over tree and plops
them into the time field of all nodes. 
DO THIS WHEN THE BLFORMAT COMMAND HAS OPTION CLOCK=YES
*/

{
void recurseBTT(NODETYPE *node,  double scalefactor, double tcurr);
double scalefactor;
scalefactor=calcMaxToTip (root);
//printf("max to tip:%f\n",scalefactor);


recurseBTT(root,  1,  scalefactor);   // this assumes the tree's lengths are in time units

//recurseBTT(root,  scalefactor,  1.0);   


}
void recurseBTT(NODETYPE *node,  double scalefactor, double tcurr)
{
	NODETYPE *child;
	if (isRoot(node))
			node->time = tcurr;
	else
			node->time = tcurr - node->length/scalefactor;
	child=node->firstdesc;
	SIBLOOP(child) 
		recurseBTT(child, scalefactor, node->time);
	return;
    
}

void scaleTree(NODETYPE * root, double calAge, NODETYPE * calNode)

/* Scale all times on a tree according to given calibration time for one node */
{
double scaleFactor;

scaleFactor=calAge/calNode->time;
(void)preOrderArg(root,SThelper,scaleFactor);
return;
}
static double SThelper(NODETYPE * node,double factor)
{
node->time *= factor;
return 0.0;  // required by the prototype...
}

/***********************************************************************************/
static int gC;

double * sort_node_time_array(NODETYPE *root)

/* returns a pointer to an array that has the sorted INTERNAL node times in increasing order
from the present backward, including the root node time.
N.B.  There will be one element of this array for each speciation event, even if there
is a polytomy, and the diversity increases several steps at one instant in time.  Therefore,
the number of elements is equal to ntaxa-1, which happens to be the number of internal nodes
in a fully bifurcating tree.  It does not matter if the cladogram is actually bifurcating
or not! (This was a headache to resolve)*/

{
void recurse_time_get(NODETYPE *node,  double times[]);
int compar(const void *v1,  const void *v2);
int n, i;
double *times;
n=numdesc(root)-1; 
gC=0;
times=(double*)myMalloc(n*sizeof(double));
recurse_time_get(root, times);
qsort((void *)times, n, sizeof(double), compar);

/*
printf("Sorted times of internal nodes\n");
printf ("%i internal times\n", n);
for (i=0;i<n;i++)
    printf("%f\n", times[i]);
*/
return times;
}

void recurse_time_get(NODETYPE *node,  double times[])
{
	NODETYPE *child;
	int i;
	if (!isTip(node))
		{
		for (i=1;i<=node_tomy(node)-1;i++) 
		    {
		    times[gC]=node->time;
		    ++gC;
		    }
		}
	child=node->firstdesc;
	SIBLOOP(child) 
		recurse_time_get(child,times);
	return;
    
}
int compar(const void *v1,  const void *v2)
{
double V;
V= *(double *)v1 - *(double *)v2;
if (V>0.0)
	return 1;
else if (V <0.0)
	return -1;
else if (V==0.0)
	return 0;

}
double get_sum_durations(NODETYPE *node)
{
	NODETYPE *child;
	double dur=0;
	if (isTip(node))
	    return 0.0;
	child=node->firstdesc;
	SIBLOOP(child)
		{
		dur+=node->time-child->time;
		/*printf("%f %f %f\n", node->time, child->time, dur);*/
		dur+=get_sum_durations(child);
		}
	return dur;
    
}
/***********************************************************************************/
void print_tree_list(PtrList treeList)
{
PtrList lnode;
TREE thisTree;
lnode=treeList;
LISTLOOP (lnode)
	{
	thisTree=lnode->item;
	printf("Tree %s\nNum taxa = %i\nNum branches = %i\nBasal tomy = %i\n\n",
		thisTree->name,thisTree->numTaxa,thisTree->numBranches,
		thisTree->basalTomy);
	DrawTree(thisTree->root,0, 0);
	}
return;
}
/***********************************************************************************/
void TreeToTaxaList(NODETYPE *node,  StrListPtr taxaList)	
{

/* Moves through clade from node, compiling list of descendants;
on entry taxaList must be a valid pointer to a possibly empty list */

	NODETYPE *child;
	if (isTip(node)) 
		{
		appendStrList(taxaList, node->taxon_name);
		return;
		}

	else	/* interior node */
		{
		child=node->firstdesc;
		SIBLOOP(child)
			{
			TreeToTaxaList(child, taxaList);
			}
		return;
		}
}
/***********************************************************************************/
void TreeToTaxaPtrList(NODETYPE *node,  PtrList NodeList)	
{

/* Moves through clade from node, compiling list of terminals NODES (!);
on entry taxaList must be a valid pointer to a possibly empty list */

	NODETYPE *child;
	if (isTip(node)) 
		{
		pListAddItem(NodeList, node);
		return;
		}

	else	/* interior node */
		{
		child=node->firstdesc;
		SIBLOOP(child)
			{
			TreeToTaxaPtrList(child, NodeList);
			}
		return;
		}
}
/***********************************************************************************/
void TreeToNodePtrList(NODETYPE *node,  PtrList NodeList)	
{

/* Moves through clade from node, compiling list of all nodes;
on entry nodeList must be a valid pointer to a possibly empty list */

	NODETYPE *child;
	pListAddItem(NodeList, node);
	child=node->firstdesc;
	SIBLOOP(child)
		{
		TreeToNodePtrList(child, NodeList);
		}
	return;
}
/***********************************************************************************/
void ABCSuperTree(NODETYPE *node, StrListPtr UniqueList, 
		    char **dataMatrix,float *wtset)	
{

/* on entry external gColumn MUST point to a valid column in dataMatrix, zero-offset;
i.e.,  start the thang at ZERO
'UniqueList' is the list of all the taxa in the analysis
 */

	extern int gColumn; /* declared in ReadNexusFile2 */
	NODETYPE *child;
	int ll, mm, j, numTaxa;
	char * taxon;
	StrListPtr aTaxaList;
	if (isTip(node)) 
		return;
	else	/* interior node */
		{
		if (!isRoot(node)) /* Don't add a char for whole tree */
		    {
		    wtset[gColumn]=node->length;
		    aTaxaList=newStrList();
		    TreeToTaxaList(node, aTaxaList);
		    numTaxa=lengthList(aTaxaList);
		    for(ll=1;ll<=numTaxa;ll++)
			    {
			    taxon=getkthStr(aTaxaList, ll);
			    mm=findMatchStr(UniqueList, taxon);
			    if (mm)
				dataMatrix[mm-1][gColumn]='1';
			    }
		    freeStrList(aTaxaList);
		    ++gColumn;
		    }
		child=node->firstdesc;
		SIBLOOP(child)
			{
			ABCSuperTree(child, UniqueList, dataMatrix,wtset);
			}
		return;
		}
}
/***********************************************************************************/
void ABCSuperTreePurvis(NODETYPE *node, StrListPtr UniqueList, 
		    char **dataMatrix,float *wtset)	
{
/* Sets up a datamatrix coded according to Purvis' supertree method.
 * Visits each internal node.  Constructs a set of characters for that node consisting of
 * 1's for a subclade and 0's for the other taxa in the clade,  not in that subclade.
 * Does a character for all the immediate subclades of that node,  then recurses in.
 * 
 * 
 */


/* on entry external gColumn MUST point to a valid column in dataMatrix, zero-offset;
i.e.,  start the thang at ZERO */

	extern int gColumn; /* declared in ReadNexusFile2 */
	NODETYPE *child;
	int ll, mm, j, numTaxa, numTaxa2;
	char * taxon;
	StrListPtr aTaxaList, aTaxaList2;
	if (isTip(node)) 
		return;
	else	/* interior node */
		{
		aTaxaList=newStrList();
		TreeToTaxaList(node, aTaxaList);
		numTaxa=lengthList(aTaxaList);
		child=node->firstdesc;
		SIBLOOP(child)
			{
			if (!isTip(child)) 
			    {
			    aTaxaList2=newStrList();
			    TreeToTaxaList(child,  aTaxaList2);
			    numTaxa2=lengthList(aTaxaList2);
			    for(ll=1;ll<=numTaxa;ll++)
				{
				 taxon=getkthStr(aTaxaList, ll);
				 mm=findMatchStr(UniqueList, taxon);
				 if (mm)
				    dataMatrix[mm-1][gColumn]='0';
				 }
			    for(ll=1;ll<=numTaxa2;ll++)
				{
				 taxon=getkthStr(aTaxaList2, ll);
				 mm=findMatchStr(UniqueList, taxon);
				 if (mm)
				    dataMatrix[mm-1][gColumn]='1';
				 }
				 
			
			    ++gColumn;
			    }
			
			ABCSuperTreePurvis(child, UniqueList, dataMatrix,wtset);
			}
		freeStrList(aTaxaList);
		return;
		}
}
/***********************************************************************************/
PtrList Tree2CladeSet(TREE thisTree, StrListPtr allTaxaList)
    {
    PtrList CladeSetList;
    CladeSetList=pNewList(); 
    Tree2CladeSets(thisTree->root, allTaxaList, thisTree->numTaxa, CladeSetList);
    return CladeSetList;
    }
void Tree2CladeSets(NODETYPE *node, StrListPtr allTaxaList, int nTaxa, 
		    PtrList SetList)	
{

/*  Recurses through a tree, obtaining a list of all the clades in the tree.
 *  A generic list,  'SetList',  of pointers to set vectors is repeatedly
 *  added to as we traverse the tree. 'allTaxaList' is a string list containing
 *  all the taxon names.  Sets are represented as integer (binary) vectors of
 *  size 'nTaxa'.  Membership in a clade for some taxon is signified by a 1 in 
 *  that position.  Note that 'allTaxaList' HAS TO BE CREATED ONCE and then used
 *  for all trees,  otherwise each set ordering will be unique and sets won't
 *  make sense.
 * 
 */


	NODETYPE *child;
	int ll, mm, j, numTaxa, i;
	char * taxon;
	Set cladeSet;
	StrListPtr aTaxaList;
	if (isTip(node)) 
		return;
	else	/* interior node */
		{
		if (!isRoot(node)) /* Don't add a char for whole tree */
		    {
		    cladeSet=newSet(nTaxa);
		    aTaxaList=newStrList();
		    TreeToTaxaList(node, aTaxaList); /* get taxa in this clade */
		    numTaxa=lengthList(aTaxaList); /* how many? */
		    for(ll=1;ll<=numTaxa;ll++)
			    {
			    taxon=getkthStr(aTaxaList, ll);
			    mm=findMatchStr(allTaxaList, taxon); /* find the position
					in the vector for this taxon */
			    if (mm)
				add_to_set(cladeSet, mm);
			    }
		    freeStrList(aTaxaList);
		    pListAddItem(SetList, cladeSet);
		    }
		child=node->firstdesc;
		SIBLOOP(child)
			{
			Tree2CladeSets(child, allTaxaList, nTaxa, SetList);
			}
		return;
		}
}
void printCladeSets(PtrList SetList)
    {
    int i, nTaxa;
    Set cladeSet2, cladeSet1;
    PtrList curP2, curP1;
    curP2=SetList;
    while(curP2)
	{
	    cladeSet2=(Set)(curP2->item);
	    print_set(cladeSet2);
	    curP2=curP2->next;
	}
#if 0	
    curP1=SetList;
    while(curP1)
      {
	cladeSet1=(Set)(curP1->item);
	curP2=curP1->next;
	while(curP2)
	    {
		cladeSet2=(Set)(curP2->item);
		test_set(cladeSet1, cladeSet2);
		curP2=curP2->next;
	    }
	curP1=curP1->next;
      }
#endif
    return;
    }	


void rootToTips(NODETYPE* node,double curLen)

/* Calculates distances from root to each tip and prints...*/

{
	NODETYPE *child;

	if (!isRoot(node))
		{
		curLen+=node->length;	/* don't add length under the root */
		}
	if (isTip(node)) 
			{
			printf("Root to tip:%s: %f\n",node->taxon_name,curLen);
			}
	child=node->firstdesc;
	SIBLOOP(child) {
			rootToTips(child,curLen);
			}
	return ;
}
/***********************************************************************************/


