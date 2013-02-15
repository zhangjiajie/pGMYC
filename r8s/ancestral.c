/*

Module to reconstruct ancestral states of continuous traits via squared change parsimony.
Traits are stored in node->nodeReal

Contains several routines from NRC and some minor modifications to same...

*/
#include "ObjFunc.h"
#include "Maximize.h"
#include "structures.h"
#include "nexus.h"
#include "ancestral.h"
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define SQR(x)        ((x)*(x))

static void pTimeArray2treeAncestral(NODETYPE * node, double lp[]);
double recurseAncestral(NODETYPE *node);



double *dvector(long nl, long nh);
void free_dvector(double *v, long nl, long nh);

int gNParm, gNT, p2tindex; 

void ancestralOptimize(TREE t,int *numIter, double ftol,double linMinDelta,int *success )
{
extern struct NexDataType *gNexDataPtr;
extern NODETYPE * gRoot;
StrListPtr DM, TL;
PtrList nodeList;
int nParm,NT,i,j, ixTL;
double *p, obj,tip_value;
char * tip_name,*tip_value_str, *dummy,*found_tip;
NODE a;
float meanTip=0.0;

DM=gNexDataPtr->DMList;
TL=gNexDataPtr->TaxaList;
gRoot = t->root;
nodeList=pNewList();
TreeToTaxaPtrList(t->root,nodeList);

nParm=numIntNodes(t->root);
NT=t->numTaxa;
p = dvector(1,nParm);
for (i=1;i<=NT;i++)
	{
	a=(NODE)(pListgetkthNode(nodeList, i)->item);
	tip_name=a->taxon_name;
	ixTL=findMatchStr(TL,tip_name); // lookup the taxon name from the matrix ordering and get the relevant data matrix ordering corresponding name
	tip_value_str=getkthStr(DM,(long)(ixTL));
	found_tip=getkthStr(TL,(long)(ixTL));
	tip_value=strtod(tip_value_str,NULL);
// printf("%i\t%s\t%s\t%s\t%f\n",i,tip_name,found_tip,tip_value_str,tip_value);
	a->nodeReal=tip_value;
	meanTip+=tip_value;
	}
meanTip/=NT;
for (i=1;i<=nParm;i++)
	p[i]=meanTip; // use mean tip value as guess for all nodes...


obj=MinND(t,0, POWELL,objAncestral,NULL,p, nParm,numIter, ftol,linMinDelta,success );



free_dvector(p,1,nParm);
}




/***********************************************************************************/

double objAncestral(double p[])


{
  extern struct NexDataType *gNexDataPtr;	
  extern NODETYPE * gRoot;    /* This global is declared when the whole algorithm is called */

  static int firstTime=1,num_branches;
  double obj;
  NODETYPE *child;
  
  p2tindex=1;
  pTimeArray2treeAncestral(gRoot,p);  // puts the array values onto the tree 


/*** Now find objective function over rest of tree ***/

  obj=recurseAncestral(gRoot);

  return obj; 
}
/**********************/

double recurseAncestral(NODETYPE *node)
{
  NODETYPE *child;
  double obj=0.0;
  if (!node) return(0.0);
  child=node->firstdesc;
  SIBLOOP(child)
    {
      obj += SQR(child->nodeReal-node->nodeReal);
      obj += recurseAncestral(child);
    }
  return obj;	
}

static void pTimeArray2treeAncestral(NODETYPE * node, double lp[])
{
NODETYPE *child;
if (!isTip(node))
	node->nodeReal=lp[p2tindex++];
child=node->firstdesc;
SIBLOOP(child)
	pTimeArray2treeAncestral(child,lp);
return;
}
/***********************************************************************************/
void printAncestral(NODETYPE *node)
{
    float diff;
    NODETYPE *child;

    if (!isTip(node))
      {
      if(*(node->taxon_name))
	    printf(" [*]  %.8s\t", node->taxon_name);
      else
	    printf("     (%i)\t", node->id);
      }
    else
	printf("      %.8s\t",node->taxon_name);

//    printf("%7.2f\t\t",node->time);
    if (!isRoot(node))
	{
	diff = node->nodeReal - node->anc->nodeReal;
	printf("%.4f\t\t%.4f\t\t%.4f",node->nodeReal,node->anc->nodeReal,diff);
	if (diff >= 0.0) printf("\t\t+\n");
	else printf("\t\t-\n");
	}
     else
	printf("%.4f\n",node->nodeReal);



    if(!isTip(node))
	{
    	child=node->firstdesc;
    	SIBLOOP(child) 
		printAncestral(child);
	}
    return;
    
}
