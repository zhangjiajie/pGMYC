#if 0
#define isEqual(a,b)		(!strcmp((a),(b)))
#endif

#include "TreeUtils.h"
#include "memory.h"
#include "root3taxa.h"
#include "nexus.h"
#include "structures.h"
#include <stdlib.h>

int		gbSLcount;
StrListPtr 	gThreeList,g3SelList;
/***********************************************************************************/
/***********************************************************************************/

StrListPtr root3taxa(StrListPtr unsortedList, NODETYPE *root)	

/* 	
	If the three taxa form a polytomy, bails. Otherwise...
	Builds a list of three taxon names using my list structure.
	First and second elements are the ingroup, third is the outgroup.
	Returns the address of the list.  NULL is the error return.
	
	The next THREE functions are all necessary for this routine.
	
	Note that the algorithm for determining which taxon of the three is an outgroup is
	not intuitive.  Each taxon in order that it is encountered on a traversal is assigned 
	either 1,2,4.  When these are added up as we proceed higher in the tree, they make for
	unique pairs of sister group numbers.  Thus, if (4,(1,2)), we have 4 versus 3 and can tell
	that the single number 4 is the outgroup.  i.e., 1,2,4 are never obtained by addition of
	themselves--this is equivalent to binary coding.
	
	NOT A GOOD X-TREE ROUTINE.  IT TRAVERSES THE ENTIRE TREE EVEN IF THREE TAXA ARE CLOSELY RELATED	

	*/



{
	extern StrListPtr g3SelList;
	int total,OG,IG1,IG2;
	char* dummy;
	StrListPtr localThreeList;
	
	g3SelList=unsortedList;
/*	if (numSelected(root)==3)  */ /* no longer confine to 3 taxa in selection */
		{
		gbSLcount=1;
		gThreeList=newStrListN(3);
		localThreeList=newStrListN(3);
		if (gThreeList && localThreeList)
			total=bSLrecurse(root);  /* should work now */
		else
			{
		/*	doGenericAlert("Allocation error in String List(build list)");*/
			return NULL;			/* Avoids allocation errors */
			}
		OG=bSLrecurse2(root);
		if ((-OG >=1) && (-OG <=3))
			{
			switch (-OG)
				{
				case 1:
					IG1=2;IG2=3;OG=1; break;
				case 2:
					IG1=1;IG2=3;OG=2; break;
				case 3:
					IG1=1;IG2=2;OG=3; break;
				}
				
			setkthNode(localThreeList,1, getkthStr(gThreeList,IG1) );
			setkthNode(localThreeList,2, getkthStr(gThreeList,IG2) );
			setkthNode(localThreeList,3, getkthStr(gThreeList,OG) );			
			freeStrList(gThreeList);
			return localThreeList;
			}
		else if (OG == -2000) /* polytomy */
			return NULL;
		}
/*	else
		return NULL;*/

}

int bSLrecurse(NODETYPE *node)

/* first pass through the tree doing some stuff when it finds any of the three taxa whose
names are contained in the global list 'g3SelList'.*/
	
{
	extern struct NexDataType *gNexDataPtr;	
	char *dummy, *taxon;
	int sum=0,k,ix;
	NODETYPE *child;
	if (!node) return(-1);
	if (isTip(node) ) 
		{
		if (gNexDataPtr->isTranslate)	/* trees stored with translation table */
			{
			ix=strtod(node->taxon_name,&dummy);/* this is the taxon code */
			taxon=getkthStr(gNexDataPtr->TransList,ix);
			}
		else
			taxon=node->taxon_name;
		if 		(
		  		(isEqual(taxon,getkthStr(g3SelList,1))) ||
		  		(isEqual(taxon,getkthStr(g3SelList,2))) ||
		  		(isEqual(taxon,getkthStr(g3SelList,3)))
		  	 	)  /**** new ****/
					{
					switch (gbSLcount)
						{
						case 1: k=1;break;
						case 2:	k=2;break;
						case 3:	k=4;break;	/* these codes get assigned to the three taxa */
						}
					setkthNode(gThreeList, gbSLcount, taxon);	/*put in list at appropriate place*/
					node->numSelectedDesc=k; 
					++gbSLcount;
					return (k);
					}
	    else
			{
			node->numSelectedDesc=0; 
			return (0);
			}
		}
	else	/* interior node */
		{
		child=node->firstdesc;
		SIBLOOP(child) 
			sum+=bSLrecurse(child);
		node->numSelectedDesc=sum;
		return (sum);
		}
}

int bSLrecurse2(NODETYPE *node)	  /* this routine should not need changing */

/* second pass through the tree making decisions about which of three taxon is
the outgroup.  Returns a code to 'root3taxa' */

{
	int sum=0,lastsum,count=0;
	NODETYPE *child;
	if (!node) return(-1);
	if (isTip(node) ) 
		return (node->numSelectedDesc);

	else	/* interior node */
		{
		child=node->firstdesc;
		SIBLOOP(child)
			{
			if (child->numSelectedDesc > 0)  /* can ignore children with no selected grandchildren*/
				{
				++count;
				lastsum=sum;		/* we need to track up to two numbers */
				sum=bSLrecurse2(child);
				if (sum<0) return (sum); /* This is how we bail out of the recursion when OG IS FOUND*/
				}
			}
		if (count >= 3)
			return (-2000); /* this a a polytomy */

		switch (count) 
			{
			case 0: return (-1000); /* this means none selected--shouldn't happen */
			case 1: break;			/* do nothing; just continue looking */
			case 2: 
				if (sum+lastsum==7)
					{
					if ((sum==1) || (lastsum==1)) return (-1);
					if ((sum==2) || (lastsum==2)) return (-2);
					if ((sum==4) || (lastsum==4)) return (-3);	/* decode the codes */
					/* abs of return value will tell us which element of list is the OG! MAGIC, huh? */
					
					}
			}
		
		return (node->numSelectedDesc);
		}
}

