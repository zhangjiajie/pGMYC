#include "nexus.h"
#include "TreeUtils.h"
#include "relativeRates.h"
#include "WuLi.h"
#include "MyUtilities.h"
#include "root3taxa.h"
#include "distance.h"

/* Global */


/***** COMMENTS *******

Note that missing data can cause problems in the following way:  Currently distances are calculated
pairwise.  This means that the relative rates tests can use sites that do not have data in all
three taxa.  Sometimes the distance between the two ingroups can therefore be zero but the two 
distances to the outgroups might not be equal!  When this occurs the variance calculations used
in WuLi and Steel et al. get mucked up and can try to take the square root of a negative number.

Long term solution is to give user a choice between (a) including only sites where data is present
in all taxa or (b) explaining the source of these errors.


*/


/**************************************************************************/
int doRelativeRates(StrListPtr aTaxaList, NODETYPE * root)
{
	extern FILE * outstream;
	extern struct NexDataType *gNexDataPtr;	
	int id[3], ix,jx,taxonID,kk,kind;
	char c, *dummy,*taxon,*tmp, *pi,*pj,*pRow1, *taxon1,*taxon2;
	StrListPtr WholeSelList,rootedList, unrootedList;
	struct MyRecType * dptr;
	struct NexDataType * nexPtr;	/* This is THE data structure for the NEXUS data */
	int i,j,k,NList;


	unrootedList=newStrListN(3);
	NList=lengthList(aTaxaList);
	WholeSelList=newStrListN(NList); /* for some reason I need to copy
			to a new list before writing to it; I suspect a problem
			in the 'setkthnode' line below but OK for now*/
	for (ix=1;ix<=NList;ix++) /* convert any taxon ids to taxon names 
								unless already stored that way*/
		{
		taxon=getkthStr(aTaxaList,ix);
		if(isStrInteger(taxon))
			{
			taxonID=strtod(taxon,&dummy);
			setkthNode(WholeSelList, ix, getkthStr(gNexDataPtr->TaxaList,taxonID));
			}
		else
			setkthNode(WholeSelList, ix, taxon);
		}




		
	nexPtr=gNexDataPtr;	

/* Preliminaries: do some checking to see if we can proceed and open an output file */

	if (nexPtr->isChars==0)
		{
		doGenericAlert("Characters not available in NEXUS file");
		return 0;	
		}
		
	if (nexPtr->isTaxa==0) 
		{
		doGenericAlert("Taxa not available in NEXUS file");
		return 0;	
		}
#if 1		
	printf("RELATIVE RATE TESTS: Method = ");
	switch (gNexDataPtr->RateBlockParms.RRtype)
		{
		case WULI:printf("Wu & Li\n"); break;
		case MIKE:printf("Mike's method\n"); break;
		case STEEL:printf("Steel et al.\n"); break;
		case TAJIMA:printf("Tajima\n"); break;
		}
	printf("(* = P<0.05; ** = P<0.01; *** = P<0.001)\n");
	printf("(Positive z-score means higher rate in first ingroup)\n\n");
	if (gNexDataPtr->RateBlockParms.isBS)
		printf("Bootstrap estimates of variance: N replicates = %li; Seed = %li\n",
		gNexDataPtr->RateBlockParms.NReps,
		gNexDataPtr->RateBlockParms.seed);
	printf("(Outgroup        (Ingroup1          Ingroup2     ))\t\tz (\"exact\")");

	if (gNexDataPtr->RateBlockParms.isBS)
		printf("\tz (bootstrap)\t[mean:an bs] (sdev: an bs)\n");	/* shift the column heading over if showing bs results */
	else printf("\n");
	printf("---------------------------------------------------------------------------------------------------\n");
#endif
(void)WuLiStub(7,3,2);
exit(1);
	i=1;j=2;k=3;
	
	for (i=1;i<j;i++)
		for (j=i+1;j<k;j++)
			for (k=j+1;k<=NList;k++)
				{ /* set up the unsorted 3-list */
				setkthNode(unrootedList,1,getkthStr(WholeSelList,(long)i) );
				setkthNode(unrootedList,2,getkthStr(WholeSelList,(long)j) );
				setkthNode(unrootedList,3,getkthStr(WholeSelList,(long)k) );
				rootedList=root3taxa(unrootedList,root); /* returns the properly sorted 3-list */
				/* Now convert the three taxon names in List to ID numbers and call WuLi */
				if (rootedList)  /* if null it means a polytomy or other error */
					{
					for (ix=0;ix<3;ix++)
						{
							taxon=getkthStr(rootedList,ix+1);
							jx = findMatchStr(nexPtr->TaxaList, taxon);
							if (jx ==0)
								doGenericAlert ("Matching taxon label not found in WuLi");
							id[ix]=jx;	/* make sure ids are on [1..ntaxa] */
						}
					(void)WuLiStub(id[0],id[1],id[2]);
					freeStrList(rootedList);
					}
				}




				
freeStrList(unrootedList);			
freeStrList(WholeSelList);	
return 1;
}
/**************************************************************************/
int doGroupRelRates(StrListPtr ig1List,StrListPtr ig2List,StrListPtr ogList)
{
long ig1Size,ig2Size, ogSize, i,j;

ig1Size=lengthList(ig1List);
ig2Size=lengthList(ig2List);
ogSize=lengthList(ogList);
if ((ig1Size==0) | (ig2Size==0) | (ogSize==0)) 
	{
	doGenericAlert("At least one of taxa lists is empty");
	return (0);
	}

} 
