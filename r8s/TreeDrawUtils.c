#include <limits.h>
#include "TreeUtils.h"
#include "MacApp.h"
#include "TreeDrawUtils.h"

double 	gMaxToTipLength;	/* used in phylogram draws */
int		gRootFlag;

/***********************************************************************************/

void  MakeTreeStruct(char *TreeDescriptionPtr,struct MyRecType *theDataPtr)
										 
	/*  Takes a string tree description and creates the linked list tree structure,
	a pointer to which is then stored in the global data structure 'theData' */										 
										 
{
	NODETYPE *root;
	Str255 temp;
	long i;


	if (TreeDescriptionPtr==NULL) {
		doGenericAlert("BAD tree description");
		theDataPtr->treeStructisMade=0; /* note that it has NOT been made */
		return;
		};

	

	root=string_to_tree(TreeDescriptionPtr);

		

	if(root ==NULL) {  	/* error: tree description bad */
		doGenericAlert("BAD tree description");
		theDataPtr->treeStructisMade=0; /* note that it has NOT been made */
		return;
		}
		
	theDataPtr->theRoot=root;	/* store the pointer to the tree structure */
	theDataPtr->displayRoot=theDataPtr->theRoot; /* initially drawn tree will draw from root */
	
	theDataPtr->treeStructisMade=1; /* note that it has been made */
	SetInternalNodeCompact(theDataPtr->theRoot,theDataPtr->collapseInternalNode);
	if (theDataPtr->allowCompact)
		SetCompactNodes(theDataPtr->theRoot);
	theDataPtr->QListHasChanged=1;
	
	gRootFlag=1; 	/* So we don't count the root */
	gMaxToTipLength=calcMaxToTip(root);
	if ((gMaxToTipLength< FLT_MAX) && (gMaxToTipLength>0.0) )
		theDataPtr->lengthsAvailable=1;
		/* previous code checks to see if ALL of the initial lengths in the tree have 
		been replaced by smaller and more reasonable numbers (they were all initialized
		to a huge value).  If something went wrong on reading a tree description list,
		particularly if soome branch did not have a length, the flag will stay zero*/
	return;
}  
/**********************************************************************************/
void MyDrawString(char *s, int flag)
{
int length;
length=strlen(s);
if (length>0)
	{
	if (!flag) 			/* draw normal */
		TextMode(srcCopy);
	else				/* invert */
		TextMode(notSrcCopy);
	DrawText(s,0,length);  
	TextMode(srcOr);
	}
return;
}
/***********************************************************************************/
void 	DrawIntName(NODETYPE *node,Rect *contRectPtr)
{
	int width,length;
	length=strlen(node->taxon_name);
	width=TextWidth(node->taxon_name,0,length);
	MoveTo(contRectPtr->right-kScrollbarAdjust-width-3,contRectPtr->bottom-3);
	DrawText(node->taxon_name,0,length);
/*	MyDrawString(node->taxon_name,0);
		 For now (!) always write this as normal uninverted text */
	TextFace(FONTSTYLE);
	return;
}
/***********************************************************************************/

void Assign_XY_Tree(NODETYPE *root,  Rect *TreeRectPtr, int treeMode)
						  
/*	Assigns X and Y display coordinates to the nodes in a tree structure.  Uses
	the  integer coordinates of the upper left and lower right of the
	box in which the tree should be displayed.  The X,Y values are stored in the
	tree structure.
	treeMode is the type of tree: 0 = cladogram, 1 = phylogram , 2= chronogram*/


{
		extern double 	gMaxToTipLength;
		int N;
		double yinc,yUpLeft;  /*have to be double for accurate placement of lines */
		if(root==NULL) 
						return;
			
		(void)maxorder(root);
		(void)numdesc(root);
		gRootFlag=1;		/* needed in recursive 'calcMaxToTip'  */
		if (treeMode == 1)
			gMaxToTipLength=calcMaxToTip(root); /*for phylogram find longest path to tip*/

		gRootFlag=1;		/* needed in recursive 'assignX'  */
		assignX(root,TreeRectPtr->left,TreeRectPtr->right,
					TreeRectPtr->right-TreeRectPtr->left+1,treeMode);
		N=root->numdesc;
		if (N==1) yinc=0;
		else yinc = (TreeRectPtr->bottom-TreeRectPtr->top)/(double)(N-1);
		yUpLeft=TreeRectPtr->top;
		assignY2(root,&yUpLeft,yinc);						
		return;
}

/***********************************************************************************/

double calcMaxToTip(NODETYPE* node)

/* Calculates maximum distance from root to tip when lengths are available */

{
	double max=0.0,temp,thisLength;
	NODETYPE *child;
	if (!node) return(0.0);

	if (gRootFlag == 1)
		{
		gRootFlag = 0;
		thisLength=0.0;
		}
	else
		{
		thisLength=node->length;	/* don't add length under the root */
		}
	if (isTip(node)  ||  (node->isCompactNode)  ) 
			{
			return (thisLength);  /* treat a tip and a compact node the same way */
			}
	child=node->firstdesc;
	SIBLOOP(child) {
			temp=calcMaxToTip(child);
			if (temp > max) max = temp;
			}
	return thisLength+max;
}
/***********************************************************************************/

int assignY2(NODETYPE *node,  double *YcurPtr,  double yinc)
{
	NODETYPE *child;
	int sum=0,count=0;
	
	if (!node) return;
	if(isTip(node)  || (node->isCompactNode) )
		{
		node->Y= *YcurPtr;
		(*YcurPtr)+=yinc;
		return(node->Y);
		}

	child=node->firstdesc;
	
	SIBLOOP(child) {
		sum+=assignY2(child,YcurPtr,yinc);
		++count;
		}
	sum/=count;
	node->Y=sum;
	return(sum);
}
/***********************************************************************************/

void assignX(NODETYPE *node,  int Xleft,  int Xright, int Xwidth, int treeMode)
{
	extern double 	gMaxToTipLength;
	if (gRootFlag) {
			node->X = Xleft;
			gRootFlag = 0;
			}
	else
			switch (treeMode)
			{
			case 0:
				node->X = Xleft + (Xright - Xleft-1)/(float)(node->order + 1);
				break;
			case 1:
				node->X = Xleft + (Xwidth-1)*node->length/gMaxToTipLength;
				break;
			case 2:
				node->X = Xleft + (Xwidth-1)*(1-node->time);
				break;
			
			}


	if (node->sib) assignX(node->sib,Xleft,Xright,Xwidth,treeMode);
	
	if (node->isCompactNode) return;  

	switch (treeMode)
		{
		case 2:
		if (node->firstdesc) 
			assignX(node->firstdesc,Xleft,Xright,Xwidth,treeMode);
		break;
		default:
		if (node->firstdesc) 
			assignX(node->firstdesc,node->X,Xright,Xwidth,treeMode);
			
		}
	return;
}
/**********************************************************************************/
void DrawHigherName(struct MyRecType *DataPtr, Rect *locContentRectPtr)
{
char *r="Root";
int length;
NODETYPE *node;
node=DataPtr->displayRoot;
if (node==DataPtr->theRoot)
	{
		MoveTo(locContentRectPtr->left+2,locContentRectPtr->bottom-3);
		MyDrawString(r,0);	
	}
else
	if (*(node->taxon_name))
		{
		MoveTo(locContentRectPtr->left+2,locContentRectPtr->bottom-3);
		MyDrawString(node->taxon_name,0);
		}
return;
}

/**********************************************************************************/
void DrawTree(WindowPtr theWindow)

										 
	/*  Plots the branches
	of it in a box in a Mac Window,the current GrafPort.  Does not yet free
	up the allocated space for the tree structure.  Nor does it write the taxon
	names!*/										 
										 
{
 	int Top,Left,Bottom,Right; 
	int x,y, windowWidth,treeAreaWidth;
	Str255 temp;
	Rect TreeRect;
	Rect *DrawRectPtr;
	int treeMode,offsetRuler,j;
	NODETYPE *root;
	extern int gMax;
	double maxLength;
	struct MyRecType * dptr;
	
	dptr=(struct MyRecType *)GetWRefCon(theWindow); /* window's data */
	root=dptr->displayRoot;
	treeMode=dptr->treeMode;
/*	ContToDrawRect(&theWindow->portRect,DrawRectPtr);*/


	switch (treeMode)	/* Checks to see if necessary info is availabe for this tree view type*/
		{
		case 0:			/* Can always show cladogram */
			break;		
		case 1:			/* phylogram*/
			if (!dptr->lengthsAvailable)
				{
				doGenericAlert("Branch lengths are not currently available");
				dptr->treeMode=0;	/* Restore default */
				return;
				}
			break;
		case 2:			/* chronogram */
			if (!dptr->timesAvailable)
				{
				doGenericAlert("Branch times are not currently available");
				dptr->treeMode=0;
				return;
				}
			break;
		
		}



	gMax=0;
	
	if (treeMode == 0) 
		offsetRuler=0;
	else
		offsetRuler=dOFFSETRULER;



	(dptr->TreeRectPtr)->top=(theWindow->portRect).top + BORDER + iTextHalfHt;
	(dptr->TreeRectPtr)->left=(theWindow->portRect).left +BORDER;
	(dptr->TreeRectPtr)->bottom=(theWindow->portRect).bottom -BORDER - iTextHalfHt - offsetRuler-kScrollbarAdjust;
	(dptr->TreeRectPtr)->right=(theWindow->portRect).right -BORDER-kScrollbarAdjust;
	
	windowWidth=dptr->TreeRectPtr->right-dptr->TreeRectPtr->left+1; /* character or pixel width of window */
				
	if (root->isCompactNode) {
		root->isCompactNode=0;
			MaxTaxLength(root);
						 /* put the width of the longest taxon name in gMax */
			treeAreaWidth=max(MINWIDTH,windowWidth-gMax-6) ;
						/* tree area is the larger of MINWIDTH and the specified
						window minus the taxon areas;
						the '-6' gives space for the space between tip and name */
			gMax=min(gMax+6,windowWidth-MINWIDTH);
					/* reset gMax to be the display width ALLOWED by the difference between
					the actual size of the window and the minimum size of the tree;  this gets
					used to determine how much of the possibly length taxon name gets displayed*/		
		dptr->TreeRectPtr->right=dptr->TreeRectPtr->left+treeAreaWidth-1;
		Assign_XY_Tree(root,dptr->TreeRectPtr,treeMode);	
		MacDrawTree(root); 
		root->isCompactNode=1;
		}
	else
		{
		MaxTaxLength(root);
		treeAreaWidth=max(MINWIDTH,windowWidth-gMax-6) ;
		gMax=min(gMax+6,windowWidth-MINWIDTH);
		dptr->TreeRectPtr->right=dptr->TreeRectPtr->left+treeAreaWidth-1;
		Assign_XY_Tree(root,dptr->TreeRectPtr,treeMode);	
		MacDrawTree(root);
		} 
		/* all this is to allow drawing of the inside of a compact node, which would
		be precluded by the fact that the root of that subtree has its flag set to compact,
		and this would force AssignXY to bail out immediately */

	if (treeMode > 0)
		{
		gRootFlag=1; 	/* So we don't count the root */
		maxLength=calcMaxToTip(root);
		drawRuler(*(dptr->TreeRectPtr),treeAreaWidth,treeMode,maxLength);
		}

	return;


}  
/***********************************************************************************/
#define dINSET_RULER 15
#define dTICKMARK_LENGTH 3
void drawRuler(Rect TreeRect,int width, int treeMode, double maxLength)
{
double length=1.0,interval;	/* factor controls roughly how long the phylogram bar is */
Str255 numStr;
int	leftRuler, rightRuler,yRuler,i,x;
switch (treeMode)
	{
	case 1:
		if (maxLength > 1.0)
			{
			while ( length< 0.1* maxLength)
				length *= 10.0;
			}
		x=length;
		length=length*width/maxLength;
		leftRuler=TreeRect.left;
		rightRuler=TreeRect.left+length;
		yRuler=TreeRect.bottom+dINSET_RULER;
		NumToString((long)x,numStr);
		MoveTo(leftRuler,yRuler);
		LineTo(rightRuler,yRuler);
		MoveTo(rightRuler + 2, yRuler);
		DrawString(numStr);
		interval=length/10.;
		for (i=0;i<=10;i++)
			{
			x=leftRuler+i*interval;
			MoveTo(x,yRuler); 
			LineTo(x, yRuler-dTICKMARK_LENGTH);
			}
		
		
		break;
	case 2:
		leftRuler=TreeRect.left;
		rightRuler=TreeRect.right;
		yRuler=TreeRect.bottom+dINSET_RULER;
		MoveTo(leftRuler,yRuler);
		LineTo(rightRuler,yRuler);
		interval=width/10.;
		for (i=0;i<=10;i++)
			{
			x=leftRuler+i*interval;
			MoveTo(x,yRuler); 
			LineTo(x, yRuler-dTICKMARK_LENGTH);
			}
		MoveTo(leftRuler-2, yRuler+FONTSIZE);
		DrawString("\p1.0");
		MoveTo(rightRuler-2, yRuler+FONTSIZE);
		DrawString("\p0.0");


		break;
	}

}
/***********************************************************************************/
void MacDrawTree(NODETYPE *node)
{
	NODETYPE *child;
	Rect theRect;
	int length,invertFlag;
	char *asterisk="*";
	char s[20];
	if (!node) return;
	if (!isTip(node)  ) 
		{
			SetRect(&theRect,node->X-NODEBOX,node->Y-NODEBOX,node->X+NODEBOX,node->Y+NODEBOX);
			FrameRect(&theRect);
			PaintRect(&theRect);
			if (node->nodeIsConstrained)
				{
				MoveTo(node->X+4,node->Y+3);
				LineTo(node->X+4,node->Y-3);
				MoveTo(node->X-5,node->Y+3);
				LineTo(node->X-5,node->Y-3);
				}
			MoveTo(node->X+3,node->Y+5);
			if  (*(node->taxon_name))  
					{
			/*		MyDrawString(asterisk,node->isQueryNode);*/ /* asterisk */
					MyDrawString(node->taxon_name,node->isQueryNode);
					
					TextFace(FONTSTYLE);
					}
		}
	if (!(node->isCompactNode) )  /* only draw branches to descendants if node is not compact */
		{
			child=node->firstdesc;
			SIBLOOP(child) 
				{
				MoveTo(node->X,child->Y);
				LineTo(child->X,child->Y);
				
				Move(TNOFFSETX,TNOFFSETY);  /*  offset from tip of branch*/
				if (isTip(child))
					{
					length=strlen(child->taxon_name);
					if (TextWidth(child->taxon_name,0,length) <= gMax)
						MyDrawString(child->taxon_name,child->isQueryNode);  /* as long as name isnt too 
																	long...*/
					TextFace(FONTSTYLE);
					}	
				MoveTo(node->X,node->Y);
				LineTo(node->X,child->Y);
				MacDrawTree(child);
				}
		}
	return;
}
/***********************************************************************************/
int maxorder(NODETYPE *node)
{
	 int max,temp;
	NODETYPE *child;
	if (!node) return(-1);
	if (isTip(node)  ||  (node->isCompactNode)  ) {node->order=0; return (0);}
			/* treat a tip and a compact node the same way */
	max=0;
	child=node->firstdesc;
	SIBLOOP(child) {
			temp=maxorder(child);
			if (temp > max) max = temp;
			}
	node->order=max+1;
	return (max+1);
}
/****************************************************/
void MaxTaxLength(NODETYPE *node)

	/* Finds the length of the longest string contained in the tree structure */
{
	NODETYPE *child;
	extern int gMax;
	int temp, length, max=0;
	if (!node) return;

	if (isTip(node))
		{
		length=strlen(node->taxon_name);
		temp=(short)TextWidth(node->taxon_name,0,(short)length);
		if (temp>gMax) gMax=temp;
		}


	if (node->isCompactNode) return;
	child=node->firstdesc;
	SIBLOOP(child) 
		{
		MaxTaxLength(child);
		}
	return;


}
/****************************************************/
NODETYPE *SearchTreeNodes(Point localPt, int radius, NODETYPE *root)

	/* Checks to see if local mouse-coordinates in 'P' are within some distance 'radius'
	from a node in the window's tree structure.  If so the function returns a pointer to
	the node; if not it returns NULL.  */
{
	NODETYPE	*foundNode;

	foundNode=TraverseScanPt(localPt, radius, root);
	return(foundNode);
}
/****************************************************/
void SetCompactNodes(NODETYPE *rootPtr)
{
extern int kMAXSCRTAXA;
int ix=0;
do	{
	numdesc(rootPtr);
	SetOneNode(rootPtr);
	++ix;
	}
	while (rootPtr->numdesc > kMAXSCRTAXA && ix< kMAXSCRTAXA);
return;
}

int SetOneNode(NODETYPE *node)
{
if (node->isCompactNode) return (0);
if ( (node->sib) && SetOneNode(node->sib) ) return (1);
if ( (node->firstdesc) && SetOneNode(node->firstdesc) ) return (1);
if (node->numdesc >kMAXSCRTAXA)
	{
	node->isCompactNode=1;
	return (1);
	}
else
	return (0);
}

void ClearCompactNodes(NODETYPE *node)
{
if (!node) return;
node->isCompactNode=0;
if (node->sib) ClearCompactNodes(node->sib);
if (node->firstdesc) ClearCompactNodes(node->firstdesc);
return;
}

void SetInternalNodeCompact(NODETYPE *node,int mode)
{
/* if the collapse internal node mode flag is set, traverse tree and set to compact
	node any node that is (a) internal, and (b) has a taxon name associated */

if (!node) return;
if 	( 
	(*(node->taxon_name))  
	&& (!isTip(node)) 
	&& mode
	)
	node->isCompactNode=1;
if (node->sib) SetInternalNodeCompact(node->sib,mode);
if (node->firstdesc) SetInternalNodeCompact(node->firstdesc,mode);
return;
}



/***********************************************************************************/

NODETYPE *TraverseScanPt(Point localPt, int radius, NODETYPE *node)

{
	NODETYPE	*foundNode;
	
	if (!node) return (NULL);
	if (  (node->X - radius) < localPt.h && (node->X + radius) > localPt.h &&
	    	  (node->Y - radius) < localPt.v && (node->Y + radius) > localPt.v ) return(node);
	else 
		{
			if (node->sib) 
				{
				foundNode=TraverseScanPt(localPt,radius, node->sib);
				if(foundNode) return (foundNode);
				}
			if (node->firstdesc) 
				{
				foundNode=TraverseScanPt(localPt,radius,node->firstdesc);
				if(foundNode) return (foundNode);
				}
			return (NULL);
		}
}
/***********************************************************************************/

NODETYPE *QTreeNodes(Point localPt, int radius, NODETYPE *node)

/* scans through the tree and returns the node corresponding to either a node that is clicked
or the node next to the taxon name that is clicked */

{
	NODETYPE *child;
	NODETYPE	*foundNode;
	Rect		taxRect;
	int			length,width,height,starth,startv;
	char		*asterisk="*";

/* check for click in node box */
	if (!node  || node->isCompactNode) return (NULL);
	if (  (node->X - radius) < localPt.h && (node->X + radius) > localPt.h &&
	    	  (node->Y - radius) < localPt.v && (node->Y + radius) > localPt.v ) return(node);

/* check for click in taxon name area, only if there is a taxon name for that node
	and ONLY for terminal taxa */
	
	if (*(node->taxon_name)  && (isTip(node)))
		{
		width=TextWidth(node->taxon_name,0,strlen(node->taxon_name));
		height=FONTSIZE/2;
		SetRect(&taxRect,node->X +4, node->Y-height,node->X +4+width, node->Y+height);
		if (PtInRect(localPt,&taxRect)) 
			{
			return (node);
		/*	doQToggleTaxon(node);
			return((NODETYPE *)-1);	*/	/* indicates a taxon was toggled */
			}
		}

/* check for click in asterisk that corresponds to higher node name, but only
for internal node that has a name */

	starth=node->X+3;
	startv=node->Y+5;
	if (*(node->taxon_name)  && (!isTip(node)))
		if (  (node->X +3) < localPt.h && (node->X +9) > localPt.h &&
	    	  (node->Y - 3) < localPt.v && (node->Y + 3) > localPt.v ) 
	    	 {
			 FLAGFLIP(node->isQueryNode);
	    	 MoveTo(starth,startv);
			 MyDrawString(asterisk,node->isQueryNode); /* asterisk */
			 return((NODETYPE *)-1);	/* just toggle the status and let 'DrawIntName handle it */
	    	 }

 
	child=node->firstdesc;
	SIBLOOP(child) 
			{
			foundNode=QTreeNodes(localPt,radius, child);
			if(foundNode) return (foundNode);
			}
	return (NULL);	/* Haven't found either a node or a taxon name so far */

}
void doQToggleTaxon(NODETYPE *node)
{
FontInfo	fi;
int 		height,width,starth,startv;
Rect		taxRect;

/*GetFontInfo( &fi );*/

starth=node->X+TNOFFSETX;
startv=node->Y+TNOFFSETY;
/*width=TextWidth(node->taxon_name,0,length);
SetRect(&taxRect,starth, startv-fi.ascent,starth+width, startv+fi.descent);
EraseRect(&taxRect);*/

MoveTo( starth,startv );

FLAGFLIP(node->isQueryNode);
MyDrawString(node->taxon_name,node->isQueryNode);  
TextFace(FONTSTYLE);	

return;

}
void QToggleClade(NODETYPE *node)

/* sets or resets the query status of all of the tips descended from node.
DOES NOT alter the internal nodes.  Decision on whether to set or reset
is based on node->toggleDesc flag. Thus, clicking on an internal node will
either set all descendants or reset all descendants depending on the history of clicks
on that internal node*/
{
	int flag;
	if (node->toggleDesc) 
		flag=0; 
	else 
		flag=1;
	QSetAllNodes(node,flag);
	node->toggleDesc=flag;
	return;
}
void QSetAllNodes(NODETYPE *node,int flag)
{

	NODETYPE *child;
	if (!node) return;
	if (isTip(node))
		node->isQueryNode=flag;		
	child=node->firstdesc;
	SIBLOOP(child) {
			QSetAllNodes(child,flag);
			}
	return;
}

/*********************************************************/
NODETYPE *ScanInternNodes(Point localPt, int radius, NODETYPE *node)


{
	NODETYPE	*foundNode;
	int			length,width,height;

/* check for click in asterisk */
	if (!node) return (NULL);
	if (!isTip(node))
		if (  (node->X +6- radius) < localPt.h && (node->X +6+ radius) > localPt.h &&
	    	  (node->Y - radius) < localPt.v && (node->Y + radius) > localPt.v ) 
	    	 {
	    	 return(node);
	    	 }



 
	if (node->sib) 
		{
		foundNode=ScanInternNodes(localPt,radius, node->sib);
		if(foundNode) return (foundNode);
		}
	if (node->firstdesc) 
		{
		foundNode=ScanInternNodes(localPt,radius,node->firstdesc);
		if(foundNode) return (foundNode);
		}
	return (NULL);

}

