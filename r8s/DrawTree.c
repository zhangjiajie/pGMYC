#include <string.h>
#include <stdio.h>
#include <ctype.h>
#include <stdlib.h>
#include "TreeUtils.h"
#include "DrawTree.h"
#include "memory.h"
#include <math.h>
#include "MyUtilities.h"

/* private functions */

void assignX(NODETYPE *node,  int Xleft,  int Xright, int Xwidth, int treeMode);
void Assign_XY_Tree(NODETYPE *root,  int xUpLeft,int yUpLeft, int xLowRight,
			int yLowRight,int treemode);
void  MakeTree(NODETYPE* root,	 int xUpLeft, int yUpLeft, int xLowRight, 
		int yLowRight, int treeWidth, int nameWidth, int treemode);
int 		assignY2(NODETYPE *node,  double *YcurPtr,  double yinc);
void		Tprint(NODETYPE *r),
		HorizLine( int x1,  int x2,  int y),
		VertLine( int y1,  int y2,  int x),
		DrawIt(NODETYPE *node),
		DrawNames(NODETYPE *node, int treemode),
		swap( int *x1,  int *x2);
void		MaxTaxLength(NODETYPE *node);
char*		TreeToString(void);

/* GLOBALS */

static double gAge;
int gHyp,gMax=0;
double gMaxToTipLength, gMaxToTipLengthRate;
char matrx[MAXHEIGHT][MAXWIDTH];

int		WINWIDTH,WINHEIGHT;

NODETYPE*	internalList[MAXINTERNALNODENAMES];	/* fixed linear array of pointers to nodes */
int		internalListix=0;
/***********************************************************************************/

void DrawTree(NODETYPE *root, int treemode, int userMaxWidth)


/*

	treemode	type of tree
	   0		   cladogram
	   1		   phylogram
	   2 		   chronogram
	   4		   ratogram
	   9		   cladogram with taxon names replaced with character state in column 1!
	   10		   phylogram with taxon names replaced with character state in column 1!
*/

{


	int numtax,treeWidth,nameWidth;

	numtax=numdesc(root);
	gAge=root->time; /* used in chronogram draws */
	if ((gAge==0.0) && (treemode==2))
	    fatal("Times do not appear to have been set on trees; try converting lengths to time");
	gMax=0;		/* width of longest taxon name */
	MaxTaxLength(root);	  /* put the width of the longest taxon name in gMax */
	gMax+=2; 	/* Add two characters for possible special characters '%' */
	gMax+=6; 	/* Add 6 characters for possible terminal character states */
	nameWidth=gMax;
	gMaxToTipLength=calcMaxToTip(root);
	gMaxToTipLengthRate=calcMaxToTipRate(root);
if (userMaxWidth !=0)
	treeWidth=userMaxWidth;
else
	treeWidth=max(max(2*numtax,MINWIDTH),userMaxWidth);
	WINWIDTH=treeWidth+nameWidth; /* window is big enough to a name area of size gMax and a tree area that
						might be as small as MINWIDTH but is larger for big trees */ 


	if (WINWIDTH > MAXWIDTH)
		{
		WINWIDTH = MAXWIDTH - nameWidth -1;
		treeWidth = WINWIDTH - nameWidth;
		printf("Tree is so large that it has been compressed horizontally and resolution may be lost: WINWIDTH=%i treeWidth=%i\n",WINWIDTH,treeWidth);
		}
	WINHEIGHT=min(2*numtax-1,MAXHEIGHT);

	if (treemode == 9 || treemode == 10) printf("Currently printing marginal probs that switch is ON at each node\n");
	MakeTree(root,0,0,WINWIDTH-1,WINHEIGHT-1,treeWidth,nameWidth,treemode);

}
/***********************************************************************************/

double calcMaxToTip(NODETYPE* node)

/* Calculates maximum distance from root to tip when lengths are available */

{
	double max=0.0,temp,thisLength;
	NODETYPE *child;
	if (!node) return(0.0);

	if (isRoot(node))
		{
		thisLength=0.0;
		}
	else
		{
		thisLength=node->length;	/* don't add length under the root */
		}
	if (isTip(node)) 
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

double calcMaxToTipRate(NODETYPE* node)

/* Calculates maximum distance from root to tip when lengths are available */

{
	double max=0.0,temp,thisLength;
	NODETYPE *child;
	if (!node) return(0.0);

	if (isRoot(node))
		{
		thisLength=0.0;
		}
	else
		{
		thisLength=node->estRate;	/* don't add length under the root */
		}
	if (isTip(node)) 
			{
			return (thisLength);  /* treat a tip and a compact node the same way */
			}
	child=node->firstdesc;
	SIBLOOP(child) {
			temp=calcMaxToTipRate(child);
			if (temp > max) max = temp;
			}
	return thisLength+max;
}
/************************************************************************************/


void  MakeTree(NODETYPE* root,	 int xUpLeft, int yUpLeft, int xLowRight, 
		int yLowRight, int treeWidth, int nameWidth, int treemode)
										 
	/*  Takes a string tree description and plots the branches
	of it in a box-- NOW ONLY USED FOR DUMB TERMINAL DRAWING */										 
										 
{
  
	NODETYPE *ax,*bx;
	int x,y, windowWidth,treeAreaWidth;
	char* s;	


	windowWidth=xLowRight-xUpLeft+1; /* character or pixel width of window */
	Assign_XY_Tree(root,xUpLeft,yUpLeft,xUpLeft+treeWidth-1,yLowRight,treemode);	
	for (y=0;y<MAXHEIGHT;y++)
		for (x=0;x<MAXWIDTH;x++)  matrx[y][x]=SPACE;

	DrawIt(root);
	DrawNames(root,treemode);
	s=TreeToString();
	printf("%s\n",s);
	myFree(s);
	return;


}  
/***********************************************************************************/
char* TreeToString(void)
{
	char *s;
	int width[MAXHEIGHT],
		TotChars,
		NChars,
		row,
		col;
	char *ptr;
	
	TotChars=0;
	for (row=0;row<=WINHEIGHT-1;row++)
		{
		NChars=WINWIDTH;
		for (col=WINWIDTH-1;col>=0;col--,NChars--)
			{
			if (matrx[row][col] != SPACE) break;
			}
		width[row]=NChars;
		TotChars+=(NChars+1); /* the extra '1' is for the CR added on each line
								4D WANTS TO SEE CR ONLY, NOT CR-LF */
		}

	s=(char*)myMalloc(sizeof(char)*(TotChars+1));
	if (!s) return NULL;


	ptr=s;
	
	for (row=0;row<=WINHEIGHT-1;row++)
		{
		for (col=0;col<=width[row]-1;col++)
			{
			*ptr=matrx[row][col];
		/*	if ((*ptr == 10) || (*ptr == 13)) don't store extraneous
						LF's or CR's that may have
						crept in via makegroup  
				*ptr=SPACE;  instead replace with space --FIXED */
			++ptr;
			}
	/*	*ptr=RETURN;
		++ptr;*/
		*ptr=LF;
		++ptr;	/* UNIX and WWW recognizes LF as CR-LF */
		}

	*ptr=0;

return (s);
}


/***********************************************************************************/
void Assign_XY_Tree(NODETYPE *root,  int xUpLeft,int yUpLeft, int xLowRight,
						 int yLowRight, int treemode	)
						  
/*	Assigns X and Y display coordinates to the nodes in a tree structure.  Uses
	the  integer coordinates of the upper left and lower right of the
	box in which the tree should be displayed.  The X,Y values are stored in the
	tree structure */


{
		 int N,yinc;
		 double yupleft;
		 yupleft=yUpLeft;
		if(root==NULL) 
						return;
			
		(void)maxorder(root);
		(void)numdesc(root);
		assignX(root,xUpLeft,xLowRight,xLowRight-xUpLeft+1,treemode);
		N=root->numdesc;
		if (N==1) yinc=0;
		else yinc = (yLowRight-yUpLeft)/(float)(N-1);
		assignY2(root,&yupleft,yinc);						
		return;
}
/***********************************************************************************/

void MaxTaxLength(NODETYPE *node)

	/* Finds the length of the longest string contained in the tree structure */
{
	NODETYPE *child;
	extern int gMax;
	int temp, length, max=0;
	if (!node) return;
	child=node->firstdesc;
	SIBLOOP(child) 
		{
		if (*(child->taxon_name) == '\0')
			length=4; /* use max likely number of digits in the id
					instead of the absent string */
		else
			length=strlen(child->taxon_name);
		if (length>gMax) 
			gMax=length;
		MaxTaxLength(child);
		}
	return;


}
/***********************************************************************************/
int assignY2(NODETYPE *node,  double *YcurPtr,  double yinc)
{
	NODETYPE *child;
	int sum=0,count=0;
	
	if (!node) return 0;
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
	if (isRoot(node)) {
			node->X = Xleft;
			}
	else
			switch (treeMode)
			{
			case 9:
			case 0:
				node->X = Xleft + (Xright - Xleft-1)/(float)(node->order + 1);
				break;
			case 10:
			case 1:
				node->X = Xleft + (Xwidth-1)*node->length/gMaxToTipLength;
				break;
			case 2:
				node->X = Xleft + (Xwidth-1)*(gAge-node->time)/gAge;
				break;
			
			case 4:
				node->X = Xleft + (Xwidth-1)*node->estRate/gMaxToTipLengthRate;
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
/***********************************************************************************/

void DrawIt(NODETYPE *node)
{
	NODETYPE *child;
	 int x,count=0;
	 char *name;
	if (!node) return;

/* 'draw' the taxon name for the node */

	child=node->firstdesc;
	SIBLOOP(child) {  /* if node is already a tip, this will just fall through */
		HorizLine(node->X,child->X,child->Y);
			/* horiz line from the node's X to its child */
		VertLine(node->Y,child->Y,node->X);
			/* vert line from the node's Y to meet the previous horiz  line */
		DrawIt(child);
		}
	return;
}

/*********************************************************************/
#define ccheck(c) (isalpha(c) || ((c)=='%'))
void DrawNames(NODETYPE *node, int treemode)
{
	extern int gLabel;
	NODETYPE *child;
	 int x,room,slength,offset,count=0,code=0,width;
	 char *name, *nameFromShell, s[15],ss[20],c;
	 double sum, switch_on;
	if (!node) return;

/* 'draw' the taxon name for the node */
	if (treemode==9 || treemode==10)  // this is the oddball case of printing the character reconstruction instead of name
		{
		if (isTip(node))
			{
			if (  (node->CL)[2]	== 1.0) c='+';
			else						c='-';
//			sprintf(s, "%*i (%c)", 1,node->opt,c);
//			name=s;
			snprintf(ss, 20,"%*i (%c) %.10s", 1,node->opt,node->state,node->taxon_name);
// had no end of trouble with abort traps until I switched over snprintf!
			name=ss;
			
			}
		else
			{
			sum = (node->CLmarg)[0] + (node->CLmarg)[1] + (node->CLmarg)[2];
			switch_on =     1 - (node->CLmarg)[0]/sum ;
			snprintf(s, 15, "%1i %4.2f (%li)", node->opt, switch_on, node->id);
			name=s;
			}
		}
	else
		name=node->taxon_name;
	if (gLabel)
	 if (!(*name)) /* if an internal node name doesn't exist ...*/
	    {
	    name=s;  /*...then set up a temp string using fixed array */
	    if(node->id == 0)
		width=1;
	    else
	        width=log10(node->id)+1;/* get the number of digits in the id */
	    sprintf(name, "%*li", width,node->id); /* just use its id number */
	    }
	slength=strlen(name)+2; /* this is how much room we need */
	x=(node->X)+1;  /* start writing  one character to the right of x */
	while(*name)   /* loop through the string */
		  {
			if (count>gMax) break;
		  	matrx[node->Y][x]=*name;
		  	++count;
		  	++x;
		  	++name;
		  }
	child=node->firstdesc;
	SIBLOOP(child) {  /* if node is already a tip, this will just fall through */
		DrawNames(child,treemode);
		}
	return;
}

/*********************************************************************/


void HorizLine( int x1,  int x2,  int y)
{
 int i;
if (x1>x2) swap(&x1,&x2);
for (i=x1;i<=x2;i++) matrx[y][i]=DASH;
return;
}
void VertLine( int y1,  int y2,  int x)
{
 int i;
if (y1>y2) swap(&y1,&y2);
for (i=y1+1;i<=y2-1;i++) matrx[i][x]=BAR;
matrx[y1][x]=PLUS;
matrx[y2][x]=PLUS;
return;
}
void swap( int *x1,  int *x2)
{
 int temp;
temp=*x1;
*x1=*x2;
*x2=temp;
return;
}

