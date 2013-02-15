#define MAXINSTRING 5000	/*max size of tree description input string */
#define MAXfSTRING 100
#define MAXINTERNALNODENAMES 100	/* needed for a global fixed length array */
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
#define	isvalidtaxchar(c)	isalnum(c) || (ispunct(c) && ((c)!=COMMA) && ((c)!=RIGHTPARENS) && ((c)!=LEFTPARENS) && ((c)!=SPACE))
#define min(a,b)			( (a)<=(b) )  ? (a):(b)
#define max(a,b)			( (a)>=(b) )  ? (a):(b)
#define SIBLOOP(c)			for (; (c); (c)=(c)->sib)
#define isTip(c)			 ( (c)->firstdesc == NULL )
#define MINWIDTH	20	/* this is the minimum width in the window 
						allowed for the tree itself (i.e, minus the taxon names) */ 
#define MAXWIDTH 150
#define MAXHEIGHT 5000



/* STRUCTURES AND PROTOTYPES */




double calcMaxToTip(NODETYPE* node);
double calcMaxToTipRate(NODETYPE* node);
void DrawTree(NODETYPE *root, int treemode,int userMaxWidth);


