/*  MODULE FOR STRING LIST STRUCTURES




StrListPtr 		newStrList(void);
StrListPtr 		lastStrNode(StrListPtr node);
StrListPtr 		kthStrNode(StrListPtr node, long k);
int			setkthNode(StrListPtr node, long k, char* s);
int			appendStrList(StrListPtr firstNode, char *s);
void			xprintStrList(StrListPtr aList);
long 			lengthList(StrListPtr node);
char*			getkthStr(StrListPtr node, long k);
void 			catkthStr(StrListPtr list, char* s, long i);
StrListPtr 		newStrListN(long numElements);
void 			freeStrList(StrListPtr node);
long 			findMatchStr(StrListPtr List, char * target)


*/



#include "structures.h"
#include "MyUtilities.h"
#include "memory.h"

StackPtr newStack(int maxElements)
{
StackPtr s;
s=(StackPtr)myMalloc(sizeof(struct Stack));
s->maxElements=maxElements;
s->numElements=0;
s->data=(double*)myMalloc(maxElements*sizeof(double));
return s;
}
void freeStack(StackPtr S)
{
myFree(S->data);
myFree (S);
return;
}
int hasElements(StackPtr S)
{
if (S->numElements>0)
	return 1;
else
	return 0;
}

#define DMIN(a,b)			( (a)<(b) )  ? (a):(b)
void pushD(StackPtr S, double x)
{
int i,top;
double *s;
s=S->data;
top=DMIN(S->numElements,S->maxElements-1);
for (i=top;i>=1;i--)
	s[i]=s[i-1];
s[0]=x;
if (S->numElements < S->maxElements)
	++(S->numElements);
return;
}
double popD(StackPtr S)
{
int i,top;
double *s,x;
s=S->data;
x=s[0];
top=S->numElements;
for (i=0;i<=top-2;i++)
	s[i]=s[i+1];
if (S->numElements > 0)
	--(S->numElements);
return x;
}

/***************************************************/
StrListPtr 
string_list_intersect(StrListPtr s1, StrListPtr s2)


{
long L1,L2;
StrListPtr snew;
snew = NULL;
L1=lengthList(s1);
L2=lengthList(s2);
while(s1)
	{
	if (findMatchStr(s2,s1->s))
		{
		if (!snew)
			snew = newStrList();
		appendStrList(snew, s1->s);
		}
	s1=s1->next;
	}
return snew;
} 

/***************************************************/
int 
string_lists_same(StrListPtr s1, StrListPtr s2)

/* returns 1 if the string lists are the same (i.e., have the same
elements in any order), otherwise returns 0 */
// 	NOT AN EFFICIENT ALGORITHM; SHOULD HASH THIS
{
long L1,L2;
L1=lengthList(s1);
L2=lengthList(s2);
if (L1 != L2)
	return 0;
while(s1)
	{
	if (!findMatchStr(s2,s1->s)) return 0;
	s1=s1->next;
	}
return 1;
} 

/***************************************************/

/* finds the first occurrence of string target in List, and returns the index for that 
string (on 1..n).  Returns 0 if not found. */

long 
findMatchStr(StrListPtr List, char * target)
{
long ix=1;
while (List)
	{
	if (isEqual(target,List->s)) return (ix);
	++ix;
	List = List ->next;
	}
	
return 0;
}

/***************************************************/
void glomStrLists(StrListPtr A,  StrListPtr B)
{
/* (destructively merges two string lists.  Finds the union of the list and
MAKES THE FIRST LIST (A) THIS UNION.  If you want a NEW list write something*/

long lengthB, ixB;
char * Belement;

lengthB = lengthList(B); 
for (ixB=1;ixB<=lengthB; ixB++)
    {
    Belement=getkthStr(B, ixB);
    if (!findMatchStr(A, Belement))
	appendStrList(A, Belement);
    }  
    
}
/***************************************************/
		
/* Prints out the contents of a string list, one element per row */

void 			
xprintStrList(StrListPtr node)
	{
	long count=0;
	
	
	while ((node != NULL)  )
		{
		if(node->s != NULL)
			printf("%s\n",node->s);
		node=node->next;
		}
		
	return;
	}

/***************************************************/
		
/* makes head node of a linked list and initializes 
string pointer and next pointer to NULL

--returns NULL on error*/

StrListPtr 		
newStrList(void)

	{
	
	StrListPtr node;
	node = (StrListPtr)myMalloc(sizeof(struct StrList));
	if (node != NULL)
		{
		node->next=NULL;
		node->s=NULL;
		}
	return node;
	}
/***************************************************/
		
/* Makes linked list of numElements elements and initializes them all to NULL

--returns NULL on error*/

StrListPtr 		
newStrListN(long numElements)

	{
	
	StrListPtr node;
	long k;
	
	if (numElements<=0 )
		return NULL;
	node = newStrList();
	for (k=1;k<=numElements-1;k++)
		appendStrList(node,NULL);
	
	return node;
	}
/***************************************************/		

/* returns last node in linked list or NULL if error */

StrListPtr 		
lastStrNode(StrListPtr node)

	{
	long count=0;
	if (node != NULL)
		while (node->next != NULL)
			{
			++count;
			node=node->next;
			}
	return node;
	}
/***************************************************/		

/* returns the kth node, or NULL if past end or bad k; k is on [1..size of list] */

StrListPtr 		
kthStrNode(StrListPtr node, long k)

{
	if ( k<0) return NULL; /* Error */
	while (  (node->next != NULL)  && (--k > 0)  )
		{
		node=node->next;	
		}
	return node;
}
/***************************************************/
		
/* Sets the kth element in the list to the value of string s;
	--returns 1 if OK, 0 if error*/

int				
setkthNode(StrListPtr list, long k, char* s)
{
	StrListPtr node;
	node = kthStrNode(list,k);
	if (node != NULL)
		{
		/** major modification in next line to allow clean overwrite of string**/
		if (node->s)
			myFree(node->s);
		node->s=DupStr(s);	/* make a persistent version of this string */
		if (node->s != NULL)
			return 1;
		else
			return 0; /* error in DupStr */
		}
	else
		return 0; /* Error */

}
/***************************************************/
		
/* Returns a ptr to the kth string element in the list, or NULL on error 

NB!  If you are going to stash this string somewhere, use DupStr to make a persistent
copy of the string itself.  Otherwise, some other routine may free the original locatiions
that this pointer points to.  */



char*			
getkthStr(StrListPtr node, long k)
{
	char* s;
	if ( k<0) return NULL; /* Error */
	while (  (node->next != NULL)  && (--k > 0)  )
		{
		node=node->next;	
		};
	if (node != NULL)
		{
		s=node->s;
		return /* DupStr */(s);	/* Return a persistent version of this or else all hell breaks
								loose when we dispose of the string list */
		}
	else
		return NULL; /* Error */

}
/***************************************************/		

/* Adds a new item to the list, but...
	If the last node in the list has a null string, it puts the item
	in that last item rather than making a new node.  This allows
	one to append to a newly created empty list.  
	A duplicate of the string is made before storing.
	If the parameter 'str' is NULL, just adds an empty node with a NULL string;
	Returns 1 if OK, 0 if error */



int		
appendStrList(StrListPtr firstNode, char *str)
{

	StrListPtr	lastNode, newNode;
	char* 		cpStr;
	
	if (firstNode == NULL ) return 0; 	/* error  */
	
	lastNode=lastStrNode(firstNode);
	if (lastNode != NULL)
		{
		if (str == NULL)	/* add a null new node */
			{
			newNode=newStrList();
			lastNode->next=newNode;
			return 1;		
			}
		else				/* add a new node that has a string */
			{
			cpStr=DupStr(str);/* make a persistent copy of the string and store */
			if (cpStr == NULL)
				return 0;	/* error */
			if (lastNode->s == NULL)  /* there's no item stored at this node, so put it here */
				{
				lastNode->s=cpStr;
				return 1;
				}
			else				/* there IS an item at this node, make a new node and store*/
				{
				newNode=newStrList();
				lastNode->next=newNode;
				newNode->s=cpStr;
				return 1;
				}
			}
		}
	else
		return 0; 	/* error */
	}
/***************************************************/		
/* returns number of elements of list */

long 		lengthList(StrListPtr node)

{
	long length=1;
	if (node == NULL) return 0; /* Error */
	while (  node->next != NULL )
		{
		node=node->next;
		++length;	
		}
	return length;
}
/***************************************************/		
/* Finds the kth string in the list and concatenates string 'ss' to it.
If that string is NULL, it just puts 'ss' there. */

void catkthStr(StrListPtr list, char* ss, long i)
	{
	StrListPtr 		ithnode;
	ithnode=kthStrNode(list,i);
	if (ithnode->s == NULL) 
		ithnode->s=DupStr(ss);		/* destination string is null, so just put string here */
	else
		concat(&(ithnode->s),ss);	/* destination exists so concatenate */
	return;
	}
/***************************************************/		
void freeStrList(StrListPtr node)
{
if (node == NULL) return;
if (node->next)
	freeStrList(node->next);
myFree(node->s);
myFree (node);
return;


}

/***************************************************/
/***************************************************/
/* GENERIC POINTER LISTS */

/* NB!  Sizeof refers to the size of the object being pointed to
    by the pointers in the list. */

/* NB!  A list is defined if it has been created by pNewList
 *	A list is defined but empty if it has no items,  which
	is indicated by having a NULL pointer as its .item
 */


/***************************************************/
/***************************************************/
		
/* makes head node of a linked list of void pointers 

--returns NULL on error*/

PtrList pNewList(void)

/*  USE THIS ONE FROM NOW ON! 
    Returns a pointer to the headnode of a generic pointer list.
 * Sets the first item to NULL,  making it an "empty" list
 */

{
	
	PtrList node;
	node = (PtrList)myMalloc(sizeof(struct PtrListStruct));
	node->next=NULL;
	node->item=NULL;
	return node;
}

PtrList pNewListAlt(size_t size) /* OLD! probably should NOT initialize item here! 
	(it makes adding items later tricky. FIX some day.  Used to initialize
	tree lists and the like,  and currently sets up their first nodes) */
{
	
	PtrList node;
	node = (PtrList)myMalloc(sizeof(struct PtrListStruct));
	node->next=NULL;
	node->item=(void *)myMalloc(size);
	
return node;
}
/***************************************************/
long pLengthList(PtrList node)

/** ! seems like this will return a length of 1 even if the first node is empty*/
{
	long length=1;
	if (node == NULL) fatal("list is empty"); /* Error */
	while (  node->next)
		{
		node=node->next;
		++length;	
		}
	return length;
}
PtrList pListLastNode(PtrList node)

	{
	long count=0;
	if (node != NULL)
		while (node->next != NULL)
			{
			++count;
			node=node->next;
			}
	return node;
	}
PtrList pListgetkthNode(PtrList node, long k)

{
	if ( k<0 ) return NULL; /* Error */
	while (  (node->next != NULL)  && (--k > 0)  )
		{
		node=node->next;	
		}
	return node;
}


PtrList pListAddNode(PtrList firstNode, size_t size)

/*** Appends a node of the given size on the list and returns a pointer
	to that node. Item is set to NULL. If firstNode is NULL, 
	create a new EMPTY list ***/

{

	PtrList	lastNode, newNode;
	
	if (firstNode == NULL ) 
		return pNewListAlt(size);
	
	lastNode=pListLastNode(firstNode);
	if (lastNode != NULL)
		{
			newNode=pNewListAlt(size);
			if (newNode)
				{
				lastNode->next=newNode;
				return newNode;	
				}
			else
				fatal("Couldn't allocate newnode in pListAddNode");	
		}
	else
		{
		fatal("Couldn't get lastNode right in pListAddNode");
		return NULL; 	/* error */
		    
		}
	}
/***************************************************/		
void pListAddItem2(PtrList firstNode, size_t size, void * ptrItem)

/** CLUNKY CODE DO NOT USE **/
/*** given a pointer to an item of size 'size', this adds this item
    to the list.  Does this either at the current last node if that has no item
    or at a new next node ***/

	{
	PtrList	newNode,  lastNode;
	if (firstNode==NULL)
	    return; /* error, no list */
	lastNode=pListLastNode(firstNode);
	if (lastNode->item == NULL)
		{
			lastNode->item=ptrItem;
		}
	else
		{
		newNode=pListAddNode(firstNode, size);
		newNode->item=ptrItem;
		}
	return;
	}
/***************************************************/		
void pListAddItem(PtrList firstNode, void * ptrItem)

/*** given a pointer to an item of size 'size', this adds this item
    to the list.  Does this either at the current last node if that has no item
    or at a new next node ***/

	{
	PtrList	newNode,  lastNode;
	if (firstNode==NULL)
	    return; /* error, no list */
	lastNode=pListLastNode(firstNode);
	if (lastNode->item == NULL)
		{
			lastNode->item=ptrItem;
		}
	else
		{
		newNode=pNewList();
		lastNode->next=newNode;
		newNode->item=ptrItem;
		}
	return;
	}
/***************************************************/		
/***************************************************/		
void DfreepList(PtrList node)

/* EXTREME free of generic list */

{
if (node == NULL) return;
if (node->next)
	freepList(node->next);
myFree(node->item);  /** NB! This may be inadequate! Because item may itself
			have pointers to dynamically created objects. OR
			IT MAY BE DESTRUCTIVE,  say,  if we just want to maintain
			a list of pointers to things.  This now ALSO DELETES
			the things pointed TOO! */
myFree (node);
return;


}
void freepList(PtrList node)

/* more GENTLE free of generic list (doesn't destroy the elements of the list,
 * just the pointers to the elements 
 */

{
if (node == NULL) return;
if (node->next)
	freepList(node->next);
myFree (node);
return;
}
/**********************************/
/**********************************/
/* Kludgy utils for set-vectors */
/* items in a set are numbered from 1..size, where size is the max possible
 * number of items.
 * NB! It is assumed that any operations on two sets forces both to have
 * the same 'size'.  Results are undefined otherwise.
 */


void test_set(Set a,  Set b)
{
    printf("\n");
    printf("set 1:       ");print_set(a);
    printf("set 2:       ");print_set(b);
    printf("intersection:");print_set(intersect_set(a, b));
    printf("union:       ");print_set(union_set(a, b));
    if (sets_Equal(a, b))
	printf("Sets EQUAL\n");
    if (is_subset(a, b))
	printf("1 is a subset of 2\n");
    if (is_superset(a, b))
	printf("1 is a superset of 2\n");
    printf("Percent overlap=%f\n", set_overlap(a, b));
    return;
	

}
double set_overlap(Set a, Set b)
    {
	int csize, dsize;
	Set c, d;
	c=intersect_set(a, b);
	d=union_set(a, b);
		
	csize=sizeof_set(c);
	dsize=sizeof_set(d);
	if (dsize==0)
	    return 0.0;
	else
	    return (double)csize/dsize;
	
    }

int is_subset(Set a,  Set b)
    { /* is a a subset of b? can also be equal to b! */
    if (sets_Equal(a, intersect_set(a, b)))
	return 1;
    else 
	return 0;
    }
int is_superset(Set a,  Set b)
    { /* is a superset of b? can also be equal to b! */
    if (sets_Equal(b, intersect_set(a, b)))
	return 1;
    else 
	return 0;
    }

Set newSet(int size)
    {
	Set theSet;
	int j;
	theSet = (Set)myMalloc(sizeof(struct SetStruct));
	theSet->size=size;
	theSet->element=(int *)myMalloc(size*sizeof(int));
	for (j=1;j<=size;j++)
			 remove_from_set(theSet,  j);
	return theSet;
    }
void add_to_set(Set theSet,  int item_id)
    {
    (theSet->element)[item_id-1]=1;
    return;
   }
void remove_from_set(Set theSet,  int item_id)
    {
    (theSet->element)[item_id-1]=0;
    return;
   }

int is_set_member(Set theSet,  int item_id)
    {
	if ((theSet->element)[item_id-1]==1)
	    return 1;
	else 
	    return 0;
    }

void print_set(Set theSet)
    {
    int i;
    for (i=1;i<=theSet->size;i++)
    if (is_set_member(theSet, i))
	printf("1");
    else
	printf("0");
    printf("\n");
    return;
    }

Set union_set(Set a,  Set b)
{
    Set c;
    int item, csize;
    c=newSet(a->size);
    for (item=1;item<=a->size;item++)
	if ((is_set_member(a, item)) || (is_set_member(b, item)))
	    add_to_set(c, item);
    return c;
}
Set intersect_set(Set a,  Set b)
{
    Set c;
    int item;
    c=newSet(a->size);
    for (item=1;item<=a->size;item++)
	if ((is_set_member(a, item)) && (is_set_member(b, item)))
	    add_to_set(c, item);
    return c;
}

int sizeof_set(Set a)
{
int n=0, item;
for (item=1;item<=a->size;item++)
    if (is_set_member(a, item))
	++n;
return n;    
}
int is_empty_set(Set a)
{
int n=0, item;
for (item=1;item<=a->size;item++)
    if (is_set_member(a, item))
	return 0;
return 1;    
}

int sets_Equal( Set a,   Set b)
/* note that two empty sets are considered equal */
/* returns 1 if two sets are equal, zero otherwise */

{
int i;
for (i=1;i<=a->size;i++)
    if (is_set_member(a, i) != is_set_member(b, i))
	return 0;   
return 1;   
}

/**********************************/
