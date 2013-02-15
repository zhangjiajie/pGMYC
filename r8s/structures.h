#ifndef _STRUCTURES
#define _STRUCTURES
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define MAX_LIST_SIZE	1000

#define STU	0	/* STRING TO UPPER? */
#define isEqual(a,b)		(!strcasecmp((a),(b))) 
#define LISTLOOP(c)		for (; (c); (c)=(c)->next)


/* My list structure */

struct Stack {
				int maxElements;
				int numElements;
				double *data;
		};

struct StrList {
				char* s;
				struct StrList* next;
				};
struct PtrListStruct {
				void * item;
				struct PtrListStruct* next;
				};

struct SetStruct {
				int size;	/* max possible elements */
				int *element;	/* pointer to an int array of
						0's and 1's */
		   };

typedef struct Stack		* StackPtr;
typedef struct PtrListStruct	* PtrList;			
typedef struct SetStruct	* Set;
typedef struct StrList		* StrListPtr;

/*************************************/

StackPtr		newStack(int maxElements);
void 			freeStack(StackPtr S);
void 			pushD(StackPtr S, double x);
double			popD(StackPtr S);
int 			hasElements(StackPtr S);

StrListPtr 		string_list_intersect(StrListPtr s1, StrListPtr s2);
int			is_subset(Set a,  Set b);
int			is_superset(Set a,  Set b);
void			test_set(Set a,  Set b);
int			is_empty_set(Set a);
int			sizeof_set(Set a);
Set			intersect_set(Set a,  Set b);
Set			union_set(Set a,  Set b);
void			print_set(Set theSet);
int			is_set_member(Set theSet,  int item_id);
void			add_to_set(Set theSet,  int item_id);
void			remove_from_set(Set theSet,  int item_id);
Set			newSet(int size);
void			glomStrLists(StrListPtr A,  StrListPtr B);
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
long			findMatchStr(StrListPtr List, char * target);
int			string_lists_same(StrListPtr s1, StrListPtr s2);

void			pListAddItem(PtrList firstNode, void * ptrItem);
void			pListAddItem2(PtrList firstNode, size_t size, void * ptrItem);
void			DfreepList(PtrList node);
void			freepList(PtrList node);
PtrList			pNewListAlt(size_t size);
PtrList			pNewList(void);
PtrList			pListLastNode(PtrList node);
PtrList			pListgetkthNode(PtrList node, long k);
PtrList			pListAddNode(PtrList firstNode, size_t size);
long 			pLengthList(PtrList node);
double			set_overlap(Set a, Set b);
int sets_Equal( Set a,   Set b);
#endif

