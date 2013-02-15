#include <stdlib.h>
#include <stdio.h>
#include "memory.h"
#include "errno.h"
/* #include <malloc.h> */

/* Modified July 1996 back to standard C without Mac toolbox calls */



void * myMalloc(size_t theSize)
{
void * p;

errno = 0;
p = (char *)malloc(theSize);
		/* print_mem_dbg(); */
if (errno)
	{
	perror("Low Level allocation error in myMalloc");
	exit(1);
	}
return p;
}

void myFree(void * p)
{
errno=0;

free(p);
if (errno)
	{
	perror("Low Level free error in myFree");
	exit(1);
	}
return;

}

void * myRealloc(void * p, size_t theSize)
{
long i,j;
void * pp;

errno=0;
pp = (char *)realloc(p,theSize);
if (errno)
	{
	perror("Low Level reallocation error in myRealloc");
	exit(1);
	}


return pp;
}

#if MEM_DBG
/* if you want to do some serious memory debugging, set this to 1 in header 
but you may have to link to library with -lmalloc in the Makefile.
PROBABLY SGI specific */


void print_mem_dbg(char *file_name,int line)
{
struct mallinfo mi;
long netSpace;
mi=mallinfo();
printf("%s:%d\n[%i  %i  %i  %i  %i  %i  %i]\n", file_name,line, 
	mi.uordblks,mi.usmblks, 
	mi.arena, mi.ordblks, mi.smblks, mi.fsmblks,mi.fordblks);

return;
}



#endif
#if 0 
     struct mallinfo  {
             int arena;         /* total space in arena */
             int ordblks;       /* number of ordinary blocks */
             int smblks;        /* number of small blocks */
             int hblkhd;        /* space in holding block headers */
             int hblks;         /* number of holding blocks */
             int usmblks;       /* space in small blocks in use */
             int fsmblks;       /* space in free small blocks */
             int uordblks;      /* space in ordinary blocks in use */

             int fordblks;      /* space in free ordinary blocks */
             int keepcost;      /* space penalty if keep option */
                                /* is used */
     }
 
#endif
