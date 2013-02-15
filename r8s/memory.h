#include <stdlib.h>
#define MEM_DBG 0	/* for memory debugging */

void * myMalloc(size_t theSize);
void myFree(void * p);
void * myRealloc(void * p, size_t theSize);
void print_mem_dbg(char *file_name,int line);
