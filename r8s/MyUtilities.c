#include "MyUtilities.h"
#include "memory.h"

/****  Miscellaneous utility commands ****/




char * slurpFile (FILE * inFileStream, long maxSize)

{
        char    *BigBuffer;
        int             c;
        long    count=0,i=0;

        
        
        BigBuffer=(char*)malloc(maxSize*sizeof(char));
        if (!BigBuffer) 
                {
                doGenericAlert("Could not allocate file buffer");
                return NULL;            
                }

        while ((c=fgetc(inFileStream)) != EOF)  /* have to define c as int so that EOF can be detected*/
                {
                if (i >= maxSize-1) /* have to save room for terminating null */
                        {
                        doGenericAlert("Slurped file exceeds allotted maximum");
                        return NULL;
                        }
                BigBuffer[i]=c;
                ++i;
                
                }
                BigBuffer[i]='\0';

return BigBuffer;       
                
}



void fatal(char *s)
{
	long i;
	printf("!FATAL ERROR!");
	doGenericAlert(s);
	for (i=1;i<100000;i++);
	exit(1);
	return;
}

void doGenericAlert(char* errorMsg)
{
char* s;
printf("********************* WARNING **********************\n\n");
printf("%s\n" ,errorMsg);
printf("\n****************************************************\n");
return;
}

void strtoupper(char *s)
{
	char *temps;
	temps=s; 
	while(  *temps ) 
		{
		*temps=toupper(*temps); 
		++temps;
		}
	/* converts string to upper case */
	return;
}
/***************************************************/

void concat(char **destHndl, char *source)

/* concatenates source into destination; have to be careful, because 'realloc' might
create a new destination pointer and we have to make sure to copy that to the old
'destHndl' */

{
char *tempPtr, *destPtr;
long lengthDest,lengthSource, length;

destPtr=*destHndl;

lengthDest=strlen(destPtr);
lengthSource=strlen(source);
length=lengthDest+lengthSource+1;

tempPtr=(char *)myRealloc(destPtr,(length*sizeof(char))); /* "myRealloc" */
if (tempPtr==NULL) 
	fatal("myReallocation error in concat");
if (tempPtr != destPtr) /* myRealloc had to move pointer */
	{
	destPtr=tempPtr;
	*destHndl=destPtr;	/* make sure to save the new pointer */
	}
strcat(destPtr,source);
return;
}

/***************************************************/
		
char*	DupStr(char* s)
{

/* makes a dynamic memory copy of a string -- returns NULL on error*/

char* sNew;
sNew=(char*)myMalloc((strlen(s)+1)*sizeof(char));
if (sNew != NULL)
	strcpy(sNew,s);
return sNew;

}
/***************************************************/

char*	makeEmptyString(void)
{
char* s;
s=(char*)myMalloc(sizeof(char));
if (s != NULL)
	*s = '\0';
return  s;
}
/***************************************************/
FILE* PromptFileName(char* promptMsg, char* mode)

{
FILE* fpntr;
char fnIn[FILENAME_MAX];	/* defined in stdio.h */
printf("%s",promptMsg);
scanf("%s",fnIn);
if (  (fpntr=fopen(fnIn,mode)) )
	return fpntr;
else
	fatal("Error in file handling");

}
/***************************************************/
int isStrInteger(char* s)

/* Checks to see if a string represents an arbitrary length integer number */

{
char * p;
p=s;
while (*p)
	{
	if (!isdigit(*p))
		return 0;
	++p;
	}
return 1;
}

/***************************************************/
#define numX	100
#define numY	60

void dumbPlot(double X[],  double Y[], int N)

/* X and Y are 0-offset arrays */


{
    double Xmax, Ymax, Xmin, Ymin, Xdif, Ydif, Xintv, Yintv;
    char m[numX+1][numY+1];
    int ix, iy, Xa, Ya;
    for (ix=0;ix<numX+1;ix++)
	for (iy=0;iy<numY+1;iy++)
	    m[ix][iy]=' ';
    array_minmax(X, N, &Xmin, &Xmax);
    array_minmax(Y, N, &Ymin, &Ymax);
    Xdif=Xmax-Xmin;
    Ydif=Ymax-Ymin;
    Xintv=Xdif/numX;
    Yintv=Ydif/numY;

    printf("Ascii Plot of %i Points\n\n", N);

#if 0
    for (ix=0;ix<N;ix++)
	{
	printf("%f\t%f\n", X[ix], Y[ix]);
	}
#endif
    for (ix=0;ix<N;ix++)
	{
	    Xa=(X[ix]-Xmin)/Xintv;
	    Ya=(Y[ix]-Ymin)/Yintv;
	    m[Xa][Ya]='*';
	}
    
    for (iy=numY;iy>=0;iy--)
	{
	if (iy==numY)
	    printf("%6.2f", Ymax);
	if (iy==0)
	    printf("%6.2f", Ymin);
	printf("\t|");
	for (ix=0;ix<=numX;ix++)
	   printf("%c",  m[ix][iy]);
	printf ("\n");
	}
	printf("\t");
	for (ix=0;ix<=numX;ix++)
	   printf("-");
	printf ("\n");
    printf("%6.2f", Xmin);
    for (ix=0;ix<=numX-12;ix++)
	   printf(" ");
    printf("%6f\n", Xmax);
	      
	
    
}

void array_minmax(double X[], int N,  double *min,  double *max)
{
    int i;
    *min=+1e100;
    *max=-1e100;
    for (i=0;i<N;i++)
	{
	    if (X[i]>*max) *max=X[i];
	    if (X[i]<*min) *min=X[i];
	}
    return;
}

void binHisto(long histo[], long N, long binSize)

/* expects a histogram array (1-off) in which element histo[k] is the count for a variate value of k.
   Puts sum of counts into bins of size binSize. N is the number of elements (array length -1)  */

{
long *binHisto, i,j,base,numBins;
numBins = N/binSize+1;
binHisto=(long *)myMalloc((numBins+1)*sizeof(long));
base=1;
for (i=1;i<=numBins;i++)
	{
	base=1+(i-1)*binSize;
	binHisto[i]=0;
	for (j=base;j<base+binSize;j++)
		binHisto[i]+=histo[j];
	}
printf("Binned histogram\n");
for (i=1;i<=numBins;i++)
	if (binHisto[i]>0)
		printf("%li\t%li\n",i*binSize,binHisto[i]);

}
