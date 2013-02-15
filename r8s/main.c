/********************************************************

		r8s 
*/

#define VERSION 1.8

/*

********************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "storeNexusFile.h"
#include "nexus.h"



int main(int argc,char * argv[])
    {
    char *p;
    extern int gInteractive;
    extern struct NexDataType *gNexDataPtr;	/* This is THE data structure for the NEXUS data */
    char *theNexusFileBuffer, fnInput[FILENAME_MAX],theArg,CLBuf[256];
    FILE * inStream =NULL;
    int cFlag=0,c=0;
    long l;
    
    gInteractive=1; /* default is interactive mode */
    gNexDataPtr=initialize_nexus();

//    fprintf(stderr,"r8s version %4.2f (compiled %s)\n",VERSION,__DATE__);
    if (argc == 1)
	{
	    doInteractive();
	    return 1;
	}
    
    else
    for (++argv, c=1;c<argc;c++)
	{
	if (**argv=='-')
		{
		p=*argv;
		++p;
		switch(tolower(*p))
				{
				case 'b':
					gInteractive=0;break; /* set to batch mode */
				case 'c':
					++argv;
					strcpy(CLBuf, *argv);
					gInteractive=0;
					cFlag=1;
					break;
			        case 'f':
					++argv;
					strcpy(fnInput, *argv);   /* set file name */
					if (  !(inStream=fopen(fnInput,"r")) )
						{
						printf("Error opening %s\n", fnInput);
						exit(1);
						}
	//				else
	//					fprintf(stderr, "[...reading file %s]\n", fnInput);
					break;
			        case 'v':
					printf("r8s version %4.2f (%s)\n",VERSION,__DATE__);
					break;
			        case 'h':
					printf("Usage: r8s [-b] [-h] [-v] [-f datafile] [-c commandstring]\n");
					printf("\t-b\tBatch process the datafile\n");
					printf("\t-h\tThis information...\n");
					printf("\t-v\tPrint version and compilation date\n");
					printf("\t-c\tOpen and execute commandstring immediately\n");
				}
		}
	++argv; if (!*argv) break;
	}
    
    if (!gNexDataPtr)
	fatal("Failure to allocate nexus data structure in main.c");
    if (inStream)
	    {
	    theNexusFileBuffer=storeNexusFile(inStream);
	    readNexusFile(theNexusFileBuffer);
	    };
    
    if(gInteractive)
	doInteractive();
    if(cFlag)
	doCommandLineControl(CLBuf);    
    
    
    return 1;
    }
