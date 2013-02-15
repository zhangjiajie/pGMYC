#include <stdio.h>
#include <stdlib.h>
#include "storeNexusFile.h"
#include "MyUtilities.h"

#define	MAX_BUFFER_SIZE	20000000  /* THIS ERROR ISN"T ALWAYS CAUGHT APPARENTLY!!*/


/***************

Prompts for a file name, reads the file into a large, fixed length buffer, and
returns a pointer to that buffer.

NEED TO SLURP ANY SIZE FILE; REWRITE WITH A SMALL BUFFER THAT KEEPS ADDING
DYNAMIC SPACE TO THE BIG BUFFER

*/ 


char * storeNexusFile (FILE * inFileStream)

{
	char	*BigBuffer;
	int		c;
	long	count=0,i=0;

	
	
	BigBuffer=(char*)malloc(MAX_BUFFER_SIZE*sizeof(char));
	if (!BigBuffer)	
		{
		doGenericAlert("Could not allocate file buffer");
		return NULL;		
		}

	while ((c=fgetc(inFileStream)) != EOF)	/* have to define c as int so that EOF can be detected*/
		{
		if (i >= MAX_BUFFER_SIZE-1) /* have to save room for terminating null */
			{
			doGenericAlert("Nexus file exceeds 500k maximum");
			return NULL;
			}
		BigBuffer[i]=c;
		++i;
		
		}
		BigBuffer[i]='\0';

return BigBuffer;	
		
}
