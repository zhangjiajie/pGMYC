#include "nexus.h"
#include "MyUtilities.h"
#include "memory.h"
#include "structures.h"



#define isNL_NEXUSwhiteSpace(c)  ( strchr(" \t\v\f", (c)) || (((c) <= 6) && ((c) >= 0)))
#define isNL_NEXUSpunct(c) ( strchr(NL_punct,(c)) )
#define NULL_RETURN {*aTokenPtr='\0';return aTokenPtr;}

#define CHECK_OVERFLOW  if (cix>=MAX_TOKEN_SIZE-1) doGenericAlert("Token Size Exceeded in nextToken")
/********************************************************/


/*	
	Gets the next token from input stream 'fpointer', and copies it onto the global
	buffer pointed to by 'aTokenPtr'.  If there is NO next token, we copy a null
	string onto that buffer.  That's a signal for the main caller routine...

	If the global variable gNewLine=1 then the newline characters, '\n' and '\r'
	ARE returned as individual tokens,  when encountered.  The normal state is
	gNewLine=0,  which treats these as white space delimiters too.  The only time
	NEXUS file needs to think about newlines is when reading interleaved matrices!

*/

char *nextToken(void)


	{
	extern char *aTokenPtr; 

	extern char * bufPtr;	/*declared and initialized in readNexusFile.c */
	extern int gNewLine;	/*declared and set in readNexusFile.c */

	char *punct="()[]{}/\\,;:=*\'\"`+";	/* these are NEXUS definitions */
	char *NL_punct="()[]{}/\\,;:=*\'\"`+\r\n";	/* NEXUS definitions plus stuff for newlines*/
	char c;
	int cix=0;	/* counter to monitor token size */
	
	*aTokenPtr='\0';
	
	if  ((c=*bufPtr++) == '\0') NULL_RETURN
	
	/* First block below handles the case where newline characters must be reported*/
	
	if (gNewLine)
	    {
	    while (( isNL_NEXUSwhiteSpace(c) ) || (c=='['))  
			    /* this whole loop is in case multiple comments separated by whitespace */ 
		    {
		    while ( isNL_NEXUSwhiteSpace(c) )  /* skip white space and land on next 'c'*/
				    {
				    c=*bufPtr++;
				    if (c=='\0')  NULL_RETURN;/* check for embedded EOF */
				    }
			    
		    if (c=='[')		/* skip the comment and land on next 'c' after comment */
			    {
			    while (c !=']')
				    {
				    c=*bufPtr++;
				    if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
				    }
			    c=*bufPtr++;	/* get next char after ']' */
			    if (c=='\0') NULL_RETURN;
			    }
		    }
	    
    
	    if (c=='\'')		/* deal with single-quoted tokens */
		    {
    
		    aTokenPtr[cix++]=c;
		    CHECK_OVERFLOW;
		    while (  (c=*bufPtr++) != '\'')
			    {
			    if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
			    aTokenPtr[cix++]=c;	/* this is a valid character in the word, add to token */
			    CHECK_OVERFLOW;
			    }
		    aTokenPtr[cix++]=c;	/* add the terminating quote too */
		    CHECK_OVERFLOW;
		    aTokenPtr[cix]='\0'; /* null terminate the string */
#if STU
		    strtoupper(aTokenPtr);
#endif
		    return(aTokenPtr);
		    }	/* return everything between single quotes, including the quotes, as a token*/		
	    aTokenPtr[cix++]=c;		/* char is either punctuation or part of word, so add it to token */
	    CHECK_OVERFLOW;
	    
	    if (!isNL_NEXUSpunct(c))	/* next char is part of word, so add all of word until white,punct,eof,
									    or Token size exceeded */
		    {
		    for (;;)
			    {
			    c=*bufPtr++;  
			    if (  isNL_NEXUSpunct(c) || isNL_NEXUSwhiteSpace(c) 
						    ||  (c == '\0') )
				    {
				    --bufPtr; /* word is terminated by some c that is not part of word;
								     push c back into stream and deal with it on
								    next call to this function; meantime, break out, 
								    and return this token*/
				    break;
				    };
			    aTokenPtr[cix++]=c;	/* this is a valid character in the word, add to token */
			    CHECK_OVERFLOW;
			    }
		    }
	    }
	else  /* identical to block above except for character test definitions! */
	    {
	    while (( isNEXUSwhiteSpace(c) ) || (c=='['))  
			    /* this whole loop is in case multiple comments separated by whitespace */ 
		    {
		    while ( isNEXUSwhiteSpace(c) )  /* skip white space and land on next 'c'*/
				    {
				    c=*bufPtr++;
				    if (c=='\0')  NULL_RETURN;/* check for embedded EOF */
				    }
			    
		    if (c=='[')		/* skip the comment and land on next 'c' after comment */
			    {
			    while (c !=']')
				    {
				    c=*bufPtr++;
				    if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
				    }
			    c=*bufPtr++;	/* get next char after ']' */
			    if (c=='\0') NULL_RETURN;
			    }
		    }
	    
    
	    if (c=='\'')		/* deal with single-quoted tokens */
		    {
		    aTokenPtr[cix++]=c;
		    CHECK_OVERFLOW;
		    while (  (c=*bufPtr++) != '\'')
			    {
			    if (c=='\0') NULL_RETURN;  /* check for embedded EOF */
			    aTokenPtr[cix++]=c;	/* this is a valid character in the word, add to token */
			    CHECK_OVERFLOW;
			    }
		    aTokenPtr[cix++]=c;	/* add the terminating quote too */
		    CHECK_OVERFLOW;
		    aTokenPtr[cix]='\0'; /* null terminate the string */
#if STU
		    strtoupper(aTokenPtr);
#endif
		    return(aTokenPtr);
		    }	/* return everything between single quotes, including the quotes, as a token*/		
	    aTokenPtr[cix++]=c;		/* char is either punctuation or part of word, so add it to token */
	    CHECK_OVERFLOW;
	    
	    if (!isNEXUSpunct(c))	/* next char is part of word, so add all of word until white,punct,eof,
									    or Token size exceeded */
		    {
		    for (;;)
			    {
			    c=*bufPtr++; 
			    if (  isNEXUSpunct(c) || isNEXUSwhiteSpace(c) 
						    ||  (c == '\0') )
				    {
				    --bufPtr; /* word is terminated by some c that is not part of word;
								     push c back into stream and deal with it on
								    next call to this function; meantime, break out, 
								    and return this token*/
				    break;
				    };
			    aTokenPtr[cix++]=c;	/* this is a valid character in the word, add to token */
			    CHECK_OVERFLOW;
			    }
		    }
	    }




	
	
		
	aTokenPtr[cix]='\0'; /* null terminate the string */
#if STU
	strtoupper(aTokenPtr);
#endif
	return(aTokenPtr);
	}
