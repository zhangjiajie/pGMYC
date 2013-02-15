#define INTEGER	    0
#define	REAL	    1
#define FLAG	    2
#define STRING	    3
#define CHARACTER   4
#define MAX_COMMANDS 25
struct cList {
		int variable_type;
		char * option_name;
		void * variable_address;
		};



/*
 * Processes a command consisting of 'option=value' strings separated by spaces and ending in ';'.
 * The array of structures 'comList' contains the syntax.  Each element of the array has the three
 * structure members corresponding to the option name string that is looked for,  the type of variable
 * it is,  and the address of a variable that the 'value' will be stored in.  For integers and reals this
 * member is just a pointer to the stored values.  For type FLAG,  the parser will look for 'option=YES' or
 * 'option=NO' and store an integer 1 or 0 respectively.  For type STRING,  the pointer points to a location
 * where a pointer to the string will be stored.  Hence to retrieve that string we have to dereference twice.
  */
  
void dummy(void)
{
char *tHndl;
double min_age, max_age; 

 
  
struct cList b[MAX_COMMANDS] =
    {
    {STRING, "TAXON", &tHndl}, 
    {REAL, "MIN_AGE", &min_age}, 
    {REAL, "MAX_AGE", &max_age} 
    };
}

void parse_command(struct cList comList[], int ncommands)
{
char * localTokenPtr,  *dummy;
extern char * aTokenPtr;
int i;
if (ncommands > MAX_COMMANDS)
    fatal("Too many commands in parse_command\n");
while (!isEqual(aTokenPtr=nextToken(),";"))
		{
		for (i=0;i<ncommands;i++)
		    if (isEqual(aTokenPtr,   (comList[i].option_name)))
			{
			if (parse_assignment(comList[i].option_name,&localTokenPtr))
			    {
			    switch (comList[i].variable_type)
				{
				 case INTEGER: *(int *)(comList[i].variable_address)=strtod(localTokenPtr,&dummy);
					myFree(localTokenPtr);
					break;
				 case REAL: *(double *)(comList[i].variable_address)=strtod(localTokenPtr,&dummy); 
					myFree(localTokenPtr);
					break;
				 case STRING: *(char**)(comList[i].variable_address)=localTokenPtr; /* don't free! */
					break;
				 case CHARACTER: *(char *)(comList[i].variable_address)=strtod(localTokenPtr,&dummy); 
					myFree(localTokenPtr);
					break;
				 case FLAG: 
					if (isEqual(localTokenPtr, "YES"))
						*(int *)(comList[i].variable_address)=1;
					else
						*(int *)(comList[i].variable_address)=0;
					myFree(localTokenPtr);
					break;
				}
			    }
			}
			

		
		}

 
 return;   
    
}
