#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include <stdlib.h>

void			 doGenericAlert(char* errorMsg);
void 			fatal(char *s);
void			strtoupper(char *s);
void			concat(char **destHndl, char *source);
char*			DupStr(char* s);
char*			makeEmptyString(void);
FILE* 			PromptFileName(char* promptMsg, char* mode);
int				isStrInteger(char* s);
void array_minmax(double X[], int N,  double *min,  double *max);
void dumbPlot(double X[],  double Y[], int N);
char * slurpFile (FILE * inFileStream, long maxSize);
void binHisto(long * histo, long N, long binSize);
