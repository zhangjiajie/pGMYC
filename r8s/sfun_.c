void sfun_(int *N,double X[],double *F, double G[])
{
int i;
double T;
*F=0;
for (i=1;i<=*N;i++)
	{

         T    = X[i-1] - i;
         *F    = *F + T*T;
         G[i-1] = 2 * T;
	}
return;
}

