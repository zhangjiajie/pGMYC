double norm(double v[],int nl, int nh);
double norm_not_active(double v[],int ac[],int nl, int nh);

double *vector(int nl, int nh);
void free_vector(double *v, int nl, int nh);
void nrerror(char error_text[]);
double **matrix(int nrl, int nrh, int ncl, int nch);
void free_matrix(double **m, int nrl, int nrh, int ncl, int nch);
