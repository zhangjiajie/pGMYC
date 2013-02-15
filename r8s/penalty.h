double penalty_rate(double p[]);
void debug_check_feasible(NODETYPE *node);
void check_feasible(NODETYPE *node);
double traversePenalty(NODETYPE *node, double p[]);
double addconstr(double x);
double penalty(double p[]);
int check_initial_point(double (*objective)(double p[]),  double p[]);
