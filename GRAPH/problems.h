#include "algos.h"

void f_ERgraph(Graph &, int, double);

void f_RRgraph(Graph &, int, double);

void f_ERCgraph(Graph &, int, double,int);

vector < vector <int> > f_SparseHopfield(Graph &, int, int);

void f_ERNSgraph(Graph&, int, double, int);

void f_toyGraph0(Graph &);

void f_toyGraph1(Graph &);

void f_toyGraph2(Graph &);

void f_toyGraph3(Graph &);

void f_toyGraph4(Graph &, int, double);

vector <int> dynamical_step(Graph&, vector <int>);

int check_cycle_condition(Graph&, vector <int>, vector <int>, int, int);

int check_L1basins_condition(Graph&, vector <int>, vector <int>, vector <int>, int, int);

int check_L2basins_condition(Graph&, vector <int>, vector <int>, vector <int>, vector <int>, int, int);

void countFixedPointsBruteForce(Graph&);

void count2CyclesBruteForce(Graph&);

void count3CyclesBruteForce(Graph&);

void count4CyclesBruteForce(Graph&);

void countL1T1BasinsBruteForce(Graph&);

void countL1T2BasinsBruteForce(Graph&);

void countL1T3BasinsBruteForce(Graph&);

void countL2T1BasinsBruteForce(Graph&);

void countL2T2BasinsBruteForce(Graph&);

void countL2T3BasinsBruteForce(Graph&);

void f_BPiterationL1(Graph &, double, int, double);

void f_BPGD_L1(Graph& G, double th, int T, double r);

void f_BPGD_L1T1(Graph& G, double th, int T, double r);

void f_BPiterationL2(Graph &, double, int, double);

void f_BPiterationL3(Graph &, double, int, double);

void f_BPiterationL4(Graph &, double, int, double);

void f_BPiterationL1T1(Graph&, double, int, double);

void f_BPiterationL1T2(Graph&, double, int, double);

void f_BPiterationL1T3(Graph&, double, int, double);

void f_BPiterationL2T1(Graph&, double, int, double);

void f_BPiterationL2T2(Graph&, double, int, double);

void f_BPiterationL2T3(Graph&, double, int, double);

