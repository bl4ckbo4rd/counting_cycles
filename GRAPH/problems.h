#include "graph.h"

void f_ERgraph(Graph &, int, double);

void f_RRgraph(Graph &, int, double);

void f_ERCgraph(Graph &, int, double,int);

vector < vector <int> > f_SparseHopfield(Graph &, int, int);

void f_ERNSgraph(Graph&, int, double, int);

void f_toyGraph0(Graph &);

void f_build_toy1(Graph& G);
void f_toyGraph1(Graph &);

void f_toyGraph2(Graph &);

void countFixedPointsBruteForce(Graph&);

void count2CyclesBruteForce(Graph&);

void count3CyclesBruteForce(Graph&);

void count4CyclesBruteForce(Graph&);

void f_BPiterationL1(Graph &, double, int, double);

void f_BPGD_L1(Graph& G, double th, int T, double r);

void f_BPiterationL2(Graph &, double, int, double);

void f_BPiterationL3(Graph &, double, int, double);

void f_BPiterationL4(Graph &, double, int, double);

void f_BPiterationL1T1(Graph&, double, int, double);

