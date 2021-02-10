#include "problems.h"

int main (int argc, char *argv[]){

    int     N;              //number of nodes
    int     M;              //number of undirected links
    double  epsilon;        //parameter that controls the asymmetry of the graph 0: symm; 1: uncorrelated
    int     seed_g;         //seed graph
    double  r;              //strength of the random part in the init cond of the BP iteration
    int     P;              //number of patterns, useful only in the case of the Hopfield model
    
    double  th = 0.0001;    //threshold for the BP algo to converge
    int     T  = 1000;      //maximum number of BP sweeps we make
    string  str;
    
    if(argc==8){
        
        int i = 1;
        N       = atoi(argv[i++]);
        M       = atoi(argv[i++]);
        epsilon = atof(argv[i++]);
        seed_g  = atoi(argv[i++]);
        r       = atof(argv[i++]);
        P       = atoi(argv[i++]);
        str     = argv[i++];
        
    }
    else{
        cout << "argument: N, M, epsilon, seed_g, r, P, str" << endl;
        return 0;
    }
    

    
    Graph G(N, seed_g);
    
    //f_toyGraph1(G);
    
    //maximum connectivity
    int z=6;
    
    vector < vector <int> > xi;

    if (str == "RR")
        f_RRgraph(G,M,epsilon);
    if (str == "ER")
        f_ERgraph(G,M,epsilon);
    if (str == "ERC")
        f_ERCgraph(G,M,epsilon,z);
    if (str == "NS")
        f_ERNSgraph(G,M,epsilon,z);
    if (str == "SH")
        xi = f_SparseHopfield(G,M,P);
    
    
    
    cout << endl;
    cout << "----------------- These are the connected components -----------------\n";
    G.connectedComponents();
    cout << "----------------------------------------------------------------------\n";
    cout << endl;
    
    clock_t start1 = clock();
    //f_BPiterationL1(G, th, T, r);
    f_BPGD_L1(G, th, T, r);
    clock_t end1 = clock();
    float sec1 = (float)(end1 - start1) / CLOCKS_PER_SEC;
    cout << "time elapsed for L=1 ----------------------> " << sec1 << endl;

    //f_toyGraph1(G);
    
    /*

    clock_t start2 = clock();
    f_BPiterationL2(G, th, T, r);
    clock_t end2 = clock();
    float sec2 = (float)(end2 - start2) / CLOCKS_PER_SEC;
    cout << "time elapsed for L=2 ----------------------> " << sec2 << endl;
    
    
    
    clock_t start21 = clock();
    f_BPiterationL1T1(G, th, T, r);
    clock_t end21 = clock();
    float sec21 = (float)(end21 - start21) / CLOCKS_PER_SEC;
    cout << "time elapsed for L=1-T=1 ----------------------> " << sec21 << endl;
    
    
  
    clock_t start3 = clock();
    f_BPiterationL3(G, th, T, r);
    clock_t end3 = clock();
    float sec3 = (float)(end3 - start3) / CLOCKS_PER_SEC;
    cout << "time elapsed for L=3 ----------------------> " << sec3 << endl;
    
    
    clock_t start4 = clock();
    f_BPiterationL4(G, th, T, r);
    clock_t end4 = clock();
    float sec4 = (float)(end4 - start4) / CLOCKS_PER_SEC;
    cout << "time elapsed for L=4 ----------------------> " << sec4 << endl;
*/
    
    return 1;

}
