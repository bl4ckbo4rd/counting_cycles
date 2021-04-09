#include "problems.h"

int main (int argc, char *argv[]){

    int     N;              //number of nodes
    int     M;              //number of undirected links
    double  epsilon;        //parameter that controls the asymmetry of the graph 0: symm; 1: uncorrelated
    int     seed_g;         //seed graph
    double  r;              //strength of the random part in the initial condition of the BP iteration
    int     P;              //number of patterns, useful only in the case of the Hopfield model
    
    double  th = 0.001;    //threshold for the BP algo to converge
    int     T  = 1000;      //maximum number of BP sweeps we make
    
    Graph G;
    
    string  str_graph_type;
    string mode;
    
    int i = 1;
    if (argc > 1)
        mode = argv[1];
    
    if (mode == "study_toy_graph" or mode == "study_random_disordered" or mode == "study_random_hopfield"){
        
        if (mode == "study_toy_graph"){
            int k;
            if (argc >= 3){
                k = atoi(argv[2]);
                if((k==0 or k==1 or k==2 or k==3) and argc==3){
                    seed_g = 0;
                    if(k==0){
                        N = 0;
                        G.initializeGraph(N, seed_g);
                        f_toyGraph0(G);
                    }
                    if(k==1){
                        N = 5;
                        G.initializeGraph(N, seed_g);
                        f_toyGraph1(G);
                    }
                    if(k==2){
                        N = 4;
                        G.initializeGraph(N, seed_g);
                        f_toyGraph2(G);
                    }
                    if(k==3){
                        N = 4;
                        G.initializeGraph(N, seed_g);
                        f_toyGraph3(G);
                    }
                }
                else if (k==4 and argc==7){
                    i = 3;
                    N = atoi(argv[i++]);
                    M = atoi(argv[i++]);
                    epsilon = atof(argv[i++]);
                    seed_g = atoi(argv[i++]);
                    G.initializeGraph(N, seed_g);
                    f_toyGraph4(G, M, epsilon);
                }
                else{
                    cout << "learn the usage launching the code with only the mode specified!" << endl;
                    return 0;
                }
                    
            }
            else{
                cout << "if mode = study_toy_graph, extra arguments are k, N, M, epsilon, seed_g. " << endl;
                cout << "if k = 1, 2, 3, other arguments are not needed. They are used only for k = 4." << endl;
                return 0;
            }
        }
        if (mode == "study_random_disordered" or mode == "study_random_hopfield"){
            if(argc==8){
                i = 2;
                N       = atoi(argv[i++]);
                M       = atoi(argv[i++]);
                epsilon = atof(argv[i++]);
                seed_g  = atoi(argv[i++]);
                r       = atof(argv[i++]);
                
                if (mode == "study_random_disordered"){
                    str_graph_type = argv[i++];
                    
                    G.initializeGraph(N, seed_g);
                    
                    //maximum connectivity
                    int z=6;
                    
                    if (str_graph_type == "RR") f_RRgraph(G,M,epsilon);
                    if (str_graph_type == "ER") f_ERgraph(G,M,epsilon);
                    if (str_graph_type == "ERC") f_ERCgraph(G,M,epsilon,z);
                    if (str_graph_type == "NS") f_ERNSgraph(G,M,epsilon,z);
                }
                else if (mode == "study_random_hopfield"){
                    P = atoi(argv[i++]);
                    
                    G.initializeGraph(N, seed_g);
                    vector < vector <int> > xi = f_SparseHopfield(G,M,P);
                    
                }
                
            }
            else{
                cout << "if mode = study_random_disordered, extra arguments are N, M, epsilon, seed_g, r, str_graph_type " << endl;
                cout << "if mode = study_random_hopfield, extra arguments are N, M, epsilon, seed_g, r, P " << endl;
                return 0;
            }
            
            
            
        }
        
    }
    else{
        cout << "specify a mode" << endl;
        cout << "mode must be \n-'study_toy_graph' or \n-'study_random_disordered' or \n-'study_random_hopfield'" << endl;
        return 0;
    }
    
    
    cout << endl;
    cout << "----------------- These are the connected components -----------------\n";
    G.connectedComponents();
    cout << "----------------------------------------------------------------------\n";
    cout << endl;
    
    
    
    clock_t start1 = clock();
    f_BPiterationL1(G, th, T, r);
    //f_BPGD_L1(G, th, T, r);
    clock_t end1 = clock();
    float sec1 = (float)(end1 - start1) / CLOCKS_PER_SEC;
    cout << "time elapsed for L=1 ----------------------> " << sec1 << endl;

    
    
    clock_t start2 = clock();
    f_BPiterationL2(G, th, T, r);
    clock_t end2 = clock();
    float sec2 = (float)(end2 - start2) / CLOCKS_PER_SEC;
    cout << "time elapsed for L=2 ----------------------> " << sec2 << endl;
    
    /*
    clock_t start21 = clock();
    f_BPiterationL1T1(G, th, T, r);
    clock_t end21 = clock();
    float sec21 = (float)(end21 - start21) / CLOCKS_PER_SEC;
    cout << "time elapsed for L=1-T=1 ----------------------> " << sec21 << endl;
    */
  
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
    
    
    return 1;

}
