#include "problems.h"

void f_ERgraph(Graph& G, int M, double epsilon){
    
    //we build an ErdosRenyi graph
    G.ErdosRenyi(M, epsilon);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
}

void f_RRgraph(Graph& G, int M, double epsilon){
    
    //we build an ErdosRenyi graph
    G.RandomRegular(M, epsilon);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
}

void f_ERCgraph(Graph& G, int M, double epsilon, int max_c){
    
    //we build a constrained ErdosRenyi graph
    G.ErdosRenyiConstrained(M, epsilon, max_c);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
}

void f_ERNSgraph(Graph& G, int M, double epsilon, int max_c){
    
    //we build a random graph with no single nodes
    G.RandomGraphNoSingle(M, epsilon, max_c);
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
}

vector < vector <int> >  f_SparseHopfield(Graph& G, int M, int P){
    
    int N = G.N;
    
    //xi_s contains all the patterns
    vector < vector <int> > xi_s;

    //here we construct the P patterns
    
    xi_s.resize(P);
    
    cout << "pattern:" << endl;
    
    for (int n = 0; n < P; ++n){
        xi_s[n].resize(N,0);
        for (int i = 0; i < N; ++i){
            if ((double)rand()/RAND_MAX < 0.5)
                xi_s[n][i] = 1;
            else
                xi_s[n][i] = -1;
            
            cout << xi_s[n][i] << " " ;
        }
    }
    cout << endl;
    
    int c = (int) (2 * M / N);
    
    vector < vector < double > > J(N, vector<double>(N,0.));
    
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            for (int m=0; m<P; m++){
                J[i][j] += xi_s[m][i] * xi_s[m][j];
            }
            
            J[i][j] /= c;
        }
    }
    
    
    //we build a symmetric RR graph
    //G.RandomRegular(M, 0);
    //f_build_toy1(G);
    
    
    
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); ++it_l){
        int l = it_l->l;
        int i = G.E[l].v_node[0];
        int j = G.E[l].v_node[1];

        G.E[l].J = J[i][j];
    }
    
    
    //here we give information on the structure of the graph
    G.graphStructure();
    
    return xi_s;
    
}


void f_BPiterationL1(Graph& G, double th, int T, double r){
    
    int L = 1;
    
    Messages mess(G,L,r);
    
    //here we iterate the BP messages for T BP-time sweeps until the largest
    //difference found between marginals is smaller or equal to th.
    
    bool verbose = 0;
    int flag_red = 0;       //set to 1 to study the reduced BP equations
    int flag_approx = 0;    //set to 1 to study approximate BP updating
    mess.BPiteration<L1>(th, flag_red, flag_approx, T, verbose);

}


void f_BPGD_L1(Graph& G, double th, int T, double r){
    
    int L = 1;
    
    Messages mess(G,L,r);
    
    //here we iterate the BP messages for T BP-time sweeps until the largest
    //difference found between marginals is smaller or equal to th.
    
    bool verbose = 0;
    int flag_red = 0;       //set to 1 to study the reduced BP equations
    int flag_approx = 0;    //set to 1 to study approximate BP updating
    

    BPGD Algo(mess);
    Algo.mess.setQ(flag_red);
    
    vector <int> v_bias;
    vector <bool> v_q;
    
    v_bias.push_back(0);
    v_q.push_back(1);
        
    Algo.initDecimation<L1>(v_bias,v_q);
    
    Algo.mess.BPiteration<L1>(th, flag_red, flag_approx, T, verbose);

    Algo.mess.Wrap_computeLastMarg<L1>();
    
    Algo.mess.linkMarginalState();
    
    Algo.mess.nodeMarginals();
    Algo.mess.nodeMarginalState();
    
};


void f_BPiterationL2(Graph& G, double th, int T, double r){
    
    int L = 2;
    
    Messages mess(G,L,r);
    
    //here we iterate the BP messages for T BP-time sweeps until the largest
    //difference found between marginals is smaller or equal to th.
    
    bool verbose = 0;
    int flag_red = 0;       //set to 1 to study the reduced BP equations
    int flag_approx = 0;    //set to 1 to study approximate BP updating
    mess.BPiteration<L2>(th, flag_red, flag_approx, T, verbose);
    
    //mess.Wrap_computeLastMarg<L2>();
    //mess.linkMarginalState();
    //mess.messState();
    
}



void f_BPiterationL3(Graph& G, double th, int T, double r){
    
    int L = 3;
    

    Messages mess(G,L,r);
    
    //here we iterate the BP messages for T BP-time sweeps until the largest
    //difference found between marginals is smaller or equal to th.
    

    bool verbose = 0;
    int flag_red = 0;       //set to 1 to study the reduced BP equations
    int flag_approx = 0;    //set to 1 to study approximate BP updating
    mess.BPiteration<L3>(th, flag_red, flag_approx, T, verbose);
    
    
    
}

void f_BPiterationL4(Graph& G, double th, int T, double r){
    
    int L = 4;
    
    Messages mess(G,L,r);
    
    //here we iterate the BP messages for T BP-time sweeps until the largest
    //difference found between marginals is smaller or equal to th.
    
    bool verbose = 0;
    int flag_red = 0;       //set to 1 to study the reduced BP equations
    int flag_approx = 0;    //set to 1 to study approximate BP updating
    mess.BPiteration<L4>(th, flag_red, flag_approx, T, verbose);
    
    
    
}


void f_BPiterationL1T1(Graph& G, double th, int T, double r){
    
    int L = 2;
    
    Messages mess(G,L,r);
    
    //here we iterate the BP messages for T BP-time sweeps until the largest
    //difference found between marginals is smaller or equal to th.
    
    bool verbose = 0;
    int flag_red = 1;       //set to 1 to study the reduced BP equations
    int flag_approx = 0;    //set to 1 to study approximate BP updating
    mess.BPiteration<L1T1>(th, flag_red, flag_approx, T, verbose);
    
    
    
}


//graph 0 has N=0, M=0
void f_toyGraph0(Graph& G){
    
    G.graphStructure();
    
    countFixedPointsBruteForce(G);
    count2CyclesBruteForce(G);
    count3CyclesBruteForce(G);
    count4CyclesBruteForce(G);
    
    return;
    
}


//graph 1 has N=5, M=4
void f_toyGraph1(Graph& G){
    
    G.addLinkWrapper(0, 1, 1., 1.);
    G.addLinkWrapper(1, 2, -1., -1.);
    G.addLinkWrapper(1, 4, -1., -1.);
    G.addLinkWrapper(2, 3, -1., -1.);

    G.graphStructure();

    countFixedPointsBruteForce(G);
    count2CyclesBruteForce(G);
    count3CyclesBruteForce(G);
    count4CyclesBruteForce(G);


    return;
    
}


//graph 2 has N=4, M=4
void f_toyGraph2(Graph& G){
    
    G.addLinkWrapper(0, 3, 0.5, 1.5);
    G.addLinkWrapper(1, 0, 0.3, 0.7);
    G.addLinkWrapper(2, 0, 1.0, 0.4);
    G.addLinkWrapper(2, 1, 1.0, 1.2);

    
    G.graphStructure();
    
    countFixedPointsBruteForce(G);
    count2CyclesBruteForce(G);
    count3CyclesBruteForce(G);
    count4CyclesBruteForce(G);

    
    return;
    
}


//graph 3 is a 1d chain of N=4 sites
void f_toyGraph3(Graph& G){
        
    G.addLinkWrapper(0, 1, 1, 1);
    G.addLinkWrapper(1, 2, 1, 1);
    G.addLinkWrapper(2, 3, 1, 1);
    
    cout << "number of directed links after loop creation " << G.numberOfTotalLinks() << endl;

    G.M = int(G.numberOfTotalLinks() / 2);
    
    G.graphStructure();
    
    countFixedPointsBruteForce(G);
    count2CyclesBruteForce(G);
    count3CyclesBruteForce(G);
    count4CyclesBruteForce(G);
    
    return;
    
}


//graph 4 is a RR graph with a frustrated loop
//with couplings so large to be decoupled from
//the rest of the graph. this motif prevents cycles
//with period smaller than 8 to be present.
void f_toyGraph4(Graph& G, int M, double epsilon){
    
    G.RandomRegular(M, epsilon);
    
    cout << "number of directed links before loop creation " << G.numberOfTotalLinks() << endl;

    
    G.addLinkWrapper(0, 1, -50, 0.0000001);
    G.addLinkWrapper(1, 2, 50, 0.0000001);
    G.addLinkWrapper(2, 3, 50, 0.0000001);
    G.addLinkWrapper(3, 0, 50, 0.0000001);
    
    
    cout << "number of directed links after loop creation " << G.numberOfTotalLinks() << endl;

    G.M = int(G.numberOfTotalLinks() / 2);
    
    G.graphStructure();
    
    countFixedPointsBruteForce(G);
    count2CyclesBruteForce(G);
    count3CyclesBruteForce(G);
    count4CyclesBruteForce(G);

    
    return;
    
}


void countFixedPointsBruteForce(Graph& G){

    int N = G.N;
    
    //this is the number of configurations
    int Nc = (int) pow(2.,N);
    //x is the vector of configurations: x[r][k] is the value of the spin k in the configuration r
    vector < vector <int> > x;
    x.resize(Nc);
    
    for (int r = 0; r < Nc; r++){
        x[r].resize(N);
        
        for (int k = 0; k < N; k++)
            x[r][k] = ( ( r & (int)pow(2.,k) ) >> k ) ;
    }
    
    int n, count = 0;
    double h, prod;
    double J_k_to_n;
    
    for (int r = 0; r < Nc; r++){
        
        vector <int> flag;
        flag.resize(N);
        
        vector <int> x1(N);
        x1 = dynamical_step(G, x[r]);

        count = check_cycle_condition(G, x[r], x1, 1, count);

        
    }

    cout << "number of fixed points: " << count << endl;
    
    
    return;
}



void count2CyclesBruteForce(Graph& G){
    
    int N = G.N;
    
    //this is the number of configurations
    int Nc = (int) pow(2.,N);
    //x is the vector of configurations: x[r][k] is the value of the spin k in the configuration r
    vector < vector <int> > x;
    x.resize(Nc);
    
    for (int r = 0; r < Nc; r++){
        x[r].resize(N);
        
        for (int k = 0; k < N; k++)
            x[r][k] = ( ( r & (int)pow(2.,k) ) >> k ) ;
    }
    
    int n, count = 0;
    double h, prod;
    double J_k_to_n;
    
    for (int r = 0; r < Nc; r++){
        
        vector <int> x1(N);
        x1 = dynamical_step(G, x[r]);
        vector <int> x2(N);
        x2 = dynamical_step(G, x1);

        count = check_cycle_condition(G, x[r], x2, 2, count);
        
    }
    
    cout << "number of cycles of lenght 2: " << count << endl;
    
    
    return;
}


void count3CyclesBruteForce(Graph& G){
    
    int N = G.N;
    
    //this is the number of configurations
    int Nc = (int) pow(2.,N);
    //x is the vector of configurations: x[r][k] is the value of the spin k in the configuration r
    vector < vector <int> > x;
    x.resize(Nc);
    
    for (int r = 0; r < Nc; r++){
        x[r].resize(N);
        
        for (int k = 0; k < N; k++)
            x[r][k] = ( ( r & (int)pow(2.,k) ) >> k ) ;
    }
    
    int n, count = 0;
    double h, prod;
    double J_k_to_n;
    
    for (int r = 0; r < Nc; r++){
        
        vector <int> x1(N);
        x1 = dynamical_step(G, x[r]);
        vector <int> x2(N);
        x2 = dynamical_step(G, x1);
        vector <int> x3(N);
        x3 = dynamical_step(G, x2);

        count = check_cycle_condition(G, x[r], x3, 3, count);
        
    }
    
    cout << "number of cycles of length 3: " << count << endl;
    
    
    return;
}


void count4CyclesBruteForce(Graph& G){
    
    int N = G.N;
    
    //this is the number of configurations
    int Nc = (int) pow(2.,N);
    //x is the vector of configurations: x[r][k] is the value of the spin k in the configuration r
    vector < vector <int> > x;
    x.resize(Nc);
    
    for (int r = 0; r < Nc; r++){
        x[r].resize(N);
        
        for (int k = 0; k < N; k++)
            x[r][k] = ( ( r & (int)pow(2.,k) ) >> k ) ;
    }
    
    int n, count = 0;
    double h, prod;
    double J_k_to_n;
    
    for (int r = 0; r < Nc; r++){
        
        vector <int> x1(N);
        x1 = dynamical_step(G, x[r]);
        vector <int> x2(N);
        x2 = dynamical_step(G, x1);
        vector <int> x3(N);
        x3 = dynamical_step(G, x2);
        vector <int> x4(N);
        x4 = dynamical_step(G, x3);

        count = check_cycle_condition(G, x[r], x4, 4, count);
        
    }
    
    cout << "number of cycles of length 4: " << count << endl;
    
    
    return;
}

void countL1T1BasinsBruteForce(Graph& G){
    
    int N = G.N;
    
    //this is the number of configurations
    int Nc = (int) pow(2.,N);
    //x is the vector of configurations: x[r][k] is the value of the spin k in the configuration r
    vector < vector <int> > x;
    x.resize(Nc);
    
    for (int r = 0; r < Nc; r++){
        x[r].resize(N);
        
        for (int k = 0; k < N; k++)
            x[r][k] = ( ( r & (int)pow(2.,k) ) >> k ) ;
    }
    
    int n, count = 0, count_fix = 0;
    double h, prod;
    double J_k_to_n;

    cout << "\n\n" << endl; 
    
    for (int r = 0; r < Nc; r++){
        
        vector <int> x1(N);
        x1 = dynamical_step(G, x[r]);
        vector <int> x2(N);
        x2 = dynamical_step(G, x1);   

        if(x[r] != x1){

            count = check_L1basins_condition(G, x[r], x1, x2, 1, count);

        }

        
    }
    
    cout << "\nNumber of confs evolving in cycles of L = 1 in T = 1 steps: " << count << endl;
    
    
    return;
}

void countL1T2BasinsBruteForce(Graph& G){
    
    int N = G.N;
    
    //this is the number of configurations
    int Nc = (int) pow(2.,N);
    //x is the vector of configurations: x[r][k] is the value of the spin k in the configuration r
    vector < vector <int> > x;
    x.resize(Nc);
    
    for (int r = 0; r < Nc; r++){
        x[r].resize(N);
        
        for (int k = 0; k < N; k++)
            x[r][k] = ( ( r & (int)pow(2.,k) ) >> k ) ;
    }
    
    int n, count = 0, count_fix = 0;
    double h, prod;
    double J_k_to_n;

    cout << "\n\n" << endl;
    
    for (int r = 0; r < Nc; r++){
        
        vector <int> x1(N);
        x1 = dynamical_step(G, x[r]);
        vector <int> x2(N);
        x2 = dynamical_step(G, x1);
        vector <int> x3(N);
        x3 = dynamical_step(G, x2);

        if(x[r] != x2 && x[r] != x3 && x1 != x2){

            count = check_L1basins_condition(G, x[r], x2, x3, 1, count);
        }    

        
    }
    
    cout << "\nNumber of confs evolving in cycles of L = 1 in T = 2 steps: " << count << endl;
    
    
    return;
}

void countL1T3BasinsBruteForce(Graph& G){
    
    int N = G.N;
    
    //this is the number of configurations
    int Nc = (int) pow(2.,N);
    //x is the vector of configurations: x[r][k] is the value of the spin k in the configuration r
    vector < vector <int> > x;
    x.resize(Nc);
    
    for (int r = 0; r < Nc; r++){
        x[r].resize(N);
        
        for (int k = 0; k < N; k++)
            x[r][k] = ( ( r & (int)pow(2.,k) ) >> k ) ;
    }
    
    int n, count = 0, count_fix = 0;
    double h, prod;
    double J_k_to_n;

    cout << "\n\n" << endl; 
    
    for (int r = 0; r < Nc; r++){
        
        vector <int> x1(N);
        x1 = dynamical_step(G, x[r]);
        vector <int> x2(N);
        x2 = dynamical_step(G, x1);
        vector <int> x3(N);
        x3 = dynamical_step(G, x2);
        vector <int> x4(N);
        x4 = dynamical_step(G, x3);

        if(x[r] != x3 && x1 != x3 && x2 != x3){

            count = check_L1basins_condition(G, x[r], x3, x4, 1, count);
        }

        
    }
    
    cout << "\nNumber of confs evolving in cycles of L = 1 in T = 3 steps: " << count << endl;
    
    
    return;
}

void countL2T1BasinsBruteForce(Graph& G){
    
    int N = G.N;
    
    //this is the number of configurations
    int Nc = (int) pow(2.,N);
    //x is the vector of configurations: x[r][k] is the value of the spin k in the configuration r
    vector < vector <int> > x;
    x.resize(Nc);
    
    for (int r = 0; r < Nc; r++){
        x[r].resize(N);
        
        for (int k = 0; k < N; k++)
            x[r][k] = ( ( r & (int)pow(2.,k) ) >> k ) ;
    }
    
    int n, count = 0, count_fix = 0;
    double h, prod;
    double J_k_to_n;

    cout << "\n\n" << endl;
    
    for (int r = 0; r < Nc; r++){
        
        vector <int> x1(N);
        x1 = dynamical_step(G, x[r]);
        vector <int> x2(N);
        x2 = dynamical_step(G, x1);
        vector <int> x3(N);
        x3 = dynamical_step(G, x2);

        if(x[r] != x2 && x[r] != x1)

            count = check_L2basins_condition(G, x[r], x1, x2, x3, 2, count);

        
    }
    
    cout << "\nNumber of confs evolving in cycles of L = 2 in T = 1 steps: " << count << endl;
    
    
    return;
}

void countL2T2BasinsBruteForce(Graph& G){
    
    int N = G.N;
    
    //this is the number of configurations
    int Nc = (int) pow(2.,N);
    //x is the vector of configurations: x[r][k] is the value of the spin k in the configuration r
    vector < vector <int> > x;
    x.resize(Nc);
    
    for (int r = 0; r < Nc; r++){
        x[r].resize(N);
        
        for (int k = 0; k < N; k++)
            x[r][k] = ( ( r & (int)pow(2.,k) ) >> k ) ;
    }
    
    int n, count = 0, count_fix = 0;
    double h, prod;
    double J_k_to_n;

    cout << "\n\n" << endl;
    
    for (int r = 0; r < Nc; r++){
        
        vector <int> x1(N);
        x1 = dynamical_step(G, x[r]);
        vector <int> x2(N);
        x2 = dynamical_step(G, x1);
        vector <int> x3(N);
        x3 = dynamical_step(G, x2);
        vector <int> x4(N);
        x4 = dynamical_step(G, x3);

        if(x1!= x2 && x2!= x3 && x1!= x3){

            count = check_L2basins_condition(G, x[r], x2, x3, x4, 2, count);
        }

        
    }
    
    cout << "\nNumber of confs evolving in cycles of L = 2 in T = 2 steps: " << count << endl;
    
    
    return;
}

vector <int> dynamical_step(Graph& G, vector <int> x){
    
    int N = G.N;
    
    vector <int> xp(N,0);
    
    for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
    
        int n = G.v[it_n->n].n;
    
        double h = 0.;
    
        for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
            vector<int>::iterator it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
            int index_k = distance (G.v[n].v_neigh.begin(), it);
            
            int index_J_k_to_n = G.v[n].v_link[index_k];
            double J_k_to_n = G.E[index_J_k_to_n].J;
            
            h += J_k_to_n  * (2 * x[*it_k] - 1);
            
        }
    
        if (h > 0)
            xp[n] = 1;
        else if (h < 0)
            xp[n] = 0;
        else if (h == 0)
            //xp[n] = x[r][n];
            xp[n] = 0;
        
        }
    
    return xp;
}


int check_cycle_condition(Graph& G, vector <int> xi, vector <int> xf, int L, int count){
    int N = G.N;

    vector <int> flag(N);
    double flag_c = 1;
    
    for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){

        int n = G.v[it_n->n].n;
        if ( ((2 * xf[n] - 1) * (2 * xi[n] -  1)) > 0)
            flag[n] = 1;
        else
            flag[n] = 0;
    }
    
    double prod = 1;
    for (int n = 0; n < N; n++)
        flag_c *= flag[n];
        
    if (flag_c == 1) {
        count++;
        cout << "********************************** cycles L = " << L << ": ******************************" << endl;
        for (int k = 0; k < N; k++)
            cout << xi[k] << " ";
        cout << endl;
            
    }
    
    return count;
}

int check_L1basins_condition(Graph& G, vector <int> xp, vector <int> xi, vector <int> xf, int L, int count){
    int N = G.N;

    vector <int> flag(N);
    double flag_c = 1;
    
    for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){

        int n = G.v[it_n->n].n;
        if ( ((2 * xf[n] - 1) * (2 * xi[n] -  1)) > 0 ){

            flag[n] = 1;
            
        }else{
            flag[n] = 0;
        }

    }

    int n;
    
    for (n = 0; n < N; n++){
        flag_c *= flag[n];
    }
    if (flag_c == 1 ) { 
        count++;
        cout << "*************************** state in a basin of a L = " << L << " cycle: **************************" << endl;
        for (int k = 0; k < N; k++)
            cout << xp[k] << " ";
        cout << "--->\t";
        for (int k = 0; k < N; k++)
            cout << xi[k] << " ";
        cout << endl;
            
    }
    
    return count;
}

int check_L2basins_condition(Graph& G, vector <int> xpp, vector <int> xp, vector <int> xi, vector <int> xf, int L, int count){
    int N = G.N;

    vector <int> flag(N);
    vector <int> flagp(N);
    double flag_p = 1;
    double flag_c = 1;
    
    for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){

        int n = G.v[it_n->n].n;
        if ( ((2 * xf[n] - 1) * (2 * xp[n] -  1)) > 0 ){

            flag[n] = 1;
            
        }else{
            flag[n] = 0;
        }
        if (((2 * xi[n] - 1) * (2 * xp[n] -  1)) > 0){

            flagp[n] = 1;

        }else{

            flagp[n] = 0;

        }
    }

    int n;
    
    for (n = 0; n < N; n++){
        flag_c *= flag[n];
        flag_p*=flagp[n];
    }    
    if (flag_c == 1 && flag_p != 1) { 
        count++;
        cout << "*************************** state in a basin of a L = " << L << " cycle: **************************" << endl;

        for (int k = 0; k < N; k++)
            cout << xpp[k] << " ";
        cout << "--->\t";
        for (int k = 0; k < N; k++)
            cout << xi[k] << " ";
        cout << endl;
            
    }
    
    return count;
}
