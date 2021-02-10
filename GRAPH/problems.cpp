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
    f_build_toy1(G);
    
    
    
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
    
    mess.Wrap_computeLastMarg<L2>();
    mess.linkMarginalState();
    
    
}



void f_BPiterationL3(Graph& G, double th, int T, double r){
    
    int L = 3;
    

    Messages mess(G,L,r);
    
    //here we iterate the BP messages for T BP-time sweeps until the largest
    //difference found between marginals is smaller or equal to th.
    

    bool verbose = 0;
    int flag_red = 1;       //set to 1 to study the reduced BP equations
    int flag_approx = 0;    //set to 1 to study approximate BP updating
    mess.BPiteration<L3>(th, flag_red, flag_approx, T, verbose);
    
    
    
}

void f_BPiterationL4(Graph& G, double th, int T, double r){
    
    int L = 4;
    
    Messages mess(G,L,r);
    
    //here we iterate the BP messages for T BP-time sweeps until the largest
    //difference found between marginals is smaller or equal to th.
    
    bool verbose = 0;
    int flag_red = 1;       //set to 1 to study the reduced BP equations
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


//graph 0 has N=5, M=0
void f_toyGraph0(Graph& G){
    
    G.graphStructure();
    
    countFixedPointsBruteForce(G);
    count2CyclesBruteForce(G);
    count3CyclesBruteForce(G);
    count4CyclesBruteForce(G);
    
    return;
    
}



//graph 1 has N=5, M=4
void f_build_toy1(Graph& G){

    int l = 0;
    int flag;
    vector<int> v;

    v = make_vector<int>() << 0 << 1;
    double J1to0 = 1;
    flag=G.addLink(l, v, J1to0);
    l++;

    v = make_vector<int>() << 1 << 0;
    double J0to1 = 1;
    flag=G.addLink(l, v, J0to1);
    l++;
    
    v = make_vector<int>() << 1 << 2;
    double J2to1 = -1;
    flag=G.addLink(l, v, J2to1);
    l++;
    
    v = make_vector<int>() << 2 << 1;
    double J1to2 = -1;
    flag=G.addLink(l, v, J1to2);
    l++;
    
    v = make_vector<int>() << 1 << 4;
    double J4to1 = -1;
    flag=G.addLink(l, v, J4to1);
    l++;
    
    v = make_vector<int>() << 4 << 1;
    double J1to4 = -1;
    flag=G.addLink(l, v, J1to4);
    l++;

    v = make_vector<int>() << 2 << 3;
    double J3to2 = -1;
    flag=G.addLink(l, v, J3to2);
    l++;
    
    v = make_vector<int>() << 3 << 2;
    double J2to3 = -1;
    flag=G.addLink(l, v, J2to3);
    
    
    
}

void f_toyGraph1(Graph& G){
    
    f_build_toy1(G);
    G.graphStructure();

    countFixedPointsBruteForce(G);
    //count2CyclesBruteForce(G);
    //count3CyclesBruteForce(G);
    //count4CyclesBruteForce(G);


    return;
    
}

//graph 1 has N=4, M=4
void f_toyGraph2(Graph& G){
    
    int l = 0;
    int flag;
    vector<int> v;
    
    v = make_vector<int>() << 0 << 3;
    double J3to0 = 1.5;
    flag=G.addLink(l, v, J3to0);
    l++;
    
    v = make_vector<int>() << 0 << 1;
    double J1to0 = 0.3;
    flag=G.addLink(l, v, J1to0);
    l++;
    
    v = make_vector<int>() << 0 << 2;
    double J2to0 = 1.0;
    flag=G.addLink(l, v, J2to0);
    l++;
    
    v = make_vector<int>() << 1 << 0;
    double J0to1 = 0.7;
    flag=G.addLink(l, v, J0to1);
    l++;
    
    v = make_vector<int>() << 1 << 2;
    double J2to1 = 1.0;
    flag=G.addLink(l, v, J2to1);
    l++;
    
    v = make_vector<int>() << 2 << 1;
    double J1to2 = 1.2;
    flag=G.addLink(l, v, J1to2);
    l++;
    
    v = make_vector<int>() << 2 << 0;
    double J0to2 = 0.4;
    flag=G.addLink(l, v, J0to2);
    l++;
    
    v = make_vector<int>() << 3 << 0;
    double J0to3 = 0.5;
    flag=G.addLink(l, v, J0to3);
    l++;
    
    
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
        
        double flag_c = 1;

        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
            
            n = G.v[it_n->n].n;
            
            h = 0.;
            
            for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
                vector<int>::iterator it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
                int index_k = distance (G.v[n].v_neigh.begin(), it);
                
                int index_J_k_to_n = G.v[n].v_link[index_k];
                J_k_to_n = G.E[index_J_k_to_n].J;
                
                h += J_k_to_n  * (2 * x[r][*it_k] - 1);
                
            }
            
            double g;
            if (h == 0)
                //g = 2 * x[r][n] -  1;
                g = -  1;
            else
                g = h;
            
            if (g * (2 * x[r][n] -  1) > 0)
                flag[n] = 1;
            else
                flag[n] = 0;
            
        }
        
        prod = 1;
        for (int n = 0; n < N; n++)
            flag_c *= flag[n];
        
        if (flag_c == 1) {
            count++;
            cout << "********************************** fixed point: ******************************" << endl;
            for (int k = 0; k < N; k++)
                cout << x[r][k] << " ";
            cout << endl;
            
        }
        
        
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
        
        vector <int> flag;
        flag.resize(N);
        
        vector <int> xp;
        xp.resize(N);
        
        vector <int> xs;
        xs.resize(N);
        
        double flag_c = 1;
        
        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
            
            n = G.v[it_n->n].n;
            
            h = 0.;
            
            for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
                vector<int>::iterator it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
                int index_k = distance (G.v[n].v_neigh.begin(), it);
                
                int index_J_k_to_n = G.v[n].v_link[index_k];
                J_k_to_n = G.E[index_J_k_to_n].J;
                
                h += J_k_to_n  * (2 * x[r][*it_k] - 1);
                
            }
            
            if (h > 0)
                xp[n] = 1;
            else if (h < 0)
                xp[n] = 0;
            else if (h == 0)
                //xp[n] = x[r][n];
                xp[n] = 0;
                
            
        }
        
        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
            
            n = G.v[it_n->n].n;
            
            h = 0.;
            
            for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
                vector<int>::iterator it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
                int index_k = distance (G.v[n].v_neigh.begin(), it);
                
                int index_J_k_to_n = G.v[n].v_link[index_k];
                J_k_to_n = G.E[index_J_k_to_n].J;
                
                h += J_k_to_n  * (2 * xp[*it_k] - 1);
                
            }
            
            if (h > 0)
                xs[n] = 1;
            else if (h < 0)
                xs[n] = 0;
            else if (h ==0)
                //xs[n] = xp[n];
                xs[n] = 0;
            
        }
        
        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){

            n = G.v[it_n->n].n;
            
            if ( ((2 * xs[n] - 1) * (2 * x[r][n] -  1)) > 0)
                flag[n] = 1;
            else
                flag[n] = 0;
            
        }
    
        prod = 1;
        for (int n = 0; n < N; n++)
            flag_c *= flag[n];
        
        if (flag_c == 1) {
            count++;
            cout << "********************************** cycles L = 2: ******************************" << endl;
            for (int k = 0; k < N; k++)
                cout << x[r][k] << " ";
            cout << endl;
            
        }
        
        
    }
    
    cout << "number of cycles: " << count << endl;
    
    
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
        
        vector <int> flag;
        flag.resize(N);
        
        vector <int> xp;
        xp.resize(N);
        
        vector <int> xs;
        xs.resize(N);
        
        vector <int> xt;
        xt.resize(N);
        
        double flag_c = 1;
        
        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
            
            n = G.v[it_n->n].n;
            
            h = 0.;
            
            for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
                vector<int>::iterator it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
                int index_k = distance (G.v[n].v_neigh.begin(), it);
                
                int index_J_k_to_n = G.v[n].v_link[index_k];
                J_k_to_n = G.E[index_J_k_to_n].J;
                
                h += J_k_to_n  * (2 * x[r][*it_k] - 1);
                
            }
            
            if (h > 0)
                xp[n] = 1;
            else if (h < 0)
                xp[n] = 0;
            else if (h ==0)
                //xp[n] = x[r][n];
                xp[n] = -1;
                
            
        }
        
        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
            
            n = G.v[it_n->n].n;
            
            h = 0.;
            
            for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
                vector<int>::iterator it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
                int index_k = distance (G.v[n].v_neigh.begin(), it);
                
                int index_J_k_to_n = G.v[n].v_link[index_k];
                J_k_to_n = G.E[index_J_k_to_n].J;
                
                h += J_k_to_n  * (2 * xp[*it_k] - 1);
                
            }
            
            if (h > 0)
                xs[n] = 1;
            else if (h < 0)
                xs[n] = 0;
            else if (h ==0)
                //xs[n] = xp[n];
                xs[n] = -1;
            
        }
        
        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
            
            n = G.v[it_n->n].n;
            
            h = 0.;
            
            for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
                vector<int>::iterator it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
                int index_k = distance (G.v[n].v_neigh.begin(), it);
                
                int index_J_k_to_n = G.v[n].v_link[index_k];
                J_k_to_n = G.E[index_J_k_to_n].J;
                
                h += J_k_to_n  * (2 * xs[*it_k] - 1);
                
            }
            
            if (h > 0)
                xt[n] = 1;
            else if (h < 0)
                xt[n] = 0;
            else if (h ==0)
                //xt[n] = xs[n];
                xt[n] = -1;
            
        }
        
        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
            
            n = G.v[it_n->n].n;
            
            if ( ((2 * xt[n] - 1) * (2 * x[r][n] -  1)) > 0)
                flag[n] = 1;
            else
                flag[n] = 0;
            
        }
        
        prod = 1;
        for (int n = 0; n < N; n++)
            flag_c *= flag[n];
        
        if (flag_c == 1) {
            count++;
            cout << "********************************** cycles L = 3: ******************************" << endl;
            for (int k = 0; k < N; k++)
                cout << x[r][k] << " ";
            cout << endl;
            
        }
        
    }
    
    cout << "number of cycles: " << count << endl;
    
    
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
        
        vector <int> flag;
        flag.resize(N);
        
        vector <int> xp;
        xp.resize(N);
        
        vector <int> xs;
        xs.resize(N);
        
        vector <int> xt;
        xt.resize(N);
        
        vector <int> xq;
        xq.resize(N);
        
        double flag_c = 1;
        
        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
            
            n = G.v[it_n->n].n;
            
            h = 0.;
            
            for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
                vector<int>::iterator it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
                int index_k = distance (G.v[n].v_neigh.begin(), it);
                
                int index_J_k_to_n = G.v[n].v_link[index_k];
                J_k_to_n = G.E[index_J_k_to_n].J;
                
                h += J_k_to_n  * (2 * x[r][*it_k] - 1);
                
            }
            
            if (h > 0)
                xp[n] = 1;
            else if (h < 0)
                xp[n] = 0;
            else if (h == 0)
                //xp[n] = x[r][n];
                xp[n] = -1;
                
            
        }
        
        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
            
            n = G.v[it_n->n].n;
            
            h = 0.;
            
            for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
                vector<int>::iterator it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
                int index_k = distance (G.v[n].v_neigh.begin(), it);
                
                int index_J_k_to_n = G.v[n].v_link[index_k];
                J_k_to_n = G.E[index_J_k_to_n].J;
                
                h += J_k_to_n  * (2 * xp[*it_k] - 1);
                
            }
            
            if (h > 0)
                xs[n] = 1;
            else if (h < 0)
                xs[n] = 0;
            else if (h==0)
                //xs[n] = xp[n];
                xs[n] = -1;

            
        }
        
        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
            
            n = G.v[it_n->n].n;
            
            h = 0.;
            
            for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
                vector<int>::iterator it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
                int index_k = distance (G.v[n].v_neigh.begin(), it);
                
                int index_J_k_to_n = G.v[n].v_link[index_k];
                J_k_to_n = G.E[index_J_k_to_n].J;
                
                h += J_k_to_n  * (2 * xs[*it_k] - 1);
                
            }
            
            if (h > 0)
                xt[n] = 1;
            else if (h < 0)
                xt[n] = 0;
            else if (h==0 )
                //xt[n] = xs[n];
                xt[n] = -1;

                
            
        }
        
        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
            
            n = G.v[it_n->n].n;
            
            h = 0.;
            
            for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
                vector<int>::iterator it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
                int index_k = distance (G.v[n].v_neigh.begin(), it);
                
                int index_J_k_to_n = G.v[n].v_link[index_k];
                J_k_to_n = G.E[index_J_k_to_n].J;
                
                h += J_k_to_n  * (2 * xt[*it_k] - 1);
                
            }
            
            if (h > 0)
                xq[n] = 1;
            else if (h < 0)
                xq[n] = 0;
            else if (h==0)
                //xq[n] = xt[n];
                xq[n] = -1;

            
        }
        
        for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
            
            n = G.v[it_n->n].n;
            
            if ( ((2 * xq[n] - 1) * (2 * x[r][n] -  1)) > 0)
                flag[n] = 1;
            else
                flag[n] = 0;
            
        }
        
        prod = 1;
        for (int n = 0; n < N; n++)
            flag_c *= flag[n];
        
        if (flag_c == 1) {
            count++;
            cout << "********************************** cycles L = 4: ******************************" << endl;
            for (int k = 0; k < N; k++)
                cout << x[r][k] << " ";
            cout << endl;
            
        }
        
    }
    
    cout << "number of cycles: " << count << endl;
    
    
    return;
}


