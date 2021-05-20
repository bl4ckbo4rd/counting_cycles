#include "messages.h"


//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Messages
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//
//-----------------------------------------------------------------------------------------------------------------------------------------------------------------//


void Messages::setQ(int flag_red){
    
    if (flag_red == 0){
        
        Q  = pow(4,L);
        QS = pow(2,L);
        
        //in the case we do not use the reduced algo, t_reduced is the identity vector v_i = i.
        t_reduced.resize(Q);
        
        for (int i=0; i<Q; i++)
            t_reduced[i] = i;
        
        factor_S_link   = 1;
        factor_S_node_1 = 1;
        factor_S_node_2 = 1;
        
    }
    
    else {
        
        Q  = pow(4,L)/2;
        QS = pow(2,L)/2;
        
        factor_S_link   = 0.5;
        factor_S_node_1 = 0.5;
        factor_S_node_2 = 2;
        
        int Q2 = 2*Q;
        //in the reduced algo, xi and xk are always ordered copying half of the possibile combinations
        //for L=1, we will thus store only
        //0 0
        //0 1
        //rather than
        //0 0
        //0 1
        //1 0
        //1 1
        
        //in cases where we encounter  1 0, we need thus to re-map it onto 0 1. To this aim, we introduce the following vector
        
        t_reduced.resize(Q2);
        
        for (int i=0; i<Q; i++)
            t_reduced[i] = i;
        for (int i=Q; i<Q2; i++)
            t_reduced[i] = Q2-i-1;
        
    }
    
    
    
    initMessages();
    

    sw = Switch();
}



Messages::Messages(Graph& p_G, int p_L, double p_r) : G (p_G), L(p_L), r(p_r) {
    
    N=G.numberOfTotalNodes();
    M=G.numberOfTotalLinks();           //in the case of directed graphs, this is the number of directed links
    
    mess.resize(N);
    update_mess.resize(N);
    
    bias.resize(N);

    int numLinkUnd = (double) M / 2;    //this is the number of undirected links
    
    super_marginal.resize(numLinkUnd);
    prev_super_marginal.resize(numLinkUnd);
    
    link_marginal.resize(numLinkUnd);
    
    
};


void Messages::initMessages(){
    
    for(vector<Node>::iterator it_i = G.v.begin(); it_i != G.v.end(); ++it_i){
        int size_di = G.v[it_i->n].numberOfLinks();  //di is the set of neighbours attached to node i. size_di is the connectivity of the node
        mess[it_i->n].resize(size_di);
        update_mess[it_i->n].resize(size_di);
        bias[it_i->n].resize(size_di);
                
    }

    
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); ++it_l){
        
        if (G.E[it_l->l].J){
            
            int i = G.E[it_l->l].v_node[0];
            int j = G.E[it_l->l].v_node[1];
            
            
            vector<int>::iterator it = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), j);
            int index_j = distance (G.v[i].v_neigh.begin(), it);
            
            //mess[i][j] will have Q components
            //if Q=4, they are indexed by t = 2 * xi + xj
            //in such a way that
            //t=0 when xi=0 and xj=0;
            //t=1 when xi=0 and xj=1;
            //t=2 when xi=1 and xj=0;
            //t=3 when xi=1 and xj=1;
            
            //if Q=16, they are indexed by t = 8 * xi + 4 * xj + 2 * xip + xjp
            
            //and so on...
            
            vector <long double> tmp_mess(Q,0.);
            vector <long double> ones(Q,1.);
            
            double sum = 0.;
            for (int i = 0; i < Q; i ++){
                tmp_mess[i] = 1./Q + r * (double)rand()/RAND_MAX;
                sum += tmp_mess[i];
            }
            
            for (int i = 0; i < Q; i ++){
                tmp_mess[i] /= Q;
            }
            
            mess[i][index_j].resize(Q);
            mess[i][index_j] = tmp_mess;
            
            update_mess[i][index_j].resize(Q);
            
            bias[i][index_j].resize(Q);
            bias[i][index_j] = ones;
            
        }
        
    }
    
    int L = super_marginal.size();
    
    for (int l = 0; l < L ; l++){
        super_marginal[l].resize(Q,0.);
        prev_super_marginal[l].resize(Q,0.);
  
        link_marginal[l].resize(2);
        link_marginal[l][0].resize(2,0.);
        link_marginal[l][1].resize(2,0.);
    }
    
    
    
}

void Messages::messUpdate(){
    
    int i, j;
    int z;
    double J_j_to_i;
    
    //for each directed link i->j
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); ++it_l){

        i = G.E[it_l->l].v_node[0];
        j = G.E[it_l->l].v_node[1];
        
        
        vector<int>::iterator it = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), j);
        int index_j = distance (G.v[i].v_neigh.begin(), it);
        
        
        //if (it_l->l == 62) cout << "message from node " << i << " to node " << j << endl;
        //compute message mess_i_to_j
        vector <vector <long double> > in_mess;
        
        //couplings from nodes k's to node i.
        vector <double> vec_J_k_to_i;
        
        z = 0; //z is the connectivity of node i less 1.
        
        //for each link k->i
        for(vector<int>::iterator it_k = G.v[i].v_neigh.begin(); it_k != G.v[i].v_neigh.end(); ++it_k){
            //where k is different from j
            if(*it_k != j){
                //if (it_l->l == 62) cout << "k: " << *it_k << endl;
                
                vector<int>::iterator it = find(G.v[*it_k].v_neigh.begin(), G.v[*it_k].v_neigh.end(), i);
                int index_i = distance (G.v[*it_k].v_neigh.begin(), it);
                
                //we store the mess_k_to_i
                in_mess.push_back(mess[*it_k][index_i]);
                
                
                it = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), *it_k);
                int index_k = distance (G.v[i].v_neigh.begin(), it);
                
                int index_J_k_to_i = G.v[i].v_link[index_k];
                
                //if (it_l->l == 62) cout << "Jk->i " << G.E[index_J_k_to_i].J << endl;
                
                vec_J_k_to_i.push_back(G.E[index_J_k_to_i].J);
                z++;
            }
            //where k is equal to j
            else{
                //if (it_l->l == 62) cout << "j: " << *it_k << endl;
                
                vector<int>::iterator it_j = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), *it_k);
                int index_j = distance (G.v[i].v_neigh.begin(), it_j);
                
                int index_J_j_to_i = G.v[i].v_link[index_j];
                
                J_j_to_i = G.E[index_J_j_to_i].J;
                
                //if (it_l->l == 62) cout << "Jj->i: " << G.E[index_J_j_to_i].J << endl;
            }
        }
        
        
        vector <long double> mess_ij(Q,0.);
        
        int xk, xi;
        int t_ki, t_ki_tmp;
        int NC;
        
        double prod_mess;
        
        for (int t_ij = 0; t_ij < Q; t_ij++){
            
            NC = allowed_conf_bp[it_l->l][t_ij].size();
            
            t_ki =0;
            
            //get the values of xi at t=0,1..,L
            for(int t=L; t>=1; t--){
                xi = (t_ij >> (2*t-1) ) & 1;
                t_ki += (1<<(2*t-2)) * xi;
            }
            
            t_ki_tmp = t_ki;

            //if(it_l->l==62)
            //    cout << "%%%%%%%%%%%%% UPDATE - t_ij: " <<t_ij << endl;
            
            for (int rr = 0; rr < NC; rr++){
                
                r = allowed_conf_bp[it_l->l][t_ij][rr];

                //if(it_l->l==62)
                //    cout << "%%%% UPDATE - r: " << r << endl;
                
                //for each value of r (i.e. for each combination of the allowed xk's)
                //we store the product of messages mess_k_to_i in prod_mess.
                prod_mess = 1.;
                
                for (int k = 0; k < z; k++){
                    
                    t_ki = t_ki_tmp;
                    
                    //if(it_l->l==62)
                    //    cout << "                    ";
                    
                    //get the allowed xk at t=0,1..,L
                    for(int t=L; t>=1; t--){
                        xk = ( r >> (t*z - 1 - k) ) & 1;
                        t_ki += (1<<(2*t-1)) * xk;
                        //if(it_l->l==62)
                        //    cout << xk << " ";
                    }
                    //if(it_l->l==62)
                    //    cout << endl;

                    
                    prod_mess *= in_mess[k][t_reduced[t_ki]];
                    
                }
            
                mess_ij[t_ij] += prod_mess * bias[i][index_j][t_ij];
            
            }
            
        }
        
        update_mess[i][index_j] = mess_ij;
        
    }
    

    
};


void Messages::messUpdate_approx(int MM){
    
    std::random_device rd;  //Will be used to obtain a seed for the random number engine
    std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
    std::uniform_int_distribution<> dis(1, RAND_MAX);
    
    
    int i, j;
    int z;
    double J_j_to_i;
    
    //for each directed link i->j
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); ++it_l){
        
        i = G.E[it_l->l].v_node[0];
        j = G.E[it_l->l].v_node[1];
        
        vector<int>::iterator it = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), j);
        int index_j = distance (G.v[i].v_neigh.begin(), it);
        
        //if (it_l->l == 62) cout << "message from node " << i << " to node " << j << endl;
        //compute message mess_i_to_j
        vector <vector <long double> > in_mess;
        
        //couplings from nodes k's to node i.
        vector <double> vec_J_k_to_i;
        
        z = 0; //z is the connectivity of node i less 1.
        
        //for each link k->i
        for(vector<int>::iterator it_k = G.v[i].v_neigh.begin(); it_k != G.v[i].v_neigh.end(); ++it_k){
            //where k is different from j
            if(*it_k != j){
                //if (it_l->l == 62) cout << "k: " << *it_k << endl;
                
                vector<int>::iterator it = find(G.v[*it_k].v_neigh.begin(), G.v[*it_k].v_neigh.end(), i);
                int index_i = distance (G.v[*it_k].v_neigh.begin(), it);
                
                //we store the mess_k_to_i
                in_mess.push_back(mess[*it_k][index_i]);
                
                
                it = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), *it_k);
                int index_k = distance (G.v[i].v_neigh.begin(), it);
                
                int index_J_k_to_i = G.v[i].v_link[index_k];
                
                //if (it_l->l == 62) cout << "Jk->i " << G.E[index_J_k_to_i].J << endl;
                
                vec_J_k_to_i.push_back(G.E[index_J_k_to_i].J);
                z++;
            }
            //where k is equal to j
            else{
                //if (it_l->l == 62) cout << "j: " << *it_k << endl;
                
                vector<int>::iterator it_j = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), *it_k);
                int index_j = distance (G.v[i].v_neigh.begin(), it_j);
                
                int index_J_j_to_i = G.v[i].v_link[index_j];
                
                J_j_to_i = G.E[index_J_j_to_i].J;
                
                //if (it_l->l == 62) cout << "Jj->i: " << G.E[index_J_j_to_i].J << endl;
            }
        }
        
        vector <long double> mess_ij(Q,0.);
        
        int xk, xi;
        int t_ki, t_ki_tmp;
        int NC;
        
        double prod_mess;
        
        for (int t_ij = 0; t_ij < Q; t_ij++){
            
            NC = allowed_conf_bp[it_l->l][t_ij].size();
            
            t_ki =0;
            
            //get the values of xi at t=0,1..,L
            for(int t=L; t>=1; t--){
                xi = (t_ij >> (2*t-1) ) & 1;
                t_ki += (1<<(2*t-2)) * xi;
            }
            
            t_ki_tmp = t_ki;
            
            if (NC > MM){
                
                for (int rr = 0; rr < MM; rr++){
                    
                    int rrr = dis(gen) % NC;
                    
                    r = allowed_conf_bp[it_l->l][t_ij][rrr];
                    
                    //for each value of r (i.e. for each combination of the allowed xk's)
                    //we store the product of messages mess_k_to_i in prod_mess.
                    prod_mess = 1.;
                    
                    for (int k = 0; k < z; k++){
                        
                        t_ki = t_ki_tmp;
                        
                        //get the allowed xk at t=0,1..,L
                        for(int t=L; t>=1; t--){
                            xk = ( r >> (t*z - 1 - k) ) & 1;
                            t_ki += (1<<(2*t-1)) * xk;
                        }
                        
                        prod_mess *= in_mess[k][t_reduced[t_ki]];
                        
                    }
                    
                    mess_ij[t_ij] += (double)NC/MM * prod_mess;
                    
                }
                
            }
            
            else{
                
                for (int rr = 0; rr < NC; rr++){
                    
                    r = allowed_conf_bp[it_l->l][t_ij][rr];
                    
                    //for each value of r (i.e. for each combination of the allowed xk's)
                    //we store the product of messages mess_k_to_i in prod_mess.
                    prod_mess = 1.;
                    
                    for (int k = 0; k < z; k++){
                        
                        t_ki = t_ki_tmp;
                        
                        //get the allowed xk at t=0,1..,L
                        for(int t=L; t>=1; t--){
                            xk = ( r >> (t*z - 1 - k) ) & 1;
                            t_ki += (1<<(2*t-1)) * xk;
                        }
                        
                        prod_mess *= in_mess[k][t_reduced[t_ki]];
                        
                    }
                    
                    mess_ij[t_ij] += prod_mess;
                    
                }
                
            }
            
        }
        
        update_mess[i][index_j] =  mess_ij;
        
    }
    
};




//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



template<typename Tag>
void Messages::look_up_table_bp(int flag_red){
    
    int l;
    int i, j;
    int z, xk;
    double hi;
    double J_j_to_i;
    
    //number of possible states of messages and marginals
    int NS = Q;
    
    vector < vector < vector < int > > > t_allowed_conf_bp(G.numberOfTotalLinks(), vector < vector <int> > (NS, vector<int>() ) );
    
    allowed_conf_bp = t_allowed_conf_bp;
    
    //for each directed link i->j
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); ++it_l){
        
        l = it_l->l;
        
        i = G.E[l].v_node[0];
        j = G.E[l].v_node[1];
        
        //couplings from nodes k's to node i.
        vector <double> vec_J_k_to_i;
        
        z = G.v[i].numberOfLinks()-1;
        
        //for each link k->i
        for(vector<int>::iterator it_k = G.v[i].v_neigh.begin(); it_k != G.v[i].v_neigh.end(); ++it_k){
            //where k is different from j
            if(*it_k != j){
                //if (it_l->l == 62) cout << "k: " << *it_k << endl;
                
                vector<int>::iterator it = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), *it_k);
                int index_k = distance (G.v[i].v_neigh.begin(), it);
                
                int index_J_k_to_i = G.v[i].v_link[index_k];
                
                //if (it_l->l == 62) cout << "Jk->i " << G.E[index_J_k_to_i].J << endl;
                
                vec_J_k_to_i.push_back(G.E[index_J_k_to_i].J);
            }
            //where k is equal to j
            else{
                //if (it_l->l == 62) cout << "j: " << *it_k << endl;
                
                vector<int>::iterator it_j = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), *it_k);
                int index_j = distance (G.v[i].v_neigh.begin(), it_j);
                
                int index_J_j_to_i = G.v[i].v_link[index_j];
                
                J_j_to_i = G.E[index_J_j_to_i].J;
                
                //if (it_l->l == 62) cout << "Jj->i: " << G.E[index_J_j_to_i].J << endl;
            }
        }
        
        constraint_bp<Tag>(z, flag_red, l, vec_J_k_to_i, J_j_to_i);
        
    }
    
}


template<typename Tag>
void Messages::look_up_table(int flag_red){
    
    int n;
    int z;
    
    vector < vector < vector < int > > > t_allowed_conf(G.numberOfTotalNodes(), vector < vector <int> > (QS, vector<int>() ) );

    allowed_conf = t_allowed_conf;
    
    
    //for each site n
    for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
        
        n = G.v[it_n->n].n;
        
        
        //couplings from node k's to node n.
        vector <double> vec_J_k_to_n;
        
        //for each link k->n
        
        z = G.v[n].numberOfLinks();

        for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
            
           
            vector<int>::iterator it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
            int index_k = distance (G.v[n].v_neigh.begin(), it);
            
            int index_J_k_to_n = G.v[n].v_link[index_k];
            
            vec_J_k_to_n.push_back(G.E[index_J_k_to_n].J);
            

            
        }
        
        
        constraint<Tag>(z, flag_red, n, vec_J_k_to_n);
        
    }

}


//--------------------------------------------------------------------------- constraint defining the cycle L=1 -------------------------------------------------------------------------//


inline void Messages::def_constraint_bp(L1, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i) {
    
    int xk;
    double hi;
    int t_ij, t_ki;
    
    L=1;
    
    int Nc = (int) pow(pow(2.,L),z);

    vector < vector <int> > x;                                  //vector of 0/1 where we store the configurations of the cavity neighbours

    compute_x(Nc,z,x);
    

    if (l ==62){
        
        cout << "Printing z and Nc: " << z << " " << Nc << endl;
        // configurations are labelled by their base 10 representation of their sequence:
        // i.e. in the case z = 3, the 3 spins may take value
        // 0 0 0  -> 0
        // 1 0 0  -> 4
        // 0 1 0  -> 2
        // 1 1 0  -> 6
        // 0 0 1  -> 1
        // 1 0 1  -> 5
        // 0 1 1  -> 3
        // 1 1 1  -> 7

        
        for (int r = 0; r < Nc; r++){
            
            cout << "configuration " << endl;
            for (int k = 0; k < z * L; k++){
                cout << x[r][k] << " ";
            }
            cout << endl;
        }
        
    }
    
    
    //this vector has a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc,0.);
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j

        if (l ==62) cout << "conf " << endl;
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_i[k] * (2 * xk - 1);
            if (l ==62){
                cout << vec_J_k_to_i[k] << " * " << (2 * xk - 1) << " plus ";
            }
        }
        if (l ==62) cout << endl;
        
    }
    
    int xi = 0;
    for (int xj = 0; xj < 2; ++xj){
        
        t_ij = 2 * xi + xj;
        
        for (int r = 0; r < Nc; r++){
            
            //local field acting on i
            double hi = J_j_to_i * (2 * xj - 1) + h[r];
            
            //if (l ==62){
                
            //    cout << "                 " << " h on i without j contribution " << h[r] << " with j contribution: " << " " << hi << endl;
            
            //}
            
            double gi;
            if (hi == 0)
                //gi = 2 * xi - 1;
                gi = -1;
            else
                gi = hi;
            
            if (gi * ( 2 * xi - 1 ) > 0){
                
                int s=0;
                for (int k=1; k <= z*L; k++)
                    s += x[r][k-1] * (1<<(L*z-k));
                
                allowed_conf_bp[l][t_ij].push_back(s);
                
                
            //    if (l ==62){
            //         cout << "allowed: " << xi << " " << xj << " " << s << endl;
            //         cout << " h on i " << hi << endl;

            //    }
                
            }
            
            
        }
    }
    
    if (flag_red == 0){
        xi = 1;
        for (int xj = 0; xj < 2; ++xj){
            
            t_ij = 2 * xi + xj;
            
            for (int r = 0; r < Nc; r++){
                
                //local field acting on i
                double hi = J_j_to_i * (2 * xj - 1) + h[r];
                
                double gi;
                if (hi == 0)
                    //gi = 2 * xi - 1;
                    gi = - 1;
                else
                    gi = hi;
                
                if (gi * ( 2 * xi - 1 ) > 0){
                    
                    int s=0;
                    for (int k=1; k <= z*L; k++)
                        s += x[r][k-1] * (1<<(L*z-k));
                    
                    allowed_conf_bp[l][t_ij].push_back(s);
                    
                    /*
                     if (l ==62)
                     cout << "allowed: " << xi << " " << xj << " " << s << endl;
                     */
                    
                }
                
                
            }
        }
    }
        
}



inline void Messages::def_constraint(L1, int flag_red, int z, int n, vector <double>& vec_J_k_to_n) {
    
    int xk;
    
    int L = 1;
    
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector < int > > x;
    
    compute_x(Nc,z,x);
    
    //this vector has a value for each configurations of the spins k's
    //local field acting on n
    vector <double> h(Nc,0.);
    
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's)
        //we compute the value of the local field given by all the neighbours of i
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_n[k] * (2 * xk - 1);
        }
        
    }
    
    
    
    //for (int xn = 0; xn < 2; ++xn){
    int xn = 0;
    for (int r = 0; r < Nc; r++){
        
        double gn;
        if (h[r] == 0)
            //gn = 2 * xn - 1;
            gn = -1;
        else
            gn = h[r];
        
        if ( gn * ( 2 * xn - 1 ) > 0){
            
            int s=0;
            for (int k=1; k <= z*L; k++)
                s += x[r][k-1] * (1<<(L*z-k));
            
            allowed_conf[n][xn].push_back(s);
            
        }
        
        
    }
    
    
    if(flag_red == 0){
        xn = 1;
        
        for (int r = 0; r < Nc; r++){
            
            double gn;
            if (h[r] == 0)
                //gn = 2 * xn - 1;
                gn = -1;
            else
                gn = h[r];
            
            if ( gn * ( 2 * xn - 1 ) > 0){
                
                int s=0;
                for (int k=1; k <= z*L; k++)
                    s += x[r][k-1] * (1<<(L*z-k));
                
                allowed_conf[n][xn].push_back(s);
                
            }
            
        }
    }

}



//--------------------------------------------------------------------------- constraint defining the cycle L=2 -------------------------------------------------------------------------//


inline void Messages::def_constraint_bp(L2, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i) {
    
    int xk, xkp;
    int z2;
    double hi, hip;
    int t_ij;
    
    L=2;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;                                  //vector of 0/1 where we store the configurations of the cavity neighbours
    
    compute_x(Nc,z,x);

    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc,0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    
    
    z2 = 2*z;

    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's and of the xkp's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_i[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's and of the xkp's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_n = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_i[k_n] * (2 * xkp - 1);
            k_n ++;
        }
        
    }
    
    int xi = 0;
    for (int xj = 0; xj < 2; ++xj){
        for (int xip = 0; xip < 2; ++xip){
            for (int xjp = 0; xjp < 2; ++xjp){
                
                t_ij = 8 * xi + 4 * xj + 2 * xip + xjp;
                
                for (int r = 0; r < Nc; r++){
                    //local field acting on i at the first time step
                    hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                    //local field acting on i at the second time step
                    hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                    
                    double gi;
                    if (hi == 0)
                        //gi = 2 * xi - 1;
                        gi = - 1;
                    else
                        gi = hi;
                    
                    double gip;
                    if (hip == 0)
                        //gip = 2 * xip - 1;
                        gip = -1;
                    else
                        gip = hip;
                    
                    if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xi - 1 ) > 0 ){
                        
                        int s=0;
                        for (int k=1; k <= z*L; k++){
                            s += x[r][k-1] * (1<<(L*z-k));
                        }
                        
                        allowed_conf_bp[l][t_ij].push_back(s);
                    }
                    
                }
                
            }
        }
        
    }
    
    if (flag_red == 0){
        
        xi = 1;
        for (int xj = 0; xj < 2; ++xj){
            for (int xip = 0; xip < 2; ++xip){
                for (int xjp = 0; xjp < 2; ++xjp){
                    
                    t_ij = 8 * xi + 4 * xj + 2 * xip + xjp;
                    
                    for (int r = 0; r < Nc; r++){
                        //local field acting on i at the first time step
                        hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                        //local field acting on i at the second time step
                        hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                        
                        double gi;
                        if (hi == 0)
                            //gi = 2 * xi - 1;
                            gi = -1;
                        else
                            gi = hi;
                        
                        double gip;
                        if (hip == 0)
                            //gip = 2 * xip - 1;
                            gip = -1;
                        else
                            gip = hip;
                        
                        if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xi - 1 ) > 0 ){
                            
                            int s=0;
                            for (int k=1; k <= z*L; k++){
                                s += x[r][k-1] * (1<<(L*z-k));
                            }
                            
                            allowed_conf_bp[l][t_ij].push_back(s);
                        }
                        
                    }
                    
                }
            }
            
        }
        
        
    }
    //}
    
}



inline void Messages::def_constraint(L2, int flag_red, int z, int n, vector <double>& vec_J_k_to_n) {
    
    
    int xk, xkp;
    int z2;
    int tn;
    
    int L = 2;
    
    // this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc,z,x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc,0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    
    
    z2 = 2*z;
    
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's and of the xkp's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_n[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's and of the xkp's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_n = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_n[k_n] * (2 * xkp - 1);
            k_n ++;
        }
        
    }
    
    //for (int xn = 0; xn < 2; ++xn){
    int xn = 0;
    for (int xnp = 0; xnp < 2; ++xnp){
        
        tn = 2 * xn + xnp;
        
        for (int r = 0; r < Nc; r++){
            
            double gn;
            if ( h[r] == 0)
                //gn = 2 * xn - 1;
                gn = -1;
            else
                gn = h[r];
            
            double gnp;
            if (hp[r] == 0)
                //gnp =  2 * xnp - 1;
                gnp = -1;
            else
                gnp = hp[r];
                
            
            if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xn - 1 ) > 0 ){
                
                int s=0;
                for (int k=1; k <= z*L; k++){
                    s += x[r][k-1] * (1<<(L*z-k));
                }
                
                allowed_conf[n][tn].push_back(s);
                
            }
            
        }
        
    }
    
    if (flag_red == 0){
        
        xn = 1;
        for (int xnp = 0; xnp < 2; ++xnp){
            
            tn = 2 * xn + xnp;
            
            for (int r = 0; r < Nc; r++){
                
                double gn;
                if ( h[r] == 0)
                    //gn = 2 * xn - 1;
                    gn = -1;
                else
                    gn = h[r];
                
                double gnp;
                if (hp[r] == 0)
                    //gnp =  2 * xnp - 1;
                    gnp = -1;
                else
                    gnp = hp[r];
                
                if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xn - 1 ) > 0 ){
                    
                    int s=0;
                    for (int k=1; k <= z*L; k++){
                        s += x[r][k-1] * (1<<(L*z-k));
                    }
                    
                    allowed_conf[n][tn].push_back(s);
                    
                }
                
            }
            
        }
        
        
    }
    
    
    
    //}
    
}


//-----------------------------------------------------------------------------constraint defining the cycle L=3-------------------------------------------------------------------------//


inline void Messages::def_constraint_bp(L3, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i) {
    
    int xk, xkp, xks;
    int z2, z3;
    double hi, hip, his;
    int t_ij;
    
    L=3;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);

    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_i[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_i[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_i[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
        
        
    }
    
    int xi = 0;
    for (int xj = 0; xj < 2; ++xj){
        for (int xip = 0; xip < 2; ++xip){
            for (int xjp = 0; xjp < 2; ++xjp){
                for (int xis = 0; xis < 2; ++xis){
                    for (int xjs = 0; xjs < 2; ++xjs){
                        
                        
                        t_ij = 32 * xi + 16 * xj + 8 * xip + 4 * xjp + 2 * xis + xjs ;
                        
                        
                        for (int r = 0; r < Nc; r++){
                            //local field acting on i at the first time step
                            hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                            //local field acting on i at the second time step
                            hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                            //local field acting on i at the third time step
                            his = J_j_to_i * (2 * xjs - 1) + hs[r];
                            
                            double gi;
                            if (hi == 0)
                                //gi = 2 * xi - 1;
                                gi = -1;
                            else
                                gi = hi;
                            
                            double gip;
                            if (hip == 0)
                                //gip = 2 * xip - 1;
                                gip = -1;
                            else
                                gip = hip;
                            
                            double gis;
                            if (his == 0)
                                //gis = 2 * xis - 1;
                                gis = -1;
                            else
                                gis = his;
                            

                            if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xi - 1 ) > 0 ){
                                
                                int s=0;
                                for (int k=1; k <= z*L; k++){
                                    s += x[r][k-1] * (1<<(L*z-k));
                                }
                                
                                allowed_conf_bp[l][t_ij].push_back(s);
                                
                            }
                            
                        }
                        
                    }
                }
                
            }
        }
        
        //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
        
    }
    
    if (flag_red == 0){
        xi = 1;
        for (int xj = 0; xj < 2; ++xj){
            for (int xip = 0; xip < 2; ++xip){
                for (int xjp = 0; xjp < 2; ++xjp){
                    for (int xis = 0; xis < 2; ++xis){
                        for (int xjs = 0; xjs < 2; ++xjs){
                            
                            
                            t_ij = 32 * xi + 16 * xj + 8 * xip + 4 * xjp + 2 * xis + xjs ;
                            
                            
                            for (int r = 0; r < Nc; r++){
                                //local field acting on i at the first time step
                                hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                                //local field acting on i at the second time step
                                hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                                //local field acting on i at the third time step
                                his = J_j_to_i * (2 * xjs - 1) + hs[r];
                                
                                double gi;
                                if (hi == 0)
                                    //gi = 2 * xi - 1;
                                    gi = -1;
                                else
                                    gi = hi;
                                
                                double gip;
                                if (hip == 0)
                                    //gip = 2 * xip - 1;
                                    gip = -1;
                                else
                                    gip = hip;
                                
                                double gis;
                                if (his == 0)
                                    //gis = 2 * xis - 1;
                                    gis = -1;
                                else
                                    gis = his;
                                
                                
                                if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xi - 1 ) > 0 ){
                                    
                                    int s=0;
                                    for (int k=1; k <= z*L; k++){
                                        s += x[r][k-1] * (1<<(L*z-k));
                                    }
                                    
                                    allowed_conf_bp[l][t_ij].push_back(s);
                                    
                                }
                                
                            }
                            
                        }
                    }
                    
                }
            }
            
            //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
            
        }
        
        
    }
    
}



inline void Messages::def_constraint(L3, int flag_red, int z, int n, vector <double>& vec_J_k_to_n) {
    
    int xk, xkp, xks;
    int z2, z3;
    int tn;
    
    L=3;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_n[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_n[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_n[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
    }
    
    int xn = 0;
    for (int xnp = 0; xnp < 2; ++xnp){
        for (int xns = 0; xns < 2; ++xns){
            
            tn = 4 * xn + 2 * xnp + xns;
            
            for (int r = 0; r < Nc; r++){
                
                double gn;
                if ( h[r] == 0)
                    //gn = 2 * xn - 1;
                    gn = 1;
                else
                    gn = h[r];
                
                double gnp;
                if (hp[r] == 0)
                    //gnp =  2 * xnp - 1;
                    gnp = -1;
                else
                    gnp = hp[r];
                
                double gns;
                if (hs[r] == 0)
                    //gns =  2 * xns - 1;
                    gns = -1;
                else
                    gns = hs[r];
                
                
                if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xn - 1 ) > 0 ){
                    int s=0;
                    for (int k=1; k <= z*L; k++){
                        s += x[r][k-1] * (1<<(L*z-k));
                    }
                    
                    allowed_conf[n][tn].push_back(s);
                    
                }
                
            }
            
        }
    }
    
    if (flag_red == 0){
        xn = 1;
        for (int xnp = 0; xnp < 2; ++xnp){
            for (int xns = 0; xns < 2; ++xns){
                
                tn = 4 * xn + 2 * xnp + xns;
                
                for (int r = 0; r < Nc; r++){
                    
                    double gn;
                    if ( h[r] == 0)
                        //gn = 2 * xn - 1;
                        gn = -1;
                    else
                        gn = h[r];
                    
                    double gnp;
                    if (hp[r] == 0)
                        //gnp =  2 * xnp - 1;
                        gnp = -1;
                    else
                        gnp = hp[r];
                    
                    double gns;
                    if (hs[r] == 0)
                        //gns =  2 * xns - 1;
                        gns = -1;
                    else
                        gns = hs[r];
                    
                    
                    if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xn - 1 ) > 0 ){
                        int s=0;
                        for (int k=1; k <= z*L; k++){
                            s += x[r][k-1] * (1<<(L*z-k));
                        }
                        
                        allowed_conf[n][tn].push_back(s);
                        
                    }
                    
                }
                
            }
        }
        
    }
    
}


//-----------------------------------------------------------------------------constraint defining the cycle L=4-------------------------------------------------------------------------//


inline void Messages::def_constraint_bp(L4, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i) {
    
    int xk, xkp, xks, xkt;
    int z2, z3, z4;
    double hi, hip, his, hit;
    int t_ij;
    
    L=4;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    //local field acting on i at the fourth time step
    vector <double> ht(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    z4 = 4*z;
    
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_i[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_i[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_i[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
        ht[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_t = 0;
        for (int k = z3; k < z4; k++){
            xkt = x[r][k];
            ht[r] += vec_J_k_to_i[k_t] * (2 * xkt - 1);
            k_t ++;
        }
        
        
        
    }
    
    int xi = 0;
    for (int xj = 0; xj < 2; ++xj){
        for (int xip = 0; xip < 2; ++xip){
            for (int xjp = 0; xjp < 2; ++xjp){
                for (int xis = 0; xis < 2; ++xis){
                    for (int xjs = 0; xjs < 2; ++xjs){
                        for (int xit = 0; xit < 2; ++xit){
                            for (int xjt = 0; xjt < 2; ++xjt){
                        
                                t_ij = 128 * xi + 64 * xj + 32 * xip + 16 * xjp + 8 * xis + 4 * xjs + 2 * xit + xjt;
                                
                                for (int r = 0; r < Nc; r++){
                                    //local field acting on i at the first time step
                                    hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                                    //local field acting on i at the second time step
                                    hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                                    //local field acting on i at the third time step
                                    his = J_j_to_i * (2 * xjs - 1) + hs[r];
                                    //local field acting on i at the fourth time step
                                    hit = J_j_to_i * (2 * xjt - 1) + ht[r];
                                    
                                    
                                    double gi;
                                    if (hi == 0)
                                        //gi = 2 * xi - 1;
                                        gi = -1;
                                    else
                                        gi = hi;
                                    
                                    double gip;
                                    if (hip == 0)
                                        //gip = 2 * xip - 1;
                                        gip = -1;
                                    else
                                        gip = hip;
                                    
                                    double gis;
                                    if (his == 0)
                                        //gis = 2 * xis - 1;
                                        gis = -1;
                                    else
                                        gis = his;
                                    
                                    double git;
                                    if (hit == 0)
                                        //git = 2 * xit - 1;
                                        git = -1;
                                    else
                                        git = hit;
                                    
                                    if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xit - 1 ) > 0 && git * ( 2 * xi - 1 ) > 0 ){
                                        
                                        int s=0;
                                        for (int k=1; k <= z*L; k++){
                                            s += x[r][k-1] * (1<<(L*z-k));
                                        }
                                        
                                        allowed_conf_bp[l][t_ij].push_back(s);
                                        
                                    }
                                    
                                }
                                
                            }
                            
                        }
                        
                    }
                }
                
            }
        }
        
        //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
        
    }
    
    if (flag_red == 0){
        xi = 1;
        for (int xj = 0; xj < 2; ++xj){
            for (int xip = 0; xip < 2; ++xip){
                for (int xjp = 0; xjp < 2; ++xjp){
                    for (int xis = 0; xis < 2; ++xis){
                        for (int xjs = 0; xjs < 2; ++xjs){
                            for (int xit = 0; xit < 2; ++xit){
                                for (int xjt = 0; xjt < 2; ++xjt){
                                    
                                    t_ij = 128 * xi + 64 * xj + 32 * xip + 16 * xjp + 8 * xis + 4 * xjs + 2 * xit + xjt;
                                    
                                    for (int r = 0; r < Nc; r++){
                                        //local field acting on i at the first time step
                                        hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                                        //local field acting on i at the second time step
                                        hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                                        //local field acting on i at the third time step
                                        his = J_j_to_i * (2 * xjs - 1) + hs[r];
                                        //local field acting on i at the fourth time step
                                        hit = J_j_to_i * (2 * xjt - 1) + ht[r];
                                        
                                        double gi;
                                        if (hi == 0)
                                            //gi = 2 * xi - 1;
                                            gi = -1;
                                        else
                                            gi = hi;
                                        
                                        double gip;
                                        if (hip == 0)
                                            //gip = 2 * xip - 1;
                                            gip = -1;
                                        else
                                            gip = hip;
                                        
                                        double gis;
                                        if (his == 0)
                                            //gis = 2 * xis - 1;
                                            gis = -1;
                                        else
                                            gis = his;
                                        
                                        double git;
                                        if (hit == 0)
                                            //git = 2 * xit - 1;
                                            git = -1;
                                        else
                                            git = hit;
                                        
                                        if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xit - 1 ) > 0 && git * ( 2 * xi - 1 ) > 0 ){
                                            
                                            int s=0;
                                            for (int k=1; k <= z*L; k++){
                                                s += x[r][k-1] * (1<<(L*z-k));
                                            }
                                            
                                            allowed_conf_bp[l][t_ij].push_back(s);
                                            
                                        }
                                        
                                    }
                                    
                                }
                                
                            }
                            
                        }
                    }
                    
                }
                
            }
            
            //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
            
        }
        
    }
    
}


inline void Messages::def_constraint(L4, int flag_red, int z, int n, vector <double>& vec_J_k_to_n){
    
    int xk, xkp, xks, xkt;
    int z2, z3, z4;
    int tn;
    
    L=4;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    //local field acting on i at the fourth time step
    vector <double> ht(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    z4 = 4*z;
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_n[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_n[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_n[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
        ht[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_t = 0;
        for (int k = z3; k < z4; k++){
            xkt = x[r][k];
            ht[r] += vec_J_k_to_n[k_t] * (2 * xkt - 1);
            k_t ++;
        }
        
    }
    
    
    
    int xn = 0;
    for (int xnp = 0; xnp < 2; ++xnp){
        for (int xns = 0; xns < 2; ++xns){
            for (int xnt = 0; xnt < 2; ++xnt){
            
                tn = 8 * xn + 4 * xnp + 2 * xns + xnt;
                
                for (int r = 0; r < Nc; r++){
                    
                    
                    double gn;
                    if ( h[r] == 0)
                        //gn = 2 * xn - 1;
                        gn = -1;
                    else
                        gn = h[r];
                    
                    double gnp;
                    if (hp[r] == 0)
                        //gnp =  2 * xnp - 1;
                        gnp = -1;
                    else
                        gnp = hp[r];
                    
                    double gns;
                    if (hs[r] == 0)
                        //gns =  2 * xns - 1;
                        gns = -1;
                    else
                        gns = hs[r];
                    
                    double gnt;
                    if (ht[r] == 0)
                        //gnt =  2 * xnt - 1;
                        gnt = -1;
                    else
                        gnt = ht[r];
                    
                    
                    
                    if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xnt - 1 ) > 0 && gnt * ( 2 * xn - 1 ) > 0 ){
                        int s=0;
                        for (int k=1; k <= z*L; k++){
                            s += x[r][k-1] * (1<<(L*z-k));
                        }
                        
                        allowed_conf[n][tn].push_back(s);
                        
                    }
                    
                }
                
            }
            
        }
    }
    
    if (flag_red == 0){
        xn = 1;
        for (int xnp = 0; xnp < 2; ++xnp){
            for (int xns = 0; xns < 2; ++xns){
                for (int xnt = 0; xnt < 2; ++xnt){
                    
                    tn = 8 * xn + 4 * xnp + 2 * xns + xnt;
                    
                    for (int r = 0; r < Nc; r++){
                        
                        
                        double gn;
                        if ( h[r] == 0)
                            //gn = 2 * xn - 1;
                            gn = -1;
                        else
                            gn = h[r];
                        
                        double gnp;
                        if (hp[r] == 0)
                            //gnp =  2 * xnp - 1;
                            gnp = -1;
                        else
                            gnp = hp[r];
                        
                        double gns;
                        if (hs[r] == 0)
                            //gns =  2 * xns - 1;
                            gns = -1;
                        else
                            gns = hs[r];
                        
                        double gnt;
                        if (ht[r] == 0)
                            //gnt =  2 * xnt - 1;
                            gnt = -1;
                        else
                            gnt = ht[r];
                        
                        
                        if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xnt - 1 ) > 0 && gnt * ( 2 * xn - 1 ) > 0 ){
                            int s=0;
                            for (int k=1; k <= z*L; k++){
                                s += x[r][k-1] * (1<<(L*z-k));
                            }
                            
                            allowed_conf[n][tn].push_back(s);
                            
                        }
                        
                    }
                    
                }
                
            }
        }
        
    }
    
}




inline void Messages::def_constraint_bp(L1T1, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i) {
    
    int xk, xkp;
    int z2;
    double hi, hip;
    int t_ij;
    
    L=2;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;                                  //vector of 0/1 where we store the configurations of the cavity neighbours
    
    compute_x(Nc,z,x);

    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc,0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    
    
    z2 = 2*z;

    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's and of the xkp's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_i[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's and of the xkp's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_n = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_i[k_n] * (2 * xkp - 1);
            k_n ++;
        }
        
    }
    
    int xi = 0;
    for (int xj = 0; xj < 2; ++xj){
        for (int xip = 0; xip < 2; ++xip){
            for (int xjp = 0; xjp < 2; ++xjp){
                
                t_ij = 8 * xi + 4 * xj + 2 * xip + xjp;
                
                for (int r = 0; r < Nc; r++){
                    //local field acting on i at the first time step
                    hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                    //local field acting on i at the second time step
                    hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                    
                    double gi;
                    if (hi == 0)
                        //gi = 2 * xi - 1;
                        gi = - 1;
                    else
                        gi = hi;
                    
                    double gip;
                    if (hip == 0)
                        //gip = 2 * xip - 1;
                        gip = -1;
                    else
                        gip = hip;
                    
                    if (gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xip - 1 ) > 0){
                            
                        int s=0;
                        for (int k=1; k <= z*L; k++){
                            s += x[r][k-1] * (1<<(L*z-k));
                        }
                        
                        allowed_conf_bp[l][t_ij].push_back(s);
                    }
                    
                    
                    
                }
                
            }
        }
        
    }
    
    if (flag_red == 0){
        
        xi = 1;
        for (int xj = 0; xj < 2; ++xj){
            for (int xip = 0; xip < 2; ++xip){
                for (int xjp = 0; xjp < 2; ++xjp){
                    
                    t_ij = 8 * xi + 4 * xj + 2 * xip + xjp;
                    
                    for (int r = 0; r < Nc; r++){
                        //local field acting on i at the first time step
                        hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                        //local field acting on i at the second time step
                        hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                        
                        double gi;
                        if (hi == 0)
                            //gi = 2 * xi - 1;
                            gi = -1;
                        else
                            gi = hi;
                        
                        double gip;
                        if (hip == 0)
                            //gip = 2 * xip - 1;
                            gip = -1;
                        else
                            gip = hip;
                    
                        if (gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xip - 1 ) > 0){
                            
                            int s=0;
                            for (int k=1; k <= z*L; k++){
                                s += x[r][k-1] * (1<<(L*z-k));
                            }
                        
                            allowed_conf_bp[l][t_ij].push_back(s);
                        }
                    
                    
                        
                        
                    }
                    
                }
            }
            
        }
        
        
    }
    //}
    
}




inline void Messages::def_constraint(L1T1, int flag_red, int z, int n, vector <double>& vec_J_k_to_n) {
    
    
    int xk, xkp;
    int z2;
    int tn;
    
    int L = 2;
    
    // this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc,z,x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc,0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    
    
    z2 = 2*z;
    
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's and of the xkp's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_n[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's and of the xkp's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_n = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_n[k_n] * (2 * xkp - 1);
            k_n ++;
        }
        
    }
    
    //for (int xn = 0; xn < 2; ++xn){
    int xn = 0;
    for (int xnp = 0; xnp < 2; ++xnp){
        
        tn = 2 * xn + xnp;
        
        for (int r = 0; r < Nc; r++){
            
            double gn;
            if ( h[r] == 0)
                //gn = 2 * xn - 1;
                gn = -1;
            else
                gn = h[r];
            
            double gnp;
            if (hp[r] == 0)
                //gnp =  2 * xnp - 1;
                gnp = -1;
            else
                gnp = hp[r];                
            
            if (gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xnp - 1 ) > 0){
                
                int s=0;
                for (int k=1; k <= z*L; k++){
                    s += x[r][k-1] * (1<<(L*z-k));
                }
                
                allowed_conf[n][tn].push_back(s);
                
            }
            
        }
        
    }
    
    if (flag_red == 0){
        
        xn = 1;
        for (int xnp = 0; xnp < 2; ++xnp){
            
            tn = 2 * xn + xnp;
            
            for (int r = 0; r < Nc; r++){
                
                double gn;
                if ( h[r] == 0)
                    //gn = 2 * xn - 1;
                    gn = -1;
                else
                    gn = h[r];
                
                double gnp;
                if (hp[r] == 0)
                    //gnp =  2 * xnp - 1;
                    gnp = -1;
                else
                    gnp = hp[r];
                
                if (gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xnp - 1 ) > 0){
                    
                    int s=0;
                    for (int k=1; k <= z*L; k++){
                        s += x[r][k-1] * (1<<(L*z-k));
                    }
                    
                    allowed_conf[n][tn].push_back(s);
                    
                }
                
            }
            
        }
        
        
    }
    
    
    
    //}
    
}

inline void Messages::def_constraint_bp(L1T2, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i) {
    
    int xk, xkp, xks;
    int z2, z3;
    double hi, hip, his;
    int t_ij;
    
    L=3;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);

    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_i[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_i[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_i[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
        
        
    }
    
    int xi = 0;
    for (int xj = 0; xj < 2; ++xj){
        for (int xip = 0; xip < 2; ++xip){
            for (int xjp = 0; xjp < 2; ++xjp){
                for (int xis = 0; xis < 2; ++xis){
                    for (int xjs = 0; xjs < 2; ++xjs){
                        
                        
                        t_ij = 32 * xi + 16 * xj + 8 * xip + 4 * xjp + 2 * xis + xjs ;
                        
                        
                        for (int r = 0; r < Nc; r++){
                            //local field acting on i at the first time step
                            hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                            //local field acting on i at the second time step
                            hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                            //local field acting on i at the third time step
                            his = J_j_to_i * (2 * xjs - 1) + hs[r];
                            
                            double gi;
                            if (hi == 0)
                                //gi = 2 * xi - 1;
                                gi = -1;
                            else
                                gi = hi;
                            
                            double gip;
                            if (hip == 0)
                                //gip = 2 * xip - 1;
                                gip = -1;
                            else
                                gip = hip;
                            
                            double gis;
                            if (his == 0)
                                //gis = 2 * xis - 1;
                                gis = -1;
                            else
                                gis = his;
                            

                            if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xis - 1 ) > 0 ){
                                
                                int s=0;
                                for (int k=1; k <= z*L; k++){
                                    s += x[r][k-1] * (1<<(L*z-k));
                                }
                                
                                allowed_conf_bp[l][t_ij].push_back(s);
                                
                            }
                            
                        }
                        
                    }
                }
                
            }
        }
        
        //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
        
    }
    
    if (flag_red == 0){
        xi = 1;
        for (int xj = 0; xj < 2; ++xj){
            for (int xip = 0; xip < 2; ++xip){
                for (int xjp = 0; xjp < 2; ++xjp){
                    for (int xis = 0; xis < 2; ++xis){
                        for (int xjs = 0; xjs < 2; ++xjs){
                            
                            
                            t_ij = 32 * xi + 16 * xj + 8 * xip + 4 * xjp + 2 * xis + xjs ;
                            
                            
                            for (int r = 0; r < Nc; r++){
                                //local field acting on i at the first time step
                                hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                                //local field acting on i at the second time step
                                hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                                //local field acting on i at the third time step
                                his = J_j_to_i * (2 * xjs - 1) + hs[r];
                                
                                double gi;
                                if (hi == 0)
                                    //gi = 2 * xi - 1;
                                    gi = -1;
                                else
                                    gi = hi;
                                
                                double gip;
                                if (hip == 0)
                                    //gip = 2 * xip - 1;
                                    gip = -1;
                                else
                                    gip = hip;
                                
                                double gis;
                                if (his == 0)
                                    //gis = 2 * xis - 1;
                                    gis = -1;
                                else
                                    gis = his;
                                
                                
                                if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xis - 1 ) > 0 ){
                                    
                                    int s=0;
                                    for (int k=1; k <= z*L; k++){
                                        s += x[r][k-1] * (1<<(L*z-k));
                                    }
                                    
                                    allowed_conf_bp[l][t_ij].push_back(s);
                                    
                                }
                                
                            }
                            
                        }
                    }
                    
                }
            }
            
            //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
            
        }
        
        
    }
    
}

inline void Messages::def_constraint(L1T2, int flag_red, int z, int n, vector <double>& vec_J_k_to_n) {
    
    int xk, xkp, xks;
    int z2, z3;
    int tn;
    
    L=3;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_n[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_n[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_n[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
    }
    
    int xn = 0;
    for (int xnp = 0; xnp < 2; ++xnp){
        for (int xns = 0; xns < 2; ++xns){
            
            tn = 4 * xn + 2 * xnp + xns;
            
            for (int r = 0; r < Nc; r++){
                
                double gn;
                if ( h[r] == 0)
                    //gn = 2 * xn - 1;
                    gn = -1;
                else
                    gn = h[r];
                
                double gnp;
                if (hp[r] == 0)
                    //gnp =  2 * xnp - 1;
                    gnp = -1;
                else
                    gnp = hp[r];
                
                double gns;
                if (hs[r] == 0)
                    //gns =  2 * xns - 1;
                    gns = -1;
                else
                    gns = hs[r];
                
                
                if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xns - 1 ) > 0 ){
                    int s=0;
                    for (int k=1; k <= z*L; k++){
                        s += x[r][k-1] * (1<<(L*z-k));
                    }
                    
                    allowed_conf[n][tn].push_back(s);
                    
                }
                
            }
            
        }
    }
    
    if (flag_red == 0){
        xn = 1;
        for (int xnp = 0; xnp < 2; ++xnp){
            for (int xns = 0; xns < 2; ++xns){
                
                tn = 4 * xn + 2 * xnp + xns;
                
                for (int r = 0; r < Nc; r++){
                    
                    double gn;
                    if ( h[r] == 0)
                        //gn = 2 * xn - 1;
                        gn = -1;
                    else
                        gn = h[r];
                    
                    double gnp;
                    if (hp[r] == 0)
                        //gnp =  2 * xnp - 1;
                        gnp = -1;
                    else
                        gnp = hp[r];
                    
                    double gns;
                    if (hs[r] == 0)
                        //gns =  2 * xns - 1;
                        gns = -1;
                    else
                        gns = hs[r];
                    
                    
                    if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xns - 1 ) > 0 ){
                        int s=0;
                        for (int k=1; k <= z*L; k++){
                            s += x[r][k-1] * (1<<(L*z-k));
                        }
                        
                        allowed_conf[n][tn].push_back(s);
                        
                    }
                    
                }
                
            }
        }
        
    }
    
}


inline void Messages::def_constraint_bp(L1T3, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i) {
    
    int xk, xkp, xks, xkt;
    int z2, z3, z4;
    double hi, hip, his, hit;
    int t_ij;
    
    L=4;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    //local field acting on i at the fourth time step
    vector <double> ht(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    z4 = 4*z;
    
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_i[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_i[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_i[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
        ht[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_t = 0;
        for (int k = z3; k < z4; k++){
            xkt = x[r][k];
            ht[r] += vec_J_k_to_i[k_t] * (2 * xkt - 1);
            k_t ++;
        }
        
        
        
    }
    
    int xi = 0;
    for (int xj = 0; xj < 2; ++xj){
        for (int xip = 0; xip < 2; ++xip){
            for (int xjp = 0; xjp < 2; ++xjp){
                for (int xis = 0; xis < 2; ++xis){
                    for (int xjs = 0; xjs < 2; ++xjs){
                        for (int xit = 0; xit < 2; ++xit){
                            for (int xjt = 0; xjt < 2; ++xjt){
                        
                                t_ij = 128 * xi + 64 * xj + 32 * xip + 16 * xjp + 8 * xis + 4 * xjs + 2 * xit + xjt;
                                
                                for (int r = 0; r < Nc; r++){
                                    //local field acting on i at the first time step
                                    hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                                    //local field acting on i at the second time step
                                    hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                                    //local field acting on i at the third time step
                                    his = J_j_to_i * (2 * xjs - 1) + hs[r];
                                    //local field acting on i at the fourth time step
                                    hit = J_j_to_i * (2 * xjt - 1) + ht[r];
                                    
                                    
                                    double gi;
                                    if (hi == 0)
                                        //gi = 2 * xi - 1;
                                        gi = -1;
                                    else
                                        gi = hi;
                                    
                                    double gip;
                                    if (hip == 0)
                                        //gip = 2 * xip - 1;
                                        gip = -1;
                                    else
                                        gip = hip;
                                    
                                    double gis;
                                    if (his == 0)
                                        //gis = 2 * xis - 1;
                                        gis = -1;
                                    else
                                        gis = his;
                                    
                                    double git;
                                    if (hit == 0)
                                        //git = 2 * xit - 1;
                                        git = -1;
                                    else
                                        git = hit;
                                    
                                    if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xit - 1 ) > 0 && git * ( 2 * xit - 1 ) > 0 ){
                                        
                                        int s=0;
                                        for (int k=1; k <= z*L; k++){
                                            s += x[r][k-1] * (1<<(L*z-k));
                                        }
                                        
                                        allowed_conf_bp[l][t_ij].push_back(s);
                                        
                                    }
                                    
                                }
                                
                            }
                            
                        }
                        
                    }
                }
                
            }
        }
        
        //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
        
    }
    
    if (flag_red == 0){
        xi = 1;
        for (int xj = 0; xj < 2; ++xj){
            for (int xip = 0; xip < 2; ++xip){
                for (int xjp = 0; xjp < 2; ++xjp){
                    for (int xis = 0; xis < 2; ++xis){
                        for (int xjs = 0; xjs < 2; ++xjs){
                            for (int xit = 0; xit < 2; ++xit){
                                for (int xjt = 0; xjt < 2; ++xjt){
                                    
                                    t_ij = 128 * xi + 64 * xj + 32 * xip + 16 * xjp + 8 * xis + 4 * xjs + 2 * xit + xjt;
                                    
                                    for (int r = 0; r < Nc; r++){
                                        //local field acting on i at the first time step
                                        hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                                        //local field acting on i at the second time step
                                        hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                                        //local field acting on i at the third time step
                                        his = J_j_to_i * (2 * xjs - 1) + hs[r];
                                        //local field acting on i at the fourth time step
                                        hit = J_j_to_i * (2 * xjt - 1) + ht[r];
                                        
                                        double gi;
                                        if (hi == 0)
                                            //gi = 2 * xi - 1;
                                            gi = -1;
                                        else
                                            gi = hi;
                                        
                                        double gip;
                                        if (hip == 0)
                                            //gip = 2 * xip - 1;
                                            gip = -1;
                                        else
                                            gip = hip;
                                        
                                        double gis;
                                        if (his == 0)
                                            //gis = 2 * xis - 1;
                                            gis = -1;
                                        else
                                            gis = his;
                                        
                                        double git;
                                        if (hit == 0)
                                            //git = 2 * xit - 1;
                                            git = -1;
                                        else
                                            git = hit;
                                        
                                        if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xit - 1 ) > 0 && git * ( 2 * xit - 1 ) > 0 ){
                                            
                                            int s=0;
                                            for (int k=1; k <= z*L; k++){
                                                s += x[r][k-1] * (1<<(L*z-k));
                                            }
                                            
                                            allowed_conf_bp[l][t_ij].push_back(s);
                                            
                                        }
                                        
                                    }
                                    
                                }
                                
                            }
                            
                        }
                    }
                    
                }
                
            }
            
            //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
            
        }
        
    }
    
}

inline void Messages::def_constraint(L1T3, int flag_red, int z, int n, vector <double>& vec_J_k_to_n){
    
    int xk, xkp, xks, xkt;
    int z2, z3, z4;
    int tn;
    
    L=4;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    //local field acting on i at the fourth time step
    vector <double> ht(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    z4 = 4*z;
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_n[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_n[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_n[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
        ht[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_t = 0;
        for (int k = z3; k < z4; k++){
            xkt = x[r][k];
            ht[r] += vec_J_k_to_n[k_t] * (2 * xkt - 1);
            k_t ++;
        }
        
    }
    
    
    
    int xn = 0;
    for (int xnp = 0; xnp < 2; ++xnp){
        for (int xns = 0; xns < 2; ++xns){
            for (int xnt = 0; xnt < 2; ++xnt){
            
                tn = 8 * xn + 4 * xnp + 2 * xns + xnt;
                
                for (int r = 0; r < Nc; r++){
                    
                    
                    double gn;
                    if ( h[r] == 0)
                        //gn = 2 * xn - 1;
                        gn = -1;
                    else
                        gn = h[r];
                    
                    double gnp;
                    if (hp[r] == 0)
                        //gnp =  2 * xnp - 1;
                        gnp = -1;
                    else
                        gnp = hp[r];
                    
                    double gns;
                    if (hs[r] == 0)
                        //gns =  2 * xns - 1;
                        gns = -1;
                    else
                        gns = hs[r];
                    
                    double gnt;
                    if (ht[r] == 0)
                        //gnt =  2 * xnt - 1;
                        gnt = -1;
                    else
                        gnt = ht[r];
                    
                    
                    
                    if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xnt - 1 ) > 0 && gnt * ( 2 * xnt - 1 ) > 0 ){
                        int s=0;
                        for (int k=1; k <= z*L; k++){
                            s += x[r][k-1] * (1<<(L*z-k));
                        }
                        
                        allowed_conf[n][tn].push_back(s);
                        
                    }
                    
                }
                
            }
            
        }
    }
    
    if (flag_red == 0){
        xn = 1;
        for (int xnp = 0; xnp < 2; ++xnp){
            for (int xns = 0; xns < 2; ++xns){
                for (int xnt = 0; xnt < 2; ++xnt){
                    
                    tn = 8 * xn + 4 * xnp + 2 * xns + xnt;
                    
                    for (int r = 0; r < Nc; r++){
                        
                        
                        double gn;
                        if ( h[r] == 0)
                            //gn = 2 * xn - 1;
                            gn = -1;
                        else
                            gn = h[r];
                        
                        double gnp;
                        if (hp[r] == 0)
                            //gnp =  2 * xnp - 1;
                            gnp = -1;
                        else
                            gnp = hp[r];
                        
                        double gns;
                        if (hs[r] == 0)
                            //gns =  2 * xns - 1;
                            gns = -1;
                        else
                            gns = hs[r];
                        
                        double gnt;
                        if (ht[r] == 0)
                            //gnt =  2 * xnt - 1;
                            gnt = -1;
                        else
                            gnt = ht[r];
                        
                        
                        if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xnt - 1 ) > 0 && gnt * ( 2 * xnt - 1 ) > 0 ){
                            int s=0;
                            for (int k=1; k <= z*L; k++){
                                s += x[r][k-1] * (1<<(L*z-k));
                            }
                            
                            allowed_conf[n][tn].push_back(s);
                            
                        }
                        
                    }
                    
                }
                
            }
        }
        
    }
    
}

inline void Messages::def_constraint_bp(L2T1, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i) {
    
    int xk, xkp, xks;
    int z2, z3;
    double hi, hip, his;
    int t_ij;
    
    L=3;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);

    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_i[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_i[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_i[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
        
        
    }
    
    int xi = 0;
    for (int xj = 0; xj < 2; ++xj){
        for (int xip = 0; xip < 2; ++xip){
            for (int xjp = 0; xjp < 2; ++xjp){
                for (int xis = 0; xis < 2; ++xis){
                    for (int xjs = 0; xjs < 2; ++xjs){
                        
                        
                        t_ij = 32 * xi + 16 * xj + 8 * xip + 4 * xjp + 2 * xis + xjs ;
                        
                        
                        for (int r = 0; r < Nc; r++){
                            //local field acting on i at the first time step
                            hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                            //local field acting on i at the second time step
                            hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                            //local field acting on i at the third time step
                            his = J_j_to_i * (2 * xjs - 1) + hs[r];
                            
                            double gi;
                            if (hi == 0)
                                //gi = 2 * xi - 1;
                                gi = -1;
                            else
                                gi = hi;
                            
                            double gip;
                            if (hip == 0)
                                //gip = 2 * xip - 1;
                                gip = -1;
                            else
                                gip = hip;
                            
                            double gis;
                            if (his == 0)
                                //gis = 2 * xis - 1;
                                gis = -1;
                            else
                                gis = his;
                            

                            if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xip - 1 ) > 0 ){
                                
                                int s=0;
                                for (int k=1; k <= z*L; k++){
                                    s += x[r][k-1] * (1<<(L*z-k));
                                }
                                
                                allowed_conf_bp[l][t_ij].push_back(s);
                                
                            }
                            
                        }
                        
                    }
                }
                
            }
        }
        
        //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
        
    }
    
    if (flag_red == 0){
        xi = 1;
        for (int xj = 0; xj < 2; ++xj){
            for (int xip = 0; xip < 2; ++xip){
                for (int xjp = 0; xjp < 2; ++xjp){
                    for (int xis = 0; xis < 2; ++xis){
                        for (int xjs = 0; xjs < 2; ++xjs){
                            
                            
                            t_ij = 32 * xi + 16 * xj + 8 * xip + 4 * xjp + 2 * xis + xjs ;
                            
                            
                            for (int r = 0; r < Nc; r++){
                                //local field acting on i at the first time step
                                hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                                //local field acting on i at the second time step
                                hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                                //local field acting on i at the third time step
                                his = J_j_to_i * (2 * xjs - 1) + hs[r];
                                
                                double gi;
                                if (hi == 0)
                                    //gi = 2 * xi - 1;
                                    gi = -1;
                                else
                                    gi = hi;
                                
                                double gip;
                                if (hip == 0)
                                    //gip = 2 * xip - 1;
                                    gip = -1;
                                else
                                    gip = hip;
                                
                                double gis;
                                if (his == 0)
                                    //gis = 2 * xis - 1;
                                    gis = -1;
                                else
                                    gis = his;
                                
                                
                                if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xip - 1 ) > 0 ){
                                    
                                    int s=0;
                                    for (int k=1; k <= z*L; k++){
                                        s += x[r][k-1] * (1<<(L*z-k));
                                    }
                                    
                                    allowed_conf_bp[l][t_ij].push_back(s);
                                    
                                }
                                
                            }
                            
                        }
                    }
                    
                }
            }
            
            //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
            
        }
        
        
    }
}

inline void Messages::def_constraint(L2T1, int flag_red, int z, int n, vector <double>& vec_J_k_to_n){
    
    int xk, xkp, xks;
    int z2, z3;
    int tn;
    
    L=3;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_n[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_n[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_n[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
    }
    
    int xn = 0;
    for (int xnp = 0; xnp < 2; ++xnp){
        for (int xns = 0; xns < 2; ++xns){
            
            tn = 4 * xn + 2 * xnp + xns;
            
            for (int r = 0; r < Nc; r++){
                
                double gn;
                if ( h[r] == 0)
                    //gn = 2 * xn - 1;
                    gn = 1;
                else
                    gn = h[r];
                
                double gnp;
                if (hp[r] == 0)
                    //gnp =  2 * xnp - 1;
                    gnp = -1;
                else
                    gnp = hp[r];
                
                double gns;
                if (hs[r] == 0)
                    //gns =  2 * xns - 1;
                    gns = -1;
                else
                    gns = hs[r];
                
                
                if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xnp - 1 ) > 0 ){
                    int s=0;
                    for (int k=1; k <= z*L; k++){
                        s += x[r][k-1] * (1<<(L*z-k));
                    }
                    
                    allowed_conf[n][tn].push_back(s);
                    
                }
                
            }
            
        }
    }
    
    if (flag_red == 0){
        xn = 1;
        for (int xnp = 0; xnp < 2; ++xnp){
            for (int xns = 0; xns < 2; ++xns){
                
                tn = 4 * xn + 2 * xnp + xns;
                
                for (int r = 0; r < Nc; r++){
                    
                    double gn;
                    if ( h[r] == 0)
                        //gn = 2 * xn - 1;
                        gn = -1;
                    else
                        gn = h[r];
                    
                    double gnp;
                    if (hp[r] == 0)
                        //gnp =  2 * xnp - 1;
                        gnp = -1;
                    else
                        gnp = hp[r];
                    
                    double gns;
                    if (hs[r] == 0)
                        //gns =  2 * xns - 1;
                        gns = -1;
                    else
                        gns = hs[r];
                    
                    
                    if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xnp - 1 ) > 0 ){
                        int s=0;
                        for (int k=1; k <= z*L; k++){
                            s += x[r][k-1] * (1<<(L*z-k));
                        }
                        
                        allowed_conf[n][tn].push_back(s);
                        
                    }
                    
                }
                
            }
        }
        
    }
    
}

inline void Messages::def_constraint_bp(L2T2, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i) {
    
    int xk, xkp, xks, xkt;
    int z2, z3, z4;
    double hi, hip, his, hit;
    int t_ij;
    
    L=4;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    //local field acting on i at the fourth time step
    vector <double> ht(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    z4 = 4*z;
    
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_i[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_i[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_i[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
        ht[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_t = 0;
        for (int k = z3; k < z4; k++){
            xkt = x[r][k];
            ht[r] += vec_J_k_to_i[k_t] * (2 * xkt - 1);
            k_t ++;
        }
        
        
        
    }
    
    int xi = 0;
    for (int xj = 0; xj < 2; ++xj){
        for (int xip = 0; xip < 2; ++xip){
            for (int xjp = 0; xjp < 2; ++xjp){
                for (int xis = 0; xis < 2; ++xis){
                    for (int xjs = 0; xjs < 2; ++xjs){
                        for (int xit = 0; xit < 2; ++xit){
                            for (int xjt = 0; xjt < 2; ++xjt){
                        
                                t_ij = 128 * xi + 64 * xj + 32 * xip + 16 * xjp + 8 * xis + 4 * xjs + 2 * xit + xjt;
                                
                                for (int r = 0; r < Nc; r++){
                                    //local field acting on i at the first time step
                                    hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                                    //local field acting on i at the second time step
                                    hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                                    //local field acting on i at the third time step
                                    his = J_j_to_i * (2 * xjs - 1) + hs[r];
                                    //local field acting on i at the fourth time step
                                    hit = J_j_to_i * (2 * xjt - 1) + ht[r];
                                    
                                    
                                    double gi;
                                    if (hi == 0)
                                        //gi = 2 * xi - 1;
                                        gi = -1;
                                    else
                                        gi = hi;
                                    
                                    double gip;
                                    if (hip == 0)
                                        //gip = 2 * xip - 1;
                                        gip = -1;
                                    else
                                        gip = hip;
                                    
                                    double gis;
                                    if (his == 0)
                                        //gis = 2 * xis - 1;
                                        gis = -1;
                                    else
                                        gis = his;
                                    
                                    double git;
                                    if (hit == 0)
                                        //git = 2 * xit - 1;
                                        git = -1;
                                    else
                                        git = hit;
                                    
                                    if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xit - 1 ) > 0 && git * ( 2 * xis - 1 ) > 0 ){
                                        
                                        int s=0;
                                        for (int k=1; k <= z*L; k++){
                                            s += x[r][k-1] * (1<<(L*z-k));
                                        }
                                        
                                        allowed_conf_bp[l][t_ij].push_back(s);
                                        
                                    }
                                    
                                }
                                
                            }
                            
                        }
                        
                    }
                }
                
            }
        }
        
        //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
        
    }
    
    if (flag_red == 0){
        xi = 1;
        for (int xj = 0; xj < 2; ++xj){
            for (int xip = 0; xip < 2; ++xip){
                for (int xjp = 0; xjp < 2; ++xjp){
                    for (int xis = 0; xis < 2; ++xis){
                        for (int xjs = 0; xjs < 2; ++xjs){
                            for (int xit = 0; xit < 2; ++xit){
                                for (int xjt = 0; xjt < 2; ++xjt){
                                    
                                    t_ij = 128 * xi + 64 * xj + 32 * xip + 16 * xjp + 8 * xis + 4 * xjs + 2 * xit + xjt;
                                    
                                    for (int r = 0; r < Nc; r++){
                                        //local field acting on i at the first time step
                                        hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                                        //local field acting on i at the second time step
                                        hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                                        //local field acting on i at the third time step
                                        his = J_j_to_i * (2 * xjs - 1) + hs[r];
                                        //local field acting on i at the fourth time step
                                        hit = J_j_to_i * (2 * xjt - 1) + ht[r];
                                        
                                        double gi;
                                        if (hi == 0)
                                            //gi = 2 * xi - 1;
                                            gi = -1;
                                        else
                                            gi = hi;
                                        
                                        double gip;
                                        if (hip == 0)
                                            //gip = 2 * xip - 1;
                                            gip = -1;
                                        else
                                            gip = hip;
                                        
                                        double gis;
                                        if (his == 0)
                                            //gis = 2 * xis - 1;
                                            gis = -1;
                                        else
                                            gis = his;
                                        
                                        double git;
                                        if (hit == 0)
                                            //git = 2 * xit - 1;
                                            git = -1;
                                        else
                                            git = hit;
                                        
                                        if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xit - 1 ) > 0 && git * ( 2 * xis - 1 ) > 0 ){
                                            
                                            int s=0;
                                            for (int k=1; k <= z*L; k++){
                                                s += x[r][k-1] * (1<<(L*z-k));
                                            }
                                            
                                            allowed_conf_bp[l][t_ij].push_back(s);
                                            
                                        }
                                        
                                    }
                                    
                                }
                                
                            }
                            
                        }
                    }
                    
                }
                
            }
            
            //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
            
        }
        
    }
    
}

inline void Messages::def_constraint(L2T2, int flag_red, int z, int n, vector <double>& vec_J_k_to_n){
    
    int xk, xkp, xks, xkt;
    int z2, z3, z4;
    int tn;
    
    L=4;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    //local field acting on i at the fourth time step
    vector <double> ht(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    z4 = 4*z;
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_n[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_n[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_n[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
        ht[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_t = 0;
        for (int k = z3; k < z4; k++){
            xkt = x[r][k];
            ht[r] += vec_J_k_to_n[k_t] * (2 * xkt - 1);
            k_t ++;
        }
        
    }
    
    
    
    int xn = 0;
    for (int xnp = 0; xnp < 2; ++xnp){
        for (int xns = 0; xns < 2; ++xns){
            for (int xnt = 0; xnt < 2; ++xnt){
            
                tn = 8 * xn + 4 * xnp + 2 * xns + xnt;
                
                for (int r = 0; r < Nc; r++){
                    
                    
                    double gn;
                    if ( h[r] == 0)
                        //gn = 2 * xn - 1;
                        gn = -1;
                    else
                        gn = h[r];
                    
                    double gnp;
                    if (hp[r] == 0)
                        //gnp =  2 * xnp - 1;
                        gnp = -1;
                    else
                        gnp = hp[r];
                    
                    double gns;
                    if (hs[r] == 0)
                        //gns =  2 * xns - 1;
                        gns = -1;
                    else
                        gns = hs[r];
                    
                    double gnt;
                    if (ht[r] == 0)
                        //gnt =  2 * xnt - 1;
                        gnt = -1;
                    else
                        gnt = ht[r];
                    
                    
                    
                    if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xnt - 1 ) > 0 && gnt * ( 2 * xns - 1 ) > 0 ){
                        int s=0;
                        for (int k=1; k <= z*L; k++){
                            s += x[r][k-1] * (1<<(L*z-k));
                        }
                        
                        allowed_conf[n][tn].push_back(s);
                        
                    }
                    
                }
                
            }
            
        }
    }
    
    if (flag_red == 0){
        xn = 1;
        for (int xnp = 0; xnp < 2; ++xnp){
            for (int xns = 0; xns < 2; ++xns){
                for (int xnt = 0; xnt < 2; ++xnt){
                    
                    tn = 8 * xn + 4 * xnp + 2 * xns + xnt;
                    
                    for (int r = 0; r < Nc; r++){
                        
                        
                        double gn;
                        if ( h[r] == 0)
                            //gn = 2 * xn - 1;
                            gn = -1;
                        else
                            gn = h[r];
                        
                        double gnp;
                        if (hp[r] == 0)
                            //gnp =  2 * xnp - 1;
                            gnp = -1;
                        else
                            gnp = hp[r];
                        
                        double gns;
                        if (hs[r] == 0)
                            //gns =  2 * xns - 1;
                            gns = -1;
                        else
                            gns = hs[r];
                        
                        double gnt;
                        if (ht[r] == 0)
                            //gnt =  2 * xnt - 1;
                            gnt = -1;
                        else
                            gnt = ht[r];
                        
                        
                        if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xnt - 1 ) > 0 && gnt * ( 2 * xns - 1 ) > 0 ){
                            int s=0;
                            for (int k=1; k <= z*L; k++){
                                s += x[r][k-1] * (1<<(L*z-k));
                            }
                            
                            allowed_conf[n][tn].push_back(s);
                            
                        }
                        
                    }
                    
                }
                
            }
        }
        
    }
    
}

inline void Messages::def_constraint_bp(L2T3, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i) {
    
    int xk, xkp, xks, xkt, xtu;
    int z2, z3, z4, z5;
    double hi, hip, his, hit, hiu;
    int t_ij;
    
    L=5;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    //local field acting on i at the fourth time step
    vector <double> ht(Nc,0.);
    //local field acting on i at the fifth time step
    vector <double> hu(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    z4 = 4*z;
    z5 = 5*z;
    
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's and xks's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_i[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_i[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_i[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
        ht[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the fourth time step
        
        int k_t = 0;
        for (int k = z3; k < z4; k++){
            xkt = x[r][k];
            ht[r] += vec_J_k_to_i[k_t] * (2 * xkt - 1);
            k_t ++;
        }

        hu[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's, xkt's and xtu's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the fifth time step
        
        int k_u = 0;
        for (int k = z4; k < z5; k++){
            xtu = x[r][k];
            hu[r] += vec_J_k_to_i[k_u] * (2 * xtu - 1);
            k_u ++;
        }
        
        
        
    }
    
    int xi = 0;
    for (int xj = 0; xj < 2; ++xj){
        for (int xip = 0; xip < 2; ++xip){
            for (int xjp = 0; xjp < 2; ++xjp){
                for (int xis = 0; xis < 2; ++xis){
                    for (int xjs = 0; xjs < 2; ++xjs){
                        for (int xit = 0; xit < 2; ++xit){
                            for (int xjt = 0; xjt < 2; ++xjt){
                                for (int xiu = 0; xiu < 2; ++xiu){
                                    for (int xju = 0; xju < 2; ++xju){
                        
                                        t_ij = 512 * xi + 256 * xj + 128 * xip + 64 * xjp + 32 * xis + 16 * xjs + 8 * xit + 4 * xjt + 2 * xiu + xju;
                                
                                        for (int r = 0; r < Nc; r++){
                                            //local field acting on i at the first time step
                                            hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                                            //local field acting on i at the second time step
                                            hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                                            //local field acting on i at the third time step
                                            his = J_j_to_i * (2 * xjs - 1) + hs[r];
                                            //local field acting on i at the fourth time step
                                            hit = J_j_to_i * (2 * xjt - 1) + ht[r];
                                            //local field acting on i at the fifth time step
                                            hiu = J_j_to_i * (2 * xju - 1) + hu[r];
                                    
                                    
                                            double gi;
                                            if (hi == 0)
                                                //gi = 2 * xi - 1;
                                                gi = -1;
                                            else
                                                gi = hi;
                                    
                                            double gip;
                                            if (hip == 0)
                                                //gip = 2 * xip - 1;
                                                gip = -1;
                                            else
                                                gip = hip;
                                    
                                            double gis;
                                            if (his == 0)
                                                //gis = 2 * xis - 1;
                                                gis = -1;
                                            else
                                                gis = his;
                                    
                                            double git;
                                            if (hit == 0)
                                                //git = 2 * xit - 1;
                                                git = -1;
                                            else
                                                git = hit;

                                            double giu;
                                            if (hiu == 0)
                                                //giu = 2 * xiu - 1;
                                                giu = -1;
                                            else
                                                giu = hiu;
                                    
                                            if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xit - 1 ) > 0 && git * ( 2 * xiu - 1 ) > 0 && giu * ( 2 * xit - 1 ) > 0 ){
                                        
                                                int s=0;
                                                for (int k=1; k <= z*L; k++){
                                                    s += x[r][k-1] * (1<<(L*z-k));
                                                }
                                        
                                                allowed_conf_bp[l][t_ij].push_back(s);
                                        
                                            }
                                    
                                        }
                                
                                    }
                            
                                }
                        
                            }
                        }
                
                    }
                }
            }
        }
        
        //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
        
    }
    
    if (flag_red == 0){
        xi = 1;
        for (int xj = 0; xj < 2; ++xj){
            for (int xip = 0; xip < 2; ++xip){
                for (int xjp = 0; xjp < 2; ++xjp){
                    for (int xis = 0; xis < 2; ++xis){
                        for (int xjs = 0; xjs < 2; ++xjs){
                            for (int xit = 0; xit < 2; ++xit){
                                for (int xjt = 0; xjt < 2; ++xjt){
                                    for (int xiu = 0; xiu < 2; ++xiu){
                                        for (int xju = 0; xju < 2; ++xju){
                        
                                            t_ij = 512 * xi + 256 * xj + 128 * xip + 64 * xjp + 32 * xis + 16 * xjs + 8 * xit + 4 * xjt + 2 * xiu + xju;
                                
                                            for (int r = 0; r < Nc; r++){
                                                //local field acting on i at the first time step
                                                hi  = J_j_to_i * (2 * xj - 1)  + h[r];
                                                //local field acting on i at the second time step
                                                hip = J_j_to_i * (2 * xjp - 1) + hp[r];
                                                //local field acting on i at the third time step
                                                his = J_j_to_i * (2 * xjs - 1) + hs[r];
                                                //local field acting on i at the fourth time step
                                                hit = J_j_to_i * (2 * xjt - 1) + ht[r];
                                                //local field acting on i at the fifth time step
                                                hiu = J_j_to_i * (2 * xju - 1) + hu[r];
                                    
                                    
                                                double gi;
                                                if (hi == 0)
                                                    //gi = 2 * xi - 1;
                                                    gi = -1;
                                                else
                                                    gi = hi;
                                    
                                                double gip;
                                                if (hip == 0)
                                                    //gip = 2 * xip - 1;
                                                    gip = -1;
                                                else
                                                    gip = hip;
                                    
                                                double gis;
                                               if (his == 0)
                                                    //gis = 2 * xis - 1;
                                                    gis = -1;
                                                else
                                                    gis = his;
                                    
                                                double git;
                                                if (hit == 0)
                                                    //git = 2 * xit - 1;
                                                    git = -1;
                                                else
                                                    git = hit;

                                                double giu;
                                                if (hiu == 0)
                                                    //giu = 2 * xiu - 1;
                                                    giu = -1;
                                                else
                                                    giu = hiu;
                                    
                                                if ( gi * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xis - 1 ) > 0 && gis * ( 2 * xit - 1 ) > 0 && git * ( 2 * xiu - 1 ) > 0 && giu * ( 2 * xit - 1 ) > 0 ){
                                            
                                                    int s=0;
                                                    for (int k=1; k <= z*L; k++){
                                                        s += x[r][k-1] * (1<<(L*z-k));
                                                    }
                                        
                                                    allowed_conf_bp[l][t_ij].push_back(s);
                                        
                                                }
                                    
                                            }
                                
                                        }
                            
                                    }
                        
                                }
                            }
                
                        }
                    }
                }
            }
        
            //if (it_l->l == 62) cout << "mess ( t_ij = " << t_ij << " ) = " << mess_ij[t_ij] << endl;
        
        }
        
    }
    
}

inline void Messages::def_constraint(L2T3, int flag_red, int z, int n, vector <double>& vec_J_k_to_n){
    
    int xk, xkp, xks, xkt, xtu;
    int z2, z3, z4, z5;
    int tn;
    
    L=5;
    
    //this is the number of combinations of the spins k's.
    int Nc = (int) pow(pow(2.,L),z);
    
    vector < vector <int> > x;
    
    compute_x(Nc, z, x);
    
    //these vectors have a value for each configurations of the spins k's
    //local field acting on i
    vector <double> h(Nc, 0.);
    //local field acting on i at the second time step
    vector <double> hp(Nc,0.);
    //local field acting on i at the third time step
    vector <double> hs(Nc,0.);
    //local field acting on i at the fourth time step
    vector <double> ht(Nc,0.);
    //local field acting on i at the fifth time step
    vector <double> hu(Nc,0.);
    
    z2 = 2*z;
    z3 = 3*z;
    z4 = 4*z;
    z5 = 5*z;
    
    for (int r = 0; r < Nc; r++){
        
        h[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's  xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the first time step
        
        for (int k = 0; k < z; k++){
            xk = x[r][k];
            h[r] += vec_J_k_to_n[k] * (2 * xk - 1);
        }
        
        hp[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the second time step
        
        int k_p = 0;
        for (int k = z; k < z2; k++){
            xkp = x[r][k];
            hp[r] += vec_J_k_to_n[k_p] * (2 * xkp - 1);
            k_p ++;
        }
        
        hs[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the third time step
        
        int k_s = 0;
        for (int k = z2; k < z3; k++){
            xks = x[r][k];
            hs[r] += vec_J_k_to_n[k_s] * (2 * xks - 1);
            k_s ++;
        }
        
        ht[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's and xkt's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the fourth time step
        
        int k_t = 0;
        for (int k = z3; k < z4; k++){
            xkt = x[r][k];
            ht[r] += vec_J_k_to_n[k_t] * (2 * xkt - 1);
            k_t ++;
        }

        hu[r] = 0;
        //for each value of r (i.e. for each combination of the xk's  xkp's xks's xkt's and xtu's)
        //we compute the value of the contribution to the local field given by the neighbours of i except j at the fifth time step
        
        int k_u = 0;
        for (int k = z4; k < z5; k++){
            xtu = x[r][k];
            hu[r] += vec_J_k_to_n[k_u] * (2 * xtu - 1);
            k_u ++;
        }
        
    }
    
    
    
    int xn = 0;
    for (int xnp = 0; xnp < 2; ++xnp){
        for (int xns = 0; xns < 2; ++xns){
            for (int xnt = 0; xnt < 2; ++xnt){
                for (int xnu = 0; xnu < 2; ++xnu){
            
                    tn = 16 * xn + 8 * xnp + 4 * xns + 2 * xnt + xnu;
                
                    for (int r = 0; r < Nc; r++){
                    
                    
                        double gn;
                        if ( h[r] == 0)
                            //gn = 2 * xn - 1;
                            gn = -1;
                        else
                            gn = h[r];
                    
                        double gnp;
                        if (hp[r] == 0)
                            //gnp =  2 * xnp - 1;
                            gnp = -1;
                        else
                            gnp = hp[r];
                    
                        double gns;
                        if (hs[r] == 0)
                            //gns =  2 * xns - 1;
                            gns = -1;
                        else
                            gns = hs[r];
                    
                        double gnt;
                        if (ht[r] == 0)
                            //gnt =  2 * xnt - 1;
                            gnt = -1;
                        else
                            gnt = ht[r];

                        double gnu;
                        if (hu[r] == 0)
                            //gnu =  2 * xnu - 1;
                            gnu = -1;
                        else
                            gnu = hu[r];
                    
                    
                    
                        if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xnt - 1 ) > 0 && gnt * ( 2 * xnu - 1 ) > 0 && gnu * ( 2 * xnt - 1 ) > 0){
                            int s=0;
                            for (int k=1; k <= z*L; k++){
                                s += x[r][k-1] * (1<<(L*z-k));
                            }
                        
                            allowed_conf[n][tn].push_back(s);
                        
                        }
                    
                    }
                
                }
            
            }
        }
    }
    
    if (flag_red == 0){
        xn = 1;
        for (int xnp = 0; xnp < 2; ++xnp){
            for (int xns = 0; xns < 2; ++xns){
                for (int xnt = 0; xnt < 2; ++xnt){
                    for (int xnu = 0; xnu < 2; ++xnu){
            
                        tn = 16 * xn + 8 * xnp + 4 * xns + 2 * xnt + xnu;
                
                        for (int r = 0; r < Nc; r++){
                    
                    
                            double gn;
                            if ( h[r] == 0)
                                //gn = 2 * xn - 1;
                                gn = -1;
                            else
                                gn = h[r];
                    
                            double gnp;
                            if (hp[r] == 0)
                                //gnp =  2 * xnp - 1;
                                gnp = -1;
                            else
                                gnp = hp[r];
                    
                            double gns;
                            if (hs[r] == 0)
                                //gns =  2 * xns - 1;
                                gns = -1;
                            else
                                gns = hs[r];
                    
                            double gnt;
                            if (ht[r] == 0)
                                //gnt =  2 * xnt - 1;
                                gnt = -1;
                            else
                                gnt = ht[r];

                            double gnu;
                            if (hu[r] == 0)
                                //gnu =  2 * xnu - 1;
                                gnu = -1;
                            else
                                gnu = hu[r];
                    
                    
                    
                            if ( gn * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xns - 1 ) > 0 && gns * ( 2 * xnt - 1 ) > 0 && gnt * ( 2 * xnu - 1 ) > 0 && gnu * ( 2 * xnt - 1 ) > 0){
                                int s=0;
                                for (int k=1; k <= z*L; k++){
                                    s += x[r][k-1] * (1<<(L*z-k));
                                }
                        
                                allowed_conf[n][tn].push_back(s);
                        
                            }
                    
                        }
                
                    }
            
                }
            }
        }
        
    }
    
}

//------------------------------------------------------------------------------ COMPUTE VECTOR OF CONFS --------------------------------------------------------------------------------//

void Messages::compute_x(int Nc, int z, vector < vector <int> >& x){
    
    //in the BP case:
    
    //for L = 1
    //in each of these combinations, x_k's may assume 0 or 1 values. we store these values in the matrix x
    //e.g. let's consider the case when i has 2 other neighbours besides j: in this case z=2 and Nc=2^2=4 and x will be the matrix
    //0 0
    //1 0
    //0 1
    //1 1
    
    //for L = 2
    //in each of these combinations, (x_k, x_kp) may assume values (00), (01), (10), (11). we store these values on the matrix x
    //e.g. let's consider the case when i has 2 other neighbours besides j: in this case z=2 and Nc=4^2=16 and x will be the matrix
    //xk1 xk2 xk1p xk2p
    //  0   0    0    0
    //  1   0    0    0
    //  0   1    0    0
    //  1   1    0    0
    //  etc...
    //  1   1    1    1
    
    //in the case when we need to compute the total marginals, the same matrix can be computed except that it refers to the total number of neighbours
    
    vector < vector <int> > t_x(Nc, vector <int>( z * L, 0));
    
    x = t_x;
    
    //cout << "Printing z and Nc: " << z << " " << Nc << endl;
    
    for (int r = 0; r < Nc; r++){
        
        //cout << "configuration " << r << " : " << endl;
        for (int k = 0; k < z * L; k++){
            x[r][k] = ( ( r & (int)pow(2.,k) ) >> k ) ;
            //cout << x[r][k] << " ";
        }
        //cout << endl;
        
    }
    
}


bool Messages::messNormalize(){
    
    long double sum;
    
    int f = 1;

    int i,j;
    
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); ++it_l){
        
        
        i = G.E[it_l->l].v_node[0];
        j = G.E[it_l->l].v_node[1];
        
        vector<int>::iterator it = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), j);
        int index_j = distance (G.v[i].v_neigh.begin(), it);
        
     

        
        sum = 0.;
        for (int t=0; t<Q; ++t){
            sum += update_mess[i][index_j][t];
        }
        
        if (sum !=0 ){
            for (int t=0; t<Q; ++t){
                update_mess[i][index_j][t]/=sum;
            }
        }
  
        /*
        if (it_l->l == 62){
          for (int t=0; t<4; ++t){
              cout << "********************" << update_mess[i][index_j][t] << endl;
          }
        }
        */
        
        else
            f=0;
        
    }
    

    return f;
    
};


void Messages::superMarginals(){
    
    int i,j;
    double sum;
    
    int l = 0;
    
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); it_l+=2){
        
        sum = 0.;
        vector <long double> marg(Q,0.);
        
        i = G.E[it_l->l].v_node[0];
        j = G.E[it_l->l].v_node[1];
        

        vector<int>::iterator it_j = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), j);
        int index_j = distance (G.v[i].v_neigh.begin(), it_j);
        
        vector<int>::iterator it_i = find(G.v[j].v_neigh.begin(), G.v[j].v_neigh.end(), i);
        int index_i = distance (G.v[j].v_neigh.begin(), it_i);
        
        //cout << "marginal of the SuperNode at link " << it_l->l << ":" << endl;
        
        for (int k=0; k<Q; k++){
            marg[k] = update_mess[i][index_j][k] * update_mess[j][index_i][t_reduced[sw[k]]];
            sum += marg[k];
            
        }
        
        for (int k=0; k<Q; k++){
            marg[k] = marg[k] / sum ;
        }
        
        super_marginal[l] = marg;
        
        l++;
        
    }
    
};


vector<int> Messages::Switch(){
    
    int Nc2 = Q;
    
    vector < vector <int> > y(Nc2, vector <int>( 2 * L, 0));
    vector <int> cp(Nc2, 0);

    for (int r = 0; r < Nc2; r++){
        
        //cout << "configuration " << r << " : " << endl;
        for (int k = 0; k < 2 * L; k++){
            y[r][k] = ( ( r & (int)pow(2.,k) ) >> k ) ;
            //cout << y[r][k] << " ";
        }
        //cout << endl;
        
        vector <int> v( 2 * L, 0);
        
        for (int i=0; i < 2 * L; i+=2){
            v[i]   = y[r][i+1];
            v[i+1] = y[r][i];
        }
        
        int s = 0;
        
        for (int i=0; i < 2 * L ; i++)
            s += pow(2,i) * v[i];
        
        cp[r] = s;
        
    }
    
    return cp;
}



void Messages::messState(){

    cout << endl;
    cout << "---*---*---*---printing messages---*---*---*---" << endl;

    int i,j;
    
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); ++it_l){
        i = G.E[it_l->l].v_node[0];
        j = G.E[it_l->l].v_node[1];
        
        vector<int>::iterator it = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), j);
        int index_j = distance (G.v[i].v_neigh.begin(), it);
        
        cout << "message from node " << i << " to node " << j << endl;
        for(int s = 0; s < Q; s ++)
            cout << mess[i][index_j][s] << " ";
        
        cout << endl;
    }
    

};



void Messages::nodeMarginals(){
  
    node_marginal.resize(N);
    int n = 0;
    int i,j, l;
    
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); ++it_l){
        
        if (it_l->l % 2 == 0)
            l = (int) (it_l->l / 2);
        else
            continue;
        
        cout << "LINK: " << it_l->l << endl;
        
        if (n<N){
        
            i = G.E[it_l->l].v_node[0];
            j = G.E[it_l->l].v_node[1];
        
            cout << " connecting nodes " << i << " " << j << endl;

            if(!node_marginal[i].size()){
                node_marginal[i].resize(2);
                node_marginal[i]=link_marginal[l][0];
                n++;
            }
        
            if(!node_marginal[j].size()){
                node_marginal[j].resize(2);
                node_marginal[j]=link_marginal[l][1];
                n++;
            }
            
            
    
        }
        else
            break;
        
    }
    
};



template<typename Tag>
void Messages::Wrap_computeLastMarg(){
    
    long double x0, x1;
       
    int l = 0;
    
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); ++it_l){
        
        if (it_l->l % 2 == 0)
            l = (int) (it_l->l / 2);
        else
            continue;
        
        //at link l, link_marginal[l] contains two distributions, one for each of the spins that are connected from l.
        //for the first spin, i, the marginal in contained in link_marginal[l][0]
        //for the second spin, j, the marginal in contained in link_marginal[l][1]
        
        int i = G.E[it_l->l].v_node[0];
        int j = G.E[it_l->l].v_node[1];
        
        
        vector <long double> marginals = computeMarginalsLastTime<Tag>(l);
                
        link_marginal[l][0][0] = marginals[0];
        link_marginal[l][0][1] = 1 - marginals[0];

        link_marginal[l][1][0] = marginals[1];
        link_marginal[l][1][1] = 1 - marginals[1];

    }
}



inline vector <long double> Messages::def_computeMarginalsLastTime(L1, int l) {
    
    
    double sum1 = 0.;

    vector<long double> v_link_marginal(2, 0.);
    
    int xi = 0;
    for (int xj = 0; xj < 2; ++xj){
        
        int t_ij = 2 * xi + xj;
        sum1 += super_marginal[l][t_ij];
    }
    
    v_link_marginal[0] = sum1;
    
    sum1 = 0.;
    
    int xj = 0;
    for (int xi = 0; xi < 2; ++xi){
        
        int t_ij = 2 * xi + xj;
        sum1 += super_marginal[l][t_ij];
    }
    
    v_link_marginal[1] = sum1;

        
    return v_link_marginal;
    
}

inline vector <long double> Messages::def_computeMarginalsLastTime(L1T1, int l) {
    
    
    double sum1 = 0.;

    vector<long double> v_link_marginal(2, 0.);
    
    int xi = 0;
    for (int xj = 0; xj < 2; ++xj){
        for(int xip = 0; xip < 2; ++xip){
            for(int xjp = 0; xjp < 2; ++xjp){
                int t_ij = 8 * xip + 4 * xjp + 2 * xi + xj;
                sum1 += super_marginal[l][t_ij];
            }
        }
    }
    
    v_link_marginal[0] = sum1;
    
    sum1 = 0.;
    
    int xj = 0;
    for (int xi = 0; xi < 2; ++xi){
        for(int xip = 0; xip < 2; ++xip){
            for(int xjp = 0; xjp < 2; ++xjp){
                int t_ij = 8 * xip + 4 * xjp + 2 * xi + xj;
                sum1 += super_marginal[l][t_ij];
            }
        }
    }
    
    v_link_marginal[1] = sum1;

        
    return v_link_marginal;
    
}



inline vector <long double> Messages::def_computeMarginalsLastTime(L2, int l) {
    
    double sum1 = 0.;

    vector<long double> v_link_marginal(2, 0.);

    
    int xip = 0;
    for (int xi = 0; xi < 2; ++xi){
        for (int xj = 0; xj < 2; ++xj){
            for (int xjp = 0; xjp < 2; ++xjp){
                
                int t_ij = 8 * xi + 4 * xj + 2 * xip + xjp;

                sum1 += super_marginal[l][t_ij];
                
            }
        }
    }

    v_link_marginal[0] = sum1;
    
    sum1=0;
    
    int xjp = 0;
    for (int xi = 0; xi < 2; ++xi){
        for (int xj = 0; xj < 2; ++xj){
            for (int xip = 0; xip < 2; ++xip){
            
                int t_ij = 8 * xi + 4 * xj + 2 * xip + xjp;

                sum1 += super_marginal[l][t_ij];
     
            }
        }
    }
    
    v_link_marginal[1] = sum1;

    
    return v_link_marginal;
    
    
}



void Messages::updateAndStore(){

    mess = update_mess;
    prev_super_marginal = super_marginal;

};



void Messages::linkMarginalState(){

    cout << endl;
    cout << "---*---*---*---printing link marginals---*---*---*---" << endl;

    int l;
    
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); ++it_l){
        
        if (it_l->l % 2 == 0)
            l = (int) (it_l->l / 2);
        else
            continue;
        
        cout << "marginal computed at link " << it_l->l << " of node " << G.E[it_l->l].v_node[0] << ":" << endl;
        cout <<  link_marginal[l][0][0] << " " << link_marginal[l][0][1] << endl;

        cout << "marginal computed at link " << it_l->l << " of node " << G.E[it_l->l].v_node[1] << ":" << endl;
        cout <<  link_marginal[l][1][0] << " " << link_marginal[l][1][1] << endl;

    }
    
};



void Messages::nodeMarginalState(){
    
    cout << endl;
    cout << "---*---*---*---printing node marginals---*---*---*---" << endl;
    
    for(vector<Node>::iterator it_i = G.v.begin(); it_i != G.v.end(); ++it_i){
        
        cout << "marginal of node " << it_i-> n << ":" << endl;

        //the following if is needed to deal with isolated nodes,

        if(node_marginal[it_i->n].size()) 
	        cout <<  node_marginal[it_i->n][0] << " " << node_marginal[it_i->n][1] << endl;
        else
		cout <<  1 << " " << 0 << endl; 
        
    }
    
    cout << "---*---*---*---printing node marginals of sigma = 1---*---*---*---" << endl;
    for(vector<Node>::iterator it_i = G.v.begin(); it_i != G.v.end(); ++it_i){
	if(node_marginal[it_i->n].size())
        	cout << node_marginal[it_i->n][1] << endl;
        else
		cout << 0 << endl;
    }
    
    
};



double Messages::compareSuperMarginals(){
    
    long double tmp, max = 0.;
    
    int L = super_marginal.size();
    
    for (int l=0; l < L; l++){
        
        for (int k = 0; k < Q -1; k++){
            tmp = abs(super_marginal[l][k] - prev_super_marginal[l][k]);
            if (tmp > max)
                max = tmp;
        }
    }
    
    return max;
    
};


void Messages::superMarginalsState(){
    
    int i,j;
    long double x0, x1;
    
    int l = 0;
    
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); it_l+=2){

        i = G.E[l].v_node[0];
        j = G.E[l].v_node[1];
        
        cout << "marginal of the SuperNode at link " << l << " between spins " << i << " " << j << endl;
        
        for (int s = 0; s < Q; s ++)
            cout << super_marginal[l][s] << " ";
        cout << endl;
        
        l++;
        
    }
    
    
};




//------------------------------------------------------------------------------------- BP iteration ---------------------------------------------------------------------------------//


    
    


template<typename Tag>
void Messages::BPiteration(double th, int flag_red, int flag_approx, int T, bool verbose){
    
    int    t = 0;
    double tmp_th = 1;

    setQ(flag_red);

    look_up_table_bp<Tag>(flag_red);
    look_up_table<Tag>(flag_red);

    //MM is the number of configurations used to compute the BP update in the approximate case
    int MM = 10000;
    
    while (tmp_th > th && t < T){
        //the problem with this approximate algo is that the random number generation is particularly slow
        if(flag_approx)
            messUpdate_approx(MM);
        else
            messUpdate();


        messNormalize();
        superMarginals();
        
        tmp_th = compareSuperMarginals();
        
        updateAndStore();
        
        if(verbose){
            messState();
            cout << endl;
            cout << "BP iteration: at time t=" << t << " the maximum error between current and previous super_marginals is " << tmp_th << endl;
        }
        
        t++;
        
    }
    
    print_BPit<Tag>(tmp_th, t);
    
    if(verbose){
        cout << "the final node_marginals are " << endl;
        cout << endl;
        superMarginalsState();
        cout << endl;
    }
    

    logPartitionFunction<Tag>();
     
    
};


//------------------------------------------------------------------------------------- COUNTING LOOPS ---------------------------------------------------------------------------------//


template <typename Tag>
void Messages::logPartitionFunction(){
    
    
    int i, j, n;
    
    double Part1 = 0.;
    double Part2 = 0.;
    double sum;
    double logz_ij, logz_hn;
    
    vector <double> t(Q,0.);
    
    //***************************** computation of z_ij (S_link) *****************************//
    for (vector<Link>::iterator it_l = G.E.begin() ; it_l != G.E.end(); ++it_l){
        i = G.E[it_l->l].v_node[0];
        j = G.E[it_l->l].v_node[1];
        
        //if there is a link whose J = 0, this link is absent and thus we do not need to sum over it
        if (G.E[it_l->l].J != 0){
            
            vector<int>::iterator it_j = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), j);
            int index_j = distance (G.v[i].v_neigh.begin(), it_j);
            
            //cout << "link " << it_l->l << endl;
            
            //cout << mess[i][index_j][0] << " " << mess[i][index_j][1] << " " << mess[i][index_j][2] << " " << mess[i][index_j][3] << endl;
            
            vector<int>::iterator it_i = find(G.v[j].v_neigh.begin(), G.v[j].v_neigh.end(), i);
            int index_i = distance (G.v[j].v_neigh.begin(), it_i);
            
            //cout << mess[j][index_i][0] << " " << mess[j][index_i][1] << " " << mess[j][index_i][2] << " " << mess[j][index_i][3] << endl;
            
            //t_ij = 2 xi + xj
            
            sum = 0.;
            for (int k=0; k<Q; k++){
                t[k] = mess[i][index_j][k] * mess[j][index_i][t_reduced[sw[k]]];
                sum += t[k];
                
                /*
                 if(it_l->l == 62)
                 cout << k << " " << sw[k] << " " << t_reduced[sw[k]] << " " << mess[i][index_j][k] << " " << mess[j][index_i][t_reduced[sw[k]]] << " " << t[k] << " " << sum <<  endl;
                 */
                
            }
            
            
            logz_ij = log(factor_S_link * sum);
            
        }
        else
            logz_ij = 0;
        
        
        
        Part1 += logz_ij;
    }
    

    
    //***************************** computation of z_hn (S_node) *****************************//
    for (vector<Node>::iterator it_n = G.v.begin() ; it_n != G.v.end(); ++it_n){
        
        n = G.v[it_n->n].n;
        
        vector <vector <long double> > in_mess;
        
        //couplings from node k's to node n.
        vector <double> vec_J_k_to_n;
        
        //for each link k->n
        
        int z =0;
        
        
        for(vector<int>::iterator it_k = G.v[n].v_neigh.begin(); it_k != G.v[n].v_neigh.end(); ++it_k){
            
            vector<int>::iterator it = find(G.v[*it_k].v_neigh.begin(), G.v[*it_k].v_neigh.end(), n);
            int index_n = distance (G.v[*it_k].v_neigh.begin(), it);
            
            //we store the mess_k_to_n
            in_mess.push_back(mess[*it_k][index_n]);
            
            
            it = find(G.v[n].v_neigh.begin(), G.v[n].v_neigh.end(), *it_k);
            int index_k = distance (G.v[n].v_neigh.begin(), it);
            
            int index_J_k_to_n = G.v[n].v_link[index_k];
            
            vec_J_k_to_n.push_back(G.E[index_J_k_to_n].J);
            z++;
            
        }
        
        
        int cnt = 0;
        for (int i = 0; i < vec_J_k_to_n.size(); i++ ){
            if (vec_J_k_to_n[i] == 0)
                cnt ++;
            else
                break;
        }
        
        double sum;
        
        
        if( cnt == vec_J_k_to_n.size() ){
            sum = 1;
        }
        
        else{
            
            int xn, xk;
            int t_kn, t_kn_tmp;
            int NC;
            
            sum = 0.;
            
            double prod_mess;
            
            
            
            for (int tn = 0; tn < QS; tn++){
                
                NC = allowed_conf[n][tn].size();
                
                /*
                 if (n==50)
                 cout << "N allowed conf " << tn << " " << NC << endl;
                 */
                
                t_kn =0;
                
                //get the values of xi at t=0,1..,L
                for(int t=L; t>=1; t--){
                    xn = (tn >> (t-1) ) & 1;
                    t_kn += (1<<(2*t-2)) * xn;
                }
                
                t_kn_tmp = t_kn;
                
                
                for (int rr = 0; rr < NC; rr++){
                    
                    r = allowed_conf[n][tn][rr];
                    
                    //for each value of r (i.e. for each combination of the allowed xk's)
                    //we store the product of messages mess_k_to_i in prod_mess.
                    prod_mess = 1.;
                    
                    for (int k = 0; k < z; k++){
                        
                        
                        t_kn = t_kn_tmp;
                        
                        //get the allowed xk at t=0,1..,L
                        for(int t=L; t>=1; t--){
                            xk = ( r >> (t*z - 1 - k) ) & 1;
                            t_kn += (1<<(2*t-1)) * xk;
                        }
                        
                        prod_mess *= in_mess[k][t_reduced[t_kn]] * factor_S_node_1;
                        
                        /*
                         if (n==50){
                         cout << 2./(1<<z) << " " <<  k << " " << t_kn << " " << t_reduced[t_kn] << " " << in_mess[k][t_reduced[t_kn]] << endl;
                         }
                         */
                        
                    }
                    
                    sum += factor_S_node_2 * prod_mess;
                    
                }
                
                
                
            }
            
        }
        
        
        logz_hn = log(sum);
        
        Part2 += logz_hn;
    }
    
    cout << "Countings loops: " << endl;
    
    cout << "S_link: "  << 0.5 * Part1 << endl;
    cout << "S_node: " << Part2 << endl;

    int print_file = 0;         //set print_file = 1/0 to print/not print the values of complexities on file .dat
    
    printEntropy<Tag>(Part1, Part2, print_file);

}


void vec_print(vector<int>& vec){
    for (int i=0; i<vec.size(); i++)
        cout << vec[i] << ' ';
    cout << endl;
}

//------------------------------------------------------------------------------------- PRINT FUNCTIONS ---------------------------------------------------------------------------------//

void Messages::def_printEntropy(L1, double Part1, double Part2, int print_file){

    
    if (std::isnan( (- 0.5 * Part1 + Part2) ) ){
        cout << "Entropy-L1: " << -1000 << endl;
        cout << "Number of cycles of length 1: " << 0 << endl;
        cout << "NUMBER-L1 " << 0 << endl;
        cout << "Entropy-L1 " << -1000 << endl;
    }
    else{
        cout << "Entropy-L1: " << (- 0.5 * Part1 + Part2) / N << endl;
        cout << "Number of cycles of length 1: " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "NUMBER-L1 " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "Entropy-L1 " << (- 0.5 * Part1 + Part2) / N << endl;
    }

    if(print_file == 1){
        ofstream writer("L1.dat" , ios::app);
        writer << (- 0.5 * Part1 + Part2) / N << endl;
    }  
}


void Messages::def_printEntropy(L2, double Part1, double Part2, int print_file){

    
    if (std::isnan( (- 0.5 * Part1 + Part2) ) ){
        cout << "Entropy-L2: " << -1000 << endl;
        cout << "Number of cycles of length 2: " << 0 << endl;
        cout << "NUMBER-L2 " << 0 << endl;
        cout << "Entropy-L2 " << -1000 << endl;
    }
    else{
        cout << "Entropy-L2: " << (- 0.5 * Part1 + Part2) / N << endl;
        cout << "Number of cycles of length 2: " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "NUMBER-L2 " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "Entropy-L2 " << (- 0.5 * Part1 + Part2) / N << endl;
    }

    if(print_file == 1){
        ofstream writer("L2.dat" , ios::app);
        writer << (- 0.5 * Part1 + Part2) / N << endl;
    }  
    
}


void Messages::def_printEntropy(L3, double Part1, double Part2, int print_file){
    
    if (std::isnan( (- 0.5 * Part1 + Part2) ) ){
        cout << "Entropy-L3: " << -1000 << endl;
        cout << "Number of cycles of length 3: " << 0 << endl;
        cout << "NUMBER-L3 " << 0 << endl;
        cout << "Entropy-L3 " << -1000 << endl;
    }
    else{
        cout << "Entropy-L3: " << (- 0.5 * Part1 + Part2) / N << endl;
        cout << "Number of cycles of length 3: " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "NUMBER-L3 " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "Entropy-L3 " << (- 0.5 * Part1 + Part2) / N << endl;
    }

    if(print_file == 1){
        ofstream writer("L3.dat" , ios::app);
        writer << (- 0.5 * Part1 + Part2) / N << endl;
    } 
}


void Messages::def_printEntropy(L4, double Part1, double Part2, int print_file){
    
    if (std::isnan( (- 0.5 * Part1 + Part2) ) ){
        cout << "Entropy-L4: " << -1000 << endl;
        cout << "Number of cycles of length 4: " << 0 << endl;
        cout << "NUMBER-L4 " << 0 << endl;
        cout << "Entropy-L4 " << -1000 << endl;
    }
    else{
        cout << "Entropy-L4: " << (- 0.5 * Part1 + Part2) / N << endl;
        cout << "Number of cycles of length 4: " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "NUMBER-L4 " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "Entropy-L4 " << (- 0.5 * Part1 + Part2) / N << endl;
    }

    if(print_file == 1){
        ofstream writer("L4.dat" , ios::app);
        writer << (- 0.5 * Part1 + Part2) / N << endl;
    } 
    
}


void Messages::def_printEntropy(L1T1, double Part1, double Part2, int print_file){
    
    if (std::isnan( (- 0.5 * Part1 + Part2) ) ){
        cout << "Entropy-L1T1: " << -1000 << endl;
        cout << "Number of confs evolving in cycles of L=1 in T=1 steps : " << 0 << endl;
        cout << "NUMBER-L1T1 " << 0 << endl;
        cout << "Entropy-L1T1 " << -1000 << endl;
    }
    else{
        cout << "Entropy-L1T1: " << (- 0.5 * Part1 + Part2) / N << endl;
        cout << "Number of confs evolving in cycles of L=1 in T=1 steps : " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "NUMBER-L1T1 " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "Entropy-L1T1 " << (- 0.5 * Part1 + Part2) / N << endl;
    }

    if(print_file == 1){
        ofstream writer("L1T1.dat" , ios::app);
        writer << (- 0.5 * Part1 + Part2) / N << endl;
    } 
    
}

void Messages::def_printEntropy(L1T2, double Part1, double Part2, int print_file){
    
    if (std::isnan( (- 0.5 * Part1 + Part2) ) ){
        cout << "Entropy-L1T2: " << -1000 << endl;
        cout << "Number of confs evolving in cycles of L=1 in T=2 steps : " << 0 << endl;
        cout << "NUMBER-L1T2 " << 0 << endl;
        cout << "Entropy-L1T2 " << -1000 << endl;
    }
    else{
        cout << "Entropy-L1T2: " << (- 0.5 * Part1 + Part2) / N << endl;
        cout << "Number of confs evolving in cycles of L=1 in T=2 steps : " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "NUMBER-L1T2 " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "Entropy-L1T2 " << (- 0.5 * Part1 + Part2) / N << endl;
    }

    if(print_file == 1){
        ofstream writer("L1T2.dat" , ios::app);
        writer << (- 0.5 * Part1 + Part2) / N << endl;
    } 
    
};

void Messages::def_printEntropy(L1T3, double Part1, double Part2, int print_file){
    
    if (std::isnan( (- 0.5 * Part1 + Part2) ) ){
        cout << "Entropy-L1T3: " << -1000 << endl;
        cout << "Number of confs evolving in cycles of L=1 in T=3 steps : " << 0 << endl;
        cout << "NUMBER-L1T3 " << 0 << endl;
        cout << "Entropy-L1T3 " << -1000 << endl;
    }
    else{
        cout << "Entropy-L1T3: " << (- 0.5 * Part1 + Part2) / N << endl;
        cout << "Number of confs evolving in cycles of L=1 in T=3 steps : " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "NUMBER-L1T3 " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "Entropy-L1T3 " << (- 0.5 * Part1 + Part2) / N << endl;
    }

    if(print_file == 1){
        ofstream writer("L1T3.dat" , ios::app);
        writer << (- 0.5 * Part1 + Part2) / N << endl;
    } 
    
};

void Messages::def_printEntropy(L2T1, double Part1, double Part2, int print_file){
    
    if (std::isnan( (- 0.5 * Part1 + Part2) ) ){
        cout << "Entropy-L2T1: " << -1000 << endl;
        cout << "Number of confs evolving in cycles of L=2 in T=1 steps : " << 0 << endl;
        cout << "NUMBER-L2T1 " << 0 << endl;
        cout << "Entropy-L2T1 " << -1000 << endl;
    }
    else{
        cout << "Entropy-L2T1: " << (- 0.5 * Part1 + Part2) / N << endl;
        cout << "Number of confs evolving in cycles of L=2 in T=1 steps : " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "NUMBER-L2T1 " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "Entropy-L2T1 " << (- 0.5 * Part1 + Part2) / N << endl;
    }

    if(print_file == 1){
        ofstream writer("L2T1.dat" , ios::app);
        writer << (- 0.5 * Part1 + Part2) / N << endl;
    } 
    
};

void Messages::def_printEntropy(L2T2, double Part1, double Part2, int print_file){
    
    if (std::isnan( (- 0.5 * Part1 + Part2) ) ){
        cout << "Entropy-L2T2: " << -1000 << endl;
        cout << "Number of confs evolving in cycles of L=2 in T=2 steps : " << 0 << endl;
        cout << "NUMBER-L2T2 " << 0 << endl;
        cout << "Entropy-L2T2 " << -1000 << endl;
    }
    else{
        cout << "Entropy-L2T2: " << (- 0.5 * Part1 + Part2) / N << endl;
        cout << "Number of confs evolving in cycles of L=2 in T=2 steps : " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "NUMBER-L2T2 " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "Entropy-L2T2 " << (- 0.5 * Part1 + Part2) / N << endl;
    }

    if(print_file == 1){
        ofstream writer("L2T2.dat" , ios::app);
        writer << (- 0.5 * Part1 + Part2) / N << endl;
    }
    
};

void Messages::def_printEntropy(L2T3, double Part1, double Part2, int print_file){
    
    if (std::isnan( (- 0.5 * Part1 + Part2) ) ){
        cout << "Entropy-L2T3: " << -1000 << endl;
        cout << "Number of confs evolving in cycles of L=2 in T=3 steps : " << 0 << endl;
        cout << "NUMBER-L2T3 " << 0 << endl;
        cout << "Entropy-L2T3 " << -1000 << endl;
    }
    else{
        cout << "Entropy-L2T3: " << (- 0.5 * Part1 + Part2) / N << endl;
        cout << "Number of confs evolving in cycles of L=2 in T=3 steps : " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "NUMBER-L2T3 " << exp (- 0.5 * Part1 + Part2) << endl;
        cout << "Entropy-L2T3 " << (- 0.5 * Part1 + Part2) / N << endl;
    }

    if(print_file == 1){
        ofstream writer("L2T3.dat" , ios::app);
        writer << (- 0.5 * Part1 + Part2) / N << endl;
    }
    
};


void Messages::def_print_BPit(L1, double tmp_th, int t){
    
    cout << endl;
    cout << "BP iteration stopped with an error equal to " << tmp_th << endl;
    cout << "ERROR-1 " << tmp_th << endl;
    cout << "BP TIME STEPS-L1 " << t << endl;
    
}


void Messages::def_print_BPit(L2, double tmp_th, int t){
    
    cout << endl;
    cout << "BP iteration stopped with an error equal to " << tmp_th << endl;
    cout << "ERROR-2 " << tmp_th << endl;
    cout << "BP TIME STEPS-L2 " << t << endl;

}


void Messages::def_print_BPit(L3, double tmp_th, int t){
    
    cout << endl;
    cout << "BP iteration stopped with an error equal to " << tmp_th << endl;
    cout << "ERROR-3 " << tmp_th << endl;
    cout << "BP TIME STEPS-L3 " << t << endl;

}


void Messages::def_print_BPit(L4, double tmp_th, int t){
    
    cout << endl;
    cout << "BP iteration stopped with an error equal to " << tmp_th << endl;
    cout << "ERROR-4 " << tmp_th << endl;
    cout << "BP TIME STEPS-L4 " << t << endl;

}


void Messages::def_print_BPit(L1T1, double tmp_th, int t){
    
    cout << endl;
    cout << "BP iteration stopped with an error equal to " << tmp_th << endl;
    cout << "ERROR-L1T1 " << tmp_th << endl;
    cout << "BP TIME STEPS-L1T1 " << t << endl;

}

void Messages::def_print_BPit(L1T2, double tmp_th, int t){
    
    cout << endl;
    cout << "BP iteration stopped with an error equal to " << tmp_th << endl;
    cout << "ERROR-L2T2 " << tmp_th << endl;
    cout << "BP TIME STEPS-L2T2 " << t << endl;

}

void Messages::def_print_BPit(L1T3, double tmp_th, int t){
    
    cout << endl;
    cout << "BP iteration stopped with an error equal to " << tmp_th << endl;
    cout << "ERROR-L1T3 " << tmp_th << endl;
    cout << "BP TIME STEPS-L1T3 " << t << endl;

}

void Messages::def_print_BPit(L2T1, double tmp_th, int t){
    
    cout << endl;
    cout << "BP iteration stopped with an error equal to " << tmp_th << endl;
    cout << "ERROR-L2T1 " << tmp_th << endl;
    cout << "BP TIME STEPS-L2T1 " << t << endl;

}

void Messages::def_print_BPit(L2T2, double tmp_th, int t){
    
    cout << endl;
    cout << "BP iteration stopped with an error equal to " << tmp_th << endl;
    cout << "ERROR-L2T2 " << tmp_th << endl;
    cout << "BP TIME STEPS-L2T2 " << t << endl;

}

void Messages::def_print_BPit(L2T3, double tmp_th, int t){
    
    cout << endl;
    cout << "BP iteration stopped with an error equal to " << tmp_th << endl;
    cout << "ERROR-L2T3 " << tmp_th << endl;
    cout << "BP TIME STEPS-L2T3 " << t << endl;

}

//these are template methods and needs to be specified for the linker

template void Messages::look_up_table_bp<L1>(int);
template void Messages::look_up_table<L1>(int);

template void Messages::look_up_table_bp<L2>(int);
template void Messages::look_up_table<L2>(int);

template void Messages::look_up_table_bp<L3>(int);
template void Messages::look_up_table<L3>(int);

template void Messages::look_up_table_bp<L4>(int);
template void Messages::look_up_table<L4>(int);

template void Messages::look_up_table_bp<L1T1>(int);
template void Messages::look_up_table<L1T1>(int);

template void Messages::look_up_table_bp<L1T2>(int);
template void Messages::look_up_table<L1T2>(int);

template void Messages::look_up_table_bp<L1T3>(int);
template void Messages::look_up_table<L1T3>(int);

template void Messages::look_up_table_bp<L2T1>(int);
template void Messages::look_up_table<L2T1>(int);

template void Messages::look_up_table_bp<L2T2>(int);
template void Messages::look_up_table<L2T2>(int);

template void Messages::look_up_table_bp<L2T3>(int);
template void Messages::look_up_table<L2T3>(int);

template void Messages::BPiteration<L1>(double, int, int, int, bool);
template void Messages::BPiteration<L2>(double, int, int, int, bool);
template void Messages::BPiteration<L3>(double, int, int, int, bool);
template void Messages::BPiteration<L4>(double, int, int, int, bool);

template void Messages::BPiteration<L1T1>(double, int, int, int, bool);
template void Messages::BPiteration<L1T2>(double, int, int, int, bool);
template void Messages::BPiteration<L1T3>(double, int, int, int, bool);
template void Messages::BPiteration<L2T1>(double, int, int, int, bool);
template void Messages::BPiteration<L2T2>(double, int, int, int, bool);
template void Messages::BPiteration<L2T3>(double, int, int, int, bool);

template void Messages::Wrap_computeLastMarg<L1>();
template void Messages::Wrap_computeLastMarg<L1T1>();
template void Messages::Wrap_computeLastMarg<L2>();

template void Messages::logPartitionFunction<L1>();
template void Messages::logPartitionFunction<L2>();
template void Messages::logPartitionFunction<L3>();
template void Messages::logPartitionFunction<L4>();

template void Messages::logPartitionFunction<L1T1>();
template void Messages::logPartitionFunction<L1T2>();
template void Messages::logPartitionFunction<L1T3>();
template void Messages::logPartitionFunction<L2T1>();
template void Messages::logPartitionFunction<L2T2>();
template void Messages::logPartitionFunction<L2T3>();
