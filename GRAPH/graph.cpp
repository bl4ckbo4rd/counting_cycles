#include "graph.h"



//---------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Node
//---------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------------//



int Node::numberOfLinks(){
  return v_link.size();
};



//---------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------------//     methods of class Factor
//---------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------------//



//input variables:
//p_l : index of the link

Link::Link(int p_l, double p_J) : l ( p_l ), J (p_J) {
    v_node.reserve(2);
};



//---------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Graph
//---------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------------//



Graph::Graph(int p_N, int p_seed_g) : N(p_N), seed_g(p_seed_g) {
    
    v.reserve(N);
    for (int i = 0; i < N; ++i)
        v.push_back (Node(i));
    
};



int Graph::numberOfTotalLinks(){
    return E.size();
};



int Graph::numberOfTotalNodes(){
    return v.size();
};



void Graph::neighsOfNode(int i){
    cout << "node " << i << " has " << v[i].numberOfLinks() << " neighbours: " << endl;
    for (vector<int>::iterator it = v[i].v_neigh.begin() ; it != v[i].v_neigh.end(); ++it)
        cout << *it << endl;
    
    cout << "and " << v[i].numberOfLinks() << " links: " << endl;
    for (vector<int>::iterator it = v[i].v_link.begin() ; it != v[i].v_link.end(); ++it)
        cout << *it << endl;
    
};



void Graph::nodesOfLink(int l){
    cout << "link " << l << " connects the two nodes : " << endl;
    for (vector<int>::iterator it = E[l].v_node.begin() ; it != E[l].v_node.end(); ++it)
        cout << *it << endl;
};





int Graph::addLink(int p_l, vector<int> v_dl, double J){
    
    vector <int> v1 = v_dl;
    
    int flag=1;
    
    //before adding a link we check that the link does not already exist.
    //if the link already exists, flag is set to 0.
    //to this aim it is sufficient to check that the first node of v_dl does not appear in a link with the other node of v_dl.
    
    for(vector<int>::iterator it_l = v[v_dl[0]].v_link.begin(); it_l != v[v_dl[0]].v_link.end(); it_l++){
        
        vector <int> v2 = E[*it_l].v_node;
        
        if (v1 == v2){
            flag=0;
        }

    }
    
    if(flag){
        
        //this is the J that n1 feels from n2: J (n2->n1)
        Link l(p_l,J);
        
        for (int i=0; i<2; i++)
            l.v_node.push_back(v_dl[i]);
        
        E.push_back(l);
        
        v[v_dl[0]].v_link.push_back(p_l);
        v[v_dl[0]].v_neigh.push_back(v_dl[1]);
        
    }
    
    return flag;
    
};



int Graph::addLinkWrapper(int i, int j, double Jitoj, double Jjtoi){
    
    int flag_ij, flag_ji;

    int l = numberOfTotalLinks();
            
    vector <int> vitoj = make_vector<int>() << j << i;
    flag_ij = addLink(l, vitoj, Jitoj);
    l++;

    vector <int> vjtoi = make_vector<int>() << i << j;
    flag_ji = addLink(l, vjtoi, Jjtoi);
    l++;
    
    //if link already exists, overwrite the values of J's with the arguments of the function
    
    int index_link;
    int k = 0;
    for (vector<int>::iterator it = v[i].v_neigh.begin() ; it != v[i].v_neigh.end(); ++it){
        if (*it==j){
            index_link = v[i].v_link[k];
            E[index_link].J = Jjtoi;
        }
        k++;
    }
    
    k = 0;
    for (vector<int>::iterator it = v[j].v_neigh.begin() ; it != v[j].v_neigh.end(); ++it){
        if (*it==i){
            index_link = v[j].v_link[k];
            E[index_link].J = Jitoj;
        }
        k++;
    }
    
    return flag_ij * flag_ji;
}

int Graph::removeLink(int l){
    
    int flag = 0;
    //check that I haven't already deleted the link l
    
    cout << "vorrei eliminare " << l << " " << E[l].J << endl;
    
    if ( E[l].J != 0){
        
        /*
        int i = E[l].v_node[0];
        int j = E[l].v_node[1];
        
        vector<int>::iterator it = find(v[i].v_neigh.begin(), v[i].v_neigh.end(), j);
        int index_j = distance (v[i].v_neigh.begin(), it);
        v[i].v_neigh.erase(v[i].v_neigh.begin() + index_j);
        v[i].v_link.erase(v[i].v_link.begin() + index_j);
        
        it = find(v[j].v_neigh.begin(), v[j].v_neigh.end(), i);
        int index_i = distance (v[j].v_neigh.begin(), it);
        v[j].v_neigh.erase(v[j].v_neigh.begin()+index_i);
        v[j].v_link.erase(v[j].v_link.begin()+index_i);
        */
        
        //E[l].v_node.clear();
        E[l].J = 0;
        
        if (l%2 == 0){
            //E[l+1].v_node.clear();
            E[l+1].J = 0;
        }
        else{
            //E[l-1].v_node.clear();
            E[l-1].J = 0;
        }
        
        flag = 1;
        
    }
    
    return flag;
    
};


void Graph::graphStructure(){
    
    cout << endl;
    cout << endl;
    
    cout << "structure of the graph: " << endl;
    cout << endl;
    
    for (int i = 0; i < N; i ++)
        neighsOfNode(i);
    
    cout << endl;
    cout << endl;

    
    for (vector<Link>::iterator it_l = E.begin() ; it_l != E.end(); ++it_l){
        cout << "coupling on link " << it_l->l << ": " << E[it_l->l].J << endl;
        nodesOfLink(it_l->l);
    }
        
    cout << endl;
    cout << endl;

};



// Method to print connected components in an
// undirected graph
void Graph::connectedComponents(){
    // Mark all the vertices as not visited
    vector<bool> visited(N,false);
    
    int k=0;
    
    for (int i=0; i<N; i++)
    {
        if (visited[i] == false)
        {
            // print all reachable vertices
            // from v
            
            vector<int> t;
            component.push_back(t);
            
            cout << "component " << k << endl;
            DFSUtil(i, visited, k);
            cout << "\n";
            k ++;
        }
    }
    
    for (int s=0; s<k; s++)
        cout << "component's size " << component[s].size() << endl;
    
}


void Graph::DFSUtil(int i, vector<bool>& visited, int k){
    // Mark the current node as visited and print it
    visited[i] = true;
    cout << i << " ";
    
    component[k].push_back(i);
    
    // Recur for all the vertices
    // adjacent to this vertex
    for(vector<int>::iterator it_j = v[i].v_neigh.begin(); it_j != v[i].v_neigh.end(); ++it_j){
        
        vector<int>::iterator it = find(v[i].v_neigh.begin(), v[i].v_neigh.end(), *it_j);
        int index_j = distance (v[i].v_neigh.begin(), it);
        
        //if the link between two nodes is set to zero, we don't count it as a real link and we treat the nodes as disconnected
        if (E[v[i].v_link[index_j]].J != 0){
            
            if(!visited[*it_j])
                DFSUtil(*it_j, visited, k);
        }
            
    }
    
}
    

void Graph::ErdosRenyi(int p_M, double epsilon){
    
    std::mt19937 gen(seed_g);
    std::normal_distribution<double> d(0.,1);
    
    srand(seed_g);
    
    double Jij, Jji;
    
    
    M = p_M;
    
    vector<int> v;
    int i,j;
    int l = 0;
    int l_und =0;
    int flag;
    
    while (l_und < M){
        i = rand() % N ;
        j = rand() % N ;
        
        if (i!=j){
            
            double S = d(gen);
            double A = d(gen);
            
            double Jij = ( 1 - epsilon/2 ) * S + epsilon/2 * A;
            double Jji = ( 1 - epsilon/2 ) * S - epsilon/2 * A;
            
            v = make_vector<int>() << i << j;
            flag=addLink(l, v, Jij);
            
            //if we hadn't already added a the couple i-j, do it for both i-j and j-i, otherwise extract another pair
            
            if (flag){
                
                l ++;
                
                v = make_vector<int>() << j << i;
                flag=addLink(l, v, Jji);
                
                l ++;
                l_und ++;
                
            }
            
        }
    }
    
    
};



void Graph::ErdosRenyiConstrained(int p_M, double epsilon, int max_c){
    
    std::mt19937 gen(seed_g);
    std::normal_distribution<double> d(0.,1);
    
    srand(seed_g);
    
    double Jij, Jji;
    
    
    M = p_M;
    
    vector<int> vv;
    int i,j;
    int l = 0;
    int l_und =0;
    int flag;
    
    while (l_und < M){
        i = rand() % N ;
        j = rand() % N ;
        
        Node n_i = v[i];
        Node n_j = v[j];

        int c_i = n_i.numberOfLinks();
        int c_j = n_j.numberOfLinks();
        
        if (c_i < max_c && c_j < max_c){
            
            if (i!=j){
                
                double S = d(gen);
                double A = d(gen);
                
                double Jij = ( 1 - epsilon/2 ) * S + epsilon/2 * A;
                double Jji = ( 1 - epsilon/2 ) * S - epsilon/2 * A;
                
                vv = make_vector<int>() << i << j;
                flag=addLink(l, vv, Jij);
                
                //if we hadn't already added a the couple i-j, do it for both i-j and j-i, otherwise extract another pair
                
                if (flag){
                    
                    l ++;
                    
                    vv = make_vector<int>() << j << i;
                    flag=addLink(l, vv, Jji);
                    
                    l ++;
                    l_und ++;
                    
                }
                
            }
            
        }
        
    }
    
    
};


void Graph::RandomGraphNoSingle(int p_M, double epsilon, int max_c){
    
    
    std::mt19937 gen(seed_g);
    std::normal_distribution<double> d(0.,1);
    
    srand(seed_g);

    double Jij, Jji;
    
    
    M = p_M;
    
    vector<int> vv;
    int i,j,k;
    int l = 0;
    int l_und =0;
    int flag;
    
    i = 0;
    

    
    //we first random match the first N/2 nodes with the other half of the system.
    
    int N2 = (double)N/2;
    vector<int> v2(N2);
    
    for (int i=0; i<N2; i++)
        v2[i]=i+N2;
    
    int M_match = N2;
    
    while (l_und < M_match){

        k = rand() % v2.size();
        
        // we pair i with j
        j = v2[k];
        
        if (i!=j){
            
            double S = d(gen);
            double A = d(gen);
            
            double Jij = ( 1 - epsilon/2 ) * S + epsilon/2 * A;
            double Jji = ( 1 - epsilon/2 ) * S - epsilon/2 * A;
            
            vv = make_vector<int>() << i << j;
            flag=addLink(l, vv, Jij);
            
            //if we hadn't already added a the couple i-j, do it for both i-j and j-i, otherwise extract another pair
            
            if (flag){
                
                l ++;
                
                vv = make_vector<int>() << j << i;
                flag=addLink(l, vv, Jji);
                
                l ++;
                l_und ++;
                
            }
            
            v2.erase (v2.begin()+k);
            
        }
        
        i++;
        
    }
    
    //then, if the number of links is greater than the ones necessary to random pair nodes, we keep adding links till we reach the desidered connectivity
    
    
    if (M > N2) {
        
        while (l_und < M){
            i = rand() % N ;
            j = rand() % N ;
            
            Node n_i = v[i];
            Node n_j = v[j];
            
            int c_i = n_i.numberOfLinks();
            int c_j = n_j.numberOfLinks();
            
            if (c_i < max_c && c_j < max_c){
                
                if (i!=j){
                    
                    double S = d(gen);
                    double A = d(gen);
                    
                    double Jij = ( 1 - epsilon/2 ) * S + epsilon/2 * A;
                    double Jji = ( 1 - epsilon/2 ) * S - epsilon/2 * A;
                    
                    vv = make_vector<int>() << i << j;
                    flag=addLink(l, vv, Jij);
                    
                    //if we hadn't already added a the couple i-j, do it for both i-j and j-i, otherwise extract another pair
                    
                    if (flag){
                        
                        l ++;
                        
                        vv = make_vector<int>() << j << i;
                        flag=addLink(l, vv, Jji);
                        
                        l ++;
                        l_und ++;
                        
                    }
                    
                }
                
            }
            
        }
        
    }
    
    // else keep removing links at random from those that we added till now
    
    else{
        
        int undLinks = 0.5 * numberOfTotalLinks();
        int dirLinks = numberOfTotalLinks();
        
        while (undLinks > M){
            

            l = rand() % dirLinks;
            

            //decrease the number of links if we managed to remove one
            undLinks -= removeLink(l);
            
        }
        
    }
    

};


void Graph::RandomRegular(int p_M, double epsilon){
    
    std::mt19937 gen(seed_g);
    std::normal_distribution<double> d(0.,1);
    
    srand(seed_g);
    
    M = p_M;
    
    int c    = 2 * M / N;
    int cN   = c * N;
    int max, min, j;
    
    int node1, node2;
    int flagGraph = 0, flag;
    int l      = 0;
    int l_und  = 0;
    int count  = 0;
    
    vector <int> v;
    vector <int> v1;
    vector <int> v2;
    
    
    while (flagGraph == 0){
        
        v1.resize(cN);
        v2.resize(cN);
        
        
        for (int i=0; i < cN; ++i){
            v1[i] = i;
            v2[i] = i;
        }
        
        
        //vec_print(v1);
        //vec_print(v2);
        
        for (int i = 0; i < cN; i++){
            
            node1 = v1[0]/c;
            
            max = v2.size() - 1;
            min = 0;
            
            int m = (int)((double)v1[0] / c) * c;
            
            v.resize(c);
            for (int i = 0; i < c; i ++)
                v[i] = m + i;
            
            
            int flag = 0;
            while (!flag){
                j = rand() % (max - min + 1) + min;
                if (max != 1){
                    if(find(v.begin(), v.end(), v2[j]) != v.end())
                        flag = 0;
                    else
                        flag = 1;
                }
                else{
                    if(find(v.begin(), v.end(), v2[j]) != v.end()){
                        l_und = 0;
                        l = 0;
                        
                        clearGraph();
                        //cout << "******************* butto tutto" << endl;
                        break;
                    }
                    else
                        flag = 1;
                    
                    
                }
                
            }
            
            node2 = v2[j]/c;
            
            //cout << "extracted nodes " << endl;
            //cout << node1 << " " << node2 << " " << v1[0] << " " << v2[j] << endl;
            
            //random matching tra i nodi di v1 e di v2
            
            int i_tmp = v1[0];
            int j_tmp = v2[j];
            
            v1.erase(std::remove(v1.begin(), v1.end(), j_tmp), v1.end());
            v1.erase(std::remove(v1.begin(), v1.end(), i_tmp), v1.end());
            
            v2.erase(std::remove(v2.begin(), v2.end(), i_tmp), v2.end());
            v2.erase(std::remove(v2.begin(), v2.end(), j_tmp), v2.end());
            
            //cout << "stampo dopo aver eliminato: " << endl;
            
            //vec_print(v1);
            //vec_print(v2);
            
            if (node1 != node2){
                
                double S = d(gen);
                double A = d(gen);
                
                double J12 = ( 1 - epsilon/2 ) * S + epsilon/2 * A;
                double J21 = ( 1 - epsilon/2 ) * S - epsilon/2 * A;
                
                v = make_vector<int>() << node1 << node2;
                flag=addLink(l, v, J12);
                
                if (flag) l++;
                
                //check that we are not adding the same link twice
                //this is a check that has to be done only on the first ordering i-j becasue if we pass this test for i-j
                //we clearly pass it also for j-i
                
                else{
                    l = 0;
                    l_und = 0;
                    
                    clearGraph();
                    //cout << "******************* butto tutto" << endl;
                    break;
                }
                
                //we reach this stage only if we could add the first link i-j
                //anytime this is possible, we add also the link j-i.
                v = make_vector<int>() << node2 << node1;
                flag=addLink(l, v, J21);
                
                if (flag) l++;
                
                //cout << J12 << " " << J21 << endl;
                //cout << "l: " << l << endl;
                
                //we increase the counter of undirected links by 1
                l_und ++;

            }
            else{
                l = 0;
                l_und = 0;
                
                clearGraph();
                //cout << "******************* butto tutto" << endl;
                break;
                
            }
            
            if (l_und >= M){
                //cout << "BASTA LINK" << endl;
                flagGraph = 1;
                break;
            }
            
        }
        
        count++;
        
    }
    
}



void Graph::clearGraph(){
    
    vector<Node>().swap(v);
    vector<Link>().swap(E);
    
    v.reserve(N);
    for (int i = 0; i < N; ++i)
        v.push_back (Node(i));
    
};



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
            bias[i][index_j] = tmp_mess;
            
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
        
        
        if (it_l->l == 62) cout << "message from node " << i << " to node " << j << endl;
        //compute message mess_i_to_j
        vector <vector <long double> > in_mess;
        
        //couplings from nodes k's to node i.
        vector <double> vec_J_k_to_i;
        
        z = 0; //z is the connectivity of node i less 1.
        
        //for each link k->i
        for(vector<int>::iterator it_k = G.v[i].v_neigh.begin(); it_k != G.v[i].v_neigh.end(); ++it_k){
            //where k is different from j
            if(*it_k != j){
                if (it_l->l == 62) cout << "k: " << *it_k << endl;
                
                vector<int>::iterator it = find(G.v[*it_k].v_neigh.begin(), G.v[*it_k].v_neigh.end(), i);
                int index_i = distance (G.v[*it_k].v_neigh.begin(), it);
                
                //we store the mess_k_to_i
                in_mess.push_back(mess[*it_k][index_i]);
                
                
                it = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), *it_k);
                int index_k = distance (G.v[i].v_neigh.begin(), it);
                
                int index_J_k_to_i = G.v[i].v_link[index_k];
                
                if (it_l->l == 62) cout << "Jk->i " << G.E[index_J_k_to_i].J << endl;
                
                vec_J_k_to_i.push_back(G.E[index_J_k_to_i].J);
                z++;
            }
            //where k is equal to j
            else{
                if (it_l->l == 62) cout << "j: " << *it_k << endl;
                
                vector<int>::iterator it_j = find(G.v[i].v_neigh.begin(), G.v[i].v_neigh.end(), *it_k);
                int index_j = distance (G.v[i].v_neigh.begin(), it_j);
                
                int index_J_j_to_i = G.v[i].v_link[index_j];
                
                J_j_to_i = G.E[index_J_j_to_i].J;
                
                if (it_l->l == 62) cout << "Jj->i: " << G.E[index_J_j_to_i].J << endl;
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
                    
                    if ( gip * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xi - 1 ) > 0 ){
                        
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
                        
                        if ( gip * ( 2 * xip - 1 ) > 0 && gip * ( 2 * xi - 1 ) > 0 ){
                            
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
                
            
            if ( gnp * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xn - 1 ) > 0 ){
                
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
                
                if ( gnp * ( 2 * xnp - 1 ) > 0 && gnp * ( 2 * xn - 1 ) > 0 ){
                    
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
    
    printEntropy<Tag>(Part1, Part2);

}


void vec_print(vector<int>& vec){
    for (int i=0; i<vec.size(); i++)
        cout << vec[i] << ' ';
    cout << endl;
}

//------------------------------------------------------------------------------------- PRINT FUNCTIONS ---------------------------------------------------------------------------------//

void Messages::def_printEntropy(L1, double Part1, double Part2){
    
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
    
}


void Messages::def_printEntropy(L2, double Part1, double Part2){
    
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
    
}


void Messages::def_printEntropy(L3, double Part1, double Part2){
    
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
}


void Messages::def_printEntropy(L4, double Part1, double Part2){
    
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
    
}


void Messages::def_printEntropy(L1T1, double Part1, double Part2){
    
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
    
}


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



//----------------------------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class Messages
//----------------------------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------------------------//

template<typename Tag>
void BPGD::initDecimation(vector<int>& v_bias, vector<bool>& v_q){
    
    
    for(int i = 0; i < mess.N; ++i)
        notFixedSpins.push_back(i);
    

    
    if(v_bias.size()){

        int size = v_bias.size();
        if (size != v_q.size()){
            cout << "error: v_q and v_bias have to have the same size!" << endl;
            return;
        }
        else{
            for (int i = 0; i < size; ++i){
                int n = v_bias[i];
                fixedSpins.push_back(n);
                fixedValues.push_back(v_q[i]);
                notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), n), notFixedSpins.end());
                
            }
        }
        
        cout << "SIZE NOT FIXED: " << notFixedSpins.size() ;
        
        setHardBias<Tag>(v_bias, v_q);
        

        
    }

};

template<typename Tag>
void BPGD::setHardBias(vector<int>& v_bias, vector<bool>& v_q){
    
    //v_bias contains the indices of nodes that have to be biased
    //v_q contains the color towards which they have to be biased
    
    int size = v_bias.size();
    if (size != v_q.size()){
        cout << "error: v_q and v_bias have to have the same size!" << endl;
        return;
    }
    else{
        for (int i = 0; i < size; ++i){
            int n = v_bias[i];
            //we set the bias towards the value specified in v_q[i]
                        
            setHardBiasSite<Tag>(n, v_q[i]);
            
        }
        
    }
    
};


inline void BPGD::def_setHardBiasSite(L1, int n, int value) {
    
    for (vector<int>::iterator it = mess.G.v[n].v_neigh.begin() ; it != mess.G.v[n].v_neigh.end(); ++it){
        int j = *it;
                           
        vector<int>::iterator it_j = find(mess.G.v[n].v_neigh.begin(), mess.G.v[n].v_neigh.end(), j);
        int index_j = distance (mess.G.v[n].v_neigh.begin(), it_j);
    

        if (value == 0){
            int xn = 1;
            for(int xj = 0; xj <1; xj++){
                int t_nj = 2 * xn + xj;
                mess.bias[n][index_j][t_nj] = 0.;
            }
        }
        else{
            int xn = 0;
            for(int xj = 0; xj <1; xj++){
                int t_nj = 2 * xn + xj;
                mess.bias[n][index_j][t_nj] = 0.;
            }
                   
        }
        
        for(int xn = 0; xn <=1; xn++){
            for(int xj = 0; xj <=1; xj++){
                int t_nj = 2 * xn + xj;
                cout << "stampo bias " <<mess.bias[n][index_j][t_nj] << endl;
            }
        }
        
    }
        
}


void BPGD::findMostBiased(vector<int> & v_bias, vector<bool>& v_q){
    
    int    i_max = notFixedSpins[0];
    bool   col;
    double tmp, max = 0.;
    
    
    for (vector<int>::iterator it_i = notFixedSpins.begin(); it_i != notFixedSpins.end(); ++it_i){
        tmp = mess.node_marginal[*it_i][0] - mess.node_marginal[*it_i][1];
        if (abs(tmp) > 0.999){
            
            if (tmp > 0)
                col = 0;
            else
                col = 1;
            
            v_bias.push_back(*it_i);
            v_q.push_back(col);
            
        }
        
    }
    
    if (v_bias.size() == 0){
        
        for (vector<int>::iterator it_i = notFixedSpins.begin(); it_i != notFixedSpins.end(); ++it_i){
            tmp = mess.node_marginal[*it_i][0] - mess.node_marginal[*it_i][1];
            if (abs(tmp) > max){
                
                if (tmp > 0)
                    col = 0;
                else
                    col = 1;
                
                max = abs(tmp);
                
                i_max = *it_i;
                
            }
            
        }
        
        v_bias.push_back(i_max);
        v_q.push_back(col);

    }
    
    
};
 


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

template void Messages::BPiteration<L1>(double, int, int, int, bool);
template void Messages::BPiteration<L2>(double, int, int, int, bool);
template void Messages::BPiteration<L3>(double, int, int, int, bool);
template void Messages::BPiteration<L4>(double, int, int, int, bool);

template void Messages::BPiteration<L1T1>(double, int, int, int, bool);

template void Messages::Wrap_computeLastMarg<L1>();
template void Messages::Wrap_computeLastMarg<L2>();


template void Messages::logPartitionFunction<L1>();
template void Messages::logPartitionFunction<L2>();
template void Messages::logPartitionFunction<L3>();
template void Messages::logPartitionFunction<L4>();

template void Messages::logPartitionFunction<L1T1>();



template void BPGD::setHardBias<L1>(vector<int>&, vector<bool>&);

template void BPGD::initDecimation<L1>(vector<int>&, vector<bool>&);
