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

Graph::Graph(){};


Graph::Graph(int p_N, int p_seed_g) : N(p_N), seed_g(p_seed_g) {
    
    v.reserve(N);
    for (int i = 0; i < N; ++i)
        v.push_back (Node(i));
    
};


Graph::initializeGraph(nt p_N, int p_seed_g): N(p_N), seed_g(p_seed_g) {
    
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

