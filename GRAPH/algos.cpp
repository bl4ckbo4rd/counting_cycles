#include "algos.h"


//----------------------------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class BPGD
//----------------------------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------------------------//

template<typename Tag>
void BPGD::initDecimation(vector<int>& v_bias, vector<bool>& v_q, vector<int>& fixedSpins, vector<bool>& fixedValues, vector<int>& notFixedSpins){
    
    
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

        vector<int>::iterator it_i = find(mess.G.v[j].v_neigh.begin(), mess.G.v[j].v_neigh.end(), n);
        int index_i = distance (mess.G.v[j].v_neigh.begin(), it_i);
    

        if (value == 0){
            int xn = 1;
            for(int xj = 0; xj <=1; xj++){
                int t_nj = 2 * xn + xj;
                mess.bias[n][index_j][t_nj] = 0.;

                int t_jn = 2 * xj + xn;
                mess.bias[j][index_i][t_jn] = 0.;
            }

        }
        else{
            int xn = 0;
            for(int xj = 0; xj <=1; xj++){
                int t_nj = 2 * xn + xj;
                mess.bias[n][index_j][t_nj] = 0.;

                int t_jn = 2 * xj + xn;
                mess.bias[j][index_i][t_jn] = 0.;
            }
                   
        }
        
        for(int xn = 0; xn <=1; xn++){
            for(int xj = 0; xj <=1; xj++){
                int t_nj = 2 * xn + xj;
                cout << "fixed spin: " << n << "and neighbour: " << j << endl;
                cout << "print bias: " <<mess.bias[n][index_j][t_nj] << endl;
            }
        }
        
    }
        
}

inline void BPGD::def_setHardBiasSite(L1T1, int n, int value) {
    
    for (vector<int>::iterator it = mess.G.v[n].v_neigh.begin() ; it != mess.G.v[n].v_neigh.end(); ++it){
        int j = *it;
                           
        vector<int>::iterator it_j = find(mess.G.v[n].v_neigh.begin(), mess.G.v[n].v_neigh.end(), j);
        int index_j = distance (mess.G.v[n].v_neigh.begin(), it_j);

        vector<int>::iterator it_i = find(mess.G.v[j].v_neigh.begin(), mess.G.v[j].v_neigh.end(), n);
        int index_i = distance (mess.G.v[j].v_neigh.begin(), it_i);
    

        if (value == 0){
            int xn = 1;
            for(int xj = 0; xj <=1; xj++){
                for(int xnp = 0; xnp <=1; xnp++){
                    for(int xjp = 0; xjp <=1; xjp++){
                        int t_nj = 8 * xnp + 4 * xjp + 2 * xn + xj;
                        mess.bias[n][index_j][t_nj] = 0.;

                        int t_jn = 8 * xnp + 4 * xjp + 2 * xj + xn;
                        mess.bias[j][index_i][t_jn] = 0.;
                    }
                }
            }

        }
        else{
            int xn = 0;
            for(int xj = 0; xj <=1; xj++){
                for(int xnp = 0; xnp <=1; xnp++){
                    for(int xjp = 0; xjp <=1; xjp++){
                        int t_nj = 8 * xnp + 4 * xjp + 2 * xn + xj;
                        mess.bias[n][index_j][t_nj] = 0.;

                        int t_jn = 8 * xnp + 4 * xjp + 2 * xj + xn;
                        mess.bias[j][index_i][t_jn] = 0.;
                    }
                }
            }
                   
        }
        
        for(int xn = 0; xn <=1; xn++){
            for(int xj = 0; xj <=1; xj++){
                for(int xnp = 0; xnp <=1; xnp++){
                    for(int xjp = 0; xjp <= 1; xjp++){
                        int t_nj = 8 * xnp + 4 * xjp + 2 * xn + xj;
                        cout << "fixed spin: " << n << "and neighbour: " << j << endl;
                        cout << "print bias: " <<mess.bias[n][index_j][t_nj] << endl;
                    }
                }
            }
        }
        
    }
        
}

inline void BPGD::def_setHardBiasSite(L1T2, int n, int value) {
    
    for (vector<int>::iterator it = mess.G.v[n].v_neigh.begin() ; it != mess.G.v[n].v_neigh.end(); ++it){
        int j = *it;
                           
        vector<int>::iterator it_j = find(mess.G.v[n].v_neigh.begin(), mess.G.v[n].v_neigh.end(), j);
        int index_j = distance (mess.G.v[n].v_neigh.begin(), it_j);

        vector<int>::iterator it_i = find(mess.G.v[j].v_neigh.begin(), mess.G.v[j].v_neigh.end(), n);
        int index_i = distance (mess.G.v[j].v_neigh.begin(), it_i);
    

        if (value == 0){
            int xn = 1;
            for(int xj = 0; xj <=1; xj++){
                for(int xnp = 0; xnp <=1; xnp++){
                    for(int xjp = 0; xjp <=1; xjp++){
                        for(int xnpp = 0; xnpp <=1; xnpp++){
                            for(int xjpp = 0; xjpp <=1; xjpp++){
                                int t_nj = 32 * xnpp + 16 * xjpp + 8 * xnp + 4 * xjp + 2 * xn + xj;
                                mess.bias[n][index_j][t_nj] = 0.;

                                int t_jn = 32 * xnpp + 16 * xjpp + 8 * xnp + 4 * xjp + 2 * xj + xn;
                                mess.bias[j][index_i][t_jn] = 0.;
                            }
                        }
                    }
                }
            }

        }
        else{
            int xn = 0;
            for(int xj = 0; xj <=1; xj++){
                for(int xnp = 0; xnp <=1; xnp++){
                    for(int xjp = 0; xjp <=1; xjp++){
                        for(int xnpp = 0; xnpp <=1; xnpp++){
                            for(int xjpp = 0; xjpp <=1; xjpp++){
                                int t_nj = 32 * xnpp + 16 * xjpp + 8 * xnp + 4 * xjp + 2 * xn + xj;
                                mess.bias[n][index_j][t_nj] = 0.;

                                int t_jn = 32 * xnpp + 16 * xjpp + 8 * xnp + 4 * xjp + 2 * xj + xn;
                                mess.bias[j][index_i][t_jn] = 0.;
                            }
                        }
                    }
                }
            }
                   
        }
        
        for(int xn = 0; xn <=1; xn++){
            for(int xj = 0; xj <=1; xj++){
                for(int xnp = 0; xnp <=1; xnp++){
                    for(int xjp = 0; xjp <= 1; xjp++){
                        for(int xnpp = 0; xnpp <=1; xnpp++){
                            for(int xjpp = 0; xjpp <= 1; xjpp++){
                                int t_nj = 32 * xnpp + 16 * xjpp + 8 * xnp + 4 * xjp + 2 * xn + xj;
                                cout << "fixed spin: " << n << "and neighbour: " << j << endl;
                                cout << "print bias: " <<mess.bias[n][index_j][t_nj] << endl;
                            }
                        }
                    }
                }
            }
        }
        
    }
        
}

template<typename Tag>
void BPGD::findMostBiased(vector<int> & v_bias, vector<bool>& v_q, vector<int> & fixedSpins, vector<bool>& fixedValues, vector<int> & notFixedSpins){
    
    int    i_max = notFixedSpins[0];
    bool   col;
    double tmp, max = 0.;
    
    
    //for (vector<int>::iterator it_i = notFixedSpins.begin(); it_i != notFixedSpins.end(); ++it_i){
    //    tmp = mess.node_marginal[*it_i][0] - mess.node_marginal[*it_i][1];
    //    if (abs(tmp) > 0.999){
            
    //        if (tmp > 0)
    //            col = 0;
    //        else
    //            col = 1;
            
    //        v_bias.push_back(*it_i);
    //        v_q.push_back(col);
            
    //    }
        
    //}
    
    //if (v_bias.size() == 0){
        
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

    //}
    

    int size = v_bias.size();
    if (size != v_q.size()){
        cout << "error: v_q and v_bias have to have the same size!" << endl;
        return;
    }else{
        int n = v_bias[size-1];
        fixedSpins.push_back(n);
        fixedValues.push_back(v_q[size-1]);
        notFixedSpins.erase(remove(notFixedSpins.begin(), notFixedSpins.end(), n), notFixedSpins.end());
        cout << n << " removed!" << endl;        
            
    }
        
    cout << "SIZE NOT FIXED: " << notFixedSpins.size() << endl;
    cout << "v_bias size: " << v_bias.size() << endl ;
    cout << "v_q size: " << v_q.size() << endl;
        
    setHardBias<Tag>(v_bias, v_q);
        

};

template<typename Tag>
void BPGD::BPGDiteration(double th, int flag_red, int flag_approx, int T, bool verbose, vector<int> & v_bias, vector<bool>& v_q, vector<int> & fixedSpins, vector<bool>& fixedValues, vector<int> & notFixedSpins){
    
    mess.look_up_table_bp<Tag>(flag_red);
    mess.look_up_table<Tag>(flag_red);

    for(int l = 0; l < 0.5*mess.G.numberOfTotalNodes()+1; l++){

        cout << endl;
        cout << "iteration: " << l << endl;
        cout << endl;

        int    t = 0;
        double tmp_th = 1;

        //MM is the number of configurations used to compute the BP update in the approximate case
        int MM = 10000;
    
        while (tmp_th > th && t < T){
            //the problem with this approximate algo is that the random number generation is particularly slow
            if(flag_approx)
                mess.messUpdate_approx(MM);
            else
                mess.messUpdate();


            mess.messNormalize();
            mess.superMarginals();
        
            tmp_th = mess.compareSuperMarginals();
        
            mess.updateAndStore();
        
            if(verbose){
                mess.messState();
                cout << endl;
                cout << "BP iteration: at time t=" << t << " the maximum error between current and previous super_marginals is " << tmp_th << endl;
            }
        
            t++;
        
        }
        mess.Wrap_computeLastMarg<Tag>();
        //mess.linkMarginalState();
        mess.nodeMarginals();
        mess.nodeMarginalState();
        findMostBiased<Tag>(v_bias, v_q, fixedSpins, fixedValues, notFixedSpins);
        mess.superMarginalsState();


        if(l == 0.5*mess.G.numberOfTotalNodes()){

            mess.print_BPit<Tag>(tmp_th, t);
    
            if(verbose){
                cout << "the final node_marginals are " << endl;
                cout << endl;
                mess.superMarginalsState();
                cout << endl;
            }
    

            mess.logPartitionFunction<Tag>();

        }
    }
     
    
};
 

//these are template methods and needs to be specified for the linker

template void BPGD::setHardBias<L1>(vector<int>&, vector<bool>&);
template void BPGD::setHardBias<L1T1>(vector<int>&, vector<bool>&);
template void BPGD::setHardBias<L1T2>(vector<int>&, vector<bool>&);

template void BPGD::initDecimation<L1>(vector<int>&, vector<bool>&, vector<int>&, vector<bool>&, vector<int>&);
template void BPGD::initDecimation<L1T1>(vector<int>&, vector<bool>&, vector<int>&, vector<bool>&, vector<int>&);
template void BPGD::initDecimation<L1T2>(vector<int>&, vector<bool>&, vector<int>&, vector<bool>&, vector<int>&);

template void BPGD::BPGDiteration<L1>(double, int, int, int, bool, vector<int>&, vector<bool>&, vector<int>&, vector<bool>&, vector<int>&);
template void BPGD::BPGDiteration<L1T1>(double, int, int, int, bool, vector<int>&, vector<bool>&, vector<int>&, vector<bool>&, vector<int>&);
template void BPGD::BPGDiteration<L1T2>(double, int, int, int, bool, vector<int>&, vector<bool>&, vector<int>&, vector<bool>&, vector<int>&);
