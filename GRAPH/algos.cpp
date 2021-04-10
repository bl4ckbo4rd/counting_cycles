#include "algos.h"


//----------------------------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------------------------// methods of class BPGD
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

template void BPGD::setHardBias<L1>(vector<int>&, vector<bool>&);

template void BPGD::initDecimation<L1>(vector<int>&, vector<bool>&);
