#include "messages.h"


using namespace std;


class BPGD {
public:
    

    Messages mess;
    
    BPGD(Messages & p_mess) : mess(p_mess) {};
    
    template<typename Tag>
    void initDecimation(vector<int> &, vector<bool> &, vector<int> &, vector<bool> &, vector<int> &);         //this method set an hard bias on some specified variables
                                                                //the first vector contains the biased variables
                                                                //the second vector contains the colours toward which the biased nodes are biased.
                                                                //it also fills the vectors fixedSpins and fixedValues and notFixedSpins.
                                                                //if no specified variable are given, the last vector is made by all the variables of the system
    
    vector <int> fixedSpins;                                    //this vector contains the spin that get frozen (decimated) along the decimation process.
                                                                //it is filled by the method initDecimation and by the method fixSpins during the decimation process.
    
    vector <bool> fixedValues;                                  //this vector contains the values of the spin that get frozen (decimated) along the decimation process.
                                                                //it is filled by the method initDecimation and by the method fixSpins during the decimation process.
    
    vector <int> notFixedSpins;                                 //this vector contains the indices of the spin that are not frozen.
                                                                //at t=0, it is formed by all the spins, unless we decide to run the decimation process
                                                                //from a specific set of fixed nodes by feeding initDecimation with it.
                                                                //As t increases and the decimation process continues, its size decreases.
                                                                //it is filled by the method initDecimation and by the method fixSpins during the decimation process.
    
    template<typename Tag>
    void setHardBias(vector<int>&, vector<bool>&);              //set an hard bias on the nodes specified by the first vector according to the value specified in the second.
                                                                //default input vectors are empty vectors.
                                                                //input variables: vector of nodes to be biased, vector of colors to towards which node biases have to be biased
    
    template<typename Tag>
    inline void setHardBiasSite(int i, int value){
        def_setHardBiasSite( Tag(), i, value );
    }
    
    
    //these are different functions, labelled by the first argument
    inline void def_setHardBiasSite(L1, int i, int value);
    
    template<typename Tag>
    void findMostBiased(vector<int>&, vector<bool>&, vector<int>&, vector<bool>&, vector<int>&);           //after having ran the BP equation till convergence, we find the most biased variables.
                                                                //by this we mean variables for which | p[0]-p[1] | > 0.999 or,if none, the variable with the largest
                                                                //absolute value of the difference p[0]-p[1].
                                                                //the search is carefully made only on the variables that have not been fixed yet.
                                                                //this method is called inside BPguidedDecimation.
    template<typename Tag>
    void BPGDiteration(double, int, int, int, bool, vector<int>&, vector<bool>&, vector<int>&, vector<bool>&, vector<int>&);     

    
};

