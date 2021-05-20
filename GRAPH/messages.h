#include "graph.h"


using namespace std;


//name the type of constraints considered
struct L1 { };     //cycle of length 1
struct L2 { };     //cycle of lenght 2
struct L3 { };     //cycle of lenght 3
struct L4 { };     //cycle of lenght 4

struct L1T1 { };  //1 time basin of attraction of cycle of lenght 1
struct L1T2 { };  //2 time basin of attraction of cycle of lenght 1
struct L1T3 { };  //3 time basin of attraction of cycle of lenght 1
struct L2T1 { };  //1 time basin of attraction of cycle of lenght 2
struct L2T2 { };  //2 time basin of attraction of cycle of lenght 2
struct L2T3 { };  //3 time basin of attraction of cycle of lenght 2

class Messages{
public:

    Messages();
    
    Messages(Graph &, int, double);
    
    int N;                                                      //number of nodes
    int M;                                                      //number of links
    int Q;                                                      //number of components of the supermarginals 4^L
    int QS;                                                     //number of component of single site marginals 2^L
    int L;                                                      //length of the cycle we consider
    
    int r;                                                      //strenght of the random term in the initial condition of the BP iteration procedure
    
    Graph G;
    
    //use tag dispatching whenever you want to define functions with the same name doing different things according to different labels
    template<typename Tag>
    inline void constraint(int z, int flag_red, int n, vector <double>& vec_J_k_to_n){
        def_constraint(Tag(), flag_red, z, n, vec_J_k_to_n);
    }

    //these are different functions, labelled by the first argument
    inline void def_constraint(L1, int flag_red, int z, int n, vector <double>& vec_J_k_to_n);
    inline void def_constraint(L2, int flag_red, int z, int n, vector <double>& vec_J_k_to_n);
    inline void def_constraint(L3, int flag_red, int z, int n, vector <double>& vec_J_k_to_n);
    inline void def_constraint(L4, int flag_red, int z, int n, vector <double>& vec_J_k_to_n);
    
    inline void def_constraint(L1T1, int flag_red, int z, int n, vector <double>& vec_J_k_to_n);
    inline void def_constraint(L1T2, int flag_red, int z, int n, vector <double>& vec_J_k_to_n);
    inline void def_constraint(L1T3, int flag_red, int z, int n, vector <double>& vec_J_k_to_n);
    inline void def_constraint(L2T1, int flag_red, int z, int n, vector <double>& vec_J_k_to_n);
    inline void def_constraint(L2T2, int flag_red, int z, int n, vector <double>& vec_J_k_to_n);
    inline void def_constraint(L2T3, int flag_red, int z, int n, vector <double>& vec_J_k_to_n);

    
    template<typename Tag>
    inline void constraint_bp(int z, int flag_red, int l, vector <double>& vec_J_k_to_i, double J_j_to_i){
        def_constraint_bp( Tag(), flag_red, z, l, vec_J_k_to_i, J_j_to_i );
    }
    //these are different functions, labelled by the first argument
    inline void def_constraint_bp(L1, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i);
    inline void def_constraint_bp(L2, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i);
    inline void def_constraint_bp(L3, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i);
    inline void def_constraint_bp(L4, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i);
    
    inline void def_constraint_bp(L1T1, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i);
    inline void def_constraint_bp(L1T2, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i);
    inline void def_constraint_bp(L1T3, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i);
    inline void def_constraint_bp(L2T1, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i);
    inline void def_constraint_bp(L2T2, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i);
    inline void def_constraint_bp(L2T3, int flag_red, int z, int l, vector <double>& vec_J_k_to_i, double J_j_to_i);
    
    
    template<typename Tag>
    inline vector<long double> computeMarginalsLastTime(int l){
        return def_computeMarginalsLastTime(Tag(), l);          //For a supervariable describing the state of several time steps, take the latest time and marginalise over previous times
    }
    
    inline vector<long double> def_computeMarginalsLastTime(L1, int l);
    inline vector<long double> def_computeMarginalsLastTime(L1T1, int l);
    inline vector<long double> def_computeMarginalsLastTime(L2, int l);

    
    
    template<typename Tag>
    void printEntropy(double Part1, double Part2, int print_file){
        def_printEntropy(Tag(), Part1, Part2, print_file);
    }
    
    void def_printEntropy(L1, double Part1, double Part2, int print_file);
    void def_printEntropy(L2, double Part1, double Part2, int print_file);
    void def_printEntropy(L3, double Part1, double Part2, int print_file);
    void def_printEntropy(L4, double Part1, double Part2, int print_file);

    void def_printEntropy(L1T1, double Part1, double Part2, int print_file);
    void def_printEntropy(L1T2, double Part1, double Part2, int print_file);
    void def_printEntropy(L1T3, double Part1, double Part2, int print_file);
    void def_printEntropy(L2T1, double Part1, double Part2, int print_file);
    void def_printEntropy(L2T2, double Part1, double Part2, int print_file);
    void def_printEntropy(L2T3, double Part1, double Part2, int print_file);

    
    
    template<typename Tag>
    void print_BPit(double tmp_th, int t){
        def_print_BPit(Tag(), tmp_th, t);
    }
    
    void def_print_BPit(L1, double tmp_th, int t);
    void def_print_BPit(L2, double tmp_th, int t);
    void def_print_BPit(L3, double tmp_th, int t);
    void def_print_BPit(L4, double tmp_th, int t);

    void def_print_BPit(L1T1, double tmp_th, int t);
    void def_print_BPit(L1T2, double tmp_th, int t);
    void def_print_BPit(L1T3, double tmp_th, int t);
    void def_print_BPit(L2T1, double tmp_th, int t);
    void def_print_BPit(L2T2, double tmp_th, int t);
    void def_print_BPit(L2T3, double tmp_th, int t);

    
    template<typename Tag>
    void look_up_table(int flag_red);                           //this function fills the vector allowed_conf for each node.
    
    template<typename Tag>
    void look_up_table_bp(int flag_red);                        //this function fills the vector allowed_conf_bp for each directed link. This will be used in each step of the BP eqs.
    
    void compute_x(int Nc, int z, vector < vector <int> >& x);  //compute the matrix of possible configurations of neighbours of a site

    vector <int> sw;                                            //auxiliary vector used to compute super marginals, that basically containes the switched i-j terms
    
    vector <int> t_reduced;                                     //this is an auxiliary vector used to map messages configurations to their simmetric ones in the reduced case
                                                                //it is defined in reducedAlgo
    
    double factor_S_link, factor_S_node_1, factor_S_node_2;     //factors set to 1 in the general case, and useful in the reduced case.
    
    vector < vector < vector <long double> > > mess;            //mess contains the messages from nodes to nodes
                                                                //Loosely speaking
                                                                //mess[i][j][t] is the t-th component of the message from node i to node j
                                                                //if Q = 4 (i.e. we study supermarginal relative to two variables):
                                                                //t = 2 * x_i + x_j, i.e. t = 0,1,2,3 while (xi xj) = [(0,0);(0, 1);(1,0);(1,1)].
                                                                //if Q = 16 (i.e. we study supermarginal relative to four variables):
                                                                //t = 8 * x_i + 4 * x_j + 2 * x_ip + x_jp, i.e. t = 0,1,2,3,...,15 while (xi xj xip xjp) = [(0,0,0,0);(0,0,0,1);(0,0,1,0);(0,0,1,1);...;(1,1,1,1)].

    vector < vector < vector <long double> > > bias;
    
    vector < vector < vector <int> > > allowed_conf_bp;         //for each i->j and (xi,xj), this vector contains the configurations of neighbours of i, indexed by
                                                                //k (diff than j), that satisfy the constraint. This configuration in encoded in a decimal form
                                                                //taking the 0/1 representation of the allowed configuration and getting its decimal form

    vector < vector < vector <int> > > allowed_conf;            //for each i and xi, this vector contains the configurations of neighbours of i, indexed by
                                                                //k  that satisfy the constraint. This configuration in encoded in a decimal form
                                                                //taking the 0/1 representation of the allowed configuration and getting its decimal form

    vector < vector < vector <long double> > > update_mess;     //update_mess contains the messages after a BP-sweep

    vector < vector < vector <long double> > > link_marginal;   //link_marginal probability
                                                                //given the particular form of the problem, where we have a super-variable for each link of the original graph
                                                                //we have the same variables in multiple super-variables and
                                                                //at link l, marginal[l] contains two marginals, one for each of the spins that are connected from l.

    vector < vector <long double> > node_marginal;              //node_marginal probability
                                                                //here we store in a less redundant way the information contained in link_marginal
                                                                //considering for each site i the marginal on that site.
                                                                //IMPORTANT: to improve performaces, it is RECOMMENDED to call this method only after the convergence of BP equations
                                                                //otherwise, if we want to call it after every BP sweep, for instance to study nodeMarginalState(),
                                                                //it is necessary to free the vector with marginal_node.swap( vector < vector<double> >() );
    
    vector < vector <long double> > super_marginal;             //super_marginal contains the marginal of the super variables
    vector < vector <long double> > prev_super_marginal;        //prev_super_marginal is used to store the values of the super_marginals at the previous BP time step
    
    
    
    void initMessages();                                        //this method initializes messages to umbiased values 1/Q.
                                                                //the marginals are set to the uniform distribution as well.
                                                                //it is called inside the constructor
                                                                 
    
    vector<int> Switch();                                       //method that switch i-j components in the vector of possible configurations x
    
    void setQ(int);                                             //this method implements all the factors and vectors that we need in the case we use the reduced algo

    void messUpdate();                                          //messUpdate implements a BP-sweep on the whole graph
                                                                //updated messages are stored in update_mess
    
    void messUpdate_approx(int);                                //implements messUpdate using a montecarlo on links for which we need to update over a number of variables
                                                                //we need to update over a number of variables larger than MM
                                                                //input variable:
                                                                //MM:                   number of configuration we sum over
    
    bool messNormalize();                                       //this method normalizes the messagess. it retuns 0 if the probability is not normalizable (all entries equal to 0)
    
    void superMarginals();                                      //this method computes (do not print) super_marginals
    
    void nodeMarginals();                                       //this method computes (do not print) node_marginals
    
    void messState();                                           //this method prints the state of the messages
    void linkMarginalState();                                   //this method prints the state of the link_marginals
    void nodeMarginalState();                                   //this method prints the state of the node_marginals
    void superMarginalsState();                                 //this method prints the state of the super marginals
    
    double compareSuperMarginals();                             //this method compare super_marginals at time step t-1 and t during BPiteration

    template<typename Tag>                                      //for each link l connected to nodes i and j, compute the probability that
    void Wrap_computeLastMarg();                                //spin i is 0 and spin j is 0
    
    template<typename Tag>
    void BPiteration(double, int, int, int, bool);              //this method iterates BP equations by calling BP_sweep until convergence
                                                                //input variables:
                                                                //th          : this value sets the convergence quality. set it to ~10^-3.
                                                                //flag_red    : flag equal to 1/0 to study reduced/non-reduced case
                                                                //flag_approx : flag to set to 1/0 to use approximate or exact BP updating
                                                                //T           : maximum iteration time. set it to ~ N.
                                                                //verbose     : set to one to have a verbose version
    
    void updateAndStore();                                      //this method sets the update_mess in mess for the next BP-sweep
                                                                //and store the current super_marginals in prev_super_marginals
    
    template<typename Tag>
    void logPartitionFunction();
    
    friend class BPGD;

};
