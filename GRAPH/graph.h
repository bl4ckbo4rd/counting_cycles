#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdlib.h>
#include <cassert>
#include <iterator>
#include <random>
#include <time.h>


using namespace std;



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//------------------------------------------------------------------------------------------------------------------------------------------------------------// class representing a node
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



class Node{
public:

    Node(int p_n) : n (p_n) {};
  
    int n;                                                //n is the index of the Node
    vector <int> v_link;                                  //v_link contains the indices of links attached to the node
    vector <int> v_neigh;                                 //v_neigh contains the indices of the nodes attached to the node
    
    int numberOfLinks();                                  //this method returns the number of links to which the Node is attached
  
};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//----------------------------------------------------------------------------------------------------------------------------------------------------------// class representing a factor
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



class Link{
public:
    
    Link(int, double);

    int l;                                                      //l is the index of the link
    vector <int> v_node;                                        //v_node contains the indices of nodes connected by l

    double J;                                                      //value of the coupling on the link
    
};



//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//
//---------------------------------------------------------------------------------------------------------------------------------------------------------// class representing the graph
//--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------//



class Graph{
public:
    
    Graph();
    Graph(int, int);
    
    
    int seed_g;
    
    int N;
    int M;
    
    vector <Node> v;                                        //v containss all the Nodes of the Graph.
                                                            //this vector is filled by the method getGraph.
    
    vector <Link> E;                                        //E contains all the directed Links of the Graph.
    
    vector< vector<int> > component;

    void initializeGraph(int, int);                         //do the job of the constructur in case the parameters N and seed_g were not specified when
                                                            //instantiating the graph object.
    
    int removeLink(int l);                                  //TODO: check if it works
    
    int addLink(int, vector<int>, double);                  //the addLink method implements the operations needed when adding a link in the graph
                                                            //namely one needs to create a Link object l, store its neighour variables and for each of them
                                                            //add the index of l as a neghbouring link.
                                                            //input variables: link index, vector of nodes connected by the link, value of the coupling:
                                                            //ex.:
                                                            //vector <int> vitoj = make_vector<int>() << j << i;
                                                            //flag_ij = addLink(l, vitoj, Jitoj);
                                                            //output: 1 if the operation can be done, 0 if not (the link already exists).
    
    int addLinkWrapper(int, int, double, double);           //wrapper of the addLink method taking as inputs the labels of the nodes one wants to connect
                                                            //and the values of the couplings that the two directed links should take. if the links (i->j and j->i)
                                                            //exist, overwrite the values of Ji->j and Jj->i.
                                                            //inputs: i, j, Ji->j, Jj->i
    
    int numberOfTotalLinks();                               //this method returns the total number of directed Links in the Graph.
    
    int numberOfTotalNodes();                               //this method returns the total number or nodes in the graph.

    void neighsOfNode(int);                                 //this method returns the links attached to the input node
    void nodesOfLink(int);                                  //this method returns the nodes connected by the link
    
    void ErdosRenyi(int, double);                           //Erdos Renyi random graph
    void RandomRegular(int, double);                        //Random regular graph
    void ErdosRenyiConstrained(int, double, int);           //Erdos Renyi random graph with the constraint that the maximum commectity is given by the third argument
    void RandomGraphNoSingle(int, double, int);             //Random Graph with no single nodes and a constraint on the maximum connectivity given by the third argument

    void connectedComponents();                             //method to count connected components
    void DFSUtil(int, vector<bool>&, int);                  //auxiliary function used to compute connected components
    
    void graphStructure();                                  //this method prints the structure of the graph
    void printConfiguration();                              //this method prints the configuration studied
    
    void clearGraph();

    
};



//---------------------------------------------------------------------------------------------------------------------------------------------------------------------// useful functions

//this function allows to fill a vector in one single line
//it has been downloaded from https://gist.github.com/pablomtz/5577626
//examples:
//vector<int> v = make_vector<int>() << 0 << 0 << 1 << 1 << 1 << 0 << 0 << 1 << 1;

template <typename T>
class make_vector {
 public:
  typedef make_vector<T> my_type;
  my_type& operator<< (const T& val) {
    data_.push_back(val);
    return *this;
  }
  operator std::vector<T>() const {
    return data_;
  }
 private:
  std::vector<T> data_;
};

//this function print all the elements of a vector.

void vec_print(vector<int>& vec);


