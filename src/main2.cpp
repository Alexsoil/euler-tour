// g++ main.cpp -lboost_system  -O3 -std=c++0x -fno-strict-aliasing -fopenmp -pedantic -Wall 
#include <boost/version.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/intrusive/rbtree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>

#include <boost/graph/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/graphviz.hpp>

#include <iostream>
#include <tuple>
#include <vector>
#include <queue>
#include <random>
#include <omp.h>
#include <stdlib.h>


#define VERTICES 10
#define EDGES 10
#define MAX_WEIGHT 200

#define THREADS 2

#define RAND_FLOAT() ((float) std::rand() / (float) INT_MAX)

#define DEBUG

#ifdef DEBUG
    #define LOG(X) (std::cout << X << std::endl)
#else
    #define LOG(X)
#endif


struct EdgeInfo{
    int weight;
};

struct VertexInfo{
    int label;
    int id;
};


typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexInfo, EdgeInfo> DGraph;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexInfo, EdgeInfo> Graph;
typedef boost::graph_traits<DGraph>::adjacency_iterator AdjacencyIterator;

typedef boost::graph_traits<DGraph>::vertex_descriptor Vertex;
typedef boost::graph_traits<DGraph>::edge_descriptor Edge;
typedef boost::graph_traits<DGraph>::vertex_iterator vertex_it;
typedef boost::graph_traits<DGraph>::edge_iterator edge_it;
typedef boost::property_map<DGraph, int EdgeInfo::*>::type WeightPrMap;
typedef boost::property_map<DGraph, int VertexInfo::*>::type VertexIdMap;
typedef boost::property_map<DGraph, int VertexInfo::*>::type VertexLabelMap;


boost::random::mt19937 rng_gen();
// Set random weights in the range [min_w, max_w] to all nodes of graph G
void graph_set_weights(Graph& G, WeightPrMap w, int min_w, int max_w);
// Set Id to each node of graph G
void graph_set_node_id(Graph& G, VertexIdMap idmap);
// Prints a list of all vertices, edges, and total number of each for graph G. Includes information stored in them too.
void print_graph(DGraph& G);

void spanning_tree(Graph& graph);

int main(void){

    #ifdef _OPENMP
        (void) omp_set_dynamic(0);
        if(omp_get_dynamic()){
            printf("warning: dynamic adjustment of threads has been set\n");
        }
        (void) omp_set_num_threads(THREADS);
    #endif

    LOG("\n: ^)");
    LOG("Using Boost " << BOOST_VERSION / 10000 << "." << BOOST_VERSION / 100 % 1000 << "." << BOOST_VERSION % 100);

    boost::random::mt19937 rng = rng_gen();
    srand(std::time(nullptr));
    DGraph randomGraph;
    generate_random_graph(randomGraph, 5, 3, rng, false, false);
    edge_it ei;
    for (ei = boost::edges(randomGraph).first; ei != boost::edges(randomGraph).second; ei++){
        std::cout << "Abnaroz" << std::endl;
        Edge e = boost::add_edge(boost::target(*ei, randomGraph), boost::source(*ei, randomGraph), randomGraph).first;
    }
    print_graph(randomGraph);
    return EXIT_SUCCESS;
}

void print_graph(DGraph& G){
    vertex_it vi;
    for (vi = boost::vertices(G).first; vi != boost::vertices(G).second; vi++){
        std::cout << "ID: " << *vi << std::endl;
    }
    edge_it ei;
    for (ei = boost::edges(G).first; ei != boost::edges(G).second; ei++){
        std::cout << "Edge: " << boost::source(*ei, G) << " -> " << boost::target(*ei, G) << std::endl;
    }
    std::cout << boost::num_vertices(G) << " vertices and " << boost::num_edges(G) << " edges" << std::endl;
    return;
}

std::vector<std::vector<int>> get_components(std::vector<int> component) {
    std::vector<std::vector<int>> result(*std::max_element(std::begin(component), std::end(component)) + 1);

    for (size_t i = 0; i < component.size(); i++) {
        result[component[i]].emplace_back(i);
    }
    
    return result;
}

boost::random::mt19937 rng_gen(){
    boost::mt19937 rng;
    rng.seed(uint32_t (std::time(0)));
    // rng.seed(uint32_t (0));
    boost::uniform_int<> u(0, 1000);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rng1(rng, u);
    return rng;
    
}
