// g++ main.cpp -lboost_system  -O3 -std=c++0x -fno-strict-aliasing -fopenmp -pedantic -Wall 
#include <boost/version.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/intrusive/rbtree.hpp>
#include <boost/graph/kruskal_min_spanning_tree.hpp>
#include <boost/intrusive/rbtree.hpp>
#include <boost/graph/random.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/graph/graphviz.hpp>

#include <typeinfo>
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



typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexInfo, EdgeInfo> Graph;

typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIterator;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::vertex_iterator vertex_it;
typedef boost::graph_traits<Graph>::edge_iterator edge_it;
typedef boost::property_map<Graph, int EdgeInfo::*>::type WeightPrMap;
typedef boost::property_map<Graph, int VertexInfo::*>::type VertexIdMap;
typedef boost::property_map<Graph, int VertexInfo::*>::type VertexLabelMap;

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexInfo, EdgeInfo> DGraph;
typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIteratorD;
typedef boost::graph_traits<DGraph>::vertex_descriptor VertexD;
typedef boost::graph_traits<DGraph>::edge_descriptor EdgeD;
typedef boost::graph_traits<DGraph>::vertex_iterator vertex_itD;
typedef boost::graph_traits<DGraph>::edge_iterator edge_itD;
typedef std::pair<edge_itD, edge_itD> EdgePairD;

typedef boost::property_map<DGraph, int EdgeInfo::*>::type WeightPrMapD;
typedef boost::property_map<DGraph, int VertexInfo::*>::type VertexIdMapD;
typedef boost::property_map<DGraph, int VertexInfo::*>::type VertexLabelMapD;





void generate_random_graph(Graph& graph);

// Set random weights in the range [min_w, max_w] to all nodes of graph G
void graph_set_weights(Graph& G, WeightPrMap w, int min_w, int max_w);
void graph_set_weightsD(DGraph& G, WeightPrMapD w);
// Set Id to each node of graph G
void graph_set_node_id(Graph& G, VertexIdMap idmap);
void graph_set_node_idD(DGraph& G, VertexIdMapD idmap);
// Prints a list of all vertices, edges, and total number of each for graph G. Includes information stored in them too.
void print_graph(Graph& G);

void spanning_tree(Graph& graph);


class Tree {
    /*
    * idmap: parrent vertex to tree vertex
    * labelmap: tree vertex to parent vertex
    */
    public:
        Tree(Graph& graph, std::vector<int> &component) {
            m_idMap = boost::get(&VertexInfo::id, m_graph);
            m_labelMap = boost::get(&VertexInfo::label, m_graph);
            m_weightMap = boost::get(&EdgeInfo::weight, m_graph);
            
            for (size_t i = 0; i < component.size(); i++) {
                boost::add_vertex(m_graph);
                m_stdIdMap[component[i]] = i;
            }
            set_label_map(component);

            for (size_t i = 0; i < component.size(); i++) {
                AdjacencyIterator ai, ae;
                for (boost::tie(ai, ae) = boost::adjacent_vertices(i, graph); ai != ae; ai++) { 
                    if (m_stdIdMap[*ai] != m_stdIdMap[i] && !boost::edge(m_stdIdMap[i], m_stdIdMap[*ai], m_graph).second) {
                        // LOG("v: " << i << " " << *ai << "|\t" << m_stdIdMap[i] << " " << m_stdIdMap[*ai]); // initial edge to tree edge corr
                        boost::add_edge(m_stdIdMap[i], m_stdIdMap[*ai], m_graph);
                    }
                }
            }

            graph_set_node_id(m_graph, m_idMap);
            graph_set_node_id(m_graph, m_labelMap);
            graph_set_weights(m_graph, m_weightMap, 1, MAX_WEIGHT);

            create_Dgraph();
        }
        ~Tree() {}


        Graph& get_graph() { return m_graph; }
    private:

        void create_Dgraph() {
            std::vector<Edge> spanning_tree;
            boost::kruskal_minimum_spanning_tree(m_graph, std::back_inserter(spanning_tree), weight_map(boost::get(&EdgeInfo::weight, m_graph)));

            m_idMapD = boost::get(&VertexInfo::id, m_Dgraph);
            m_labelMapD = boost::get(&VertexInfo::label, m_Dgraph);
            m_weightMapD = boost::get(&EdgeInfo::weight, m_Dgraph);

            for (vertex_it vi = boost::vertices(m_graph).first; vi != boost::vertices(m_graph).second; vi++) {
                boost::add_vertex(m_Dgraph);
            }

            for (vertex_it vi = boost::vertices(m_graph).first; vi != boost::vertices(m_graph).second; vi++) {
                AdjacencyIterator ai, ae;
                for (boost::tie(ai, ae) = boost::adjacent_vertices(*vi, m_graph); ai != ae; ai++) { 
                    if (*vi != *ai && !boost::edge(*vi, *ai, m_Dgraph).second) {
                        boost::add_edge(*vi, *ai, m_Dgraph);
                        boost::add_edge(*ai, *vi, m_Dgraph);
                    }

                    if (std::find(spanning_tree.begin(), spanning_tree.end(), boost::edge(*vi, *ai, m_graph).first) == spanning_tree.end()) {
                        auxiliaryEdges.push_back(boost::edge(*vi, *ai, m_graph).first);
                    }
                }
                boost::add_edge(*vi, *vi, m_Dgraph);
            }

            LOG("Auxiliary Edges");
            for (size_t i = 0; i < auxiliaryEdges.size(); i++) {
                LOG(boost::source(auxiliaryEdges[i], m_Dgraph) << " " << boost::target(auxiliaryEdges[i], m_Dgraph));
            }

            LOG("\n Vertices:");
            for (vertex_itD it = boost::vertices(m_Dgraph).first; it != boost::vertices(m_Dgraph).second; it++) {
                LOG(m_stdLabelMap[*it]);
            }

            LOG("\n Edges:");
            VertexD u, v;
            for (EdgePairD ep = boost::edges(m_Dgraph); ep.first != ep.second; ++ep.first) {
                u = boost::source(*ep.first, m_Dgraph);
                v = boost::target(*ep.first, m_Dgraph);
                LOG(m_stdLabelMap[u] << " => " << m_stdLabelMap[v]);
            }

            graph_set_node_idD(m_Dgraph, m_idMapD);
            graph_set_node_idD(m_Dgraph, m_labelMapD);
            graph_set_weightsD(m_Dgraph, m_weightMapD);

        }

        void set_label_map(std::vector<int>& vec) {
            LOG("label map");
            for (vertex_it vi = boost::vertices(m_graph).first; vi != boost::vertices(m_graph).second; vi++) {
                m_stdLabelMap[*vi] = vec[*vi];
            }
        }

    public:
        Graph m_graph;
        WeightPrMap m_weightMap;
        VertexIdMap m_idMap;
        VertexLabelMap m_labelMap;

        std::map<int, int> m_stdIdMap;
        std::map<int, int> m_stdLabelMap;
        
        DGraph m_Dgraph;
        WeightPrMapD m_weightMapD;
        VertexIdMapD m_idMapD;
        VertexLabelMapD m_labelMapD;
  
        std::vector<Edge> auxiliaryEdges;
};





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

    srand(std::time(nullptr));
    Graph randomGraph;

    WeightPrMap weightMap = get(&EdgeInfo::weight, randomGraph);
    VertexIdMap idMap = get(&VertexInfo::id, randomGraph);
    VertexLabelMap labelMap = get(&VertexInfo::label, randomGraph);
    
    generate_random_graph(randomGraph);

    graph_set_node_id(randomGraph, idMap);
    graph_set_node_id(randomGraph, labelMap);
    graph_set_weights(randomGraph, weightMap, 1, MAX_WEIGHT);
    

    spanning_tree(randomGraph);    


    return EXIT_SUCCESS;
}

std::vector<std::vector<int>> get_components(std::vector<int> component) {
    std::vector<std::vector<int>> result(*std::max_element(std::begin(component), std::end(component)) + 1);

    for (size_t i = 0; i < component.size(); i++) {
        result[component[i]].emplace_back(i);
    }
    
    return result;
}

void spanning_tree(Graph& graph) {
    std::vector<int> component(boost::num_vertices(graph));
    int num = boost::connected_components(graph, &component[0]);

    std::vector<std::vector<int>> component_vertex = get_components(component);

    for (size_t i = 0; i < component_vertex.size(); i++) {
        std::cout << i << ": ";
        for (size_t j = 0; j < component_vertex[i].size(); j++) {
            std::cout << component_vertex[i][j] << " ";
        }
        std::cout << std::endl;
    }

    Tree* a = new Tree(graph, component_vertex[0]);
    LOG("tree created");
    print_graph(a->get_graph());
    LOG("Abnaroz");
    delete a;

    // Tree a(graph, component_vertex[0]);

    print_graph(graph);

    return; 
}


void print_graph(Graph& G) {
    std::ofstream stream;
    stream.open("outputfile");

    boost::dynamic_properties dp;

    dp.property("node_id", boost::get(&VertexInfo::id, G));
    dp.property("label", boost::get(&VertexInfo::label, G));
    dp.property("weight", boost::get(&EdgeInfo::weight, G));
    dp.property("label", boost::get(&EdgeInfo::weight, G));
    
    write_graphviz_dp(stream, G, dp);
    
    stream.close();

    system("dot -Tsvg ./outputfile > test.svg");
    // system("eog test.svg");
    return;
}


void generate_random_graph(Graph& graph) {
    std::srand(std::time(nullptr));

    for (int i = 0; i < VERTICES; i++) {
        boost::add_vertex(graph);
    }

    size_t edge_num = 0;
    volatile bool flag = false;


    #pragma omp parallel for shared(flag)
    for (vertex_it it = boost::vertices(graph).first; it != boost::vertices(graph).second; it++) {
        if ( RAND_FLOAT() > 0.4f ) {
            for (vertex_it vi = boost::vertices(graph).first; vi != boost::vertices(graph).second; vi++) {
                if (flag) {
                    continue;
                }

                if (*it != *vi && !boost::edge(*it, *vi, graph).second && RAND_FLOAT() > 0.5f) {
                    #pragma omp critical (add_edge)
                    {   
                        edge_num++;
                        if (edge_num >= EDGES) {
                            flag = true;
                        }

                        boost::add_edge(*it, *vi, graph);
                        // boost::add_edge(*vi, *it, graph);
                    }    
                }
            }
        }
        
        #pragma omp critical (add_edge)
        boost::add_edge(*it, *it, graph);
    }
}

boost::random::mt19937 rng_gen(){
    boost::mt19937 rng;
    rng.seed(uint32_t (std::time(0)));
    // rng.seed(uint32_t (0));
    boost::uniform_int<> u(0, 1000);
    boost::variate_generator<boost::mt19937&, boost::uniform_int<> > rng1(rng, u);
    return rng;
    
}

void graph_set_weights(Graph& G, WeightPrMap w, int min_w, int max_w){
    std::time_t now = std::time(0);
    boost::random::mt19937 gen{static_cast<std::uint32_t>(now)};
    boost::random::uniform_int_distribution<> dist{min_w, max_w};
    edge_it ei;
    for (ei = boost::edges(G).first; ei != boost::edges(G).second; ei++){
        // w[*ei] = dist(gen);
        w[*ei] = 1;
    }
    return;
}

void graph_set_weightsD(DGraph& G, WeightPrMapD w) {
    
    edge_itD ei;
    for (ei = boost::edges(G).first; ei != boost::edges(G).second; ei++){
        // w[*ei] = dist(gen);
        w[*ei] = 1;
    }
    return;
}


void graph_set_node_id(Graph& G, VertexIdMap idmap) {
    for (vertex_it vi = boost::vertices(G).first; vi != boost::vertices(G).second; vi++) {
        idmap[*vi] = *vi;
    }
}

void graph_set_node_idD(DGraph& G, VertexIdMapD idmap) {
    for (vertex_itD vi = boost::vertices(G).first; vi != boost::vertices(G).second; vi++) {
        idmap[*vi] = *vi;
    }
}