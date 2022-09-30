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

// #include <LEDA/core/dynamic_trees.h>

#include <typeinfo>
#include <iostream>
#include <tuple>
#include <vector>
#include <stack>
#include <random>
#include <omp.h>
#include <stdlib.h>
#include <string>


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


// Undirected Graph/Vertices/Edges/Iterators/Maps
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexInfo, EdgeInfo> Graph;
typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIterator;
typedef boost::graph_traits<Graph>::vertex_descriptor Vertex;
typedef boost::graph_traits<Graph>::edge_descriptor Edge;
typedef boost::graph_traits<Graph>::vertex_iterator vertex_it;
typedef boost::graph_traits<Graph>::edge_iterator edge_it;
typedef boost::property_map<Graph, int EdgeInfo::*>::type WeightPrMap;
typedef boost::property_map<Graph, int VertexInfo::*>::type VertexIdMap;
typedef boost::property_map<Graph, int VertexInfo::*>::type VertexLabelMap;

// Directed Graph/Vertices/Edges/Iterators/Maps
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, VertexInfo, EdgeInfo> DGraph;
typedef boost::graph_traits<Graph>::adjacency_iterator AdjacencyIteratorD;
typedef boost::graph_traits<DGraph>::vertex_descriptor VertexD;
typedef boost::graph_traits<DGraph>::edge_descriptor EdgeD;
typedef boost::graph_traits<DGraph>::vertex_iterator vertex_itD;
typedef boost::graph_traits<DGraph>::edge_iterator edge_itD;
typedef boost::graph_traits<DGraph>::out_edge_iterator out_edge_itD;
typedef std::pair<edge_itD, edge_itD> EdgePairD;
typedef boost::property_map<DGraph, int EdgeInfo::*>::type WeightPrMapD;
typedef boost::property_map<DGraph, int VertexInfo::*>::type VertexIdMapD;
typedef boost::property_map<DGraph, int VertexInfo::*>::type VertexLabelMapD;

// Used in Euler Path Creation
int getLastElementNotEqualValue(std::vector<int>& vec, int value) {
    if (!vec.empty()) {
        return vec.back() != value;
    } else {
        return true;
    }
}

//Generate the Undirected graph
void generate_random_graph(Graph& graph);

// Set random weights in the range [min_w, max_w] to all nodes of graph G
void graph_set_weights(Graph& G, WeightPrMap w, int min_w, int max_w);
void graph_set_weightsD(DGraph& G, WeightPrMapD w);
// Set Id to each node of graph G
void graph_set_node_id(Graph& G, VertexIdMap idmap);
void graph_set_node_idD(DGraph& G, VertexIdMapD idmap);
// Prints a list of all vertices, edges, and total number of each for graph G. Includes information stored in them too.
void print_graph(Graph& G);

template <typename T>
void print_vector(std::vector<T>& vec, const std::string& prefix, const std::string& postfix) {
    std::cout << prefix;
    for (size_t i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << postfix;
    }
    std::cout << std::endl;
}


class Tree {
    /*
    * 
    */
    public:
        //Constructor (Undirected graph, integer vector of graph connected components)
        Tree(Graph& graph, std::vector<int> &component) {
            //Maps for Vertices/Edges of Undirected graph
            m_idMap = boost::get(&VertexInfo::id, m_graph);
            m_labelMap = boost::get(&VertexInfo::label, m_graph);
            m_weightMap = boost::get(&EdgeInfo::weight, m_graph);
            
            // Map graph vertices to m_graph vertices
            for (size_t i = 0; i < component.size(); i++) {
                boost::add_vertex(m_graph);
                m_stdIdMap[component[i]] = i;
            }
            // Set std_labelMap
            set_label_map(component);

            for (size_t i = 0; i < component.size(); i++) {
                //Iterate adjacent edges to vertex
                AdjacencyIterator ai, ae;
                //Iterate through each vertex in this component's edges
                for (boost::tie(ai, ae) = boost::adjacent_vertices(i, graph); ai != ae; ai++) {
                    //Insert all edges that belong to m_graph
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


            create_euler_path();
        }
        // Destructor
        ~Tree() {}

        // Prints a vector
        template <typename T>
        void print_tree_vector(std::vector<T>& vec, const std::string& prefix, const std::string& postfix) {
            std::cout << prefix;
            for (size_t i = 0; i < vec.size(); i++) {
                std::cout << m_stdLabelMap[vec[i]] << postfix;
            }
            std::cout << std::endl;
        }


        //eulerPath is a vector of vertices
        void printEulerPath() {
            print_tree_vector<int>(m_eulerPath, "Euler Path", " ");
        }

        Graph& get_graph() { return m_graph; }
    private:

        //Generate Directed graph from Undirected
        void create_Dgraph() {
            //Undirected edge vector to store edges of spanning tree
            std::vector<Edge> spanning_tree;
            //Kruskal spanning tree finder (back_inserter for simplicity)
            boost::kruskal_minimum_spanning_tree(m_graph, std::back_inserter(spanning_tree), weight_map(boost::get(&EdgeInfo::weight, m_graph)));

            //Maps for Directed graph
            m_idMapD = boost::get(&VertexInfo::id, m_Dgraph);
            m_labelMapD = boost::get(&VertexInfo::label, m_Dgraph);
            m_weightMapD = boost::get(&EdgeInfo::weight, m_Dgraph);

            //Create a directed graph vertex for each undirected graph vertex
            for (vertex_it vi = boost::vertices(m_graph).first; vi != boost::vertices(m_graph).second; vi++) {
                boost::add_vertex(m_Dgraph);
            }
            //Doing the same thing as for Directed but 2 edges back and forth
            for (vertex_it vi = boost::vertices(m_graph).first; vi != boost::vertices(m_graph).second; vi++) {
                AdjacencyIterator ai, ae;
                for (boost::tie(ai, ae) = boost::adjacent_vertices(*vi, m_graph); ai != ae; ai++) { 
                    if (*vi != *ai && !boost::edge(*vi, *ai, m_Dgraph).second) {
                        boost::add_edge(*vi, *ai, m_Dgraph);
                        boost::add_edge(*ai, *vi, m_Dgraph);
                    }
                    //If edge doesn't exist in spanning tree, it's auxiliary and is stored elsewhere
                    if (std::find(spanning_tree.begin(), spanning_tree.end(), boost::edge(*vi, *ai, m_graph).first) == spanning_tree.end()) {
                        auxiliaryEdges.push_back(boost::edge(*vi, *ai, m_graph).first);
                        auxiliaryEdgesD.push_back(boost::edge(*vi, *ai, m_Dgraph).first);
                    }
                }
                //Self edge for Directed
                boost::add_edge(*vi, *vi, m_Dgraph);
            }

            LOG("Auxiliary Edges");
            //Print auxiliary Edges
            for (size_t i = 0; i < auxiliaryEdges.size(); i++) {
                LOG(boost::source(auxiliaryEdges[i], m_Dgraph) << " " << boost::target(auxiliaryEdges[i], m_Dgraph));
            }

            LOG("\n Vertices:");
            //Print vertices
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

        void create_euler_path() {
            
            std::stack<int> s;
            std::vector<EdgeD> visitedEdges;
            std::vector<int> path;

            // Stack for vertices to be examined
            s.push(0);
            // Always start at the first vertex of the tree (arbitrary      )
            int root = 0;
            out_edge_itD ef, ee;

            std::map<int, int> node_roots;

            while(!s.empty()) {
                // Take last inserted vertex
                int vertex = s.top();
                s.pop();
                std::vector<EdgeD> temp;
                
                // Iterate through its edges
                for (boost::tie(ef, ee) = out_edges(vertex, m_Dgraph);  ef != ee; ef++) {
                    // If the edge is not auxiliary, it has not already been visited, and it's not an edge going to the parent vertex, insert to Euler Tour
                    if ((std::find(auxiliaryEdgesD.begin(), auxiliaryEdgesD.end(), *ef) == auxiliaryEdgesD.end()) && 
                        (std::find(visitedEdges.begin(), visitedEdges.end(), *ef) == visitedEdges.end()) && 
                        getLastElementNotEqualValue(path, boost::target(*ef, m_Dgraph)) ) {
                        temp.push_back(*ef);
                    }
                }

                // Do not add the parent if you are root, or if the parent has already been added to the stack
                if (!path.empty() && node_roots.find(vertex) == node_roots.end() && vertex != root) {
                    // LOG(path.back());
                    node_roots[vertex] = path.back();
                    s.push(path.back());
                }

                // If self edge has not beend accounted for, do so
                if (std::find(visitedEdges.begin(), visitedEdges.end(), boost::edge(vertex, vertex, m_Dgraph).first) == visitedEdges.end()) {
                    s.push(vertex);
                }

                //If adjacent vertices are not ourself or parent, add vertices to the stack
                for (size_t i = 0; i < temp.size(); i++) {
                    // LOG(m_stdLabelMap[boost::source(temp[i], m_Dgraph)] << " " << m_stdLabelMap[boost::target(temp[i], m_Dgraph)]);
                    // LOG("=");
                    int v = boost::target(temp[i], m_Dgraph);
                    if (v != vertex && v != node_roots[vertex]) {
                        s.push(v);
                    }
                    visitedEdges.push_back(temp[i]);
                }
                
                // Vertex checks done, add to Euler Path
                path.push_back(vertex);

                // std::cout << "Stack: ";
                // printstack(s, m_stdLabelMap);
                        
            }

            std::cout << "Path: ";
            for (size_t i = 0; i < path.size(); i++) {
                std::cout << m_stdLabelMap[path[i]] << " ";
            }
            std::cout << std::endl;

            std::copy(path.begin(), path.end(), std::back_inserter(m_eulerPath));
            return;
        }

    public:
        // Undirected data
        Graph m_graph;
        WeightPrMap m_weightMap;
        VertexIdMap m_idMap;
        VertexLabelMap m_labelMap;

        // idmap: parrent vertex to tree vertex
        std::map<int, int> m_stdIdMap;
        // labelmap: tree vertex to parent vertex
        std::map<int, int> m_stdLabelMap;
        
        // Directed data
        DGraph m_Dgraph;
        WeightPrMapD m_weightMapD;
        VertexIdMapD m_idMapD;
        VertexLabelMapD m_labelMapD;
        
        // Aux edges (undirected/directed)
        std::vector<Edge> auxiliaryEdges;
        std::vector<EdgeD> auxiliaryEdgesD;

        // Euler path (vertices!)
        std::vector<int> m_eulerPath;
};


int reRootEuler(int vertex, std::vector<int>& eulerPath);
int link(std::vector<std::vector<int>>& dynamicTreeVec, int u, int v);
int cut(std::vector<std::vector<int>>& dynamicTreeVec, int u, int v);

std::vector<Tree*> spanning_tree(Graph& graph);
Tree* getTreeFromVertex(std::vector<Tree*>* treeVec, int vertex);
int getTreeIndexFromVertex(std::vector<std::vector<int>>& dynamicTreeVec, int vertex);




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
    //Generate Graph
    generate_random_graph(randomGraph);

    graph_set_node_id(randomGraph, idMap);
    graph_set_node_id(randomGraph, labelMap);
    graph_set_weights(randomGraph, weightMap, 1, MAX_WEIGHT);
    

    print_graph(randomGraph);


    spanning_tree(randomGraph);    


    return EXIT_SUCCESS;
}

// Return as many vectors as are components, with each component's vertices each
std::vector<std::vector<int>> get_components(std::vector<int> component) {
    // 2d vector
    std::vector<std::vector<int>> result(*std::max_element(std::begin(component), std::end(component)) + 1);

    for (size_t i = 0; i < component.size(); i++) {
        result[component[i]].emplace_back(i);
    }
    
    return result;
}

void printstack(std::stack<int> s, std::map<int, int> lmap) {
    while (!s.empty()) {
        std::cout << lmap[s.top()] << " ";
        s.pop();
    }
    std::cout << std::endl;
}

std::vector<Tree*> spanning_tree(Graph& graph) {
    // Get ids of connected components for each vertex
    std::vector<int> component(boost::num_vertices(graph));
    int num = boost::connected_components(graph, &component[0]);
    
    // Vector of component vectors with each having vertices of a specific component
    std::vector<std::vector<int>> component_vertex = get_components(component);

    // Vector of Tree obejcts
    std::vector<Tree*> tree_vector;

    for (size_t i = 0; i < component_vertex.size(); i++) {
        std::cout << i << ": ";
        for (size_t j = 0; j < component_vertex[i].size(); j++) {
            std::cout << component_vertex[i][j] << " ";
        }
        std::cout << std::endl;
        
        // Create trees
        tree_vector.push_back(new Tree(graph, component_vertex[i]));
        
    }

    std::vector<std::vector<int>> dynamicTrees;

    // For each tree
    for (Tree* a : tree_vector) {
        a->printEulerPath();
        
        std::vector<int> eulerpath;
        // Iterate through tree's euler tour and insert the original graph's vertices related to the dynamic euler tour
        for (std::vector<int>::iterator it = a->m_eulerPath.begin(); it != a->m_eulerPath.end(); ++it) {
            eulerpath.push_back(a->m_stdLabelMap[*it]);
        }
        
        dynamicTrees.push_back(eulerpath);

        std::cout << "____________" << std::endl;
    }


    LOG(" ");
    print_vector<int>(dynamicTrees[0], "dt0: ", " ");
    // Reroot to the 2nd vertex
    reRootEuler(2, dynamicTrees[0]);
    print_vector<int>(dynamicTrees[0], "dt0: ", " ");


    // link
    if (dynamicTrees.size() >= 2) {
        // Try to link dynamic trees 1 and 2
        LOG("Dynamic Trees size: " << dynamicTrees.size());
        print_vector<int>(dynamicTrees[0], "dt0: ", " ");
        print_vector<int>(dynamicTrees[1], "dt1: ", " ");
        if (!link(dynamicTrees, dynamicTrees[0][0], dynamicTrees[1][0])) {
            print_vector<int>(dynamicTrees[0], "dt0: ", " ");
            LOG("Dynamic Trees size: " << dynamicTrees.size());
        } else {
            LOG("ERROR on link");
        }
        // Try to cut dynamic tree
        if (!cut(dynamicTrees, dynamicTrees[0][1], dynamicTrees[0][2])) {
            LOG("Dynamic Trees size: " << dynamicTrees.size());
            for (size_t i = 0; i < dynamicTrees.size(); i++) {
                print_vector<int>(dynamicTrees[i], "dt [i] ", " ");
            }
        } else {
            LOG("ERROR on cut");
        }

    }
    // Free Memory
    for (Tree* a : tree_vector) {
        delete a;
    }

    return tree_vector; 
}

// Print vector image of original graph
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
    //Add vertices as defined
    for (int i = 0; i < VERTICES; i++) {
        boost::add_vertex(graph);
    }

    size_t edge_num = 0;
    //volatile for openMP
    volatile bool flag = false;
    //Parallelized for more speed
    #pragma omp parallel for shared(flag)
    for (vertex_it it = boost::vertices(graph).first; it != boost::vertices(graph).second; it++) {
        if ( RAND_FLOAT() > 0.4f ) {
            //If this is true, this vertex will have edges
            for (vertex_it vi = boost::vertices(graph).first; vi != boost::vertices(graph).second; vi++) {
                if (flag) { // If flag is true, we stop adding edges
                    continue;
                }
                
                //If this vertex is not pointing at itself and edge doesn't already exist
                if (*it != *vi && !boost::edge(*it, *vi, graph).second && RAND_FLOAT() > 0.5f) {
                    #pragma omp critical (add_edge)
                    {   
                        edge_num++;
                        if (edge_num >= EDGES) {
                            flag = true;
                        }

                        boost::add_edge(*it, *vi, graph);
                    }    
                }
            }
        }
        
        //Generate Self Edge
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


//All weights are 1 and are only used in the Minimum Spanning Tree calculation algorithm
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

//Not Relevant
void graph_set_weightsD(DGraph& G, WeightPrMapD w) {
    edge_itD ei;
    for (ei = boost::edges(G).first; ei != boost::edges(G).second; ei++){
        // w[*ei] = dist(gen);
        w[*ei] = 1;
    }
    return;
}

// Sets idmap
void graph_set_node_id(Graph& G, VertexIdMap idmap) {
    for (vertex_it vi = boost::vertices(G).first; vi != boost::vertices(G).second; vi++) {
        idmap[*vi] = *vi;
    }
}

// Sets idmap for directed graphs
void graph_set_node_idD(DGraph& G, VertexIdMapD idmap) {
    for (vertex_itD vi = boost::vertices(G).first; vi != boost::vertices(G).second; vi++) {
        idmap[*vi] = *vi;
    }
}

// Reroots the Euler Tour at specified vertex
int reRootEuler(int vertex, std::vector<int>& eulerPath) {
    std::vector<int>::iterator newRoot = std::find(eulerPath.begin(), eulerPath.end(), vertex);

    // it's has already vertex as root
    if(vertex == eulerPath[0]) {
        return 0;
    }

    if (newRoot != eulerPath.end()) {
        // Cut at the root
        std::vector<int> A (eulerPath.begin(), newRoot);
        std::vector<int> B (newRoot, eulerPath.end());

        A.erase(A.begin());// remove first element of A
        A.push_back(vertex);
        
        print_vector<int>(A, "A: ", " ");
        print_vector<int>(B, "B: ", " ");

        eulerPath.clear();
        // Concatenate
        eulerPath.insert(eulerPath.end(), B.begin(), B.end());
        eulerPath.insert(eulerPath.end(), A.begin(), A.end());

        //success
        return 0;
    } else {
        //fail (vertex not in spanning tree)
        return 1;
    }
}

// Link two vertices and Euler Tours
int link(std::vector<std::vector<int>>& dynamicTreeVec, int u, int v) {
    int uIndex = getTreeIndexFromVertex(dynamicTreeVec, u);
    int vIndex = getTreeIndexFromVertex(dynamicTreeVec, v);
    // If vertex not in Spanning Tree
    if (uIndex == -1 || vIndex == -1) {
        return 1; // error
    }

    std::vector<int>::iterator uIt = std::find(dynamicTreeVec[uIndex].begin(), dynamicTreeVec[uIndex].end(), u);
    std::vector<int>::iterator vIt = std::find(dynamicTreeVec[vIndex].begin(), dynamicTreeVec[vIndex].end(), v);

    // Possibly unnecessary but too afraid to remove it
    if (uIt != dynamicTreeVec[uIndex].end() && vIt != dynamicTreeVec[vIndex].end()) {
        // Throw error if cannot reroot
        if(reRootEuler(u, dynamicTreeVec[uIndex])) return 1;
        if(reRootEuler(v, dynamicTreeVec[vIndex])) return 1;

        // Add Euler Tour the v belongs to, to Euler Tour that u belongs to
        dynamicTreeVec[uIndex].insert(dynamicTreeVec[uIndex].end(), dynamicTreeVec[vIndex].begin(), dynamicTreeVec[vIndex].end());
        dynamicTreeVec[uIndex].push_back(u);
        
        dynamicTreeVec.erase(dynamicTreeVec.begin() + vIndex);

        return 0; // success
    } else {
        return 1;
    }
}

// Cut edge between u and v
int cut(std::vector<std::vector<int>>& dynamicTreeVec, int u, int v) {
    int uIndex = getTreeIndexFromVertex(dynamicTreeVec, u);
    int vIndex = getTreeIndexFromVertex(dynamicTreeVec, v);
    
    if (uIndex != vIndex || u == v) { // each vertex is in diff tree or u=v
        return 1;
    }

    // Search for u in tree, while u is found but u+1 is not v, search for next occurance of u
    std::vector<int>::iterator uvEdge = std::find(dynamicTreeVec[uIndex].begin(), dynamicTreeVec[uIndex].end(), u);
    while(uvEdge != dynamicTreeVec[uIndex].end() && *(uvEdge + 1) != v) {
        uvEdge = std::find(uvEdge + 1, dynamicTreeVec[uIndex].end(), u);
    }

    // Search for v in tree, while v is found but v+1 is not u, search for next occurance of v
    std::vector<int>::iterator vuEdge = std::find(dynamicTreeVec[uIndex].begin(), dynamicTreeVec[uIndex].end(), v);
    while(vuEdge != dynamicTreeVec[uIndex].end() && *(vuEdge + 1) != u) {
        vuEdge = std::find(vuEdge + 1, dynamicTreeVec[uIndex].end(), v);
    }

    // If either edge is not found, error out
    if (vuEdge == dynamicTreeVec[uIndex].end() || uvEdge == dynamicTreeVec[uIndex].end() ) {
        return 1;
    }

    std::vector<int> J(dynamicTreeVec[uIndex].begin(), uvEdge + 1);
    std::vector<int> K(uvEdge + 1, vuEdge + 1);
    std::vector<int> L(vuEdge + 1, dynamicTreeVec[uIndex].end());
    print_vector(J, "J: ", " ");
    print_vector(K, "K: ", " ");
    print_vector(L, "L: ", " ");

    J.pop_back();

    std::vector<int> E1 = K;
    std::vector<int> E2 = J;

    // Add L
    E2.insert(E2.end(), L.begin(), L.end());

    
    // Remove old tree
    dynamicTreeVec.erase(dynamicTreeVec.begin() + uIndex);
    // Add resulting trees
    dynamicTreeVec.push_back(E1);
    dynamicTreeVec.push_back(E2);

    return 0;
}

// Returns The tree object that vertex belongs to
Tree* getTreeFromVertex(std::vector<Tree*>& treeVec, int vertex) {
    for(Tree* tree : treeVec) {
        if (tree->m_stdIdMap.find(vertex) != tree->m_stdIdMap.end()) {
            return tree;
        }
    }
    return nullptr;
}

// Returns the index of Tree that vertex belongs to
int getTreeIndexFromVertex(std::vector<std::vector<int>>& dynamicTreeVec, int vertex) {
    for (size_t i = 0; i < dynamicTreeVec.size(); i++) {
        if (std::find(dynamicTreeVec[i].begin(), dynamicTreeVec[i].end(), vertex) != dynamicTreeVec[i].end()) {
            return i;
        }
    }
    return -1;
}