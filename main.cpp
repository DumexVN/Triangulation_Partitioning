//-----------------------------------------------------------------------------
// Author   : N Vu
// Program  : Random Aggregation Clustering Type III.a
// File: Main
//-----------------------------------------------------------------------------
// see readme.txt for more details

#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include <limits.h>
#include <random>
#include <time.h>

#include "util.h"

#include <chrono>
#include <process.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

//define typedef
typedef unsigned long int uint;
typedef std::pair<uint, uint> edge;
typedef std::vector<uint> adj_list;
typedef std::vector<adj_list> Graph;

char * outPath = NULL;

//initilise random device
 //using mersene twister
std::random_device rd;
std::mt19937 gen(rd());

/// SAVE RESULT AS: VertexID\tCommunity (seperated by tab)
/// \brief saveResult
/// \param component
///
void saveResult(const std::vector<uint> &component)
{
    std::ofstream myfile;
    if (outPath == NULL)
    {
        printf("Output Path is not defined!");
        return;
    }
    myfile.open (outPath);
    for (size_t i = 0; i < component.size(); i++)
    {
        myfile << i << '\t' << component[i] << '\n';
    }
    myfile.close();
}

///
/// \brief postprocess
/// \param G: original graph G
/// \param hierarchy: points_to edge
///
void postprocess(const std::vector<edge> &hierarchy, const unsigned int &V)
{
    typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
    BoostGraph BG;
    for(size_t i = 0; i < V; i++)   boost::add_vertex(BG);

    for(auto &e : hierarchy)
    {
        boost::add_edge(e.first, e.second, BG);
    }
    //result if the connected component
    printf("Boost Graph V-%d, E-%d", boost::num_vertices(BG), boost::num_edges(BG));
    std::vector<uint> component(boost::num_vertices(BG));
    uint num = boost::connected_components(BG, &component[0]);
    printf("\nNumber of Community: %lu \n", num);
    saveResult(component);
}

///
/// \brief process: do aggregation III.a
/// \param G
///
bool process(const Graph &G, std::vector<edge> &points_to)
{
    printf("\n Processing ... \n");

    for (size_t i = 0; i < G.size(); i++)
    {
        //get main vertex v adjacency
        adj_list v_adj = G[i];
        uint max = 0;
        std::vector<uint> dup_max;
        for (size_t j = 0; j < v_adj.size(); j++)
        {
            //get neighbour vertex u adjacency
            adj_list u_adj = G[v_adj[j]];
            unsigned int comm = Util::NumberOfCommonElement(v_adj, u_adj); //number of common
            if (comm == max)    dup_max.push_back(v_adj[j]);
            else if (comm > max)
            {
                max = comm;
                dup_max.clear();
                dup_max.push_back(v_adj[j]);
            }
        }
        //points to the highest common neighbours
        // if there are more than 1, pick 1 at random
        if (dup_max.size() == 1)    points_to.push_back(edge(i, dup_max[0]));
        else
        {
            std::uniform_int_distribution<uint> dis(0, dup_max.size()-1);
            uint ran = dis(gen);
            points_to.push_back(edge(i,dup_max[ran])); //return the index of chosen vertex
        }
    }
    //DONE
    return true;
}


///
/// \brief build_graph: build graph from input
/// \param raw_edge: std::vector of lines read from input
/// \param myGraph: empty Graph to write
///
bool build_graph(const std::vector<edge> &raw_edge, const uint &MAX_V, Graph &myGraph)
{
    printf("Building Graph ... \n");
    for(uint i = 0; i <= MAX_V; i++)
    {
        adj_list empty;
        myGraph.push_back(empty);
    }

    for(uint i = 0; i < raw_edge.size(); i++)
    {
        edge e = raw_edge.at(i);
        uint from = e.first, to = e.second;
        assert(from < (MAX_V+1) && "ASSERTION FAILED! INDEX OUT OF BOUND");
        myGraph[from].push_back(to); //get adj_list of source vertex
        //check for duplicate neighbour
        //if (std::find(v.begin(), v.end(),value)!=v.end())
    }
    //sort each adj_list
    uint sum_d = 0;
    for (auto &adj : myGraph)
    {
        std::sort(adj.begin(), adj.end());
        sum_d += adj.size();
    }
    printf("Graph done ... \n");
    printf("Vertices: %d \t Sum Degree: %lu \n", myGraph.size(), sum_d);
    printf("Proceed to processing ... \n");
    //process the graph
    return true;
}


///
/// \brief read_graph: Read in edge.txt file
/// it is assumed that each line contains source'\t'target (seperated by tab)
/// it is assumed that index starts from 0
///
bool read_input(const char * filepath, std::vector<edge> &raw_edge, uint &MAX_V)
{
    uint MIN_V = 0;
    std::ifstream infile(filepath);
    if (!infile.good())
    {
        printf("File Erorr! Terminating ...");
        return false;
    }
    std::string line;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        uint a, b;
        if (!(iss >> a >> b) || (a == b))
        {
            printf("Error While Reading Input! Terminating ...");
            return false;
        } // error
        // process pair (a,b)
        if (std::max(a,b) > MAX_V)  MAX_V = std::max(a,b);
        if (std::min(a,b) < MIN_V)  MIN_V = std::min(a,b);
        edge e(a,b);
        raw_edge.push_back(e);
    }
    printf("Finished Reading Edge File!\n");
    printf("Number of Edges: %d \n",raw_edge.size());
    printf("(Recheck) Vertex with smallest index: %lu; Largest index: %lu \n", MIN_V, MAX_V);
    if (MIN_V != 0){    printf("Index does not start with 0. Recheck! \n"); return false;}
    else    return true;
}


///
/// IT IS ASSUMED THAT COMMUNITY FILE IS NOT NECCESSARILY SORTED I.E. 0,1,2,3....
/// COMMUNITY FILE FORMAT IS SIMILAR TO EDGES i.e. 0 \t 1 meaning vertex 0 belongs in community 1, and so on
/// \brief readCommunityFile
/// \param community: empty community to write to
/// \param filepath
///
uint readCommunityFile(std::vector<std::vector<uint> > &community, const char * filepath)
{
    printf("Parsing Community File: %s ... \n", filepath);
    //an exception is added here for LFR indices, since LFR index starts from 0
    std::set<std::string> exception;
    exception.insert("community.dat");

    //do reading
    std::ifstream infile(filepath);
    if (!infile.good())
    {
        printf("File Erorr! Terminating ... \n");
        return 0;
    }
    bool except = (exception.find(std::string(filepath)) != exception.end());
    std::string line;
    uint MAX_C = 0, MAX_V = 0;
    std::vector<edge> membership;
    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        uint a, b;
        if (!(iss >> a >> b))
        {
            printf("Error While Reading Input! Terminating ... \n");
            return 0;
        } // error
        // process pair (a,b)
        if (except){ a-=1; b-=1;}//reduce index by 1 in cases of i.e. LFR
        if (b > MAX_C)  MAX_C = b;
        if (a > MAX_V)  MAX_V = a;
        edge e(a,b);
        membership.push_back(e);
    }

    for (uint i = 0; i <= MAX_C; i++)
    {
        std::vector<uint> empty;
        community.push_back(empty);
    }

    for (uint i = 0; i < membership.size(); i++)
    {
        edge e = membership.at(i);
        uint v = e.first, c = e.second;
        assert((c <= MAX_C) && "While Building Community For Quality: Index Out Of Bounds!");
        community[c].push_back(v);
    }

    for(auto & comm: community)
        std::sort(comm.begin(), comm.end()); //sort for clustering matching

    return MAX_V;
}

///
/// QUALITY i.e. MODULARITY OR PAIRWISE INDICES MATCHING
/// \brief parseResultForQuality
///
void CalculateClusterComparison(const char * truthfile, const char * clusterfile)
{
    std::vector<std::vector<uint> > truth, result;
    uint MAX_V_t = readCommunityFile(truth, truthfile),
         MAX_V_r = readCommunityFile(result, clusterfile);
    if (MAX_V_t == MAX_V_r && MAX_V_t != 0)
    {
        printf("Calculating Clustering Quality ... \n");
        std::vector<double> quality = Util::ComputePairWiseMatchingIndex(MAX_V_t, truth, result);
        printf("RAND: %f\tJACCARD: %f\tARI: %f\n", quality[0], quality[1], quality[2]);
    }
    else
    {
        printf("There was some problem reading community files!"
               "\n Either The Indices did not match or file is empty!"
               "\n Indices Matching Were Not Calculated ... \n");
    }
}



int main(int argc, char *argv[])
{
    if (argc == 4 && argv[1][1] == 'c') // ./ragg -c filepath outpath
    {
        gen.seed(time(0) + _getpid());
        outPath = argv[3];
        //read file
        std::string filepath(argv[2]);
        std::vector<edge> raw_edge; //edge file
        uint MAX_V = 0; //highest index vertex
        //start loading
        bool load = read_input(&filepath[0], raw_edge, MAX_V);
        if (!load)   {printf("Problem While Reading Input! Terminating ... \n"); return 0;}
        //load success
        Graph G;
        bool build = build_graph(raw_edge, MAX_V, G);
        if (!build) {printf("Problem While Building Graph ! Terminating ... \n"); return 0;}
        //build sucess
        raw_edge.clear(); //in case the graph input is large
        std::vector<edge> hierarchy; //result is stored as edge (selected->neighbour)
        //start processing
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); //start timer
        bool success = process(G, hierarchy);
        if (!success)    {printf("Some problem encountered while clustering the graph! Terminating ... \n"); return 0;}
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now(); //end timer
        postprocess(hierarchy,(MAX_V+1));
        std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        printf("\n ** Finished! Elapsed time: %f sec \n", (elapsed));
        printf("Results are written in the " );
        printf(outPath);
        printf("\s file! \n");
    }
    else if (argc == 4 && argv[1][1] == 'q')    // ./ragg -q community1 community2
    {
        std::string truthfile("community.dat"), clusterfile("cluster.txt");
        CalculateClusterComparison(&truthfile[0], &clusterfile[0]);
    }
    else
        printf("\nBad input format! \n"
               "For clustering try: ./ragg -c filename.txt\n"
               "For calculating quality try: ./ragg -q community1 community2\n");
    return 0;

}
