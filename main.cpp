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
#include <queue>

#include "util.h"
#include "graphgen.h"
#include "graphutil.h"
#include "graphgeo.h"

#include <chrono>
#include <process.h>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/clustering_coefficient.hpp>
#include <boost/graph/erdos_renyi_generator.hpp>
#include <boost/random/linear_congruential.hpp>

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

///
/// \brief writeStringToFile: writting a string to file
/// \param str: first char in string
/// \param outPath: outpath
///
void writeStringToFile(const char * str, const char * outPath)
{
    std::ofstream myfile;
    if (outPath == NULL)
    {
        printf("Output Path is not defined!\n");
        return;
    }
    myfile.open(outPath);
    myfile << str;
    myfile.close();
}

/// SAVE RESULT AS: VertexID\tCommunity (seperated by tab)
/// \brief saveResult
/// \param component
///
void saveResult(const std::vector<uint> &component)
{
    std::ofstream myfile;
    if (outPath == NULL)
    {
        printf("Output Path is not defined!\n");
        return;
    }
    myfile.open (outPath);
    myfile << "Id\tCluster\n";
    for (size_t i = 0; i < component.size(); i++)
    {
        myfile << i << '\t' << component[i] << '\n';
    }
    myfile.close();
}

/// SAVE EDGE/HIERARCHY RESULT AS: VertexID_FROM\VERTEXID_TO (seperated by tab)
/// \brief saveResult
/// \param component
///
void saveHierResult(const std::vector<edge> &hier,
                    const char * hierPath,
                    const char &delim)
{
    std::ofstream myfile;
    if (hierPath == NULL)
    {
        printf("Output Path is not defined! \n");
        return;
    }
    myfile.open (hierPath);
    for (size_t i = 0; i < hier.size(); i++)
    {
        myfile << hier[i].first << delim << hier[i].second << '\n';
    }
    myfile.close();
}

///
/// \brief postprocess
/// \param G: original graph G
/// \param hierarchy: points_to edge
/// \param communities: empty holder for result communities
/// \return: number of components
///
int postprocess(const std::vector<edge> &hierarchy,
                std::vector<std::vector<uint> > &communities,
                const unsigned int &V)
{
    typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
    BoostGraph BG;
    for(size_t i = 0; i <= V; i++)   boost::add_vertex(BG);

    for(auto &e : hierarchy)
    {
        boost::add_edge(e.first, e.second, BG);
    }
    //result if the connected component
    printf("SubGraph H V-%d, E-%d", boost::num_vertices(BG), boost::num_edges(BG));
    std::vector<uint> component(boost::num_vertices(BG));
    //running connected component
    uint num = boost::connected_components(BG, &component[0]);
    GraphUtil::sortVerticesIntoCommunities(num, component, communities);
    printf("\nNumber of Community: %lu \n", num);
    saveResult(component);
    saveHierResult(hierarchy,
                   "C:/Users/Dumex/Desktop/Social_Final/karate/frags.txt",
                    '\t');
    return num;
}


/*-------------------------------------------------------------------------------------------------/
 *
 *                                          K-MXT Functions
 *
 *-------------------------------------------------------------------------------------------------*/
///
/// \brief kMT: in this version, the min weight is set as a function of a vertex
/// i.e. minWeight(v) = d(v)/a where a > 0
/// \param G: input Graph
/// \param points_to: the oriented edge
/// \param isolated: vertices which have no oriented edge
/// \param k: number of oriented edge to include
/// \param frac: fraction of degree of every vertex
/// \return
///
bool kMT_f(const Graph &G,
         std::vector<edge> &points_to,
         std::vector<uint> &isolated,
         const int &k,
         const unsigned int &frac)
{
    if (frac == 0){printf("Fraction of Vertex's Degree Neeeds to be > 0! Terminating"); return 0;} //if frac == 0
    printf("\n Now Doing k-MT(f) ... \n");
    //declaring comparator for min heap
    auto comparator =
            [](const std::pair<uint,uint> &p1, const std::pair<uint,uint> &p2)
            { return p1.first > p2.first; }; //lambda function, testing
    //starting loop
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); //start timer
    for (size_t i = 0; i < G.size(); i++)
    {
        //init
        adj_list v_adj = G[i]; // v adjacency list
        uint minWeight = ceil( (float) v_adj.size()/frac);
        //k-min-heap

        std::priority_queue<std::pair<uint,uint>,
                            std::vector<std::pair<uint,uint> >,
                            decltype(comparator)> q(comparator);
        //starting for vertex
        for (size_t j = 0; j < v_adj.size(); j++)
        {
            //get neighbour vertex u adjacency
            adj_list u_adj = G[v_adj[j]];
            unsigned int comm = Util::NumberOfCommonElement(v_adj, u_adj); //number of common
            if (comm < minWeight)   continue;   //if the edge weight is less than minWeight, skip
            std::pair<uint,uint> p(comm, v_adj[j]); //storing as pair: <comm,index>
            //get min of the heap
            if (q.size() < k){q.push(p);continue;}
            if (comm > q.top().first) //if number of comms is greater than top of min heap
            {
                    q.pop();
                    q.push(p);
            }
            else if (comm == q.top().first) //if equal, take 1 at ran
            {
                std::uniform_real_distribution<double> dist(0, 1);
                if(dist(gen) <= 0.5)
                {
                    q.pop();
                    q.push(p);
                }
            }
        }
        if (q.empty())  isolated.push_back(i); //if q is empty then V[i] is isolated
        // points to the highest common neighbours
        while(!q.empty())
        {
            std::pair<uint,uint> p = q.top();
            points_to.push_back(edge(i,p.second)); //return the index of chosen vertex
            q.pop();
        }
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now(); //end timer
    std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    printf("\n ** K-MT FINISHED! Elapsed time: %f sec \n", (elapsed));
    //assert((G.size()*k == points_to.size()) && "FATAL ERROR WHILE FINISHING kMT: Number of edges does not match!\n");
    //DONE
    return true;
}

//-------------------------------------//
///
/// \brief kMT: do kMXT
/// \param G: input graph g
/// \param points_to: components
/// \param k: number of neighbours to select
/// \param minWeight: edge minimum weight to be selected
/// \return
///
bool kMT(const Graph &G,
         std::vector<edge> &points_to,
         std::vector<uint> &isolated,
         const int &k,
         const int &minWeight)
{
    printf("\n Now Doing k-MT(w) ... \n");
    //declaring comparator for min heap
    auto comparator =
            [](const std::pair<uint,uint> &p1, const std::pair<uint,uint> &p2)
            { return p1.first > p2.first; }; //lambda function, testing
    //starting loop
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); //start timer
    for (size_t i = 0; i < G.size(); i++)
    {
        //init
        adj_list v_adj = G[i]; // v adjacency list
        //k-min-heap

        std::priority_queue<std::pair<uint,uint>,
                            std::vector<std::pair<uint,uint> >,
                            decltype(comparator)> q(comparator);
        //starting for vertex
        for (size_t j = 0; j < v_adj.size(); j++)
        {
            //get neighbour vertex u adjacency
            adj_list u_adj = G[v_adj[j]];
            unsigned int comm = Util::NumberOfCommonElement(v_adj, u_adj); //number of common
            if (comm < minWeight)   continue;   //if the edge weight is less than minWeight, skip
            std::pair<uint,uint> p(comm, v_adj[j]); //storing as pair: <comm,index>
            //get min of the heap
            if (q.size() < k){q.push(p);continue;}
            if (comm > q.top().first) //if number of comms is greater than top of min heap
            {
                    q.pop();
                    q.push(p);
            }
            else if (comm == q.top().first) //if equal, take 1 at ran
            {
                std::uniform_real_distribution<double> dist(0, 1);
                if(dist(gen) <= 0.5)
                {
                    q.pop();
                    q.push(p);
                }
            }
        }
        if (q.empty())  isolated.push_back(i); //if q is empty then V[i] is isolated
        // points to the highest common neighbours
        while(!q.empty())
        {
            std::pair<uint,uint> p = q.top();
            points_to.push_back(edge(i,p.second)); //return the index of chosen vertex
            q.pop();
        }
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now(); //end timer
    std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    printf("\n ** K-MT FINISHED! Elapsed time: %f sec \n", (elapsed));
    //assert((G.size()*k == points_to.size()) && "FATAL ERROR WHILE FINISHING kMT: Number of edges does not match!\n");
    //DONE
    return true;
}


///
/// \brief kMT: k-Max-Triangulated
/// \param G: graph G
/// \param k: number of k to point to
/// \param points_to: empty hierarchy to write to
/// \return: true or false for success/fail execution
///
bool kMT(const Graph &G, std::vector<edge> &points_to, const int &k)
{
    printf("\n Now Doing k-MT ... \n");
    //declaring comparator for min heap
    auto comparator =
            [](const std::pair<uint,uint> &p1, const std::pair<uint,uint> &p2)
            { return p1.first > p2.first; }; //lambda function, testing
    //starting loop
    std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); //start timer
    for (size_t i = 0; i < G.size(); i++)
    {
        //init
        adj_list v_adj = G[i]; // v adjacency list
        //k-min-heap

        std::priority_queue<std::pair<uint,uint>,
                            std::vector<std::pair<uint,uint> >,
                            decltype(comparator)> q(comparator);
        //starting for vertex
        for (size_t j = 0; j < v_adj.size(); j++)
        {
            //get neighbour vertex u adjacency
            adj_list u_adj = G[v_adj[j]];
            unsigned int comm = Util::NumberOfCommonElement(v_adj, u_adj); //number of common
            std::pair<uint,uint> p(comm, v_adj[j]); //storing as pair: <comm,index>
            //get min of the heap
            if (q.size() < k){q.push(p);continue;}
            if (comm > q.top().first) //if number of comms is greater than top of min heap
            {
                if (q.size() < k)   q.push(p); //if size is < k
                else if (q.size() == k) //if it is full, pop
                {
                    q.pop();
                    q.push(p);
                }
            }
            else if (comm == q.top().first) //if equal, take 1 at ran
            {
                std::uniform_real_distribution<double> dist(0, 1);
                if(dist(gen) <= 0.5)
                {
                    q.pop();
                    q.push(p);
                }
            }
        }
        // points to the highest common neighbours

        while(!q.empty())
        {
            std::pair<uint,uint> p = q.top();
            points_to.push_back(edge(i,p.second)); //return the index of chosen vertex
            q.pop();
        }
    }
    std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now(); //end timer
    std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
    printf("\n ** K-MT FINISHED! Elapsed time: %f sec \n", (elapsed));
    //assert((G.size()*k == points_to.size()) && "FATAL ERROR WHILE FINISHING kMT: Number of edges does not match!\n");
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
        adj_list v = myGraph[from], u = myGraph[to];
        if (std::find(v.begin(), v.end(),to)==v.end())
            myGraph[from].push_back(to); //get adj_list of source vertex
        if (std::find(u.begin(), u.end(),from)==u.end())
            myGraph[to].push_back(from); // not adding double edge
        //check for duplicate neighbour

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
/// \brief edgeComparator: custom edge comparator
/// \param a
/// \param b
/// \return
///
static bool edgeComparator(const edge & a, const edge & b)
{
    return (std::minmax(a.first, a.second) <
                std::minmax(b.first, b.second));
}

///
/// \brief read_graph: Read in edge.txt file
/// it is assumed that each line contains source'\t'target (seperated by tab)
/// it is assumed that index starts from 0
///
bool read_input(const char * filepath, std::vector<edge> &raw_edge, uint &MAX_V)
{
    std::string exceptionFile("network.dat");
    //do reading
    std::string path(filepath);
    bool except = (path.find(exceptionFile)!=std::string::npos);
    uint MIN_V = 0;
    std::ifstream infile(filepath);
    std::set<edge, bool(*)(const edge &a, const edge &b)> s(&edgeComparator);
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
        if(except){ a-=1;   b-=1;}
        // process pair (a,b)
        if (std::max(a,b) > MAX_V)  MAX_V = std::max(a,b);
        if (std::min(a,b) < MIN_V)  MIN_V = std::min(a,b);
        s.insert(edge(a,b));
    }
    raw_edge.assign(s.begin(),s.end());
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
    std::string exceptionFile("community.dat");
    //do reading
    std::ifstream infile(filepath);
    if (!infile.good())
    {
        printf("File Erorr! Terminating ... \n");
        return 0;
    }
    //bool except = (exception.find(std::string(filepath)) != exception.end());
    std::string path(filepath);
    bool except = (path.find(exceptionFile)!=std::string::npos);
    std::string line;
    uint MAX_C = 0, MAX_V = 0;
    std::vector<edge> membership;
    std::getline(infile,line); //skip the first line
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
        if (b > MAX_C && b != INT_MAX)  MAX_C = b;
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
        uint N = MAX_V_t;
        //uint N = MAX_V_t + 1; //index starts from 0
        std::vector<double> quality = Util::ComputePairWiseMatchingIndex(N, truth, result);
        printf("RAND: %f\tJACCARD: %f\tARI: %f\n", quality[0], quality[1], quality[2]);
    }
    else
    {
        printf("There was some problem reading community files!"
               "\n Either The Indices did not match or file is empty!"
               "\n Indices Matching Were Not Calculated ... \n");
    }
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//----------------------------------- HIDDEN BLOCK EXP -------------------------------------------
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
///
/// \brief singleRunHiddenBlockModel: run hidden block model a single time
/// \param l: number of vertex per block
/// \param no_runs: number of execution (each time a different Gnp is generated)
/// \param k: number of max-triangulated neighbours to be selected
/// \param q: within edge probability
/// \param p: global edge probability
/// \return: a string contains all stats
///
std::string singleRunHiddenBlockModel(const int& l, const int &no_runs, const int&k,
                                      const double &q, const double &p)
{
    //init
    GraphGen graphGen;
    uint n = l*4;
    std::string record;
    //init the truth community
    std::vector<std::vector<uint> > truth_community = {std::vector<uint>(),
                                                       std::vector<uint>(),
                                                       std::vector<uint>(),
                                                       std::vector<uint>()};
    for(int i = 0; i < n; i++)
    {
        int v_community = i/l;
        truth_community[v_community].push_back(i);
    }
    //statistics/indicices to keep track of
    double sum_frac = 0;
    unsigned int num_comps = 0;
    std::vector<double> quality_indices = {0.0,0.0,0.0};
    //starting main loop
    for (int i = 0; i < no_runs; i++)
    {
        Graph G;
        std::vector<edge> hierarchy; //hierarchy tree
        std::vector<std::vector<uint> > communities;
        graphGen.generateGnp_BlockModel(G, l, q, p); //generate Block Model
        bool success = kMT(G, hierarchy, k); // else call kMT
        if (!success)    {printf("Some problem encountered while clustering the graph! Terminating ... \n"); return "?!";}
        //save the graph
        if (i == 0)
        {
            std::stringstream outpathss;
            outpathss << "C:/Users/Dumex/Desktop/HiddenBlockModel/GraphFile/" << n << "_" << q << "_" << k
                      << "_hier.graphml" ;
            GraphUtil::saveOutputAsGraphML(hierarchy, l, &outpathss.str()[0]);
        }

        //get the statistics
        num_comps += postprocess(hierarchy,communities,l*4);    //number of components
        sum_frac += Util::fractionOfIntraEdge(hierarchy, l, k); //fraction of correct edge
        std::vector<double> quality = Util::ComputePairWiseMatchingIndex(n, truth_community, communities); //RAND, JACCARD, ARI
        for (int i = 0; i < 3; i++) quality_indices[i] += quality[i];
        //done
    }
    //record
    std::stringstream ss;
    ss << n << '\t' << q << '\t' << p << '\t' << k << '\t' << (num_comps/no_runs) << '\t' << (sum_frac/no_runs)
       << '\t';
    for (int i = 0; i < 3; i++) ss << (quality_indices[i]/no_runs) << '\t';
    ss << '\n';
    record += ss.str();
    return record;
}

///
/// \brief autoRunHiddenBlockModel_FixedDegree: auto execution of hidden block model
/// \param min_l: starting number of vertex per block
/// \param max_l: stoping number of vertex per block
/// \param no_runs: number of runs
/// \param k: number of max-triangulated-neighbours
///
void autoRunHiddenBlockModel_FixedDegree_FixedPrThreshold(const int& min_l, const int& max_l,
                                         const int &no_runs, const int &k)
{
    //init
    GraphGen graphGen;
    //stats
    /*  n\tq\tp\tNumberOfComponents\tfractionofWrongEdge
     */
    std::string record = "*AutoRunHiddenBlockModel: Fixed Prob Threshold";
    record += "Number of k-Max-NeighbourL: k= "+ k + '\n';
    record += "\nn\tq\tp\t#Comp\tfractionofWrongEdge\tRAND\tJACCARD\tARI\n";
    //start
    for (int l = min_l; l <= max_l; l+= 500)
    {
        //graph atts
        uint n = l*4;
        double p = sqrt(n)/(n); // p is kept at the threshold
        double min_q = sqrt(l)/l,
               max_q = min_q*4,
               step_q = (max_q-min_q)/5; //playing around q
        //generate truth community
        std::vector<std::vector<uint> > truth_community = {std::vector<uint>(),
                                                           std::vector<uint>(),
                                                           std::vector<uint>(),
                                                           std::vector<uint>()};
        for(int i = 0; i < n; i++)
        {
            int v_community = i/l;
            truth_community[v_community].push_back(i);
        }
        //doing work
        for (double q = min_q; q <= max_q; q += step_q)
        {
            //run atts
            double sum_frac = 0;
            unsigned int num_comps = 0;
            std::vector<double> quality_indices = {0.0,0.0,0.0};
            for (int i = 0; i < no_runs; i++)
            {
                Graph G;
                std::vector<edge> hierarchy; //hierarchy tree
                std::vector<std::vector<uint> > communities;
                graphGen.generateGnp_BlockModel(G, l, q, p); //generate Block Model
                bool success = kMT(G, hierarchy, k); // else call kMT
                if (!success)    {printf("Some problem encountered while clustering the graph! Terminating ... \n"); return;}

                //get the statistics
                num_comps += postprocess(hierarchy,communities,l*4);
                sum_frac += Util::fractionOfIntraEdge(hierarchy, l, k);
                std::vector<double> quality = Util::ComputePairWiseMatchingIndex(n, truth_community, communities);
                for (int i = 0; i < 3; i++) quality_indices[i] += quality[i];
                //done
            }
            //record
            std::stringstream ss;
            ss << n << '\t' << q << '\t' << p << '\t' << (num_comps/no_runs) << '\t' << (sum_frac/no_runs)
               << '\t';
            for (int i = 0; i < 3; i++) ss << (quality_indices[i]/no_runs) << '\t';
            ss << '\n';
            record += ss.str();
        }
    }
    std::string outPath = "C:/Users/Dumex/Desktop/stats.txt";
    writeStringToFile(&record[0], &outPath[0]);
}


///
/// \brief autoRunHiddenBlockModel: caller of the singlerun method,
/// use this to playaround with parameters q,p,k
///
void autoRunHiddenBlockModel()
{
    //stats
    /*  n\tq\tp\tNumberOfComponents\tfractionofWrongEdge
     */
    std::string outPath = "C:/Users/Dumex/Desktop/HiddenBlockModel/exp_log.txt"; //modify this path
    std::string record = "*AutoRunHiddenBlockModel: Within Edge = 2 Inter-edge";
    record += "\nn\tq\tp\tk\t#Comp\tfractionofWrongEdge\tRAND\tJACCARD\tARI\n";
    //params
    int l_min = 5000,
        l_max = 5000,
        l_step = 500,
        no_runs = 5;
    for (int l = l_min; l <= l_max; l += l_step)
    {
        int n = l*4;
        double c = 0.45;
        double p_threshold = 1/(pow(n,c)), // p is kept at the threshold
               max_q = 10*p_threshold, //twice as many inter-cluster edges
               q_stepsize = p_threshold;
               //q_stepsize = (q_threshold*4)/10;
        for (double q = 0.0; q < max_q ; q+=q_stepsize)
        {
            for(int k = 1; k <= 5; k++)
            {
                record += singleRunHiddenBlockModel(l, no_runs, k , q , p_threshold);
            }
        }
    }

    writeStringToFile(&record[0], &outPath[0]);
}

///
/// \brief HiddenBlockModel_ClusteringCoefficient_EXP: watching the clustering coefficient as q grows
///
void HiddenBlockModel_ClusteringCoefficient_EXP(const uint &l)
{
    uint n = l*4;
    double p = sqrt(n)/(n), // p is kept at the threshold
           max_q = 2.6*2.6*2*p, //twice as many inter-cluster edges
           q_stepsize = p;
           //q_stepsize = (q_threshold*4)/10;
    std::stringstream outputss;
    outputss << "**** Fixed N = "<< n << "\tp = " << p << "\n";
    outputss << "q\tcc\tintra_cc\tinter_cc\n";
    std::vector<std::vector<uint> > truth_community = {std::vector<uint>(),
                                                       std::vector<uint>(),
                                                       std::vector<uint>(),
                                                       std::vector<uint>()};
    for(int i = 0; i < n; i++)
    {
        int v_community = i/l;
        truth_community[v_community].push_back(i);
    }

    std::vector<std::vector<double> > ratioMatrix;
    for (double q = 0.0; q < max_q ; q+=q_stepsize)
    {
        GraphGen graphGen;
        //init the truth community

        //starting main loop
        Graph G;
        graphGen.generateGnp_BlockModel(G, l, q, p); //generate Block Model

        //get the statistics
        std::vector<double> ratio;
        for (int i = 0; i < G.size(); i++)  ratio.push_back(GraphUtil::vertexHiddenBlockWeightRatio(G, i, l));
        ratioMatrix.push_back(ratio);
        /*
        std::vector<double> cc, intra_cc, inter_cc;
        double mean_cc = 0.0 , mean_intra_cc = 0.0, mean_inter_cc = 0.0;
        mean_cc = GraphUtil::clusteringCoefficient(G,cc);
        mean_intra_cc = GraphUtil::intraClusteringCoefficient(G, truth_community, intra_cc);
        mean_inter_cc = GraphUtil::interClusteringCoeffcient(G, truth_community, inter_cc);
        outputss << q << "\t" << mean_cc << "\t" << mean_intra_cc << "\t" << mean_inter_cc << "\n";
        //done*/
    }
    std::string out = "C:/Users/Dumex/Desktop/HiddenBlockModel/cc_log.txt";
    writeStringToFile(&outputss.str()[0], &out[0]);
    Util::writeMatrixToFile(ratioMatrix,"C:/Users/Dumex/Desktop/HiddenBlockModel/ratio_log.txt", ';');
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//----------------------------------- HIDDEN BLOCK EXP - END--------------------------------------
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//************************************************************************************************
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//----------------------------------- GEOGRAPHICAL GRAPH EXP -------------------------------------
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
///
/// \brief saveConvexHulls: output the convex hulls of the geometric graphs
/// the convex is seperated into two files:
/// 1. vertex   lat;long;id
///             xx;yy;i
/// 2. edge     id_from,id_to
/// \param hulls: set of hulls,
/// for each hull, the points in the exterior ring is drawn progressively i.e. [0] -> [1] -> [2] -> [0]
/// \param outDir: output directory
///
void saveConvexHullsPolygon(const std::vector<std::vector<std::pair<double,double> > > &hulls,
                           const char * outDir)
{
    //init
    std::string v_path = outDir + std::string("poly.txt");
    std::ofstream v_file;
    v_file.open(v_path);
    v_file << "lat;long;id;com\n";
    uint num_v = 0;
    //starting
    for (int i = 0; i < hulls.size(); i++)
    {
        std::vector<std::pair<double, double> > convex_hull = hulls.at(i);
        if (convex_hull.empty())    continue; //how can one be empty?
        for (int j = 0; j < convex_hull.size(); j++)
        {
            std::pair<double,double> point;
            point = convex_hull.at(j);
            v_file << std::setprecision(8) << point.first << ";"
                                           << point.second << ";"
                                           << num_v  << ";"
                                           << (i+1) << "\n";
            num_v++;
        }
    }
}


///
/// \brief saveConvexHulls: output the convex hulls of the geometric graphs
/// the convex is seperated into two files:
/// 1. vertex   lat;long;id
///             xx;yy;i
/// 2. edge     id_from,id_to
/// \param hulls: set of hulls,
/// for each hull, the points in the exterior ring is drawn progressively i.e. [0] -> [1] -> [2] -> [0]
/// \param outDir: output directory
///
void saveConvexHulls(const std::vector<std::vector<std::pair<double,double> > > &hulls,
                     const char * outDir)
{
    //init
    std::string v_path = outDir + std::string("hull_vertex.txt"),
                e_path = outDir + std::string("hull_edge.txt");
    std::ofstream v_file, e_file;
    v_file.open(v_path);
    v_file << "lat;long;id\n";
    e_file.open(e_path);
    uint num_v = 0;
    //starting
    for (int i = 0; i < hulls.size(); i++)
    {
        std::vector<std::pair<double, double> > convex_hull = hulls.at(i);
        uint starting_v = num_v; //record the starting vertex
        for (int j = 0; j < convex_hull.size(); j++)
        {
            std::pair<double,double> point = convex_hull.at(j);
            v_file << std::setprecision(8) << point.first << ";" << point.second << ";" << num_v <<'\n';
            if (j == convex_hull.size() -1) e_file << num_v << "," << starting_v << '\n';
            else    e_file << num_v << "," << (num_v+1) << '\n';
            num_v++;
        }
    }
}

///
/// \brief geoGraphExp: experiment with geo graph
/// \param filePath: inpath
///
void geoGraphExp(const char * filePath, const double &limitingDistance)
{
    //load file
    GraphGeo geoG;
    bool success = geoG.readSpatialData(filePath);
    if (!success){printf("Loading File: Failed!\n"); return;}
    //init
    Graph G;
    success = geoG.generateGeoDiscGraph(G, limitingDistance);
    if (!success){printf("GeoGraph Generating: Failed!\n"); return;}
    //save the graph
    //GraphUtil::saveGraphAsUniqueEdgeList(G, "C:/Users/Dumex/Desktop/Small_Set/KMXT/Tralfagar/geograph.txt", ',');
    //do work
    //
    for (int k = 1; k <= 4; k++)
    {
        std::vector<edge> hier;
        std::vector<uint> isolated;
        //kMT(G,hier,k);
        kMT(G,hier,isolated,k,80);
        std::stringstream ss;
        /* Saving the current kMXT
         *
        ss << "C:/Users/Dumex/Desktop/GeoEXp/testout_" << k << ".txt";
        saveHierResult(hier, &ss.str()[0], ',');    //cluster
        /* END
         * Now Calculating Connected Components
         */
        std::vector<std::vector<uint> > communities;
        uint num = postprocess(hier, communities, G.size());
        /* END
         * Saving Size Distribution
         */
        double average_size = 0.0;
        std::ostringstream sizedis, sizedis_file, ashapeinput;
        //FILE HEADERS
        sizedis << "comm;size\n";
        sizedis_file << "C:/Users/Dumex/Desktop/Small_Set/KMXT/Tralfagar/clusterSize_d" << limitingDistance <<"_k"  << k << ".txt";
        // CALCULATE AVERAGE SIZE
        for (int i = 0; i < communities.size(); i++)
        {
            average_size += (double) communities[i].size()/num;
        }

        std::vector<std::vector<std::pair<double,double> > > allpos;
        /* END
         * Calculating Convex Hulls and Alphashape
         */
        ashapeinput.precision(10);
        ashapeinput << "lat;long;id;com\n";
        uint no_unique_coor = 0,
             alphaShapeUniqueCount = 1,
             uniqueClusterCount = 1;
        for (int i = 0; i < num; i++)
        {
            std::vector<uint> comm = communities.at(i);
            //if (comm.size() < 100) continue; // skip if size is less than average
            if (comm.size() < average_size) continue;
            // SAVE CLUSTER SIZE
            sizedis << uniqueClusterCount << ";" << communities[i].size() << '\n';
            uniqueClusterCount++;
            // SAVE CONVEX HULL
            std::vector<std::pair<double,double> > pos,
                    hull;
            geoG.getPos(comm, pos); //recovering original positions
            geoG.convexHull(pos, hull); //calculating convex hulls
            allpos.push_back(hull); //store
            //getting unique position for ashape
           // std::set<std::pair<double,double> > uniqueCoor(pos.begin(),pos.end());
            std::vector<std::pair<double,double> > uniqueCoor = pos;
            //if (uniqueCoor.size() <= 3) continue;
            for (auto it = uniqueCoor.begin();
                 it != uniqueCoor.end(); it++) {
                 std::pair<double,double> p = *it;
                 ashapeinput << p.first << ";"
                             << p.second << ";"
                             << no_unique_coor << ";"
                             << alphaShapeUniqueCount << '\n';
                 no_unique_coor++;
            }
            alphaShapeUniqueCount++;
        }
        //WRITE SIZE DISTRIBUTION TO FILE
        writeStringToFile(&sizedis.str()[0], &sizedis_file.str()[0]);
        //WRITE CONVEX HULL TO FILE
        std::stringstream oo;
        oo << "C:/Users/Dumex/Desktop/Small_Set/KMXT/Tralfagar/ConvexHull_" << k << "/";
        geoG.removeWithinPolygons(allpos); //uncomment this line to remove all polygons reside within others
        saveConvexHulls(allpos, &oo.str()[0]); //CONVEX HULL FOR DRAWING ON QGIS
        saveConvexHullsPolygon(allpos, &oo.str()[0]); //CONVEX HULL FOR DENSITY CALCULATION
        oo << "ashape_in.txt";
        writeStringToFile(&ashapeinput.str()[0],&oo.str()[0]);
    }
}

void geoGraphExp(const Graph &G)
{
    const float limitingDistance = 0.0001;
    GraphGeo geoG;
    bool success = geoG.readSpatialData("C:/Users/Dumex/Desktop/GeoEXp/innerLondon.txt");
    if (!success){printf("Loading File: Failed!\n"); return;}

    for (int k = 1; k <= 4; k++)
    {
        std::vector<edge> hier;
        kMT(G,hier,k);
        //kMT(G,hier,k,0);
        std::stringstream ss;
        /* Saving the current kMXT
         *
        ss << "C:/Users/Dumex/Desktop/GeoEXp/testout_" << k << ".txt";
        saveHierResult(hier, &ss.str()[0], ',');    //cluster
        /* END
         * Now Calculating Connected Components
         */
        std::vector<std::vector<uint> > communities;
        uint num = postprocess(hier, communities, G.size());
        /* END
         * Saving Size Distribution
         */
        double average_size = 0.0;
        std::ostringstream sizedis, sizedis_file, ashapeinput;
        //FILE HEADERS
        sizedis << "comm;size\n";
        sizedis_file << "C:/Users/Dumex/Desktop/GeoEXp/clusterSize_d" << limitingDistance <<"_k"  << k << ".txt";
        // CALCULATE AVERAGE SIZE
        for (int i = 0; i < communities.size(); i++)
        {
            average_size += (double) communities[i].size()/num;
        }

        std::vector<std::vector<std::pair<double,double> > > allpos;
        /* END
         * Calculating Convex Hulls and Alphashape
         */
        ashapeinput.precision(10);
        ashapeinput << "lat;long;id;com\n";
        uint no_unique_coor = 0,
             alphaShapeUniqueCount = 1,
             uniqueClusterCount = 1;
        for (int i = 0; i < num; i++)
        {
            std::vector<uint> comm = communities.at(i);
            //if (comm.size() < 100) continue; // skip if size is less than average
            if (comm.size() < average_size) continue;
            // SAVE CLUSTER SIZE
            sizedis << uniqueClusterCount << ";" << communities[i].size() << '\n';
            uniqueClusterCount++;
            // SAVE CONVEX HULL
            std::vector<std::pair<double,double> > pos,
                    hull;
            geoG.getPos(comm, pos); //recovering original positions
            geoG.convexHull(pos, hull); //calculating convex hulls
            allpos.push_back(hull); //store
            //getting unique position for ashape
           // std::set<std::pair<double,double> > uniqueCoor(pos.begin(),pos.end());
            std::vector<std::pair<double,double> > uniqueCoor = pos;
            //if (uniqueCoor.size() <= 3) continue;
            for (auto it = uniqueCoor.begin();
                 it != uniqueCoor.end(); it++) {
                 std::pair<double,double> p = *it;
                 ashapeinput << p.first << ";"
                             << p.second << ";"
                             << no_unique_coor << ";"
                             << alphaShapeUniqueCount << '\n';
                 no_unique_coor++;
            }
            alphaShapeUniqueCount++;
        }
        //WRITE SIZE DISTRIBUTION TO FILE
        writeStringToFile(&sizedis.str()[0], &sizedis_file.str()[0]);
        //WRITE CONVEX HULL TO FILE
        std::stringstream oo;
        oo << "C:/Users/Dumex/Desktop/GeoEXp/ConvexHull_" << k << "/";
        geoG.removeWithinPolygons(allpos); //uncomment this line to remove all polygons reside within others
        saveConvexHulls(allpos, &oo.str()[0]); //CONVEX HULL FOR DRAWING ON QGIS
        saveConvexHullsPolygon(allpos, &oo.str()[0]); //CONVEX HULL FOR DENSITY CALCULATION
        oo << "ashape_in.txt";
        writeStringToFile(&ashapeinput.str()[0],&oo.str()[0]);
    }
}



///
/// \brief inspectingCommunity: inspecting the given graph for clustering coeffcient
/// \param G
///
float community_clustering_coefficient(const Graph &in,
                                      const std::vector<std::vector<uint> > &communities)
{
    printf("\nCalulating clustering coefficient of induced subgraphs...");
    //get the new graph
    Graph G;
    GraphUtil::isolate_communities(in, communities, G);
    //get the unique edge list
    std::vector<edge> E, newE;
    GraphUtil::get_unique_edge_from_graph(G, E);
    GraphUtil::reindex_edge(E, newE);
    //init BoostGraph
    typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
    BoostGraph BG;

    for(auto &e : newE)
    {
        boost::add_edge(e.first, e.second, BG);
    }
    printf("Fragment Subgraph H V-%d, E-%d\n", boost::num_vertices(BG), boost::num_edges(BG));

    // The clustering property, container, and map define the containment
    // and abstract accessor for the clustering coefficients of vertices.
    typedef boost::exterior_vertex_property<BoostGraph, float> ClusteringProperty;
    typedef ClusteringProperty::container_type ClusteringContainer;
    typedef ClusteringProperty::map_type ClusteringMap;

    ClusteringContainer coefs(boost::num_vertices(BG));
    ClusteringMap cm(coefs, BG);
    float cc = boost::all_clustering_coefficients(BG, cm); //ignore this, isolated vertices are considered discarded
    /*
    boost::graph_traits<BoostGraph>::vertex_iterator i, end;
    /*
    for(boost::tie(i, end) = boost::vertices(BG); i != end; ++i) {
            std::cout << boost::get(cm, *i) << '\n';
        }
        */
    printf("!\nDone !!!");
   return cc;
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//----------------------------------- GEOGRAPHICAL GRAPH END -------------------------------------
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------


int main(int argc, char *argv[])
{
    //geo graph but load edge
    //geoGraphExp("C:/Users/Dumex/Desktop/Small_Set/tralfagar.txt", 0.01);
    //return 0;
    /* BLOCK MODEL EXP*/
  // autoRunHiddenBlockModel_FixedDegree_FixedPrThreshold(1000,2000,10,2);
    autoRunHiddenBlockModel();
  // HiddenBlockModel_ClusteringCoefficient_EXP(2500);
    return 0;
    gen.seed(time(0) + _getpid());
    /*
    for (int k = 10; k <= 15; k++)
    {
        uint n = 100000;
        float p = (float) k/n;
        GraphGen gen;
        Graph G;
        gen.generateGnp(G,n,p);
        std::stringstream path;
        path << "C:/Users/Dumex/Desktop/Social_Final/Gnp_speedtest/n_" << n
             << "_p_" << p << ".txt";
        GraphUtil::saveGraphAsUniqueEdgeList(G, &path.str()[0], '\t');
    }*/


   //  Main Terminal
    if (argc == 6 && argv[1][1] == 'c') // ./ragg -c filepath outpath truthpath
    {
        /*
         *  BUILD GRAPH AND Ground truth communities
         */
        gen.seed(time(0) + _getpid());
        outPath = argv[3];
        //read file
        std::string filepath(argv[2]);
        std::vector<edge> rawEdge; //edge file
        uint MAX_V = 0; //highest index vertex
        //start loading
        bool load = read_input(&filepath[0], rawEdge, MAX_V);
        if (!load)   {printf("Problem While Reading Input! Terminating ... \n"); return 0;}

        //load success
        Graph G;
        bool build = build_graph(rawEdge, MAX_V, G);
        if (!build) {printf("Problem While Building Graph ! Terminating ... \n"); return 0;}
        GraphUtil::saveGraphAsUniqueEdgeList(G,"C:/Users/Dumex/Desktop/Social_Final/LFR/smaller/gephi_in.csv", ',');
        //return 0;
        /*
        //get the edge weight distribution
        printf("\n ** Geting Edge W **\n");
        std::vector<uint> edgeW;
        GraphUtil::getEdgeTriangleDistribution(G, rawEdge, edgeW);
        std::ostringstream ws;
        for(size_t i = 0; i < edgeW.size(); i++)    ws << i << "\t" << edgeW[i] << '\n';
        writeStringToFile(&ws.str()[0], "C:/Users/Dumex/Desktop/Social_Final/kmxt_w/lfr_n_10k_k50_maxk_100_size_10_1000/edgeWeight.txt");
        return 0;
        //*/
        float C = GraphUtil::graph_clustering_coefficient(rawEdge);
        std::vector<std::vector<uint> > truthCommunities;
        readCommunityFile(truthCommunities, argv[4]);
        //float groundTruthCC = community_clustering_coefficient(G, truthCommunities);
        //printf("\n ------ CC of induced ground truth communities:\t%f------\n",groundTruthCC);
        /*
         *  DOING WORK
         */
        //build sucess
        rawEdge.clear(); //in case the graph input is large
        std::vector<edge> hierarchy; //result is stored as edge (selected->neighbour)
        std::vector<std::vector<uint> > communities; //community vector
        std::vector<uint> isolated;

        //start processing
        std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); //start timer
        //for(uint minWeight = 0; minWeight <= 220; minWeight+=20)
        //for(uint k = 1; k <= 5; k++)
        //{
        uint k = 2;
        uint minWeight = 3;
        bool success = kMT(G,
                           hierarchy,
                           isolated,
                           k,
                           minWeight);
        if (!success)    {printf("Some problem encountered while clustering the graph! Terminating ... \n"); return 0;}
        std::chrono::high_resolution_clock::time_point end = std::chrono::high_resolution_clock::now(); //end timer
        postprocess(hierarchy,communities,MAX_V);
        std::chrono::duration<double> elapsed = std::chrono::duration_cast<std::chrono::duration<double>>(end - start);
        /*
         * FINISHING
         */
        //finishing
        //Cleaning up: i.e. remove isolated vertices etc.
        std::vector<std::vector<uint> > cleanedCommunities, modifiedCommunities;
        for(int i = 0; i < communities.size(); i++)
        {
            if (communities[i].size() > 1)  cleanedCommunities.push_back(communities[i]);
        }

        GraphUtil::modifyCommunity(truthCommunities,
                                   isolated,
                                   modifiedCommunities);
        //Results
        printf("\n ** Finished! Elapsed time: %f sec \n", (elapsed));
        printf("Results are written in the " );
        printf(outPath);
        printf("\s file! \n");
        //calculate modularity
        uint nonIsolatedVertices = G.size() - isolated.size();
        float Q = GraphUtil::modularity(G, communities),
              MC =  community_clustering_coefficient(G, communities),
              mARI = Util::ComputePairWiseMatchingIndex(nonIsolatedVertices,
                                                        modifiedCommunities,
                                                        cleanedCommunities)[2], //*/,
              averageDiam = GraphUtil::get_average_diameter_of_subgraphs(G, communities);

        FILE * statfile;
        statfile = fopen(argv[5], "a");
        int averageSize = 0;
        if (cleanedCommunities.size() > 0) averageSize = nonIsolatedVertices/(cleanedCommunities.size());
        printf("SIZE: %d", averageSize);
        fprintf(statfile,
               "\n ----------------- RESULTS ----------------"
               "\n ---- Param: k-%d; minWeight-%d -----------"
               "\n ** Modularity Q:\t%f "
               "\n ** Global Clustering Coefficient:\t%f "
               "\n ** Subgraphs clustering coefficient:\t%f"
               "\n ** Modified ARI (subtracting isolated vertices):\t%f"
               "\n ** Average Communitiy size:\t%d"
               "\n ** Average Diamater:\t%f"
               "\n ** No of Communities:\t%d"
               "\n ** No of isolated:\t%d"
               "\n ------------------------------------------",
                k, minWeight,
                Q, C, MC,
                mARI,
                averageSize,
                averageDiam,
                cleanedCommunities.size(),
                isolated.size());
        //glolba clustering coefficient
        //calculate connected component clustering coefficient
        hierarchy.clear(); communities.clear(); isolated.clear(); //remove later
        //}
    }
    else if (argc == 4 && argv[1][1] == 'q')    // ./ragg -q community1 community2
    {
        //std::string truthfile("community.dat"), clusterfile("cluster.txt");
        CalculateClusterComparison(argv[2], argv[3]);
    }
    else
        printf("\nBad input format! \n"
               "For clustering try: ./ragg -c filename.txt\n"
               "For calculating quality try: ./ragg -q community1 community2\n");
    return 0;

}
