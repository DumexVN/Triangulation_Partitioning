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
    for(size_t i = 0; i < V; i++)   boost::add_vertex(BG);

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
                   "C:/Users/Dumex/Qt_WorkingSpace/build-IIIa_Optimise-Desktop_Qt_5_7_0_MinGW_32bit-Debug/hier.txt",
                    ';');
    return num;
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
    //assert((G.size()*k == points_to.size()) && "FATAL ERROR WHILE FINISHING kMT: Number of edges does not match!\n");
    //DONE
    return true;
}

///
/// \brief process: do aggregation III.a
/// \param G
///
/*
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
        // points to the highest common neighbours
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
*/

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
        myGraph[to].push_back(from); // hmm, some files have only single edge
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
    int l_min = 500,
        l_max = 3500,
        l_step = 500,
        no_runs = 5;
    for (int l = l_min; l <= l_max; l += l_step)
    {
        int n = l*4;
        double p_threshold = sqrt(n)/(n), // p is kept at the threshold
               max_q = 2.6*2.6*2*p_threshold, //twice as many inter-cluster edges
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
            v_file << point.first << ";" << point.second << ";" << num_v <<'\n';
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
    //GraphUtil::saveGraphAsUniqueEdgeList(G, "C:/Users/Dumex/Desktop/GeoEXp/geograph.txt", ',');
    //do work
    std::vector<std::vector<double> > CC; //container for cluster coeffcient
    std::vector<double> cc; //input graph clustering coefficient
    GraphUtil::clusteringCoefficient(G, cc);
    CC.push_back(cc);
    cc.clear();
    //
    for (int k = 1; k <= 5; k++)
    {
        std::vector<edge> hier;
        kMT(G,hier,k);
        std::stringstream ss;
        ss << "C:/Users/Dumex/Desktop/GeoEXp/testout_" << k << ".txt";
        saveHierResult(hier, &ss.str()[0], ',');    //cluster
        //producing convex hulls
        std::vector<std::vector<uint> > communities;
        uint num = postprocess(hier, communities, G.size());
        //saving size distribution
        double average_size = 0.0;
        std::ostringstream sizedis, sizedis_file, ashapeinput;
        sizedis_file << "C:/Users/Dumex/Desktop/GeoEXp/sizedis_d" << limitingDistance <<"_k"  << k << ".txt";
        for (int i = 0; i < communities.size(); i++)
        {
            sizedis << communities[i].size() << '\n';
            average_size += (double) communities[i].size()/num;
        }
        writeStringToFile(&sizedis.str()[0], &sizedis_file.str()[0]);
        //getting clustering coefficient
        /*
        std::vector<double> subcc(G.size(), 0.0);
        for (int i = 0; i < communities.size(); i++)
        {
            std::vector<uint> comm = communities[i];
            std::vector<double> hcc;
            Graph H;
            geoG.generateGeoDiscSubGraph(H, comm, limitingDistance);
            GraphUtil::clusteringCoefficient(H, hcc);
            for(int j = 0; j < comm.size(); j++)
            {
                subcc[comm[j]] = hcc[j];
            }
        }
        CC.push_back(subcc);
                */
        //
        std::vector<std::vector<std::pair<double,double> > > allpos;
        //prequisite for alphashape
        ashapeinput.precision(10);
        ashapeinput << "lat;long;id;com\n";
        uint no_unique_coor = 0,
                no_comm = 1;
        for (int i = 0; i < num; i++)
        {
            std::vector<uint> comm = communities.at(i);
            if (comm.size() < 3 || comm.size() < average_size) continue; // skip if size is less than average
            std::vector<std::pair<double,double> > pos,
                    hull;
            geoG.getPos(comm, pos); //recovering original positions
            geoG.convexHull(pos, hull); //calculating convex hulls
            allpos.push_back(hull); //store
            //getting unique position for ashape
            std::set<std::pair<double,double> > uniqueCoor(pos.begin(),
                                                           pos.end());
            if (uniqueCoor.size() <= 3) continue;
            for (auto it = uniqueCoor.begin();
                 it != uniqueCoor.end(); it++) {
                 std::pair<double,double> p = *it;
                 ashapeinput << p.first << ";"
                             << p.second << ";"
                             << no_unique_coor << ";"
                             << no_comm << '\n';
                 no_unique_coor++;
            }
            no_comm++;

        }
        std::stringstream oo;
        oo << "C:/Users/Dumex/Desktop/GeoEXp/ConvexHull_" << k << "/";
        saveConvexHulls(allpos, &oo.str()[0]);
        oo << "ashape_in.txt";
        writeStringToFile(&ashapeinput.str()[0],&oo.str()[0]);
    }
    Util::writeMatrixToFile(CC, "C:/Users/Dumex/Desktop/GeoEXp/cc.txt", ';');
}

//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------
//----------------------------------- GEOGRAPHICAL GRAPH END -------------------------------------
//------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------


int main(int argc, char *argv[])
{
    /* GEO EXP
     */
    geoGraphExp("C:/Users/Dumex/Desktop/GeoEXp/innerLondon.txt", 0.025);
    /* BLOCK MODEL EXP*/
  //  autoRunHiddenBlockModel_FixedDegree_FixedPrThreshold(1000,2000,10,2);
  // autoRunHiddenBlockModel();
  //  HiddenBlockModel_ClusteringCoefficient_EXP(2500);
    /* Main Terminal
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
        bool success = kMT(G, hierarchy, 1);
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
    */
    return 0;

}
