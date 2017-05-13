#ifndef GRAPHUTIL_H
#define GRAPHUTIL_H

#include <vector>
#include <array>
#include <random>
#include <chrono>
#include <cmath>

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>
#include <boost/graph/exterior_property.hpp>
#include <boost/graph/clustering_coefficient.hpp>

#include "graphgeo.h"
#include "util.h"

class GraphUtil{
    typedef unsigned long int uint;
    typedef std::pair<uint, uint> edge;
    typedef std::vector<uint> adj_list;
    typedef std::vector<adj_list> Graph;


public:
    ///
    /// \brief degreeDistribution: get a graph degree distribution
    /// \param G: graph G
    /// \param deg: empty holder for degree
    ///
    static void degreeDistribution(const Graph &G,
                                   std::vector<uint> &deg)
    {
        for (int i = 0; i < G.size(); i++)
        {
            deg.push_back(G[i].size());
        }
    }

    ///
    /// \brief Edge Weight Ratio of the hidden block model
    /// Ratio(v)
    /// \param G: input graph G
    /// \param v_id: index of queried vertex
    /// \param comm: v's community
    /// \return: ratio weight_intra/weight_inter
    ///
    static double vertexHiddenBlockWeightRatio(const Graph &G,
                                           const uint &v_id,
                                           const uint &l)
    {
        assert(v_id < G.size());
        adj_list v_adj = G.at(v_id);
        uint dv = v_adj.size(),
             intra_triangle = 0,
             inter_triangle = 0;
        if (v_adj.size() <= 1)  return 0;
        for (int i = 0; i < dv; i++)
        {
            uint u_id = v_adj[i]; //adjacent index
            adj_list u_adj = G.at(u_id);    //u adjacent
            if (u_id/l == v_id/l) //same community
            {
                intra_triangle += Util::NumberOfCommonElement(v_adj, u_adj);
            }
            else
            {
                inter_triangle += Util::NumberOfCommonElement(v_adj, u_adj);
            }
        }
        return (double) intra_triangle/inter_triangle;
    }


    ///
    /// \brief vertexBlockClusteringCoefficient: get a clustering coefficnet for a specific vertex
    /// \param G: input graph G
    /// \param v_id: index of queried vertex
    /// \param comm: v's community
    /// \return: ratio weight_intra/weight_inter
    ///
    static double vertexFractionOfTriangle(const Graph &G,
                                           const uint &v_id,
                                           const std::vector<uint> &comm)
    {
        assert(v_id < G.size());
        adj_list v_adj = G.at(v_id);
        uint dv = v_adj.size(),
             intra_triangle = 0,
             inter_triangle = 0;
        if (v_adj.size() <= 1)  return 0;
        for (int i = 0; i < dv; i++)
        {
            uint u_id = v_adj[i]; //adjacent index
            adj_list u_adj = G.at(u_id);    //u adjacent
            if (std::find(comm.begin(), comm.end(), u_id) != comm.end()) //same community
            {
                intra_triangle += Util::NumberOfCommonElement(v_adj, u_adj);
            }
            else
            {
                inter_triangle += Util::NumberOfCommonElement(v_adj, u_adj);
            }
        }
        return (double) intra_triangle/inter_triangle;
    }

    /// *** UNOPTIMISED ***
    /// \brief interClusteringCoeffcient: clustering coeffcient of a vertex with inter-edges
    /// \param G: graph G
    /// \param cc: place holder for clustering coefficient
    /// \return: mean
    ///
    static double interClusteringCoeffcient(const Graph &G,
                                            const std::vector<std::vector<uint> > communities,
                                            std::vector<double> &cc)
    {
        //init
        cc = std::vector<double>(G.size(), 0.0);
        double cc_sum = 0.0;
        //start
        for (int i = 0; i < communities.size(); i++)
        {
            //init some useful data structures
            std::vector<uint> a_community = communities.at(i);
            std::sort(a_community.begin(), a_community.end());
            std::set<uint> a_community_unique(a_community.begin(), a_community.end());
            //working on a community
            for (int j = 0; j < a_community.size(); j++)
            {
                uint v_id = a_community.at(j); //query a vertex
                adj_list v_adj = G.at(v_id); //get v adjacency list
                if (v_adj.size() <= 1){ cc[v_id] = 0; continue;}
                uint closed_triangles = 0;
                for (int k = 0; k < v_adj.size(); k++)
                {
                    if (a_community_unique.find(v_adj[k]) == a_community_unique.end()) //if a neighbour is inside the community
                        closed_triangles += Util::NumberOfCommonElement(v_adj, G.at(v_adj.at(k)));
                }
                double clustering_coeffcient = (double) closed_triangles/(v_adj.size()*(v_adj.size()-1));
                cc[v_id] = clustering_coeffcient;
                cc_sum += clustering_coeffcient;
            }
        }
        return (cc_sum/G.size());
    }

    /// *** UNOPTIMISED ***
    /// \brief intraClusteringCoefficient: clustering coeffcient of a vertex with intra-edge
    /// \param G: graph G
    /// \param cc: place holder for clustering coefficient
    /// \return: means
    ///
    static double intraClusteringCoefficient(const Graph &G,
                                             const std::vector<std::vector<uint> > communities,
                                             std::vector<double> &cc)
    {
        //init
        cc = std::vector<double>(G.size(), 0.0);
        double cc_sum = 0.0;
        //start
        for (int i = 0; i < communities.size(); i++)
        {
            //init some useful data structures
            std::vector<uint> a_community = communities.at(i);
            std::sort(a_community.begin(), a_community.end());
            std::set<uint> a_community_unique(a_community.begin(), a_community.end());
            //working on a community
            for (int j = 0; j < a_community.size(); j++)
            {
                uint v_id = a_community.at(j); //query a vertex
                adj_list v_adj = G.at(v_id); //get v adjacency list
                if (v_adj.size() <= 1){ cc[v_id] = 0; continue;}
                uint closed_triangles = 0;
                for (int k = 0; k < v_adj.size(); k++)
                {
                    if (a_community_unique.find(v_adj[k]) != a_community_unique.end()) //if a neighbour is inside the community
                        closed_triangles += Util::NumberOfCommonElement(v_adj, G.at(v_adj.at(k)));
                }
                double clustering_coeffcient = (double) closed_triangles/(v_adj.size()*(v_adj.size()-1));
                cc[v_id] = clustering_coeffcient;
                cc_sum += clustering_coeffcient;
            }
        }
        return (cc_sum/G.size());
    }


    ///
    /// \brief clusteringCoefficient: calculate the local clustering coefficient for each vertex
    /// \param G: graph G
    /// \param cc: vector for scores
    /// \return: the average clustering coefficient i.e. sum_{cc(v)}/|V|
    ///
    static double clusteringCoefficient(const Graph &G, std::vector<double> &cc)
    {
        //init
        uint V = G.size();
        std::vector<edge> uniqueE;
        uniqueEdgeList(G, uniqueE); //get the unique edge list
        //create boost graph
        typedef boost::adjacency_list<boost::setS, boost::vecS, boost::undirectedS> BoostGraph;
        BoostGraph BG;
        for(size_t i = 0; i < V; i++)   boost::add_vertex(BG);

        for(auto &e : uniqueE)
        {
            boost::add_edge(e.first, e.second, BG);
        }
        assert((boost::num_edges(BG) == uniqueE.size()) && "ERROR: WHILE CALCULATING COEFFICIENTS: NUM EDGES DOES NOT MATCH");
        // The clustering property, container, and map define the containment
        // and abstract accessor for the clustering coefficients of vertices.
        typedef boost::exterior_vertex_property<BoostGraph, float> ClusteringProperty;
        typedef ClusteringProperty::container_type ClusteringContainer;
        typedef ClusteringProperty::map_type ClusteringMap;

        // Compute the clustering coefficients of each vertex in the graph
        ClusteringContainer coefs(V);
        ClusteringMap cm(coefs, BG);
        float average = boost::all_clustering_coefficients(BG, cm);

        // Print the clustering coefficient of each vertex.
        boost::graph_traits<BoostGraph>::vertex_iterator i, end;
        for(boost::tie(i, end) = boost::vertices(BG); i != end; ++i) {
           cc.push_back(boost::get(cm, *i));
        }
        return average;
    }

    ///
    /// \brief uniqueEdgeList: produce a unique edge list of G (since G is adjacency list)
    /// \param G: input graph G
    /// \param E: empty result holder
    ///
    static void uniqueEdgeList(const Graph &G, std::vector<edge> &uniqueE)
    {
        //init
        uint V = G.size(), //number of vertices
             E = 0;//number of edges
        //starting
        for (int i = 0; i < V; i++)
        {
            adj_list v_adj = G.at(i); //v adjacency list
            E += v_adj.size(); // E += d(v)
            for (int j = 0; j < v_adj.size(); j++)
            {
                uint u_index = v_adj.at(j); //index of the neighbour u
                if (i < u_index)    uniqueE.push_back(edge(i, u_index));  //since the adjacency list is sorted
            }
        }
        assert((uniqueE.size() == E/2 ) && "ERROR: WHILE GENERATING UNIQUE EDGES: NUMBER OF EDGES DOES NOT MATCH");
    }

    ///
    /// \brief sortVerticesIntoCommunities
    /// \param num_comm: number of communities
    /// \param membership: the complete list of vertices membership
    /// i.e. v[0] = 1 means vertex 0 belongs in community 1
    /// \param communities: empty result holder
    /// communities[1] = 0: vertex 0 now is put in communities[1]
    ///
    static void sortVerticesIntoCommunities(const int &num_comm,
                                            const std::vector<uint> &membership,
                                            std::vector<std::vector<uint> > &communities)
    {
        //init empty comm
        for (int i = 0; i < num_comm; i++)
        {
            std::vector<uint> comm;
            communities.push_back(comm);
        }
        //putting vertex in correct comm
        for (int i = 0; i < membership.size(); i++)
        {
            uint comm_index = membership[i]; //community of vertex i
            communities[comm_index].push_back(i); //putting i into its designated comm
        }
    }

    ///
    /// \brief saveGraphAsEdgeList: save graph as edge list
    /// *NOTE* this method produces double edge i.e. both edge uv and vu are present
    /// \param G: graph G
    /// \param outPath: output path
    /// \param delim: delimiter
    ///
    static void saveGraphAsEdgeList(const Graph &G, const char *outPath,
                                    const char &delim)
    {
        if (outPath == NULL)
        {
            printf("Output Path is not defined!\n");
            return;
        }
        printf("- Writing Output File");
        std::ofstream outfile;
        outfile.open(outPath);
        for (int i = 0; i < G.size(); i++)
        {
            adj_list v_adj = G.at(i); //ith vertex adjacency list
            std::stringstream out;
            for (int j = 0; j < v_adj.size(); j++)
            {
                uint v_neighbour = v_adj.at(j);
                out << i << delim << v_neighbour << '\n'; //concat
            }
            outfile << out.str();
        }
    }

    ///
    /// \brief saveGraphAsEdgeList: only unique edge is present
    /// \param G: graph G
    /// \param outPath: output path
    /// \param delim: delimiter
    ///
    static void saveGraphAsUniqueEdgeList(const Graph &G, const char *outPath,
                                    const char &delim)
    {
        if (outPath == NULL)
        {
            printf("Output Path is not defined!\n");
            return;
        }
        printf("- Writing Unique Edge List File");
        std::ofstream outfile;
        outfile.open(outPath);
        std::vector<edge> uniqueEdge;
        uniqueEdgeList(G, uniqueEdge);
        for (int i = 0; i < uniqueEdge.size(); i++)
        {
            outfile << uniqueEdge.at(i).first << delim << uniqueEdge.at(i).second << '\n';
        }
    }

    ///
    /// \brief saveGraphAsCSV: saving graph as adjacency list as .csv
    /// \param G
    ///
    static void saveGraphAsAdjacencyListCSV(const Graph &G, const char * outPath)
    {
        if (outPath == NULL)
        {
            printf("Output Path is not defined!\n");
            return;
        }
        printf("- Writing Output File");
        std::ofstream outfile;
        outfile.open(outPath);
        for (int i = 0; i < G.size(); i++)
        {
            adj_list v_adj = G.at(i); //ith vertex adjacency list
            std::stringstream out;
            out << i << ";"; //prepare outstring
            for (int j = 0; j < v_adj.size(); j++)
            {
                uint v_neighbour = v_adj.at(j);
                out << v_neighbour << ";"; //concat
            }
            outfile << out.str() << '\n';
        }
    }

    ///
    /// \brief saveGraphAsGraphML: saving graph as graphML
    /// graphML is in form of XML, separating vertex and edge attributes
    /// use this when attributes are required for i.e. Gephi
    /// \param l: number of vertex per block
    /// \param outPath: output path
    /// \param hier: hierarchy edge
    ///
    static void saveOutputAsGraphML(const std::vector<edge> &hier, const int &l, const char * outPath){
        std::ofstream outfile;
        if (outPath == NULL)
        {
            printf("Output Path is not definedn");
            return;
        }
        outfile.open (outPath);
        printf("- Writing Output File");
        outfile << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>" << '\n'
            <<   "<graphml>" << '\n';
        //sample
        /*
         * <key id="d0" for="node" attr.name="color" attr.type="string">
         * <default>yellow</default>
         * </key>
         */
        outfile << "<key attr.name=\"x\" attr.type=\"float\" for=\"node\" id=\"x\"/>\n";
        outfile << "<key attr.name=\"y\" attr.type=\"float\" for=\"node\" id=\"y\"/>\n";
        outfile << "<graph id=\"G\" edgedefault=\"undirected\">" << '\n';
        //begin writing vertex's position
        unsigned int n = l*4;
        //getting layout
        std::vector<std::pair<double, double> > v_pos = layoutBlockModel(l);
        for (int i = 0; i < n; i++ )
        {
            //sample
            /*
             * </node>
             * <node id="1">
             * <data key="x">341.79462</data>
             * <data key="y">457.72717</data>
             * </node>
             */
            outfile << "<node id=\""<<i<<"\">\n";
            outfile << "<data key=\"x\">"<<v_pos.at(i).first<<"</data>\n";
            outfile << "<data key=\"y\">"<<v_pos.at(i).second<<"</data>\n";
            outfile << "</node>\n";
        }
        // now writing edge

        for (size_t i = 0; i < hier.size(); i++)
        {
            //sample
            /*
             * <edge id="e1" source="n0" target="n1"/>
             */
            outfile << "<edge id=\""<<i<<"\" source=\""<<hier[i].first<<
                                    "\" target=\""<<hier[i].second<<"\"/>\n" ;
        }
        outfile << "</graph>\n</graphml>";
        outfile.close();
        printf("Output Hierarchy Saved!");
    }

    ///
    /// \brief layoutBlockModel: layout block model, presumably 4 blocks
    /// \param l: number of vertices per block
    /// \return
    ///
    static std::vector<std::pair<double, double> > layoutBlockModel(const int &l){
        //init data structures
        typedef std::pair<double,double> pos;
        std::vector<pos> vec_pos;

        //seed generator
        std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
        auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);
        auto epoch = now_ms.time_since_epoch();
        auto value = std::chrono::duration_cast<std::chrono::milliseconds>(epoch);
        long seed1 = value.count();
        std::mt19937 gen(seed1);
        std::uniform_real_distribution<double> rad(0, 1);
        std::uniform_real_distribution<double> angle(0, 360);

        pos centres [4] = {pos(2000.0,2000.0),
                           pos(2000.0,0.0),
                           pos(0.0,2000.0),
                           pos(0.0,0.0)};
        for(int i = 0; i < 4; i++)
        {
            //centre of the balls
            unsigned int UNIT_BALL_RAD = 500;
            pos centre = centres[i];
            for (int j = 0; j < l; j++)
            {
                double r = rad(gen);
                double phi = angle(gen);
                double x = sqrt(r)*cos(phi)*UNIT_BALL_RAD + centre.first; //rescale and relocate to middle of screen
                double y = sqrt(r)*sin(phi)*UNIT_BALL_RAD + centre.second; //rescale and relocate to middle of screen
                pos v(x,y);
                vec_pos.push_back(v);
            }
        }
        return vec_pos;
    }

};

#endif // GRAPHUTIL_H
