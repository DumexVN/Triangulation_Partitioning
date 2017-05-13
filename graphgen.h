
//-----------------------------------------------------------------------------
// Author   : N Vu
// Program  : Random Aggregation Clustering Type III.a
// File: Utility
//-----------------------------------------------------------------------------
// see readme.txt for more details

#ifndef GRAPHGEN_H
#define GRAPHGEN_H

#include <random>
#include <chrono>
#include <vector>
#include <assert.h>

class GraphGen
{
    typedef unsigned long int uint;
    //typedef std::pair<uint, uint> edge;
    typedef std::vector<uint> adj_list;
    typedef std::vector<adj_list> Graph;

private:
    std::mt19937 gen;

public:

    GraphGen()
    {
        //seed generator
        std::chrono::time_point<std::chrono::system_clock> now = std::chrono::system_clock::now();
        auto now_ms = std::chrono::time_point_cast<std::chrono::milliseconds>(now);
        auto epoch = now_ms.time_since_epoch();
        auto value = std::chrono::duration_cast<std::chrono::milliseconds>(epoch);
        long seed1 = value.count();
        std::mt19937 generator(seed1);
        gen = generator;
    }

    ///
    /// \brief generateGnp: generate Gnp
    /// \param G: empty Graph
    /// \param n: number of vertices
    /// \param p: probability of success
    ///
    void generateGnp(Graph & G, const uint &n, const double &p)
    {
        //initialise  full adjacency list
        for (int i = 0; i < n; i++)
        {
            adj_list empty;
            G.push_back(empty);
        }
        //generate
        std::uniform_real_distribution<> dist(0, 1);
        for (int i = 0; i < n; i++)
        {
            for (int j = i+1; j < n; j++)
            {
                double ran = dist(gen);
                if (ran <= p){ //success
                    G[i].push_back(j);
                    G[j].push_back(i);
                }
            }
        }
    }

    ///
    /// \brief generateGnp_BlockModel: generate hidden Gnp block model
    /// \param G: empty Graph G
    /// \param l: number of vertices per block
    /// \param q: probability of edges inside block
    /// \param globalp: probability of edges between block
    ///
    void generateGnp_BlockModel(Graph & G, const uint &l,
                                const double &q, const double &globalp)
    {
        uint n = l*4;
        double r = q + globalp;
        //initialise  full adjacency list
        for (int i = 0; i < n; i++)
        {
            adj_list empty;
            G.push_back(empty);
        }
        //generating edges
        std::uniform_real_distribution<> dist(0, 1);
        for (int i = 0; i < n; i++)
        {
            int cv = i/l; //block of v
            for (int j = i+1; j < n; j++)
            {
                bool connect = false;
                int cj = j/l; //block of j
                double ran = dist(gen);
                if (cv == cj && ran <= r)   connect = true;
                else if (cv != cj && ran <= globalp)    connect = true;
                /*
                if (ran <= globalp)  connect = true; //global edge
                if (!connect && (cv == cj) && (ran <= q))   connect = true; //within-block edge*/
                if (connect){ //success
                    G[i].push_back(j);
                    G[j].push_back(i);
                }
            }
        }
        printf("Gnp Block Model Generated: n = %d, p = %f, q = %f \n ", n, globalp, q);
    }
};



#endif // GRAPHGEN_H
