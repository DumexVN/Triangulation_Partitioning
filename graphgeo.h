#ifndef GRAPHGEO_H
#define GRAPHGEO_H

#include <cmath>
#include <vector>
#include <cstdlib>

#include <boost/geometry.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/adapted/boost_tuple.hpp>

#include "stdlib.h"

class GraphGeo{
    typedef unsigned long int uint;
    typedef std::pair<double,double> geo;
    typedef std::pair<uint, geo> vertex;
public:

    GraphGeo(){ V = std::vector<vertex>();}

    ///
    /// \brief getPos: return the position for a collection of vertices, which is done through index look up
    /// \param vertex: indicies of queried vertices
    /// \param pos: its original geometrical position
    ///
    ///
    void getPos(const std::vector<uint> &vertex,
                std::vector<geo> &pos)
    {
        uint precision = 10000;
        //do work
        for (int i = 0; i < vertex.size(); i++)
        {
            uint v_id = vertex.at(i);
            assert((v_id < V.size()) && "INDEX OUT OF BOUND WHILE GETTING GEO POS\n");
            geo unrounded = V.at(v_id).second;
            double roundlat = roundf(unrounded.first * precision) / precision,
                   roundlong = roundf(unrounded.second * precision) / precision;
            pos.push_back(std::pair<double,double>(roundlat,roundlong));
        }
    }

    ///
    /// \brief convexHull: generate a convex hull covering all points given in &in
    /// \param in: input-sets of point
    /// \param hull: output, the convex hull
    ///
    void convexHull(const std::vector<geo> &in,
                    std::vector<geo> &hull)
    {
        //init boost geometric
        namespace bg = boost::geometry;
        typedef bg::model::point<double, 2, bg::cs::cartesian> point_t;
        typedef bg::model::linestring<point_t> linestring_t; //linestring is a collection of points
        typedef bg::model::polygon<point_t> poly;

        //load data into linestring
        linestring_t ls1;
        for (int i = 0; i < in.size(); i++) bg::append(ls1, point_t(in[i].first, in[i].second));
        //execute convex hull
        poly convexHull;
        bg::convex_hull(ls1, convexHull);
        //clearing
        for(auto it1 = boost::begin(boost::geometry::exterior_ring(convexHull));
            it1 != boost::end(boost::geometry::exterior_ring(convexHull));
            ++it1)
        {
            hull.push_back(geo(bg::get<0>(*it1),
                               bg::get<1>(*it1)));
        }
    }

    ///
    /// \brief generateGeoDiscGraph: generate a disc graph
    /// \param G: empty place holder for graph
    /// \param d: limiting distance i.e. 50m
    /// \return
    ///
    bool generateGeoDiscGraph(std::vector<std::vector<uint> > &G,
                              const double &d )
    {   //init
        uint n = V.size(),
             e = 0;
        for (int i = 0; i < n; i++) G.push_back(std::vector<uint>());
        //starting
        for (int i = 0; i < n; i++)
        {
            vertex v = V.at(i); //get v
            for (int j = i+1; j < n; j++)
            {
                vertex u = V.at(j); //get u
                double distance = distanceHarvensine(v.second.first, u.second.first, //look quite bad
                                                     v.second.second, u.second.second); //calculate distance
                if (distance <= d)
                {
                    G[i].push_back(j);
                    G[j].push_back(i);
                    e++;
                }
            }
        }
        printf("Geometric Graph Generated: V %d - E %d", n, e);
        return true;
    }

    ///
    /// \brief generateGeoDiscSubGraph: generate a subgraph H of G
    /// \param H:place holder for subgraph H
    /// \param subV: array of vertices included in the subgraph
    /// \param d: limiting distance i.e. 50m
    /// \return
    ///
    bool generateGeoDiscSubGraph(std::vector<std::vector<uint> > &H,
                                 const std::vector<uint> &subV,
                                 const double &d)
    {
        H = std::vector<std::vector<uint> >(subV.size(), std::vector<uint>());
        for (int i = 0; i < subV.size(); i++)
        {
            uint v_index = subV.at(i);
            vertex v = V.at(v_index);
            for (int j = i+1; j < subV.size(); j++)
            {
                uint u_index = subV.at(j);
                vertex u = V.at(u_index); //get u
                double distance = distanceHarvensine(v.second.first, u.second.first, //look quite bad
                                                     v.second.second, u.second.second); //calculate distance
                if (distance <= d)
                {
                    H[i].push_back(j);
                    H[j].push_back(i);
                }
            }
        }
    }

    ///
    /// \brief readSpatialData: read spatial data input
    ///
    bool readSpatialData(const char * inPath)
    {
        std::ifstream infile(inPath);
        if (!infile.good())
        {
            printf("File Erorr! Terminating ...");
            return false;
        }
        std::string line;
        std::vector<vertex> vList;
        int line_count = 0;
        while (std::getline(infile, line))
        {
            if (line_count == 0){line_count++;continue;}
            std::istringstream iss(line);
            std::vector<std::string> tokenVec;
            std::string token;
            while(std::getline(iss, token, ';'))//can be replaced by line, putting iss for clearer intention
            {
                tokenVec.push_back(token);
            }
            uint id = std::atoi(tokenVec.at(2).c_str());
            double lat = std::atof(tokenVec.at(0).c_str()),
                   lon = std::atof(tokenVec.at(1).c_str());
            vList.push_back(vertex(id, geo(lat,lon))); //create new vertex and append
        }
        std::sort(vList.begin(), vList.end(), [](const vertex &v, const vertex &u) {
               return v.first < u.first;
           }); //finally, sort the vector according to the order of vertex' indices
        //check sum
        uint maxVID = vList.at(vList.size()-1).first + 1; //if index starts from 0
        assert((maxVID == vList.size()) && "Vertex indices does not correspond to the number of vertex!");
        V = vList;
        printf("Finished Reading Spatial Graph File!\n");
        return true;
    }

    ///
    /// \brief degreeFromDecimal
    /// \param x
    /// \return
    ///
    double degreeFromDecimal(const double &x)
    {
        const double Pi = 3.141592653589793;
        return  (x * Pi / 180);
    }

    ///
    /// \brief distanceHarvensine: calculate earth surface distance
    /// \param lat1:    latitude
    /// \param lat2:    lattitude
    /// \param lon1:    longitude
    /// \param lon2:    longitude
    /// \return
    ///
    double distanceHarvensine(const double &lat1,
                              const double &lat2,
                              const double &lon1,
                              const double &lon2)
    {
        const double R = 6371.0;

        double dlon = degreeFromDecimal(lon2 - lon1);
        double dlat = degreeFromDecimal(lat2 - lat1);
        double deglat1 = degreeFromDecimal(lat1);
        double deglat2 = degreeFromDecimal(lat2);

        double a = sin(dlat/2)*sin(dlat/2) + cos(deglat1) * cos(deglat2) * sin(dlon/2)*sin(dlon/2);
        double c = 2 * atan2( sqrt(a), sqrt(1-a) );
        double d = R * c;

        return d;
    }

private:
    std::vector<vertex> V;
};
#endif // GRAPHGEO_H
