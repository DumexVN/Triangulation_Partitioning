//-----------------------------------------------------------------------------
// Author   : N Vu
// Program  : Random Aggregation Clustering Type III.a
// File: Utility
//-----------------------------------------------------------------------------
// see readme.txt for more details


#ifndef UTIL_H
#define UTIL_H

#include <algorithm>
#include <vector>
#include <assert.h>

class Util
{
    typedef unsigned long long int uint64;
    typedef unsigned long int uint32;

private:

    static uint64 calA(std::vector<uint64> param)
    {
        uint64 a = param[0];
        a = a/2;
      //  qDebug() << "a:" << a;
        return a;
    }

    static uint64 calB(std::vector<uint64> param)
    {
        uint64 ai = param[0], nij_square = param[1];
        uint64 b = (ai-nij_square)/2;
     //   qDebug() << "b:" << b;
        return b;
    }

    static uint64 calC(std::vector<uint64> param)
    {
        uint64 ni_square = param[0], nij_square = param[1];
        uint64 c = (ni_square - nij_square)/2;
      //  qDebug() << "c:" << c;
        return c;
    }

    static uint64 calD(std::vector<uint64> param)
    {
        uint64 n_square = param[0],
                nij_square = param[1],
                ni_square = param[2],
                nj_square = param[3];
        uint64 d = n_square + nij_square - ni_square - nj_square; //type ii: diff diff
        d /= 2;
     //   qDebug() << "d:" << d;
        return d;
    }

    /// Calculate Adjusted Rand Index
    /// \brief calAdRand
    /// \param param
    /// \return
    ///
    static double calAdRand(std::vector<uint64> param)
    {
        uint64 ni_choose_2 = param[0],
                nj_choose_2 = param[1],
                n_choose_2 = param[2],
                nij_choose_2 = param[3];
        double nc = (double)ni_choose_2*nj_choose_2/n_choose_2;
        double nom = (double)(nij_choose_2 - nc);
        double sum = (double) (ni_choose_2 + nj_choose_2)/2;
        double denom = sum - nc;
        double ARI = nom/denom;
        return ARI;
    }

public:

    /// Write a matrix to file
    /// The matrix is in unconventional form
    /// each vector is the height, so it is a rather 'vertical' matrix instead of a horizontal
    ///
    ///
    template<typename T>
    static void writeMatrixToFile(const std::vector<std::vector<T> > &matrix,
                                  const char * outPath,
                                  const char &delim)
    {
        std::ofstream outfile;
        outfile.open(outPath);
        long int w = matrix.size(), h = matrix[0].size();
        for (int j = 0; j < h; j++)
        {
            for (int i = 0; i < w; i++) outfile << matrix[i][j] << delim;
            outfile << '\n';
        }
        outfile.close();
    }

    ///
    /// \brief fractionOfIntraEdge: calculate fraction of incorrect edge i.e. intra edge
    /// \param hier: the hierarchy
    /// \param l: number of edge per block
    /// \param k: number of edge per vertex
    /// \return: fraction
    ///
    static double fractionOfIntraEdge(const std::vector<std::pair<uint32, uint32> > &hier,
                                      const int &l, const int &k)
    {
        uint32 e = l*k*4; //total number of vertex
        uint32 incorrect = 0;
        for (int i = 0; i < hier.size(); i++)
        {
            std::pair<uint32,uint32> edge = hier.at(i);
            uint32 from = edge.first, to = edge.second;
            if ( (from/l) != (to/l) ) incorrect++;
        }
        return (double)incorrect/e;
    }


    ///
    /// \brief Util::NumberOfCommonElement
    /// \param vec1
    /// \param vec2
    /// \return the number of common element (or set intersection),
    /// IT IS ASSUMED THAT BOTH VECTORS ARE SORTED
    ///
    static unsigned int NumberOfCommonElement(const std::vector<unsigned long int> &vec1,
                                const std::vector<unsigned long int> &vec2)
    {
        unsigned int comm = 0;
        for(size_t i = 0, j = 0; i < vec1.size() && j < vec2.size();)
        {
            if (vec1[i] == vec2[j])
            {
                comm++;
                i++;
                j++;
            }
            else if (vec1[i] < vec2[j]) i++;
            else j++;
        }
        return comm;
    }


    typedef std::vector<uint32> community;
    typedef std::vector<community> AllCommunity;

    static std::vector<double> ComputePairWiseMatchingIndex(const long unsigned int &n, const AllCommunity &ground_truth_communities,
                                                                     const AllCommunity &large_result)
    {
         uint64 nij = 0, nij_minus = 0, nij_square = 0, nij_choose_2 = 0, n_choose_2 = 0;
         uint64 n_minus_1 = n-1, n_times_n_minus = n*n_minus_1;
         n_choose_2 = n_times_n_minus/2;

         size_t row = ground_truth_communities.size(), column = large_result.size();
         std::vector<uint64> ni, nj; //ni: sum row, nj: sum column
         for(size_t j = 0; j < column; j++)
             nj.push_back(0);

         for (size_t i = 0; i < row; i++)
         {
             uint64 sum_row = 0;
             community X = ground_truth_communities[i];
             for (size_t j = 0; j < column; j++)
             {
                 community Y = large_result[j];
                 uint32 entry = NumberOfCommonElement(X,Y); //sets are assumed to be sorted
                 uint64 entry_square = entry*entry,
                        entry_entryminus = entry*(entry-1);
                 nij_minus += entry_entryminus; // nij(nij-1)
                 nij_square += entry_square; // nij^2
                 nij_choose_2 += entry_entryminus/2; //(nij choose 2) for adjust rand
                 nij += entry;
                 sum_row += entry;
                 //
                 uint64 sum_col = nj[j];
                 uint64 new_sum = entry + sum_col;
                 nj[j] = new_sum;
             }
             ni.push_back(sum_row);
         }

         uint64 n_square = 0,
                 ni_sum = 0, //sum row
                 nj_sum =0, // sum column
                 ni_choose_2 = 0, //bionomial row
                 nj_choose_2 = 0, //bionomial column
                 ni_square = 0,  // sum each row square
                 nj_square = 0; // sum each column square

         for (size_t i = 0; i < ni.size(); i++)
         {
             uint64 entry = ni[i];
             uint64 entry_square = std::pow(entry,2);
             uint64 entry_choose_2 = entry*(entry-1)/2;
             ni_square += entry_square;
             ni_choose_2 += entry_choose_2;
             ni_sum+=  entry;
         }

         for (size_t i = 0; i < nj.size(); i++)
         {
             uint64 entry = nj[i];
             uint64 entry_square = entry*entry;
             uint64 entry_choose_2 = entry*(entry-1)/2;
             nj_square += entry_square;
             nj_choose_2 += entry_choose_2;
             nj_sum+=  entry;
         }

         //check sum
         uint64 sum_i = 0, sum_j = 0;
         for (size_t i = 0; i < ni.size(); i++) sum_i += ni[i];
         for (size_t j = 0; j < nj.size(); j++) sum_j += nj[j];
         assert((sum_i == sum_j && sum_i == n) && ("FATAL ERROR WHILE CALCULATING PAIRWISE MATCHING"));
         ni.clear();
         nj.clear();

         n_square = std::pow(n,2);

         std::vector<uint64> param_a,param_b,param_c,param_d;
         param_a.push_back(nij_minus);
         param_b.push_back(ni_square);
         param_b.push_back(nij_square);
         param_c.push_back(nj_square);
         param_c.push_back(nij_square);
         param_d.push_back(n_square);
         param_d.push_back(nij_square);
         param_d.push_back(ni_square);
         param_d.push_back(nj_square);

         uint64 a = calA(param_a);
         uint64 d = calD(param_d);
         uint64 c = calC(param_c); // type iii: same diff
         uint64 b = calB(param_b); //type iv: diff same

         double RAND = (double) (a+d)/(a+b+c+d);
         double Jaccard = (double) a/(a+b+c);

         //check sum
         assert((a+b+c+d == n_choose_2) && "FATAL ERROR WHILE CALCULATING PAIRWISE MATCHING");
      //   qDebug() << a/(a+b+c);
         //
         std::vector<uint64> param_ARI;
         param_ARI.push_back(ni_choose_2);
         param_ARI.push_back(nj_choose_2);
         param_ARI.push_back(n_choose_2);
         param_ARI.push_back(nij_choose_2);
         double ARI = calAdRand(param_ARI);

         std::vector<double> result;
         result.push_back(RAND);
         result.push_back(Jaccard);
         result.push_back(ARI);
         return result;
     }
};




#endif // UTIL_H
