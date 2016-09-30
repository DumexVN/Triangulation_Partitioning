# Triangulation_Partitioning
An Optimised Version of Triangulation Partitioning. 
1. input is assumed to be in edge_file format i.e. source_id[tab]target_id 
2. index assumed to starts from 0
3. to partition: -c input_file_path output_file_path
****************************
# To compare similarities between two partitions T & C (|T| = |C| = n)
-q partion1 partion2
Results are returned in a tuple (RAND_index, Jaccard_index, ARI_index)

Computation is done using a contigency table:
![Alt text](http://imgur.com/a/YEPwX "Contigency Table")

The Adjusted Rand Index is then given by:
![Alt text](http://imgur.com/a/lW6E9 "ARI")

For RAND, Jaccard, count:
- a, the number of pairs of elements in T that are in the same set in C and in the same set in T
- b, the number of pairs of elements in T that are in different sets in C and in different sets T
- c, the number of pairs of elements in T that are in the same set in C and in different sets in T
- d, the number of pairs of elements in T that are in different sets in C and in the same set in T

then Rand=(a+b)/(a+b+c+d); Jaccard=a/(a+b+c)
a,b,c,d can be calculated with equations from [1]
[1] Lawrence Hubert and Phipps Arabie (1985). "Comparing partitions". Journal of Classification. 2 (1): 193–218.
