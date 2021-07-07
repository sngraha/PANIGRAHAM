PANIGRAHAM: Practical Non-blocking Graph Algorithms. 
==================

What is PANIGRAHAM?
-----------------


PANIGRAHAM is a package designed to support fully dynamic graph analytics. PANIGRAHAM  is composed of concurrent graph data structure that provides breadth-first search, single-source shortest-path, and betweenness centrality with concurrent dynamic updates of both edges and vertices.  The directory structure of the components is as follows:

    PGCn-ICn/ - Source code for the PGCn and PGICn along with make file.
    Ligra/ - Source code for the Ligra along with make file.
    Stinger/ - Source code for the BFS.
    datasets/ - Graph datasets and operation distribution files.

Compilation
--------

Recommended environment for PGCn-ICn

* g++ &gt;= 5.3.0

To compile PGCn-ICn

```
$ make 
```

To run
```
./main-<OP>-PGCn/PGICn <number of threads> <Threshod value> <input Graph> <operation distribution file>
```
Run Examples

for BFS-PGCn

```
$ ./main-BFS-PGCn 14 10 ../datasets/rMatGraph-1k ../datasets/gds9802-1M 1000
```
for BFS-PGICn
```
$ ./main-BFS-PGICn 14 10 ../datasets/rMatGraph-1k ../datasets/gds9802-1M 1000
```
for BC-PGCn
```
$ ./main-BC-PGCn 14 10 ../datasets/rMatGraph-1k ../datasets/gds9802-1M 1000
```
for SSSP-PGCn <input weighted graph> 
```
$ ./main-SSSP-PGCn 14 10 ../datasets/rMatWtGraph-1k ../datasets/gds9802-1M 1000
```


For Ligra
===========

Recommended environment
---------

* Intel icpc compiler
* g++ &gt;= 4.8.0 with support for Cilk+, 

To compile with g++ using Cilk, define the environment variable CILK. Then set the number of threads. More information can be found (https://github.com/jshun/ligra).

To run
```
export CILK_NWORKERS=<number of threads>

./<OP>Ligra -n <number of operations> -f <operation distribution file> -s <input Graph>
```
Run Examples

```
$ export CILK_NWORKERS=14

$ ./BFSLigra -n 1000 -f ../datasets/gds9802-1M -s ../datasets/rMat-1k
$ ./SSSPLigra -n 1000 -f ../datasets/gds9010-1M -s ../datasets/rMatWt-8k 
```

For Stinger
===========
Place bfs.c file from Stinger/ to /src/standalone/breadth_first_search then build the Stinger based on instructions given (https://github.com/stingergraph/stinger)

To run, place the operation distribution file in side bin/ directory then 
```
$ export OMP_NUM_THREADS=<number of threads>

./stinger_breadth_first_search <input Graph> <Op file> 
```
Run Examples 

```
$ export OMP_NUM_THREADS=14

$ ./stinger_breadth_first_search rMat-1k.bin a.1k.bin 
```

Input Format for PGCn and PGICn applications
-----------

The input format of unweighted graphs should be in edge list formats 

&lt;n>   # number of vertices
&lt;m>   # number of edges
&lt;src0 des0>  
&lt;src2 des2>  
...  
&lt;src(m-1) des(m-1)>  

The input format of weighted graphs should be in edge list formats 

&lt;n>   # number of vertices
&lt;m>   # number of edges
&lt;src0 des0 wt0>  
&lt;src2 des2 wt1>  
...  
&lt;src(m-1) des(m-1) wt(m-1)>  


Input Format for Ligra applications
-----------
AdjacencyGraph  
&lt;n>  
&lt;m>  
&lt;o0>  
&lt;o1>  
...  
&lt;o(n-1)>  
&lt;e0>  
&lt;e1>  
...  
&lt;e(m-1)>  


Input Format for Stinger application is same as Ligra in binary format.
-----------
&lt;endian check value>  
&lt;n>  
&lt;m>  
&lt;o0>  
&lt;o1>  
...  
&lt;o(n-1)>  
&lt;e0>  
&lt;e1>  
...  
&lt;e(m-1)>  




#Resources  
#-------- 
#Bapi Chatterjee, Sathya Peri and Muktikanta Sa. [Non-blocking Dynamic Unbounded Graphs with Worst-case Amortized Bounds](https://arxiv.org/abs/2003.01697). 


#If you have any questions, please contact: sngrahagraph@gmail.com, cs15resch11012@iith.ac.in

