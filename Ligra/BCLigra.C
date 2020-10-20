        // This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
//#include <chrono>
//#include <unistd.h>
#include "ligra.h"
#include <vector>
#include "math.h"
#include "parallel.h"
#include "parseCommandLine.h"
#include <chrono>
#include <unistd.h>
#include <bits/stdc++.h> 
#include <ctime>
#include <time.h>
#include <stdio.h>
using namespace std;
#define BILLION  1000000000.0;
//time_t start1,end1;
double seconds;
int ops;
int opsUpdate;
int opsLookup;
int opsapp;
long sumedgeCount;
  typedef struct 
{
    int     secs;
    int     usecs;
}TIME_DIFF;

TIME_DIFF * my_difftime (struct timeval * start, struct timeval * end)
{
	TIME_DIFF * diff = (TIME_DIFF *) malloc ( sizeof (TIME_DIFF) );
 
	if (start->tv_sec == end->tv_sec) 
	{
        	diff->secs = 0;
        	diff->usecs = end->tv_usec - start->tv_usec;
    	}
   	else 
	{
        	diff->usecs = 1000000 - start->tv_usec;
        	diff->secs = end->tv_sec - (start->tv_sec + 1);
        	diff->usecs += end->tv_usec;
        	if (diff->usecs >= 1000000) 
		{
        	    diff->usecs -= 1000000;
	            diff->secs += 1;
	        }
	}
        return diff;
}
  struct timeval tv1, tv2;
  TIME_DIFF * difference1, * difference2 ;
  int dig1,temp1; 
  int dig2,temp2; 
  double duration1 = 0.0,duration2 = 0.0;
  
  double *TimePerOps;
struct timeval tv3, tv4, tv5;
TIME_DIFF * difference;
enum type {APP=0, ADDE, REME,ADDV, REMV};
typedef double fType;
struct BC_F {
  fType* NumPaths;
  bool* Visited;

  BC_F(fType* _NumPaths, bool* _Visited) : 
    NumPaths(_NumPaths), Visited(_Visited) {}
  inline bool update(intT s, intT d){ //Update function for forward phase
    fType oldV = NumPaths[d];
    NumPaths[d] += NumPaths[s];
    return oldV == 0.0;
  }
  inline bool updateAtomic (intT s, intT d) { //atomic Update, basically an add
    volatile fType oldV, newV; 
    do { 
      oldV = NumPaths[d]; newV = oldV + NumPaths[s];
    } while(!CAS(&NumPaths[d],oldV,newV));
    return oldV == 0.0;
  }
  inline bool cond (intT d) { return Visited[d] == 0; } //check if visited
};

struct BC_Back_F {
  fType* Dependencies;
  bool* Visited;
  BC_Back_F(fType* _Dependencies, bool* _Visited) : 
    Dependencies(_Dependencies), Visited(_Visited) {}
  inline bool update(intT s, intT d){ //Update function for backwards phase
    fType oldV = Dependencies[d];
    Dependencies[d] += Dependencies[s];
    return oldV == 0.0;
  }
  inline bool updateAtomic (intT s, intT d) { //atomic Update
    volatile fType oldV, newV;
    do {
      oldV = Dependencies[d];
      newV = oldV + Dependencies[s];
    } while(!CAS(&Dependencies[d],oldV,newV));
    return oldV == 0.0;
  }
  inline bool cond (intT d) { return Visited[d] == 0; } //check if visited
};

//vertex map function to mark visited vertexSubset
struct BC_Vertex_F {
  bool* Visited;
  BC_Vertex_F(bool* _Visited) : Visited(_Visited) {}
  inline bool operator() (intT i) {
    Visited[i] = 1;
    return 1;
  }
};

//vertex map function (used on backwards phase) to mark visited vertexSubset
//and add to Dependencies score
struct BC_Back_Vertex_F {
  bool* Visited;
  fType* Dependencies, *inverseNumPaths;
  BC_Back_Vertex_F(bool* _Visited, fType* _Dependencies, fType* _inverseNumPaths) : 
    Visited(_Visited), Dependencies(_Dependencies), inverseNumPaths(_inverseNumPaths) {}
  inline bool operator() (intT i) {
    Visited[i] = 1;
    Dependencies[i] += inverseNumPaths[i];
    return 1;
  }
};


template <class vertex>
void printGraph(graph<vertex> GA){
 vertex * v = GA.V;
 for (long i = 0; i < GA.n; i++){
  cout<<"V["<<i<<"] & "<<v[i].getOutDegree()<<" ->";
  for (long j = 0; j<v[i].getOutDegree(); j++)
   cout<<v[i].getOutNeighbor(j)<<" ";
  cout<<endl;
  } 
}
 ofstream fcoutcsv;
ofstream fcout;
time_t start1,end2, end1;
template <class vertex>
void BC(intT start, graph<vertex> GA, int nops, int wp,  string opsfile) {
uintT u,u_,wt;
int itr = 1;
intT n,m;
int x;
intE* edges;
vertex* vnew;
int *Opstime = new int[nops];
cout<<"Vertecies:"<<GA.n<<" Edges:"<<GA.m<<endl;
fcoutcsv<<GA.n<<","<<GA.m<<",";
ifstream cinops(opsfile);
 int nop;
 cinops>>nop;
 int optype;
 int ii;
 gettimeofday(&tv3,NULL);
 while(itr<=wp)
   { 
   n =GA.n;
   m = GA.m;
   vertex * v = GA.V;
   uintT j=0;
   cinops>>optype;
   if(optype == APP){
   uintT no = n;
   cinops>>u; 
   if(u>=n)  goto last1;
   uintT start =  u;
  intT threshold = GA.m/20;

  fType* NumPaths = newA(fType,n);
  {parallel_for(intT i=0;i<n;i++) NumPaths[i] = 0.0;}
  NumPaths[start] = 1.0;

  bool* Visited = newA(bool,n);
  {parallel_for(intT i=0;i<n;i++) Visited[i] = 0;}
  Visited[start] = 1;
  vertexSubset Frontier(n,start);
 
  vector<vertexSubset> Levels;
  Levels.push_back(Frontier);

  intT round = 0;
  timer t1,t2;
  while(!Frontier.isEmpty()){ //first phase
    round++;
    vertexSubset output = edgeMap(GA, Frontier, BC_F(NumPaths,Visited),threshold);
    vertexMap(output, BC_Vertex_F(Visited)); //mark visited
    Levels.push_back(output); //save frontier onto Levels
    Frontier = output;
  }

  fType* Dependencies = newA(fType,n);
  {parallel_for(intT i=0;i<n;i++) Dependencies[i] = 0.0;}

  //invert numpaths
  fType* inverseNumPaths = NumPaths;
  {parallel_for(intT i=0;i<n;i++) inverseNumPaths[i] = 1/inverseNumPaths[i];}

  Levels[round].del();
  //reuse Visited
  {parallel_for(intT i=0;i<n;i++) Visited[i] = 0;}
  Frontier = Levels[round-1];
  vertexMap(Frontier,BC_Back_Vertex_F(Visited,Dependencies,inverseNumPaths));

  //tranpose graph
  GA.transpose();
 
  for(intT r=round-2;r>=0;r--) { //backwards phase
    vertexSubset output = edgeMap(GA, Frontier, 
			      BC_Back_F(Dependencies,Visited),threshold);
    output.del();
    Frontier.del();
    Frontier = Levels[r]; //gets frontier from Levels array
    //vertex map to mark visited and update Dependencies scores
    vertexMap(Frontier,BC_Back_Vertex_F(Visited,Dependencies,inverseNumPaths));
  }
  
  Frontier.del();

  //Update dependencies scores
  parallel_for(intT i=0;i<n;i++) {
    Dependencies[i]=(Dependencies[i]-inverseNumPaths[i])/inverseNumPaths[i];
  }
  free(inverseNumPaths);
  free(Visited);
  free(Dependencies);
  opsapp++;
  ops++;
  } 
  else  if(optype == ADDV){
           cinops>>u; 
            if(u<n) goto last1;
            n = n + 1;
            vnew = newA(vertex,n); // new v
            edges = newA(intE,m); // new edge
            uintT j = 0;
            for (uintT i = 0; i < n-1; i++){
                vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex
	        for(intE k=0; k<v[i].getOutDegree(); k++, j++) { 
         	 edges[j] = v[i].getOutNeighbor(k); // copy out Neighbor
	        }
	        vnew[i].setOutNeighbors((edges+ (j - v[i].getOutDegree()))); // set out neighbors
           }
          vnew[n-1].setOutDegree(0); // set out degree to zero
          vnew[n-1].setOutNeighbors((edges+j));// set out neighbors j-1
          GA.del(); // free the old graph
          GA = graph<vertex>(vnew,n,m,edges); // create new graph with modified one
          opsUpdate++;
          ops++;
         }
     else if(optype == REMV){
           cinops>>u;
           if(u>=n) goto last1;
            vertex * v = GA.V;
           uintT outdegreeu = v[u].getOutDegree();
           bool flag = false;
           if(outdegreeu != 0){
             vnew = newA(vertex,n);
             edges = newA(intE,m-outdegreeu);
             uintT j = 0;
             for (uintT i = 0; i < n; i++){
	        if( i == u){ // removeV is equal to v[i], delete all its outgoing edges
         	   vnew[i].setOutDegree(0); // set the out degree to zero
	           vnew[i].setOutNeighbors(edges + j);// set out neighbors to zero
	           flag = true;
	        }        
	        else{
	          vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex out degree
	          for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
         	      edges[j] = v[i].getOutNeighbor(k); // copy the out going edges
	          }
	          vnew[i].setOutNeighbors((edges+ (j - v[i].getOutDegree()))); // copy the out neighbors
                 }
              }
              GA.del();
              m = m - outdegreeu;
              GA = graph<vertex>(vnew,n,m,edges); // create new graph with modified one
              }
          opsUpdate++;
          ops++;
        }

      else if(optype == ADDE){
           cinops>>u>>u_>>wt; 
           if(u >= n || u_ >=n) goto last1;
           vertex * v = GA.V;
            int loc = -1;
            for(uintT ii=0; ii<v[u].getOutDegree(); ii++) { // find the location
	             if( u_ <= v[u].getOutNeighbor(ii)){
	              loc = ii;
	              break;
	             }
             }
             if(loc >=0 && u_ == v[u].getOutNeighbor(loc)){// edge already present
         	   goto gotoend4;	
         	}
             else{
                
                  vnew = newA(vertex,n);
                  edges = newA(intE,m+1);
                  m = m + 1;
                  uintT j = 0;
                  for (uintT i = 0; i < n; i++){
                        if( i == u ){ 
                 	  vnew[i].setOutDegree(v[i].getOutDegree()+1); //increase the out degree by one  
	                  for(uintT k=0, k1=0; k1<v[i].getOutDegree(); k1++,j++) { 
                 	       edges[j] = v[i].getOutNeighbor(k); // copy the out going edges
                 	       k++;
                 	       
	                  }
	                  vnew[i].setOutNeighbors((edges+ (j - (v[i].getOutDegree()+1)))); // copy the out neighbors
	                }        
	                else{
	                  vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex out degree
	                  for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
                 	      edges[j] = v[i].getOutNeighbor(k); // copy the out going edges
	                  }
	                  vnew[i].setOutNeighbors((edges+ (j - v[i].getOutDegree()))); // copy the out neighbors
                         }
                  }
             GA.del();
             GA = graph<vertex>(vnew,n,m,edges);
             }	      
         
         gotoend4:
         opsUpdate++;
         ops++;
        } //End putE operation
     else if(optype == REME){
           cinops>>u>>u_; 
           if(u >= n || u_ >=n) goto last1;
           vertex * v = GA.V;
           uintT outdegreeu = v[u].getOutDegree();
            int loc = -1;
            for(uintT ii=0; ii<v[u].getOutDegree(); ii++) { // find the location
	             if( u_ == v[u].getOutNeighbor(ii)){
	              loc = ii;
	              break;
	             }
             }
             if(loc == -1 ){// edge not present
         	   goto gotoend3;	
         	}
             else{
                  vnew = newA(vertex,n);
                  edges = newA(intE,m-1);
                  m = m - 1;
                  uintT j = 0;
                  for (uintT i = 0; i < n; i++){
                        if( i == u ){ 
                 	  vnew[i].setOutDegree(v[i].getOutDegree()-1); //increase the out degree by one  
	                  for(uintT k=0; k<v[i].getOutDegree(); k++) { 
                 	      if(k != loc)
                 	          edges[j++] = v[i].getOutNeighbor(k); // copy the out going edges
	                  }
	                  vnew[i].setOutNeighbors((edges+ (j - (v[i].getOutDegree()-1)))); // copy the out neighbors
	                }        
	                else{
	                  vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex out degree
	                  for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
                 	      edges[j] = v[i].getOutNeighbor(k); // copy the out going edges
	                  }
	                  vnew[i].setOutNeighbors((edges+ (j - v[i].getOutDegree()))); // copy the out neighbors
                         }
                  
                  }
             }	      
             
         GA.del();
         GA = graph<vertex>(vnew,n,m,edges);
         gotoend3:
         opsUpdate++;
         ops++;
        } //End RemE operation
        
    last1:
    itr++;
 }
 gettimeofday(&tv4,NULL);
  //cin>>ii;
  difference1 = my_difftime (&tv3, &tv4);
  dig1 = 1;
  temp1 = difference1->usecs;
  while(temp1>=10)
  {	
	dig1++;
	temp1 = temp1/10;
  }
  temp1 =1;
  for(int i=1;i<=dig1;i++)
    temp1 = temp1 * 10;
  duration1 = (double) difference1->secs + ((double)difference1->usecs / (double)temp1);
  fcout << "Warm time: " << duration1 <<" secs."<<endl;
  fcoutcsv<<duration1<<",";
  
 gettimeofday(&tv4,NULL);
 while(itr<=nops)
   { 
   n =GA.n;
   m = GA.m;
   vertex * v = GA.V;
   uintT j=0;
   cinops>>optype;
   if(optype == APP){
   uintT no = n;
   cinops>>u;// = (rand() % (no-1)); 
   if(u>=n)  goto last2;
   uintT start =  u;
  intT threshold = GA.m/20;

  fType* NumPaths = newA(fType,n);
  {parallel_for(intT i=0;i<n;i++) NumPaths[i] = 0.0;}
  NumPaths[start] = 1.0;

  bool* Visited = newA(bool,n);
  {parallel_for(intT i=0;i<n;i++) Visited[i] = 0;}
  Visited[start] = 1;
  vertexSubset Frontier(n,start);
 
  vector<vertexSubset> Levels;
  Levels.push_back(Frontier);

  intT round = 0;
  timer t1,t2;
  while(!Frontier.isEmpty()){ //first phase
    round++;
    vertexSubset output = edgeMap(GA, Frontier, BC_F(NumPaths,Visited),threshold);
    vertexMap(output, BC_Vertex_F(Visited)); //mark visited
    Levels.push_back(output); //save frontier onto Levels
    Frontier = output;
  }

  fType* Dependencies = newA(fType,n);
  {parallel_for(intT i=0;i<n;i++) Dependencies[i] = 0.0;}

  //invert numpaths
  fType* inverseNumPaths = NumPaths;
  {parallel_for(intT i=0;i<n;i++) inverseNumPaths[i] = 1/inverseNumPaths[i];}

  Levels[round].del();
  //reuse Visited
  {parallel_for(intT i=0;i<n;i++) Visited[i] = 0;}
  Frontier = Levels[round-1];
  vertexMap(Frontier,BC_Back_Vertex_F(Visited,Dependencies,inverseNumPaths));

  //tranpose graph
  GA.transpose();
 
  for(intT r=round-2;r>=0;r--) { //backwards phase
    vertexSubset output = edgeMap(GA, Frontier, 
			      BC_Back_F(Dependencies,Visited),threshold);
    output.del();
    Frontier.del();
    Frontier = Levels[r]; //gets frontier from Levels array
    //vertex map to mark visited and update Dependencies scores
    vertexMap(Frontier,BC_Back_Vertex_F(Visited,Dependencies,inverseNumPaths));
  }
  
  Frontier.del();

  //Update dependencies scores
  parallel_for(intT i=0;i<n;i++) {
    Dependencies[i]=(Dependencies[i]-inverseNumPaths[i])/inverseNumPaths[i];
  }
  free(inverseNumPaths);
  free(Visited);
  free(Dependencies);
  opsapp++;
  ops++;
  } 
  else  if(optype == ADDV){
           cinops>>u; 
            if(u<n) goto last2;
            n = n + 1;
            vnew = newA(vertex,n); // new v
            edges = newA(intE,m); // new edge
            uintT j = 0;
            for (uintT i = 0; i < n-1; i++){
                vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex
	        for(intE k=0; k<v[i].getOutDegree(); k++, j++) { 
         	 edges[j] = v[i].getOutNeighbor(k); // copy out Neighbor
	        }
	        vnew[i].setOutNeighbors((edges+ (j - v[i].getOutDegree()))); // set out neighbors
           }
          vnew[n-1].setOutDegree(0); // set out degree to zero
          vnew[n-1].setOutNeighbors((edges+j));// set out neighbors j-1
          GA.del(); // free the old graph
          GA = graph<vertex>(vnew,n,m,edges); // create new graph with modified one
          opsUpdate++;
          ops++;
         }
     else if(optype == REMV){
           cinops>>u;
           if(u>=n) goto last2;
            vertex * v = GA.V;
           uintT outdegreeu = v[u].getOutDegree();
           bool flag = false;
           if(outdegreeu != 0){
             vnew = newA(vertex,n);
             edges = newA(intE,m-outdegreeu);
             uintT j = 0;
             for (uintT i = 0; i < n; i++){
	        if( i == u){ // removeV is equal to v[i], delete all its outgoing edges
         	   vnew[i].setOutDegree(0); // set the out degree to zero
	           vnew[i].setOutNeighbors(edges + j);// set out neighbors to zero
	           flag = true;
	        }        
	        else{
	          vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex out degree
	          for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
         	      edges[j] = v[i].getOutNeighbor(k); // copy the out going edges
	          }
	          vnew[i].setOutNeighbors((edges+ (j - v[i].getOutDegree()))); // copy the out neighbors
                 }
              }
              GA.del();
              m = m - outdegreeu;
              GA = graph<vertex>(vnew,n,m,edges); // create new graph with modified one
              }
          opsUpdate++;
          ops++;
        }

      else if(optype == ADDE){
           cinops>>u>>u_>>wt; 
           if(u >= n || u_ >=n) goto last2;
           vertex * v = GA.V;
            int loc = -1;
            for(uintT ii=0; ii<v[u].getOutDegree(); ii++) { // find the location
	             if( u_ <= v[u].getOutNeighbor(ii)){
	              loc = ii;
	              break;
	             }
             }

             if(loc >=0 && u_ == v[u].getOutNeighbor(loc)){// edge already present
         	   goto gotoend1;	
         	}
             else{
                
                  vnew = newA(vertex,n);
                  edges = newA(intE,m+1);
                  m = m + 1;
                  uintT j = 0;
                  for (uintT i = 0; i < n; i++){
                        if( i == u ){ 
                 	  vnew[i].setOutDegree(v[i].getOutDegree()+1); //increase the out degree by one  
	                  for(uintT k=0, k1=0; k1<v[i].getOutDegree(); k1++,j++) { 
                 	       edges[j] = v[i].getOutNeighbor(k); // copy the out going edges
                 	       k++;
                 	       
	                  }
	                  vnew[i].setOutNeighbors((edges+ (j - (v[i].getOutDegree()+1)))); // copy the out neighbors
	                }        
	                else{
	                  vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex out degree
	                  for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
                 	      edges[j] = v[i].getOutNeighbor(k); // copy the out going edges
	                  }
	                  vnew[i].setOutNeighbors((edges+ (j - v[i].getOutDegree()))); // copy the out neighbors
                         }
                  }
             GA.del();
             GA = graph<vertex>(vnew,n,m,edges);
             }	      
         
         gotoend1:
         opsUpdate++;
         ops++;
        } //End putE operation
     else if(optype == REME){
           cinops>>u>>u_; 
           if(u >= n || u_ >=n) goto last2;
           vertex * v = GA.V;
           uintT outdegreeu = v[u].getOutDegree();
            int loc = -1;
            for(uintT ii=0; ii<v[u].getOutDegree(); ii++) { // find the location
	             if( u_ == v[u].getOutNeighbor(ii)){
	              loc = ii;
	              break;
	             }
             }
             if(loc == -1 ){// edge not present
         	   goto gotoend2;	
         	}
             else{
                  vnew = newA(vertex,n);
                  edges = newA(intE,m-1);
                  m = m - 1;
                  uintT j = 0;
                  for (uintT i = 0; i < n; i++){
                        if( i == u ){ 
                 	  vnew[i].setOutDegree(v[i].getOutDegree()-1); //increase the out degree by one  
	                  for(uintT k=0; k<v[i].getOutDegree(); k++) { 
                 	      if(k != loc)
                 	          edges[j++] = v[i].getOutNeighbor(k); // copy the out going edges
	                  }
	                  vnew[i].setOutNeighbors((edges+ (j - (v[i].getOutDegree()-1)))); // copy the out neighbors
	                }        
	                else{
	                  vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex out degree
	                  for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
                 	      edges[j] = v[i].getOutNeighbor(k); // copy the out going edges
	                  }
	                  vnew[i].setOutNeighbors((edges+ (j - v[i].getOutDegree()))); // copy the out neighbors
                         }
                  
                  }
             }	      
             
         GA.del();
         GA = graph<vertex>(vnew,n,m,edges);
         gotoend2:
         opsUpdate++;
         ops++;
        } //End RemE operation
        
    last2:
    itr++;
 }
 gettimeofday(&tv5,NULL);
  difference2 = my_difftime (&tv4, &tv5);
  dig2 = 1;
  temp2 = difference2->usecs;
  while(temp2>=10)
  {	
	dig2++;
	temp2 = temp2/10;
  }
  temp2 =1;
  for(int i=1;i<=dig2;i++)
    temp2 = temp2 * 10;
  duration2 = (double) difference2->secs + ((double)difference2->usecs / (double)temp2);
  cout << "Execution time: " << duration2 <<" secs."<<endl;  
 }
 
 
int parallel_main(int argc, char* argv[]) {  
  commandLine P(argc,argv," [-s] <inFile>");
  char* iFile = P.getArgument(0);
  bool symmetric = P.getOptionValue("-s");
  bool binary = P.getOptionValue("-b");
  long start = P.getOptionLongValue("-r",0);
  long rounds = P.getOptionLongValue("-rounds",1);
  long nops=P.getOptionLongValue("-n",10);
  struct timeval tv1, tv2;
  TIME_DIFF * difference;
  int dig,temp; 
  double duration = 0.0;
  string opsfile=P.getOperatioFile("-f");
  //int nth = P.getOptionValue1("-nt");
  //string foutname = P.getOperatioFile("-fout");
  //fcout << "Number of Threads: " << nth << endl;
  int wp = nops *(0.05);
  graph<symmetricVertex> G = readGraph<symmetricVertex>(iFile,symmetric,binary); //symmetric graph
  BC((intT)start,G,nops, wp, opsfile);
}


