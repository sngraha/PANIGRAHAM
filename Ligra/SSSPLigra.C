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
#include "ligra.h"
#include "gettime.h"
#include "parseCommandLine.h"
#define BILLION  1000000000.0;
using namespace std;
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
struct timeval tv3, tv4, tv5;
TIME_DIFF * difference;
enum type {APP=0,ADDE, REME,ADDV, REMV};
struct BF_F {
  int* ShortestPathLen;
  int* Visited;
  BF_F(int* _ShortestPathLen, int* _Visited) : 
    ShortestPathLen(_ShortestPathLen), Visited(_Visited) {}
  inline bool update (intT s, intT d, intE edgeLen) { //Update ShortestPathLen if found a shorter path
    intE newDist = ShortestPathLen[s] + edgeLen;
    if(ShortestPathLen[d] > newDist) {
      ShortestPathLen[d] = newDist;
      if(Visited[d] == 0) { Visited[d] = 1 ; return 1;}
    }
    return 0;
  }
  inline bool updateAtomic (intT s, intT d, intE edgeLen){ //atomic Update
    intE newDist = ShortestPathLen[s] + edgeLen;
    return (writeMin(&ShortestPathLen[d],newDist) &&
	    CAS(&Visited[d],0,1));
  }
  inline bool cond (intT d) { return cond_true(d); } //does nothing
};

//reset visited vertices
struct BF_Vertex_F {
  int* Visited;
  BF_Vertex_F(int* _Visited) : Visited(_Visited) {}
  inline bool operator() (intT i){
    Visited[i] = 0;
    return 1;
  }
};

template <class vertex>
void BellmanFord(intT start, wghGraph<vertex> GA, int nops, string opsfile) {
  uintT u,u_;
int itr = 1;
intT n,m;
int x;
intE* edges;
vertex* vnew;
ops = 0;
cout<<"Vertecies:"<<GA.n<<" Edges:"<<GA.m<<endl;
ifstream cinops(opsfile);
int nop;
 cinops>>nop;
 int optype;
//printGraph(GA);
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
  cinops>>u;
  if(u>=n)  goto last; 
   uintT start =  u;
  int* ShortestPathLen = newA(int,n);
  {parallel_for(intT i=0;i<n;i++) ShortestPathLen[i] = INT_MAX/2;}
  ShortestPathLen[start] = 0;

  int* Visited = newA(int,n);
  {parallel_for(intT i=0;i<n;i++) Visited[i] = 0;}

  vertexSubset Frontier(n,start); //initial frontier

  intT round = 0;
  intT numVisited = 0;
  while(!Frontier.isEmpty()){
    round++;
    if(round == n) {
      //negative weight cycle
      {parallel_for(intT i=0;i<n;i++) ShortestPathLen[i] = -(INT_MAX/2);}
      break;
    }
     numVisited+=Frontier.numNonzeros();
            sumedgeCount += numVisited;
    vertexSubset output = edgeMap(GA, Frontier, BF_F(ShortestPathLen,Visited), GA.m/20,DENSE_FORWARD);
    vertexMap(output,BF_Vertex_F(Visited));
    Frontier.del();
    Frontier = output;
  } 

  Frontier.del();
  free(Visited);
  opsapp++;
  ops++;
  } 
  else  if(optype == ADDV){
           cinops>>u; 
           if(u>=n) goto last;
      n = n + 1;
    vnew = newA(vertex,n); // new v
    edges = newA(intE,2*m); // new edge
    uintT j = 0;
    for (uintT i = 0; i < n-1; i++){
        vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex
	for(intE k=0; k<v[i].getOutDegree(); k++, j++) { 
 	 edges[2*j] = v[i].getOutNeighbor(k); // copy out Neighbor
 	 edges[2*j+1] = v[i].getOutWeight(k); // copy out Neighbor
	}
	vnew[i].setOutNeighbors((edges+ 2*(j - v[i].getOutDegree()))); // set out neighbors
   }
  vnew[n-1].setOutDegree(0); // set out degree to zero
  vnew[n-1].setOutNeighbors((edges+2*j-1));// set out neighbors
  GA.del(); // free the old graph
  GA = wghGraph<vertex>(vnew,n,m,edges); // create new graph with modified one
  opsUpdate++;
  ops++;
 } // End of if
else if(optype == REMV){
           cinops>>u;
           if(u>=n) goto last;
    vertex * v = GA.V;
   uintT outdegreeu = v[u].getOutDegree();
   if(outdegreeu != 0){
     //n = n - 1;
     vnew = newA(vertex,n);
     edges = newA(intE,2*m);
     uintT j = 0;
     for (uintT i = 0; i < n; i++){
	if( i == u){ // removeV is equal to v[i], delete all its outgoing edges
 	   vnew[i].setOutDegree(0); // set the out degree to zero
	   vnew[i].setOutNeighbors(edges + 2*j);// set out neighbors to zero
	}        
	else{
	  vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex out degree
	  for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
 	      edges[2*j] = v[i].getOutNeighbor(k); // copy the out going edges
 	      edges[2*j+1] = v[i].getOutWeight(k); // copy out Neighbor
	  }
	  vnew[i].setOutNeighbors((edges+ 2*(j - v[i].getOutDegree()))); // copy the out neighbors
         }
      }
      GA.del();
      GA = wghGraph<vertex>(vnew,n,m,edges); // create new graph with modified one
      }
  opsUpdate++;
  ops++;
  
 }
 else if(optype == ADDE){
 uintT wt;
           cinops>>u>>u_>>wt;
           if(u >= n || u_ >=n)
            goto last;
  vertex * v = GA.V;
   uintT outdegreeu = v[u].getOutDegree();
 
   if(outdegreeu != 0){
     m = m + 1;
     vnew = newA(vertex,n);
     edges = newA(intE,2*m);
     uintT j = 0;
     for (uintT i = 0; i < n; i++){
	if( i == u){ 
	  int loc = -1;
	  for(uintT ii=0; ii<v[i].getOutDegree(); ii++) { // find the location
	     if( u_ <= v[i].getOutNeighbor(ii)){
	      loc = ii;
	      break;
	     }
 	   }
 	   if(loc >=0 && u_ == v[i].getOutNeighbor(loc)){// edge already present
 	    for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
	         edges[2*j] = v[i].getOutNeighbor(k); // copy the out going edges
	         edges[2*j+1] = v[i].getOutWeight(k);
	        }
	    vnew[i].setOutNeighbors((edges+ 2*(j - v[i].getOutDegree()))); // copy the out neighbors   
	    vnew[i].setOutDegree(outdegreeu); // copy the out degree
	    m = m - 1; // decrease size of m by one
	  }
	  else if(loc == -1){ ///new edge loc is at the end
	   for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
 	     edges[2*j] = v[i].getOutNeighbor(k); // copy the out going edges
 	     edges[2*j+1] = v[i].getOutWeight(k);
	  }
	  edges[2*j] = u_;// set the new edge u_ at the end
	  edges[2*j+1] = wt;//rand()%10+1;// set the new weight
	  
	  vnew[i].setOutDegree(outdegreeu + 1); // increase the out degree by one  
	  vnew[i].setOutNeighbors((edges+ 2*(j - v[i].getOutDegree()))); // copy the out neighbors
	  }
	  else{ 
	   for(uintT k=0; k<v[i].getOutDegree()+1; j++) { 
	     if( k == loc){
 	      edges[2*j] = u_; // set the new edge u_
 	      edges[2*j+1] = wt;//rand()%10+1;// set the new weight
 	      loc = -1; // reset loc
 	      }
 	     else{
 	      edges[2*j] = v[i].getOutNeighbor(k); // copy the out going edges
 	      edges[2*j+1] = v[i].getOutWeight(k++);
 	      }
	    }
	  vnew[i].setOutDegree(outdegreeu + 1); // increase the out degree by one  
	  vnew[i].setOutNeighbors((edges+ 2*(j - v[i].getOutDegree())-2)); // copy the out neighbors
	
	 }
	} // end of if i == u        
	else{
	  
	  for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
 	     edges[2*j] = v[i].getOutNeighbor(k); // copy the out going edges
 	     edges[2*j+1] = v[i].getOutWeight(k);
	  }
	  vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex out degree
	  vnew[i].setOutNeighbors((edges+ 2*(j - v[i].getOutDegree()))); // copy the out neighbors
         
      }
   }
 GA.del();
 GA = wghGraph<vertex>(vnew,n,m,edges); 
 } // 
 else{
     m = m + 1;
     vnew = newA(vertex,n);
     edges = newA(intE,2*m);
     uintT j = 0;
     for (uintT i = 0; i < n; i++){
      if( i == u){ // addE  is equal to v[i], set its outgoing edges
 	     edges[2*j] = u_; // set the edge value
 	     edges[2*j+1] = wt;//rand()%10+1;// set the new weight
 	     vnew[i].setOutDegree(v[i].getOutDegree()+1); // copy each vertex out degree
	     vnew[i].setOutNeighbors((edges+ 2*(j - v[i].getOutDegree()))); // copy the out neighbors
	     j++;
	  }
	  else{
	   for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
 	     edges[2*j] = v[i].getOutNeighbor(k); // copy the out going edges
 	     edges[2*j+1] = v[i].getOutWeight(k);
	  }
	  vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex out degree
	  vnew[i].setOutNeighbors((edges+ 2*(j - v[i].getOutDegree()))); // copy the out neighbors

      
    }
    
 }
  GA.del();
  GA = wghGraph<vertex>(vnew,n,m,edges);
 }
 opsUpdate++;
 ops++;
} //End putE operation

else if(optype == REME){
          cinops>>u>>u_; 
   if(u >= n || u_ >=n)
            goto last;
   vertex * v = GA.V;
   uintT outdegreeu = v[u].getOutDegree();
   
   if(outdegreeu != 0){
     vnew = newA(vertex,n);
     edges = newA(intE,2*m);
     uintT j = 0;
     for (uintT i = 0; i < n; i++){
	if( i == u){ // addE  is equal to v[i], set its outgoing edges
	  int loc = -1;
	  for(uintT ii=0; ii<v[i].getOutDegree(); ii++) { // find the location
	     if( u_ <= v[i].getOutNeighbor(ii)){
	      loc = ii;
	      break;
	     }
 	   }
 	   if(loc >=0 && u_ == v[i].getOutNeighbor(loc)){// edge  present
 	   for(uintT k=0; k<v[i].getOutDegree(); k++ ) { 
	     if( k == loc){
	       loc = -1;
 	      continue;//don't copy the edge, just ingore
 	      }
 	     else{
 	      edges[2*j] = v[i].getOutNeighbor(k); // copy the out going edges
 	      edges[2*j+1] = v[i].getOutWeight(k);
 	      j++;
 	      }
	    }
	  vnew[i].setOutDegree(v[i].getOutDegree() - 1); // increase the out degree by one  
	  vnew[i].setOutNeighbors((edges+ 2*(j - v[i].getOutDegree()+1))); // copy the out neighbors
	  m = m - 1;
	  }
 	   else{// edge not present
 	        for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
	         edges[2*j] = v[i].getOutNeighbor(k); // copy the out going edges
	         edges[2*j+1] = v[i].getOutWeight(k);
	        }
	    vnew[i].setOutNeighbors((edges+ 2*(j - v[i].getOutDegree()))); // copy the out neighbors   
	    vnew[i].setOutDegree(v[i].getOutDegree()); // copy the out degree
	    
	  }
	} // end of if i == u        
	else{
	  
	  for(uintT k=0; k<v[i].getOutDegree(); k++, j++) { 
 	     edges[2*j] = v[i].getOutNeighbor(k); // copy the out going edges
 	     edges[2*j+1] = v[i].getOutWeight(k);
	  }
	  vnew[i].setOutDegree(v[i].getOutDegree()); // copy each vertex out degree
	  vnew[i].setOutNeighbors((edges+ 2*(j - v[i].getOutDegree()))); // copy the out neighbors
         
      }
   }
  GA.del();
  GA = wghGraph<vertex>(vnew,n,m,edges);
 } 
 opsUpdate++;
 ops++;
} //End RemoveV operation



 last:
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
  long rounds = P.getOptionLongValue("-rounds",3);
  long nops=P.getOptionLongValue("-n",10);
  string opsfile=P.getOperatioFile("-f");
  //int nth = P.getOptionValue1("-nt");
  // string foutname = P.getOperatioFile("-fout");
  //cout << "Number of Threads: " << nth << endl;
  int wp = nops *(0.05);
   wghGraph<symmetricWghVertex> WG = 
   readWghGraph<symmetricWghVertex>(iFile,symmetric,binary);
   BellmanFord((intT)start,WG,nops,opsfile);
}   
  
