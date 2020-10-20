/*
 * File:
 *  
 *
 * Author(s):
 *     Muktikanta Sa   <cs15resch11012@iith.ac.in>
 *   
 *   
 * Description:
 *   
 *
*/
#include <bits/stdc++.h> 
#include <chrono>
#include <unistd.h>
#include <sys/resource.h>
#include <thread>
#include <pthread.h>
#include"BFS-PGCn.cpp"
#define Warmup 5
#define BILLION  1000000000.0;
time_t start1,end1;
atomic<long> vertexID;
double seconds;
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
enum type {APP=0, ADDE, REME, ADDV, REMV };
struct timeval tv1, tv2, tv3;
TIME_DIFF * difference1, *difference2;
int NTHREADS;
atomic <int> ops;
atomic <int> opsUpdate;
atomic <int> opsLookup;
atomic <int> opsapp;
int ***dynamicops; // total number of ops, # rows same as # threads
double ***updateTime; // Update time with # of other operations
typedef struct infothread{
  int tid; //thread ID
  int core_id; // core id
  int cols; // column size
  int wp; // warm up value
  Hash G;
}tinfo;

long n,m;
int  nops;
int itr = 0;


void loadOps(string file, int nocols){
int i,j,k;
dynamicops = (int ***)malloc(NTHREADS*sizeof(int**));
  for (i = 0; i< NTHREADS; i++) {
     dynamicops[i] = (int **) malloc(nocols*sizeof(int *));
        for (j = 0; j < nocols; j++) {
            dynamicops[i][j] = (int *)malloc(4*sizeof(int));
       }
    }
 ifstream cinn(file);   
 long n,m;
 int u, v, v_, wt;
 int nop;
 cinn>>nop;
 int optype;
 for (i = 0; i< NTHREADS; i++) {
   for (j = 0; j < nocols; j++) {
        //for (k = 0; j < 4; k++) {
            cinn>> optype;
            switch(optype){
                   
               
                case ADDE : 
                      cinn>>u>>v>>wt;  
                      dynamicops[i][j][0] = optype;
                      dynamicops[i][j][1] = u;
                      dynamicops[i][j][2] = v;
                      dynamicops[i][j][3] = wt;
                      
                    break;
              
                case REME : 
                      cinn>>u>>v;  
                      dynamicops[i][j][0] = optype;
                      dynamicops[i][j][1] = u;
                      dynamicops[i][j][2] = v;
                      dynamicops[i][j][3] = -1;
                break;                              
                
                case APP : 
                      cinn>>u;  
                      dynamicops[i][j][0] = optype;
                      dynamicops[i][j][1] = u;
                      dynamicops[i][j][2] = -1;
                      dynamicops[i][j][3] = -1;
                break; 
                case ADDV : 
                      cinn>>u;  
                      dynamicops[i][j][0] = optype;
                      dynamicops[i][j][1] = u;
                      dynamicops[i][j][2] = -1;
                      dynamicops[i][j][3] = -1;
                break;   
                case REMV : 
                      cinn>>u;  
                      dynamicops[i][j][0] = optype;
                      dynamicops[i][j][1] = u;
                      dynamicops[i][j][2] = -1;
                      dynamicops[i][j][3] = -1;
                break;                  
                default: break;                                          
                }          
            
       }
    }
   

}

void* pthread_call(void* t)
{
        tinfo *ti=(tinfo*)t;
        int tid = ti->tid;
        int core_id = ti->core_id;
        const pthread_t pid = pthread_self();
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(core_id, &cpuset);
        int rc = pthread_setaffinity_np(pid,sizeof(cpu_set_t), &cpuset);
        Hash G1=ti->G;
        int u, v;
        int optype;
        int cols = ti->cols;
        int wp = ti->wp;
        struct timespec starttime, endtime;
        for(int i=wp; i<cols; i++){
            optype = dynamicops[tid][i][0];

            switch(optype){

                case APP : 
                      u = dynamicops[tid][i][1];             
                      G1.GetBFS(u, tid);
                      opsapp++;
                      ops++;
                      break;

                case ADDE : 
                      u = dynamicops[tid][i][1];              
                      v = dynamicops[tid][i][2];     
                      G1.PutE(u,v, tid);
                      opsUpdate++;
                      ops++;
                      break;
                                 
                 case REME : 
                      u = dynamicops[tid][i][1];              
                      v = dynamicops[tid][i][2];  
                      G1.RemoveE(u,v, tid);
                      opsUpdate++;
                      ops++;
                      break;
                 case ADDV : 
                      u = dynamicops[tid][i][1];             
                      G1.insert(u, NTHREADS, tid);
                      opsapp++;
                      ops++;
                      break;
                 case REMV : 
                      u = dynamicops[tid][i][1];             
                      G1.remove(u, NTHREADS, tid);
                      opsapp++;
                      ops++;
                      break;     
                 default: break;     
        }
      }                
      
}

void* pthread_callwp(void* t)
{
       tinfo *ti=(tinfo*)t;
        int tid = ti->tid;
        int core_id = ti->core_id;
        const pthread_t pid = pthread_self();
        cpu_set_t cpuset;
        CPU_ZERO(&cpuset);
        CPU_SET(core_id, &cpuset);
        int rc = pthread_setaffinity_np(pid,sizeof(cpu_set_t), &cpuset);
        Hash G1=ti->G;
        int u, v;
        int optype;
        int cols = ti->wp;
        for(int i=0; i<cols; i++){
            optype = dynamicops[tid][i][0];
            switch(optype){

                case APP : 
                      u = dynamicops[tid][i][1];             
                      G1.GetBFS(u, tid);
                      break;
      
                case ADDE : 
                      u = dynamicops[tid][i][1];              
                      v = dynamicops[tid][i][2];     
                      G1.PutE(u,v, tid);
                      break;
                                 
                 case REME : 
                      u = dynamicops[tid][i][1];              
                      v = dynamicops[tid][i][2];  
                      G1.RemoveE(u,v, tid);
                      break;
                case ADDV : 
                      u = dynamicops[tid][i][1];             
                      G1.insert(u, NTHREADS, tid);
                      break;
                 case REMV : 
                      u = dynamicops[tid][i][1];             
                      G1.remove(u, NTHREADS, tid);
                      break; 
                 default: break;     
        }

      }                
}

int main(int argc, char*argv[]) 
{
        Hash sg;
        
        int i;
        if(argc < 3)
        {
                cout << "Enter 3 command line arguments - #threads, init graph and ops file" << endl;
                return 0;
        }
       
       ifstream cingph(argv[3]); // Input graph
 	cingph>>n>>m;
	cingph.close();
	th = atoi(argv[2]); // threshold value
        NTHREADS = atoi(argv[1]);  // # of threads 
        nops =   atoi(argv[5]); // # of operation to be perform
        int cols_size = nops/NTHREADS;
        int wp = cols_size *(0.05);// 5% warm up
        sg.initHNode(n);
        loadOps(argv[4], cols_size);// graph ops file
        sg.initGraphFromFile(argv[3],NTHREADS,1); // load graph from file
        pthread_t *thr = new pthread_t[NTHREADS];
        
        // Make threads Joinable for sure.
        pthread_attr_t attr;
        pthread_attr_init (&attr);
        pthread_attr_setdetachstate (&attr, PTHREAD_CREATE_JOINABLE);
        
        ops = 0;
        opsUpdate = 0;
        opsLookup = 0;;
        opsapp = 0;
        gettimeofday(&tv1,NULL);
         for (i=0;i < NTHREADS;i++)
        {
              tinfo *t =(tinfo*) malloc(sizeof(tinfo));
               t->tid = i;
                if(i >= 14 && i<= 27 && NTHREADS == 28)
                  t->core_id = i + 14; 
                else
                  t->core_id = i; 
                
                t->G = sg;
                t->wp = wp;
                t->cols = cols_size;
                pthread_create(&thr[i], &attr, pthread_callwp, (void*)t);
        }
      
        for (i = 0; i < NTHREADS; i++)
        {
                pthread_join(thr[i], NULL);
        }

       gettimeofday(&tv2,NULL);
          
        for (i=0;i < NTHREADS;i++)
        {
              tinfo *t =(tinfo*) malloc(sizeof(tinfo));
                t->tid = i;
                if(i >= 14 && i<= 27 && NTHREADS == 28)
                  t->core_id = i + 14; 
                else
                  t->core_id = i; 
                
                t->G = sg;
                t->wp = wp;
                t->cols = cols_size;
                pthread_create(&thr[i], &attr, pthread_call, (void*)t);
        }
        for (i = 0; i < NTHREADS; i++)
        {
                pthread_join(thr[i], NULL);
        }
        gettimeofday(&tv3,NULL);
        delete []thr;
        difference1 = my_difftime (&tv1, &tv2);
        difference2 = my_difftime (&tv2, &tv3);
	int dig1 = 1;
	int temp1 = difference1->usecs;
	while(temp1>=10)
	{	
		dig1++;
		temp1 = temp1/10;
	}
	temp1 =1;
	for(i=1;i<=dig1;i++)
		temp1 = temp1 * 10;
	double duration1 = (double) difference1->secs + ((double)difference1->usecs / (double)temp1);
        int dig2 = 1;
	int temp2 = difference2->usecs;
	while(temp2>=10)
	{	
		dig2++;
		temp2 = temp2/10;
	}
	temp2 =1;
	for(i=1;i<=dig2;i++)
		temp2 = temp2 * 10;
	double duration2 = (double) difference2->secs + ((double)difference2->usecs / (double)temp2);
	cout<<"Number of Threads:"<<NTHREADS<<endl;
    	cout << "Graph Size V: " << n <<" and E:"<<m<<endl;
    	cout << "Execution Time : " << duration2 <<" secs."<<endl;
        return 0;
}


