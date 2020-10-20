/* -*- mode: C; mode: folding; fill-column: 70; -*- */
#define _XOPEN_SOURCE 600
#define _LARGEFILE64_SOURCE 1
#define _FILE_OFFSET_BITS 64

#if defined(_OPENMP)
#include "omp.h"
#endif
#include<unistd.h>
#include "stinger_core/stinger_atomics.h"
#include "stinger_core/stinger.h"
#include "stinger_core/xmalloc.h"

#include "stinger_utils/stinger_utils.h"
#include "stinger_utils/timer.h"
//#include <thread>
#include "bfs.h"
//#include <chrono.h>
#include <unistd.h>
#include<string.h>
#include <sys/resource.h>

time_t start1,end1;
//atomic<long> vertexID;
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

struct timeval tv1, tv2;
TIME_DIFF * difference;


#define INC 1

const int64_t endian_check = 0x1234ABCD;
int64_t stinger_remove_vertex1 (struct stinger *G, int64_t vtx_id);
int stinger_save_to_file1 (struct stinger * S, uint64_t maxVtx, const char * stingerfile);
int stinger_open_from_file1 (const char * stingerfile, struct stinger * S, uint64_t * maxVtx);
int64_t stinger_insert_vertex1 (struct stinger *G, uint64_t maxVtx, int64_t vtx_id);

static int64_t nv, ne, naction;
static int64_t *restrict off;
static int64_t *restrict from;
static int64_t *restrict ind;
static int64_t *restrict weight;
static int64_t *restrict action;

/* handles for I/O memory */
static int64_t *restrict graphmem;
static int64_t *restrict actionmem;

static char *initial_graph_name = INITIAL_GRAPH_NAME_DEFAULT;
static char *action_stream_name = ACTION_STREAM_NAME_DEFAULT;
static char *dataset_file_name = NULL;

static long batch_size = BATCH_SIZE_DEFAULT;
static long nbatch = NBATCH_DEFAULT;

static struct stinger *S, *S1, *S2;

static double * update_time_trace;

int
main (const int argc, char *argv[])
{
  parse_args (argc, argv, &initial_graph_name, &action_stream_name, 
              &batch_size, &nbatch);
  STATS_INIT ();

  load_graph_and_action_stream (initial_graph_name, &nv, &ne,
                                (int64_t **) & off, (int64_t **) & ind,
                                (int64_t **) & weight,
                                (int64_t **) & graphmem, action_stream_name,
                                &naction, (int64_t **) & action,
                                (int64_t **) & actionmem);

  print_initial_graph_stats (nv, ne, batch_size, nbatch, naction);
  BATCH_SIZE_CHECK ();

//start: new added 1 
#if defined(_OPENMP)
  OMP("omp parallel")
  {
  OMP("omp master")
  PRINT_STAT_INT64 ("num_threads", (long int) omp_get_num_threads());
  }
#endif

  update_time_trace = xmalloc (nbatch * sizeof(*update_time_trace));
 //end: new added 1 


  /* Convert to STINGER */
  tic ();
  S = stinger_new ();
  //S1 = stinger_new ();
  //printf("\noffsets:");
  //for(int64_t i=0; i<nv; i++){
  //   printf(" %lld",(long long)off[i]);
  //}
  stinger_set_initial_edges (S, nv, 0, off, ind, weight, NULL, NULL, 0);
  PRINT_STAT_DOUBLE ("time_stinger", toc ());
  fflush (stdout);

  int64_t numSteps = 1;
  //int64_t src_dest_pair[2] = { 124, 381 };

  //tic ();
  //int64_t size_intersect =    st_conn_stinger (S, nv, ne, src_dest_pair, 1, numSteps);
  //PRINT_STAT_DOUBLE ("time_st_conn_stinger", toc ());
  //PRINT_STAT_INT64 ("size_intersect", size_intersect);
  
 //char * fname="../dataset1/gds9802-1M-1024";
 //char * fname="../dataset1/gds9010-1M-1024";
 char *  fname = "../dataset1/gds9505-1M-1024";

 
 FILE *fp = fopen(fname, "r"); // file to read all the opeartions 
 //FILE *fp = fopen("../dataset1/gds5050-1M-1024", "r"); // file to read all the opeartions 
 long nop=0;
 int type, type1;
 int i,j,wt;
 int64_t nv1;
 tic ();
 fscanf(fp,"%d",&nop); 
 
 //FILE *fp1 = fopen("sleeptime", "r");
 //int sleeptime;
 //MTA("mta assert parallel")
 //MTA("mta block dynamic schedule")
 //OMP("omp parallel for")

gettimeofday(&tv1,NULL);
 for(int k = 0; k < nop/1000; k++) {
   
      fscanf(fp,"%d",&type);
      //fscanf(fp1,"%d",&sleeptime);
      //printf("type: %d",type);
      if(type == 0){
          fscanf(fp,"%d",&i);
          int64_t src_dest_pair[2];
          src_dest_pair[0] = i;// source vertex, bfs starts from i
          
          
          src_dest_pair[1] = 0;
          if( i >= stinger_max_nv(S)) continue;
          st_conn_stinger(S, nv, ne, src_dest_pair, 1, numSteps);
          }
     else if(type == 1){
          fscanf(fp,"%d %d %d",&i, &j, &wt);
          if( i >= stinger_max_nv(S) || j >= stinger_max_nv(S)) continue;
         //usleep(sleeptime);
          stinger_insert_edge (S, 0, i, j, wt, k+1);
          }
      else if(type == 2){
          fscanf(fp,"%d %d",&i, &j);
          if( i >= stinger_max_nv(S) || j >= stinger_max_nv(S)) continue;
          //usleep(sleeptime);
          stinger_remove_edge(S, 0, i, j); 
          }
      else if(type == 3){
          fscanf(fp,"%d",&i);
          if( i >= stinger_max_nv(S)) continue;
         //usleep(sleeptime);
          nv = nv+1;
          stinger_insert_vertex1(S, nv, i); 
          }
       else if(type == 4){
          fscanf(fp,"%d",&i);
          if( i >= stinger_max_nv(S)) continue;
          //usleep(sleeptime);
          stinger_remove_vertex1(S, (int64_t)i); 
          }        
        
     }
      gettimeofday(&tv2,NULL);
      difference = my_difftime (&tv1, &tv2);
	int dig = 1;
	int temp = difference->usecs;
	while(temp>=10)
	{	
		dig++;
		temp = temp/10;
	}
	temp =1;
	for(i=1;i<=dig;i++)
		temp = temp * 10;
	double duration = (double) difference->secs + ((double)difference->usecs / (double)temp);
	int num_cpus = omp_get_num_threads();
	//fprintf(ftime,"Thread:%d\n",num_cpus);
	//fprintf(ftime,"Time:%lf\n \n",duration);
	//fprintf(ftime,"%d,",num_cpus);
	//fprintf(ftime,"%lf\n",duration);
   PRINT_STAT_DOUBLE ("time_update_stinger", toc ());     

  stinger_free_all (S);
  free (graphmem);
  free (actionmem);
  STATS_END ();
}


void
bfs_stinger (const struct stinger *G, const int64_t nv, const int64_t ne,
             const int64_t startV,
             int64_t * marks, const int64_t numSteps,
             int64_t * Q, int64_t * QHead, int64_t * neighbors)
{
  int64_t j, k, s;
  int64_t nQ, Qnext, Qstart, Qend;
  int64_t w_start;


  MTA ("mta assert nodep")
    for (j = 0; j < nv; j++) {
      marks[j] = 0;
    }

  s = startV;
  /* Push node s onto Q and set bounds for first Q sublist */
  Q[0] = s;
  Qnext = 1;
  nQ = 1;
  QHead[0] = 0;
  QHead[1] = 1;
  marks[s] = 1;

 PushOnStack:                   /* Push nodes onto Q */

  /* Execute the nested loop for each node v on the Q AND
     for each neighbor w of v  */
  Qstart = QHead[nQ - 1];
  Qend = QHead[nQ];
  w_start = 0;

  MTA ("mta assert no dependence")
    MTA ("mta block dynamic schedule")
    for (j = Qstart; j < Qend; j++) {
      int64_t v = Q[j];
      size_t d;
      size_t deg_v = stinger_outdegree (G, v);
      int64_t my_w = stinger_int64_fetch_add (&w_start, deg_v);
      stinger_gather_typed_successors (G, 0, v, &d, &neighbors[my_w], deg_v);
      assert (d == deg_v);

      MTA ("mta assert nodep")
        for (k = 0; k < deg_v; k++) {
          int64_t d, w, l;
          w = neighbors[my_w + k];
          /* If node has not been visited, set distance and push on Q */
          if (stinger_int64_fetch_add (marks + w, 1) == 0) {
            Q[stinger_int64_fetch_add (&Qnext, 1)] = w;
          }
        }
    }

  if (Qnext != QHead[nQ] && nQ < numSteps) {
    nQ++;
    QHead[nQ] = Qnext;
    goto PushOnStack;
  }
}


int64_t
st_conn_stinger (const struct stinger *G, const int64_t nv, const int64_t ne,
                 const int64_t * sources, const int64_t num,
                 const int64_t numSteps)
{
  int64_t k, x;

  int64_t *Q_big = (int64_t *) xmalloc (INC * nv * sizeof (int64_t));
  int64_t *marks_s_big = (int64_t *) xmalloc (INC * nv * sizeof (int64_t));
  int64_t *marks_t_big = (int64_t *) xmalloc (INC * nv * sizeof (int64_t));
  int64_t *QHead_big =
    (int64_t *) xmalloc (INC * 2 * numSteps * sizeof (int64_t));
  int64_t *neighbors_big = (int64_t *) xmalloc (INC * ne * sizeof (int64_t));

  int64_t count = 0;

  k = 0;
  x = 0;

  OMP ("omp parallel for")
    MTA ("mta assert parallel")
    MTA ("mta loop future")
    MTA ("mta assert nodep")
    MTA ("mta assert no alias")
    for (x = 0; x < INC; x++) {
      int64_t *Q = Q_big + x * nv;
      int64_t *marks_s = marks_s_big + x * nv;
      int64_t *marks_t = marks_t_big + x * nv;
      int64_t *QHead = QHead_big + x * 2 * numSteps;
      int64_t *neighbors = neighbors_big + x * ne;

      for (int64_t claimedk = stinger_int64_fetch_add (&k, 2);
           claimedk < 2 * num; claimedk = stinger_int64_fetch_add (&k, 2)) {
        int64_t s = sources[claimedk];
        int64_t deg_s = stinger_outdegree (G, s);
        int64_t t = sources[claimedk + 1];
        int64_t deg_t = stinger_outdegree (G, t);

        if (deg_s == 0 || deg_t == 0) {
          stinger_int64_fetch_add (&count, 1);
        } else {
          bfs_stinger (G, nv, ne, s, marks_s, numSteps, Q, QHead, neighbors);
          //bfs_stinger (G, nv, ne, t, marks_t, numSteps, Q, QHead, neighbors);
          int64_t local_count = 0;

          MTA ("mta assert nodep")
            for (int64_t j = 0; j < nv; j++) {
              if (marks_s[j] && marks_t[j])
                stinger_int64_fetch_add (&local_count, 1);
            }

          if (local_count == 0)
            stinger_int64_fetch_add (&count, 1);
        }
      }
    }

  free (neighbors_big);
  free (QHead_big);
  free (marks_t_big);
  free (marks_s_big);
  free (Q_big);

  return count;

}


int64_t
st_conn_stinger_source (const struct stinger * G, const int64_t nv,
                        const int64_t ne, const int64_t from,
                        const int64_t * sources, const int64_t num,
                        const int64_t numSteps)
{
  int64_t k;

  int64_t *Q = (int64_t *) xmalloc (nv * sizeof (int64_t));
  int64_t *marks_s = (int64_t *) xmalloc (nv * sizeof (int64_t));
  int64_t *marks_t = (int64_t *) xmalloc (nv * sizeof (int64_t));
  int64_t *QHead = (int64_t *) xmalloc (2 * numSteps * sizeof (int64_t));
  int64_t *neighbors = (int64_t *) xmalloc (ne * sizeof (int64_t));

  int64_t count = 0;

  int64_t deg_s = stinger_outdegree (G, from);
  if (deg_s == 0) {
    free (neighbors);
    free (QHead);
    free (marks_t);
    free (marks_s);
    free (Q);
    return num;
  }

  bfs_stinger (G, nv, ne, from, marks_s, numSteps, Q, QHead, neighbors);

  for (k = 0; k < num; k++) {
    int64_t t = sources[k];
    int64_t deg_t = stinger_outdegree (G, t);

    if (deg_t == 0) {
      stinger_int64_fetch_add (&count, 1);
    } else {
      bfs_stinger (G, nv, ne, t, marks_t, numSteps, Q, QHead, neighbors);

      int64_t local_count = 0;

      MTA ("mta assert nodep")
        for (int64_t j = 0; j < nv; j++) {
          if (marks_s[j] && marks_t[j])
            stinger_int64_fetch_add (&local_count, 1);
        }

      if (local_count == 0)
        stinger_int64_fetch_add (&count, 1);
    }
  }

  free (neighbors);
  free (QHead);
  free (marks_t);
  free (marks_s);
  free (Q);

  return count;

}


void
bfs_csr (const int64_t nv, const int64_t ne, const int64_t * off,
         const int64_t * ind, const int64_t startV, int64_t * marks,
         const int64_t numSteps, int64_t * Q, int64_t * QHead)
{
  int64_t j, k, s;
  int64_t nQ, Qnext, Qstart, Qend;


  MTA ("mta assert nodep")
    for (j = 0; j < nv; j++) {
      marks[j] = 0;
    }

  s = startV;
  /* Push node s onto Q and set bounds for first Q sublist */
  Q[0] = s;
  Qnext = 1;
  nQ = 1;
  QHead[0] = 0;
  QHead[1] = 1;
  marks[s] = 1;

 PushOnStack:                   /* Push nodes onto Q */

  /* Execute the nested loop for each node v on the Q AND
     for each neighbor w of v  */
  Qstart = QHead[nQ - 1];
  Qend = QHead[nQ];

  MTA ("mta assert no dependence")
    MTA ("mta block dynamic schedule")
    for (j = Qstart; j < Qend; j++) {
      int64_t v = Q[j];
      int64_t myStart = off[v];
      int64_t myEnd = off[v + 1];

      MTA ("mta assert nodep")
        for (k = myStart; k < myEnd; k++) {
          int64_t d, w, l;
          w = ind[k];
          /* If node has not been visited, set distance and push on Q */
          if (stinger_int64_fetch_add (marks + w, 1) == 0) {
            Q[stinger_int64_fetch_add (&Qnext, 1)] = w;
          }
        }
    }

  if (Qnext != QHead[nQ] && nQ < numSteps) {
    nQ++;
    QHead[nQ] = Qnext;
    goto PushOnStack;
  }
}


int64_t
st_conn_csr (const int64_t nv, const int64_t ne, const int64_t * off,
             const int64_t * ind, const int64_t * sources, const int64_t num,
             const int64_t numSteps)
{
  int64_t k, x;

  int64_t *Q_big = (int64_t *) xmalloc (INC * nv * sizeof (int64_t));
  int64_t *marks_s_big = (int64_t *) xmalloc (INC * nv * sizeof (int64_t));
  int64_t *marks_t_big = (int64_t *) xmalloc (INC * nv * sizeof (int64_t));
  int64_t *QHead_big =
    (int64_t *) xmalloc (INC * 2 * numSteps * sizeof (int64_t));
  int64_t *neighbors_big = (int64_t *) xmalloc (INC * ne * sizeof (int64_t));

  int64_t count = 0;

  k = 0;

  MTA ("mta assert parallel")
    MTA ("mta loop future")
    for (x = 0; x < INC; x++) {
      int64_t *Q = Q_big + x * nv;
      int64_t *marks_s = marks_s_big + x * nv;
      int64_t *marks_t = marks_t_big + x * nv;
      int64_t *QHead = QHead_big + x * 2 * numSteps;
      int64_t *neighbors = neighbors_big + x * ne;

      for (int64_t claimedk = stinger_int64_fetch_add (&k, 2);
           claimedk < 2 * num; claimedk = stinger_int64_fetch_add (&k, 2)) {

        int64_t s = sources[claimedk];
        int64_t t = sources[claimedk + 1];
        bfs_csr (nv, ne, off, ind, s, marks_s, numSteps, Q, QHead);
        bfs_csr (nv, ne, off, ind, t, marks_t, numSteps, Q, QHead);

        int64_t local_count = 0;

        MTA ("mta assert nodep")
          for (int64_t j = 0; j < nv; j++) {
            if (marks_s[j] && marks_t[j])
              stinger_int64_fetch_add (&local_count, 1);
          }

        if (local_count == 0)
          stinger_int64_fetch_add (&count, 1);
      }
    }

  free (neighbors_big);
  free (QHead_big);
  free (marks_t_big);
  free (marks_s_big);
  free (Q_big);

  return count;

}

int64_t
stinger_insert_vertex1 (struct stinger *G, uint64_t maxVtx, int64_t vtx_id)
{
  #define tdeg(X,Y) offsets[((X) * (maxVtx+2)) + (Y+2)]

  /* TODO TODO TODO fix this function and the load function to work wihtout static sizes */

  MAP_STING(S);
  uint64_t * restrict offsets = xcalloc((maxVtx+2) * (S->max_netypes), sizeof(uint64_t));

  for(int64_t type = 0; type < (S->max_netypes); type++) {
    struct stinger_eb * local_ebpool = ebpool->ebpool;
    OMP("omp parallel for")
    MTA("mta assert parallel")
    for(uint64_t block = 0; block < ETA(S,type)->high; block++) {
      struct stinger_eb * cureb = local_ebpool + ETA(S, type)->blocks[block];
      int64_t num = cureb->numEdges;
      if (num) {
	stinger_int64_fetch_add(&tdeg(type,cureb->vertexID), num);
      }
    }
  }

  for(int64_t type = 0; type < (S->max_netypes); type++) {
    prefix_sum(maxVtx+2, offsets + (type * (maxVtx+2)));
  }

  uint64_t * restrict type_offsets = xmalloc(((S->max_netypes)+1) * sizeof(uint64_t));

  type_offsets[0] = 0;
  for(int64_t type = 0; type < (S->max_netypes); type++) {
    type_offsets[type+1] = offsets[type * (maxVtx + 2) + (maxVtx + 1)] + type_offsets[type];
  }

  int64_t * restrict ind = xmalloc(4 * type_offsets[(S->max_netypes)] * sizeof(int64_t));
  int64_t * restrict weight = ind + type_offsets[(S->max_netypes)];
  int64_t * restrict timefirst = weight + type_offsets[(S->max_netypes)];
  int64_t * restrict timerecent = timefirst + type_offsets[(S->max_netypes)];

#undef tdeg
#define tdeg(X,Y) offsets[((X) * (maxVtx+2)) + (Y+1)]

  for(int64_t type = 0; type < (S->max_netypes); type++) {
    struct stinger_eb * local_ebpool = ebpool->ebpool;
    OMP("omp parallel for")
    MTA("mta assert parallel")
    for(uint64_t block = 0; block < ETA(S,type)->high; block++) {
      struct stinger_eb * cureb = local_ebpool + ETA(S,type)->blocks[block];
      int64_t num = cureb->numEdges;
      if (num) {
	int64_t my_off = stinger_int64_fetch_add(&tdeg(type,cureb->vertexID),num) + type_offsets[type];
	int64_t stop = my_off + num;
	for(uint64_t k = 0; k < stinger_eb_high(cureb) && my_off < stop; k++) {
	  if(!stinger_eb_is_blank(cureb, k)) {
	    ind[my_off] = stinger_eb_adjvtx(cureb,k);
	    weight[my_off] = stinger_eb_weight(cureb,k);
	    timefirst[my_off] = stinger_eb_first_ts(cureb,k);
	    timerecent[my_off] = stinger_eb_ts(cureb,k);
	    my_off++;
	  }
	}
      }
    }
  }
#undef tdeg




  int64_t etypes = (S->max_netypes);
  int64_t written = 0;
  

  //written += fwrite(&endian_check, sizeof(int64_t), 1, fp);
  //written += fwrite(&maxVtx, sizeof(int64_t), 1, fp);
  //written += fwrite(&etypes, sizeof(int64_t), 1, fp);
  //written += fwrite(type_offsets, sizeof(int64_t), etypes+1, fp);
  //written += fwrite(offsets, sizeof(int64_t), (maxVtx+2) * etypes, fp);
  //written += fwrite(ind, sizeof(int64_t), 4 * type_offsets[etypes], fp);
//printf("\nendian_check:%lld, maxVtx:%lld, etypes:%lld \n",(long long)endian_check,(long long)maxVtx, (long long)etypes);
  //int64_t * restrict ind = xmalloc(4 * type_offsets[etypes] * sizeof(int64_t));
  //int64_t * restrict weight = ind + type_offsets[etypes];
  //int64_t * restrict timefirst = weight + type_offsets[etypes];
  //int64_t * restrict timerecent = timefirst + type_offsets[etypes];
  //result = fread(ind, sizeof(int64_t), 4 * type_offsets[etypes], fp);
  
 // printf("\nind(%lld) :", 4 * type_offsets[etypes]);
  stinger_free_all (S);
  
   S = stinger_new ();
   
  for(int64_t type = 0; type < etypes; type++) {
    stinger_set_initial_edges(S, maxVtx, type,
	offsets + type * (maxVtx+2),
	ind + type_offsets[type],
	weight + type_offsets[type],
	timerecent + type_offsets[type],
	timefirst + type_offsets[type],
	0);
  }
  
   //printf("\n After offsets(%lld) :", (maxVtx+2) * etypes);
  //for(int64_t i=0; i<(maxVtx+2) * etypes; i++){
  //   printf(" %lld",(long long)offsets[i]);
  //}
}



/** @brief Checkpoint a STINGER data structure to disk.
 *  Format (64-bit words):
 *    - endianness check
 *    - max NV (max vertex ID + 1)
 *    - number of edge types
 *    - type offsets (length = etypes + 1) offsets of the beginning on each type in the CSR structure
 *    - offsets (length = etypes * (maxNV+2)) concatenated CSR offsets, one NV+2-array per edge type
 *    - ind/adj (legnth = type_offsets[etypes] = total number of edges) concatenated CSR adjacency arrays,
 *		      offsets into this are type_offsets[etype] + offets[etype*(maxNV+2) + sourcevertex]
 *    - weight (legnth = type_offsets[etypes] = total number of edges) concatenated csr weight arrays,
 *		      offsets into this are type_offsets[etype] + offets[etype*(maxnv+2) + sourcevertex]
 *    - time first (legnth = type_offsets[etypes] = total number of edges) concatenated csr timefirst arrays,
 *		      offsets into this are type_offsets[etype] + offets[etype*(maxnv+2) + sourcevertex]
 *    - time recent (legnth = type_offsets[etypes] = total number of edges) concatenated csr time recent arrays,
 *		      offsets into this are type_offsets[etype] + offets[etype*(maxnv+2) + sourcevertex]
 * @param S A pointer to the Stinger to be saved.
 * @param maxVtx The maximum vertex ID + 1.
 * @param stingerfile The path and name of the output file.
 * @return 0 on success, -1 on failure.
 */
int
stinger_save_to_file1 (struct stinger * S, uint64_t maxVtx, const char * stingerfile)
{
#define tdeg(X,Y) offsets[((X) * (maxVtx+2)) + (Y+2)]

  /* TODO TODO TODO fix this function and the load function to work wihtout static sizes */

  MAP_STING(S);
  uint64_t * restrict offsets = xcalloc((maxVtx+2) * (S->max_netypes), sizeof(uint64_t));

  for(int64_t type = 0; type < (S->max_netypes); type++) {
    struct stinger_eb * local_ebpool = ebpool->ebpool;
    OMP("omp parallel for")
    MTA("mta assert parallel")
    for(uint64_t block = 0; block < ETA(S,type)->high; block++) {
      struct stinger_eb * cureb = local_ebpool + ETA(S, type)->blocks[block];
      int64_t num = cureb->numEdges;
      if (num) {
	stinger_int64_fetch_add(&tdeg(type,cureb->vertexID), num);
      }
    }
  }

  for(int64_t type = 0; type < (S->max_netypes); type++) {
    prefix_sum(maxVtx+2, offsets + (type * (maxVtx+2)));
  }

  uint64_t * restrict type_offsets = xmalloc(((S->max_netypes)+1) * sizeof(uint64_t));

  type_offsets[0] = 0;
  for(int64_t type = 0; type < (S->max_netypes); type++) {
    type_offsets[type+1] = offsets[type * (maxVtx + 2) + (maxVtx + 1)] + type_offsets[type];
  }

  int64_t * restrict ind = xmalloc(4 * type_offsets[(S->max_netypes)] * sizeof(int64_t));
  int64_t * restrict weight = ind + type_offsets[(S->max_netypes)];
  int64_t * restrict timefirst = weight + type_offsets[(S->max_netypes)];
  int64_t * restrict timerecent = timefirst + type_offsets[(S->max_netypes)];

#undef tdeg
#define tdeg(X,Y) offsets[((X) * (maxVtx+2)) + (Y+1)]

  for(int64_t type = 0; type < (S->max_netypes); type++) {
    struct stinger_eb * local_ebpool = ebpool->ebpool;
    OMP("omp parallel for")
    MTA("mta assert parallel")
    for(uint64_t block = 0; block < ETA(S,type)->high; block++) {
      struct stinger_eb * cureb = local_ebpool + ETA(S,type)->blocks[block];
      int64_t num = cureb->numEdges;
      if (num) {
	int64_t my_off = stinger_int64_fetch_add(&tdeg(type,cureb->vertexID),num) + type_offsets[type];
	int64_t stop = my_off + num;
	for(uint64_t k = 0; k < stinger_eb_high(cureb) && my_off < stop; k++) {
	  if(!stinger_eb_is_blank(cureb, k)) {
	    ind[my_off] = stinger_eb_adjvtx(cureb,k);
	    weight[my_off] = stinger_eb_weight(cureb,k);
	    timefirst[my_off] = stinger_eb_first_ts(cureb,k);
	    timerecent[my_off] = stinger_eb_ts(cureb,k);
	    my_off++;
	  }
	}
      }
    }
  }
#undef tdeg

#if !defined(__MTA__)
  FILE * fp = fopen(stingerfile, "wb");

  if(!fp) {
    fprintf (stderr, "%s %d: Can't open file \"%s\" for writing\n", __func__, __LINE__, stingerfile);
    free(offsets); free(type_offsets); free(ind);
    return -1;
  }

  int64_t etypes = (S->max_netypes);
  int64_t written = 0;

  written += fwrite(&endian_check, sizeof(int64_t), 1, fp);
  written += fwrite(&maxVtx, sizeof(int64_t), 1, fp);
  written += fwrite(&etypes, sizeof(int64_t), 1, fp);
  written += fwrite(type_offsets, sizeof(int64_t), etypes+1, fp);
  written += fwrite(offsets, sizeof(int64_t), (maxVtx+2) * etypes, fp);
  written += fwrite(ind, sizeof(int64_t), 4 * type_offsets[etypes], fp);
printf("\nendian_check:%lld, maxVtx:%lld, etypes:%lld \n",(long long)endian_check,(long long)maxVtx, (long long)etypes);
 
 printf("\ntype_offsets:");
  for(int64_t i=0; i<etypes; i++){
     printf(" %lld",(long long)type_offsets[i]);
  }
   printf("\noffsets(%lld) :", (maxVtx+2) * etypes);
  for(int64_t i=0; i<(maxVtx+2) * etypes; i++){
     printf(" %lld",(long long)offsets[i]);
  }
  /*
  printf("\nind(%lld) :", 4 * type_offsets[etypes]);
  for(int64_t i=0; i<4 * type_offsets[etypes]; i++){
     printf(" %lld",(long long)ind[i]);
     if(i%type_offsets[etypes] == 0) printf("\n");
  }
  */
  //printf("\n vdata(%lld):",maxVtx);
  for(uint64_t v = 0; v < maxVtx; v++) {
    int64_t vdata[2];
    vdata[0] = stinger_vtype_get(S,v);
    vdata[1] = stinger_vweight_get(S,v);
    fwrite(vdata, sizeof(int64_t), 2, fp);
    //printf(" [%lld,%lld]",(long long)vdata[0], (long long)vdata[1]);
  }

  stinger_names_save(stinger_physmap_get(S), fp);
  stinger_names_save(stinger_vtype_names_get(S), fp);
  stinger_names_save(stinger_etype_names_get(S), fp);

  fclose(fp);
#else
  xmt_luc_io_init();
  int64_t etypes = (S->max_netypes);
  int64_t written = 0;
  size_t sz = (3 + etypes + 1 + (maxVtx+2) * etypes + 4 * type_offsets[etypes]) * sizeof(int64_t);
  int64_t * xmt_buf = xmalloc (sz);
  xmt_buf[0] = endian_check;
  xmt_buf[1] = maxVtx;
  xmt_buf[2] = etypes;
  for (int64_t i = 0; i < etypes+1; i++) {
    xmt_buf[3+i] = type_offsets[i];
  }
  for (int64_t i = 0; i < (maxVtx+2) * etypes; i++) {
    xmt_buf[3+etypes+1+i] = offsets[i];
  }
  for (int64_t i = 0; i < 4 * type_offsets[etypes]; i++) {
    xmt_buf[3+etypes+1+(maxVtx+2)*etypes+i] = ind[i];
  }
  xmt_luc_snapout(stingerfile, xmt_buf, sz);
  written = sz;
  free(xmt_buf);
#endif

  if(written != (3 + etypes + 1 + (maxVtx+2) * etypes + 4 * type_offsets[etypes])) {
    free(offsets); free(type_offsets); free(ind);
    return -1;
  } else {
    free(offsets); free(type_offsets); free(ind);
    return 0;
  }
}

/** @brief Restores a STINGER checkpoint from disk.
 * @param stingerfile The path and name of the input file.
 * @param S A pointer to an empty STINGER structure.
 * @param maxVtx Output pointer for the the maximum vertex ID + 1.
 * @return 0 on success, -1 on failure.
*/
int
stinger_open_from_file1 (const char * stingerfile, struct stinger * S, uint64_t * maxVtx)
{
#if !defined(__MTA__)
  FILE * fp = fopen(stingerfile, "rb");

  if(!fp) {
    fprintf (stderr, "%s %d: Can't open file \"%s\" for reading\n", __func__, __LINE__, stingerfile);
    return -1;
  }
#endif

  if (!S) {
    return -1;
  }

  int64_t local_endian;
  int64_t etypes = 0;
#if !defined(__MTA__)
  int result = 0;
  result += fread(&local_endian, sizeof(int64_t), 1, fp);
  result += fread(maxVtx, sizeof(int64_t), 1, fp);
  result += fread(&etypes, sizeof(int64_t), 1, fp);
printf("\nstinger_open_from_file1 local_endian:%lld, maxVtx:%lld, etypes:%lld \n",(long long)local_endian,(long long)maxVtx, (long long)etypes);
  if(result != 3) {
    fprintf (stderr, "%s %d: Fread of file \"%s\" failed.\n", __func__, __LINE__, stingerfile);
    return -1;
  }
#else
  xmt_luc_io_init();
  size_t sz;
  xmt_luc_stat(stingerfile, &sz);
  int64_t * xmt_buf = (int64_t *) xmalloc (sz);
  xmt_luc_snapin(stingerfile, xmt_buf, sz);
  local_endian = xmt_buf[0];
  *maxVtx = xmt_buf[1];
  etypes = xmt_buf[2];
#endif
//int ap;scanf("%d",&ap);
  if(local_endian != endian_check) {
    *maxVtx = bs64(*maxVtx);
    etypes = bs64(etypes);
  }

  if(*maxVtx > (S->max_nv)) {
    fprintf (stderr, "%s %d: Vertices in file \"%s\" larger than the maximum number of vertices.\n", __func__, __LINE__, stingerfile);
    return -1;
  }
  if(etypes > (S->max_netypes)) {
    fprintf (stderr, "%s %d: Edge types in file \"%s\" larger than the maximum number of edge types.\n", __func__, __LINE__, stingerfile);
    return -1;
  }

#if !defined(__MTA__)
  uint64_t * restrict offsets = xcalloc((*maxVtx+2) * etypes, sizeof(uint64_t));
  uint64_t * restrict type_offsets = xmalloc((etypes+1) * sizeof(uint64_t));

  result = fread(type_offsets, sizeof(int64_t), etypes+1, fp);
  
 // printf("\ntype_offsets:");
  //for(int64_t i=0; i<etypes; i++){
  //   printf(" %lld",(long long)type_offsets[i]);
 // }
 //scanf("%d",&ap);
  result += fread(offsets, sizeof(int64_t), (*maxVtx+2) * etypes, fp);
  printf("\noffsets(%lld) :", (long long)(*maxVtx+2) * etypes);
  //for(int64_t i=0; i<(*maxVtx+2) * etypes; i++){
  //   printf(" %lld",(long long)offsets[i]);
  //}
  if(result != etypes + 1 + (((*maxVtx) + 2) * etypes)) {
    fprintf (stderr, "%s %d: Fread of file \"%s\" failed. Types and type offsets failed.\n", __func__, __LINE__, stingerfile);
    return -1;
  }
#else
  uint64_t * restrict type_offsets = &xmt_buf[3];
  uint64_t * restrict offsets = type_offsets + etypes+1;
#endif

  if(local_endian != endian_check) {
    bs64_n(etypes + 1, type_offsets);
    bs64_n(*maxVtx + 2, offsets);
  }

#if !defined(__MTA__)
  int64_t * restrict ind = xmalloc(4 * type_offsets[etypes] * sizeof(int64_t));
  int64_t * restrict weight = ind + type_offsets[etypes];
  int64_t * restrict timefirst = weight + type_offsets[etypes];
  int64_t * restrict timerecent = timefirst + type_offsets[etypes];
  result = fread(ind, sizeof(int64_t), 4 * type_offsets[etypes], fp);
  
  printf("\nind(%lld) :", 4 * type_offsets[etypes]);
 // for(int64_t i=0; i<4 * type_offsets[etypes]; i++){
  //   printf(" %lld",(long long)ind[i]);
 //    if(i%type_offsets[etypes] == 0) printf("\n");
 // }

  if(result != 4 * type_offsets[etypes]) {
    fprintf (stderr, "%s %d: Fread of file \"%s\" failed. Edges and types failed.\n", __func__, __LINE__, stingerfile);
    return -1;
  }
#else
  int64_t * restrict ind = offsets + (*maxVtx+2)*etypes;
  int64_t * restrict weight = ind + type_offsets[etypes];
  int64_t * restrict timefirst = weight + type_offsets[etypes];
  int64_t * restrict timerecent = timefirst + type_offsets[etypes];
#endif

  if(local_endian != endian_check) {
    bs64_n(4 * type_offsets[etypes], ind);
  }

  for(int64_t type = 0; type < etypes; type++) {
    stinger_set_initial_edges(S, *maxVtx, type,
	offsets + type * (*maxVtx+2),
	ind + type_offsets[type],
	weight + type_offsets[type],
	timerecent + type_offsets[type],
	timefirst + type_offsets[type],
	0);
  }

  int ignore = 0;
  printf("\n vdata(%lld):",*maxVtx);
  for(uint64_t v = 0; v < *maxVtx; v++) {
    int64_t vdata[2];
    ignore = fread(vdata, sizeof(int64_t), 2, fp);
    //printf(" [%lld,%lld]",(long long)vdata[0], (long long)vdata[1]);
    stinger_vtype_set(S, v, vdata[0]);
    stinger_vweight_set(S, v, vdata[1]);
  }

  stinger_names_load(stinger_physmap_get(S), fp);
  stinger_names_load(stinger_vtype_names_get(S), fp);
  stinger_names_load(stinger_etype_names_get(S), fp);

#if !defined(__MTA__)
  fclose(fp);
  free(offsets); free(type_offsets); free(ind);
#else
  free(xmt_buf);
#endif
  return 0;
}

/** @brief Removes a vertex and all incident edges from the graph
 *
 *  Removes all edges incident to the vertex and unmaps the vertex
 *  in the physmap.  This algorithm first removes all out edges. During
 *  this pass it will assume an undirected graph and attempt to remove
 *  all reverse edges if they exist.  In the case of an undirected graph
 *  the algorithm will be complete and will terminate quickly.
 *
 *  If, however, the graph is a directed graph and there are incident edges
 *  still pointing to the vertex (stinger_indegree_get for the vtx > 0).  Then
 *  all edges are traversed looking for incident edges until the indegree is 0.
 *
 *  As a final step, the vertex is removed from the physmap, freeing this vertex
 *  to be re-assigned to a new vertex string.
 *
 *  @param G The STINGER data structure
 *  @param vtx_id The vertex to remove
 *  @return 0 if deletion is successful, -1 if it fails
 */
int64_t
stinger_remove_vertex1 (struct stinger *G, int64_t vtx_id)
{
  // Remove out edges
  STINGER_PARALLEL_FORALL_EDGES_OF_VTX_BEGIN(G,vtx_id) {
    stinger_remove_edge_pair(G,STINGER_EDGE_TYPE,STINGER_EDGE_SOURCE,STINGER_EDGE_DEST);
  } STINGER_PARALLEL_FORALL_EDGES_OF_VTX_END();

  // Remove remaining in edges
  if (stinger_indegree_get(G,vtx_id) > 0) {
    STINGER_FORALL_EDGES_OF_ALL_TYPES_BEGIN(G) {
      int64_t from = STINGER_EDGE_SOURCE;
      int64_t to = STINGER_EDGE_DEST;
      int64_t type = STINGER_EDGE_TYPE;
      if (to == vtx_id) {
        stinger_remove_edge(G,type,from,to);
        if (stinger_indegree_get(G,vtx_id) == 0) {
          break;
        }
      }
    } STINGER_FORALL_EDGES_OF_ALL_TYPES_END();
  }

  // Remove from physmap
  return stinger_physmap_vtx_remove_id(stinger_physmap_get(G),stinger_vertices_get(G),vtx_id);
}

