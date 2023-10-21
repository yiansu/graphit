#include <iostream> 
#include <vector>
#include <algorithm>
#include "intrinsics.h"
#ifdef GEN_PYBIND_WRAPPERS
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
namespace py = pybind11;
#endif

#if defined(USE_BASELINE) || defined(USE_FORKJOIN)
#include <functional>
#include <taskparts/benchmark.hpp>
#elif defined(USE_OPENMP)
#include "utility.hpp"
#elif defined(USE_HB_COMPILER) || defined(USE_HB_MANUAL)
#include "loop_handler.hpp"
#endif

Graph edges;
int  * __restrict parent;
template <typename TO_FUNC , typename APPLY_FUNC> VertexSubset<NodeID>* edgeset_apply_pull_parallel_from_vertexset_to_filter_func_with_frontier2(Graph & g , VertexSubset<NodeID>* from_vertexset, TO_FUNC to_func, APPLY_FUNC apply_func) 
{ 
    int64_t numVertices = g.num_nodes(), numEdges = g.num_edges();
  VertexSubset<NodeID> *next_frontier = new VertexSubset<NodeID>(g.num_nodes(), 0);
  bool * next = newA(bool, g.num_nodes());
  ligra::parallel_for_lambda((int)0, (int)numVertices, [&] (int i) { next[i] = 0; });
  from_vertexset->toDense();
  ligra::parallel_for_lambda((NodeID)0, (NodeID)g.num_nodes(), [&] (NodeID d) {
    if (to_func(d)){ 
      for(NodeID s : g.in_neigh(d)){
        if (from_vertexset->bool_map_[s] ) { 
          if( apply_func ( s , d ) ) { 
            next[d] = 1; 
            if (!to_func(d)) break; 
          }
        }
      } //end of loop on in neighbors
    } //end of to filtering 
  }); //end of outer for loop
  next_frontier->num_vertices_ = sequence::sum(next, numVertices);
  next_frontier->bool_map_ = next;
  next_frontier->is_dense = true;
  return next_frontier;
} //end of edgeset apply function 
struct parent_generated_vector_op_apply_func_0
{
void operator() (NodeID v) 
  {
    parent[v] =  -(1) ;
  };
};
struct updateEdge
{
bool operator() (NodeID src, NodeID dst) 
  {
    bool output1 ;
    parent[dst] = src;
    output1 = (bool) 1;
    return output1;
  };
};
struct toFilter
{
bool operator() (NodeID v) 
  {
    bool output ;
    output = (parent[v]) == ( -(1) );
    return output;
  };
};
struct reset
{
void operator() (NodeID v) 
  {
    parent[v] =  -(1) ;
  };
};

#if defined(USE_BASELINE) || defined(TEST_CORRECTNESS)
void BFS_baseline(Graph &g, int start_vertex) {
  // ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
  //   reset()(vertexsetapply_iter);
  // });;
  // VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , (0) );
  // builtin_addVertex(frontier, start_vertex) ;
  // parent[start_vertex] = start_vertex;
  // while ( (builtin_getVertexSetSize(frontier) ) != ((0) ))
  // {
  //   VertexSubset<int> *  output = edgeset_apply_pull_parallel_from_vertexset_to_filter_func_with_frontier2(edges, frontier, toFilter(), updateEdge()); 
  //   deleteObject(frontier) ;
  //   frontier = output;
  // }
  // deleteObject(frontier) ;

  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    reset()(i);
  }
  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , (0) );
  builtin_addVertex(frontier, start_vertex) ;
  parent[start_vertex] = start_vertex;
  while ( (builtin_getVertexSetSize(frontier) ) != ((0) ))
  {
    int64_t numVertices = g.num_nodes(), numEdges = g.num_edges();
    VertexSubset<NodeID> *output = new VertexSubset<NodeID>(g.num_nodes(), 0);
    bool * next = newA(bool, g.num_nodes());
    for (uint64_t i = 0; i < numVertices; i++) {
      next[i] = 0;
    }
    frontier->toDense();
    for (uint64_t d = 0; d < g.num_nodes(); d++) {
      if (toFilter()(d)){ 
        for(NodeID s : g.in_neigh(d)){
          if (frontier->bool_map_[s] ) { 
            if( updateEdge() ( s , d ) ) { 
              next[d] = 1; 
              if (!toFilter()(d)) break; 
            }
          }
        } //end of loop on in neighbors
      } //end of to filtering 
    } //end of outer for loop
    output->num_vertices_ = sequence::sum(next, numVertices);
    output->bool_map_ = next;
    output->is_dense = true;

    deleteObject(frontier) ;
    frontier = output;
  }
  deleteObject(frontier) ;
}

#if defined(TEST_CORRECTNESS)
void test_correctness(int start_vertex) {
  // save previous output
  int  * __restrict parent_output = new int [ builtin_getVertices(edges) ];
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    parent_output[i] = parent[i];
  }
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    reset()(vertexsetapply_iter);
  });;

  // run serial and generate output
  BFS_baseline(edges, start_vertex);

  // compare previous output with serial-version output
  uint64_t num_diffs = 0;
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    if (parent_output[i] != parent[i]) {
      std::cout << " parent_output[i]=" << parent_output[i] << " parent[i]=" << parent[i] << " at i=" << i << std::endl;
      num_diffs++;
    }
  }
  if (num_diffs > 0) {
    printf("\033[0;31mINCORRECT!\033[0m");
    printf("  num_diffs = %lu\n", num_diffs);
  } else {
    printf("\033[0;32mCORRECT!\033[0m\n");
  }
}
#endif
#endif

#if defined(USE_OPENMP)
void BFS_openmp(Graph &g, int start_vertex) {
  #if defined(OMP_SCHEDULE_STATIC)
    #pragma omp parallel for schedule(static)
  #elif defined(OMP_SCHEDULE_DYNAMIC)
    #pragma omp parallel for schedule(dynamic)
  #elif defined(OMP_SCHEDULE_GUIDED)
    #pragma omp parallel for schedule(guided)
  #endif
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    reset()(i);
  }
  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , (0) );
  builtin_addVertex(frontier, start_vertex) ;
  parent[start_vertex] = start_vertex;
  while ( (builtin_getVertexSetSize(frontier) ) != ((0) ))
  {
    int64_t numVertices = g.num_nodes(), numEdges = g.num_edges();
    VertexSubset<NodeID> *output = new VertexSubset<NodeID>(g.num_nodes(), 0);
    bool * next = newA(bool, g.num_nodes());
    #if defined(OMP_SCHEDULE_STATIC)
      #pragma omp parallel for schedule(static)
    #elif defined(OMP_SCHEDULE_DYNAMIC)
      #pragma omp parallel for schedule(dynamic)
    #elif defined(OMP_SCHEDULE_GUIDED)
      #pragma omp parallel for schedule(guided)
    #endif
    for (uint64_t i = 0; i < numVertices; i++) {
      next[i] = 0;
    }
    frontier->toDense();
    #if defined(OMP_SCHEDULE_STATIC)
      #pragma omp parallel for schedule(static)
    #elif defined(OMP_SCHEDULE_DYNAMIC)
      #pragma omp parallel for schedule(dynamic)
    #elif defined(OMP_SCHEDULE_GUIDED)
      #pragma omp parallel for schedule(guided)
    #endif
    for (uint64_t d = 0; d < g.num_nodes(); d++) {
      if (toFilter()(d)){ 
        for(NodeID s : g.in_neigh(d)){
          if (frontier->bool_map_[s] ) { 
            if( updateEdge() ( s , d ) ) { 
              next[d] = 1; 
              if (!toFilter()(d)) break; 
            }
          }
        } //end of loop on in neighbors
      } //end of to filtering 
    } //end of outer for loop
    output->num_vertices_ = sequence::sum(next, numVertices);
    output->bool_map_ = next;
    output->is_dense = true;

    deleteObject(frontier) ;
    frontier = output;
  }
  deleteObject(frontier) ;
}

#elif defined(USE_HB_COMPILER)
#if defined(RUN_HEARTBEAT)
  bool run_heartbeat = true;
#else
  bool run_heartbeat = false;
#endif

void HEARTBEAT_loop0(uint64_t maxIter) {
  for (uint64_t i = 0; i < maxIter; i++) {
    reset()(i);
  }
}

void HEARTBEAT_loop1(uint64_t maxIter, bool *next) {
  for (uint64_t i = 0; i < maxIter; i++) {
    next[i] = 0;
  }
}

void HEARTBEAT_loop2(uint64_t maxIter, Graph &g, VertexSubset<int> *frontier, bool *next) {
  for (uint64_t d = 0; d < maxIter; d++) {
    if (toFilter()(d)){ 
      for(NodeID s : g.in_neigh(d)){
        if (frontier->bool_map_[s] ) { 
          if( updateEdge() ( s , d ) ) { 
            next[d] = 1; 
            if (!toFilter()(d)) break; 
          }
        }
      } //end of loop on in neighbors
    } //end of to filtering 
  } //end of outer for loop
}

void BFS_hbc(Graph &g, int start_vertex) {
  HEARTBEAT_loop0(builtin_getVertices(edges));
  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , (0) );
  builtin_addVertex(frontier, start_vertex) ;
  parent[start_vertex] = start_vertex;
  while ( (builtin_getVertexSetSize(frontier) ) != ((0) ))
  {
    int64_t numVertices = g.num_nodes(), numEdges = g.num_edges();
    VertexSubset<NodeID> *output = new VertexSubset<NodeID>(g.num_nodes(), 0);
    bool * next = newA(bool, g.num_nodes());
    HEARTBEAT_loop1(numVertices, next);
    frontier->toDense();
    HEARTBEAT_loop2(g.num_nodes(), edges, frontier, next);
    output->num_vertices_ = sequence::sum(next, numVertices);
    output->bool_map_ = next;
    output->is_dense = true;

    deleteObject(frontier) ;
    frontier = output;
  }
  deleteObject(frontier) ;
}

#endif

int main(int argc, char * argv[])
{
  auto graph_file_name = "../Twitter.el";
  if (const auto env_p = std::getenv("INPUT_GRAPH")) {
    graph_file_name = env_p;
  }
  auto runs = 1;
  if (const auto env_p = std::getenv("RUNS")) {
    if (const auto env_p_p = std::atol(env_p)) {
      runs = env_p_p;
    }
  }
  edges = builtin_loadEdgesFromFile ( graph_file_name ) ;
  int start_vertex = atoi( "100" ) ;
  parent = new int [ builtin_getVertices(edges) ];
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    parent_generated_vector_op_apply_func_0()(vertexsetapply_iter);
  });;
  for ( int trail = (0) ; trail < (runs) ; trail++ )
  {
  //   startTimer() ;
// ####################
#if defined(USE_BASELINE) || defined(USE_FORKJOIN)
  taskparts::benchmark_nativeforkjoin([&] (auto sched) {
#elif defined(USE_OPENMP)
  utility::run([&] {
#elif defined(USE_HB_COMPILER) || defined(USE_HB_MANUAL)
  run_bench([&] {
#endif

#if defined(USE_BASELINE)
    BFS_baseline(edges, start_vertex);
#elif defined(USE_OPENMP)
    BFS_openmp(edges, start_vertex);
#elif defined(USE_HB_COMPILER)
    BFS_hbc(edges, start_vertex);
#endif

#if defined(TEST_CORRECTNESS)
    test_correctness(start_vertex);
#endif

#if defined(USE_BASELINE) || defined(USE_FORKJOIN)
  }, [&] (auto sched) {
  }, [&] (auto sched) {
  });
#else
  }, [&] {}, [&] {});
#endif
// ####################
  //   float elapsed_time = stopTimer() ;
  //   std::cout << "elapsed time: "<< std::endl;
  //   std::cout << elapsed_time<< std::endl;
  }
};
#ifdef GEN_PYBIND_WRAPPERS
PYBIND11_MODULE(, m) {
}
#endif

