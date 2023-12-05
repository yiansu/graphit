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
int  * __restrict IDs;
template <typename APPLY_FUNC > VertexSubset<NodeID>* edgeset_apply_pull_parallel_deduplicatied_from_vertexset_with_frontier3(Graph & g , VertexSubset<NodeID>* from_vertexset, APPLY_FUNC apply_func) 
{ 
    int64_t numVertices = g.num_nodes(), numEdges = g.num_edges();
  VertexSubset<NodeID> *next_frontier = new VertexSubset<NodeID>(g.num_nodes(), 0);
  bool * next = newA(bool, g.num_nodes());
  ligra::parallel_for_lambda((int)0, (int)numVertices, [&] (int i) { next[i] = 0; });
  from_vertexset->toDense();
  ligra::parallel_for_lambda((NodeID)0, (NodeID)g.num_nodes(), [&] (NodeID d) {
    for(NodeID s : g.in_neigh(d)){
      if (from_vertexset->bool_map_[s] ) { 
        if( apply_func ( s , d ) ) { 
          next[d] = 1; 
        }
      }
    } //end of loop on in neighbors
  }); //end of outer for loop
  next_frontier->num_vertices_ = sequence::sum(next, numVertices);
  next_frontier->bool_map_ = next;
  next_frontier->is_dense = true;
  return next_frontier;
} //end of edgeset apply function 
struct IDs_generated_vector_op_apply_func_0
{
void operator() (NodeID v) 
  {
    IDs[v] = (1) ;
  };
};
struct updateEdge
{
bool operator() (NodeID src, NodeID dst) 
  {
    bool output2 ;
    bool IDs_trackving_var_1 = (bool) 0;
    if ( ( IDs[dst]) > ( IDs[src]) ) { 
      IDs[dst]= IDs[src]; 
      IDs_trackving_var_1 = true ; 
    } 
    output2 = IDs_trackving_var_1;
    return output2;
  };
};
struct init
{
void operator() (NodeID v) 
  {
    IDs[v] = v;
  };
};

#if defined(USE_BASELINE) || defined(TEST_CORRECTNESS)
void CC_baseline(Graph &g, int n) {
  // VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , n);
  // ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
  //   init()(vertexsetapply_iter);
  // });;
  // while ( (builtin_getVertexSetSize(frontier) ) != ((0) ))
  // {
  //   VertexSubset<int> *  output = edgeset_apply_pull_parallel_deduplicatied_from_vertexset_with_frontier3(edges, frontier, updateEdge()); 
  //   deleteObject(frontier) ;
  //   frontier = output;
  // }
  // deleteObject(frontier) ;

  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , n);
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    init()(i);
  }
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
      for (uint64_t i = g.get_in_neighbors_begin_index_(d); i < g.get_in_neighbors_end_index_(d); i++) {
        if (frontier->bool_map_[g.get_in_neighbors_()[i]] ) { 
          if( updateEdge() ( g.get_in_neighbors_()[i] , d ) ) { 
            next[d] = 1; 
          }
        }
      } //end of loop on in neighbors
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
void test_correctness() {
  // save previous output
  int  * __restrict IDs_output = new int [ builtin_getVertices(edges) ];
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    IDs_output[i] = IDs[i];
  }
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    IDs_generated_vector_op_apply_func_0()(vertexsetapply_iter);
  });;

  // run serial and generate output
  CC_baseline(edges, builtin_getVertices(edges));

  // compare previous output with serial-version output
  uint64_t num_diffs = 0;
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    if (IDs_output[i] != IDs[i]) {
      std::cout << " IDs_output[i]=" << IDs_output[i] << " IDs[i]=" << IDs[i] << " at i=" << i << std::endl;
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

#include <omp.h>

void CC_openmp(Graph &g, int n) {
#if defined(OMP_NESTED_PARALLELISM)
  omp_set_max_active_levels(2);
#endif
  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , n);
#if !defined(OMP_CHUNKSIZE)
#if defined(OMP_SCHEDULE_STATIC)
  #pragma omp parallel for schedule(static)
#elif defined(OMP_SCHEDULE_DYNAMIC)
  #pragma omp parallel for schedule(dynamic)
#elif defined(OMP_SCHEDULE_GUIDED)
  #pragma omp parallel for schedule(guided)
#endif
#else
#if defined(OMP_SCHEDULE_STATIC)
  #pragma omp parallel for schedule(static, OMP_CHUNKSIZE)
#elif defined(OMP_SCHEDULE_DYNAMIC)
  #pragma omp parallel for schedule(dynamic, OMP_CHUNKSIZE)
#elif defined(OMP_SCHEDULE_GUIDED)
  #pragma omp parallel for schedule(guided, OMP_CHUNKSIZE)
#endif
#endif
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    init()(i);
  }
  while ( (builtin_getVertexSetSize(frontier) ) != ((0) ))
  {
    int64_t numVertices = g.num_nodes(), numEdges = g.num_edges();
    VertexSubset<NodeID> *output = new VertexSubset<NodeID>(g.num_nodes(), 0);
    bool * next = newA(bool, g.num_nodes());
#if !defined(OMP_CHUNKSIZE)
#if defined(OMP_SCHEDULE_STATIC)
    #pragma omp parallel for schedule(static)
#elif defined(OMP_SCHEDULE_DYNAMIC)
    #pragma omp parallel for schedule(dynamic)
#elif defined(OMP_SCHEDULE_GUIDED)
    #pragma omp parallel for schedule(guided)
#endif
#else
#if defined(OMP_SCHEDULE_STATIC)
    #pragma omp parallel for schedule(static, OMP_CHUNKSIZE)
#elif defined(OMP_SCHEDULE_DYNAMIC)
    #pragma omp parallel for schedule(dynamic, OMP_CHUNKSIZE)
#elif defined(OMP_SCHEDULE_GUIDED)
    #pragma omp parallel for schedule(guided, OMP_CHUNKSIZE)
#endif
#endif
    for (uint64_t i = 0; i < numVertices; i++) {
      next[i] = 0;
    }
    frontier->toDense();
#if !defined(OMP_CHUNKSIZE)
#if defined(OMP_SCHEDULE_STATIC)
    #pragma omp parallel for schedule(static)
#elif defined(OMP_SCHEDULE_DYNAMIC)
    #pragma omp parallel for schedule(dynamic)
#elif defined(OMP_SCHEDULE_GUIDED)
    #pragma omp parallel for schedule(guided)
#endif
#else
#if defined(OMP_SCHEDULE_STATIC)
    #pragma omp parallel for schedule(static, OMP_CHUNKSIZE)
#elif defined(OMP_SCHEDULE_DYNAMIC)
    #pragma omp parallel for schedule(dynamic, OMP_CHUNKSIZE)
#elif defined(OMP_SCHEDULE_GUIDED)
    #pragma omp parallel for schedule(guided, OMP_CHUNKSIZE)
#endif
#endif
    for (uint64_t d = 0; d < g.num_nodes(); d++) {
#if defined(OMP_NESTED_PARALLELISM)
#if !defined(OMP_CHUNKSIZE)
#if defined(OMP_SCHEDULE_STATIC)
      #pragma omp parallel for schedule(static)
#elif defined(OMP_SCHEDULE_DYNAMIC)
      #pragma omp parallel for schedule(dynamic)
#elif defined(OMP_SCHEDULE_GUIDED)
      #pragma omp parallel for schedule(guided)
#endif
#else
#if defined(OMP_SCHEDULE_STATIC)
      #pragma omp parallel for schedule(static, OMP_CHUNKSIZE)
#elif defined(OMP_SCHEDULE_DYNAMIC)
      #pragma omp parallel for schedule(dynamic, OMP_CHUNKSIZE)
#elif defined(OMP_SCHEDULE_GUIDED)
      #pragma omp parallel for schedule(guided, OMP_CHUNKSIZE)
#endif
#endif
#endif
      for (uint64_t i = g.get_in_neighbors_begin_index_(d); i < g.get_in_neighbors_end_index_(d); i++) {
        if (frontier->bool_map_[g.get_in_neighbors_()[i]] ) { 
          if( updateEdge() ( g.get_in_neighbors_()[i] , d ) ) { 
            next[d] = 1; 
          }
        }
      } //end of loop on in neighbors
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
    init()(i);
  }
}

void HEARTBEAT_loop1(uint64_t maxIter, bool *next) {
  for (uint64_t i = 0; i < maxIter; i++) {
    next[i] = 0;
  }
}

void HEARTBEAT_loop3(uint64_t startIter, uint64_t maxIter, Graph &g, VertexSubset<NodeID> *frontier, bool *next, uint64_t d) {
  for (uint64_t i = startIter; i < maxIter; i++) {
    if (frontier->bool_map_[g.get_in_neighbors_()[i]] ) { 
      if( updateEdge() ( g.get_in_neighbors_()[i] , d ) ) { 
        next[d] = 1; 
      }
    }
  } //end of loop on in neighbors
}

void HEARTBEAT_loop2(uint64_t maxIter, Graph &g, VertexSubset<NodeID> *frontier, bool *next) {
  for (uint64_t d = 0; d < maxIter; d++) {
    HEARTBEAT_loop3(g.get_in_neighbors_begin_index_(d), g.get_in_neighbors_end_index_(d), g, frontier, next, d);
  } //end of outer for loop
}

void CC_hbc(Graph &g, int n) {
  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , n);
  HEARTBEAT_loop0(builtin_getVertices(edges));
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
  auto graph_file_name = "inputs/Twitter.el";
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
  IDs = new int [ builtin_getVertices(edges) ];
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    IDs_generated_vector_op_apply_func_0()(vertexsetapply_iter);
  });;
  int n = builtin_getVertices(edges) ;
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
    CC_baseline(edges, n);
#elif defined(USE_OPENMP)
    CC_openmp(edges, n);
#elif defined(USE_HB_COMPILER)
    CC_hbc(edges, n);
#endif

#if defined(TEST_CORRECTNESS)
    test_correctness();
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

