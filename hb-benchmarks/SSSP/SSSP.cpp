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

WGraph edges;
int  * __restrict SP;
template <typename APPLY_FUNC > VertexSubset<NodeID>* edgeset_apply_pull_parallel_weighted_deduplicatied_from_vertexset_with_frontier3(WGraph & g , VertexSubset<NodeID>* from_vertexset, APPLY_FUNC apply_func) 
{ 
    int64_t numVertices = g.num_nodes(), numEdges = g.num_edges();
  VertexSubset<NodeID> *next_frontier = new VertexSubset<NodeID>(g.num_nodes(), 0);
  bool * next = newA(bool, g.num_nodes());
  ligra::parallel_for_lambda((int)0, (int)numVertices, [&] (int i) { next[i] = 0; });
  from_vertexset->toDense();
  ligra::parallel_for_lambda((NodeID)0, (NodeID)g.num_nodes(), [&] (NodeID d) {
    for(WNode s : g.in_neigh(d)){
      if (from_vertexset->bool_map_[s.v] ) { 
        if( apply_func ( s.v , d, s.w ) ) { 
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
struct SP_generated_vector_op_apply_func_0
{
void operator() (NodeID v) 
  {
    SP[v] = (2147483647) ;
  };
};
struct updateEdge
{
bool operator() (NodeID src, NodeID dst, int weight) 
  {
    bool output2 ;
    bool SP_trackving_var_1 = (bool) 0;
    if ( ( SP[dst]) > ( (SP[src] + weight)) ) { 
      SP[dst]= (SP[src] + weight); 
      SP_trackving_var_1 = true ; 
    } 
    output2 = SP_trackving_var_1;
    return output2;
  };
};
struct reset
{
void operator() (NodeID v) 
  {
    SP[v] = (2147483647) ;
  };
};

#if defined(USE_BASELINE) || defined(TEST_CORRECTNESS)
void SSSP_baseline(WGraph &g, int start_vertex) {
  // ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
  //   reset()(vertexsetapply_iter);
  // });;
  // int n = builtin_getVertices(edges) ;
  // VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , (0) );
  // builtin_addVertex(frontier, start_vertex) ;
  // SP[start_vertex] = (0) ;
  // int rounds = (0) ;
  // while ( (builtin_getVertexSetSize(frontier) ) != ((0) ))
  // {
  //   VertexSubset<int> *  output = edgeset_apply_pull_parallel_weighted_deduplicatied_from_vertexset_with_frontier3(edges, frontier, updateEdge()); 
  //   deleteObject(frontier) ;
  //   frontier = output;
  //   rounds = (rounds + (1) );
  //   if ((rounds) == (n))
  //     { 
  //     std::cout << "negative cycle"<< std::endl;
  //     break;
  //     } 
  // }
  // deleteObject(frontier) ;

  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    reset()(i);
  }
  int n = builtin_getVertices(edges) ;
  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , (0) );
  builtin_addVertex(frontier, start_vertex) ;
  SP[start_vertex] = (0) ;
  int rounds = (0) ;
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
          if( updateEdge() ( g.get_in_neighbors_()[i] , d, 1 ) ) { 
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
    rounds = (rounds + (1) );
    if ((rounds) == (n))
      { 
      std::cout << "negative cycle"<< std::endl;
      break;
      } 
  }
  deleteObject(frontier) ;
}

#if defined(TEST_CORRECTNESS)
void test_correctness(int start_vertex) {
  // save previous output
  int  * __restrict SP_output = new int [ builtin_getVertices(edges) ];
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    SP_output[i] = SP[i];
  }
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    reset()(vertexsetapply_iter);
  });;

  // run serial and generate output
  SSSP_baseline(edges, start_vertex);

  // compare previous output with serial-version output
  uint64_t num_diffs = 0;
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    if (SP_output[i] != SP[i]) {
      std::cout << " SP_output[i]=" << SP_output[i] << " SP[i]=" << SP[i] << " at i=" << i << std::endl;
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
void SSSP_openmp(WGraph &g, int start_vertex) {
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
  int n = builtin_getVertices(edges) ;
  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , (0) );
  builtin_addVertex(frontier, start_vertex) ;
  SP[start_vertex] = (0) ;
  int rounds = (0) ;
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
      for (uint64_t i = g.get_in_neighbors_begin_index_(d); i < g.get_in_neighbors_end_index_(d); i++) {
        if (frontier->bool_map_[g.get_in_neighbors_()[i]] ) { 
          if( updateEdge() ( g.get_in_neighbors_()[i] , d, 1 ) ) { 
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
    rounds = (rounds + (1) );
    if ((rounds) == (n))
      { 
      std::cout << "negative cycle"<< std::endl;
      break;
      } 
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

void HEARTBEAT_loop3(uint64_t startIter, uint64_t maxIter, WGraph &g, VertexSubset<int> *frontier, bool *next, uint64_t d) {
  for (uint64_t i = startIter; i < maxIter; i++) {
    if (frontier->bool_map_[g.get_in_neighbors_()[i]] ) { 
      if( updateEdge() ( g.get_in_neighbors_()[i] , d, 1 ) ) { 
        next[d] = 1; 
      }
    }
  } //end of loop on in neighbors
}

void HEARTBEAT_loop2(uint64_t maxIter, WGraph &g, VertexSubset<int> *frontier, bool *next) {
  for (uint64_t d = 0; d < maxIter; d++) {
    HEARTBEAT_loop3(g.get_in_neighbors_begin_index_(d), g.get_in_neighbors_end_index_(d), g, frontier, next, d);
  } //end of outer for loop
}

void SSSP_hbc(WGraph &g, int start_vertex) {
  HEARTBEAT_loop0(builtin_getVertices(edges));
  int n = builtin_getVertices(edges) ;
  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , (0) );
  builtin_addVertex(frontier, start_vertex) ;
  SP[start_vertex] = (0) ;
  int rounds = (0) ;
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
    rounds = (rounds + (1) );
    if ((rounds) == (n))
      { 
      std::cout << "negative cycle"<< std::endl;
      break;
      } 
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
  edges = builtin_loadWeightedEdgesFromFile ( graph_file_name ) ;
  int start_vertex = atoi( "100" ) ;
  SP = new int [ builtin_getVertices(edges) ];
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    SP_generated_vector_op_apply_func_0()(vertexsetapply_iter);
  });;
  for ( int trail = (0) ; trail < (runs) ; trail++ )
  {
    // startTimer() ;
// ####################
#if defined(USE_BASELINE) || defined(USE_FORKJOIN)
  taskparts::benchmark_nativeforkjoin([&] (auto sched) {
#elif defined(USE_OPENMP)
  utility::run([&] {
#elif defined(USE_HB_COMPILER) || defined(USE_HB_MANUAL)
  run_bench([&] {
#endif

#if defined(USE_BASELINE)
    SSSP_baseline(edges, start_vertex);
#elif defined(USE_OPENMP)
    SSSP_openmp(edges, start_vertex);
#elif defined(USE_HB_COMPILER)
    SSSP_hbc(edges, start_vertex);
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
    // float elapsed_time = stopTimer() ;
    // std::cout << "elapsed time: "<< std::endl;
    // std::cout << elapsed_time<< std::endl;
    // std::cout << "rounds"<< std::endl;
    // std::cout << rounds<< std::endl;
  }
};
#ifdef GEN_PYBIND_WRAPPERS
PYBIND11_MODULE(, m) {
}
#endif

