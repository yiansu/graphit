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
double  * __restrict old_rank;
double  * __restrict new_rank;
int  * __restrict out_degree;
double  * __restrict contrib;
double  * __restrict error;
int  * __restrict generated_tmp_vector_2;
double damp; 
double beta_score; 
template <typename APPLY_FUNC > void edgeset_apply_pull_parallel6(Graph & g , APPLY_FUNC apply_func) 
{ 
    int64_t numVertices = g.num_nodes(), numEdges = g.num_edges();
  ligra::parallel_for_lambda((NodeID)0, (NodeID)g.num_nodes(), [&] (NodeID d) {
    for(NodeID s : g.in_neigh(d)){
      apply_func ( s , d );
    } //end of loop on in neighbors
  }); //end of outer for loop
} //end of edgeset apply function 
struct error_generated_vector_op_apply_func_5
{
void operator() (NodeID v) 
  {
    error[v] = ((float) 0) ;
  };
};
struct contrib_generated_vector_op_apply_func_4
{
void operator() (NodeID v) 
  {
    contrib[v] = ((float) 0) ;
  };
};
struct generated_vector_op_apply_func_3
{
void operator() (NodeID v) 
  {
    out_degree[v] = generated_tmp_vector_2[v];
  };
};
struct new_rank_generated_vector_op_apply_func_1
{
void operator() (NodeID v) 
  {
    new_rank[v] = ((float) 0) ;
  };
};
struct old_rank_generated_vector_op_apply_func_0
{
void operator() (NodeID v) 
  {
    old_rank[v] = (((float) 1)  / builtin_getVertices(edges) );
  };
};
struct computeContrib
{
void operator() (NodeID v) 
  {
    contrib[v] = (old_rank[v] / out_degree[v]);
  };
};
struct updateEdge
{
void operator() (NodeID src, NodeID dst) 
  {
    new_rank[dst] += contrib[src];
  };
};
struct updateVertex
{
void operator() (NodeID v) 
  {
    double old_score = old_rank[v];
    new_rank[v] = (beta_score + (damp * new_rank[v]));
    error[v] = fabs((new_rank[v] - old_rank[v])) ;
    old_rank[v] = new_rank[v];
    new_rank[v] = ((float) 0) ;
  };
};
struct printRank
{
void operator() (NodeID v) 
  {
    std::cout << old_rank[v]<< std::endl;
  };
};
struct reset
{
void operator() (NodeID v) 
  {
    old_rank[v] = (((float) 1)  / builtin_getVertices(edges) );
    new_rank[v] = ((float) 0) ;
  };
};

#if defined(USE_BASELINE) || defined(TEST_CORRECTNESS)
void PR_baseline(Graph &g) {
  // for ( int i = (0) ; i < (20) ; i++ )
  // {
  //   ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
  //     computeContrib()(vertexsetapply_iter);
  //   });;
  //   edgeset_apply_pull_parallel6(edges, updateEdge()); 
  //   ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
  //     updateVertex()(vertexsetapply_iter);
  //   });;
  // }
  for ( int i = (0) ; i < (20) ; i++ )
  {
    for (uint64_t nodeID = 0; nodeID < builtin_getVertices(edges); nodeID++) {
      computeContrib()(nodeID);
    }

    for (uint64_t d = 0; d < g.num_nodes(); d++) {
      for (uint64_t i = g.get_in_neighbors_begin_index_(d); i < g.get_in_neighbors_end_index_(d); i++) {
        updateEdge()(g.get_in_neighbors_()[i], d);
      }
    }

    for (uint64_t nodeID = 0; nodeID < builtin_getVertices(edges); nodeID++) {
      updateVertex()(nodeID);
    }
  }
}

#if defined(TEST_CORRECTNESS)
void test_correctness() {
  // save previous output
  double  * __restrict old_rank_output = new double [ builtin_getVertices(edges) ];
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    old_rank_output[i] = old_rank[i];
  }
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    reset()(vertexsetapply_iter);
  });;

  // run serial and generate output
  PR_baseline(edges);

  // compare previous ourput with serial-version output
  uint64_t num_diffs = 0;
  double epsilon = 0.01;
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    auto diff = std::abs(old_rank_output[i] - old_rank[i]);
    if (diff > epsilon) {
      printf("diff=%f old_rank_output[i]=%f old_rank[i]=%f at i=%lu\n", diff, old_rank_output[i], old_rank[i], i);
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
void PR_openmp(Graph &g) {
  for ( int i = (0) ; i < (20) ; i++ )
  {
    #pragma omp parallel for schedule(static)
    for (uint64_t nodeID = 0; nodeID < builtin_getVertices(edges); nodeID++) {
      computeContrib()(nodeID);
    }

    #if defined(OMP_SCHEDULE_STATIC)
      #pragma omp parallel for schedule(static)
    #elif defined(OMP_SCHEDULE_DYNAMIC)
      #pragma omp parallel for schedule(dynamic)
    #elif defined(OMP_SCHEDULE_GUIDED)
      #pragma omp parallel for schedule(guided)
    #endif
    for (uint64_t d = 0; d < g.num_nodes(); d++) {
      for (uint64_t i = g.get_in_neighbors_begin_index_(d); i < g.get_in_neighbors_end_index_(d); i++) {
        updateEdge()(g.get_in_neighbors_()[i], d);
      }
    }

    #pragma omp parallel for schedule(static)
    for (uint64_t nodeID = 0; nodeID < builtin_getVertices(edges); nodeID++) {
      updateVertex()(nodeID);
    }
  }
}

#elif defined(USE_HB_COMPILER)
#if defined(RUN_HEARTBEAT)
  bool run_heartbeat = true;
#else
  bool run_heartbeat = false;
#endif

void HEARTBEAT_loop0(uint64_t maxIter) {
  for (uint64_t nodeID = 0; nodeID < maxIter; nodeID++) {
    computeContrib()(nodeID);
  }
}

void HEARTBEAT_loop2(
  uint64_t startIter,
  uint64_t maxIter,
  Graph &g,
  uint64_t d
) {
  for (uint64_t i = startIter; i < maxIter; i++) {
    updateEdge()(g.get_in_neighbors_()[i], d);
  }
}

void HEARTBEAT_loop1(Graph &g, uint64_t maxIter) {
  for (uint64_t d = 0; d < maxIter; d++) {
    HEARTBEAT_loop2(g.get_in_neighbors_begin_index_(d), g.get_in_neighbors_end_index_(d), g, d);
  }
}

void HEARTBEAT_loop3(uint64_t maxIter) {
  for (uint64_t nodeID = 0; nodeID < maxIter; nodeID++) {
    updateVertex()(nodeID);
  }
}

void PR_hbc(Graph &g) {
  for ( int i = (0) ; i < (20) ; i++ )
  {
    HEARTBEAT_loop0(builtin_getVertices(edges));
    HEARTBEAT_loop1(g, g.num_nodes());
    HEARTBEAT_loop3(builtin_getVertices(edges));
  }
}
#endif


int main(int argc, char * argv[])
{
  edges = builtin_loadEdgesFromFile ("USAroad.el") ;
  old_rank = new double [ builtin_getVertices(edges) ];
  new_rank = new double [ builtin_getVertices(edges) ];
  out_degree = new int [ builtin_getVertices(edges) ];
  contrib = new double [ builtin_getVertices(edges) ];
  error = new double [ builtin_getVertices(edges) ];
  damp = ((float) 0.85) ;
  beta_score = ((((float) 1)  - damp) / builtin_getVertices(edges) );
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    old_rank_generated_vector_op_apply_func_0()(vertexsetapply_iter);
  });;
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    new_rank_generated_vector_op_apply_func_1()(vertexsetapply_iter);
  });;
  generated_tmp_vector_2 = builtin_getOutDegrees(edges) ;
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    generated_vector_op_apply_func_3()(vertexsetapply_iter);
  });;
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    contrib_generated_vector_op_apply_func_4()(vertexsetapply_iter);
  });;
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    error_generated_vector_op_apply_func_5()(vertexsetapply_iter);
  });;
  // for ( int trail = (0) ; trail < (10) ; trail++ )
  // {
    // startTimer() ;
    ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
      reset()(vertexsetapply_iter);
    });;
// ####################
#if defined(USE_BASELINE) || defined(USE_FORKJOIN)
  taskparts::benchmark_nativeforkjoin([&] (auto sched) {
#elif defined(USE_OPENMP)
  utility::run([&] {
#elif defined(USE_HB_COMPILER) || defined(USE_HB_MANUAL)
  run_bench([&] {
#endif

#if defined(USE_BASELINE)
    PR_baseline(edges);
#elif defined(USE_OPENMP)
    PR_openmp(edges);
#elif defined(USE_HB_COMPILER)
    PR_hbc(edges);
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
    // double elapsed_time = stopTimer() ;
    // std::cout << "elapsed time: "<< std::endl;
    // std::cout << elapsed_time<< std::endl;
  }
// };
#ifdef GEN_PYBIND_WRAPPERS
PYBIND11_MODULE(, m) {
}
#endif

