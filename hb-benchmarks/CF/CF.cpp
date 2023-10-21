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
typedef double double20 [ 20]; 
double20 * __restrict  latent_vec;
double20 * __restrict  error_vec;
double step; 
double lambda; 
int K; 
template <typename APPLY_FUNC > void edgeset_apply_pull_parallel_weighted0(WGraph & g , APPLY_FUNC apply_func) 
{ 
    int64_t numVertices = g.num_nodes(), numEdges = g.num_edges();
  ligra::parallel_for_lambda((NodeID)0, (NodeID)g.num_nodes(), [&] (NodeID d) {
    for(WNode s : g.in_neigh(d)){
      apply_func ( s.v , d, s.w );
    } //end of loop on in neighbors
  }); //end of outer for loop
} //end of edgeset apply function 
struct updateEdge
{
void operator() (NodeID src, NodeID dst, int rating) 
  {
    double estimate = (0) ;
    for ( int i = (0) ; i < K; i++ )
    {
      estimate += (latent_vec[src][i] * latent_vec[dst][i]);
    }
    double err = (rating - estimate);
    for ( int i = (0) ; i < K; i++ )
    {
      error_vec[dst][i] += (latent_vec[src][i] * err);
    }
  };
};
struct updateVertex
{
void operator() (NodeID v) 
  {
    for ( int i = (0) ; i < K; i++ )
    {
      latent_vec[v][i] += (step * (( -lambda * latent_vec[v][i]) + error_vec[v][i]));
      error_vec[v][i] = (0) ;
    }
  };
};
struct initVertex
{
void operator() (NodeID v) 
  {
    for ( int i = (0) ; i < K; i++ )
    {
      latent_vec[v][i] = ((float) 0.5) ;
      error_vec[v][i] = (0) ;
    }
  };
};

#if defined(USE_BASELINE) || defined(TEST_CORRECTNESS)
void CF_baseline(WGraph &g) {
  // ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
  //   initVertex()(vertexsetapply_iter);
  // });;
  // for ( int i = (0) ; i < (10) ; i++ )
  // {
  //   edgeset_apply_pull_parallel_weighted0(edges, updateEdge()); 
  //   ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
  //     updateVertex()(vertexsetapply_iter);
  //   });;
  // }

  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    initVertex()(i);
  }
  for ( int i = (0) ; i < (10) ; i++ )
  {
    for (uint64_t d = 0; d < g.num_nodes(); d++) {
      for (uint64_t i = g.get_in_neighbors_begin_index_(d); i < g.get_in_neighbors_end_index_(d); i++) {
        updateEdge() ( g.get_in_neighbors_()[i] , d, 1 );
      } //end of loop on in neighbors
    } //end of outer for loop
    for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
      updateVertex()(i);
    }
  }
}

#if defined(TEST_CORRECTNESS)
void test_correctness() {
  // save previous output
  double20 * __restrict  latent_vec_output = new double20 [ builtin_getVertices(edges) ];
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    for (uint64_t k = 0; k < K; k++) {
      latent_vec_output[i][k] = latent_vec[i][k];
    }
  }
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    initVertex()(vertexsetapply_iter);
  });;

  // run serial and generate output
  CF_baseline(edges);

  // compare previous ourput with serial-version output
  uint64_t num_diffs = 0;
  double epsilon = 0.01;
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    for (uint64_t k = 0; k < K; k++) {
      auto diff = std::abs(latent_vec_output[i][k] - latent_vec[i][k]);
      if (diff > epsilon) {
        std::cout << "diff=" << diff << " latent_vec_output[i][k]=" << latent_vec_output[i][k] << " latent_vec[i][k]=" << latent_vec[i][k] << " at i=" << i << ", k=" << k << std::endl;
        num_diffs++;
      }
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
void CF_openmp(WGraph &g) {
  #if defined(OMP_SCHEDULE_STATIC)
    #pragma omp parallel for schedule(static)
  #elif defined(OMP_SCHEDULE_DYNAMIC)
    #pragma omp parallel for schedule(dynamic)
  #elif defined(OMP_SCHEDULE_GUIDED)
    #pragma omp parallel for schedule(guided)
  #endif
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    initVertex()(i);
  }
  for ( int i = (0) ; i < (10) ; i++ )
  {
    #if defined(OMP_SCHEDULE_STATIC)
      #pragma omp parallel for schedule(static)
    #elif defined(OMP_SCHEDULE_DYNAMIC)
      #pragma omp parallel for schedule(dynamic)
    #elif defined(OMP_SCHEDULE_GUIDED)
      #pragma omp parallel for schedule(guided)
    #endif
    for (uint64_t d = 0; d < g.num_nodes(); d++) {
      for (uint64_t i = g.get_in_neighbors_begin_index_(d); i < g.get_in_neighbors_end_index_(d); i++) {
        updateEdge() ( g.get_in_neighbors_()[i] , d, 1 );
      } //end of loop on in neighbors
    } //end of outer for loop
    #if defined(OMP_SCHEDULE_STATIC)
      #pragma omp parallel for schedule(static)
    #elif defined(OMP_SCHEDULE_DYNAMIC)
      #pragma omp parallel for schedule(dynamic)
    #elif defined(OMP_SCHEDULE_GUIDED)
      #pragma omp parallel for schedule(guided)
    #endif
    for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
      updateVertex()(i);
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
  for (uint64_t i = 0; i < maxIter; i++) {
    initVertex()(i);
  }
}

void HEARTBEAT_loop2(uint64_t startIter, uint64_t maxIter, WGraph &g, uint64_t d) {
  for (uint64_t i = startIter; i < maxIter; i++) {
    updateEdge() ( g.get_in_neighbors_()[i] , d, 1 );
  } //end of loop on in neighbors
}

void HEARTBEAT_loop1(uint64_t maxIter, WGraph &g) {
  for (uint64_t d = 0; d < maxIter; d++) {
    HEARTBEAT_loop2(g.get_in_neighbors_begin_index_(d), g.get_in_neighbors_end_index_(d), g, d);
  } //end of outer for loop
}

void HEARTBEAT_loop3(uint64_t maxIter) {
  for (uint64_t i = 0; i < maxIter; i++) {
    updateVertex()(i);
  }
}

void CF_hbc(WGraph &g) {
  HEARTBEAT_loop0(builtin_getVertices(edges));
  for ( int i = (0) ; i < (10) ; i++ )
  {
    HEARTBEAT_loop1(g.num_nodes(), edges);
    HEARTBEAT_loop3(builtin_getVertices(edges));
  }
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
  latent_vec = new double20 [ builtin_getVertices(edges) ];
  error_vec = new double20 [ builtin_getVertices(edges) ];
  step = ((float) 3.5e-07) ;
  lambda = ((float) 0.001) ;
  K = (20) ;
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
    CF_baseline(edges);
#elif defined(USE_OPENMP)
    CF_openmp(edges);
#elif defined(USE_HB_COMPILER)
    CF_hbc(edges);
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
};
#ifdef GEN_PYBIND_WRAPPERS
PYBIND11_MODULE(, m) {
}
#endif

