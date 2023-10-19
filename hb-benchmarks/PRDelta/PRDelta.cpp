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
typedef struct struct_delta_out_degree { 
  double delta;
  int out_degree;
} struct_delta_out_degree;
double  * __restrict cur_rank;
double  * __restrict ngh_sum;
struct_delta_out_degree  * __restrict array_of_struct_delta_out_degree;
int  * __restrict generated_tmp_vector_3;
double damp; 
double beta_score; 
double epsilon2; 
double epsilon; 
template <typename APPLY_FUNC > void edgeset_apply_pull_parallel_from_vertexset5(Graph & g , VertexSubset<NodeID>* from_vertexset, APPLY_FUNC apply_func) 
{ 
    int64_t numVertices = g.num_nodes(), numEdges = g.num_edges();
  from_vertexset->toDense();
  ligra::parallel_for_lambda((NodeID)0, (NodeID)g.num_nodes(), [&] (NodeID d) {
    for(NodeID s : g.in_neigh(d)){
      if (from_vertexset->bool_map_[s] ) { 
        apply_func ( s , d );
      }
    } //end of loop on in neighbors
  }); //end of outer for loop
} //end of edgeset apply function 
struct generated_vector_op_apply_func_4
{
void operator() (NodeID v) 
  {
    array_of_struct_delta_out_degree[v].out_degree  = generated_tmp_vector_3[v];
  };
};
struct delta_generated_vector_op_apply_func_2
{
void operator() (NodeID v) 
  {
    array_of_struct_delta_out_degree[v].delta  = (((float) 1)  / builtin_getVertices(edges) );
  };
};
struct ngh_sum_generated_vector_op_apply_func_1
{
void operator() (NodeID v) 
  {
    ngh_sum[v] = ((float) 0) ;
  };
};
struct cur_rank_generated_vector_op_apply_func_0
{
void operator() (NodeID v) 
  {
    cur_rank[v] = (0) ;
  };
};
struct updateEdge
{
void operator() (NodeID src, NodeID dst) 
  {
    ngh_sum[dst] += (array_of_struct_delta_out_degree[src].delta  / array_of_struct_delta_out_degree[src].out_degree );
  };
};
struct updateVertexFirstRound
{
bool operator() (NodeID v) 
  {
    bool output ;
    array_of_struct_delta_out_degree[v].delta  = ((damp * ngh_sum[v]) + beta_score);
    cur_rank[v] += array_of_struct_delta_out_degree[v].delta ;
    array_of_struct_delta_out_degree[v].delta  = (array_of_struct_delta_out_degree[v].delta  - (((float) 1)  / builtin_getVertices(edges) ));
    output = (fabs(array_of_struct_delta_out_degree[v].delta ) ) > ((epsilon2 * cur_rank[v]));
    ngh_sum[v] = (0) ;
    return output;
  };
};
struct updateVertex
{
bool operator() (NodeID v) 
  {
    bool output ;
    array_of_struct_delta_out_degree[v].delta  = (ngh_sum[v] * damp);
    cur_rank[v] += array_of_struct_delta_out_degree[v].delta ;
    ngh_sum[v] = (0) ;
    output = (fabs(array_of_struct_delta_out_degree[v].delta ) ) > ((epsilon2 * cur_rank[v]));
    return output;
  };
};
struct printRank
{
void operator() (NodeID v) 
  {
    std::cout << cur_rank[v]<< std::endl;
  };
};
struct reset
{
void operator() (NodeID v) 
  {
    cur_rank[v] = (0) ;
    ngh_sum[v] = ((float) 0) ;
    array_of_struct_delta_out_degree[v].delta  = (((float) 1)  / builtin_getVertices(edges) );
  };
};

#if defined(USE_BASELINE) || defined(TEST_CORRECTNESS)
void PRDelta_baseline(Graph &g, int n) {
  // VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , n);
  // ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
  //   reset()(vertexsetapply_iter);
  // });;
  // for ( int i = (1) ; i < (11) ; i++ )
  // {
  //   edgeset_apply_pull_parallel_from_vertexset5(edges, frontier, updateEdge()); 
  //   VertexSubset<int> *  output ;
  //   if ((i) == ((1) ))
  //     { 
  //     output = builtin_const_vertexset_filter <updateVertexFirstRound>(updateVertexFirstRound(), builtin_getVertices(edges) );
  //     } 
  //   else
  //     { 
  //     output = builtin_const_vertexset_filter <updateVertex>(updateVertex(), builtin_getVertices(edges) );

  //     } 
  //   deleteObject(frontier) ;
  //   frontier = output;
  // }
  // deleteObject(frontier) ;

  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , n);
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    reset()(i);
  }

  for ( int i = (1) ; i < (11) ; i++ )
  {
    frontier->toDense();
    for (uint64_t d = 0; d < g.num_nodes(); d++) {
      for (uint64_t i = g.get_in_neighbors_begin_index_(d); i < g.get_in_neighbors_end_index_(d); i++) {
        if (frontier->bool_map_[g.get_in_neighbors_()[i]] ) { 
          updateEdge()( g.get_in_neighbors_()[i], d );
        }
      }
    };

    VertexSubset<int> *  output ;
    if ((i) == ((1) ))
      { 
      output = new VertexSubset<NodeID>( builtin_getVertices(edges), 0);
      bool * next0 = newA(bool, builtin_getVertices(edges));
      for (uint64_t v = 0; v < builtin_getVertices(edges); v++) {
        next0[v] = 0;
        if (updateVertexFirstRound()(v))
            next0[v] = 1;
      }
      output->num_vertices_ = sequence::sum(next0, builtin_getVertices(edges));
      output->bool_map_ = next0;
      output->is_dense = true;
      } 
    else
      { 
      output = new VertexSubset<NodeID>( builtin_getVertices(edges), 0);
      bool * next0 = newA(bool, builtin_getVertices(edges));
      for(uint64_t v = 0; v < builtin_getVertices(edges); v++) {
        next0[v] = 0;
        if (updateVertex()(v))
            next0[v] = 1;
      }
      output->num_vertices_ = sequence::sum(next0, builtin_getVertices(edges));
      output->bool_map_ = next0;
      output->is_dense = true;
      } 
    deleteObject(frontier) ;
    frontier = output;
  }
  deleteObject(frontier) ;
}

#if defined(TEST_CORRECTNESS)
void test_correctness() {
  // save previous output
  double  * __restrict cur_rank_output = new double [ builtin_getVertices(edges) ];
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    cur_rank_output[i] = cur_rank[i];
  }
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    reset()(vertexsetapply_iter);
  });;

  // run serial and generate output
  PRDelta_baseline(edges, builtin_getVertices(edges));

  // compare previous ourput with serial-version output
  uint64_t num_diffs = 0;
  double epsilon = 0.01;
  for (uint64_t i = 0; i < builtin_getVertices(edges); i++) {
    auto diff = std::abs(cur_rank_output[i] - cur_rank[i]);
    if (diff > epsilon) {
      std::cout << "diff=" << diff << " cur_rank_output[i]=" << cur_rank_output[i] << " cur_rank[i]=" << cur_rank[i] << " at i=" << i << std::endl;
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
void PRDelta_openmp(Graph &g, int n) {
  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , n);
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

  for ( int i = (1) ; i < (11) ; i++ )
  {
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
          updateEdge()( g.get_in_neighbors_()[i], d );
        }
      }
    };

    VertexSubset<int> *  output ;
    if ((i) == ((1) ))
      { 
      output = new VertexSubset<NodeID>( builtin_getVertices(edges), 0);
      bool * next0 = newA(bool, builtin_getVertices(edges));
      #if defined(OMP_SCHEDULE_STATIC)
        #pragma omp parallel for schedule(static)
      #elif defined(OMP_SCHEDULE_DYNAMIC)
        #pragma omp parallel for schedule(dynamic)
      #elif defined(OMP_SCHEDULE_GUIDED)
        #pragma omp parallel for schedule(guided)
      #endif
      for (uint64_t v = 0; v < builtin_getVertices(edges); v++) {
        next0[v] = 0;
        if (updateVertexFirstRound()(v))
            next0[v] = 1;
      }
      output->num_vertices_ = sequence::sum(next0, builtin_getVertices(edges));
      output->bool_map_ = next0;
      output->is_dense = true;
      } 
    else
      { 
      output = new VertexSubset<NodeID>( builtin_getVertices(edges), 0);
      bool * next0 = newA(bool, builtin_getVertices(edges));
      #if defined(OMP_SCHEDULE_STATIC)
        #pragma omp parallel for schedule(static)
      #elif defined(OMP_SCHEDULE_DYNAMIC)
        #pragma omp parallel for schedule(dynamic)
      #elif defined(OMP_SCHEDULE_GUIDED)
        #pragma omp parallel for schedule(guided)
      #endif
      for(uint64_t v = 0; v < builtin_getVertices(edges); v++) {
        next0[v] = 0;
        if (updateVertex()(v))
            next0[v] = 1;
      }
      output->num_vertices_ = sequence::sum(next0, builtin_getVertices(edges));
      output->bool_map_ = next0;
      output->is_dense = true;
      } 
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

void HEARTBEAT_loop2(
  uint64_t startIter,
  uint64_t maxIter,
  Graph &g,
  VertexSubset<int> *frontier,
  uint64_t d
) {
  for (uint64_t i = startIter; i < maxIter; i++) {
    if (frontier->bool_map_[g.get_in_neighbors_()[i]] ) { 
      updateEdge()( g.get_in_neighbors_()[i], d );
    }
  }
}

void HEARTBEAT_loop1(uint64_t maxIter, Graph &g, VertexSubset<int> *frontier) {
  for (uint64_t d = 0; d < maxIter; d++) {
    HEARTBEAT_loop2(g.get_in_neighbors_begin_index_(d), g.get_in_neighbors_end_index_(d), g, frontier, d);
  };
}

void HEARTBEAT_loop4(uint64_t maxIter, bool *next0) {
  for(uint64_t v = 0; v < maxIter; v++) {
    next0[v] = 0;
    if (updateVertex()(v))
        next0[v] = 1;
  }
}

void HEARTBEAT_loop3(uint64_t maxIter, bool *next0) {
  for (uint64_t v = 0; v < maxIter; v++) {
    next0[v] = 0;
    if (updateVertexFirstRound()(v))
        next0[v] = 1;
  }
}

void PRDelta_hbc(Graph &g, int n) {
  VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , n);
  HEARTBEAT_loop0(builtin_getVertices(edges));

  for ( int i = (1) ; i < (11) ; i++ )
  {
    frontier->toDense();
    HEARTBEAT_loop1(g.num_nodes(), edges, frontier);

    VertexSubset<int> *  output ;
    if ((i) == ((1) ))
      { 
      output = new VertexSubset<NodeID>( builtin_getVertices(edges), 0);
      bool * next0 = newA(bool, builtin_getVertices(edges));
      HEARTBEAT_loop3(builtin_getVertices(edges), next0);
      output->num_vertices_ = sequence::sum(next0, builtin_getVertices(edges));
      output->bool_map_ = next0;
      output->is_dense = true;
      } 
    else
      { 
      { 
      output = new VertexSubset<NodeID>( builtin_getVertices(edges), 0);
      bool * next0 = newA(bool, builtin_getVertices(edges));
      HEARTBEAT_loop4(builtin_getVertices(edges), next0);
      output->num_vertices_ = sequence::sum(next0, builtin_getVertices(edges));
      output->bool_map_ = next0;
      output->is_dense = true;
      } 
      } 
    deleteObject(frontier) ;
    frontier = output;
  }
  deleteObject(frontier) ;
}
#endif

int main(int argc, char * argv[])
{
  edges = builtin_loadEdgesFromFile ( "Twitter.el" ) ;
  cur_rank = new double [ builtin_getVertices(edges) ];
  ngh_sum = new double [ builtin_getVertices(edges) ];
  array_of_struct_delta_out_degree = new struct_delta_out_degree [ builtin_getVertices(edges) ];
  damp = ((float) 0.85) ;
  beta_score = ((((float) 1)  - damp) / builtin_getVertices(edges) );
  epsilon2 = ((float) 0.1) ;
  epsilon = ((float) 1e-07) ;
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    cur_rank_generated_vector_op_apply_func_0()(vertexsetapply_iter);
  });;
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    ngh_sum_generated_vector_op_apply_func_1()(vertexsetapply_iter);
  });;
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    delta_generated_vector_op_apply_func_2()(vertexsetapply_iter);
  });;
  generated_tmp_vector_3 = builtin_getOutDegrees(edges) ;
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    generated_vector_op_apply_func_4()(vertexsetapply_iter);
  });;
  int n = builtin_getVertices(edges) ;
  // for ( int trail = (0) ; trail < (10) ; trail++ )
  // {
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
    PRDelta_baseline(edges, n);
#elif defined(USE_OPENMP)
    PRDelta_openmp(edges, n);
#elif defined(USE_HB_COMPILER)
    PRDelta_hbc(edges, n);
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
  // }
};
#ifdef GEN_PYBIND_WRAPPERS
PYBIND11_MODULE(, m) {
}
#endif
