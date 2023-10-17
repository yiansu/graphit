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
int main(int argc, char * argv[])
{
  edges = builtin_loadEdgesFromFile ( argv_safe((1) , argv, argc)) ;
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
  for ( int trail = (0) ; trail < (10) ; trail++ )
  {
    startTimer() ;
    ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
      reset()(vertexsetapply_iter);
    });;
    for ( int i = (0) ; i < (20) ; i++ )
    {
      ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
        computeContrib()(vertexsetapply_iter);
      });;
      edgeset_apply_pull_parallel6(edges, updateEdge()); 
      ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
        updateVertex()(vertexsetapply_iter);
      });;
    }
    double elapsed_time = stopTimer() ;
    std::cout << "elapsed time: "<< std::endl;
    std::cout << elapsed_time<< std::endl;
  }
};
#ifdef GEN_PYBIND_WRAPPERS
PYBIND11_MODULE(, m) {
}
#endif

