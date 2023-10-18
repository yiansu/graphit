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
int main(int argc, char * argv[])
{
  edges = builtin_loadEdgesFromFile ("USAroad.el") ;
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
  for ( int trail = (0) ; trail < (10) ; trail++ )
  {
    startTimer() ;
    VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , n);
    ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
      reset()(vertexsetapply_iter);
    });;
    for ( int i = (1) ; i < (11) ; i++ )
    {
      edgeset_apply_pull_parallel_from_vertexset5(edges, frontier, updateEdge()); 
      VertexSubset<int> *  output ;
      if ((i) == ((1) ))
       { 
        output = builtin_const_vertexset_filter <updateVertexFirstRound>(updateVertexFirstRound(), builtin_getVertices(edges) );
       } 
      else
       { 
        output = builtin_const_vertexset_filter <updateVertex>(updateVertex(), builtin_getVertices(edges) );

       } 
      deleteObject(frontier) ;
      frontier = output;
    }
    deleteObject(frontier) ;
    double elapsed_time = stopTimer() ;
    std::cout << "elapsed time: "<< std::endl;
    std::cout << elapsed_time<< std::endl;
  }
};
#ifdef GEN_PYBIND_WRAPPERS
PYBIND11_MODULE(, m) {
}
#endif

