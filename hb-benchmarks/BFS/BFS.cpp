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
int main(int argc, char * argv[])
{
  edges = builtin_loadEdgesFromFile ( "USAroad.el" ) ;
  parent = new int [ builtin_getVertices(edges) ];
  ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
    parent_generated_vector_op_apply_func_0()(vertexsetapply_iter);
  });;
  for ( int trail = (0) ; trail < (10) ; trail++ )
  {
    startTimer() ;
    ligra::parallel_for_lambda((int)0, (int)builtin_getVertices(edges) , [&] (int vertexsetapply_iter) {
      reset()(vertexsetapply_iter);
    });;
    VertexSubset<int> *  frontier = new VertexSubset<int> ( builtin_getVertices(edges)  , (0) );
    int start_vertex = atoi(argv_safe((2) , argv, argc)) ;
    builtin_addVertex(frontier, start_vertex) ;
    parent[start_vertex] = start_vertex;
    while ( (builtin_getVertexSetSize(frontier) ) != ((0) ))
    {
      VertexSubset<int> *  output = edgeset_apply_pull_parallel_from_vertexset_to_filter_func_with_frontier2(edges, frontier, toFilter(), updateEdge()); 
      deleteObject(frontier) ;
      frontier = output;
    }
    deleteObject(frontier) ;
    float elapsed_time = stopTimer() ;
    std::cout << "elapsed time: "<< std::endl;
    std::cout << elapsed_time<< std::endl;
  }
};
#ifdef GEN_PYBIND_WRAPPERS
PYBIND11_MODULE(, m) {
}
#endif

