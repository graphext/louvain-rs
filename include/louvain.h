#include <cstdint>

using NodeID = uint32_t;

struct Edge {
  NodeID source;
  NodeID target;
  float weight;
};

extern "C" {

void louvain(const NodeID *nodes_ptr,
             uintptr_t nodes_len,
             const Edge *edges_ptr,
             uintptr_t edges_len,
             double resolution,
             uint32_t noise,
             int32_t *communities_ptr,
             uintptr_t communities_len);

} // extern "C"
