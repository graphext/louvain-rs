use std::{collections::HashMap, slice::from_raw_parts, slice::from_raw_parts_mut};

use petgraph::graph::NodeIndex;
use petgraph::Graph;

use crate::io::{Edge, NodeID};
use crate::Modularity;

#[no_mangle]
pub extern "C" fn louvain(
    nodes_ptr: *const NodeID,
    nodes_len: usize,
    edges_ptr: *const Edge,
    edges_len: usize,
    resolution: f64,
    noise: u32,
    communities_ptr: *mut i32,
    communities_len: usize,
) {
    let mut start_time = chrono::Utc::now();

    let nodes = unsafe { from_raw_parts(nodes_ptr, nodes_len) };
    let edges = unsafe { from_raw_parts(edges_ptr, edges_len) };
    let communities = unsafe { from_raw_parts_mut(communities_ptr, communities_len) };

    println!("Nodes: {}", nodes.len());
    println!("Edges: {}", edges.len());

    let mut inv_map: HashMap<NodeID, NodeIndex> = HashMap::new();
    let mut graph = Graph::new_undirected();
    for ref node in nodes {
        let i = graph.add_node(**node);
        inv_map.insert(**node, i);
    }
    for edge in edges {
        let a = *inv_map.get(&edge.source).unwrap();
        let b = *inv_map.get(&edge.target).unwrap();
        if let Some(edge_idx) = graph.find_edge(a, b) {
            let weight = graph[edge_idx];
            graph.update_edge(a, b, weight + edge.weight);
        } else {
            graph.add_edge(a, b, edge.weight);
        }
    }
    println!(
        "Construct the graph computed in {}",
        chrono::Utc::now().signed_duration_since(start_time)
    );

    start_time = chrono::Utc::now();
    let mut modularity = Modularity::new(resolution, noise);
    let results = modularity.execute(&graph);
    println!(
        "Louvain Clusters computed in {}",
        chrono::Utc::now().signed_duration_since(start_time)
    );

    let num_of_communities = *modularity.communityByNode.iter().max().unwrap_or(&0) as usize + 1; // Parent Community included
    println!(
        "Number of Clusters: {} -  with resolution {}",
        num_of_communities, resolution
    );
    println!("Final Modularity: {:?}", results);

    communities.copy_from_slice(&modularity.communityByNode);
}
