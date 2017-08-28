#![recursion_limit="128"]

#[macro_use]
extern crate stdweb;

#[macro_use]
extern crate serde_derive;
extern crate serde_json;

extern crate petgraph;
extern crate louvain;
extern crate chrono;

use std::collections::{HashMap};
use petgraph::{Graph};
use petgraph::graph::{NodeIndex};

use louvain::io::{Node, NodeID, Edge, Community};
use louvain::{Modularity};

fn compute_louvain(nodes_json: String, edges_json: String, resolution: f64) -> String {
    println!("Computing Louvain");
    let mut start_time = chrono::Utc::now();
    let nodes: Vec<Node> = serde_json::from_str(&nodes_json).unwrap();
    println!("Read nodes {:?}", chrono::Utc::now().signed_duration_since(start_time));
    start_time = chrono::Utc::now();
    let edges: Vec<Edge> = serde_json::from_str(&edges_json).unwrap();
    println!("Read links {:?}", chrono::Utc::now().signed_duration_since(start_time));

    println!("Nodes: {}", nodes.len());
    println!("Edges: {}", edges.len());

    start_time = chrono::Utc::now();
    let mut inv_map: HashMap<NodeID, NodeIndex> = HashMap::new();
    let mut graph = Graph::new_undirected();
    for ref node in &nodes {
        let i = graph.add_node(node.id.clone());
        inv_map.insert(node.id.clone(), i);
    }
    for edge in edges {
        graph.add_edge(
            *inv_map.get(&edge.source).unwrap(),
            *inv_map.get(&edge.target).unwrap(),
            edge.weight);
    }
    println!("Construct the graph {:?}", chrono::Utc::now().signed_duration_since(start_time));

    start_time = chrono::Utc::now();
    let mut modularity = Modularity::new(resolution);
    let results = modularity.execute(& graph);
    println!("Louvain Clusters computed in {:?}", chrono::Utc::now().signed_duration_since(start_time));

    start_time = chrono::Utc::now(); 
    let num_of_communities = modularity.communityByNode.iter().max().unwrap_or(&0) + 2; // Parent Community included
    println!("Number of Clusters: {:?} -  with resolution {:?}", num_of_communities-1, resolution);
    println!("Final Modularity: {:?}", results);
    let mut communities: Vec<Community> = vec![Default::default(); num_of_communities];
    communities[0] = Community{
        id:0,
        parent:-1,
        nodes: graph.node_indices().map(|i| i.index()).collect(),
        children: (1.. num_of_communities as i32).collect()
    };

    for (node, &com_id) in modularity.communityByNode.iter().enumerate() {
        let com = com_id + 1;
        communities[com].id = (com) as i32;
        communities[com].parent = 0;
        communities[com].nodes.push(node);
    }

    let communities_json = serde_json::to_string(&communities).unwrap();
    println!("Written communities {:?}", chrono::Utc::now().signed_duration_since(start_time));

    communities_json
}



fn main() {
    stdweb::initialize();
    println!("Hello, world!");

    js! { @(no_return)
        window.compute_louvain = @{compute_louvain};
    };
    
    stdweb::event_loop();
}
