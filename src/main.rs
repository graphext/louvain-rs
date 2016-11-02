#![feature(proc_macro)]

#[macro_use]
extern crate serde_derive;
extern crate petgraph;
extern crate serde_json;
extern crate louvain;
extern crate chrono;

use std::env;
use std::collections::{HashMap};

use petgraph::{Graph};
use petgraph::graph::{NodeIndex};

use louvain::{Modularity};
use louvain::io::{Node, NodeID, Edge, Community, read_json_file, write_json_file};


fn main() {
    // println!("Reading files...");
    // let nodes: Vec<Node> = read_json_file("../miserables_nodes.json");
    // let edges: Vec<Edge> = read_json_file("../miserables_links.json");
    let args : Vec<String> = env::args().collect();
    let nodes: Vec<Node> = read_json_file(&args[1]);
    let edges: Vec<Edge> = read_json_file(&args[2]);

    //compute_louvain( &nodes, &edges);
    // println!("Nodes: {}", nodes.len());
    // println!("Edges: {}", edges.len());

    //let graph = Graph {nodes: nodes, edges: edges};

    let mut inv_map: HashMap<NodeID, NodeIndex> = HashMap::new();
    let mut graph = Graph::new_undirected();
    for ref node in &nodes {
        let i = graph.add_node(node.id.clone());
        inv_map.insert(node.id.clone(), i);
    }
    for edge in edges {
        // match graph.find_edge(*inv_map.get(&edge.source).unwrap(), *inv_map.get(&edge.target).unwrap()) {
        //     Some(e) => println!("----> There is a parallel edge!!! {:?} -> {:?}", *inv_map.get(&edge.source).unwrap(), *inv_map.get(&edge.target).unwrap()),
        //     None => ()
        // }
        graph.add_edge(
            *inv_map.get(&edge.source).unwrap(),
            *inv_map.get(&edge.target).unwrap(),
            edge.weight);
    }

    let mut modularity = Modularity::new();
    let start_time = chrono::UTC::now();
    modularity.execute(& graph);
    println!("Clusters compued in {:?}", chrono::UTC::now() - start_time);

    let mut communities: Vec<Community> = vec![Default::default(); modularity.communityByNode.iter().max().unwrap_or(&0)+1];
    communities.push(Community{
        id:0,
        parent:-1,
        nodes: graph.node_indices().map(|i| i.index()).collect(),
        children: Vec::new()
    });

    for (node, &com_id) in modularity.communityByNode.iter().enumerate() {
        let com = com_id + 1;
        communities[com].id = (com) as i32;
        communities[com].parent = 0;
        communities[com].nodes.push(node);
    }

    // println!("{:?}", communities);

    write_json_file("communities.json", &communities);
}
