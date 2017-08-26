extern crate petgraph;
extern crate chrono;
extern crate louvain;

use std::env;
use std::collections::{HashMap};

use petgraph::{Graph};
use petgraph::graph::{NodeIndex};

use louvain::{Modularity};
use louvain::io::{Node, NodeID, Edge, Community, read_json_file, write_json_file};

fn main() { 
    let mut start_time = chrono::Utc::now();
    // println!("Reading files...");
    // let nodes: Vec<Node> = read_json_file("../miserables_nodes.json");
    // let edges: Vec<Edge> = read_json_file("../miserables_links.json");
    let args : Vec<String> = env::args().collect();
    let nodes: Vec<Node> = read_json_file(&args[1]);
    println!("Read nodes {:?}", chrono::Utc::now().signed_duration_since(start_time));
    start_time = chrono::Utc::now();
    let edges: Vec<Edge> = read_json_file(&args[2]);
    println!("Read links {:?}", chrono::Utc::now().signed_duration_since(start_time));

    let resolution : f64 = match args.get(3) {
        Some(res) => res.parse().expect("Resolution should be a float. Default = 1.0"),
        None => 1.0
    };
    //compute_louvain( &nodes, &edges);
    println!("Nodes: {}", nodes.len());
    println!("Edges: {}", edges.len());

    //let graph = Graph {nodes: nodes, edges: edges};

    start_time = chrono::Utc::now();
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
    println!("Construct the graph computed in {:?}", chrono::Utc::now().signed_duration_since(start_time));


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

    // println!("{:?}", communities);

    write_json_file("communities.json", &communities);
    println!("Written communities {:?}", chrono::Utc::now().signed_duration_since(start_time));
}
