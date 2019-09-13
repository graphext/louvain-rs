use std::collections::{HashMap};

use petgraph::{Graph};
use petgraph::graph::{NodeIndex};
use clap::{Arg, App};

use serde_json::{Value};

use louvain::{Modularity};
use louvain::io::{Node, NodeID, Edge, read_json_file, write_json_file};

fn main() {
    let mut start_time = chrono::Utc::now();

    let matches = App::new("Louvain")
        .version("1.0")
        .author("Juan Morales <crispamares@gmail.com>")
        .about("Compute louvain clustering")
        .arg(Arg::with_name("nodes")
            .short("n")
            .long("nodes")
            .value_name("FILE")
            .help("Sets the nodes input file")
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("links")
            .short("l")
            .long("links")
            .value_name("FILE")
            .help("Sets the links input file")
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("output")
            .short("o")
            .long("output")
            .value_name("FILE")
            .help("Sets the output nodes file")
            .required(true)
            .takes_value(true))
        .arg(Arg::with_name("resolution")
            .short("r")
            .long("resolution")
            .value_name("FLOAT")
            .help("Sets the resolution. Higher this value the bigger the clusters")
            .default_value("0.1")
            .takes_value(true))
        .arg(Arg::with_name("noise")
            .long("noise")
            .value_name("INT")
            .help("Sets the minimum of nodes in a cluster to not be considered as noise")
            .default_value("1")
            .takes_value(true))
        .get_matches();


    let nodes: Vec<Node> = read_json_file(matches.value_of("nodes").unwrap());
    println!("Read nodes {}", chrono::Utc::now().signed_duration_since(start_time));
    start_time = chrono::Utc::now();
    let edges: Vec<Edge> = read_json_file(matches.value_of("links").unwrap());
    println!("Read links {}", chrono::Utc::now().signed_duration_since(start_time));
    
    let resolution : f64 = match matches.value_of("resolution") {
        Some(res) => res.parse().expect("Resolution should be a float. Default = 1.0"),
        None => 1.0
    };
    let noise : u32 = match matches.value_of("noise") {
        Some(res) => res.parse().expect("Noise should be a u32. Default = 1"),
        None => 1
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
    println!("Construct the graph computed in {}", chrono::Utc::now().signed_duration_since(start_time));


    start_time = chrono::Utc::now();
    let mut modularity = Modularity::new(resolution, noise);
    let results = modularity.execute(& graph);
    println!("Louvain Clusters computed in {}", chrono::Utc::now().signed_duration_since(start_time));

    start_time = chrono::Utc::now();
    let num_of_communities = *modularity.communityByNode.iter().max().unwrap_or(&0) as usize + 2; // Parent Community included
    println!("Number of Clusters: {} -  with resolution {}", num_of_communities-1, resolution);
    println!("Final Modularity: {:?}", results);

    let mut all_nodes: Vec<Value> = read_json_file(matches.value_of("nodes").unwrap());
    for (i, node) in all_nodes.iter_mut().enumerate() {
        let map = node.as_object_mut().unwrap();
        map.insert("_cluster".into(), modularity.communityByNode[i].into());
    }
    println!("Read and modify nodes {}", chrono::Utc::now().signed_duration_since(start_time));

    start_time = chrono::Utc::now();
    write_json_file(matches.value_of("output").unwrap(), &all_nodes);
    println!("Written output {}", chrono::Utc::now().signed_duration_since(start_time));

}
