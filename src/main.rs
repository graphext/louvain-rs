extern crate serde;
extern crate serde_json;
extern crate petgraph;

extern crate louvain;

use std::path::Path;
use std::error::Error;
use std::fs::File;
use std::io::Read;
use std::collections::{HashMap};

use petgraph::{Graph};
use petgraph::graph::{NodeIndex};

use louvain::{Node, Edge, Modularity};



fn read_json<T>(file_path: &str) -> Vec<T>
    where T: serde::Deserialize
{
    let path = Path::new(file_path);

    let mut file = match File::open(&path) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.description()),
        Ok(file) => file,
    };

    let mut buffer = String::new();
    file.read_to_string(&mut buffer).unwrap();

    let array: Vec<T> = serde_json::from_str(&buffer).unwrap();
    return array;
}


fn main() {
    println!("Reading files...");
    // let nodes: Vec<Node> = read_json("../miserables_nodes.json");
    // let edges: Vec<Edge> = read_json("../miserables_links.json");
    let nodes: Vec<Node> = read_json("../nodes.json");
    let edges: Vec<Edge> = read_json("../links.json");

    //compute_louvain( &nodes, &edges);
    println!("Nodes: {}", nodes.len());
    println!("Edges: {}", edges.len());

    //let graph = Graph {nodes: nodes, edges: edges};

    let mut inv_map: HashMap<String, NodeIndex> = HashMap::new();
    let mut graph = Graph::new_undirected();
    for node in nodes {
        let i = graph.add_node(node.id.clone());
        inv_map.insert(node.id, i);
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

    // println!("Creating CommunityStructure...");
    //let cs = CommunityStructure::new( & graph, &mut modularity);

    let result = modularity.execute(& graph);
    println!("{:?}", result);
}
