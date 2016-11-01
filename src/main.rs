#![feature(proc_macro)]

#[macro_use]
extern crate serde_derive;

extern crate serde;
extern crate serde_json;
extern crate petgraph;

extern crate louvain;

use std::env;
use std::path::Path;
use std::error::Error;
use std::fs::File;
use std::io::Read;
use std::collections::{HashMap};

use petgraph::{Graph};
use petgraph::graph::{NodeIndex};

use louvain::{Node, Edge, Modularity};


#[derive(Serialize, Debug, Default, Clone)]
pub struct Community {
    id: i32,
    parent: i32,
    nodes: Vec<usize>,
    children: Vec<i32>,
}


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
    // println!("Reading files...");
    // let nodes: Vec<Node> = read_json("../miserables_nodes.json");
    // let edges: Vec<Edge> = read_json("../miserables_links.json");
    let args : Vec<String> = env::args().collect();
    let nodes: Vec<Node> = read_json(&args[1]);
    let edges: Vec<Edge> = read_json(&args[2]);

    //compute_louvain( &nodes, &edges);
    // println!("Nodes: {}", nodes.len());
    // println!("Edges: {}", edges.len());

    //let graph = Graph {nodes: nodes, edges: edges};

    let mut inv_map: HashMap<String, NodeIndex> = HashMap::new();
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
    let result = modularity.execute(& graph);
    // println!("Modularity: {:?}", result);

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

    println!("{}", serde_json::to_string(&communities).unwrap());
}
