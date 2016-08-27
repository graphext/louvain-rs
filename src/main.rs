extern crate serde;
extern crate serde_json;

extern crate louvain;

use std::path::Path;
use std::error::Error;
use std::fs::File;
use std::io::Read;

use louvain::{Node, Edge, compute_louvain};

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
    let nodes: Vec<Node> = read_json("../nodes.json");
    let edges: Vec<Edge> = read_json("../links.json");

    compute_louvain( &nodes, &edges);
    println!("{}", nodes.len());
    println!("{}", edges.len());
}
