use serde;
use serde_json;
use std::path::Path;
use std::io::Read;
use std::fs::File;
use std::error::Error;

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Hash, Clone)]
pub struct NodeID(String);

#[derive(Serialize, Deserialize, Debug)]
pub struct Node {
    pub id: NodeID,
    community: Option<CommunityTag>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Edge {
    pub source: NodeID,
    pub target: NodeID,
    pub weight: f32,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct CommunityTag {
    id: u32,
    color: u32,
}

#[derive(Serialize, Debug, Default, Clone)]
pub struct Community {
    pub id: i32,
    pub parent: i32,
    pub nodes: Vec<usize>,
    pub children: Vec<i32>,
}

pub fn read_json_file<T>(file_path: &str) -> Vec<T>
    where T: serde::de::DeserializeOwned
{
    let path = Path::new(file_path);

    let mut s = String::new();
    match File::open(&path).unwrap().read_to_string(&mut s) {
        Err(why) => panic!("couldn't open {}: {}", path.display(), why.description()),
        Ok(file) => file,
    };

    let array: Vec<T> = serde_json::from_str(&s).unwrap();
    
    return array;
}


pub fn write_json_file<T>(file_path: &str, array: &Vec<T>)
    where T: serde::Serialize
{
    let path = Path::new(file_path);

    let file = match File::create(&path) {
        Err(why) => panic!("couldn't create {}: {}", path.display(), why.description()),
        Ok(file) => file,
    };

    serde_json::to_writer(file, array).unwrap();
}
