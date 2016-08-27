#![feature(custom_derive, plugin)]
#![plugin(serde_macros)]

#![allow(non_snake_case)]  // FIXME: Good while porting

use std::collections::HashMap;

struct Graph {
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Node {
    id: String,
    community: Option<CommunityTag>,
    // name: String,
    // _index: String,
    // _type: String,
    // date: String,
    // url: String,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct Edge {
    source: String,
    target: String,
    weight: f32,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct CommunityTag {
    id: u32,
    color: u32,
}

#[derive(Default, Clone)]
struct Community {
    id: i32,
    weightSum: f64,
//    structure: CommunityStructure,  Not here because of borrowing issues
    nodes: Vec<usize>,
    connectionsWeight: HashMap<i32, f32>,
    connectionsCount: HashMap<i32, i32>,
}

impl Community {

    fn new() -> Community {
        Default::default()
    }

    fn seed(&mut self, node: usize, weight: f64) {
        self.nodes.push(node);
        self.weightSum += weight;
    }

}

struct CommunityStructure {
    nodeConnectionsWeight: Vec<HashMap<i32, f32>>,
    nodeConnectionsCount: Vec<HashMap<i32, i32>>,
    nodeCommunities: Vec<Community>,
    map: HashMap<String, usize>,
    weights: Vec<f64>,
    graphWeightSum: f64,
    topology: Vec<Vec<Edge>>,
    communities: Vec<Community>,
    N: usize,
    invMap: HashMap<usize, Community>,
}

impl CommunityStructure {

    fn new(graph: &Graph) -> CommunityStructure {
        let N = graph.nodes.len();
        let mut cs = CommunityStructure {
            nodeConnectionsWeight: Vec::with_capacity(N),
            nodeConnectionsCount: Vec::with_capacity(N),
            nodeCommunities: Vec::with_capacity(N),
            map: HashMap::new(),
            weights: Vec::with_capacity(N),
            graphWeightSum: 0.0,
            topology: Vec::with_capacity(N),
            communities: Vec::new(),
            N: N,
            invMap: HashMap::new(),
        };

        let mut index: usize = 0;
        for node in & graph.nodes {
            cs.map.insert(node.id.clone(), index);
            cs.nodeCommunities[index] = Community::new();

            cs.nodeConnectionsWeight[index] = HashMap::new();
            cs.nodeConnectionsCount[index] = HashMap::new();
            cs.weights[index] = 0.0;
            cs.nodeCommunities[index].seed(index, cs.weights[index]);

            let mut hidden = Community::new();
            hidden.nodes.push(index);
            cs.invMap.insert(index, hidden);

            cs.communities.push( cs.nodeCommunities[index].clone() );
            index += 1;
        }

        cs
    }
}


struct Modularity {
    pub MODULARITY_CLASS: &'static str,
//    progress: ProgressTicket,
    isCanceled: bool,
//    structure: CommunityStructure,
    modularity: f64,
    modularityResolution: f64,
    isRandomized: bool,
    useWeight: bool,
    resolution: f64,
}

impl Modularity {

    fn new() -> Modularity {
        Modularity {
            MODULARITY_CLASS: "modularity_class",
            isCanceled: false,
            modularity: 0.,
            modularityResolution: 0.,
            isRandomized: false,
            useWeight: true,
            resolution: 1.0,
        }
    }
}




pub fn compute_louvain(nodes: & Vec<Node>, edges: & Vec<Edge>) {

}




#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
    }
}
