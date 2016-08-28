#![feature(custom_derive, plugin)]
#![plugin(serde_macros)]

#![allow(non_snake_case)]  // FIXME: Good while porting

use std::collections::{HashMap, HashSet};
use std::cell::RefCell;
use std::rc::Rc;
use std::iter::FromIterator;
use std::ops::{AddAssign, SubAssign};

struct Graph {
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
}

impl Graph {

    fn getNeighbors(& self, node_id: & String) -> Vec<&String> {
        let mut neighbours = Vec::new();
        for edge in &self.edges {
            if edge.source == *node_id {
                neighbours.push(&edge.target);
            }
        }
        neighbours
    }

    fn getEdges(& self, source_id: &String, target_id: &String) -> Vec<&Edge> {
        let mut parallel_edges = Vec::new();
        for edge in &self.edges {
            if edge.source == *source_id && edge.target == *target_id {
                parallel_edges.push(edge);
            }
        }
        parallel_edges
    }
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

#[derive(Debug)]
pub struct ModEdge {
    source: usize,
    target: usize,
    weight: f32,
}

#[derive(Serialize, Deserialize, Debug)]
pub struct CommunityTag {
    id: u32,
    color: u32,
}

#[derive(Default, PartialEq)]
struct Community {
    id: i32,
    weightSum: f32,
//    structure: CommunityStructure,  Not here because of borrowing issues
    nodes: Vec<usize>,
    connectionsWeight: RefCell<HashMap<i32, f32>>,
    connectionsCount: RefCell<HashMap<i32, i32>>,
}

impl Community {

    fn new() -> Community {
        Default::default()
    }

    fn seed(&mut self, node: usize, cs: & CommunityStructure) {
        self.nodes.push(node);
        self.weightSum += cs.weights[node];
    }

    fn add(&mut self, node: usize, cs: & CommunityStructure) -> bool {
        self.nodes.push(node);
        self.weightSum += cs.weights[node];
        true
    }

    fn remove(&mut self, node: usize, cs: &mut CommunityStructure) -> bool {
        self.nodes.remove(node);
        self.weightSum -= cs.weights[node];
        if self.nodes.is_empty() {
            cs.communities.retain(|ref c| *c.borrow() != *self);
        }
        true
    }


}

struct CommunityStructure {
    nodeConnectionsWeight: Vec<HashMap<i32, f32>>,
    nodeConnectionsCount: Vec<HashMap<i32, i32>>,
    nodeCommunities: Vec<Rc<RefCell<Community>>>,
    map: HashMap<String, usize>,
    weights: Vec<f32>,
    graphWeightSum: f32,
    topology: Vec<Vec<ModEdge>>,
    communities: Vec<Rc<RefCell<Community>>>,
    N: usize,
    invMap: HashMap<usize, Community>,
}

impl CommunityStructure {

    fn new(graph: &Graph, modularity: &Modularity) -> CommunityStructure {
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
            cs.nodeCommunities[index] = Rc::new(RefCell::new(Community::new()));

            cs.nodeConnectionsWeight[index] = HashMap::new();
            cs.nodeConnectionsCount[index] = HashMap::new();
            cs.weights[index] = 0.0;
            cs.nodeCommunities[index].borrow_mut().seed(index, &cs);

            let mut hidden = Community::new();
            hidden.nodes.push(index);
            cs.invMap.insert(index, hidden);

            cs.communities.push( cs.nodeCommunities[index].clone() );
            index += 1;
        }

        for node in &graph.nodes {
            let node_index = cs.map.get(&node.id).unwrap();
            cs.topology[*node_index] = Vec::new();

            let uniqueNeighbors : HashSet<&String> = HashSet::from_iter(graph.getNeighbors(&node.id));
            for neighbor in uniqueNeighbors {
                if node.id == *neighbor { continue }
                let neighbor_index = cs.map.get(neighbor).unwrap();

                //Sum all parallel edges weight:
                let mut weight : f32  = 0.0;
                for edge in graph.getEdges(&node.id, neighbor) {
                    weight += if modularity.useWeight {edge.weight} else {1.0};
                }

                //Finally add a single edge with the summed weight of all parallel edges:
                cs.weights[*node_index] += weight;
                let modularity_edge = ModEdge {source: *node_index, target: *neighbor_index, weight: weight};
                cs.topology[*node_index].push(modularity_edge);
                let adjCom = &cs.nodeCommunities[*neighbor_index];

                cs.nodeConnectionsWeight[*node_index].insert(adjCom.borrow().id, weight);
                cs.nodeConnectionsCount[*node_index].insert(adjCom.borrow().id, 1);

                let nodeCom = &cs.nodeCommunities[*node_index];
                nodeCom.borrow().connectionsWeight.borrow_mut().insert(adjCom.borrow().id, weight);
                nodeCom.borrow().connectionsCount.borrow_mut().insert(adjCom.borrow().id, 1);

                cs.nodeConnectionsWeight[*neighbor_index].insert(nodeCom.borrow().id, weight);
                cs.nodeConnectionsCount[*neighbor_index].insert(nodeCom.borrow().id, 1);

                adjCom.borrow_mut().connectionsWeight.borrow_mut().insert(nodeCom.borrow().id, weight);
                adjCom.borrow_mut().connectionsCount.borrow_mut().insert(nodeCom.borrow().id, 1);

                cs.graphWeightSum += weight;
            }
        }

        cs.graphWeightSum /= 2.0;
        cs
    }

    fn addNodeTo(&mut self, node : usize, to: Rc<RefCell<Community>>) {

        fn add_node<V: AddAssign + Copy>(map : &mut HashMap<i32,V>, key: i32, weight: V) {
            let w = map.entry(key).or_insert(weight);
            *w += weight;
        }

        self.nodeCommunities[node] = to.clone();

        let mut com = to.borrow_mut();
        com.add(node, & self);

        for e in &self.topology[node] {
            let neighbor = &e.target;

            add_node(&mut self.nodeConnectionsWeight[*neighbor], com.id, e.weight);
            add_node(&mut self.nodeConnectionsCount[*neighbor], com.id, 1);

            let mut adjCom = self.nodeCommunities[*neighbor].borrow_mut();

            add_node(&mut adjCom.connectionsWeight.borrow_mut(), com.id, e.weight);
            add_node(&mut adjCom.connectionsCount.borrow_mut(), com.id, 1);

            add_node(&mut self.nodeConnectionsWeight[node], adjCom.id, e.weight);
            add_node(&mut self.nodeConnectionsCount[node], adjCom.id, 1);

            if com.id != adjCom.id {
                add_node(&mut com.connectionsWeight.borrow_mut(), adjCom.id, e.weight);
                add_node(&mut com.connectionsCount.borrow_mut(), adjCom.id, 1);
            }

        }
    }


    fn removeNodeFromItsCommunity(&mut self, node: usize) {

        fn remove_node( weights_map : &mut HashMap<i32,f32>,
            count_map : &mut HashMap<i32,i32>,
            key: i32, weight: f32
        )
        {
            let count = count_map.get(&key).unwrap().clone();

            if count -1 == 0 {
                weights_map.remove(&key);
                count_map.remove(&key);
            } else {
                let count = count_map.get_mut(&key).unwrap();
                *count -= 1;
                let w = weights_map.get_mut(&key).unwrap();
                *w -= weight;
            }
        }

        let mut community = self.nodeCommunities[node].borrow_mut();

        for e in &self.topology[node] {
            let neighbor = &e.target;

            ////////
            //Remove Node Connection to this community
            remove_node( &mut self.nodeConnectionsWeight[*neighbor],
                &mut self.nodeConnectionsCount[*neighbor],
                community.id, e.weight );

            ///////////////////
            //Remove Adjacency Community's connection to this community
            let mut adjCom = self.nodeCommunities[*neighbor].borrow_mut();
            //let connectionsWeight = &adjCom.connectionsWeight;
            //let connectionsCount = &adjCom.connectionsCount;
            remove_node( &mut adjCom.connectionsWeight.borrow_mut(),
                &mut adjCom.connectionsCount.borrow_mut(),
                community.id, e.weight );

        }
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
