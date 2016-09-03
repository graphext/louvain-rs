#![feature(custom_derive, plugin)]
#![plugin(serde_macros)]

#![allow(non_snake_case)]  // FIXME: Good while porting

use std::collections::{HashMap, HashSet};
use std::cell::{RefCell, RefMut};
use std::rc::{Rc, Weak};
use std::iter::FromIterator;
use std::cmp::PartialEq;
use std::ops::{AddAssign, SubAssign};

pub struct Graph {
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

#[derive(Debug, PartialEq)]
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


type CommunityId = i32;

#[derive(Debug, Default, PartialEq)]
pub struct Community {
    id: CommunityId,
    weightSum: f32,
    nodes: Vec<usize>,
    connectionsWeight: RefCell<HashMap<CommunityId, f32>>,
    connectionsCount: RefCell<HashMap<CommunityId, i32>>,
}

impl Community {

    fn new(id: CommunityId) -> Community {
        Community { id: id, ..Default::default() }
    }

    fn seed(&mut self, node: usize, cs: & CommunityStructure) {
        self.nodes.push(node);
        self.weightSum += cs.weights[node];
    }

    fn add(&mut self, node: usize,  cs: & CommunityStructure) -> bool {
        self.nodes.push(node);
        self.weightSum += cs.weights[node];
        true
    }

    fn remove(&mut self, node: usize, cs: &mut CommunityStructure) -> bool {
        self.nodes.remove(node);
        self.weightSum -= cs.weights[node];
        if self.nodes.is_empty() {
            cs.communities.retain(|&c| c != self.id);
        }
        true
    }


}

#[derive(Default)]
pub struct CommunityCatalog {
    lastId: CommunityId,
    map: HashMap<CommunityId, Community>
}

impl CommunityCatalog {
    pub fn createNew(&mut self) -> CommunityId {
        self.lastId += 1;
        self.map.insert(self.lastId, Community::new(self.lastId));
        self.lastId
    }
}


#[derive(Default, PartialEq)]
pub struct CommunityStructure {
    nodeConnectionsWeight: Vec<HashMap<CommunityId, f32>>,
    nodeConnectionsCount: Vec<HashMap<CommunityId, i32>>,
    nodeCommunities: Vec<CommunityId>,
    map: HashMap<String, usize>,
    weights: Vec<f32>,
    graphWeightSum: f32,
    topology: Vec<Vec<ModEdge>>,
    communities: Vec<CommunityId>,
    N: usize,
    invMap: HashMap<usize, CommunityId>,
}

impl CommunityStructure {

    pub fn new(graph: &Graph, modularity: &Modularity,
            cc: &mut CommunityCatalog) -> CommunityStructure {
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
            cs.nodeCommunities.push( cc.createNew() );

            cs.nodeConnectionsWeight.push(HashMap::new());
            cs.nodeConnectionsCount.push(HashMap::new());
            cs.weights.push(0.0);
            cc.map.get_mut(& cs.nodeCommunities[index]).unwrap().seed(index, & cs);

            let hidden = cc.createNew();
            cc.map.get_mut(& hidden).unwrap().nodes.push(index);
            cs.invMap.insert(index, hidden);

            cs.communities.push( cs.nodeCommunities[index] );
            index += 1;
        }

        unsafe {
            cs.topology.set_len(graph.nodes.len());
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
                let mut adjCom = cc.map.get(&cs.nodeCommunities[*neighbor_index]).unwrap();

                cs.nodeConnectionsWeight[*node_index].insert(adjCom.id, weight);
                cs.nodeConnectionsCount[*node_index].insert(adjCom.id, 1);

                let nodeCom = cc.map.get(&cs.nodeCommunities[*node_index]).unwrap();
                nodeCom.connectionsWeight.borrow_mut().insert(adjCom.id, weight);
                nodeCom.connectionsCount.borrow_mut().insert(adjCom.id, 1);

                cs.nodeConnectionsWeight[*neighbor_index].insert(nodeCom.id, weight);
                cs.nodeConnectionsCount[*neighbor_index].insert(nodeCom.id, 1);

                adjCom.connectionsWeight.borrow_mut().insert(nodeCom.id, weight);
                adjCom.connectionsCount.borrow_mut().insert(nodeCom.id, 1);

                cs.graphWeightSum += weight;
            }
        }

        cs.graphWeightSum /= 2.0;
        cs
    }

    fn addNodeTo(&mut self, node : usize, com: &mut Community,
                cc: &mut CommunityCatalog) {

        fn add_node<V: AddAssign + Copy>(map : &mut HashMap<i32,V>, key: i32, weight: V) {
            let w = map.entry(key).or_insert(weight);
            *w += weight;
        }

        self.nodeCommunities[node] = com.id.clone();

        com.add(node, self);

        for e in &self.topology[node] {
            let neighbor = &e.target;

            add_node(&mut self.nodeConnectionsWeight[*neighbor], com.id, e.weight);
            add_node(&mut self.nodeConnectionsCount[*neighbor], com.id, 1);

            let mut adjCom = cc.map.get_mut(& self.nodeCommunities[*neighbor]).unwrap();

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


    pub fn removeNodeFromItsCommunity(&mut self, node: usize,
            cc: &mut CommunityCatalog) {

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

        {
            let community = cc.map.get(& self.nodeCommunities[node]).unwrap();

            for e in &self.topology[node] {
                let neighbor = &e.target;

                ////////
                //Remove Node Connection to this community
                remove_node( &mut self.nodeConnectionsWeight[*neighbor],
                    &mut self.nodeConnectionsCount[*neighbor],
                    community.id.clone(), e.weight );

                ///////////////////
                //Remove Adjacency Community's connection to this community
                let adjCom = cc.map.get(& self.nodeCommunities[*neighbor]).unwrap();
                remove_node( &mut adjCom.connectionsWeight.borrow_mut(),
                    &mut adjCom.connectionsCount.borrow_mut(),
                    community.id.clone(), e.weight );

                if node == *neighbor {
                    continue;
                }

                if adjCom.id != community.id {
                    remove_node( &mut community.connectionsWeight.borrow_mut(),
                        &mut community.connectionsCount.borrow_mut(),
                        adjCom.id, e.weight);
                }

                remove_node( &mut self.nodeConnectionsWeight[node],
                    &mut self.nodeConnectionsCount[node],
                    adjCom.id, e.weight );
            }
        }
        let mut community = cc.map.get_mut(& self.nodeCommunities[node]).unwrap();
        community.remove(node, self);
    }


    fn moveNodeTo(&mut self, node: usize, to: &mut Community,
            cc: &mut CommunityCatalog) {
        self.removeNodeFromItsCommunity(node, cc);
        self.addNodeTo(node, to, cc);
    }

    fn zoomOut(&mut self,
            cc: &mut CommunityCatalog) {
        let M = self.communities.len();
        self.communities.sort();
        let mut newTopology: Vec<Vec<ModEdge>> = Vec::with_capacity(M);
        let mut index : usize = 0;
        let mut nodeCommunities: Vec<CommunityId> = Vec::with_capacity(M);
        let mut nodeConnectionsWeight: Vec<HashMap<CommunityId, f32>> = Vec::with_capacity(M);
        let mut nodeConnectionsCount: Vec<HashMap<CommunityId, i32>> = Vec::with_capacity(M);
        let mut newInvMap: HashMap<usize, CommunityId> = HashMap::new();

        for com_id in & self.communities {
            nodeConnectionsWeight.push(HashMap::new());
            nodeConnectionsCount.push(HashMap::new());

            newTopology.push( Vec::new() );
            nodeCommunities.push(cc.createNew());

            let mut weightSum : f32 = 0.0;
            let hidden_id = cc.createNew();

            /*
             * FIXME: two unnecessary clones of node vectors to content the Borrow Checker
             */
            let nodes = cc.map.get(com_id).unwrap().nodes.clone();
            for nodeInt in nodes {
                let oldHiddenNodes = {
                    let oldHidden = cc.map.get(& self.invMap.get(&nodeInt).unwrap()).unwrap();
                    oldHidden.nodes.clone()
                };
                let mut hidden = cc.map.get_mut(& hidden_id).unwrap();
                hidden.nodes.extend(oldHiddenNodes);
            }

            newInvMap.insert(index, hidden_id);
            {
                let com = cc.map.get(com_id).unwrap();
                for adjCom_id in com.connectionsWeight.borrow().keys() {
                    let target = self.communities.binary_search(adjCom_id).unwrap();
                    let weight = com.connectionsWeight.borrow().get(adjCom_id).unwrap().clone();
                    weightSum += if target == index { 2.0 * weight} else { weight };

                    let e = ModEdge {source: index, target: target, weight: weight};
                    newTopology[index].push(e);
                }
            }
            self.weights[index] = weightSum;
            let com = cc.map.get_mut(& nodeCommunities[index]).unwrap();
            com.seed(index, & self);

            index += 1;
        }
        self.communities.clear(); // FIXME: Free communities in cc too

        for i in 0..M {
            let com = cc.map.get(& nodeCommunities[i]).unwrap();
            self.communities.push(com.id);
            for e in & newTopology[i] {
                nodeConnectionsWeight[i].insert(nodeCommunities[e.target], e.weight);
                nodeConnectionsCount[i].insert(nodeCommunities[e.target], 1);
                com.connectionsWeight.borrow_mut().insert(nodeCommunities[e.target], e.weight);
                com.connectionsCount.borrow_mut().insert(nodeCommunities[e.target], 1);
            }
        }

        self.N = M;
        self.topology = newTopology;
        self.invMap = newInvMap;
    }

}


pub struct Modularity {
    pub MODULARITY_CLASS: &'static str,
//    progress: ProgressTicket,
//    structure: CommunityStructure,
    modularity: f64,
    modularityResolution: f64,
    isRandomized: bool,
    useWeight: bool,
    resolution: f64,
    cc : CommunityCatalog
}

impl Modularity {

    pub fn new() -> Modularity {
        Modularity {
            MODULARITY_CLASS: "modularity_class",
            modularity: 0.0,
            modularityResolution: 0.0,
            isRandomized: false,
            useWeight: true,
            resolution: 1.0,
            cc: Default::default(),
        }
    }

    pub fn execute(&mut self, graph: & Graph) {

        //let structure = CommunityStructure::new( & graph, & self, self.cc);
        // structure = new Modularity.CommunityStructure(graph);
        // int[] comStructure = new int[graph.getNodeCount()];
        //
        // if (graph.getNodeCount() > 0) {//Fixes issue #713 Modularity Calculation Throws Exception On Empty Graph
        //     HashMap<String, Double> computedModularityMetrics = computeModularity(graph, structure, comStructure, resolution, isRandomized, useWeight);
        //     modularity = computedModularityMetrics.get("modularity");
        //     modularityResolution = computedModularityMetrics.get("modularityResolution");
        // } else {
        //     modularity = 0;
        //     modularityResolution = 0;
        // }

        // saveValues(comStructure, graph, structure);

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
