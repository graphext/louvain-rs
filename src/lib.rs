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
type CommunityStructureId = i32;

#[derive(Debug, Default, PartialEq)]
struct Community {
    id: CommunityId,
    weightSum: f32,
    structure: CommunityStructureId,
    nodes: Vec<usize>,
    connectionsWeight: RefCell<HashMap<CommunityId, f32>>,
    connectionsCount: RefCell<HashMap<CommunityId, i32>>,
}

impl Community {

    fn new(id: CommunityId, cs_id: CommunityStructureId) -> Community {
        Community { id: id, structure: cs_id, ..Default::default() }
    }

    fn seed(&mut self, node: usize, csc: & CommunityStructureCatalog) {
        self.nodes.push(node);
        let cs = csc.map.get(&self.structure).unwrap();
        self.weightSum += cs.weights[node];
    }

    fn add(&mut self, node: usize, csc: & CommunityStructureCatalog) -> bool {
        self.nodes.push(node);
        let cs = csc.map.get(&self.structure).unwrap();
        self.weightSum += cs.weights[node];
        true
    }

    fn remove(&mut self, node: usize, csc: &mut CommunityStructureCatalog) -> bool {
        self.nodes.remove(node);
        let mut cs = csc.map.get_mut(&self.structure).unwrap();
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
    fn createNew(&mut self, cs_id: CommunityStructureId) -> CommunityId {
        self.lastId += 1;
        self.map.insert(self.lastId, Community::new(self.lastId, cs_id));
        self.lastId
    }
}

#[derive(Default)]
pub struct CommunityStructureCatalog {
    lastId: CommunityStructureId,
    map: HashMap<CommunityStructureId, CommunityStructure>
}

impl CommunityStructureCatalog {
    fn createNew(&mut self, graph: &Graph, modularity: &Modularity,
            cc: &mut CommunityCatalog) -> CommunityStructureId {
        self.lastId += 1;
        let cs = CommunityStructure::new(self.lastId, graph, modularity, cc, self);
        self.map.insert(self.lastId, cs);
        self.lastId
    }
}


#[derive(Default, PartialEq)]
pub struct CommunityStructure {
    id: CommunityStructureId,
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

    pub fn new(cs_id: CommunityStructureId, graph: &Graph, modularity: &Modularity,
            cc: &mut CommunityCatalog,
            csc: &mut CommunityStructureCatalog) -> CommunityStructure {
        let N = graph.nodes.len();
        let mut cs = CommunityStructure {
            id: cs_id,
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
            cs.nodeCommunities.push( cc.createNew(cs_id) );

            cs.nodeConnectionsWeight.push(HashMap::new());
            cs.nodeConnectionsCount.push(HashMap::new());
            cs.weights.push(0.0);
            cc.map.get_mut(& cs.nodeCommunities[index]).unwrap().seed(index, csc);

            let hidden = cc.createNew(cs_id);
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
                cc: &mut CommunityCatalog,
                csc: &mut CommunityStructureCatalog) {

        fn add_node<V: AddAssign + Copy>(map : &mut HashMap<i32,V>, key: i32, weight: V) {
            let w = map.entry(key).or_insert(weight);
            *w += weight;
        }

        self.nodeCommunities[node] = com.id.clone();

        com.add(node, csc);

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
            cc: &mut CommunityCatalog,
            csc: &mut CommunityStructureCatalog) {

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
        community.remove(node, csc);
    }


    fn moveNodeTo(&mut self, node: usize, to: &mut Community,
            cc: &mut CommunityCatalog,
            csc: &mut CommunityStructureCatalog) {
        self.removeNodeFromItsCommunity(node, cc, csc);
        self.addNodeTo(node, to, cc, csc);
    }

    fn zoomOut(&mut self,
            cc: &mut CommunityCatalog,
            csc: &mut CommunityStructureCatalog) {
        let M = self.communities.len();
        let mut newTopology: Vec<Vec<ModEdge>> = Vec::with_capacity(M);
        let mut index : usize = 0;
        let mut nodeCommunities: Vec<CommunityId> = Vec::with_capacity(M);
        let mut nodeConnectionsWeight: Vec<HashMap<CommunityId, f32>> = Vec::with_capacity(M);
        let mut nodeConnectionsCount: Vec<HashMap<CommunityId, i32>> = Vec::with_capacity(M);
        let mut newInvMap: HashMap<usize, CommunityId> = HashMap::new();

        for (i, com) in self.communities.iter().enumerate() {
            nodeConnectionsWeight.push(HashMap::new());
            nodeConnectionsCount.push(HashMap::new());

            newTopology.push( Vec::new() );
            nodeCommunities.push(cc.createNew(self.id));

            // Set<Community> iter = com.connectionsWeight.keySet();
            // double weightSum = 0;
            //
            // Community hidden = new Community(structure);
            // for (Integer nodeInt : com.nodes) {
            //     Community oldHidden = invMap.get(nodeInt);
            //     hidden.nodes.addAll(oldHidden.nodes);
            // }

        }
    }

}


pub struct Modularity {
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

    pub fn new() -> Modularity {
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
