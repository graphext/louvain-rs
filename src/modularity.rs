#![allow(non_snake_case)]  // FIXME: Good while porting

use std::collections::{HashMap, HashSet};
use std::cell::{RefCell};
use std::iter::FromIterator;
use std::ops::{AddAssign};
use rand::{thread_rng, sample};
use petgraph::{Graph as PetGraph, Undirected};
use petgraph::graph::{NodeIndex, EdgeIndex};
use petgraph::visit::EdgeRef;

use io::NodeID;

type Graph = PetGraph<NodeID, f32, Undirected, u32>;

trait EasyGraph {
    fn opposite(&self, node: NodeIndex, edge :EdgeIndex) -> NodeIndex;
    fn getParallelEdges(&self, node: NodeIndex, neighbor: NodeIndex) -> Vec<&f32>;
}

impl EasyGraph for Graph {
    fn opposite(&self, node: NodeIndex, edge :EdgeIndex) -> NodeIndex {
        let (a, b) = self.edge_endpoints(edge).unwrap();
        if a == node { b } else { a }
    }

    fn getParallelEdges(&self, node: NodeIndex, neighbor: NodeIndex) -> Vec<&f32> {
        let mut weights = vec![];
        for edge in self.edges(node) {
            if edge.target() == neighbor {
                weights.push(edge.weight());
            }
        }
        weights
    }
}

#[derive(Debug, PartialEq)]
pub struct ModEdge {
    source: usize,
    target: usize,
    weight: f32,
}


type CommunityId = i32;

#[derive(Debug, Default, PartialEq)]
pub struct Community {
    id: CommunityId,
    weightSum: f64,
    nodes: Vec<usize>,
    connectionsWeight: RefCell<HashMap<CommunityId, f32>>,
    connectionsCount: RefCell<HashMap<CommunityId, i32>>,
}

impl Community {

    fn new(id: CommunityId) -> Community {
        Community { id: id, ..Default::default() }
    }

    fn size(&self) -> usize {
        self.nodes.len()
    }

    fn seed(&mut self, node: usize, cs: & CommunityStructure) {
        self.nodes.push(node);
        self.weightSum += cs.weights[node] as f64;
    }

    fn add(&mut self, node: usize,  cs: & CommunityStructure) -> bool {
        self.nodes.push(node);
        self.weightSum += cs.weights[node] as f64;
        true
    }

    fn remove(&mut self, node: usize, cs: &mut CommunityStructure) -> bool {
        self.nodes.retain(|&n| n != node);
        self.weightSum -= cs.weights[node] as f64;
        if self.nodes.is_empty() {
            cs.communities.retain(|&c| c != self.id); // FIXME: Maybe remove from cc too
        }
        true
    }


}

#[derive(Default, Debug)]
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


#[derive(Debug, Default, PartialEq)]
pub struct CommunityStructure {
    nodeConnectionsWeight: Vec<HashMap<CommunityId, f32>>,
    nodeConnectionsCount: Vec<HashMap<CommunityId, i32>>,
    nodeCommunities: Vec<CommunityId>,
    map: HashMap<NodeIndex, usize>,
    weights: Vec<f64>,
    graphWeightSum: f64,
    topology: Vec<Vec<ModEdge>>,
    communities: Vec<CommunityId>,
    N: usize,
    invMap: HashMap<usize, CommunityId>,
}

impl CommunityStructure {

    pub fn new(graph: &Graph, modularity: &mut Modularity) -> CommunityStructure {
        let N = graph.node_count();
        let mut cc = &mut modularity.cc;
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
        for node in graph.node_indices() {
            cs.map.insert(node.clone(), index);
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

        for node in graph.node_indices() {
            let node_index = cs.map.get(&node).unwrap();

            cs.topology.push(Vec::new());

            let uniqueNeighbors : HashSet<NodeIndex> = HashSet::from_iter(graph.neighbors(node));
            for neighbor in uniqueNeighbors {
                if node == neighbor { continue }
                let neighbor_index = cs.map.get(& neighbor).unwrap();

                //Sum all parallel edges weight:
                let mut weight : f32  = 0.0;
                for &w in graph.getParallelEdges(node, neighbor) {
                    weight += if modularity.useWeight {w} else {1.0};
                }

                //Finally add a single edge with the summed weight of all parallel edges:
                cs.weights[*node_index] += weight as f64;
                let modularity_edge = ModEdge {source: *node_index, target: *neighbor_index, weight: weight};
                cs.topology[*node_index].push(modularity_edge);
                let adjCom = cc.map.get(&cs.nodeCommunities[*neighbor_index]).unwrap();

                cs.nodeConnectionsWeight[*node_index].insert(adjCom.id, weight);
                cs.nodeConnectionsCount[*node_index].insert(adjCom.id, 1);

                let nodeCom = cc.map.get(&cs.nodeCommunities[*node_index]).unwrap();
                nodeCom.connectionsWeight.borrow_mut().insert(adjCom.id, weight);
                nodeCom.connectionsCount.borrow_mut().insert(adjCom.id, 1);

                cs.nodeConnectionsWeight[*neighbor_index].insert(nodeCom.id, weight);
                cs.nodeConnectionsCount[*neighbor_index].insert(nodeCom.id, 1);

                adjCom.connectionsWeight.borrow_mut().insert(nodeCom.id, weight);
                adjCom.connectionsCount.borrow_mut().insert(nodeCom.id, 1);

                cs.graphWeightSum += weight as f64;
            }
        }

        cs.graphWeightSum /= 2.0;
        cs
    }

    fn addNodeTo(&mut self, node : usize, com_id: CommunityId, cc: &mut CommunityCatalog) {

        fn add_node<V:AddAssign + Copy + From<u8> >(map : &mut HashMap<i32,V>, key: i32, weight: V) {
            let w = map.entry(key).or_insert(V::from(0));
            *w += weight;
        }
        self.nodeCommunities[node] = com_id;

        {
            let mut com = cc.map.get_mut(& com_id).unwrap();
            com.add(node, self);
        }

        for e in &self.topology[node] {
            let neighbor = &e.target;

            add_node(&mut self.nodeConnectionsWeight[*neighbor], com_id, e.weight);
            add_node(&mut self.nodeConnectionsCount[*neighbor], com_id, 1);

            let adjCom = cc.map.get(& self.nodeCommunities[*neighbor]).unwrap();

            add_node(&mut adjCom.connectionsWeight.borrow_mut(), com_id, e.weight);
            add_node(&mut adjCom.connectionsCount.borrow_mut(), com_id, 1);

            add_node(&mut self.nodeConnectionsWeight[node], adjCom.id, e.weight);
            add_node(&mut self.nodeConnectionsCount[node], adjCom.id, 1);

            if com_id != adjCom.id {
                let com = cc.map.get(& com_id).unwrap();
                add_node(&mut com.connectionsWeight.borrow_mut(), adjCom.id, e.weight);
                add_node(&mut com.connectionsCount.borrow_mut(), adjCom.id, 1);
            }

        }


    }


    pub fn removeNodeFromItsCommunity(&mut self, node: usize, cc: &mut CommunityCatalog) {

        fn remove_node( weights_map : &mut HashMap<i32,f32>,
            count_map : &mut HashMap<i32,i32>,
            key: i32, weight: f32
        )
        {

            // println!("*** remove community {} from {:?}", key, count_map);
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
        {
            let mut community = cc.map.get_mut(& self.nodeCommunities[node]).unwrap();
            community.remove(node, self);
        }
    }


    fn moveNodeTo(&mut self, node: usize, to: CommunityId, cc: &mut CommunityCatalog) {
        self.removeNodeFromItsCommunity(node, cc);
        self.addNodeTo(node, to, cc);
    }

    fn zoomOut(&mut self, cc: &mut CommunityCatalog) {
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
            self.weights[index] = weightSum as f64;
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
        self.nodeCommunities = nodeCommunities;
        self.nodeConnectionsWeight = nodeConnectionsWeight;
        self.nodeConnectionsCount = nodeConnectionsCount;
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
    cc : CommunityCatalog,
    pub communityByNode: Vec<usize>,
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
            communityByNode: Default::default()
        }
    }

    pub fn execute(&mut self, graph: & Graph) -> (f64, f64){

        let mut structure = CommunityStructure::new( & graph, self );
        let mut comStructure : Vec<usize>= vec![0; graph.node_count()];

        if graph.node_count() > 0 {
            let (modularity, modularityResolution) = self.computeModularity(graph, &mut structure, &mut comStructure);
            self.modularity = modularity;
            self.modularityResolution = modularityResolution;
        } else {
            self.modularity = 0.0;
            self.modularityResolution = 0.0;
        }
        // saveValues(comStructure, graph, structure);
        self.communityByNode = comStructure;
        (self.modularity, self.modularityResolution)
    }

    fn computeModularity(&mut self, graph: & Graph, cs : &mut CommunityStructure, comStructure: &mut Vec<usize>) -> (f64, f64) {

        let nodeDegrees = cs.weights.clone();
        let totalWeight = cs.graphWeightSum;

        let mut someChange = true;
        while someChange {
            someChange = false;
            let mut localChange = true;
            while localChange {
                localChange = false;
                let mut start: usize = 0;
                if self.isRandomized {
                    let mut rng = thread_rng();
                    start = sample(&mut rng, 1..cs.N, 1)[0];
                }
                for step in 0..cs.N {
                    let i = (step + start) % cs.N;
                    let resolution = self.resolution;
                    if let Some(bestCommunity) = self.updateBestCommunity(cs, i, resolution) {
                        if cs.nodeCommunities[i] != bestCommunity {
                            cs.moveNodeTo(i, bestCommunity, &mut self.cc);
                            localChange = true;
                        }
                    }
                }
                someChange = localChange || someChange;
                //println!("localChange: {}", localChange);
            }
            if someChange {
                cs.zoomOut(&mut self.cc);
                // println!("Zooming Out: {} communities left", cs.N);
            }
        }
        self.fillComStructure(cs, comStructure);
        let degreeCount = self.fillDegreeCount(graph, cs, comStructure, & nodeDegrees);

        let computedModularity = self.finalQ(comStructure, & degreeCount, graph, cs, totalWeight, 1.);
        let computedModularityResolution = self.finalQ(comStructure, & degreeCount, graph, cs, totalWeight, self.resolution);

        // let communities: Vec<&Community> = cs.communities.iter().map(|c| self.cc.map.get(c).unwrap()).collect();
        // println!("{:?}", comStructure);

        (computedModularity, computedModularityResolution)
    }

    fn fillComStructure(&self, cs : & CommunityStructure, comStructure: &mut Vec<usize>) {
        for (count, com_id) in cs.communities.iter().enumerate() {
            let com = self.cc.map.get(com_id).unwrap();
            for node in & com.nodes {
                let hidden_id = cs.invMap.get(node).unwrap();
                let hidden = self.cc.map.get(hidden_id).unwrap();
                for nodeInt in & hidden.nodes {
                    comStructure[* nodeInt] = count;
                }
            }
        }
    }

    fn fillDegreeCount(&self, graph: & Graph, cs : & CommunityStructure,
        comStructure: & Vec<usize>, nodeDegrees: & Vec<f64>) -> Vec<f64> {

        let mut degreeCount : Vec<f64> = vec![0.0; cs.communities.len()];

        for node in graph.node_indices() {
            let index = * cs.map.get(& node).unwrap();
            if self.useWeight {
                degreeCount[comStructure[index]] += nodeDegrees[index];
            } else {
                degreeCount[comStructure[index]] += graph.edges(node).count() as f64;
            }

        }
        degreeCount
    }

    fn updateBestCommunity(&mut self, cs: & CommunityStructure, node_id: usize, currentResolution: f64) -> Option<CommunityId> {
        let mut best = 0.0;
        let mut bestCommunity: Option<CommunityId> = None;
        // println!("{:?} !!!!{:?}", node_id, cs.nodeConnectionsWeight[node_id]);
        for com_id in cs.nodeConnectionsWeight[node_id].keys() {
            let qValue = self.q(node_id, com_id, cs, currentResolution);
            if qValue > best {
                best = qValue;
                bestCommunity = Some(*com_id);
            }
        }
        bestCommunity
    }

    fn q(& self, node: usize, com_id: & CommunityId, cs: & CommunityStructure, currentResolution: f64) -> f64{
        let mut edgesTo: f64 = 0.0;
        if let Some(edgesToFloat) = cs.nodeConnectionsWeight[node].get(com_id) {
            edgesTo = *edgesToFloat as f64;
        }
        let weightSum: f64 = self.cc.map.get(com_id).unwrap().weightSum;
        let nodeWeight: f64 = cs.weights[node];
        let mut qValue : f64 = currentResolution * edgesTo - (nodeWeight * weightSum) / (2.0 * cs.graphWeightSum);

        let nodeCommunities_len = self.cc.map.get(& cs.nodeCommunities[node]).unwrap().size();
        if (cs.nodeCommunities[node] == *com_id) && (nodeCommunities_len > 1) {
            qValue = currentResolution * edgesTo - (nodeWeight * (weightSum - nodeWeight)) / (2.0 * cs.graphWeightSum);
        }
        if (cs.nodeCommunities[node] == *com_id) && (nodeCommunities_len == 1) {
            qValue = 0.0;
        }
        qValue
    }

    fn finalQ(&self, comStructure: & Vec<usize>, degrees: & Vec<f64>, graph: & Graph,
         cs: & CommunityStructure, totalWeight: f64, usedResolution: f64) -> f64 {

        let mut res: f64 = 0.0;
        let mut internal: Vec<f64> = vec![0.0; degrees.len()];

        for node in graph.node_indices() {
            let n_index = * cs.map.get(& node).unwrap();
            for edge in graph.edges(node) { // FIXME: check if getEdges() is inner edges
                if node == edge.target() {
                    continue;
                }
                let neigh_index = * cs.map.get(& edge.target()).unwrap();
                if comStructure[neigh_index] == comStructure[n_index] {
                    if self.useWeight {
                        internal[comStructure[neigh_index]] += * edge.weight() as f64;
                    } else {
                        internal[comStructure[neigh_index]] += 1.0;
                    }
                }
            }
        }
        for i in 0..degrees.len() {
            internal[i] /= 2.0;
            res += usedResolution * (internal[i] / totalWeight) - (degrees[i] / (2.0 * totalWeight)).powi(2);
        }
        res
    }

}




#[cfg(test)]
mod tests {
    #[test]
    fn it_works() {
    }
}
