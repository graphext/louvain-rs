#![allow(non_snake_case)]  // FIXME: Good while porting

use std::collections::{HashSet}; // , BTreeMap as HashMap};
use fnv::{FnvHashMap as HashMap};
use std::cell::{RefCell};
use std::iter::FromIterator;
use std::ops::{AddAssign};
use rand::{thread_rng, seq};
use petgraph::{Graph as PetGraph, Undirected};
use petgraph::graph::{NodeIndex, EdgeIndex};
use petgraph::visit::EdgeRef;
use slotmap::DenseSlotMap;
use std::f32;

type Graph<T> = PetGraph<T, f32, Undirected, u32>;

trait EasyGraph {
    fn opposite(&self, node: NodeIndex, edge :EdgeIndex) -> NodeIndex;
    fn getParallelEdges(&self, node: NodeIndex, neighbor: NodeIndex) -> Vec<&f32>;
}

impl<T> EasyGraph for Graph<T> {
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


type CommunityId = slotmap::DefaultKey;

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
    map: DenseSlotMap<CommunityId, Community>
}

impl CommunityCatalog {
    pub fn createNew(&mut self) -> CommunityId {
        self.map.insert_with_key(|id| Community::new(id))
    }

    #[inline]
    pub fn get(&self, key: &CommunityId) -> Option<&Community> {
        self.map.get(*key)
    }

    #[inline]
    pub fn get_mut(&mut self, key: &CommunityId) -> Option<&mut Community> {
        self.map.get_mut(*key)
    }

    #[inline]
    pub fn remove(&mut self, key: &CommunityId) {
        self.map.remove(*key);
    }
}


#[derive(Debug, Default, PartialEq)]
pub struct CommunityStructure {
    N: usize,                                              // Number of nodes / communities
    communities: Vec<CommunityId>,                         // Direct communities
    nodeConnectionsWeight: Vec<HashMap<CommunityId, f32>>, 
    nodeConnectionsCount: Vec<HashMap<CommunityId, i32>>,  
    nodeCommunities: Vec<CommunityId>,                     
    weights: Vec<f64>,                                     // One per node
    topology: Vec<Vec<ModEdge>>,                           // One per node
    //map: HashMap<NodeIndex, usize>,
    invMap: HashMap<usize, CommunityId>,
    graphWeightSum: f64,
}

impl CommunityStructure {

    pub fn new<T>(graph: &Graph<T>, modularity: &mut Modularity) -> CommunityStructure {
        let N = graph.node_count();
        let cc = &mut modularity.cc;
        let mut cs = CommunityStructure {
            N: N,
            communities: Vec::new(),
            nodeConnectionsWeight: Vec::with_capacity(N),
            nodeConnectionsCount: Vec::with_capacity(N),
            nodeCommunities: Vec::with_capacity(N),
            weights: Vec::with_capacity(N),
            topology: Vec::with_capacity(N),
            //map: HashMap::default(),
            invMap: HashMap::default(),
            graphWeightSum: 0.0,
        };

        let mut index: usize = 0;
        // Create one community and one inverse community per node
        // All weights to 0.0
        for _node in graph.node_indices() {
            //cs.map.insert(node.clone(), index);
            cs.nodeCommunities.push( cc.createNew() );

            cs.nodeConnectionsWeight.push(HashMap::default());
            cs.nodeConnectionsCount.push(HashMap::default());
            cs.weights.push(0.0);
            cc.get_mut(& cs.nodeCommunities[index]).unwrap().seed(index, & cs);

            // New hidden community
            let hidden = cc.createNew();
            cc.get_mut(& hidden).unwrap().nodes.push(index);
            cs.invMap.insert(index, hidden);

            cs.communities.push( cs.nodeCommunities[index] );
            index += 1;
        }

        for node in graph.node_indices() {
            //let node_index = cs.map.get(&node).unwrap();

            cs.topology.push(Vec::new());

            let uniqueNeighbors : HashSet<NodeIndex> = HashSet::from_iter(graph.neighbors(node));
            for neighbor in uniqueNeighbors {
                if node == neighbor { continue }
                //let neighbor_index = cs.map.get(& neighbor).unwrap();

                //Sum all parallel edges weight:
                let mut weight : f32  = 0.0;
                for &w in graph.getParallelEdges(node, neighbor) {
                    weight += if modularity.useWeight {w} else {1.0};
                }

                //Finally add a single edge with the summed weight of all parallel edges:
                cs.weights[node.index()] += weight as f64;
                let modularity_edge = ModEdge {source: node.index(), target: neighbor.index(), weight: weight};
                cs.topology[node.index()].push(modularity_edge);
                let adjCom = cc.get(&cs.nodeCommunities[neighbor.index()]).unwrap();

                cs.nodeConnectionsWeight[node.index()].insert(adjCom.id, weight);
                cs.nodeConnectionsCount[node.index()].insert(adjCom.id, 1);

                let nodeCom = cc.get(&cs.nodeCommunities[node.index()]).unwrap();
                nodeCom.connectionsWeight.borrow_mut().insert(adjCom.id, weight);
                nodeCom.connectionsCount.borrow_mut().insert(adjCom.id, 1);

                cs.nodeConnectionsWeight[neighbor.index()].insert(nodeCom.id, weight);
                cs.nodeConnectionsCount[neighbor.index()].insert(nodeCom.id, 1);

                adjCom.connectionsWeight.borrow_mut().insert(nodeCom.id, weight);
                adjCom.connectionsCount.borrow_mut().insert(nodeCom.id, 1);

                cs.graphWeightSum += weight as f64;
            }
        }

        cs.graphWeightSum /= 2.0;
        cs
    }

    fn addNodeTo(&mut self, node : usize, com_id: CommunityId, cc: &mut CommunityCatalog) {

        #[inline]
        fn add_node<V:AddAssign + Copy + From<u8> >(map : &mut HashMap<CommunityId,V>, key: CommunityId, weight: V) {
            let w = map.entry(key).or_insert(V::from(0));
            *w += weight;
        }
        
        self.nodeCommunities[node] = com_id;

        {
            let com = cc.get_mut(& com_id).unwrap();
            com.add(node, self);
        }

        for e in &self.topology[node] {
            let neighbor = &e.target;

            add_node(&mut self.nodeConnectionsWeight[*neighbor], com_id, e.weight);
            add_node(&mut self.nodeConnectionsCount[*neighbor], com_id, 1);

            let adjCom = cc.get(& self.nodeCommunities[*neighbor]).unwrap();

            add_node(&mut adjCom.connectionsWeight.borrow_mut(), com_id, e.weight);
            add_node(&mut adjCom.connectionsCount.borrow_mut(), com_id, 1);

            add_node(&mut self.nodeConnectionsWeight[node], adjCom.id, e.weight);
            add_node(&mut self.nodeConnectionsCount[node], adjCom.id, 1);

            if com_id != adjCom.id {
                let com = cc.get(& com_id).unwrap();
                add_node(&mut com.connectionsWeight.borrow_mut(), adjCom.id, e.weight);
                add_node(&mut com.connectionsCount.borrow_mut(), adjCom.id, 1);
            }

        }


    }


    pub fn removeNodeFromItsCommunity(&mut self, node: usize, cc: &mut CommunityCatalog) {

        #[inline]
        fn remove_node( weights_map : &mut HashMap<CommunityId,f32>,
            count_map : &mut HashMap<CommunityId,i32>,
            key: CommunityId, weight: f32
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
            let community = cc.get(& self.nodeCommunities[node]).unwrap();

            for e in &self.topology[node] {
                let neighbor = &e.target;

                ////////
                //Remove Node Connection to this community
                remove_node( &mut self.nodeConnectionsWeight[*neighbor],
                    &mut self.nodeConnectionsCount[*neighbor],
                    community.id.clone(), e.weight );

                ///////////////////
                //Remove Adjacency Community's connection to this community
                let adjCom = cc.get(& self.nodeCommunities[*neighbor]).unwrap();
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
            let community = cc.get_mut(& self.nodeCommunities[node]).unwrap();
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
        let mut newInvMap: HashMap<usize, CommunityId> = HashMap::default();

        for com_id in & self.communities {
            nodeConnectionsWeight.push(HashMap::default());
            nodeConnectionsCount.push(HashMap::default());

            newTopology.push( Vec::new() );
            nodeCommunities.push(cc.createNew());

            let mut weightSum : f32 = 0.0;
            let hidden_id = cc.createNew();
            
            /*
             * FIXME: two unnecessary clones of node vectors to content the Borrow Checker
             */
            let nodes = cc.get(com_id).unwrap().nodes.clone();
            for nodeInt in nodes {
                let oldHiddenNodes = {
                    let oldHidden = cc.get(& self.invMap.get(&nodeInt).unwrap()).unwrap();
                    oldHidden.nodes.clone()
                };
                let hidden = cc.get_mut(& hidden_id).unwrap();
                hidden.nodes.extend(oldHiddenNodes);
            }

            newInvMap.insert(index, hidden_id);
            {
                let com = cc.get(com_id).unwrap();
                for adjCom_id in com.connectionsWeight.borrow().keys() {
                    let target = self.communities.binary_search(adjCom_id).unwrap();
                    let weight = com.connectionsWeight.borrow().get(adjCom_id).unwrap().clone();
                    weightSum += if target == index { 2.0 * weight} else { weight };

                    let e = ModEdge {source: index, target: target, weight: weight};
                    newTopology[index].push(e);
                }
            }
            self.weights[index] = weightSum as f64;
            let com = cc.get_mut(& nodeCommunities[index]).unwrap();
            com.seed(index, & self);

            index += 1;
        }

        for com_id in & self.communities { cc.remove(com_id); }
        self.communities.clear();

        for i in 0..M {
            let com = cc.get(& nodeCommunities[i]).unwrap();
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
    modularity: f64,
    modularityResolution: f64,
    isRandomized: bool,
    useWeight: bool,
    resolution: f64,
    noise: u32,
    cc : CommunityCatalog,
    pub communityByNode: Vec<i32>,
}

impl Modularity {
    
    ///
    /// resolution:  1.0 More resolution produce less Clusters
    /// noise:  0 clusters with number_of_nodes <= noise are considered noise (-1 tag)
    /// 
    pub fn new(resolution: f64, noise: u32) -> Modularity {
        Modularity {
            modularity: 0.0,
            modularityResolution: 0.0,
            isRandomized: false,
            useWeight: true,
            resolution: resolution,
            noise: noise,
            cc: Default::default(),
            communityByNode: Default::default()
        }
    }

    pub fn execute<T>(&mut self, graph: & Graph<T>) -> (f64, f64){

        let mut structure = CommunityStructure::new( & graph, self );
        let mut communities : Vec<i32>= vec![0; graph.node_count()];

        if graph.node_count() > 0 {
            let (modularity, modularityResolution) = self.computeModularity(graph, &mut structure, &mut communities);
            self.modularity = modularity;
            self.modularityResolution = modularityResolution;
        } else {
            self.modularity = 0.0;
            self.modularityResolution = 0.0;
        }
        // saveValues(comStructure, graph, structure);
        self.communityByNode = communities;

        (self.modularity, self.modularityResolution)
    }

    fn computeModularity<T>(&mut self, graph: & Graph<T>, cs : &mut CommunityStructure, communities: &mut Vec<i32>) -> (f64, f64) {

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
                    start = seq::sample_iter(&mut rng, 1..cs.N, 1).unwrap()[0];
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
        let mut comStructure : Vec<usize>= vec![0; graph.node_count()];
        let noiseMap = self.fillComStructure(cs, &mut comStructure);
        let degreeCount = self.fillDegreeCount(graph, cs, &mut comStructure, & nodeDegrees);

        let computedModularity = self.finalQ(&mut comStructure, & degreeCount, graph, totalWeight, 1.);
        let computedModularityResolution = self.finalQ(&mut comStructure, & degreeCount, graph, totalWeight, self.resolution);

        for i in 0..communities.len() {
            communities[i] = if noiseMap[i] { -1 } else { comStructure[i] as i32 }
        }
        // let communities: Vec<&Community> = cs.communities.iter().map(|c| self.cc.get(c).unwrap()).collect();
        // println!("{:?}", comStructure);
       
        (computedModularity, computedModularityResolution)
    }

    fn fillComStructure(&self, cs : & CommunityStructure, comStructure: &mut Vec<usize>) -> Vec<bool> {
        let mut noiseMap: Vec<bool> = vec![false; comStructure.len()];
        for (count, com_id) in cs.communities.iter().enumerate() {
            let com = self.cc.get(com_id).unwrap();
            for node in & com.nodes {
                let hidden_id = cs.invMap.get(node).unwrap();
                let hidden = self.cc.get(hidden_id).unwrap();
                let isNoise = hidden.nodes.len() as u32 <= self.noise;
                for nodeInt in & hidden.nodes {
                    comStructure[* nodeInt] = count;
                    noiseMap[* nodeInt] = isNoise;
                }
            }
        }
        noiseMap
    }

    fn fillDegreeCount<T>(&self, graph: & Graph<T>, cs : & CommunityStructure,
        comStructure: & Vec<usize>, nodeDegrees: & Vec<f64>) -> Vec<f64> {

        let mut degreeCount : Vec<f64> = vec![0.0; cs.communities.len()];

        for node in graph.node_indices() {
            
            if self.useWeight {
                degreeCount[comStructure[node.index()]] += nodeDegrees[node.index()];
            } else {
                degreeCount[comStructure[node.index()]] += graph.edges(node).count() as f64;
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

    fn q(& self, 
        node: usize, 
        com_id: & CommunityId, 
        cs: & CommunityStructure, 
        currentResolution: f64
    ) -> f64
    {
        let mut edgesTo: f64 = 0.0;
        if let Some(edgesToFloat) = cs.nodeConnectionsWeight[node].get(com_id) {
            edgesTo = *edgesToFloat as f64;
        }
        let weightSum: f64 = self.cc.get(com_id).unwrap().weightSum;
        let nodeWeight: f64 = cs.weights[node];
        let mut qValue : f64 = currentResolution * edgesTo - (nodeWeight * weightSum) / (2.0 * cs.graphWeightSum);

        let nodeCommunities_len = self.cc.get(& cs.nodeCommunities[node]).unwrap().size();
        if (cs.nodeCommunities[node] == *com_id) && (nodeCommunities_len > 1) {
            qValue = currentResolution * edgesTo - (nodeWeight * (weightSum - nodeWeight)) / (2.0 * cs.graphWeightSum);
        }
        if (cs.nodeCommunities[node] == *com_id) && (nodeCommunities_len == 1) {
            qValue = 0.0;
        }
        qValue
    }

    fn finalQ<T>(&self, comStructure: & Vec<usize>, degrees: & Vec<f64>, graph: & Graph<T>, totalWeight: f64, usedResolution: f64) -> f64 {

        let mut res: f64 = 0.0;
        let mut internal: Vec<f64> = vec![0.0; degrees.len()];

        for node in graph.node_indices() {
            for edge in graph.edges(node) { // FIXME: check if getEdges() is inner edges
                if node == edge.target() {
                    continue;
                }
                let neighbor = edge.target().index();
                if comStructure[neighbor] == comStructure[node.index()] {
                    if self.useWeight {
                        internal[comStructure[neighbor]] += * edge.weight() as f64;
                    } else {
                        internal[comStructure[neighbor]] += 1.0;
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

// AUX functions

#[allow(dead_code)]
fn check_weights(cc: &CommunityCatalog, cs: &CommunityStructure) -> bool {
    for i in 0..cs.nodeCommunities.len() - 1 {
        let com0 = cc.get(&cs.nodeCommunities[i]).unwrap();
        if let Some(w0_1) = com0.connectionsWeight.borrow().get(&cs.nodeCommunities[i+1]) {
            let nw0_1 = cs.nodeConnectionsWeight[i].get(&cs.nodeCommunities[i+1]).unwrap();
            // println!("{}: {:?} = {:?} --> {:?}", i, w0_1, nw0_1, f32::abs(w0_1 - nw0_1) < f32::EPSILON );
            if f32::abs(w0_1 - nw0_1) > f32::EPSILON {
                return false;
            }
        }
    }
    true
}