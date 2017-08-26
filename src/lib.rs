#[macro_use]
extern crate serde_derive;
extern crate serde;
extern crate serde_json;

extern crate rand;
extern crate petgraph;

pub mod modularity;
pub use self::modularity::Modularity;

pub mod io;
