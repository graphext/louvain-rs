# Rust implementation of Louvain

Louvain is a clustering algorithm that can run on big networks efficiently.

~~Nightly is needed.~~ Now use DenseSlotMap instead of SlotMap... Some further profiling is needed to measure the performance regression.

## How to compile it ?

You can use the docker image "rust:1.42-slim-stretch", has everything installed.

Run in the shell:

```sh
cargo build --release
```

Cargo will download all dependencies and generate in the target directory some assets:

- `target/release/louvain`: An executable that reads nodes.json and links.json
    Type `louvain --help` to show the documentation for all the params

- `target/release/liblouvain.a` A C library that expose only one function. Use it in combination with `include/louvain.h`