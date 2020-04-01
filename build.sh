# This script builds louvain with statical linking using musl instead of gnu's libc
# Execute this in the image: "rust:1.42-slim-stretch"
#   docker run -it --rm -v "`pwd`":/hostdata:cached rust:1.42-slim-stretch sh
#
rustup target add x86_64-unknown-linux-musl
cargo build --release --target=x86_64-unknown-linux-musl
echo "---------------------------------------------------------------------------"
echo "The output is in target/x86_64-unknown-linux-musl/release/louvain"
echo "---------------------------------------------------------------------------"
