# This script builds louvain with statical linking using musl instead of gnu's libc
# Execute this in the image: "rust:1.42-slim-stretch"
#   docker run -it --rm -v "`pwd`":/hostdata:cached rust:1.42-slim-stretch sh
#
RUSTFLAGS="-C target-feature=-crt-static" cargo build --release
echo "---------------------------------------------------------------------------"
echo "The output is in target/release/louvain"
echo "---------------------------------------------------------------------------"
