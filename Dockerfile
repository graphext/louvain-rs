#
# This Dockerfile acts as the louvain executable, don't use this to build the library
#  
# If you want to build the library use the same image "rust:1.42-slim-stretch"
#   and run `cargo build --release`
#
FROM rust:1.42-slim-stretch AS builder
WORKDIR /build_app/
COPY . .
RUN ./build.sh

FROM scratch
COPY --from=builder /build_app/target/x86_64-unknown-linux-musl/release/louvain /louvain
ENTRYPOINT ["/louvain"]
