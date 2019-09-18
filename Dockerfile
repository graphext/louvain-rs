FROM rust:1.37-slim-stretch AS builder
WORKDIR /build_app/
COPY . .
RUN ./build.sh

FROM scratch
COPY --from=builder /build_app/target/x86_64-unknown-linux-musl/release/louvain /louvain
ENTRYPOINT ["/louvain"]
