#
# This Dockerfile acts as the louvain executable, don't use this to build the library
#
# If you want to build the library use the same image "rust:1.42-slim-stretch"
#   and run `cargo build --release`
#
FROM rust:1.42-alpine3.11 AS builder
WORKDIR /build_app/
COPY . .
RUN apk --update add musl-dev && ./build.sh

FROM scratch
COPY --from=builder /build_app/target/release/louvain /louvain
COPY --from=builder /build_app/target/release/liblouvain.a /lib/
COPY --from=builder /build_app/include /lib/
ENTRYPOINT ["/louvain"]
