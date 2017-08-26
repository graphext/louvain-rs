## The total memory is fixed... must be multiple of 16MB
EMMAKEN_CFLAGS="-s TOTAL_MEMORY=167772160" cargo web start --host 0.0.0.0 --release

## This flag truns down some optimizations
# EMMAKEN_CFLAGS="-s ALLOW_MEMORY_GROWTH=1" cargo web start --host 0.0.0.0 --release           
