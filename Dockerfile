# Stringmol (perf branch) — pinned build environment.
# CUDA devel base so the CPU reference and the future GPU build share one toolchain
# (g++ 13.3, CUDA 13.2). Pinned by digest: upgrades are deliberate, tested changes.
FROM nvidia/cuda:13.2.0-devel-ubuntu24.04@sha256:f9492f2eea77fbc3d0c14fa8738f35946b42da72917bf5959d284ca39b4f209a

WORKDIR /opt/stringmol
COPY src/ src/
COPY config/ config/
RUN mkdir -p release && cd src && make all

# Runs execute in a mounted directory containing run.conf; all output lands there.
WORKDIR /runs
ENTRYPOINT ["/opt/stringmol/src/smspatial", "x"]
CMD ["run.conf"]
