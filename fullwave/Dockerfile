# Full-wave reference solver for Impedora: AWS Palace (3D FEM, Apache-2.0) + Gmsh.
# Builds natively for the host CPU (arm64 on Apple Silicon, amd64 elsewhere).
#   Built and run by:  npm run fullwave        (see fullwave/README.md)
FROM ubuntu:24.04

ARG PALACE_VERSION=v0.19.0
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential gfortran cmake git ca-certificates curl pkg-config python3 python3-gmsh \
      libopenmpi-dev openmpi-bin libopenblas-dev liblapack-dev zlib1g-dev \
    && rm -rf /var/lib/apt/lists/*

# Palace superbuild: downloads and builds its own MFEM, hypre, libCEED, SuperLU_DIST, METIS.
# ARPACK replaces SLEPc (Palace needs one eigensolver); optional accelerators are off.
# Palace's unit tests are not needed, but their CMake fetches Catch2 from GitHub at configure time,
# which can fail inside the build (it did on the first Mac build). The sed drops that one line and
# grep stops the build if a future Palace version moves it.
RUN git clone --depth 1 --branch "${PALACE_VERSION}" https://github.com/awslabs/palace.git /opt/palace-src \
    && sed -i '/add_subdirectory(..\/test\/unit/d' /opt/palace-src/palace/CMakeLists.txt \
    && ! grep -q 'test/unit' /opt/palace-src/palace/CMakeLists.txt \
    && cmake -S /opt/palace-src -B /opt/palace-build \
         -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/opt/palace \
         -DPALACE_WITH_SLEPC=OFF -DPALACE_WITH_ARPACK=ON -DPALACE_WITH_LIBXSMM=OFF -DPALACE_WITH_MAGMA=OFF \
         -DPALACE_WITH_SUNDIALS=OFF -DPALACE_WITH_GSLIB=OFF -DPALACE_WITH_SUPERLU=ON \
    && cmake --build /opt/palace-build -j"$(nproc)" \
    && rm -rf /opt/palace-src /opt/palace-build

ENV PATH=/opt/palace/bin:$PATH \
    PALACE_VERSION=${PALACE_VERSION} \
    OMPI_ALLOW_RUN_AS_ROOT=1 OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1 \
    OMPI_MCA_rmaps_base_oversubscribe=1 \
    HOME=/tmp

WORKDIR /work
