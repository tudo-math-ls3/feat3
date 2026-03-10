#!/usr/bin/env bash

set -e # Exit on error
set -o pipefail # Exit if any part of pipeline fails

USE_OPENMP=OFF
USE_OPENMP_NUM=0
USE_OPENMP_FLAG=""

NP=1
BUILD_TESTS=OFF

POSITIONAL_ARGS=()

while [[ $# -gt 0 ]]; do
  case $1 in
    -h|--help)
      echo "Install script for Trilinos and its dependencies."
      echo "Supported options:"
      echo ""
      echo "-w <PATH>, --working-dir <PATH>"
      echo -e "\tWorking directory for this script. Source and build directories will be created in the working directory. Required."
      echo "-i <PATH>, --install-dir <PATH>"
      echo -e "\tInstall directory for this script. Dependencies will be installed into this directory. Required."
      echo "-t <PATH>, --trilinos-dir <PATH>"
      echo -e "\tPath to root of trilinos sources."
      echo "-o,--openmp"
      echo -e "\tEnable OpenMP support for dependencies and Trilinos."
      echo "-j"
      echo -e "\tSet number of build processes"
      exit 0
      ;;
    -w|--working-dir)
      WORKING_DIR=$(realpath "$2")
      shift # past argument
      shift # past value
      ;;
    -i|--install-dir)
      INSTALL_DIR=$(realpath "$2")
      shift # past argument
      shift # past value
      ;;
    -t|--trilinos-dir)
      TRILINOS_DIR=$(realpath "$2")
      shift # past argument
      shift # past value
      ;;
    -j)
      NP="$2"
      shift
      shift
      ;;
    -o|--openmp)
      USE_OPENMP=ON
      USE_OPENMP_NUM=1
      USE_OPENMP_FLAG=-fopenmp
      shift
      ;;
    -*)
      echo "Unknown option $1"
      exit 1
      ;;
    *)
      POSITIONAL_ARGS+=("$1") # save positional arg
      shift # past argument
      ;;
  esac
done

set -- "${POSITIONAL_ARGS[@]}" # restore positional parameters

BUILD_ID=$1
hasTag() { case $BUILD_ID in *$1* ) return 0;; *) return 1;; esac ;}

BAD_ARGS=false

if [ -z "$WORKING_DIR" ]; then echo "No working directory set"; BAD_ARGS=true; fi
if [ -z "$INSTALL_DIR" ]; then echo "No install directory set"; BAD_ARGS=true; fi

if [ $BAD_ARGS = true ]; then echo "Missing required parameters. Aborting!"; exit 1; fi


# Setup folder structure
mkdir -p "$WORKING_DIR"
mkdir -p "$WORKING_DIR"/logs
mkdir -p "$INSTALL_DIR"

cd "$WORKING_DIR"

export CMAKE_PREFIX_PATH=$INSTALL_DIR
export CMAKE_INSTALL_PREFIX=$INSTALL_DIR
export CC=gcc
export prefix=$INSTALL_DIR
export CFLAGS=-O3

if hasTag "trilinos" || hasTag "umfpack" || hasTag "mumps" || hasTag "scalapack"; then
    if [[ ! -d "$WORKING_DIR"/openblas-src ]]; then
        echo "[OpenBLAS] Download step..."
        git clone --depth=1 --branch  v0.3.31 https://github.com/OpenMathLib/OpenBLAS openblas-src
    else
        echo "[OpenBLAS] Skipping download step: $WORKING_DIR/openblas-src already exists"
    fi

    if [[ ! -d "$WORKING_DIR"/openblas-build ]]; then
        echo "[OpenBLAS] Build step..."

        mkdir -p openblas-build
        cd openblas-build
        cp -r "$WORKING_DIR"/openblas-src/* .

        make  USE_OPENMP=$USE_OPENMP_NUM NO_SHARED=1 -j $NP 2>&1 | tee "$WORKING_DIR"/logs/openblas-build.log


        cd ..
    else
        echo "[OpenBLAS] Skipping build step: $WORKING_DIR/openblas-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib/libopenblas.a ]; then
        echo "[OpenBLAS] Install step..."

        cd "$WORKING_DIR"/openblas-build
        make install USE_OPENMP=$USE_OPENMP_NUM NO_SHARED=1 PREFIX="$INSTALL_DIR" 2>&1 | tee "$WORKING_DIR"/logs/openblas-install.log
        cd ..
    else
        echo "[OpenBLAS] Skipping install step: $INSTALL_DIR/lib/libopenblas.a already exists"
    fi
fi

if hasTag "boost"; then
    if [ ! -d "$WORKING_DIR"/boost-src ]; then
        echo "[Boost] Download step..."
        git clone --depth=1 --branch boost-1.88.0 --recursive https://github.com/boostorg/boost.git boost-src
    else
        echo "[Boost] Skipping download step: $WORKING_DIR/boost-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/boost-build ]; then
        echo "[Boost] Build step..."
        mkdir -p boost-build
        cd boost-build
        cmake ../boost-src \
            -D CMAKE_BUILD_TYPE="Release" 2>&1 | tee "$WORKING_DIR"/logs/boost-configure.log

        cmake --build . -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/boost-build.log
        cd ..
    else
        echo "[Boost] Skipping build step: $WORKING_DIR/boost-build already exists"
    fi

    if [ ! -d "$INSTALL_DIR"/include/boost ]; then
        cd boost-build
        cmake --install . --prefix="$INSTALL_DIR" 2>&1 | tee "$WORKING_DIR"/logs/boost-install-log
        cd ..
    else
        echo "[Boost] Skipping install step: $INSTALL_DIR/include/boost already exists"
    fi
fi

if hasTag "cgal"; then
    if [ ! -d "$WORKING_DIR"/cgal-src ]; then
        echo "[CGAL] Download step..."
        git clone --depth=1 --branch v6.0.1 --recursive https://github.com/CGAL/cgal.git cgal-src
    else
        echo "[CGAL] Skipping download step: $WORKING_DIR/cgal-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/cgal-build ]; then
        echo "[CGAL] Build step..."
        mkdir -p cgal-build
        cd cgal-build
        cmake ../cgal-src \
            -D CMAKE_BUILD_TYPE=Release \
            -D EIGEN3_INCLUDE_DIR="/usr/include/eigen3/" 2>&1 | tee "$WORKING_DIR"/logs/cgal-configure.log

        cmake --build . -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/cgal-build.log
        cd ..
    else
        echo "[CGAL] Skipping build step: $WORKING_DIR/cgal-build already exists"
    fi

    if [ ! -d "$INSTALL_DIR"/include/CGAL ]; then
        echo "[CGAL] Install step..."
        cd cgal-build
        cmake --install . --prefix="$INSTALL_DIR" 2>&1 | tee "$WORKING_DIR"/logs/cgal-install.log
        cd ..
    else
        echo "[CGAL] Skipping install step: $INSTALL_DIR/include/CGAL already exists"
    fi
fi

if hasTag "pmp"; then
    if [ ! -d "$WORKING_DIR"/pmp-src ]; then
        echo "[PMP] Download step..."
        git clone --depth=1 --branch 3.0.0 --recursive https://github.com/pmp-library/pmp-library.git pmp-src
        cd pmp-src
        # Patch pmp to correct version
        sed -i -e 's/VERSION 2.0.1/VERSION 3.0.0/' CMakeLists.txt
        cd ..
    else
        echo "[PMP] Skipping download step: $WORKING_DIR/pmp-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/pmp-build ]; then
        echo "[PMP] Build step..."
        mkdir -p pmp-build
        cd pmp-build
        cmake ../pmp-src \
            -D CMAKE_BUILD_TYPE=Release \
            -D PMP_BUILD_EXAMPLES=OFF \
            -D PMP_BUILD_TESTS=OFF \
            -D PMP_BUILD_DOCS=OFF \
            -D PMP_BUILD_VIS=OFF \
            -D BUILD_SHARED_LIBS=OFF \
            2>&1 | tee "$WORKING_DIR"/logs/pmp-configure.log

        cmake --build . -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/pmp-build.log
        cd ..
    else
        echo "[PMP] Skipping build step: $WORKING_DIR/pmp-build already exists"
    fi

    if [ ! -d "$INSTALL_DIR"/include/pmp ]; then
        echo "[PMP] Install step..."
        cd pmp-build
        cmake --install . --prefix="$INSTALL_DIR" 2>&1 | tee "$WORKING_DIR"/logs/pmp-install.log
        cd ..
    else
        echo "[PMP] Skipping install step: $INSTALL_DIR/include/pmp already exists"
    fi
fi

if hasTag "meshhexer"; then
    if [ ! -d "$WORKING_DIR"/meshhexer-src ]; then
        echo "[MesHexer] Download step..."
        git clone https://github.com/tudo-math-ls3/MeshHexer.git meshhexer-src
    else
        echo "[MesHexer] Skipping download step: $WORKING_DIR/meshhexer-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/meshhexer-build ]; then
        echo "[MesHexer] Build step..."
        mkdir -p meshhexer-build
        cd meshhexer-build
        cmake ../meshhexer-src \
            -D CMAKE_BUILD_TYPE=Release \
            -D CMAKE_PREFIX_PATH="$INSTALL_DIR" \
            -D MESHHEXER_ALLOW_EXTERNAL_DOWNLOAD=OFF \
            -D MESHHEXER_BUILD_APPLICATIONS=OFF \
            -D MESHHEXER_DOCUMENTATION=OFF \
            -D MESHHEXER_HAVE_OMP=${USE_OPENMP} \
            -D MESHHEXER_TESTING=OFF \
            2>&1 | tee "$WORKING_DIR"/logs/meshhexer-configure.log

        cmake --build . -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/meshhexer-build.log
        cd ..
    else
        echo "[MesHexer] Skipping build step: $WORKING_DIR/meshhexer-build already exists"
    fi

    if [ ! -d "$INSTALL_DIR"/include/meshhexer ]; then
        echo "[MeshHexer] Install step..."
        cd meshhexer-build
        cmake --install . --prefix="$INSTALL_DIR" 2>&1 | tee "$WORKING_DIR"/logs/meshhexer-install.log
        cd ..
    else
        echo "[MesHexer] Skipping install step: $INSTALL_DIR/include/meshhexer already exists"
    fi
fi

if hasTag "hypre"; then
    if [ ! -d "$WORKING_DIR"/hypre-src ]; then
        echo "[Hypre] Download step..."
        git clone --branch v2.33.0 --depth=1 https://github.com/hypre-space/hypre.git hypre-src
    else
        echo "[Hypre] Skipping download step: $WORKING_DIR/hypre-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/hypre-build ]; then
        echo "[Hypre] Build step..."
        mkdir -p hypre-build
        cd hypre-build
        cmake ../hypre-src/src \
            -D CMAKE_BUILD_TYPE=Release \
            -D HYPRE_ENABLE_PERSISTENT_COMM=ON \
            -D HYPRE_USE_PERSISTENT_COMM=ON \
            -D HYPRE_ENABLE_OPENMP=${USE_OPENMP} \
            -D HYPRE_ENABLE_LTO=ON \
            2>&1 | tee "$WORKING_DIR"/logs/hypre-configure.log

        cmake --build . -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/meshhexer-build.log
        cd ..
    else
        echo "[Hypre] Skipping build step: $WORKING_DIR/hypre-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib64/libHYPRE.a ]; then
        echo "[Hypre] Install step..."
        cd hypre-build
        cmake --install . --prefix="$INSTALL_DIR" 2>&1 | tee "$WORKING_DIR"/logs/hypre-install.log
        cd ..
    else
        echo "[Hypre] Skipping install step: $INSTALL_DIR/lib64/libHYPRE.a already exists"
    fi
fi

if hasTag "alglib"; then
    if [ ! -d "$WORKING_DIR"/alglib-src ]; then
        echo "[Alglib] Download step..."
        wget http://www.alglib.net/translator/re/alglib-4.04.0.cpp.gpl.zip
        unzip alglib-4.04.0.cpp.gpl.zip
        rm alglib-4.04.0.cpp.gpl.zip
        mv alglib-cpp alglib-src

        cd alglib-src

        # Add build system
        cat > src/CMakeLists.txt <<'EOF'
cmake_minimum_required(VERSION 3.12)

project(Alglib VERSION 1.0 LANGUAGES CXX)

set(CMAKE_BUILD_TYPE Release)

set(HEADERS
alglibinternal.h
alglibmisc.h
ap.h
dataanalysis.h
diffequations.h
fasttransforms.h
integration.h
interpolation.h
kernels_avx2.h
kernels_fma.h
kernels_sse2.h
linalg.h
optimization.h
solvers.h
specialfunctions.h
statistics.h
stdafx.h
)

set(SOURCES
alglibinternal.cpp
alglibmisc.cpp
ap.cpp
dataanalysis.cpp
diffequations.cpp
fasttransforms.cpp
integration.cpp
interpolation.cpp
kernels_avx2.cpp
kernels_fma.cpp
kernels_sse2.cpp
linalg.cpp
optimization.cpp
solvers.cpp
specialfunctions.cpp
statistics.cpp
)

add_library(alglib STATIC)
target_sources(alglib PRIVATE ${SOURCES})
target_sources(alglib PUBLIC FILE_SET HEADERS FILES ${HEADERS})

include(GNUInstallDirs)

install(TARGETS alglib FILE_SET HEADERS)
EOF
        cd ..
    else
        echo "[Alglib] Skipping download step: $WORKING_DIR/hypre-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/alglib-build ]; then
        echo "[Alglib] Build step..."
        mkdir -p alglib-build
        cd alglib-build
        cmake ../alglib-src/src \
            -D CMAKE_BUILD_TYPE=Release \
            2>&1 | tee "$WORKING_DIR"/logs/alglib-configure.log

        cmake --build . -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/alglib-build.log
        cd ..
    else
        echo "[Alglib] Skipping build step: $WORKING_DIR/alglib-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib64/libalglib.a ]; then
        echo "[Alglib] Install step..."
        cd alglib-build
        cmake --install . --prefix="$INSTALL_DIR" 2>&1 | tee "$WORKING_DIR"/logs/alglib-install.log
        cd ..
    else
        echo "[Alglib] Skipping install step: $INSTALL_DIR/lib64/libalglib.a already exists"
    fi
fi

if hasTag "cudss"; then
    if [ ! -d "$WORKING_DIR"/cudss-src ]; then
        echo "[cuDSS] Download step..."
        wget https://developer.download.nvidia.com/compute/cudss/redist/libcudss/linux-x86_64/libcudss-linux-x86_64-0.7.1.4_cuda12-archive.tar.xz
        tar xvf libcudss-linux-x86_64-0.7.1.4_cuda12-archive.tar.xz
        mv libcudss-linux-x86_64-0.7.1.4_cuda12-archive cudss-src
    else
        echo "[cuDSS] Skipping download step: $WORKING_DIR/cudss-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/cudss-build ]; then
        echo "[cuDSS] Build step..."
        mkdir -p cudss-build
        # Nothing to do. cuDSS is distributed as a binary blob
    else
        echo "[cuDSS] Skipping build step: $WORKING_DIR/cudss-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib/libcudss_static.a ]; then
        echo "[cuDSS] Install step..."
        cd cudss-src
        cp -ra ./include/* "$INSTALL_DIR"/include
        cp -ra ./lib/* "$INSTALL_DIR"/lib
        cd ..
    else
        echo "[cuDSS] Skipping install step: $INSTALL_DIR/lib/libcudss_static.a already exists"
    fi
fi

if hasTag "floatx"; then
    if [ ! -d "$WORKING_DIR"/floatx-src ]; then
        echo "[FloatX] Download step..."
        git clone https://github.com/oprecomp/FloatX.git floatx-src
    else
        echo "[FloatX] Skipping download step: $WORKING_DIR/floatx-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/floatx-build ]; then
        echo "[FloatX] Build step..."
        mkdir -p floatx-build
        # Nothing to do. FloatX is header only
    else
        echo "[FloatX] Skipping build step: $WORKING_DIR/floatx-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/include/floatx.hpp ]; then
        echo "[FloatX] Install step..."
        cd floatx-src
        cp -ra ./src/floatx.hpp $INSTALL_DIR/include/floatx.hpp
        cd ..
    else
        echo "[FloatX] Skipping install step: $INSTALL_DIR/include/floatx.hpp already exists"
    fi
fi

if hasTag "umfpack" || hasTag "trilinos"; then
    if [[ ! -d "$WORKING_DIR"/suitesparse-src ]]; then
        echo "[SuiteSparse] Download step..."
        git clone https://github.com/DrTimothyAldenDavis/SuiteSparse suitesparse-src
        cd suitesparse-src
        git checkout 8ac3f515ad91ae3d0137fe98239e52d2a689eac3
        cd ..
    else
        echo "[SuiteSparse] Skipping download step: $WORKING_DIR/suitesparse-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/suitesparse-build ]; then
        echo "[SuiteSparse] Build step..."
        mkdir -p suitesparse-build
        cd suitesparse-build

        # NOTE: CMakes checks for valid BLAS libraries by looking for the sgemm_ symbol.
        # For some reason our openblas.a fails that check, despite containing a sgemm_ symbol.
        # Directly setting the BLAS via BLAS_LIBRARIES circumvents this check.
        # Start looking here if any UMFPACK/BLAS weirdness starts turning up
        cmake \
        -D BLA_STATIC:BOOL=ON \
        -D BLAS_LIBRARIES:PATH="$INSTALL_DIR/lib/libopenblas.a" \
        -D LAPACK_LIBRARIES:PATH="$INSTALL_DIR/lib/libopenblas.a" \
        -D BUILD_SHARED_LIBS:BOOL=OFF \
        -D CMAKE_BUILD_TYPE:STRING=Release \
        -D CMAKE_PREFIX_PATH:PATH="$INSTALL_DIR" \
        -D SUITESPARSE_ENABLE_PROJECTS:STRING="umfpack;amd;suitesparse_config;cholmod" \
        -D SUITESPARSE_USE_FORTRAN:BOOL=OFF \
        -D SUITESPARSE_USE_OPENMP:BOOL=${USE_OPENMP} \
        ../suitesparse-src 2>&1 | tee "$WORKING_DIR"/logs/suitesparse-configure.log

        cmake --build . -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/suitesparse-build.log

        cd ..
    else
        echo "[SuiteSparse] Skipping build step: $WORKING_DIR/suitesparse-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib64/libumfpack.a ]; then
        echo "[SuiteSparse] Install step..."
        cd suitesparse-build
        cmake --install . --prefix "$INSTALL_DIR" 2>&1 | tee "$WORKING_DIR"/logs/suitesparse-install.log
        cd ..
    else
        echo "[SuiteSparse] Skipping install step: $INSTALL_DIR/lib64/libumfpack.a already exists"
    fi
fi

if hasTag "zfp"; then
    if [ ! -d "$WORKING_DIR"/zfp-src ]; then
        echo "[ZFP] Download step..."
        git clone --depth=1 --branch 1.0.0 https://github.com/LLNL/zfp zfp-src
    else
        echo "[ZFP] Skipping install step: $WORKING_DIR/zfp-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/zfp-build ]; then
        echo "[ZFP] Build step..."
        mkdir -p zfp-build
        cd zfp-build

        cmake ../zfp-src \
            -D CMAKE_BUILD_TYPE=Release \
            -D CMAKE_POSITION_INDEPENDENT_CODE=ON \
            -D BUILD_SHARED_LIBS=OFF 2>&1 | tee "$WORKING_DIR"/logs/zfp-configure.log

        cmake --build . -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/zfp-build.log

        cd ..
    else
        echo "[ZFP] Skipping build step: $WORKING_DIR/zfp-build already exists"
    fi

    if [ ! -d "$INSTALL_DIR"/include/zfp ]; then
        cd zfp-build
        cmake --install . --prefix="$INSTALL_DIR"
        cd ..
    else
        echo "[ZFP] Skipping install step: $INSTALL_DIR/include/zfp already exists"
    fi
fi

if hasTag "zlib" || hasTag "trilinos"; then
    if [[ ! -d "$WORKING_DIR"/zlib-src ]]; then
        echo "[Zlib] Download step..."
        git clone --depth 1 --branch v1.3.1 https://github.com/madler/zlib zlib-src
    else
        echo "[Zlib] Skipping download step: $WORKING_DIR/zlib-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/zlib-build ]; then
        echo "[Zlib] Build step..."
        mkdir -p zlib-build
        cd zlib-build

        ../zlib-src/configure --static 2>&1 | tee "$WORKING_DIR"/logs/zlib-configure.log
        make -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/zlib-build.log

        cd ..
    else
        echo "[Zlib] Skipping build step: $WORKING_DIR/zlib-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib/libz.a ]; then
        echo "[Zlib] Install step..."
        cd zlib-build

        make install 2>&1 | tee "$WORKING_DIR"/logs/zlib-install.log

        cd ..
    else
        echo "[Zlib] Skipping install step: $INSTALL_DIR/lib/libz.a already exists"
    fi
fi

if hasTag "hdf5" || hasTag "trilinos"; then
    if [[ ! -d "$WORKING_DIR"/hdf5-src ]]; then
        echo "[HDF5] Download step..."
        git clone --depth 1 --branch hdf5-1.14.6 https://github.com/HDFGroup/hdf5 hdf5-src
    else
        echo "[HDF5] Skipping download step: $WORKING_DIR/hdf5-src already exists"
    fi

    if [[ ! -d $WORKING_DIR/hdf5-build ]]; then
        echo "[HDF5] Build step..."
        mkdir -p hdf5-build
        cd hdf5-build

        CC=mpicc ../hdf5-src/configure \
            --enable-parallel \
            --prefix="$INSTALL_DIR" \
            --with-zlib="$INSTALL_DIR"/include,"$INSTALL_DIR"/lib \
            --disable-shared \
            --enable-tests="" \
            --enable-build-mode=production 2>&1 | tee "$WORKING_DIR"/logs/hdf5-configure.log

        make -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/hdf5-build.log
        cd ..
    else
        echo "[HDF5] Skipping build step: $WORKING_DIR/hdf5-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib/libhdf5.a ]; then
        echo "[HDF5] Install step..."
        cd hdf5-build

        make install 2>&1 | tee "$WORKING_DIR"/logs/hdf5-install.log

        cd ..
    else
        echo "[HDF5] Skipping install step: $INSTALL_DIR/lib/libhdf5.a already exists"
    fi
fi

if hasTag "metis" || hasTag "trilinos"; then
    if [[ ! -d "$WORKING_DIR"/parmetis-src ]]; then
        echo "[ParMETIS] Download step..."
        git clone https://github.com/tudo-math-ls3/ParMETIS.git parmetis-src
    else
        echo "[ParMETIS] Skipping download step: $WORKING_DIR/parmetis-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/parmetis-build ]; then
        echo "[ParMETIS] Build step..."
        mkdir -p parmetis-build
        cd parmetis-build

        cmake ../parmetis-src -DINDEXTYPEWIDTH=32 -DREALTYPEWIDTH=64 2>&1 | tee "$WORKING_DIR"/logs/parmetis-configure.log
        cmake --build . -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/parmetis-build.log

        cd ..
    else
        echo "[ParMETIS] Skipping build step: $WORKING_DIR/parmetis-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib64/libparmetis.a ]; then
        echo "[ParMETIS] Install step..."

        cd parmetis-build
        cmake --install . --prefix="$INSTALL_DIR" 2>&1 | tee "$WORKING_DIR"/logs/parmetis-install.log
        cd ..
    else
        echo "[ParMETIS] Skipping install step: $INSTALL_DIR/lib64/libparmetis.a already exists"
    fi
fi

if hasTag "superlu"; then
    if [ ! -d "$WORKING_DIR"/superludist-src ]; then
        echo "[SuperLU_DIST] Download step..."
        git clone --branch v9.1.0 --depth=1 https://github.com/xiaoyeli/superlu_dist.git superludist-src
    else
        echo "[SuperLU_DIST] Skipping download step: $WORKING_DIR/superludist-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/superludist-build ]; then
        echo "[SuperLU_DIST] Build step..."
        mkdir -p superludist-build

        cd superludist-build
        cmake ../superludist-src \
            -D TPL_PARMETIS_INCLUDE_DIRS="$INSTALL_DIR/include" \
            -D TPL_PARMETIS_LIBRARIES="$INSTALL_DIR/lib64/libparmetis.a;$INSTALL_DIR/lib64/libmetis.a;$INSTALL_DIR/lib64/libgklib.a" \
            -D TPL_BLAS_INCLUDE_DIRS="$INSTALL_DIR/include" \
            -D TPL_BLAS_LIBRARIES="$INSTALL_DIR/lib/libopenblas.a" \
            -D BUILD_SHARED_LIBS=OFF \
            -D CMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
            -D XSDK_INDEX_SIZE=32 \
            2>&1 | tee "$WORKING_DIR"/logs/superludist-configure.log

        cmake --build . -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/superludist-build.log
        cd ..
    else
        echo "[SuperLU_DIST] Skipping build step: $WORKING_DIR/superludist-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib64/libsuperlu_dist.a ]; then
        echo "[SuperLU_DIST] Install step..."
        cd superludist-build
        cmake --install . --prefix="$INSTALL_DIR" 2>&1 | tee "$WORKING_DIR"/superludist-install.log
        cd ..
    else
        echo "[SuperLU_DIST] Skipping install step: $INSTALL_DIR/lib64/libsuperlu_dist.a already exists"
    fi
fi

if hasTag "trilinos"; then
    if [[ ! -d "$WORKING_DIR"/scalapack-src ]]; then
        echo "[ScaLAPACK] Download step..."
        git clone https://github.com/Reference-ScaLAPACK/scalapack.git scalapack-src
        cd scalapack-src
        git checkout c856b245a5f4285500ddc39894132adbb98e6c3f
        cd ..
    else
        echo "[ScaLAPACK] Skipping download step: $WORKING_DIR/scalapack-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/scalapack-build ]; then
        echo "[ScaLAPACK] Build step..."
        mkdir -p scalapack-build
        cd scalapack-build

        # NOTE: The CMAKE_INSTALL_PREFIX is required here,
        # despite setting the prefix for the installation later as well
        # If this is not set here, scalapack will try to install its pkgconfig
        # files into /usr/local/lib/pkgconfig
        cmake ../scalapack-src \
            -DCMAKE_INSTALL_PREFIX="$INSTALL_DIR" \
            -DLAPACK_LIBRARIES=$INSTALL_DIR/lib/libopenblas.a \
            -DLAPACK_LIBRARIES=$INSTALL_DIR/lib/libopenblas.a \
            2>&1 | tee "$WORKING_DIR"/logs/scalapack-configure.log
        cmake --build . --target scalapack --config Release -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/scalapack-build.log

        cd ..
    else
        echo "[ScaLAPACK] Skipping build step: $WORKING_DIR/scalapack-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib64/libscalapack.a ]; then
        echo "[ScaLAPACK] Install step..."
        cd scalapack-build
        cmake --install . --prefix="$INSTALL_DIR" 2>&1 | tee "$WORKING_DIR"/logs/scalapack-install.log
        cd ..
    else
        echo "[ScaLAPACK] Skipping install step: $INSTALL_DIR/lib64/libscalapack.a already exists"
    fi
fi

if hasTag "zoltan" || hasTag "trilinos" || hasTag "mumps"; then
    if [ ! -d "$WORKING_DIR"/scotch-src ]; then
        echo "[scotch] Download step..."
        git clone https://gitlab.inria.fr/scotch/scotch.git scotch-src
    else
        echo "[scotch] Skipping download step: $WORKING_DIR/scotch-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/scotch-build ]; then
        echo "[scotch] Build step..."
        mkdir -p scotch-build
        cd scotch-build

        cmake ../scotch-src \
            -DCMAKE_BUILD_TYPE:STRING=Release \
            -DINSTALL_METIS_HEADERS:BOOL=OFF 2>&1 | tee "$WORKING_DIR"/logs/scotch-configure.log
        cmake --build . -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/scotch-build.log


        cd ..
    else
        echo "[scotch] Skipping build step: $WORKING_DIR/scotch-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib64/libscotch.a ]; then
        echo "[scotch] Install step..."
        cd scotch-build
        cmake --install . --prefix="$INSTALL_DIR" 2>&1 | tee "$WORKING_DIR"/logs/scotch-install.log
        cd ..
    else
        echo "[scotch] Skipping install step: $INSTALL_DIR/lib/libscotch.a already exists"
    fi
fi

if hasTag "mumps" || hasTag "trilinos"; then
    if [[ ! -d "$WORKING_DIR"/mumps-src ]]; then
        echo "[MUMPS] Download step..."
        wget https://mumps-solver.org/MUMPS_5.7.3.tar.gz
        mkdir mumps-src
        tar xf MUMPS_5.7.3.tar.gz --strip-components 1 -C mumps-src
        rm MUMPS_5.7.3.tar.gz
    else
        echo "[MUMPS] Skipping download step: $WORKING_DIR/mumps-src already exists"
    fi

    if [[ ! -d "$WORKING_DIR"/mumps-build ]]; then
        echo "[MUMPS] Build step..."
        cp -ra mumps-src mumps-build

        cd mumps-build

        sed -e "s,@INSTALL_DIR@,${INSTALL_DIR},g" -e "s,@OMP_FLAG@,${USE_OPENMP_FLAG},g" <<'EOF' > Makefile.inc
#
#  This file is part of MUMPS 5.7.3, released
#  on Mon Jul 15 11:44:21 UTC 2024
#
#Begin orderings

# NOTE that PORD is distributed within MUMPS by default. It is recommended to
# install other orderings. For that, you need to obtain the corresponding package
# and modify the variables below accordingly.
# For example, to have Metis available within MUMPS:
#          1/ download Metis and compile it
#          2/ uncomment (suppress # in first column) lines
#             starting with LMETISDIR,  LMETIS
#          3/ add -Dmetis in line ORDERINGSF
#             ORDERINGSF  = -Dpord -Dmetis
#          4/ Compile and install MUMPS
#             make clean; make   (to clean up previous installation)
#
#          Metis/ParMetis and SCOTCH/PT-SCOTCH (ver 6.0 and later) orderings are recommended.
#

SCOTCHDIR  = @INSTALL_DIR@
ISCOTCH    = -I$(SCOTCHDIR)/include
LSCOTCH    = -L$(SCOTCHDIR)/lib -lptesmumps -lptscotch -lptscotcherr -lscotch

LPORDDIR = $(topdir)/PORD/lib/
IPORD    = -I$(topdir)/PORD/include/
LPORD    = -L$(LPORDDIR) -lpord$(PLAT)

LMETISDIR = @INSTALL_DIR@/lib64
IMETIS    = -I@INSTALL_DIR@/include
LMETIS    = -L$(LMETISDIR) -lparmetis -lmetis -lgklib

ORDERINGSF = -Dscotch -Dmetis -Dpord -Dptscotch -Dparmetis
#ORDERINGSF  = -Dpord -Dparmetis -Dptscotch
#ORDERINGSF  = -Dpord -Dparmetis
ORDERINGSC  = $(ORDERINGSF)

LORDERINGS = $(LMETIS) $(LPORD) $(LSCOTCH)
IORDERINGSF = $(ISCOTCH)
IORDERINGSC = $(IMETIS) $(IPORD) $(ISCOTCH)

#End orderings
########################################################################
################################################################################

PLAT    =
LIBEXT  = .a
LIBEXT_SHARED  = .so
SONAME = -soname
SHARED_OPT = -shared
FPIC_OPT = -fPIC
# Adapt/uncomment RPATH_OPT to avoid modifying
# LD_LIBRARY_PATH in case of shared libraries
# RPATH_OPT = -Wl,-rpath,/path/to/MUMPS_x.y.z/lib/
OUTC    = -o
OUTF    = -o
RM = /bin/rm -f
CC = mpicc
FC = mpifort
FL = mpifort
AR = ar vr # WARNING: This comment is load-bearing. Removing it will cause MUMPS builds to fail, due to a missing trailing space behind vr
#RANLIB = ranlib
RANLIB  = echo
# Make this variable point to the path where the Intel MKL library is
# installed. It is set to the default install directory for Intel MKL.
#MKLROOT=/sfw/intel/2024.0.1/mkl/latest/lib/intel64
#LAPACK = -L$(MKLROOT) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
#SCALAP = -L$(MKLROOT) -lmkl_scalapack_lp64 -lmkl_blacs_intelmpi_lp64
LAPACK = @INSTALL_DIR@/lib/libopenblas.a
SCALAP  = -L@INSTALL_DIR@/lib -lscalapack

LIBPAR = $(SCALAP) $(LAPACK)

INCSEQ = -I$(topdir)/libseq
LIBSEQ  = $(LAPACK) -L$(topdir)/libseq -lmpiseq$(PLAT)

#LIBBLAS = -L$(MKLROOT) -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core
LIBBLAS = @INSTALL_DIR@/lib/libopenblas.a
LIBOTHERS = -lpthread

#Preprocessor defs for calling Fortran from C (-DAdd_ or -DAdd__ or -DUPPER)
CDEFS   = -DAdd_

#Begin Optimized options
OPTF    = -O @OMP_FLAG@ -fallow-argument-mismatch -DGEMMT_AVAILABLE
OPTL    = -O @OMP_FLAG@
OPTC    = -O @OMP_FLAG@
#End Optimized options

INCS = $(INCPAR)
LIBS = $(LIBPAR)
LIBSEQNEEDED =

EOF

        make d -j "${NP}" 2>&1 | tee "$WORKING_DIR"/logs/mumps-build.log
        cd ..
    else
        echo "[MUMPS] Skipping build step: $WORKING_DIR/mumps-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib/libmumps_common.a ]; then
        echo "[MUMPS] Install step..."
        cd mumps-build
        cp -r lib "$INSTALL_DIR"
        cp -r include "$INSTALL_DIR"
        cd ..
    else
        echo "[MUMPS] Skipping install step: $INSTALL_DIR/lib/libmumps_common.a already exists"
    fi
fi

if hasTag "zoltan"; then
    if [ ! -d "$WORKING_DIR"/zoltan-src ]; then
        echo "[Zoltan] Download step..."
        git clone --branch v3.901 --depth=1 https://github.com/sandialabs/Zoltan.git zoltan-src
    else
        echo "[Zoltan] Skipping download step: $WORKING_DIR/zoltan-src already exists"
    fi

    if [ ! -d "$WORKING_DIR"/zoltan-build ]; then
        echo "[Zoltan] Build step..."
        mkdir -p zoltan-build

        cd zoltan-build
        CC=mpicc ../zoltan-src/configure \
            --prefix="$INSTALL_DIR" \
            --with-gnumake \
            --with-mpi="$MPI_HOME" \
            --with-scotch \
            --with-scotch-incdir="$INSTALL_DIR"/include/ \
            --with-scotch-libdir="$INSTALL_DIR"/lib64 \
            --with-parmetis \
            --with-parmetis-incdir="$INSTALL_DIR"/include \
            --with-parmetis-libdir="$INSTALL_DIR"/lib64 \
            --with-libs="-lgklib" \
            2>&1 | tee "$WORKING_DIR"/logs/zoltan-configure.log

        make everything 2>&1 | tee "$WORKING_DIR"/logs/zoltan-build.log
        cd ..
    else
        echo "[Zoltan] Skipping build step: $WORKING_DIR/superludist-build already exists"
    fi

    if [ ! -f "$INSTALL_DIR"/lib/libzoltan.a ]; then
        echo "[Zoltan] Install step..."
        cd zoltan-build
        make install 2>&1 | tee "$WORKING_DIR"/logs/zoltan-install.log
        cd ..
    else
        echo "[Zoltan] Skipping install step: $INSTALL_DIR/lib/libzoltan.a already exists"
    fi
fi

if hasTag "trilinos" && [ -d "$TRILINOS_DIR" ]; then
    if [ ! -d "$WORKING_DIR"/trilinos-build ]; then
        echo "[Trilinos] Build step..."
        mkdir -p trilinos-build
        cd trilinos-build
        sed -i 's/REQUIRED_LIBS_NAMES umfpack amd/REQUIRED_LIBS_NAMES umfpack amd cholmod camd ccolamd colamd suitesparseconfig/' "$TRILINOS_DIR"/cmake/TPLs/FindTPLUMFPACK.cmake

        export CC=mpicc
        export CXX=mpicxx
        export FC=mpif90

        # build tests(should be activated for testing purposes.... duhh)

        HDF5_INCLUDE="${INSTALL_DIR}"/include
        HDF5_LIB="${INSTALL_DIR}"/lib
        UMFPACK_ROOT="${INSTALL_DIR}"
        MUMPS_ROOT="${INSTALL_DIR}"

        cmake "$TRILINOS_DIR"  \
            -D Trilinos_ASSERT_MISSING_PACKAGES=OFF \
            -D CMAKE_VERBOSE_MAKEFILE:BOOL=OFF \
            -D "CMAKE_C_FLAGS:STRING=-DH5_HAVE_PARALLEL -DSHYLU_NODEBASKER -DFROSCH_TIMER_DETAILS=2" \
            -D "CMAKE_CXX_FLAGS:STRING=-O3 -fPIC -fno-var-tracking -DH5_HAVE_PARALLEL -DSHYLU_NODEBASKER -DFROSCH_TIMER_DETAILS=2" \
            -D CMAKE_CXX_STANDARD:STRING=20 \
            -D MPI_EXEC_MAX_NUMPROCS:STRING=4 \
            -D Trilinos_ENABLE_Fortran:BOOL=ON \
            -D Trilinos_ENABLE_OpenMP:BOOL=${USE_OPENMP} \
            -D Trilinos_ENABLE_Teuchos:BOOL=ON \
            -D Trilinos_ENABLE_Epetra:BOOL=OFF \
            -D Trilinos_ENABLE_EpetraExt:BOOL=OFF \
            -D Trilinos_ENABLE_AztecOO:BOOL=OFF \
            -D Trilinos_ENABLE_Belos:BOOL=ON \
            -D Trilinos_ENABLE_Anasazi:BOOL=ON \
            -D Trilinos_ENABLE_Amesos:BOOL=OFF \
            -D Trilinos_ENABLE_Amesos2:BOOL=ON \
            -D Trilinos_ENABLE_NOX:BOOL=ON \
            -D Trilinos_ENABLE_Zoltan:BOOL=ON \
            -D Trilinos_ENABLE_Zoltan2:BOOL=ON \
            -D Trilinos_ENABLE_Ifpack2:BOOL=ON \
            -D Trilinos_ENABLE_ML:BOOL=OFF \
            -D Trilinos_ENABLE_Thyra:BOOL=ON \
            -D Trilinos_ENABLE_ShyLU_DDFROSch:BOOL=ON \
            -D ShyLU_DDFROSch_ENABLE_Xpetra:BOOL=ON \
            -D ShyLU_DD_ENABLE_TESTS:STRING="${BUILD_TESTS}" \
            -D Trilinos_ENABLE_Stratimikos:BOOL=ON \
            -D Trilinos_ENABLE_MueLu:BOOL=OFF \
            -D Trilinos_ENABLE_Tpetra:BOOL=ON \
            -D Trilinos_ENABLE_Xpetra:BOOL=ON \
            -D Tpetra_INST_INT_INT:BOOL=ON \
            -D Tpetra_INST_INT_LONG_LONG:BOOL=OFF\
            -D TPL_FIND_SHARED_LIBS:BOOL=OFF \
            -D TPL_ENABLE_MPI:BOOL=ON \
            -D TPL_ENABLE_Boost:BOOL=ON \
            -D TPL_ENABLE_LAPACK:BOOL=ON \
            -D TPL_ENABLE_Matio:BOOL=OFF \
            -D TPL_ENABLE_METIS:BOOL=ON \
            -D METIS_LIBRARY_DIRS:PATH="$INSTALL_DIR/lib64" \
            -D METIS_INCLUDE_DIRS:PATH="$INSTALL_DIR/include" \
            -D METIS_LIBRARY_NAMES:STRING="metis;gklib" \
            -D TPL_ENABLE_BLAS:BOOL=ON \
            -D TPL_ENABLE_ParMETIS:BOOL=ON \
            -D ParMETIS_LIBRARY_DIRS:PATH="$INSTALL_DIR/lib64"  \
            -D ParMETIS_INCLUDE_DIRS:PATH="$INSTALL_DIR/include" \
            -D ParMETIS_LIBRARY_NAMES:STRING="parmetis" \
            -D TPL_ENABLE_Scotch:BOOL=ON \
            -D Scotch_LIBRARY_DIRS:PATH="$INSTALL_DIR/lib64"  \
            -D Scotch_INCLUDE_DIRS:PATH="$INSTALL_DIR/include" \
            -D Scotch_LIBRARY_NAMES:STRING="ptscotch;ptscotcherr;scotch;scotcherr" \
            -D Amesos_ENABLE_MUMPS:BOOL=OFF \
            -D Amesos2_ENABLE_MUMPS:BOOL=ON \
            -D Amesos2_ENABLE_SuperLUDist:BOOL=OFF \
            -D Amesos_ENABLE_SuperLUDist:BOOL=OFF \
            -D Amesos2_ENABLE_UMFPACK:BOOL=ON \
            -D Amesos2_ENABLE_PARDISO_MKL:BOOL=OFF \
            -D TPL_ENABLE_SuperLUDist=OFF \
            -D TPL_ENABLE_PARDISO_MKL=OFF \
            -D TPL_ENABLE_HDF5:BOOL=ON \
            -D Trilinos_ENABLE_TESTS:BOOL="${BUILD_TESTS}" \
            -D NOX_ENABLE_TESTS:BOOL="${BUILD_TESTS}" \
            -D BUILD_SHARED_LIBS=OFF                           \
            -D CMAKE_BUILD_TYPE=Release \
            -D HDF5_LIBRARY_DIRS:PATH="$INSTALL_DIR/lib" \
            -D HDF5_INCLUDE_DIRS:PATH="$INSTALL_DIR/include" \
            -D TPL_BLAS_LIBRARIES="$INSTALL_DIR"/lib/libopenblas.a \
            -D TPL_LAPACK_LIBRARIES="$INSTALL_DIR"/lib/libopenblas.a \
            -D BLAS_INCLUDE_DIRS:PATH="$INSTALL_DIR/include" \
            -D LAPACK_INCLUDE_DIRS:PATH="$INSTALL_DIR/include" \
            -D TPL_ENABLE_SCALAPACK:BOOL=ON \
            -D SCALAPACK_LIBRARY_DIRS:STRING="${INSTALL_DIR}/lib" \
            -D SCALAPACK_LIBRARY_NAMES:STRING="libscalapack.a" \
            -D UMFPACK_LIBRARIES="${INSTALL_DIR}/lib64/libamd.a;${INSTALL_DIR}/lib64/libcamd.a;${INSTALL_DIR}/lib64/libccolamd.a;${INSTALL_DIR}/lib64/libcholdmod.a;${INSTALL_DIR}/lib64/libcolamd.a;${INSTALL_DIR}/lib64/libsuitesparseconfig.a;${INSTALL_DIR}/lib64/libumfpack.a;$INSTALL_DIR/lib/libopenblas.a" \
            -D UMFPACK_INCLUDE_DIRS:STRING="${INSTALL_DIR}/include/suitesparse" \
            -D MUMPS_LIBRARY_DIRS:STRING="${INSTALL_DIR}/lib" \
            -D MUMPS_INCLUDE_DIRS:STRING="${INSTALL_DIR}/include" \
            -D MUMPS_LIBRARY_NAMES:STRING="dmumps;mumps_common;pord;parmetis;metis;esmumps;ptscotch;ptscotcherr;scotch;scotcherr;scalapack;openblas" \
            -D TRIBITS_TPL_FIND_INCLUDE_DIRS_AND_LIBRARIES_VERBOSE:BOOL=ON \
            -D CMAKE_PREFIX_PATH:PATH="$INSTALL_DIR" \
            -D CMAKE_INSTALL_PREFIX:PATH="$INSTALL_DIR" 2>&1 |
            tee "$WORKING_DIR"/logs/trilinos-configure.log

        make all -j4 2>&1 | tee "$WORKING_DIR"/logs/trilinos-build.log || exit
        make install 2>&1 | tee "$WORKING_DIR"/logs/trilinos-install.log || exit

            #-D Trilinos_EXTRA_LINK_FLAGS:STRING="-lgfortran;-lomp;-lmpi_mpifh;-L${INSTALL_DIR}/lib;-L${INSTALL_DIR}/lib64;-lscalapack;-lopenblas;-lpord;-lGKlib" \


        cd ..
    else
        echo "[Trilinos] Skipping build step: $WORKING_DIR/trilinos-build already exists"
    fi
fi
