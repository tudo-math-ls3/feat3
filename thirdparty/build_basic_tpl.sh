tpl_dir = "${TPL_DIR}"
install_dir = "${INSTALL_DIR}"
export PATH=${install_dir}:$PATH

cd $tpl_dir

module purge
module load gcc/13.2.0 cmake/3.28.3 openmpi/4.1.6 binutils/2.43 openblas/0.3.29

export CC=gcc
export CXX=g++


## install umfpack
if [[ ! -d SuiteSparse ]]; then
        git clone --depth 1 --branch v7.10.2 https://github.com/DrTimothyAldenDavis/SuiteSparse.git
fi

cd SuiteSparse && mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX="${install_dir}" -DSUITESPARSE_ENABLE_PROJECTS="umfpack" -DCMAKE_BUILD_TYPE="Release" \
      -DSUITESPARSE_USE_OPENMP="ON" -DBUILD_SHARED_LIBS="OFF" -DSUITESPARSE_USE_64BIT_BLAS="ON" -DSUITESPARSE_USE_FORTRAN="OFF" ..
make -j4
make install

export CC=mpicc
export CXX=mpicxx

## install hypre
if [[ ! -d hypre ]]; then
        git clone --depth 1 -branch v2.25.0 https://github.com/hypre-space/hypre
fi

cd hypre && mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX="${install_dir}" -DCMAKE_BUILD_TYPE="Release" \
      -DHYPRE_WITH_OPENMP="ON" -DHYPRE_ENABLE_CUDA_STREAM="OFF" -DHYPRE_ENABLE_CURAND="OFF" -DHYPRE_ENABLE_CUSPARSE="OFF" ../src/
make -j4
make install

## install boost

if [[ ! -d boost ]]; then
        git clone --depth 1 --branch boost-1.88.0 --recursive https://github.com/boostorg/boost.git
fi

export CC=gcc
export CXX=g++

cd boost && mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX="${install_dir}" -DCMAKE_BUILD_TYPE="Release" \
      -DBUILD_SHARED_LIBS="OFF" ..
make -j4
make install


## install cgal
if [[ ! -d cgal ]]; then
        git clone --depth 1 --branch v6.0.1 --recursive https://github.com/CGAL/cgal.git
fi

export CC=gcc
export CXX=g++

cd cgal && mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX="${install_dir}" -DCMAKE_BUILD_TYPE="Release" ..
make -j4
make install

## install parmetis
if [[ ! -d ParMETIS ]]; then
        git clone --depth 1 --recursive https://github.com/otkmesser/ParMETIS
fi

export CC=mpicc
export CXX=mpicxx

cd ParMETIS && mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX="${install_dir}" -DCMAKE_BUILD_TYPE="Release" ..
make -j4
make install

## install zlib
if [[ ! -d zlib ]]; then
        git clone --depth 1 --branch v1.3.1 https://github.com/madler/zlib.git
fi

export CC=gcc
export CXX=g++

cd zlib && mkdir -p build && cd build
cmake -DCMAKE_INSTALL_PREFIX="${install_dir}" -DCMAKE_BUILD_TYPE="Release"
      -DCMAKE_BUILD_TESTING=OFF -DZLIB_BUILD_SHARED=OFF -DZLIB_INSTALL_COMPAT_DLL=OFF -DZLIB_BUILD_MINIZIP=OFF ..
make -j4
make install
