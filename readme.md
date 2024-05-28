Finite Element Analysis Toolbox 3
=================================

Welcome to the Finite Element Analysis Toolbox 3 (FEAT3) library, which is developed and maintained
by the FEAT group under Stefan Turek from the chair of Applied Mathematics and Numerics (LS3) at the
Department of Mathematics of the Technische Universitaet Dortmund, which is located in Germany.

FEAT3 is an open-source finite element software library developed primarily in standard-conforming
C++ and it is a part of the FEATLFOW software family and it is the common successor to both the
FEATLFOW2 and FEAST software projects, which are both descendants of the original FEAT 1 libraries:
http://www.featflow.de

The primary application field of FEAT3 are the incompressible Navier-Stokes equations in two or
three dimensions on complex domains with focus on providing fast and robust solution methods for
the simulation of highly viscous non-Newtonian fluids. FEAT3 aims to combine the flexibility and
robustness of both conforming and non-conforming Finite Element Methods with fast and robust
iterative solvers based on geometric multigrid methods by providing efficient implementations of
the discretization and solution steps, which utilize both the traditional MPI-based "scale-out"
parallelization as well as offloading of crucial code paths to modern accelerator hardware via CUDA.

License
-------
FEAT3 is released under the GNU General Public License version 3;
please refer to the 'copyright.txt' in the top level directory for details.

GitHub Repository
-----------------
The source code of FEAT3 is maintained on an internal git server, but we keep a publicly accessible
mirror of our master branch at a GitHub repository, which can be found here:
https://github.com/tudo-math-ls3/feat3

Doxygen Code Documentation
--------------------------
The source code of FEAT3 is documented via [doxygen](https://www.doxygen.nl) and the documentation
for the current version of the master branch is publicly available on the GitHub page under
https://tudo-math-ls3.github.io/feat3

The code documentation can be created in your local git clone checkout via the command `make doc`
once the build system has been configured (see below).

Dependencies
------------
FEAT3 is written in standard-conforming C\+\+17 and therefore it can (in theory) be compiled by any
modern C++ compiler independent of the underlying operating system and hardware without the need to
install any additional third-party libraries, since all other dependencies (including MPI and CUDA)
are treated as opt-in, however, the functionality of FEAT3 will obviously be severely limited when
compiled without any third-party libraries.

FEAT3 is regularly being built and tested with various compilers and is confirmed to compile and
work properly with at least the following C++ compilers:
* GNU g++ compiler, version 7.0 or newer
* Intel C/C++ compiler, version 19.0 or newer
* Intel oneAPI DPC++/C++ compiler, version 2021 or newer
* Clang LLVM compiler, version 8.0 or newer
* Microsoft Visual C++ 2019 compiler or newer

In addition to a modern C++ compiler, you also need the following tools to set up the build system
of FEAT3:
* CMake version 3.18.0 or newer
* Python version 3.7 or newer
* gnumake or ninja or some other build system supported by CMake

All other dependencies (e.g. MPI, CUDA, MKL, etc.) are opt-in and have to be enabled by supplying
their corresponding _build-id tag_ to our configure script to be used by FEAT3, see below.

Configure and Build
-------------------
The build system of FEAT3 is designed for _out-of-source_ builds, i.e. one has to create a new
build directory outside of the source directory and call the configure script from within that
build directory, which will then call CMake to generate the build system files.

Assuming your local git clone of the FEAT3 repository is located under `~/feat3.git`, you can simply
run the following commands to create a build directory, configure and build FEAT3 in its simplest
form:

    mkdir ~/feat3.build
    cd ~/feat3.build
    ~/feat3.git/configure
    make

To enable more sophisticated features, you have to supply a _build-id_ and optionally some
additional flags to the configure scripts. In its simplest form, the build-id contains the name
of the compiler that you want to use as well as either the **debug** or **opt** tag, depending on
whether you want to compile in debug mode or in optimized mode. In addition to that, the build-id
may also contain the tags for the additional features that you want to enable, e.g. **mpi** or
**cuda**:

    ~/feat3.git/configure gcc-opt-mpi-cuda

To obtain a list of all supported build-id tags, including the tags for the compiler that you want
to use, as well as all other additional options supported by the configure script, simply call

    ~/feat3.git/configure help

Git-Commit Hooks
----------------
This section is only important for FEAT3 developers: please execute the following four commands
after you have cloned the repository to install our commit hooks **before** you make any changes to
the code and **before** you make any commits, because our internal git repository will reject your
commits if the formatting of your code is incorrect, which can be a pain to fix afterwards.

Again, assuming that your git clone of the FEAT3 repository is found under `~/feat3.git`, execute
the following commands to set up the git hooks in your local repository:

    cd ~/feat3.git/.git/hooks
    git init
    git pull .. remotes/origin/hooks
    cd ../..
