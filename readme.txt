FEAT 3 README

=== DEPENDECIES ===
FEAT depends on:
cmake 2.8 or greater
python 2.6 or greater
gcc 4.8+ or icc 15+ or clang 3.6+ or MSVC++ 2015+

Additionally you will need a proper buildsystem, e.g. gnumake.

=== BUILDING FEAT ===
In its simplest shape, building FEAT is just as easy as
./configure
make all
This will automagically detect a proper tool/setting setup and build all FEAT components.

'./configure help' reveals a lot more configuration possiblities.

=== TEST RUN ===
To check if your system created a correct FEAT binary, you can run a set of unittests with the additional command
ctest

=== DOCUMENATION ===
A comprehensive doxygen documentation is available in the doc/html folder, after issueing the
make doc
command ( do not forget the ./configure command, if not already configured).

The tutorial folder contains some sample applications with extensive in code documentation.

=== ADDITIONAL SOFTWARE SUPPORT ===
OpenMPI 1.8.5+
MPICH 3.2+
Microsoft MPI 7.1
CUDA 6.5+
CCache 3.2+
Ninja 1.6+
Valgrind 3.10+
Score-P 1.4.2+
Intel MKL 15+

Alglib, Fparser, Umfpack and Parmetis are downloaded and build automatically, when included in the configure flags.

=== git commit hooks ===
When working with the git repository, the following instructions will install proper pre commit hooks.
In the FEAT source root directory, execute:
 cd .git/hooks
 git init
 git pull .. remotes/origin/hooks
