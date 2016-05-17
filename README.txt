FEAST README

=== DEPENDECIES ===
FEAST depends on:
cmake 2.8 or greater
python 2.6 or greater
gcc 4.8+ or icc 15+ or clang 3.6+

Additionally you will need a proper buildsystem, e.g. gnumake.

=== BUILDING FEAST ===
In its simplest shape, building FEAST is just as easy as
./configure-feast
make all
This will automagically detect a proper tool/setting setup and build all FEAST components.

'./configure-feast help' reveals a lot more configuration possiblities.

=== TEST RUN ===
To check if your system created a correct FEAST binary, you can run a set of unittests with the additional command
ctest

=== ADDITIONAL SOFTWARE SUPPORT ===
OpenMPI 1.8.5+
MPICH 3.2+
CUDA 6.5+
CCache 3.2+
Ninja 1.6+
Valgrind 3.10+
Score-P 1.4.2+

Alglib, Fparser, Umfpack and Parmetis are downloaded and build automatically.
