FEAST README

=== DEPENDECIES ===
FEAST depends on:
cmake 2.8 or greater
python 2.6 or greater
gcc 4.9+ or icc 14+ or clang 3.6+

Additionally you will need a proper c++ compiler and buildsystem.
In most cases, this might be g++ and gnumake.

=== BUILDING FEAST ===
In its simplest shape, building FEAST is just as easy as
./configure-feast
make all
This will automagically detect a proper tool/setting setup and build all FEAST components.

'./configure-feast help' reveals a lot more configuration possiblities.

=== TEST RUN ===
To check if your system created a correct FEAST binary, you can run a set of unittests with the additional command
ctest
