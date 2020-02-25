FEAT 3 README

=== LICENSE ===
FEAT3 is released under the GNU General Public License version 3,
see the file 'copyright.txt' in the top level directory for details.


=== DEPENDECIES ===
FEAT depends on:
cmake 3.9 or greater
python 2.6 or greater
gcc 4.8+ or icc 15+ or clang 3.6+ or MSVC++ 2015+ or AOCC 1.3+

Additionally you will need a proper buildsystem, e.g. gnumake.

=== BUILDING FEAT ===
The usual workflow to build FEAT starts with a second independent directory, e.g. called feat-build.
This is the folder where all binaries and object files during the build process will be stored.
You need to create it on your own (and not in your FEAT source directory).
From the newly created feat-build folder, call
# PATH_TO_FEAT_SOURCE/configure
to trigger the configuration process with default parameters
Afterwards
make all

In its simplest shape, building FEAT is just as easy as
./configure
# make all
will build the FEAT kernel and some sample applications and tutorials in the ./tutorial directory.

'./configure help' reveals a lot more configuration possiblities.

=== TEST RUN ===
To check if your system created a correct FEAT binary, you can run a set of unittests with the additional command
# ctest

=== DOCUMENATION ===
A comprehensive doxygen documentation is available in the doc/html folder, after issueing the
# make doc
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
Intel MKL 16+ / 11.3.4+

Alglib, Fparser, Umfpack and Parmetis are downloaded and build automatically, when included in the configure flags.

=== git commit hooks ===
When working with the git repository, the following instructions will install proper pre commit hooks.
In the FEAT source root directory, execute:
 cd .git/hooks
 git init
 git pull .. remotes/origin/hooks
