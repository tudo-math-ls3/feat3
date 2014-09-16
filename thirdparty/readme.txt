This file keeps track of the third-party library versions that have been
tested with FEAST.

ALGLIB                  3.8.2
AMD                     2.3.1, 2.4.0
SuiteSparse_config      4.2.1, 4.3.1
UMFPACK                 5.6.2, 5.7.0

As of 16.09.2014, recent versions of AMD, SuiteSparse_config and UMFPACK
are provided at their usual URLs. If you want to always use the current
versions of these packages, simply remove the version number from the
filename variables in cmake_modules/umfpack.py. You might end up with an
untested and incompatible version, though.
