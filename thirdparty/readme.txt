This file keeps track of the third-party library versions that have been
tested with FEAT.

ALGLIB                  3.10.0, 3.13.0
fparser                 4.5.2
half                    1.12.0
ParMETIS                4.0.3
SuiteSparse             4.4.3, 5.2.0
Triangle                1.6
zlib                    1.2.11



Notes: ParMETIS for Visual Studio
---------------------------------
If you want to use ParMETIS with Visual Studio, you will need to perform
two manual steps after unpacking the downloaded parmetis package:

1. Rename the directory from 'parmetis-4.0.3' to 'parmetis'

2. Run the 'parmetis_patch_win.py' script, which will patch one header
   file of the GKlib used by ParMETIS. Without this patch, you will
   not be able to compile ParMETIS.
