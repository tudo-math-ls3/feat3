This file keeps track of the third-party library versions that have been
tested with FEAST.

ALGLIB                  3.8.2, 3.9.0
ParMETIS                4.0.3
SuiteSparse             4.4.3



Notes: ParMETIS for Visual Studio
---------------------------------
If you want to use ParMETIS with Visual Studio, you will need to perform
two manual steps after unpacking the downloaded parmetis package:

1. Rename the directory from 'parmetis-4.0.3' to 'parmetis'

2. Run the 'parmetis_patch_win.py' script, which will patch one header
   file of the GKlib used by ParMETIS. Without this patch, you will
   not be able to compile ParMETIS.
