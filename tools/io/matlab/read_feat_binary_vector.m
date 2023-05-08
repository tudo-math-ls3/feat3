%% FEAT3: Finite Element Analysis Toolbox, Version 3
%% Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
%% FEAT3 is released under the GNU General Public License version 3,
%% see the file 'copyright.txt' in the top level directory for details.
function u = read_feat_binary_vector(filename)
  fileID = fopen(filename);
  feat_header   = fread(fileID, [14, 1], '*uint64');
  size_total = feat_header(12);
  size_native = feat_header(14);
  bd = size_total / size_native;
  u   = fread(fileID, [bd, size_native], 'double');
  fclose(fileID);
  u = u';
end