%% FEAT3: Finite Element Analysis Toolbox, Version 3
%% Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
%% FEAT3 is released under the GNU General Public License version 3,
%% see the file 'copyright.txt' in the top level directory for details.
function write_feat_binary_vector(u, filename)
% WRITE_FEAT_BINARY_VECTOR export vecor to FEAT3-compatible binary file
% input arguments:
% u the vector that should be exported. if it is a blockvector it has to be
% in the format u = [u_1 u_2 ... u_n] and will be saved as
% DenseVectorBlocker
  vector    = u';
  lv        = numel(u);
  block_dim = size(u, 2);
  if (block_dim == 1)
    idx = 1;
  else
    idx = 10;
  end
  gsize     = 8*(lv+14)+16;
  header    = uint64([gsize,        ... %% outputstream size
                      idx,          ... %% magic number -> BlockedVector
                      25769803784,  ... %% double "feature hash"
                      4294967304,   ... %% uint64 "feature hash"
                      1,            ... %% _elements size
                      0,            ... %% _indices size
                      1,            ... %% _elements_size size
                      0,            ... %% _indices_size size
                      1,            ... %% _scalar_index size
                      0,            ... %% _scalar_dt size
                      17,           ... %% compression magic number -> no compression
                      lv,           ... %% length of _elements[0]
                      lv*8,         ... %% size (bytes) of _elements[0]
                      lv/block_dim]);   %% _scalar_index[0] -> native size
  fileID = fopen(filename, 'w');
  fwrite(fileID, header, 'uint64');
  fwrite(fileID, vector, 'double');
  fwrite(fileID, [0 0], 'uint64');
  fclose(fileID);
end