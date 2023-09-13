%% FEAT3: Finite Element Analysis Toolbox, Version 3
%% Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
%% FEAT3 is released under the GNU General Public License version 3,
%% see the file 'copyright.txt' in the top level directory for details.
function write_feat_binary_matrix(mat, filename)
if (issparse(mat))
    write_feat_binary_sparse_matrix(mat, filename);
else
    write_feat_binary_dense_matrix(mat, filename);
end
end

function write_feat_binary_sparse_matrix(mat, filename)
lm        = nnz(mat);
[m, n]    = size(mat);
[col_ind, row_ind, vals] = find(transpose(mat));
row_ptr = cumsum(histcounts(row_ind, 0:(m+1)));
col_ind = col_ind - 1; % 0-based indexing

gsize     = 8 *(21 + 2*lm + m + 1) + 16;
header    = uint64([gsize,        ... %% outputstream size
                    4,            ... %% magic number -> SparseMtrixCSR /fm_csr
                    25769803784,  ... %% double "feature hash"
                    4294967304,   ... %% uint64 "feature hash"
                    1,            ... %% _elements size
                    2,            ... %% _indices size
                    1,            ... %% _elements_size size
                    2,            ... %% _indices_size size
                    4,            ... %% _scalar_index size
                    0,            ... %% _scalar_dt size
                    17,           ... %% compression magic number -> no compression
                    lm,           ... %% length of _elements[0]
                    lm*8,         ... %% size (bytes) of _elements[0]
                    lm,           ... %% length of _indices[0]
                    (m+1),        ... %% length of _indices[1]
                    lm*8,         ... %% size (bytes) of _indices[0]
                    (m+1)*8,      ... %% size (bytes) of _indices[1]
                    m*n,          ... %% _scalar_index[0] -> size
                    m,            ... %% _scalar_index[1] -> rows
                    n,            ... %% _scalar_index[2] -> columns
                    lm]);         ... %% _scalar_index[3] -> used elements

fileID = fopen(filename, 'w');
fwrite(fileID, header, 'uint64');
fwrite(fileID, vals,   'double');
fwrite(fileID, col_ind,'uint64');
fwrite(fileID, row_ptr,'uint64');
fwrite(fileID, [0 0],  'uint64');
fclose(fileID);
end

function write_feat_binary_dense_matrix(mat, filename)
lv        = numel(mat);
[m, n]    = size(mat);
gsize     = 8*(lv+16)+16;
header    = uint64([gsize,        ... %% outputstream size
                    7,            ... %% magic number -> DenseMatrix /fm_dm
                    25769803784,  ... %% double "feature hash"
                    4294967304,   ... %% uint64 "feature hash"
                    1,            ... %% _elements size
                    0,            ... %% _indices size
                    1,            ... %% _elements_size size
                    0,            ... %% _indices_size size
                    3,            ... %% _scalar_index size
                    0,            ... %% _scalar_dt size
                    17,           ... %% compression magic number -> no compression
                    lv,           ... %% length of _elements[0]
                    lv*8,         ... %% size (bytes) of _elements[0]
                    lv,           ... %% _scalar_index[0] -> native size
                    m,            ... %% _scalar_index[1] -> rows
                    n]);              %% _scalar_index[2] -> columns

fileID = fopen(filename, 'w');
fwrite(fileID, header, 'uint64');
fwrite(fileID, transpose(mat), 'double');
fwrite(fileID, [0 0], 'uint64');
fclose(fileID);
end