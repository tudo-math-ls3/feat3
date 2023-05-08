%% FEAT3: Finite Element Analysis Toolbox, Version 3
%% Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
%% FEAT3 is released under the GNU General Public License version 3,
%% see the file 'copyright.txt' in the top level directory for details.
function mat = read_feat_binary_matrix(filename)
fileID = fopen(filename);
fseek(fileID, 8,"bof");
type_id  = fread(fileID, 1, '*uint64');
if (type_id == 4)
    mat = read_feat_binary_sparse_matrix(fileID);
elseif (type_id == 7)
    mat = read_feat_binary_dense_matrix(fileID);
else
    error("Incompatible type ID in input file.")
end
fclose(fileID);
end

function mat = read_feat_binary_sparse_matrix(fid)
frewind(fid);
feat_header   = fread(fid, [21, 1], '*uint64');
m = feat_header(19);
n = feat_header(20);
nel = feat_header(21);
vals = fread(fid, nel, '*double');
col_ind = fread(fid, nel, '*uint64');
row_ptr = fread(fid, m+1, '*uint64') + 1;
row_ind = ones(nel, 1, 'uint64');
for i=1:m
    row_ind( row_ptr(i):(row_ptr(i+1)-1) ) = i;
end
mat = sparse(row_ind, col_ind+1, vals, m, n);
end

function mat = read_feat_binary_dense_matrix(fid)
frewind(fid);
feat_header   = fread(fid, [16, 1], '*uint64');
m = feat_header(15);
n = feat_header(16);
mat = fread(fid, [n, m], '*double');
mat = transpose(mat);
end