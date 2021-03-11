// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/adjacency/export_tga.hpp>
#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Adjacency;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'ell2tga ell-file tga-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    String input(argv[1]);
    String output(argv[2]);

    SparseMatrixELL<Mem::Main, double> matrix(FileMode::fm_ell, input);
    ExportTGA::write(output, matrix);
}
