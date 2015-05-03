#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/adjacency/export_tga.hpp>
#include <iostream>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::Adjacency;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'csr2tga csr-file tga-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    String input(argv[1]);
    String output(argv[2]);

    SparseMatrixCSR<Mem::Main, double> matrix(FileMode::fm_csr, input);
    ExportTGA::write(output, matrix);
}
