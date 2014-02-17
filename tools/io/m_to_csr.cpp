#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <iostream>

using namespace FEAST;
using namespace FEAST::LAFEM;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'm2csr m-file csr-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    String input(argv[1]);
    String output(argv[2]);

    SparseMatrixCSR<Mem::Main, double> csr(FileMode::fm_m, input);
    csr.write_out(FileMode::fm_csr, output);
}
