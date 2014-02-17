#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <iostream>

using namespace FEAST;
using namespace FEAST::LAFEM;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'ell2csr ell-file csr-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    String input(argv[1]);
    String output(argv[2]);

    SparseMatrixELL<Mem::Main, double> ell(FileMode::fm_ell, input);
    SparseMatrixCSR<Mem::Main, double> csr(ell);
    csr.write_out(FileMode::fm_csr, output);
}
