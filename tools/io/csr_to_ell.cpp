#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <iostream>

using namespace FEAST;
using namespace FEAST::LAFEM;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'csr2ell csr-file ell-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    String input(argv[1]);
    String output(argv[2]);

    SparseMatrixCSR<Mem::Main, double> csr(FileMode::fm_csr, input);
    SparseMatrixELL<Mem::Main, double> ell(csr);
    ell.write_out(FileMode::fm_ell, output);
}
