#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <iostream>

using namespace FEAST;
using namespace FEAST::LAFEM;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'm2ell m-file ell-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    String input(argv[1]);
    String output(argv[2]);

    SparseMatrixELL<Mem::Main, double> ell(FileMode::fm_m, input);
    ell.write_out(FileMode::fm_ell, output);
}
