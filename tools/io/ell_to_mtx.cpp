#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <iostream>

using namespace FEAST;
using namespace FEAST::LAFEM;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'ell2mtx ell-file mtx-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    String input(argv[1]);
    String output(argv[2]);

    SparseMatrixELL<Mem::Main, double> ell(input);
    ell.write_out(fm_mtx, output);
}
