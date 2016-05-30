#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'exp2dv exp-file dv-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    String input(argv[1]);
    String output(argv[2]);

    DenseVector<Mem::Main, double> dv(FileMode::fm_exp, input);
    dv.write_out(FileMode::fm_dv, output);
}
