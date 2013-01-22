#include <iostream>
#include <kernel/lafem/dense_vector.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'dv2exp dv-file exp-file'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    String input(argv[1]);
    String output(argv[2]);

    DenseVector<Mem::Main, double> dv(fm_dv, input);
    dv.write_out(fm_exp, output);
}
