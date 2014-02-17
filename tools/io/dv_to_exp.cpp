#include <kernel/lafem/dense_vector.hpp>
#include <iostream>

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

    DenseVector<Mem::Main, double> dv(FileMode::fm_dv, input);
    dv.write_out(FileMode::fm_exp, output);
}
