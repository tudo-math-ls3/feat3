// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;

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
