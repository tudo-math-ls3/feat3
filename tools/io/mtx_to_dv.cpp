// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;

int main(int argc, char ** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  if (argc != 3)
  {
    std::cout<<"Usage 'mtx2dv mtx-file dv-file'"<<"\n";
    FEAT::Runtime::abort();
  }

  String input(argv[1]);
  String output(argv[2]);

  DenseVector<double, Index> dv(FileMode::fm_mtx, input);
  dv.write_out(FileMode::fm_dv, output);
  return 0;
}
