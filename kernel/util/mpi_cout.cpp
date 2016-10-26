#include<kernel/util/mpi_cout.hpp>

#include <iostream>

using namespace FEAT;
using namespace Util;

void Util::mpi_cout(FEAT::String string, std::function<bool (int, int)> func)
{
  Dist::Comm comm(Dist::Comm::world());

  if (func(comm.rank(), comm.size()))
  {
    std::cout << string;
  }
}
