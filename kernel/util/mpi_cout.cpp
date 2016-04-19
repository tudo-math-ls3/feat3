#include<kernel/util/mpi_cout.hpp>

#include <iostream>

using namespace FEAST;
using namespace Util;

void Util::mpi_cout(FEAST::String string,
std::function<bool (Index, Index)> func)
{
  Index rank(Comm::rank());
  Index ranks(Comm::size());

  if (func(rank, ranks))
  {
    std::cout<<string;
  }
}
