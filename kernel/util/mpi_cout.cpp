#include<kernel/util/mpi_cout.hpp>

#include <iostream>

using namespace FEAT;
using namespace Util;

void Util::mpi_cout(FEAT::String string,
std::function<bool (Index, Index)> func)
{
  Index rank(Comm::rank());
  Index ranks(Comm::size());

  if (func(rank, ranks))
  {
    std::cout<<string;
  }
}
