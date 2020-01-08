// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/archs.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/exception.hpp>

using namespace FEAT;

std::ostream & FEAT::operator<< (std::ostream & left, Archs::TagValue value)
{
  switch (value)
  {
    case Archs::tv_none:
      left << Archs::None::name();
      break;

    case Archs::tv_serial:
      left << Archs::Serial::name();
      break;

    case Archs::tv_parallel:
      left << Archs::Parallel::name();
      break;

    default:
      left << "Unknown Backend";
      break;
  }

  return left;
}

std::ostream & FEAT::operator<< (std::ostream & left, Mem::TagValue value)
{
  switch (value)
  {
    case Mem::tv_main:
      left << Mem::Main::name();
      break;

    case Mem::tv_cuda:
      left << Mem::CUDA::name();
      break;

    default:
      left << "Unknown Backend";
      break;
  }

  return left;
}

std::ostream & FEAT::operator<< (std::ostream & left, Algo::TagValue value)
{
  switch (value)
  {
    case Algo::tv_generic:
      left << Algo::Generic::name();
      break;

    case Algo::tv_mkl:
      left << Algo::MKL::name();
      break;

    case Algo::tv_cuda:
      left << Algo::CUDA::name();
      break;

    default:
      left << "Unknown Backend";
      break;
  }

  return left;
}
