#include <kernel/archs.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/exception.hpp>

using namespace FEAST;

std::ostream & FEAST::operator<< (std::ostream & left, Archs::TagValue value)
{
  switch (value)
  {
    case Archs::tv_nil:
      left << Archs::None::name();
      break;

    case Archs::tv_cpu:
      left << Archs::CPU::name();
      break;

    case Archs::tv_generic:
      left << Archs::Generic::name();
      break;

    case Archs::tv_gpu:
      left << Archs::GPU::name();
      break;

    default:
      left << "Unknown Backend";
      break;
  }

  //throw InternalError("Unexpected value for Arch::TagValue '" + stringify(long(value)) + "'");

  return left;
}
