#include <kernel/hornet/archs.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/exception.hpp>

using namespace FEAST;
using namespace FEAST::HOrNET;

const String Archs::CPU::name = "cpu";

std::ostream & FEAST::HOrNET::operator<< (std::ostream & left, Archs::TagValue value)
{
  do
  {
    switch (value)
    {
      case Archs::tv_cpu:
        left << Archs::CPU::name;
        continue;

      default:
        left << "Uknown Backend";
        continue;
    }

    throw InternalError("Unexpected value for Arch::TagValue '" + stringify(long(value)) + "'");
  } while (false);

  return left;
}
