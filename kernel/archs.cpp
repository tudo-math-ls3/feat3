#include <kernel/archs.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/exception.hpp>

using namespace FEAST;

const String Archs::Nil::name = "Nil";
const String Archs::CPU::name = "cpu";

std::ostream & FEAST::operator<< (std::ostream & left, Archs::TagValue value)
{
  do
  {
    switch (value)
    {
      case Archs::tv_nil:
        left << Archs::Nil::name;
        continue;

      case Archs::tv_cpu:
        left << Archs::CPU::name;
        continue;

      default:
        left << "Unknown Backend";
        continue;
    }

    throw InternalError("Unexpected value for Arch::TagValue '" + stringify(long(value)) + "'");
  } while (false);

  return left;
}
