// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/difference.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

void Difference<Mem::Main, Algo::MKL>::value(float * r, const float * const x, const float * const y, const Index size)
{
  vsSub((MKL_INT)size, x, y, r);
}

void Difference<Mem::Main, Algo::MKL>::value(double * r, const double * const x, const double * const y, const Index size)
{
  vdSub((MKL_INT)size, x, y, r);
}
