// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/max_element.hpp>

#include <mkl.h>


using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

Index MaxElement<Mem::Main>::value_mkl(const float * const x, const Index size)
{
  return cblas_isamax((MKL_INT)size, x, 1);
}

Index MaxElement<Mem::Main>::value_mkl(const double * const x, const Index size)
{
  return cblas_idamax((MKL_INT)size, x, 1);
}
