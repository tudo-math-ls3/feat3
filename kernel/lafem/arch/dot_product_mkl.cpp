// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/dot_product.hpp>

#include <mkl.h>


using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

float DotProduct<Mem::Main, Algo::MKL>::value(const float * const x, const float * const y, const Index size)
{
    return cblas_sdot((MKL_INT)size, x, 1, y, 1);
}

double DotProduct<Mem::Main, Algo::MKL>::value(const double * const x, const double * const y, const Index size)
{
    return cblas_ddot((MKL_INT)size, x, 1, y, 1);
}
