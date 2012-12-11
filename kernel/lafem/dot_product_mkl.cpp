// includes, FEAST
#include <kernel/lafem/dot_product.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;

float DotProduct<Algo::MKL>::value(const DenseVector<Mem::Main, float> & x, const DenseVector<Mem::Main, float> & y)
{
    return cblas_sdot(x.size(), x.elements(), 1, y.elements(), 1);
}

double DotProduct<Algo::MKL>::value(const DenseVector<Mem::Main, double> & x, const DenseVector<Mem::Main, double> & y)
{
    return cblas_ddot(x.size(), x.elements(), 1, y.elements(), 1);
}
