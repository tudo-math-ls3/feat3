// includes, FEAST
#include <kernel/lafem/norm.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;

float Norm2<Algo::MKL>::value(const DenseVector<Mem::Main, float> & x)
{
    return cblas_snrm2(x.size(), x.elements(), 1);
}

double Norm2<Algo::MKL>::value(const DenseVector<Mem::Main, double> & x)
{
    return cblas_dnrm2(x.size(), x.elements(), 1);
}
