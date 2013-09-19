// includes, FEAST
#include <kernel/lafem/difference.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;

void Difference<Algo::MKL>::value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & x, const DenseVector<Mem::Main, float> & y)
{
    vsSub((MKL_INT)x.size(), x.elements(), y.elements(), r.elements());
}

void Difference<Algo::MKL>::value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & x, const DenseVector<Mem::Main, double> & y)
{
    vdSub((MKL_INT)x.size(), x.elements(), y.elements(), r.elements());
}
