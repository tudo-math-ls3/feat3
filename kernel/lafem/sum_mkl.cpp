// includes, FEAST
#include <kernel/lafem/sum.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;

void Sum<Algo::MKL>::value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & x, const DenseVector<Mem::Main, float> & y)
{
    vsAdd(x.size(), x.elements(), y.elements(), r.elements());
}

void Sum<Algo::MKL>::value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & x, const DenseVector<Mem::Main, double> & y)
{
    vdAdd(x.size(), x.elements(), y.elements(), r.elements());
}
