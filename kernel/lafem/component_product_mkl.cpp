// includes, FEAST
#include <kernel/lafem/component_product.hpp>

#include <mkl.h>

using namespace FEAST;
using namespace FEAST::LAFEM;

void ComponentProduct<Algo::MKL>::value(DenseVector<Mem::Main, float> & r, const DenseVector<Mem::Main, float> & x, const DenseVector<Mem::Main, float> & y)
{
    vsMul(x.size(), x.elements(), y.elements(), r.elements());
}

void ComponentProduct<Algo::MKL>::value(DenseVector<Mem::Main, double> & r, const DenseVector<Mem::Main, double> & x, const DenseVector<Mem::Main, double> & y)
{
    vdMul(x.size(), x.elements(), y.elements(), r.elements());
}
