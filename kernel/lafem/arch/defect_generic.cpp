// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/lafem/arch/difference.hpp>

#include <cstring>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;


template void Defect<Mem::Main, Algo::Generic>::csr(float *, const float * const, const float * const, const Index * const, const Index * const, const float * const, const Index);
template void Defect<Mem::Main, Algo::Generic>::csr(double *, const double * const, const double * const, const Index * const, const Index * const, const double * const, const Index);

template void Defect<Mem::Main, Algo::Generic>::ell(float *, const float * const, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index);
template void Defect<Mem::Main, Algo::Generic>::ell(double *,const double * const,  const double * const, const Index * const, const Index * const, const double * const, const Index, const Index);

template void Defect<Mem::Main, Algo::Generic>::coo(float *, const float * const, const float * const, const Index * const, const Index * const, const float * const, const Index, const Index);
template void Defect<Mem::Main, Algo::Generic>::coo(double *, const double * const, const double * const, const Index * const, const Index * const, const double * const, const Index, const Index);
