// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/diagonal.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template void Diagonal<Mem::Main>::csr_generic(float *, const float * const, const Index * const, const Index * const, const Index);
template void Diagonal<Mem::Main>::csr_generic(double *, const double * const, const Index * const, const Index * const, const Index);

template void Diagonal<Mem::Main>::ell_generic(float *, const float * const, const Index * const,
  const Index * const, const Index * const, const Index, const Index);
template void Diagonal<Mem::Main>::ell_generic(double *, const double * const, const Index * const,
  const Index * const, const Index * const, const Index, const Index);
