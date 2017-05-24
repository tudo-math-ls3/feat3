// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/lumping.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template void Lumping<Mem::Main>::csr_generic(float *, const float * const, const Index * const, const Index * const, const Index);
template void Lumping<Mem::Main>::csr_generic(double *, const double * const, const Index * const, const Index * const, const Index);

template void Lumping<Mem::Main>::bcsr_generic(float *, const float * const, const Index * const, const Index * const, const Index, const int, const int);
template void Lumping<Mem::Main>::bcsr_generic(double *, const double * const, const Index * const, const Index * const, const Index, const int, const int);

template void Lumping<Mem::Main>::ell_generic(float *, const float * const, const Index * const,
  const Index * const, const Index * const, const Index, const Index);
template void Lumping<Mem::Main>::ell_generic(double *, const double * const, const Index * const,
  const Index * const, const Index * const, const Index, const Index);
