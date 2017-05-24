// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/row_norm.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template void RowNorm<Mem::Main>::csr_generic_norm2(float*, const float* const, const Index* const,
const Index* const, const Index);
template void RowNorm<Mem::Main>::csr_generic_norm2(double*, const double* const, const Index* const,
const Index* const, const Index);

template void RowNorm<Mem::Main>::csr_generic_norm2sqr(float*, const float* const, const Index* const,
const Index* const, const Index);
template void RowNorm<Mem::Main>::csr_generic_norm2sqr(double*, const double* const, const Index* const,
const Index * const, const Index);

template void RowNorm<Mem::Main>::csr_generic_scaled_norm2sqr(float*, const float* const, const float* const, const Index* const,
const Index* const, const Index);
template void RowNorm<Mem::Main>::csr_generic_scaled_norm2sqr(double*, const double* const, const double* const, const Index* const,
const Index* const, const Index);

template void RowNorm<Mem::Main>::bcsr_generic_norm2(float *, const float * const, const Index * const, const Index * const, const Index, const int, const int);
template void RowNorm<Mem::Main>::bcsr_generic_norm2(double *, const double * const, const Index * const, const Index * const, const Index, const int, const int);

template void RowNorm<Mem::Main>::bcsr_generic_norm2sqr(float*, const float* const, const Index* const,
const Index* const, const Index, const int, const int);
template void RowNorm<Mem::Main>::bcsr_generic_norm2sqr(double*, const double* const, const Index* const,
const Index * const, const Index, const int, const int);

template void RowNorm<Mem::Main>::bcsr_generic_scaled_norm2sqr(float*, const float* const, const float* const,
const Index* const, const Index* const, const Index, const int, const int);
template void RowNorm<Mem::Main>::bcsr_generic_scaled_norm2sqr(double*, const double* const, const double* const,
const Index* const, const Index* const, const Index, const int, const int);

template void RowNorm<Mem::Main>::ell_generic_norm2(float*, const float* const, const Index* const,
  const Index* const, const Index* const, const Index, const Index);
template void RowNorm<Mem::Main>::ell_generic_norm2(double*, const double* const, const Index* const,
  const Index* const, const Index* const, const Index, const Index);

template void RowNorm<Mem::Main>::ell_generic_norm2sqr(float*, const float* const, const Index* const,
  const Index* const, const Index* const, const Index, const Index);
template void RowNorm<Mem::Main>::ell_generic_norm2sqr(double*, const double* const, const Index* const,
  const Index* const, const Index* const, const Index, const Index);

template void RowNorm<Mem::Main>::ell_generic_scaled_norm2sqr(float*, const float* const, const float* const,
const Index* const, const Index* const, const Index* const, const Index, const Index);
template void RowNorm<Mem::Main>::ell_generic_scaled_norm2sqr(double*, const double* const, const double* const,
const Index* const, const Index* const, const Index* const, const Index, const Index);
