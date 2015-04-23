// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/defect.hpp>
#include <kernel/lafem/arch/difference.hpp>

#include <cstring>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template void Defect<Mem::Main>::ell_generic(float *, const float * const, const float * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const float * const, const Index, const Index);
template void Defect<Mem::Main>::ell_generic(double *, const double * const, const double * const, const unsigned long * const, const unsigned long * const, const unsigned long * const, const double * const, const Index, const Index);
template void Defect<Mem::Main>::ell_generic(float *, const float * const, const float * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const float * const, const Index, const Index);
template void Defect<Mem::Main>::ell_generic(double *, const double * const, const double * const, const unsigned int * const, const unsigned int * const, const unsigned int * const, const double * const, const Index, const Index);
