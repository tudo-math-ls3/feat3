// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/dot_product.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::LAFEM::Arch;

template float DotProduct<Mem::Main>::value_generic(const float * const, const float * const, const Index);
template double DotProduct<Mem::Main>::value_generic(const double * const, const double * const, const Index);

template float TripleDotProduct<Mem::Main>::value_generic(const float * const, const float * const, const float * const, const Index);
template double TripleDotProduct<Mem::Main>::value_generic(const double * const, const double * const, const double * const, const Index);
