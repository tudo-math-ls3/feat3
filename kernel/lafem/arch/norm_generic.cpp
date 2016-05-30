// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/norm.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template float Norm2<Mem::Main>::value_generic(const float * const, const Index);
template double Norm2<Mem::Main>::value_generic(const double * const, const Index);
