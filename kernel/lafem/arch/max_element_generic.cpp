// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/arch/max_element.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template Index MaxElement<Mem::Main>::value_generic(const float * const, const Index);
template Index MaxElement<Mem::Main>::value_generic(const double * const, const Index);
