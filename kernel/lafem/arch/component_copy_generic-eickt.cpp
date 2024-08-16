// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/component_copy.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template void ComponentCopy::value_generic(float *, const float * const, const int, const int, const Index);
template void ComponentCopy::value_generic(double *, const double * const, const int, const int, const Index);
template void ComponentCopy::value_to_generic(const float * const, float *, const int, const int, const Index);
template void ComponentCopy::value_to_generic(const double * const, double *, const int, const int, const Index);
