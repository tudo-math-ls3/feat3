// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/apply.hpp>

#include <cstring>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template void Apply::csr_generic(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
template void Apply::csr_generic(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);
template void Apply::csr_generic(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const bool);
template void Apply::csr_generic(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const bool);

template void Apply::cscr_generic(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint64_t * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const Index, const bool);
template void Apply::cscr_generic(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint64_t * const, const std::uint64_t * const, const std::uint64_t * const, const Index, const Index, const Index, const Index, const bool);
template void Apply::cscr_generic(float *, const float, const float * const, const float, const float * const, const float * const, const std::uint32_t * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const Index, const bool);
template void Apply::cscr_generic(double *, const double, const double * const, const double, const double * const, const double * const, const std::uint32_t * const, const std::uint32_t * const, const std::uint32_t * const, const Index, const Index, const Index, const Index, const bool);

template void Apply::dense_generic(float *, const float, const float, const float * const, const float * const, const float * const, const Index, const Index);
template void Apply::dense_generic(double *, const double, const double, const double * const, const double * const, const double * const, const Index, const Index);
