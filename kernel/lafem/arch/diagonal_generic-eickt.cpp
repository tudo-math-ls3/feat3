// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/diagonal.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::LAFEM::Arch;

template void Diagonal::csr_generic(std::uint64_t *, const std::uint64_t * const, const std::uint64_t * const, const Index);
template void Diagonal::csr_generic(std::uint32_t *, const std::uint32_t * const, const std::uint32_t * const, const Index);
