// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/memory_pool.hpp>

// static member initialisation
std::map<void*, FEAT::Util::Intern::MemoryInfo> FEAT::MemoryPool<FEAT::Mem::Main>::_pool;
std::map<void*, FEAT::Util::Intern::MemoryInfo> FEAT::MemoryPool<FEAT::Mem::Main>::_pinned_pool;
