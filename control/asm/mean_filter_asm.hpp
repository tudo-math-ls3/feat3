// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/global/mean_filter.hpp>
#include <kernel/global/gate.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Asm
    {
      /**
       * \brief Assembles a global mean filter
       *
       * \param[in] gate
       * A \transient reference to the gate of the space to assemble the filter for
       *
       * \param[in] space
       * A \transient reference to the space to assemble the filter for
       *
       * \param[in] cubature_name
       * The name of the cubature formula to use.
       */
      template<typename Gate_, typename Space_>
      Global::MeanFilter<typename Gate_::DataType, typename Gate_::IndexType> asm_mean_filter(Gate_& gate, const Space_& space, const String& cubature_name)
      {
        // get a clone of the frequencies vector
        auto vec_f = gate.get_freqs().clone(LAFEM::CloneMode::Deep);
        auto vec_v = gate.get_freqs().clone(LAFEM::CloneMode::Deep);
        auto vec_w = gate.get_freqs().clone(LAFEM::CloneMode::Deep);

        // assemble the mean filter
        Assembly::MeanFilterAssembler::assemble(vec_v, vec_w, space, cubature_name);

        // synchronize the vectors
        gate.sync_1(vec_v);
        gate.sync_0(vec_w);

        // build the mean filter
        return Global::MeanFilter<typename Gate_::DataType, typename Gate_::IndexType>(
          std::move(vec_v), std::move(vec_w), std::move(vec_f), gate.get_comm());
      }
    } // namespace Asm
  } // namespace Control
} // namespace FEAT
