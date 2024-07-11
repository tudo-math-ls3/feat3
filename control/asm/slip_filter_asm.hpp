// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef CONTROL_ASSEMBLY_SLIP_FILTER_SYNC_HPP
#define CONTROL_ASSEMBLY_SLIP_FILTER_SYNC_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/global/gate.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Asm
    {
      /**
       * \brief Assembles a slip filter for a set of mesh-parts given by their names
       *
       * \attention
       * The assembled slip filter still has to be synchronized by using the #sync_slip_filter() function
       *
       * \param[inout] filter
       * A \transient reference to the filter to be assembled
       *
       * \param[in] dom_level
       * A \transient reference to the domain level on which is to be assembled
       *
       * \param[in] space
       * A \transient reference to the space to be assemble on
       *
       * \param[in] mesh_part_names
       * A string containing all mesh part names on which the slip filter is to be assembled.
       * If set to "*", all mesh-parts with be added to the filter.
       */
      template<typename DataType_, typename IndexType_, int block_size_, typename DomainLevel_, typename Space_>
      void asm_slip_filter(LAFEM::SlipFilter<DataType_, IndexType_, block_size_>& filter,
        const DomainLevel_& dom_level, const Space_& space, const String& mesh_part_names)
      {
        // create slip-filter assembler
        Assembly::SlipFilterAssembler<typename DomainLevel_::TrafoType> slip_filter_asm(dom_level.trafo);

        std::deque<String> mp_names;
        if(mesh_part_names == "*")
          mp_names = dom_level.get_mesh_node()->get_mesh_part_names(true);
        else
          mp_names = mesh_part_names.split_by_whitespaces();

        // loop over all mesh parts
        for(const auto& mp_name : mp_names)
        {
          auto* mesh_part_node = dom_level.get_mesh_node()->find_mesh_part_node(mp_name);
          XASSERT(mesh_part_node != nullptr);

          // let's see if we have that mesh part
          // if it is nullptr, then our patch is not adjacent to that boundary part
          auto* mesh_part = mesh_part_node->get_mesh();
          if (mesh_part != nullptr)
          {
            // add to boundary assembler
            slip_filter_asm.add_mesh_part(*mesh_part);
          }
        }

        // finally, assemble the filter
        slip_filter_asm.assemble(filter, space);
      }

      /**
       * \brief Synchronizes a slip filter with a gate
       *
       * \param[in] gate
       * A \transient reference to the gate to be used for synchronization
       *
       * \param[inout] slip_filter
       * A \transient reference to the slip filter that is to be synchronized
       */
      template<typename Gate_, typename SlipFilter_>
      void sync_slip_filter(Gate_& gate, SlipFilter_& slip_filter)
      {
        auto& slip_filter_vector = slip_filter.get_filter_vector();
        auto tmp = gate.get_freqs().clone(LAFEM::CloneMode::Layout);
        tmp.format();

        if(slip_filter_vector.used_elements() > 0)
        {
          auto* tmp_elements = tmp.template elements<LAFEM::Perspective::native>();
          auto* sfv_elements = slip_filter_vector.template elements<LAFEM::Perspective::native>();

          // Copy sparse filter vector contents to DenseVector
          for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
          {
            Index idense(slip_filter_vector.indices()[isparse]);
            tmp_elements[idense] = sfv_elements[isparse];
          }

          gate.sync_0(tmp);

          // Copy sparse filter vector contents to DenseVector
          for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
          {
            Index idense(slip_filter_vector.indices()[isparse]);
            tmp_elements[idense].normalize();
            sfv_elements[isparse] = tmp_elements[idense];
          }
        }
        else
        {
          gate.sync_0(tmp);
        }
      }
    } // namespace Asm
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_ASSEMBLY_SLIP_FILTER_SYNC_HPP
