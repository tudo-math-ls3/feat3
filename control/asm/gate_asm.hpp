// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef CONTROL_ASSEMBLY_GATE_ASM_HPP
#define CONTROL_ASSEMBLY_GATE_ASM_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/assembly/mirror_assembler.hpp>

#include <control/domain/domain_control.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Asm
    {
      /**
       * \brief Assembles a scalar or blocked gate for a system level
       *
       * \attention
       * The underlying vector type of the gate to be assembled may be a LAFEM::DenseVector
       * or a LAFEM::DenseVectorBlocked, but nothing else!
       *
       * \param[in] virt_lvl
       * The virtual domain level on which the gate is to be assembled
       *
       * \param[in] space
       * The space for which the gate is to be assembled
       *
       * \param[out] gate
       * The gate that is to be assembled
       *
       * \param[in] skip_empty_halos
       * Specifies whether empty halos should be skipped; if set to false the program execution
       * will be aborted if an empty halo is encountered.
       */
      template<typename DomainLevel_, typename Space_, typename Gate_>
      void asm_gate(const Domain::VirtualLevel<DomainLevel_>& virt_lvl,
        const Space_& space, Gate_& gate, bool skip_empty_halos)
      {
        typedef typename Gate_::LocalVectorType VectorType;
        typedef typename Gate_::MirrorType MirrorType;

        const auto& dom_level = virt_lvl.level();
        const auto& dom_layer = virt_lvl.layer();

        // set the gate comm
        gate.set_comm(dom_layer.comm_ptr());

        // loop over all ranks
        for(Index i(0); i < dom_layer.neighbor_count(); ++i)
        {
          int rank = dom_layer.neighbor_rank(i);

          // try to find our halo
          auto* halo = dom_level.find_halo_part(rank);
          if(halo == nullptr)
          {
            if(skip_empty_halos)
              continue;
            XABORTM("empty halo for neighbor rank " + stringify(rank) + " detected");
          }

          // assemble the mirror
          MirrorType mirror;
          Assembly::MirrorAssembler::assemble_mirror(mirror, space, *halo);

          // push mirror into gate
          if(!mirror.empty())
            gate.push(rank, std::move(mirror));
        }

        // create local template vector
        VectorType tmpl_s(space.get_num_dofs());

        // compile gate
        gate.compile(std::move(tmpl_s));
      }

      /**
       * \brief Builds a system gate based on TupleMirror from 2 separate gates
       *
       * \param[out] gate_sys
       * The 2-component system gate that is to be build
       *
       * \param[in] gate_0
       * A \transient reference to the first component gate
       *
       * \param[in] gate_1
       * A \transient reference to the second component gate
       *
       * \note This function is copy-&-paste-friendly, so if you require an overload that will build
       * a gate from three or more components, simply add another overload of this function and add
       * another component gate.
       */
      template<typename SysGate_, typename Gate0_, typename Gate1_>
      void build_gate_tuple(SysGate_& gate_sys, const Gate0_& gate_0, const Gate1_& gate_1)
      {
        typedef typename SysGate_::LocalVectorType LocalSystemVectorType;
        typedef typename SysGate_::MirrorType SystemMirrorType;

        // set the communicator; this one must be identical for all gates
        XASSERTM(gate_0.get_comm() == gate_1.get_comm(), "gates have incompatible communicators");
        gate_sys.set_comm(gate_0.get_comm());

        // collect all neighbor ranks; these may differ for each gate
        std::set<int> neighbors;
        const std::vector<int>& ranks_0 = gate_0.get_ranks();
        const std::vector<int>& ranks_1 = gate_1.get_ranks();
        for(auto r : ranks_0)
          neighbors.insert(r);
        for(auto r : ranks_1)
          neighbors.insert(r);

        // get the mirrors
        const auto& mirrors_0 = gate_0.get_mirrors();
        const auto& mirrors_1 = gate_1.get_mirrors();

        // build a temporary vector from the gate frequency vectors
        LocalSystemVectorType sys_loc;
        sys_loc.template at<0>().clone(gate_0.get_freqs(), LAFEM::CloneMode::Layout);
        sys_loc.template at<1>().clone(gate_1.get_freqs(), LAFEM::CloneMode::Layout);

        // loop over all neighbor ranks
        for(auto rank : neighbors)
        {
          // create a mirror with empty sub-mirrors of correct size first
          SystemMirrorType sys_mir = SystemMirrorType::make_empty(sys_loc);
          for(std::size_t i(0); i < ranks_0.size(); ++i)
          {
            if(ranks_0[i] == rank)
            {
              sys_mir.template at<0>().clone(mirrors_0[i], LAFEM::CloneMode::Shallow);
              break;
            }
          }
          for(std::size_t i(0); i < ranks_1.size(); ++i)
          {
            if(ranks_1[i] == rank)
            {
              sys_mir.template at<1>().clone(mirrors_1[i], LAFEM::CloneMode::Shallow);
              break;
            }
          }
          // sys_mir cannot be empty, since at least one of the above gates must have a mirror for this neighbor
          XASSERTM(!sys_mir.empty(), "invalid empty mirror");
          gate_sys.push(rank, std::move(sys_mir));
        }

        // compile the system gate
        gate_sys.compile(std::move(sys_loc));
      }

      /**
       * \brief Builds a system gate based on TupleMirror from 3 separate gates
       *
       * \param[out] gate_sys
       * The 2-component system gate that is to be build
       *
       * \param[in] gate_0
       * A \transient reference to the first component gate
       *
       * \param[in] gate_1
       * A \transient reference to the second component gate
       *
       * \param[in] gate_2
       * A \transient reference to the third component gate
       *
       * \note This function is copy-&-paste-friendly, so if you require an overload that will build
       * a gate from three or more components, simply add another overload of this function and add
       * another component gate.
       */
      template<typename SysGate_, typename Gate0_, typename Gate1_, typename Gate2_>
      void build_gate_tuple(SysGate_& gate_sys, const Gate0_& gate_0, const Gate1_& gate_1, const Gate2_& gate_2)
      {
        typedef typename SysGate_::LocalVectorType LocalSystemVectorType;
        typedef typename SysGate_::MirrorType SystemMirrorType;

        // set the communicator; this one must be identical for all gates
        XASSERTM(gate_0.get_comm() == gate_1.get_comm(), "gates have incompatible communicators");
        XASSERTM(gate_0.get_comm() == gate_2.get_comm(), "gates have incompatible communicators");
        gate_sys.set_comm(gate_0.get_comm());

        // collect all neighbor ranks; these may differ for each gate
        std::set<int> neighbors;
        const std::vector<int>& ranks_0 = gate_0.get_ranks();
        const std::vector<int>& ranks_1 = gate_1.get_ranks();
        const std::vector<int>& ranks_2 = gate_2.get_ranks();
        for(auto r : ranks_0)
          neighbors.insert(r);
        for(auto r : ranks_1)
          neighbors.insert(r);
        for(auto r : ranks_2)
          neighbors.insert(r);

        // get the mirrors
        const auto& mirrors_0 = gate_0.get_mirrors();
        const auto& mirrors_1 = gate_1.get_mirrors();
        const auto& mirrors_2 = gate_2.get_mirrors();

        // build a temporary vector from the gate frequency vectors
        LocalSystemVectorType sys_loc;
        sys_loc.template at<0>().clone(gate_0.get_freqs(), LAFEM::CloneMode::Layout);
        sys_loc.template at<1>().clone(gate_1.get_freqs(), LAFEM::CloneMode::Layout);
        sys_loc.template at<2>().clone(gate_2.get_freqs(), LAFEM::CloneMode::Layout);

        // loop over all neighbor ranks
        for(auto rank : neighbors)
        {
          // create a mirror with empty sub-mirrors of correct size first
          SystemMirrorType sys_mir = SystemMirrorType::make_empty(sys_loc);
          for(std::size_t i(0); i < ranks_0.size(); ++i)
          {
            if(ranks_0[i] == rank)
            {
              sys_mir.template at<0>().clone(mirrors_0[i], LAFEM::CloneMode::Shallow);
              break;
            }
          }
          for(std::size_t i(0); i < ranks_1.size(); ++i)
          {
            if(ranks_1[i] == rank)
            {
              sys_mir.template at<1>().clone(mirrors_1[i], LAFEM::CloneMode::Shallow);
              break;
            }
          }
          for(std::size_t i(0); i < ranks_2.size(); ++i)
          {
            if(ranks_2[i] == rank)
            {
              sys_mir.template at<2>().clone(mirrors_2[i], LAFEM::CloneMode::Shallow);
              break;
            }
          }
          // sys_mir cannot be empty, since at least one of the above gates must have a mirror for this neighbor
          XASSERTM(!sys_mir.empty(), "invalid empty mirror");
          gate_sys.push(rank, std::move(sys_mir));
        }

        // compile the system gate
        gate_sys.compile(std::move(sys_loc));
      }
    } // namespace Asm
  } // namespace Control
} // namespace FEAT

#endif // CONTROL_ASSEMBLY_GATE_ASM_HPP
