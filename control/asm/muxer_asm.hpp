// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2024 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/muxer.hpp>
#include <kernel/assembly/mirror_assembler.hpp>

#include <control/domain/domain_control.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Asm
    {
      /**
       * \brief Assembles a scalar or blocked muxer for a system level
       *
       * \attention
       * The underlying vector type of the muxer to be assembled may be a LAFEM::DenseVector
       * or a LAFEM::DenseVectorBlocked, but nothing else!
       *
       * \param[in] virt_lvl_c
       * The virtual coarse domain level for which the muxer is to be assembled
       *
       * \param[in] space
       * The space for which the muxer is to be assembled
       *
       * \param[out] muxer
       * The muxer that is to be assembled
       */
      template<typename DomainLevel_, typename SpaceLambda_, typename Muxer_>
      void asm_muxer(const Domain::VirtualLevel<DomainLevel_>& virt_lvl_c, SpaceLambda_&& space_lambda, Muxer_& muxer)
      {
        typedef typename Muxer_::LocalVectorType VectorType;
        typedef typename Muxer_::MirrorType MirrorType;

        // assemble muxer parent
        if(virt_lvl_c.is_parent())
        {
          XASSERT(virt_lvl_c.is_child());

          const auto& layer_c = virt_lvl_c.layer_c();
          const DomainLevel_& level_p = virt_lvl_c.level_p();

          // loop over all children
          for(Index i(0); i < layer_c.child_count(); ++i)
          {
            const auto* child = level_p.find_patch_part(int(i));
            XASSERT(child != nullptr);
            MirrorType child_mirror;
            Assembly::MirrorAssembler::assemble_mirror(child_mirror, *space_lambda(level_p), *child);
            muxer.push_child(std::move(child_mirror));
          }
        }

        // assemble muxer child
        if(virt_lvl_c.is_child())
        {
          const auto& layer_c = virt_lvl_c.layer_c();
          const DomainLevel_& level_c = virt_lvl_c.level_c();

          Index num_dofs = space_lambda(level_c)->get_num_dofs();

          // set parent and sibling comms
          muxer.set_parent(
            layer_c.sibling_comm_ptr(),
            layer_c.get_parent_rank(),
            MirrorType::make_identity(num_dofs)
          );

          // compile muxer
          muxer.compile(VectorType(num_dofs));
        }
      }

      /**
       * \brief Builds a system muxer based on TupleMirror from 2 separate muxers
       *
       * \attention
       * The system muxer still has to be compiled after this function returns!
       *
       * \param[out] muxer_sys
       * The 2-component system muxer that is to be build
       *
       * \param[in] tmpl_vec
       * A \transient reference to a template vector that is allocated to correct size, but its contents are irrelevant.
       * Note: You can use the system gate's frequency vector as a template vector here.
       *
       * \param[in] muxer_0
       * A \transient reference to the first component muxer
       *
       * \param[in] muxer_1
       * A \transient reference to the second component muxer
       *
       * \note This function is copy-&-paste-friendly, so if you require an overload that will build
       * a gate from three or more components, simply add another overload of this function and add
       * another component gate.
       */
      template<typename SysMuxer_, typename Muxer0_, typename Muxer1_>
      void build_muxer_tuple(SysMuxer_& muxer_sys, const typename SysMuxer_::LocalVectorType& tmpl_vec,
        const Muxer0_& muxer_0, const Muxer1_& muxer_1)
      {
        typedef typename SysMuxer_::MirrorType SystemMirrorType;

        // the sibling communicator must be equal for all muxers
        XASSERTM(muxer_0.get_sibling_comm() == muxer_1.get_sibling_comm(), "muxers have incompatible sibling communicators");

        // do we actually have a parent?
        if(muxer_0.get_sibling_comm() != nullptr)
        {
          // build the parent mirror
          SystemMirrorType parent_mirror;
          parent_mirror.template at<0>().clone(muxer_0.get_parent_mirror(), LAFEM::CloneMode::Shallow);
          parent_mirror.template at<1>().clone(muxer_1.get_parent_mirror(), LAFEM::CloneMode::Shallow);

          // set the parent
          muxer_sys.set_parent(muxer_0.get_sibling_comm(), muxer_0.get_parent_rank(), std::move(parent_mirror));
        }

        // get the child mirrors
        const auto& child_mirrors_0 = muxer_0.get_child_mirrors();
        const auto& child_mirrors_1 = muxer_1.get_child_mirrors();

        // all child mirror vectors must have the same size
        XASSERT(child_mirrors_0.size() == child_mirrors_1.size());

        // loop over all child ranks
        for(std::size_t i(0); i < child_mirrors_0.size(); ++i)
        {
          // build the child mirror
          SystemMirrorType child_mirror;
          child_mirror.template at<0>().clone(child_mirrors_0[i], LAFEM::CloneMode::Shallow);
          child_mirror.template at<1>().clone(child_mirrors_1[i], LAFEM::CloneMode::Shallow);

          // child mirror cannot be empty
          XASSERTM(!child_mirror.empty(), "invalid empty mirror");
          muxer_sys.push_child(std::move(child_mirror));
        }

        // compile
        muxer_sys.compile(tmpl_vec);
      }

      /**
       * \brief Builds a system muxer based on TupleMirror from 3 separate muxers
       *
       * \attention
       * The system muxer still has to be compiled after this function returns!
       *
       * \param[out] muxer_sys
       * The 2-component system muxer that is to be build
       *
       * \param[in] tmpl_vec
       * A \transient reference to a template vector that is allocated to correct size, but its contents are irrelevant.
       * Note: You can use the system gate's frequency vector as a template vector here.
       *
       * \param[in] muxer_0
       * A \transient reference to the first component muxer
       *
       * \param[in] muxer_1
       * A \transient reference to the second component muxer
       *
       * \param[in] muxer_2
       * A \transient reference to the third component muxer
       *
       * \note This function is copy-&-paste-friendly, so if you require an overload that will build
       * a gate from three or more components, simply add another overload of this function and add
       * another component gate.
       */
      template<typename SysMuxer_, typename Muxer0_, typename Muxer1_, typename Muxer2_>
      void build_muxer_tuple(SysMuxer_& muxer_sys, const typename SysMuxer_::LocalVectorType& tmpl_vec,
        const Muxer0_& muxer_0, const Muxer1_& muxer_1, const Muxer2_& muxer_2)
      {
        typedef typename SysMuxer_::MirrorType SystemMirrorType;

        // the sibling communicator must be equal for all muxers
        XASSERTM(muxer_0.get_sibling_comm() == muxer_1.get_sibling_comm(), "muxers have incompatible sibling communicators");
        XASSERTM(muxer_0.get_sibling_comm() == muxer_2.get_sibling_comm(), "muxers have incompatible sibling communicators");

        // do we actually have a parent?
        if(muxer_0.get_sibling_comm() != nullptr)
        {
          // build the parent mirror
          SystemMirrorType parent_mirror;
          parent_mirror.template at<0>().clone(muxer_0.get_parent_mirror(), LAFEM::CloneMode::Shallow);
          parent_mirror.template at<1>().clone(muxer_1.get_parent_mirror(), LAFEM::CloneMode::Shallow);
          parent_mirror.template at<2>().clone(muxer_2.get_parent_mirror(), LAFEM::CloneMode::Shallow);

          // set the parent
          muxer_sys.set_parent(muxer_0.get_sibling_comm(), muxer_0.get_parent_rank(), std::move(parent_mirror));
        }

        // get the child mirrors
        const auto& child_mirrors_0 = muxer_0.get_child_mirrors();
        const auto& child_mirrors_1 = muxer_1.get_child_mirrors();
        const auto& child_mirrors_2 = muxer_2.get_child_mirrors();

        // all child mirror vectors must have the same size
        XASSERT(child_mirrors_0.size() == child_mirrors_1.size());
        XASSERT(child_mirrors_0.size() == child_mirrors_2.size());

        // loop over all child ranks
        for(std::size_t i(0); i < child_mirrors_0.size(); ++i)
        {
          // build the child mirror
          SystemMirrorType child_mirror;
          child_mirror.template at<0>().clone(child_mirrors_0[i], LAFEM::CloneMode::Shallow);
          child_mirror.template at<1>().clone(child_mirrors_1[i], LAFEM::CloneMode::Shallow);
          child_mirror.template at<2>().clone(child_mirrors_2[i], LAFEM::CloneMode::Shallow);

          // child mirror cannot be empty
          XASSERTM(!child_mirror.empty(), "invalid empty mirror");
          muxer_sys.push_child(std::move(child_mirror));
        }

        // compile
        muxer_sys.compile(tmpl_vec);
      }
    } // namespace Asm
  } // namespace Control
} // namespace FEAT
