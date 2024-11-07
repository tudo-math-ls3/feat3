// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/splitter.hpp>
#include <kernel/assembly/mirror_assembler.hpp>

#include <control/domain/domain_control.hpp>
#include <control/asm/muxer_asm.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Asm
    {
      /**
       * \brief Assembles a base splitter for a system level
       *
       * \tparam VectorType_
       * The underlying vector type: may be a LAFEM::DenseVector or a LAFEM::DenseVectorBlocked, but nothing else!
       *
       * \tparam MirrorType_
       * The underlying vector type: may be a LAFEM::VectorMirror compatible to VectorType_
       *
       * \param[inout] sys_lvl
       * The system level in which the base splitter is to be assembled
       *
       * \param[in] virt_lvl_c
       * The virtual domain level for which the base splitter is to be assembled
       *
       * \param[in] space
       * The space for which the base splitter is to be assembled
       *
       * \param[out] base_splitter
       * The base splitter that is to be assembled
       */
      template<typename DomainLevel_, typename SpaceLambda_, typename BaseSplitter_>
      void asm_splitter(const Domain::VirtualLevel<DomainLevel_>& virt_lvl,
        SpaceLambda_&& space_lambda, BaseSplitter_& base_splitter)
      {
        typedef typename BaseSplitter_::LocalVectorType VectorType;
        typedef typename BaseSplitter_::MirrorType MirrorType;

        // get the layer
        const auto& layer = virt_lvl.layer();

        // nothing to assemble?
        if(layer.comm().size() <= 1)
          return;

        // ensure that this virtual level has a base-mesh
        XASSERTM(virt_lvl.has_base(), "cannot assemble base splitter because level has no base");

        // is this the root process?
        if(layer.comm().rank() == 0)
        {
          const DomainLevel_& level_b = virt_lvl.level_b();

          // set parent vector template
          base_splitter.set_base_vector_template(VectorType(space_lambda(level_b)->get_num_dofs()));

          // assemble patch mirrors on root process
          for(int i(0); i < layer.comm().size(); ++i)
          {
            const auto* patch = level_b.find_patch_part(i);
            XASSERT(patch != nullptr);
            MirrorType patch_mirror;
            Assembly::MirrorAssembler::assemble_mirror(patch_mirror, *space_lambda(level_b), *patch);
            base_splitter.push_patch(std::move(patch_mirror));
          }
        }

        // assemble base splitter child
        const DomainLevel_& level = virt_lvl.level();

        Index num_dofs = space_lambda(level)->get_num_dofs();

        // set parent and sibling comms
        base_splitter.set_root(layer.comm_ptr(), 0, MirrorType::make_identity(num_dofs));

        // compile base splitter
        base_splitter.compile(VectorType(num_dofs));
      }

      /**
       * \brief Builds a system splitter based on TupleMirror from 2 separate splitters
       *
       * \attention
       * The system splitter still has to be compiled after this function returns!
       *
       * \param[out] splitter_sys
       * The 2-component system splitter that is to be build
       *
       * \param[in] tmpl_vec
       * A \transient reference to a template vector that is allocated to correct size, but its contents are irrelevant.
       * Note: You can use the system gate's frequency vector as a template vector here.
       *
       * \param[in] splitter_0
       * A \transient reference to the first component splitter
       *
       * \param[in] splitter_1
       * A \transient reference to the second component splitter
       *
       * \note This function is copy-&-paste-friendly, so if you require an overload that will build
       * a gate from three or more components, simply add another overload of this function and add
       * another component gate.
       */
      template<typename SysSplitter_, typename Splitter0_, typename Splitter1_>
      void build_splitter_tuple(SysSplitter_& splitter_sys, const typename SysSplitter_::LocalVectorType& tmpl_vec, const Splitter0_& splitter_0, const Splitter1_& splitter_1)
      {
        // build the internal tuple muxer
        typename SysSplitter_::MuxerType sys_muxer;
        build_muxer_tuple(sys_muxer, tmpl_vec, splitter_0.get_muxer(), splitter_1.get_muxer());

        // build the internal base vector template
        typename SysSplitter_::LocalVectorType base_vector_tmpl;
        base_vector_tmpl.template at<0>().clone(splitter_0.get_base_vector_template(), LAFEM::CloneMode::Shallow);
        base_vector_tmpl.template at<1>().clone(splitter_1.get_base_vector_template(), LAFEM::CloneMode::Shallow);

        // create the splitter
        splitter_sys = SysSplitter_(std::move(sys_muxer), std::move(base_vector_tmpl));
      }

      /**
       * \brief Builds a system splitter based on TupleMirror from 3 separate splitters
       *
       * \attention
       * The system splitter still has to be compiled after this function returns!
       *
       * \param[out] splitter_sys
       * The 2-component system splitter that is to be build
       *
       * \param[in] tmpl_vec
       * A \transient reference to a template vector that is allocated to correct size, but its contents are irrelevant.
       * Note: You can use the system gate's frequency vector as a template vector here.
       *
       * \param[in] splitter_0
       * A \transient reference to the first component splitter
       *
       * \param[in] splitter_1
       * A \transient reference to the second component splitter
       *
       * \param[in] splitter_2
       * A \transient reference to the third component splitter
       *
       * \note This function is copy-&-paste-friendly, so if you require an overload that will build
       * a gate from three or more components, simply add another overload of this function and add
       * another component gate.
       */
      template<typename SysSplitter_, typename Splitter0_, typename Splitter1_, typename Splitter2_>
      void build_splitter_tuple(SysSplitter_& splitter_sys, const typename SysSplitter_::LocalVectorType& tmpl_vec,
        const Splitter0_& splitter_0, const Splitter1_& splitter_1, const Splitter2_& splitter_2)
      {
        // build the internal tuple muxer
        typename SysSplitter_::MuxerType sys_muxer;
        build_muxer_tuple(sys_muxer, tmpl_vec, splitter_0.get_muxer(), splitter_1.get_muxer());

        // build the internal base vector template
        typename SysSplitter_::LocalVectorType base_vector_tmpl;
        base_vector_tmpl.template at<0>().clone(splitter_0.get_base_vector_template(), LAFEM::CloneMode::Shallow);
        base_vector_tmpl.template at<1>().clone(splitter_1.get_base_vector_template(), LAFEM::CloneMode::Shallow);
        base_vector_tmpl.template at<2>().clone(splitter_2.get_base_vector_template(), LAFEM::CloneMode::Shallow);

        // create the splitter
        splitter_sys = SysSplitter_(std::move(sys_muxer), std::move(base_vector_tmpl));
      }
    } // namespace Asm
  } // namespace Control
} // namespace FEAT
