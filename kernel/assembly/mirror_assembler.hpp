// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_MIRROR_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_MIRROR_ASSEMBLER_HPP 1

// includes, FEAT
#include <kernel/assembly/symbolic_assembler.hpp>

// includes, FEAT-LAFEM
#include <kernel/lafem/matrix_mirror.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/power_mirror.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/lafem/sparse_matrix_banded.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Space_,
        typename MeshPart_,
        int shape_dim_>
      struct DofMirrorHelper
      {
        static Index count(const Space_& space, const MeshPart_& mesh_part)
        {
          // fetch the target set for this dimension
          const typename MeshPart_::template TargetSet<shape_dim_>::Type&
            target_set(mesh_part.template get_target_set<shape_dim_>());
          if(target_set.get_num_entities() <= 0)
            return 0;

          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_>::Type dof_assign(space);
          if(dof_assign.get_max_assigned_dofs() <= 0)
            return 0;

          // loop over all target indices
          Index num_entities = target_set.get_num_entities();
          Index count(0);
          for(Index i(0); i < num_entities; ++i)
          {
            dof_assign.prepare(target_set[i]);
            count += Index(dof_assign.get_num_assigned_dofs());
            dof_assign.finish();
          }

          return count;
        }

        template<typename IT_>
        static Index fill(IT_ idx[], Index offset, const Space_& space, const MeshPart_& mesh_part)
        {
          // fetch the target set for this dimension
          const typename MeshPart_::template TargetSet<shape_dim_>::Type&
            target_set(mesh_part.template get_target_set<shape_dim_>());
          if(target_set.get_num_entities() <= 0)
            return offset;

          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_>::Type dof_assign(space);
          if(dof_assign.get_max_assigned_dofs() <= 0)
            return offset;

          // loop over all target indices
          Index num_entities = target_set.get_num_entities();
          for(Index i(0); i < num_entities; ++i)
          {
            dof_assign.prepare(target_set[i]);
            int num_assign(dof_assign.get_num_assigned_dofs());
            for(int j(0); j < num_assign; ++j, ++offset)
            {
              idx[offset] = IT_(dof_assign.get_index(j));
            }
            dof_assign.finish();
          }

          return offset;
        }
      };

      template<
        typename Space_,
        typename MeshPart_,
        int shape_dim_ = Space_::shape_dim>
      struct DofMirrorHelpWrapper
      {
        static Index count(const Space_& space, const MeshPart_& mesh_part)
        {
          // recursive call
          return DofMirrorHelpWrapper<Space_, MeshPart_, shape_dim_ - 1>::count(space, mesh_part) +
            DofMirrorHelper<Space_, MeshPart_, shape_dim_>::count(space, mesh_part);
        }

        template<typename IT_>
        static Index fill(IT_ idx[], const Space_& space, const MeshPart_& mesh_part)
        {
          Index offset =  DofMirrorHelpWrapper<Space_, MeshPart_, shape_dim_ - 1>::fill(idx, space, mesh_part);
          return DofMirrorHelper<Space_, MeshPart_, shape_dim_>::fill(idx, offset, space, mesh_part);
        }
      };

      template<
        typename Space_,
        typename MeshPart_>
      struct DofMirrorHelpWrapper<Space_, MeshPart_, 0>
      {
        static Index count(const Space_& space, const MeshPart_& mesh_part)
        {
          return DofMirrorHelper<Space_, MeshPart_, 0>::count(space, mesh_part);
        }

        template<typename IT_>
        static Index fill(IT_ idx[], const Space_& space, const MeshPart_& mesh_part)
        {
          return DofMirrorHelper<Space_, MeshPart_, 0>::fill(idx, 0, space, mesh_part);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Dof-Mirror assembler class template.
     *
     * \author Peter Zajac
     */
    class MirrorAssembler
    {
    public:
      /**
       * \brief Assembles a VectorMirror from a space and a mesh-part.
       *
       * \param[out] vec_mirror
       * A \transient reference to the vector mirror that is to be assembled.
       *
       * \param[in] space
       * A \transient reference to the finite element space to be used.
       *
       * \param[in] mesh_part
       * A \transient reference to the mesh-part that is to be mirrored.
       */
      template<
        typename MemType_,
        typename DataType_,
        typename IndexType_,
        typename Space_,
        typename MeshPart_>
      static void assemble_mirror(
        LAFEM::VectorMirror<MemType_, DataType_, IndexType_>& vec_mirror,
        const Space_& space, const MeshPart_& mesh_part)
      {
        // count number of dofs in mirror
        const Index count = Intern::DofMirrorHelpWrapper<Space_, MeshPart_>::count(space, mesh_part);

        // allocate mirror in main memory
        LAFEM::VectorMirror<Mem::Main, DataType_, IndexType_> vmir(space.get_num_dofs(), count);

        // fill mirror indices
        if(count > Index(0))
        {
          Intern::DofMirrorHelpWrapper<Space_, MeshPart_>::fill(vmir.indices(), space, mesh_part);
        }

        // convert mirror
        vec_mirror.convert(vmir);
      }

      template<
        typename MemType_,
        typename DataType_,
        typename IndexType_,
        Index blocks_,
        typename Space_,
        typename MeshPart_>
      static void assemble_mirror(
        LAFEM::PowerMirror<LAFEM::VectorMirror<MemType_, DataType_, IndexType_>, blocks_>& vec_mirror,
        const Space_& space, const MeshPart_& mesh_part)
      {
        assemble_mirror(vec_mirror._sub_mirror, space, mesh_part);
      }
    }; // class MirrorAssembler
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_MIRROR_ASSEMBLER_HPP
