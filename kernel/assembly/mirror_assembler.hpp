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
        static Index count(Index& contribs, const Space_& space, const MeshPart_& mesh_part)
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
            int num_assign(dof_assign.get_num_assigned_dofs());
            for(int j(0); j < num_assign; ++j)
            {
              contribs += Index(dof_assign.get_num_contribs(j));
            }
            count += Index(num_assign);
            dof_assign.finish();
          }

          return count;
        }

        static Index fill(Index*& ptr, Index idx[], Index offset, const Space_& space, const MeshPart_& mesh_part)
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
            for(int j(0); j < num_assign; ++j, ++ptr)
            {
              *ptr = offset;
              int num_contribs(dof_assign.get_num_contribs(j));
              for(int k(0); k < num_contribs; ++k, ++offset)
              {
                idx[offset] = dof_assign.get_index(j, k);
              }
            }
            *ptr = offset;
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
        static Index count(Index& contribs, const Space_& space, const MeshPart_& mesh_part)
        {
          // recursive call
          return DofMirrorHelpWrapper<Space_, MeshPart_, shape_dim_ - 1>::count(contribs, space, mesh_part) +
            DofMirrorHelper<Space_, MeshPart_, shape_dim_>::count(contribs, space, mesh_part);
        }

        static Index fill(Index*& ptr, Index idx[], const Space_& space, const MeshPart_& mesh_part)
        {
          Index offset =  DofMirrorHelpWrapper<Space_, MeshPart_, shape_dim_ - 1>::fill(ptr, idx, space, mesh_part);
          return DofMirrorHelper<Space_, MeshPart_, shape_dim_>::fill(ptr, idx, offset, space, mesh_part);
        }
      };

      template<
        typename Space_,
        typename MeshPart_>
      struct DofMirrorHelpWrapper<Space_, MeshPart_, 0>
      {
        static Index count(Index& contribs, const Space_& space, const MeshPart_& mesh_part)
        {
          return DofMirrorHelper<Space_, MeshPart_, 0>::count(contribs, space, mesh_part);
        }

        static Index fill(Index*& ptr, Index idx[], const Space_& space, const MeshPart_& mesh_part)
        {
          return DofMirrorHelper<Space_, MeshPart_, 0>::fill(ptr, idx, 0, space, mesh_part);
        }
      };

      template<
        typename Space_,
        int shape_dim_ = Space_::shape_dim>
      struct MaxDofContrib
      {
        static int value(const Space_& space)
        {
          int max_contrib(MaxDofContrib<Space_, shape_dim_-1>::value(space));
          typename Space_::template DofAssignment<shape_dim_>::Type dof_assign(space);
          return std::max(max_contrib, dof_assign.get_max_contribs());
        }
      };

      template<typename Space_>
      struct MaxDofContrib<Space_, 0>
      {
        static int value(const Space_& space)
        {
          typename Space_::template DofAssignment<0>::Type dof_assign(space);
          return dof_assign.get_max_contribs();
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
       * \brief Assembles the Dof-Mirror adjacency graph.
       */
      template<
        typename Space_,
        typename MeshPart_>
      static Adjacency::Graph assemble_mirror_graph(
        const Space_& space,
        const MeshPart_& mesh_part)
      {
        // ensure that the space has not more than 1 Dof contribution...
        if(Intern::MaxDofContrib<Space_>::value(space) != 1)
          throw InternalError("Cannot compute Dof-Mirror graph: multiple DOF contributions!");

        Index contribs(0);
        Index count = Intern::DofMirrorHelpWrapper<Space_, MeshPart_>::count(contribs, space, mesh_part);
        Adjacency::Graph graph(count, space.get_num_dofs(), contribs);
        Index* ptr = graph.get_domain_ptr();
        Index* idx = graph.get_image_idx();
        Intern::DofMirrorHelpWrapper<Space_, MeshPart_>::fill(ptr, idx, space, mesh_part);
        return graph;
      }

      /**
       * \brief Assembles a Vector-Mirror from graph.
       *
       * \param[out] vec_mirror
       * The vector mirror that is to be assembled.
       *
       * \param[in] graph
       * The mirror adjacency graph.
       */
      template<typename MemType_, typename DataType_, typename IndexType_>
      static void assemble_mirror(LAFEM::VectorMirror<MemType_, DataType_, IndexType_>& vec_mirror, const Adjacency::Graph& graph)
      {
        typedef typename LAFEM::VectorMirror<MemType_, DataType_, IndexType_>::MirrorMatrixType MirrorMatrixType;
        vec_mirror.get_gather() = MirrorMatrixType(graph);
        vec_mirror.get_scatter() = vec_mirror.get_gather().transpose();

        vec_mirror.get_gather().format(DataType_(1));
        vec_mirror.get_scatter().format(DataType_(1));
      }

      /**
       * \brief Assembles a Vector-Mirror from a space and a mesh-part.
       *
       * \param[out] vec_mirror
       * The vector mirror that is to be assembled.
       *
       * \param[in] space
       * A reference to the finite element space to be used.
       *
       * \param[in] mesh_part
       * A reference to the mesh-part that is to be mirrored.
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
        assemble_mirror(vec_mirror, assemble_mirror_graph(space, mesh_part));
      }

      template<
        typename MemType_,
        typename DataType_,
        typename IndexType_,
        Index count_,
        typename Space_,
        typename MeshPart_>
      static void assemble_mirror(
        LAFEM::PowerMirror<LAFEM::VectorMirror<MemType_, DataType_, IndexType_>, count_>& vec_mirror,
        const Space_& space, const MeshPart_& mesh_part)
      {
        LAFEM::VectorMirror<MemType_, DataType_, IndexType_> submirror;
        assemble_mirror(submirror, assemble_mirror_graph(space, mesh_part));
        vec_mirror = LAFEM::PowerMirror<LAFEM::VectorMirror<MemType_, DataType_, IndexType_>, count_>(std::move(submirror));
      }
    }; // class MirrorAssembler
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_MIRROR_ASSEMBLER_HPP
