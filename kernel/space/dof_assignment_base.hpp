// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_DOF_ASSIGNMENT_BASE_HPP
#define KERNEL_SPACE_DOF_ASSIGNMENT_BASE_HPP 1

// includes, FEAT
#include <kernel/space/base.hpp>

namespace FEAT
{
  namespace Space
  {
    /**
     * \brief Mapping of mesh entities to FE basis functions associated with them.
     *
     * In FEAT, every FE basis function is associated with mesh entities of exactly one shape dimension, i.e.
     * vertices, facets or cells. The DofMapping maps a mesh cell (meaning only the highest-dimensional entities of
     * a mesh) to a set of basis functions with nonempty support on that cell. It does, however, not contain any
     * information about which entities of that cell are associated with which basis function. This is done by the
     * DofAssignment.
     *
     * For the DofAssignment, mesh entities are open. This means i.e. the DofAssignment for an edge only knows about
     * basis functions directly associated with the edge and NOT about basis functions associated with the vertices
     * that represent the edge's end points.
     *
     * \see DofMappingBase
     *
     * \author Peter Zajac
     */
    template<
      typename Space_,
      int shape_dim_,
      typename DataType_>
    class DofAssignmentBase
    {
    public:
      /// space typedef
      typedef Space_ SpaceType;

    protected:
      /// space reference
      const SpaceType& _space;
      /// currently active cell index
      Index _cell_index;

      /// constructor
      explicit DofAssignmentBase(const SpaceType& space) :
        _space(space),
        _cell_index(~Index(0))
      {
      }

    public:
      void prepare(Index cell_index)
      {
        _cell_index = cell_index;
      }

      void finish()
      {
        _cell_index = ~Index(0);
      }

#ifdef DOXYGEN
      int get_max_assigned_dofs() const;

      int get_num_assigned_dofs() const;

      Index get_index(int assign_idx) const;
#endif // DOXYGEN
    }; // class DofAssignmentBase

    /// \cond internal
    namespace Intern
    {
      template<typename Tag_, template<typename, int> class Traits_, int cell_dim_>
      struct UniformDofAssignHelper
      {
        template<typename Mesh_>
        static Index dof_offset(const Mesh_& mesh)
        {
          return UniformDofAssignHelper<Tag_, Traits_, cell_dim_ - 1>::dof_offset(mesh)
            + mesh.get_num_entities(cell_dim_-1) * Index(Traits_<Tag_, cell_dim_ - 1>::count);
        }
      };

      template<typename Tag_, template<typename, int> class Traits_>
      struct UniformDofAssignHelper<Tag_, Traits_, 0>
      {
        template<typename Mesh_>
        static Index dof_offset(const Mesh_&)
        {
          return Index(0);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Uniform Dof-Assignment class template.
     *
     * \tparam Space_
     * The finite element space this dof-assignment is used by.
     *
     * \tparam shape_dim_
     * The dimension of the shape that this dof-assignment refers to.
     *
     * \tparam DataType_
     * The data-type to be used by the dof-assignment.
     *
     * \tparam DofTraits_
     * A dof-traits class template that defines the number of dofs per dimension.
     *
     * \tparam DofTag_
     * A tag class that is passed as a first parameter to the DofTraits_ class template.
     *
     * \author Peter Zajac
     */
    template<
      typename Space_,
      int shape_dim_,
      typename DataType_,
      template<typename Tag_, int dim_> class DofTraits_,
      typename DofTag_>
    class DofAssignmentUniform
      : public DofAssignmentBase<Space_, shape_dim_, DataType_>
    {
    public:
      typedef DofAssignmentBase<Space_, shape_dim_, DataType_> BaseClass;
      typedef Space_ SpaceType;
      typedef typename SpaceType::ShapeType ShapeType;

      // number of dofs for each cell
      static constexpr int dofs_per_cell = DofTraits_<DofTag_, shape_dim_>::count;

    protected:
      /// the offset of the first dof of this cell dimension
      Index _dof_offset;

    public:
      explicit DofAssignmentUniform(const SpaceType& space) :
        BaseClass(space),
        _dof_offset(Intern::UniformDofAssignHelper<DofTag_, DofTraits_, shape_dim_>::dof_offset(space.get_mesh()))
      {
      }

      /** \copydoc DofAssignmentBase::get_max_assigned_dofs() */
      int get_max_assigned_dofs() const
      {
        return dofs_per_cell;
      }

      /** \copydoc DofAssignmentBase::get_num_assigned_dofs() */
      int get_num_assigned_dofs() const
      {
        return dofs_per_cell;
      }

      /** \copydoc DofAssignmentBase::get_index() */
      Index get_index(int assign_idx) const
      {
        return _dof_offset + Index(dofs_per_cell) * this->_cell_index + Index(assign_idx);
      }
    };
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_DOF_ASSIGNMENT_BASE_HPP
