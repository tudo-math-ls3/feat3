#pragma once
#ifndef KERNEL_SPACE_DOF_MAPPING_BASE_HPP
#define KERNEL_SPACE_DOF_MAPPING_BASE_HPP 1

// includes, FEAT
#include <kernel/space/base.hpp>

namespace FEAT
{
  namespace Space
  {
    /**
     * \brief Mapping from mesh cells to Dof that have support on them
     *
     * This class acts as a base-class and interface documentation for Finite Element Dof mapping implementations.
     *
     * In FEAT, every FE basis function is associated with a mesh entity of a certain shape dimension, i.e.
     * vertices, facets or cells. The DofMapping maps a mesh cell (meaning only the highest-dimensional entities of
     * a mesh) to a set of basis functions with nonempty support on the closure of the cell. It does, however,
     * not contain any information about which entities of that cell are associated with which basis function. This
     * is done by the DofAssignment.
     *
     * \see DofAssignmentBase
     *
     * \tparam Space_
     * The finite element space that this Dof mapping is used by.
     *
     * \author Peter Zajac
     */
    template<typename Space_>
    class DofMappingBase
    {
    public:
      /// space typedef
      typedef Space_ SpaceType;
      /// shape typedef
      typedef typename SpaceType::ShapeType ShapeType;

    protected:
      /// space reference
      const SpaceType& _space;
      /// currently active cell index
      Index _cell_index;

      /// constructor
      explicit DofMappingBase(const SpaceType& space) :
        _space(space),
        _cell_index(~Index(0))
      {
      }

    public:

      /**
       * \brief Prepares the dof-mapping for a given cell.
       *
       * \param[in] cell_index
       * The index of the cell that is to be used by the dof-mapping.
       */
      void prepare(Index cell_index)
      {
        _cell_index = cell_index;
      }

      /**
       * \brief Releases the dof-mapping from the current cell.
       */
      void finish()
      {
        // reset cell index
        _cell_index = ~Index(0);
      }

#ifdef DOXYGEN
      /**
       * \brief Returns the number of local dofs.
       *
       * \returns
       * The total number of local degrees of freedom for the currently active cell.
       */
      int get_num_local_dofs() const;

      /**
       * \brief Return the number of global dofs.
       *
       * \returns
       * The total number of global degress of freedom for the finite element space.
       */
      Index get_num_global_dofs() const;
#endif // DOXYGEN

      /**
       * \brief Returns the maximum number of dof contributions.
       *
       * In the case of hanging vertices, values to a dof associated with such a vertex can get distributed to
       * adjacent dof. This is the maximum number of of such dof.
       *
       * As nonconforming meshes are not implemented yet, this is always 1.
       */
      int get_max_contribs() const
      {
        return 1;
      }

      /**
       * \brief Returns the number of dof contributions.
       *
       * \param[in] local_dof_idx
       * The index of the local dof whose contribution count is to be returned.
       *
       * \returns
       * The number of global dof contributions for the specified local dof.
       */
      int get_num_contribs(int /*local_dof_idx*/) const
      {
        return 1;
      }

#ifdef DOXYGEN
      /**
       * \brief Returns the mapped dof contribution index.
       *
       * \param[in] local_dof_idx
       * The index of the local dof whose mapped contribution index is to be returned.
       *
       * \param[in] contrib_idx
       * The contribution index for the local dof.
       *
       * \returns
       * The mapped dof contribution index.
       */
      Index get_index(int local_dof_idx, int contrib_idx = 0) const;
#endif // DOXYGEN

      /**
       * \brief Returns the mapped dof contribution weight.
       *
       * \param[in] local_dof_idx
       * The index of the local dof whose mapped contribution index is to be returned.
       *
       * \param[in] contrib_idx
       * The contribution index for the local dof.
       *
       * \returns
       * The mapped dof contribution index.
       */
      Real get_weight(int /*local_dof_idx*/, int /*contrib_idx*/ = 0) const
      {
        return Real(1.0);
      }
    }; // class DofMappingBase<...>
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_DOF_MAPPING_BASE_HPP
