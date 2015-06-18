#pragma once
#ifndef KERNEL_SPACE_DOF_MAPPING_BASE_HPP
#define KERNEL_SPACE_DOF_MAPPING_BASE_HPP 1

// includes, FEAST
#include <kernel/space/base.hpp>

namespace FEAST
{
  namespace Space
  {
    /**
     * \brief Finite-Element Dof-Mapping base-class template
     *
     * This class acts as a base-class and interface documentation for Finite-Element Dof-Mapping implementations.
     *
     * \tparam Space_
     * The finite-element space that this dof-mapping is used by.
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
} // namespace FEAST

#endif // KERNEL_SPACE_DOF_MAPPING_BASE_HPP
