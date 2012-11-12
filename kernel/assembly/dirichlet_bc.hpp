#pragma once
#ifndef KERNEL_ASSEMBLY_DIRICHLET_BC_HPP
#define KERNEL_ASSEMBLY_DIRICHLET_BC_HPP 1

// includes, FEAST
#include <kernel/assembly/intern/dof_set.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Homogene Dirichlet boundary condition assembly class template.
     *
     * \tparam Space_
     * The type of the trial space for which the boundary conditions are to be assembled.
     *
     * \author Peter Zajac
     */
    template<typename Space_>
    class HomogeneDirichletBC
    {
    public:
      /// space type
      typedef Space_ SpaceType;

    private:
      /// dof-index set typedef
      typedef std::set<Index> IdxSet;

      /// space reference
      const SpaceType& _space;

      /// dof-index set
      IdxSet _idx_set;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] space
       * A reference to the space for which the boundary conditions are to be assembled.
       */
      explicit HomogeneDirichletBC(const SpaceType& space) :
        _space(space),
        _idx_set()
      {
      }

      /// virtual destructor
      virtual ~HomogeneDirichletBC()
      {
      }

      /**
       * \brief Adds the dofs on a cell-set or submesh to the dof-set.
       *
       * \param[in] cell_set
       * A reference to a cell-set object.
       */
      template<typename CellSet_>
      void add_cell_set(const CellSet_& cell_set)
      {
        Intern::DofSetAddWrapper<Space_, CellSet_>::apply(_idx_set, _space, cell_set);
      }

      /**
       * \brief Builds a unit-filter.
       *
       * This function assembles a unit-filter implementing the homogene Dirichlet boundary conditions.
       *
       * \param[out] filter
       * A reference to the unit-filter to be assembled.
       */
      template<typename MemType_, typename DataType_>
      void build_filter(LAFEM::UnitFilter<MemType_, DataType_>& filter)
      {
        // allocate a filter of the correct size
        filter = LAFEM::UnitFilter<MemType_, DataType_>(Index(_idx_set.size()));

        // fetch the filter's arrays
        Index* indices(filter.get_indices());
        DataType_* values(filter.get_values());

        // loop over all dof-indices
        typename IdxSet::const_iterator it(_idx_set.begin());
        typename IdxSet::const_iterator jt(_idx_set.end());
        for(Index i(0); it != jt; ++it, ++i)
        {
          // store the dof-index and set its value to zero
          indices[i] = *it;
          values[i] = DataType_(0);
        }
      }
    }; // class HomogeneDirichletBC
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_DIRICHLET_BC_HPP
