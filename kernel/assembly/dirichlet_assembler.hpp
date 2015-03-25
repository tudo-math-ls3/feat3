#pragma once
#ifndef KERNEL_ASSEMBLY_DIRICHLET_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_DIRICHLET_ASSEMBLER_HPP 1

// includes, FEAST
#include <kernel/assembly/base.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/dense_vector.hpp>

// includes, system
#include <set>
#include <map>
#include <vector>

namespace FEAST
{
  namespace Assembly
  {
    /// \cond internal
    namespace Intern
    {
      // forward declarations
      template<int shape_dim_>
      struct DirichletWrapper;
    } // namespace Intern
    /// \endcond

    /**
     * \brief Dirichlet boundary condition assembly class template.
     *
     * \tparam Space_
     * The type of the trial space for which the boundary conditions are to be assembled.
     *
     * \author Peter Zajac
     */
    template<typename Space_>
    class DirichletAssembler
    {
    public:
      /// space type
      typedef Space_ SpaceType;

    private:
      /// shape dimension
      static constexpr int shape_dim = SpaceType::shape_dim;

      /// dof-index set typedef
      typedef std::set<Index> IdxSet;

      /// space reference
      const SpaceType& _space;

      /// dof-index set
      IdxSet _cells[shape_dim + 1];

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] space
       * A reference to the space for which the boundary conditions are to be assembled.
       */
      explicit DirichletAssembler(const SpaceType& space) :
        _space(space)
      {
      }

      /// virtual destructor
      virtual ~DirichletAssembler()
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
        Intern::DirichletWrapper<shape_dim>::merge(_cells, cell_set);
      }

      /**
       * \brief Builds a homogeneous unit-filter.
       *
       * This function assembles a unit-filter implementing the homogeneous Dirichlet boundary conditions.
       *
       * \param[in,out] filter
       * A reference to the unit-filter where the entries are to be added.
       */
      template<typename MemType_, typename DataType_, typename IndexType_>
      void assemble(LAFEM::UnitFilter<MemType_, DataType_, IndexType_>& filter) const
      {
        // build index set
        std::set<Index> idx_set;
        Intern::DirichletWrapper<shape_dim>::assemble(idx_set, _space, _cells);

        // loop over all dof-indices
        typename std::set<Index>::const_iterator it(idx_set.begin());
        typename std::set<Index>::const_iterator jt(idx_set.end());
        for(Index i(0); it != jt; ++it, ++i)
        {
          // store the dof-index
          filter.add(IndexType_(*it), DataType_(0));
        }
      }

      /**
       * \brief Builds an inhomogeneous unit-filter.
       *
       * This function assembles a unit-filter implementing the inhomogeneous Dirichlet boundary conditions given by
       * the vector vector_. In general, this vector will be a representation of an FE function.
       *
       * \param[in,out] filter
       * A reference to the unit-filter where the entries are to be added.
       *
       * \param[in] vector_
       * A LAFEM::DenseVector containing (among other stuff) the entries to be added.
       *
       */
      template<typename MemType_, typename DataType_, typename IndexType_>
      void assemble(LAFEM::UnitFilter<MemType_, DataType_, IndexType_>& filter,
      const LAFEM::DenseVector<MemType_, DataType_, IndexType_>& vector_) const
      {
        // build index set
        std::set<Index> idx_set;
        Intern::DirichletWrapper<shape_dim>::assemble(idx_set, _space, _cells);

        // loop over all dof-indices
        typename std::set<Index>::const_iterator it(idx_set.begin());
        typename std::set<Index>::const_iterator jt(idx_set.end());
        for(Index i(0); it != jt; ++it, ++i)
        {
          // store the dof-index
          filter.add(IndexType_(*it), vector_(IndexType_(*it)));
        }
      }

      /**
       * \brief Builds an (inhomogeneous) unit-filter.
       *
       * This function assembles a unit-filter implementing Dirichlet boundary conditions based on
       * a functor.
       *
       * \param[in,out] filter
       * A reference to the unit-filter where the entries are to be added.
       *
       * \param[in] functor
       * An object implementing the Analytic::Functor interface representing the boundary value
       * function.
       */
      template<
        typename MemType_,
        typename DataType_,
        typename IndexType_,
        typename Functor_>
      void assemble(LAFEM::UnitFilter<MemType_, DataType_, IndexType_>& filter, const Functor_& functor) const
      {
        // build index-value map
        std::map<Index, DataType_> idx_map;
        Intern::DirichletWrapper<shape_dim>::assemble(idx_map, _space, _cells, functor);

        // loop over all dof-indices
        typename std::map<Index, DataType_>::const_iterator it(idx_map.begin());
        typename std::map<Index, DataType_>::const_iterator jt(idx_map.end());
        for(Index i(0); it != jt; ++it, ++i)
        {
          // store the dof-index
          filter.add(IndexType_(it->first), it->second);
        }
      }
      /**
       * \brief Builds a homogeneous blocked unit filter.
       *
       * This function assembles a blocked unit filter implementing the homogeneous Dirichlet boundary conditions.
       *
       * \param[in,out] filter
       * A reference to the blocked unit filter where the entries are to be added.
       */
      template<typename MemType_, typename DataType_, typename IndexType_, Index BlockSize_>
      void assemble(LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, BlockSize_>& filter) const
      {
        // build index set
        std::set<Index> idx_set;
        Intern::DirichletWrapper<shape_dim>::assemble(idx_set, _space, _cells);

        // loop over all dof-indices
        typename std::set<Index>::const_iterator it(idx_set.begin());
        typename std::set<Index>::const_iterator jt(idx_set.end());
        typename LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, BlockSize_>::ValueType tmp(DataType_(0));
        for(Index i(0); it != jt; ++it, ++i)
        {
          // store the dof-index
          filter.add(IndexType_(*it), tmp);
        }
      }

      /**
       * \brief Builds an inhomogeneous blocked unit filter.
       *
       * This function assembles a blocked unit filter implementing the inhomogeneous Dirichlet boundary conditions
       * given by the vector vector_. In general, this vector will be a representation of an FE function.
       *
       * \param[in,out] filter
       * A reference to the unit-filter where the entries are to be added.
       *
       * \param[in] vector_
       * A LAFEM::DenseVector containing (among other stuff) the entries to be added.
       *
       */
      template<typename MemType_, typename DataType_, typename IndexType_, Index BlockSize_>
      void assemble(LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, BlockSize_>& filter,
      const LAFEM::DenseVectorBlocked<MemType_, DataType_, IndexType_, BlockSize_>& vector_) const
      {
        // build index set
        std::set<Index> idx_set;
        Intern::DirichletWrapper<shape_dim>::assemble(idx_set, _space, _cells);

        // loop over all dof-indices
        typename std::set<Index>::const_iterator it(idx_set.begin());
        typename std::set<Index>::const_iterator jt(idx_set.end());
        for(Index i(0); it != jt; ++it, ++i)
        {
          // store the dof-index
          filter.add(IndexType_(*it), vector_(IndexType_(*it)));
        }
      }

    }; // class DirichletAssembler

    /// \cond internal
    namespace Intern
    {
      /// Dirichlet-BC assembly helper class; this is where the magic happens
      template<int shape_dim_>
      struct DirichletHelper
      {
        template<typename CellSet_>
        static void merge(std::set<Index>& idx, const CellSet_& cell_set)
        {
          // fetch the target set for this dimension
          const typename CellSet_::template TargetSet<shape_dim_>::Type&
            target_set(cell_set.template get_target_set<shape_dim_>());

          // merge
          const Index num_entities = target_set.get_num_entities();
          for(Index i(0); i < num_entities; ++i)
          {
            idx.insert(target_set[i]);
          }
        }

        /// homogeneous assembly
        template<typename Space_>
        static void assemble(std::set<Index>& idx, const Space_& space, const std::set<Index>& cells)
        {
          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_>::Type dof_assign(space);

          // loop over all target indices
          std::set<Index>::const_iterator it(cells.begin());
          std::set<Index>::const_iterator jt(cells.end());
          for(; it != jt; ++it)
          {
            dof_assign.prepare(*it);
            const Index num_assign(dof_assign.get_num_assigned_dofs());
            for(Index j(0); j < num_assign; ++j)
            {
              if(dof_assign.get_derive_order(j) > Index(0))
                continue;

              const Index num_contribs(dof_assign.get_num_contribs(j));
              for(Index k(0); k < num_contribs; ++k)
              {
                idx.insert(dof_assign.get_index(j, k));
              }
            }
            dof_assign.finish();
          }
        }

        /// functor-based assembly
        template<
          typename DataType_,
          typename Space_,
          typename Functor_>
        static void assemble(
          std::map<Index, DataType_>& idx,
          const Space_& space,
          const std::set<Index>& cells,
          const Functor_& functor)
        {
          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_, DataType_>::Type dof_assign(space);

          // create a node-functional object
          typename Space_::template NodeFunctional<Functor_, shape_dim_, DataType_>::Type node_func(space, functor);

          // check for empty node functional set
          if(node_func.get_max_assigned_dofs() <= 0)
            return;

          // loop over all target indices
          std::set<Index>::const_iterator it(cells.begin());
          std::set<Index>::const_iterator jt(cells.end());
          for(; it != jt; ++it)
          {
            dof_assign.prepare(*it);
            node_func.prepare(*it);
            const Index num_assign(dof_assign.get_num_assigned_dofs());
            for(Index j(0); j < num_assign; ++j)
            {
              if(dof_assign.get_derive_order(j) > Index(0))
                continue;

              // evaluate node functional
              DataType_ value(node_func(j));
              const Index num_contribs(dof_assign.get_num_contribs(j));
              for(Index k(0); k < num_contribs; ++k)
              {
                Index index(dof_assign.get_index(j, k));
                DataType_ weight(dof_assign.get_weight(j, k));
                idx.insert(std::make_pair(index, (weight*value)));
              }
            }
            node_func.finish();
            dof_assign.finish();
          }
        }
      };

      /// Dirichlet-BC assembly wrapper; calls the helper for all shape dimensions
      template<int shape_dim_>
      struct DirichletWrapper
      {
        template<typename CellSet_>
        static void merge(std::set<Index>* idx, const CellSet_& cell_set)
        {
          DirichletWrapper<shape_dim_ - 1>::merge(idx, cell_set);
          DirichletHelper<shape_dim_>::merge(idx[shape_dim_], cell_set);
        }

        template<typename Space_>
        static void assemble(std::set<Index>& idx, const Space_& space, const std::set<Index>* cells)
        {
          DirichletWrapper<shape_dim_ - 1>::assemble(idx, space, cells);
          DirichletHelper<shape_dim_>::assemble(idx, space, cells[shape_dim_]);
        }

        template<
          typename DataType_,
          typename Space_,
          typename Functor_>
        static void assemble(
          std::map<Index, DataType_>& idx,
          const Space_& space,
          const std::set<Index>* cells,
          const Functor_& functor)
        {
          DirichletWrapper<shape_dim_ - 1>::assemble(idx, space, cells, functor);
          DirichletHelper<shape_dim_>::assemble(idx, space, cells[shape_dim_], functor);
        }
      };

      template<>
      struct DirichletWrapper<0>
      {
        template<typename CellSet_>
        static void merge(std::set<Index>* idx, const CellSet_& cell_set)
        {
          DirichletHelper<0>::merge(idx[0], cell_set);
        }

        template<typename Space_>
        static void assemble(std::set<Index>& idx, const Space_& space, const std::set<Index>* cells)
        {
          DirichletHelper<0>::assemble(idx, space, cells[0]);
        }

        template<
          typename DataType_,
          typename Space_,
          typename Functor_>
        static void assemble(
          std::map<Index, DataType_>& idx,
          const Space_& space,
          const std::set<Index>* cells,
          const Functor_& functor)
        {
          DirichletHelper<0>::assemble(idx, space, cells[0], functor);
        }
      };
    } // namespace Intern
    /// \endcond

  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_DIRICHLET_ASSEMBLER_HPP
