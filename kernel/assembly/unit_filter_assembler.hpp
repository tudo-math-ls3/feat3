#pragma once
#ifndef KERNEL_ASSEMBLY_UNIT_FILTER_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_UNIT_FILTER_ASSEMBLER_HPP 1

// includes, FEAST
#include <kernel/assembly/base.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/geometry/mesh_part.hpp>

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
      struct UnitAsmWrapper;
    } // namespace Intern
    /// \endcond

    /**
     * \brief Unit-Filter assembly class template.
     *
     * \tparam Mesh_
     * The type of the mesh on which the unit filter is to be assembled.
     *
     * \author Peter Zajac
     */
    template<typename Mesh_>
    class UnitFilterAssembler
    {
    public:
      /// space type
      typedef Mesh_ MeshType;

    private:
      /// shape dimension
      static constexpr int shape_dim = MeshType::shape_dim;

      /// dof-index set typedef
      typedef std::set<Index> IdxSet;

      /// dof-index set
      IdxSet _cells[shape_dim + 1];

    public:
      /**
       * \brief Constructor.
       */
      explicit UnitFilterAssembler()
      {
      }

      /// virtual destructor
      virtual ~UnitFilterAssembler()
      {
      }

      /**
       * \brief Adds the dofs on a mesh-part to the dof-set.
       *
       * \param[in] mesh_part
       * A reference to a mesh part object.
       */
      void add_mesh_part(const Geometry::MeshPart<MeshType>& mesh_part)
      {
        Intern::UnitAsmWrapper<shape_dim>::merge(_cells, mesh_part);
      }

      /**
       * \brief Builds a homogeneous unit-filter.
       *
       * This function assembles a unit-filter implementing the homogeneous Dirichlet boundary conditions.
       *
       * \param[in,out] filter
       * A reference to the unit-filter where the entries are to be added.
       * The filter will be allocated with the correct size if it is empty.
       *
       * \param[in] space
       * The finite-element space for which the filter is to be assembled.
       */
      template<typename MemType_, typename DataType_, typename IndexType_, typename Space_>
      void assemble(LAFEM::UnitFilter<MemType_, DataType_, IndexType_>& filter, const Space_& space) const
      {
        // build index set
        std::set<Index> idx_set;
        Intern::UnitAsmWrapper<shape_dim>::assemble(idx_set, space, _cells);

        // allocate filter if necessary
        if(filter.size() == Index(0))
        {
          filter = LAFEM::UnitFilter<MemType_, DataType_, IndexType_>(space.get_num_dofs());
        }

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
       * \param[in] space
       * The finite-element space for which the filter is to be assembled.
       *
       * \param[in] vector_
       * A LAFEM::DenseVector containing (among other stuff) the entries to be added.
       *
       */
      template<typename MemType_, typename DataType_, typename IndexType_, typename Space_>
      void assemble(
        LAFEM::UnitFilter<MemType_, DataType_, IndexType_>& filter,
        const Space_& space,
        const LAFEM::DenseVector<MemType_, DataType_, IndexType_>& vector_) const
      {
        // build index set
        std::set<Index> idx_set;
        Intern::UnitAsmWrapper<shape_dim>::assemble(idx_set, space, _cells);

        // allocate filter if necessary
        if(filter.size() == Index(0))
        {
          filter = LAFEM::UnitFilter<MemType_, DataType_, IndexType_>(space.get_num_dofs());
        }

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
       * \param[in] space
       * The finite-element space for which the filter is to be assembled.
       *
       * \param[in] function
       * An object implementing the Analytic::Function interface representing the boundary value
       * function.
       */
      template<
        typename MemType_,
        typename DataType_,
        typename IndexType_,
        typename Space_,
        typename Function_>
      void assemble(
        LAFEM::UnitFilter<MemType_, DataType_, IndexType_>& filter,
        const Space_& space,
        const Function_& function) const
      {
        // build index-value map
        std::map<Index, DataType_> idx_map;
        Intern::UnitAsmWrapper<shape_dim>::assemble(idx_map, space, _cells, function);

        // allocate filter if necessary
        if(filter.size() == Index(0))
        {
          filter = LAFEM::UnitFilter<MemType_, DataType_, IndexType_>(space.get_num_dofs());
        }

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
       *
       * \param[in] space
       * The finite-element space for which the filter is to be assembled.
       */
      template<typename MemType_, typename DataType_, typename IndexType_, int BlockSize_, typename Space_>
      void assemble(
        LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, BlockSize_>& filter,
        const Space_& space) const
      {
        // build index set
        std::set<Index> idx_set;
        Intern::UnitAsmWrapper<shape_dim>::assemble(idx_set, space, _cells);

        // allocate filter if necessary
        if (filter.size() == Index(0))
        {
          filter = LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, BlockSize_>(space.get_num_dofs());
        }

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
       * \param[in] space
       * The finite-element space for which the filter is to be assembled.
       */
      template<typename MemType_, typename DataType_, typename IndexType_, int BlockSize_, typename Space_>
      void assemble(
        LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, BlockSize_>& filter,
        const Space_& space,
        const LAFEM::DenseVectorBlocked<MemType_, DataType_, IndexType_, BlockSize_>& vector_) const
      {
        // build index set
        std::set<Index> idx_set;
        Intern::UnitAsmWrapper<shape_dim>::assemble(idx_set, space, _cells);

        // allocate filter if necessary
        if (filter.size() == Index(0))
        {
          filter = LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, BlockSize_>(space.get_num_dofs());
        }

        // loop over all dof-indices
        typename std::set<Index>::const_iterator it(idx_set.begin());
        typename std::set<Index>::const_iterator jt(idx_set.end());
        for(Index i(0); it != jt; ++it, ++i)
        {
          // store the dof-index
          filter.add(IndexType_(*it), vector_(IndexType_(*it)));
        }
      }

    }; // class UnitFilterAssembler

    /// \cond internal
    namespace Intern
    {
      /// Dirichlet-BC assembly helper class; this is where the magic happens
      template<int shape_dim_>
      struct UnitAsmHelper
      {
        template<typename MeshPart_>
        static void merge(std::set<Index>& idx, const MeshPart_& mesh_part)
        {
          // fetch the target set for this dimension
          const typename MeshPart_::template TargetSet<shape_dim_>::Type&
            target_set(mesh_part.template get_target_set<shape_dim_>());

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
            const int num_assign(dof_assign.get_num_assigned_dofs());
            for(int j(0); j < num_assign; ++j)
            {
              const int num_contribs(dof_assign.get_num_contribs(j));
              for(int k(0); k < num_contribs; ++k)
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
          typename Function_>
        static void assemble(
          std::map<Index, DataType_>& idx,
          const Space_& space,
          const std::set<Index>& cells,
          const Function_& function)
        {
          // create a node-functional object
          typedef typename Space_::template NodeFunctional<shape_dim_, DataType_>::Type NodeFunc;

          // check for empty node functional set
          static constexpr Index max_dofs = NodeFunc::max_assigned_dofs;
          if(max_dofs <= 0)
            return;

          // create node functional
          NodeFunc node_func(space);

          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_, DataType_>::Type dof_assign(space);

          // create node data; avoid zero-length vectors
          Tiny::Vector<DataType_, max_dofs+1> node_data;

          // loop over all target indices
          std::set<Index>::const_iterator it(cells.begin());
          std::set<Index>::const_iterator jt(cells.end());
          for(; it != jt; ++it)
          {
            node_func.prepare(*it);
            node_func(node_data, function);
            node_func.finish();

            dof_assign.prepare(*it);

            const int num_assign(dof_assign.get_num_assigned_dofs());
            for(int j(0); j < num_assign; ++j)
            {
              const int num_contribs(dof_assign.get_num_contribs(j));
              for(int k(0); k < num_contribs; ++k)
              {
                Index index(dof_assign.get_index(j, k));
                DataType_ weight(dof_assign.get_weight(j, k));
                idx.insert(std::make_pair(index, (weight*node_data[j])));
              }
            }
            dof_assign.finish();
          }
        }
      };

      /// Dirichlet-BC assembly wrapper; calls the helper for all shape dimensions
      template<int shape_dim_>
      struct UnitAsmWrapper
      {
        template<typename MeshPart_>
        static void merge(std::set<Index>* idx, const MeshPart_& mesh_part)
        {
          UnitAsmWrapper<shape_dim_ - 1>::merge(idx, mesh_part);
          UnitAsmHelper<shape_dim_>::merge(idx[shape_dim_], mesh_part);
        }

        template<typename Space_>
        static void assemble(std::set<Index>& idx, const Space_& space, const std::set<Index>* cells)
        {
          UnitAsmWrapper<shape_dim_ - 1>::assemble(idx, space, cells);
          UnitAsmHelper<shape_dim_>::assemble(idx, space, cells[shape_dim_]);
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
          UnitAsmWrapper<shape_dim_ - 1>::assemble(idx, space, cells, functor);
          UnitAsmHelper<shape_dim_>::assemble(idx, space, cells[shape_dim_], functor);
        }
      };

      template<>
      struct UnitAsmWrapper<0>
      {
        template<typename MeshPart_>
        static void merge(std::set<Index>* idx, const MeshPart_& mesh_part)
        {
          UnitAsmHelper<0>::merge(idx[0], mesh_part);
        }

        template<typename Space_>
        static void assemble(std::set<Index>& idx, const Space_& space, const std::set<Index>* cells)
        {
          UnitAsmHelper<0>::assemble(idx, space, cells[0]);
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
          UnitAsmHelper<0>::assemble(idx, space, cells[0], functor);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_UNIT_FILTER_ASSEMBLER_HPP