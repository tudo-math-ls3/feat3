// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/assembly/base.hpp>
#include <kernel/analytic/function.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_vector_factory.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/geometry/mesh_part.hpp>

// includes, system
#include <set>
#include <map>
#include <vector>

namespace FEAT
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
      UnitFilterAssembler() = default;
      UnitFilterAssembler(UnitFilterAssembler&&) = default;
      UnitFilterAssembler& operator=(UnitFilterAssembler&&) = default;
      UnitFilterAssembler(const UnitFilterAssembler&) = delete;
      UnitFilterAssembler& operator=(const UnitFilterAssembler&) = delete;
      virtual ~UnitFilterAssembler() = default;

      /**
       * \brief Adds the dofs on a mesh-part to the dof-set.
       *
       * \param[in] mesh_part
       * A \transient reference to a mesh part object.
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
       * A \transient reference to the unit-filter where the entries are to be added.
       * The filter will be allocated with the correct size if it is empty.
       *
       * \param[in] space
       * The \transient finite-element space for which the filter is to be assembled.
       */
      template<typename DataType_, typename IndexType_, typename Space_>
      void assemble(LAFEM::UnitFilter<DataType_, IndexType_>& filter, const Space_& space) const
      {
        // create a sparse vector factory of correct size
        LAFEM::SparseVectorFactory<DataType_, IndexType_> svf(space.get_num_dofs());

        // add existing filter entries to factory
        if(!filter.empty())
        {
          XASSERT(filter.size() == svf.size());
          if(!filter.hollow())
            svf.add(filter.get_filter_vector());
        }

        // assemble filter vector
        Intern::UnitAsmWrapper<shape_dim>::assemble(svf, space, _cells);

        // assign filter vector
        filter.get_filter_vector() = svf.make_sv();
      }

      /**
       * \brief Builds an inhomogeneous unit-filter.
       *
       * This function assembles a unit-filter implementing the inhomogeneous Dirichlet boundary conditions given by
       * the vector vector_. In general, this vector will be a representation of an FE function.
       *
       * \param[in,out] filter
       * A \transient reference to the unit-filter where the entries are to be added.
       *
       * \param[in] space
       * The \transient finite-element space for which the filter is to be assembled.
       *
       * \param[in] vector
       * A \transient LAFEM::DenseVector containing (among other stuff) the entries to be added.
       */
      template<typename DTF_, typename ITF_, typename DTV_, typename ITV_, typename Space_>
      void assemble(
        LAFEM::UnitFilter<DTF_, ITF_>& filter,
        const Space_& space,
        const LAFEM::DenseVector<DTV_, ITV_>& vector) const
      {
        // create a sparse vector factory of correct size
        LAFEM::SparseVectorFactory<DTF_, ITF_> svf(space.get_num_dofs());

        // add existing filter entries to factory
        if(!filter.empty())
        {
          XASSERT(filter.size() == svf.size());
          if(!filter.hollow())
            svf.add(filter.get_filter_vector());
        }

        // ensure the vector is of appropriate size
        XASSERT(svf.size() <= vector.size());

        // assemble filter vector by vector
        Intern::UnitAsmWrapper<shape_dim>::assemble(svf, space, _cells, vector);

        // assign filter vector
        filter.get_filter_vector() = svf.make_sv();
      }

      /**
       * \brief Builds an (inhomogeneous) unit-filter.
       *
       * This function assembles a unit-filter implementing Dirichlet boundary conditions based on
       * a function.
       *
       * \param[in,out] filter
       * A \transient reference to the unit-filter where the entries are to be added.
       *
       * \param[in] space
       * The \transient finite-element space for which the filter is to be assembled.
       *
       * \param[in] function
       * An \transient object implementing the Analytic::Function interface representing the scalar boundary value
       * function.
       */
      template<
        typename DataType_,
        typename IndexType_,
        typename Space_,
        typename Function_>
      void assemble(
        LAFEM::UnitFilter<DataType_, IndexType_>& filter,
        const Space_& space,
        const Function_& function) const
      {
        // ensure that the function is scalar
        typedef typename Function_::ImageType FuncImageType;
        static_assert(FuncImageType::is_scalar, "only scalar functions are supported");

        // create a sparse vector factory of correct size
        LAFEM::SparseVectorFactory<DataType_, IndexType_> svf(space.get_num_dofs());

        // add existing filter entries to factory
        if(!filter.empty())
        {
          XASSERT(filter.size() == svf.size());
          if(!filter.hollow())
            svf.add(filter.get_filter_vector());
        }

        // assemble filter vector by function
        Intern::UnitAsmWrapper<shape_dim>::assemble(svf, space, _cells, function);

        // assign filter vector
        filter.get_filter_vector() = svf.make_sv();
      }

      /**
       * \brief Builds a homogeneous blocked unit filter.
       *
       * This function assembles a blocked unit filter implementing the homogeneous Dirichlet boundary conditions.
       *
       * \param[in,out] filter
       * A \transient reference to the blocked unit filter where the entries are to be added.
       *
       * \param[in] space
       * The \transient finite-element space for which the filter is to be assembled.
       */
      template<typename DataType_, typename IndexType_, int block_size_, typename Space_>
      void assemble(
        LAFEM::UnitFilterBlocked<DataType_, IndexType_, block_size_>& filter,
        const Space_& space) const
      {
        // create a sparse vector factory of correct size
        LAFEM::SparseVectorFactory<DataType_, IndexType_, block_size_> svf(space.get_num_dofs());

        // add existing filter entries to factory
        if(!filter.empty())
        {
          XASSERT(filter.size() == svf.size());
          if(!filter.hollow())
            svf.add(filter.get_filter_vector());
        }

        // assemble filter vector
        Intern::UnitAsmWrapper<shape_dim>::assemble(svf, space, _cells);

        // assign filter vector
        filter.get_filter_vector() = svf.make_svb();
      }

      /**
       * \brief Builds an inhomogeneous blocked unit filter.
       *
       * This function assembles a blocked unit filter implementing the inhomogeneous Dirichlet boundary conditions.
       *
       * \param[in,out] filter
       * A \transient reference to the blocked unit filter where the entries are to be added.
       *
       * \param[in] space
       * The \transient finite-element space for which the filter is to be assembled.
       *
       * \param[in] function
       * An \transient object implementing the Analytic::Function interface representing the vector-valued boundary value
       * function.
       */
      template<typename DataType_, typename IndexType_, int block_size_, typename Space_, typename Function_>
      void assemble(
        LAFEM::UnitFilterBlocked<DataType_, IndexType_, block_size_>& filter,
        const Space_& space,
        const Function_& function) const
      {
        // ensure that the function is a vector field of correct dimension
        typedef typename Function_::ImageType FuncImageType;
        static_assert(FuncImageType::is_vector, "only vector fields are supported");
        static_assert(FuncImageType::image_dim == block_size_, "invalid filter block size");

        // create a sparse vector factory of correct size
        LAFEM::SparseVectorFactory<DataType_, IndexType_, block_size_> svf(space.get_num_dofs());

        // add existing filter entries to factory
        if(!filter.empty())
        {
          XASSERT(filter.size() == svf.size());
          if(!filter.hollow())
            svf.add(filter.get_filter_vector());
        }

        // assemble filter vector by function
        Intern::UnitAsmWrapper<shape_dim>::assemble(svf, space, _cells, function);

        // assign filter vector
        filter.get_filter_vector() = svf.make_svb();
      }

      /**
       * \brief Builds an inhomogeneous blocked unit filter.
       *
       * This function assembles a blocked unit filter implementing the inhomogeneous Dirichlet boundary conditions
       * given by the vector vector_. In general, this vector will be a representation of an FE function.
       *
       * \param[in,out] filter
       * A \transient reference to the unit-filter where the entries are to be added.
       *
       * \param[in] vector_
       * A \transient LAFEM::DenseVector containing (among other stuff) the entries to be added.
       *
       * \param[in] space
       * The \transient finite-element space for which the filter is to be assembled.
       */
      template<typename DTF_, typename ITF_, typename DTV_, typename ITV_, int block_size_, typename Space_>
      void assemble(
        LAFEM::UnitFilterBlocked<DTF_, ITF_, block_size_>& filter,
        const Space_& space,
        const LAFEM::DenseVectorBlocked<DTV_, ITV_, block_size_>& vector) const
      {
        // create a sparse vector factory of correct size
        LAFEM::SparseVectorFactory<DTF_, ITF_, block_size_> svf(space.get_num_dofs());

        // add existing filter entries to factory
        if(!filter.empty())
        {
          XASSERT(filter.size() == svf.size());
          if(!filter.hollow())
            svf.add(filter.get_filter_vector());
        }

        // ensure the vector is of appropriate size
        XASSERT(svf.size() <= vector.size());

        // assemble filter vector by function
        Intern::UnitAsmWrapper<shape_dim>::assemble(svf, space, _cells, vector);

        // assign filter vector
        filter.get_filter_vector() = svf.make_svb();
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
        template<typename DT_, typename IT_, int block_size_, typename Space_>
        static void assemble(LAFEM::SparseVectorFactory<DT_, IT_, block_size_>& svf, const Space_& space, const std::set<Index>& cells)
        {
          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_>::Type dof_assign(space);

          Tiny::Vector<DT_, block_size_> null_val;
          null_val.format();

          // loop over all target indices
          std::set<Index>::const_iterator it(cells.begin());
          std::set<Index>::const_iterator jt(cells.end());
          for(; it != jt; ++it)
          {
            dof_assign.prepare(*it);
            const int num_assign(dof_assign.get_num_assigned_dofs());
            for(int j(0); j < num_assign; ++j)
            {
              svf.add(IT_(dof_assign.get_index(j)), null_val);
            }
            dof_assign.finish();
          }
        }

        /// function-based assembly
        template<typename DT_, typename IT_, typename Space_, typename Function_>
        static void assemble(
          LAFEM::SparseVectorFactory<DT_, IT_, 1>& svf,
          const Space_& space,
          const std::set<Index>& cells,
          const Function_& function)
        {
          // create a node-functional object
          typedef typename Space_::template NodeFunctional<shape_dim_, DT_>::Type NodeFunc;

          // check for empty node functional set
          static constexpr Index max_dofs = NodeFunc::max_assigned_dofs;
          if(max_dofs <= 0)
            return;

          // create node functional
          NodeFunc node_func(space);

          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_, DT_>::Type dof_assign(space);

          // create node data; avoid zero-length vectors
          Tiny::Vector<DT_, max_dofs+1> node_data;

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
              svf.add(dof_assign.get_index(j), node_data[j]);
            }
            dof_assign.finish();
          }
        }

        /// function-based assembly
        template<typename DT_, typename IT_, int dim_, typename Space_, typename Function_>
        static void assemble(
          LAFEM::SparseVectorFactory<DT_, IT_, dim_>& svf,
          const Space_& space,
          const std::set<Index>& cells,
          const Function_& function)
        {
          // create a node-functional object
          typedef typename Space_::template NodeFunctional<shape_dim_, DT_>::Type NodeFunc;

          // check for empty node functional set
          static constexpr Index max_dofs = NodeFunc::max_assigned_dofs;
          if(max_dofs <= 0)
            return;

          // create node functional
          NodeFunc node_func(space);

          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_, DT_>::Type dof_assign(space);

          typedef Tiny::Vector<DT_, dim_> ValueType;

          // create node data; avoid zero-length vectors
          Tiny::Vector<ValueType, max_dofs+1> node_data;

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
              svf.add(dof_assign.get_index(j), node_data[j]);
            }
            dof_assign.finish();
          }
        }

        /// vector-based assembly
        template<typename DT_, typename IT_, typename Space_, typename DTV_, typename ITV_>
        static void assemble(
          LAFEM::SparseVectorFactory<DT_, IT_, 1>& svf,
          const Space_& space,
          const std::set<Index>& cells,
          const LAFEM::DenseVector<DTV_, ITV_>& vector)
        {
          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_, DT_>::Type dof_assign(space);

          const Memory::TypedView<DTV_> val_view = vector.elements_view_r();

          // loop over all target indices
          std::set<Index>::const_iterator it(cells.begin());
          std::set<Index>::const_iterator jt(cells.end());
          for(; it != jt; ++it)
          {
            dof_assign.prepare(*it);

            const int num_assign(dof_assign.get_num_assigned_dofs());
            for(int j(0); j < num_assign; ++j)
            {
              const Index idx = dof_assign.get_index(j);
              svf.add(idx, val_view(idx));
            }
            dof_assign.finish();
          }
        }

        /// vector-based assembly
        template<typename DT_, typename IT_, int dim_, typename Space_, typename DTV_, typename ITV_>
        static void assemble(
          LAFEM::SparseVectorFactory<DT_, IT_, dim_>& svf,
          const Space_& space,
          const std::set<Index>& cells,
          const LAFEM::DenseVectorBlocked<DTV_, ITV_, dim_>& vector)
        {
          // create a dof-assignment object
          typename Space_::template DofAssignment<shape_dim_, DT_>::Type dof_assign(space);

          typedef Tiny::Vector<DTV_, dim_> ValueType;
          const Memory::TypedView<ValueType> val_view = vector.elements_view_r();

          // loop over all target indices
          std::set<Index>::const_iterator it(cells.begin());
          std::set<Index>::const_iterator jt(cells.end());
          for(; it != jt; ++it)
          {
            dof_assign.prepare(*it);

            const int num_assign(dof_assign.get_num_assigned_dofs());
            for(int j(0); j < num_assign; ++j)
            {
              const Index idx = dof_assign.get_index(j);
              svf.add(idx, val_view(idx));
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

        template<typename DT_, typename IT_, int dim_, typename Space_>
        static void assemble(LAFEM::SparseVectorFactory<DT_, IT_, dim_>& svf, const Space_& space, const std::set<Index>* cells)
        {
          UnitAsmWrapper<shape_dim_ - 1>::assemble(svf, space, cells);
          UnitAsmHelper<shape_dim_>::assemble(svf, space, cells[shape_dim_]);
        }

        template<typename DT_, typename IT_, int dim_, typename Space_, typename Whatever_>
        static void assemble(
          LAFEM::SparseVectorFactory<DT_, IT_, dim_>& svf,
          const Space_& space,
          const std::set<Index>* cells,
          const Whatever_& whatever)
        {
          UnitAsmWrapper<shape_dim_ - 1>::assemble(svf, space, cells, whatever);
          UnitAsmHelper<shape_dim_>::assemble(svf, space, cells[shape_dim_], whatever);
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

        template<typename DT_, typename IT_, int dim_, typename Space_>
        static void assemble(LAFEM::SparseVectorFactory<DT_, IT_, dim_>& svf, const Space_& space, const std::set<Index>* cells)
        {
          UnitAsmHelper<0>::assemble(svf, space, cells[0]);
        }

        template<typename DT_, typename IT_, int dim_, typename Space_, typename Whatever_>
        static void assemble(
          LAFEM::SparseVectorFactory<DT_, IT_, dim_>& svf,
          const Space_& space,
          const std::set<Index>* cells,
          const Whatever_& whatever)
        {
          UnitAsmHelper<0>::assemble(svf, space, cells[0], whatever);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Assembly
} // namespace FEAT
