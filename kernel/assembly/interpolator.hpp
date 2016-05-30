#pragma once
#ifndef KERNEL_ASSEMBLY_INTERPOLATOR_HPP
#define KERNEL_ASSEMBLY_INTERPOLATOR_HPP 1

// includes, FEAT
#include <kernel/assembly/base.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>

namespace FEAT
{
  namespace Assembly
  {
    /// \cond internal
    namespace Intern
    {
      template<
        typename Vector_,
        typename Function_,
        typename Space_,
        int shape_dim_,
        bool enable_ =
          (Space_::template NodeFunctional<shape_dim_, typename Vector_::DataType>::Type::max_assigned_dofs > 0)>
      struct InterpolatorCore
      {
        static void project(Vector_&, const Function_&, const Space_&)
        {
          // do nothing
        };
      };

      template<
        typename Vector_,
        typename Function_,
        typename Space_,
        int shape_dim_>
      struct InterpolatorCore<Vector_, Function_, Space_, shape_dim_, true>
      {
      public:
        static void project(Vector_& vector, const Function_& function, const Space_& space)
        {
          typedef typename Vector_::DataType DataType;

          // create a node-functional object
          typedef typename Space_::template NodeFunctional<shape_dim_, DataType>::Type NodeFunc;

          // check for empty node functional set
          static constexpr int max_dofs = NodeFunc::max_assigned_dofs;

          typedef typename NodeFunc::template Value<Function_>::Type ValueType;

          // create node functional
          NodeFunc node_func(space);

          // create node data; avoid zero-length vectors
          Tiny::Vector<ValueType, max_dofs> node_data;

          // define dof assignment
          typedef typename Space_::template DofAssignment<shape_dim_, DataType>::Type DofAssignType;
          DofAssignType dof_assign(space);

          // loop over all entities
          const Index num_entities = space.get_mesh().get_num_entities(shape_dim_);
          for(Index i(0); i < num_entities; ++i)
          {
            // evaluate the node functional
            node_func.prepare(i);
            node_func(node_data, function);
            node_func.finish();

            // prepare dof-assignment
            dof_assign.prepare(i);

            // loop over all assigned DOFs
            const int num_dofs = dof_assign.get_num_assigned_dofs();
            for(int j(0); j < num_dofs; ++j)
            {
              // loop over all contributions
              const int num_contribs = dof_assign.get_num_contribs(j);
              for(int k(0); k < num_contribs; ++k)
              {
                Index idx(dof_assign.get_index(j,k));
                vector(idx, vector(idx) + dof_assign.get_weight(j,k) * node_data[j]);
              }
            }

            // finish
            dof_assign.finish();
          }
        }
      };

      template<typename Space_, int shape_dim_ = Space_::shape_dim>
      class InterpolatorWrapper
      {
      public:
        template<typename Vector_, typename Function_>
        static void project(Vector_& vector, const Function_& function, const Space_& space)
        {
          // recurse down
          InterpolatorWrapper<Space_, shape_dim_ - 1>::project(vector, function, space);

          // call interpolator core
          InterpolatorCore<Vector_, Function_, Space_, shape_dim_>::project(vector, function, space);
        }
      };

      template<typename Space_>
      class InterpolatorWrapper<Space_, 0>
      {
      public:
        template<typename Vector_, typename Function_>
        static void project(Vector_& vector, const Function_& function, const Space_& space)
        {
          // call interpolator core
          InterpolatorCore<Vector_, Function_, Space_, 0>::project(vector, function, space);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Interpolator class template
     *
     * This class template implements an interpolator which projects an (analytic) function into a
     * Finite-Element space by evaluating its node functionals.
     *
     * \author Peter Zajac
     */
    class Interpolator
    {
    public:
      /**
       * \brief Interpolates a scalar function.
       *
       * \param[out] vector
       * A vector that shall receive the interpolation coefficients.
       *
       * \param[in] function
       * An object implementing the Analytic::Function interface.
       *
       * \param[in] space
       * The Finite-Element space into which the function is to be projected.
       */
      template<
        typename Mem_,
        typename DT_,
        typename IT_,
        typename Function_,
        typename Space_>
      static void project(
        LAFEM::DenseVector<Mem_, DT_, IT_>& vector,
        const Function_& function,
        const Space_& space)
      {
        typedef typename Space_::TrafoType TrafoType;
        typedef typename Function_::ImageType ImageType;

        // as a very first step, ensure that the function has the right dimensions
        static_assert(Function_::domain_dim == TrafoType::world_dim, "invalid function domain dimension");
        static_assert(ImageType::is_scalar, "only scalar functions can be interpolated to scalar vectors");

        // create vector in Mem::Main
        LAFEM::DenseVector<Mem::Main, DT_, IT_> vec(space.get_num_dofs(), DT_(0));

        // project function
        Intern::InterpolatorWrapper<Space_>::project(vec, function, space);

        // convert to output vector
        vector.convert(vec);
      }

      /**
       * \brief Interpolates a vector field.
       *
       * \param[out] vector
       * A vector that shall receive the interpolation coefficients.
       *
       * \param[in] function
       * An object implementing the Analytic::Function interface.
       *
       * \param[in] space
       * The Finite-Element space into which the function is to be projected.
       */
      template<
        typename Mem_,
        typename DT_,
        typename IT_,
        int block_size_,
        typename Function_,
        typename Space_>
      static void project(
        LAFEM::DenseVectorBlocked<Mem_, DT_, IT_, block_size_>& vector,
        const Function_& function,
        const Space_& space)
      {
        typedef typename Space_::TrafoType TrafoType;
        typedef typename Function_::ImageType ImageType;

        // as a very first step, ensure that the function has the right dimensions
        static_assert(Function_::domain_dim == TrafoType::world_dim, "invalid function domain dimension");
        static_assert(ImageType::is_vector, "only vector fields can be interpolated to blocked vectors");
        static_assert(ImageType::image_dim == block_size_, "invalid vector field size");

        // create vector in Mem::Main
        LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, block_size_> vec(space.get_num_dofs(), DT_(0));

        // project function
        Intern::InterpolatorWrapper<Space_>::project(vec, function, space);

        // convert to output vector
        vector.convert(vec);
      }
    }; // class Interpolator
  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_INTERPOLATOR_HPP
