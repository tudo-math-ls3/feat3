#pragma once
#ifndef KERNEL_ASSEMBLY_INTERPOLATOR_HPP
#define KERNEL_ASSEMBLY_INTERPOLATOR_HPP 1

// includes, FEAST
#include <kernel/assembly/base.hpp>

namespace FEAST
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
        int shape_dim_>
      class InterpolatorCore
      {
      public:
        static void project(Vector_& vector, const Function_& function, const Space_& space)
        {
          typedef typename Vector_::DataType DataType;

          // create a node-functional object
          typedef typename Space_::template NodeFunctional<shape_dim_, DataType>::Type NodeFunc;

          // check for empty node functional set
          static constexpr int max_dofs = NodeFunc::max_assigned_dofs;
          if(max_dofs <= 0)
            return;

          // create node functional
          NodeFunc node_func(space);

          // create node data; avoid zero-length vectors
          Tiny::Vector<DataType, max_dofs+1> node_data;

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
              // evaluate node functional
              DataType v = node_data[j];

              // loop over all contributions
              const int num_contribs = dof_assign.get_num_contribs(j);
              for(int k(0); k < num_contribs; ++k)
              {
                Index idx(dof_assign.get_index(j,k));
                vector(idx, vector(idx) + dof_assign.get_weight(j,k) * v);
              }
            }

            // finish
            dof_assign.finish();
          }
        }
      };

      template<
        typename Vector_,
        typename Function_,
        typename Space_,
        int shape_dim_ = Space_::shape_dim>
      class InterpolatorWrapper
      {
      public:
        static void project(Vector_& vector, const Function_& function, const Space_& space)
        {
          // recurse down
          InterpolatorWrapper<Vector_, Function_, Space_, shape_dim_ - 1>::project(vector, function, space);

          // call interpolator core
          InterpolatorCore<Vector_, Function_, Space_, shape_dim_>::project(vector, function, space);
        }
      };

      template<
        typename Vector_,
        typename Function_,
        typename Space_>
      class InterpolatorWrapper<Vector_, Function_, Space_, 0>
      {
      public:
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
       * \brief Interpolates a functor.
       *
       * \param[out] vector
       * A vector that shall receive the interpolation coefficients.
       *
       * \param[in] function
       * An object implementing the AnalyticFunction interface.
       *
       * \param[in] space
       * The Finite-Element space into which the functor is to be projected.
       */
      template<
        typename Vector_,
        typename Function_,
        typename Space_>
      static void project(
        Vector_& vector,
        const Function_& function,
        const Space_& space)
      {
        vector = Vector_(space.get_num_dofs(), typename Vector_::DataType(0));
        Intern::InterpolatorWrapper<Vector_, Function_, Space_>::project(vector, function, space);
      }
    }; // class Interpolator
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_INTERPOLATOR_HPP
