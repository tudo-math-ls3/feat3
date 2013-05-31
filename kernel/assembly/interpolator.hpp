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
        typename Functor_,
        typename Space_,
        int shape_dim_>
      class InterpolatorCore
      {
      public:
        static void project(Vector_& vector, const Functor_& functor, const Space_& space)
        {
          typedef typename Vector_::DataType DataType;

          // define node functional
          typedef typename Space_::template NodeFunctional<Functor_, shape_dim_, DataType>::Type NodeFuncType;
          NodeFuncType node_func(space, functor);

          // skip empty node functional sets
          if(node_func.get_max_assigned_dofs() <= 0)
            return;

          // define dof assignment
          typedef typename Space_::template DofAssignment<shape_dim_, DataType>::Type DofAssignType;
          DofAssignType dof_assign(space);

          // loop over all entities
          const Index num_entities = space.get_mesh().get_num_entities(shape_dim_);
          for(Index i(0); i < num_entities; ++i)
          {
            // prepare node functional
            node_func.prepare(i);

            // prepare dof-assignment
            dof_assign.prepare(i);

            // loop over all assigned DOFs
            const Index num_dofs = node_func.get_num_assigned_dofs();
            for(Index j(0); j < num_dofs; ++j)
            {
              // evaluate node functional
              DataType v = node_func(j);

              // loop over all contributions
              const Index num_contribs = dof_assign.get_num_contribs(j);
              for(Index k(0); k < num_contribs; ++k)
              {
                Index idx(dof_assign.get_index(j,k));
                vector(idx, vector(idx) + dof_assign.get_weight(j,k) * v);
              }
            }

            // finish
            dof_assign.finish();
            node_func.finish();
          }
        }
      };

      template<
        typename Vector_,
        typename Functor_,
        typename Space_,
        int shape_dim_ = Space_::shape_dim>
      class InterpolatorWrapper
      {
      public:
        static void project(Vector_& vector, const Functor_& functor, const Space_& space)
        {
          // recurse down
          InterpolatorWrapper<Vector_, Functor_, Space_, shape_dim_ - 1>::project(vector, functor, space);

          // call interpolator core
          InterpolatorCore<Vector_, Functor_, Space_, shape_dim_>::project(vector, functor, space);
        }
      };

      template<
        typename Vector_,
        typename Functor_,
        typename Space_>
      class InterpolatorWrapper<Vector_, Functor_, Space_, 0>
      {
      public:
        static void project(Vector_& vector, const Functor_& functor, const Space_& space)
        {
          // call interpolator core
          InterpolatorCore<Vector_, Functor_, Space_, 0>::project(vector, functor, space);
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
       * \param[in] functor
       * The functor that is to be interpolated. Has to implement the Space::DeriveFunctor interface.
       *
       * \param[in] space
       * The Finite-Element space into which the functor is to be projected.
       */
      template<
        typename Vector_,
        typename Functor_,
        typename Space_>
      static void project(
        Vector_& vector,
        const Functor_& functor,
        const Space_& space)
      {
        vector = Vector_(space.get_num_dofs(), typename Vector_::DataType(0));
        Intern::InterpolatorWrapper<Vector_, Functor_, Space_>::project(vector, functor, space);
      }
    }; // class Interpolator
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_INTERPOLATOR_HPP
