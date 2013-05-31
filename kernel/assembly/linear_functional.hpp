#pragma once
#ifndef KERNEL_ASSEMBLY_LINEAR_FUNCTIONAL_HPP
#define KERNEL_ASSEMBLY_LINEAR_FUNCTIONAL_HPP 1

// includes, FEAST
#include <kernel/assembly/asm_traits.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Linear Functional Assembly class template
     *
     * This class template implements the assembly of a standard linear functional into a vector.
     *
     * \tparam Functor_
     * The linear functor that is to be assembled. \see LinearFunctorBase
     *
     * \author Peter Zajac
     */
    template<typename Functor_>
    class LinearFunctional
    {
    public:
      /// linear functor type
      typedef Functor_ FunctorType;

    public:
      /**
       * \brief Assembles a linear functional into a vector.
       *
       * \param[in,out] vector
       * The vector that is to be assembled.
       *
       * \param[in] functor
       * A reference to the functor defining the linearform.
       *
       * \param[in] cubature_factory
       * A reference to the cubature factory to be used for integration.
       *
       * \param[in] space
       * A reference to the finite-element space to be used.
       */
      template<
        typename Vector_,
        typename CubatureFactory_,
        typename Space_>
      static void assemble_vector(
        Vector_& vector,
        const FunctorType& functor,
        const CubatureFactory_& cubature_factory,
        const Space_& space,
        typename Vector_::DataType alpha = typename Vector_::DataType(1))
      {
        /// vector type
        typedef Vector_ VectorType;

        /// space type
        typedef Space_ SpaceType;

        /// assembly traits
        typedef AsmTraits1<
          typename Vector_::DataType,
          Space_,
          typename Functor_::TrafoConfig,
          typename Functor_::SpaceConfig> AsmTraits;

        /// data type
        typedef typename AsmTraits::DataType DataType;

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create a functor evaluator
        typename FunctorType::template Evaluator<AsmTraits> func_eval(functor);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create local vector data
        typename AsmTraits::LocalVectorDataType lvad(dof_mapping);

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        VectorScatterAxpy<VectorType> scatter_axpy(vector);

        // loop over all cells of the mesh
        for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // prepare functor evaluator
          func_eval.prepare(trafo_eval);

          // fetch number of local dofs
          Index num_loc_dofs = space_eval.get_num_local_dofs();

          // clear local matrix
          lvad.clear();

          // loop over all quadrature points and integrate
          for(Index k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_data(trafo_eval, cubature_rule.get_point(k));

            // compute basis function data
            space_data(space_eval, trafo_data);

            // test function loop
            for(Index i(0); i < num_loc_dofs; ++i)
            {
              // evaluate functor and integrate
              lvad(i) += trafo_data.jac_det * cubature_rule.get_weight(k) *
                func_eval(trafo_data, space_data.phi[i]);
              // continue with next trial function
            }
            // continue with next test function
          }

          // finish functor evaluator
          func_eval.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // incorporate local matrix
          scatter_axpy(lvad, alpha);

          // finish dof-mapping
          dof_mapping.finish();

          // continue with next cell
        }

        // okay, that's it
      }
    }; // class LinearFunctional
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_LINEAR_FUNCTIONAL_HPP
