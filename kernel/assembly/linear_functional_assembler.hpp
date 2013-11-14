#pragma once
#ifndef KERNEL_ASSEMBLY_LINEAR_FUNCTIONAL_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_LINEAR_FUNCTIONAL_ASSEMBLER_HPP 1

// includes, FEAST
#include <kernel/assembly/asm_traits.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Linear Functional Assembly class template
     *
     * This class template implements the assembly of a linear functional into a vector.
     *
     * \author Peter Zajac
     */
    class LinearFunctionalAssembler
    {
    public:
      /**
       * \brief Assembles a linear functional into a vector.
       *
       * \param[in,out] vector
       * The vector that is to be assembled.
       *
       * \param[in] functional
       * A reference to the linear functional implementing the LinearFunctional interface to be assembled.
       *
       * \param[in] space
       * A reference to the finite-element (test) space to be used.
       *
       * \param[in] cubature_factory
       * A reference to the cubature factory to be used for integration.
       */
      template<
        typename Vector_,
        typename Functional_,
        typename CubatureFactory_,
        typename Space_>
      static void assemble_vector(
        Vector_& vector,
        const Functional_& functional,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        typename Vector_::DataType alpha = typename Vector_::DataType(1))
      {
        // vector type
        typedef Vector_ VectorType;
        // functor type
        typedef Functional_ FunctionalType;
        // space type
        typedef Space_ SpaceType;

        // assembly traits
        typedef AsmTraits1<
          typename VectorType::DataType,
          SpaceType,
          typename FunctionalType::TrafoConfig,
          typename FunctionalType::TestConfig> AsmTraits;

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::TestEvaluator test_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create a functor evaluator
        typename FunctionalType::template Evaluator<AsmTraits> func_eval(functional);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData test_data;

        // create local vector data
        typename AsmTraits::LocalVectorDataType lvad(dof_mapping);

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        LAFEM::ScatterAxpy<VectorType> scatter_axpy(vector);

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare test-space evaluator
          test_eval.prepare(trafo_eval);

          // prepare functor evaluator
          func_eval.prepare(trafo_eval);

          // fetch number of local dofs
          Index num_loc_dofs = test_eval.get_num_local_dofs();

          // clear local matrix
          lvad.clear();

          // loop over all quadrature points and integrate
          for(Index k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute test basis function data
            test_eval(test_data, trafo_data);

            // test function loop
            for(Index i(0); i < num_loc_dofs; ++i)
            {
              // evaluate functor and integrate
              lvad(i) += trafo_data.jac_det * cubature_rule.get_weight(k) *
                func_eval(trafo_data, test_data.phi[i]);
              // continue with next trial function
            }
            // continue with next test function
          }

          // finish functor evaluator
          func_eval.finish();

          // finish evaluators
          test_eval.finish();
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
    }; // class LinearFunctionalAssembler
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_LINEAR_FUNCTIONAL_ASSEMBLER_HPP
