// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/assembly/asm_traits.hpp>

namespace FEAT
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
       * A \transient reference to the vector that is to be assembled.
       *
       * \param[in] functional
       * A \transient reference to the linear functional implementing the LinearFunctional interface to be assembled.
       *
       * \param[in] space
       * A \transient reference to the finite-element (test) space to be used.
       *
       * \param[in] cubature_factory
       * A reference to the cubature factory to be used for integration.
       *
       * \param[in] alpha
       * The scaling factor for the linear functional.
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
        // validate vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");

        // vector type
        typedef Vector_ VectorType;
        // functional type
        typedef Functional_ FunctionalType;
        // space type
        typedef Space_ SpaceType;

        // assembly traits
        typedef AsmTraits1<
          typename VectorType::DataType,
          SpaceType,
          FunctionalType::trafo_config,
          FunctionalType::test_config> AsmTraits;

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::TestEvaluator test_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create a functional evaluator
        typename FunctionalType::template Evaluator<AsmTraits> func_eval(functional);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::TestEvalData test_data;

        // the value type of the functional
        typedef typename FunctionalType::template Evaluator<AsmTraits>::ValueType ValueType;

        // ensure that the functional and vector value types are compatible
        static_assert(std::is_same<ValueType, typename VectorType::ValueType>::value,
          "vector and linear functional have different value types!");

        // create local vector data
        typename AsmTraits::template TLocalVector<ValueType> loc_vec;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename VectorType::ScatterAxpy scatter_axpy(vector);

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare test-space evaluator
          test_eval.prepare(trafo_eval);

          // prepare functional evaluator
          func_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_dofs = test_eval.get_num_local_dofs();

          // format local matrix
          loc_vec.format();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute test basis function data
            test_eval(test_data, trafo_data);

            // prepare functional
            func_eval.set_point(trafo_data);

            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // evaluate functional and integrate
              Tiny::axpy(loc_vec(i), func_eval.eval(test_data.phi[i]),
                 trafo_data.jac_det * cubature_rule.get_weight(k));
              // continue with next trial function
            }
            // continue with next test function
          }

          // finish functional evaluator
          func_eval.finish();

          // finish evaluators
          test_eval.finish();
          trafo_eval.finish();

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // incorporate local matrix
          scatter_axpy(loc_vec, dof_mapping, alpha);

          // finish dof-mapping
          dof_mapping.finish();

          // continue with next cell
        }

        // okay, that's it
      }
    }; // class LinearFunctionalAssembler
  } // namespace Assembly
} // namespace FEAT
