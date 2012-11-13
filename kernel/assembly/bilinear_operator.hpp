#pragma once
#ifndef KERNEL_ASSEMBLY_BILINEAR_OPERATOR_HPP
#define KERNEL_ASSEMBLY_BILINEAR_OPERATOR_HPP 1

// includes, FEAST
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/assembly/local_system_data.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Bilinear Operator Assembly class template
     *
     * This class template implements the assembly of a standard bilinear operator into a matrix.
     *
     * \tparam Matrix_
     * The type of the matrix that is to be assembled.
     *
     * \tparam Functor_
     * The bilinear functor that is to be assembled. \see BilinearFunctorBase
     *
     * \tparam TestSpace_
     * The type of the Finite Element that is to be used as the test space.
     *
     * \tparam TrialSpace_
     * The type of the Finite Element that is to be used as the trial space. If set to Archs::None, the test
     * space is also used as the trial space for the assembly.
     *
     * \author Peter Zajac
     */
    template<
      typename Matrix_,
      typename Functor_,
      typename TestSpace_,
      typename TrialSpace_ = Archs::None>
    class BilinearOperator
    {
    public:
      /// matrix type
      typedef Matrix_ MatrixType;
      /// bilinear functor type
      typedef Functor_ FunctorType;
      /// test-space type
      typedef TestSpace_ TestSpaceType;
      /// trial-space type
      typedef TrialSpace_ TrialSpaceType;

      /// assembly traits
      typedef AsmTraits2<
        typename MatrixType::data_type,
        TestSpace_,
        TrialSpace_,
        typename FunctorType::TrafoConfig,
        typename FunctorType::TestConfig,
        typename FunctorType::TrialConfig> AsmTraits;

      /// data type
      typedef typename AsmTraits::DataType DataType;

      /// cubature rule type
      typedef typename AsmTraits::CubatureRuleType CubatureRuleType;

    public:
      /**
       * \brief Assembles a bilinear operator.
       *
       * \param[in,out] matrix
       * The matrix that is to be assembled.
       *
       * \param[in] functor
       * A reference to the functor defining the bilinearform.
       *
       * \param[in] test_space
       * A reference to the finite-element test-space to be used.
       *
       * \param[in] trial_space
       * A reference to the finite-element trial-space to be used.
       *
       * \param[in] cubature_name
       * A string containing the name of the cubature rule to be used for integration.
       *
       * \param[in] alpha
       * The scaling factor for the bilinear operator.
       */
      static void assemble(
        MatrixType& matrix,
        FunctorType& functor,
        const TestSpaceType& test_space,
        const TrialSpaceType& trial_space,
        const String& cubature_name,
        DataType alpha = DataType(1))
      {
        typedef typename AsmTraits::CubatureDynamicFactoryType CubatureDynamicFactoryType;
        CubatureRuleType cubature_rule(CubatureDynamicFactoryType::create(cubature_name));
        assemble(matrix, functor, test_space, trial_space, cubature_rule, alpha);
      }

      /**
       * \brief Assembles a bilinear operator.
       *
       * \param[in,out] matrix
       * The matrix that is to be assembled.
       *
       * \param[in] functor
       * A reference to the functor defining the bilinearform.
       *
       * \param[in] test_space
       * A reference to the finite-element test-space to be used.
       *
       * \param[in] trial_space
       * A reference to the finite-element trial-space to be used.
       *
       * \param[in] cubature_rule
       * A reference to the cubature rule to be used for integration.
       *
       * \param[in] alpha
       * The scaling factor for the bilinear operator.
       */
      static void assemble(
        MatrixType& matrix,
        FunctorType& functor,
        const TestSpaceType& test_space,
        const TrialSpaceType& trial_space,
        const CubatureRuleType& cubature_rule,
        DataType alpha = DataType(1))
      {
        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = test_space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create space evaluators
        typename AsmTraits::TestSpaceEvaluator test_eval(test_space);
        typename AsmTraits::TrialSpaceEvaluator trial_eval(trial_space);

        // create dof-mappings
        typename AsmTraits::TestDofMapping test_dof_mapping(test_space);
        typename AsmTraits::TrialDofMapping trial_dof_mapping(trial_space);

        // create a functor evaluator
        typename FunctorType::template Evaluator<AsmTraits> func_eval(functor);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::TestSpaceEvalData test_data;
        typename AsmTraits::TrialSpaceEvalData trial_data;

        // create local matrix data
        typename AsmTraits::LocalMatrixDataType lmd(test_dof_mapping, trial_dof_mapping);

        // create matrix scatter-axpy
        MatrixScatterAxpy<MatrixType> scatter_axpy(matrix);

        // loop over all cells of the mesh
        for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluators
          test_eval.prepare(trafo_eval);
          trial_eval.prepare(trafo_eval);

          // prepare functor evaluator
          func_eval.prepare(trafo_eval);

          // fetch number of local dofs
          Index num_loc_test_dofs = test_eval.get_num_local_dofs();
          Index num_loc_trial_dofs = trial_eval.get_num_local_dofs();

          // clear local matrix
          lmd.clear();

          // loop over all quadrature points and integrate
          for(Index k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_data(trafo_eval, cubature_rule.get_point(k));

            // compute basis function data
            test_data(test_eval, trafo_data);
            trial_data(trial_eval, trafo_data);

            // test function loop
            for(Index i(0); i < num_loc_test_dofs; ++i)
            {
              // trial function loop
              for(Index j(0); j < num_loc_trial_dofs; ++j)
              {
                // evaluate functor and integrate
                lmd(i,j) += trafo_data.jac_det * cubature_rule.get_weight(k) *
                  func_eval(trafo_data,
                    typename AsmTraits::TrialFuncData(trial_data, j),
                    typename AsmTraits::TestFuncData(test_data, i));
              }
              // continue with next trial function
            }
            // continue with next test function
          }

          // finish functor evaluator
          func_eval.finish();

          // finish evaluators
          trial_eval.finish();
          test_eval.finish();
          trafo_eval.finish();

          // initialise dof-mappings
          test_dof_mapping.prepare(cell);
          trial_dof_mapping.prepare(cell);

          // incorporate local matrix
          scatter_axpy(lmd, alpha);

          // finish dof-mapping
          trial_dof_mapping.finish();
          test_dof_mapping.finish();

          // continue with next cell
        }

        // okay, that's it
      }
    }; // class BilinearOperator<...>

    /**
     * \brief BilinearOperator specialisation for identical test and trial spaces.
     *
     * \author Peter Zajac
     */
    template<
      typename Matrix_,
      typename Functor_,
      typename Space_>
    class BilinearOperator<Matrix_, Functor_, Space_, Archs::None>
    {
    public:
      /// matrix type
      typedef Matrix_ MatrixType;
      /// bilinear functor type
      typedef Functor_ FunctorType;
      /// space type
      typedef Space_ SpaceType;
      /// test-space type
      typedef Space_ TestSpaceType;
      /// trial-space type
      typedef Space_ TrialSpaceType;

      /// assembly traits
      typedef AsmTraits1<
        typename Matrix_::data_type,
        Space_,
        typename Functor_::TrafoConfig,
        Space::ConfigOr<
          typename Functor_::TestConfig,
          typename Functor_::TrialConfig> > AsmTraits;

      /// data type
      typedef typename AsmTraits::DataType DataType;

      /// cubature rule type
      typedef typename AsmTraits::CubatureRuleType CubatureRuleType;

    public:
      /**
       * \brief Assembles a bilinear operator.
       *
       * \param[in,out] matrix
       * The matrix that is to be assembled.
       *
       * \param[in] functor
       * A reference to the functor defining the bilinearform.
       *
       * \param[in] space
       * A reference to the finite-element to be used as the test- and trial-space.
       *
       * \param[in] cubature_name
       * A string containing the name of the cubature rule to be used for integration.
       *
       * \param[in] alpha
       * The scaling factor for the bilinear operator.
       */
      static void assemble(
        MatrixType& matrix,
        FunctorType& functor,
        const SpaceType& space,
        const String& cubature_name,
        DataType alpha = DataType(1))
      {
        typedef typename AsmTraits::CubatureDynamicFactoryType CubatureDynamicFactoryType;
        CubatureRuleType cubature_rule(CubatureDynamicFactoryType::create(cubature_name));
        assemble(matrix, functor, space, cubature_rule, alpha);
      }

      /**
       * \brief Assembles a bilinear operator.
       *
       * \param[in,out] matrix
       * The matrix that is to be assembled.
       *
       * \param[in] functor
       * A reference to the functor defining the bilinearform.
       *
       * \param[in] space
       * A reference to the finite-element to be used as the test- and trial-space.
       *
       * \param[in] cubature_rule
       * A reference to the cubature rule to be used for integration.
       *
       * \param[in] alpha
       * The scaling factor for the bilinear operator.
       */
      static void assemble(
        MatrixType& matrix,
        FunctorType& functor,
        const SpaceType& space,
        const CubatureRuleType& cubature_rule,
        DataType alpha = DataType(1))
      {
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

        // create local matrix data
        typename AsmTraits::LocalMatrixDataType lmd(dof_mapping, dof_mapping);

        // create matrix scatter-axpy
        MatrixScatterAxpy<MatrixType> scatter_axpy(matrix);

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
          lmd.clear();

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
              // trial function loop
              for(Index j(0); j < num_loc_dofs; ++j)
              {
                // evaluate functor and integrate
                lmd(i,j) += trafo_data.jac_det * cubature_rule.get_weight(k) *
                  func_eval(trafo_data,
                    typename AsmTraits::FuncData(space_data, j),
                    typename AsmTraits::FuncData(space_data, i));
              }
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
          scatter_axpy(lmd, alpha);

          // finish dof-mapping
          dof_mapping.finish();

          // continue with next cell
        }

        // okay, that's it
      }
    }; // class BilinearOperator<..., Archs::None>
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_BILINEAR_OPERATOR_HPP
