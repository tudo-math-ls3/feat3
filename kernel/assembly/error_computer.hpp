#pragma once
#ifndef KERNEL_ASSEMBLY_ERROR_COMPUTER_HPP
#define KERNEL_ASSEMBLY_ERROR_COMPUTER_HPP 1

// includes, FEAST
#include <kernel/assembly/asm_traits.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /**
     * \brief Scalar L2-error computer class
     *
     * \author Peter Zajac
     */
    class ScalarErrorComputerL2
    {
    private:
      /**
       * \brief Trafo configuration tag class
       *
       * \see Trafo::ConfigBase
       */
      struct TrafoConfig : public Trafo::ConfigBase
      {
        static constexpr bool need_img_point = true;
        static constexpr bool need_jac_det = true;
      };

      /**
       * \brief Space configuration tag class
       *
       * \see Space::ConfigBase
       */
      struct SpaceConfig : public Space::ConfigBase
      {
        static constexpr bool need_value = true;
      };
      /// \endcond

    public:
      /**
       * \brief Computes the L2-error.
       *
       * This function computes the L2-norm of the difference of a discrete Finite-Element function
       * and an analytic function, i.e. it computes
       *   \f[ \|u - u_h\|_{L^2}. \f]
       *
       * \param[in] vector
       * The coefficient vector of the FE function.
       *
       * \param[in] function
       * An object implementing the AnalyticFunction interface representing the analytic function to be
       * tested against.
       *
       * \param[in] space
       * The Finite-Element space that the coefficient vector belongs to.
       *
       * \param[in] cubature_factory
       * A reference to the cubature factory to be used for integration.
       *
       * \returns
       * The L2-error.
       */
      template<
        typename Vector_,
        typename Function_,
        typename Space_,
        typename CubatureFactory_>
      static typename Vector_::DataType compute(
        const Vector_& vector,
        const Function_& function,
        const Space_& space,
        const CubatureFactory_& cubature_factory)
      {
        // ensure the functor offers function values
        static_assert(Function_::can_value != 0, "analytic function can't compute function values");

        /// vector type
        typedef Vector_ VectorType;
        /// analytic function type
        typedef Function_ FunctionType;
        /// space type
        typedef Space_ SpaceType;
        /// assembly traits
        typedef AsmTraits1<typename Vector_::DataType, SpaceType, TrafoConfig, SpaceConfig> AsmTraits;
        /// data type
        typedef typename AsmTraits::DataType DataType;

        // create the cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a function evaluator
        typename FunctionType::template Evaluator<typename AsmTraits::AnalyticEvalTraits> func_eval(function);

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create local vector data
        typename AsmTraits::LocalVectorType lvad;

        // create matrix scatter-axpy
        LAFEM::GatherAxpy<VectorType> gather_axpy(vector);

        // initialise result
        DataType result(DataType(0));

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // format local vector
          lvad.format();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // incorporate local matrix
          gather_axpy(lvad, dof_mapping);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // prepare functor evaluator
          func_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // evaluate function value
            typename AsmTraits::ValueType value(func_eval.value(trafo_data));

            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // subtract basis function value
              value -= lvad[i] * space_data.phi[i].value;
              // continue with next trial function
            }

            // update result
            result += trafo_data.jac_det * cubature_rule.get_weight(k) * value * value;

            // continue with next cubature point
          }

          // finish functor evaluator
          func_eval.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();


          // continue with next cell
        }

        // okay, that's it
        return Math::sqrt(result);
      }
    }; // class ScalarErrorComputerL2

    /**
     * \brief Scalar H1-error computer class
     *
     * \author Peter Zajac
     */
    class ScalarErrorComputerH1
    {
    private:
      /**
       * \brief Trafo configuration tag class
       *
       * \see Trafo::ConfigBase
       */
      struct TrafoConfig : public Trafo::ConfigBase
      {
        static constexpr bool need_img_point = true;
        static constexpr bool need_jac_det = true;
      };

      /**
       * \brief Space configuration tag class
       *
       * \see Space::ConfigBase
       */
      struct SpaceConfig : public Space::ConfigBase
      {
        static constexpr bool need_grad = true;
      };
      /// \endcond

    public:
      /**
       * \brief Computes the H1-error.
       *
       * This function computes the H1-semi-norm of the difference of a discrete Finite-Element function
       * and an analytic function, i.e. it computes
       *   \f[ |u - u_h|_{H^1}. \f]
       *
       *
       * \param[in] function
       * An object implementing the AnalyticFunction interface representing the analytic function to be
       * tested against.
       *
       * \param[in] space
       * The Finite-Element space that the coefficient vector belongs to.
       *
       * \param[in] cubature_factory
       * A reference to the cubature factory to be used for integration.
       *
       * \returns
       * The H1-error.
       */
      template<
        typename Vector_,
        typename Function_,
        typename Space_,
        typename CubatureFactory_>
      static typename Vector_::DataType compute(
        const Vector_& vector,
        const Function_& function,
        const Space_& space,
        const CubatureFactory_& cubature_factory)
      {
        // ensure the functor offers gradients
        static_assert(Function_::can_grad != 0, "analytic function can't compute gradients");

        /// vector type
        typedef Vector_ VectorType;
        /// analytic function type
        typedef Function_ FunctionType;
        /// space type
        typedef Space_ SpaceType;
        /// assembly traits
        typedef AsmTraits1<typename Vector_::DataType, SpaceType, TrafoConfig, SpaceConfig> AsmTraits;
        /// data type
        typedef typename AsmTraits::DataType DataType;

        // create the cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a function evaluator
        typename FunctionType::template Evaluator<typename AsmTraits::AnalyticEvalTraits> func_eval(function);

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create local vector data
        typename AsmTraits::LocalVectorType lvad;

        // create matrix scatter-axpy
        LAFEM::GatherAxpy<VectorType> gather_axpy(vector);

        // initialise result
        DataType result(DataType(0));

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // format local vector
          lvad.format();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // incorporate local matrix
          gather_axpy(lvad, dof_mapping);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // prepare functor evaluator
          func_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // evaluate function gradient
            typename AsmTraits::GradientType value(func_eval.gradient(trafo_data));

            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // subtract basis function gradient
              value -= lvad[i] * space_data.phi[i].grad;
              // continue with next trial function
            }

            // update result
            result += trafo_data.jac_det * cubature_rule.get_weight(k) * Tiny::dot(value, value);

            // continue with next cubature point
          }

          // finish functor evaluator
          func_eval.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }

        // okay, that's it
        return Math::sqrt(result);
      }
    }; // class ScalarErrorComputerH1

    /**
     * \brief Scalar H2-error computer class
     *
     * \author Peter Zajac
     */
    class ScalarErrorComputerH2
    {
      /**
       * \brief Trafo configuration tag class
       *
       * \see Trafo::ConfigBase
       */
      struct TrafoConfig : public Trafo::ConfigBase
      {
        static constexpr bool need_img_point = true;
        static constexpr bool need_jac_det = true;
      };

      /**
       * \brief Space configuration tag class
       *
       * \see Space::ConfigBase
       */
      struct SpaceConfig : public Space::ConfigBase
      {
        static constexpr bool need_hess = true;
      };
      /// \endcond

    public:
      /**
       * \brief Computes the H2-error.
       *
       * This function computes the H2-semi-norm of the difference of a discrete Finite-Element function
       * and an analytic function, i.e. it computes
       *   \f[ |u - u_h|_{H^2}. \f]
       *
       * \param[in] vector
       * The coefficient vector of the FE function.
       *
       * \param[in] function
       * An object implementing the AnalyticFunction interface representing the analytic function to be
       * tested against.
       *
       * \param[in] space
       * The Finite-Element space that the coefficient vector belongs to.
       *
       * \param[in] cubature_factory
       * A reference to the cubature factory to be used for integration.
       *
       * \returns
       * The H2-error.
       */
      template<
        typename Vector_,
        typename Function_,
        typename Space_,
        typename CubatureFactory_>
      static typename Vector_::DataType compute(
        const Vector_& vector,
        const Function_& function,
        const Space_& space,
        const CubatureFactory_& cubature_factory)
      {
        // ensure the function offers hessians
        static_assert(Function_::can_hess != 0, "analytic function can't compute hessians");

        /// vector type
        typedef Vector_ VectorType;
        /// analytic function type
        typedef Function_ FunctionType;
        /// space type
        typedef Space_ SpaceType;
        /// assembly traits
        typedef AsmTraits1<typename Vector_::DataType, SpaceType, TrafoConfig, SpaceConfig> AsmTraits;
        /// data type
        typedef typename AsmTraits::DataType DataType;

        // create the cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a function evaluator
        typename FunctionType::template Evaluator<typename AsmTraits::AnalyticEvalTraits> func_eval(function);

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create local vector data
        typename AsmTraits::LocalVectorType lvad;

        // create matrix scatter-axpy
        LAFEM::GatherAxpy<VectorType> gather_axpy(vector);

        // initialise result
        DataType result(DataType(0));

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // format local vector
          lvad.format();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // incorporate local matrix
          gather_axpy(lvad, dof_mapping);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // prepare functor evaluator
          func_eval.prepare(trafo_eval);

          // fetch number of local dofs
          int num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // evaluate functor
            typename AsmTraits::HessianType value(func_eval.hessian(trafo_data));

            // test function loop
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // subtract basis function gradient
              value -= lvad[i] * space_data.phi[i].hess;
              // continue with next trial function
            }

            // update result
            result += trafo_data.jac_det * cubature_rule.get_weight(k) * value.hessian_sqr_norm();

            // continue with next cubature point
          }

          // finish functor evaluator
          func_eval.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }

        // okay, that's it
        return Math::sqrt(result);
      }
    }; // class ScalarErrorComputerH2
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_ERROR_COMPUTER_HPP
