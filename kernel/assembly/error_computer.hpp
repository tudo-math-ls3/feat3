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
      /// \cond internal
      template<typename AsmTraits_>
      struct EvalTraits
      {
        typedef typename AsmTraits_::TrafoEvaluator TrafoEvaluator;
        typedef typename AsmTraits_::TrafoData TrafoData;
        typedef typename AsmTraits_::DataType CoeffType;
        typedef typename AsmTraits_::DataType ValueType;
        enum
        {
          domain_dim = AsmTraits_::domain_dim,
          image_dim = AsmTraits_::image_dim
        };
      };

      struct TrafoConfig : public Trafo::ConfigBase
      {
        enum
        {
          need_img_point = 1,
          need_jac_det = 1
        };
      };

      struct SpaceConfig : public Space::ConfigBase
      {
        enum
        {
          need_value = 1
        };
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
       * \param[in] functor
       * An object implementing the Analytic::Functor interface representing the analytic function to be
       * tested against.
       *
       * \param[in] vector
       * The coefficient vector of the FE function.
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
        typename Functor_,
        typename Vector_,
        typename Space_,
        typename CubatureFactory_>
      static typename Vector_::DataType compute(
        const Functor_& functor,
        const Vector_& vector,
        const Space_& space,
        const CubatureFactory_& cubature_factory)
      {
        // ensure the functor offers function values
        static_assert(Functor_::can_value != 0, "analytic functor can't compute function values");

        /// vector type
        typedef Vector_ VectorType;
        /// linear functor type
        typedef Functor_ FunctorType;
        /// space type
        typedef Space_ SpaceType;
        /// assembly traits
        typedef AsmTraits1<typename Vector_::DataType, Space_, TrafoConfig, SpaceConfig> AsmTraits;
        /// data type
        typedef typename AsmTraits::DataType DataType;

        // create the cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a functor value evaluator
        typename Functor_::template ValueEvaluator<EvalTraits<AsmTraits> > func_eval(functor);

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
        typename AsmTraits::LocalVectorDataType lvad(dof_mapping);

        // create matrix scatter-axpy
        LAFEM::GatherAxpy<VectorType> gather_axpy(vector);

        // initialise result
        DataType result(DataType(0));

        // loop over all cells of the mesh
        for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
        {
          // clear local vector
          lvad.clear();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // incorporate local matrix
          gather_axpy(lvad);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // prepare functor evaluator
          func_eval.prepare(trafo_eval);

          // fetch number of local dofs
          Index num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(Index k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_data(trafo_eval, cubature_rule.get_point(k));

            // compute basis function data
            space_data(space_eval, trafo_data);

            // evaluate functor
            DataType value(0);
            func_eval(value, trafo_data);

            // test function loop
            for(Index i(0); i < num_loc_dofs; ++i)
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
      /// \cond internal
      template<typename AsmTraits_>
      struct EvalTraits
      {
        typedef typename AsmTraits_::TrafoEvaluator TrafoEvaluator;
        typedef typename AsmTraits_::TrafoData TrafoData;
        typedef typename AsmTraits_::DataType CoeffType;
        enum
        {
          domain_dim = AsmTraits_::domain_dim,
          image_dim = AsmTraits_::image_dim
        };
        typedef Tiny::Vector<CoeffType, image_dim> ValueType;
      };

      struct TrafoConfig : public Trafo::ConfigBase
      {
        enum
        {
          need_img_point = 1,
          need_jac_det = 1
        };
      };

      struct SpaceConfig : public Space::ConfigBase
      {
        enum
        {
          need_grad = 1
        };
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
       * \param[in] functor
       * An object implementing the Analytic::Functor interface representing the analytic function to be
       * tested against.
       *
       * \param[in] vector
       * The coefficient vector of the FE function.
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
        typename Functor_,
        typename Vector_,
        typename Space_,
        typename CubatureFactory_>
      static typename Vector_::DataType compute(
        const Functor_& functor,
        const Vector_& vector,
        const Space_& space,
        const CubatureFactory_& cubature_factory)
      {
        // ensure the functor offers gradients
        static_assert(Functor_::can_grad != 0, "analytic functor can't compute gradients");

        /// vector type
        typedef Vector_ VectorType;
        /// linear functor type
        typedef Functor_ FunctorType;
        /// space type
        typedef Space_ SpaceType;
        /// assembly traits
        typedef AsmTraits1<typename Vector_::DataType, Space_, TrafoConfig, SpaceConfig> AsmTraits;
        /// data type
        typedef typename AsmTraits::DataType DataType;

        // create the cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a functor gradient evaluator
        typename Functor_::template GradientEvaluator<EvalTraits<AsmTraits> > func_eval(functor);

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
        typename AsmTraits::LocalVectorDataType lvad(dof_mapping);

        // create matrix scatter-axpy
        LAFEM::GatherAxpy<VectorType> gather_axpy(vector);

        // initialise result
        DataType result(DataType(0));

        // loop over all cells of the mesh
        for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
        {
          // clear local vector
          lvad.clear();

          // initialise dof-mapping
          dof_mapping.prepare(cell);

          // incorporate local matrix
          gather_axpy(lvad);

          // finish dof-mapping
          dof_mapping.finish();

          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // prepare functor evaluator
          func_eval.prepare(trafo_eval);

          // fetch number of local dofs
          Index num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(Index k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_data(trafo_eval, cubature_rule.get_point(k));

            // compute basis function data
            space_data(space_eval, trafo_data);

            // evaluate functor
            Tiny::Vector<DataType, AsmTraits::image_dim> value(DataType(0));
            func_eval(value, trafo_data);

            // test function loop
            for(Index i(0); i < num_loc_dofs; ++i)
            {
              // subtract basis function gradient
              value -= lvad[i] * space_data.phi[i].grad;
              // continue with next trial function
            }

            // update result
            result += trafo_data.jac_det * cubature_rule.get_weight(k) * dot(value, value);

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
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_ERROR_COMPUTER_HPP
