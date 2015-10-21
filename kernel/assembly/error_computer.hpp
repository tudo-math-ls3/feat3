#pragma once
#ifndef KERNEL_ASSEMBLY_ERROR_COMPUTER_HPP
#define KERNEL_ASSEMBLY_ERROR_COMPUTER_HPP 1

// includes, FEAST
#include <kernel/analytic/function.hpp>
#include <kernel/assembly/asm_traits.hpp>

namespace FEAST
{
  namespace Assembly
  {
    /// \cond internal
    namespace Intern
    {
      template<bool h0_, typename AsmTraits_, typename AnaTraits_>
      struct SecHelperH0
      {
        typedef typename AsmTraits_::DataType DataType;

        template<typename FuncEval_, typename LocalVec_>
        static DataType eval(
          FuncEval_&,
          const typename AsmTraits_::TrafoEvalData&,
          const typename AsmTraits_::SpaceEvalData&,
          const LocalVec_&,
          const int
          )
        {
          return DataType(0);
        }
      };

      template<typename AsmTraits_, typename AnaTraits_>
      struct SecHelperH0<true, AsmTraits_, AnaTraits_>
      {
        typedef typename AsmTraits_::DataType DataType;

        template<typename FuncEval_, typename LocalVec_>
        static DataType eval(
          FuncEval_& func_eval,
          const typename AsmTraits_::TrafoEvalData& trafo_data,
          const typename AsmTraits_::SpaceEvalData& space_data,
          const LocalVec_& lvad,
          const int num_loc_dofs
          )
        {
          // evaluate function value
          typename AnaTraits_::ValueType value(DataType(0));
          func_eval.value(value, trafo_data.img_point);

          // test function loop
          for(int i(0); i < num_loc_dofs; ++i)
          {
            // subtract basis function value
            value -= lvad[i] * space_data.phi[i].value;
            // continue with next trial function
          }

          // update result
          return value * value;
        }
      };

      template<bool h1_, typename AsmTraits_, typename AnaTraits_>
      struct SecHelperH1
      {
        typedef typename AsmTraits_::DataType DataType;

        template<typename FuncEval_, typename LocalVec_>
        static DataType eval(
          FuncEval_&,
          const typename AsmTraits_::TrafoEvalData&,
          const typename AsmTraits_::SpaceEvalData&,
          const LocalVec_&,
          const int
          )
        {
          return DataType(0);
        }
      };

      template<typename AsmTraits_, typename AnaTraits_>
      struct SecHelperH1<true, AsmTraits_, AnaTraits_>
      {
        typedef typename AsmTraits_::DataType DataType;

        template<typename FuncEval_, typename LocalVec_>
        static DataType eval(
          FuncEval_& func_eval,
          const typename AsmTraits_::TrafoEvalData& trafo_data,
          const typename AsmTraits_::SpaceEvalData& space_data,
          const LocalVec_& lvad,
          const int num_loc_dofs
          )
        {
          // evaluate function gradient
          typename AnaTraits_::GradientType value(DataType(0));
          func_eval.gradient(value, trafo_data.img_point);

          // test function loop
          for(int i(0); i < num_loc_dofs; ++i)
          {
            // subtract basis function gradient
            value -= lvad[i] * space_data.phi[i].grad;
            // continue with next trial function
          }

          // update result
          return Tiny::dot(value, value);
        }
      };

      template<bool h2_, typename AsmTraits_, typename AnaTraits_>
      struct SecHelperH2
      {
        typedef typename AsmTraits_::DataType DataType;

        template<typename FuncEval_, typename LocalVec_>
        static DataType eval(
          FuncEval_&,
          const typename AsmTraits_::TrafoEvalData&,
          const typename AsmTraits_::SpaceEvalData&,
          const LocalVec_&,
          const int
          )
        {
          return DataType(0);
        }
      };

      template<typename AsmTraits_, typename AnaTraits_>
      struct SecHelperH2<true, AsmTraits_, AnaTraits_>
      {
        typedef typename AsmTraits_::DataType DataType;

        template<typename FuncEval_, typename LocalVec_>
        static DataType eval(
          FuncEval_& func_eval,
          const typename AsmTraits_::TrafoEvalData& trafo_data,
          const typename AsmTraits_::SpaceEvalData& space_data,
          const LocalVec_& lvad,
          const int num_loc_dofs
          )
        {
          // evaluate functor
          typename AnaTraits_::HessianType value(DataType(0));
          func_eval.hessian(value, trafo_data.img_point);

          // test function loop
          for(int i(0); i < num_loc_dofs; ++i)
          {
            // subtract basis function gradient
            value -= lvad[i] * space_data.phi[i].hess;
            // continue with next trial function
          }

          // update result
          return value.norm_hessian_sqr();
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Scalar Error information structure
     *
     * This structure encapsulated the up to three error norms computed by the ScalarErrorComputer
     * class template.
     */
    template<typename DataType_>
    struct ScalarErrorInfo
    {
      /// Specifies whether the H0-error was computed.
      bool have_h0;
      /// Specifies whether the H1-error was computed.
      bool have_h1;
      /// Specifies whether the H2-error was computed.
      bool have_h2;

      /// The computed H0-error.
      DataType_ norm_h0;
      /// The computed H1-error.
      DataType_ norm_h1;
      /// The computed H2-error.
      DataType_ norm_h2;

      /// standard constructor
      ScalarErrorInfo() :
        have_h0(false),
        have_h1(false),
        have_h2(false),
        norm_h0(DataType_(0)),
        norm_h1(DataType_(0)),
        norm_h2(DataType_(0))
      {
      }

      /// conversion constructor
      template<typename DT2_>
      ScalarErrorInfo(const ScalarErrorInfo<DT2_>& other) :
        have_h0(other.have_h0),
        have_h1(other.have_h1),
        have_h2(other.have_h2),
        norm_h0(DataType_(other.norm_h0)),
        norm_h1(DataType_(other.norm_h1)),
        norm_h2(DataType_(other.norm_h2))
      {
      }

      /// conversion assignment operator
      template<typename DT2_>
      ScalarErrorInfo& operator=(const ScalarErrorInfo<DT2_>& other)
      {
        have_h0 = other.have_h0;
        have_h1 = other.have_h1;
        have_h2 = other.have_h2;
        norm_h0 = DataType_(other.norm_h0);
        norm_h1 = DataType_(other.norm_h1);
        norm_h2 = DataType_(other.norm_h2);
      }

      /// prints the info to an output stream
      friend std::ostream& operator<<(std::ostream& os, const ScalarErrorInfo& sei)
      {
        // print errors
        if(sei.have_h0)
          os << "H0-Error: " << stringify_fp_sci(sei.norm_h0) << std::endl;
        if(sei.have_h1)
          os << "H1-Error: " << stringify_fp_sci(sei.norm_h1) << std::endl;
        if(sei.have_h2)
          os << "H2-Error: " << stringify_fp_sci(sei.norm_h2) << std::endl;
        return os;
      }
    }; // struct ScalarErrorInfo

    /**
     * \brief Scalar H0/H1/H2-error computer class template
     *
     * This class template implements the computation of errors between a finite element
     * function and an analytic reference function in the H0-norm as well as the H1- and
     * H2-semi-norms, if desired.
     *
     * \tparam max_norm_
     * Specifies the maximal H^k-semi-norm in which the error is to be measured:
     * - max_norm_ = 0: compute only H0-norm of error
     * - max_norm_ = 1: compute H0-norm and H1-semi-norm of error
     * - max_norm_ = 2: compute H0-norm, H1- and H2-semi-norms of error
     *
     * \attention
     * To compute the H1- and H2-semi-norms of the error, both the analytic reference function
     * as well as the finite element space need to be able to compute gradients and hessians,
     * respectively.
     *
     * \author Peter Zajac
     */
    template<int max_norm_ = 0>
    class ScalarErrorComputer
    {
      static_assert(max_norm_ >= 0, "invalid max_norm_ parameter");

    private:
      /// \cond internal
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
        static constexpr bool need_value = (max_norm_ >= 0);
        static constexpr bool need_grad  = (max_norm_ >= 1);
        static constexpr bool need_hess  = (max_norm_ >= 2);
      };
      /// \endcond

    public:
      /**
       * \brief Computes the H0/H1/H2-error.
       *
       * This function computes the norm of the difference of a discrete Finite-Element function
       * and an analytic function in the H0, H1 and H2 semi-norms, depending on the template argument
       * of this class.
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
       * A ScalarErrorInfo object which contains the H0-, H1- and H2-error (semi-)norms.
       */
      template<
        typename Vector_,
        typename Function_,
        typename Space_,
        typename CubatureFactory_>
      static ScalarErrorInfo<typename Vector_::DataType> compute(
        const Vector_& vector,
        const Function_& function,
        const Space_& space,
        const CubatureFactory_& cubature_factory)
      {
        // make sure the function has everything we need
        static_assert(Function_::can_value, "function values are required for H0-Error");
        static_assert(Function_::can_grad || (max_norm_ < 1), "function gradients are required for H1-Error");
        static_assert(Function_::can_hess || (max_norm_ < 2), "function hessians are required for H2-Error");

        // vector type
        typedef Vector_ VectorType;
        // analytic function type
        typedef Function_ FunctionType;
        // space type
        typedef Space_ SpaceType;
        // assembly traits
        typedef AsmTraits1<typename Vector_::DataType, SpaceType, TrafoConfig, SpaceConfig> AsmTraits;
        // data type
        typedef typename AsmTraits::DataType DataType;

        // create the cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // fetch the trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

        // create a function evaluator
        typename FunctionType::template Evaluator<AnalyticEvalTraits> func_eval(function);

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
        typename VectorType::GatherAxpy gather_axpy(vector);

        // initialise result
        ScalarErrorInfo<DataType> result;
        result.have_h0 = (max_norm_ >= 0);
        result.have_h1 = (max_norm_ >= 1);
        result.have_h2 = (max_norm_ >= 2);

        // H0/H1/H2 helpers
        typedef Intern::SecHelperH0<(max_norm_ >= 0), AsmTraits, AnalyticEvalTraits> SecH0;
        typedef Intern::SecHelperH1<(max_norm_ >= 1), AsmTraits, AnalyticEvalTraits> SecH1;
        typedef Intern::SecHelperH2<(max_norm_ >= 2), AsmTraits, AnalyticEvalTraits> SecH2;

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

          // fetch number of local dofs
          int num_loc_dofs = space_eval.get_num_local_dofs();

          // loop over all quadrature points and integrate
          for(int k(0); k < cubature_rule.get_num_points(); ++k)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(k));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // compute integration weight
            DataType omega = trafo_data.jac_det * cubature_rule.get_weight(k);

            // update results
            result.norm_h0 += omega * SecH0::eval(func_eval, trafo_data, space_data, lvad, num_loc_dofs);
            result.norm_h1 += omega * SecH1::eval(func_eval, trafo_data, space_data, lvad, num_loc_dofs);
            result.norm_h2 += omega * SecH2::eval(func_eval, trafo_data, space_data, lvad, num_loc_dofs);

            // continue with next cubature point
          }

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }

        // take square-roots
        result.norm_h0 = Math::sqrt(result.norm_h0);
        result.norm_h1 = Math::sqrt(result.norm_h1);
        result.norm_h2 = Math::sqrt(result.norm_h2);

        // okay, that's it
        return result;
      }
    }; // class ScalarErrorComputer
  } // namespace Assembly
} // namespace FEAST

#endif // KERNEL_ASSEMBLY_ERROR_COMPUTER_HPP
