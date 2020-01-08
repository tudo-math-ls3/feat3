// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_ERROR_COMPUTER_HPP
#define KERNEL_ASSEMBLY_ERROR_COMPUTER_HPP 1

// includes, FEAT
#include <kernel/analytic/function.hpp>
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/util/dist.hpp>

namespace FEAT
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
          typename AnaTraits_::ValueType value = func_eval.value(trafo_data.img_point);

          // subtract FE function value
          for(int i(0); i < num_loc_dofs; ++i)
            value -= lvad[i] * space_data.phi[i].value;

          // return result
          return Math::sqr(value);
        }
      };

      template<bool h1_, typename AsmTraits_, typename AnaTraits_, bool sub_dim_>
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

      // full-dimension variant for H1
      template<typename AsmTraits_, typename AnaTraits_>
      struct SecHelperH1<true, AsmTraits_, AnaTraits_, false>
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
          typename AnaTraits_::GradientType grad = func_eval.gradient(trafo_data.img_point);

          // subtract FE function gradient
          for(int i(0); i < num_loc_dofs; ++i)
            grad.axpy(-lvad[i], space_data.phi[i].grad);

          // return result
          return grad.norm_euclid_sqr();
        }
      };

      // sub-dimensional variant for H1
      template<typename AsmTraits_, typename AnaTraits_>
      struct SecHelperH1<true, AsmTraits_, AnaTraits_, true>
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
          // get the dimensions
          static constexpr int dom_dim = AsmTraits_::domain_dim;

          // evaluate function gradient
          typename AnaTraits_::GradientType grad = func_eval.gradient(trafo_data.img_point);

          // transform back onto reference element by applying chain rule
          Tiny::Vector<DataType, dom_dim> ref_grad;
          ref_grad.set_vec_mat_mult(grad, trafo_data.jac_mat);

          // subtract FE function reference gradient
          for(int i(0); i < num_loc_dofs; ++i)
            ref_grad.axpy(-lvad[i], space_data.phi[i].ref_grad);

          // compute gram matrix of transformation
          Tiny::Matrix<DataType, dom_dim, dom_dim> gram_mat, gram_inv;
          gram_mat.set_gram(trafo_data.jac_mat);
          gram_inv.set_inverse(gram_mat);

          // return result
          return gram_inv.scalar_product(ref_grad, ref_grad);
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
          // evaluate function hessian
          typename AnaTraits_::HessianType hess = func_eval.hessian(trafo_data.img_point);

          // subtract FE function hessian
          for(int i(0); i < num_loc_dofs; ++i)
            hess.axpy(-lvad[i], space_data.phi[i].hess);

          // return result
          return hess.norm_hessian_sqr();
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

      /**
       * \brief Synchronises the error information over a communicator
       *
       * This function sums up the error information of all patches in a
       * parallel simulation to obtain the information for the global mesh.
       *
       * \param[in] comm
       * The communication over which to synchronise.
       */
      void synchronise(const Dist::Comm& comm)
      {
        DataType_ verr[3] =
        {
          have_h0 ? Math::sqr(norm_h0) : DataType_(0),
          have_h1 ? Math::sqr(norm_h1) : DataType_(0),
          have_h2 ? Math::sqr(norm_h2) : DataType_(0)
        };

        comm.allreduce(verr, verr, std::size_t(3), Dist::op_sum);

        if(have_h0) norm_h0 = Math::sqrt(verr[0]);
        if(have_h1) norm_h1 = Math::sqrt(verr[1]);
        if(have_h2) norm_h2 = Math::sqrt(verr[2]);
      }

      /**
       * \brief Formats the error information as a string.
       *
       * \param[in] precision
       * The precision for floating point values.
       *
       * \param[in] pad_size
       * The leading string padding size. Should be >= 8 to achieve vertical alignment.
       *
       * \param[in] pad_char
       * The leading string padding character.
       *
       * \returns
       * A string containing the formatted error info.
       */
      String format_string(int precision = 0, std::size_t pad_size = 8u, char pad_char = '.') const
      {
        String s;
        if(have_h0)
          s += String("H0-Error").pad_back(pad_size, pad_char) + ": " + stringify_fp_sci(norm_h0, precision) + "\n";
        if(have_h1)
          s += String("H1-Error").pad_back(pad_size, pad_char) + ": " + stringify_fp_sci(norm_h1, precision) + "\n";
        if(have_h2)
          s += String("H2-Error").pad_back(pad_size, pad_char) + ": " + stringify_fp_sci(norm_h2, precision) + "\n";
        return s;
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
     * \tparam sub_dimensional_
     * Specifies whether the errors are to be computed on a sub-dimensional mesh,
     * i.e. on a surface mesh representing a manifold. It is recommended to set this
     * to \c false if working on a fully-dimensional mesh.
     *
     * \attention
     * To compute the H1- and H2-semi-norms of the error, both the analytic reference function
     * as well as the finite element space need to be able to compute gradients and hessians,
     * respectively.
     *
     * \todo implement H2-error on sub-dimensional surface computation
     *
     * \author Peter Zajac
     */
    template<int max_norm_ = 0, bool sub_dimensional_ = false>
    class ScalarErrorComputer
    {
      static_assert(max_norm_ >= 0, "invalid max_norm_ parameter");

    private:
      /// \cond internal
      /**
       * \brief Trafo configuration tag class
       */
      static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_mat | TrafoTags::jac_det;

      /**
       * \brief Space configuration tag class
       */
      static constexpr SpaceTags space_config =
        ((max_norm_ >= 0) ? SpaceTags::value : SpaceTags::none) |
        ((max_norm_ >= 1) ? (sub_dimensional_ ? SpaceTags::ref_grad : SpaceTags::grad) : SpaceTags::none) |
        ((max_norm_ >= 2) ? SpaceTags::hess : SpaceTags::none);
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

        // validate vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");

        // vector type
        typedef Vector_ VectorType;
        // analytic function type
        typedef Function_ FunctionType;
        // space type
        typedef Space_ SpaceType;
        // assembly traits
        typedef AsmTraits1<typename Vector_::DataType, SpaceType, trafo_config, space_config> AsmTraits;
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
        typedef Intern::SecHelperH1<(max_norm_ >= 1), AsmTraits, AnalyticEvalTraits, sub_dimensional_> SecH1;
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

    /// \cond internal
    namespace Intern
    {
      template<bool h0_, typename AsmTraits_, typename AnaTraits_>
      struct VecHelperH0
      {
        typedef typename AsmTraits_::DataType DataType;
        typedef Tiny::Vector<DataType, AnaTraits_::image_dim> ResultType;

        template<typename FuncEval_, typename LocalVec_>
        static ResultType eval(
          FuncEval_&,
          const typename AsmTraits_::TrafoEvalData&,
          const typename AsmTraits_::SpaceEvalData&,
          const LocalVec_&,
          const int
          )
        {
          return ResultType::null();
        }
      };

      template<typename AsmTraits_, typename AnaTraits_>
      struct VecHelperH0<true, AsmTraits_, AnaTraits_>
      {
        typedef typename AsmTraits_::DataType DataType;
        typedef Tiny::Vector<DataType, AnaTraits_::image_dim> ResultType;

        template<typename FuncEval_, typename LocalVec_>
        static ResultType eval(
          FuncEval_& func_eval,
          const typename AsmTraits_::TrafoEvalData& trafo_data,
          const typename AsmTraits_::SpaceEvalData& space_data,
          const LocalVec_& lvad,
          const int num_loc_dofs
          )
        {
          // evaluate function value
          typename AnaTraits_::ValueType value = func_eval.value(trafo_data.img_point);

          // subtract FE function value
          for(int i(0); i < num_loc_dofs; ++i)
            value.axpy(-space_data.phi[i].value, lvad[i]);

          // compute norm of each component
          ResultType result;
          for(int k(0); k < AnaTraits_::image_dim; ++k)
            result[k] = Math::sqr(value[k]);

          return result;
        }
      };

      template<bool h1_, typename AsmTraits_, typename AnaTraits_>
      struct VecHelperH1
      {
        typedef typename AsmTraits_::DataType DataType;
        typedef Tiny::Vector<DataType, AnaTraits_::image_dim> ResultType;

        template<typename FuncEval_, typename LocalVec_>
        static ResultType eval(
          FuncEval_&,
          const typename AsmTraits_::TrafoEvalData&,
          const typename AsmTraits_::SpaceEvalData&,
          const LocalVec_&,
          const int
          )
        {
          return ResultType::null();
        }
      };

      template<typename AsmTraits_, typename AnaTraits_>
      struct VecHelperH1<true, AsmTraits_, AnaTraits_>
      {
        typedef typename AsmTraits_::DataType DataType;
        typedef Tiny::Vector<DataType, AnaTraits_::image_dim> ResultType;

        template<typename FuncEval_, typename LocalVec_>
        static ResultType eval(
          FuncEval_& func_eval,
          const typename AsmTraits_::TrafoEvalData& trafo_data,
          const typename AsmTraits_::SpaceEvalData& space_data,
          const LocalVec_& lvad,
          const int num_loc_dofs
          )
        {
          // evaluate reference function gradient
          typename AnaTraits_::GradientType grad = func_eval.gradient(trafo_data.img_point);

          // subtract FE function gradient
          for(int i(0); i < num_loc_dofs; ++i)
            grad.add_outer_product(lvad[i], space_data.phi[i].grad, -DataType(1));

          // compute norm of each component
          ResultType result;
          for(int k(0); k < AnaTraits_::image_dim; ++k)
            result[k] = grad[k].norm_euclid_sqr();

          return result;
        }
      };

      template<bool h2_, typename AsmTraits_, typename AnaTraits_>
      struct VecHelperH2
      {
        typedef typename AsmTraits_::DataType DataType;
        typedef Tiny::Vector<DataType, AnaTraits_::image_dim> ResultType;

        template<typename FuncEval_, typename LocalVec_>
        static ResultType eval(
          FuncEval_&,
          const typename AsmTraits_::TrafoEvalData&,
          const typename AsmTraits_::SpaceEvalData&,
          const LocalVec_&,
          const int
          )
        {
          return ResultType::null();
        }
      };

      template<typename AsmTraits_, typename AnaTraits_>
      struct VecHelperH2<true, AsmTraits_, AnaTraits_>
      {
        typedef typename AsmTraits_::DataType DataType;
        typedef Tiny::Vector<DataType, AnaTraits_::image_dim> ResultType;

        template<typename FuncEval_, typename LocalVec_>
        static ResultType eval(
          FuncEval_& func_eval,
          const typename AsmTraits_::TrafoEvalData& trafo_data,
          const typename AsmTraits_::SpaceEvalData& space_data,
          const LocalVec_& lvad,
          const int num_loc_dofs
          )
        {
          // evaluate reference function hessian
          typename AnaTraits_::HessianType hess = func_eval.hessian(trafo_data.img_point);

          // subtract FE function hessian
          for(int i(0); i < num_loc_dofs; ++i)
          {
            hess.add_vec_mat_outer_product(lvad[i], space_data.phi[i].hess, -DataType(1));
          }

          // compute norm of each component
          ResultType result;
          for(int k(0); k < AnaTraits_::image_dim; ++k)
            result[k] = hess[k].norm_hessian_sqr();

          return result;
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Vector Error information structure
     *
     * This structure encapsulated the up to three error norms computed by the VectorErrorComputer
     * class template.
     */
    template<typename DataType_, int dim_>
    struct VectorErrorInfo
    {
      /// The dimension of the analysed vector field
      static constexpr int dim = dim_;

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

      /**
       * \brief Vector field components H0-errors
       *
       * This entry contains the H0-errors of all components of the vector field.
       */
      Tiny::Vector<DataType_, dim_> norm_h0_comp;

      /**
       * \brief Vector field components H1-errors
       *
       * This entry contains the H1-errors of all components of the vector field.
       */
      Tiny::Vector<DataType_, dim_> norm_h1_comp;

      /**
       * \brief Vector field components H2-errors
       *
       * This entry contains the H2-errors of all components of the vector field.
       */
      Tiny::Vector<DataType_, dim_> norm_h2_comp;


      /// standard constructor
      VectorErrorInfo() :
        have_h0(false),
        have_h1(false),
        have_h2(false),
        norm_h0(DataType_(0)),
        norm_h1(DataType_(0)),
        norm_h2(DataType_(0)),
        norm_h0_comp(DataType_(0)),
        norm_h1_comp(DataType_(0)),
        norm_h2_comp(DataType_(0))
      {
      }

      /// conversion constructor
      template<typename DT2_>
      VectorErrorInfo(const VectorErrorInfo<DT2_, dim_>& other) :
        have_h0(other.have_h0),
        have_h1(other.have_h1),
        have_h2(other.have_h2),
        norm_h0(DataType_(other.norm_h0)),
        norm_h1(DataType_(other.norm_h1)),
        norm_h2(DataType_(other.norm_h2)),
        norm_h0_comp(other.norm_h0_comp),
        norm_h1_comp(other.norm_h1_comp),
        norm_h2_comp(other.norm_h2_comp)
      {
      }

      /// conversion assignment operator
      template<typename DT2_>
      VectorErrorInfo& operator=(const VectorErrorInfo<DT2_, dim_>& other)
      {
        have_h0 = other.have_h0;
        have_h1 = other.have_h1;
        have_h2 = other.have_h2;
        norm_h0 = DataType_(other.norm_h0);
        norm_h1 = DataType_(other.norm_h1);
        norm_h2 = DataType_(other.norm_h2);
        norm_h0_comp = other.norm_h0_comp;
        norm_h1_comp = other.norm_h1_comp;
        norm_h2_comp = other.norm_h2_comp;
      }

      /**
       * \brief Synchronises the error information over a communicator
       *
       * This function sums up the error information of all patches in a
       * parallel simulation to obtain the information for the global mesh.
       *
       * \param[in] comm
       * The communication over which to synchronise.
       */
      void synchronise(const Dist::Comm& comm)
      {
        DataType_ verr[3+3*dim_] =
        {
          have_h0 ? Math::sqr(norm_h0) : DataType_(0),
          have_h1 ? Math::sqr(norm_h1) : DataType_(0),
          have_h2 ? Math::sqr(norm_h2) : DataType_(0)
        };
        for(int i(0); i < dim_; ++i)
        {
          verr[3+0*dim_+i] = (have_h0 ? Math::sqr(norm_h0_comp[i]) : DataType_(0));
          verr[3+1*dim_+i] = (have_h0 ? Math::sqr(norm_h1_comp[i]) : DataType_(0));
          verr[3+2*dim_+i] = (have_h0 ? Math::sqr(norm_h2_comp[i]) : DataType_(0));
        }

        comm.allreduce(verr, verr, std::size_t(3*dim+3), Dist::op_sum);

        if(have_h0)
        {
          norm_h0 = Math::sqrt(verr[0]);
          for(int i(0); i < dim_; ++i)
            norm_h0_comp[i] = Math::sqrt(verr[3+i]);
        }
        if(have_h1)
        {
          norm_h1 = Math::sqrt(verr[1]);
          for(int i(0); i < dim_; ++i)
            norm_h1_comp[i] = Math::sqrt(verr[3+dim_+i]);
        }
        if(have_h2)
        {
          norm_h2 = Math::sqrt(verr[2]);
          for(int i(0); i < dim_; ++i)
            norm_h2_comp[i] = Math::sqrt(verr[3+2*dim_+i]);
        }
      }

      /**
       * \brief Formats the error information as a string.
       *
       * \param[in] precision
       * The precision for floating point values.
       *
       * \param[in] pad_size
       * The leading string padding size. Should be >= 8 to achieve vertical alignment.
       *
       * \param[in] pad_char
       * The leading string padding character.
       *
       * \returns
       * A string containing the formatted error info.
       */
      String format_string(int precision = 0, std::size_t pad_size = 8u, char pad_char = '.') const
      {
        String s;
        if(have_h0)
        {
          s += String("H0-Error").pad_back(pad_size, pad_char) + ": " + stringify_fp_sci(norm_h0, precision) + " [";
          for(int i(0); i < dim_; ++i)
            (s += " ") += stringify_fp_sci(norm_h0_comp[i], precision);
          s += " ]\n";
        }
        if(have_h1)
        {
          s += String("H1-Error").pad_back(pad_size, pad_char) + ": " + stringify_fp_sci(norm_h1, precision) + " [";
          for(int i(0); i < dim_; ++i)
            (s += " ") += stringify_fp_sci(norm_h1_comp[i], precision);
          s += " ]\n";
        }
        if(have_h2)
        {
          s += String("H2-Error").pad_back(pad_size, pad_char) + ": " + stringify_fp_sci(norm_h2, precision) + " [";
          for(int i(0); i < dim_; ++i)
            (s += " ") += stringify_fp_sci(norm_h2_comp[i], precision);
          s += " ]\n";
        }
        return s;
      }

      /// prints the info to an output stream
      friend std::ostream& operator<<(std::ostream& os, const VectorErrorInfo& sei)
      {
        // print errors
        if(sei.have_h0)
        {
          os << "H0-Error: " << stringify_fp_sci(sei.norm_h0) << " [";
          for(int i(0); i < dim_; ++i)
            os << " " << stringify_fp_sci(sei.norm_h0_comp[i]);
          os << " ]" << std::endl;
        }
        if(sei.have_h1)
        {
          os << "H1-Error: " << stringify_fp_sci(sei.norm_h1) << " [";
          for(int i(0); i < dim_; ++i)
            os << " " << stringify_fp_sci(sei.norm_h1_comp[i]);
          os << " ]" << std::endl;
        }
        if(sei.have_h2)
        {
          os << "H2-Error: " << stringify_fp_sci(sei.norm_h2) << " [";
          for(int i(0); i < dim_; ++i)
            os << " " << stringify_fp_sci(sei.norm_h2_comp[i]);
          os << " ]" << std::endl;
        }
        return os;
      }
    }; // struct VectorErrorInfo

    /**
     * \brief Vector H0/H1/H2-error computer class template
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
    class VectorErrorComputer
    {
      static_assert(max_norm_ >= 0, "invalid max_norm_ parameter");

    private:
      /// \cond internal
      /**
       * \brief Trafo configuration tag class
       */
      static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_det;

      /**
       * \brief Space configuration tag class
       */
      static constexpr SpaceTags space_config =
        ((max_norm_ >= 0) ? SpaceTags::value : SpaceTags::none) |
        ((max_norm_ >= 1) ? SpaceTags::grad : SpaceTags::none) |
        ((max_norm_ >= 2) ? SpaceTags::hess : SpaceTags::none);
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
       * A VectorErrorInfo object which contains the H0-, H1- and H2-error (semi-)norms.
       */
      template<
        typename Vector_,
        typename Function_,
        typename Space_,
        typename CubatureFactory_>
      static VectorErrorInfo<typename Vector_::DataType, Function_::ImageType::scalar_components> compute(
        const Vector_& vector,
        const Function_& function,
        const Space_& space,
        const CubatureFactory_& cubature_factory)
      {
        // make sure the function has everything we need
        static_assert(Function_::can_value, "function values are required for H0-Error");
        static_assert(Function_::can_grad || (max_norm_ < 1), "function gradients are required for H1-Error");
        static_assert(Function_::can_hess || (max_norm_ < 2), "function hessians are required for H2-Error");

        // validate vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");

        // vector type
        typedef Vector_ VectorType;
        // analytic function type
        typedef Function_ FunctionType;
        // space type
        typedef Space_ SpaceType;
        /// The domain dim
        static constexpr int dim = Function_::ImageType::scalar_components;
        /// We use the scalar assembly traits and just define our own LocalVectorData below
        typedef AsmTraits1
        <
          typename Vector_::DataType,
          SpaceType,
          trafo_config,
          space_config
        > AsmTraits;
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
        typedef Tiny::Vector<DataType, dim> VectorValue;
        typedef Tiny::Vector<VectorValue, AsmTraits::max_local_test_dofs> LocalVectorType;
        LocalVectorType lvad;

        // create matrix scatter-axpy
        typename VectorType::GatherAxpy gather_axpy(vector);

        // initialise result
        VectorErrorInfo<DataType, dim> info;
        info.have_h0 = (max_norm_ >= 0);
        info.have_h1 = (max_norm_ >= 1);
        info.have_h2 = (max_norm_ >= 2);

        // H0/H1/H2 helpers
        typedef Intern::VecHelperH0<(max_norm_ >= 0), AsmTraits, AnalyticEvalTraits> VecH0;
        typedef Intern::VecHelperH1<(max_norm_ >= 1), AsmTraits, AnalyticEvalTraits> VecH1;
        typedef Intern::VecHelperH2<(max_norm_ >= 2), AsmTraits, AnalyticEvalTraits> VecH2;

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
            info.norm_h0_comp.axpy(omega, VecH0::eval(func_eval, trafo_data, space_data, lvad, num_loc_dofs));
            info.norm_h1_comp.axpy(omega, VecH1::eval(func_eval, trafo_data, space_data, lvad, num_loc_dofs));
            info.norm_h2_comp.axpy(omega, VecH2::eval(func_eval, trafo_data, space_data, lvad, num_loc_dofs));

            // continue with next cubature point
          }

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();

          // continue with next cell
        }

        // compute the rest
        for(int i(0); i < dim; ++i)
        {
          info.norm_h0 += info.norm_h0_comp[i];
          info.norm_h1 += info.norm_h1_comp[i];
          info.norm_h2 += info.norm_h2_comp[i];
          info.norm_h0_comp[i] = Math::sqrt(info.norm_h0_comp[i]);
          info.norm_h1_comp[i] = Math::sqrt(info.norm_h1_comp[i]);
          info.norm_h2_comp[i] = Math::sqrt(info.norm_h2_comp[i]);
        }

        // take square-roots
        info.norm_h0 = Math::sqrt(info.norm_h0);
        info.norm_h1 = Math::sqrt(info.norm_h1);
        info.norm_h2 = Math::sqrt(info.norm_h2);

        // okay, that's it
        return info;
      }
    }; // class VectorErrorComputer

  } // namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_ERROR_COMPUTER_HPP
