#pragma once
#ifndef KERNEL_TRAFO_EVAL_DATA_HPP
#define KERNEL_TRAFO_EVAL_DATA_HPP 1

// includes, FEAST
#include <kernel/trafo/base.hpp>

namespace FEAST
{
  namespace Trafo
  {
    /// \cond internal
    namespace Intern
    {
      template<typename EvalTraits_, bool need_img_point_>
      struct ImgPointData
      {
        template<typename Evaluator_>
        void eval(const Evaluator_&, typename EvalTraits_::DomainPointConstRef) {}
      };

      template<typename EvalTraits_, bool need_jac_mat_>
      struct JacMatData
      {
        template<typename Evaluator_>
        void eval(const Evaluator_&, typename EvalTraits_::DomainPointConstRef) {}
      };

      template<typename EvalTraits_, bool need_jac_inv_>
      struct JacInvData
      {
        template<bool need_jac_mat_>
        void eval(JacMatData<EvalTraits_, need_jac_mat_>&) {}
      };

      template<typename EvalTraits_, bool need_jac_det_>
      struct JacDetData
      {
        template<bool need_jac_mat_>
        void eval(JacMatData<EvalTraits_, need_jac_mat_>&) {}
      };

      template<typename EvalTraits_>
      struct ImgPointData<EvalTraits_, true>
      {
        /// image point
        typename EvalTraits_::ImagePointType img_point;

        template<typename Evaluator_>
        void eval(const Evaluator_& evaluator, typename EvalTraits_::DomainPointConstRef dom_point)
        {
          static_assert(Evaluator_::can_img_point != 0, "trafo evaluator can't compute image point coordinates");
          // let the evaluator map the point
          evaluator.map_point(img_point, dom_point);
        }
      };

      template<typename EvalTraits_>
      struct JacMatData<EvalTraits_, true>
      {
        /// jacobian matrix
        typename EvalTraits_::JacMatType jac_mat;

        template<typename Evaluator_>
        void eval(const Evaluator_& evaluator, typename EvalTraits_::DomainPointConstRef dom_point)
        {
          static_assert(Evaluator_::can_jac_mat != 0, "trafo evaluator can't compute jacobian matrices");
          // let the evaluator compute the jacobian matrix
          evaluator.calc_jac_mat(jac_mat, dom_point);
        }
      };

      template<typename EvalTraits_>
      struct JacInvData<EvalTraits_, true>
      {
        /// jacobian inverse matrix
        typename EvalTraits_::JacMatType jac_inv;

        void eval(JacMatData<EvalTraits_, true>& jmd)
        {
          //static_assert(Evaluator_::can_jac_inv != 0, "trafo evaluator can't compute jacobian inverse matrices");
          // invert the jacobian matrix
          jac_inv.set_inverse(jmd.jac_mat);
        }
      };

      template<typename EvalTraits_>
      struct JacDetData<EvalTraits_, true>
      {
        /// jacobian determinant
        typename EvalTraits_::JacDetType jac_det;

        void eval(JacMatData<EvalTraits_, true>& jmd)
        {
          //static_assert(Evaluator_::can_jac_det != 0, "trafo evaluator can't compute jacobian determinants");
          // compute volume of jacobian matrix
          jac_det = jmd.jac_mat.vol();
        }
      };

      /// config helper class template
      template<typename Cfg_>
      struct CfgHelper
      {
        enum
        {
          /// jacobian matrices are needed for jacobian inverse matrices and jacobian determinants
          need_jac_mat = int(Cfg_::need_jac_mat) | int(Cfg_::need_jac_inv) | int(Cfg_::need_jac_det)
        };
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Trafo evaluation data structure.
     *
     * \tparam EvalTraits_
     * The trafo evaluator traits that this evaluation data shall use.
     *
     * \tparam Cfg_
     * A trafo config class that specifies what data shall be supplied. See Trafo::ConfigBase for details.
     *
     * \author Peter Zajac
     */
    template<
      typename EvalTraits_,
      typename Cfg_>
    class EvalData :
      public Intern::ImgPointData<EvalTraits_, Cfg_::need_img_point != 0>,
      public Intern::JacMatData<EvalTraits_, Intern::CfgHelper<Cfg_>::need_jac_mat != 0>,
      public Intern::JacInvData<EvalTraits_, Cfg_::need_jac_inv != 0>,
      public Intern::JacDetData<EvalTraits_, Cfg_::need_jac_det != 0>
    {
    public:
      /// support enumeration
      enum
      {
        /// specifies whether domain point coordinates are given
        have_dom_point = 1,
        /// specifies whether image point coordinates are given
        have_img_point = Cfg_::need_img_point,
        /// specifies whether the jacobian matrix is given
        have_jac_mat = Intern::CfgHelper<Cfg_>::need_jac_mat,
        /// specifies whether the jacobian inverse matrix is given
        have_jac_inv = Cfg_::need_jac_inv,
        /// specifies whether the jacobian determinant is given
        have_jac_det = Cfg_::need_jac_det
      };

      /// \cond internal
      typedef Intern::ImgPointData<EvalTraits_, have_img_point != 0> ImgPointBase;
      typedef Intern::JacMatData<EvalTraits_, have_jac_mat != 0> JacMatBase;
      typedef Intern::JacInvData<EvalTraits_, have_jac_inv != 0> JacInvBase;
      typedef Intern::JacDetData<EvalTraits_, have_jac_det != 0> JacDetBase;
      /// \endcond

    public:
      /// trafo evaluation traits
      typedef EvalTraits_ EvalTraits;

      /// domain point
      typename EvalTraits::DomainPointType dom_point;

      /**
       * \brief Evaluation operator
       *
       * \param[in] evaluator
       * The trafo evaluator that is to be used for evaluation.
       *
       * \param[in] dom_point_
       * A reference to the domain point in which the evaluation shall take place.
       */
      template<typename Evaluator_>
      void operator()(const Evaluator_& evaluator, typename EvalTraits::DomainPointConstRef dom_point_)
      {
        // store domain point
        dom_point = dom_point_;

        // compute other data
        ImgPointBase::eval(evaluator, dom_point);
        JacMatBase::eval(evaluator, dom_point);
        JacInvBase::eval(static_cast<JacMatBase&>(*this));
        JacDetBase::eval(static_cast<JacMatBase&>(*this));
      }
    }; // class EvalData<...>
  } // namespace Trafo
} // namespace FEAST

#endif // KERNEL_TRAFO_EVAL_DATA_HPP
