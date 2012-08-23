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
      template<typename Evaluator_, bool need_img_point_>
      struct ImgPointData
      {
        void eval(const Evaluator_&, typename Evaluator_::DomainPointConstRef) {}
      };

      template<typename Evaluator_, bool need_jac_mat_>
      struct JacMatData
      {
        void eval(const Evaluator_&, typename Evaluator_::DomainPointConstRef) {}
      };

      template<typename Evaluator_, bool need_jac_inv_>
      struct JacInvData
      {
        void eval(typename Evaluator_::JacMatConstRef) {}
      };

      template<typename Evaluator_, bool need_jac_det_>
      struct JacDetData
      {
        void eval(typename Evaluator_::JacMatConstRef) {}
      };

      template<typename Evaluator_>
      struct ImgPointData<Evaluator_, true>
      {
        static_assert(Evaluator_::can_img_point != 0, "trafo evaluator can't compute image point coordinates");

        /// image point
        typename Evaluator_::ImagePointType img_point;

        void eval(const Evaluator_& evaluator, typename Evaluator_::DomainPointConstRef dom_point)
        {
          // let the evaluator map the point
          evaluator.map_point(img_point, dom_point);
        }
      };

      template<typename Evaluator_>
      struct JacMatData<Evaluator_, true>
      {
        static_assert(Evaluator_::can_jac_mat != 0, "trafo evaluator can't compute jacobian matrices");

        /// jacobian matrix
        typename Evaluator_::JacMatType jac_mat;

        void eval(const Evaluator_& evaluator, typename Evaluator_::DomainPointConstRef dom_point)
        {
          // let the evaluator compute the jacobian matrix
          evaluator.calc_jac_mat(jac_mat, dom_point);
        }
      };

      template<typename Evaluator_>
      struct JacInvData<Evaluator_, true>
      {
        static_assert(Evaluator_::can_jac_inv != 0, "trafo evaluator can't compute jacobian inverse matrices");

        /// jacobian inverse matrix
        typename Evaluator_::JacMatType jac_inv;

        void eval(typename Evaluator_::JacMatConstRef jac_mat)
        {
          // invert the jacobian matrix
          jac_inv.set_inverse(jac_mat);
        }
      };

      template<typename Evaluator_>
      struct JacDetData<Evaluator_, true>
      {
        static_assert(Evaluator_::can_jac_det != 0, "trafo evaluator can't compute jacobian determinants");

        /// jacobian determinant
        typename Evaluator_::JacDetType jac_det;

        void eval(typename Evaluator_::JacMatConstRef jac_mat)
        {
          // compute volume of jacobian matrix
          jac_det = jac_mat.vol();
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
     * \tparam Evaluator_
     * The trafo evaluator that this evaluation data shall use.
     *
     * \tparam Cfg_
     * A trafo config class that specifies what data shall be supplied. See Trafo::ConfigBase for details.
     *
     * \author Peter Zajac
     */
    template<
      typename Evaluator_,
      typename Cfg_>
    class EvalData :
      public Intern::ImgPointData<Evaluator_, Cfg_::need_img_point != 0>,
      public Intern::JacMatData<Evaluator_, Intern::CfgHelper<Cfg_>::need_jac_mat != 0>,
      public Intern::JacInvData<Evaluator_, Cfg_::need_jac_inv != 0>,
      public Intern::JacDetData<Evaluator_, Cfg_::need_jac_det != 0>
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
      typedef Intern::ImgPointData<Evaluator_, have_img_point != 0> ImgPointBase;
      typedef Intern::JacMatData<Evaluator_, have_jac_mat != 0> JacMatBase;
      typedef Intern::JacInvData<Evaluator_, have_jac_inv != 0> JacInvBase;
      typedef Intern::JacDetData<Evaluator_, have_jac_det != 0> JacDetBase;
      /// \endcond

    public:
      /// trafo evaluation policy
      typedef typename Evaluator_::EvalPolicy EvalPolicy;

      /// domain point
      typename EvalPolicy::DomainPointType dom_point;

      /**
       * \brief Evaluation operator
       *
       * \param[in] evaluator_
       * The trafo evaluator that is to be used for evaluation.
       *
       * \param[in] dom_point_
       * A reference to the domain point in which the evaluation shall take place.
       */
      void operator()(const Evaluator_& evaluator, typename EvalPolicy::DomainPointConstRef dom_point_)
      {
        // store domain point
        dom_point = dom_point_;

        // compute other data
        ImgPointBase::eval(evaluator, dom_point);
        JacMatBase::eval(evaluator, dom_point);
        JacInvBase::eval(JacMatBase::jac_mat);
        JacDetBase::eval(JacMatBase::jac_mat);
      }
    }; // class EvalData<...>
  } // namespace Trafo
} // namespace FEAST

#endif // KERNEL_TRAFO_EVAL_DATA_HPP
