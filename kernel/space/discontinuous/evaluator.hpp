// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_DISCONTINUOUS_EVALUATOR_HPP
#define KERNEL_SPACE_DISCONTINUOUS_EVALUATOR_HPP 1

// includes, FEAT
#include <kernel/space/parametric_evaluator.hpp>
#include <kernel/space/discontinuous/variant.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Discontinuous
    {
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_,
        typename VariantTag_,
        typename Shape_ = typename Space_::ShapeType>
      class Evaluator DOXY({});

      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_,
        typename Shape_>
      class Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Variant::StdPolyP<0>, Shape_> :
        public EvaluatorBase<
          Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Variant::StdPolyP<0>, Shape_>,
          TrafoEvaluator_,
          SpaceEvalTraits_>
      {
      public:
        /// space evaluation traits
        typedef SpaceEvalTraits_ SpaceEvalTraits;

        /// trafo evaluation traits
        typedef typename TrafoEvaluator_::EvalTraits TrafoEvalTraits;

        /// basis value coefficient type
        typedef typename SpaceEvalTraits::DataType DataType;

        /// can compute function values
        static constexpr SpaceTags eval_caps = SpaceTags::value;

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] space
         * A reference to the Element using this evaluator.
         */
        explicit Evaluator(const Space_& DOXY(space))
        {
        }

        /// virtual destructor
        virtual ~Evaluator()
        {
        }

        /**
         * \brief Returns the number of local DOFs.
         *
         * \returns
         * The number of local dofs.
         */
        int get_num_local_dofs() const
        {
          return Index(1);
        }

        /** \copydoc Space::EvaluatorBase::eval_values() */
        template<SpaceTags space_cfg_, TrafoTags trafo_cfg_>
        void eval_values(
          EvalData<SpaceEvalTraits, space_cfg_>& data,
          const Trafo::EvalData<TrafoEvalTraits, trafo_cfg_>& DOXY(trafo_data)) const
        {
          data.phi[0].value = DataType(1);
        }
      };

      static constexpr SpaceTags ref_caps_p1 = SpaceTags::value | SpaceTags::grad;

      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_,
        int shape_dim_>
      class Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Variant::StdPolyP<1>, Shape::Simplex<shape_dim_> > :
        public ParametricEvaluator<
          Evaluator<
            Space_,
            TrafoEvaluator_,
            SpaceEvalTraits_,
            Variant::StdPolyP<1>,
            Shape::Simplex<shape_dim_> >,
          TrafoEvaluator_,
          SpaceEvalTraits_,
          ref_caps_p1>
      {
      public:
        /// base-class typedef
        typedef ParametricEvaluator<Evaluator, TrafoEvaluator_, SpaceEvalTraits_, ref_caps_p1> BaseClass;

        /// space type
        typedef Space_ SpaceType;

        /// space evaluation traits
        typedef SpaceEvalTraits_ SpaceEvalTraits;

        /// evaluation policy
        typedef typename SpaceEvalTraits::EvalPolicy EvalPolicy;

        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;

        /// data type
        typedef typename SpaceEvalTraits::DataType DataType;

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] space
         * A reference to the Element using this evaluator.
         */
        explicit Evaluator(const SpaceType& DOXY(space))
        {
        }

        /**
         * \brief Returns the number of local DOFs.
         *
         * \returns
         * The number of local dofs.
         */
        int get_num_local_dofs() const
        {
          return (shape_dim_ + 1);
        }

        /**
         * \brief Evaluates the basis function values on the reference cell.
         *
         * \param[out] data
         * A reference to a basis value vector receiving the result.
         *
         * \param[in] point
         * A reference to the point on the reference cell where to evaluate.
         */
        template<typename EvalData_>
        void eval_ref_values(
          EvalData_& data,
          const DomainPointType& point) const
        {
          data.phi[0].ref_value = DataType(1);
          for(int i(0); i < shape_dim_; ++i)
            data.phi[0].ref_value -= (data.phi[i+1].ref_value = point[i]);
        }

        /**
         * \brief Evaluates the basis function gradients on the reference cell.
         *
         * \param[out] data
         * A reference to a basis gradient vector receiveing the result.
         *
         * \param[in] point
         * A reference to the point on the reference cell where to evaluate.
         */
        template<typename EvalData_>
        void eval_ref_gradients(
          EvalData_& data,
          const DomainPointType& DOXY(point)) const
        {
          for(int i(0); i < shape_dim_; ++i)
          {
            data.phi[0].ref_grad[i] = -DataType(1);
            for(int j(0); j < shape_dim_; ++j)
              data.phi[i+1].ref_grad[j] = DataType(i == j ? 1 : 0);
          }
        }
      }; // class Evaluator<...,Simplex<.>>

      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_,
        int shape_dim_>
      class Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Variant::StdPolyP<1>, Shape::Hypercube<shape_dim_>> :
        public EvaluatorBase<
          Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Variant::StdPolyP<1>, Shape::Hypercube<shape_dim_>>,
          TrafoEvaluator_,
          SpaceEvalTraits_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Evaluator, TrafoEvaluator_, SpaceEvalTraits_> BaseClass;

        /// space type
        typedef Space_ SpaceType;

        /// space evaluation traits
        typedef SpaceEvalTraits_ SpaceEvalTraits;

        /// trafo evaluator type
        typedef TrafoEvaluator_ TrafoEvaluator;

        /// trafo evaluator traits
        typedef typename TrafoEvaluator::EvalTraits TrafoEvalTraits;

        /// trafo type
        typedef typename TrafoEvaluator::TrafoType TrafoType;

        /// mesh type
        typedef typename TrafoType::MeshType MeshType;

        /// data type
        typedef typename SpaceEvalTraits::DataType DataType;

        /// evaluation policy
        typedef typename SpaceEvalTraits::EvalPolicy EvalPolicy;

        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;
        /// image point type
        typedef typename EvalPolicy::ImagePointType ImagePointType;

        /// jacobian matrix type
        typedef typename EvalPolicy::JacobianMatrixType JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef typename EvalPolicy::JacobianInverseType JacobianInverseType;

        static constexpr SpaceTags eval_caps = SpaceTags::value | SpaceTags::grad;

        template<SpaceTags cfg_>
        struct ConfigTraits
        {
          /// evaluation data configuration
          static constexpr SpaceTags config = cfg_;

          /// trafo configuration
          static constexpr TrafoTags trafo_config = TrafoTags::img_point;

          /// evaluation data typedef
          typedef Space::EvalData<SpaceEvalTraits, config> EvalDataType;
        };

      protected:
        /// inverse linearized trafo config
        static constexpr TrafoTags inv_lin_trafo_config = TrafoTags::dom_point | TrafoTags::img_point | TrafoTags::jac_inv;

        /// inverse linearized trafo data
        typedef typename TrafoEvaluator::template ConfigTraits<inv_lin_trafo_config>::EvalDataType InvLinTrafoData;

        /// inverse linearized trafo matrix
        JacobianInverseType _inv_lin_mat;

        // inverse linearized trafo vector
        ImagePointType _inv_lin_vec;

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] space
         * A reference to the Element using this evaluator.
         */
        explicit Evaluator(const SpaceType& DOXY(space))
        {
        }

        /// virtual destructor
        virtual ~Evaluator()
        {
        }

        /**
         * \brief Returns the number of local DOFs.
         *
         * \returns
         * The number of local dofs.
         */
        int get_num_local_dofs() const
        {
          return Index(shape_dim_ + 1);
        }

        /**
         * \brief Prepares the evaluator for a given cell.
         *
         * \param[in] trafo_eval
         * A reference to the trafo evaluator containing the cell information.
         */
        void prepare(const TrafoEvaluator& trafo_eval)
        {
          // create a domain point in the barycentre of the cell
          DomainPointType dom_point(DataType(0));

          // create the trafo data
          InvLinTrafoData trafo_data;
          trafo_eval(trafo_data, dom_point);

          // store inverse trafo linearization
          _inv_lin_mat = trafo_data.jac_inv;
          _inv_lin_vec = trafo_data.img_point;
        }

        /** \copydoc Space::EvaluatorBase::eval_values() */
        template<SpaceTags space_cfg_, TrafoTags trafo_cfg_>
        void eval_values(
          EvalData<SpaceEvalTraits, space_cfg_>& data,
          const Trafo::EvalData<TrafoEvalTraits, trafo_cfg_>& trafo_data) const
        {
          // transform image point
          DomainPointType pt;
          pt.set_mat_vec_mult(_inv_lin_mat, trafo_data.img_point - _inv_lin_vec);

          // evaluate basis functions
          data.phi[0].value = DataType(1);
          for(int i(0); i < shape_dim_; ++i)
            data.phi[i+1].value = pt[i];
        }

        /** \copydoc Space::EvaluatorBase::Eval_gradients */
        template<SpaceTags space_cfg_, TrafoTags trafo_cfg_>
        void eval_gradients(
          EvalData<SpaceEvalTraits, space_cfg_>& data,
          const Trafo::EvalData<TrafoEvalTraits, trafo_cfg_>& trafo_data) const
        {
          // transform image point to local coordinate system
          DomainPointType pt, loc_grad;
          pt.set_mat_vec_mult(_inv_lin_mat, trafo_data.img_point - _inv_lin_vec);

          // loop over all basis functions
          for(int i(0); i < (shape_dim_+1); ++i)
          {
            // compute local gradient
            for(int j(0); j < shape_dim_; ++j)
              loc_grad[j] = DataType(i == (j+1) ? 1 : 0);

            // multiply by transpose inverse linearized trafo matrix for chain rule
            data.phi[i].grad.set_vec_mat_mult(loc_grad, _inv_lin_mat);
          }
        }
      }; // Evaluator<..., Variant::StdNonPar, Shape::Quadrilateral>
    } // namespace Discontinuous
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_DISCONTINUOUS_EVALUATOR_HPP
