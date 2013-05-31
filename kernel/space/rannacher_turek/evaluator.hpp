#pragma once
#ifndef KERNEL_SPACE_RANNACHER_TUREK_EVALUATOR_HPP
#define KERNEL_SPACE_RANNACHER_TUREK_EVALUATOR_HPP 1

// includes, FEAST
#include <kernel/space/evaluator_base.hpp>
#include <kernel/trafo/eval_data.hpp>
#include <kernel/space/rannacher_turek/variant.hpp>

namespace FEAST
{
  namespace Space
  {
    namespace RannacherTurek
    {
      /**
       * \brief Rannacher-Turek Element Evaluator class template declaration.
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename SpaceEvalTraits_,
        typename TrafoEvaluator_,
        typename VariantTag_,
        typename Shape_ = typename Space_::ShapeType>
      class Evaluator DOXY({});

      /**
       * \brief Implementation of standard non-parametric quadrilateral Rannacher-Turek evaluator
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename SpaceEvalTraits_,
        typename TrafoEvaluator_>
      class Evaluator<Space_, SpaceEvalTraits_, TrafoEvaluator_, Variant::StdNonPar, Shape::Quadrilateral> :
        public EvaluatorBase<
          Evaluator<Space_, SpaceEvalTraits_, TrafoEvaluator_, Variant::StdNonPar, Shape::Quadrilateral>,
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

        /// facet-index-set type
        typedef typename MeshType::template IndexSet<2,1>::Type FacetIndexSetType;

        /// data type
        typedef typename SpaceEvalTraits::DataType DataType;

        /// evaluation policy
        typedef typename SpaceEvalTraits::EvalPolicy EvalPolicy;

        /// domain point coefficient
        typedef typename EvalPolicy::DomainCoordType DomainCoordType;
        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;
        /// image point type
        typedef typename EvalPolicy::ImagePointType ImagePointType;
        /// const image point reference
        typedef typename EvalPolicy::ImagePointConstRef ImagePointConstRef;

        /// jacobian matrix type
        typedef typename EvalPolicy::JacMatType JacMatType;
        /// jacobian inverse matrix type
        typedef typename EvalPolicy::JacInvType JacInvType;

        /// basis value coefficient type
        typedef typename SpaceEvalTraits::BasisValueCoeff BasisValueCoeff;

        /// basis gradient coefficient type
        typedef typename SpaceEvalTraits::BasisGradientCoeff BasisGradientCoeff;

        /** \copydoc EvaluatorBase::EvaluatorCapabilities */
        enum EvaluatorCapabilities
        {
          /// can compute function values
          can_value = 1,

          /// can compute gradients
          can_grad = 1
        };

      protected:
        /// inverse linearised trafo config
        struct InvLinTrafoConfig :
          public Trafo::ConfigBase
        {
          enum
          {
            need_img_point = 1,
            need_jac_inv = 1
          };
        };

        /// inverse linearised trafo data
        typedef Trafo::EvalData<TrafoEvalTraits, InvLinTrafoConfig> InvLinTrafoData;


        /// trafo config for facet trafo
        struct FacetTrafoConfig :
          public Trafo::ConfigBase
        {
          enum
          {
            need_img_point = 1,
            need_jac_det = 1
          };
        };

        /// trafo evaluator for facets (=edges)
        typedef typename TrafoType::template Evaluator<Shape::Hypercube<1>, DataType>::Type FacetTrafoEvaluator;
        /// facet evaluation policy
        typedef typename FacetTrafoEvaluator::EvalPolicy FacetEvalPolicy;
        /// facet evaluation traits
        typedef typename FacetTrafoEvaluator::EvalTraits FacetEvalTraits;
        /// facet trafo data
        typedef Trafo::EvalData<FacetEvalTraits, FacetTrafoConfig> FacetTrafoData;

        /// basis function coefficient matrix
        typedef Tiny::Matrix<DataType, 4, 4> CoeffMatrixType;

        /// inverse linearised trafo matrix
        JacInvType _inv_lin_mat;

        // inverse linearised trafo vector
        ImagePointType _inv_lin_vec;

        /// basis function coefficient matrix
        CoeffMatrixType _coeff_mat;

        void _build_inv_lin_trafo(const TrafoEvaluator& trafo_eval)
        {
          // create a domain point in the barycentre of the cell
          DomainPointType dom_point(DomainCoordType(0));

          // create the trafo data
          InvLinTrafoData trafo_data;
          trafo_data(trafo_eval, dom_point);

          // store inverse trafo linearisation
          _inv_lin_mat = trafo_data.jac_inv;
          _inv_lin_vec = trafo_data.img_point;
        }

        /// computes the basis function coefficient matrix for the current cell
        void _build_coeff_matrix(const TrafoEvaluator& trafo_eval)
        {
          // fetch the trafo
          const TrafoType& trafo = trafo_eval.get_trafo();

          // fetch the mesh
          const MeshType& mesh = trafo.get_mesh();

          // fetch the facet-index-set of the mesh
          const FacetIndexSetType& facet_index_set = mesh.template get_index_set<2,1>();

          // fetch the cell index of the currently active cell
          const Index cell = trafo_eval.get_cell_index();

          // create a nodal matrix
          CoeffMatrixType _nodal_mat(DataType(0));

          // create a facet trafo evaluator
          FacetTrafoEvaluator facet_eval(trafo);

          // create facet evaluation data
          FacetTrafoData facet_data;

          // define 2-point Gauss cubature point coordinate
          static const DataType g = std::sqrt(DataType(1) / DataType(3));
          const typename FacetEvalPolicy::DomainPointType g1(-g), g2(+g);
          DomainPointType q1, q2;

          // loop over all 4 edges of the quad
          for(int i(0); i < 4; ++i)
          {
            // initialise facet evaluator
            facet_eval.prepare(facet_index_set(cell, i));

            // map first cubature point
            facet_data(facet_eval, g1);
            DataType w1(facet_data.jac_det);
            q1.set_mat_vec_mult(_inv_lin_mat, facet_data.img_point - _inv_lin_vec);

            // map second cubature point
            facet_data(facet_eval, g2);
            DataType w2(facet_data.jac_det);
            q2.set_mat_vec_mult(_inv_lin_mat, facet_data.img_point - _inv_lin_vec);

            // release facet evaluator
            facet_eval.finish();

            // compute reciprocal edge length
            DataType v = DataType(1) / (w1 + w2);

            // evaluate node functionals
            _nodal_mat(0,i) = DataType(1); // = v*(w1 + w2)
            _nodal_mat(1,i) = v*(w1*q1(0) + w2*q2(0));
            _nodal_mat(2,i) = v*(w1*q1(1) + w2*q2(1));
            _nodal_mat(3,i) = v*(w1*(q1(0)*q1(0) - q1(1)*q1(1)) + w2*(q2(0)*q2(0) - q2(1)*q2(1)));
          }

          // invert nodal matrix to obtain coefficient matrix
          _coeff_mat.set_inverse(_nodal_mat);
        }

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] space
         * A reference to the Element using this evaluator.
         */
        explicit Evaluator(const SpaceType& /*space*/)
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
        Index get_num_local_dofs() const
        {
          return 4;
        }

        /**
         * \brief Prepares the evaluator for a given cell.
         *
         * \param[in] trafo_eval
         * A reference to the trafo evaluator containing the cell information.
         */
        void prepare(const TrafoEvaluator& trafo_eval)
        {
          // prepare inverse linearised trafo
          _build_inv_lin_trafo(trafo_eval);

          // build coefficient matrix
          _build_coeff_matrix(trafo_eval);
        }

        /** \copydoc Space::EvaluatorBase::eval_values() */
        template<typename SpaceCfg_, typename TrafoCfg_>
        void eval_values(
          EvalData<SpaceEvalTraits, SpaceCfg_>& data,
          const Trafo::EvalData<TrafoEvalTraits, TrafoCfg_>& trafo_data) const
        {
          // transform image point
          DomainPointType pt;
          pt.set_mat_vec_mult(_inv_lin_mat, trafo_data.img_point - _inv_lin_vec);

          // "evaluate" monomials
          DataType x(pt(0));
          DataType y(pt(1));
          DataType r(x*x - y*y);

          // loop over all basis functions
          for(int i(0); i < 4; ++i)
          {
            data.phi[i].value = _coeff_mat(i,0) + _coeff_mat(i,1)*x + _coeff_mat(i,2)*y + _coeff_mat(i,3)*r;
          }
        }

        /** \copydoc Space::EvaluatorBase::Eval_gradients */
        template<typename SpaceCfg_, typename TrafoCfg_>
        void eval_gradients(
          EvalData<SpaceEvalTraits, SpaceCfg_>& data,
          const Trafo::EvalData<TrafoEvalTraits, TrafoCfg_>& trafo_data) const
        {
          // transform image point to local coordinate system
          DomainPointType pt, loc_grad;
          pt.set_mat_vec_mult(_inv_lin_mat, trafo_data.img_point - _inv_lin_vec);

          // loop over all basis functions
          for(int i(0); i < 4; ++i)
          {
            // calculate local gradients
            loc_grad(0) = _coeff_mat(i,1) + DataType(2) * _coeff_mat(i,3) * pt(0);
            loc_grad(1) = _coeff_mat(i,2) - DataType(2) * _coeff_mat(i,3) * pt(1);

            // multiply by transpose inverse linearised trafo matrix for chain rule
            data.phi[i].grad.set_vec_mat_mult(loc_grad, _inv_lin_mat);
          }
        }

      }; // Evaluator<..., Variant::StdNonPar, Shape::Quadrilateral>
    } // namespace RannacherTurek
  } // namespace Space
} // namespace FEAST

#endif // KERNEL_SPACE_RANNACHER_TUREK_EVALUATOR_HPP
