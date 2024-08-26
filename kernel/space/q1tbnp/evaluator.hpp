// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_Q1TBNP_EVALUATOR_HPP
#define KERNEL_SPACE_Q1TBNP_EVALUATOR_HPP 1

// includes, FEAT
#include <kernel/space/parametric_evaluator.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Q1TBNP
    {
      /**
       * \brief Q1TBNP Evaluator reference capabilities
       *
       * \author Peter Zajac
       */
      static constexpr SpaceTags ref_caps = SpaceTags::ref_value | SpaceTags::ref_grad;

      /**
       * \brief Q1TBNP Element Evaluator class template declaration.
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_,
        typename Shape_ = typename Space_::ShapeType>
      class Evaluator DOXY({});

      /**
       * \brief Implementation of 2D Q1TBNP evaluator
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_>
      class Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Shape::Quadrilateral> :
        public EvaluatorBase<
          Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Shape::Quadrilateral>,
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
        static constexpr TrafoTags inv_lin_trafo_config = TrafoTags::dom_point | TrafoTags::img_point | TrafoTags::jac_inv | TrafoTags::jac_det;

        /// inverse linearized trafo data
        typedef typename TrafoEvaluator::template ConfigTraits<inv_lin_trafo_config>::EvalDataType InvLinTrafoData;

        /// trafo evaluator for facets (=edges)
        typedef typename TrafoType::template Evaluator<Shape::Hypercube<1>, DataType>::Type FacetTrafoEvaluator;
        /// facet evaluation traits
        typedef typename FacetTrafoEvaluator::EvalTraits FacetEvalTraits;

        /// trafo config for facet trafo
        static constexpr TrafoTags facet_trafo_config = TrafoTags::img_point | TrafoTags::jac_det;

        /// facet trafo data
        typedef typename FacetTrafoEvaluator::template ConfigTraits<facet_trafo_config>::EvalDataType FacetTrafoData;

        /// basis function coefficient matrix
        typedef Tiny::Matrix<DataType, 5, 5> CoeffMatrixType;

        /// inverse linearized trafo matrix
        JacobianInverseType _inv_lin_mat;

        // inverse linearized trafo vector
        ImagePointType _inv_lin_vec;

        /// basis function coefficient matrix
        CoeffMatrixType _coeff_mat;

        void _build_inv_lin_trafo(const TrafoEvaluator& trafo_eval)
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
          const DataType g = Math::sqrt(DataType(1) / DataType(3));
          const typename FacetEvalTraits::DomainPointType g1(-g), g2(+g);
          DomainPointType q1, q2;

          // loop over all 4 edges of the quad
          for(int i(0); i < 4; ++i)
          {
            // initialize facet evaluator
            facet_eval.prepare(Index(facet_index_set(cell, i)));

            // map first cubature point
            facet_eval(facet_data, g1);
            DataType w1(facet_data.jac_det);
            q1.set_mat_vec_mult(_inv_lin_mat, facet_data.img_point - _inv_lin_vec);

            // map second cubature point
            facet_eval(facet_data, g2);
            DataType w2(facet_data.jac_det);
            q2.set_mat_vec_mult(_inv_lin_mat, facet_data.img_point - _inv_lin_vec);

            // release facet evaluator
            facet_eval.finish();

            // compute reciprocal edge length
            DataType v = DataType(1) / (w1 + w2);

            // evaluate node functionals
            _nodal_mat(0,i) = DataType(1); // = v*(w1 + w2)
            _nodal_mat(1,i) = v*(w1*q1[0] + w2*q2[0]);
            _nodal_mat(2,i) = v*(w1*q1[1] + w2*q2[1]);
            _nodal_mat(3,i) = v*(w1*((q1[0]+q1[1])*(q1[0]-q1[1])) + w2*((q2[0]+q2[1])*(q2[0]-q2[1])));
            _nodal_mat(4,i) = v*(w1*q1[0]*q1[1] + w2*q2[0]*q2[1]);
          }

          // loop over all 4 cubature points
          InvLinTrafoData trafo_data;
          DataType vol_t = DataType(0);
          for(int i(0); i < 4; ++i)
          {
            DomainPointType dom_point(DataType(0));
            dom_point[0] = DataType(1 - 2*(i&1)) * g;
            dom_point[1] = DataType(1 - 1*(i&2)) * g;
            trafo_eval(trafo_data, dom_point);
            dom_point.set_mat_vec_mult(_inv_lin_mat, trafo_data.img_point - _inv_lin_vec);
            const DataType x = dom_point[0];
            const DataType y = dom_point[1];
            const DataType w = trafo_data.jac_det * x*y;
            vol_t += trafo_data.jac_det;
            _nodal_mat(0,4) += w;
            _nodal_mat(1,4) += w*x;
            _nodal_mat(2,4) += w*y;
            _nodal_mat(3,4) += w*(x+y)*(x-y);
            _nodal_mat(4,4) += w*x*y;
          }
          vol_t = DataType(1) / vol_t;
          for(int i(0); i < 5; ++i)
            _nodal_mat(i,4) *= vol_t;

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
          return 5;
        }

        /**
         * \brief Prepares the evaluator for a given cell.
         *
         * \param[in] trafo_eval
         * A reference to the trafo evaluator containing the cell information.
         */
        void prepare(const TrafoEvaluator& trafo_eval)
        {
          // prepare inverse linearized trafo
          _build_inv_lin_trafo(trafo_eval);

          // build coefficient matrix
          _build_coeff_matrix(trafo_eval);
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

          // "evaluate" monomials
          const DataType x(pt[0]);
          const DataType y(pt[1]);
          const DataType r((x+y)*(x-y));
          const DataType b(x*y);

          // loop over all basis functions
          for(int i(0); i < 5; ++i)
          {
            data.phi[i].value = _coeff_mat(i,0) + _coeff_mat(i,1)*x + _coeff_mat(i,2)*y + _coeff_mat(i,3)*r + _coeff_mat(i,4)*b;
          }
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
          for(int i(0); i < 5; ++i)
          {
            // calculate local gradients
            loc_grad(0) = _coeff_mat(i,1) + DataType(2) * _coeff_mat(i,3) * pt[0] + _coeff_mat(i,4) * pt[1];
            loc_grad(1) = _coeff_mat(i,2) - DataType(2) * _coeff_mat(i,3) * pt[1] + _coeff_mat(i,4) * pt[0];

            // multiply by transpose inverse linearized trafo matrix for chain rule
            data.phi[i].grad.set_vec_mat_mult(loc_grad, _inv_lin_mat);
          }
        }
      }; // Evaluator<..., Shape::Quadrilateral>

      /**
       * \brief Implementation of 3D Q1TBNP evaluator
       *
       * \author Peter Zajac
       */
      template<
        typename Space_,
        typename TrafoEvaluator_,
        typename SpaceEvalTraits_>
      class Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Shape::Hexahedron> :
        public EvaluatorBase<
          Evaluator<Space_, TrafoEvaluator_, SpaceEvalTraits_, Shape::Hexahedron>,
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
        typedef typename MeshType::template IndexSet<3,2>::Type FacetIndexSetType;

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
        static constexpr TrafoTags inv_lin_trafo_config = TrafoTags::dom_point | TrafoTags::img_point | TrafoTags::jac_inv | TrafoTags::jac_det;

        /// inverse linearized trafo data
        typedef typename TrafoEvaluator::template ConfigTraits<inv_lin_trafo_config>::EvalDataType InvLinTrafoData;

        /// trafo evaluator for facets
        typedef typename TrafoType::template Evaluator<Shape::Hypercube<2>, DataType>::Type FacetTrafoEvaluator;
        /// facet evaluation traits
        typedef typename FacetTrafoEvaluator::EvalTraits FacetEvalTraits;

        /// trafo config for facet trafo
        static constexpr TrafoTags facet_trafo_config = TrafoTags::img_point | TrafoTags::jac_det;

        /// facet trafo data
        typedef typename FacetTrafoEvaluator::template ConfigTraits<facet_trafo_config>::EvalDataType FacetTrafoData;

        /// basis function coefficient matrix
        typedef Tiny::Matrix<DataType, 9, 9> CoeffMatrixType;

        /// inverse linearized trafo matrix
        JacobianInverseType _inv_lin_mat;

        // inverse linearized trafo vector
        ImagePointType _inv_lin_vec;

        /// basis function coefficient matrix
        CoeffMatrixType _coeff_mat;

        void _build_inv_lin_trafo(const TrafoEvaluator& trafo_eval)
        {
          // create a domain point in the barycenter of the cell
          DomainPointType dom_point(DataType(0));

          // create the trafo data
          InvLinTrafoData trafo_data;
          trafo_eval(trafo_data, dom_point);

          // store inverse trafo linearization
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
          const FacetIndexSetType& facet_index_set = mesh.template get_index_set<3,2>();

          // fetch the cell index of the currently active cell
          const Index cell = trafo_eval.get_cell_index();

          // create a nodal matrix
          CoeffMatrixType _nodal_mat(DataType(0));

          // create a facet trafo evaluator
          FacetTrafoEvaluator facet_eval(trafo);

          // create facet evaluation data
          FacetTrafoData facet_data;

          // define 2-point Gauss cubature point coordinate
          const DataType g = Math::sqrt(DataType(1) / DataType(3));
          typename FacetEvalTraits::DomainPointType g1, g2, g3, g4;
          DomainPointType q1, q2, q3, q4;

          // set up cubature points
          g1[0] = g2[0] = -g;
          g3[0] = g4[0] = +g;
          g1[1] = g3[1] = -g;
          g2[1] = g4[1] = +g;

          // loop over all 6 faces of the hexa
          for(int i(0); i < 6; ++i)
          {
            // initialize facet evaluator
            facet_eval.prepare(Index(facet_index_set(cell, i)));

            // map first cubature point
            facet_eval(facet_data, g1);
            DataType w1(facet_data.jac_det);
            q1.set_mat_vec_mult(_inv_lin_mat, facet_data.img_point - _inv_lin_vec);

            // map second cubature point
            facet_eval(facet_data, g2);
            DataType w2(facet_data.jac_det);
            q2.set_mat_vec_mult(_inv_lin_mat, facet_data.img_point - _inv_lin_vec);

            // map third cubature point
            facet_eval(facet_data, g3);
            DataType w3(facet_data.jac_det);
            q3.set_mat_vec_mult(_inv_lin_mat, facet_data.img_point - _inv_lin_vec);

            // map fourth cubature point
            facet_eval(facet_data, g4);
            DataType w4(facet_data.jac_det);
            q4.set_mat_vec_mult(_inv_lin_mat, facet_data.img_point - _inv_lin_vec);

            // release facet evaluator
            facet_eval.finish();

            // compute reciprocal quad area
            DataType v = DataType(1) / (w1 + w2 + w3 + w4);

            // evaluate node functionals
            _nodal_mat(0,i) = DataType(1); // = v*(w1 + w2 + w3 + w4)
            _nodal_mat(1,i) = v*(w1*q1[0] + w2*q2[0] + w3*q3[0] + w4*q4[0]);
            _nodal_mat(2,i) = v*(w1*q1[1] + w2*q2[1] + w3*q3[1] + w4*q4[1]);
            _nodal_mat(3,i) = v*(w1*q1[2] + w2*q2[2] + w3*q3[2] + w4*q4[2]);
            _nodal_mat(4,i) = v*( // x^2 - y^2
              w1*(q1[0] + q1[1])*(q1[0] - q1[1]) + w2*(q2[0] + q2[1])*(q2[0] - q2[1]) +
              w3*(q3[0] + q3[1])*(q3[0] - q3[1]) + w4*(q4[0] + q4[1])*(q4[0] - q4[1]));
            _nodal_mat(5,i) = v*( // y^2 - z^2
              w1*(q1[1] + q1[2])*(q1[1] - q1[2]) + w2*(q2[1] + q2[2])*(q2[1] - q2[2]) +
              w3*(q3[1] + q3[2])*(q3[1] - q3[2]) + w4*(q4[1] + q4[2])*(q4[1] - q4[2]));
            _nodal_mat(6,i) = v*(w1*q1[0]*q1[1] + w2*q2[0]*q2[1] + w3*q3[0]*q3[1] + w4*q4[0]*q4[1]);
            _nodal_mat(7,i) = v*(w1*q1[1]*q1[2] + w2*q2[1]*q2[2] + w3*q3[1]*q3[2] + w4*q4[1]*q4[2]);
            _nodal_mat(8,i) = v*(w1*q1[2]*q1[0] + w2*q2[2]*q2[0] + w3*q3[2]*q3[0] + w4*q4[2]*q4[0]);
          }

          // loop over all 8 cubature points
          InvLinTrafoData trafo_data;
          DataType vol_t = DataType(0);
          for(int i(0); i < 8; ++i)
          {
            DomainPointType dom_point(DataType(0));
            dom_point[0] = DataType(1 - ((i&1) << 1)) * g;
            dom_point[1] = DataType(1 - ((i&2)     )) * g;
            dom_point[2] = DataType(1 - ((i&4) >> 1)) * g;
            trafo_eval(trafo_data, dom_point);
            dom_point.set_mat_vec_mult(_inv_lin_mat, trafo_data.img_point - _inv_lin_vec);
            const DataType x = dom_point[0];
            const DataType y = dom_point[1];
            const DataType z = dom_point[2];
            const DataType w = trafo_data.jac_det;
            vol_t += trafo_data.jac_det;
            _nodal_mat(0,6) += w*x*y;
            _nodal_mat(1,6) += w*x*x*y;
            _nodal_mat(2,6) += w*y*x*y;
            _nodal_mat(3,6) += w*z*x*y;
            _nodal_mat(4,6) += w*(x+y)*(x-y)*x*y;
            _nodal_mat(5,6) += w*(y+z)*(y-z)*x*y;
            _nodal_mat(6,6) += w*x*y*x*y;
            _nodal_mat(7,6) += w*y*z*x*y;
            _nodal_mat(8,6) += w*z*x*x*y;

            _nodal_mat(0,7) += w*y*z;
            _nodal_mat(1,7) += w*x*y*z;
            _nodal_mat(2,7) += w*y*y*z;
            _nodal_mat(3,7) += w*z*y*z;
            _nodal_mat(4,7) += w*(x+y)*(x-y)*y*z;
            _nodal_mat(5,7) += w*(y+z)*(y-z)*y*z;
            _nodal_mat(6,7) += w*x*y*y*z;
            _nodal_mat(7,7) += w*y*z*y*z;
            _nodal_mat(8,7) += w*z*x*y*z;

            _nodal_mat(0,8) += w*z*x;
            _nodal_mat(1,8) += w*x*z*x;
            _nodal_mat(2,8) += w*y*z*x;
            _nodal_mat(3,8) += w*z*z*x;
            _nodal_mat(4,8) += w*(x+y)*(x-y)*z*x;
            _nodal_mat(5,8) += w*(y+z)*(y-z)*z*x;
            _nodal_mat(6,8) += w*x*y*z*x;
            _nodal_mat(7,8) += w*y*z*z*x;
            _nodal_mat(8,8) += w*z*x*z*x;
          }
          vol_t = DataType(1) / vol_t;
          for(int i(0); i < 9; ++i)
          {
            _nodal_mat(i,6) *= vol_t;
            _nodal_mat(i,7) *= vol_t;
            _nodal_mat(i,8) *= vol_t;
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
          return 9;
        }

        /**
        * \brief Prepares the evaluator for a given cell.
        *
        * \param[in] trafo_eval
        * A reference to the trafo evaluator containing the cell information.
        */
        void prepare(const TrafoEvaluator& trafo_eval)
        {
          // prepare inverse linearized trafo
          _build_inv_lin_trafo(trafo_eval);

          // build coefficient matrix
          _build_coeff_matrix(trafo_eval);
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

          // "evaluate" monomials
          DataType x(pt[0]);
          DataType y(pt[1]);
          DataType z(pt[2]);
          DataType rxy((x+y)*(x-y)); // = x^2 - y^2
          DataType ryz((y+z)*(y-z)); // = y^2 - z^2

          // loop over all basis functions
          for(int i(0); i < 9; ++i)
          {
            data.phi[i].value = _coeff_mat(i,0) + _coeff_mat(i,4)*rxy + _coeff_mat(i,5)*ryz
              + x*(_coeff_mat(i,1) + y*_coeff_mat(i,6))
              + y*(_coeff_mat(i,2) + z*_coeff_mat(i,7))
              + z*(_coeff_mat(i,3) + x*_coeff_mat(i,8));
          }
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
          for(int i(0); i < 9; ++i)
          {
            // calculate local gradients
            loc_grad[0] = _coeff_mat(i,1) + pt[1]*_coeff_mat(i,6) + DataType(2) * pt[0] * _coeff_mat(i,4);
            loc_grad[1] = _coeff_mat(i,2) + pt[2]*_coeff_mat(i,7) + DataType(2) * pt[1] * (_coeff_mat(i,5) - _coeff_mat(i,4));
            loc_grad[2] = _coeff_mat(i,3) + pt[0]*_coeff_mat(i,8) - DataType(2) * pt[2] * _coeff_mat(i,5);

            // multiply by transpose inverse linearized trafo matrix for chain rule
            data.phi[i].grad.set_vec_mat_mult(loc_grad, _inv_lin_mat);
          }
        }
      }; // Evaluator<..., Shape::Hexahedron>
    } // namespace Q1TBNP
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_Q1TBNP_EVALUATOR_HPP
