// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_TRAFO_STANDARD_EVALUATOR_HPP
#define KERNEL_TRAFO_STANDARD_EVALUATOR_HPP 1

// includes, FEAT
#include <kernel/trafo/evaluator_base.hpp>
#include "details.hpp"

namespace FEAT
{
  namespace Trafo
  {
    namespace Standard
    {
      /**
       * \brief Standard trafo evaluator class template
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename EvalPolicy_,
        typename Shape_ = typename EvalPolicy_::ShapeType>
      class Evaluator:
        public EvaluatorBase<Trafo_, Evaluator<Trafo_, EvalPolicy_, Shape_ >, EvalPolicy_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Trafo_, Evaluator, EvalPolicy_> BaseClass;
        /// shape type
        typedef Shape_ ShapeType;
        /// trafo type using this evaluator
        typedef Trafo_ TrafoType;
        /// trafo evaluation traits
        typedef EvalPolicy_ EvalPolicy;

        /// type of the underlying mesh
        typedef typename TrafoType::MeshType MeshType;

        /// evaluation data type
        typedef typename EvalPolicy::DataType DataType;
        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;
        /// image point type
        typedef typename EvalPolicy::ImagePointType ImagePointType;
        /// jacobian matrix type
        typedef typename EvalPolicy::JacobianMatrixType JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef typename EvalPolicy::JacobianInverseType JacobianInverseType;
        /// jacobian determinant type
        typedef typename EvalPolicy::JacobianDeterminantType JacobianDeterminantType;
        /// hessian tensor type
        typedef typename EvalPolicy::HessianTensorType HessianTensorType;
        /// hessian inverse tensor type
        typedef typename EvalPolicy::HessianInverseType HessianInverseType;

        /// domain dimension
        static constexpr int domain_dim = EvalPolicy::domain_dim;
        /// image dimension
        static constexpr int image_dim = EvalPolicy::image_dim;
        /// number of verts
        static constexpr int num_verts = Shape::FaceTraits<ShapeType, 0>::count;

        /// evaluation helper, see details.hpp
        typedef EvalHelper<DataType, DomainPointType, ImagePointType, JacobianMatrixType, JacobianInverseType, JacobianDeterminantType,
                          HessianTensorType, HessianInverseType ,ShapeType, image_dim> EvalHp;

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten |
          (domain_dim == image_dim ? (TrafoTags::jac_inv | TrafoTags::hess_inv) : TrafoTags::none);

      protected:
        /// the coefficients of the trafo
        DataType _coeff[image_dim][num_verts];

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] trafo
         * A reference to the trafo using this evaluator.
         */
        explicit Evaluator(const TrafoType& trafo) :
          BaseClass(trafo)
        {
        }

        /**
         * \brief Prepares the evaluator for a given cell.
         *
         * \param[in] cell_index
         * The index of the cell for which the evaluator is to be prepared.
         */
        void prepare(Index cell_index)
        {
          // prepare base-class
          BaseClass::prepare(cell_index);

          // fetch the mesh from the trafo
          const MeshType& mesh = this->_trafo.get_mesh();

          // fetch the vertex set from the mesh
          typedef typename MeshType::VertexSetType VertexSetType;
          const VertexSetType& vertex_set = mesh.get_vertex_set();

          // fetch the index set
          typedef typename MeshType::template IndexSet<domain_dim, 0>::Type IndexSetType;
          const IndexSetType& index_set = mesh.template get_index_set<domain_dim, 0>();

          EvalHp::set_coefficients(_coeff, index_set[cell_index], vertex_set.begin());
        }

        /**
         * \brief Maps a point from the reference cell to the selected cell.
         *
         * \param[out] img_point
         * A reference to the point on the selected cell that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell that is to be mapped.
         */
        void map_point(ImagePointType& img_point, const DomainPointType& dom_point) const
        {
          EvalHp::map_point(img_point, dom_point, _coeff);
        }

        /**
         * \brief Calculates the jacobian matrix for a given point.
         *
         * \param[out] jac_mat
         * A reference to the jacobian matrix that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell where the jacobian matrix is to be computed.
         */
        void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& dom_point) const
        {
          EvalHp::calc_jac_mat(jac_mat, dom_point, _coeff);
        }

        /**
         * \brief Computes the hessian tensor for a given domain point.
         *
         * \param[out] hess_ten
         * A reference to the hessian tensor that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the domain point on the reference cell for which the hessian tensor is to be computed.
         */
        void calc_hess_ten(HessianTensorType& hess_ten, const DomainPointType& dom_point) const
        {
          EvalHp::calc_hess_ten(hess_ten, dom_point, _coeff);
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \returns
         * The volume of the current cell.
         */
        DataType volume() const
        {
          return EvalHp::volume(_coeff);
        }

        /**
         * \brief Computes and returns the directed mesh width.
         *
         * This function approximates the cell width along a given normalized ray direction vector.
         *
         * \param[in] ray
         * A (normalized) direction vector. Must not be a null vector.
         *
         * \returns
         * The mesh width in direction of the input ray vector.
         */
        DataType width_directed(const ImagePointType& ray) const
        {
          return EvalHp::width_directed(ray, _coeff);
        }
      }; // class Evaluator<Hypercube<1>,...>

      // /* ************************************************************************************* */

      /**
       * \brief Specialization of standard trafo evaluator for Vertex shape
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename EvalPolicy_>
      class Evaluator<Trafo_, EvalPolicy_, Shape::Vertex > :
        public EvaluatorBase<Trafo_, Evaluator<Trafo_, EvalPolicy_, Shape::Vertex >, EvalPolicy_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Trafo_, Evaluator, EvalPolicy_> BaseClass;
        /// shape type
        typedef Shape::Vertex ShapeType;
        /// trafo type using this evaluator
        typedef Trafo_ TrafoType;
        /// trafo evaluation traits
        typedef EvalPolicy_ EvalPolicy;

        /// type of the underlying mesh
        typedef typename TrafoType::MeshType MeshType;

        /// evaluation data type
        typedef typename EvalPolicy::DataType DataType;
        /// domain point type
        typedef typename EvalPolicy::DomainPointType DomainPointType;
        /// image point type
        typedef typename EvalPolicy::ImagePointType ImagePointType;
        /// jacobian matrix type
        typedef typename EvalPolicy::JacobianMatrixType JacobianMatrixType;
        /// jacobian inverse matrix type
        typedef typename EvalPolicy::JacobianInverseType JacobianInverseType;
        /// jacobian determinant type
        typedef typename EvalPolicy::JacobianDeterminantType JacobianDeterminantType;
        /// hessian tensor type
        typedef typename EvalPolicy::HessianTensorType HessianTensorType;
        /// hessian inverse tensor type
        typedef typename EvalPolicy::HessianInverseType HessianInverseType;

        /// domain dimension
        static constexpr int domain_dim = EvalPolicy::domain_dim;
        /// image dimension
        static constexpr int image_dim = EvalPolicy::image_dim;

        /// we can compute domain and image points
        static constexpr TrafoTags eval_caps = TrafoTags::dom_point | TrafoTags::img_point;

      protected:
        /// the coefficients of the trafo
        DataType _coeff[image_dim];

      public:
        /**
         * \brief Constructor.
         *
         * \param[in] trafo
         * A reference to the trafo using this evaluator.
         */
        explicit Evaluator(const TrafoType& trafo) :
          BaseClass(trafo)
        {
        }

        /**
         * \brief Prepares the evaluator for a given cell.
         *
         * \param[in] cell_index
         * The index of the cell for which the evaluator is to be prepared.
         */
        void prepare(Index cell_index)
        {
          // prepare base-class
          BaseClass::prepare(cell_index);

          // fetch the mesh from the trafo
          const MeshType& mesh = this->_trafo.get_mesh();

          // fetch the vertex set from the mesh
          typedef typename MeshType::VertexSetType VertexSetType;
          const VertexSetType& vertex_set = mesh.get_vertex_set();

          // fetch the vertex
          typedef typename VertexSetType::VertexType VertexType;
          const VertexType& vtx = vertex_set[cell_index];

          // calculate transformation coefficients
          for(int i(0); i < image_dim; ++i)
          {
            _coeff[i] = DataType(vtx[i]);
          }
        }

        /**
         * \brief Maps a point from the reference cell to the selected cell.
         *
         * \param[out] img_point
         * A reference to the point on the selected cell that is to be computed.
         *
         * \param[in] dom_point
         * A reference to the point on the reference cell that is to be mapped.
         */
        void map_point(ImagePointType& img_point, const DomainPointType& DOXY(dom_point)) const
        {
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] = _coeff[i];
          }
        }
      }; // class Evaluator<Vertex,...>
    } // namespace Standard
  } // namespace Trafo
} // namespace FEAT

#endif // KERNEL_TRAFO_STANDARD_EVALUATOR_HPP
