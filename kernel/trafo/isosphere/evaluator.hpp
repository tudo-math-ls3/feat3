// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_TRAFO_ISOSPHERE_EVALUATOR_HPP
#define KERNEL_TRAFO_ISOSPHERE_EVALUATOR_HPP 1

// includes, FEAT
#include <kernel/trafo/evaluator_base.hpp>

namespace FEAT
{
  namespace Trafo
  {
    namespace IsoSphere
    {
      /// \cond internal
      namespace Intern
      {
        // p0, p1 and p2 are the 1D basis functions on the reference interval [-1,+1].
        // These are used for the tensor-product approach in the Hypercube evaluators.

        // basis function for left vertex
        template<typename T_>
        inline T_ p0(T_ x)
        {
          return T_(0.5) * x * (x - T_(1));
        }

        // basis function for right vertex
        template<typename T_>
        inline T_ p1(T_ x)
        {
          return T_(0.5) * x * (x + T_(1));
        }

        // basis function for edge midpoint
        template<typename T_>
        inline T_ p2(T_ x)
        {
          return (T_(1) - x) * (T_(1) + x);
        }

        // first order derivatives

        template<typename T_>
        inline T_ d1p0(T_ x)
        {
          return x - T_(0.5);
        }

        template<typename T_>
        inline T_ d1p1(T_ x)
        {
          return x + T_(0.5);
        }

        template<typename T_>
        inline T_ d1p2(T_ x)
        {
          return -T_(2) * x;
        }

        // second order derivatives

        template<typename T_>
        inline T_ d2p0(T_)
        {
          return T_(1);
        }

        template<typename T_>
        inline T_ d2p1(T_)
        {
          return T_(1);
        }

        template<typename T_>
        inline T_ d2p2(T_)
        {
          return -T_(2);
        }
      } // namespace Intern
      /// \endcond

      /**
       * \brief Iso-Sphere trafo evaluator class template
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename EvalPolicy_,
        typename Shape_ = typename EvalPolicy_::ShapeType>
      class Evaluator DOXY({});

      /* ************************************************************************************* */

      /**
       * \brief Specialisation of Iso-Sphere trafo evaluator for Vertex shape
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
        Tiny::Vector<DataType, image_dim> _coeff;

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

          // fetch the vertex and normalise
          _coeff = this->_trafo.get_mesh().get_vertex_set()[cell_index];
          _coeff.normalise();
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
          img_point = _coeff;
        }
      }; // class Evaluator<Vertex,...>


      /* ************************************************************************************* */

      /**
       * \brief Specialisation of Iso-Sphere trafo evaluator for Simplex<1> shape
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename EvalPolicy_>
      class Evaluator<Trafo_, EvalPolicy_, Shape::Simplex<1> > :
        public EvaluatorBase<Trafo_, Evaluator<Trafo_, EvalPolicy_, Shape::Simplex<1> >, EvalPolicy_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Trafo_, Evaluator, EvalPolicy_> BaseClass;
        /// shape type
        typedef Shape::Hypercube<1> ShapeType;
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

        static_assert(domain_dim < image_dim, "this trafo is for surfaces only");

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten;

      protected:
        /// the coefficients of the trafo
        Tiny::Matrix<DataType, 3, image_dim> _coeffs;

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

          // get vertex coords
          _coeffs[0] = vertex_set[index_set(cell_index, 0)];
          _coeffs[1] = vertex_set[index_set(cell_index, 1)];

          // get edge midpoint
          _coeffs[2] = _coeffs[0] + _coeffs[1];

          // normalise all coefficients
          for(int i(0); i < 3; ++i)
            _coeffs[i].normalise();
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
          const DataType x = dom_point[0];

          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] =
              _coeffs[0][i] * DataType(0.5) * x * (x - DataType(1)) +
              _coeffs[1][i] * DataType(0.5) * x * (x + DataType(1)) +
              _coeffs[2][i] * (DataType(1) - x) * (DataType(1) - x);
          }
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
          const DataType x = dom_point[0];

          for(int i(0); i < image_dim; ++i)
          {
            jac_mat[i][0] =
              _coeffs[0][i] * (x - DataType(0.5)) +
              _coeffs[1][i] * (x + DataType(0.5)) -
              _coeffs[2][i] * DataType(2) * x;
          }
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
        void calc_hess_ten(HessianTensorType& hess_ten, const DomainPointType& DOXY(dom_point)) const
        {
          for(int i(0); i < image_dim; ++i)
          {
            hess_ten[i][0] = _coeffs[0][i] + _coeffs[1][i] - DataType(2) *_coeffs[2][i];
          }
        }
      }; // class Evaluator<Simplex<1>,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialisation of Iso-Sphere trafo evaluator for Simplex<2> shape
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename EvalPolicy_>
      class Evaluator<Trafo_, EvalPolicy_, Shape::Simplex<2> > :
        public EvaluatorBase<Trafo_, Evaluator<Trafo_, EvalPolicy_, Shape::Simplex<2> >, EvalPolicy_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Trafo_, Evaluator, EvalPolicy_> BaseClass;
        /// shape type
        typedef Shape::Simplex<2> ShapeType;
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

        static_assert(domain_dim < image_dim, "this trafo is for surfaces only");

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten;

      protected:
        /// the coefficients of the trafo
        Tiny::Matrix<DataType, 6, image_dim> _coeffs;

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

          // get vertex coords
          _coeffs[0] = vertex_set[index_set(cell_index, 0)];
          _coeffs[1] = vertex_set[index_set(cell_index, 1)];
          _coeffs[2] = vertex_set[index_set(cell_index, 2)];

          // get edge midpoints
          _coeffs[3] = _coeffs[1] + _coeffs[2];
          _coeffs[4] = _coeffs[0] + _coeffs[2];
          _coeffs[5] = _coeffs[0] + _coeffs[1];

          // normalise all coefficients
          for(int i(0); i < 6; ++i)
            _coeffs[i].normalise();
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
          const DataType x = dom_point[0];
          const DataType y = dom_point[1];

          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] =
              _coeffs[0][i] * (DataType(2) * (x + y - DataType(0.5)) * (x + y - DataType(1))) +
              _coeffs[1][i] * (DataType(2) * x * (x - DataType(0.5))) +
              _coeffs[2][i] * (DataType(2) * y * (y - DataType(0.5))) +
              _coeffs[3][i] * (DataType(4) * x * y) +
              _coeffs[4][i] * (DataType(4) * y * (DataType(1) - x - y)) +
              _coeffs[5][i] * (DataType(4) * x * (DataType(1) - x - y));
          }
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
          const DataType x = dom_point[0];
          const DataType y = dom_point[1];

          for(int i(0); i < image_dim; ++i)
          {
            jac_mat[i][0] =
              _coeffs[0][i] * (DataType(4) * (x + y) - DataType(3)) +
              _coeffs[1][i] * (DataType(4) * x - DataType(1)) +
            //_coeffs[2][i] * (DataType(0)) +
              _coeffs[3][i] * (DataType(4) * y) -
              _coeffs[4][i] * (DataType(4) * y) -
              _coeffs[5][i] * (DataType(4) * (DataType(2)*x + y - DataType(1)));
            jac_mat[i][1] =
              _coeffs[0][i] * (DataType(4) * (x + y) - DataType(3)) +
            //_coeffs[1][i] * (DataType(0)) +
              _coeffs[2][i] * (DataType(4) * y - DataType(1)) +
              _coeffs[3][i] * (DataType(4) * x) -
              _coeffs[4][i] * (DataType(4) * (DataType(2)*y + x - DataType(1))) -
              _coeffs[5][i] * (DataType(4) * x);
          }
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
        void calc_hess_ten(HessianTensorType& hess_ten, const DomainPointType& DOXY(dom_point)) const
        {
          for(int i(0); i < image_dim; ++i)
          {
            hess_ten[i][0][0] = DataType(4) * (_coeffs[0][i] + _coeffs[1][i] - DataType(2) * _coeffs[5][i]);
            hess_ten[i][1][1] = DataType(4) * (_coeffs[0][i] + _coeffs[2][i] - DataType(2) * _coeffs[4][i]);
            hess_ten[i][0][1] =
            hess_ten[i][1][0] = DataType(4) * (_coeffs[0][i] + _coeffs[3][i] + _coeffs[4][i] + _coeffs[5][i]);
          }
        }
      }; // class Evaluator<Simplex<2>,...>


      /* ************************************************************************************* */

      /**
       * \brief Specialisation of Iso-Sphere trafo evaluator for Hypercube<1> shape
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename EvalPolicy_>
      class Evaluator<Trafo_, EvalPolicy_, Shape::Hypercube<1> > :
        public EvaluatorBase<Trafo_, Evaluator<Trafo_, EvalPolicy_, Shape::Hypercube<1> >, EvalPolicy_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Trafo_, Evaluator, EvalPolicy_> BaseClass;
        /// shape type
        typedef Shape::Hypercube<1> ShapeType;
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

        static_assert(domain_dim < image_dim, "this trafo is for surfaces only");

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten;

      protected:
        /// the coefficients of the trafo
        Tiny::Matrix<DataType, 3, image_dim> _coeffs;

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

          // get vertex coords
          _coeffs[0] = vertex_set[index_set(cell_index, 0)];
          _coeffs[1] = vertex_set[index_set(cell_index, 1)];

          // compute edge midpoints
          _coeffs[2] = _coeffs[0] + _coeffs[1];

          // normalise all coefficients
          for(int i(0); i < 3; ++i)
            _coeffs[i].normalise();
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
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] =
              _coeffs[0][i] * Intern::p0(dom_point[0]) +
              _coeffs[1][i] * Intern::p1(dom_point[0]) +
              _coeffs[2][i] * Intern::p2(dom_point[0]);
          }
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
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat[i][0] =
              _coeffs[0][i] * Intern::d1p0(dom_point[0]) +
              _coeffs[1][i] * Intern::d1p1(dom_point[0]) +
              _coeffs[2][i] * Intern::d1p2(dom_point[0]);
          }
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
          for(int i(0); i < image_dim; ++i)
          {
            hess_ten[i][0][0] =
              _coeffs[0][i] * Intern::d2p0(dom_point[0]) +
              _coeffs[1][i] * Intern::d2p1(dom_point[0]) +
              _coeffs[2][i] * Intern::d2p2(dom_point[0]);
          }
        }
      }; // class Evaluator<Hypercube<1>,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialisation of Iso-Sphere trafo evaluator for Hypercube<2> shape
       *
       * \author Peter Zajac
       */
      template<
        typename Trafo_,
        typename EvalPolicy_>
      class Evaluator<Trafo_, EvalPolicy_, Shape::Hypercube<2> > :
        public EvaluatorBase<Trafo_, Evaluator<Trafo_, EvalPolicy_, Shape::Hypercube<2> >, EvalPolicy_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Trafo_, Evaluator, EvalPolicy_> BaseClass;
        /// shape type
        typedef Shape::Hypercube<2> ShapeType;
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

        static_assert(domain_dim < image_dim, "this trafo is for surfaces only");

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten;

      protected:
        /// the coefficients of the trafo
        Tiny::Matrix<DataType, 9, image_dim> _coeffs;

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

          // get vertex coords
          _coeffs[0] = vertex_set[index_set(cell_index, 0)];
          _coeffs[1] = vertex_set[index_set(cell_index, 1)];
          _coeffs[2] = vertex_set[index_set(cell_index, 2)];
          _coeffs[3] = vertex_set[index_set(cell_index, 3)];

          // compute edge midpoints
          _coeffs[4] = _coeffs[0] + _coeffs[1];
          _coeffs[5] = _coeffs[2] + _coeffs[3];
          _coeffs[6] = _coeffs[0] + _coeffs[2];
          _coeffs[7] = _coeffs[1] + _coeffs[3];

          // compute quad midpoint
          _coeffs[8] = _coeffs[0] + _coeffs[1] + _coeffs[2] + _coeffs[3];

          // normalise all coefficients
          for(int i(0); i < 9; ++i)
            _coeffs[i].normalise();
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
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] =
              _coeffs[0][i] * Intern::p0(dom_point[0]) * Intern::p0(dom_point[1]) +
              _coeffs[1][i] * Intern::p1(dom_point[0]) * Intern::p0(dom_point[1]) +
              _coeffs[2][i] * Intern::p0(dom_point[0]) * Intern::p1(dom_point[1]) +
              _coeffs[3][i] * Intern::p1(dom_point[0]) * Intern::p1(dom_point[1]) +
              _coeffs[4][i] * Intern::p2(dom_point[0]) * Intern::p0(dom_point[1]) +
              _coeffs[5][i] * Intern::p2(dom_point[0]) * Intern::p1(dom_point[1]) +
              _coeffs[6][i] * Intern::p0(dom_point[0]) * Intern::p2(dom_point[1]) +
              _coeffs[7][i] * Intern::p1(dom_point[0]) * Intern::p2(dom_point[1]) +
              _coeffs[8][i] * Intern::p2(dom_point[0]) * Intern::p2(dom_point[1]);
          }
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
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat[i][0] =
              _coeffs[0][i] * Intern::d1p0(dom_point[0]) * Intern::p0(dom_point[1]) +
              _coeffs[1][i] * Intern::d1p1(dom_point[0]) * Intern::p0(dom_point[1]) +
              _coeffs[2][i] * Intern::d1p0(dom_point[0]) * Intern::p1(dom_point[1]) +
              _coeffs[3][i] * Intern::d1p1(dom_point[0]) * Intern::p1(dom_point[1]) +
              _coeffs[4][i] * Intern::d1p2(dom_point[0]) * Intern::p0(dom_point[1]) +
              _coeffs[5][i] * Intern::d1p2(dom_point[0]) * Intern::p1(dom_point[1]) +
              _coeffs[6][i] * Intern::d1p0(dom_point[0]) * Intern::p2(dom_point[1]) +
              _coeffs[7][i] * Intern::d1p1(dom_point[0]) * Intern::p2(dom_point[1]) +
              _coeffs[8][i] * Intern::d1p2(dom_point[0]) * Intern::p2(dom_point[1]);
            jac_mat[i][1] =
              _coeffs[0][i] * Intern::p0(dom_point[0]) * Intern::d1p0(dom_point[1]) +
              _coeffs[1][i] * Intern::p1(dom_point[0]) * Intern::d1p0(dom_point[1]) +
              _coeffs[2][i] * Intern::p0(dom_point[0]) * Intern::d1p1(dom_point[1]) +
              _coeffs[3][i] * Intern::p1(dom_point[0]) * Intern::d1p1(dom_point[1]) +
              _coeffs[4][i] * Intern::p2(dom_point[0]) * Intern::d1p0(dom_point[1]) +
              _coeffs[5][i] * Intern::p2(dom_point[0]) * Intern::d1p1(dom_point[1]) +
              _coeffs[6][i] * Intern::p0(dom_point[0]) * Intern::d1p2(dom_point[1]) +
              _coeffs[7][i] * Intern::p1(dom_point[0]) * Intern::d1p2(dom_point[1]) +
              _coeffs[8][i] * Intern::p2(dom_point[0]) * Intern::d1p2(dom_point[1]);
          }
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
          for(int i(0); i < image_dim; ++i)
          {
            hess_ten[i][0][0] =
              _coeffs[0][i] * Intern::d2p0(dom_point[0]) * Intern::p0(dom_point[1]) +
              _coeffs[1][i] * Intern::d2p1(dom_point[0]) * Intern::p0(dom_point[1]) +
              _coeffs[2][i] * Intern::d2p0(dom_point[0]) * Intern::p1(dom_point[1]) +
              _coeffs[3][i] * Intern::d2p1(dom_point[0]) * Intern::p1(dom_point[1]) +
              _coeffs[4][i] * Intern::d2p2(dom_point[0]) * Intern::p0(dom_point[1]) +
              _coeffs[5][i] * Intern::d2p2(dom_point[0]) * Intern::p1(dom_point[1]) +
              _coeffs[6][i] * Intern::d2p0(dom_point[0]) * Intern::p2(dom_point[1]) +
              _coeffs[7][i] * Intern::d2p1(dom_point[0]) * Intern::p2(dom_point[1]) +
              _coeffs[8][i] * Intern::d2p2(dom_point[0]) * Intern::p2(dom_point[1]);
            hess_ten[i][1][1] =
              _coeffs[0][i] * Intern::p0(dom_point[0]) * Intern::d2p0(dom_point[1]) +
              _coeffs[1][i] * Intern::p1(dom_point[0]) * Intern::d2p0(dom_point[1]) +
              _coeffs[2][i] * Intern::p0(dom_point[0]) * Intern::d2p1(dom_point[1]) +
              _coeffs[3][i] * Intern::p1(dom_point[0]) * Intern::d2p1(dom_point[1]) +
              _coeffs[4][i] * Intern::p2(dom_point[0]) * Intern::d2p0(dom_point[1]) +
              _coeffs[5][i] * Intern::p2(dom_point[0]) * Intern::d2p1(dom_point[1]) +
              _coeffs[6][i] * Intern::p0(dom_point[0]) * Intern::d2p2(dom_point[1]) +
              _coeffs[7][i] * Intern::p1(dom_point[0]) * Intern::d2p2(dom_point[1]) +
              _coeffs[8][i] * Intern::p2(dom_point[0]) * Intern::d2p2(dom_point[1]);
            hess_ten[i][0][1] =
            hess_ten[i][1][0] =
              _coeffs[0][i] * Intern::d1p0(dom_point[0]) * Intern::d1p0(dom_point[1]) +
              _coeffs[1][i] * Intern::d1p1(dom_point[0]) * Intern::d1p0(dom_point[1]) +
              _coeffs[2][i] * Intern::d1p0(dom_point[0]) * Intern::d1p1(dom_point[1]) +
              _coeffs[3][i] * Intern::d1p1(dom_point[0]) * Intern::d1p1(dom_point[1]) +
              _coeffs[4][i] * Intern::d1p2(dom_point[0]) * Intern::d1p0(dom_point[1]) +
              _coeffs[5][i] * Intern::d1p2(dom_point[0]) * Intern::d1p1(dom_point[1]) +
              _coeffs[6][i] * Intern::d1p0(dom_point[0]) * Intern::d1p2(dom_point[1]) +
              _coeffs[7][i] * Intern::d1p1(dom_point[0]) * Intern::d1p2(dom_point[1]) +
              _coeffs[8][i] * Intern::d1p2(dom_point[0]) * Intern::d1p2(dom_point[1]);
          }
        }
      }; // class Evaluator<Hypercube<2>,...>
    } // namespace IsoSphere
  } // namespace Trafo
} // namespace FEAT

#endif // KERNEL_TRAFO_ISOSPHERE_EVALUATOR_HPP
