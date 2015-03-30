#pragma once
#ifndef KERNEL_TRAFO_STANDARD_EVALUATOR_HPP
#define KERNEL_TRAFO_STANDARD_EVALUATOR_HPP 1

// includes, FEAST
#include <kernel/trafo/evaluator_base.hpp>

namespace FEAST
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
      class Evaluator DOXY({});

      /* ************************************************************************************* */

      /**
       * \brief Specialisation of standard trafo evaluator for Vertex shape
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

        /// can compute domain points
        static constexpr bool can_dom_point = true;
        /// can compute image points
        static constexpr bool can_img_point = true;
        /// can compute jacobian matrices
        static constexpr bool can_jac_mat = false;
        /// can compute jacobian inverse matrices if domain and image dimensions coincide
        static constexpr bool can_jac_inv = false;
        /// can compute jacobian determinants
        static constexpr bool can_jac_det = false;
        /// can compute hessian tensors
        static constexpr bool can_hess_ten = false;
        /// can compute inverse hessian tensors if domain and image dimensions coincide
        static constexpr bool can_hess_inv = false;

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
          typedef typename VertexSetType::ConstVertexReference ConstVertexReference;
          ConstVertexReference vtx = vertex_set[cell_index];

          // calculate transformation coefficients
          for(Index i(0); i < Index(image_dim); ++i)
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
          for(Index i(0); i < Index(image_dim); ++i)
          {
            img_point[i] = _coeff[i];
          }
        }
      }; // class Evaluator<Vertex,...>


      /* ************************************************************************************* */

      /**
       * \brief Specialisation of standard trafo evaluator for Simplex<1> shape
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

        /// can compute domain points
        static constexpr bool can_dom_point = true;
        /// can compute image points
        static constexpr bool can_img_point = true;
        /// can compute jacobian matrices
        static constexpr bool can_jac_mat = true;
        /// can compute jacobian inverse matrices if domain and image dimensions coincide
        static constexpr bool can_jac_inv = (domain_dim == image_dim);
        /// can compute jacobian determinants
        static constexpr bool can_jac_det = true;
        /// can compute hessian tensors
        static constexpr bool can_hess_ten = true;
        /// can compute inverse hessian tensors if domain and image dimensions coincide
        static constexpr bool can_hess_inv = (domain_dim == image_dim);

      protected:
        /// the coefficients of the trafo
        DataType _coeff[image_dim][2];

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

          // fetch the vertices of the edge
          typedef typename VertexSetType::ConstVertexReference ConstVertexReference;
          ConstVertexReference v0 = vertex_set[index_set(cell_index, 0)];
          ConstVertexReference v1 = vertex_set[index_set(cell_index, 1)];

          // calculate transformation coefficients
          for(Index i(0); i < Index(image_dim); ++i)
          {
            _coeff[i][0] = DataType(v0[i]);
            _coeff[i][1] = DataType(v1[i] - v0[i]);
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
        void map_point(ImagePointType& img_point, const DomainPointType& dom_point) const
        {
          for(Index i(0); i < Index(image_dim); ++i)
          {
            img_point[i] = _coeff[i][0] + _coeff[i][1] * dom_point[0];
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
        void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& DOXY(dom_point)) const
        {
          for(Index i(0); i < Index(image_dim); ++i)
          {
            jac_mat(i,0) = _coeff[i][1];
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
          for(Index i(0); i < Index(image_dim); ++i)
          {
            hess_ten(i,0,0) = DataType(0);
          }
        }
      }; // class Evaluator<Simplex<1>,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialisation of standard trafo evaluator for Simplex<2> shape
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

        /// can compute domain points
        static constexpr bool can_dom_point = true;
        /// can compute image points
        static constexpr bool can_img_point = true;
        /// can compute jacobian matrices
        static constexpr bool can_jac_mat = true;
        /// can compute jacobian inverse matrices if domain and image dimensions coincide
        static constexpr bool can_jac_inv = (domain_dim == image_dim);
        /// can compute jacobian determinants
        static constexpr bool can_jac_det = true;
        /// can compute hessian tensors
        static constexpr bool can_hess_ten = true;
        /// can compute inverse hessian tensors if domain and image dimensions coincide
        static constexpr bool can_hess_inv = (domain_dim == image_dim);

      protected:
        /// the coefficients of the trafo
        DataType _coeff[image_dim][3];

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

          // fetch the vertices of the edge
          typedef typename VertexSetType::ConstVertexReference ConstVertexReference;
          ConstVertexReference v0 = vertex_set[index_set(cell_index, 0)];
          ConstVertexReference v1 = vertex_set[index_set(cell_index, 1)];
          ConstVertexReference v2 = vertex_set[index_set(cell_index, 2)];

          // calculate transformation coefficients
          for(Index i(0); i < Index(image_dim); ++i)
          {
            _coeff[i][0] = DataType(v0[i]);
            _coeff[i][1] = DataType(v1[i] - v0[i]);
            _coeff[i][2] = DataType(v2[i] - v0[i]);
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
        void map_point(ImagePointType& img_point, const DomainPointType& dom_point) const
        {
          for(Index i(0); i < Index(image_dim); ++i)
          {
            img_point[i] = _coeff[i][0] + _coeff[i][1] * dom_point[0] + _coeff[i][2] * dom_point[1];
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
        void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& DOXY(dom_point)) const
        {
          for(Index i(0); i < Index(image_dim); ++i)
          {
            jac_mat(i,0) = _coeff[i][1];
            jac_mat(i,1) = _coeff[i][2];
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
          for(Index i(0); i < Index(image_dim); ++i)
          {
            hess_ten.format(DataType(0));
          }
        }
      }; // class Evaluator<Simplex<2>,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialisation of standard trafo evaluator for Simplex<3> shape
       *
       * \author Stefan Wahlers
       */
      template<
        typename Trafo_,
        typename EvalPolicy_>
      class Evaluator<Trafo_, EvalPolicy_, Shape::Simplex<3> > :
        public EvaluatorBase<Trafo_, Evaluator<Trafo_, EvalPolicy_, Shape::Simplex<3> >, EvalPolicy_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Trafo_, Evaluator, EvalPolicy_> BaseClass;
        /// shape type
        typedef Shape::Simplex<3> ShapeType;
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

        /// can compute domain points
        static constexpr bool can_dom_point = true;
        /// can compute image points
        static constexpr bool can_img_point = true;
        /// can compute jacobian matrices
        static constexpr bool can_jac_mat = true;
        /// can compute jacobian inverse matrices if domain and image dimensions coincide
        static constexpr bool can_jac_inv = (domain_dim == image_dim);
        /// can compute jacobian determinants
        static constexpr bool can_jac_det = true;
        /// can compute hessian tensors
        static constexpr bool can_hess_ten = true;
        /// can compute inverse hessian tensors if domain and image dimensions coincide
        static constexpr bool can_hess_inv = (domain_dim == image_dim);

      protected:
        /// the coefficients of the trafo
        DataType _coeff[image_dim][4];

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

          // fetch the vertices of the edge
          typedef typename VertexSetType::ConstVertexReference ConstVertexReference;
          ConstVertexReference v0 = vertex_set[index_set(cell_index, 0)];
          ConstVertexReference v1 = vertex_set[index_set(cell_index, 1)];
          ConstVertexReference v2 = vertex_set[index_set(cell_index, 2)];
          ConstVertexReference v3 = vertex_set[index_set(cell_index, 3)];

          // calculate transformation coefficients
          for(Index i(0); i < Index(image_dim); ++i)
          {
            _coeff[i][0] = DataType(v0[i]);
            _coeff[i][1] = DataType(v1[i] - v0[i]);
            _coeff[i][2] = DataType(v2[i] - v0[i]);
            _coeff[i][3] = DataType(v3[i] - v0[i]);
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
        void map_point(ImagePointType& img_point, const DomainPointType& dom_point) const
        {
          for(Index i(0); i < Index(image_dim); ++i)
          {
            img_point[i] = _coeff[i][0] + _coeff[i][1] * dom_point[0]
                                        + _coeff[i][2] * dom_point[1]
                                        + _coeff[i][3] * dom_point[2];
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
        void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& DOXY(dom_point)) const
        {
          for(Index i(0); i < Index(image_dim); ++i)
          {
            jac_mat(i,0) = _coeff[i][1];
            jac_mat(i,1) = _coeff[i][2];
            jac_mat(i,2) = _coeff[i][3];
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
          for(Index i(0); i < Index(image_dim); ++i)
          {
            hess_ten.format(DataType(0));
          }
        }
      }; // class Evaluator<Simplex<3>,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialisation of standard trafo evaluator for Hypercube<1> shape
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

        /// can compute domain points
        static constexpr bool can_dom_point = true;
        /// can compute image points
        static constexpr bool can_img_point = true;
        /// can compute jacobian matrices
        static constexpr bool can_jac_mat = true;
        /// can compute jacobian inverse matrices if domain and image dimensions coincide
        static constexpr bool can_jac_inv = (domain_dim == image_dim);
        /// can compute jacobian determinants
        static constexpr bool can_jac_det = true;
        /// can compute hessian tensors
        static constexpr bool can_hess_ten = true;
        /// can compute inverse hessian tensors if domain and image dimensions coincide
        static constexpr bool can_hess_inv = (domain_dim == image_dim);

      protected:
        /// the coefficients of the trafo
        DataType _coeff[image_dim][2];

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

          // fetch the vertices of the edge
          typedef typename VertexSetType::ConstVertexReference ConstVertexReference;
          ConstVertexReference v0 = vertex_set[index_set(cell_index, 0)];
          ConstVertexReference v1 = vertex_set[index_set(cell_index, 1)];

          // calculate transformation coefficients
          for(Index i(0); i < Index(image_dim); ++i)
          {
            _coeff[i][0] = DataType(0.5) * DataType( v0[i] + v1[i]);
            _coeff[i][1] = DataType(0.5) * DataType(-v0[i] + v1[i]);
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
        void map_point(ImagePointType& img_point, const DomainPointType& dom_point) const
        {
          for(Index i(0); i < Index(image_dim); ++i)
          {
            img_point[i] = _coeff[i][0] + _coeff[i][1] * dom_point[0];
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
        void calc_jac_mat(JacobianMatrixType& jac_mat, const DomainPointType& DOXY(dom_point)) const
        {
          for(Index i(0); i < Index(image_dim); ++i)
          {
            jac_mat(i,0) = _coeff[i][1];
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
          for(Index i(0); i < Index(image_dim); ++i)
          {
            hess_ten(i,0,0) = DataType(0);
          }
        }
      }; // class Evaluator<Hypercube<1>,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialisation of standard trafo evaluator for Hypercube<2> shape
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

        /// can compute domain points
        static constexpr bool can_dom_point = true;
        /// can compute image points
        static constexpr bool can_img_point = true;
        /// can compute jacobian matrices
        static constexpr bool can_jac_mat = true;
        /// can compute jacobian inverse matrices if domain and image dimensions coincide
        static constexpr bool can_jac_inv = (domain_dim == image_dim);
        /// can compute jacobian determinants
        static constexpr bool can_jac_det = true;
        /// can compute hessian tensors
        static constexpr bool can_hess_ten = true;
        /// can compute inverse hessian tensors if domain and image dimensions coincide
        static constexpr bool can_hess_inv = (domain_dim == image_dim);

      protected:
        /// the coefficients of the trafo
        DataType _coeff[image_dim][4];

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

          // fetch the vertices of the edge
          typedef typename VertexSetType::ConstVertexReference ConstVertexReference;
          ConstVertexReference v0 = vertex_set[index_set(cell_index, 0)];
          ConstVertexReference v1 = vertex_set[index_set(cell_index, 1)];
          ConstVertexReference v2 = vertex_set[index_set(cell_index, 2)];
          ConstVertexReference v3 = vertex_set[index_set(cell_index, 3)];

          // calculate transformation coefficients
          for(Index i(0); i < Index(image_dim); ++i)
          {
            _coeff[i][0] = DataType(0.25) * DataType( v0[i] + v1[i] + v2[i] + v3[i]);
            _coeff[i][1] = DataType(0.25) * DataType(-v0[i] + v1[i] - v2[i] + v3[i]);
            _coeff[i][2] = DataType(0.25) * DataType(-v0[i] - v1[i] + v2[i] + v3[i]);
            _coeff[i][3] = DataType(0.25) * DataType( v0[i] - v1[i] - v2[i] + v3[i]);
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
        void map_point(ImagePointType& img_point, const DomainPointType& dom_point) const
        {
          for(Index i(0); i < Index(image_dim); ++i)
          {
            img_point[i] = _coeff[i][0] + _coeff[i][1] * dom_point[0] +
              (_coeff[i][2] + _coeff[i][3] * dom_point[0]) * dom_point[1];
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
          for(Index i(0); i < Index(image_dim); ++i)
          {
            jac_mat(i,0) = _coeff[i][1] + dom_point[1] * _coeff[i][3];
            jac_mat(i,1) = _coeff[i][2] + dom_point[0] * _coeff[i][3];
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
          for(Index i(0); i < Index(image_dim); ++i)
          {
            hess_ten(i,0,0) = hess_ten(i,1,1) = DataType(0);
            hess_ten(i,0,1) = hess_ten(i,1,0) = _coeff[i][3];
          }
        }
      }; // class Evaluator<Hypercube<2>,...>

      /**
       * \brief Specialisation of standard trafo evaluator for Hypercube<3> shape
       *
       * \author Stefan Wahlers
       */
      template<
        typename Trafo_,
        typename EvalPolicy_>
      class Evaluator<Trafo_, EvalPolicy_, Shape::Hypercube<3> > :
        public EvaluatorBase<Trafo_, Evaluator<Trafo_, EvalPolicy_, Shape::Hypercube<3> >, EvalPolicy_>
      {
      public:
        /// base-class typedef
        typedef EvaluatorBase<Trafo_, Evaluator, EvalPolicy_> BaseClass;
        /// shape type
        typedef Shape::Hypercube<3> ShapeType;
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

        /// can compute domain points
        static constexpr bool can_dom_point = true;
        /// can compute image points
        static constexpr bool can_img_point = true;
        /// can compute jacobian matrices
        static constexpr bool can_jac_mat = true;
        /// can compute jacobian inverse matrices if domain and image dimensions coincide
        static constexpr bool can_jac_inv = (domain_dim == image_dim);
        /// can compute jacobian determinants
        static constexpr bool can_jac_det = true;
        /// can compute hessian tensors
        static constexpr bool can_hess_ten = true;
        /// can compute inverse hessian tensors if domain and image dimensions coincide
        static constexpr bool can_hess_inv = (domain_dim == image_dim);

      protected:
        /// the coefficients of the trafo
        DataType _coeff[image_dim][8];

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

          // fetch the vertices of the edge
          typedef typename VertexSetType::ConstVertexReference ConstVertexReference;
          ConstVertexReference v0 = vertex_set[index_set(cell_index, 0)];
          ConstVertexReference v1 = vertex_set[index_set(cell_index, 1)];
          ConstVertexReference v2 = vertex_set[index_set(cell_index, 2)];
          ConstVertexReference v3 = vertex_set[index_set(cell_index, 3)];
          ConstVertexReference v4 = vertex_set[index_set(cell_index, 4)];
          ConstVertexReference v5 = vertex_set[index_set(cell_index, 5)];
          ConstVertexReference v6 = vertex_set[index_set(cell_index, 6)];
          ConstVertexReference v7 = vertex_set[index_set(cell_index, 7)];

          // calculate transformation coefficients
          // j = _coeff[i][j] for all j = 0....7
          // v = 0 + 1*x + 2*y + 3*z + x*y*4 + x*z*5 + y*z*6 + x*y*z*7
          for(Index i(0); i < Index(image_dim); ++i)
          {
            _coeff[i][0] = DataType(0.125) * DataType( + v0[i] + v1[i] + v2[i] + v3[i] + v4[i] + v5[i] + v6[i] + v7[i]);
            _coeff[i][1] = DataType(0.125) * DataType( - v0[i] + v1[i] - v2[i] + v3[i] - v4[i] + v5[i] - v6[i] + v7[i]);
            _coeff[i][2] = DataType(0.125) * DataType( - v0[i] - v1[i] + v2[i] + v3[i] - v4[i] - v5[i] + v6[i] + v7[i]);
            _coeff[i][3] = DataType(0.125) * DataType( - v0[i] - v1[i] - v2[i] - v3[i] + v4[i] + v5[i] + v6[i] + v7[i]);
            _coeff[i][4] = DataType(0.125) * DataType( + v0[i] - v1[i] - v2[i] + v3[i] + v4[i] - v5[i] - v6[i] + v7[i]);
            _coeff[i][5] = DataType(0.125) * DataType( + v0[i] - v1[i] + v2[i] - v3[i] - v4[i] + v5[i] - v6[i] + v7[i]);
            _coeff[i][6] = DataType(0.125) * DataType( + v0[i] + v1[i] - v2[i] - v3[i] - v4[i] - v5[i] + v6[i] + v7[i]);
            _coeff[i][7] = DataType(0.125) * DataType( - v0[i] + v1[i] + v2[i] - v3[i] + v4[i] - v5[i] - v6[i] + v7[i]);
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
        void map_point(ImagePointType& img_point, const DomainPointType& dom_point) const
        {
          for(Index i(0); i < Index(image_dim); ++i)
          {
            img_point[i] = _coeff[i][0]
              + dom_point[0] * (_coeff[i][1])
              + dom_point[1] * (_coeff[i][2] + dom_point[0]*_coeff[i][4])
              + dom_point[2] * (_coeff[i][3] + dom_point[0]*_coeff[i][5]
                + dom_point[1]*(_coeff[i][6] + dom_point[0]*_coeff[i][7]));
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
          for(Index i(0); i < Index(image_dim); ++i)
          {
            jac_mat(i,0) = _coeff[i][1] + dom_point[1] * _coeff[i][4] + dom_point[2] * (_coeff[i][5] + dom_point[1] * _coeff[i][7]);
            jac_mat(i,1) = _coeff[i][2] + dom_point[0] * _coeff[i][4] + dom_point[2] * (_coeff[i][6] + dom_point[0] * _coeff[i][7]);
            jac_mat(i,2) = _coeff[i][3] + dom_point[0] * _coeff[i][5] + dom_point[1] * (_coeff[i][6] + dom_point[0] * _coeff[i][7]);
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
          for(Index i(0); i < Index(image_dim); ++i)
          {
            hess_ten(i,0,0) = hess_ten(i,1,1) = hess_ten(i,2,2) = DataType(0);
            hess_ten(i,0,1) = hess_ten(i,1,0) = _coeff[i][4] + _coeff[i][7] * dom_point[2];
            hess_ten(i,0,2) = hess_ten(i,2,0) = _coeff[i][5] + _coeff[i][7] * dom_point[1];
            hess_ten(i,1,2) = hess_ten(i,2,1) = _coeff[i][6] + _coeff[i][7] * dom_point[0];
          }
        }
      }; // class Evaluator<Hypercube<3>,...>

    } // namespace Standard
  } // namespace Trafo
} // namespace FEAST

#endif // KERNEL_TRAFO_STANDARD_EVALUATOR_HPP
