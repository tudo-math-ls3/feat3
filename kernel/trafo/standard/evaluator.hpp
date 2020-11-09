// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_TRAFO_STANDARD_EVALUATOR_HPP
#define KERNEL_TRAFO_STANDARD_EVALUATOR_HPP 1

// includes, FEAT
#include <kernel/trafo/evaluator_base.hpp>

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
      class Evaluator DOXY({});

      /* ************************************************************************************* */

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


      /* ************************************************************************************* */

      /**
       * \brief Specialization of standard trafo evaluator for Simplex<1> shape
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

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten |
          (domain_dim == image_dim ? (TrafoTags::jac_inv | TrafoTags::hess_inv) : TrafoTags::none);

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
          typedef typename VertexSetType::VertexType VertexType;
          const VertexType& v0 = vertex_set[index_set(cell_index, 0)];
          const VertexType& v1 = vertex_set[index_set(cell_index, 1)];

          // calculate transformation coefficients
          for(int i(0); i < image_dim; ++i)
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
          for(int i(0); i < image_dim; ++i)
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
          for(int i(0); i < image_dim; ++i)
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
          hess_ten.format();
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \returns
         * The volume of the current cell.
         */
        DataType volume() const
        {
          DataType v = DataType(0);
          for(int i(0); i < image_dim; ++i)
            v += Math::sqr(_coeff[i][1]);
          return Math::sqrt(v);
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
        DataType width_directed(const ImagePointType& DOXY(ray)) const
        {
          // in 1D, the width is always equal to the volume
          return volume();
        }
      }; // class Evaluator<Simplex<1>,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialization of standard trafo evaluator for Simplex<2> shape
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

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten |
          (domain_dim == image_dim ? (TrafoTags::jac_inv | TrafoTags::hess_inv) : TrafoTags::none);

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
          typedef typename VertexSetType::VertexType VertexType;
          const VertexType& v0 = vertex_set[index_set(cell_index, 0)];
          const VertexType& v1 = vertex_set[index_set(cell_index, 1)];
          const VertexType& v2 = vertex_set[index_set(cell_index, 2)];

          // calculate transformation coefficients
          for(int i(0); i < image_dim; ++i)
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
          for(int i(0); i < image_dim; ++i)
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
          for(int i(0); i < image_dim; ++i)
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
          hess_ten.format();
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \returns
         * The volume of the current cell.
         */
        DataType volume() const
        {
          JacobianMatrixType jac_mat;
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = _coeff[i][1];
            jac_mat(i,1) = _coeff[i][2];
          }
          return jac_mat.vol() / DataType(2);
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
          // This one is a little tricky:
          // We follow a similar approach as on quadrilaterals here, i.e. we first map the ray
          // vector onto the reference cells by multiplying it to the inverse jacobian matrix.
          // However, our reference cell is the "right triangle" spanned by the three points
          // (0,0) (1,0) (0,1) and therefore its edges have different lengths, which would
          // result in a "skew" width. To circumvent this, we will afterwards map the ray
          // vector onto a "equilateral unit triangle" (i.e. the triangle where all edges
          // have length = 1) and we compute the norm of that mapped vector then.
          JacobianMatrixType jac_mat;
          JacobianInverseType jac_inv;
          DomainPointType ref_ray;

          // compute jacobian matrix
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = _coeff[i][1];
            jac_mat(i,1) = _coeff[i][2];
          }

          // invert jacobian matrix and multiply by ray vector
          jac_inv.set_inverse(jac_mat);
          ref_ray.set_mat_vec_mult(jac_inv, ray);

          // We have mapped the ray onto our reference triangle. Now let's map that one
          // onto a equilateral unit triangle, which can be performed by multiplying
          // the following matrix by the reference ray vector:
          //
          //  eut_ray := [ -1       +1      ] * ref_ray
          //             [ sqrt(3)  sqrt(3) ]
          //
          // Finally, we have to compute the norm of that vector, which is explicitly
          // written down in the denominator of the following expression:
          return DataType(2) / Math::sqrt(Math::sqr(ref_ray[1]-ref_ray[0]) + DataType(3) * Math::sqr(ref_ray[0]+ref_ray[1]));
        }
      }; // class Evaluator<Simplex<2>,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialization of standard trafo evaluator for Simplex<3> shape
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

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten |
          (domain_dim == image_dim ? (TrafoTags::jac_inv | TrafoTags::hess_inv) : TrafoTags::none);

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
          typedef typename VertexSetType::VertexType VertexType;
          const VertexType& v0 = vertex_set[index_set(cell_index, 0)];
          const VertexType& v1 = vertex_set[index_set(cell_index, 1)];
          const VertexType& v2 = vertex_set[index_set(cell_index, 2)];
          const VertexType& v3 = vertex_set[index_set(cell_index, 3)];

          // calculate transformation coefficients
          for(int i(0); i < image_dim; ++i)
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
          for(int i(0); i < image_dim; ++i)
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
          for(int i(0); i < image_dim; ++i)
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
          hess_ten.format();
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \returns
         * The volume of the current cell.
         */
        DataType volume() const
        {
          JacobianMatrixType jac_mat;
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = _coeff[i][1];
            jac_mat(i,1) = _coeff[i][2];
            jac_mat(i,2) = _coeff[i][3];
          }
          return jac_mat.vol() / DataType(6);
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
          // We follow the same approach as for triangles here:
          // First, map the ray onto the reference element and then
          // map it to a regular tetrahedron to compute its norm.
          JacobianMatrixType jac_mat;
          JacobianInverseType jac_inv;
          DomainPointType ref_ray;

          // compute jacobian matrix
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = _coeff[i][1];
            jac_mat(i,1) = _coeff[i][2];
            jac_mat(i,2) = _coeff[i][3];
          }

          // invert jacobian matrix and multiply by ray vector
          jac_inv.set_inverse(jac_mat);
          ref_ray.set_mat_vec_mult(jac_inv, ray);

          // The transformation from the reference tetrahedron to the regular unit tetrahedron
          // can be performed by applying the following matrix-vector product:
          //
          //             [ 2  1         1      ]
          //  eut_ray := [ 0  1        -1      ] * ref_ray
          //             [ 0  sqrt(2)  sqrt(2) ]
          //
          // Finally, we have to compute the norm of that vector, which is explicitly
          // written down in the denominator of the following expression:
          return DataType(2) / Math::sqrt(Math::sqr(DataType(2) * ref_ray[0] + ref_ray[1] + ref_ray[2]) +
            Math::sqr(ref_ray[1] - ref_ray[2]) + DataType(2) * Math::sqr(ref_ray[1] + ref_ray[2]));
        }
      }; // class Evaluator<Simplex<3>,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialization of standard trafo evaluator for Hypercube<1> shape
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

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten |
          (domain_dim == image_dim ? (TrafoTags::jac_inv | TrafoTags::hess_inv) : TrafoTags::none);

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
          typedef typename VertexSetType::VertexType VertexType;
          const VertexType& v0 = vertex_set[index_set(cell_index, 0)];
          const VertexType& v1 = vertex_set[index_set(cell_index, 1)];

          // calculate transformation coefficients
          for(int i(0); i < image_dim; ++i)
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
          for(int i(0); i < image_dim; ++i)
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
          for(int i(0); i < image_dim; ++i)
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
          hess_ten.format();
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \returns
         * The volume of the current cell.
         */
        DataType volume() const
        {
          DataType v = DataType(0);
          for(int i(0); i < image_dim; ++i)
            v += Math::sqr(_coeff[i][1]);
          return DataType(2) * Math::sqrt(v);
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
        DataType width_directed(const ImagePointType& DOXY(ray)) const
        {
          // in 1D, the width is always equal to the volume
          return volume();
        }
      }; // class Evaluator<Hypercube<1>,...>

      /* ************************************************************************************* */

      /**
       * \brief Specialization of standard trafo evaluator for Hypercube<2> shape
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

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten |
          (domain_dim == image_dim ? (TrafoTags::jac_inv | TrafoTags::hess_inv) : TrafoTags::none);

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
          typedef typename VertexSetType::VertexType VertexType;
          const VertexType& v0 = vertex_set[index_set(cell_index, 0)];
          const VertexType& v1 = vertex_set[index_set(cell_index, 1)];
          const VertexType& v2 = vertex_set[index_set(cell_index, 2)];
          const VertexType& v3 = vertex_set[index_set(cell_index, 3)];

          // calculate transformation coefficients
          for(int i(0); i < image_dim; ++i)
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
          for(int i(0); i < image_dim; ++i)
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
          for(int i(0); i < image_dim; ++i)
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
          for(int i(0); i < image_dim; ++i)
          {
            hess_ten(i,0,0) = hess_ten(i,1,1) = DataType(0);
            hess_ten(i,0,1) = hess_ten(i,1,0) = _coeff[i][3];
          }
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \returns
         * The volume of the current cell.
         */
        DataType volume() const
        {
          // According to Varignon's theorem, the area/volume of a quadrilateral is
          // equal to twice the area of the dual parallelogram of the quadrilateral,
          // which is spanned by the four edge midpoints of the quadrilateral.
          // The Jacobian matrix of this transformation evaluated at the barycentre
          // of the reference element spans a parallelogram, which intersects with
          // our original quadrilateral in the edge midpoints and therefore (again
          // using Varignon's theorem) has the same area as the original quadrilateral.
          // Now the area of the "Jacobian parallelogram" is equal to four times
          // the determinant of its Jacobian determinant, which finally gives us a
          // formula for our quadrilateral area: 4*det(Jacobian(0,0)).
          // Note that this is not a lousy approximation, but a true identity :)

          // compute jacobian matrix at barycentre
          JacobianMatrixType jac_mat;
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = _coeff[i][1];
            jac_mat(i,1) = _coeff[i][2];
          }

          // return scaled volume
          return DataType(4) * jac_mat.vol();
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
          JacobianMatrixType jac_mat;
          JacobianInverseType jac_inv;
          DomainPointType ref_ray;

          // compute jacobian matrix at barycentre
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = _coeff[i][1];
            jac_mat(i,1) = _coeff[i][2];
          }

          // invert jacobian matrix and multiply by ray vector
          jac_inv.set_inverse(jac_mat);
          ref_ray.set_mat_vec_mult(jac_inv, ray);

          // return scaled inverse ray norm
          return DataType(2) / ref_ray.norm_euclid();
        }
      }; // class Evaluator<Hypercube<2>,...>

      /**
       * \brief Specialization of standard trafo evaluator for Hypercube<3> shape
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

        /// we can compute domain and image points,
        /// jacobian matrices and determinants as well as hessian tensors.
        static constexpr TrafoTags eval_caps =
          TrafoTags::dom_point | TrafoTags::img_point |
          TrafoTags::jac_mat | TrafoTags::jac_det | TrafoTags::hess_ten |
          (domain_dim == image_dim ? (TrafoTags::jac_inv | TrafoTags::hess_inv) : TrafoTags::none);

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
          typedef typename VertexSetType::VertexType VertexType;
          const VertexType& v0 = vertex_set[index_set(cell_index, 0)];
          const VertexType& v1 = vertex_set[index_set(cell_index, 1)];
          const VertexType& v2 = vertex_set[index_set(cell_index, 2)];
          const VertexType& v3 = vertex_set[index_set(cell_index, 3)];
          const VertexType& v4 = vertex_set[index_set(cell_index, 4)];
          const VertexType& v5 = vertex_set[index_set(cell_index, 5)];
          const VertexType& v6 = vertex_set[index_set(cell_index, 6)];
          const VertexType& v7 = vertex_set[index_set(cell_index, 7)];

          // calculate transformation coefficients
          // j = _coeff[i][j] for all j = 0....7
          // v = 0 + 1*x + 2*y + 3*z + x*y*4 + x*z*5 + y*z*6 + x*y*z*7
          for(int i(0); i < image_dim; ++i)
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
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] = _coeff[i][0]
              + dom_point[0] * (_coeff[i][1])
              + dom_point[1] * (_coeff[i][2] + dom_point[0]*_coeff[i][4])
              + dom_point[2] * (_coeff[i][3] + dom_point[0]*_coeff[i][5]
              + dom_point[1] * (_coeff[i][6] + dom_point[0]*_coeff[i][7]));
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
          for(int i(0); i < image_dim; ++i)
          {
            hess_ten(i,0,0) = hess_ten(i,1,1) = hess_ten(i,2,2) = DataType(0);
            hess_ten(i,0,1) = hess_ten(i,1,0) = _coeff[i][4] + _coeff[i][7] * dom_point[2];
            hess_ten(i,0,2) = hess_ten(i,2,0) = _coeff[i][5] + _coeff[i][7] * dom_point[1];
            hess_ten(i,1,2) = hess_ten(i,2,1) = _coeff[i][6] + _coeff[i][7] * dom_point[0];
          }
        }

        /**
         * \brief Computes and returns the volume of the current cell.
         *
         * \returns
         * The volume of the current cell.
         */
        DataType volume() const
        {
          // In contrast to 2D, it is not sufficient to evaluate the Jacobian determinant
          // in the barycentre of the cell to compute the cell's volume, as this will give
          // very inaccurate results for non-parallelepiped cells. Instead, we have to use
          // the 2x2x2 Gauss-Legendre cubature rule to integrate the volume of the cell,
          // which is hard-coded in the code below.

          // compute 1D Gauss-Legendre root
          const DataType cx = DataType(1) / Math::sqrt(DataType(3));

          JacobianMatrixType jac_mat;
          DomainPointType cub_pt;
          DataType vol = DataType(0);

          // loop over all 8 cubature points
          for(int i(0); i < 8; ++i)
          {
            // set cubature point coords by magic bitshifts
            for(int j(0); j < 3; ++j)
              cub_pt[j] = DataType((((i >> j) & 1) << 1) - 1) * cx;

            // compute jacobian matrix and add its volume
            calc_jac_mat(jac_mat, cub_pt);
            vol += jac_mat.vol();
          }

          return vol;
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
          JacobianMatrixType jac_mat;
          JacobianInverseType jac_inv;
          DomainPointType ref_ray;

          // compute jacobian matrix at barycentre
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = _coeff[i][1];
            jac_mat(i,1) = _coeff[i][2];
            jac_mat(i,2) = _coeff[i][3];
          }

          // invert jacobian matrix and multiply by ray vector
          jac_inv.set_inverse(jac_mat);
          ref_ray.set_mat_vec_mult(jac_inv, ray);

          // return scaled inverse ray norm
          return DataType(2) / ref_ray.norm_euclid();
        }
      }; // class Evaluator<Hypercube<3>,...>
    } // namespace Standard
  } // namespace Trafo
} // namespace FEAT

#endif // KERNEL_TRAFO_STANDARD_EVALUATOR_HPP
