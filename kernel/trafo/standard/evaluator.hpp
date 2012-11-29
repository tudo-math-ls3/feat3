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

        /// trafo coefficient type
        typedef typename EvalPolicy::TrafoCoeffType CoeffType;

        /// const domain point reference
        typedef typename EvalPolicy::DomainPointConstRef DomainPointConstRef;

        /// image coordinate type
        typedef typename EvalPolicy::ImageCoordType ImageCoordType;
        /// image point reference
        typedef typename EvalPolicy::ImagePointRef ImagePointRef;

        /// jacobian matrix coefficient type
        typedef typename EvalPolicy::JacMatCoeff JacMatCoeff;
        /// jacobian matrix reference
        typedef typename EvalPolicy::JacMatRef JacMatRef;

        /// dummy enumeration
        enum
        {
          /// domain dimension
          domain_dim = EvalPolicy::domain_dim,
          /// image dimension
          image_dim = EvalPolicy::image_dim
        };

        /// capability enumeration
        enum EvaluatorCapabilities
        {
          /// can compute domain points
          can_dom_point = 1,
          /// can compute image points
          can_img_point = 1,
          /// can compute jacobian matrices
          can_jac_mat = 1,
          /// can compute jacobian inverse matrices if domain and image dimensions coincide
          can_jac_inv = (int(domain_dim) == int(image_dim)) ? 1 : 0,
          /// can compute jacobian determinants
          can_jac_det = 1
        };

      protected:
        /// the coefficients of the trafo
        CoeffType _coeff[image_dim][2];

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
          for(int i(0); i < image_dim; ++i)
          {
            _coeff[i][0] = CoeffType(0.5) * CoeffType( v0[i] + v1[i]);
            _coeff[i][1] = CoeffType(0.5) * CoeffType(-v0[i] + v1[i]);
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
        void map_point(ImagePointRef img_point, DomainPointConstRef dom_point) const
        {
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] =
              ImageCoordType(_coeff[i][0]) + ImageCoordType(_coeff[i][1]) * ImageCoordType(dom_point[0]);
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
        void calc_jac_mat(JacMatRef jac_mat, DomainPointConstRef /*dom_point*/) const
        {
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = JacMatCoeff(_coeff[i][1]);
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

        /// trafo coefficient type
        typedef typename EvalPolicy::TrafoCoeffType CoeffType;

        /// const domain point reference
        typedef typename EvalPolicy::DomainPointConstRef DomainPointConstRef;

        /// image coordinate type
        typedef typename EvalPolicy::ImageCoordType ImageCoordType;
        /// image point reference
        typedef typename EvalPolicy::ImagePointRef ImagePointRef;

        /// jacobian matrix coefficient type
        typedef typename EvalPolicy::JacMatCoeff JacMatCoeff;
        /// jacobian matrix reference
        typedef typename EvalPolicy::JacMatRef JacMatRef;

        /// dummy enumeration
        enum
        {
          /// domain dimension
          domain_dim = EvalPolicy::domain_dim,
          /// image dimension
          image_dim = EvalPolicy::image_dim
        };

        /// dummy enumeration
        enum
        {
          /// can compute domain points
          can_dom_point = 1,
          /// can compute image points
          can_img_point = 1,
          /// can compute jacobian matrices
          can_jac_mat = 1,
          /// can compute jacobian inverse matrices if domain and image dimensions coincide
          can_jac_inv = (int(domain_dim) == int(image_dim)) ? 1 : 0,
          /// can compute jacobian determinants
          can_jac_det = 1
        };

      protected:
        /// the coefficients of the trafo
        CoeffType _coeff[image_dim][4];

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
          for(int i(0); i < image_dim; ++i)
          {
            _coeff[i][0] = CoeffType(0.25) * CoeffType( v0[i] + v1[i] + v2[i] + v3[i]);
            _coeff[i][1] = CoeffType(0.25) * CoeffType(-v0[i] + v1[i] - v2[i] + v3[i]);
            _coeff[i][2] = CoeffType(0.25) * CoeffType(-v0[i] - v1[i] + v2[i] + v3[i]);
            _coeff[i][3] = CoeffType(0.25) * CoeffType( v0[i] - v1[i] - v2[i] + v3[i]);
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
        void map_point(ImagePointRef img_point, DomainPointConstRef dom_point) const
        {
          for(int i(0); i < image_dim; ++i)
          {
            img_point[i] =
              ImageCoordType(_coeff[i][0]) + ImageCoordType(_coeff[i][1]) * ImageCoordType(dom_point[0]) +
              (ImageCoordType(_coeff[i][2]) + ImageCoordType(_coeff[i][3]) * ImageCoordType(dom_point[0])) *
              ImageCoordType(dom_point[1]);
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
        void calc_jac_mat(JacMatRef jac_mat, DomainPointConstRef dom_point) const
        {
          for(int i(0); i < image_dim; ++i)
          {
            jac_mat(i,0) = JacMatCoeff(_coeff[i][1]) + JacMatCoeff(dom_point[1]) * JacMatCoeff(_coeff[i][3]);
            jac_mat(i,1) = JacMatCoeff(_coeff[i][2]) + JacMatCoeff(dom_point[0]) * JacMatCoeff(_coeff[i][3]);
          }
        }
      }; // class Evaluator<Hypercube<2>,...>
    } // namespace Standard
  } // namespace Trafo
} // namespace FEAST

#endif // KERNEL_TRAFO_STANDARD_EVALUATOR_HPP
