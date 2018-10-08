#pragma once
#ifndef KERNEL_TRAFO_STANDARD_MAPPING_HPP
#define KERNEL_TRAFO_STANDARD_MAPPING_HPP 1

// includes, FEAT
#include <kernel/trafo/mapping_base.hpp>
#include <kernel/trafo/standard/evaluator.hpp>

namespace FEAT
{
  namespace Trafo
  {
    /**
     * \brief Standard Transformation namespace
     *
     * This namespace encapsulates all classes related to the implementation of the standard first-order
     * (i.e. P1/Q1) transformation mapping.
     */
    namespace Standard
    {
      /**
       * \brief Standard transformation mapping class template
       *
       * This class implements the standard first-order transformation mapping for any sort of mesh.
       *
       * \tparam Mesh_
       * The mesh class that this transformation is to be defined on.
       *
       * \author Peter Zajac
       */
      template<typename Mesh_>
      class Mapping :
        public MappingBase<Mesh_>
      {
      public:
        /// base-class typedef
        typedef MappingBase<Mesh_> BaseClass;
        /// mesh type
        typedef Mesh_ MeshType;
        /// data type
        typedef typename MeshType::VertexSetType::CoordType CoordType;
        /// shape type
        typedef typename MeshType::ShapeType ShapeType;
        /// Shape of the facets that make up the boundary we can i.e. compute the normals for,
        /// i.e. Simplex<2> faces for a Simplex<3> mesh
        typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension - 1>::ShapeType FacetType;
        /// Number of coefficients of the reference cell transformation is the number of vertices
        static constexpr int num_coeff = Shape::FaceTraits<ShapeType, 0>::count;
        /// Number of coefficients of the reference cell transformation for faces is the number of vertices
        static constexpr int num_coeff_facet = Shape::FaceTraits<FacetType, 0>::count;

      public:
        /** \copydoc MappingBase::Evaluator */
        template<
          typename Shape_ = ShapeType,
          typename CoordType_ = Real>
        class Evaluator
        {
        private:
          /// evaluation policy
          typedef Trafo::StandardEvalPolicy<Shape_, CoordType_, MeshType::world_dim> EvalPolicy;

        public:
          /// evaluator type
          typedef Trafo::Standard::Evaluator<Mapping, EvalPolicy> Type;
        };

      public:
        /**
         * \brief Constructor
         *
         * \param[in] mesh
         * A reference to the mesh that this trafo mapping is to be defined on.
         */
        explicit Mapping(MeshType& mesh) :
          BaseClass(mesh)
        {
        }

        /**
         * \brief Computes the volume of one cell
         *
         * \tparam ShapeType_
         * Shape type of the underlying mesh.
         *
         * \tparam CoordType_
         * Precision, by default the same precision as for the coordinates.
         *
         * \param[in] cell
         * Index for which the volume is computed.
         *
         * \returns The volume of cell.
         */
        template<typename ShapeType_ = ShapeType, typename CoordType_ = CoordType>
        CoordType_ compute_vol(const Index cell) const
        {
          typename Evaluator<ShapeType_, CoordType_>::Type evaluator(*this);

          evaluator.prepare(cell);
          CoordType_ vol = evaluator.volume();
          evaluator.finish();
          return vol;
        }

        /**
         * \brief Computes the volume of the whole mesh given by this trafo
         *
         * \returns The shape_dim dimensional volume
         */
        CoordType compute_vol()
        {
          typename Evaluator<ShapeType, CoordType>::Type evaluator(*this);

          CoordType vol(0);

          Index num_cells = this->get_mesh().get_num_entities(ShapeType::dimension);
          for(Index cell(0); cell < num_cells; ++cell)
          {
            evaluator.prepare(cell);
            vol += evaluator.volume();
            evaluator.finish();
          }

          return vol;
        }

        /**
         * \brief Computes all oriented normals on one face
         *
         * \param[in] k
         * Number of the face to compute the normals at
         *
         * \returns
         * A matrix containing the oriented normals in all Lagrange points row-wise
         *
         */
        Tiny::Matrix<CoordType, num_coeff_facet, MeshType::world_dim> compute_oriented_normals(Index k) const
        {
          // Evaluator on Faces
          typedef typename Evaluator<FacetType, CoordType>::Type FaceEvaluator;
          // Type of the Jacobian it returns
          typedef typename FaceEvaluator::JacobianMatrixType JacobianMatrixType;

          // Evaluator for the trafo
          FaceEvaluator trafo_eval(*this);
          trafo_eval.prepare(k);

          // This will be the coordinates of the Lagrange points on the reference cell the face is parametrised over
          Tiny::Matrix<CoordType, num_coeff_facet, ShapeType::dimension - 1> xloc(CoordType(0));
          for(int i(0); i < num_coeff_facet; ++i)
          {
            for(int d(0); d < ShapeType::dimension - 1 ; ++d)
              xloc(i,d) = Shape::ReferenceCell<FacetType>::template vertex<CoordType>(i,d);
          }

          // The Jacobian matrix of the transformation
          JacobianMatrixType jac_mat;
          // This will hold the local tangential surface coordinate system at one Lagrange point of the local face
          Tiny::Matrix<CoordType, MeshType::world_dim-1, MeshType::world_dim> tau(CoordType(0));

          // This will hold the oriented unit normal vector for all Lagrange points of the face
          Tiny::Matrix<CoordType, num_coeff_facet, MeshType::world_dim> nu(CoordType(0));
          for(int i(0); i < num_coeff_facet; ++i)
          {
            tau.format();

            // Compute the tangential vectors at the Lagrange points: They are the gradients wrt. the parametrisation
            // of the transformation evaluated in the Lagrange points
            trafo_eval.calc_jac_mat(jac_mat, xloc[i]);

            // Now write the cross product of the tangentials into nu[i]
            nu[i] = Tiny::orthogonal(jac_mat);
            // Normalise nu
            nu[i].normalise();
          }

          return nu;

        } // compute_oriented_normals

      }; // class Mapping<...>
    } // namespace Standard
  } // namespace Trafo
} // namespace FEAT

#endif // KERNEL_TRAFO_STANDARD_MAPPING_HPP
