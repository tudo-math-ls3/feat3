#pragma once
#ifndef KERNEL_TRAFO_STANDARD_MAPPING_HPP
#define KERNEL_TRAFO_STANDARD_MAPPING_HPP 1

// includes, FEAST
#include <kernel/trafo/mapping_base.hpp>
#include <kernel/trafo/standard/evaluator.hpp>
#include <kernel/trafo/standard/volume.hpp>

namespace FEAST
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

        /** \copydoc MappingBase::Evaluator */
        template<
          typename Shape_ = ShapeType,
          typename DataType_ = Real>
        class Evaluator
        {
        private:
          /// evaluation policy
          typedef Trafo::StandardEvalPolicy<Shape_, DataType_, MeshType::world_dim> EvalPolicy;

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
         * \tparam DataType_
         * Precision, by default the same precision as for the coordinates.
         *
         * \param[in] cell
         * Index for which the volume is computed.
         *
         * \returns The volume of cell.
         *
         * \author Jordi Paul
         */
        template<typename ShapeType_, typename DataType_ = CoordType>
          DataType_ compute_vol(const Index cell)
          {
            // Extract the transformation, the underlying mesh's index and vertex sets and stuff 'em into the
            // CellVolumeEvaluator who does the actual work.
            return CellVolumeEvaluator<ShapeType_>::template compute_vol<DataType_>(
                *this,
                this->get_mesh().template get_index_set<ShapeType_::dimension,0>(),
                this->get_mesh().get_vertex_set(),
                cell);
          }
      }; // class Mapping<...>
    } // namespace Standard
  } // namespace Trafo
} // namespace FEAST

#endif // KERNEL_TRAFO_STANDARD_MAPPING_HPP
