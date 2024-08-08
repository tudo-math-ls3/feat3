// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_TRAFO_STANDARD_MAPPING_HPP
#define KERNEL_TRAFO_STANDARD_MAPPING_HPP 1

// includes, FEAT
#include <kernel/trafo/mapping_base.hpp>
#include <kernel/trafo/standard/evaluator.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/atlas/chart.hpp>

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
         * A \resident reference to the mesh that this trafo mapping is to be defined on.
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
         * \brief Adds a mesh-part and its associated chart to the trafo.
         *
         * \param[in] mesh_part
         * The \transient mesh-part that is to be added.
         *
         * \param[in] chart
         * The \transient chart that the mesh part is associated with.
         */
        void add_meshpart_chart(const Geometry::MeshPart<MeshType>& DOXY(mesh_part), const Geometry::Atlas::ChartBase<MeshType>& DOXY(chart))
        {
          // nothing to do here; the only purpose of this function is to supply
          // the same interface as the Trafo::Isoparam::Mapping class
        }
      }; // class Mapping<...>
    } // namespace Standard
  } // namespace Trafo
} // namespace FEAT

#endif // KERNEL_TRAFO_STANDARD_MAPPING_HPP
