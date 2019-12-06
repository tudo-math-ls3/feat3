// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_TRAFO_ISOPARAM_MAPPING_HPP
#define KERNEL_TRAFO_ISOPARAM_MAPPING_HPP 1

// includes, FEAT
#include <kernel/trafo/mapping_base.hpp>
#include <kernel/trafo/isoparam/evaluator.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/atlas/chart.hpp>

namespace FEAT
{
  namespace Trafo
  {
    /**
     * \brief Iso-Parametric Transformation namespace
     *
     * This namespace encapsulates all classes related to the implementation of the higher-order
     * "iso-parametric" transformation mappings.
     */
    namespace Isoparam
    {
      /// \cond internal
      namespace Intern
      {
        template<typename Shape_, int dim_ = Shape_::dimension-1>
        struct PartChartHelper
        {
          template<typename AVC_, typename Chart_>
          static void emplace(const Geometry::TargetSetHolder<Shape_>& tsh, AVC_& v, const Chart_& chart)
          {
            PartChartHelper<Shape_, dim_-1>::emplace(tsh, v, chart);

            const auto& trg = tsh.template get_target_set<dim_>();
            for(Index i(0); i < trg.get_num_entities(); ++i)
            {
              v[dim_].at(trg[i]) = &chart;
            }
          }
        };

        template<typename Shape_>
        struct PartChartHelper<Shape_,0>
        {
          template<typename AVC_, typename Chart_>
          static void emplace(const Geometry::TargetSetHolder<Shape_>&, AVC_&, const Chart_&)
          {
            // nothing to do here
          }
        };
      } // namespace Intern
      /// \endcond

      /**
       * \brief Standard transformation mapping class template
       *
       * This class implements the standard first-order transformation mapping for any sort of mesh.
       *
       * \attention
       * This mapping implementation is <b>highly</b> experimental.
       * Do not use it unless you really know what you are doing!
       *
       * \warning
       * It seems that the implementation for any degree > 2 is still buggy.
       * Do not use any higher degree except for debugging this implementation!
       *
       * \tparam Mesh_
       * The mesh class that this transformation is to be defined on.
       *
       * \tparam degree_
       * The desired polynomial degree of the mapping. Must be 1 <= degree <= 3.
       *
       * \todo implement Simplex evaluators (all)
       * \todo implement volume() function of all evaluators
       * \todo implement directed_width() function of all evaluators
       *
       * \author Peter Zajac
       */
      template<typename Mesh_, int degree_>
      class Mapping :
        public MappingBase<Mesh_>
      {
      public:
        static_assert(degree_ > 0, "invalid mapping degree");
        static_assert(degree_ < 4, "invalid mapping degree");

        /// base-class typedef
        typedef MappingBase<Mesh_> BaseClass;
        /// mesh type
        typedef Mesh_ MeshType;
        /// chart type
        typedef Geometry::Atlas::ChartBase<MeshType> ChartType;
        /// data type
        typedef typename MeshType::VertexSetType::CoordType CoordType;
        /// shape type
        typedef typename MeshType::ShapeType ShapeType;
        /// Shape of the facets that make up the boundary we can i.e. compute the normals for,
        /// i.e. Simplex<2> faces for a Simplex<3> mesh
        typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension - 1>::ShapeType FacetType;

        /// our shape dimension
        static constexpr int shape_dim = ShapeType::dimension;

      protected:
        // for each shape dimension, a vector of charts
        std::array<std::vector<const ChartType*>, shape_dim+1> _shape_charts;

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
          typedef Trafo::Isoparam::Evaluator<Mapping, EvalPolicy, degree_> Type;
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
          // allocate chart vectors
          for(int dim(1); dim <= shape_dim; ++dim)
          {
            _shape_charts.at(std::size_t(dim)).resize(mesh.get_num_entities(dim), nullptr);
          }
        }

        /**
         * \brief Adds a mesh-part and its associated chart to the trafo.
         *
         * \param[in] mesh_part
         * The mesh-part that is to be added.
         *
         * \param[in] chart
         * The chart that the mesh part is associated with.
         */
        void add_meshpart_chart(const Geometry::MeshPart<MeshType>& mesh_part, const Geometry::Atlas::ChartBase<MeshType>& chart)
        {
          Intern::PartChartHelper<ShapeType>::emplace(mesh_part.get_target_set_holder(), _shape_charts, chart);
        }

        /**
         * \brief Clears all added mesh-part/chart pairs.
         */
        void clear()
        {
          for(int dim(1); dim <= shape_dim; ++dim)
          {
            for(auto& x : _shape_charts.at(std::size_t(dim)))
              x = nullptr;
          }
        }

        /**
         * \brief Return the shape-charts vector.
         *
         * This function is intended for internal use by the evaluator classes only.
         */
        const std::vector<const ChartType*>& get_charts_vector(int dim) const
        {
          XASSERTM((dim >= 0) && (dim <= shape_dim), "invalid dimension");
          return _shape_charts.at(std::size_t(dim));
        }
      }; // class Mapping<...>
    } // namespace Isoparam
  } // namespace Trafo
} // namespace FEAT

#endif // KERNEL_TRAFO_ISOPARAM_MAPPING_HPP
