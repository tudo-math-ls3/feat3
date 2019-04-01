// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_CHART_HPP
#define KERNEL_GEOMETRY_ATLAS_CHART_HPP 1

// includes, FEAT
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/xml_scanner.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Atlas namespace
     */
    namespace Atlas
    {
      /**
       * \brief Chart base-class
       *
       * \tparam Mesh_
       * The type of the mesh that is to be parameterised.
       *
       * \author Peter Zajac
       */
      template<typename Mesh_>
      class ChartBase
      {
      public:
        /// our mesh type
        typedef Mesh_ MeshType;
        /// our mesh part type
        typedef MeshPart<Mesh_> PartType;
        /// our vertex set type
        typedef typename MeshType::VertexSetType VertexSetType;
        /// Type of a single vertex
        typedef typename VertexSetType::VertexType WorldPoint;
        /// out coordinate type
        typedef typename VertexSetType::CoordType CoordType;

      public:
        /// virtual DTOR
        virtual ~ChartBase() {}

        /// \returns The size of dynamically allocated memory in bytes.
        virtual std::size_t bytes() const
        {
          return std::size_t(0);
        }

        /**
         * \brief Specifies whether the chart can perform explicit projection.
         *
         * \returns
         * \c true if explicit projection is possible, otherwise \c false.
         */
        virtual bool can_explicit() const = 0;

        /**
         * \brief Specifies whether the chart can perform implicit projection.
         *
         * \returns
         * \c true if implicit projection is possible, otherwise \c false.
         */
        virtual bool can_implicit() const = 0;

        /**
         * \brief Adapts a mesh using this chart.
         *
         * \param[inout] mesh
         * The mesh that is to be adapted.
         *
         * \param[in] meshpart
         * The mesh part that describes the part to adapt.
         */
        virtual void adapt(MeshType& mesh, const PartType& meshpart) const = 0;

        /**
         * \brief Adapts a mesh part using this chart.
         *
         * \param[inout] mesh
         * The mesh part that is to be adapted.
         *
         * \param[in] meshpart
         * The mesh part that describes the part to adapt.
         */
        virtual void adapt(PartType& mesh, const PartType& meshpart) const = 0;

        /**
         * \brief Applies a "proper rigid" transformation onto the chart.
         *
         * Let \e v denote the \p origin world point, \e w the \p offset world point and \e R
         * the rotation matrix corresponding to the \p angles, then this function applies the
         * following transformation for any point \e x of the chart:
         *
         *   \f[ x \mapsto w + R\cdot (x - v) \f]
         *
         * \param[in] origin
         * The origin of the transformation. This is subtracted from any point before applying the
         * rotation.
         *
         * \param[in] angles
         * The angles of the rotation matrix.
         * - 2D: the rotation angle in radians is stored as:
         *   - angles(0): rotation angle
         *   - angles(1): -ignored-
         * - 3D: the rotation angles in radians stored as:
         *   - angles(0): yaw angle
         *   - angles(1): pitch angle
         *   - angles(2): roll angle
         *
         * \param[in] offset
         * The offset of the transformation. This is added to any point after applying the rotation.
         */
        virtual void transform(const WorldPoint& origin, const WorldPoint& angles, const WorldPoint& offset) = 0;

        /**
         * \brief Maps a parameter to a world point
         *
         * \param[in] param
         * The parameter to be mapped.
         *
         * \returns
         * The mapped world point.
         */
        virtual WorldPoint map(const WorldPoint& param) const = 0;

        /**
         * \brief Projects a point onto the chart
         *
         * \param[in] point
         * The point to be projected.
         *
         * \returns
         * The projected point.
         */
        virtual WorldPoint project(const WorldPoint& point) const = 0;

        /**
         * \brief Computes the distance of a point to this chart
         *
         * \param[in] point
         * The world point to compute the distance for
         *
         * \returns The distance to this chart
         */
        virtual CoordType dist(const WorldPoint& point) const = 0;

        /**
         * \brief Computes the distance of a point to this chart
         *
         * \param[in] point
         * The world point to compute the distance for
         *
         * \param[in] grad_dist
         * The gradient of the distance function in point.
         *
         * \returns The distance to this chart
         */
        virtual CoordType dist(const WorldPoint& point, WorldPoint& grad_dist) const = 0;

        /**
         * \brief Computes the signed distance of a point to this chart
         *
         * \param[in] point
         * The world point to compute the distance for
         *
         * \returns The signed distance to this chart
         */
        virtual CoordType signed_dist(const WorldPoint& point) const = 0;

        /**
         * \brief Computes the signed distance of a point to this chart
         *
         * \param[in] point
         * The world point to compute the distance for
         *
         * \param[in] grad_dist
         * The gradient of the distance function in point.
         *
         * \returns The signed distance to this chart
         */
        virtual CoordType signed_dist(const WorldPoint& point, WorldPoint& grad_dist) const = 0;

        /**
         * \brief Writes the type as String
         *
         * Needed for writing this to mesh files.
         *
         * \returns The class name as String.
         */
        virtual String get_type() const = 0;

        /**
         * \brief Writes the Chart into a stream in XML format.
         *
         * \param[in,out] os
         * The output stream to write into.
         *
         * \param[in] sindent
         * The indentation string.
         */
        virtual void write(std::ostream& os, const String& sindent) const = 0;
      }; // class ChartBase<...>

      /// \cond internal
      namespace Intern
      {
        template<bool enable_>
        struct ImplicitChartHelper
        {
          template<typename CT_, typename MT_, typename PT_>
          static bool adapt(const CT_&, MT_&, const PT_&)
          {
            return false;
          }

          template<typename CT_, typename WP_>
          static bool project(const CT_&, WP_&)
          {
            return false;
          }
        };

        template<>
        struct ImplicitChartHelper<true>
        {
          template<typename CT_, typename MT_, typename PT_>
          static bool adapt(const CT_& chart, MT_& mesh, const PT_& part)
          {
            // First of all, check whether the chart can really perform
            // implicit adaption
            if(!chart.can_implicit())
              return false;

            // Try to project the whole meshpart
            chart.project_meshpart(mesh, part);

            // okay
            return true;
          }

          template<typename CT_, typename WP_>
          static bool project(const CT_& chart, WP_& wp)
          {
            chart.project_point(wp);
            return true;
          }
        };

        template<bool enable_>
        struct ExplicitChartHelper
        {
          template<typename CT_, typename MT_, typename PT_>
          static bool adapt(const CT_&, MT_&, const PT_&)
          {
            return false;
          }

          template<typename CT_, typename WP_, typename PP_>
          static bool map(const CT_&, WP_&, const PP_&)
          {
            return false;
          }
        };

        template<>
        struct ExplicitChartHelper<true>
        {
          template<typename CT_, typename MT_, typename PT_>
          static bool adapt(const CT_& chart, MT_& mesh, const PT_& part)
          {
            // a world point in the mesh
            typedef typename CT_::WorldPoint WorldPoint;
            // a parameter point in the part
            typedef typename CT_::ParamPoint ParamPoint;

            // vertex set type of our mesh
            typedef typename MT_::VertexSetType VertexSetType;

            // attribute type of our mesh part
            typedef typename PT_::AttributeSetType AttributeSetType;

            // First of all, check whether the chart can really perform
            // explicit adaption
            if(!chart.can_explicit())
              return false;

            // Try to fetch the parametrisation attribute.
            const AttributeSetType* attrib = part.find_attribute("param");
            if(attrib == nullptr)
              return false;

            // We have the attribute; check whether it matches our chart
            int attrib_dim = attrib->get_dimension();
            XASSERTM(attrib_dim == CT_::param_dim, "Invalid chart attribute dimension");

            // Get the vertex set of the mesh
            VertexSetType& vtx = mesh.get_vertex_set();

            // Get the vertex target set of the part
            const TargetSet& vidx = part.template get_target_set<0>();

            // loop over all vertices in the mesh part
            Index num_vtx = vidx.get_num_entities();
            for(Index i(0); i < num_vtx; ++i)
            {
              // apply the chart's map function
              chart.map_param(
                 reinterpret_cast<      WorldPoint&>(vtx[vidx[i]]),
                *reinterpret_cast<const ParamPoint*>(attrib->raw_at(i))
                );
            }

            // okay
            return true;
          }

          template<typename CT_, typename WP_, typename PP_>
          static bool map(const CT_& chart, WP_& wp, const PP_& pp)
          {
            chart.map_param(wp, pp);
            return true;
          }
        };
      } // namespace Intern
      /// \endcond

      /**
       * \brief Chart CRTP base-class template
       *
       * This class template acts as a CRTP base-class for the actual chart implementations.
       *
       * \tparam Derived_
       * The derived class.
       *
       * \tparam Mesh_
       * The mesh class to be parameterised by this chart class.
       *
       * \tparam Traits_
       * A traits class that specifies the chart's properties.
       *
       * \author Peter Zajac
       */
      template<typename Derived_, typename Mesh_, typename Traits_>
      class ChartCRTP :
        public ChartBase<Mesh_>
      {
      public:
        /// base-class type
        typedef ChartBase<Mesh_> BaseClass;
        /// traits type
        typedef Traits_ TraitsType;
        /// mesh type
        typedef typename BaseClass::MeshType MeshType;
        /// mesh-part type
        typedef typename BaseClass::PartType PartType;
        /// attribute type of our mesh-part
        typedef typename PartType::AttributeSetType AttributeSetType;
        /// vertex set type of our mesh
        typedef typename MeshType::VertexSetType VertexSetType;
        /// coordinate type
        typedef typename VertexSetType::CoordType CoordType;

        /// specifies whether this chart is explicit
        static constexpr bool is_explicit = TraitsType::is_explicit;
        /// specifies whether this chart is implicit
        static constexpr bool is_implicit = TraitsType::is_implicit;

        /// the world dimension of this chart
        static constexpr int world_dim = TraitsType::world_dim;
        /// the parameter dimension of this chart
        static constexpr int param_dim = TraitsType::param_dim;

        /// our world point type
        //typedef Tiny::Vector<CoordType, world_dim> WorldPoint;
        typedef typename BaseClass::WorldPoint WorldPoint;
        /// out parameter type
        typedef Tiny::Vector<CoordType, param_dim> ParamPoint;

      protected:
        /**
         * \brief Casts \c this to its true type
         *
         * \returns A Derived_ reference to \c this
         */
        Derived_& cast() {return static_cast<Derived_&>(*this);}

        /// \copydoc cast()
        const Derived_& cast() const {return static_cast<const Derived_&>(*this);}

      public:
        /// \copydoc BaseClass::can_explicit()
        virtual bool can_explicit() const override
        {
          return is_explicit;
        }

        /// \copydoc BaseClass::can_implicit()
        virtual bool can_implicit() const override
        {
          return is_implicit;
        }

        /**
         * \brief Adapts a whole MeshPart
         *
         * \param[in] mesh
         * Mesh to be adapted
         *
         * \param[in] part
         * MeshPart identifying the region to be adapted
         */
        virtual void adapt(MeshType& mesh, const PartType& part) const override
        {
          // ensure that the mesh world dimension is compatible
          XASSERTM(MeshType::world_dim == world_dim, "Mesh/Chart world dimension mismatch");

          // Try to adapt explicity
          if(Intern::ExplicitChartHelper<is_explicit>::adapt(cast(), mesh, part))
            return;

          // Try to adapt implicitly
          if(Intern::ImplicitChartHelper<is_implicit>::adapt(cast(), mesh, part))
            return;

          // If we come out here, we have no way of adaption...
          throw InternalError(__func__,__FILE__,__LINE__,"No adaption possible!");
        }

        /**
         * \brief Adapts a whole MeshPart referring to another MeshPart
         *
         * \param[in] parent_meshpart
         * MeshPart to be adapted
         *
         * \param[in] meshpart
         * MeshPart identifying the region to be adapted
         *
         * \todo: Implement this
         *
         * There is currently no code that uses MeshParts referring to other MeshParts instead of a RootMesh.
         */
        virtual void adapt(PartType& DOXY(parent_meshpart), const PartType& DOXY(meshpart)) const override
        {
          throw InternalError(__func__,__FILE__,__LINE__,"Adaption of MeshPart not possible yet");
        }

        /// \copydoc BaseClass::map()
        virtual WorldPoint map(const WorldPoint& param) const override
        {
          ASSERTM(is_explicit, "cannot map point: Chart is not explicit");
          WorldPoint wp;
          ParamPoint pp(param.template size_cast<param_dim>());
          Intern::ExplicitChartHelper<is_explicit>::map(cast(), wp, pp);
          return wp;
        }

        /// \copydoc BaseClass::project()
        virtual WorldPoint project(const WorldPoint& point) const override
        {
          ASSERTM(is_implicit, "cannot project point: Chart is not implicit");
          WorldPoint wp(point);
          Intern::ImplicitChartHelper<is_implicit>::project(cast(), wp);
          return wp;
        }

        /// \copydoc BaseClass::dist()
        virtual CoordType dist(const WorldPoint& point) const override
        {
          return (this->cast()).compute_dist(point);
        }

        /// \copydoc BaseClass::dist(const WorldPoint&,WorldPoint&)
        virtual CoordType dist(const WorldPoint& point, WorldPoint& grad_dist) const override
        {
          return (this->cast()).compute_dist(point, grad_dist);
        }

        /// \copydoc BaseClass::signed_dist()
        virtual CoordType signed_dist(const WorldPoint& point) const override
        {
          return (this->cast()).compute_signed_dist(point);
        }

        /// \copydoc BaseClass::signed_dist(const WorldPoint&,WorldPoint&)
        virtual CoordType signed_dist(const WorldPoint& point, WorldPoint& grad_dist) const override
        {
          return (this->cast()).compute_signed_dist(point, grad_dist);
        }

      }; // class ChartCRTP<...>

    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAT
#endif // KERNEL_GEOMETRY_ATLAS_CHART_HPP
