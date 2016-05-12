#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_CHART_HPP
#define KERNEL_GEOMETRY_ATLAS_CHART_HPP 1

// includes, FEAST
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/xml_scanner.hpp>

namespace FEAST
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
         * \brief Adapts a mesh using this chart.
         *
         * \param[inout] mesh
         * The mesh that is to be adapted.
         *
         * \param[in] part
         * The mesh part that describes the part to adapt.
         */
        virtual void adapt(MeshType& mesh, const PartType& part) const = 0;

        /**
         * \brief Adapts a mesh part using this chart.
         *
         * \param[inout] mesh
         * The mesh part that is to be adapted.
         *
         * \param[in] part
         * The mesh part that describes the part to adapt.
         */
        virtual void adapt(PartType& mesh, const PartType& part) const = 0;

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
        };

        template<bool enable_>
        struct ExplicitChartHelper
        {
          template<typename CT_, typename MT_, typename PT_>
          static bool adapt(const CT_&, MT_&, const PT_&)
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
            typedef typename PT_::MeshAttributeType AttributeType;

            // First of all, check whether the chart can really perform
            // explicit adaption
            if(!chart.can_explicit())
              return false;

            // Try to fetch the parametrisation attribute.
            const AttributeType* attrib = part.find_attribute("param");
            if(attrib == nullptr)
              return false;

            // We have the attribute; check whether it matches our chart
            int attrib_dim = attrib->get_num_coords();
            if(attrib_dim != CT_::param_dim)
            {
              // invalid attribute dimension
              throw InternalError("Invalid chart attribute dimension");
            }

            // Get the vertex set of the mesh
            VertexSetType& vtx = mesh.get_vertex_set();

            // Get the vertex target set of the part
            const TargetSet& vidx = part.template get_target_set<0>();

            // loop over all vertices in the mesh part
            Index num_vtx = vidx.get_num_entities();
            for(Index i(0); i < num_vtx; ++i)
            {
              // apply the chart's map function
              chart.map(
                 reinterpret_cast<      WorldPoint&>(vtx[vidx[i]]),
                *reinterpret_cast<const ParamPoint*>((*attrib)[i])
                );
            }

            // okay
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
        typedef typename PartType::MeshAttributeType AttributeType;
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
        typedef Tiny::Vector<CoordType, world_dim> WorldPoint;
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
        /**
         * \brief Specifies whether the chart can perform explicit projection.
         *
         * This function returns #is_explicit by default, but it may be
         * overridden by the derived class in case that explicit projection
         * can be disabled at runtime (e.g. due to missing parameters).
         */
        bool can_explicit() const
        {
          return is_explicit;
        }

        /**
         * \brief Specifies whether the chart can perform implicit projection.
         *
         * This function returns #is_implicit by default, but it may be
         * overridden by the derived class in case that implicit projection
         * can be disabled at runtime.
         */
        bool can_implicit() const
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
          if(MeshType::world_dim != world_dim)
            throw InternalError("Mesh/Chart world dimension mismatch");

          // Try to adapt explicity
          if(Intern::ExplicitChartHelper<is_explicit>::adapt(cast(), mesh, part))
            return;

          // Try to adapt implicitly
          if(Intern::ImplicitChartHelper<is_implicit>::adapt(cast(), mesh, part))
            return;

          // If we come out here, we have no way of adaption...
          throw InternalError("No adaption possible");
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
          throw InternalError("Adaption of MeshPart not possible yet");
        }

      }; // class ChartCRTP<...>

    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_ATLAS_CHART_HPP
