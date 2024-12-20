// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#include "kernel/util/assertion.hpp"
#include <kernel/geometry/atlas/chart.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/string.hpp>
#include <kernel/geometry/cgal.hpp>
#include <kernel/runtime.hpp>

#include <memory>

namespace FEAT
{
  namespace Geometry
  {
    namespace Atlas
    {
      /// CGAL chart traits
      struct CGALSurfaceMeshTraits
      {
        /// No explicit map is available in general
        static constexpr bool is_explicit = false;
        /// We support implicit projection
        static constexpr bool is_implicit = true;
        /// This is a world_dim dimensional object
        static constexpr int world_dim = 3;
        /// If there was a parametrization, it would be the object's shape dim
        static constexpr int param_dim = 2;
      };

      /**
       * \brief Boundary description by a cgal based surface in 3d
       *
       * \tparam Mesh_
       * Type for the mesh this boundary description refers to
       *
       * \author Maximilian Esser
       *         Based on SurfaceMesh, so credit to Jordi Paul and Peter Zajac.
       *
       */
      template<typename Mesh_>
      class CGALSurfaceMesh :
#ifdef FEAT_HAVE_CGAL
        public ChartCRTP<CGALSurfaceMesh<Mesh_>, Mesh_, CGALSurfaceMeshTraits>
      {
      public:
        /// Dimension of the CGAL mesh
        static constexpr int shape_dim = 2;
        /// The CRTP base class
        typedef ChartCRTP<CGALSurfaceMesh<Mesh_>, Mesh_, CGALSurfaceMeshTraits> BaseClass;
        /// Floating point type for coordinates
        typedef typename BaseClass::CoordType CoordType;
        /// Vector type for world points
        typedef typename BaseClass::WorldPoint WorldPoint;


      private:
        ///Our cgal handler
        Geometry::CGALWrapper<CoordType> cgal;

      public:
        /// delete default constructor explicitly
        CGALSurfaceMesh() = delete;
        /// For now, no move constructor
        CGALSurfaceMesh(CGALSurfaceMesh&&) = delete;
        /// And no move assign
        CGALSurfaceMesh& operator=(CGALSurfaceMesh&&) = delete;
        ///Also, no copy operator
        CGALSurfaceMesh(const CGALSurfaceMesh&) = delete;
        ///Rule of 5
        CGALSurfaceMesh& operator=(const CGALSurfaceMesh&) = delete;


        CGALSurfaceMesh(const String& filename, Geometry::CGALFileMode file_mode) :
        cgal{filename, file_mode}
        {}

        CGALSurfaceMesh(std::istream & file, CGALFileMode file_mode) :
        cgal{file, file_mode}
        {}

        static std::unique_ptr<CGALSurfaceMesh<Mesh_>> create_cgal_surface_mesh(std::istream& file, Geometry::CGALFileMode file_mode)
        {
          return std::make_unique<CGALSurfaceMesh<Mesh_>>(file, file_mode);
        }


        ~CGALSurfaceMesh() = default;

        /// \copydoc ChartBase::bytes
        std::size_t bytes() const override
        {
          return cgal.bytes();
        }

        /** \copydoc ChartBase::get_type() */
        String get_type() const override
        {
          return "CGALSurfaceMesh";
        }

        /// \copydoc ChartBase:transform()
        void transform(const WorldPoint& origin, const WorldPoint& angles, const WorldPoint& offset) override
        {
          // create rotation matrix
          Tiny::Matrix<CoordType, 3, 3> rot;
          rot.set_rotation_3d(angles[0], angles[1], angles[2]);
          Tiny::Vector<CoordType, 3> trans = (-1) * rot * origin + offset;
          cgal.transform(rot, trans);
        }

        /** \copydoc ChartBase::write */
        void write(std::ostream& /*os*/, const String& /*sindent*/) const override
        {
          XABORTM("Not implemented yet");
        }

        /**
         * \brief Projects a single point to the surface given by the surface mesh
         *
         * \param[in,out] point
         * The point that gets projected.
         */
        void project_point(WorldPoint& point) const
        {
          point = cgal.closest_point(point);
        }

        /**
         * \brief Projects a single point to the surface given by the surface mesh
         *
         * \param[in,out] point
         * The point that gets projected.
         *
         * \param[in,out] surface_normal
         * The normal on the connected surface element.
         *
         */
        void project_point(WorldPoint& point, WorldPoint& normal) const
        {
          WorldPoint tmp_norm;
          point = cgal.closest_point(point, tmp_norm);
          normal = tmp_norm;
        }


        /**
         * \brief Orthogonally projects all vertices at facets of a MeshPart
         *
         * \param[in,out] mesh
         * The mesh the Meshpart refers to and whose vertices are to be projected
         *
         * \param[in] meshpart
         * The meshpart identifying the boundary of the mesh that is to be projected
         *
         * The MeshPart has to contain facets and a topology for them. To keep the search as local as possible,
         * the MeshPart is traversed facet-wise and after finishing the current facet, its neighbors are added to
         * the search stack.
         *
         */
        void project_meshpart(Mesh_& mesh, const MeshPart<Mesh_>& meshpart) const
        {
          Index num_facets(meshpart.get_num_entities(shape_dim));

          // There is nothing to do if the meshpart does not have any cells (which refer to facets of the original
          // mesh)
          if(num_facets == Index(0))
            return;

          // The number of vertices that need to be projected
          const Index num_verts(meshpart.get_num_entities(0));
          // Mapping of vertices from the meshpart to the real mesh
          const auto& ts_verts(meshpart.template get_target_set<0>());

          //get vertex set of mesh
          auto& vtx(mesh.get_vertex_set());

          //we could simply go through the vertices of the meshpart, project and hope for the best...
          // but we should and can at least check, if the new facet is at about the same distance away
          for(Index i = 0; i < num_verts; ++i)
          {
            project_point(vtx[ts_verts[i]]);
          }
          //already done
        }


        /// \copydoc ChartBase::dist()
        CoordType compute_dist(const WorldPoint& point) const
        {
          return Math::sqrt(cgal.squared_distance(point[0], point[1], point[2]));
        }

        /// \copydoc ChartBase::dist()
        /// This makes
        CoordType compute_dist(const WorldPoint& point, WorldPoint& grad_distance) const
        {
          // WorldPoint normal;
          WorldPoint projected_point(point);
          project_point(projected_point);
          // project_point(projected_point, normal);

          grad_distance = (projected_point - point);
          CoordType dist = grad_distance.norm_euclid();

          // If the distance is too small, we set the gradient vector to zero
          if(grad_distance.norm_euclid() < Math::eps<CoordType>())
          {
            grad_distance.format(CoordType(0));
          }
          else
          {
            grad_distance.normalize();
          }

          return dist;
        }

        /// \copydoc ChartBase::signed_dist()
        CoordType compute_signed_dist(const WorldPoint& point) const
        {
          return cgal.point_inside(point[0], point[1], point[2]) ? -compute_dist(point) : compute_dist(point);
        }

        /// \copydoc ChartBase::signed_dist()
        CoordType compute_signed_dist(const WorldPoint& point, WorldPoint& grad_distance) const
        {
          return cgal.point_inside(point[0], point[1], point[2]) ? -compute_dist(point, grad_distance) : compute_dist(point, grad_distance);
        }
#else
      public ChartBase<Mesh_>  //dummy implementation... you should not create a CGALSurfaceMesh chart, if you do not have cgal loaded...
      {
      public:
        /// Dimension of the CGAL mesh
        static constexpr int shape_dim = 2;
        /// The CRTP base class
        typedef ChartBase<Mesh_> BaseClass;
        /// Floating point type for coordinates
        typedef typename BaseClass::CoordType CoordType;
        /// Vector type for world points
        typedef typename BaseClass::WorldPoint WorldPoint;
        /// Meshtype of Baseclass
        typedef typename BaseClass::MeshType MeshType;
        /// Parttype of Baseclass
        typedef typename BaseClass::PartType PartType;

        CGALSurfaceMesh()
        {
          XABORTM("ERROR: Can not create a CGALSurfaceMesh without CGAL enabled");
        }

        bool can_explicit() const override
        {
          return false;
        }

        bool can_implicit() const override
        {
          return false;
        }

        void adapt(MeshType&, const PartType&) const override
        {
          XABORTM("Thou shalt not arrive here!");
        }

        void adapt(PartType&, const PartType&) const override
        {
          XABORTM("Thou shalt not arrive here!");
        }

        void transform(const WorldPoint&, const WorldPoint&, const WorldPoint&) override
        {
          XABORTM("Thou shalt not arrive here!");
        }

        WorldPoint map(const WorldPoint&) const override
        {
          XABORTM("Thou shalt not arrive here!");
          return WorldPoint{};
        }

        WorldPoint project(const WorldPoint&) const override
        {
          XABORTM("Thou shalt not arrive here!");
          return WorldPoint{};
        }

        CoordType dist(const WorldPoint&) const override
        {
          XABORTM("Thou shalt not arrive here!");
          return CoordType{};
        }

        CoordType dist(const WorldPoint&, WorldPoint&) const override
        {
          XABORTM("Thou shalt not arrive here!");
          return CoordType{};
        }

        CoordType signed_dist(const WorldPoint&) const override
        {
          XABORTM("Thou shalt not arrive here!");
          return CoordType{};
        }

        CoordType signed_dist(const WorldPoint&, WorldPoint&) const override
        {
          XABORTM("Thou shalt not arrive here!");
          return CoordType{};
        }

        String get_type() const override
        {
          return "CGALSurfaceMeshDummy";
        }

        void write(std::ostream&, const String&) const override
        {
          return;
        }

#endif
      }; //class CGALSurfaceMesh


/*  For now, CGALSurfaceMeshParser should not be used... in any case, lets keep the implementation if it is required someday...
      template<typename Mesh_, typename ChartReturn_ = ChartBase<Mesh_>, bool enable_ = (Mesh_::shape_dim > 2)>
      class CGALSurfaceMeshChartParser :
        public Xml::DummyParser
      {
      public:
        explicit CGALSurfaceMeshChartParser(std::unique_ptr<ChartReturn_>&)
        {
          XABORTM("Thou shall not arrive here");
        }
      };

      template<typename Mesh_, typename ChartReturn_>
      class CGALSurfaceMeshChartParser<Mesh_, ChartReturn_, true> :
        public Xml::MarkupParser
      {
      private:
        typedef CGALSurfaceMesh<Mesh_> ChartType;
        typedef typename ChartType::CoordType CoordType;
        std::unique_ptr<ChartBase<Mesh_>>& _chart;

      public:
        explicit CGALSurfaceMeshChartParser(std::unique_ptr<ChartReturn_>& chart) :
          _chart(chart)
        {
        }

        bool attribs(std::map<String,bool>& attrs) const override
        {
          attrs.emplace("filename", true);
          return true;
        }

        void create(
          int iline,
          const String& sline,
          const String&,
          const std::map<String, String>& attrs,
          bool) override
        {
#ifdef FEAT_HAVE_CGAL
          String filename;

          // try to parse the filename
          if(!attrs.find("filename")->second.parse(filename))
            throw Xml::GrammarError(iline, sline, "Failed to parse filename");

          // get filemode
          CGALFileMode mode;
          auto split = filename.split_by_charset(".");
          if(split.back().compare_no_case("off") == 0)
          {
            mode = CGALFileMode::fm_off;
          }
          else if (split.back().compare_no_case("obj"))
          {
            mode = CGALFileMode::fm_obj;
          }
          else
          {
            XABORTM("ERROR: File extension ." + split.back() + " is not valid/implemented yet.");
          }

          ///TODO: Add options for affine transformation

          // everything seems fine, let's create the chart then
          _chart.reset(new ChartType(filename, mode));
#else
          std::cout << "ERROR: Trying to parse a CGALMeshpart chart without CGAL enabled\n";
          FEAT::Runtime::abort();
          //no unused warnings by casting to void:
          (void)iline;
          (void)sline;
          (void)attrs;
#endif
        }

        void close(int, const String&) override
        {
        }

        bool content(int, const String&) override
        {
          return false;
        }

        std::shared_ptr<Xml::MarkupParser> markup(int, const String&, const String&) override
        {
          return nullptr;
        }
      }; // class CGALSurfaceMeshChartParser<...>
*/
    } //namespace Atlas
  } //namespace Geometry
} //namespace FEAT
