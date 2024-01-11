// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#include "feat_config.hpp"
#include "kernel/util/assertion.hpp"
#ifndef KERNEL_GEOMETRY_ATLAS_CGAL_CHART_HPP
#define KERNEL_GEOMETRY_ATLAS_CGAL_CHART_HPP 1

#include <kernel/geometry/atlas/chart.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/string.hpp>
#include <kernel/geometry/cgal.hpp>


namespace FEAT{
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
        Geometry::CGALWrapper cgal;

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


        virtual ~CGALSurfaceMesh() = default;

        /// \copydoc ChartBase::bytes
        virtual std::size_t bytes() const override
        {
          return cgal.bytes();
        }

        /** \copydoc ChartBase::get_type() */
        virtual String get_type() const override
        {
          return "CGALSurfaceMesh";
        }

        /// \copydoc ChartBase:transform()
        virtual void transform(const WorldPoint& origin, const WorldPoint& angles, const WorldPoint& offset) override
        {
          // create rotation matrix as double
          Tiny::Matrix<double, 3, 3> rot;
          rot.set_rotation_3d(double(angles[0]), double(angles[1]), double(angles[2]));
          Tiny::Vector<double, 3> trans = (-1) * rot * origin + offset;
          cgal.transform(rot, trans);
        }

        /** \copydoc ChartBase::write */
        virtual void write(std::ostream& /*os*/, const String& /*sindent*/) const override
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
          point = WorldPoint(cgal.closest_point(CGALWrapper::PointType(point)));
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
          CGALWrapper::PointType tmp_norm;
          point = WorldPoint(cgal.closest_point(CGALWrapper::PointType(point), tmp_norm));
          normal = WorldPoint(tmp_norm);
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
          return CoordType(Math::sqrt(cgal.squared_distance(double(point[0]), double(point[1]), double(point[2]))));
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
          return cgal.point_inside(double(point[0]), double(point[1]), double(point[2])) ? -compute_dist(point) : compute_dist(point);
        }

        /// \copydoc ChartBase::signed_dist()
        CoordType compute_signed_dist(const WorldPoint& point, WorldPoint& grad_distance) const
        {
          return cgal.point_inside(double(point[0]), double(point[1]), double(point[2])) ? -compute_dist(point, grad_distance) : compute_dist(point, grad_distance);
        }


      }; //class CGALSurfaceMesh

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

        virtual bool attribs(std::map<String,bool>& attrs) const override
        {
          attrs.emplace("filename", true);
          return true;
        }

        virtual void create(
          int iline,
          const String& sline,
          const String&,
          const std::map<String, String>& attrs,
          bool) override
        {
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

          //check if filename begins with env variable accesser
          if(filename.starts_with("FEAT_SOURCE_DIR"))
            filename.replace_all("FEAT_SOURCE_DIR", FEAT_SOURCE_DIR);

          if(filename.starts_with("FEAT_BINARY_DIR"))
            filename.replace_all("FEAT_BINARY_DIR", FEAT_BINARY_DIR);

          ///TODO: Add options for affine transformation

          // everything seems fine, let's create the chart then
          _chart.reset(new ChartType(filename, mode));
        }

        virtual void close(int, const String&) override
        {
        }

        virtual bool content(int, const String&) override
        {
          return false;
        }

        virtual std::shared_ptr<Xml::MarkupParser> markup(int, const String&, const String&) override
        {
          return nullptr;
        }
      }; // class CGALSurfaceMeshChartParser<...>
    } //namespace Atlas
  } //namespace Geometry
} //namespace FEAT
#endif
