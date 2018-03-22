#pragma once
#ifndef KERNEL_GEOMETRY_EXPORT_EPS_HPP
#define KERNEL_GEOMETRY_EXPORT_EPS_HPP 1

// includes, FEAT
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/math.hpp>

// includes, system
#include <iostream>
#include <fstream>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief EPS exporter class
     *
     * This class implements an exporter for the Encapsulated Post-Script (EPS)
     * file format, which can be used to generate mesh plots for LaTeX documents.
     *
     * \note
     * This class can only export 2D meshes.
     *
     * \author Peter Zajac
     */
    class ExportEPS
    {
    public:
      /**
       * \brief Exports a mesh to an EPS file.
       *
       * \param[in,out] os
       * The output stream to which to write to.
       *
       * \param[in] mesh
       * The 2D conformal mesh that is to be exported.
       *
       * \param[in] width, height
       * The maximum dimensions of the bounding box in millimeters.
       *
       * \param[in] stroke
       * The stroke width in millimeters.
       *
       * \param[in] extra
       * The extra bounding box offset in millimeters.
       */
      template<typename Shape_, typename DataType_>
      static void write(
        std::ostream& os,
        const ConformalMesh<Shape_, 2, DataType_>& mesh,
        const double width = 100.0,
        const double height = 100.0,
        const double stroke = 0.1,
        const double extra = 0.0)
      {
        XASSERTM(width  >= 1.0, "invalid bounding box width");
        XASSERTM(height >= 1.0, "invalid bounding box height");
        XASSERTM(stroke > 0.01, "invalid stroke width");
        XASSERTM(extra >= 0.0, "invalid bounding box extra");

        // get vertices-at-edge
        const auto& vtx = mesh.get_vertex_set();
        const auto& idx = mesh.template get_index_set<1,0>();

        // compute minimal and maximal x/y coords
        DataType_ x_min, x_max, y_min, y_max;
        x_min = x_max = vtx[0][0];
        y_min = y_max = vtx[0][1];
        for(Index i(1); i < vtx.get_num_vertices(); ++i)
        {
          const auto& v = vtx[i];
          Math::minimax(v[0], x_min, x_max);
          Math::minimax(v[1], y_min, y_max);
        }

        // compute x/y scale:
        // Note: PostScript units are "points"
        // 1 inch = 25.4 millimeters
        // 1 inch = 72 points
        // ==> 1 millimeter = 72 / 25.4 points
        const double x_sc = (72.0 / 25.4) * width  / double(x_max - x_min);
        const double y_sc = (72.0 / 25.4) * height / double(y_max - y_min);

        // compute scaling factor
        const double scale = Math::min(x_sc, y_sc);

        // compute bounding box (rounding up)
        const int box_x = int(scale * double(x_max - x_min) + 2.0*extra) + 1;
        const int box_y = int(scale * double(y_max - y_min) + 2.0*extra) + 1;

        // compute offset from extra
        const double box_o = 72.0 * extra / 25.4;

        // write EPS header
        os << "%!!PS-Adobe-3.0 EPSF-3.0" << std::endl;

        // write bounding box
        os << "%%BoundingBox: 0 0 " << box_x << " " << box_y << std::endl;

        // set stroke width
        os << (72.0 * stroke / 25.4) << " setlinewidth" << std::endl;

        // begin a new path
        os << "newpath" << std::endl;

        // write all edges
        for(Index i(0); i < idx.get_num_entities(); ++i)
        {
          // get the vertices
          const auto& v0 = vtx[idx[i][0]];
          const auto& v1 = vtx[idx[i][1]];

          // transform vertices
          double v0x = box_o + (double(v0[0]) - x_min) * scale;
          double v0y = box_o + (double(v0[1]) - y_min) * scale;
          double v1x = box_o + (double(v1[0]) - x_min) * scale;
          double v1y = box_o + (double(v1[1]) - y_min) * scale;

          // write vertices
          os << v0x << " " << v0y << " moveto" << std::endl;
          os << v1x << " " << v1y << " lineto" << std::endl;
        }

        // draw path and show page
        os << "stroke" << std::endl;
        os << "showpage" << std::endl;
      }

      /**
       * \brief Exports a mesh to an EPS file.
       *
       * \param[in] filename
       * The name of the output file.
       *
       * \param[in] mesh
       * The 2D conformal mesh that is to be exported.
       *
       * \param[in] width, height
       * The maximum dimensions of the bounding box in millimeters.
       *
       * \param[in] stroke
       * The stroke width in millimeters.
       *
       * \param[in] extra
       * The extra bounding box offset in millimeters.
       */
      template<typename Shape_, typename DataType_>
      static void write(
        const String& filename,
        const ConformalMesh<Shape_, 2, DataType_>& mesh,
        const double width = 100.0,
        const double height = 100.0,
        const double stroke = 0.1,
        const double extra = 0.0)
      {
        std::ofstream ofs(filename);
        if(!ofs.is_open() || !ofs.good())
          throw FileError(String("Failed to open output file: ") + filename);
        write(ofs, mesh, width, height, stroke, extra);
        ofs.close();
      }
    };
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_EXPORT_EPS_HPP
