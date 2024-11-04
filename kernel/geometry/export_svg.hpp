// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

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
     * \brief SVG exporter class
     *
     * This class implements an exporter for the Scalable Vector Graphics (SVG) file format,
     * which can be displayed directly in any modern internet browser or converted to another
     * vector graphics format via some external tool (e.g. Inkscape) to be used in LaTeX documents.
     *
     * \note
     * This class can only export 2D meshes.
     *
     * \author Pia Ritter
     */
    class ExportSVG
    {
    public:
      /**
       * \brief Exports a mesh to a SVG file.
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

        const double x_sc =  width  / double(x_max - x_min);
        const double y_sc = height / double(y_max - y_min);

        // compute scaling factor
        const double scale = Math::min(x_sc, y_sc);

        // write SVG header
        os << "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n";
        os << "<svg version=\"1.1\" id=\"svg2\" ";

        const double X = double(x_max - x_min);
        const double Y = double(y_max - y_min);

        const double adjusted_width = scale * X + 2.0 * extra;
        const double adjusted_height = scale * Y + 2.0 * extra;

        // write dimensions
        os << "width=\"" << adjusted_width << "mm\" height=\"" << adjusted_height << "mm\" ";

        // write viewbox
        os << "viewBox=\"" << (-extra/scale) << " " << (-extra/scale) << " ";
        os << (x_max - x_min + 2.0*extra/scale) << " " << (y_max - y_min + 2.0*extra/scale);
        os << "\" xmlns=\"http://www.w3.org/2000/svg\">\n";

        // start path and set stroke width
        os << "<path style=\"fill:none;stroke:#000000;";
        os << "stroke-width:" << stroke/scale << ";stroke-opacity:1\" id=\"mesh\" d=\" \n";

        // write all edges
        for(Index i(0); i < idx.get_num_entities(); ++i)
        {
          // get the vertices
          const auto& v0 = vtx[idx[i][0]];
          const auto& v1 = vtx[idx[i][1]];

          // transform vertices
          double v0x = (double(v0[0]) - x_min);
          double v0y = (y_max - double(v0[1]));
          double v1x = (double(v1[0]) - x_min);
          double v1y = (y_max - double(v1[1]));

          // write vertices
          os << "M " << v0x << " " << v0y;
          os << " L " << v1x << " " << v1y << "\n";
        }

        // finish
        os << "\" />\n";
        os << "</svg>";
      }

      /**
       * \brief Exports a mesh to a SVG file.
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
