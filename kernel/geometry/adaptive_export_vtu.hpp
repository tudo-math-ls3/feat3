// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GEOMETRY_ADAPTIVE_EXPORT_VTU_HPP
#define KERNEL_GEOMETRY_ADAPTIVE_EXPORT_VTU_HPP 1

#include <kernel/geometry/adaptive_mesh.hpp>
#include <kernel/geometry/adaptive_mesh_layer.hpp>
#include <kernel/geometry/export_vtk.hpp>

namespace FEAT
{
  namespace Geometry
  {
    template<typename MeshType_>
    class AdaptiveExportVTU;

    template<typename TemplateSet_, typename Shape_>
    class AdaptiveExportVTU<Geometry::AdaptiveMesh<TemplateSet_, Shape_>>
    {
    public:
      typedef Shape_ ShapeType;
      typedef Geometry::AdaptiveMesh<TemplateSet_, Shape_> MeshType;

      /// our VTK shape type
      typedef Intern::VTKShape<ShapeType> VTKShapeType;

      static constexpr int shape_dim = ShapeType::dimension;

    protected:
      /// the underlying regular mesh
      const typename MeshType::FoundationMeshType& _mesh_reg;
      /// the chosen adaptive mesh layer
      const Geometry::AdaptiveMeshLayer<MeshType>& _mesh_ada;

      // number of vertices in regular/adaptive mesh
      Index _num_verts_reg, _num_verts_ada, _num_verts_total;
      // number of elements in regular/adaptive mesh
      Index _num_elems_reg, _num_elems_ada, _num_elems_total;

      // regular element mask
      std::vector<char> _reg_elem_mask;

      std::deque<std::pair<String, std::vector<double>>> _vertex_scalars;
      std::deque<std::pair<String, std::vector<double>>> _cell_scalars;

    public:
      explicit AdaptiveExportVTU(Geometry::AdaptiveMeshLayer<MeshType>& mesh) :
        _mesh_reg(mesh.foundation_mesh()),
        _mesh_ada(mesh),
        _num_verts_reg(_mesh_reg.get_num_entities(0)),
        _num_verts_ada(mesh.get_num_entities(0)),
        _num_verts_total(_num_verts_reg + _num_verts_ada),
        _num_elems_reg(_mesh_reg.get_num_entities(shape_dim)),
        _num_elems_ada(mesh.get_num_entities(shape_dim)),
        _num_elems_total(
          _num_elems_reg + _num_elems_ada - mesh.adaptive_mesh().get_num_entities(Geometry::Layer{0}, shape_dim)),
        _reg_elem_mask(_num_elems_reg, 1)
      {
        // figure out which regular elements are refined adaptively and mask them out
        Index count = _num_elems_ada;
        for(Index i(0); i < _num_elems_reg; ++i)
          count += Index(_reg_elem_mask[i] = (mesh.adaptive_mesh().get_overlap_cell(i).has_value() ? 0 : 1));
        XASSERTM(_num_elems_total == count, "invalid element count");
      }

      virtual ~AdaptiveExportVTU()
      {
      }

      // no copy, no problems
      AdaptiveExportVTU(const AdaptiveExportVTU&) = delete;
      AdaptiveExportVTU& operator=(const AdaptiveExportVTU&) = delete;

      template<typename BridgeVector_>
      void add_vertex_scalar(const String& name, const BridgeVector_& bridge_vector)
      {
        XASSERTM(bridge_vector.regular().size() == _num_verts_reg, "invalid number of regular values in bridge vector");
        XASSERTM(
          bridge_vector.adaptive().size() == _num_verts_ada,
          "invalid number of adaptive values in bridge vector");

        const auto* val_reg = bridge_vector.regular().elements();
        const auto* val_ada = bridge_vector.adaptive().elements();

        std::vector<double> values(_num_verts_total);

        for(Index i(0); i < _num_verts_reg; ++i)
          values[i] = double(val_reg[i]);
        for(Index i(0); i < _num_verts_ada; ++i)
          values[_num_verts_reg + i] = double(val_ada[i]);

        _vertex_scalars.push_back(std::make_pair(name, std::move(values)));
      }

      template<typename BridgeVector_>
      void add_cell_scalar(const String& name, const BridgeVector_& bridge_vector)
      {
        XASSERTM(bridge_vector.regular().size() == _num_elems_reg, "invalid number of regular values in bridge vector");
        XASSERTM(
          bridge_vector.adaptive().size() == _num_elems_ada,
          "invalid number of adaptive values in bridge vector");

        const auto* val_reg = bridge_vector.regular().elements();
        const auto* val_ada = bridge_vector.adaptive().elements();

        std::vector<double> values(_num_elems_total);

        Index j(0);
        for(Index i(0); i < _num_elems_reg; ++i)
        {
          if(_reg_elem_mask[i] > 0)
          {
            values[j] = double(val_reg[i]);
            ++j;
          }
        }
        for(Index i(0); i < _num_elems_ada; ++i, ++j)
          values[j] = double(val_ada[i]);
        XASSERT(j == _num_elems_total);

        _cell_scalars.push_back(std::make_pair(name, std::move(values)));
      }

      void write(const String& filename) const
      {
        String vtu_name(filename + ".vtu");
        std::ofstream ofs(vtu_name.c_str());
        if(!(ofs.is_open() && ofs.good()))
          throw FileError("Failed to create '" + vtu_name + "'");

        // write
        write_vtu(ofs);

        // and close
        ofs.close();
      }

      void write_vtu(std::ostream& os) const
      {
        // fetch basic information
        const int verts_per_cell = Shape::FaceTraits<ShapeType, 0>::count;

        // write VTK header
        os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">\n";

        // write mesh header
        os << "<UnstructuredGrid>\n";
        os << "<Piece NumberOfPoints=\"" << _num_verts_total << "\" NumberOfCells=\"" << _num_elems_total << "\">\n";

        // write point data
        if((!_vertex_scalars.empty()) /*|| (!_vertex_vectors.empty())*/)
        {
          os << "<PointData>\n";

          // write vertex variables
          for(Index i(0); i < Index(_vertex_scalars.size()); ++i)
          {
            const auto& var(_vertex_scalars[i]);
            os << "<DataArray type=\"Float64\" Name=\"" << var.first << "\" Format=\"ascii\">\n";
            for(Index j(0); j < _num_verts_total; ++j)
            {
              os << stringify_fp_sci(var.second[j]) << "\n";
            }
            os << "</DataArray>\n";
          }
          // write vertex fields
          /*for(Index i(0); i < Index(_vertex_vectors.size()); ++i)
          {
            const auto& var(_vertex_vectors[i]);
            os << "<DataArray type=\"Float64\" Name=\"" << var.first;
            os <<"\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
            for(Index j(0); j < _num_verts_total; ++j)
            {
              os << stringify_fp_sci(var.second[3*j+0], _var_prec) << " ";
              os << stringify_fp_sci(var.second[3*j+1], _var_prec) << " ";
              os << stringify_fp_sci(var.second[3*j+2], _var_prec) << "\n";
            }
            os << "</DataArray>\n";
          }*/

          os << "</PointData>\n";
        }

        // write cell data
        if(!_cell_scalars.empty() /*|| !_cell_vectors.empty()*/)
        {
          os << "<CellData>\n";
          if(!_cell_scalars.empty())
          {
            for(Index i(0); i < Index(_cell_scalars.size()); ++i)
            {
              const auto& var(_cell_scalars[i]);
              os << "<DataArray type=\"Float64\" Name=\"" << var.first << "\" Format=\"ascii\">\n";
              for(Index j(0); j < _num_elems_total; ++j)
              {
                os << stringify_fp_sci(var.second[j]) << "\n";
              }
              os << "</DataArray>\n";
            }
          }

          /*if(!_cell_vectors.empty())
          {
            // write cell fields
            for(Index i(0); i < Index(_cell_vectors.size()); ++i)
            {
              const auto& var(_cell_vectors[i]);
              os << "<DataArray type=\"Float64\" Name=\"" << var.first;
              os <<"\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
              for(Index j(0); j < _num_elems_total; ++j)
              {
                os << stringify_fp_sci(var.second[3*j+0], _var_prec) << " ";
                os << stringify_fp_sci(var.second[3*j+1], _var_prec) << " ";
                os << stringify_fp_sci(var.second[3*j+2], _var_prec) << "\n";
              }
              os << "</DataArray>\n";
            }
          }*/
          os << "</CellData>\n";
        }

        // write vertices
        const auto& vtx_reg = _mesh_reg.get_vertex_set();
        const auto& vtx_ada = _mesh_ada.get_vertex_set();

        os << "<Points>\n";
        os << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">\n";
        // write regular vertices
        for(Index i(0); i < _num_verts_reg; ++i)
        {
          os << vtx_reg[i][0];
          for(int j(1); j < shape_dim; ++j)
          {
            os << " " << vtx_reg[i][j];
          }
          for(int j(shape_dim); j < 3; ++j)
          {
            os << " 0";
          }
          os << "\n";
        }
        // write adaptive vertices
        for(Index i(0); i < _num_verts_ada; ++i)
        {
          os << vtx_ada[i][0];
          for(int j(1); j < shape_dim; ++j)
          {
            os << " " << vtx_ada[i][j];
          }
          for(int j(shape_dim); j < 3; ++j)
          {
            os << " 0";
          }
          os << "\n";
        }
        os << "</DataArray>\n";
        os << "</Points>\n";

        // write cells
        const auto& idx_reg = _mesh_reg.template get_index_set<shape_dim, 0>();
        const auto& idx_ada = _mesh_ada.template get_index_set<shape_dim, 0>();
        os << "<Cells>\n";
        os << "<DataArray type=\"UInt32\" Name=\"connectivity\">\n";
        // write regular vertices
        for(Index i(0); i < _num_elems_reg; ++i)
        {
          // is this element masked?
          if(_reg_elem_mask[i] > 0)
          {
            os << idx_reg(i, VTKShapeType::map(0));
            for(int j(1); j < verts_per_cell; ++j)
            {
              os << " " << idx_reg(i, VTKShapeType::map(j));
            }
            os << "\n";
          }
        }
        // write adaptive vertices
        for(Index i(0); i < _num_elems_ada; ++i)
        {
          os << (_num_verts_reg + idx_ada(i, VTKShapeType::map(0)));
          for(int j(1); j < verts_per_cell; ++j)
          {
            os << " " << (_num_verts_reg + idx_ada(i, VTKShapeType::map(j)));
          }
          os << "\n";
        }
        os << "</DataArray>\n";
        os << "<DataArray type=\"UInt32\" Name=\"offsets\">\n";
        for(Index i(0); i < _num_elems_total; ++i)
        {
          os << ((i + 1) * verts_per_cell) << "\n";
        }
        os << "</DataArray>\n";
        os << "<DataArray type=\"UInt32\" Name=\"types\">\n";
        for(Index i(0); i < _num_elems_total; ++i)
        {
          os << VTKShapeType::type << "\n";
        }
        os << "</DataArray>\n";
        os << "</Cells>\n";

        // finish
        os << "</Piece>\n";
        os << "</UnstructuredGrid>\n";
        os << "</VTKFile>\n";
      }
    }; // class AdaptiveExportVTU<Geometry::AdaptiveMesh<...>>
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_ADAPTIVE_EXPORT_VTU_HPP
