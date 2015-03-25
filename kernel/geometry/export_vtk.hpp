#pragma once
#ifndef KERNEL_GEOMETRY_EXPORT_VTK_HPP
#define KERNEL_GEOMETRY_EXPORT_VTK_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/structured_mesh.hpp>

// includes, STL
#include <fstream>
#include <vector>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      struct VTKHelper;

      template<>
      struct VTKHelper< Shape::Simplex<1> >
      {
        static constexpr int type = 3; // VTK_LINE
        static inline int map(int i)
        {
          return i;
        }
      };

      template<>
      struct VTKHelper< Shape::Simplex<2> >
      {
        static constexpr int type = 5; // VTK_TRIANGLE
        static inline int map(int i)
        {
          return i;
        }
      };

      template<>
      struct VTKHelper< Shape::Simplex<3> >
      {
        static constexpr int type = 10; // VTK_TETRA
        static inline int map(int i)
        {
          return i;
        }
      };

      template<>
      struct VTKHelper< Shape::Hypercube<1> >
      {
        static constexpr int type = 3; // VTK_LINE
        static inline int map(int i)
        {
          return i;
        }
      };

      template<>
      struct VTKHelper< Shape::Hypercube<2> >
      {
        static constexpr int type = 9; // VTK_QUAD
        static inline int map(int i)
        {
          static int v[] = {0, 1, 3, 2};
          return v[i];
        }
      };

      template<>
      struct VTKHelper< Shape::Hypercube<3> >
      {
        static constexpr int type = 12; // VTK_HEXAHEDRON
        static inline int map(int i)
        {
          static int v[] = {0, 1, 3, 2, 4, 5, 7, 6};
          return v[i];
        }
      };

      template<typename Mesh_>
      struct VTKHeader;

      template<
        typename Shape_,
        int num_coords_,
        int stride_,
        typename Coord_>
      struct VTKHeader< ConformalMesh<Shape_, num_coords_, stride_, Coord_> >
      {
        typedef ConformalMesh<Shape_, num_coords_, stride_, Coord_> MeshType;
        static void write(std::ostream& os, const MeshType& mesh)
        {
          // write mesh type
          os << "DATASET UNSTRUCTURED_GRID" << std::endl;

          // write vertex coordinates
          const typename MeshType::VertexSetType& vtx = mesh.get_vertex_set();
          Index num_verts = vtx.get_num_vertices();
          int num_coords = vtx.get_num_coords();
          ASSERT_((num_coords >= 1) && (num_coords <= 3));
          os << "POINTS " << num_verts << " double" << std::endl;
          for(Index i(0); i < num_verts; ++i)
          {
            os << vtx[i][0];
            for(int j(1); j < num_coords; ++j)
            {
              os << " " << vtx[i][j];
            }
            for(int j(num_coords); j < 3; ++j)
            {
              os << " 0.0";
            }
            os << std::endl;
          }

          typedef VTKHelper<typename MeshType::ShapeType> VTKHelperType;

          // fetch index set
          const typename MeshType::template IndexSet<MeshType::shape_dim,0>::Type& idx =
            mesh.template get_index_set<MeshType::shape_dim, 0>();
          Index num_cells = mesh.get_num_entities(MeshType::shape_dim);
          int num_idx = idx.get_num_indices();

          // write cells
          os << "CELLS " << num_cells << " " << (Index(num_idx+1)*num_cells) << std::endl;
          for(Index i(0); i < num_cells; ++i)
          {
            os << num_idx;
            for(int j(0); j < num_idx; ++j)
            {
              os << " " << idx[i][VTKHelperType::map(j)];
            }
            os << std::endl;
          }

          // write cell types
          os << "CELL_TYPES " << num_cells << std::endl;
          for(Index i(0); i < num_cells; ++i)
          {
            os << VTKHelperType::type << std::endl;
          }
        }
      };

      template<
        int shape_dim_,
        int num_coords_,
        int stride_,
        typename Coord_>
      struct VTKHeader< StructuredMesh<shape_dim_, num_coords_, stride_, Coord_> >
      {
        typedef StructuredMesh<shape_dim_, num_coords_, stride_, Coord_> MeshType;
        static void write(std::ostream& os, const MeshType& mesh)
        {
          // write mesh type
          os << "DATASET STRUCTURED_GRID" << std::endl;

          // write dimensions
          os << "DIMENSIONS";
          for(int i(0); i < shape_dim_; ++i)
          {
            os << " " << (mesh.get_num_slices(i) + 1);
          }
          for(int i(shape_dim_); i < 3; ++i)
          {
            os << " 1";
          }
          os << std::endl;

          // write vertex coordinates
          const typename MeshType::VertexSetType& vtx = mesh.get_vertex_set();
          Index num_verts = vtx.get_num_vertices();
          int num_coords = vtx.get_num_coords();
          ASSERT_((num_coords >= 1) && (num_coords <= 3));
          os << "POINTS " << num_verts << " double" << std::endl;
          for(Index i(0); i < num_verts; ++i)
          {
            os << vtx[i][0];
            for(int j(1); j < num_coords; ++j)
            {
              os << " " << vtx[i][j];
            }
            for(int j(num_coords); j < 3; ++j)
            {
              os << " 0.0";
            }
            os << std::endl;
          }
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Provisional VTK exporter class template
     *
     * This class template is a provisional VTK exporter which will be replaced by a more mature one later.
     * \author Peter Zajac
     */
    template<typename Mesh_>
    class ExportVTK
    {
    public:
      /// mesh type
      typedef Mesh_ MeshType;
      /// vertex set type
      typedef typename MeshType::VertexSetType VertexSetType;

    protected:
      typedef std::pair<String, double*> VarPair;
      typedef std::vector<VarPair> VarVector;

      /// reference to mesh to be exported
      const MeshType& _mesh;
      /// number of vertices in mesh
      Index _num_verts;
      /// number of cells in mesh
      Index _num_cells;
      /// vertex variable list
      VarVector _vars_vertex;
      /// cell variable list
      VarVector _vars_cell;

    public:
      explicit ExportVTK(const MeshType& mesh) :
        _mesh(mesh),
        _num_verts(mesh.get_num_entities(0)),
        _num_cells(mesh.get_num_entities(MeshType::shape_dim))
      {
        CONTEXT("ExportVTK::ExportVTK()");
      }

      virtual ~ExportVTK()
      {
        CONTEXT("ExportVTK::~ExportVTK()");
        while(!_vars_cell.empty())
        {
          delete [] _vars_cell.back().second;
          _vars_cell.pop_back();
        }
        while(!_vars_vertex.empty())
        {
          delete [] _vars_vertex.back().second;
          _vars_vertex.pop_back();
        }
      }

      template<typename T_>
      void add_scalar_vertex(const String& name, const T_* data)
      {
        ASSERT_(data != nullptr);
        double* d = new double[_num_verts];
        for(Index i(0); i < _num_verts; ++i)
        {
          d[i] = double(data[i]);
        }
        _vars_vertex.push_back(VarPair(name, d));
      }

      template<typename T_>
      void add_scalar_cell(const String& name, const T_* data)
      {
        ASSERT_(data != nullptr);
        double* d = new double[_num_cells];
        for(Index i(0); i < _num_cells; ++i)
        {
          d[i] = double(data[i]);
        }
        _vars_cell.push_back(VarPair(name, d));
      }

      bool write(const String& filename) const
      {
        CONTEXT("ExportVTK::begin()");

        // try to open a file
        std::ofstream ofs(filename.c_str());
        if(!(ofs.is_open() && ofs.good()))
          return false;

        // write VTK header
        ofs << "# vtk DataFile Version 2.0" << std::endl;
        ofs << "Generated by FEAST v" << version_major << "." << version_minor << "." << version_patch << std::endl;
        ofs << "ASCII" << std::endl;

        // write mesh header
        Intern::VTKHeader<MeshType>::write(ofs, _mesh);

        // write vertex variables
        if(!_vars_vertex.empty())
        {
          ofs << "POINT_DATA " << _num_verts << std::endl;
          for(Index i(0); i < Index(_vars_vertex.size()); ++i)
          {
            const VarPair& var(_vars_vertex[i]);
            ofs << "SCALARS " << var.first << " double 1" << std::endl;
            ofs << "LOOKUP_TABLE default" << std::endl;
            for(Index j(0); j < _num_verts; ++j)
            {
              ofs << var.second[j] << std::endl;
            }
          }
        }

        // write cell variables
        if(!_vars_cell.empty())
        {
          ofs << "CELL_DATA " << _num_cells << std::endl;
          for(Index i(0); i < Index(_vars_cell.size()); ++i)
          {
            const VarPair& var(_vars_cell[i]);
            ofs << "SCALARS " << var.first << " double 1" << std::endl;
            ofs << "LOOKUP_TABLE default" << std::endl;
            for(Index j(0); j < _num_cells; ++j)
            {
              ofs << var.second[j] << std::endl;
            }
          }
        }

        // close output stream
        ofs.close();

        // okay
        return true;
      }
    }; // class ExportVTK
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_EXPORT_VTK_HPP
