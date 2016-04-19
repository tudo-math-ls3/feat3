#pragma once
#ifndef KERNEL_GEOMETRY_EXPORT_VTK_HPP
#define KERNEL_GEOMETRY_EXPORT_VTK_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/structured_mesh.hpp>
#include <kernel/util/file_error.hpp>
#include <kernel/util/function_scheduler.hpp>

// includes, STL
#include <fstream>
#include <vector>
#include <deque>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      struct VTKShape;

      template<>
      struct VTKShape< Shape::Simplex<1> >
      {
        static constexpr int type = 3; // VTK_LINE
        static inline int map(int i)
        {
          return i;
        }
      };

      template<>
      struct VTKShape< Shape::Simplex<2> >
      {
        static constexpr int type = 5; // VTK_TRIANGLE
        static inline int map(int i)
        {
          return i;
        }
      };

      template<>
      struct VTKShape< Shape::Simplex<3> >
      {
        static constexpr int type = 10; // VTK_TETRA
        static inline int map(int i)
        {
          return i;
        }
      };

      template<>
      struct VTKShape< Shape::Hypercube<1> >
      {
        static constexpr int type = 3; // VTK_LINE
        static inline int map(int i)
        {
          return i;
        }
      };

      template<>
      struct VTKShape< Shape::Hypercube<2> >
      {
        static constexpr int type = 9; // VTK_QUAD
        static inline int map(int i)
        {
          // bit-stunt: {0, 1, 2, 3} -> {0, 1, 3, 2}
          return (i ^ ((i >> 1) & 1));
        }
      };

      template<>
      struct VTKShape< Shape::Hypercube<3> >
      {
        static constexpr int type = 12; // VTK_HEXAHEDRON
        static inline int map(int i)
        {
          // bit-stunt: {0, 1, ..., 7} -> {0, 1, 3, 2, 4, 5, 7, 6}
          return (i ^ ((i >> 1) & 1));
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief VTK exporter class template
     *
     * This class templates implements an exporter for the XML-based VTK file formats.
     * This exporter is capable of writing (stand-alone) serial VTU files representing
     * unstructured grids as well as parallel PVTU files representing partitionings.
     *
     * \note This class template supports both the Geometry::ConformalMesh and
     * Geometry::StructuredMesh classes as input, however, both types of meshes are
     * exported as unstructured meshes in the sense of VTK.
     *
     * \tparam Mesh_
     * The type of the mesh to be exported.
     *
     * \author Peter Zajac
     */
    template<typename Mesh_>
    class ExportVTK
    {
    public:
      /// mesh type
      typedef Mesh_ MeshType;
      /// our shape type
      typedef typename MeshType::ShapeType ShapeType;
      /// our VTK shape type
      typedef Intern::VTKShape<ShapeType> VTKShapeType;

    protected:
      /// our variable container
      typedef std::deque<std::pair<String, std::vector<double>>> VarDeque;

      /// reference to mesh to be exported
      const MeshType& _mesh;
      /// number of vertices in mesh
      Index _num_verts;
      /// number of cells in mesh
      Index _num_cells;
      /// vertex variable list
      VarDeque _vertex_scalars;
      /// vertex field list
      VarDeque _vertex_vectors;
      /// cell variable list
      VarDeque _cell_scalars;
      /// cell field list
      VarDeque _cell_vectors;
      /// precision of variables
      int _var_prec;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] mesh
       * A reference to the mesh that is to be exported. Must remain unchanged for the lifetime of this exporter.
       *
       * \param[in] var_prec
       * Specifies the precision of the variable entries. If set to 0, the runtime's default precision is used.
       * See the \p precision parameter of the #stringify_fp_sci() function for details.
       */
      explicit ExportVTK(const MeshType& mesh, int var_prec = 0) :
        _mesh(mesh),
        _num_verts(mesh.get_num_entities(0)),
        _num_cells(mesh.get_num_entities(MeshType::shape_dim)),
        _var_prec(Math::max(0, var_prec))
      {
      }

      /// destructor
      virtual ~ExportVTK()
      {
      }

      /**
       * \brief Clears all vertex and cell variables in the exporter.
       */
      void clear()
      {
        _cell_scalars.clear();
        _vertex_scalars.clear();
        _vertex_vectors.clear();
        _cell_vectors.clear();
      }

      /**
       * \brief Adds a scalar vertex variable to the exporter.
       *
       * \param[in] name
       * The name of the variable to be exported.
       *
       * \param[in] data
       * An array of floating point values. Its length is assumed to correspond to the number of
       * vertices of the mesh. Must not be \p nullptr.
       *
       * \note
       * This function creates a (deep) copy of the data array, so the \p data array
       * can be deleted or overwritten after the return of this function.
       */
      template<typename T_>
      void add_vertex_scalar(const String& name, const T_* data)
      {
        ASSERT_(data != nullptr);
        std::vector<double> d(_num_verts);
        for(Index i(0); i < _num_verts; ++i)
        {
          d[i] = double(data[i]);
        }
        _vertex_scalars.push_back(std::make_pair(name, std::move(d)));
      }

      /**
       * \brief Adds a vector-field vertex variable to the exporter.
       *
       * This functions adds a 1D, 2D or 3D vector-field variable to the exporter.
       *
       * \param[in] name
       * The name of the variable to be exported.
       *
       * \param[in] x, y, z
       * The three arrays of floating point values. Their lengths are assumed to correspond to the number of
       * vertices of the mesh. The first array \p x must not be \p nullptr.
       *
       * \note
       * This function creates a (deep) copy of the three data arrays, so the data arrays
       * can be deleted or overwritten after the return of this function.
       */
      template<typename T_>
      void add_vertex_vector(const String& name, const T_* x, const T_* y = nullptr, const T_* z = nullptr)
      {
        ASSERT_(x != nullptr);
        std::vector<double> d(3*_num_verts);

        if(z != nullptr)
        {
          for(Index i(0); i < _num_verts; ++i)
          {
            d[3*i+0] = double(x[i]);
            d[3*i+1] = double(y[i]);
            d[3*i+2] = double(z[i]);
          }
        }
        else if(y != nullptr)
        {
          for(Index i(0); i < _num_verts; ++i)
          {
            d[3*i+0] = double(x[i]);
            d[3*i+1] = double(y[i]);
            d[3*i+2] = 0.0;
          }
        }
        else
        {
          for(Index i(0); i < _num_verts; ++i)
          {
            d[3*i+0] = double(x[i]);
            d[3*i+1] = 0.0;
            d[3*i+2] = 0.0;
          }
        }
        _vertex_vectors.push_back(std::make_pair(name, std::move(d)));
      }

      /**
       * \brief Adds a vector-field vertex variable given by a LAFEM::***VectorBlocked to the exporter.
       *
       * This functions adds a 1D, 2D or 3D vector-field variable to the exporter.
       *
       * \tparam VectorType_
       * Type of the input. Only DenseVectorBlocked and SparseVectorBlocked have the right interface at the moment.
       *
       * \param[in] name
       * The name of the variable to be exported.
       *
       * \param[in] v
       * The vector containing the field values in the mesh vertices.
       *
       * \note
       * This function creates a (deep) copy of the vector data.
       */
      template<typename VectorType_>
      void add_vertex_vector(const String& name, const VectorType_& v)
      {
        std::vector<double> d(3*_num_verts);

        for(Index i(0); i < _num_verts; ++i)
        {
          for(Index j(0); j < 3; ++j)
            d[Index(3)*i+j] = double(0);

          for(int j(0); j < v.BlockSize; ++j)
            d[Index(3)*i+Index(j)] = v(i)[j];

        }
        _vertex_vectors.push_back(std::make_pair(name, std::move(d)));
      }

      /**
       * \brief Adds a scalar cell variable to the exporter.
       *
       * \param[in] name
       * The name of the variable to be exported.
       *
       * \param[in] data
       * An array of floating point values. Its length is assumed to correspond to the number of
       * cells of the mesh. Must not be \p nullptr.
       *
       * \note
       * This function creates a (deep) copy of the data array, so the \p data array
       * can be deleted or overwritten after the return of this function.
       */
      template<typename T_>
      void add_cell_scalar(const String& name, const T_* data)
      {
        ASSERT_(data != nullptr);
        std::vector<double> d(_num_cells);
        for(Index i(0); i < _num_cells; ++i)
        {
          d[i] = double(data[i]);
        }
        _cell_scalars.push_back(std::make_pair(name, std::move(d)));
      }

      /**
       * \brief Adds a vector-field cell variable to the exporter.
       *
       * This functions adds a 1D, 2D or 3D vector-field variable to the exporter.
       *
       * \param[in] name
       * The name of the variable to be exported.
       *
       * \param[in] x, y, z
       * The three arrays of floating point values. Their lengths are assumed to correspond to the number of
       * cells of the mesh. The first array \p x must not be \p nullptr.
       *
       * \note
       * This function creates a (deep) copy of the three data arrays, so the data arrays
       * can be deleted or overwritten after the return of this function.
       */
      template<typename T_>
      void add_cell_vector(const String& name, const T_* x, const T_* y = nullptr, const T_* z = nullptr)
      {
        ASSERT_(x != nullptr);
        std::vector<double> d(3*_num_cells);

        if(z != nullptr)
        {
          for(Index i(0); i < _num_cells; ++i)
          {
            d[3*i+0] = double(x[i]);
            d[3*i+1] = double(y[i]);
            d[3*i+2] = double(z[i]);
          }
        }
        else if(y != nullptr)
        {
          for(Index i(0); i < _num_cells; ++i)
          {
            d[3*i+0] = double(x[i]);
            d[3*i+1] = double(y[i]);
            d[3*i+2] = 0.0;
          }
        }
        else
        {
          for(Index i(0); i < _num_cells; ++i)
          {
            d[3*i+0] = double(x[i]);
            d[3*i+1] = 0.0;
            d[3*i+2] = 0.0;
          }
        }
        _cell_vectors.push_back(std::make_pair(name, std::move(d)));
      }

      /**
       * \brief Adds a vector-field cell variable given by a LAFEM::***VectorBlocked to the exporter.
       *
       * This functions adds a 1D, 2D or 3D vector-field variable to the exporter.
       *
       * \tparam VectorType_
       * Type of the input. Only DenseVectorBlocked and SparseVectorBlocked have the right interface at the moment.
       *
       * \param[in] name
       * The name of the variable to be exported.
       *
       * \param[in] v
       * The vector containing the field values in the mesh cell.
       *
       * \note
       * This function creates a (deep) copy of the vector data.
       */
      template<typename VectorType_>
      void add_cell_vector(const String& name, const VectorType_& v)
      {
        std::vector<double> d(3*_num_cells);

        for(Index i(0); i < _num_cells; ++i)
        {
          for(Index j(0); j < 3; ++j)
            d[Index(3)*i+j] = double(0);

          for(int j(0); j < v.BlockSize; ++j)
            d[Index(3)*i+Index(j)] = v(i)[j];

        }
        _cell_vectors.push_back(std::make_pair(name, std::move(d)));
      }

      /**
       * \brief Writes out the data to a serial XML-VTU file.
       *
       * \param[in] filename
       * The filename to which to export to. The extension ".vtu" is automatically appended to the filename.
       */
      void write(const String& filename) const
      {
        // try to open the output file
        String vtu_name(filename + ".vtu");
        std::ofstream ofs(vtu_name.c_str());
        if(!(ofs.is_open() && ofs.good()))
          throw FileError("Failed to create '" + vtu_name + "'");

        // write
        write_vtu(ofs);

        // and close
        ofs.close();
      }

      /**
       * \brief Writes out the data to a parallel XML-PVTU file.
       *
       * This function writes out the serial VTU file of the partition represented by this exporter object
       * and, if the \p rank parameter is 0, also the corresponding parallel PVTU file.
       *
       * \param[in] filename
       * The filename to which to export to. The extension is automatically appended to the filename.
       *
       * \param[in] rank
       * The rank of this exporter. Must be 0 <= \p rank < \p nparts.
       * This usually coincides with the MPI-rank of the current process.
       *
       * \param[in] nparts
       * The total number of partitions. This usually coincides with the number of MPI-processes.
       *
       * \note
       * Specifying \p nparts < 1 is equivalent to calling <c>write(filename)</c>.
       */
      void write(const String& filename, const int rank, const int nparts)
      {
        // call the standard serial write version if nparts < 1
        if(nparts < 1)
        {
          write(filename);
          return;
        }

        // verify rank
        if((rank < 0) || (rank >= nparts))
          throw InternalError("Invalid rank '" + stringify(rank) + "'");

        // Add rank cell array since we're parallel if we come to here
        double* rank_array = new double[_num_cells];
        for(Index i(0); i < _num_cells; ++i)
          rank_array[i] = double(rank);

        add_cell_scalar("rank", rank_array);
        delete[] rank_array;

        // compute number of non-zero digits in (nparts-1) for padding
        const std::size_t ndigits = Math::ilog10(std::size_t(nparts-1));

        // write serial VTU file: "filename.#rank.vtu"
        write(filename + "." + stringify(rank).pad_front(ndigits, '0'));

        // we're done unless we have rank = 0
        if(rank != 0)
          return;

        // try to open our output file
        String pvtu_name(filename + ".pvtu");
        std::ofstream ofs(pvtu_name.c_str());
        if(!(ofs.is_open() && ofs.good()))
          throw FileError("Failed to create '" + pvtu_name + "'");

        // extract the file title from our filename
        std::size_t p = filename.find_last_of("\\/");
        String file_title = filename.substr(p == filename.npos ? 0 : ++p);

        // write PVTU file
        write_pvtu(ofs, file_title, nparts);

        // and close
        ofs.close();
      }

      void write_scheduled(const String& filename, const int rank, const int nparts)
      {
        auto func = [&] () { write_scheduled(filename, rank, nparts); };
        Util::schedule_function(func, Util::ScheduleMode::clustered);
      }

      /**
       * \brief Writes out the mesh and variable data in serial XML-VTU format.
       *
       * \param[in] os
       * The output stream to which to write to.
       */
      void write_vtu(std::ostream& os) const
      {
        // fetch basic infomation
        const int num_coords = MeshType::world_dim;
        const int verts_per_cell = Shape::FaceTraits<ShapeType,0>::count;

        // write VTK header
        os << "<VTKFile type=\"UnstructuredGrid\" version=\"0.1\">" << std::endl;
        os << "<!-- Generated by FEAT v" << version_major << "." << version_minor;
        os << "." << version_patch << " -->" << std::endl;

        // write mesh header
        os << "<UnstructuredGrid>" << std::endl;
        os << "<Piece NumberOfPoints=\"" << _num_verts << "\" NumberOfCells=\"" << _num_cells << "\">" << std::endl;

        // write point data
        if((!_vertex_scalars.empty()) || (!_vertex_vectors.empty()))
        {
          os << "<PointData>" << std::endl;

          // write vertex variables
          for(Index i(0); i < Index(_vertex_scalars.size()); ++i)
          {
            const auto& var(_vertex_scalars[i]);
            os << "<DataArray type=\"Float64\" Name=\"" << var.first <<"\" Format=\"ascii\">" << std::endl;
            for(Index j(0); j < _num_verts; ++j)
            {
              os << stringify_fp_sci(var.second[j], _var_prec) << std::endl;
            }
            os << "</DataArray>" << std::endl;
          }
          // write vertex fields
          for(Index i(0); i < Index(_vertex_vectors.size()); ++i)
          {
            const auto& var(_vertex_vectors[i]);
            os << "<DataArray type=\"Float64\" Name=\"" << var.first;
            os <<"\" NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;
            for(Index j(0); j < _num_verts; ++j)
            {
              os << stringify_fp_sci(var.second[3*j+0], _var_prec) << " ";
              os << stringify_fp_sci(var.second[3*j+1], _var_prec) << " ";
              os << stringify_fp_sci(var.second[3*j+2], _var_prec) << std::endl;
            }
            os << "</DataArray>" << std::endl;
          }

          os << "</PointData>" << std::endl;
        }

        // write cell data
        if(!_cell_scalars.empty() || !_cell_vectors.empty())
        {
          os << "<CellData>" << std::endl;
          if(!_cell_scalars.empty())
          {
            for(Index i(0); i < Index(_cell_scalars.size()); ++i)
            {
              const auto& var(_cell_scalars[i]);
              os << "<DataArray type=\"Float64\" Name=\"" << var.first <<"\" Format=\"ascii\">" << std::endl;
              for(Index j(0); j < _num_cells; ++j)
              {
                os << stringify_fp_sci(var.second[j], _var_prec) << std::endl;
              }
              os << "</DataArray>" << std::endl;
            }
          }

          if(!_cell_vectors.empty())
          {
            // write cell fields
            for(Index i(0); i < Index(_cell_vectors.size()); ++i)
            {
              const auto& var(_cell_vectors[i]);
              os << "<DataArray type=\"Float64\" Name=\"" << var.first;
              os <<"\" NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;
              for(Index j(0); j < _num_cells; ++j)
              {
                os << stringify_fp_sci(var.second[3*j+0], _var_prec) << " ";
                os << stringify_fp_sci(var.second[3*j+1], _var_prec) << " ";
                os << stringify_fp_sci(var.second[3*j+2], _var_prec) << std::endl;
              }
              os << "</DataArray>" << std::endl;
            }
          }
          os << "</CellData>" << std::endl;
        }

        // write vertices
        const auto& vtx = _mesh.get_vertex_set();
        os << "<Points>" << std::endl;
        os << "<DataArray type=\"Float32\" NumberOfComponents=\"3\" Format=\"ascii\">" << std::endl;
        for(Index i(0); i < _num_verts; ++i)
        {
          os << vtx[i][0];
          for(int j(1); j < num_coords; ++j)
          {
            os << " " << vtx[i][j];
          }
          for(int j(num_coords); j < 3; ++j)
          {
            os << " 0";
          }
          os << std::endl;
        }
        os << "</DataArray>" << std::endl;
        os << "</Points>" << std::endl;

        // write cells
        const auto& idx = _mesh.template get_index_set<MeshType::shape_dim, 0>();
        os << "<Cells>" << std::endl;
        os << "<DataArray type=\"UInt32\" Name=\"connectivity\">" << std::endl;
        for(Index i(0); i < _num_cells; ++i)
        {
          os << idx(i, Index(VTKShapeType::map(0)));
          for(int j(1); j < verts_per_cell; ++j)
          {
            os << " " << idx(i, Index(VTKShapeType::map(j)));
          }
          os << std::endl;
        }
        os << "</DataArray>" << std::endl;
        os << "<DataArray type=\"UInt32\" Name=\"offsets\">" << std::endl;
        for(Index i(0); i < _num_cells; ++i)
        {
          os << ((i+1) * verts_per_cell) << std::endl;
        }
        os << "</DataArray>" << std::endl;
        os << "<DataArray type=\"UInt32\" Name=\"types\">" << std::endl;
        for(Index i(0); i < _num_cells; ++i)
        {
          os << VTKShapeType::type << std::endl;
        }
        os << "</DataArray>" << std::endl;
        os << "</Cells>" << std::endl;

        // finish
        os << "</Piece>" << std::endl;
        os << "</UnstructuredGrid>" << std::endl;
        os << "</VTKFile>" << std::endl;
      }

      /**
       * \brief Writes out the partition data in parallel XML-PVTU format.
       *
       * \param[in] os
       * The output stream to which to write to.
       *
       * \param[in] file_title
       * The file title of the serial VTU files.
       *
       * \param[in] nparts
       * The total number of partitions.
       */
      void write_pvtu(std::ostream& os, const String& file_title, const int nparts) const
      {
        // write VTK header
        os << "<VTKFile type=\"PUnstructuredGrid\" version=\"0.1\">" << std::endl;
        os << "<!-- Generated by FEAT v" << version_major << "." << version_minor;
        os << "." << version_patch << " -->" << std::endl;
        os << "<PUnstructuredGrid GhostLevel=\"0\">" << std::endl;

        // write vertex data
        if((!_vertex_scalars.empty()) || (!_vertex_vectors.empty()))
        {
          os << "<PPointData>" << std::endl;

          // write vertex variables
          for(Index i(0); i < Index(_vertex_scalars.size()); ++i)
          {
            os << "<PDataArray type=\"Float64\" Name=\"" << _vertex_scalars[i].first <<"\" />" << std::endl;
          }
          // write vertex fields
          for(Index i(0); i < Index(_vertex_vectors.size()); ++i)
          {
            os << "<PDataArray type=\"Float64\" Name=\"" << _vertex_vectors[i].first;
            os <<"\" NumberOfComponents=\"3\" />" << std::endl;
          }

          os << "</PPointData>" << std::endl;
        }

        // write cell variables
        if(!_cell_scalars.empty() || !_cell_vectors.empty())
        {
          os << "<PCellData>" << std::endl;
          // write cell scalars
          for(Index i(0); i < Index(_cell_scalars.size()); ++i)
          {
            os << "<PDataArray type=\"Float64\" Name=\"" << _cell_scalars[i].first <<"\" />" << std::endl;
          }
          // write cell fields
          for(Index i(0); i < Index(_cell_vectors.size()); ++i)
          {
            os << "<PDataArray type=\"Float64\" Name=\"" << _cell_vectors[i].first;
            os <<"\" NumberOfComponents=\"3\" />" << std::endl;
          }
          os << "</PCellData>" << std::endl;
        }

        // write vertices
        os << "<PPoints>" << std::endl;
        os << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\" />" << std::endl;
        os << "</PPoints>" << std::endl;

        // compute number of non-zero digits in (nparts-1) for padding
        const std::size_t ndigits = Math::ilog10(std::size_t(nparts-1));

        // now let's write our piece data
        for(int i(0); i < nparts; ++i)
        {
          os << "<Piece Source=\"" << file_title << "." << stringify(i).pad_front(ndigits, '0');
          os << ".vtu\" />" << std::endl;
        }

        // finish
        os << "</PUnstructuredGrid>" << std::endl;
        os << "</VTKFile>" << std::endl;
      }
    }; // class ExportVTK
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_EXPORT_VTK_HPP
