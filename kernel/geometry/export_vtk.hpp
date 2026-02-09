// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/index_set.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/structured_mesh.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/dist_file_io.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/pack.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/tiny_algebra.hpp>

// includes, STL
#include <deque>
#include <fstream>
#include <functional>
#include <string_view>
#include <vector>

#ifdef FEAT_HAVE_CGAL
#include <kernel/geometry/cgal.hpp>
#endif


namespace FEAT::Geometry
{
  /// \cond internal
  namespace Intern
  {
    template<typename Shape_>
    struct VTKShape;

    template<>
    struct VTKShape<Shape::Simplex<1>>
    {
      static constexpr std::uint32_t type = 3; // VTK_LINE
      static inline int map(int i)
      {
        return i;
      }
    };

    template<>
    struct VTKShape<Shape::Simplex<2>>
    {
      static constexpr std::uint32_t type = 5; // VTK_TRIANGLE
      static inline int map(int i)
      {
        return i;
      }
    };

    template<>
    struct VTKShape<Shape::Simplex<3>>
    {
      static constexpr std::uint32_t type = 10; // VTK_TETRA
      static inline int map(int i)
      {
        return i;
      }
    };

    template<>
    struct VTKShape<Shape::Hypercube<1>>
    {
      static constexpr std::uint32_t type = 3; // VTK_LINE
      static inline int map(int i)
      {
        return i;
      }
    };

    template<>
    struct VTKShape<Shape::Hypercube<2>>
    {
      static constexpr std::uint32_t type = 9; // VTK_QUAD
      static inline int map(int i)
      {
        // bit-stunt: {0, 1, 2, 3} -> {0, 1, 3, 2}
        return (i ^ ((i >> 1) & 1));
      }
    };

    template<>
    struct VTKShape<Shape::Hypercube<3>>
    {
      static constexpr std::uint32_t type = 12; // VTK_HEXAHEDRON
      static inline int map(int i)
      {
        // bit-stunt: {0, 1, ..., 7} -> {0, 1, 3, 2, 4, 5, 7, 6}
        return (i ^ ((i >> 1) & 1));
      }
    };

    template<typename T_>
    struct VTKType;

    template<>
    struct VTKType<std::uint8_t>
    {
      static constexpr std::string_view type = "UInt8";
    };

    template<>
    struct VTKType<std::uint16_t>
    {
      static constexpr std::string_view type = "UInt16";
    };

    template<>
    struct VTKType<std::uint32_t>
    {
      static constexpr std::string_view type = "UInt32";
    };

    template<>
    struct VTKType<std::uint64_t>
    {
      static constexpr std::string_view type = "UInt64";
    };

    template<>
    struct VTKType<std::int8_t>
    {
      static constexpr std::string_view type = "Int8";
    };

    template<>
    struct VTKType<std::int16_t>
    {
      static constexpr std::string_view type = "Int16";
    };

    template<>
    struct VTKType<std::int32_t>
    {
      static constexpr std::string_view type = "Int32";
    };

    template<>
    struct VTKType<std::int64_t>
    {
      static constexpr std::string_view type = "Int64";
    };

    template<>
    struct VTKType<float>
    {
      static constexpr std::string_view type = "Float32";
    };

    template<>
    struct VTKType<double>
    {
      static constexpr std::string_view type = "Float64";
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
    * This class templates supports exporting in ascii, binary, and compressed
    * binary formats. The default format is binary, if zlib is not available,
    * and compressed binary, if zlib is available.
    *
    * \note This class template supports both the Geometry::ConformalMesh and
    * Geometry::StructuredMesh classes as input, however, both types of meshes are
    * exported as unstructured meshes in the sense of VTK.
    *
    * \tparam Mesh_
    * The type of the mesh to be exported.
    *
    * \tparam cell_dim_
    * The dimension of the mesh entities to be exported.
    *
    * \author Peter Zajac
    */
  template<typename Mesh_, int cell_dim_ = Mesh_::ShapeType::dimension>
  class ExportVTK
  {
  public:
    /// mesh type
    using MeshType = Mesh_;
    /// our shape type
    using ShapeType = typename MeshType::ShapeType;
    /// shape type of exported cells
    using CellShapeType = typename Shape::FaceTraits<ShapeType, cell_dim_>::ShapeType;
    /// our VTK shape type
    using VTKShapeType = Intern::VTKShape<CellShapeType>;

    /// Coordinates per vertex. The VTK documentation states that this is always three
    static constexpr int num_coords = 3;

    /// Vertices per cell
    static constexpr int verts_per_cell = Shape::FaceTraits<CellShapeType, 0>::count;

#ifdef FEAT_HAVE_ZLIB
    static constexpr bool can_compress = true;
#else
    static constexpr bool can_compress = false;
#endif

  protected:
    /// our variable container
    using VarDeque = std::deque<std::pair<String, std::vector<double>>>;

    // NOTE(mmuegge): Depending on wether we are writing a mesh or a meshpart,
    // we may need to do some additional handling of indices.
    // The accessor functions allow us to do so,
    // while presenting a consistent interface to the actual file writing logic.

    /// Accessor function for vertices
    std::function<void(Index, Tiny::Vector<float, num_coords>&)> _vertices;
    /// Accessor function for cells
    std::function<void(Index, IndexTuple<verts_per_cell>&)> _cells;
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
    /// ascii output flag
    bool _ascii_output = false;
    /// compressed output flag. Only applicable if _binary_output is true
    bool _use_compression = can_compress;

  public:
    /**
      * \brief Constructor
      *
      * \param[in] mesh
      * A \resident reference to the mesh that is to be exported.
      * Must remain unchanged for the lifetime of this exporter.
      *
      * \param[in] var_prec
      * Specifies the precision of the variable entries. If set to 0, the runtime's default precision is used.
      * See the \p precision parameter of the #stringify_fp_sci() function for details.
      */
    explicit ExportVTK(const MeshType& mesh, int var_prec = 0) :
      _num_verts(mesh.get_num_entities(0)),
      _num_cells(mesh.get_num_entities(cell_dim_)),
      _var_prec(Math::max(0, var_prec))
    {
      const auto& vertex_set = mesh.get_vertex_set();
      const auto& index_set = mesh.template get_index_set<cell_dim_, 0>();

      _vertices = [&](Index i, auto& vert) {
        // Actual vertex might not be 3D
        // vert is zeroed before this is called
        for(int j(0); j < vertex_set.num_coords; j++)
        {
          vert[j] = float(vertex_set[i][j]);
        }
      };
      _cells = [&](Index i, auto& tuple) {
        for(int j(0); j < verts_per_cell; j++)
        {
          tuple[j] = index_set(i, j);
        }
      };
    }

    /**
      * \brief Constructor
      *
      * \param[in] mesh
      * A \resident reference to the mesh that \c part was created from.
      * Must remain unchanged for the lifetime of this exporter.
      *
      * \param[in] part
      * A \resident reference to the meshpart that is to be exported.
      * Must remain unchanged for the lifetime of this exporter.
      *
      * \param[in] var_prec
      * Specifies the precision of the variable entries. If set to 0, the runtime's default precision is used.
      * See the \p precision parameter of the #stringify_fp_sci() function for details.
      */
    explicit ExportVTK(const MeshType& mesh, const MeshPart<MeshType>& part, int var_prec = 0) :
      _num_verts(part.get_num_entities(0)),
      _num_cells(part.get_num_entities(cell_dim_)),
      _var_prec(Math::max(0, var_prec))
    {
      const auto& vertex_set = mesh.get_vertex_set();
      const auto& vertex_target_set = part.template get_target_set<0>();
      const auto& cell_target_set = part.template get_target_set<cell_dim_>();

      _vertices = [&](Index i, auto& vert) {
        // Actual vertex might not be 3D
        // vert is zeroed before this is called
        for(int j(0); j < vertex_set.num_coords; j++)
        {
          vert[j] = float(vertex_set[i][j]);
        }
      };

      if(part.has_topology())
      {
        // The meshpart has its own topology. Use that.
        const auto& index_set = part.template get_index_set<cell_dim_, 0>();
        _cells = [&](Index i, auto& tuple) {
          for(int j(0); j < verts_per_cell; j++)
          {
            tuple[j] = index_set(i, j);
          }
        };
      }
      else
      {
        // The meshpart has no own topology. Use the mesh topology.
        const auto& index_set = mesh.template get_index_set<cell_dim_, 0>();

        // The mesh index-set returns indices of vertices the mesh's index-space,
        // but the .vtu will be written in the meshparts index-space.
        // We thus need to search for the corresponding index on the meshpart for any vertex index.
        _cells = [&](Index i, auto& tuple)
        {
          for(int j(0); j < verts_per_cell; j++)
          {
            const Index vertex = index_set(cell_target_set[i], j);
            for(Index k(0); k < vertex_target_set.get_num_entities(); k++)
            {
              if(vertex_target_set[k] == vertex)
              {
                tuple[j] = k;
              }
            }
          }
        };
      }
    }

#if defined(FEAT_HAVE_CGAL) || defined(DOXYGEN)
    /**
      * \brief Constructor
      *
      * \param[in] cw
      * A \resident reference to the cgal surface mesh.
      * Must remain unchanged for the lifetime of this exporter.
      *
      * \param[in] var_prec
      * Specifies the precision of the variable entries. If set to 0, the runtime's default precision is used.
      * See the \p precision parameter of the #stringify_fp_sci() function for details.
      */
    template<typename DT_>
    explicit ExportVTK(const Geometry::CGALWrapper<DT_>& cw, int var_prec = 0) :
      _num_verts(cw.get_num_entities(0)),
      _num_cells(cw.get_num_entities(cell_dim_)),
      _var_prec(Math::max(0, var_prec))
    {
      XASSERTM(cell_dim_ == 2, "CGAL Export only possible with surface mesh");

      _vertices = [&](Index i, auto& vert) {
        auto point = cw.point(i);
        for(int j(0); j < verts_per_cell; j++)
        {
          vert[j] = float(point[j]);
        }
      };
      _cells = [&](Index i, auto& tuple) {
        auto it = cw.vertices_around_face().image_begin(i);
        auto end = cw.vertices_around_face().image_end(i);
        Index j(0);
        for(;it != end; ++it)
        {
          tuple[j++] = *it;
        }
      };
    }
#endif

    /// destructor
    virtual ~ExportVTK() = default;

    /// Copy constructor
    ExportVTK(const ExportVTK& other) = default;

    /// Move constructor
    ExportVTK(ExportVTK&& other) noexcept = default;

    /// Copy-assignment constructor
    ExportVTK& operator=(const ExportVTK& other) = default;

    /// Move-assignment constructor
    ExportVTK& operator=(ExportVTK&& other) = default;

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
     * \brief Enable or disable binary output
     *
     * \param[in] flag Boolean flag indicating whether binary output should be enabled
     */
    void enable_ascii_output(bool flag)
    {
      _ascii_output = flag;
    }

    /**
     * \brief Enable or disable compression
     *
     * \param[in] flag Boolean flag indicating whether compression should be enabled
     *
     * \note Compression is only supported for binary output and only used if FEAT3 has been built with zlib support.
     */
    void enable_compression(bool flag)
    {
      _use_compression = flag;
    }

    /**
      * \brief Adds a scalar vertex variable to the exporter.
      *
      * \param[in] name
      * The name of the variable to be exported.
      *
      * \param[in] data
      * An \transient array of floating point values. Its length is assumed to correspond to the number of
      * vertices of the mesh. Must not be \p nullptr.
      *
      * \param[in] scaling_factor
      * A scaling factor applied to each datapoint.
      *
      * \note
      * This function creates a (deep) copy of the data array, so the \p data array
      * can be deleted or overwritten after the return of this function.
      */
    template<typename T_>
    void add_vertex_scalar(const String& name, const T_* data, double scaling_factor = 1.0)
    {
      XASSERTM(data != nullptr, "data array is nullptr");
      std::vector<double> d(_num_verts);
      for(Index i(0); i < _num_verts; ++i)
      {
        d[i] = scaling_factor * double(data[i]);
      }
      _vertex_scalars.emplace_back(name, std::move(d));
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
      * The three \transient arrays of floating point values. Their lengths are assumed to correspond
      * to the number of vertices of the mesh. The first array \p x must not be \p nullptr.
      *
      * \param[in] scaling_factor
      * A scaling factor applied to each datapoint.
      *
      * \note
      * This function creates a (deep) copy of the three data arrays, so the data arrays
      * can be deleted or overwritten after the return of this function.
      */
    template<typename T_>
    void add_vertex_vector(
      const String& name,
      const T_* x,
      const T_* y = nullptr,
      const T_* z = nullptr,
      double scaling_factor = 1.0)
    {
      XASSERTM(x != nullptr, "x-data array is nullptr");
      std::vector<double> d(3 * _num_verts);

      if(y != nullptr && z != nullptr)
      {
        for(Index i(0); i < _num_verts; ++i)
        {
          d[(3 * i) + 0] = scaling_factor * double(x[i]);
          d[(3 * i) + 1] = scaling_factor * double(y[i]);
          d[(3 * i) + 2] = scaling_factor * double(z[i]);
        }
      }
      else if(y != nullptr)
      {
        for(Index i(0); i < _num_verts; ++i)
        {
          d[(3 * i) + 0] = scaling_factor * double(x[i]);
          d[(3 * i) + 1] = scaling_factor * double(y[i]);
          d[(3 * i) + 2] = 0.0;
        }
      }
      else
      {
        for(Index i(0); i < _num_verts; ++i)
        {
          d[(3 * i) + 0] = scaling_factor * double(x[i]);
          d[(3 * i) + 1] = 0.0;
          d[(3 * i) + 2] = 0.0;
        }
      }
      _vertex_vectors.emplace_back(name, std::move(d));
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
      * A \transient reference to the vector containing the field values in the mesh vertices.
      *
      * \param[in] scaling_factor
      * A scaling factor applied to each datapoint.
      *
      * \note
      * This function creates a (deep) copy of the vector data.
      */
    template<typename VectorType_>
    void add_vertex_vector(const String& name, const VectorType_& v, double scaling_factor = 1.0)
    {
      std::vector<double> d(3 * _num_verts);

      for(Index i(0); i < _num_verts; ++i)
      {
        for(Index j(0); j < 3; ++j)
        {
          d[(Index(3) * i) + j] = double(0);
        }

        for(int j(0); j < v.BlockSize; ++j)
        {
          d[(Index(3) * i) + Index(j)] = scaling_factor * double(v(i)[j]);
        }
      }
      _vertex_vectors.emplace_back(name, std::move(d));
    }

    /**
      * \brief Adds a vector-field vertex variable given by a std::vector of Tiny::Vector to the exporter.
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
      * A \transient reference to the vector containing the field values in the mesh vertices.
      *
      * \param[in] scaling_factor
      * A scaling factor applied to each datapoint.
      *
      * \note
      * This function creates a (deep) copy of the vector data.
      */
    template<typename DT_, int dim_>
    void add_vertex_vector(const String& name, const std::vector<Tiny::Vector<DT_, dim_>>& v, double scaling_factor = 1.0)
    {
      std::vector<double> d(3 * _num_verts);

      for(Index i(0); i < _num_verts; ++i)
      {
        for(Index j(0); j < 3; ++j)
        {
          d[(Index(3) * i) + j] = double(0);
        }

        for(int j(0); j < dim_; ++j)
        {
          d[(Index(3) * i) + Index(j)] = scaling_factor * double(v[i][j]);
        }
      }
      _vertex_vectors.emplace_back(name, std::move(d));
    }

    /**
      * \brief Adds a scalar cell variable to the exporter.
      *
      * \param[in] name
      * The name of the variable to be exported.
      *
      * \param[in] data
      * An \transient array of floating point values. Its length is assumed to correspond to the number of
      * cells of the mesh. Must not be \p nullptr.
      *
      * \param[in] scaling_factor
      * A scaling factor applied to each datapoint.
      *
      * \note
      * This function creates a (deep) copy of the data array, so the \p data array
      * can be deleted or overwritten after the return of this function.
      */
    template<typename T_>
    void add_cell_scalar(const String& name, const T_* data, double scaling_factor = 1.0)
    {
      XASSERTM(data != nullptr, "data array is nullptr");
      std::vector<double> d(_num_cells);
      for(Index i(0); i < _num_cells; ++i)
      {
        d[i] = scaling_factor * double(data[i]);
      }
      _cell_scalars.emplace_back(name, std::move(d));
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
      * The three \transient arrays of floating point values. Their lengths are assumed to correspond
      * to the number of cells of the mesh. The first array \p x must not be \p nullptr.
      *
      * \param[in] scaling_factor
      * A scaling factor applied to each datapoint.
      *
      * \note
      * This function creates a (deep) copy of the three data arrays, so the data arrays
      * can be deleted or overwritten after the return of this function.
      */
    template<typename T_>
    void add_cell_vector(
      const String& name,
      const T_* x,
      const T_* y = nullptr,
      const T_* z = nullptr,
      double scaling_factor = 1.0)
    {
      XASSERTM(x != nullptr, "x-data array is nullptr");
      std::vector<double> d(3 * _num_cells);

      if(y != nullptr && z != nullptr)
      {
        for(Index i(0); i < _num_cells; ++i)
        {
          d[(3 * i) + 0] = scaling_factor * double(x[i]);
          d[(3 * i) + 1] = scaling_factor * double(y[i]);
          d[(3 * i) + 2] = scaling_factor * double(z[i]);
        }
      }
      else if(y != nullptr)
      {
        for(Index i(0); i < _num_cells; ++i)
        {
          d[(3 * i) + 0] = scaling_factor * double(x[i]);
          d[(3 * i) + 1] = scaling_factor * double(y[i]);
          d[(3 * i) + 2] = 0.0;
        }
      }
      else
      {
        for(Index i(0); i < _num_cells; ++i)
        {
          d[(3 * i) + 0] = scaling_factor * double(x[i]);
          d[(3 * i) + 1] = 0.0;
          d[(3 * i) + 2] = 0.0;
        }
      }
      _cell_vectors.emplace_back(name, std::move(d));
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
      * A \transient reference to the vector containing the field values in the mesh cell.
      *
      * \param[in] scaling_factor
      * A scaling factor applied to each datapoint.
      *
      * \note
      * This function creates a (deep) copy of the vector data.
      */
    template<typename VectorType_>
    void add_cell_vector(const String& name, const VectorType_& v, double scaling_factor = 1.0)
    {
      std::vector<double> d(3 * _num_cells);

      for(Index i(0); i < _num_cells; ++i)
      {
        for(Index j(0); j < 3; ++j)
        {
          d[(Index(3) * i) + j] = double(0);
        }

        for(int j(0); j < v.BlockSize; ++j)
        {
          d[(Index(3) * i) + Index(j)] = scaling_factor * double(v(i)[j]);
        }
      }
      _cell_vectors.emplace_back(name, std::move(d));
    }

    /**
      * \brief Adds a vector-field cell variable given by a std::vector of Tiny::Vector to the exporter.
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
      * A \transient reference to the vector containing the field values in the mesh cell.
      *
      * \param[in] scaling_factor
      * A scaling factor applied to each datapoint.
      *
      * \note
      * This function creates a (deep) copy of the vector data.
      */
    template<typename DT_, int dim_>
    void add_cell_vector(const String& name, const std::vector<Tiny::Vector<DT_, dim_>>& v, double scaling_factor = 1.0)
    {
      std::vector<double> d(3 * _num_cells);

      for(Index i(0); i < _num_cells; ++i)
      {
        for(Index j(0); j < 3; ++j)
        {
          d[(Index(3) * i) + j] = double(0);
        }

        for(int j(0); j < dim_; ++j)
        {
          d[(Index(3) * i) + Index(j)] = scaling_factor * double(v[i][j]);
        }
      }
      _cell_vectors.emplace_back(name, std::move(d));
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
      {
        throw FileError("Failed to create '" + vtu_name + "'");
      }

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
      if(nparts <= 1)
      {
        write(filename);
        return;
      }

      // verify rank
      XASSERTM((rank >= 0) && (rank < nparts), "Invalid rank");

      // Add rank cell array since we're parallel if we come to here
      std::vector<double> rank_array((std::size_t(_num_cells)), double(rank));
      add_cell_scalar("rank", rank_array.data());

      // compute number of non-zero digits in (nparts-1) for padding
      const std::size_t ndigits = Math::ilog10(std::size_t(nparts - 1));

      // write serial VTU file: "filename.#rank.vtu"
      write(filename + "." + stringify(rank).pad_front(ndigits, '0'));

      // we're done unless we have rank = 0
      if(rank != 0)
      {
        return;
      }

      // try to open our output file
      String pvtu_name(filename + ".pvtu");
      std::ofstream ofs(pvtu_name.c_str());
      if(!(ofs.is_open() && ofs.good()))
      {
        throw FileError("Failed to create '" + pvtu_name + "'");
      }

      // extract the file title from our filename
      std::size_t p = filename.find_last_of("\\/");
      String file_title = filename.substr(p == FEAT::String::npos ? 0 : ++p);

      // write PVTU file
      write_pvtu(ofs, file_title, nparts);

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
      * \param[in] comm
      * The communicator of this exporter.
      */
    void write(const String& filename, const Dist::Comm& comm)
    {
      // Add rank cell array since we're parallel if we come to here
      std::vector<double> rank_array(std::size_t(_num_cells), double(comm.rank()));
      add_cell_scalar("rank", rank_array.data());

      // compute number of non-zero digits in (nparts-1) for padding
      const auto ndigits = std::size_t(Math::max(Math::ilog10(comm.size() - 1), 1));

      // write serial VTU file into a stringstream
      std::stringstream stream;
      write_vtu(stream);

      // generate pattern for filename: "filename.#rank.vtu"
      String pattern = filename + "." + String(ndigits, '*') + ".vtu";

      // write distributed VTU files
      DistFileIO::write_sequence(stream, pattern, comm);

      // we're done unless we have rank = 0
      if(comm.rank() != 0)
      {
        return;
      }

      // try to open our output file
      String pvtu_name(filename + ".pvtu");
      std::ofstream ofs(pvtu_name.c_str());
      if(!(ofs.is_open() && ofs.good()))
      {
        throw FileError("Failed to create '" + pvtu_name + "'");
      }

      // extract the file title from our filename
      std::size_t p = filename.find_last_of("\\/");
      String file_title = filename.substr(p == FEAT::String::npos ? 0 : ++p);

      // write PVTU file
      write_pvtu(ofs, file_title, comm.size());

      // and close
      ofs.close();
    }

    /**
      * \brief Writes out the mesh and variable data in serial XML-VTU format.
      *
      * \param[in] os
      * The output stream to which to write to.
      */
    void write_vtu(std::ostream& os) const
    {
      // fetch basic information

      // write VTK header
      os << R"(<VTKFile type="UnstructuredGrid" version="0.1" byte_order=")";
      os << (_little_endian() ? "LittleEndian" : "BigEndian") << "\"";

      if(!_ascii_output && _use_compression && can_compress)
      {
#ifndef FEAT_HAVE_ZLIB
        XABORTM("Requested compression in ExportVTK, but ZLib is not available!");
#endif
        os << " compressor=\"vtkZLibDataCompressor\"";
      }

      os << ">\n";
      os << "<!-- Generated by FEAT v" << version_major << "." << version_minor;
      os << "." << version_patch << " -->\n";

      // write mesh header
      os << "<UnstructuredGrid>\n";
      os << "<Piece NumberOfPoints=\"" << _num_verts << "\" NumberOfCells=\"" << _num_cells << "\">\n";

      // write point data
      if((!_vertex_scalars.empty()) || (!_vertex_vectors.empty()))
      {
        os << "<PointData>\n";
        for(const auto& var : _vertex_scalars)
        {
          _write_data_array(os, var, 1);
        }

        for(const auto& var : _vertex_vectors)
        {
          _write_data_array(os, var, 3);
        }
        os << "</PointData>\n";
      }

      // write cell data
      if(!_cell_scalars.empty() || !_cell_vectors.empty())
      {
        os << "<CellData>\n";

        for(const auto& var : _cell_scalars)
        {
          _write_data_array(os, var, 1);
        }

        for(const auto& var : _cell_vectors)
        {
          _write_data_array(os, var, 3);
        }

        os << "</CellData>\n";
      }

      // write vertices
      os << "<Points>\n";
      if(_ascii_output)
      {
        // Ascii output. Write data array adhoc without copying all vertices.
        os << R"(<DataArray type="Float32" NumberOfComponents="3" format="ascii" Name="test">)" << "\n";
        Tiny::Vector<float, num_coords> vertex;
        for(Index i(0); i < _num_verts; ++i)
        {
          vertex.format(0);
          _vertices(i, vertex);
          os << vertex[0];
          for(int j(1); j < num_coords; ++j)
          {
            os << " " << vertex[j];
          }
          os << "\n";
        }
        os << "</DataArray>\n";
      }
      else
      {
        // Binary output. Copy vertices to buffer for base64 encoding
        // This might not be needed, but ensures we are independent of the internal representation of the vertices
        std::pair<String, std::vector<float>> vertices("", {});
        vertices.second.reserve(3 * _num_cells);

        Tiny::Vector<float, num_coords> vertex;
        for(Index i(0); i < _num_verts; ++i)
        {
          vertex.format(0);
          _vertices(i, vertex);
          for(int j(0); j < num_coords; ++j)
          {
            vertices.second.push_back(static_cast<float>(vertex[j]));
          }
        }

        _write_data_array(os, vertices, num_coords);
      }
      os << "</Points>\n";

      // write cells
      os << "<Cells>\n";
      if(_ascii_output)
      {
        os << R"(<DataArray type="UInt32" Name="connectivity">)" << "\n";
        IndexTuple<verts_per_cell> cell;
        for(Index i(0); i < _num_cells; ++i)
        {
          _cells(i, cell);
          os << cell[VTKShapeType::map(0)];
          for(int j(1); j < verts_per_cell; ++j)
          {
            os << " " << cell[VTKShapeType::map(j)];
          }
          os << "\n";
        }
        os << "</DataArray>\n";
      }
      else
      {
        // Same as for vertices. Copy for base64 encoding
        std::pair<String, std::vector<std::uint32_t>> connectivity("connectivity", {});
        connectivity.second.reserve(_num_cells * verts_per_cell);

        IndexTuple<verts_per_cell> cell;
        for(Index i(0); i < _num_cells; ++i)
        {
          _cells(i, cell);
          for(int j(0); j < verts_per_cell; ++j)
          {
            connectivity.second.push_back(static_cast<std::uint32_t>(cell[VTKShapeType::map(j)]));
          }
        }

        _write_data_array(os, connectivity);
      }

      if(_ascii_output)
      {
        os << R"(<DataArray type="UInt32" Name="offsets">)" << "\n";
        for(Index i(0); i < _num_cells; ++i)
        {
          os << ((i + 1) * verts_per_cell) << "\n";
        }
        os << "</DataArray>\n";
      }
      else
      {
        // Same as for vertices. Copy for base64 encoding
        std::pair<String, std::vector<std::uint32_t>> offsets("offsets", {});
        offsets.second.reserve(_num_cells);

        for(Index i(0); i < _num_cells; ++i)
        {
          offsets.second.push_back(static_cast<std::uint32_t>((i + 1) * verts_per_cell));
        }

        _write_data_array(os, offsets);
      }

      if(_ascii_output)
      {
        os << R"(<DataArray type="UInt32" Name="types">)" << "\n";
        for(Index i(0); i < _num_cells; ++i)
        {
          os << VTKShapeType::type << "\n";
        }
        os << "</DataArray>\n";
      }
      else
      {
        // Same as for vertices. Copy for base64 encoding
        std::pair<String, std::vector<std::uint32_t>> types("types", {});
        types.second.reserve(_num_cells);

        for(Index i(0); i < _num_cells; ++i)
        {
          types.second.push_back(VTKShapeType::type);
        }

        _write_data_array(os, types);
      }
      os << "</Cells>\n";

      // finish
      os << "</Piece>\n";
      os << "</UnstructuredGrid>\n";
      os << "</VTKFile>\n";
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
      os << R"(<VTKFile type="PUnstructuredGrid" version="0.1" byte_order=")";
      os << (_little_endian() ? "LittleEndian" : "BigEndian") << "\"";

      if(!_ascii_output && _use_compression && can_compress)
      {
#ifndef FEAT_HAVE_ZLIB
        XABORTM("Requested compression in ExportVTK, but ZLib is not available!");
#endif
        os << " compressor=\"vtkZLibDataCompressor\"";
      }

      os << ">\n";
      os << "<!-- Generated by FEAT v" << version_major << "." << version_minor;
      os << "." << version_patch << " -->\n";
      os << "<PUnstructuredGrid GhostLevel=\"0\">\n";

      // write vertex data
      if((!_vertex_scalars.empty()) || (!_vertex_vectors.empty()))
      {
        os << "<PPointData>\n";

        // write vertex variables
        for(const auto& vertex_scalar : _vertex_scalars)
        {
          os << "<PDataArray type=\"Float64\" Name=\"" << vertex_scalar.first << "\" />\n";
        }
        // write vertex fields
        for(const auto& vertex_vector : _vertex_vectors)
        {
          os << "<PDataArray type=\"Float64\" Name=\"" << vertex_vector.first;
          os << "\" NumberOfComponents=\"3\" />\n";
        }

        os << "</PPointData>\n";
      }

      // write cell variables
      if(!_cell_scalars.empty() || !_cell_vectors.empty())
      {
        os << "<PCellData>\n";
        // write cell scalars
        for(const auto& cell_scalar : _cell_scalars)
        {
          os << "<PDataArray type=\"Float64\" Name=\"" << cell_scalar.first << "\" />\n";
        }
        // write cell fields
        for(const auto& cell_vector : _cell_vectors)
        {
          os << "<PDataArray type=\"Float64\" Name=\"" << cell_vector.first;
          os << "\" NumberOfComponents=\"3\" />\n";
        }
        os << "</PCellData>\n";
      }

      // write vertices
      os << "<PPoints>\n";
      os << "<PDataArray type=\"Float32\" NumberOfComponents=\"3\" />\n";
      os << "</PPoints>\n";

      // compute number of non-zero digits in (nparts-1) for padding
      const std::size_t ndigits = Math::ilog10(std::size_t(nparts - 1));

      // now let's write our piece data
      for(int i(0); i < nparts; ++i)
      {
        os << "<Piece Source=\"" << file_title << "." << stringify(i).pad_front(ndigits, '0');
        os << ".vtu\" />\n";
      }

      // finish
      os << "</PUnstructuredGrid>\n";
      os << "</VTKFile>\n";
    }
  private:

    /**
     * \brief Header for uncompressed binary data
     *
     * Note that this header is mentioned nowhere in the VTK documentation, but is definitely required.
     * For uncompressed data it contains the number of bytes of the following data (excluding the header itself)
     * stored as a unsigned 32bit integer.
     */
    struct UncompressedHeader
    {
      std::uint32_t num_bytes;

      explicit UncompressedHeader(std::uint32_t number_of_bytes) : num_bytes(number_of_bytes)
      {
      }
    };

    /**
     * \brief Header for compressed binary data
     *
     * See https://mathema.tician.de/what-they-dont-tell-you-about-vtk-xml-binary-formats/
     */
    struct CompressedHeader
    {
      std::uint32_t blocks = 1;
      std::uint32_t blocksize;
      std::uint32_t last_blocksize;
      std::uint32_t compressed_blocksizes[1];

      explicit CompressedHeader(std::uint32_t uncompressed_bytes, std::uint32_t compressed_bytes) :
        blocksize(uncompressed_bytes),
        last_blocksize(uncompressed_bytes),
        compressed_blocksizes{compressed_bytes}
      {
      }
    };

    /**
     * \brief Writes an DataArray XML-node to the give output stream
     *
     * \tparam T_ Type of data to write, supported are std::uint{8,16,32,64}_t, std::int{8,16,32,64}_t, float, and double.
     *
     * \param[in] os Outputstream to write to
     * \param[in] data Pair of array name and data to write
     * \param[in] number_of_components Number of components per entry, i.e. 1 for scalar arrays, 3 for 3D vector arrays
     *
     * If the name of the data array is empty, no Name attribute will be added to the emitted DataArray tag
     * If \c number_of_components is zero, no NumberOfComponents attribute will be added to the emitted DataArray tag
     *
     * \c data is automatically compressed and base64 encoded, depending on settings for binary output and compression.
     */
    template<typename T_>
    void _write_data_array(std::ostream& os, const std::pair<String, std::vector<T_>>& data, Index number_of_components=0) const
    {
      // Ensure there is no padding in headers.
      // Otherwise the reinterpret_casts below would produce garbage data
      static_assert(sizeof(UncompressedHeader) == sizeof(std::uint32_t));
      static_assert(sizeof(CompressedHeader) == 4 * sizeof(std::uint32_t));

      std::string_view format = _ascii_output ? "ascii" : "binary";
      std::string_view type = Intern::VTKType<T_>::type;

      // Write XML-tag
      os << "<DataArray type=\"" << type << "\" format=\"" << format << "\"";
      if(data.first != "")
      {
        os << " Name = \"" << data.first << "\"";
      }
      if(number_of_components != 0)
      {
       os << " NumberOfComponents=\"" << number_of_components << "\"";
      }
      os << ">\n";

      if(_ascii_output)
      {
        Index num_blocks = data.second.size() / number_of_components;
        for(Index i(0); i < num_blocks; ++i)
        {
          for(Index j(0); j < number_of_components; j++)
          {
            os << stringify_fp_sci(data.second[(i * number_of_components) + j], _var_prec);
            os << ((j == number_of_components - 1) ? "\n" : " ");
          }
        }
      }
      else if(_use_compression && can_compress)
      {
        auto uncompressed_size = static_cast<std::uint32_t>(data.second.size() * sizeof(T_));

        // Encode with compression
        Pack::Type pack_type = Pack::deduct_type<T_>() | Pack::Type::Mask_Z;
        std::size_t buffer_size = Pack::estimate_size(data.second.size(), pack_type);

        std::vector<Pack::u8> buffer(buffer_size, 0);
        std::size_t compressed_size = Pack::encode(buffer.data(), data.second.data(), buffer_size, data.second.size(), pack_type, false);

        // Set up header and data pointers
        CompressedHeader header(uncompressed_size, static_cast<std::uint32_t>(compressed_size));
        const auto* header_begin = reinterpret_cast<const Pack::u8*>(&header);
        const auto* header_end = header_begin + sizeof(CompressedHeader);

        const auto* data_begin = reinterpret_cast<const Pack::u8*>(buffer.data());
        const auto* data_end = data_begin + compressed_size;

        // For compressed data, the header and the data are encoded separately
        os << Pack::base64_encode(header_begin, header_end) << Pack::base64_encode(data_begin, data_end) << "\n";
      }
      else
      {
        UncompressedHeader header(static_cast<std::uint32_t>(data.second.size() * sizeof(T_)));

        const auto* header_begin = reinterpret_cast<const Pack::u8*>(&header);
        const auto* header_end = header_begin + sizeof(header);

        const auto* data_begin = reinterpret_cast<const Pack::u8*>(data.second.data());
        const auto* data_end = data_begin + header.num_bytes;

        os << Pack::base64_encode(header_begin, header_end) << Pack::base64_encode(data_begin, data_end) << "\n";
      }

      // Close XML-tag
      os << "</DataArray>\n";
    }

    // Copied from binary_stream-test.cpp
    bool _little_endian() const
    {
      std::int32_t x = 1;
      void* xv = (void*)&x;
      auto* x16 = (int16_t*)xv;
      return *x16 == 1;
    }
  }; // class ExportVTK
} // namespace FEAT::Geometry
