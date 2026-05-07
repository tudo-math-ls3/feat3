// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <cstring>
#include <filesystem>
#include <fstream>
#include <ios>
#include <istream>
#include <string_view>

#include <kernel/base_header.hpp>
#include <kernel/geometry/io_common.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/shape.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT::Geometry
{
  /**
   * \brief Gmsh file reader class
   *
   * This class implements a reader which parses Gmesh's .msh files into
   * a RootMeshNode.
   *
   * The basic usage is as follows:
   * -# Create an object of the GmshFileReader
   * -# Scan through the file to ensure it contains a mesh that FEAT can represent using read_root_markup()
   * -# Call the parse() function with the appropriate mesh type
   *
   * For more details on meshes, see the related doxygen page \ref mesh_file_format.
   *
   * Limitations:
   * -# Currently only nodes and elements are parsed. No charts, meshparts, or partitions will be created.
   * -# Only the version 2.2 legacy format is supported.
   * -# Only meshes with contiguous indices for vertices and elements are supported.
   *
   * See Chapter 10 of https://gmsh.info/dev/doc/texinfo/gmsh.pdf for detail on the mesh format.
   *
   * \author Markus Muegge
   */
  class GmshFileReader
  {
    /// Utility class for walking through potentially mixed ascii and binary .msh files
    struct FileView
    {
      /// End of buffer (for bounds checks)
      const char* end;
      /// Current position in buffer
      const char* current;

      /// Constructor
      explicit FileView(const String& bytes) : FileView(bytes, bytes.data())
      {
      }

      /// Constructor
      FileView(const String& bytes, const char* pos) :
        end(bytes.data() + bytes.size()),
        current(pos)
      {
      }

      /// Current position in buffer
      const char* position() const
      {
        return current;
      }

      bool done() const
      {
        return current == end;
      }

      /// Return a view of the data up to (but excluding) the next '\n' char, then advance the view beyond that '\n'
      /// char
      std::string_view next_line()
      {
        // Calling this at the end of the buffer is an error
        XASSERT(current < end);

        // If this is called, we assume that the next part of the data is valid ascii
        // We can thus safely skip to the next non-whitespace char
        while(current < end && (std::isspace(static_cast<unsigned char>(*current)) != 0))
        {
          current++;
        }

        // Search for next newline
        const char* eol = static_cast<const char*>(std::memchr(current, '\n', std::size_t(end - current)));

        if(eol == nullptr)
        {
          eol = end;
        }

        // String views are constructed via a (start, len) pair
        // So convert from end-of-line pointer to length
        const auto line_len = std::size_t(eol - current);

        const std::string_view line{current, line_len};

        // Advance past line. +1 to skip over newline
        current = eol + 1 < end ? eol + 1 : end;

        // Check we are still in bounds
        XASSERT(current <= end);

        return line;
      }

      /// Advance the view \c len chars
      void advance(std::size_t len)
      {
        XASSERT(current + len <= end);
        current += len;
      }

      /// Parses a line (split by whitespace) into the given parameters
      template<typename... Args>
      void parse_line(Args&... args)
      {
        const std::string_view line = next_line();
        std::istringstream stream((std::string(line)));

        ((stream >> args), ...);

        XASSERT(!stream.fail());
      }

      /// Copy a value from the buffer and then advance that many bytes
      template<typename T>
      void memcpy(T* t)
      {
        XASSERT(current + sizeof(T) <= end);
        std::memcpy(t, current, sizeof(T));
        advance(sizeof(T));
      }
    };

    /// Header for a chunk of elements in binary .msh files
    struct BinaryElementsHeader
    {
      /// Type of elements in chunk
      int type;
      /// Number of elements in chunk
      int num_elements;
      /// Tags per element of the chunk
      int tags_per_element;

      /// Start of elements of chunk in mesh (excluding the chunks header)
      const char* elements_start;
    };

    /// Contents of the .msh file
    String _file_contents;

    bool _have_root_markup = false;

    // The following is all the data that gets parsed from the file
    // during read_root_markup

    /// .msh file version
    double _version = 0.0;
    /// .msh ascii/binary indicator
    bool _is_binary = false;
    /// .msh size for reals
    std::size_t _data_size = 8;
    /// Number of nodes in the file
    std::size_t _num_nodes = 0;
    /// Number of elements in the file
    std::size_t _num_elements = 0;
    /// .msh element type
    int _element_type = -1;
    /// Start of nodes in the file
    const char* _nodes_start = nullptr;
    /// Start of elements in the file, for ascii files
    const char* _elements_start = nullptr;
    /// Element chunks, for binary files
    std::vector<BinaryElementsHeader> _element_headers;

    // FEAT specific data

    /// Mesh type string
    String _mesh_type_string = "";
    /// Shape type
    ParsedShapeType _shape_type = ParsedShapeType::unknown;
    /// Shape dimension
    int _shape_dim = 0;

  public:
    /// Input stream constructor
    explicit GmshFileReader(const std::istream& stream)
    {
      // NOTE: We are reading the whole file into memory here.
      // Because .msh files may contain meshes with mixed elements
      // we have to scan through all elements anyway to determine if
      // FEAT can even represent the mesh.
      // Binary .msh files also contains a mix of ASCII and binary data
      // which makes reading the file line by line annoying.
      read_data(stream);
    }

    /// File path constructor
    explicit GmshFileReader(const std::filesystem::path& path) : GmshFileReader(std::ifstream(path, std::ios::binary))
    {
      XASSERTM(std::filesystem::exists(path), "GmshFileReader: File " + stringify(path) + " does not exist");
      std::ifstream file(path, std::ios::binary);
      XASSERTM(file.good(), "GmshFileReader: Bad stream");
      read_data(file);
    }

    /**
     * \brief Returns the mesh-type string of the parsed mesh
     *
     * \note
     * This function can only be called after calling the #read_root_markup()
     * function as it will always return an empty string otherwise.
     *
     * \returns
     * The mesh-type string of the parsed mesh
     */
    const String& get_meshtype_string() const
    {
      return _mesh_type_string;
    }

    /**
     * \brief Returns the mesh-type of the parsed mesh
     *
     * \note
     * The .msh format has no ability to describe structured meshes,
     * so the mesh type is always MeshType::conformal.
     *
     * \returns
     * The mesh-type of the parsed mesh
     */
    ParsedMeshType get_mesh_type() const
    {
      return ParsedMeshType::conformal;
    }

    /**
     * \brief Returns the shape-type of the parsed mesh
     *
     * \note
     * This function can only be called after calling the #read_root_markup()
     * function as it will always return ShapeType::unknown otherwise.
     *
     * \returns
     * The shape-type of the parsed mesh
     */
    ParsedShapeType get_shape_type() const
    {
      return _shape_type;
    }

    /**
     * \brief Returns the shape-dimension of the parsed mesh
     *
     * \note
     * This function can only be called after calling the #read_root_markup()
     * function as it will always return 0 otherwise.
     *
     * \returns
     * The shape-dimension of the parsed mesh
     */
    int get_shape_dim() const
    {
      return _shape_dim;
    }

    /**
     * \brief Returns the world-dimension of the parsed mesh
     *
     * \note
     * Nodes in .msh files always contain x, y, and z coordinates,
     * thus the world dimension is always three.
     *
     * \returns
     * The world-dimension from the root markup node.
     */
    int get_world_dim() const
    {
      return 3;
    }

    /**
     * \brief Reads the root markup from the input stream and analyses it.
     */
    void read_root_markup()
    {
      // NOTE(mmuegge): Unlike FEATs .xml format, .msh allows meshes with mixed
      // elements. In ASCII meshes each element is tagged with its type, in
      // binary meshes elements are split into chunks of elements of the same
      // type. We thus need to scan through the whole file to ensure it contains
      // a single element type that is representable by FEAT.

      FileView view(_file_contents);

      XASSERT(view.next_line() == "$MeshFormat");

      // Next line is:
      // version-number file-type data-size
      view.parse_line(_version, _is_binary, _data_size);
      XASSERTM(_version == 2.2, "Unsupported .msh file format version!");
      XASSERTM(_data_size == 8, "Unsupported .msh data size!");

      if(_is_binary)
      {
        int one(0);
        view.memcpy(&one);

        XASSERTM(one == 1, ".msh files with non-native endianess are not supported!");
      }

      XASSERT(view.next_line() == "$EndMeshFormat");

      while(!view.done())
      {
        std::string_view line = view.next_line();

        if(line == "$Nodes")
        {
          scan_nodes(view);
        }
        else if(line == "$Elements")
        {
          scan_elements(view);
        }
      }

      XASSERTM(_element_type != -1, ".msh file with no elements");
      _shape_type = element_type_to_shape_type(_element_type);
      _shape_dim = element_type_to_shape_dim(_element_type);
      _mesh_type_string = "conformal:" + stringify(_shape_type) + ":" + stringify(_shape_dim) + ":3";
      _have_root_markup = true;
    }

    /**
     * \brief Parses the mesh file into a mesh node and a mesh atlas.
     *
     * \param[in,out] mesh_atlas
     * Currently unused, only accepted for compatibility
     *
     * \param[in,out] part_set
     * Currently unused, only accepted for compatibility
     *
     * \returns
     * A unique pointer to the root mesh node containing the mesh.
     * This mesh node is automatically linked to the given mesh atlas.
     */
    template<typename MeshType_>
    std::unique_ptr<RootMeshNode<MeshType_>> parse(MeshAtlas<MeshType_>& mesh_atlas, PartitionSet* part_set = nullptr)
    {
      using MeshShapeType = typename MeshType_::ShapeType;
      using DataType = typename MeshType_::CoordType;
      using VertexType = typename MeshType_::VertexType;

      static constexpr std::size_t num_verts = Shape::FaceTraits<MeshShapeType, 0>::count;
      static constexpr int shape_dim = MeshType_::shape_dim;

      if(!_have_root_markup)
      {
        read_root_markup();
      }

      XASSERT(num_verts == nodes_per_element(_element_type));

      // Validate that the chosen mesh type fits to the data
      if(sizeof(DataType) != _data_size)
      {
        XABORTM("CoordType of MeshType_ does not match data type of .msh file");
      }

      XASSERT(MeshShapeType::dimension == _shape_dim);

      if(_shape_type == ParsedShapeType::hypercube && _shape_dim == 2)
      {
        XASSERT((std::is_same_v<MeshShapeType, Shape::Hypercube<2>>));
      }
      else if(_shape_type == ParsedShapeType::hypercube && _shape_dim == 3)
      {
        XASSERT((std::is_same_v<MeshShapeType, Shape::Hypercube<3>>));
      }
      else if(_shape_type == ParsedShapeType::simplex && _shape_dim == 2)
      {
        XASSERT((std::is_same_v<MeshShapeType, Shape::Simplex<2>>));
      }
      else if(_shape_type == ParsedShapeType::simplex && _shape_dim == 3)
      {
        XASSERT((std::is_same_v<MeshShapeType, Shape::Simplex<3>>));
      }

      // part_set is only for compatibility
      (void)part_set;

      // Create mesh
      std::array<Index, shape_dim + 1> size;
      size.fill(0);
      size[0] = _num_nodes;
      size[shape_dim] = _num_elements;

      auto mesh = std::make_unique<MeshType_>(size.data());

      // Parse nodes
      auto& vertex_set = mesh->get_vertex_set();
      FileView node_view(_file_contents, _nodes_start);
      if(_is_binary)
      {
        // Safety: If we reached this point, then scan_nodes succeeded,
        // which means there are the expected number of bytes between $Nodes and $EndNodes.
        // The following memory accesses are thus safe.
        for(std::size_t node(0); node < _num_nodes; node++)
        {
          // Read node index
          int vertex_index(0);
          node_view.memcpy(&vertex_index);

          // This might occur in meshes with non-contiguous indices for vertices.
          // We currently do not support those.
          XASSERT(std::size_t(vertex_index - 1) < _num_nodes);

          // Copy coords
          for(int j(0); j < get_world_dim(); j++)
          {
            // NOTE: .msh files are 1-indexed, hence vertex_index - 1 to transform to 0-based indices
            node_view.memcpy(&vertex_set[Index(vertex_index - 1)][j]);
          }
        }
      }
      else
      {
        // Safety: Cast to int is safe because _num_nodes is stored as an int as per .msh docs.
        for(int i(0); i < (int)_num_nodes; i++)
        {
          Index vertex_index(~Index(0));
          double x(0.0);
          double y(0.0);
          double z(0.0);

          node_view.parse_line(vertex_index, x, y, z);

          // This might occur in meshes with non-contiguous indices for vertices.
          // We currently do not support those.
          XASSERT(std::size_t(vertex_index - 1) < _num_nodes);

          // NOTE: .msh files are 1-indexed, hence vertex_index - 1 to transform to 0-based indices
          vertex_set[vertex_index - 1] = VertexType{x, y, z};
        }
      }

      // Parse elements
      auto& element_index_set = mesh->template get_index_set<shape_dim, 0>();
      if(_is_binary)
      {
        for(BinaryElementsHeader& header : _element_headers)
        {
          FileView chunk_view(_file_contents, header.elements_start);

          // Safety: If we reached this point, then scan_elements succeeded,
          // which means there are the expected number of bytes between $Elements and $EndElements.
          // The following memory accesses are thus safe.
          for(std::size_t e(0); e < std::size_t(header.num_elements); e++)
          {
            // Determine element index
            int element_index(0);
            chunk_view.memcpy(&element_index);

            // This might occur in meshes with non-contiguous indices for elements.
            // We currently do not support those.
            XASSERT(std::size_t(element_index - 1) < _num_elements);

            // Skip past tags
            chunk_view.advance(std::size_t(header.tags_per_element) * sizeof(int));

            for(int n(0); n < (int)num_verts; n++)
            {
              int node_index(0);
              chunk_view.memcpy(&node_index);

              // NOTE: .msh indexing is one-based
              element_index_set(Index(element_index - 1), msh_to_feat_order(header.type, n)) = Index(node_index - 1);
            }
          }
        }
      }
      else
      {
        FileView view(_file_contents, _elements_start);

        for(int i(0); i < (int)_num_elements; i++)
        {
          std::string_view line = view.next_line();
          std::istringstream stream((std::string(line)));

          int element_index(0);
          int type(0);
          int num_tags(0);

          stream >> element_index >> type >> num_tags;

          // This might occur in meshes with non-contiguous indices for elements.
          // We currently do not support those.
          XASSERT(std::size_t(element_index - 1) < _num_elements);

          // Skip over tags
          for(int tag(0); tag < num_tags; tag++)
          {
            int dummy(0);
            stream >> dummy;
          }

          for(int n(0); n < (int)num_verts; n++)
          {
            int node_number(0);
            stream >> node_number;

            element_index_set(Index(element_index - 1), msh_to_feat_order(type, n)) = Index(node_number - 1);
          }

          XASSERT(!stream.fail());
        }
      }

      // Validate element orientation

      mesh->deduct_topology_from_top();
      return std::make_unique<RootMeshNode<MeshType_>>(std::move(mesh), &mesh_atlas);
    }

  private:
    /// Size (in bytes) of one chunk of elements in a binary .msh file
    static std::size_t binary_element_chunk_size(const BinaryElementsHeader& header)
    {
      // 4 bytes for element index, 4 bytes per tag, 4 bytes per vertex index
      return std::size_t(header.num_elements) * (sizeof(int) + (sizeof(int) * std::size_t(header.tags_per_element)) +
                                                 (sizeof(int) * std::size_t(nodes_per_element(header.type))));
    }

    /// .msh element type to ShapeType mapping
    static ParsedShapeType element_type_to_shape_type(int type)
    {
      switch(type)
      {
      case 2:  return ParsedShapeType::simplex;   // Triangle
      case 3:  return ParsedShapeType::hypercube; // Quadrangle
      case 4:  return ParsedShapeType::simplex;   // 4-point tetrahedron
      case 5:  return ParsedShapeType::hypercube; // 8-point hexahedron
      default: XABORTM("Unsupported .msh element type!");
      }
    }

    /// .msh element type to shape dimension mapping
    static int element_type_to_shape_dim(int type)
    {
      switch(type)
      {
      case 2:  return 2; // Triangle
      case 3:  return 2; // Quadrangle
      case 4:  return 3; // 4-point tetrahedron
      case 5:  return 3; // 8-point hexahedron
      default: XABORTM("Unsupported .msh element type!");
      }
    }

    /// .msh element type to number of nodes per element
    static int nodes_per_element(int type)
    {
      switch(type)
      {
      case 2:  return 3; // Triangle
      case 3:  return 4; // Quadrangle
      case 4:  return 4; // 4-point tetrahedron
      case 5:  return 8; // 8-point hexahedron
      default: XABORTM("Unsupported .msh element type!");
      }
    }

    /**
     * \brief Re-ordering for cell vertices
     *
     * .msh file and FEAT use different conventions for the ordering of local
     * vertices of a cell. This function provides a mapping to correct for these
     * different conventions
     */
    static int msh_to_feat_order(int type, int node)
    {
      // Re-ordering for quadrangles
      // .msh are ordered counter-clockwise, FEAT in a zig-zag pattern
      static constexpr std::array<int, 4> quadrangle_ordering{0, 1, 3, 2};

      // Re-ordering for hexahedrons
      static constexpr std::array<int, 8> hexahedron_ordering{0, 1, 3, 2, 4, 5, 7, 6};

      // Re-ordering
      switch(type)
      {
      case 2:  return node; // triangles are ordered the same
      case 3:  return quadrangle_ordering.at(std::size_t(node));
      case 4:  return node; // tetrahedrons are ordered the same
      case 5:  return hexahedron_ordering.at(std::size_t(node));
      default: XABORTM("Unsupported .msh element type!");
      }
    }

    /// Gather number of nodes and start position of node data
    void scan_nodes(FileView& view)
    {
      // Next line contains number of nodes
      view.parse_line(_num_nodes);

      // Skip past nodes
      _nodes_start = view.position();
      if(_is_binary)
      {
        // Each node consists of an integer index and 3 coordinates
        view.advance(_num_nodes * (sizeof(int) + (3 * _data_size)));
      }
      else
      {
        for(Index i(0); i < _num_nodes; i++)
        {
          view.next_line();
        }
      }

      XASSERT(view.next_line() == "$EndNodes");
    }

    /// Gather start position(s) of element data and ensure the mesh contains only a single element type
    void scan_elements(FileView& view)
    {
      view.parse_line(_num_elements);

      _elements_start = view.position();
      if(_is_binary)
      {
        // Each chunk has a header containing:
        // - the element type,
        // - the number of elements of the chunk,
        // - the number of tags per element of the chunk
        // We handle each chunk one by one until we have seen all elements
        std::size_t elements_handled(0);
        while(elements_handled < _num_elements)
        {
          BinaryElementsHeader header{};

          // Read header values
          view.memcpy(&header.type);
          view.memcpy(&header.num_elements);
          view.memcpy(&header.tags_per_element);

          // Record current positions as start of elements for this header
          header.elements_start = view.position();

          // Skip past actual element entries
          view.advance(binary_element_chunk_size(header));

          if(_element_type == -1)
          {
            _element_type = header.type;
          }
          else
          {
            XASSERT(_element_type == header.type);
          }

          elements_handled += std::size_t(header.num_elements);
          _element_headers.push_back(header);
        }
      }
      else
      {
        for(Index i(0); i < _num_elements; i++)
        {
          int number(0);
          int type(0);

          view.parse_line(number, type);

          if(_element_type == -1)
          {
            _element_type = type;
          }
          else
          {
            XASSERT(_element_type == type);
          }
        }
      }

      XASSERT(view.next_line() == "$EndElements");
    }

    /// Read file into string
    void read_data(const std::istream& stream)
    {
      XASSERT(!stream.fail());

      // TODO: Is this an unnecessary copy?
      std::stringstream buffer;
      buffer << stream.rdbuf();

      _file_contents = buffer.str();
    }
  }; // class MshFileReader
} // namespace FEAT::Geometry
