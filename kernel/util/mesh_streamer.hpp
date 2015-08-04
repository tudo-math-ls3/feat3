#pragma once
#ifndef KERNEL_UTIL_MESH_STREAMER_HPP
#define KERNEL_UTIL_MESH_STREAMER_HPP 1

// includes, FEAST
#include <kernel/util/assertion.hpp>
#include <kernel/util/file_error.hpp>

// includes, system
#include <deque>
#include <iostream>
#include <vector>
#include <map>
#include <fstream>

namespace FEAST
{
  /**
   * \brief MeshStreamer class
   *
   * This class enables the storage and analysis of data given as a FEAST mesh/adjacency/coord file.
   * (See the external documentation for further information.)
   *
   * \author Constantin Christof
   * \author Stefan Wahlers
   */
  class MeshStreamer
  {
  public:
    /**
     * \brief Class for saving mesh attributes belonging to one shape dimension
     *
     * \author Jordi Paul
     */
    class AttributeContainer
    {
      public:
        /// Only the default floating point type is supported for values
        typedef Real ValueType;
        /// A value can be vector-valued and this is the type for that
        typedef std::vector<ValueType> ValueVec;

        /// Unique identifier String
        String identifier;
        /// Dimension of the values
        int value_dim;
        /// Number of values
        Index value_count;
        /// Container holding all vector-valued values
        std::vector<ValueVec> values;

        /// \brief Default constructor
        explicit AttributeContainer() :
          identifier(""),
          value_dim(0),
          value_count(0)
        {
          CONTEXT("AttributeContainer::AttributeContainer()");
        }

        /**
         * \brief Useful constructor
         *
         * \param[in] identifier_
         * The identifier for this object.
         *
         * \param[in] value_dim_
         * Dimension of the values.
         *
         * \param[in] value_count_
         * The number of values.
         *
         */
        explicit AttributeContainer(String identifier_, int value_dim_, Index value_count_) :
          identifier(identifier_),
          value_dim(value_dim_),
          value_count(value_count_)
        {
          CONTEXT("AttributeContainer::AttributeContainer(String, Index, Index)");
        }

        /// Default DTOR
        ~AttributeContainer()
        {
          CONTEXT("AttributeContainer::~AttributeContainer()");
        }
    }; // class AttributeContainer

    /**
     * \brief Baseclass for MeshDataContainer containing all data for MeshParts that is parsed from files
     */
    class BaseContainer
    {
    public:
      /// Basic information
      String name;
      /// Name of the parent
      String parent;
      /// Info string
      String info;

      /// Number of entities
      Index num_entities[4];
      /// Slices
      std::vector<Index> slices;
      /// Number of attributes
      Index attribute_count;
      /// Number of coordinates
      Index coord_per_vertex;
      /// parent_indices[d](i) = j means that entity i of dimension d in this container corresponds to entity j in
      // the parent
      std::vector<Index> parent_indices[4];
      /// One set of AttributeContainers for every dimension including dim = 0
      std::vector<AttributeContainer> attributes[4];

    public:
      /// Default CTOR
      BaseContainer()
      {
        CONTEXT("BaseContainer::BaseContainer()");
        for(int i(0); i < 4; ++i)
          num_entities[i] = 0;
      }

      /// Default DTOR
      virtual ~BaseContainer()
      {
        CONTEXT("BaseContainer::~BaseContainer()");
      }

      /**
       * \brief Parses a mesh-info-data-input stream.
       *
       * \param[in] ifs
       * A reference to the input stream to be parsed.
       *
       * \param[in] cur_line
       * A counter that specifies the number of lines read so far.
       *
       * \returns
       * The string containing the information.
       */
      String _parse_info_section(Index& cur_line, std::istream& ifs);

      /**
       * \brief Parses a counts subchunk
       *
       * \param[in] ifs
       * A reference to the input stream to be parsed.
       *
       * \param[in] cur_line
       * A counter that specifies the number of lines read so far.
       *
       * \returns
       * The number of lines that has been read when the routine is done.
       */
      Index _parse_counts_chunk(Index cur_line, std::istream& ifs);

      /**
       * \brief Parses a parent-indices subchunk
       *
       * \param[in] cur_line
       * A counter that specifies the number of lines read so far.
       *
       * \param[in] ifs
       * A reference to the input stream to be parsed.
       *
       * \param[in] line
       * The last line that was read (and therefore the first line of the subchunk)
       *
       * \returns
       * The number of lines that has been read when the routine is done.
       */
      Index _parse_parents_chunk(Index cur_line, std::istream& ifs, String line);
    }; // class MeshStreamer::BaseContainer

    /**
     * \brief Mesh data container class
     *
     * This class stores the data related to a (sub-) mesh.
     *
     * \author Constantin Christof
     * \author Stefan Wahlers
     */
    class MeshDataContainer :
      public BaseContainer
    {
    public:
      // coordinate- and adjacency vectors
      typedef std::vector<double> CoordVec;
      typedef std::vector<CoordVec> CoordStack;
      typedef std::vector<Index> AdjVec;
      typedef std::vector<AdjVec> AdjStack;

      /// mesh-type enumeration
      enum MeshType
      {
        /// unknown mesh type
        mt_unknown = 0,
        /// conformal mesh type
        mt_conformal,
        /// structured mesh type
        mt_structured
      };

      /// shape-type enumerations
      enum ShapeType
      {
        /// unknown shape type
        st_unknown = 0,
        /// 0d vertex set
        st_vert,
        /// 1D edge mesh
        st_edge,
        /// 2D triangular mesh
        st_tria,
        /// 2D quadrilateral mesh
        st_quad,
        /// 3D tetrahedral mesh
        st_tetra,
        /// 3D hexahedral mesh
        st_hexa
      };

    public:
      /// Name of the chart that contains an analytic description of the geometry
      String chart;
      /// Mesh data
      MeshType mesh_type;
      /// Shape of the mesh
      ShapeType shape_type;
      /// Coordinates
      CoordStack coords;
      /// Adjacencies: Adjacency[i][j] contains which entities of dimension i lie at entities of dimension j
      AdjStack adjacencies[4][4];
      /// Should parent information be deducted?
      String deduct_target_sets;
      /// Should the topology be deducted from the parent information?
      bool deduct_topology;

    public:
      // default CTOR
      MeshDataContainer()
      {
        CONTEXT("MeshDataContainer::MeshDataContainer()");
      }

      // default DTOR
      virtual ~MeshDataContainer()
      {
        CONTEXT("MeshDataContainer::~MeshDataContainer()");
      }

      /**
      * \brief Converts a MeshType to a String.
      *
      * \param[in] mesh_type
      * A mesh_type object to be converted to a String
      *
      * \returns
      * The mesh_type as a String
      */
      static String convert_mesh_type(const MeshType mesh_type);

       /**
      * \brief Converts a String to a MeshType.
      *
      * \param[in] mesh_type
      * A String to be converted to a MeshType
      *
      * \returns
      * The mesh_type as a MeshType
      */
      static MeshType convert_mesh_type(const String mesh_type);

      /**
      * \brief Converts a ShapeType to a String.
      *
      * \param[in] shape_type
      * A shape_type object to be converted to a String
      *
      * \returns
      * The shape_type as a String
      */
      static String convert_shape_type(const ShapeType shape_type);

      /**
       * \brief Converts a String to a ShapeType
       *
       * \param[in] shape_type
       * A String to be converted to a ShapeType
       *
       * \returns
       * The shape_type as a ShapeType
       */
      static ShapeType convert_shape_type(const String shape_type);

      /**
       * \brief Parses an attribute subchunk
       *
       * \param[in] ifs
       * A reference to the input stream to be parsed.
       *
       * \param[in] cur_line
       * A counter that specifies the number of lines read so far.
       *
       * \returns
       * The number of lines that has been read when the routine is done.
       */
      Index _parse_attribute_chunk(Index cur_line, std::istream& ifs);

      /**
      * \brief Parses a mesh-section-data-input stream.
      *
      * \param[in] cur_line
      * A counter that specifies the number of lines read so far.
      *
      * \param[in] is_meshpart
      * Specifies whether the data belongs to a sub-mesh (true) or the root mesh (false).
      *
      * \param[in] ifs
      * A reference to the input stream to be parsed.
      *
      * \returns
      * The number of lines that has been read when the routine is done.
      */
      Index _parse_mesh_section(Index cur_line, bool is_meshpart, std::istream& ifs);

      /**
       * \brief Parses the header of a mesh section
       *
       * \param[in] cur_line
       * A counter that specifies the number of lines read so far.
       *
       * \param[in] is_meshpart
       * Specifies whether the data belongs to a sub-mesh (true) or the root mesh (false).
       *
       * \param[in] ifs
       * A reference to the input stream to be parsed.
       *
       * \returns
       * The number of lines that has been read when the routine is done.
       */
      Index _parse_mesh_header_section(Index cur_line, bool is_meshpart, std::istream& ifs);

      /**
       * \brief Parses a given FEAST- coord file.
       *
       * This function parses the FEAST- coord file given by "filepath" and saves the data in
       * the MeshDataContainer that is specified by the pointer "mesh".
       *
       * \param[in] filename
       * A String containing the path to the coord file.
       *
       */
      void parse_coord_file(String filename);

      /**
       * \brief Parses a coord-file-data-input stream.
       *
       * \param[in] ifs
       * A reference to the input stream to be parsed.
       *
       */
      void parse_coord_file(std::istream& ifs);

      /**
       * \brief Parses a given FEAST- adjacency file.
       *
       * This function parses the FEAST- adjacency file given by "filepath" and saves the data in
       * the MeshDataContainer that is specified by the pointer "mesh".
       *
       * \param[in] filename
       *
       * A String containing the path to the adjacency file.
       *
       */
      void parse_adjacency_file(String filename);

      /**
       * \brief Parses an adjacency-file-data-input stream.
       *
       * \param[in] ifs
       * A reference to the input stream to be parsed.
       *
       */
      void parse_adjacency_file(std::istream& ifs);

      /**
       * \brief Parses a coords subchunk
       *
       * \param[in] cur_line
       * A counter that specifies the number of lines read so far.
       *
       * \param[in] ifs
       * A reference to the input stream to be parsed.
       *
       * \returns
       * The number of lines that has been read when the routine is done.
       */
      Index _parse_coords_chunk(Index cur_line, std::istream& ifs);

      /**
       * \brief Parses an adjacency subchunk
       *
       * \param[in] ifs
       * A reference to the input stream to be parsed.
       *
       * \param[in] cur_line
       * A counter that specifies the number of lines read so far.
       *
       * \param[in] line
       * The last line that was read (and therefore the first line of the subchunk)
       *
       * \returns
       * The number of lines that has been read when the routine is done.
       */
      Index _parse_adjacency_chunk(Index cur_line, std::istream& ifs, String line);

    }; // MeshDataContainer

    /**
     * \brief Container for saving chart information parsed from files
     *
     * As opposed to the other containers, only meta information gets parsed into this container (i.e. the name).
     * Because charts need very different parameters (i.e. an analytic description of a circle only needs midpoint
     * and radius, but a polygon line needs several points), everything in the data chunk gets parsed into a
     * std::deque of Strings and passed to the factory that is used to construct the Chart objects.
     *
     * \author Jordi Paul
     *
     */
    class ChartContainer
    {
      public:
        /// Description
        String info;
        /// Name the chart gets referenced by, i.e. by MeshParts
        String name;
        /// Type as a String, i.e. circle
        String type;
        /// Container holding the data the streamer reads from the file and then passes to a factory
        std::deque<String> data;
        /// Line in the file the data section starts at
        Index start_of_data;

        MeshDataContainer mesh_data;

        /**
         * \brief Default constructor
         */
        explicit ChartContainer() :
          info(""),
          name(""),
          type(""),
          data(),
          start_of_data(0),
          mesh_data()
        {
        }

        /**
         * \brief Copy constructor
         *
         * \param[in] other
         * The other ChartContainer to be copied.
         */
        ChartContainer(const ChartContainer& other) :
          info(other.info),
          name(other.name),
          type(other.type),
          data(other.data),
          start_of_data(other.start_of_data),
          mesh_data(other.mesh_data)
          {
          }

        /// Destructor with nothing to do
        ~ChartContainer()
        {
        }

        /**
         * \brief Parses a chart section from a stream
         *
         * \param[in] cur_line
         * The starting line number in the input stream.
         *
         * \param[in] ifs
         * The input stream to parse from.
         *
         * \returns
         * The number of lines that has been read when the routine is done.
         *
         */
        Index _parse_chart_section(Index cur_line, std::istream& ifs);

        /**
         * \brief Parses a chart's header section from a stream
         *
         * \param[in] cur_line
         * The starting line number in the input stream.
         *
         * \param[in] ifs
         * The input stream to parse from.
         *
         * \returns
         * The number of lines that has been read when the routine is done.
         *
         */
        Index _parse_chart_header_section(Index cur_line, std::istream& ifs);

    }; // class ChartContainer

    /**
     * \brief Mesh node class
     *
     * This class implements a tree node that contains a MeshDataContainer as well as its
     * child cell-set and sub-mesh nodes.
     *
     * \author Peter Zajac
     */
    class MeshNode
    {
    public:
      /// a map of meshpart nodes
      typedef std::map<String, MeshNode*, String::NoCaseLess> MeshpartMap;

      /// the mesh data container of this node
      MeshDataContainer mesh_data;
      /// the sub-mesh map of this node
      MeshpartMap meshpart_map;

    public:
      /// Virtual DTOR
      virtual ~MeshNode()
      {
        MeshpartMap::iterator it(meshpart_map.begin()), jt(meshpart_map.end());
        for(; it != jt; ++it)
        {
          delete it->second;
        }
      }

      /**
       * \brief Finds a sub-mesh within this sub-tree.
       */
      MeshNode* find_meshpart(String name)
      {
        // perform depth-first-search
        MeshpartMap::iterator it(meshpart_map.begin()), jt(meshpart_map.end());
        for(; it != jt; ++it)
        {
          if(it->first.compare_no_case(name) == 0)
            return it->second;

          MeshNode* node = it->second->find_meshpart(name);
          if(node != nullptr)
            return node;
        }

        // parent not found
        return nullptr;
      }

      /**
       * \brief returns the total number of meshparts in the subtree
       * which are related to this MeshNode
       */
      Index get_num_meshparts_below()
      {
        Index count(Index(meshpart_map.size()));
        // perform depth-first-search
        MeshpartMap::iterator it(meshpart_map.begin()), jt(meshpart_map.end());
        for(; it != jt; ++it)
        {
          count += it->second->get_num_meshparts_below();
        }
        return count;
      }

      /**
       * \brief Writes the mesh data of this mesh and all is_meshpartes related to this mesh
       *  into the output stream.
       */
      void write(std::ostream& ofs, bool is_meshpart) const;
    };

  private:
    /// Number of charts
    Index _num_charts;
    /// Number of mesh parts
    Index _num_meshparts;
    /// Description
    String _info;
    /// root mesh node
    MeshNode* _root_mesh_node;
    /// For error checking: Number of charts parsed
    Index _num_parsed_charts;
    /// For error checking: Number of meshparts parsed
    Index _num_parsed_meshparts;

  public:
    /// All Charts this Mesh and all of its MeshParts refer to
    std::deque<ChartContainer> charts;

  public:
    /// Default Constructor
    MeshStreamer();

    /**
     * \brief Read-From-File Constructor
     *
     * This constructor automatically parses the mesh file.
     *
     * \param[in] filename
     * The name of the mesh file to be parsed.
     */
    explicit MeshStreamer(const String& filename);

    /// Virtual Destructor
    virtual ~MeshStreamer();

    /**
     * \brief Parses multiple files given in one String
     *
     * \param[in] filenames
     * All filenames seperated by whitespace.
     */
    void parse_multiple_files(std::deque<String>& filenames);

    /**
     * \brief Parses a given FEAST- mesh file.
     *
     * This function parses the FEAST- mesh file given by "filepath" and saves the data.
     *
     * \param[in] filename
     * A String containing the path to the mesh file.
     */
    void parse_mesh_file(String filename);

    /**
     * \brief Parses a mesh-data-input stream.
     *
     * \param[in] ifs
     * A reference to the input stream to be parsed.
     *
     */
    void parse_mesh_file(std::istream& ifs);

    /**
     * \brief Inserts mesh_node into the tree structure which depends on root
     *
     * \param[in] mesh_node
     * The MeshNode to be inserted.
     */
    void _insert_meshpart(MeshStreamer::MeshNode* mesh_node);

    /**
     * \brief deletes mesh_node from the tree structure
     *
     * \param[in] mesh_node
     * The MeshNode to be deleted.
     */
    void _delete_meshpart(MeshNode* mesh_node);

    /**
     * \brief deletes mesh_node by name from the tree structure
     *
     * \param[in] name
     * The name of the MeshNode to be deleted.
     */
    void _delete_meshpart(String name);

    /**
     * \brief Returns the number of is_meshpartes.
     */
    Index get_num_meshparts() const;

    /**
     * \brief Returns the file's information.
     */
    String get_info() const;

    /**
     * \brief Returns a pointer to the root-mesh node.
     */
    MeshStreamer::MeshNode* get_root_mesh_node()
    {
      return _root_mesh_node;
    }

    /**
     * \brief Returns a mesh data container.
     *
     * \param[in] name
     * The name of the (sub-)mesh that is to be returned.
     */
    MeshStreamer::MeshDataContainer* get_mesh(String name = "root");

    /**
     * \brief Writes the mesh data into a file
     *
     * \param[in] filename
     * The name of the file the data has to be written in
     */
    void write_mesh_file(String filename) const;

    /**
     * \brief Writes the mesh data into a stream
     *
     * \param[in] ofs
     * The output stream that is to be written to.
     */
    void write_mesh_file(std::ostream& ofs) const;


  private:
    /// \brief Finds a MeshNode by name String
    MeshStreamer::MeshNode* _find_meshpart_parent(String parent_name);

    /**
     * \brief Parses a mesh-header-data-input stream.
     *
     * \param[in] ifs
     * A reference to the input stream to be parsed.
     *
     * \param[in] cur_line
     * A counter that specifies the number of lines read so far.
     *
     * \returns
     * The number of lines that has been read when the routine is done.
     */
     Index _parse_header_section(Index cur_line, std::istream& ifs);

    /**
     * \brief Parses a mesh-info-data-input stream.
     *
     * \param[in] ifs
     * A reference to the input stream to be parsed.
     *
     * \param[in] cur_line
     * A counter that specifies the number of lines read so far.
     *
     * \returns
     * The string containing the information.
     */
    String _parse_info_section(Index& cur_line, std::istream& ifs);


    /**
     * \brief Parses a mesh-section-data-input stream.
     *
     * \param[in] ifs
     * A reference to the input stream to be parsed.
     *
     * \param[in] cur_line
     * A counter that specifies the number of lines read so far.
     *
     * \param[in] is_meshpart
     * Specifies whether the data belongs to a sub-mesh (true) or the root mesh (false).
     *
     * \returns
     * The number of lines that has been read when the routine is done.
     */
    Index _parse_mesh_section(Index cur_line, bool is_meshpart, std::istream& ifs);

  }; // class MeshStreamer

  /// \cond internal
  namespace Intern
  {
    /**
     * \brief Wrapper class for writing back information to a MeshDataContainer
     *
     * \tparam dim_
     * Dimension for the highest dimensional shape, i.e. 3 for a Hypercube<3> mesh
     *
     * \author Jordi Paul
     *
     */
    template<int dim_>
    struct MeshDataContainerUpdater
    {
      /**
       * \brief Updates data in a MeshDataContainer from an IndexSetHolder
       *
       * If a streamed mesh only has vertex@cell information, the other information like edge@face etc. is
       * generated by the RedundantIndexSetBuilder and can be written back to the MeshDataContainer so it can be
       * used by other methods.
       *
       * \tparam IndexSetHolderType_
       * Type for the IndexSetHolder to copy the new information from
       *
       * \param[in,out] mesh_data
       * MeshDataContainer to be updated
       *
       * \param[in] ish
       * IndexSetHolder containing the new information
       *
       *
       */
      template<typename IndexSetHolderType_>
      static void update_from_ish(MeshStreamer::MeshDataContainer* DOXY(mesh_data),
      const IndexSetHolderType_& DOXY(ish))
      {
      }

    };

    template<>
    struct MeshDataContainerUpdater<2>
    {
      template<typename IndexSetHolderType_>
      static void update_from_ish( MeshStreamer::MeshDataContainer* mesh_data, const IndexSetHolderType_& ish)
      {
        // The first relevant data is vertex@edge, meaning dimension 1
        if(mesh_data->num_entities[1] == 0)
        {
          ASSERT(mesh_data->adjacencies[0][1].size() == 0, "_mesh_data contains no edges, but adjacency vector is not empty!");
          // Update vertex@edge information
          auto& vert_at_edge = ish.template get_index_set<1,0>();
          // Set number of edges
          mesh_data->num_entities[1] = vert_at_edge.get_num_entities();

          for(Index edge(0); edge < mesh_data->num_entities[1]; ++edge)
          {
            // idx will contain all indices for edge
            std::vector<Index> idx;
            for(int i(0); i < vert_at_edge.num_indices; ++i)
              idx.push_back(vert_at_edge(edge,Index(i)));

            // Add this edge to the global adjacency vector
            mesh_data->adjacencies[0][1].push_back(idx);
          }

          // Next: edge@cell
          ASSERT(mesh_data->adjacencies[1][2].size() == 0, "_mesh_data has no edges, but edge@cell adjacency information!");
          // Update edge@cell information
          auto& edge_at_cell = ish.template get_index_set<1,0>();

          // Do not set the number of shapes of dimension 2 here. Either this method is called in a 2d context
          // (then that number is correct anyway) or in a 3d context (then that number is set in the higher
          // dimensional) version

          for(Index cell(0); cell < mesh_data->num_entities[2]; ++cell)
          {
            // idx will contain all indices for edge
            std::vector<Index> idx;
            // The _index_bound for the cell to edges mapping is the number of edges
            for(int edge(0); edge < edge_at_cell.num_indices; ++edge)
              idx.push_back(edge_at_cell(cell,Index(edge)));
            // Add this edge to the global adjacency vector
            mesh_data->adjacencies[1][2].push_back(idx);
          }
        }
      } // MeshDataContainerUpdater<2>::update_from_ish()

    }; //MeshDataContainerUpdater<2>

    template<>
    struct MeshDataContainerUpdater<3>
    {
      template<typename IndexSetHolderType_>
      static void update_from_ish(MeshStreamer::MeshDataContainer* mesh_data, const IndexSetHolderType_& ish)
      {
        // Update 2d data first
        MeshDataContainerUpdater<2>::update_from_ish(mesh_data, ish);

        // Next: Everything face related
        // Check if the data is missing
        if(mesh_data->num_entities[2] == 0)
        {
          // vertex@face
          ASSERT(mesh_data->adjacencies[0][2].size() == 0, "_mesh_data has no faces, but vertex adjacency vector is not empty!");
          auto& vert_at_face = ish.template get_index_set<2,0>();

          // Set number of faces
          mesh_data->num_entities[2] = vert_at_face.get_num_entities();

          for(Index face(0); face < mesh_data->num_entities[2]; ++face)
          {
            // idx will contain all indices for face
            std::vector<Index> idx;
            for(int i(0); i < vert_at_face.num_indices; ++i)
              idx.push_back(vert_at_face(face,Index(i)));

            // Add this edge to the global adjacency vector
            mesh_data->adjacencies[0][2].push_back(idx);
          }

          // edge@face was set by the 2d routine called earlier

          // Next: edge@cell
          ASSERT(mesh_data->adjacencies[1][3].size() == 0, "_mesh_data has no edges, but adjacency vector for faces is not empty!");

          auto& edge_at_cell= ish.template get_index_set<2,1>();

          // Update edge@cell information
          for(Index cell(0); cell < mesh_data->num_entities[3]; ++cell)
          {
            // idx will contain all indices for edge
            std::vector<Index> idx;
            // num_indices is the (local!) number of edges per cell
            for(int edge(0); edge < edge_at_cell.num_indices; ++edge)
              idx.push_back(edge_at_cell(cell,Index(edge)));
            // Add this edge to the global adjacency vector
            mesh_data->adjacencies[1][3].push_back(idx);
          }

          // Last: face@cell
          ASSERT(mesh_data->adjacencies[2][3].size() == 0, "_mesh_data has no faces, but adjacency vector for faces is not empty!");
          auto& face_at_cell= ish.template get_index_set<3,2>();

          // Update face@cell information
          for(Index cell(0); cell < mesh_data->num_entities[3]; ++cell)
          {
            // idx will contain all indices for face
            std::vector<Index> idx;
            // num_indices is the (local!) number of faces per cell
            for(int face(0); face < face_at_cell.num_indices; ++face)
              idx.push_back(face_at_cell(cell,Index(face)));
            // Add this face to the global adjacency vector
            mesh_data->adjacencies[2][3].push_back(idx);
          }

        } // handling faces

      } // MeshDataContainerUpdater<3>::update_from_ish()

    }; // MeshDataContainerUpdater<3>
  } // namespace Intern
  /// \endcond
} // namespace FEAST

#endif // KERNEL_UTIL_MESH_STREAMER_HPP
