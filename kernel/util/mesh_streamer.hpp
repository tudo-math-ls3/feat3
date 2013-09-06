#pragma once
#ifndef KERNEL_UTIL_MESH_STREAMER_HPP
#define KERNEL_UTIL_MESH_STREAMER_HPP 1

// includes, FEAST
#include <kernel/util/file_error.hpp>

// includes, system
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
    class BaseContainer
    {
    public:
      // basic information
      String name;
      String parent;
      String info;

      // number of entities
      Index vertex_count;
      Index edge_count;
      Index tria_count;
      Index quad_count;
      Index tetra_count;
      Index hexa_count;
      std::vector<Index> slices;

      Index coord_per_vertex;

      // parent indices ([0] -> vertices etc.)
      std::vector<Index> parent_indices [4];

    public:
      // default CTOR
      BaseContainer()
      {
        CONTEXT("BaseContainer::BaseContainer()");
      }

      // default DTOR
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
       * The number of lines that has been read when the programme is done.
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
       * The number of lines that has been read when the programme is done.
       */
      Index _parse_parents_chunk(Index cur_line, std::istream& ifs, String line);
    };

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
        /// 1D edge mesh
        st_edge,
        /// 2D triangular mesh
        st_tria,
        /// 2D quadrilateral mesh
        st_quad,
        /// 3D tetrahedral mesh
        st_tetra,
        /// 3D hexahedral mesh
        st_hexa,
        /// 2D mixed triangular/quadrilateral mesh
        st_tria_quad,
        /// 3D mixed tetrahedral/hexahedral mesh
        st_tetra_hexa
      };

    public:

      // basic information
      String chart;

      // mesh data
      MeshType mesh_type;
      ShapeType shape_type;

      // coordinates
      CoordStack coords;

      // adjacencies
      AdjStack adjacencies[4][4];

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
      * \brief Parses a mesh-section-data-input stream.
      *
      * \param[in] cur_line
      * A counter that specifies the number of lines read so far.
      *
      * \param[in] submesh
      * Specifies whether the data belongs to a sub-mesh (true) or the root mesh (false).
      *
      * \param[in] ifs
      * A reference to the input stream to be parsed.
      *
      * \returns
      * The number of lines that has been read when the programme is done.
      */
      Index _parse_mesh_section(Index cur_line, bool submesh, std::istream& ifs);

      /**
       * \brief Parses a given FEAST- coord file.
       *
       * This function parses the FEAST- coord file given by "filepath" and saves the data in
       * the MeshDataContainer that is specified by the pointer "mesh".
       *
       * \param[in] filepath
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
       * \param[in] file
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
       * The number of lines that has been read when the programme is done.
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
       * The number of lines that has been read when the programme is done.
       */
      Index _parse_adjacency_chunk(Index cur_line, std::istream& ifs, String line);

    }; // MeshDataContainer

    /**
     * \brief Cell set container class
     *
     * This class stores the data related to a cell set.
     *
     * \author Constantin Christof
     */
    class CellSetContainer :
      public BaseContainer
    {
    public:
      // default CTOR
      CellSetContainer()
      {
        CONTEXT("CellSetContainer::CellSetContainer()");
      }

      // default DTOR
      ~CellSetContainer()
      {
        CONTEXT("CellSetContainer::~CellSetContainer()");
      }

      /**
       * \brief Parses a mesh-cellset-data-input stream.
       *
       * \param[in] ifs
       * A reference to the input stream to be parsed.
       *
       * \param[in] cur_line
       * A counter that specifies the number of lines read so far.
       *
       * \returns
       * The number of lines that has been read when the programme is done.
       */
      Index _parse_cellset_section(Index cur_line, std::istream& ifs);

    }; // CellSetContainer

    class CellSetNode;

    /**
     * \brief Base-Class for Cell-Set node parents
     *
     * This class acts as a base class for the CellSetNode and MeshNode classes.
     *
     * \author Peter Zajac
     */
    class CellSetParent
    {
    public:
      // a map of cell-set nodes
      typedef std::map<String, CellSetNode*, String::NoCaseLess> CellSetMap;

      // the cell-set map of this node
      CellSetMap cell_set_map;

    public:
      virtual ~CellSetParent()
      {
        CellSetMap::iterator it(cell_set_map.begin()), jt(cell_set_map.end());
        for(; it != jt; ++it)
        {
          delete it->second;
        }
      }

      /**
       * \brief Finds a cell-set node within this sub-tree.
       */
      CellSetNode* find_cell_set(String name)
      {
        // perform depth-first-search
        CellSetMap::iterator it(cell_set_map.begin()), jt(cell_set_map.end());
        for(; it != jt; ++it)
        {
          if(it->first.compare_no_case(name) == 0)
            return it->second;

          CellSetNode* node = it->second->find_cell_set(name);
          if(node != nullptr)
            return node;
        }

        // parent not found
        return nullptr;
      }
    };

    /**
     * \brief Cell-Set node class
     *
     * This class implements a tree node that contains a CellSetContainer as well as
     * its child cell-set nodes.
     *
     * \author Peter Zajac
     */
    class CellSetNode :
      public CellSetParent
    {
    public:
      CellSetContainer cell_set;

      /**
       * \brief Writes the stored cell set data into the output stream.
       */
      void write(std::ostream &ofs) const;

    };

    /**
     * \brief Mesh node class
     *
     * This class implements a tree node that contains a MeshDataContainer as well as its
     * child cell-set and sub-mesh nodes.
     *
     * \author Peter Zajac
     */
    class MeshNode :
      public CellSetParent
    {
    public:
      // a map of sub-mesh nodes
      typedef std::map<String, MeshNode*, String::NoCaseLess> SubMeshMap;

      // the mesh data container of this node
      MeshDataContainer mesh_data;
      // the sub-mesh map of this node
      SubMeshMap sub_mesh_map;

    public:
      virtual ~MeshNode()
      {
        SubMeshMap::iterator it(sub_mesh_map.begin()), jt(sub_mesh_map.end());
        for(; it != jt; ++it)
        {
          delete it->second;
        }
      }

      /**
       * \brief Finds a sub-mesh within this sub-tree.
       */
      MeshNode* find_sub_mesh(String name)
      {
        // perform depth-first-search
        SubMeshMap::iterator it(sub_mesh_map.begin()), jt(sub_mesh_map.end());
        for(; it != jt; ++it)
        {
          if(it->first.compare_no_case(name) == 0)
            return it->second;

          MeshNode* node = it->second->find_sub_mesh(name);
          if(node != nullptr)
            return node;
        }

        // parent not found
        return nullptr;
      }

      /**
       * \brief returns the total number of submeshes of the subtree
       * which are related to this MeshNode
       */
      Index get_num_sub_meshes_below()
      {
        Index count(Index(sub_mesh_map.size()));
        // perform depth-first-search
        SubMeshMap::iterator it(sub_mesh_map.begin()), jt(sub_mesh_map.end());
        for(; it != jt; ++it)
        {
          count += it->second->get_num_sub_meshes_below();
        }
        return count;
      }

      /**
       * \brief Finds a cell-set parent node within this sub-tree.
       *
       * The returned node may be either a CellSetNode or a MeshNode.
       */
      CellSetParent* find_cell_set_parent(String parent_name)
      {
        // search mesh nodes
        MeshNode* node = find_sub_mesh(parent_name);
        if(node != nullptr)
          return node;

        // search cell sets
        return find_cell_set(parent_name);
      }

      /**
       * \brief Writes the mesh data of this mesh and all submeshes related to this mesh
       *  into the output stream.
       */
      void write(std::ostream& ofs, bool submesh) const;
    };

  private:

    // general information
    Index _num_submeshes;
    Index _num_cellsets;

    // file paths
    String _chart_path;
    // info
    String _info;

    // root mesh node
    MeshNode* _root_mesh_node;

  public:
    /// Default Constructor
    MeshStreamer();

    /// Virtual Destructor
    virtual ~MeshStreamer();

    /**
     * \brief Parses a given FEAST- mesh file.
     *
     * This function parses the FEAST- mesh file given by "filepath" and saves the data.
     *
     * \param[in] filepath
     * A String containing the path to the mesh file.
     *
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
     * \brief Returns the chart path.
     */
    String get_chart_path() const;

    /**
     * \brief Inserts mesh_node into the tree structure which depends on root
     *
     * \param[in] mesh_node
     * The MeshNode to be inserted.
     */
    void _insert_sub_mesh(MeshStreamer::MeshNode* mesh_node);

    /**
     * \brief deletes mesh_node from the tree structure
     *
     * \param[in] mesh_node
     * The MeshNode to be deleted.
     */
    void _delete_sub_mesh(MeshNode* mesh_node);

    /**
     * \brief deletes mesh_node by name from the tree structure
     *
     * \param[in] name
     * The name of the MeshNode to be deleted.
     */
    void _delete_sub_mesh(String name);

    /**
     * \brief Returns the number of submeshes.
     */
    Index get_num_submeshes() const;

    /**
     * \brief Returns the number of cellsets.
     */
    Index get_num_cellsets() const;

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
     * \brief Returns a cell-set container.
     *
     * \param[in] name
     * The name of the cell-set that is to be returned.
     */
    MeshStreamer::CellSetContainer* get_cell_set(String name);

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

    MeshStreamer::CellSetParent* _find_cell_set_parent(String parent_name);
    MeshStreamer::MeshNode* _find_sub_mesh_parent(String parent_name);

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
     * The number of lines that has been read when the programme is done.
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
     * \param[in] submesh
     * Specifies whether the data belongs to a sub-mesh (true) or the root mesh (false).
     *
     * \returns
     * The number of lines that has been read when the programme is done.
     */
    Index _parse_mesh_section(Index cur_line, bool submesh, std::istream& ifs);

    /**
     * \brief Parses a mesh-cellset-data-input stream.
     *
     * \param[in] ifs
     * A reference to the input stream to be parsed.
     *
     * \param[in] cur_line
     * A counter that specifies the number of lines read so far.
     *
     * \returns
     * The number of lines that has been read when the programme is done.
     */
    Index _parse_cellset_section(Index cur_line, std::istream& ifs);
  }; // class MeshStreamer
} // namespace FEAST

#endif // KERNEL_UTIL_MESH_STREAMER_HPP
