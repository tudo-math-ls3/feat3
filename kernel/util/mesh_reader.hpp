#pragma once
#ifndef KERNEL_UTIL_MESH_READER_HPP
#define KERNEL_UTIL_MESH_READER_HPP 1

// includes, FEAST
#include <kernel/util/file_error.hpp>

// includes, system
#include <iostream>
#include <vector>
#include <map>

namespace FEAST
{
  /**
   * \brief Meshreader class
   *
   * This class enables the storage and analysis of data given as a FEAST mesh/adjacency/coord file.
   * (See the external documentation for further information.)
   *
   * \author Constantin Christof
   */
  class MeshReader
  {
  public:
    class BaseContainer
    {
    public:
      // basic information
      String name;
      String parent;

      // number of entities
      Index vertex_number;
      Index edge_number;
      Index tria_number;
      Index quad_number;
      Index tetra_number;
      Index hexa_number;

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
    };

    /**
     * \brief Mesh data container class
     *
     * This class stores the data related to a (sub-) mesh.
     *
     * \author Constantin Christof
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

    public:

      // basic information
      String chart;
      String coord_version;
      String adjacency_version;

      // mesh data
      String mesh_type;
      String shape_type;
      Index coord_per_vertex;

      // file paths
      String coord_path;
      String adj_path;

      // coordinates
      CoordStack coords;

      // adjacencies
      AdjStack adjacencies [4][4];

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
    };

  private:
    // general information
    Index _version;
    Index _num_submeshes;
    Index _num_cellsets;

    // file paths
    String _chart_path;

    // root mesh node
    MeshNode* _root_mesh_node;

  public:
    /// Default Constructor
    MeshReader();

    /// Virtual Destructor
    virtual ~MeshReader();

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
     * \brief Parses a given FEAST- coord file.
     *
     * This function parses the FEAST- coord file given by "filepath" and saves the data within
     * the MeshDataContainer that is specified by the pointer "mesh".
     *
     * \param[in] filepath
     * A String containing the path to the coord file.
     *
     * \param[in] mesh
     * A pointer to the MeshDataContainer the data shall be saved in.
     *
     */
    void parse_coord_file(String filename, MeshDataContainer *mesh);

    /**
     * \brief Parses a coord-file-data-input stream.
     *
     * \param[in] ifs
     * A reference to the input stream to be parsed.
     *
     * \param[in] mesh
     * A pointer to the MeshDataContainer the data shall be saved in.
     *
     */
    void parse_coord_file(std::istream& ifs, MeshDataContainer *mesh);

    /**
     * \brief Parses a given FEAST- adjacency file.
     *
     * This function parses the FEAST- adjacency file given by "filepath" and saves the data within
     * the MeshDataContainer that is specified by the pointer "mesh".
     *
     * \param[in] filepath
     * A String containing the path to the adjacency file.
     *
     * \param[in] mesh
     * A pointer to the MeshDataContainer the data shall be saved in.
     *
     */
    void parse_adjacency_file(String filename, MeshDataContainer *mesh);

    /**
     * \brief Parses an adjacency-file-data-input stream.
     *
     * \param[in] ifs
     * A reference to the input stream to be parsed.
     *
     * \param[in] mesh
     * A pointer to the MeshDataContainer the data shall be saved in.
     *
     */
    void parse_adjacency_file(std::istream& ifs, MeshDataContainer *mesh);

    /**
     * \brief Returns the version.
     */
    Index get_version() const;

    /**
     * \brief Returns the chart path.
     */
    String get_chart_path() const;

    /**
     * \brief Returns the number of submeshes.
     */
    Index get_num_submeshes() const;

    /**
     * \brief Returns the number of cellsets.
     */
    Index get_num_cellsets() const;

    /**
     * \brief Returns a pointer to the root-mesh node.
     */
    MeshReader::MeshNode* get_root_mesh_node()
    {
      return _root_mesh_node;
    }

    /**
     * \brief Returns a mesh data container.
     *
     * \param[in] name
     * The name of the (sub-)mesh that is to be returned.
     */
    MeshReader::MeshDataContainer* get_mesh(String name = "root");

    /**
     * \brief Returns a cell-set container.
     *
     * \param[in] name
     * The name of the cell-set that is to be returned.
     */
    MeshReader::CellSetContainer* get_cell_set(String name);

    /*
     * \brief Returns the mesh data container that belongs to the given name.
     * The returned bool is false if no such mesh was found.
     *
     * \param[in] name
     * The name of the mesh that is needed.
     */
    //std::pair<MeshReader::MeshDataContainer, bool> get_mesh(String mesh_name);

    /*
     * \brief Returns the cellset container that belongs to the given name.
     * The returned bool is false if no such set was found.
     *
     * \param[in] name
     * The name of the cellset that is needed.
     */
    //std::pair<MeshReader::CellSetContainer, bool> get_cellset(String cellset_name);

    /*
     * \brief Returns true if "meshes" is empty.
     */
    //bool no_meshes() const;

    /*
     * \brief Returns true if "cellsets" is empty.
     */
    //bool no_cellsets() const;

  private:

    MeshReader::CellSetParent* _find_cell_set_parent(String parent_name);
    MeshReader::MeshNode* _find_sub_mesh_parent(String parent_name);

    /**
     * \brief Parses a mesh-header-data-input stream.
     *
     * \param[in] ifs
     * A reference to the input stream to be parsed.
     *
     * \param[in] cur_line
     * A counter that specifies the number of "lines" read so far.
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
     * A counter that specifies the number of "lines" read so far.
     *
     * \returns
     * The number of lines that has been read when the programme is done.
     */
    Index _parse_info_section(Index cur_line, std::istream& ifs);

    /**
     * \brief Parses a mesh-section-data-input stream.
     *
     * \param[in] ifs
     * A reference to the input stream to be parsed.
     *
     * \param[in] cur_line
     * A counter that specifies the number of "lines" read so far.
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
     * A counter that specifies the number of "lines" read so far.
     *
     * \returns
     * The number of lines that has been read when the programme is done.
     */
    Index _parse_cellset_section(Index cur_line, std::istream& ifs);

  }; // class MeshReader
} // namespace FEAST

#endif // KERNEL_UTIL_MESH_READER_HPP
