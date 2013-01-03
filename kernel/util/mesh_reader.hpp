#pragma once
#ifndef KERNEL_UTIL_MESH_READER_HPP
#define KERNEL_UTIL_MESH_READER_HPP 1

// includes, FEAST
#include <kernel/util/file_error.hpp>

// includes, system
#include <iostream>
#include <vector>

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

    /**
     * \brief Mesh data container class
     *
     * This class stores the data related to a (sub-) mesh.
     *
     * \author Constantin Christof
     */
    class MeshDataContainer
    {

      public:
        // coordinate- and adjacency vectors
        typedef std::vector<double> CoordVec;
        typedef std::vector<CoordVec> CoordStack;
        typedef std::vector<Index> AdjVec;
        typedef std::vector<AdjVec> AdjStack;

      public:

        // basic information
        String name;
        String parent;
        String chart;
        String coord_version;
        String adjacency_version;

        // mesh data
        String mesh_type;
        String shape_type;
        Index coord_per_vertex;
        Index vertex_number;
        Index edge_number;
        Index tria_number;
        Index quad_number;
        Index tetra_number;
        Index hexa_number;

        // file paths
        String coord_path;
        String adj_path;

        // coordinates
        CoordStack coords;

        // adjacencies
        AdjStack adjacencies [4][4];

        // parent indices
        std::vector<Index> parent_indices [4];

      public:
        // default CTOR
        MeshDataContainer()
        {
          CONTEXT("MeshDataContainer::MeshDataContainer()");
        }

        // default DTOR
        MeshDataContainer::~MeshDataContainer()
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
    class CellSetContainer
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
        CellSetContainer()
        {
          CONTEXT("CellSetContainer::CellSetContainer()");
        }

        // default DTOR
        CellSetContainer::~CellSetContainer()
        {
          CONTEXT("CellSetContainer::~CellSetContainer()");
        }
    }; // CellSetContainer

  private:
    // general information
    String version;
    Index number_of_submeshes;
    Index number_of_cellsets;

    // file paths
    String chart_path;

    // mesh data
    std::vector<MeshDataContainer> meshes;

    // cellset data
    std::vector<CellSetContainer> cell_sets;

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
     * \brief Returns the version-string.
     */
    String get_version();

    /**
     * \brief Returns the chart path.
     */
    String get_chart_path();

    /**
     * \brief Returns the number of submeshes.
     */
    Index get_number_of_submeshes();

    /**
     * \brief Returns the number of cellsets.
     */
    Index get_number_of_cellsets();

    /**
     * \brief Returns the mesh data container that belongs to the given name.
     * The returned bool is false if no such mesh was found.
     *
     * \param[in] name
     * The name of the mesh that is needed.
     */
    std::pair<MeshReader::MeshDataContainer, bool> get_mesh(String mesh_name);

    /**
     * \brief Returns the cellset container that belongs to the given name.
     * The returned bool is false if no such set was found.
     *
     * \param[in] name
     * The name of the cellset that is needed.
     */
    std::pair<MeshReader::CellSetContainer, bool> get_cellset(String cellset_name);

    /**
     * \brief Returns true if "meshes" is empty.
     */
    bool no_meshes();

    /**
     * \brief Returns true if "cellsets" is empty.
     */
    bool no_cellsets();

  private:

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
    Index parse_header_section(Index cur_line, std::istream& ifs);

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
    Index parse_info_section(Index cur_line, std::istream& ifs);

    /**
     * \brief Parses a mesh-section-data-input stream.
     *
     * \param[in] ifs
     * A reference to the input stream to be parsed.
     *
     * \param[in] cur_line
     * A counter that specifies the number of "lines" read so far.
     *
     * \param[in] flag
     * A string that specifies if the data belongs to a submesh (flag = "sub")
     * or the root mesh (flag = "root")
     *
     * \returns
     * The number of lines that has been read when the programme is done.
     */
    Index parse_mesh_section(Index cur_line, String flag, std::istream& ifs);

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
    Index parse_cellset_section(Index cur_line, std::istream& ifs);

  }; // class MeshReader
} // namespace FEAST

#endif // KERNEL_UTIL_MESH_READER_HPP