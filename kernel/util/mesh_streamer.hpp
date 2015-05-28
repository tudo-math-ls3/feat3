#pragma once
#ifndef KERNEL_UTIL_MESH_STREAMER_HPP
#define KERNEL_UTIL_MESH_STREAMER_HPP 1

// includes, FEAST
#include <kernel/util/assertion.hpp>
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
      /// Basic information
      String name;
      /// Name of the parent
      String parent;
      /// Info string
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

      /// parent_indices[d](i) = j means that entity i of dimension d in this container corresponds to entity j in the
      /// parent
      std::vector<Index> parent_indices[4];

    public:
      // Default CTOR
      BaseContainer()
      {
        CONTEXT("BaseContainer::BaseContainer()");
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
      /// Name of the chart file that contains an analytic description of the geometry
      String chart;
      /// Mesh data
      MeshType mesh_type;
      /// Shape of the mesh
      ShapeType shape_type;
      /// Coordinates
      CoordStack coords;
      /// Adjacencies: Adjacency[i][j] contains which entities of dimension i lie at entities of dimension j
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
    template<Index dim_>
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

      /**
       * \brief Returns a reference to a given target_data's count of entities of a given dimension
       *
       * Because of different shapes and dimensions, this count has different names according to the situation.
       * I.e. if the mesh is of Hypercube<3> shape, the entities of dimension 2 are quads and their count is
       * quad_count.
       *
       * \warning Mixed meshes are not supported because they are not implemented (yet).
       *
       * \param[in] target_dim
       * Dimension of the entities, i.e. 1 for edges.
       *
       * \param[in,out] target_data
       * Data refering to another MeshDataContainer whose entity count is to be returned, i.e. a CellSubSet
       * refering to its Mesh.
       *
       * \param[in] mesh_data
       * MeshDataContainer that target_data refers to.
       *
       * \returns A reference to the number of entities of dimension target_dim entry in target_data.
       *
       */
      static Index& get_entity_count(Index target_dim, MeshStreamer::BaseContainer* target_data,
      const MeshStreamer::MeshDataContainer* const mesh_data)
      {
        // Determine what the dimension target_dim entities are
        switch(target_dim)
        {
          case 1 :
            return target_data->edge_count;
          case 2 :
            // We have simplices
            if(mesh_data->shape_type == mesh_data->st_tria || mesh_data->shape_type == mesh_data->st_tetra)
              return target_data->tria_count;
            // We have hypercubes
            if(mesh_data->shape_type == mesh_data->st_quad || mesh_data->shape_type == mesh_data->st_hexa)
              return target_data->quad_count;
            // We have a mixed mesh
            if(mesh_data->shape_type == mesh_data->st_tria_quad || mesh_data->shape_type == mesh_data->st_tetra_hexa)
              throw InternalError("No mixed meshes. Please.");
            break;
          case 3 :
            // We have simplices
            if(mesh_data->shape_type == mesh_data->st_tria || mesh_data->shape_type == mesh_data->st_tetra)
              return target_data->tetra_count;
            // We have hypercubes
            if(mesh_data->shape_type == mesh_data->st_quad || mesh_data->shape_type == mesh_data->st_hexa)
              return target_data->hexa_count;
            // We have a mixed mesh
            if(mesh_data->shape_type == mesh_data->st_tria_quad || mesh_data->shape_type == mesh_data->st_tetra_hexa)
              throw InternalError("No mixed meshes. Please.");
            break;
          default:
            throw InternalError("Unknown target_dim!");
            break;
        }
        throw InternalError("Could not determine target entity count!");
      }

      /**
       * \brief Updated parent_indices of a BaseContainer refering to another MeshDataContainer
       *
       * If the BaseContainer contains only vertex information (namely parent_indices[0]) but the dimension of that
       * container is > 0 (i.e. if the dimension was 2, it could hold edge and triangle/quad parent relations),
       * this missing information can be extracted from the MeshDataContainer the BaseContainer refers to.
       *
       * This is implemented in the generic template and referenced from the specialisations below.
       *
       * \param[in,out] target_data
       * Data refering to another MeshDataContainer whose parent_indices are to be updated, i.e. a CellSubSet
       * refering to its Mesh.
       *
       * \param[in] mesh_data
       * MeshDataContainer that target_data refers to.
       *
       * \param[in] target_dim
       * Dimension of the entities, i.e. 1 for edges.
       *
       */
      static void update_parent_data(MeshStreamer::BaseContainer* target_data,
      const MeshStreamer::MeshDataContainer* const mesh_data, Index target_dim = dim_)
      {
        if(target_dim < 1)
          // Nothing could be missing
          return;
        if(target_dim == 2 || target_dim == 3)
          // First call the lower dimensional version
          update_parent_data(target_data, mesh_data, target_dim-1);
        if(target_dim > 3)
          throw InternalError("CellSubSets exist only up to dimension 3!");

        // Get a reference to the number of entities of target_dim in target_data, the shape is read from mesh_data
        Index& entity_count(get_entity_count(target_dim, target_data, mesh_data));

        // If there are no entities of that dimension in target data, generate this information
        if(entity_count == 0)
        {
          ASSERT(target_data->parent_indices[target_dim].size() == 0, "CellSubSet does not have entities of target_dim, but a parent mapping for them!");
          // For each vertex in the parent, this will contain its index in the CellSubSet, or the number of
          // vertices if it is not present
          Index* parent_vertex_to_sub(new Index[mesh_data->vertex_count]);
          for(Index i(0); i < mesh_data->vertex_count; ++i)
            parent_vertex_to_sub[i] = mesh_data->vertex_count;

          for(Index i(0); i < target_data->vertex_count; ++i)
            parent_vertex_to_sub[(target_data->parent_indices[0])[i]] = i;

          // For every entity in the parent, check if all its vertices are in the CellSubSet
          for(Index entity(0); entity < mesh_data->adjacencies[0][target_dim].size(); ++entity)
          {
            bool is_in_sub(true);
            // Check all local vertices
            for(Index j(0); j < ((mesh_data->adjacencies[0][target_dim])[0]).size(); ++j)
            {
              // This is the global index of local vertex j
              Index i = ((mesh_data->adjacencies[0][target_dim])[entity])[j];

              if(parent_vertex_to_sub[i] == mesh_data->vertex_count)
                is_in_sub = false;
            }

            // Add this entity if all vertices were in the CellSubSet
            if(is_in_sub)
              target_data->parent_indices[target_dim].push_back(entity);
          }
          // Update the entity count
          entity_count = target_data->parent_indices[target_dim].size();

          // Clean up
          delete[] parent_vertex_to_sub;
        }
      }// MeshDataContainerUpdater<dim>::update_parent_data()
    };

    template<>
    struct MeshDataContainerUpdater<2>
    {
      template<typename IndexSetHolderType_>
      static void update_from_ish( MeshStreamer::MeshDataContainer* mesh_data, const IndexSetHolderType_& ish)
      {
        // Mixed meshes are not supported yet
        if(mesh_data->shape_type == mesh_data->st_tria_quad || mesh_data->shape_type == mesh_data->st_tetra_hexa)
          throw InternalError("No mixed meshes, please.");

        // The first relevant data is vertex@edge, meaning dimension 1
        if(mesh_data->edge_count == 0)
        {
          ASSERT(mesh_data->adjacencies[0][1].size() == 0, "_mesh_data contains no edges, but adjacency vector is not empty!");
          // Update vertex@edge information
          auto& vert_at_edge = ish.template get_index_set<1,0>();
          // Set number of edges
          mesh_data->edge_count = vert_at_edge.get_num_entities();

          for(Index edge(0); edge < mesh_data->edge_count; ++edge)
          {
            // idx will contain all indices for edge
            std::vector<Index> idx;
            for(Index i(0); i < vert_at_edge.num_indices; ++i)
              idx.push_back(vert_at_edge(edge,i));

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

          for(Index cell(0); cell < mesh_data->tria_count + mesh_data->quad_count; ++cell)
          {
            // idx will contain all indices for edge
            std::vector<Index> idx;
            // The _index_bound for the cell to edges mapping is the number of edges
            for(Index edge(0); edge < edge_at_cell.num_indices; ++edge)
              idx.push_back(edge_at_cell(cell,edge));
            // Add this edge to the global adjacency vector
            mesh_data->adjacencies[1][2].push_back(idx);
          }
        }
      } // MeshDataContainerUpdater<2>::update_from_ish()

      static void update_parent_data(MeshStreamer::BaseContainer* target_data,
      const MeshStreamer::MeshDataContainer* const mesh_data, Index target_dim = Index(2))
      {
        // Just call the generic version
        MeshDataContainerUpdater<0>::update_parent_data(target_data, mesh_data, target_dim);
      }


    }; //MeshDataContainerUpdater<2>

    template<>
    struct MeshDataContainerUpdater<3>
    {
      template<typename IndexSetHolderType_>
      static void update_from_ish(MeshStreamer::MeshDataContainer* mesh_data, const IndexSetHolderType_& ish)
      {
        if(mesh_data->shape_type == mesh_data->st_tria_quad || mesh_data->shape_type == mesh_data->st_tetra_hexa)
          throw InternalError("No mixed meshes, please.");

        // Update 2d data first
        MeshDataContainerUpdater<2>::update_from_ish(mesh_data, ish);

        // Next: Everything face related
        // Check if the data is missing
        if(mesh_data->tria_count + mesh_data->quad_count == 0)
        {
          // vertex@face
          ASSERT(mesh_data->adjacencies[0][2].size() == 0, "_mesh_data has no faces, but vertex adjacency vector is not empty!");
          auto& vert_at_face = ish.template get_index_set<2,0>();

          // Set number of faces
          if( mesh_data->shape_type == mesh_data->st_tetra)
            mesh_data->tria_count = vert_at_face.get_num_entities();

          if( mesh_data->shape_type == mesh_data->st_hexa)
            mesh_data->quad_count = vert_at_face.get_num_entities();

          for(Index face(0); face < mesh_data->tria_count + mesh_data->quad_count; ++face)
          {
            // idx will contain all indices for face
            std::vector<Index> idx;
            for(Index i(0); i < vert_at_face.num_indices; ++i)
              idx.push_back(vert_at_face(face,i));

            // Add this edge to the global adjacency vector
            mesh_data->adjacencies[0][2].push_back(idx);
          }

          // edge@face was set by the 2d routine called earlier

          // Next: edge@cell
          ASSERT(mesh_data->adjacencies[1][3].size() == 0, "_mesh_data has no edges, but adjacency vector for faces is not empty!");

          auto& edge_at_cell= ish.template get_index_set<2,1>();

          // Update edge@cell information
          for(Index cell(0); cell < mesh_data->tetra_count + mesh_data->hexa_count; ++cell)
          {
            // idx will contain all indices for edge
            std::vector<Index> idx;
            // num_indices is the (local!) number of edges per cell
            for(Index edge(0); edge < edge_at_cell.num_indices; ++edge)
              idx.push_back(edge_at_cell(cell,edge));
            // Add this edge to the global adjacency vector
            mesh_data->adjacencies[1][3].push_back(idx);
          }

          // Last: face@cell
          ASSERT(mesh_data->adjacencies[2][3].size() == 0, "_mesh_data has no faces, but adjacency vector for faces is not empty!");
          auto& face_at_cell= ish.template get_index_set<3,2>();

          // Update face@cell information
          for(Index cell(0); cell < mesh_data->tetra_count + mesh_data->hexa_count; ++cell)
          {
            // idx will contain all indices for face
            std::vector<Index> idx;
            // num_indices is the (local!) number of faces per cell
            for(Index face(0); face < face_at_cell.num_indices; ++face)
              idx.push_back(face_at_cell(cell,face));
            // Add this face to the global adjacency vector
            mesh_data->adjacencies[2][3].push_back(idx);
          }

        } // handling faces

      } // MeshDataContainerUpdater<3>::update_from_ish()

      static void update_parent_data(MeshStreamer::BaseContainer* target_data,
      const MeshStreamer::MeshDataContainer* const mesh_data, Index target_dim = Index(3))
      {
        // Just call the generic version
        MeshDataContainerUpdater<0>::update_parent_data(target_data, mesh_data, target_dim);
      }
    }; // MeshDataContainerUpdater<3>
  } // namespace Intern
  /// \cond internal
} // namespace FEAST

#endif // KERNEL_UTIL_MESH_STREAMER_HPP
