#pragma once
#ifndef KERNEL_GEOMETRY_MESH_NODE_HPP
#define KERNEL_GEOMETRY_MESH_NODE_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_streamer_factory.hpp>
#include <kernel/geometry/macro_factory.hpp>
#include <kernel/geometry/patch_factory.hpp>
#include <kernel/geometry/patch_halo_builder.hpp>
#include <kernel/geometry/patch_meshpart_factory.hpp>
#include <kernel/geometry/patch_meshpart_splitter.hpp>
#include <kernel/geometry/intern/dual_adaptor.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/util/mesh_streamer.hpp>

// includes, STL
#include <map>
#include <deque>
#include <vector>

namespace FEAST
{
  namespace Geometry
  {
    // forward declarations
    /// \cond internal
    namespace Intern
    {
      template<typename MeshType_>
      struct TypeConverter;

      template<int dim_>
      struct AdjacenciesFiller;
    }

    template<typename Policy_>
    class MeshPartNode DOXY({});
    /// \endcond

    /**
     * \brief Adapt mode enumeration
     */
    enum class AdaptMode
    {
      none  = 0x0,
      chart = 0x1,
      dual  = 0x2
    };

    /// \cond internal
    AdaptMode operator|(AdaptMode x, AdaptMode y)
    {
      return static_cast<AdaptMode>(int(x) | int(y));
    }

    AdaptMode operator&(AdaptMode x, AdaptMode y)
    {
      return static_cast<AdaptMode>(int(x) & int(y));
    }
    /// \endcond

    /**
     * \brief Mesh Node base class
     *
     * A MeshNode is a container for bundling a mesh with MeshParts referring to it.
     *
     * \tparam Policy_
     * Bundle of type names for meshes and partial meshes.
     *
     * \tparam MeshNodePolicy_
     * Bundle of type names for other MeshNodes containing meshes and MeshParts.
     *
     * \author Peter Zajac
     */
    template<
      typename RootMesh_,
      typename ThisMesh_>
    class MeshNode
    {
    public:
      /// the root mesh type, has to be a Mesh
      typedef RootMesh_ RootMeshType;
      /// mesh type of this node, can be a Mesh or MeshPart
      typedef ThisMesh_ MeshType;
      /// the mesh part type
      typedef MeshPart<RootMeshType> MeshPartType;
      /// the mesh part node type
      typedef MeshPartNode<RootMeshType> MeshPartNodeType;
      /// the mesh atlas type
      typedef MeshAtlas<RootMeshType> MeshAtlasType;
      /// the mesh chart type
      typedef Atlas::ChartBase<RootMeshType> MeshChartType;

    protected:
      /**
       * \brief Container class for bundling MeshPartNodes with their corresponding charts
       */
      class MeshPartNodeBin
      {
      public:
        /// This container's MeshPartNode
        MeshPartNodeType* node;
        /// Chart belonging to node
        const MeshChartType* chart;

      public:
        /**
         * \brief Constructor
         *
         * \param node_
         * MeshPartNode for this container.
         *
         * \param chart_
         * Chart for this container.
         *
         */
        explicit MeshPartNodeBin(MeshPartNodeType* node_, const MeshChartType* chart_) :
          node(node_),
          chart(chart_)
        {
        }
      }; // class MeshPartNodeBin

      /// submesh node bin container type
      typedef std::map<String, MeshPartNodeBin> MeshPartNodeContainer;
      /// submesh node iterator type
      typedef typename MeshPartNodeContainer::iterator MeshPartNodeIterator;
      /// submesh node const-iterator type
      typedef typename MeshPartNodeContainer::const_iterator MeshPartNodeConstIterator;
      /// submesh node reverse-iterator type
      typedef typename MeshPartNodeContainer::reverse_iterator MeshPartNodeReverseIterator;

    protected:
      /// a pointer to the mesh of this node
      MeshType* _mesh;
      /// child submesh nodes
      MeshPartNodeContainer _mesh_part_nodes;

    protected:
      /**
       * \brief Constructor.
       *
       * \param[in] mesh
       * A pointer to the mesh for this node.
       */
      explicit MeshNode(MeshType* mesh) :
        _mesh(mesh),
        _mesh_part_nodes()
      {
        CONTEXT(name() + "::MeshNode()");
      }

    public:
      /**
       * \brief Virtual destructor
       *
       * \note _mesh gets deleted because passing it to the MeshNode's constructor passes ownership to it.
       */
      virtual ~MeshNode()
      {
        CONTEXT(name() + "::~MeshNode()");

        // loop over all submesh nodes in reverse order and delete them
        MeshPartNodeReverseIterator it(_mesh_part_nodes.rbegin());
        MeshPartNodeReverseIterator jt(_mesh_part_nodes.rend());
        for(; it != jt; ++it)
        {
          if(it->second.node != nullptr)
          {
            delete it->second.node;
          }
        }

        // delete mesh
        if(_mesh != nullptr)
        {
          delete _mesh;
        }
      }

      /// \returns The size of dynamically allocated memory in bytes.
      std::size_t bytes() const
      {
        std::size_t s(0);
        if(_mesh != nullptr)
          s += _mesh->bytes();
        for(auto it = _mesh_part_nodes.begin(); it != _mesh_part_nodes.end(); ++it)
          s += it->second.node->bytes();
        return s;
      }

      /**
       * \brief Returns the mesh of this node.
       * \returns
       * A pointer to the mesh contained in this node.
       */
      MeshType* get_mesh()
      {
        return _mesh;
      }

      /**
       * \brief Gets the identifier
       *
       * \returns The name of the mesh or meshpart contained in this node.
       */
      String get_identifier() const
      {
        if(_mesh != nullptr)
          return _mesh->get_identifier();
        return "";
      }

      /** \copydoc get_mesh() */
      const MeshType* get_mesh() const
      {
        return _mesh;
      }

      /**
       * \brief Returns the names of all mesh parts of this node.
       */
      std::deque<String> get_mesh_part_names() const
      {
        std::deque<String> names;
        for(auto it = _mesh_part_nodes.begin(); it != _mesh_part_nodes.end(); ++it)
        {
          names.push_back((*it).first);
        }
        return names;
      }

      /**
       * \brief Adds a new mesh-part child node.
       *
       * \param[in] part_name
       * The name of the child node.
       *
       * \param[in] mesh_part_node
       * A pointer to the mesh_part node to be added.
       *
       * \param[in] chart
       * A pointer to the chart that the subnode is to be associated with. May be \c nullptr.
       *
       * \returns
       * \p mesh_part_node if the insertion was successful, otherwise \c nullptr.
       *
       * \note Because mesh_part_node gets added to (this) with identifier part_name, the identifier and
       * parent_identifier get overwritten with part_name and this->get_identifier() respectively.
       */
      MeshPartNodeType* add_mesh_part_node(
        const String& part_name,
        MeshPartNodeType* mesh_part_node,
        const MeshChartType* chart = nullptr)
      {
        CONTEXT(name() + "::add_mesh_part_node()");
        if(mesh_part_node != nullptr)
        {
          if(_mesh_part_nodes.insert(std::make_pair(part_name, MeshPartNodeBin(mesh_part_node, chart))).second)
          {
            // Set identifier in the MeshPartNode
            mesh_part_node->set_identifier(part_name);

            // Set parent_identifier in the MeshPartNode
            mesh_part_node->set_parent_identifier(this->get_identifier());

            return mesh_part_node;
          }
        }
        return nullptr;
      }

      /**
       * \brief Adds a new mesh-part child node.
       *
       * \param[in] part_name
       * The name of the child node.
       *
       * \param[in] mesh_part
       * A pointer to the mesh_part to be added.
       *
       * \param[in] chart
       * A pointer to the chart that the subnode is to be associated with. May be \c nullptr.
       *
       * \returns
       * A pointer to the newly created mesh-part node if the insertion was successful, otherwise \c nullptr.
       */
      MeshPartNodeType* add_mesh_part(
        const String& part_name,
        MeshPartType* mesh_part,
        const MeshChartType* chart = nullptr)
      {
        CONTEXT(name() + "::add_mesh_part()");
        MeshPartNodeType* part_node = new MeshPartNodeType(mesh_part);
        if(add_mesh_part_node(part_name, part_node, chart) == nullptr)
        {
          delete part_node;
          return nullptr;
        }

        part_node->set_parent_identifier(this->get_identifier());

        return part_node;
      }
      /**
       * \brief Sets the chart for a particular mesh part.
       *
       * \param[in] part_name
       * The name of the mesh part that the chart is to be assigned to.
       *
       * \param[in] chart
       * The chart that is to be assigned. May also be \c nullptr.
       *
       * \returns
       * \c true if the chart was assigned successfully or \c false if no mesh part
       * with the specified name was found.
       */
      bool set_mesh_part_chart(const String& part_name, const MeshChartType* chart)
      {
        MeshPartNodeIterator it(_mesh_part_nodes.find(part_name));
        if(it == _mesh_part_nodes.end())
          return false;
        it->second.chart = chart;
        return true;
      }

      /**
       * \brief Searches this container for a MeshPartNode
       *
       * \param[in] part_name
       * The name of the node to be found.
       *
       * \returns
       * A pointer to the mesh_part node associated with \p part_name or \c nullptr if no such node was found.
       */
      MeshPartNodeType* find_mesh_part_node(const String& part_name)
      {
        CONTEXT(name() + "::find_mesh_part_node()");
        MeshPartNodeIterator it(_mesh_part_nodes.find(part_name));
        return (it != _mesh_part_nodes.end()) ? it->second.node : nullptr;
      }

      /** \copydoc find_mesh_part_node() */
      const MeshPartNodeType* find_mesh_part_node(const String& part_name) const
      {
        CONTEXT(name() + "::find_mesh_part_node() [const]");
        MeshPartNodeConstIterator it(_mesh_part_nodes.find(part_name));
        return (it != _mesh_part_nodes.end()) ? it->second.node : nullptr;
      }

      /**
       * \brief Searches this container for a MeshPart
       *
       * \param[in] part_name
       * The name of the node associated with the mesh_part to be found.
       *
       * \returns
       * A pointer to the mesh_part associated with \p part_name or \c nullptr if no such node was found.
       */
      MeshPartType* find_mesh_part(const String& part_name)
      {
        CONTEXT(name() + "::find_mesh_part()");
        MeshPartNodeType* node = find_mesh_part_node(part_name);
        return (node != nullptr) ? node->get_mesh() : nullptr;
      }

      /** \copydoc find_mesh_part() */
      const MeshPartType* find_mesh_part(const String& part_name) const
      {
        CONTEXT(name() + "::find_mesh_part() [const]");
        const MeshPartNodeType* node = find_mesh_part_node(part_name);
        return (node != nullptr) ? node->get_mesh() : nullptr;
      }

      /**
       * \brief Searches for a chart belonging to a MeshPart by name
       *
       * \param[in] part_name
       * The name of the node associated with the chart to be found.
       *
       * \returns
       * A pointer to the chart associated with \p part_name of \c nullptr if no such node was found or if
       * the corresponding node did not have a chart.
       */
      MeshChartType* find_mesh_part_chart(const String& part_name)
      {
        CONTEXT(name() + "::find_mesh_part_chart()");
        MeshPartNodeIterator it(_mesh_part_nodes.find(part_name));
        return (it != _mesh_part_nodes.end()) ? it->second.chart : nullptr;
      }

      /** \copydoc find_mesh_part_chart() */
      const MeshChartType* find_mesh_part_chart(const String& part_name) const
      {
        CONTEXT(name() + "::find_mesh_part_chart() [const]");
        MeshPartNodeConstIterator it(_mesh_part_nodes.find(part_name));
        return (it != _mesh_part_nodes.end()) ? it->second.chart : nullptr;
      }

      /**
       * \brief Adapts this mesh node.
       *
       * This function loops over all MeshPart nodes and uses their associated charts (if given)
       * to adapt the mesh in this node.
       *
       * \param[in] recursive
       * If set to \c true, all mesh_part nodes are adapted prior to adapting this node.
       */
      void adapt(bool recursive = true)
      {
        CONTEXT(name() + "::adapt()");

        // loop over all mesh_part nodes
        MeshPartNodeIterator it(_mesh_part_nodes.begin());
        MeshPartNodeIterator jt(_mesh_part_nodes.end());
        for(; it != jt; ++it)
        {
          // adapt child node
          if(recursive)
          {
            it->second.node->adapt(true);
          }

          // adapt this node if a chart is given
          if(it->second.chart != nullptr)
          {
            it->second.chart->adapt(*_mesh, *(it->second.node->get_mesh()));
          }
        }
      }

      /**
       * \brief Adapts this mesh node.
       *
       * This function adapts this node by a specific chart whose name is given.
       *
       * \param[in] part_name
       * The name of the mesh_part node that is to be used for adaption.
       *
       * \param[in] recursive
       * If set to \c true, the mesh_part node associated with \p part_name will be adapted prior to adapting this node.
       *
       * \returns
       * \c true if this node was adapted successfully or \c false if no node is associated with \p part_name or if
       * the node did not contain any chart.
       */
      bool adapt_by_name(const String& part_name, bool recursive = false)
      {
        CONTEXT(name() + "::adapt_by_name()");

        // try to find the corresponding mesh_part node
        MeshPartNodeIterator it(_mesh_part_nodes.find(part_name));
        if(it == _mesh_part_nodes.end())
          return false;

        // adapt child node
        if(recursive)
        {
          it->second.node->adapt(true);
        }

        // adapt this node
        if(it->second.chart != nullptr)
        {
          it->second.chart->adapt(*_mesh, *(it->second.node->get_mesh()));
          return true;
        }

        // no chart associated
        return false;
      }

      /**
       * \brief Returns the name of the class.
       * \returns
       * The name of the class as a String.
       */
      static String name()
      {
        return "MeshNode<...>";
      }

    protected:
      /**
       * \brief Refines all child nodes of this node.
       *
       * \param[in,out] refined_node
       * A reference to the node generated by refining this node.
       */
      void refine_children(MeshNode& refined_node) const
      {
        // refine mesh_partes
        refine_mesh_parts(refined_node);
      }

      /**
       * \brief Refines all child MeshPart nodes of this node.
       *
       * \param[in,out] refined_node
       * A reference to the node generated by refining this node.
       */
      void refine_mesh_parts(MeshNode& refined_node) const
      {
        MeshPartNodeConstIterator it(_mesh_part_nodes.begin());
        MeshPartNodeConstIterator jt(_mesh_part_nodes.end());

        for(; it != jt; ++it)
        {
          refined_node.add_mesh_part_node(it->first, it->second.node->refine(*_mesh), it->second.chart);
        }
      }

    }; // class MeshNode

    /* ***************************************************************************************** */

    /**
     * \brief MeshPart mesh tree node class template
     *
     * This class template implements a mesh tree node containing a MeshPart.
     *
     * \author Peter Zajac
     */
    template<typename RootMesh_>
    class MeshPartNode
      : public MeshNode<RootMesh_, MeshPart<RootMesh_>>
    {
    public:
      /// base class typedef
      typedef MeshNode<RootMesh_, MeshPart<RootMesh_>> BaseClass;

      /// this mesh type
      typedef typename BaseClass::MeshType MeshType;
      /// mesh part type
      typedef typename BaseClass::MeshPartType MeshPartType;
      /// mesh atlas type
      typedef typename BaseClass::MeshAtlasType MeshAtlasType;
      /// mesh chart type
      typedef typename BaseClass::MeshChartType MeshChartType;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] mesh_part
       * A pointer to the MeshPart for this node.
       */
      explicit MeshPartNode(MeshPartType* mesh_part) :
        BaseClass(mesh_part)
      {
        CONTEXT(name() + "::MeshPartNode()");
      }

      /**
       * \brief Constructs a MeshPartNode from a MeshStreamer::MeshNode and adds all MeshParts therein
       *
       * \param streamer_node
       * Container holding streamed mesh data.
       *
       * \param atlas
       * MeshAtlas holding all charts.
       *
       * If a mesh part from the streamer node specifies a chart, it is meant to be contained in atlas and added to
       * the (sub) MeshPartNode that gets added to (this).
       *
       */
      explicit MeshPartNode(MeshStreamer::MeshNode& streamer_node, MeshAtlasType* atlas = nullptr) :
        BaseClass(nullptr)
      {
        // create a factory for our mesh
        MeshStreamerFactory<MeshPartType> my_factory(&streamer_node.mesh_data);

        // create our mesh
        this->_mesh = new MeshType(my_factory);

        // now loop over all mesh parts
        for(auto& it : streamer_node.meshpart_map)
        {
          // get the node of this mesh part
          MeshStreamer::MeshNode& streamer_subnode = *it.second;

          // create its node
          MeshPartNode<RootMesh_>* my_mesh_part_node = new MeshPartNode<RootMesh_>(streamer_subnode, atlas);

          // check for a chart to assign
          MeshChartType* chart = nullptr;
          if(atlas != nullptr)
          {
            String chart_name = streamer_subnode.mesh_data.chart;
            if(!chart_name.empty())
            {
              chart = atlas->find_mesh_chart(chart_name);
              if(chart == nullptr)
                throw InternalError("Chart not found!");
            }
          }

          // Add the new MeshPartNode to this node
          this->add_mesh_part_node(it.first, my_mesh_part_node, chart);
        }
      }

      /// virtual destructor
      virtual ~MeshPartNode()
      {
        CONTEXT(name() + "::~MeshPartNode()");
      }

      void set_identifier(const String& id)
      {
        if (this->_mesh != nullptr)
          this->_mesh->set_identifier(id);
      }

      /**
      * \brief Gets the parent's identifier
      *
      * \returns The parent's identifer
      */
      String get_parent_identifier() const
      {
        if(this->_mesh != nullptr)
          return this->_mesh->get_parent_identifier();
        return "";
      }

      /**
      * \brief Sets the parent_identifier
      *
      * \param[in] parent_identifier
      * New parent identifier.
      *
      * This is used to keep the identifiers consistent when adding MeshNodes as children to other MeshNodes.
      *
      */
      void set_parent_identifier(String parent_identifier)
      {
        if(this->_mesh != nullptr)
          this->_mesh->set_parent_identifier(parent_identifier);
      }

      /**
       * \brief Refines this node and its sub-tree.
       *
       * \param[in] parent
       * A reference to the parent mesh/cell subset of this node's cell subset.
       *
       * \returns
       * A pointer to a MeshPartNode containing the refined cell subset tree.
       */
      template<typename ParentType_>
      MeshPartNode* refine(const ParentType_& parent) const
      {
        CONTEXT(name() + "::refine()");

        // the mesh part may be a nullptr; in this case also return a node containing a nullptr
        if(this->_mesh == nullptr)
        {
          return new MeshPartNode(nullptr);
        }

        // create a refinery
        StandardRefinery<MeshPartType> refinery(*this->_mesh, parent);

        // create a new MeshPartNode
        MeshPartNode* fine_node = new MeshPartNode(new MeshPartType(refinery));

        // refine our children
        refine_mesh_parts(*fine_node);

        // okay
        return fine_node;
      }

      /**
       * \brief Returns the name of the class.
       * \returns
       * The name of the class as a String.
       */
      static String name()
      {
        return "MeshPartNode<...>";
      }

    protected:
      /**
       * \brief Refines this node's child nodes.
       *
       * \note This function is called by this node's refine() function, therefore there is no need
       * for the user to call this function.
       *
       * \param[in,out] refined_node
       * A reference to the node generated from refinement of this node.
       */
      void refine_mesh_parts(MeshPartNode& refined_node) const
      {
        CONTEXT(name() + "::refine_mesh_parts()");

        typename BaseClass::MeshPartNodeConstIterator it(this->_mesh_part_nodes.begin());
        typename BaseClass::MeshPartNodeConstIterator jt(this->_mesh_part_nodes.end());
        for(; it != jt; ++it)
        {
          refined_node.add_mesh_part_node(it->first, it->second.node->refine(*this->_mesh));
        }
      }

    }; // class MeshPartNode<...>

    /* ***************************************************************************************** */

    /**
     * \brief Root mesh node class template
     *
     * This class template is used for the root node of a mesh tree.
     *
     * \author Peter Zajac
     */
    template<typename RootMesh_>
    class RootMeshNode
      : public MeshNode<RootMesh_, RootMesh_>
    {
    public:
      /// base class typedef
      typedef MeshNode<RootMesh_, RootMesh_> BaseClass;

      /// this mesh type
      typedef typename BaseClass::MeshType MeshType;
      /// mesh part type
      typedef typename BaseClass::MeshPartType MeshPartType;
      /// mesh atlas type
      typedef typename BaseClass::MeshAtlasType MeshAtlasType;
      /// mesh chart type
      typedef typename BaseClass::MeshChartType MeshChartType;

    protected:
      /// our atlas
      MeshAtlasType* _atlas;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] mesh
       * A pointer to the mesh for this node.
       *
       * \param[in] atlas
       * A pointer to the atlas for this mesh tree.
       */
      explicit RootMeshNode(MeshType* mesh, MeshAtlasType* atlas = nullptr) :
        BaseClass(mesh),
        _atlas(atlas)
      {
      }

      /**
       * \brief Constructs a RootMeshNode from a streamed mesh
       *
       * \param[in] mesh_reader
       * MeshStreamer that contains the information from the streamed mesh.
       *
       * \param atlas
       * MeshAtlas to be added to the new RootMeshNode, can be nullptr.
       *
       */
      explicit RootMeshNode(MeshStreamer& mesh_reader, MeshAtlasType* atlas = nullptr) :
        BaseClass(nullptr),
        _atlas(atlas)
      {
        MeshStreamer::MeshNode& streamer_node = *mesh_reader.get_root_mesh_node();

        // create a factory for our mesh
        MeshStreamerFactory<MeshType> my_factory(&streamer_node.mesh_data);

        // create our mesh
        this->_mesh = new MeshType(my_factory);

        // now loop over all mesh parts
        for(auto& it : streamer_node.meshpart_map)
        {
          // get the node of this mesh part
          MeshStreamer::MeshNode& streamer_subnode = *it.second;

          // create its node
          MeshPartNode<RootMesh_>* my_mesh_part_node = new MeshPartNode<RootMesh_>(streamer_subnode, atlas);

          // check for a chart to assign
          MeshChartType* chart = nullptr;
          if(atlas != nullptr)
          {
            String chart_name = streamer_subnode.mesh_data.chart;
            if(!chart_name.empty())
            {
              chart = atlas->find_mesh_chart(chart_name);
              if(chart == nullptr)
                throw InternalError("Chart not found!");
            }
          }

          // Add the new MeshPartNode to this node
          this->add_mesh_part_node(it.first, my_mesh_part_node, chart);
        }
      }

      /**
       * \brief Creates a MeshStreamer object for writing
       *
       * \returns MeshStreamer containing all data to be written.
       *
       */
      MeshStreamer* create_mesh_writer()
      {
        MeshStreamer* mesh_writer(new MeshStreamer);

        // Write MeshStreamer part
        // Set these to 0 because they increase when we add MeshParts/MeshCharts
        mesh_writer->_num_charts = 0;
        mesh_writer->_num_meshparts = 0;
        mesh_writer->_info = "Written by MeshNode::create_mesh_writer";
        // Create MeshStreamer::MeshNode for holding the root mesh
        mesh_writer->_root_mesh_node = new MeshStreamer::MeshNode;
        MeshWriter<RootMesh_>::write_mesh(mesh_writer->_root_mesh_node->mesh_data, *(this->get_mesh()));

        // Add MeshStreamer::MeshNodes for holding the MeshParts
        for(auto& it : this->_mesh_part_nodes)
        {
          auto* mpart = it.second.node->get_mesh();
          if(mpart == nullptr)
            continue;

          MeshStreamer::MeshNode* mesh_part(new MeshStreamer::MeshNode);
          MeshWriter<RootMesh_>::write_mesh_part(mesh_part->mesh_data, *(mpart));

          // Add the chart (if any)
          auto* chart(it.second.chart);
          if(chart != nullptr)
          {
            // Get chart
            mesh_part->mesh_data.chart = it.first;

            MeshStreamer::ChartContainer* chart_data(new MeshStreamer::ChartContainer);

            MeshWriter<RootMesh_>::write_chart(*chart_data, it.first, *(it.second.chart));

            // Insert copies to an std::vector...
            mesh_writer->insert_chart(*chart_data);
            // ... so we need to clean up afterwards
            delete chart_data;
          }

          // Add the new mesh part to the MeshStreamer
          mesh_writer->_insert_meshpart(mesh_part);
        }

        return mesh_writer;
      }

      /**
       * \brief Virtual destructor
       *
       * \note _atlas does not get deleted because it is not owned by the RootMeshNode, as several MeshNodes may refer
       * to the same atlas.
       */
      virtual ~RootMeshNode()
      {
      }

      /**
       * \brief Creates mesh-parts for all base cells of this root mesh.
       *
       * This function must be called for the unrefined base root mesh, if one
       * intends to use the Assembly::BaseSplitter class for the redistribution
       * of data in dynamic load balancing.
       */
      void create_base_splitting()
      {
        // make sure we have a mesh
        if(this->_mesh == nullptr)
          throw InternalError("No mesh assigned");

        // get number of cells
        Index num_cells = this->_mesh->get_num_entities(MeshType::shape_dim);

        // loop over all cell indices
        for(Index cell(0); cell < num_cells; ++cell)
        {
          // create a macro factory
          MacroFactory<MeshPartType> factory(*this->_mesh, cell);
          // create a new mesh part and add to this mesh node
          this->add_mesh_part("_base:" + stringify(cell), new MeshPartType(factory));
        }
      }

      /**
       * \brief Refines this node and its sub-tree.
       *
       * \param[in] adapt_mode
       * Mode for adaption, defaults to chart.
       *
       * \returns
       * A pointer to a RootMeshNode containing the refined mesh tree.
       */
      RootMeshNode* refine(AdaptMode adapt_mode = AdaptMode::chart) const
      {
        // create a refinery
        StandardRefinery<MeshType> refinery(*this->_mesh);

        // create a new root mesh node
        RootMeshNode* fine_node = new RootMeshNode(new MeshType(refinery), this->_atlas);

        // refine our children
        this->refine_children(*fine_node);

        // adapt by chart?
        if((adapt_mode & AdaptMode::chart) != AdaptMode::none)
        {
          fine_node->adapt(true);
        }

        // adapt dual?
        if((adapt_mode & AdaptMode::dual) != AdaptMode::none)
        {
          Intern::DualAdaptor<MeshType>::adapt(*fine_node->get_mesh(), *this->get_mesh());
        }

        // okay
        return fine_node;
      }

      /**
       * \brief Extracts a patch from the root mesh as a new mesh node
       *
       * \param[in] rank
       * The rank of the patch to be created.
       *
       * \param[in] ranks_at_elem
       * A graph specifies the ranks for each element of the base mesh
       *
       * \param[in] comm_ranks
       * A vector of all ranks adjacent to the current rank
       *
       * \returns
       * A new mesh node representing the extracted patch.
       */
      RootMeshNode* extract_patch(
        const Index rank,
        const Adjacency::Graph& ranks_at_elem,
        const std::vector<Index>& comm_ranks)
      {
        // get base root mesh
        MeshType& base_root_mesh = *this->get_mesh();

        // get all mesh part names
        std::deque<String> part_names = this->get_mesh_part_names();

        // Step 1: create mesh part of the patch
        MeshPartType* patch_mesh_part = nullptr;
        {
          PatchMeshPartFactory<MeshType> part_factory(rank, base_root_mesh, ranks_at_elem);

          // create patch mesh part
          patch_mesh_part = new MeshPartType(part_factory);
          patch_mesh_part->template deduct_target_sets_from_top<MeshType::shape_dim>(
            base_root_mesh.get_index_set_holder());

          /// \todo really add to this mesh node?
          this->add_mesh_part(String("_patch:") + stringify(rank), patch_mesh_part);
        }

        // Step 2: Create root mesh of partition by using PatchFactory
        MeshType* patch_root_mesh = nullptr;
        {
          // create patch root mesh
          PatchFactory<MeshType> patch_factory(base_root_mesh, *patch_mesh_part);
          patch_root_mesh = new MeshType(patch_factory);
        }

        // Step 3: Create root mesh node for mesh
        RootMeshNode* patch_node = new RootMeshNode(patch_root_mesh, this->_atlas);

        // Step 4: intersect boundary and other base mesh parts
        {
          // create mesh part splitter
          PatchMeshPartSplitter<MeshType> part_splitter(base_root_mesh, *patch_mesh_part);

          // loop over all base-mesh mesh parts
          for(auto it = part_names.begin(); it != part_names.end(); ++it)
          {
            // get base mesh part node
            auto* base_part_node = this->find_mesh_part_node(*it);
            if(base_part_node == nullptr)
              throw InternalError(String("base-mesh part '") + (*it) + "' not found!");

            // our split mesh part
            MeshPartType* split_part = nullptr;

            // get base-mesh part
            MeshPartType* base_part = base_part_node->get_mesh();

            // build
            if((base_part != nullptr) && part_splitter.build(*base_part))
            {
              // create our mesh part
              split_part = new MeshPartType(part_splitter);
            }

            // insert patch mesh part
            patch_node->add_mesh_part(*it, split_part);
          }
        }

        // Step 5: Create halos
        {
          // create halo builder
          PatchHaloBuilder<MeshType> halo_builder(ranks_at_elem, base_root_mesh, *patch_mesh_part);

          // loop over all comm ranks
          for(auto it = comm_ranks.begin(); it != comm_ranks.end(); ++it)
          {
            // build halo
            halo_builder.build(*it);

            // create halo mesh part
            MeshPartType* halo_part = new MeshPartType(halo_builder);

            // insert into patch mesh node
            patch_node->add_mesh_part(String("_halo:") + stringify(*it), halo_part);
          }
        }

        // return patch mesh part
        return patch_node;
      }

      /**
       * \brief Returns the name of the class.
       * \returns
       * The name of the class as a String.
       */
      static String name()
      {
        return "RootMeshNode<...>";
      }
    }; // class RootMeshNode
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_MESH_NODE_HPP
