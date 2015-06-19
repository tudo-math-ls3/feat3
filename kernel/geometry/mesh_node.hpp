#pragma once
#ifndef KERNEL_GEOMETRY_MESH_NODE_HPP
#define KERNEL_GEOMETRY_MESH_NODE_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_streamer_factory.hpp>
#include <kernel/util/mesh_streamer.hpp>

// includes, STL
#include <map>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    // dummy chart class; will be replaced by 'real' charts later...
    class DummyChart
    {
    public:
      template<typename A_, typename B_>
      void adapt(A_&, const B_&) const
      {
      }
    };

    /// \endcond

    /**
     * \brief Standard MeshNode Policy for ConformalMesh
     *
     * \tparam Shape_
     * The shape tag class for the mesh.
     *
     * This contains type names for meshes, partial meshes that refer to them, and charts for those meshes.
     */
    template<typename Shape_>
    struct StandardConformalMeshNodePolicy
    {
      /// root mesh type
      typedef ConformalMesh<Shape_> RootMeshType;
      /// chart for the root mesh
      typedef DummyChart RootMeshChartType;

      /// Type for partial meshes refering to a root mesh
      typedef MeshPart<RootMeshType> MeshPartType;
      /// Type for charts of partial meshes
      typedef DummyChart MeshPartChartType;

    };

    /// \cond internal
    // helper policy template for RootMeshNode class template
    template<typename MeshNodePolicy_>
    struct RootMeshNodePolicy
    {
      typedef typename MeshNodePolicy_::RootMeshType MeshType;
      typedef typename MeshNodePolicy_::RootMeshChartType ChartType;
    };

    // helper policy template MeshPartNode class template
    template<typename MeshNodePolicy_>
    struct MeshPartNodePolicy
    {
      typedef typename MeshNodePolicy_::MeshPartType MeshType;
      typedef typename MeshNodePolicy_::MeshPartChartType ChartType;
    };

    // forward declarations
    template<typename Policy_>
    class MeshPartNode DOXY({});
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
      typename Policy_,
      typename MeshNodePolicy_>
    class MeshNode
    {
    public:
      /// submesh type
      typedef typename Policy_::MeshPartType MeshPartType;
      /// submesh node type
      typedef MeshPartNode<Policy_> MeshPartNodeType;

      /// mesh type of this node
      typedef typename MeshNodePolicy_::MeshType MeshType;
      /// mesh chart type of this node
      typedef typename MeshNodePolicy_::ChartType MeshChartType;

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
        _mesh(mesh)
      {
        CONTEXT(name() + "::MeshNode()");
      }

    public:
      /// virtual destructor
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

      /**
       * \brief Returns the mesh of this node.
       * \returns
       * A pointer to the mesh contained in this node.
       */
      MeshType* get_mesh()
      {
        return _mesh;
      }

      /** \copydoc get_mesh() */
      const MeshType* get_mesh() const
      {
        return _mesh;
      }

      /**
       * \brief Adds a new mesh_part child node.
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
            return mesh_part_node;
          }
        }
        return nullptr;
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
       * \brief Searches for a chart belonging to a MeshPart by index
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
       * This function loops over all mesh_part nodes and uses their associated charts (if given)
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
            it->second.chart->adapt(*_mesh, *(it->second.node->_mesh));
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
          it->second.chart->adapt(*_mesh, *(it->second.node->_mesh));
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
       * \brief Refines all child mesh_part nodes of this node.
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
     * \brief Root mesh node class template
     *
     * This class template is used for the root node of a mesh tree.
     *
     * \author Peter Zajac
     */
    template<typename Policy_>
    class RootMeshNode
      : public MeshNode<Policy_, RootMeshNodePolicy<Policy_> >
    {
    public:
      /// base class typedef
      typedef MeshNode<Policy_, RootMeshNodePolicy<Policy_> > BaseClass;

      /// the mesh type of this node
      typedef typename BaseClass::MeshType MeshType;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] mesh
       * A pointer to the mesh for this node.
       */
      explicit RootMeshNode(MeshType* mesh) :
        BaseClass(mesh)
      {
      }

      /**
       * \brief Constructs a RootMeshNode from a streamed mesh
       *
       * \param[in] mesh_reader
       * MeshStreamer that contains the information from the streamed mesh.
       *
       */
      explicit RootMeshNode(MeshStreamer& mesh_reader) :
        BaseClass(nullptr)
      {
        // Construct a new Geometry::Mesh using the MeshStreamer and a MeshStreamerFactory
        MeshStreamerFactory<MeshType> my_factory(mesh_reader);
        this->_mesh = new MeshType(my_factory);

        // Generate a MeshStreamer::MeshNode, as this then contains the information about the MeshPartes and
        // MeshParts present in the MeshStreamer
        MeshStreamer::MeshNode* root(mesh_reader.get_root_mesh_node());
        ASSERT_(root != nullptr);

        // Add MeshPartNodes and MeshPartNodes to the new RootMeshNode by iterating over the MeshStreamer::MeshNode
        for(auto& it:root->sub_mesh_map)
        {
          // Create a factory for the MeshPart
          MeshStreamerFactory<typename Policy_::MeshPartType> mesh_part_factory(mesh_reader, it.first);
          // Construct the MeshPart using that factory
          typename Policy_::MeshPartType* my_mesh_part(new typename Policy_::MeshPartType(mesh_part_factory));
          // Create the MeshPartNode using that MeshPart
          MeshPartNode<Policy_>* my_mesh_part_node(new MeshPartNode<Policy_>(my_mesh_part));
          // Add the new MeshPartNode to the RootMeshNode
          this->add_mesh_part_node(it.first, my_mesh_part_node);
        }
      }

      /**
       * \brief Creates a MeshStreamer object for writing
       *
       * \todo Implement this
       *
       * \returns MeshStreamer containing all data to be written.
       *
       */
      MeshStreamer* create_mesh_writer()
      {
        MeshStreamer* mesh_writer(new MeshStreamer);

        return mesh_writer;
      }

      /// virtual destructor
      virtual ~RootMeshNode()
      {
      }

      /**
       * \brief Refines this node and its sub-tree.
       *
       * \returns
       * A pointer to a RootMeshNode containing the refined mesh tree.
       */
      RootMeshNode* refine() const
      {
        // create a refinery
        StandardRefinery<MeshType> refinery(*this->_mesh);

        // create a new root mesh node
        RootMeshNode* fine_node = new RootMeshNode(new MeshType(refinery));

        // refine our children
        this->refine_children(*fine_node);

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
        return "RootMeshNode<...>";
      }
    }; // class RootMeshNode

    /* ***************************************************************************************** */

    /**
     * \brief MeshPart mesh tree node class template
     *
     * This class template implements a mesh tree node containing a MeshPart.
     *
     * \author Peter Zajac
     */
    template<typename Policy_>
    class MeshPartNode
      : public MeshNode<Policy_, MeshPartNodePolicy<Policy_>>
    {
    public:
      /// base class typedef
      typedef MeshNode<Policy_, MeshPartNodePolicy<Policy_>> BaseClass;
      /// policy type
      typedef Policy_ Policy;
      /// cell subset type
      typedef typename Policy_::MeshPartType MeshPartType;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] subset
       * A pointer to the cell subset for this node.
       */
      explicit MeshPartNode(MeshPartType* mesh_part) :
        BaseClass(mesh_part)
      {
        CONTEXT(name() + "::MeshPartNode()");
      }

      /// virtual destructor
      virtual ~MeshPartNode()
      {
        CONTEXT(name() + "::~MeshPartNode()");
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

        // create a refinery
        StandardRefinery<MeshPartType, ParentType_> refinery(*this->_mesh, parent);

        // create a new cell subset node
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
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_MESH_NODE_HPP
