#pragma once
#ifndef KERNEL_GEOMETRY_MESH_NODE_HPP
#define KERNEL_GEOMETRY_MESH_NODE_HPP 1

// includes, FEAST
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_part_node.hpp>
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
     */
    template<typename Shape_>
    struct StandardConformalMeshNodePolicy
    {
      /// root mesh type
      typedef ConformalMesh<Shape_> RootMeshType;
      /// chart for the root mesh
      typedef DummyChart RootMeshChartType;

      /// chart for the sub-mesh
      /// MeshPart type
      typedef MeshPart<RootMeshType> MeshPartType;
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

    // helper policy template for SubMeshNode class template
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
     * \author Peter Zajac
     */
    template<
      typename Policy_,
      typename MeshNodePolicy_>
    class MeshNode
      : public MeshPartParent<Policy_>
    {
    public:
      /// base class typedef
      typedef MeshPartParent<Policy_> BaseClass;

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
       * \brief MeshPartNode bin class
       */
      class MeshPartNodeBin
      {
      public:
        MeshPartNodeType* node;
        const MeshChartType* chart;

      public:
        explicit MeshPartNodeBin(MeshPartNodeType* node_, const MeshChartType* chart_) :
          node(node_),
          chart(chart_)
        {
        }
      };

      /// submesh node bin container type
      typedef std::map<Index,MeshPartNodeBin> MeshPartNodeContainer;
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
        BaseClass(),
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
       * \param[in] id
       * The id of the child node.
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
        Index id,
        MeshPartNodeType* mesh_part_node,
        const MeshChartType* chart = nullptr)
      {
        CONTEXT(name() + "::add_mesh_part_node()");
        if(mesh_part_node != nullptr)
        {
          if(_mesh_part_nodes.insert(std::make_pair(id,MeshPartNodeBin(mesh_part_node, chart))).second)
          {
            return mesh_part_node;
          }
        }
        return nullptr;
      }

      /**
       * \brief Searches for a mesh_part node.
       *
       * \param[in] id
       * The id of the node to be found.
       *
       * \returns
       * A pointer to the mesh_part node associated with \p id or \c nullptr if no such node was found.
       */
      MeshPartNodeType* find_mesh_part_node(Index id)
      {
        CONTEXT(name() + "::find_mesh_part_node()");
        MeshPartNodeIterator it(_mesh_part_nodes.find(id));
        return (it != _mesh_part_nodes.end()) ? it->second.node : nullptr;
      }

      /** \copydoc find_mesh_part_node() */
      const MeshPartNodeType* find_mesh_part_node(Index id) const
      {
        CONTEXT(name() + "::find_mesh_part_node() [const]");
        MeshPartNodeConstIterator it(_mesh_part_nodes.find(id));
        return (it != _mesh_part_nodes.end()) ? it->second.node : nullptr;
      }

      /**
       * \brief Searches for a mesh_part.
       *
       * \param[in] id
       * The id of the node associated with the mesh_part to be found.
       *
       * \returns
       * A pointer to the mesh_part associated with \p id or \c nullptr if no such node was found.
       */
      MeshPartType* find_mesh_part(Index id)
      {
        CONTEXT(name() + "::find_mesh_part()");
        MeshPartNodeType* node = find_mesh_part_node(id);
        return (node != nullptr) ? node->get_mesh() : nullptr;
      }

      /** \copydoc find_mesh_part() */
      const MeshPartType* find_mesh_part(Index id) const
      {
        CONTEXT(name() + "::find_mesh_part() [const]");
        const MeshPartNodeType* node = find_mesh_part_node(id);
        return (node != nullptr) ? node->get_mesh() : nullptr;
      }

      /**
       * \brief Searches for a mesh_part chart.
       *
       * \param[in] id
       * The id of the node associated with the chart to be found.
       *
       * \returns
       * A pointer to the chart associated with \p id of \c nullptr if no such node was found or if
       * the corresponding node did not have a chart.
       */
      MeshChartType* find_mesh_part_chart(Index id)
      {
        CONTEXT(name() + "::find_mesh_part_chart()");
        MeshPartNodeIterator it(_mesh_part_nodes.find(id));
        return (it != _mesh_part_nodes.end()) ? it->second.chart : nullptr;
      }

      /** \copydoc find_mesh_part_chart() */
      const MeshChartType* find_mesh_part_chart(Index id) const
      {
        CONTEXT(name() + "::find_mesh_part_chart() [const]");
        MeshPartNodeConstIterator it(_mesh_part_nodes.find(id));
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
       * This function adapts this node by a specific chart whose id is given.
       *
       * \param[in] id
       * The id of the mesh_part node that is to be used for adaption.
       *
       * \param[in] recursive
       * If set to \c true, the mesh_part node associated with \p id will be adapted prior to adapting this node.
       *
       * \returns
       * \c true if this node was adapted successfully or \c false if no node is associated with \p id or if
       * the node did not contain any chart.
       */
      bool adapt_by_id(Index id, bool recursive = false)
      {
        CONTEXT(name() + "::adapt_by_id()");

        // try to find the corresponding mesh_part node
        MeshPartNodeIterator it(_mesh_part_nodes.find(id));
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
        // Careful: MeshStreamer::MeshPartNodes and MeshStreamer::MeshPartNodes use Strings as their names and
        // identifiers, whereas the RootMeshNode uses an Index for this
        Index i(0);
        for(auto& it:root->sub_mesh_map)
        {
          // Create a factory for the MeshPart
          MeshStreamerFactory<typename Policy_::MeshPartType> mesh_part_factory(mesh_reader, it.first);
          // Construct the MeshPart using that factory
          typename Policy_::MeshPartType* my_mesh_part(new typename Policy_::MeshPartType(mesh_part_factory));
          // Create the MeshPartNode using that MeshPart
          MeshPartNode<Policy_>* my_mesh_part_node(new MeshPartNode<Policy_>(my_mesh_part));
          // Add the new MeshPartNode to the RootMeshNode
          this->add_mesh_part_node(i++, my_mesh_part_node);
        }

      }

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

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_MESH_NODE_HPP
