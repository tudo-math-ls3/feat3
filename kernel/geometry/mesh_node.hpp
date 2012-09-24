#pragma once
#ifndef KERNEL_GEOMETRY_MESH_NODE_HPP
#define KERNEL_GEOMETRY_MESH_NODE_HPP 1

// includes, FEAST
#include <kernel/geometry/cell_sub_set_node.hpp>

// includes, STL
#include <map>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond interal
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
      /// mesh policy for the root mesh
      typedef ConformalMeshPolicy<Shape_> RootMeshPolicy;
      /// root mesh type
      typedef ConformalMesh<RootMeshPolicy> RootMeshType;
      /// standard refinery for the root mesh
      typedef StandardRefinery<RootMeshType> RootMeshRefineryType;
      /// chart for the root mesh
      typedef DummyChart RootMeshChartType;

      /// mesh policy for the sub-mesh
      typedef ConformalSubMeshPolicy<Shape_> SubMeshPolicy;
      /// sub-mesh type
      typedef ConformalSubMesh<SubMeshPolicy> SubMeshType;
      /// standard refinery for the sub-mesh
      typedef StandardRefinery<SubMeshType> SubMeshRefineryType;
      /// chart for the sub-mesh
      typedef DummyChart SubMeshChartType;

      /// cell subset type
      typedef CellSubSet<Shape_> CellSubSetType;
      /// standard refinery for the cell subset
      typedef StandardRefinery<CellSubSetType> CellSubSetRefineryType;
    };

    /// \cond internal
    /// helper policy template for RootMeshNode class template
    template<typename MeshNodePolicy_>
    struct RootMeshNodePolicy
    {
      typedef typename MeshNodePolicy_::RootMeshType MeshType;
      typedef typename MeshNodePolicy_::RootMeshRefineryType RefineryType;
      typedef typename MeshNodePolicy_::RootMeshChartType ChartType;
    };

    /// helper policy template for SubMeshNode class template
    template<typename MeshNodePolicy_>
    struct SubMeshNodePolicy
    {
      typedef typename MeshNodePolicy_::SubMeshType MeshType;
      typedef typename MeshNodePolicy_::SubMeshRefineryType RefineryType;
      typedef typename MeshNodePolicy_::SubMeshChartType ChartType;
    };
    /// \endcond

    // forward declarations
    template<typename Policy_>
    class SubMeshNode;

    /**
     * \brief Mesh Node base class
     *
     * \author Peter Zajac
     */
    template<
      typename Policy_,
      typename MeshNodePolicy_>
    class MeshNode
      : public CellSubSetParent<Policy_>
    {
    public:
      /// base class typedef
      typedef CellSubSetParent<Policy_> BaseClass;

      /// submesh type
      typedef typename Policy_::SubMeshType SubMeshType;
      /// submesh node type
      typedef SubMeshNode<Policy_> SubMeshNodeType;

      /// mesh type of this node
      typedef typename MeshNodePolicy_::MeshType MeshType;
      /// mesh refinery type of this node
      typedef typename MeshNodePolicy_::RefineryType RefineryType;
      /// mesh chart type of this node
      typedef typename MeshNodePolicy_::ChartType MeshChartType;

    protected:
      /**
       * \brief SubMeshNode bin class
       */
      class SubMeshNodeBin
      {
      public:
        SubMeshNodeType* node;
        const MeshChartType* chart;

      public:
        explicit SubMeshNodeBin(
          SubMeshNodeType* node_,
          const MeshChartType* chart_)
           :
          node(node_),
          chart(chart_)
        {
        }
      };

      /// submesh node bin container type
      typedef std::map<Index,SubMeshNodeBin> SubMeshNodeContainer;
      /// submesh node iterator type
      typedef typename SubMeshNodeContainer::iterator SubMeshNodeIterator;
      /// submesh node const-iterator type
      typedef typename SubMeshNodeContainer::const_iterator SubMeshNodeConstIterator;
      /// submesh node reverse-iterator type
      typedef typename SubMeshNodeContainer::reverse_iterator SubMeshNodeReverseIterator;

    protected:
      /// a pointer to the mesh of this node
      MeshType* _mesh;
      /// a pointer to the refinery this node's mesh was created by
      RefineryType* _refinery;
      /// child submesh nodes
      SubMeshNodeContainer _submesh_nodes;

    protected:
      /**
       * \brief Constructor.
       *
       * \param[in] mesh
       * A pointer to the mesh for this node.
       *
       * \param[in] refinery
       * A pointer to the refinery that created \p mesh. May be \c nullptr.
       */
      explicit MeshNode(
        MeshType* mesh,
        RefineryType* refinery = nullptr)
         :
        CellSubSetParent<Policy_>(),
        _mesh(mesh),
        _refinery(refinery)
      {
        CONTEXT(name() + "::MeshNode()");
      }

    public:
      /// virtual destructor
      virtual ~MeshNode()
      {
        CONTEXT(name() + "::~MeshNode()");

        // loop over all submesh nodes in reverse order and delete them
        SubMeshNodeReverseIterator it(_submesh_nodes.rbegin());
        SubMeshNodeReverseIterator jt(_submesh_nodes.rend());
        for(; it != jt; ++it)
        {
          if(it->second.node != nullptr)
          {
            delete it->second.node;
          }
        }

        // delete refinery
        if(_refinery != nullptr)
        {
          delete _refinery;
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
       * \brief Adds a new submesh child node.
       *
       * \param[in] id
       * The id of the child node.
       *
       * \param[in] submesh_node
       * A pointer to the submesh node to be added.
       *
       * \param[in] chart
       * A pointer to the chart that the subnode is to be associated with. May be \c nullptr.
       *
       * \returns
       * \p submesh_node if the insertion was successful, otherwise \c nullptr.
       */
      SubMeshNodeType* add_submesh_node(
        Index id,
        SubMeshNodeType* submesh_node,
        const MeshChartType* chart = nullptr)
      {
        CONTEXT(name() + "::add_submesh_node()");
        if(submesh_node != nullptr)
        {
          if(_submesh_nodes.insert(std::make_pair(id,SubMeshNodeBin(submesh_node, chart))).second)
          {
            return submesh_node;
          }
        }
        return nullptr;
      }

      /**
       * \brief Searches for a submesh node.
       *
       * \param[in] id
       * The id of the node to be found.
       *
       * \returns
       * A pointer to the submesh node associated with \p id or \c nullptr if no such node was found.
       */
      SubMeshNodeType* find_submesh_node(Index id)
      {
        CONTEXT(name() + "::find_submesh_node()");
        SubMeshNodeIterator it(_submesh_nodes.find(id));
        return (it != _submesh_nodes.end()) ? it->second.node : nullptr;
      }

      /** \copydoc find_submesh_node() */
      const SubMeshNodeType* find_submesh_node(Index id) const
      {
        CONTEXT(name() + "::find_submesh_node() [const]");
        SubMeshNodeConstIterator it(_submesh_nodes.find(id));
        return (it != _submesh_nodes.end()) ? it->second.node : nullptr;
      }

      /**
       * \brief Searches for a submesh.
       *
       * \param[in] id
       * The id of the node associated with the submesh to be found.
       *
       * \returns
       * A pointer to the submesh associated with \p id or \c nullptr if no such node was found.
       */
      SubMeshType* find_submesh(Index id)
      {
        CONTEXT(name() + "::find_submesh()");
        SubMeshNodeType* node = find_submesh_node(id);
        return (node != nullptr) ? node->get_mesh() : nullptr;
      }

      /** \copydoc find_submesh() */
      const SubMeshType* find_submesh(Index id) const
      {
        CONTEXT(name() + "::find_submesh() [const]");
        const SubMeshNodeType* node = find_submesh_node(id);
        return (node != nullptr) ? node->get_mesh() : nullptr;
      }

      /**
       * \brief Searches for a submesh chart.
       *
       * \param[in] id
       * The id of the node associated with the chart to be found.
       *
       * \returns
       * A pointer to the chart associated with \p id of \c nullptr if no such node was found or if
       * the corresponding node did not have a chart.
       */
      MeshChartType* find_submesh_chart(Index id)
      {
        CONTEXT(name() + "::find_submesh_chart()");
        SubMeshNodeIterator it(_submesh_nodes.find(id));
        return (it != _submesh_nodes.end()) ? it->second.chart : nullptr;
      }

      /** \copydoc find_submesh_chart() */
      const MeshChartType* find_submesh_chart(Index id) const
      {
        CONTEXT(name() + "::find_submesh_chart() [const]");
        SubMeshNodeConstIterator it(_submesh_nodes.find(id));
        return (it != _submesh_nodes.end()) ? it->second.chart : nullptr;
      }

      /**
       * \brief Adapts this mesh node.
       *
       * This function loops over all submesh nodes and uses their associated charts (if given)
       * to adapt the mesh in this node.
       *
       * \param[in] recursive
       * If set to \c true, all submesh nodes are adapted prior to adapting this node.
       */
      void adapt(bool recursive = true)
      {
        CONTEXT(name() + "::adapt()");

        // loop over all submesh nodes
        SubMeshNodeIterator it(_submesh_nodes.begin());
        SubMeshNodeIterator jt(_submesh_nodes.end());
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
       * The id of the submesh node that is to be used for adaption.
       *
       * \param[in] recursive
       * If set to \c true, the submesh node associated with \p id will be adapted prior to adapting this node.
       *
       * \returns
       * \c true if this node was adapted successfully or \c false if no node is associated with \p id or if
       * the node did not contain any chart.
       */
      bool adapt_by_id(Index id, bool recursive = false)
      {
        CONTEXT(name() + "::adapt_by_id()");

        // try to find the corresponding submesh node
        SubMeshNodeIterator it(_submesh_nodes.find(id));
        if(it == _submesh_nodes.end())
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
        // refine submeshes
        refine_submeshes(refined_node);

        // refine subsets
        refine_subsets(refined_node);
      }

      /**
       * \brief Refines all child submesh nodes of this node.
       *
       * \param[in,out] refined_node
       * A reference to the node generated by refining this node.
       */
      void refine_submeshes(MeshNode& refined_node) const
      {
        SubMeshNodeConstIterator it(_submesh_nodes.begin());
        SubMeshNodeConstIterator jt(_submesh_nodes.end());
        for(; it != jt; ++it)
        {
          refined_node.add_submesh_node(it->first, it->second.node->refine(*_mesh), it->second.chart);
        }
      }

      /**
       * \brief Refines all child cell subset nodes of this node.
       *
       * \param[in,out] refined_node
       * A reference to the node generated by refining this node.
       */
      void refine_subsets(MeshNode& refined_node) const
      {
        typename BaseClass::CellSubSetNodeConstIterator it(this->_subset_nodes.begin());
        typename BaseClass::CellSubSetNodeConstIterator jt(this->_subset_nodes.end());
        for(; it != jt; ++it)
        {
          refined_node.add_subset_node(it->first, it->second->refine(*_mesh));
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
      /// the refinery type of this node
      typedef typename BaseClass::RefineryType RefineryType;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] mesh
       * A pointer to the mesh for this node.
       *
       * \param[in] refinery
       * A pointer to the refinery that created \p mesh. May be \c nullptr.
       */
      explicit RootMeshNode(
        MeshType* mesh,
        RefineryType* refinery = nullptr)
         :
        BaseClass(mesh, refinery)
      {
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
        RefineryType* refinery = new RefineryType(*this->_mesh);

        // refine the mesh
        MeshType* fine_mesh = refinery->refine();

        // create a new root mesh node
        RootMeshNode* fine_node = new RootMeshNode(fine_mesh, refinery);

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
     * \brief Sub-Mesh node class template
     *
     * This class template is used for all mesh nodes of a mesh tree except for the root node.
     *
     * \author Peter Zajac
     */
    template<typename Policy_>
    class SubMeshNode
      : public MeshNode<Policy_, SubMeshNodePolicy<Policy_> >
    {
    public:
      /// base class typedef
      typedef MeshNode<Policy_, SubMeshNodePolicy<Policy_> > BaseClass;

      /// the mesh type of this node
      typedef typename BaseClass::MeshType MeshType;
      /// the refinery type of this node
      typedef typename BaseClass::RefineryType RefineryType;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] mesh
       * A pointer to the mesh for this node.
       *
       * \param[in] refinery
       * A pointer to the refinery that created \p mesh. May be \c nullptr.
       */
      explicit SubMeshNode(
        MeshType* mesh,
        RefineryType* refinery = nullptr)
         :
        BaseClass(mesh, refinery)
      {
      }

      /// virtual destructor
      virtual ~SubMeshNode()
      {
      }

      /**
       * \brief Refines this node and its sub-tree.
       *
       * \param[in] parent_mesh
       * A reference to the parent mesh of this node's submesh.
       *
       * \returns
       * A pointer to a SubMeshNode containing the refined mesh tree.
       */
      template<typename ParentMeshType_>
      SubMeshNode* refine(const ParentMeshType_& parent_mesh) const
      {
        // create a refinery
        RefineryType* refinery = new RefineryType(*this->_mesh);

        // refine the mesh
        MeshType* fine_mesh = refinery->refine(parent_mesh);

        // create a new sub-mesh node
        SubMeshNode* fine_node = new SubMeshNode(fine_mesh, refinery);

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
        return "SubMeshNode<...>";
      }
    }; // class SubMeshNode

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_MESH_NODE_HPP
