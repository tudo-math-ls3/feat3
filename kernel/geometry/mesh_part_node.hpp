#pragma once
#ifndef KERNEL_GEOMETRY_MESH_PART_NODE_HPP
#define KERNEL_GEOMETRY_MESH_PART_NODE_HPP 1

// includes, FEAST
#include <kernel/geometry/mesh_part.hpp>

// includes, STL
#include <map>

namespace FEAST
{
  namespace Geometry
  {
    // forward declarations
    template<typename Policy_>
    class MeshPartNode DOXY({});

    /* ***************************************************************************************** */

    /**
     * \brief MeshPartNode parent base class
     *
     * This class is used by all mesh tree node classes which can have cell subsets as child nodes.
     *
     * \author Peter Zajac
     */
    template<typename Policy_>
    class MeshPartParent
    {
    public:
      /// policy typedef
      typedef Policy_ Policy;
      /// cell subset type
      typedef typename Policy_::MeshPartType MeshPartType;
      /// cell subset node type
      typedef MeshPartNode<Policy_> MeshPartNodeType;

    protected:
      typedef std::map<Index, MeshPartNodeType*> MeshPartNodeContainer;
      typedef typename MeshPartNodeContainer::iterator MeshPartNodeIterator;
      typedef typename MeshPartNodeContainer::const_iterator MeshPartNodeConstIterator;
      typedef typename MeshPartNodeContainer::reverse_iterator MeshPartNodeReverseIterator;

    protected:
      /// child cell subset nodes
      MeshPartNodeContainer _mesh_part_nodes;

      /// protected default constructor
      MeshPartParent()
      {
        CONTEXT(name() + "::MeshPartParent()");
      }

    public:
      /// virtual destructor
      virtual ~MeshPartParent()
      {
        CONTEXT(name() + "::~MeshPartParent()");

        // loop over all child cell subset nodes in reverse order and delete them
        MeshPartNodeReverseIterator it(_mesh_part_nodes.rbegin());
        MeshPartNodeReverseIterator jt(_mesh_part_nodes.rend());
        for(; it != jt; ++it)
        {
          if((it->second) != nullptr)
          {
            delete it->second;
          }
        }
      }

      /**
       * \brief Adds a new mesh part child node.
       *
       * \param[in] id
       * The id of the child node.
       *
       * \param[in] mesh_part_node
       * A pointer to the mesh part node to be added.
       *
       * \returns
       * \p mesh_part_node if the insertion was successful, otherwise \c nullptr.
       */
      MeshPartNodeType* add_mesh_part_node(
        Index id,
        MeshPartNodeType* mesh_part_node)
      {
        CONTEXT(name() + "::add_mesh_part_node()");
        if(mesh_part_node != nullptr)
        {
          if(_mesh_part_nodes.insert(std::make_pair(id,mesh_part_node)).second)
          {
            return mesh_part_node;
          }
        }
        return nullptr;
      }

      /**
       * \brief Searches for a cell subset node.
       *
       * \param[in] id
       * The id of the node to be found.
       *
       * \returns
       * A pointer to the cell subset node associated with \p id or \c nullptr if no such node was found.
       */
      MeshPartNodeType* find_mesh_part_node(Index id)
      {
        CONTEXT(name() + "::find_mesh_part_node()");
        MeshPartNodeIterator it(_mesh_part_nodes.find(id));
        return (it != _mesh_part_nodes.end()) ? it->second : nullptr;
      }

      /** \copydoc find_mesh_part_node() */
      const MeshPartNodeType* find_mesh_part_node(Index id) const
      {
        CONTEXT(name() + "::find_mesh_part_node() [const]");
        MeshPartNodeConstIterator it(_mesh_part_nodes.find(id));
        return (it != _mesh_part_nodes.end()) ? it->second : nullptr;
      }

      /**
       * \brief Searches for a cell subset.
       *
       * \param[in] id
       * The id of the node associated with the cell subset to be found.
       *
       * \returns
       * A pointer to the cell subset associated with \p id or \c nullptr if no such node was found.
       */
      MeshPartType* find_mesh_part(Index id)
      {
        CONTEXT(name() + "::find_mesh_part()");
        MeshPartNodeType* node = find_mesh_part_node(id);
        return (node != nullptr) ? node->get_mesh_part() : nullptr;
      }

      /** \copydoc find_mesh_part() */
      const MeshPartType* find_mesh_part(Index id) const
      {
        CONTEXT(name() + "::find_mesh_part() [const]");
        const MeshPartNodeType* node = find_mesh_part_node(id);
        return (node != nullptr) ? node->get_mesh_part() : nullptr;
      }

      /**
       * \brief Returns the name of the class.
       * \returns
       * The name of the class as a String.
       */
      static String name()
      {
        return "MeshPartParent<...>";
      }
    }; // class MeshPartParent<...>

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
      : public MeshPartParent<Policy_>
    {
    public:
      /// base class typedef
      typedef MeshPartParent<Policy_> BaseClass;
      /// policy type
      typedef Policy_ Policy;
      /// cell subset type
      typedef typename Policy_::MeshPartType MeshPartType;

    protected:
      /// a pointer to the cell subset of this node
      MeshPartType* _mesh_part;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] subset
       * A pointer to the cell subset for this node.
       */
      explicit MeshPartNode(MeshPartType* subset) :
        BaseClass(),
        _mesh_part(subset)
      {
        CONTEXT(name() + "::MeshPartNode()");
      }

      /// virtual destructor
      virtual ~MeshPartNode()
      {
        CONTEXT(name() + "::~MeshPartNode()");
        if(_mesh_part != nullptr)
        {
          delete _mesh_part;
        }
      }

      /**
       * \brief Returns the cell subset of this node.
       * \returns
       * A pointer to the MeshPart contained in this node.
       */
      MeshPartType* get_mesh_part()
      {
        return _mesh_part;
      }

      /** \copydoc get_mesh_part() */
      const MeshPartType* get_mesh_part() const
      {
        return _mesh_part;
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
        StandardRefinery<MeshPartType, ParentType_> refinery(*_mesh_part, parent);

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
          refined_node.add_mesh_part_node(it->first, it->second->refine(*_mesh_part));
        }
      }
    }; // class MeshPartNode<...>

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_MESH_PART_NODE_HPP
