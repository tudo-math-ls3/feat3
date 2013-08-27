#pragma once
#ifndef KERNEL_GEOMETRY_CELL_SUB_SET_NODE_HPP
#define KERNEL_GEOMETRY_CELL_SUB_SET_NODE_HPP 1

// includes, FEAST
#include <kernel/geometry/cell_sub_set.hpp>

// includes, STL
#include <map>

namespace FEAST
{
  namespace Geometry
  {
    // forward declarations
    template<typename Policy_>
    class CellSubSetNode DOXY({});

    /* ***************************************************************************************** */

    /**
     * \brief CellSubSetNode parent base class
     *
     * This class is used by all mesh tree node classes which can have cell subsets as child nodes.
     *
     * \author Peter Zajac
     */
    template<typename Policy_>
    class CellSubSetParent
    {
    public:
      /// policy typedef
      typedef Policy_ Policy;
      /// cell subset type
      typedef typename Policy_::CellSubSetType CellSubSetType;
      /// cell subset node type
      typedef CellSubSetNode<Policy_> CellSubSetNodeType;

    protected:
      typedef std::map<Index, CellSubSetNodeType*> CellSubSetNodeContainer;
      typedef typename CellSubSetNodeContainer::iterator CellSubSetNodeIterator;
      typedef typename CellSubSetNodeContainer::const_iterator CellSubSetNodeConstIterator;
      typedef typename CellSubSetNodeContainer::reverse_iterator CellSubSetNodeReverseIterator;

    protected:
      /// child cell subset nodes
      CellSubSetNodeContainer _subset_nodes;

      /// protected default constructor
      CellSubSetParent()
      {
        CONTEXT(name() + "::CellSubSetParent()");
      }

    public:
      /// virtual destructor
      virtual ~CellSubSetParent()
      {
        CONTEXT(name() + "::~CellSubSetParent()");

        // loop over all child cell subset nodes in reverse order and delete them
        CellSubSetNodeReverseIterator it(_subset_nodes.rbegin());
        CellSubSetNodeReverseIterator jt(_subset_nodes.rend());
        for(; it != jt; ++it)
        {
          if((it->second) != nullptr)
          {
            delete it->second;
          }
        }
      }

      /**
       * \brief Adds a new cell subset child node.
       *
       * \param[in] id
       * The id of the child node.
       *
       * \param[in] subset_node
       * A pointer to the subset node to be added.
       *
       * \returns
       * \p subset_node if the insertion was successful, otherwise \c nullptr.
       */
      CellSubSetNodeType* add_subset_node(
        Index id,
        CellSubSetNodeType* subset_node)
      {
        CONTEXT(name() + "::add_subset_node()");
        if(subset_node != nullptr)
        {
          if(_subset_nodes.insert(std::make_pair(id,subset_node)).second)
          {
            return subset_node;
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
      CellSubSetNodeType* find_subset_node(Index id)
      {
        CONTEXT(name() + "::find_subset_node()");
        CellSubSetNodeIterator it(_subset_nodes.find(id));
        return (it != _subset_nodes.end()) ? it->second : nullptr;
      }

      /** \copydoc find_subset_node() */
      const CellSubSetNodeType* find_subset_node(Index id) const
      {
        CONTEXT(name() + "::find_subset_node() [const]");
        CellSubSetNodeConstIterator it(_subset_nodes.find(id));
        return (it != _subset_nodes.end()) ? it->second : nullptr;
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
      CellSubSetType* find_subset(Index id)
      {
        CONTEXT(name() + "::find_subset()");
        CellSubSetNodeType* node = find_subset_node(id);
        return (node != nullptr) ? node->get_subset() : nullptr;
      }

      /** \copydoc find_subset() */
      const CellSubSetType* find_subset(Index id) const
      {
        CONTEXT(name() + "::find_subset() [const]");
        const CellSubSetNodeType* node = find_subset_node(id);
        return (node != nullptr) ? node->get_subset() : nullptr;
      }

      /**
       * \brief Returns the name of the class.
       * \returns
       * The name of the class as a String.
       */
      static String name()
      {
        return "CellSubSetParent<...>";
      }
    }; // class CellSubSetParent<...>

    /* ***************************************************************************************** */

    /**
     * \brief CellSubSet mesh tree node class template
     *
     * This class template implements a mesh tree node containing a CellSubSet.
     *
     * \author Peter Zajac
     */
    template<typename Policy_>
    class CellSubSetNode
      : public CellSubSetParent<Policy_>
    {
    public:
      /// base class typedef
      typedef CellSubSetParent<Policy_> BaseClass;
      /// policy type
      typedef Policy_ Policy;
      /// cell subset type
      typedef typename Policy_::CellSubSetType CellSubSetType;

    protected:
      /// a pointer to the cell subset of this node
      CellSubSetType* _subset;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] subset
       * A pointer to the cell subset for this node.
       */
      explicit CellSubSetNode(CellSubSetType* subset) :
        BaseClass(),
        _subset(subset)
      {
        CONTEXT(name() + "::CellSubSetNode()");
      }

      /// virtual destructor
      virtual ~CellSubSetNode()
      {
        CONTEXT(name() + "::~CellSubSetNode()");
        if(_subset != nullptr)
        {
          delete _subset;
        }
      }

      /**
       * \brief Returns the cell subset of this node.
       * \returns
       * A pointer to the CellSubSet contained in this node.
       */
      CellSubSetType* get_subset()
      {
        return _subset;
      }

      /** \copydoc get_subset() */
      const CellSubSetType* get_subset() const
      {
        return _subset;
      }

      /**
       * \brief Refines this node and its sub-tree.
       *
       * \param[in] parent
       * A reference to the parent mesh/cell subset of this node's cell subset.
       *
       * \returns
       * A pointer to a CellSubSetNode containing the refined cell subset tree.
       */
      template<typename ParentType_>
      CellSubSetNode* refine(const ParentType_& parent) const
      {
        CONTEXT(name() + "::refine()");

        // create a refinery
        StandardRefinery<CellSubSetType, ParentType_> refinery(*_subset, parent);

        // create a new cell subset node
        CellSubSetNode* fine_node = new CellSubSetNode(new CellSubSetType(refinery));

        // refine our children
        refine_subsets(*fine_node);

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
        return "CellSubSetNode<...>";
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
      void refine_subsets(CellSubSetNode& refined_node) const
      {
        CONTEXT(name() + "::refine_subsets()");

        typename BaseClass::CellSubSetNodeConstIterator it(this->_subset_nodes.begin());
        typename BaseClass::CellSubSetNodeConstIterator jt(this->_subset_nodes.end());
        for(; it != jt; ++it)
        {
          refined_node.add_subset_node(it->first, it->second->refine(*_subset));
        }
      }
    }; // class CellSubSetNode<...>

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CELL_SUB_SET_NODE_HPP
