#pragma once
#ifndef KERNEL_UTIL_DYNAMIC_GRAPH_HPP
#define KERNEL_UTIL_DYNAMIC_GRAPH_HPP 1

// includes, FEAST
#include <kernel/util/adjactor.hpp>

// includes, system
#include <vector>
#include <set>

namespace FEAST
{
  /**
   * \brief Dynamic Adjacency Graph implementation
   *
   * This class implements a dynamic counterpart of the Graph class.
   *
   * \author Peter Zajac
   */
  class DynamicGraph
  {
  protected:
    /// \cond internal
    typedef std::set<Index> IdxSet;
    typedef std::vector<IdxSet> IdxSetVec;
    /// \endcond

  public:
    /// ImageIterator for Adjactor interface implementation
    typedef IdxSet::const_iterator ImageIterator;

    /**
     * \brief Render type enumeration
     *
     * \see Graph::RenderType
     */
    enum RenderType
    {
      /**
       * \brief Render-As-Is mode
       *
       * \see Graph::rt_as_is
       */
      rt_as_is = 0,

      /**
       * \brief Render-Transpose mode
       *
       * \see Graph::rt_transpose
       */
      rt_transpose = 2
    };

  protected:
    /// total number of domain nodes
    Index _num_nodes_domain;
    /// total number of image nodes
    Index _num_nodes_image;
    /// index-set-vector
    IdxSetVec _indices;

  public:
    /**
     * \brief Constructor
     *
     * This constructor creates an empty dynamic adjacency graph.
     *
     * \param[in] num_nodes_domain
     * The total number of domain nodes for the graph.
     *
     * \param[in] num_nodes_image
     * The total number of image nodes for the graph.
     */
    explicit DynamicGraph(
      Index num_nodes_domain,
      Index num_nodes_image)
       :
      _num_nodes_domain(num_nodes_domain),
      _num_nodes_image(num_nodes_image),
      _indices(num_nodes_domain)
    {
      CONTEXT("DynamicGraph::DynamicGraph()");
    }

    /**
     * \brief Render constructor
     *
     * This constructor renders an object implementing the Adjactor interface into a graph.
     *
     * \param[in] render_type
     * The render type. See #RenderType for details.
     *
     * \param[in] adjactor
     * A reference to the adjactor that is to be rendered.
     */
    template<typename Adjactor_>
    explicit DynamicGraph(
      DynamicGraph::RenderType render_type,
      const Adjactor_& adjactor)
       :
      _num_nodes_domain(0),
      _num_nodes_image(0),
      _indices()
    {
      CONTEXT("DynamicGraph::DynamicGraph() [render]");
      switch(render_type)
      {
      case rt_as_is:
        _render_as_is(adjactor);
        break;

      case rt_transpose:
        _render_transpose(adjactor);
        break;
      }
    }

    /**
     * \brief Composite-Render constructor
     *
     * This constructor renders a composition of two objects implementing the Adjactor interface into a graph.
     *
     * \param[in] render_type
     * The render type. See #RenderType for details.
     *
     * \param[in] adjactor1
     * A reference to the first adjactor in the composition that is to be rendered.
     *
     * \param[in] adjactor2
     * A reference to the second adjactor in the composition that is to be rendered.
     */
    template<
      typename Adjactor1_,
      typename Adjactor2_>
    explicit DynamicGraph(
      DynamicGraph::RenderType render_type,
      const Adjactor1_& adjactor1,
      const Adjactor2_& adjactor2)
       :
      _num_nodes_domain(0),
      _num_nodes_image(0),
      _indices()
    {
      CONTEXT("DynamicGraph::DynamicGraph() [composite render]");
      switch(render_type)
      {
      case rt_as_is:
        _render_as_is(adjactor1, adjactor2);
        break;

      case rt_transpose:
        _render_transpose(adjactor1, adjactor2);
        break;
      }
    }

    /// virtual destructor
    virtual ~DynamicGraph()
    {
      CONTEXT("DynamicGraph::~DynamicGraph()");
    }

    /**
     * \brief Returns the degree of a domain node.
     *
     * This function returns the degree of a domain node, i.e. the number of image nodes which are adjacent to the
     * specified domain node.
     *
     * \param[in] domain_node
     * The domain node index whose degree is to be returned.
     *
     * \returns
     * The degree of the specified domain node.
     */
    Index degree(Index domain_node) const
    {
      CONTEXT("DynamicGraph::degree()");
      ASSERT(domain_node < _num_nodes_domain, "Domain node index out of range");
      return Index(_indices[domain_node].size());
    }

    /**
     * \brief Returns the degree of the graph.
     *
     * This function returns the degree of the graph, i.e. the maximum of the degrees of all domain nodes.
     *
     * \attention This function performs a loop over all domain nodes and therefore has linear runtime. It is
     * recommended to store the result of this function in a variable if the degree is needed multiple times.
     *
     * \returns
     * The degree of the graph.
     */
    Index degree() const
    {
      CONTEXT("DynamicGraph::degree()");
      Index deg(0);
      for(Index i(0); i < _num_nodes_domain; ++i)
      {
        deg = std::max(deg, Index(_indices[i].size()));
      }
      return deg;
    }

    /**
     * \brief Clears the graph.
     *
     * This function erases all adjacency entries in the graph.
     */
    void clear()
    {
      CONTEXT("DynamicGraph::clear()");
      for(Index i(0); i < _num_nodes_domain; ++i)
      {
        _indices[i].clear();
      }
    }

    /**
     * \brief Checks whether an adjacency entry exists.
     *
     * \param[in] domain_node
     * The index of the domain node of the adjacency entry.
     *
     * \param[in] image_node
     * The index of the image node of the adjacency entry.
     *
     * \returns
     * \c true, if the adjacency exists, otherwise \c false.
     */
    bool exists(Index domain_node, Index image_node) const
    {
      CONTEXT("DynamicGraph::exists()");
      ASSERT(domain_node < _num_nodes_domain, "Domain node index out of range");
      ASSERT(image_node < _num_nodes_image, "Image node index out of range");
      return _indices[domain_node].find(image_node) != _indices[domain_node].end();
    }

    /**
     * \brief Inserts a new adjacency entry to the graph.
     *
     * \param[in] domain_node
     * The index of the domain node of the adjacency entry.
     *
     * \param[in] image_node
     * The index of the image node of the adjacency entry.
     *
     * \returns
     * \c true, if the adjacency was inserted or \c false, if the adjacency already existed.
     */
    bool insert(Index domain_node, Index image_node)
    {
      CONTEXT("DynamicGraph::insert()");
      ASSERT(domain_node < _num_nodes_domain, "Domain node index out of range");
      ASSERT(image_node < _num_nodes_image, "Image node index out of range");
      return _indices[domain_node].insert(image_node).second;
    }

    /**
     * \brief Erases an adjacency entry from the graph.
     *
     * \param[in] domain_node
     * The index of the domain node of the adjacency entry.
     *
     * \param[in] image_node
     * The index of the image node of the adjacency entry.
     *
     * \returns
     * \c true, if the adjacency was erased or \c false, if the adjacency did not exist.
     */
    bool erase(Index domain_node, Index image_node)
    {
      CONTEXT("DynamicGraph::erase()");
      ASSERT(domain_node < _num_nodes_domain, "Domain node index out of range");
      ASSERT(image_node < _num_nodes_image, "Image node index out of range");
      // note: std::set.erase() returns the number of deleted elements
      return _indices[domain_node].erase(image_node) > 0;
    }

    /* *************************************************** */
    /*  R E N D E R   F U N C T I O N   T E M P L A T E S  */
    /* *************************************************** */
  private:
    /// \cond internal
    /// renders adjactor
    template<typename Adjactor_>
    void _render_as_is(const Adjactor_& adj)
    {
      CONTEXT("DynamicGraph::_render_as_is()");

      typedef typename Adjactor_::ImageIterator AImIt;

      // set counts
      _num_nodes_domain = adj.get_num_nodes_domain();
      _num_nodes_image = adj.get_num_nodes_image();
      _indices.resize(_num_nodes_domain);

      // loop over all domain nodes
      for(Index i(0); i < _num_nodes_domain; ++i)
      {
        AImIt cur(adj.image_begin(i));
        AImIt end(adj.image_end(i));
        for(; cur != end; ++cur)
        {
          insert(i, *cur);
        }
      }
    }

    /// renders transposed adjactor
    template<typename Adjactor_>
    void _render_transpose(const Adjactor_& adj)
    {
      CONTEXT("DynamicGraph::_render_transpose()");

      typedef typename Adjactor_::ImageIterator AImIt;

      // set counts
      _num_nodes_domain = adj.get_num_nodes_image();
      _num_nodes_image = adj.get_num_nodes_domain();
      _indices.resize(_num_nodes_domain);

      // loop over all image nodes
      for(Index i(0); i < _num_nodes_image; ++i)
      {
        AImIt cur(adj.image_begin(i));
        AImIt end(adj.image_end(i));
        for(; cur != end; ++cur)
        {
          insert(*cur, i);
        }
      }
    }

    /// renders adjactor composition
    template<
      typename Adjactor1_,
      typename Adjactor2_>
    void _render_as_is(
      const Adjactor1_& adj1,
      const Adjactor2_& adj2)
    {
      CONTEXT("DynamicGraph::_render_as_is() [composite]");

      // validate adjactor dimensions
      ASSERT(adj1.get_num_nodes_image() <= adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

      typedef typename Adjactor1_::ImageIterator AImIt1;
      typedef typename Adjactor2_::ImageIterator AImIt2;

      // set counts
      _num_nodes_domain = adj1.get_num_nodes_domain();
      _num_nodes_image = adj2.get_num_nodes_image();
      _indices.resize(_num_nodes_domain);

      // loop over all domain nodes
      for(Index i(0); i < _num_nodes_domain; ++i)
      {
        AImIt1 cur1(adj1.image_begin(i));
        AImIt1 end1(adj1.image_end(i));
        for(; cur != end; ++cur)
        {
          AImIt2 cur2(adj2.image_begin(*cur1));
          AImIt2 end2(adj2.image_end(*cur1));
          for(; cur2 != end2; ++cur2)
          {
            insert(i, *cur2);
          }
        }
      }
    }

    template<
      typename Adjactor1_,
      typename Adjactor2_>
    void _render_transpose(
      const Adjactor1_& adj1,
      const Adjactor2_& adj2)
    {
      CONTEXT("DynamicGraph::_render_transpose() [composite]");

      // validate adjactor dimensions
      ASSERT(adj1.get_num_nodes_image() <= adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

      typedef typename Adjactor1_::ImageIterator AImIt1;
      typedef typename Adjactor2_::ImageIterator AImIt2;

      // set counts
      _num_nodes_domain = adj2.get_num_nodes_image();
      _num_nodes_image = adj1.get_num_nodes_domain();
      _indices.resize(_num_nodes_domain);

      // loop over all domain nodes
      for(Index i(0); i < _num_nodes_image; ++i)
      {
        AImIt1 cur1(adj1.image_begin(i));
        AImIt1 end1(adj1.image_end(i));
        for(; cur != end; ++cur)
        {
          AImIt2 cur2(adj2.image_begin(*cur1));
          AImIt2 end2(adj2.image_end(*cur1));
          for(; cur2 != end2; ++cur2)
          {
            insert(*cur2, i);
          }
        }
      }
    }
    /// \endcond

    /* ******************************************************************* */
    /*  A D J A C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
    /* ******************************************************************* */
  public:
    /** \copydoc Adjactor::get_num_nodes_domain() */
    inline Index get_num_nodes_domain() const
    {
      CONTEXT("DynamicGraph::get_num_nodes_domain()");
      return _num_nodes_domain;
    }

    /** \copydoc Adjactor::get_num_nodes_image() */
    inline Index get_num_nodes_image() const
    {
      CONTEXT("DynamicGraph::get_num_nodes_image()");
      return _num_nodes_image;
    }

    /** \copydoc Adjactor::image_begin() */
    inline ImageIterator image_begin(Index domain_node) const
    {
      CONTEXT("DynamicGraph::image_begin()");
      ASSERT(domain_node < _num_nodes_domain, "Domain node index out of range");
      return _indices[domain_node].begin();
    }

    /** \copydoc Adjactor::image_end() */
    inline ImageIterator image_end(Index domain_node) const
    {
      CONTEXT("DynamicGraph::image_end()");
      ASSERT(domain_node < _num_nodes_domain, "Domain node index out of range");
      return _indices[domain_node].end();
    }
  }; // class DynamicGraph
} // namespace FEAST

#endif // KERNEL_UTIL_DYNAMIC_GRAPH_HPP
