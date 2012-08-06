#pragma once
#ifndef KERNEL_UTIL_GRAPH_HPP
#define KERNEL_UTIL_GRAPH_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/adjactor.hpp>
#include <kernel/util/assertion.hpp>

// includes, system
#include <iostream>
#include <stdlib.h>

namespace FEAST
{
  // forward declaration
  template<typename Adjactor_>
  class GraphRenderer;

  /**
   * \brief Adjacency Graph implementation
   *
   * \todo detailed description
   *
   * \author Peter Zajac
   */
  class Graph
  {
  public:
    /**
     * \brief ImageIterator for Graph class
     *
     * For the purpose of the Adjactor interface implementation a const Index-pointer is the optimal choice for
     * the image node iterator.
     */
    typedef const Index* ImageIterator;

    /**
     * \brief Share enumeration
     *
     * This enumeration is used by the constructor to determine how the arrays passed to the constructor are to be
     * treated.
     */
    enum Share
    {
      /**
       * \brief Share input arrays
       *
       * This value specifies that the Graph object should point directly to the input arrays passed to the
       * constructor. The Graph object will not delete the arrays upon destruction. It is the responsibility
       * of the caller to ensure that the input arrays remain valid for the lifetime of the Graph object.
       */
      share,

      /**
       * \brief Own input arrays
       *
       * This value specifies that the Graph object should point directly to the input arrays passed to the
       * constructor. In contrast to #share, the arrays will be deleted by the destructor upon destruction of
       * the Graph object.
       */
      own,

      /**
       * \brief Copy input arrays
       *
       * This value specifies that the Graph object should allocate and maintain its own copy of the input
       * arrays passed to the constructor.
       */
      copy
    };

  protected:
    /// total number of domain nodes
    Index _num_nodes_domain;
    /// total number of image nodes
    Index _num_nodes_image;
    /// total number of image node indices
    Index _num_indices_image;

    /**
     * \brief Domain pointer array
     *
     * Dimension: #_num_nodes_domain+1
     */
    Index* _domain_ptr;

    /**
     * \brief Domain end-pointer array
     *
     * Dimension: #_num_nodes_domain
     */
    Index* _domain_end;

    /**
     * \brief Image node index array
     *
     * Dimension: #_num_indices_image
     */
    Index* _image_idx;

    /**
     * \brief Specifies whether the graph's arrays are shared or not
     *
     * This value specifies whether the Graph object will delete the #_domain_ptr, #_domain_end and #_image_idx
     * arrays within the destructor.
     */
    bool _shared;

  public:

    /**
     * \brief Default constructor.
     *
     * This constructor creates a new empty graph, but does not allocate any arrays.
     */
    Graph() :
      _num_nodes_domain(0),
      _num_nodes_image(0),
      _num_indices_image(0),
      _domain_ptr(nullptr),
      _domain_end(nullptr),
      _image_idx(nullptr),
      _shared(false)
    {
    }

    /**
     * \brief Constructor.
     *
     * This constructor creates a new graph and allocates the Graph's arrays to the corresponding lengths.
     *
     * \warning This constructor does not initialise the allocated arrays -- they have to be modified by the user
     * after construction.
     *
     * \param[in] num_nodes_domain
     * The total number of domain nodes for the graph.
     *
     * \param[in] num_nodes_image
     * The total number of image nodes for the graph.
     *
     * \param[in] num_indices_image
     * The total number of image node indices for the graph.
     *
     * \param[in] alloc_domain_end
     * If \c true, then the #_domain_end array will be allocated, otherwise #_domain_end is set to \c nullptr.
     */
    Graph(
      Index num_nodes_domain,
      Index num_nodes_image,
      Index num_indices_image,
      bool alloc_domain_end = false)
       :
      _num_nodes_domain(num_nodes_domain),
      _num_nodes_image(num_nodes_image),
      _num_indices_image(num_indices_image),
      _domain_ptr(nullptr),
      _domain_end(nullptr),
      _image_idx(nullptr),
      _shared(false)
    {
      _domain_ptr = new Index[_num_nodes_domain+1];
      if(alloc_domain_end)
      {
        _domain_end = new Index[_num_nodes_domain];
      }
      _image_idx = new Index[_num_indices_image];
    }

    /**
     * \brief Constructor
     *
     * This constructor creates a new graph based on the arrays given.
     *
     * \param[in] num_nodes_domain
     * The total number of domain nodes for the graph.
     *
     * \param[in] num_nodes_image
     * The total number of image nodes for the graph.
     *
     * \param[in] num_indices_image
     * The total number of image node indices for the graph.
     *
     * \param[in] domain_ptr
     * The domain pointer array for the graph. Must not be \c nullptr.
     *
     * \param[in] domain_end
     * The domain end-pointer array for the graph. May be \c nullptr.
     *
     * \param[in] image_idx
     * The image node index array for the graph. Must not be \c nullptr.
     *
     * \param[in] share
     * Specifies how the input arrays are treated. See Graph::Share for details.
     */
    Graph(
      Index num_nodes_domain,
      Index num_nodes_image,
      Index num_indices_image,
      Index* domain_ptr,
      Index* domain_end,
      Index* image_idx,
      Share share = Graph::copy)
       :
      _num_nodes_domain(num_nodes_domain),
      _num_nodes_image(num_nodes_image),
      _num_indices_image(num_indices_image),
      _domain_ptr(domain_ptr),
      _domain_end(domain_end),
      _image_idx(image_idx),
      _shared(share == Graph::share)
    {
      CONTEXT("Graph::Graph()");
      if(share == Graph::copy)
      {
        // we need to make copies of the arrays
        _domain_ptr = new Index[num_nodes_domain+1];
        for(Index i(0); i <= num_nodes_domain; ++i)
        {
          _domain_ptr[i] = domain_ptr[i];
        }
        if(domain_end != nullptr)
        {
          _domain_end = new Index[num_nodes_domain];
          for(Index i(0); i < num_nodes_domain; ++i)
          {
            _domain_end[i] = domain_end[i];
          }
        }
        _image_idx = new Index[num_indices_image];
        for(Index i(0); i < num_indices_image; ++i)
        {
          _image_idx[i] = image_idx[i];
        }
      }
    }

    /// DTOR
    virtual ~Graph()
    {
      CONTEXT("Graph::destroy()");
      if(!_shared)
      {
        if(_image_idx != nullptr)
          delete [] _image_idx;
        if(_domain_ptr != nullptr)
          delete [] _domain_ptr;
        if(_domain_end != nullptr)
          delete [] _domain_end;
      }
    }

    /**
     * \brief Returns the degree of a domain node.
     *
     * This function returns the degree of a domain node, i.e. the number of image nodes which are adjacent to the
     * specified domain node.
     *
     * \note The degree might be greater than the total number of image nodes if the graph is not injective.
     *
     * \param[in] domain_node
     * The domain node index whose degree is to be returned.
     *
     * \returns
     * The degree of the specified domain node.
     */
    inline Index degree(Index domain_node) const
    {
      CONTEXT("Graph::degree()");
      ASSERT(domain_node < _num_nodes_domain, "Primal node index out of range");
      if(_domain_end != nullptr)
        return _domain_end[domain_node] - _domain_ptr[domain_node];
      else
        return _domain_ptr[domain_node+1] - _domain_ptr[domain_node];
    }

    /**
     * \brief Returns the degree of the graph.
     *
     * This function returns the degree of the graph, i.e. the maximum of the degrees of all domain nodes.
     *
     * \note The degree of the graph might be greater than the total number of image nodes if the graph is not
     * injective.
     *
     * \note This function performs a loop over all domain nodes and therefore has linear runtime. It is
     * recommended to store the result of this function in a variable if the degree is needed multiple times.
     *
     * \returns
     * The degree of the graph.
     */
    inline Index degree() const
    {
      CONTEXT("Graph::degree()");
      Index deg = 0;
      if(_domain_end != nullptr)
      {
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          deg = std::max(deg, _domain_end[i] - _domain_ptr[i]);
        }
      }
      else
      {
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          deg = std::max(deg, _domain_ptr[i+1] - _domain_ptr[i]);
        }
      }
      return deg;
    }

    /**
     * \brief Returns the domain pointer array.
     * \returns The domain pointer array.
     */
    inline Index* get_domain_ptr()
    {
      CONTEXT("Graph::get_domain_ptr()");
      return _domain_ptr;
    }

    /** \copydoc get_domain_ptr() */
    inline const Index* get_domain_ptr() const
    {
      CONTEXT("Graph::get_domain_ptr()");
      return _domain_ptr;
    }

    /**
     * \brief Returns the domain end-pointer array.
     *
     * \warning Please note that a graph does not necessarily have a domain end-pointer array, thus this function
     * might return a \c nullptr even if the graph is not empty.
     *
     * \returns The domain end-pointer array.
     */
    inline Index* get_domain_end()
    {
      CONTEXT("Graph::get_domain_end()");
      return _domain_end;
    }

    /** \copydoc get_domain_end() */
    inline const Index* get_domain_end() const
    {
      CONTEXT("Graph::get_domain_end()");
      return _domain_end;
    }

    /**
     * \brief Returns the image node index array.
     * \returns The image node index array.
     */
    inline Index* get_image_idx()
    {
      CONTEXT("Graph::get_image_idx()");
      return _image_idx;
    }

    /** \copydoc get_image_idx() */
    inline const Index* get_image_idx() const
    {
      CONTEXT("Graph::get_image_idx()");
      return _image_idx;
    }

    /**
     * \brief Specifies whether the graph's arrays are shared or not.
     *
     * \returns \c true, if the graph's arrays are shared, otherwise \c false.
     */
    inline bool is_shared() const
    {
      CONTEXT("Graph::is_shared()");
      return _shared;
    }

    /**
     * \brief Dumps the graph into a stream.
     *
     * \param[in] stream
     * A stream to which the graph is to be dumped into.
     */
    template<typename Stream_>
    void dump(Stream_& stream) const
    {
      CONTEXT("Graph::dump()");
      stream << "number of domain nodes: " << _num_nodes_domain << std::endl;
      stream << "number of image nodes  : " << _num_nodes_image << std::endl;
      stream << "number of image indices: " << _num_indices_image << std::endl;

      const Index* _ptr = _domain_ptr;
      const Index* _end = nullptr;
      if(_domain_end != nullptr)
        _end = _domain_end;
      else
        _end = &_domain_ptr[1];

      if (_num_nodes_domain > 0)
      {
        stream << "node | degree | neighbours: " << std::endl;
        for(Index i(0) ; i < _num_nodes_domain ; ++i)
        {

          stream << i << " | " << (_end[i] - _ptr[i]);
          if (_end[i] - _ptr[i] > 0)
          {
            stream << " | " << _image_idx[_ptr[i]];
            for(Index j(_ptr[i]+1) ; j < _end[i] ; ++j)
            {
              stream << ", " << _image_idx[j];
            }
          }
          stream << std::endl;
        }
      }
    }

    /**
     * \brief Dumps the graph into a string.
     *
     * \returns A string containing the dumped graph.
     */
    inline String dump() const
    {
      CONTEXT("Graph::dump()");
      std::ostringstream oss;
      dump(oss);
      return oss.str();
    }

    /**
     * \brief Renders an adjactor into a graph.
     *
     * \tparam Adjactor_
     * A class implementing the Adjactor interface. This parameter is determined automatically.
     *
     * \param[in] adj
     * A reference to the adjactor which is to be rendered.
     *
     * \param[in] transpose
     * If \c true, then the transposed adjactor will be rendered into the graph.
     *
     * \param[in] injectify
     * If \c true, then the injectified adjactor will be rendered into the graph.
     *
     * \returns
     * A pointer to a Graph containing the rendered adjactor.
     */
    template<typename Adjactor_>
    static Graph* render(
      const Adjactor_& adj,
      bool injectify = true,
      bool transpose = false)
    {
      CONTEXT("Graph::render()");
      return GraphRenderer<Adjactor_>::render(adj, injectify, transpose);
    }

    /**
     * \brief Renders a composition of two Adjactors into a Graph.
     *
     * \todo detailed description
     *
     * \param[in] adj1
     * A reference to the first adjactor of the composition which is to be rendered.
     *
     * \param[in] adj2
     * A reference to the second adjactor of the composition which is to be rendered.
     *
     * \param[in] transpose
     * If \c true, then the transposed adjactor will be rendered into the graph.
     *
     * \param[in] injectify
     * If \c true, then the injectified adjactor will be rendered into the graph.
     *
     * \returns
     * A pointer to a Graph containing the rendered adjactor composition.
     */
    template<
      typename Adjactor1_,
      typename Adjactor2_>
    static Graph* render_composite(
      const Adjactor1_& adj1,
      const Adjactor2_& adj2,
      bool injectify = true,
      bool transpose = false)
    {
      CONTEXT("Graph::render_composite()");
      typedef CompositeAdjactor<Adjactor1_, Adjactor2_> CA12;
      return GraphRenderer<CA12>::render_composite(adj1, adj2, injectify, transpose);
    }


    /* ********************************************************************* */
    /*  A D J U N C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
    /* ********************************************************************* */
    /** \copydoc Adjactor::get_num_nodes_domain() */
    inline Index get_num_nodes_domain() const
    {
      CONTEXT("Graph::get_num_nodes_domain()");
      return _num_nodes_domain;
    }

    /** \copydoc Adjactor::get_num_nodes_image() */
    inline Index get_num_nodes_image() const
    {
      CONTEXT("Graph::get_num_nodes_image()");
      return _num_nodes_image;
    }

    /** \copydoc Adjactor::image_begin() */
    inline ImageIterator image_begin(Index domain_node) const
    {
      CONTEXT("Graph::image_begin()");
      ASSERT(domain_node < _num_nodes_domain, "Primal node index out of range");
      return &_image_idx[_domain_ptr[domain_node]];
    }

    /** \copydoc Adjactor::image_end() */
    inline ImageIterator image_end(Index domain_node) const
    {
      CONTEXT("Graph::image_end()");
      ASSERT(domain_node < _num_nodes_domain, "Primal node index out of range");

      if(_domain_end != nullptr)
        return &_image_idx[_domain_end[domain_node]];
      else
        return &_image_idx[_domain_ptr[domain_node+1]];
    }
  }; // class Graph

  /**
   * \brief Graph renderer implementation
   *
   * \author Peter Zajac
   */
  template<typename Adjactor_>
  class GraphRenderer
  {
  protected:
    /// \cond internal
    static Index _aux_inj(const Adjactor_& adj, Index i, Index* idx)
    {
      CONTEXT("GraphRenderer::_aux_inj()");

      Index num_idx = 0;
      typename Adjactor_::ImageIterator cur(adj.image_begin(i));
      typename Adjactor_::ImageIterator end(adj.image_end(i));
      for(; cur != end; ++cur)
      {
        Index jdx = *cur;
        bool found = false;
        for(Index k(0); k < num_idx; ++k)
        {
          if(idx[k] == jdx)
          {
            found = true;
            break;
          }
        }
        if(!found)
        {
          idx[num_idx] = jdx;
          ++num_idx;
        }
      }
      return num_idx;
    }

    /// renders adjactor
    static Graph* _render(const Adjactor_& adj)
    {
      CONTEXT("GraphRenderer::_render()");

      typedef typename Adjactor_::ImageIterator AImIt;

      // get counts
      Index num_nodes_domain = adj.get_num_nodes_domain();
      Index num_nodes_image = adj.get_num_nodes_image();
      Index num_indices_image = 0;

      // allocate pointer array
      Index* domain_ptr = new Index[num_nodes_domain+1];

      // count number of adjacencies and build pointer array
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        domain_ptr[i] = num_indices_image;
        AImIt cur(adj.image_begin(i));
        AImIt end(adj.image_end(i));
        for(; cur != end; ++cur)
        {
          ++num_indices_image;
        }
      }
      domain_ptr[num_nodes_domain] = num_indices_image;

      // allocate and build index array
      Index* image_idx = new Index[num_indices_image];
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        Index* idx = &image_idx[domain_ptr[i]];
        AImIt cur(adj.image_begin(i));
        AImIt end(adj.image_end(i));
        for(; cur != end; ++cur, ++idx)
        {
          *idx = *cur;
        }
      }

      // create the Graph
      return new Graph(num_nodes_domain, num_nodes_image, num_indices_image,
        domain_ptr, nullptr, image_idx, Graph::own);
    }

    /// renders injectified adjactor
    static Graph* _render_inj(const Adjactor_& adj)
    {
      CONTEXT("GraphRenderer::_render_inj()");

      // get counts
      Index num_nodes_domain = adj.get_num_nodes_domain();
      Index num_nodes_image = adj.get_num_nodes_image();
      Index num_indices_image = 0;

      // allocate pointer array
      Index* domain_ptr = new Index[num_nodes_domain+1];

      // allocate auxiliary index array
      Index* aux = new Index[num_nodes_image];

      // count number of adjacencies and build pointer array
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        domain_ptr[i] = num_indices_image;
        num_indices_image += _aux_inj(adj, i, aux);
      }
      domain_ptr[num_nodes_domain] = num_indices_image;

      // delete auxiliary array
      delete [] aux;

      // allocate and build index array
      Index* image_idx = new Index[num_indices_image];
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        _aux_inj(adj, i, &image_idx[domain_ptr[i]]);
      }

      // create graph
      return new Graph(num_nodes_domain, num_nodes_image, num_indices_image,
        domain_ptr, nullptr, image_idx, Graph::own);
    }

    /// renders transposed adjactor
    static Graph* _render_trans(const Adjactor_& adj)
    {
      CONTEXT("GraphRenderer::_render_trans()");

      typedef typename Adjactor_::ImageIterator AImIt;

      // get counts
      Index num_nodes_domain = adj.get_num_nodes_image();
      Index num_nodes_image = adj.get_num_nodes_domain();
      Index num_indices_image = 0;

      // allocate and format pointer array
      Index *domain_ptr = new Index[num_nodes_domain+1];
      for(Index i(0); i <= num_nodes_domain; ++i)
      {
        domain_ptr[i] = 0;
      }

      // count number of adjacencies
      for(Index j(0); j < num_nodes_image; ++j)
      {
        AImIt cur(adj.image_begin(j));
        AImIt end(adj.image_end(j));
        for(; cur != end; ++cur)
        {
          ++domain_ptr[(*cur)+1];
        }
      }

      // build pointer array
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        domain_ptr[i+1] += domain_ptr[i];
      }
      num_indices_image = domain_ptr[num_nodes_domain];

      // allocate and build index array
      Index* image_idx = new Index[num_indices_image];
      Index** image_ptr = new Index*[num_nodes_domain];
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        image_ptr[i] = &image_idx[domain_ptr[i]];
      }

      for(Index j(0); j < num_nodes_image; ++j)
      {
        AImIt cur(adj.image_begin(j));
        AImIt end(adj.image_end(j));
        for(; cur != end; ++cur)
        {
          Index*& idx = image_ptr[*cur];
          *idx = j;
          ++idx;
        }
      }

      delete [] image_ptr;

      // create graph
      return new Graph(num_nodes_domain, num_nodes_image, num_indices_image,
        domain_ptr, nullptr, image_idx, Graph::own);
    }

    /// renders transposed injectified adjactor
    static Graph* _render_trans_inj(const Adjactor_& adj)
    {
      CONTEXT("GraphRenderer::_render_trans_inj()");

      // get counts
      Index num_nodes_domain = adj.get_num_nodes_image();
      Index num_nodes_image = adj.get_num_nodes_domain();
      Index num_indices_image = 0;

      // allocate pointer array
      Index* domain_ptr = new Index[num_nodes_domain+1];
      for(Index i(0); i <= num_nodes_domain; ++i)
      {
        domain_ptr[i] = 0;
      }

      // allocate auxiliary index array
      Index* aux = new Index[num_nodes_domain];

      // loop over all image nodes
      for(Index j(0); j < num_nodes_image; ++j)
      {
        Index num_aux = _aux_inj(adj, j, aux);
        for(Index k(0); k < num_aux; ++k)
        {
          ++domain_ptr[aux[k]+1];
        }
        num_indices_image += num_aux;
      }

      Index* image_idx = new Index[num_indices_image];
      Index** image_ptr = new Index*[num_nodes_domain];

      // build pointer arrays
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        domain_ptr[i+1] += domain_ptr[i];
        image_ptr[i] = &image_idx[domain_ptr[i]];
      }

      // build image index array
      for(Index j(0); j < num_nodes_image; ++j)
      {
        Index num_aux = _aux_inj(adj, j, aux);
        for(Index k(0); k < num_aux; ++k)
        {
          Index*& idx = image_ptr[aux[k]];
          *idx = j;
          ++idx;
        }
      }

      // delete auxiliary arrays
      delete [] image_ptr;
      delete [] aux;

      // create graph
      return new Graph(num_nodes_domain, num_nodes_image, num_indices_image,
        domain_ptr, nullptr, image_idx, Graph::own);
    }

    /// \endcond
  public:
    /** \copydoc Graph::render() */
    static Graph* render(
      const Adjactor_& adj,
      bool injectify = true,
      bool transpose = false)
    {
      CONTEXT("GraphRenderer::render()");
      if(transpose)
      {
        if(injectify)
        {
          return _render_trans_inj(adj);
        }
        else
        {
          return _render_trans(adj);
        }
      }
      else
      {
        if(injectify)
        {
          return _render_inj(adj);
        }
        else
        {
          return _render(adj);
        }
      }
    }
  }; // class GraphRenderer

  /**
   * \brief Partial specialisation of GraphRenderer for CompositeAdjactor
   *
   * \author Peter Zajac
   */
  template<
    typename Adjactor1_,
    typename Adjactor2_>
  class GraphRenderer< CompositeAdjactor<Adjactor1_, Adjactor2_> >
  {
  protected:
    /// \cond internal
    static Index _aux_inj(const Adjactor1_& adj1, const Adjactor2_& adj2, Index i, Index* idx)
    {
      CONTEXT("GraphRenderer<CompositeAdjactor>::_aux_inj()");

      Index num_idx = 0;
      typename Adjactor1_::ImageIterator cur1(adj1.image_begin(i));
      typename Adjactor1_::ImageIterator end1(adj1.image_end(i));
      for(; cur1 != end1; ++cur1)
      {
        typename Adjactor2_::ImageIterator cur2(adj2.image_begin(*cur1));
        typename Adjactor2_::ImageIterator end2(adj2.image_end(*cur1));

        // loop over all image node indices
        for(; cur2 != end2; ++cur2)
        {
          Index jdx = *cur2;

          // check if we already have that image node index
          bool found = false;
          for(Index k(0); k < num_idx; ++k)
          {
            if(idx[k] == jdx)
            {
              found = true;
              break;
            }
          }
          // add the index if we don't have it already
          if(!found)
          {
            idx[num_idx] = jdx;
            ++num_idx;
          }
        }
      }
      return num_idx;
    }

    /// renders adjactor composition
    static Graph* _render(
      const Adjactor1_& adj1,
      const Adjactor2_& adj2)
    {
      CONTEXT("GraphRenderer<CompositeAdjactor>::_render()");

      typedef typename Adjactor1_::ImageIterator AImIt1;
      typedef typename Adjactor2_::ImageIterator AImIt2;

      // get counts
      Index num_nodes_domain = adj1.get_num_nodes_domain();
      Index num_nodes_image = adj2.get_num_nodes_image();
      Index num_indices_image = 0;

      // allocate pointer array
      Index* domain_ptr = new Index[num_nodes_domain+1];

      // count number of adjacencies and build pointer array
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        domain_ptr[i] = num_indices_image;
        AImIt1 cur1(adj1.image_begin(i));
        AImIt1 end1(adj1.image_end(i));
        for(; cur1 != end1; ++cur1)
        {
          AImIt2 cur2(adj2.image_begin(*cur1));
          AImIt2 end2(adj2.image_end(*cur1));
          for(; cur2 != end2; ++cur2)
          {
            ++num_indices_image;
          }
        }
      }
      domain_ptr[num_nodes_domain] = num_indices_image;

      // allocate and build index array
      Index* image_idx = new Index[num_indices_image];
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        Index* idx = &image_idx[domain_ptr[i]];
        AImIt1 cur1(adj1.image_begin(i));
        AImIt1 end1(adj1.image_end(i));
        for(; cur1 != end1; ++cur1)
        {
          AImIt2 cur2(adj2.image_begin(*cur1));
          AImIt2 end2(adj2.image_end(*cur1));
          for(; cur2 != end2; ++cur2, ++idx)
          {
            *idx = *cur2;
          }
        }
      }

      // create graph
      return new Graph(num_nodes_domain, num_nodes_image, num_indices_image,
        domain_ptr, nullptr, image_idx, Graph::own);
    }

    /// renders injectified adjactor composition
    static Graph* _render_inj(
      const Adjactor1_& adj1,
      const Adjactor2_& adj2)
    {
      CONTEXT("GraphRenderer<CompositeAdjactor>::_render_inj()");

      typedef typename Adjactor1_::ImageIterator AImIt1;
      typedef typename Adjactor2_::ImageIterator AImIt2;

      // get counts
      Index num_nodes_domain = adj1.get_num_nodes_domain();
      Index num_nodes_image = adj2.get_num_nodes_image();
      Index num_indices_image = 0;

      // allocate pointer array
      Index* domain_ptr = new Index[num_nodes_domain+1];

      // allocate auxiliary array
      Index* aux = new Index[num_nodes_image];

      // count number of adjacencies and build pointer array
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        domain_ptr[i] = num_indices_image;
        num_indices_image += _aux_inj(adj1, adj2, i, aux);
      }
      domain_ptr[num_nodes_domain] = num_indices_image;

      // delete auxiliary array
      delete [] aux;

      // allocate and build index array
      Index* image_idx = new Index[num_indices_image];
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        _aux_inj(adj1, adj2, i, &image_idx[domain_ptr[i]]);
      }

      // create graph
      return new Graph(num_nodes_domain, num_nodes_image, num_indices_image,
        domain_ptr, nullptr, image_idx, Graph::own);
    }

    /// renders transposed adjactor composition
    static Graph* _render_trans(
      const Adjactor1_& adj1,
      const Adjactor2_& adj2)
    {
      CONTEXT("GraphRenderer<CompositeAdjactor>::_render_trans()");

      typedef typename Adjactor1_::ImageIterator AImIt1;
      typedef typename Adjactor2_::ImageIterator AImIt2;

      // get counts
      Index num_nodes_domain = adj2.get_num_nodes_image();
      Index num_nodes_image = adj1.get_num_nodes_domain();
      Index num_indices_image = 0;

      // allocate and format pointer array
      Index *domain_ptr = new Index[num_nodes_domain+1];
      for(Index i(0); i <= num_nodes_domain; ++i)
      {
        domain_ptr[i] = 0;
      }

      // count number of adjacencies
      for(Index j(0); j < num_nodes_image; ++j)
      {
        AImIt1 cur1(adj1.image_begin(j));
        AImIt1 end1(adj1.image_end(j));
        for(; cur1 != end1; ++cur1)
        {
          AImIt2 cur2(adj2.image_begin(*cur1));
          AImIt2 end2(adj2.image_end(*cur1));
          for(; cur2 != end2; ++cur2)
          {
            ++domain_ptr[(*cur2)+1];
          }
        }
      }

      // build pointer array
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        domain_ptr[i+1] += domain_ptr[i];
      }
      num_indices_image = domain_ptr[num_nodes_domain];

      // allocate and build index array
      Index* image_idx = new Index[num_indices_image];
      Index** image_ptr = new Index*[num_nodes_domain];
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        image_ptr[i] = &image_idx[domain_ptr[i]];
      }

      for(Index j(0); j < num_nodes_image; ++j)
      {
        AImIt1 cur1(adj1.image_begin(j));
        AImIt1 end1(adj1.image_end(j));
        for(; cur1 != end1; ++cur1)
        {
          AImIt2 cur2(adj2.image_begin(*cur1));
          AImIt2 end2(adj2.image_end(*cur1));
          for(; cur2 != end2; ++cur2)
          {
            Index*& idx = image_ptr[*cur2];
            *idx = j;
            ++idx;
          }
        }
      }

      delete [] image_ptr;

      // create graph
      return new Graph(num_nodes_domain, num_nodes_image, num_indices_image,
        domain_ptr, nullptr, image_idx, Graph::own);
    }

    /// renders transposed injectified adjactor composition
    static Graph* _render_trans_inj(
      const Adjactor1_& adj1,
      const Adjactor2_& adj2)
    {
      CONTEXT("GraphRenderer<CompositeAdjactor>::_render_trans_inj()");

      // get counts
      Index num_nodes_domain = adj2.get_num_nodes_image();
      Index num_nodes_image = adj1.get_num_nodes_domain();
      Index num_indices_image = 0;

      // allocate pointer array
      Index* domain_ptr = new Index[num_nodes_domain+1];
      for(Index i(0); i <= num_nodes_domain; ++i)
      {
        domain_ptr[i] = 0;
      }

      // allocate auxiliary index array
      Index* aux = new Index[num_nodes_domain];

      // loop over all image nodes
      for(Index j(0); j < num_nodes_image; ++j)
      {
        Index num_aux = _aux_inj(adj1, adj2, j, aux);
        for(Index k(0); k < num_aux; ++k)
        {
          ++domain_ptr[aux[k]+1];
        }
        num_indices_image += num_aux;
      }

      Index* image_idx = new Index[num_indices_image];
      Index** image_ptr = new Index*[num_nodes_domain];

      // build pointer arrays
      for(Index i(0); i < num_nodes_domain; ++i)
      {
        domain_ptr[i+1] += domain_ptr[i];
        image_ptr[i] = &image_idx[domain_ptr[i]];
      }

      // build image index array
      for(Index j(0); j < num_nodes_image; ++j)
      {
        Index num_aux = _aux_inj(adj1, adj2, j, aux);
        for(Index k(0); k < num_aux; ++k)
        {
          Index*& idx = image_ptr[aux[k]];
          *idx = j;
          ++idx;
        }
      }

      // delete auxiliary arrays
      delete [] image_ptr;
      delete [] aux;

      // create graph
      return new Graph(num_nodes_domain, num_nodes_image, num_indices_image,
        domain_ptr, nullptr, image_idx, Graph::own);
    }
    /// \endcond

  public:
    /** \copydoc Graph::render_composite() */
    static Graph* render_composite(
      const Adjactor1_& adj1,
      const Adjactor2_& adj2,
      bool injectify = true,
      bool transpose = false)
    {
      CONTEXT("GraphRenderer<CompositeAdjactor>::render()");
      if(transpose)
      {
        if(injectify)
        {
          return _render_trans_inj(adj1, adj2);
        }
        else
        {
          return _render_trans(adj1, adj2);
        }
      }
      else
      {
        if(injectify)
        {
          return _render_inj(adj1, adj2);
        }
        else
        {
          return _render(adj1, adj2);
        }
      }
    }

    /**
     * \brief Renders a CompositeAdjactor into a Graph.
     *
     * \todo detailed description
     *
     * \param[in] adj
     * A reference to the CompositeAdjactor which is to be rendered.
     *
     * \param[in] transpose
     * If \c true, then the transposed adjactor will be rendered into the graph.
     *
     * \param[in] injectify
     * If \c true, then the injectified adjactor will be rendered into the graph.
     *
     * \returns
     * A pointer to a Graph containing the rendered CompositeAdjactor.
     */
    static Graph* render(
      const CompositeAdjactor<Adjactor1_,Adjactor2_>& adj,
      bool injectify = true,
      bool transpose = false)
    {
      CONTEXT("GraphRenderer<CompositeAdjactor>::render()");
      return render_composite(adj.get_adjactor1(), adj.get_adjactor2(), injectify, transpose);
    }
  }; // class GraphRenderer< CompositeAdjactor<...> >
} // namespace FEAST

#endif // KERNEL_UTIL_GRAPH_HPP
