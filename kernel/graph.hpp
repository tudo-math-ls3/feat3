/* GENERAL_REMARK_BY_HILMAR:
 * This class is really just a rudimentary graph implementation which I needed for testing MPI stuff.
 * As far as I remember Peter wants to take a closer look and merge it with his more sophisticated graph
 * implementation. (Feel free to throw my implementation away completely, Peter!) Also see my comments
 * below.
 *
 * HILMAR WON'T TOUCH THIS FILE ANYMORE! Please remove this comment-block as soon as possible... :-)
 */
#pragma once
#ifndef KERNEL_GRAPH_HPP
#define KERNEL_GRAPH_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/adjunctor.hpp>
#include <kernel/util/assertion.hpp>

// includes, system
#include <iostream>
#include <stdlib.h>

namespace FEAST
{
#ifdef OLD_GRAPH

  /*
  * \brief class providing a graph data structure for defining connectivity of subdomains / matrix patches / processes
  *
  * This data structure is very similar to that described in the MPI-2.2 standard (p. 250ff) for creating process
  * topologies (for the difference, see description of member #_index).
  *
  * In the case this data structure is used for storing the MPI communication graph, it is a \em global representation,
  * i.e. a process storing this graph knows the \em complete communication graph. (For a distributed representation, see
  * data structure GraphDistributed.) To construct a corresponding MPI communicator, one can use the function
  *   int MPI_Graph_create(MPI_Comm comm_old, int num_nodes, int *_index, int *_neighbours, int reorder,
  *                        MPI_Comm *comm_graph)
  * (see MPI-2.2 standard, p. 250). When passing the array #_index to this routine, remember to omit the first entry
  * (see description of the array #_index).
  *
  * The data structure can also be used to create distributed graph data structures via the MPI routine
  * MPI_Dist_graph_create(...), (see MPI-2.2 standard, example 7.3 on page 256).
  *
  * COMMENT_HILMAR: This is only a very rough first version, which will be surely adapted to our needs...
  *
  * COMMENT_HILMAR:
  * This graph structure is the most general one. When it comes to MPI communication, we surely have to
  * distinguish face, edge and vertex neighbours, i.e., making this distinction within an enhanced graph class
  * or using more than one graph.
  *
  * \author Hilmar Wobker
  */
  class Graph
  {

  private:

    /* *****************
    * member variables *
    *******************/
    /// number of nodes in the graph, which are numbered from 0 to num_nodes-1
    Index const _num_nodes;

    /**
    * \brief access information for the array #_neighbours
    *
    * The i-th entry of this array stores the total number of neighbours of the graph nodes 0, ..., i-1, where
    * index[0] = 0. So, the i-th entry gives the position in the array #_neighbours, in which the subarray storing the
    * node neighbours of node i starts. The degree (=number of neighbours) of node i is given by
    * _index[i+1] - _index[i]. (Difference to the data structure described in the MPI-2.2 standard: There, the 0-th
    * entry is omitted. However, adding this entry eases array access since the first node does not have to be handled
    * in a special way.) For an example see #_neighbours.
    *
    * Dimension: [#_num_nodes+1]
    */
    Index* _index;

    /**
    * \brief node neighbours within the graph (i.e. edges of the graph), represented as a list of node numbers
    *
    * The neighbours of node i are stored in the subarray _neighbours[#_index[i]], ..., _neighbours[#_index[i+1]-1].
    * The order of the neighbours within the subarrays is arbitrary.
    *
    * Example: domain consisting of 7 subdomains, each subdomain is a node in the graph. Graph edges represent
    * neighbouring subdomains, including diagonal neighbours.
    * \verbatim
    *   -----------------
    *   | 0 | 1 | 2 | 3 |
    *   -----------------
    *   | 4 | 5 |
    *   ---------
    *   | 6 |
    *   -----
    *   nodes     neighbours
    *   0         4,5,1
    *   1         0,4,5,2
    *   2         1,5,3
    *   3         2
    *   4         0,1,5,6
    *   5         0,1,2,4,6
    *   6         4,5
    * _num_nodes = 7
    * _index      = [0,        3,           7,      10, 11,          15,              20,   22]
    * _neighbours = [4, 5, 1,  0, 4, 5, 2,  1, 5, 3, 2,  0, 1, 5, 6,  0, 1, 2, 4,  6,  4, 5]
    * \endverbatim
    *
    * Dimension: [total number of neighbours] = [number of edges] = [#_index[#_num_nodes]]
    */
    Index* _neighbours;


  public:

    /* *****************
    * member variables *
    *******************/

    /* *************************
    * constructor & destructor *
    ***************************/
    /// CTOR
    Graph(
      Index const num_nodes,
      Index const* index,
      Index const* neighbours)
      : _num_nodes(num_nodes),
        _index(nullptr),
        _neighbours(nullptr)
    {
      CONTEXT("Graph::Graph()");
      // copy index array
      _index = new Index[num_nodes+1];
      for(Index i(0) ; i < num_nodes+1 ; ++i)
      {
        _index[i] = index[i];
      }

      // copy neighbour array (its size is given by _index[_num_nodes])
      _neighbours = new Index[_index[_num_nodes]];
      for(Index i(0) ; i < _index[_num_nodes] ; ++i)
      {
        _neighbours[i] = neighbours[i];
      }
    }

    /// DTOR
    ~Graph()
    {
      CONTEXT("Graph::~Graph()");
      delete [] _index;
      _index = nullptr;
      delete [] _neighbours;
      _neighbours = nullptr;
    }

    /* ******************
    * getters & setters *
    ********************/
    /**
    * \brief getter for the number of nodes
    *
    * \return number of nodes #_num_nodes
    */
    inline Index num_nodes() const
    {
      CONTEXT("Graph::num_nodes()");
      return _num_nodes;
    }

    /**
    * \brief getter for the index array
    *
    * \return pointer to the index array #_index
    */
    inline Index* index() const
    {
      CONTEXT("Graph::index()");
      return _index;
    }

    /**
    * \brief getter for the neighbours array
    *
    * \return pointer to the edge array #_neighbours
    */
    inline Index* neighbours() const
    {
      CONTEXT("Graph::neighbours()");
      return _neighbours;
    }

    /* *****************
    * member functions *
    *******************/
    /// print the graph to the given stream
    void print(std::ostream& stream) const
    {
      CONTEXT("Graph::print()");
      stream << "number of nodes: " << _num_nodes << std::endl;
      if (_num_nodes > 0)
      {
        stream << "node | degree | neighbours: " << std::endl;
        for(Index i(0) ; i < _num_nodes ; ++i)
        {
          stream << i << " | " << _index[i+1] - _index[i];
          if (_index[i+1] - _index[i] > 0)
          {
            stream << " | " << _neighbours[_index[i]];
            for(Index j(_index[i]+1) ; j < _index[i+1] ; ++j)
            {
              stream << ", " << _neighbours[j];
            }
          }
          stream << std::endl;
        }
      }
    }

    /// returns the graph print as string
    inline String print() const
    {
      CONTEXT("Graph::print()");
      std::ostringstream oss;
      print(oss);
      return oss.str();
    }
  }; // class Graph
#else // OLD_GRAPH

  /**
   * \brief Adjacency Graph implementation
   * \todo detailed description
   * \author Peter Zajac
   */
  class Graph
  {
  public:
    /**
     * \brief Dual node iterator for Graph class
     */
    typedef const Index* DualIterator;

    enum Share
    {
      share,
      copy,
      own
    };


  protected:
    /// total number of primal nodes
    Index _num_nodes_primal;
    /// total number of dual nodes
    Index _num_nodes_dual;
    /// total number of dual node indices
    Index _num_indices_dual;

    /**
     * \brief Primal pointer array
     *
     * Dimension: #_num_nodes_primal+1
     */
    Index* _primal_ptr;

    /**
     * \brief Primal end-pointer array
     *
     * Dimension: #_num_nodes_primal
     */
    Index* _primal_end;

    /**
     * \brief Dual node index array
     *
     * Dimension: #_num_indices_dual
     */
    Index* _dual_idx;

    /**
     * \brief Specifies whether the arrays are shared or not
     */
    bool _shared;

    /// \cond internal

    /// \endcond

  public:

    /**
     * \brief Default constructor.
     *
     * This constructor creates a new empty graph, but does not allocate any arrays.
     */
    Graph() :
      _num_nodes_primal(0),
      _num_nodes_dual(0),
      _num_indices_dual(0),
      _primal_ptr(nullptr),
      _primal_end(nullptr),
      _dual_idx(nullptr),
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
     * \param[in] num_nodes_primal
     * The total number of primal nodes for the graph.
     *
     * \param[in] num_nodes_dual
     * The total number of dual nodes for the graph.
     *
     * \param[in] num_indices_dual
     * The total number of dual node indices for the graph.
     *
     * \param[in] alloc_primal_end
     * If \c true, then the #_primal_end array will be allocated, otherwise #_primal_end is set to \c nullptr.
     */
    Graph(
      Index num_nodes_primal,
      Index num_nodes_dual,
      Index num_indices_dual,
      bool alloc_primal_end = false)
       :
      _num_nodes_primal(num_nodes_primal),
      _num_nodes_dual(num_nodes_dual),
      _num_indices_dual(num_indices_dual),
      _primal_ptr(nullptr),
      _primal_end(nullptr),
      _dual_idx(nullptr),
      _shared(false)
    {
      _primal_ptr = new Index[_num_nodes_primal+1];
      if(alloc_primal_end)
      {
        _primal_end = new Index[_num_nodes_primal];
      }
      _dual_idx = new Index[_num_indices_dual];
    }

    /**
     * \brief Constructor
     *
     * This constructor creates a new graph based on the arrays given.
     *
     * \param[in] num_nodes_primal
     * The total number of primal nodes for the graph.
     *
     * \param[in] num_nodes_dual
     * The total number of dual nodes for the graph.
     *
     * \param[in] num_indices_dual
     * The total number of dual node indices for the graph.
     *
     * \param[in] primal_ptr
     * The primal pointer array for the graph. Must not be \c nullptr.
     *
     * \param[in] primal_end
     * The primal end-pointer array for the graph. May be \c nullptr.
     *
     * \param[in] dual_idx
     * The dual node index array for the graph. Must not be \c nullptr.
     *
     * \param[in] share
     * If \c false, then the graph will allocate its own #_primal_ptr, #_primal_end and #_dual_idx arrays
     * and copy the arrays passed to the constructor.\n
     * If \c true, the graph will use the arrays passed to this constructor.
     */
    Graph(
      Index num_nodes_primal,
      Index num_nodes_dual,
      Index num_indices_dual,
      Index* primal_ptr,
      Index* primal_end,
      Index* dual_idx,
      Share share = Graph::copy)
       :
      _num_nodes_primal(num_nodes_primal),
      _num_nodes_dual(num_nodes_dual),
      _num_indices_dual(num_indices_dual),
      _primal_ptr(primal_ptr),
      _primal_end(primal_end),
      _dual_idx(dual_idx),
      _shared(share == Graph::share)
    {
      CONTEXT("Graph::Graph()");
      if(share == Graph::copy)
      {
        // we need to make copies of the arrays
        _primal_ptr = new Index[num_nodes_primal+1];
        for(Index i(0); i <= num_nodes_primal; ++i)
        {
          _primal_ptr[i] = primal_ptr[i];
        }
        if(primal_end != nullptr)
        {
          _primal_end = new Index[num_nodes_primal];
          for(Index i(0); i < num_nodes_primal; ++i)
          {
            _primal_end[i] = primal_end[i];
          }
        }
        _dual_idx = new Index[num_indices_dual];
        for(Index i(0); i < num_indices_dual; ++i)
        {
          _dual_idx[i] = dual_idx[i];
        }
      }
    }

    /// DTOR
    virtual ~Graph()
    {
      CONTEXT("Graph::destroy()");
      if(!_shared)
      {
        if(_dual_idx != nullptr)
          delete [] _dual_idx;
        if(_primal_ptr != nullptr)
          delete [] _primal_ptr;
        if(_primal_end != nullptr)
          delete [] _primal_end;
      }
      /*_num_nodes_primal = 0;
      _num_nodes_dual = 0;
      _num_indices_dual = 0;
      _primal_ptr = nullptr;
      _primal_end = nullptr;
      _dual_idx = nullptr;
      _shared = false;*/
    }

    /**
     * \brief Returns the degree of a primal node.
     *
     * This function returns the degree of a primal node, i.e. the number of dual nodes which are adjacent to the
     * specified primal node.
     *
     * \note The degree might be greater than the total number of dual nodes if the graph is not injective.
     *
     * \param[in] primal_node
     * The primal node index whose degree is to be returned.
     *
     * \returns
     * The degree of the specified primal node.
     */
    inline Index degree(Index primal_node) const
    {
      CONTEXT("Graph::degree()");
      ASSERT(primal_node < _num_nodes_primal, "Primal node index out of range");
      if(_primal_end != nullptr)
        return _primal_end[primal_node] - _primal_ptr[primal_node];
      else
        return _primal_ptr[primal_node+1] - _primal_ptr[primal_node];
    }

    /**
     * \brief Returns the degree of the graph.
     *
     * This function returns the degree of the graph, i.e. the maximum of the degrees of all primal nodes.
     *
     * \note The degree of the graph might be greater than the total number of dual nodes if the graph is not
     * injective.
     *
     * \note This function performs a loop over all primal nodes and therefore has linear runtime. It is therefore
     * recommended to store the result of this function in a variable if the degree is needed multiple times.
     *
     * \returns
     * The degree of the graph.
     */
    inline Index degree() const
    {
      CONTEXT("Graph::degree()");
      Index deg = 0;
      if(_primal_end != nullptr)
      {
        for(Index i(0); i < _num_nodes_primal; ++i)
        {
          deg = std::max(deg, _primal_end[i] - _primal_ptr[i]);
        }
      }
      else
      {
        for(Index i(0); i < _num_nodes_primal; ++i)
        {
          deg = std::max(deg, _primal_ptr[i+1] - _primal_ptr[i]);
        }
      }
      return deg;
    }

    /**
     * \brief Returns the primal pointer array.
     * \returns The primal pointer array.
     */
    inline Index* get_primal_ptr()
    {
      return _primal_ptr;
    }

    /** \copydoc get_primal_ptr() */
    inline const Index* get_primal_ptr() const
    {
      return _primal_ptr;
    }

    /**
     * \brief Returns the primal end-pointer array.
     *
     * \warning Please note that a graph does not necessarily have a primal end-pointer array, thus this function
     * might return a \c nullptr even if the graph is not empty.
     *
     * \returns The primal end-pointer array.
     */
    inline Index* get_primal_end()
    {
      return _primal_end;
    }

    /** \copydoc get_primal_end() */
    inline const Index* get_primal_end() const
    {
      return _primal_end;
    }

    /**
     * \brief Returns the dual node index array.
     * \returns The dual node index array.
     */
    inline Index* get_dual_idx()
    {
      return _dual_idx;
    }

    /** \copydoc get_dual_idx() */
    inline const Index* get_dual_idx() const
    {
      return _dual_idx;
    }

    /**
     * \brief Specifies whether the graph's arrays are shared or not.
     *
     * \returns \c true, if the graph's arrays are shared, otherwise \c false.
     */
    inline bool is_shared() const
    {
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
      stream << "number of primal nodes: " << _num_nodes_primal << std::endl;
      stream << "number of dual nodes  : " << _num_nodes_dual << std::endl;
      stream << "number of dual indices: " << _num_indices_dual << std::endl;

      const Index* _ptr = _primal_ptr;
      const Index* _end = nullptr;
      if(_primal_end != nullptr)
        _end = _primal_end;
      else
        _end = &_primal_ptr[1];

      if (_num_nodes_primal > 0)
      {
        stream << "node | degree | neighbours: " << std::endl;
        for(Index i(0) ; i < _num_nodes_primal ; ++i)
        {

          stream << i << " | " << (_end[i] - _ptr[i]);
          if (_end[i] - _ptr[i] > 0)
          {
            stream << " | " << _dual_idx[_ptr[i]];
            for(Index j(_ptr[i]+1) ; j < _end[i] ; ++j)
            {
              stream << ", " << _dual_idx[j];
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

    /* ********************************************************************* */
    /*  A D J U N C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
    /* ********************************************************************* */
    inline Index num_nodes_primal() const
    {
      CONTEXT("Graph::num_nodes_primal()");
      return _num_nodes_primal;
    }

    inline Index num_nodes_dual() const
    {
      CONTEXT("Graph::num_nodes_dual()");
      return _num_nodes_dual;
    }

    inline DualIterator dual_begin(Index primal_node) const
    {
      CONTEXT("Graph::dual_begin()");
      ASSERT(primal_node < _num_nodes_primal, "Primal node index out of range");
      return &_dual_idx[_primal_ptr[primal_node]];
    }

    inline DualIterator dual_end(Index primal_node) const
    {
      CONTEXT("Graph::dual_end()");
      ASSERT(primal_node < _num_nodes_primal, "Primal node index out of range");

      if(_primal_end != nullptr)
        return &_dual_idx[_primal_end[primal_node]];
      else
        return &_dual_idx[_primal_ptr[primal_node+1]];
    }
  }; // class Graph

  /**
   * \brief Graph renderer implementation
   */
  template<typename Adjunctor_>
  class GraphRenderer
  {
  protected:
    /// \cond internal
    template<typename Adjunctor_>
    static Index _aux_inj(const Adjunctor_& adj, Index i, Index* idx)
    {
      CONTEXT("GraphRenderer::_aux_inj()");

      Index num_idx = 0;
      typename Adjunctor_::DualIterator cur(adj.dual_begin(i));
      typename Adjunctor_::DualIterator end(adj.dual_end(i));
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

    /// renders adjunctor
    template<typename Adjunctor_>
    static Graph* _render(const Adjunctor_& adj)
    {
      CONTEXT("GraphRenderer::_render()");

      typedef typename Adjunctor_::DualIterator ADuIt;

      // get counts
      Index num_nodes_primal = adj.num_nodes_primal();
      Index num_nodes_dual = adj.num_nodes_dual();
      Index num_indices_dual = 0;

      // allocate pointer array
      Index* primal_ptr = new Index[num_nodes_primal+1];

      // count number of adjacencies and build pointer array
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        primal_ptr[i] = num_indices_dual;
        ADuIt cur(adj.dual_begin(i));
        ADuIt end(adj.dual_end(i));
        for(; cur != end; ++cur)
        {
          ++num_indices_dual;
        }
      }
      primal_ptr[num_nodes_primal] = num_indices_dual;

      // allocate and build index array
      Index* dual_idx = new Index[num_indices_dual];
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        Index* idx = &dual_idx[primal_ptr[i]];
        ADuIt cur(adj.dual_begin(i));
        ADuIt end(adj.dual_end(i));
        for(; cur != end; ++cur, ++idx)
        {
          *idx = *cur;
        }
      }

      // create the Graph
      return new Graph(num_nodes_primal, num_nodes_dual, num_indices_dual,
        primal_ptr, nullptr, dual_idx, Graph::own);
    }

    /// renders injectified adjunctor
    template<typename Adjunctor_>
    static Graph* _render_inj(const Adjunctor_& adj)
    {
      CONTEXT("GraphRenderer::_render_inj()");

      // get counts
      Index num_nodes_primal = adj.num_nodes_primal();
      Index num_nodes_dual = adj.num_nodes_dual();
      Index num_indices_dual = 0;

      // allocate pointer array
      Index* primal_ptr = new Index[num_nodes_primal+1];

      // allocate auxiliary index array
      Index* aux = new Index[num_nodes_dual];

      // count number of adjacencies and build pointer array
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        primal_ptr[i] = num_indices_dual;
        num_indices_dual += _aux_inj(adj, i, aux);
      }
      primal_ptr[num_nodes_primal] = num_indices_dual;

      // delete auxiliary array
      delete [] aux;

      // allocate and build index array
      Index* dual_idx = new Index[num_indices_dual];
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        _aux_inj(adj, i, &dual_idx[primal_ptr[i]]);
      }

      // create graph
      return new Graph(num_nodes_primal, num_nodes_dual, num_indices_dual,
        primal_ptr, nullptr, dual_idx, Graph::own);
    }

    /// renders transposed adjunctor
    template<typename Adjunctor_>
    static Graph* _render_trans(const Adjunctor_& adj)
    {
      CONTEXT("GraphRenderer::_render_trans()");

      typedef typename Adjunctor_::DualIterator ADuIt;

      // get counts
      Index num_nodes_primal = adj.num_nodes_dual();
      Index num_nodes_dual = adj.num_nodes_primal();
      Index num_indices_dual = 0;

      // allocate and format pointer array
      Index *primal_ptr = new Index[num_nodes_primal+1];
      for(Index i(0); i <= num_nodes_primal; ++i)
      {
        primal_ptr[i] = 0;
      }

      // count number of adjacencies
      for(Index j(0); j < num_nodes_dual; ++j)
      {
        ADuIt cur(adj.dual_begin(j));
        ADuIt end(adj.dual_end(j));
        for(; cur != end; ++cur)
        {
          ++primal_ptr[(*cur)+1];
        }
      }

      // build pointer array
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        primal_ptr[i+1] += primal_ptr[i];
      }
      num_indices_dual = primal_ptr[num_nodes_primal];

      // allocate and build index array
      Index* dual_idx = new Index[num_indices_dual];
      Index** dual_ptr = new Index*[num_nodes_primal];
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        dual_ptr[i] = &dual_idx[primal_ptr[i]];
      }

      for(Index j(0); j < num_nodes_dual; ++j)
      {
        ADuIt cur(adj.dual_begin(j));
        ADuIt end(adj.dual_end(j));
        for(; cur != end; ++cur)
        {
          Index*& idx = dual_ptr[*cur];
          *idx = j;
          ++idx;
        }
      }

      delete [] dual_ptr;

      // create graph
      return new Graph(num_nodes_primal, num_nodes_dual, num_indices_dual,
        primal_ptr, nullptr, dual_idx, Graph::own);
    }

    /// renders transposed injectified adjunctor
    template<typename Adjunctor_>
    static Graph* _render_trans_inj(const Adjunctor_& adj)
    {
      CONTEXT("GraphRenderer::_render_trans_inj()");

      // get counts
      Index num_nodes_primal = adj.num_nodes_dual();
      Index num_nodes_dual = adj.num_nodes_primal();
      Index num_indices_dual = 0;

      // allocate pointer array
      Index* primal_ptr = new Index[num_nodes_primal+1];
      for(Index i(0); i <= num_nodes_primal; ++i)
      {
        primal_ptr[i] = 0;
      }

      // allocate auxiliary index array
      Index* aux = new Index[num_nodes_primal];

      // loop over all dual nodes
      for(Index j(0); j < num_nodes_dual; ++j)
      {
        Index num_aux = _aux_inj(adj, j, aux);
        for(Index k(0); k < num_aux; ++k)
        {
          ++primal_ptr[aux[k]+1];
        }
        num_indices_dual += num_aux;
      }

      Index* dual_idx = new Index[num_indices_dual];
      Index** dual_ptr = new Index*[num_nodes_primal];

      // build pointer arrays
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        primal_ptr[i+1] += primal_ptr[i];
        dual_ptr[i] = &dual_idx[primal_ptr[i]];
      }

      // build dual index array
      for(Index j(0); j < num_nodes_dual; ++j)
      {
        Index num_aux = _aux_inj(adj, j, aux);
        for(Index k(0); k < num_aux; ++k)
        {
          Index*& idx = dual_ptr[aux[k]];
          *idx = j;
          ++idx;
        }
      }

      // delete auxiliary arrays
      delete [] dual_ptr;
      delete [] aux;

      // create graph
      return new Graph(num_nodes_primal, num_nodes_dual, num_indices_dual,
        primal_ptr, nullptr, dual_idx, Graph::own);
    }

    /// \endcond
  public:
    /**
     * \brief Renders an Adjunctor into a Graph.
     *
     * \todo detailed description
     *
     * \tparam Adjunctor_
     * A class implementing the Adjunctor interface. This parameter is determined automatically.
     *
     * \param[in] adj
     * A reference to the adjunctor which is to be rendered.
     *
     * \param[in] transpose
     * If \c true, then the transposed adjunctor will be rendered into the graph.
     *
     * \param[in] injectify
     * If \c true, then the injectified adjunctor will be rendered into the graph.
     *
     * \returns
     * A pointer to a Graph containing the rendered adjunctor.
     */
    static Graph* render(
      const Adjunctor_& adj,
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
   * \brief Partial specialisation of GraphRenderer for CompositeAdjunctor
   */
  template<
    typename Adjunctor1_,
    typename Adjunctor2_>
  class GraphRenderer< CompositeAdjunctor<Adjunctor1_, Adjunctor2_> >
  {
  protected:
    /// \cond internal
    static Index _aux_inj(const Adjunctor1_& adj1, const Adjunctor2_& adj2, Index i, Index* idx)
    {
      CONTEXT("GraphRenderer<CompositeAdjunctor>::_aux_inj()");

      Index num_idx = 0;
      typename Adjunctor1_::DualIterator cur1(adj1.dual_begin(i));
      typename Adjunctor1_::DualIterator end1(adj1.dual_end(i));
      for(; cur1 != end1; ++cur1)
      {
        typename Adjunctor2_::DualIterator cur2(adj2.dual_begin(*cur1));
        typename Adjunctor2_::DualIterator end2(adj2.dual_end(*cur1));

        // loop over all dual node indices
        for(; cur2 != end2; ++cur2)
        {
          Index jdx = *cur2;

          // check if we already have that dual node index
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

    /// renders adjunctor composition
    static Graph* _render(
      const Adjunctor1_& adj1,
      const Adjunctor2_& adj2)
    {
      CONTEXT("GraphRenderer<CompositeAdjunctor>::_render()");

      typedef typename Adjunctor1_::DualIterator ADuIt1;
      typedef typename Adjunctor2_::DualIterator ADuIt2;

      // get counts
      Index num_nodes_primal = adj1.num_nodes_primal();
      Index num_nodes_dual = adj2.num_nodes_dual();
      Index num_indices_dual = 0;

      // allocate pointer array
      Index* primal_ptr = new Index[num_nodes_primal+1];

      // count number of adjacencies and build pointer array
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        primal_ptr[i] = num_indices_dual;
        ADuIt1 cur1(adj1.dual_begin(i));
        ADuIt1 end1(adj1.dual_end(i));
        for(; cur1 != end1; ++cur1)
        {
          ADuIt2 cur2(adj2.dual_begin(*cur1));
          ADuIt2 end2(adj2.dual_end(*cur1));
          for(; cur2 != end2; ++cur2)
          {
            ++num_indices_dual;
          }
        }
      }
      primal_ptr[num_nodes_primal] = num_indices_dual;

      // allocate and build index array
      Index* dual_idx = new Index[num_indices_dual];
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        Index* idx = &dual_idx[primal_ptr[i]];
        ADuIt1 cur1(adj1.dual_begin(i));
        ADuIt1 end1(adj1.dual_end(i));
        for(; cur1 != end1; ++cur1)
        {
          ADuIt2 cur2(adj2.dual_begin(*cur1));
          ADuIt2 end2(adj2.dual_end(*cur1));
          for(; cur2 != end2; ++cur2, ++idx)
          {
            *idx = *cur2;
          }
        }
      }

      // create graph
      return new Graph(num_nodes_primal, num_nodes_dual, num_indices_dual,
        primal_ptr, nullptr, dual_idx, Graph::own);
    }

    /// renders injectified adjunctor composition
    static Graph* _render_inj(
      const Adjunctor1_& adj1,
      const Adjunctor2_& adj2)
    {
      CONTEXT("GraphRenderer<CompositeAdjunctor>::_render_inj()");

      typedef typename Adjunctor1_::DualIterator ADuIt1;
      typedef typename Adjunctor2_::DualIterator ADuIt2;

      // get counts
      Index num_nodes_primal = adj1.num_nodes_primal();
      Index num_nodes_dual = adj2.num_nodes_dual();
      Index num_indices_dual = 0;

      // allocate pointer array
      Index* primal_ptr = new Index[num_nodes_primal+1];

      // allocate auxiliary array
      Index* aux = new Index[num_nodes_dual];

      // count number of adjacencies and build pointer array
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        primal_ptr[i] = num_indices_dual;
        num_indices_dual += _aux_inj(adj1, adj2, i, aux);
      }
      primal_ptr[num_nodes_primal] = num_indices_dual;

      // delete auxiliary array
      delete [] aux;

      // allocate and build index array
      Index* dual_idx = new Index[num_indices_dual];
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        _aux_inj(adj1, adj2, i, &dual_idx[primal_ptr[i]]);
      }

      // create graph
      return new Graph(num_nodes_primal, num_nodes_dual, num_indices_dual,
        primal_ptr, nullptr, dual_idx, Graph::own);
    }

    /// renders transposed adjunctor composition
    static Graph* _render_trans(
      const Adjunctor1_& adj1,
      const Adjunctor2_& adj2)
    {
      CONTEXT("GraphRenderer<CompositeAdjunctor>::_render_trans()");

      typedef typename Adjunctor1_::DualIterator ADuIt1;
      typedef typename Adjunctor2_::DualIterator ADuIt2;

      // get counts
      Index num_nodes_primal = adj2.num_nodes_dual();
      Index num_nodes_dual = adj1.num_nodes_primal();
      Index num_indices_dual = 0;

      // allocate and format pointer array
      Index *primal_ptr = new Index[num_nodes_primal+1];
      for(Index i(0); i <= num_nodes_primal; ++i)
      {
        primal_ptr[i] = 0;
      }

      // count number of adjacencies
      for(Index j(0); j < num_nodes_dual; ++j)
      {
        ADuIt1 cur1(adj1.dual_begin(j));
        ADuIt1 end1(adj1.dual_end(j));
        for(; cur1 != end1; ++cur1)
        {
          ADuIt2 cur2(adj2.dual_begin(*cur1));
          ADuIt2 end2(adj2.dual_end(*cur1));
          for(; cur2 != end2; ++cur2)
          {
            ++primal_ptr[(*cur2)+1];
          }
        }
      }

      // build pointer array
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        primal_ptr[i+1] += primal_ptr[i];
      }
      num_indices_dual = primal_ptr[num_nodes_primal];

      // allocate and build index array
      Index* dual_idx = new Index[num_indices_dual];
      Index** dual_ptr = new Index*[num_nodes_primal];
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        dual_ptr[i] = &dual_idx[primal_ptr[i]];
      }

      for(Index j(0); j < num_nodes_dual; ++j)
      {
        ADuIt1 cur1(adj1.dual_begin(j));
        ADuIt1 end1(adj1.dual_end(j));
        for(; cur1 != end1; ++cur1)
        {
          ADuIt2 cur2(adj2.dual_begin(*cur1));
          ADuIt2 end2(adj2.dual_end(*cur1));
          for(; cur2 != end2; ++cur2)
          {
            Index*& idx = dual_ptr[*cur2];
            *idx = j;
            ++idx;
          }
        }
      }

      delete [] dual_ptr;

      // create graph
      return new Graph(num_nodes_primal, num_nodes_dual, num_indices_dual,
        primal_ptr, nullptr, dual_idx, Graph::own);
    }

    /// renders transposed injectified adjunctor composition
    static Graph* _render_trans_inj(
      const Adjunctor1_& adj1,
      const Adjunctor2_& adj2)
    {
      CONTEXT("GraphRenderer<CompositeAdjunctor>::_render_trans_inj()");

      // get counts
      Index num_nodes_primal = adj2.num_nodes_dual();
      Index num_nodes_dual = adj1.num_nodes_primal();
      Index num_indices_dual = 0;

      // allocate pointer array
      Index* primal_ptr = new Index[num_nodes_primal+1];
      for(Index i(0); i <= num_nodes_primal; ++i)
      {
        primal_ptr[i] = 0;
      }

      // allocate auxiliary index array
      Index* aux = new Index[num_nodes_primal];

      // loop over all dual nodes
      for(Index j(0); j < num_nodes_dual; ++j)
      {
        Index num_aux = _aux_inj(adj1, adj2, j, aux);
        for(Index k(0); k < num_aux; ++k)
        {
          ++primal_ptr[aux[k]+1];
        }
        num_indices_dual += num_aux;
      }

      Index* dual_idx = new Index[num_indices_dual];
      Index** dual_ptr = new Index*[num_nodes_primal];

      // build pointer arrays
      for(Index i(0); i < num_nodes_primal; ++i)
      {
        primal_ptr[i+1] += primal_ptr[i];
        dual_ptr[i] = &dual_idx[primal_ptr[i]];
      }

      // build dual index array
      for(Index j(0); j < num_nodes_dual; ++j)
      {
        Index num_aux = _aux_inj(adj1, adj2, j, aux);
        for(Index k(0); k < num_aux; ++k)
        {
          Index*& idx = dual_ptr[aux[k]];
          *idx = j;
          ++idx;
        }
      }

      // delete auxiliary arrays
      delete [] dual_ptr;
      delete [] aux;

      // create graph
      return new Graph(num_nodes_primal, num_nodes_dual, num_indices_dual,
        primal_ptr, nullptr, dual_idx, Graph::own);
    }
    /// \endcond

  public:
    /**
     * \brief Renders a composition of two Adjunctors into a Graph.
     *
     * \todo detailed description
     *
     * \param[in] adj1
     * A reference to the first adjunctor of the composition which is to be rendered.
     *
     * \param[in] adj2
     * A reference to the second adjunctor of the composition which is to be rendered.
     *
     * \param[in] transpose
     * If \c true, then the transposed adjunctor will be rendered into the graph.
     *
     * \param[in] injectify
     * If \c true, then the injectified adjunctor will be rendered into the graph.
     *
     * \returns
     * A pointer to a Graph containing the rendered adjunctor composition.
     */
    static Graph* render(
      const Adjunctor1_& adj1,
      const Adjunctor2_& adj2,
      bool injectify = true,
      bool transpose = false)
    {
      CONTEXT("GraphRenderer<CompositeAdjunctor>::render()");
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
     * \brief Renders a CompositeAdjunctor into a Graph.
     *
     * \todo detailed description
     *
     * \param[in] adj
     * A reference to the CompositeAdjunctor which is to be rendered.
     *
     * \param[in] transpose
     * If \c true, then the transposed adjunctor will be rendered into the graph.
     *
     * \param[in] injectify
     * If \c true, then the injectified adjunctor will be rendered into the graph.
     *
     * \returns
     * A pointer to a Graph containing the rendered CompositeAdjunctor.
     */
    static Graph* render(
      const CompositeAdjunctor<Adjunctor1_,Adjunctor2_>& adj,
      bool injectify = true,
      bool transpose = false)
    {
      CONTEXT("GraphRenderer<CompositeAdjunctor>::render()");
      return render(adj.get_adjunctor1(), adj.get_adjunctor2(), injectify, transpose);
    }
  }; // class GraphRenderer< CompositeAdjunctor<...> >

#endif // OLD_GRAPH
} // namespace FEAST

#endif // KERNEL_GRAPH_HPP
