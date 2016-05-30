#pragma once
#ifndef KERNEL_ADJACENCY_GRAPH_HPP
#define KERNEL_ADJACENCY_GRAPH_HPP 1

// includes, FEAT
#include <kernel/adjacency/base.hpp>
#include <kernel/adjacency/adjactor.hpp>
#include <kernel/util/assertion.hpp>

namespace FEAT
{
  namespace Adjacency
  {
    // forward declaration
    class Permutation;

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
       * For the purpose of the Adjactor interface implementation a <c>const Index</c>-pointer is the optimal choice
       * for the image node iterator.
       */
      typedef const Index* ImageIterator;

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
       * This member specifies whether the Graph object will delete the #_domain_ptr, #_domain_end and #_image_idx
       * arrays within the destructor.
       */
      bool _shared;

    public:

      /**
       * \brief Default constructor.
       *
       * This constructor creates a new empty graph, but does not allocate any arrays.
       */
      Graph();

      /**
       * \brief Allocation Constructor.
       *
       * This constructor creates a new graph and allocates the Graph's arrays to the corresponding lengths.
       *
       * \note This constructor does not initialise the allocated arrays -- they have to be initialised by the user
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
        bool alloc_domain_end = false);

      /**
       * \brief "Using-Arrays" Constructor
       *
       * This constructor creates a new graph using the arrays given to this constructor.
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
       * \param[in] shared
       * Specifies whether the graph's arrays are shared or not.
       *  - If set to \c false then the graph will deallocate all arrays passed to this constructor
       *    upon destruction.
       *  - If set to \c true then the caller remains responsible for the deallocation of the arrays.
       */
      Graph(
        Index num_nodes_domain,
        Index num_nodes_image,
        Index num_indices_image,
        Index* domain_ptr,
        Index* domain_end,
        Index* image_idx,
        bool shared);

      /**
       * \brief "Copy-Arrays" Constructor
       *
       * This constructor creates a new graph using copies of the arrays passed to this function.
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
       */
      Graph(
        Index num_nodes_domain,
        Index num_nodes_image,
        Index num_indices_image,
        const Index* domain_ptr,
        const Index* domain_end,
        const Index* image_idx);

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
      explicit Graph(
        RenderType render_type,
        const Adjactor_& adjactor)
         :
        _num_nodes_domain(0),
        _num_nodes_image(0),
        _num_indices_image(0),
        _domain_ptr(nullptr),
        _domain_end(nullptr),
        _image_idx(nullptr),
        _shared(false)
      {
        switch(render_type)
        {
        case rt_as_is:
          _render_as_is(adjactor);
          break;

        case rt_injectify:
          _render_injectify(adjactor);
          break;

        case rt_transpose:
          _render_transpose(adjactor);
          break;

        case rt_injectify_transpose:
          _render_injectify_transpose(adjactor);
          break;

        default:
          throw InternalError("Invalid render_type parameter!");
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
      explicit Graph(
        RenderType render_type,
        const Adjactor1_& adjactor1,
        const Adjactor2_& adjactor2)
         :
        _num_nodes_domain(0),
        _num_nodes_image(0),
        _num_indices_image(0),
        _domain_ptr(nullptr),
        _domain_end(nullptr),
        _image_idx(nullptr),
        _shared(false)
      {
        switch(render_type)
        {
        case rt_as_is:
          _render_as_is(adjactor1, adjactor2);
          break;

        case rt_injectify:
          _render_injectify(adjactor1, adjactor2);
          break;

        case rt_transpose:
          _render_transpose(adjactor1, adjactor2);
          break;

        case rt_injectify_transpose:
          _render_injectify_transpose(adjactor1, adjactor2);
          break;

        default:
          throw InternalError("Invalid render_type parameter!");
        }
      }

      /// move CTOR
      Graph(Graph&& other);

      /// move-assign operator
      Graph& operator=(Graph&& other);

      /**
       * \brief "Permutation" copy CTOR
       *
       * This constructor creates a new graph by permuting the domain and image set
       * of the graph that is passed to this function.
       *
       * \param[in] other
       * The graph that has to be permuted.
       *
       * \param[in] domain_perm
       * The permutation of the domain set.
       *
       * \param[in] image_perm
       * The permutation of the image set.
       */
      Graph(const Graph& other, const Permutation& domain_perm, const Permutation& image_perm);

      /// virtual destructor
      virtual ~Graph();

      /**
       * \brief Clones this graph.
       *
       * \returns A deep-copy of this graph.
       */
      Graph clone() const
      {
        if(_domain_ptr != nullptr)
          return Graph(_num_nodes_domain, _num_nodes_image, _num_indices_image, _domain_ptr, _domain_end, _image_idx);
        else
          return Graph();
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
      Index degree(Index domain_node) const
      {
        ASSERT(domain_node < _num_nodes_domain, "Domain node index out of range");
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
       * \attention This function performs a loop over all domain nodes and therefore has linear runtime. It is
       * recommended to store the result of this function in a variable if the degree is needed multiple times.
       *
       * \returns
       * The degree of the graph.
       */
      Index degree() const;

      /**
       * \brief Returns the domain pointer array.
       * \returns The domain pointer array.
       */
      Index* get_domain_ptr()
      {
        return _domain_ptr;
      }

      /** \copydoc get_domain_ptr() */
      const Index* get_domain_ptr() const
      {
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
      Index* get_domain_end()
      {
        return _domain_end;
      }

      /** \copydoc get_domain_end() */
      const Index* get_domain_end() const
      {
        return _domain_end;
      }

      /**
       * \brief Returns the image node index array.
       * \returns The image node index array.
       */
      Index* get_image_idx()
      {
        return _image_idx;
      }

      /** \copydoc get_image_idx() */
      const Index* get_image_idx() const
      {
        return _image_idx;
      }

      /**
       * \brief Returns the total number indices.
       *
       * \returns The total number of indices in the graph.
       */
      Index get_num_indices() const
      {
        return _num_indices_image;
      }

      /**
       * \brief Specifies whether the graph's arrays are shared or not.
       *
       * \returns \c true, if the graph's arrays are shared, otherwise \c false.
       */
      bool is_shared() const
      {
        return _shared;
      }

      /**
       * \brief Sorts the image indices to non-descending order.
       */
      void sort_indices();

      /* *************************************************** */
      /*  R E N D E R   F U N C T I O N   T E M P L A T E S  */
      /* *************************************************** */
    private:
      /// \cond internal
      template<typename Adjactor_>
      static Index _aux_inj(const Adjactor_& adj, Index i, Index idx[])
      {
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
      template<typename Adjactor_>
      void _render_as_is(const Adjactor_& adj)
      {
        typedef typename Adjactor_::ImageIterator AImIt;

        // get counts
        _num_nodes_domain = adj.get_num_nodes_domain();
        _num_nodes_image = adj.get_num_nodes_image();
        _num_indices_image = 0;

        // allocate pointer array
        _domain_ptr = new Index[_num_nodes_domain+1];

        // count number of adjacencies and build pointer array
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          _domain_ptr[i] = _num_indices_image;
          AImIt cur(adj.image_begin(i));
          AImIt end(adj.image_end(i));
          for(; cur != end; ++cur)
          {
            ++_num_indices_image;
          }
        }
        _domain_ptr[_num_nodes_domain] = _num_indices_image;

        // allocate and build index array
        _image_idx = new Index[_num_indices_image];
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          Index* idx = &_image_idx[_domain_ptr[i]];
          AImIt cur(adj.image_begin(i));
          AImIt end(adj.image_end(i));
          for(; cur != end; ++cur, ++idx)
          {
            *idx = *cur;
          }
        }

        // okay
      }

      /// renders injectified adjactor
      template<typename Adjactor_>
      void _render_injectify(const Adjactor_& adj)
      {
        // get counts
        _num_nodes_domain = adj.get_num_nodes_domain();
        _num_nodes_image = adj.get_num_nodes_image();
        _num_indices_image = 0;

        // allocate pointer array
        _domain_ptr = new Index[_num_nodes_domain + 1];

        // allocate auxiliary index array
        Index* aux = new Index[_num_nodes_image];

        // count number of adjacencies and build pointer array
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          _domain_ptr[i] = _num_indices_image;
          _num_indices_image += _aux_inj(adj, i, aux);
        }
        _domain_ptr[_num_nodes_domain] = _num_indices_image;

        // delete auxiliary array
        delete [] aux;

        // allocate and build index array
        _image_idx = new Index[_num_indices_image];
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          _aux_inj(adj, i, &_image_idx[_domain_ptr[i]]);
        }

        // okay
      }

      /// renders transposed adjactor
      template<typename Adjactor_>
      void _render_transpose(const Adjactor_& adj)
      {
        typedef typename Adjactor_::ImageIterator AImIt;

        // get counts
        _num_nodes_domain = adj.get_num_nodes_image();
        _num_nodes_image = adj.get_num_nodes_domain();
        _num_indices_image = 0;

        // allocate and format pointer array
        _domain_ptr = new Index[_num_nodes_domain + 1];
        for(Index i(0); i <= _num_nodes_domain; ++i)
        {
          _domain_ptr[i] = 0;
        }

        // count number of adjacencies
        for(Index j(0); j < _num_nodes_image; ++j)
        {
          AImIt cur(adj.image_begin(j));
          AImIt end(adj.image_end(j));
          for(; cur != end; ++cur)
          {
            ++_domain_ptr[(*cur) + 1];
          }
        }

        // build pointer array
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          _domain_ptr[i+1] += _domain_ptr[i];
        }
        _num_indices_image = _domain_ptr[_num_nodes_domain];

        // allocate and build index array
        _image_idx = new Index[_num_indices_image];
        Index** image_ptr = new Index*[_num_nodes_domain];
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          image_ptr[i] = &_image_idx[_domain_ptr[i]];
        }

        for(Index j(0); j < _num_nodes_image; ++j)
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

        // okay
      }

      /// renders transposed injectified adjactor
      template<typename Adjactor_>
      void _render_injectify_transpose(const Adjactor_& adj)
      {
        // get counts
        _num_nodes_domain = adj.get_num_nodes_image();
        _num_nodes_image = adj.get_num_nodes_domain();
        _num_indices_image = 0;

        // allocate pointer array
        _domain_ptr = new Index[_num_nodes_domain + 1];
        for(Index i(0); i <= _num_nodes_domain; ++i)
        {
          _domain_ptr[i] = 0;
        }

        // allocate auxiliary index array
        Index* aux = new Index[_num_nodes_domain];

        // loop over all image nodes
        for(Index j(0); j < _num_nodes_image; ++j)
        {
          Index num_aux = _aux_inj(adj, j, aux);
          for(Index k(0); k < num_aux; ++k)
          {
            ++_domain_ptr[aux[k]+1];
          }
          _num_indices_image += num_aux;
        }

        _image_idx = new Index[_num_indices_image];
        Index** image_ptr = new Index*[_num_nodes_domain];

        // build pointer arrays
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          _domain_ptr[i+1] += _domain_ptr[i];
          image_ptr[i] = &_image_idx[_domain_ptr[i]];
        }

        // build image index array
        for(Index j(0); j < _num_nodes_image; ++j)
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

        // okay
      }

      template<
        typename Adjactor1_,
        typename Adjactor2_>
      static Index _aux_inj(const Adjactor1_& adj1, const Adjactor2_& adj2, Index i, Index idx[])
      {
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
      template<
        typename Adjactor1_,
        typename Adjactor2_>
      void _render_as_is(
        const Adjactor1_& adj1,
        const Adjactor2_& adj2)
      {
        // validate adjactor dimensions
        ASSERT(adj1.get_num_nodes_image() <= adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

        typedef typename Adjactor1_::ImageIterator AImIt1;
        typedef typename Adjactor2_::ImageIterator AImIt2;

        // get counts
        _num_nodes_domain = adj1.get_num_nodes_domain();
        _num_nodes_image = adj2.get_num_nodes_image();
        _num_indices_image = 0;

        // allocate pointer array
        _domain_ptr = new Index[_num_nodes_domain + 1];

        // count number of adjacencies and build pointer array
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          _domain_ptr[i] = _num_indices_image;
          AImIt1 cur1(adj1.image_begin(i));
          AImIt1 end1(adj1.image_end(i));
          for(; cur1 != end1; ++cur1)
          {
            AImIt2 cur2(adj2.image_begin(*cur1));
            AImIt2 end2(adj2.image_end(*cur1));
            for(; cur2 != end2; ++cur2)
            {
              ++_num_indices_image;
            }
          }
        }
        _domain_ptr[_num_nodes_domain] = _num_indices_image;

        // allocate and build index array
        _image_idx = new Index[_num_indices_image];
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          Index* idx = &_image_idx[_domain_ptr[i]];
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

        // okay
      }

      /// renders injectified adjactor composition
      template<
        typename Adjactor1_,
        typename Adjactor2_>
      void _render_injectify(
        const Adjactor1_& adj1,
        const Adjactor2_& adj2)
      {
        // validate adjactor dimensions
        ASSERT(adj1.get_num_nodes_image() <= adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

        // get counts
        _num_nodes_domain = adj1.get_num_nodes_domain();
        _num_nodes_image = adj2.get_num_nodes_image();
        _num_indices_image = 0;

        // allocate pointer array
        _domain_ptr = new Index[_num_nodes_domain + 1];

        // allocate auxiliary array
        Index* aux = new Index[_num_nodes_image];

        // count number of adjacencies and build pointer array
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          _domain_ptr[i] = _num_indices_image;
          _num_indices_image += _aux_inj(adj1, adj2, i, aux);
        }
        _domain_ptr[_num_nodes_domain] = _num_indices_image;

        // delete auxiliary array
        delete [] aux;

        // allocate and build index array
        _image_idx = new Index[_num_indices_image];
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          _aux_inj(adj1, adj2, i, &_image_idx[_domain_ptr[i]]);
        }

        // okay
      }

      /// renders transposed adjactor composition
      template<
        typename Adjactor1_,
        typename Adjactor2_>
      void _render_transpose(
        const Adjactor1_& adj1,
        const Adjactor2_& adj2)
      {
        // validate adjactor dimensions
        ASSERT(adj1.get_num_nodes_image() <= adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

        typedef typename Adjactor1_::ImageIterator AImIt1;
        typedef typename Adjactor2_::ImageIterator AImIt2;

        // get counts
        _num_nodes_domain = adj2.get_num_nodes_image();
        _num_nodes_image = adj1.get_num_nodes_domain();
        _num_indices_image = 0;

        // allocate and format pointer array
        _domain_ptr = new Index[_num_nodes_domain + 1];
        for(Index i(0); i <= _num_nodes_domain; ++i)
        {
          _domain_ptr[i] = 0;
        }

        // count number of adjacencies
        for(Index j(0); j < _num_nodes_image; ++j)
        {
          AImIt1 cur1(adj1.image_begin(j));
          AImIt1 end1(adj1.image_end(j));
          for(; cur1 != end1; ++cur1)
          {
            AImIt2 cur2(adj2.image_begin(*cur1));
            AImIt2 end2(adj2.image_end(*cur1));
            for(; cur2 != end2; ++cur2)
            {
              ++_domain_ptr[(*cur2) + 1];
            }
          }
        }

        // build pointer array
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          _domain_ptr[i+1] += _domain_ptr[i];
        }
        _num_indices_image = _domain_ptr[_num_nodes_domain];

        // allocate and build index array
        _image_idx = new Index[_num_indices_image];
        Index** image_ptr = new Index*[_num_nodes_domain];
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          image_ptr[i] = &_image_idx[_domain_ptr[i]];
        }

        for(Index j(0); j < _num_nodes_image; ++j)
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

        // okay
      }

      /// renders transposed injectified adjactor composition
      template<
        typename Adjactor1_,
        typename Adjactor2_>
      void _render_injectify_transpose(
        const Adjactor1_& adj1,
        const Adjactor2_& adj2)
      {
        // validate adjactor dimensions
        ASSERT(adj1.get_num_nodes_image() <= adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

        // get counts
        _num_nodes_domain = adj2.get_num_nodes_image();
        _num_nodes_image = adj1.get_num_nodes_domain();
        _num_indices_image = 0;

        // allocate pointer array
        _domain_ptr = new Index[_num_nodes_domain+1];
        for(Index i(0); i <= _num_nodes_domain; ++i)
        {
          _domain_ptr[i] = 0;
        }

        // allocate auxiliary index array
        Index* aux = new Index[_num_nodes_domain];

        // loop over all image nodes
        for(Index j(0); j < _num_nodes_image; ++j)
        {
          Index num_aux = _aux_inj(adj1, adj2, j, aux);
          for(Index k(0); k < num_aux; ++k)
          {
            ++_domain_ptr[aux[k]+1];
          }
          _num_indices_image += num_aux;
        }

        _image_idx = new Index[_num_indices_image];
        Index** image_ptr = new Index*[_num_nodes_domain];

        // build pointer arrays
        for(Index i(0); i < _num_nodes_domain; ++i)
        {
          _domain_ptr[i+1] += _domain_ptr[i];
          image_ptr[i] = &_image_idx[_domain_ptr[i]];
        }

        // build image index array
        for(Index j(0); j < _num_nodes_image; ++j)
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

        // okay
      }
      /// \endcond
      /* ******************************************************************* */
      /*  A D J A C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
      /* ******************************************************************* */
    public:
      /** \copydoc Adjactor::get_num_nodes_domain() */
      inline Index get_num_nodes_domain() const
      {
        return _num_nodes_domain;
      }

      /** \copydoc Adjactor::get_num_nodes_image() */
      inline Index get_num_nodes_image() const
      {
        return _num_nodes_image;
      }

      /** \copydoc Adjactor::image_begin() */
      inline ImageIterator image_begin(Index domain_node) const
      {
        ASSERT(domain_node < _num_nodes_domain, "Domain node index out of range");
        return &_image_idx[_domain_ptr[domain_node]];
      }

      /** \copydoc Adjactor::image_end() */
      inline ImageIterator image_end(Index domain_node) const
      {
        ASSERT(domain_node < _num_nodes_domain, "Domain node index out of range");

        if(_domain_end != nullptr)
          return &_image_idx[_domain_end[domain_node]];
        else
          return &_image_idx[_domain_ptr[domain_node+1]];
      }
    }; // class Graph
  } // namespace Adjacency
} // namespace FEAT

#endif // KERNEL_ADJACENCY_GRAPH_HPP