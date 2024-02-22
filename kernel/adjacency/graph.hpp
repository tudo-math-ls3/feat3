// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ADJACENCY_GRAPH_HPP
#define KERNEL_ADJACENCY_GRAPH_HPP 1

// includes, FEAT
#include <kernel/adjacency/base.hpp>
#include <kernel/adjacency/adjactor.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/exception.hpp>

// includes, system
#include <vector>

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
      using IndexVector = std::vector < Index>;

      /**
       * \brief ImageIterator for Graph class
       *
       * For the purpose of the Adjactor interface implementation a <c>const Index</c>-pointer is the optimal choice
       * for the image node iterator.
       */
      typedef IndexVector::const_iterator ImageIterator;

      // magic number for Graph serialization
      static constexpr std::uint64_t magic = 0x5052474A44413346ull; // "F3ADJGRP"

    protected:
      /// total number of image nodes
      Index _num_nodes_image;

      /**
       * \brief Domain pointer Vector
       *
       * Dimension: #_num_nodes_domain+1
       */
      IndexVector _domain_ptr;

      /**
       * \brief Image node index Vector
       *
       * Dimension: #_num_indices_image
       */
      IndexVector _image_idx;


    public:

      /**
       * \brief Default constructor.
       *
       * This constructor creates a new empty graph, but does not allocate any vectors.
       */
      Graph();

      /**
       * \brief Allocation Constructor.
       *
       * This constructor creates a new graph and allocates the Graph's vectors to the corresponding lengths.
       *
       * \note This constructor does not initialize the allocated vectors -- they have to be initialized by the user
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
       */
      explicit Graph(
        Index num_nodes_domain,
        Index num_nodes_image,
        Index num_indices_image);

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
       * The \transient domain pointer array for the graph. Must not be \c nullptr.
       *
       * \param[in] image_idx
       * The \transient image node index array for the graph. Must not be \c nullptr.
       */
      explicit Graph(
        Index num_nodes_domain,
        Index num_nodes_image,
        Index num_indices_image,
        const Index* domain_ptr,
        const Index* image_idx);

      /**
       * \brief "Copy-Vector" Constructor
       *
       * This constructor creates a new graph using copies of the vectors passed to this function.
       *
       * \param[in] num_nodes_image
       * The total number of image nodes for the graph.
       *
       * \param[in] domain_ptr
       * The \transient domain pointer vector for the graph. Must not be \c nullptr.
       *
       * \param[in] image_idx
       * The \transient image node index vector for the graph. Must not be \c nullptr.
       */
      explicit Graph(
        Index num_nodes_image,
        const IndexVector& domain_ptr,
        const IndexVector& image_idx);

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
      explicit Graph(RenderType render_type, const Adjactor_& adjactor) :
        _num_nodes_image(0),
        _domain_ptr(),
        _image_idx()
      {
        switch(render_type)
        {
        case RenderType::as_is:
        case RenderType::as_is_sorted:
          _render_as_is(adjactor);
          if(render_type == RenderType::as_is_sorted)
            this->sort_indices();
          break;

        case RenderType::injectify:
        case RenderType::injectify_sorted:
          _render_injectify(adjactor);
          if(render_type == RenderType::injectify_sorted)
            this->sort_indices();
          break;

        case RenderType::transpose:
        case RenderType::transpose_sorted:
          _render_transpose(adjactor);
          // transpose is automatically sorted
          break;

        case RenderType::injectify_transpose:
        case RenderType::injectify_transpose_sorted:
          _render_injectify_transpose(adjactor);
          // transpose is automatically sorted
          break;

        default:
          XABORTM("Invalid render_type parameter!");
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
       * A \transient reference to the first adjactor in the composition that is to be rendered.
       *
       * \param[in] adjactor2
       * A \transient reference to the second adjactor in the composition that is to be rendered.
       */
      template<typename Adjactor1_, typename Adjactor2_>
      explicit Graph(RenderType render_type, const Adjactor1_& adjactor1, const Adjactor2_& adjactor2) :
        _num_nodes_image(0),
        _domain_ptr(),
       _image_idx()
      {
        switch(render_type)
        {
        case RenderType::as_is:
        case RenderType::as_is_sorted:
          _render_as_is(adjactor1, adjactor2);
          if(render_type == RenderType::as_is_sorted)
            this->sort_indices();
          break;

        case RenderType::injectify:
        case RenderType::injectify_sorted:
          _render_injectify(adjactor1, adjactor2);
          if(render_type == RenderType::injectify_sorted)
            this->sort_indices();
          break;

        case RenderType::transpose:
        case RenderType::transpose_sorted:
          _render_transpose(adjactor1, adjactor2);
          // transpose is automatically sorted
          break;

        case RenderType::injectify_transpose:
        case RenderType::injectify_transpose_sorted:
          _render_injectify_transpose(adjactor1, adjactor2);
          // transpose is automatically sorted
          break;

        default:
          XABORTM("Invalid render_type parameter!");
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
       * The \transient graph that has to be permuted.
       *
       * \param[in] domain_perm
       * The \transient permutation of the domain set.
       *
       * \param[in] image_perm
       * The \transient permutation of the image set.
       */
      explicit Graph(const Graph& other, const Permutation& domain_perm, const Permutation& image_perm);

      /**
       * \brief Deserialization CTOR
       *
       * This constructor creates the graph from a buffer that has been created by using
       * the serialize function.
       *
       * \param[in] buffer
       * The \transient serialization buffer to create the graph from.
       */
      explicit Graph(const std::vector<char>& buffer);

      /// virtual destructor
      virtual ~Graph();

      /// Clears the graph
      void clear();

      /**
       * \brief Clones this graph.
       *
       * \returns A deep-copy of this graph.
       */
      Graph clone() const
      {
        if(!_domain_ptr.empty())
          return Graph(_num_nodes_image, _domain_ptr, _image_idx);
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
        ASSERTM(domain_node < get_num_nodes_domain(), "Domain node index out of range");
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
        return _domain_ptr.data();
      }

      /** \copydoc get_domain_ptr() */
      const Index* get_domain_ptr() const
      {
        return _domain_ptr.data();
      }

      /**
       * \brief Returns the image node index array.
       * \returns The image node index array.
       */
      Index* get_image_idx()
      {
        return _image_idx.data();
      }

      /** \copydoc get_image_idx() */
      const Index* get_image_idx() const
      {
        return _image_idx.data();
      }

      /**
       * \brief Returns the total number indices.
       *
       * \returns The total number of indices in the graph.
       */
      Index get_num_indices() const
      {
        return Index(_image_idx.size());
      }

      /**
       * \brief Sorts the image indices to non-descending order.
       */
      void sort_indices();

      /**
       * \brief Permutes the image indices.
       *
       * \param[in] inv_perm
       * The \transient inverse image indices permutation.
       */
      void permute_indices(const Adjacency::Permutation& inv_perm);

      /**
       * \brief Serializes the graph into a buffer.
       */
      std::vector<char> serialize() const;

      /* *************************************************** */
      /*  R E N D E R   F U N C T I O N   T E M P L A T E S  */
      /* *************************************************** */
    private:
      /// \cond internal

      /// renders adjactor
      template<typename Adjactor_>
      void _render_as_is(const Adjactor_& adj)
      {
        typedef typename Adjactor_::ImageIterator AImIt;

        // get counts
        //_num_nodes_domain = adj.get_num_nodes_domain();
        _num_nodes_image = adj.get_num_nodes_image();
        Index num_indices_image = 0;

        // allocate pointer vector
        _domain_ptr = IndexVector(adj.get_num_nodes_domain() +1);

        // count number of adjacencies and build pointer vector
        for(Index i(0); i < adj.get_num_nodes_domain(); ++i)
        {
          _domain_ptr[i] = num_indices_image;
          AImIt cur(adj.image_begin(i));
          AImIt end(adj.image_end(i));
          for(; cur != end; ++cur)
          {
            ++num_indices_image;
          }
        }
        _domain_ptr[adj.get_num_nodes_domain()] = num_indices_image;

        // allocate and build index vector
        _image_idx = IndexVector(num_indices_image);
        for(Index i(0); i < adj.get_num_nodes_domain(); ++i)
        {
          Index* idx = &_image_idx[_domain_ptr[i]];
          AImIt cur(adj.image_begin(i));
          AImIt end(adj.image_end(i));
          for(; cur != end; ++cur, ++idx)
          {
            *idx = *cur;
          }
        }
      }

      /// renders injectified adjactor
      template<typename Adjactor_>
      void _render_injectify(const Adjactor_& adj)
      {
        // get counts
        _num_nodes_image = adj.get_num_nodes_image();
        Index num_indices_image = 0;

        // allocate pointer vector
        _domain_ptr = IndexVector(adj.get_num_nodes_domain() + 1);

        // allocate auxiliary mask vector
        std::vector<char> vidx_mask(adj.get_num_nodes_image(),0);
        char* idx_mask = vidx_mask.data();

        // count number of adjacencies and build pointer vector
        for(Index i(0); i < adj.get_num_nodes_domain(); ++i)
        {
          _domain_ptr[i] = num_indices_image;
          for(auto it = adj.image_begin(i); it != adj.image_end(i); ++it)
          {
            if(idx_mask[*it] == 0)
            {
              ++num_indices_image;
              idx_mask[*it] = 1;
            }
          }
          for(auto it = adj.image_begin(i); it != adj.image_end(i); ++it)
            idx_mask[*it] = 0;
        }
        _domain_ptr[adj.get_num_nodes_domain()] = num_indices_image;

        // allocate and build index vector
        _image_idx = IndexVector(num_indices_image);
        for(Index i(0); i < adj.get_num_nodes_domain(); ++i)
        {
          Index k = _domain_ptr[i];
          for(auto it = adj.image_begin(i); it != adj.image_end(i); ++it)
          {
            if(idx_mask[*it] == 0)
            {
              _image_idx[k] = *it;
              ++k;
              idx_mask[*it] = 1;
            }
          }
          for(auto it = adj.image_begin(i); it != adj.image_end(i); ++it)
            idx_mask[*it] = 0;
        }
      }

      /// renders transposed adjactor
      template<typename Adjactor_>
      void _render_transpose(const Adjactor_& adj)
      {
        typedef typename Adjactor_::ImageIterator AImIt;

        // get counts
        _num_nodes_image = adj.get_num_nodes_domain();
        Index num_indices_image = 0;

        // allocate and format pointer vector
        _domain_ptr = IndexVector(adj.get_num_nodes_image() + 1,Index(0));

        // count number of adjacencies
        for(Index j(0); j < adj.get_num_nodes_domain(); ++j)
        {
          AImIt cur(adj.image_begin(j));
          AImIt end(adj.image_end(j));
          for(; cur != end; ++cur)
          {
            ++_domain_ptr[(*cur) + 1];
          }
        }

        // build pointer vector
        for(Index i(0); i < adj.get_num_nodes_image(); ++i)
        {
          _domain_ptr[i+1] += _domain_ptr[i];
        }
        num_indices_image = _domain_ptr[adj.get_num_nodes_image()];

        // allocate and build index vector
        _image_idx = IndexVector(num_indices_image);
        std::vector<Index*> vimg_ptr(adj.get_num_nodes_image(), nullptr);
        Index** image_ptr = vimg_ptr.data();
        Index* image_idx = _image_idx.data();
        for(Index i(0); i < adj.get_num_nodes_image(); ++i)
        {
          image_ptr[i] = &image_idx[_domain_ptr[i]];
        }

        for(Index j(0); j < adj.get_num_nodes_domain(); ++j)
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
      }

      /// renders transposed injectified adjactor
      template<typename Adjactor_>
      void _render_injectify_transpose(const Adjactor_& adj)
      {
        // get counts
        _num_nodes_image = adj.get_num_nodes_domain();
        Index num_indices_image = 0;

        // allocate pointer vector
        _domain_ptr =IndexVector(adj.get_num_nodes_image() + 1,0);
        // allocate auxiliary mask vector
        std::vector<char> vidx_mask(adj.get_num_nodes_image(), 0);
        char* idx_mask = vidx_mask.data();

        // loop over all image nodes
        for(Index j(0); j < adj.get_num_nodes_domain(); ++j)
        {
          for(auto it = adj.image_begin(j); it != adj.image_end(j); ++it)
          {
            if(idx_mask[*it] == 0)
            {
              ++num_indices_image;
              ++_domain_ptr[(*it)+1];
              idx_mask[*it] = 1;
            }
          }
          for(auto it = adj.image_begin(j); it != adj.image_end(j); ++it)
            idx_mask[*it] = 0;
        }

        _image_idx = IndexVector(num_indices_image);
        std::vector<Index*> vimg_ptr(adj.get_num_nodes_image(), nullptr);
        Index** image_ptr = vimg_ptr.data();
        Index* image_idx = _image_idx.data();
        // build pointer vectors
        for(Index i(0); i <adj.get_num_nodes_image(); ++i)
        {
          _domain_ptr[i+1] += _domain_ptr[i];
          image_ptr[i] = &image_idx[_domain_ptr[i]];
        }

        // build image index vector
        for(Index j(0); j < adj.get_num_nodes_domain(); ++j)
        {
          for(auto it = adj.image_begin(j); it != adj.image_end(j); ++it)
          {
            if(idx_mask[*it] == 0)
            {
              Index*& idx = image_ptr[*it];
              *idx = j;
              ++idx;
              idx_mask[*it] = 1;
            }
          }
          for(auto it = adj.image_begin(j); it != adj.image_end(j); ++it)
            idx_mask[*it] = 0;
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
        // validate adjactor dimensions
        XASSERTM(adj1.get_num_nodes_image() == adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

        typedef typename Adjactor1_::ImageIterator AImIt1;
        typedef typename Adjactor2_::ImageIterator AImIt2;

        // get counts
        _num_nodes_image = adj2.get_num_nodes_image();
        Index num_indices_image = 0;

        // allocate pointer vector
        _domain_ptr = IndexVector(adj1.get_num_nodes_domain() + 1);

        // count number of adjacencies and build pointer vector
        for(Index i(0); i < adj1.get_num_nodes_domain(); ++i)
        {
          _domain_ptr[i] = num_indices_image;
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
        _domain_ptr[adj1.get_num_nodes_domain()] = num_indices_image;

        // allocate and build index vector
        _image_idx =IndexVector(num_indices_image);
        for(Index i(0); i < adj1.get_num_nodes_domain(); ++i)
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
        XASSERTM(adj1.get_num_nodes_image() == adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

        // get counts
        _num_nodes_image = adj2.get_num_nodes_image();
        Index num_indices_image = 0;

        // allocate pointer vector
        _domain_ptr =IndexVector(adj1.get_num_nodes_domain() + 1);

        // allocate auxiliary mask vector
        std::vector<char> vidx_mask(adj2.get_num_nodes_image(), 0);
        char* idx_mask = vidx_mask.data();

        // count number of adjacencies and build pointer vector
        for(Index i(0); i < adj1.get_num_nodes_domain(); ++i)
        {
          _domain_ptr[i] = num_indices_image;
          for(auto it = adj1.image_begin(i); it != adj1.image_end(i); ++it)
          {
            for(auto jt = adj2.image_begin(*it); jt != adj2.image_end(*it); ++jt)
            {
              if(idx_mask[*jt] == 0)
              {
                ++num_indices_image;
                idx_mask[*jt] = 1;
              }
            }
          }
          // reset mask
          for(auto it = adj1.image_begin(i); it != adj1.image_end(i); ++it)
            for(auto jt = adj2.image_begin(*it); jt != adj2.image_end(*it); ++jt)
              idx_mask[*jt] = 0;
        }
        _domain_ptr[adj1.get_num_nodes_domain()] = num_indices_image;

        // allocate and build index vector
        _image_idx = IndexVector(num_indices_image);
        for(Index i(0); i < adj1.get_num_nodes_domain(); ++i)
        {
          Index k = _domain_ptr[i];
          for(auto it = adj1.image_begin(i); it != adj1.image_end(i); ++it)
          {
            for(auto jt = adj2.image_begin(*it); jt != adj2.image_end(*it); ++jt)
            {
              if(idx_mask[*jt] == 0)
              {
                _image_idx[k] = *jt;
                ++k;
                idx_mask[*jt] = 1;
              }
            }
          }
          // reset mask
          for(auto it = adj1.image_begin(i); it != adj1.image_end(i); ++it)
            for(auto jt = adj2.image_begin(*it); jt != adj2.image_end(*it); ++jt)
              idx_mask[*jt] = 0;
        }
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
        XASSERTM(adj1.get_num_nodes_image() == adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

        typedef typename Adjactor1_::ImageIterator AImIt1;
        typedef typename Adjactor2_::ImageIterator AImIt2;

        // get counts
        _num_nodes_image = adj1.get_num_nodes_domain();
        Index num_indices_image = 0;

        // allocate and format pointer vector
        _domain_ptr =IndexVector(adj2.get_num_nodes_image() + 1,Index(0));

        // count number of adjacencies
        for(Index j(0); j < adj1.get_num_nodes_domain(); ++j)
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

        // build pointer vector
        for(Index i(0); i < adj2.get_num_nodes_image(); ++i)
        {
          _domain_ptr[i+1] += _domain_ptr[i];
        }
        num_indices_image = _domain_ptr[adj2.get_num_nodes_image()];

        // allocate and build index vector
        _image_idx =IndexVector(num_indices_image);
        std::vector<Index*> vimg_ptr(adj2.get_num_nodes_image(), nullptr);
        Index** image_ptr = vimg_ptr.data();
        Index* image_idx = _image_idx.data();

        for(Index i(0); i < adj2.get_num_nodes_image(); ++i)
        {
          image_ptr[i] = &image_idx[_domain_ptr[i]];
        }

        for(Index j(0); j < adj1.get_num_nodes_domain(); ++j)
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
        XASSERTM(adj1.get_num_nodes_image() == adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

        // get counts
        _num_nodes_image = adj1.get_num_nodes_domain();
        Index num_indices_image = 0;

        // allocate pointer vector

        _domain_ptr = IndexVector(adj2.get_num_nodes_image() +1,Index(0));

        // allocate auxiliary mask vector
        std::vector<char> vidx_mask(adj2.get_num_nodes_image(), 0);
        char* idx_mask = vidx_mask.data();

        // loop over all image nodes
        for(Index j(0); j < adj1.get_num_nodes_domain(); ++j)
        {
          for(auto it = adj1.image_begin(j); it != adj1.image_end(j); ++it)
          {
            for(auto jt = adj2.image_begin(*it); jt != adj2.image_end(*it); ++jt)
            {
              if(idx_mask[*jt] == 0)
              {
                ++num_indices_image;
                ++_domain_ptr[(*jt)+1];
                idx_mask[*jt] = 1;
              }
            }
          }
          // reset mask
          for(auto it = adj1.image_begin(j); it != adj1.image_end(j); ++it)
            for(auto jt = adj2.image_begin(*it); jt != adj2.image_end(*it); ++jt)
              idx_mask[*jt] = 0;
        }

        _image_idx = IndexVector(num_indices_image);
        std::vector<Index*> vimg_ptr(adj2.get_num_nodes_image(), nullptr);
        Index** image_ptr = vimg_ptr.data();
        Index* image_idx = _image_idx.data();

        // build pointer vector
        for(Index i(0); i < adj2.get_num_nodes_image(); ++i)
        {
          _domain_ptr[i+1] += _domain_ptr[i];
          image_ptr[i] = &image_idx[_domain_ptr[i]];
        }

        // build image index vector
        for(Index j(0); j < adj1.get_num_nodes_domain(); ++j)
        {
          for(auto it = adj1.image_begin(j); it != adj1.image_end(j); ++it)
          {
            for(auto jt = adj2.image_begin(*it); jt != adj2.image_end(*it); ++jt)
            {
              if(idx_mask[*jt] == 0)
              {
                Index*& idx = image_ptr[*jt];
                *idx = j;
                ++idx;
                idx_mask[*jt] = 1;
              }
            }
          }
          // reset mask
          for(auto it = adj1.image_begin(j); it != adj1.image_end(j); ++it)
            for(auto jt = adj2.image_begin(*it); jt != adj2.image_end(*it); ++jt)
              idx_mask[*jt] = 0;
        }
      }

      /// \endcond
      /* ******************************************************************* */
      /*  A D J A C T O R   I N T E R F A C E   I M P L E M E N T A T I O N  */
      /* ******************************************************************* */
    public:


      inline Index get_num_nodes_domain() const
      {
        return (_domain_ptr.empty() ? Index(0) : Index(_domain_ptr.size() - 1));

      }


      inline Index get_num_nodes_image() const
      {
        return _num_nodes_image;
      }


      /** \copydoc Adjactor::image_begin() */
      inline ImageIterator image_begin(Index domain_node) const
      {
        ASSERTM(domain_node +1 < _domain_ptr.size(), "Domain node index out of range");

        return _image_idx.begin() + IndexVector::difference_type(_domain_ptr[domain_node]);
      }

      /** \copydoc Adjactor::image_end() */
      inline ImageIterator image_end(Index domain_node) const
      {
        ASSERTM(domain_node +1< _domain_ptr.size(), "Domain node index out of range");

        return _image_idx.begin() + IndexVector::difference_type(_domain_ptr[domain_node + Index(1)]);
      }
    }; // class Graph
  } // namespace Adjacency
} // namespace FEAT
#endif // KERNEL_ADJACENCY_GRAPH_HPP
