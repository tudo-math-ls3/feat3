// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/adjacency/graph.hpp>
#include <kernel/adjacency/permutation.hpp>

#include <algorithm> // for std::max/min/sort

namespace FEAT
{
  namespace Adjacency
  {
    Graph::Graph() :
      _num_nodes_domain(0),
      _num_nodes_image(0),
      _num_indices_image(0),
      _domain_ptr(nullptr),
      _image_idx(nullptr),
      _shared(false)
    {
    }

    // allocation constructor
    Graph::Graph(
      Index num_nodes_domain,
      Index num_nodes_image,
      Index num_indices_image)
        :
      _num_nodes_domain(num_nodes_domain),
      _num_nodes_image(num_nodes_image),
      _num_indices_image(num_indices_image),
      _domain_ptr(nullptr),
      _image_idx(nullptr),
      _shared(false)
    {
      _domain_ptr = new Index[_num_nodes_domain+1];
      _image_idx = new Index[_num_indices_image];
    }

    // "Using-Arrays" Constructor
    Graph::Graph(
      Index num_nodes_domain,
      Index num_nodes_image,
      Index num_indices_image,
      Index* domain_ptr,
      Index* image_idx,
      bool shared)
        :
      _num_nodes_domain(num_nodes_domain),
      _num_nodes_image(num_nodes_image),
      _num_indices_image(num_indices_image),
      _domain_ptr(domain_ptr),
      _image_idx(image_idx),
      _shared(shared)
    {
    }

    // "Copy-Arrays" Constructor
    Graph::Graph(
      Index num_nodes_domain,
      Index num_nodes_image,
      Index num_indices_image,
      const Index* domain_ptr,
      const Index* image_idx)
        :
      _num_nodes_domain(num_nodes_domain),
      _num_nodes_image(num_nodes_image),
      _num_indices_image(num_indices_image),
      _domain_ptr(nullptr),
      _image_idx(nullptr),
      _shared(false)
    {
      _domain_ptr = new Index[num_nodes_domain+1];
      for(Index i(0); i <= num_nodes_domain; ++i)
      {
        _domain_ptr[i] = domain_ptr[i];
      }
      _image_idx = new Index[num_indices_image];
      for(Index i(0); i < num_indices_image; ++i)
      {
        _image_idx[i] = image_idx[i];
      }
    }

    // move ctor
    Graph::Graph(Graph&& other) :
      _num_nodes_domain(other._num_nodes_domain),
      _num_nodes_image(other._num_nodes_image),
      _num_indices_image(other._num_indices_image),
      _domain_ptr(other._domain_ptr),
      _image_idx(other._image_idx),
      _shared(other._shared)
    {
      other._num_nodes_domain = other._num_nodes_image = other._num_indices_image = Index(0);
      other._domain_ptr = nullptr;
      other._image_idx = nullptr;
      other._shared = false;
    }

    /// move assignment
    Graph& Graph::operator=(Graph&& other)
    {
      // avoid self-move
      if(this == &other)
        return *this;

      if(!_shared)
      {
        if(_image_idx != nullptr)
          delete [] _image_idx;
        if(_domain_ptr != nullptr)
          delete [] _domain_ptr;
      }

      _num_nodes_domain = other._num_nodes_domain;
      _num_nodes_image = other._num_nodes_image;
      _num_indices_image = other._num_indices_image;
      _domain_ptr = other._domain_ptr;
      _image_idx = other._image_idx;
      _shared = other._shared;

      other._num_nodes_domain = other._num_nodes_image = other._num_indices_image = Index(0);
      other._domain_ptr = nullptr;
      other._image_idx = nullptr;
      other._shared = false;

      return *this;
    }

    // "Permutation" copy CTOR
    Graph::Graph(const Graph& other, const Permutation& domain_perm, const Permutation& image_perm) :
      _num_nodes_domain(other.get_num_nodes_domain()),
      _num_nodes_image(other.get_num_nodes_image()),
      _num_indices_image(other.get_num_indices()),
      _domain_ptr(new Index[_num_nodes_domain+1]),
      _image_idx(new Index[_num_indices_image]),
      _shared(false)
    {
      // get domain permutation
      const Index* domain_perm_pos = domain_perm.get_perm_pos();

      // calculate new domain array
      _domain_ptr[0] = other._domain_ptr[0];
      for(Index i(0); i < _num_nodes_domain; ++i)
      {
        _domain_ptr[i+1] = other._domain_ptr[domain_perm_pos[i]+1]
                        - other._domain_ptr[domain_perm_pos[i]]
                        + _domain_ptr[i];
      }

      // get image permutation
      const Index* image_perm_pos = image_perm.get_perm_pos();

      // calculate new image array
      Index count = 0;
      for(Index i(0); i < _num_nodes_domain; ++i)
      {
        Index _offset  = other._domain_ptr[domain_perm_pos[i]];
        Index _row_end = other._domain_ptr[domain_perm_pos[i]+1];
        for(Index j(_offset); j < _row_end; ++j)
        {
          _image_idx[count] = image_perm_pos[other._image_idx[j]];
          ++count;
        }
      }
    }

    Graph::Graph(const std::vector<char>& buffer) :
      _num_nodes_domain(0),
      _num_nodes_image(0),
      _num_indices_image(0),
      _domain_ptr(nullptr),
      _image_idx(nullptr),
      _shared(false)
    {
      // 40 bytes required by header only
      XASSERTM(buffer.size() >= std::size_t(40u), "invalid buffer size");

      typedef std::uint64_t u64;
      const u64* v = reinterpret_cast<const u64*>(buffer.data());

      // check magic
      XASSERTM(v[0] == Graph::magic, "invalid magic number");

      // check buffer size
      XASSERTM(v[1] == u64(buffer.size()), "invalid buffer size");

      this->_num_nodes_domain  = Index(v[2]);
      this->_num_nodes_image   = Index(v[3]);
      this->_num_indices_image = Index(v[4]);

      const u64* x = &v[5];

      // allocate and copy domain pointer array
      if(this->_num_nodes_domain > Index(0))
      {
        this->_domain_ptr = new Index[this->_num_nodes_domain+1u];
        for(Index i(0); i <= this->_num_nodes_domain; ++i)
          this->_domain_ptr[i] = Index(x[i]);
        x += std::streamoff(this->_num_nodes_domain+1u);
      }

      // allocate and copy image index array
      if(this->_num_indices_image > Index(0))
      {
        this->_image_idx =  new Index[this->_num_indices_image];
        for(Index i(0); i < this->_num_indices_image; ++i)
          this->_image_idx[i] = Index(x[i]);
      }
    }

    // destructor
    Graph::~Graph()
    {
      if(!_shared)
      {
        if(_image_idx != nullptr)
          delete [] _image_idx;
        if(_domain_ptr != nullptr)
          delete [] _domain_ptr;
      }
    }

    Index Graph::degree() const
    {
      Index deg = 0;
      for(Index i(0); i < _num_nodes_domain; ++i)
      {
        deg = std::max(deg, _domain_ptr[i+1] - _domain_ptr[i]);
      }
      return deg;
    }

    void Graph::sort_indices()
    {
      XASSERTM(_domain_ptr != nullptr, "domain pointer array is missing");
      XASSERTM(_image_idx != nullptr, "image index array is missing");

      // loop over all domain nodes
      for(Index i(0); i < _num_nodes_domain; ++i)
      {
        // let the STL do the sorting work
        std::sort(&_image_idx[_domain_ptr[i]], &_image_idx[_domain_ptr[i+1]]);
      }
    }

    std::vector<char> Graph::serialize() const
    {
      // compute buffer size
      typedef std::uint64_t u64;
      u64 s = u64(_num_nodes_domain > 0u ? 6u : 5u) + u64(_num_nodes_domain) + u64(_num_indices_image);

      // allocate buffer
      std::vector<char> buffer(s * 8u);
      u64* v = reinterpret_cast<u64*>(buffer.data());

      // create header
      v[0] = Graph::magic;
      v[1] = s * u64(8); // buffer size
      v[2] = u64(this->_num_nodes_domain);
      v[3] = u64(this->_num_nodes_image);
      v[4] = u64(this->_num_indices_image);

      // empty graph?
      if(this->_num_nodes_domain <= Index(0))
        return buffer;

      XASSERT(this->_domain_ptr != nullptr);

      // copy domain pointer
      u64* x = &v[5];
      for(Index i(0); i <= this->_num_nodes_domain; ++i)
        x[i] = u64(this->_domain_ptr[i]);

      // copy image indices
      x += std::streamoff(this->_num_nodes_domain+1u);
      for(Index i(0); i < this->_num_indices_image; ++i)
        x[i] = u64(this->_image_idx[i]);

      // okay
      return buffer;
    }
  } // namespace Adjacency
} // namespace FEAT
