// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
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
      _num_nodes_image(0),
      _domain_ptr(),
      _image_idx()
    {}

    // allocation constructor
    Graph::Graph(
      Index num_nodes_domain,
      Index num_nodes_image,
      Index num_indices_image)
        :
      _num_nodes_image(num_nodes_image),
      _domain_ptr(num_nodes_domain + 1),
      _image_idx(num_indices_image)
    {}

    // "Copy-Array" Constructor
    Graph::Graph(
      Index num_nodes_domain,
      Index num_nodes_image,
      Index num_indices_image,
      const Index* domain_ptr,
      const Index* image_idx)
        :
      _num_nodes_image(num_nodes_image),
      _domain_ptr(num_nodes_domain + 1),
      _image_idx(num_indices_image)
    {
      for(Index i(0); i <= num_nodes_domain; ++i)
      {
        _domain_ptr[i] = domain_ptr[i];
      }
      for(Index i(0); i < num_indices_image; ++i)
      {
        _image_idx[i] = image_idx[i];
      }
    }
    // "Copy-Vector" Constructor
    Graph::Graph(
      Index num_nodes_image,
      const IndexVector& domain_ptr,
      const IndexVector& image_idx)
      :
      _num_nodes_image(num_nodes_image),
      _domain_ptr(domain_ptr),
      _image_idx(image_idx)

    {}


    // move ctor
    Graph::Graph(Graph&& other) :
      _num_nodes_image(other._num_nodes_image),
      _domain_ptr(std::forward<IndexVector>(other._domain_ptr)),
      _image_idx(std::forward<IndexVector>(other._image_idx))
    {
      other._num_nodes_image = Index(0);
      other._domain_ptr.clear();
      other._image_idx.clear();
    }

    /// move assignment
    Graph& Graph::operator=(Graph&& other)
    {
      // avoid self-move
      if(this == &other)
        return *this;

      _num_nodes_image = other._num_nodes_image;
      _domain_ptr = std::forward<IndexVector>(other._domain_ptr);
      _image_idx = std::forward<IndexVector>(other._image_idx);

      other._num_nodes_image = Index(0);
      other._domain_ptr.clear();
      other._image_idx.clear();

      return *this;
    }

    // "Permutation" copy CTOR
    Graph::Graph(const Graph& other, const Permutation& domain_perm, const Permutation& image_perm) :
      _num_nodes_image(other.get_num_nodes_image()),
      _domain_ptr(other.get_num_nodes_domain()+1),
      _image_idx(other.get_num_indices())
    {
      // get domain permutation
      const Index* domain_perm_pos = domain_perm.get_perm_pos();

      // calculate new domain vector
      _domain_ptr[0] = other._domain_ptr[0];
      for(Index i(0); i < other.get_num_nodes_domain(); ++i)
      {
        _domain_ptr[i+1] = other._domain_ptr[domain_perm_pos[i]+1]
                        - other._domain_ptr[domain_perm_pos[i]]
                        + _domain_ptr[i];
      }

      // get image permutation
      const Index* image_perm_pos = image_perm.get_perm_pos();

      // calculate new image vector
      Index count = 0;
      for(Index i(0); i < other.get_num_nodes_domain(); ++i)
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
      _num_nodes_image(0),
      _domain_ptr(),
      _image_idx()
    {

      // 40 bytes required by header only
      XASSERTM(buffer.size() >= std::size_t(40u), "invalid buffer size");

      typedef std::uint64_t u64;
      const u64* v = reinterpret_cast<const u64*>(buffer.data());

      // check magic
      XASSERTM(v[0] == Graph::magic, "invalid magic number");

      // check buffer size
      XASSERTM(v[1] == u64(buffer.size()), "invalid buffer size");

     Index domain_ptr_size  = Index(v[2]);
     this->_num_nodes_image = Index(v[3]);
     Index num_indices_image = Index(v[4]);

      const u64* x = &v[5];

      // allocate and copy domain pointer vector
      if(domain_ptr_size > Index(0))
      {
        this->_domain_ptr = IndexVector(domain_ptr_size+1u);
        for(Index i(0); i <= this->get_num_nodes_domain(); ++i)
          this->_domain_ptr[i] = Index(x[i]);
        x += std::streamoff(this->get_num_nodes_domain() +1u);
      }

      // allocate and copy image index vector
      if(num_indices_image > Index(0))
      {
        this->_image_idx = IndexVector(num_indices_image);
        for(Index i(0); i < num_indices_image; ++i)
          this->_image_idx[i] = Index(x[i]);
      }
    }

    // destructor
    Graph::~Graph()
    {}

    Index Graph::degree() const
    {
      Index deg = 0;
      for(Index i(0); i+1 < _domain_ptr.size(); ++i)
      {
        deg = std::max(deg, _domain_ptr[i+1] - _domain_ptr[i]);
      }
      return deg;
    }

    void Graph::sort_indices()
    {
      XASSERTM(!_domain_ptr.empty(), "domain pointer vector is missing");
      XASSERTM(!_image_idx.empty(), "image index vector is missing");

      // loop over all domain nodes
      for(Index i(0); i+1 < _domain_ptr.size(); ++i)
      {
        // let the STL do the sorting work
        std::sort(_image_idx.begin() + IndexVector::difference_type(_domain_ptr[i]),
          _image_idx.begin() + IndexVector::difference_type(_domain_ptr[i+1]));
      }
    }

    std::vector<char> Graph::serialize() const
    {
      // compute buffer size
      typedef std::uint64_t u64;
      u64 s = 5u + u64(_domain_ptr.size()) + u64(_image_idx.size());
      // allocate buffer
      std::vector<char> buffer(s * 8u);
      u64* v = reinterpret_cast<u64*>(buffer.data());

      // create header
      v[0] = Graph::magic;
      v[1] = s * u64(8); // buffer size
      v[2] = u64(this->_domain_ptr.size()-1);
      v[3] = u64(this->_num_nodes_image);
      v[4] = u64(this->_image_idx.size());


      // empty graph?
      if(this->_domain_ptr.size() <= Index(1))
        return buffer;

      XASSERT(!(this->_domain_ptr.empty()));

      // copy domain pointer
      u64* x = &v[5];
      for(Index i(0); i < this->_domain_ptr.size(); ++i)
        x[i] = u64(this->_domain_ptr[i]);

      // copy image indices
      x += std::streamoff(this->_domain_ptr.size());
      for(Index i(0); i < this->_image_idx.size(); ++i)
        x[i] = u64(this->_image_idx[i]);

      // okay
      return buffer;
    }
  } // namespace Adjacency
} // namespace FEAT
