#pragma once
#ifndef KERNEL_LAFEM_SPARSE_LAYOUT_HPP
#define KERNEL_LAFEM_SPARSE_LAYOUT_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/lafem/memory_pool.hpp>

#include <vector>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Layout scheme for sparse matrix containers.
     *
     * \tparam SM_ The Container, the layout belongs to
     * Note, that for technical reasons, all Container types must be typed bool.
     *
     * This class acts as an data wrapper for all index arrays, describing a specific sparse matrix layout.
     * It enables FEAST to store layout related data only once per layout per matrix type.
     * \TODO Enable layout conversion between matrix types
     * In addition, one is able to create a new matrix with a given layout without assembling it a second time.
     *
     * \author Dirk Ribbrock
     */
    template <typename SM_>
    class SparseLayout
    {
    public:
      std::vector<Index*> _indices;
      std::vector<Index> _indices_size;
      std::vector<Index> _scalar_index;

      SparseLayout(const std::vector<Index *> & indices, const std::vector<Index> & indices_size, const std::vector<Index> & scalar_index) :
        _indices(indices),
        _indices_size(indices_size),
        _scalar_index(scalar_index)
      {
        for(auto it(this->_indices.begin()); it != this->_indices.end(); ++it)
          MemoryPool<typename SM_::MemType>::instance()->increase_memory(*it);
      }

      /// copy constructor
      SparseLayout(const SparseLayout& other) :
        _indices(other._indices),
        _indices_size(other._indices_size),
        _scalar_index(other._scalar_index)
      {
        for(auto it(this->_indices.begin()); it != this->_indices.end(); ++it)
          MemoryPool<typename SM_::MemType>::instance()->increase_memory(*it);
      }

      /// virtual destructor
      virtual ~SparseLayout()
      {
        for(auto it(this->_indices.begin()); it != this->_indices.end(); ++it)
          MemoryPool<typename SM_::MemType>::instance()->release_memory(*it);
      }
    };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SPARSE_LAYOUT_HPP
