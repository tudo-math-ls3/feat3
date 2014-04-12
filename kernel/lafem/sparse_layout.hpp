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
     * Supported File modes.
     */
    enum class SparseLayoutId
    {
      lt_csr = 0,
      lt_coo = 1,
      lt_ell
    };

    namespace Intern
    {
      template <SparseLayoutId LT_>
      struct LayoutId;

      template <>
      struct LayoutId<SparseLayoutId::lt_csr>
      {
        template<typename Mem_, typename DT_, typename IT_>
        using MatrixType = SparseMatrixCSR<Mem_, DT_, IT_>;
      };

      template <>
      struct LayoutId<SparseLayoutId::lt_ell>
      {
        template<typename Mem_, typename DT_, typename IT_>
        using MatrixType = SparseMatrixELL<Mem_, DT_, IT_>;
      };

      template <>
      struct LayoutId<SparseLayoutId::lt_coo>
      {
        template<typename Mem_, typename DT_, typename IT_>
        using MatrixType = SparseMatrixCOO<Mem_, DT_, IT_>;
      };
    }

    /**
     * \brief Layout scheme for sparse matrix containers.
     *
     * \tparam Mem_ The memory where the layout data is lying.
     * \tparam Layout_ The Matrix Type, which represented by the layout.
     *
     * This class acts as an data wrapper for all index arrays, describing a specific sparse matrix layout.
     * It enables FEAST to store layout related data only once per layout per matrix type.
     * \TODO Enable layout conversion between matrix types
     * In addition, one is able to create a new matrix with a given layout without assembling it a second time.
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename IT_, SparseLayoutId Layout_>
    class SparseLayout
    {
    public:
      std::vector<IT_*> _indices;
      std::vector<Index> _indices_size;
      std::vector<Index> _scalar_index;

      template<typename DT_>
      using MatrixType = typename Intern::LayoutId<Layout_>::template MatrixType<Mem_, DT_, IT_>;

      SparseLayout(const std::vector<IT_ *> & indices, const std::vector<Index> & indices_size, const std::vector<Index> & scalar_index) :
        _indices(indices),
        _indices_size(indices_size),
        _scalar_index(scalar_index)
      {
        for(auto i : this->_indices)
          MemoryPool<Mem_>::instance()->increase_memory(i);
      }

      /// copy constructor
      SparseLayout(const SparseLayout & other) :
        _indices(other._indices),
        _indices_size(other._indices_size),
        _scalar_index(other._scalar_index)
      {
        for(auto i : this->_indices)
          MemoryPool<Mem_>::instance()->increase_memory(i);
      }

      /// move constructor
      SparseLayout(SparseLayout && other) :
        _indices(std::move(other._indices)),
        _indices_size(std::move(other._indices_size)),
        _scalar_index(std::move(other._scalar_index))
      {
      }

      /// virtual destructor
      virtual ~SparseLayout()
      {
        for(auto i : this->_indices)
          MemoryPool<Mem_>::instance()->release_memory(i);
      }

      /// operator=
      SparseLayout & operator= (const SparseLayout & other)
      {
        if(this == &other)
          return *this;

        _indices = other._indices;
        _indices_size = other._indices_size;
        _scalar_index = other._scalar_index;

        for(auto i : this->_indices)
          MemoryPool<Mem_>::instance()->increase_memory(i);

        return *this;
      }

      /// move operator=
      SparseLayout & operator= (SparseLayout && other)
      {
        _indices = std::move(other._indices);
        _indices_size = std::move(other._indices_size);
        _scalar_index = std::move(other._scalar_index);

        return *this;
      }

    /**
    * \brief Returns a list of all Index arrays.
    *
    * \returns A list of all Index arrays.
    */
    const std::vector<IT_*> & get_indices() const
    {
      return _indices;
    }

    /**
    * \brief Returns a list of all Index array sizes.
    *
    * \returns A list of all Index array sizes.
    */
    const std::vector<Index> & get_indices_size() const
    {
      return _indices_size;
    }

    /**
    * \brief Returns a list of all scalar values with datatype index.
    *
    * \returns A list of all scalars with datatype index.
    */
    const std::vector<Index> & get_scalar_index() const
    {
      return _scalar_index;
    }

  };
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SPARSE_LAYOUT_HPP
