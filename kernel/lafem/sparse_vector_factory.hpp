// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/sparse_vector.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>

// includes, system
#include <map>
#include <iterator>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Factory for SparseVector and SparseVectorBlocked construction.
     *
     * \tparam DT_
     * The data type to be used for the sparse vector.
     *
     * \tparam IT_
     * The index type to be used for the sparse vector.
     *
     * \author Peter Zajac
     */
    template<typename DT_, typename IT_, int block_size_ = 1>
    class SparseVectorFactory
    {
    public:
      typedef DT_ ScalarDataType;
      typedef Tiny::Vector<DT_, block_size_> BlockedValueType;

    private:
      /// Total size of vector
      Index _size;
      /// entry map
      std::map<IT_, BlockedValueType> _entries;

    public:
      SparseVectorFactory() :
        _size(0)
      {
      }

      explicit SparseVectorFactory(Index size_) :
        _size(size_)
      {
      }

      explicit SparseVectorFactory(const SparseVector<DT_, IT_>& vec) :
        _size(vec.size())
      {
        this->add(vec);
      }

      explicit SparseVectorFactory(const SparseVectorBlocked<DT_, IT_, block_size_>& vec) :
        _size(vec.size())
      {
        this->add(vec);
      }

      /**
       * \brief Adds a new vector entry to the factory.
       *
       * \param[in] i The vector-index of the entry to be added.
       * \param[in] v_i The value of the entry to be added.
       */
      void add(Index i, DT_ v_i)
      {
        XASSERTM(i < _size, "invalid sparse vector index");
        auto ib = _entries.emplace(IT_(i), BlockedValueType(v_i));
        if(!ib.second)
          ib.first->second = BlockedValueType(v_i);
      }

      /**
       * \brief Adds a new vector entry to the factory.
       *
       * \param[in] i The vector-index of the entry to be added.
       * \param[in] v_i The value of the entry to be added.
       */
      void add(Index i, BlockedValueType v_i)
      {
        XASSERTM(i < _size, "invalid sparse vector index");
        auto ib = _entries.emplace(IT_(i), v_i);
        if(!ib.second)
          ib.first->second = v_i;
      }

      void add(const SparseVector<DT_, IT_>& vec)
      {
        XASSERTM(vec.size() == this->size(), "invalid vector size");
        const Memory::TypedView<IT_> v_idx = vec.indices_view_r();
        const Memory::TypedView<DT_> v_val = vec.elements_view_r();
        Index nze = vec.num_nzes();
        for(Index i = 0; i < nze; ++i)
        {
          auto ib = _entries.emplace(v_idx[i], BlockedValueType(v_val[i]));
          if(!ib.second)
            ib.first->second = BlockedValueType(v_val[i]);
        }
      }

      void add(const SparseVectorBlocked<DT_, IT_, block_size_>& vec)
      {
        XASSERTM(vec.size() == this->size(), "invalid vector size");
        const Memory::TypedView<IT_> v_idx = vec.indices_view_r();
        const Memory::TypedView<Tiny::Vector<DT_, block_size_>> v_val = vec.elements_view_r();
        Index nze = vec.num_nzes();
        for(Index i = 0; i < nze; ++i)
        {
          auto ib = _entries.emplace(v_idx[i], v_val[i]);
          if(!ib.second)
            ib.first->second = v_val[i];
        }
      }

      /// Returns the total vector size
      Index size() const
      {
        return _size;
      }

      /// Returns the number of nonzero entries
      Index num_nzes() const
      {
        return Index(_entries.size());
      }

      SparseVector<DT_, IT_> make_sv() const
      {
        SparseVector<DT_, IT_> vec(this->size(), this->num_nzes());
        Memory::TypedView<IT_> v_idx = vec.indices_view_w();
        Memory::TypedView<DT_> v_val = vec.elements_view_w();
        auto it = _entries.begin();
        auto it_end = _entries.end();
        for(Index i = 0; it != it_end; ++i, ++it)
        {
          v_idx[i] = it->first;
          v_val[i] = it->second[0];
        }
        return vec;
      }

      SparseVectorBlocked<DT_, IT_, block_size_> make_svb() const
      {
        SparseVectorBlocked<DT_, IT_, block_size_> vec(this->size(), this->num_nzes());
        Memory::TypedView<IT_> v_idx = vec.indices_view_w();
        Memory::TypedView<Tiny::Vector<DT_, block_size_>> v_val = vec.elements_view_w();
        auto it = _entries.begin();
        auto it_end = _entries.end();
        for(Index i = 0; it != it_end; ++i, ++it)
        {
          v_idx[i] = it->first;
          v_val[i] = it->second;
        }
        return vec;
      }
    }; // class SparseVectorFactory
  } // namespace LAFEM
} // namespace FEAT
