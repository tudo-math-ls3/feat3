#pragma once
#ifndef KERNEL_LAFEM_SPARSE_VECTOR_HPP
#define KERNEL_LAFEM_SPARSE_VECTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/vector_base.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Sparse vector class template.
     *
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a vector with non zero elements in a sparse layout. \n \n
     * Data survey: \n
     * _elements[0]: raw number values \n
     * _indices[0]: non zero indices \n
     * _scalar_index[0]: container size
     * _scalar_index[1]: non zero element count (used elements) \n
     * _scalar_index[2]: allocated elements \n
     * _scalar_index[3]: allocation size increment \n
     * _scalar_index[4]: boolean flag, if container is sorted \n
     * _scalar_dt[0]: zero element
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_>
    class SparseVector : public Container<Mem_, DT_>, public VectorBase
    {
      private:
        template <typename T1_, typename T2_>
        static void _insertion_sort(T1_ * key, T2_ * val1, Index size)
        {
          T1_ swap_key;
          T2_ swap1;
          for (Index i(1), j ; i < size ; ++i)
          {
            swap_key = MemoryPool<Mem_>::get_element(key, i);
            swap1 = MemoryPool<Mem_>::get_element(val1, i);
            j = i;
            while (j > 0 && MemoryPool<Mem_>::get_element(key, j - 1) > swap_key)
            {
              MemoryPool<Mem_>::copy(key + j, key + j - 1, 1);
              MemoryPool<Mem_>::copy(val1 + j, val1 + j - 1, 1);
              --j;
            }
            MemoryPool<Mem_>::set_memory(key + j, swap_key);
            MemoryPool<Mem_>::set_memory(val1 + j, swap1);
          }
        }

        Index & _used_elements()
        {
          return this->_scalar_index.at(1);
        }

        Index & _allocated_elements()
        {
          return this->_scalar_index.at(2);
        }

        Index & _sorted()
        {
          return this->_scalar_index.at(4);
        }

      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our memory architecture type
        typedef Mem_ MemType;

        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional vector.
         */
        explicit SparseVector() :
          Container<Mem_, DT_> (0)
        {
          CONTEXT("When creating SparseVector");
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
          this->_scalar_dt.push_back(DT_(0));
        }

        /**
         * \brief Constructor
         *
         * \param[in] size The size of the created vector.
         *
         * Creates a vector with a given size.
         */
        explicit SparseVector(Index size) :
          Container<Mem_, DT_>(size)
        {
          CONTEXT("When creating SparseVector");
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(0);
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
          this->_scalar_dt.push_back(DT_(0));
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source vector.
         *
         * Creates a shallow copy of a given vector.
         */
        SparseVector(const SparseVector<Mem_, DT_> & other) :
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying SparseVector");
        }

        /**
         * \brief Move Constructor
         *
         * \param[in] other The source vector.
         *
         * Moves another vector to this vector.
         */
        SparseVector(SparseVector<Mem_, DT_> && other) :
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When moving SparseVector");
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source vector.
         *
         * Creates a copy of a given vector from another memory architecture.
         */
        template <typename Mem2_, typename DT2_>
        explicit SparseVector(const SparseVector<Mem2_, DT2_> & other) :
            Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying SparseVector");
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this vector.
         */
        SparseVector<Mem_, DT_> clone() const
        {
          CONTEXT("When cloning SparseVector");

          SparseVector<Mem_, DT_> t;
          ((Container<Mem_, DT_>&)t).clone(*this);
          return t;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source vector.
         *
         * Assigns another vector to the target vector.
         */
        SparseVector<Mem_, DT_> & operator= (const SparseVector<Mem_, DT_> & other)
        {
          CONTEXT("When assigning SparseVector");

          this->assign(other);

          return *this;
        }

        /**
         * \brief Assignment move operator
         *
         * \param[in] other The source vector.
         *
         * Moves another vector to the target vector.
         */
        SparseVector<Mem_, DT_> & operator= (SparseVector<Mem_, DT_> && other)
        {
          CONTEXT("When moving SparseVector");

          this->move(std::move(other));

          return *this;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source vector.
         *
         * Assigns a vector from another memory architecture to the target vector.
         */
        template <typename Mem2_, typename DT2_>
        SparseVector<Mem_, DT_> & operator= (const SparseVector<Mem2_, DT2_> & other)
        {
          CONTEXT("When assigning SparseVector");

         this->assign(other);

          return *this;
        }

        /**
         * \brief Get a pointer to the data array.
         *
         * \returns Pointer to the data array.
         */
        DT_ * elements()
        {
          if (sorted() == 0)
            const_cast<SparseVector<Mem_, DT_> *>(this)->sort();
          return this->_elements.at(0);
        }

        DT_ const * elements() const
        {
          if (sorted() == 0)
            const_cast<SparseVector<Mem_, DT_> *>(this)->sort();
          return this->_elements.at(0);
        }

        /**
         * \brief Get a pointer to the non zero indices array.
         *
         * \returns Pointer to the indices array.
         */
        Index * indices()
        {
          if (sorted() == 0)
            const_cast<SparseVector<Mem_, DT_> *>(this)->sort();
          return this->_indices.at(0);
        }

        Index const * indices() const
        {
          if (sorted() == 0)
            const_cast<SparseVector<Mem_, DT_> *>(this)->sort();
          return this->_indices.at(0);
        }

        /**
         * \brief Set specific vector element.
         *
         * \param[in] index The index of the vector element.
         * \param[in] val The val to be set.
         */
        void operator()(Index index, DT_ val)
        {
          CONTEXT("When setting SparseVector element");

          ASSERT(index < this->_scalar_index.at(0), "Error: " + stringify(index) + " exceeds sparse vector size " + stringify(this->_scalar_index.at(0)) + " !");

          _sorted() = 0;

          // vector is empty, no arrays allocated
          if (this->_elements.size() == 0)
          {
            this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(alloc_increment()));
            this->_elements_size.push_back(alloc_increment());
            this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(alloc_increment()));
            this->_indices_size.push_back(alloc_increment());
            _allocated_elements() = alloc_increment();
            MemoryPool<Mem_>::set_memory(elements(), val);
            MemoryPool<Mem_>::set_memory(indices(), index);
            _used_elements() = 1;
          }

          // append element in already allocated arrays
          else if(used_elements() < allocated_elements())
          {
            MemoryPool<Mem_>::set_memory(elements() + used_elements(), val);
            MemoryPool<Mem_>::set_memory(indices() + used_elements(), index);
            ++_used_elements();
          }

          // reallocate arrays, append element
          else
          {
            _allocated_elements() += alloc_increment();

            DT_ * elements_new(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(allocated_elements()));
            Index * indices_new(MemoryPool<Mem_>::instance()->template allocate_memory<Index>(allocated_elements()));

            MemoryPool<Mem_>::copy(elements_new, elements(), used_elements());
            MemoryPool<Mem_>::copy(indices_new, indices(), used_elements());

            MemoryPool<Mem_>::instance()->release_memory(elements());
            MemoryPool<Mem_>::instance()->release_memory(indices());

            this->_elements.at(0) = elements_new;
            this->_indices.at(0) = indices_new;

            MemoryPool<Mem_>::set_memory(elements() + used_elements(), val);
            MemoryPool<Mem_>::set_memory(indices() + used_elements(), index);

            ++_used_elements();
            this->_elements_size.at(0) = allocated_elements();
            this->_indices_size.at(0) = allocated_elements();
          }
        }

        void sort()
        {
          if (sorted() == 0)
          {
            //first of all, mark vector as sorted, because otherwise we would call ourselves inifite times
            _sorted() = 1;

            // check if there is anything to be sorted
            if(used_elements() <= Index(0))
              return;

            _insertion_sort(indices(), elements(), used_elements());

            // find and mark duplicate entries
            for (Index i(1) ; i < used_elements() ; ++i)
            {
              if (MemoryPool<Mem_>::get_element(indices(), i - 1) == MemoryPool<Mem_>::get_element(indices(), i))
              {
                MemoryPool<Mem_>::set_memory(indices(), std::numeric_limits<Index>::max());
              }
            }

            // sort out marked duplicated elements
            _insertion_sort(indices(), elements(), used_elements());
            Index junk(0);
            while (MemoryPool<Mem_>::get_element(indices(), used_elements() - 1 - junk) == std::numeric_limits<Index>::max()
                && junk < used_elements())
              ++junk;
            _used_elements() -= junk;
          }
        }

        /**
         * \brief Retrieve specific vector element.
         *
         * \param[in] index The index of the vector element.
         *
         * \returns Specific vector element.
         */
        DT_ operator()(Index index) const
        {
          CONTEXT("When retrieving SparseVector element");

          ASSERT(index < this->_scalar_index.at(0), "Error: " + stringify(index) + " exceeds sparse vector size " + stringify(this->_scalar_index.at(0)) + " !");

          if (this->_elements.size() == 0)
            return zero_element();

          if (sorted() == 0)
            const_cast<SparseVector<Mem_, DT_> *>(this)->sort();

          Index i(0);
          while (i < used_elements())
          {
            if (MemoryPool<Mem_>::get_element(indices(), i) >= index)
              break;
            ++i;
          }

          if (MemoryPool<Mem_>::get_element(indices(), i) == index)
            return MemoryPool<Mem_>::get_element(elements(), i);
          else
            return zero_element();
        }

        /**
         * \brief Reset all elements to zero.
         */
        void clear()
        {
          CONTEXT("When clearing SparseVector");
          MemoryPool<Mem_>::instance()->release_memory(elements());
          MemoryPool<Mem_>::instance()->release_memory(indices());

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();
          _used_elements() = 0;
          _allocated_elements() = 0;
        }

        /**
         * \brief Retrieve non zero element count.
         *
         * \returns Non zero element count.
         */
        Index used_elements() const
        {
          if (sorted() == 0)
            const_cast<SparseVector<Mem_, DT_> *>(this)->sort();
          return this->_scalar_index.at(1);
        }

        /**
         * \brief Retrieve non zero element.
         *
         * \returns Non zero element.
         */
        DT_ zero_element() const
        {
          return this->_scalar_dt.at(0);
        }

        /**
         * \brief Retrieve amount of allocated elements.
         *
         * \return Allocated element count.
         */
        Index allocated_elements() const
        {
          return this->_scalar_index.at(2);
        }

        /**
         * \brief Retrieve allocation incrementation value.
         *
         * \return Allocation increment.
         */
        Index alloc_increment() const
        {
          return this->_scalar_index.at(3);
        }

        /**
         * \brief Retrieve status of element sorting.
         *
         * \return Sorting status.
         */
        Index sorted() const
        {
          return this->_scalar_index.at(4);
        }

        /**
         * \brief Returns a descriptive string.
         *
         * \returns A string describing the container.
         */
        static String type_name()
        {
          return "SparseVector";
        }
    }; // class SparseVector<...>


    /**
     * \brief SparseVector comparison operator
     *
     * \param[in] a A vector to compare with.
     * \param[in] b A vector to compare with.
     */
    template <typename Mem_, typename Mem2_, typename DT_> bool operator== (const SparseVector<Mem_, DT_> & a, const SparseVector<Mem2_, DT_> & b)
    {
      CONTEXT("When comparing SparseVectors");

      if (a.size() != b.size())
        return false;
      if (a.get_elements().size() != b.get_elements().size())
        return false;
      if (a.get_indices().size() != b.get_indices().size())
        return false;

      for (Index i(0) ; i < a.size() ; ++i)
        if (a(i) != b(i))
          return false;

      return true;
    }

    /**
     * \brief SparseVector streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The vector to be streamed.
     */
    template <typename Mem_, typename DT_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseVector<Mem_, DT_> & b)
    {
      lhs << "[";
      for (Index i(0) ; i < b.size() ; ++i)
      {
        lhs << "  " << b(i);
      }
      lhs << "]";

      return lhs;
    }

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_DENSE_VECTOR_HPP
