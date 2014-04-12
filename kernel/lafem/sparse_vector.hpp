#pragma once
#ifndef KERNEL_LAFEM_SPARSE_VECTOR_HPP
#define KERNEL_LAFEM_SPARSE_VECTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/vector_base.hpp>
#include <kernel/lafem/dense_vector.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Sparse vector class template.
     *
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
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
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class SparseVector : public Container<Mem_, DT_, IT_>, public VectorBase
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
        /// Our indextype
        typedef IT_ IndexType;
        /// Our memory architecture type
        typedef Mem_ MemType;

        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional vector.
         */
        explicit SparseVector() :
          Container<Mem_, DT_, IT_> (0)
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
        explicit SparseVector(Index size_in) :
          Container<Mem_, DT_, IT_>(size_in)
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
         * \param[in] elements A list of non zero elements.
         * \param[in] indices A list of non zero element indices.
         *
         * \note The elements must be sorted by their indices.
         *
         * Creates a vector with a given size.
         */
        explicit SparseVector(Index size_in, DenseVector<Mem_, DT_, IT_> & elements_in, DenseVector<Mem_, IT_, IT_> & indices_in) :
          Container<Mem_, DT_, IT_>(size_in)
        {
          CONTEXT("When creating SparseVector");

          if (indices_in.size() != elements_in.size())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size missmatch!");

          this->_scalar_index.push_back(elements_in.size());
          this->_scalar_index.push_back(elements_in.size());
          this->_scalar_index.push_back(1000);
          this->_scalar_index.push_back(1);
          this->_scalar_dt.push_back(DT_(0));

          this->_elements.push_back(elements_in.elements());
          this->_elements_size.push_back(elements_in.size());
          this->_indices.push_back(indices_in.elements());
          this->_indices_size.push_back(indices_in.size());

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
        }

        /**
         * \brief Move Constructor
         *
         * \param[in] other The source vector.
         *
         * Moves another vector to this vector.
         */
        SparseVector(SparseVector && other) :
          Container<Mem_, DT_, IT_>(std::forward<SparseVector>(other))
        {
          CONTEXT("When moving SparseVector");
        }

        /**
         * \brief Assignment move operator
         *
         * \param[in] other The source vector.
         *
         * Moves another vector to the target vector.
         */
        SparseVector & operator= (SparseVector && other)
        {
          CONTEXT("When moving SparseVector");

          this->move(std::forward<SparseVector>(other));

          return *this;
        }

        /** \brief Clone operation
         *
         * Create a deep copy of itself.
         *
         */
        SparseVector clone() const
        {
          CONTEXT("When cloning SparseVector");
          SparseVector t;
          t.clone(*this, true);
          return t;
        }

        /** \brief Clone operation
         *
         * Become a deep copy of a given vector.
         *
         * \param[in] other The source container.
         * \param[in] clone_indices Should we create a deep copy of the index arrays, too ?
         *
         */
        void clone(const SparseVector & other, bool clone_indices = true)
        {
          CONTEXT("When cloning SparseVector");
          Container<Mem_, DT_, IT_>::clone(other, clone_indices);
        }

        /** \brief Clone operation
         *
         * Become a deep copy of a given vector.
         *
         * \param[in] other The source container.
         * \param[in] clone_indices Should we create a deep copy of the index arrays, too ?
         *
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void clone(const SparseVector<Mem2_, DT2_, IT2_> & other, bool clone_indices = true)
        {
          CONTEXT("When cloning SparseVector");
          SparseVector t;
          t.assign(other);
          Container<Mem_, DT_, IT_>::clone(t, clone_indices);
        }


        /**
         * \brief Convertion method
         *
         * Use source vector content as content of current vector
         *
         * \param[in] other The source container.
         *
         * \note This creates a deep copy in any case!
         *
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void convert(const SparseVector<Mem2_, DT2_, IT2_> & other)
        {
          CONTEXT("When converting SparseVector");
          this->clone(other);
        }

        /**
         * \brief Get a pointer to the data array.
         *
         * \returns Pointer to the data array.
         */
        DT_ * elements()
        {
          if (sorted() == 0)
            const_cast<SparseVector *>(this)->sort();
          return this->_elements.at(0);
        }

        DT_ const * elements() const
        {
          if (sorted() == 0)
            const_cast<SparseVector *>(this)->sort();
          return this->_elements.at(0);
        }

        /**
         * \brief Get a pointer to the non zero indices array.
         *
         * \returns Pointer to the indices array.
         */
        IT_ * indices()
        {
          if (sorted() == 0)
            const_cast<SparseVector *>(this)->sort();
          return this->_indices.at(0);
        }

        IT_ const * indices() const
        {
          if (sorted() == 0)
            const_cast<SparseVector *>(this)->sort();
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
            MemoryPool<Mem_>::instance()->template set_memory<DT_>(this->_elements.back(), DT_(4711), alloc_increment());
            this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(alloc_increment()));
            this->_indices_size.push_back(alloc_increment());
            MemoryPool<Mem_>::instance()->template set_memory<IT_>(this->_indices.back(), IT_(4711), alloc_increment());
            _allocated_elements() = alloc_increment();
            MemoryPool<Mem_>::set_memory(elements(), val);
            MemoryPool<Mem_>::set_memory(indices(), IT_(index));
            _used_elements() = 1;
          }

          // append element in already allocated arrays
          else if(used_elements() < allocated_elements())
          {
            MemoryPool<Mem_>::set_memory(elements() + used_elements(), val);
            MemoryPool<Mem_>::set_memory(indices() + used_elements(), IT_(index));
            ++_used_elements();
          }

          // reallocate arrays, append element
          else
          {
            _allocated_elements() += alloc_increment();

            DT_ * elements_new(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(allocated_elements()));
            MemoryPool<Mem_>::instance()->template set_memory<DT_>(elements_new, DT_(4711), allocated_elements());
            IT_ * indices_new(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(allocated_elements()));
            MemoryPool<Mem_>::instance()->template set_memory<IT_>(indices_new, IT_(4711), allocated_elements());

            MemoryPool<Mem_>::copy(elements_new, elements(), used_elements());
            MemoryPool<Mem_>::copy(indices_new, indices(), used_elements());

            MemoryPool<Mem_>::instance()->release_memory(elements());
            MemoryPool<Mem_>::instance()->release_memory(indices());

            this->_elements.at(0) = elements_new;
            this->_indices.at(0) = indices_new;

            MemoryPool<Mem_>::set_memory(elements() + used_elements(), val);
            MemoryPool<Mem_>::set_memory(indices() + used_elements(), IT_(index));

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
                MemoryPool<Mem_>::set_memory(indices(), std::numeric_limits<IT_>::max());
              }
            }

            // sort out marked duplicated elements
            _insertion_sort(indices(), elements(), used_elements());
            Index junk(0);
            while (MemoryPool<Mem_>::get_element(indices(), used_elements() - 1 - junk) == std::numeric_limits<IT_>::max()
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
            const_cast<SparseVector *>(this)->sort();

          Index i(0);
          while (i < used_elements())
          {
            if (MemoryPool<Mem_>::get_element(indices(), i) >= index)
              break;
            ++i;
          }

          if (i < used_elements() && MemoryPool<Mem_>::get_element(indices(), i) == index)
            return MemoryPool<Mem_>::get_element(elements(), i);
          else
            return zero_element();
        }

        /**
         * \brief Reset all elements of the container to a given value or zero if missing.
         *
         * \param[in] value The value to be set (defaults to 0)
         *
         */
        void format(DT_ value = DT_(0))
        {
          CONTEXT("When clearing SparseVector");

          if (value == DT_(0))
          {
            MemoryPool<Mem_>::instance()->release_memory(elements());
            MemoryPool<Mem_>::instance()->release_memory(indices());

            this->_elements.clear();
            this->_indices.clear();
            this->_elements_size.clear();
            this->_indices_size.clear();
            _used_elements() = 0;
            _allocated_elements() = 0;
          }
          else
          {
            ((Container<Mem_, DT_, IT_> &)*this).format(value);
          }
        }

        /**
         * \brief Retrieve non zero element count.
         *
         * \returns Non zero element count.
         */
        const Index & used_elements() const override
        {
          if (sorted() == 0)
            const_cast<SparseVector *>(this)->sort();
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
        const Index & allocated_elements() const
        {
          return this->_scalar_index.at(2);
        }

        /**
         * \brief Retrieve allocation incrementation value.
         *
         * \return Allocation increment.
         */
        const Index & alloc_increment() const
        {
          return this->_scalar_index.at(3);
        }

        /**
         * \brief Retrieve status of element sorting.
         *
         * \return Sorting status.
         */
        const Index & sorted() const
        {
          return this->_scalar_index.at(4);
        }

        /**
         * \brief Returns a descriptive string.
         *
         * \returns A string describing the container.
         */
        static String name()
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
    template <typename Mem_, typename Mem2_, typename DT_, typename IT_> bool operator== (const SparseVector<Mem_, DT_, IT_> & a, const SparseVector<Mem2_, DT_, IT_> & b)
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
    template <typename Mem_, typename DT_, typename IT_>
    std::ostream &
    operator<< (std::ostream & lhs, const SparseVector<Mem_, DT_, IT_> & b)
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
