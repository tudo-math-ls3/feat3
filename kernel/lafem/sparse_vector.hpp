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
#include <kernel/util/math.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Sparse vector class template.
     *
     * \tparam Mem_ The \ref FEAST::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     *
     * This class represents a vector with non zero elements in a sparse layout. \n \n
     * Data survey: \n
     * _elements[0]: raw number values \n
     * _indices[0]: non zero indices \n
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: non zero element count (used elements) \n
     * _scalar_index[2]: allocated elements \n
     * _scalar_index[3]: allocation size increment \n
     * _scalar_index[4]: boolean flag, if container is sorted \n
     * _scalar_dt[0]: zero element
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class SparseVector : public Container<Mem_, DT_, IT_>, public VectorBase
    {
      private:
        void _read_from_mtx(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in);
          if (! file.is_open())
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
          _read_from_mtx(file);
          file.close();
        }

        void _read_from_mtx(std::istream& file)
        {
          Index rows;
          Index nnz;
          String line;
          std::getline(file, line);
          if (line.find("%%MatrixMarket matrix coordinate real general") == String::npos)
          {
            throw InternalError(__func__, __FILE__, __LINE__, "Input-file is not a compatible mtx-file");
          }
          while(!file.eof())
          {
            std::getline(file,line);
            if (file.eof())
              throw InternalError(__func__, __FILE__, __LINE__, "Input-file is empty");

            String::size_type begin(line.find_first_not_of(" "));
            if (line.at(begin) != '%')
              break;
          }
          {
            String::size_type begin(line.find_first_not_of(" "));
            line.erase(0, begin);
            String::size_type end(line.find_first_of(" "));
            String srows(line, 0, end);
            rows = (Index)atol(srows.c_str());
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String scols(line, 0, end);
            Index cols((Index)atol(scols.c_str()));
            line.erase(0, end);
            if (cols != 1)
              throw InternalError(__func__, __FILE__, __LINE__, "Input-file is no sparse-vector-file");

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String snnz(line, 0, end);
            nnz = (Index)atol(snnz.c_str());
            line.erase(0, end);
          }

          DenseVector<Mem::Main, IT_, IT_> ind(nnz);
          DenseVector<Mem::Main, DT_, IT_> val(nnz);

          IT_ * pind(ind.elements());
          DT_ * pval(val.elements());

          while(!file.eof())
          {
            std::getline(file, line);
            if (file.eof())
              break;

            String::size_type begin(line.find_first_not_of(" "));
            line.erase(0, begin);
            String::size_type end(line.find_first_of(" "));
            String srow(line, 0, end);
            IT_ row((IT_)atol(srow.c_str()));
            --row;
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            line.erase(0, end);

            begin = line.find_first_not_of(" ");
            line.erase(0, begin);
            end = line.find_first_of(" ");
            String sval(line, 0, end);
            DT_ tval((DT_)atof(sval.c_str()));

            *pval = tval;
            *pind = row;
            ++pval;
            ++pind;
          }
          SparseVector<Mem::Main, DT_, IT_> tmp(rows, val, ind, false);
          this->assign(tmp);
        }

        template <typename T1_, typename T2_>
        static void _insertion_sort(T1_ * key, T2_ * val1, Index size)
        {
          T1_ swap_key;
          T2_ swap1;
          for (Index i(1), j ; i < size ; ++i)
          {
            swap_key = Util::MemoryPool<Mem_>::get_element(key, i);
            swap1 = Util::MemoryPool<Mem_>::get_element(val1, i);
            j = i;
            while (j > 0 && Util::MemoryPool<Mem_>::get_element(key, j - 1) > swap_key)
            {
              Util::MemoryPool<Mem_>::copy(key + j, key + j - 1, 1);
              Util::MemoryPool<Mem_>::copy(val1 + j, val1 + j - 1, 1);
              --j;
            }
            Util::MemoryPool<Mem_>::set_memory(key + j, swap_key);
            Util::MemoryPool<Mem_>::set_memory(val1 + j, swap1);
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
          this->_scalar_index.push_back(Math::min<Index>(0, 1000));
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
          this->_scalar_index.push_back(Math::min<Index>(size_in, 1000));
          this->_scalar_index.push_back(1);
          this->_scalar_dt.push_back(DT_(0));
        }

        /**
         * \brief Constructor
         *
         * \param[in] size The size of the created vector.
         * \param[in] elements A list of non zero elements.
         * \param[in] indices A list of non zero element indices.
         * \param[in] is_sorted Indicates, if the elements are sorted by their indices: is_sorted = true (default)
         *
         * Creates a vector with a given size.
         */
        explicit SparseVector(Index size_in, DenseVector<Mem_, DT_, IT_> & elements_in,
                              DenseVector<Mem_, IT_, IT_> & indices_in, bool is_sorted = true) :
          Container<Mem_, DT_, IT_>(size_in)
        {
          CONTEXT("When creating SparseVector");

          if (indices_in.size() != elements_in.size())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size missmatch!");

          this->_scalar_index.push_back(elements_in.size());
          this->_scalar_index.push_back(elements_in.size());
          this->_scalar_index.push_back(Math::min<Index>(size_in, 1000));
          this->_scalar_index.push_back(Index(is_sorted));
          this->_scalar_dt.push_back(DT_(0));

          this->_elements.push_back(elements_in.elements());
          this->_elements_size.push_back(elements_in.size());
          this->_indices.push_back(indices_in.elements());
          this->_indices_size.push_back(indices_in.size());

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            Util::MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            Util::MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
        }

        /**
         * \brief Constructor
         *
         * \param[in] std::vector<char> A std::vector, containing the byte array.
         *
         * Creates a vector from the given byte array.
         */
        template <typename DT2_ = DT_, typename IT2_ = IT_>
        explicit SparseVector(std::vector<char> input) :
          Container<Mem_, DT_, IT_>(0)
        {
          CONTEXT("When creating SparseVector");
          deserialise<DT2_, IT2_>(input);
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
         * \brief Constructor
         *
         * \param[in] mode The used file format.
         * \param[in] filename The source file.
         *
         * Creates a vector based on the source file.
         */
        explicit SparseVector(FileMode mode, String filename) :
          Container<Mem_, DT_, IT_>(0)
        {
          CONTEXT("When creating SparseVector");

          switch(mode)
          {
            case FileMode::fm_mtx:
              _read_from_mtx(filename);
              break;
            default:
              throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
          }
        }

        /**
         * \brief Constructor
         *
         * \param[in] mode The used file format.
         * \param[in] file The source filestream.
         *
         * Creates a vector based on the source filestream.
         */
        explicit SparseVector(FileMode mode, std::istream& file) :
          Container<Mem_, DT_, IT_>(0)
        {
          CONTEXT("When creating SparseVector");

          switch(mode)
          {
            case FileMode::fm_mtx:
              _read_from_mtx(file);
              break;
            default:
              throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
          }
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
         * \param[in] clone_indices Should we create a deep copy of the index arrays, too ?
         *
         * \return A deep copy of itself.
         *
         */
        SparseVector clone(bool clone_indices = true) const
        {
          CONTEXT("When cloning SparseVector");
          SparseVector t;
          t.clone(*this, clone_indices);
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

        /** \brief Shallow copy operation
         *
         * Create a shallow copy of itself.
         *
         */
        SparseVector shared() const
        {
          CONTEXT("When sharing SparseVector");
          SparseVector r;
          r.assign(*this);
          return r;
        }

        /**
         * \brief Conversion method
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
         * \brief Retrieve specific vector element.
         *
         * \param[in] index The index of the vector element.
         *
         * \returns Specific vector element.
         */
        const DT_ operator()(Index index) const
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
            if (Util::MemoryPool<Mem_>::get_element(indices(), i) >= index)
              break;
            ++i;
          }

          if (i < used_elements() && Util::MemoryPool<Mem_>::get_element(indices(), i) == index)
            return Util::MemoryPool<Mem_>::get_element(elements(), i);
          else
            return zero_element();
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

          // flag container as not sorted anymore
          // CAUTION: do not use any method triggering resorting until we are finished
          _sorted() = 0;

          // vector is empty, no arrays allocated
          if (this->_elements.size() == 0)
          {
            this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(alloc_increment()));
            this->_elements_size.push_back(alloc_increment());
            Util::MemoryPool<Mem_>::instance()->template set_memory<DT_>(this->_elements.back(), DT_(4711), alloc_increment());
            this->_indices.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(alloc_increment()));
            this->_indices_size.push_back(alloc_increment());
            Util::MemoryPool<Mem_>::instance()->template set_memory<IT_>(this->_indices.back(), IT_(4711), alloc_increment());
            _allocated_elements() = alloc_increment();
            Util::MemoryPool<Mem_>::set_memory(this->_elements.at(0), val);
            Util::MemoryPool<Mem_>::set_memory(this->_indices.at(0), IT_(index));
            _used_elements() = 1;
          }

          // append element in already allocated arrays
          else if(_used_elements() < allocated_elements())
          {
            Util::MemoryPool<Mem_>::set_memory(this->_elements.at(0) + _used_elements(), val);
            Util::MemoryPool<Mem_>::set_memory(this->_indices.at(0) + _used_elements(), IT_(index));
            ++_used_elements();
          }

          // reallocate arrays, append element
          else
          {
            _allocated_elements() += alloc_increment();

            DT_ * elements_new(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(allocated_elements()));
            Util::MemoryPool<Mem_>::instance()->template set_memory<DT_>(elements_new, DT_(4711), allocated_elements());
            IT_ * indices_new(Util::MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(allocated_elements()));
            Util::MemoryPool<Mem_>::instance()->template set_memory<IT_>(indices_new, IT_(4711), allocated_elements());

            Util::MemoryPool<Mem_>::copy(elements_new, this->_elements.at(0), _used_elements());
            Util::MemoryPool<Mem_>::copy(indices_new, this->_indices.at(0), _used_elements());

            Util::MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(0));
            Util::MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(0));

            this->_elements.at(0) = elements_new;
            this->_indices.at(0) = indices_new;

            Util::MemoryPool<Mem_>::set_memory(this->_elements.at(0) + _used_elements(), val);
            Util::MemoryPool<Mem_>::set_memory(this->_indices.at(0) + _used_elements(), IT_(index));

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
            // CAUTION: do not use any method triggering resorting until we are finished
            _sorted() = 1;

            // check if there is anything to be sorted
            if(_used_elements() == Index(0))
              return;

            _insertion_sort(this->_indices.at(0), this->_elements.at(0), _used_elements());

            // find and mark duplicate entries
            for (Index i(1) ; i < _used_elements() ; ++i)
            {
              if (Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), i - 1) == Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), i))
              {
                Util::MemoryPool<Mem_>::set_memory(this->_indices.at(0) + i - 1, std::numeric_limits<IT_>::max());
              }
            }

            // sort out marked duplicated elements
            _insertion_sort(this->_indices.at(0), this->_elements.at(0), _used_elements());
            Index junk(0);
            while (Util::MemoryPool<Mem_>::get_element(this->_indices.at(0), _used_elements() - 1 - junk) == std::numeric_limits<IT_>::max()
                && junk < _used_elements())
              ++junk;
            _used_elements() -= junk;
          }
        }

        /**
         * \brief Deserialisation of complete container entity.
         *
         * \param[in] std::vector<char> A std::vector, containing the byte array.
         *
         * Recreate a complete container entity by a single binary array.
         */
        template <typename DT2_ = DT_, typename IT2_ = IT_>
        void deserialise(std::vector<char> input)
        {
          this->template _deserialise<DT2_, IT2_>(FileMode::fm_sv, input);
        }

        /**
         * \brief Serialisation of complete container entity.
         *
         * \param[in] mode FileMode enum, describing the actual container specialisation.
         * \param[out] std::vector<char> A std::vector, containing the byte array.
         *
         * Serialize a complete container entity into a single binary array.
         *
         * See \ref FEAST::LAFEM::Container::_serialise for details.
         */
        template <typename DT2_ = DT_, typename IT2_ = IT_>
        std::vector<char> serialise()
        {
          return this->template _serialise<DT2_, IT2_>(FileMode::fm_sv);
        }

        /**
         * \brief Write out vector to file.
         *
         * \param[in] mode The used file format.
         * \param[in] filename The file where the matrix shall be stored.
         */
        void write_out(FileMode mode, String filename) const
        {
          CONTEXT("When writing out SparseVector");

          switch(mode)
          {
            case FileMode::fm_mtx:
              write_out_mtx(filename);
              break;
            default:
                throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
          }
        }

        /**
         * \brief Write out vector to file.
         *
         * \param[in] mode The used file format.
         * \param[in] file The stream that shall be written to.
         */
        void write_out(FileMode mode, std::ostream& file) const
        {
          CONTEXT("When writing out SparseVector");

          switch(mode)
          {
            case FileMode::fm_mtx:
              write_out_mtx(file);
              break;
            default:
                throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
          }
        }

        /**
         * \brief Write out matrix to MatrixMarktet mtx file.
         *
         * \param[in] filename The file where the vector shall be stored.
         */
        void write_out_mtx(String filename) const
        {
          std::ofstream file(filename.c_str(), std::ofstream::out);
          if (! file.is_open())
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
          write_out_mtx(file);
          file.close();
        }

        /**
         * \brief Write out matrix to MatrixMarktet mtx file.
         *
         * \param[in] file The stream that shall be written to.
         */
        void write_out_mtx(std::ostream& file) const
        {
          SparseVector<Mem::Main, DT_, IT_> temp;
          temp.convert(*this);

          file << "%%MatrixMarket matrix coordinate real general" << std::endl;
          file << temp.size() << " " << 1 << " " << temp.used_elements() << std::endl;

          const Index u_elem(temp.used_elements());
          const IT_ * pind(temp.indices());
          const DT_ * pval(temp.elements());
          for (Index i(0) ; i < u_elem ; ++i, ++pind, ++pval)
          {
              file << *pind+1 << " " << 1 << " " << std::scientific << *pval << std::endl;
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


        /**
         * \brief SparseVector comparison operator
         *
         * \param[in] a A vector to compare with.
         * \param[in] b A vector to compare with.
         */
        template <typename Mem2_> friend bool operator== (const SparseVector & a, const SparseVector<Mem2_, DT_, IT_> & b)
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
        friend std::ostream & operator<< (std::ostream & lhs, const SparseVector & b)
        {
          lhs << "[";
          for (Index i(0) ; i < b.size() ; ++i)
          {
            lhs << "  " << b(i);
          }
          lhs << "]";

          return lhs;
        }
    }; // class SparseVector<...>


  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_SPARSE_VECTOR_HPP
