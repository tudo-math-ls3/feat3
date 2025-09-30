// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/math.hpp>
#include <kernel/adjacency/permutation.hpp>
#include <kernel/lafem/arch/max_abs_index.hpp>
#include <kernel/lafem/arch/min_abs_index.hpp>
#include <kernel/lafem/arch/max_index.hpp>
#include <kernel/lafem/arch/min_index.hpp>
#include <kernel/lafem/arch/max_rel_diff.hpp>


namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Sparse vector class template.
     *
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
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename DT_, typename IT_ = Index>
    class SparseVector : public Container<DT_, IT_>
    {
    private:
      template <typename T1_, typename T2_>
      static void _insertion_sort(T1_ * key, T2_ * val1, Index size)
      {
        T1_ swap_key;
        T2_ swap1;
        for (Index i(1), j ; i < size ; ++i)
        {
          swap_key = key[i];
          swap1 = val1[i];
          j = i;
          while (j > 0 && key[j-1] > swap_key)
          {
            key[j] = key[j-1];
            val1[j] = val1[j-1];
            --j;
          }
          key[j] = swap_key;
          val1[j] = swap1;
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

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional vector.
       */
      explicit SparseVector() :
        Container<DT_, IT_> (0)
      {
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(Math::min<Index>(0, 1000));
        this->_scalar_index.push_back(1);
      }

      /**
       * \brief Constructor
       *
       * \param[in] size_in The size of the created vector.
       *
       * Creates a vector with a given size.
       */
      explicit SparseVector(Index size_in) :
        Container<DT_, IT_>(size_in)
      {
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(Math::min<Index>(size_in, 1000));
        this->_scalar_index.push_back(1);
      }

      /**
       * \brief Constructor
       *
       * \param[in] size_in The size of the created vector.
       * \param[in] elements_in A list of non zero elements.
       * \param[in] indices_in A list of non zero element indices.
       * \param[in] is_sorted Indicates, if the elements are sorted by their indices: is_sorted = true (default)
       *
       * Creates a vector with a given size.
       */
      explicit SparseVector(Index size_in, DenseVector<DT_, IT_> & elements_in,
                            DenseVector<IT_, IT_> & indices_in, bool is_sorted = true) :
        Container<DT_, IT_>(size_in)
      {
        XASSERT(size_in != Index(0));
        XASSERTM(indices_in.size() == elements_in.size(), "Vector size mismatch!");

        this->_scalar_index.push_back(elements_in.size());
        this->_scalar_index.push_back(elements_in.size());
        this->_scalar_index.push_back(Math::min<Index>(size_in, 1000));
        this->_scalar_index.push_back(Index(is_sorted));

        this->_elements.push_back(elements_in.elements());
        this->_elements_size.push_back(elements_in.size());
        this->_indices.push_back(indices_in.elements());
        this->_indices_size.push_back(indices_in.size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool::increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool::increase_memory(this->_indices.at(i));

        this->sort();
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Creates a vector from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit SparseVector(std::vector<char> input) :
        Container<DT_, IT_>(0)
      {
        deserialize<DT2_, IT2_>(input);
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source vector.
       *
       * Moves another vector to this vector.
       */
      SparseVector(SparseVector && other) :
        Container<DT_, IT_>(std::forward<SparseVector>(other))
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a vector based on the source file.
       */
      explicit SparseVector(FileMode mode, const String& filename) :
        Container<DT_, IT_>(0)
      {
        read_from(mode, filename);
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
        Container<DT_, IT_>(0)
      {
        read_from(mode, file);
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
        this->move(std::forward<SparseVector>(other));

        return *this;
      }

      /** \brief Clone operation
       *
       * Create a clone of this container.
       *
       * \param[in] clone_mode The actual cloning procedure.
       * \returns The created clone.
       *
       */
      SparseVector clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        SparseVector t;
        t.clone(*this, clone_mode);
        return t;
      }

      /** \brief Clone operation
       *
       * Create a clone of another container.
       *
       * \param[in] other The source container to create the clone from.
       * \param[in] clone_mode The actual cloning procedure.
       *
       */
      template<typename DT2_, typename IT2_>
      void clone(const SparseVector<DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Deep)
      {
        Container<DT_, IT_>::clone(other, clone_mode);
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
      template <typename DT2_, typename IT2_>
      void convert(const SparseVector<DT2_, IT2_> & other)
      {
        // clean up previous mess before creating a new vector
        this->sort();
        this->clone(other);
      }

      /**
       * \brief Get a pointer to the data array.
       *
       * \returns Pointer to the data array.
       */
      template <Perspective = Perspective::native>
      DT_ * elements()
      {
        if(this->_elements.empty())
          return nullptr;
        if (sorted() == 0)
          const_cast<SparseVector *>(this)->sort();
        return this->_elements.at(0);
      }

      template <Perspective = Perspective::native>
      DT_ const * elements() const
      {
        if(this->_elements.empty())
          return nullptr;
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
        if(this->_indices.empty())
          return nullptr;
        if (sorted() == 0)
          const_cast<SparseVector *>(this)->sort();
        return this->_indices.at(0);
      }

      IT_ const * indices() const
      {
        if(this->_indices.empty())
          return nullptr;
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
        ASSERTM(index < this->_scalar_index.at(0), "index exceeds sparse vector size");

        MemoryPool::synchronize();

        if (this->_elements.empty())
          return DT_(0.);

        if (sorted() == 0)
          const_cast<SparseVector *>(this)->sort();

        Index i(0);
        while (i < used_elements())
        {
          if (indices()[i] >= index)
            break;
          ++i;
        }

        if (i < used_elements() && indices()[i] == index)
          return elements()[i];
        else
          return DT_(0.);
      }


      /**
       * \brief Set specific vector element.
       *
       * \param[in] index The index of the vector element.
       * \param[in] val The val to be set.
       */
      void operator()(Index index, DT_ val)
      {
        ASSERTM(index < this->_scalar_index.at(0), "index exceeds sparse vector size");

        // flag container as not sorted anymore
        // CAUTION: do not use any method triggering resorting until we are finished
        _sorted() = 0;

        // vector is empty, no arrays allocated
        if (this->_elements.empty())
        {
          this->_elements.push_back(MemoryPool::template allocate_memory<DT_>(alloc_increment()));
          this->_elements_size.push_back(alloc_increment());
          MemoryPool::template set_memory<DT_>(this->_elements.back(), DT_(4711), alloc_increment());
          this->_indices.push_back(MemoryPool::template allocate_memory<IT_>(alloc_increment()));
          this->_indices_size.push_back(alloc_increment());
          MemoryPool::template set_memory<IT_>(this->_indices.back(), IT_(4711), alloc_increment());
          this->_allocated_elements() = alloc_increment();
          this->_elements.at(0)[0] = val;
          this->_indices.at(0)[0] = IT_(index);
          this->_used_elements() = 1;
        }

        // append element in already allocated arrays
        else if(_used_elements() < allocated_elements())
        {
          this->_elements.at(0)[this->_used_elements()] = val;
          this->_indices.at(0)[this->_used_elements()] = IT_(index);
          ++this->_used_elements();
        }

        // reallocate arrays, append element
        else
        {
          _allocated_elements() += alloc_increment();

          DT_ * elements_new(MemoryPool::template allocate_memory<DT_>(allocated_elements()));
          MemoryPool::template set_memory<DT_>(elements_new, DT_(4711), allocated_elements());
          IT_ * indices_new(MemoryPool::template allocate_memory<IT_>(allocated_elements()));
          MemoryPool::template set_memory<IT_>(indices_new, IT_(4711), allocated_elements());

          MemoryPool::copy(elements_new, this->_elements.at(0), _used_elements());
          MemoryPool::copy(indices_new, this->_indices.at(0), _used_elements());

          MemoryPool::release_memory(this->_elements.at(0));
          MemoryPool::release_memory(this->_indices.at(0));

          this->_elements.at(0) = elements_new;
          this->_indices.at(0) = indices_new;

          this->_elements.at(0)[this->_used_elements()] = val;
          this->_indices.at(0)[this->_used_elements()] = IT_(index);

          ++this->_used_elements();
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

          IT_ * pindices = this->_indices.at(0);
          DT_ * pelements = this->_elements.at(0);

          _insertion_sort(pindices, pelements, _used_elements());

          // find and mark duplicate entries
          for (Index i(1) ; i < _used_elements() ; ++i)
          {
            if (pindices[i-1] == pindices[i])
            {
              pindices[i-1] = std::numeric_limits<IT_>::max();
            }
          }

          // sort out marked duplicated elements
          _insertion_sort(pindices, pelements, _used_elements());
          Index junk(0);
          while (pindices[_used_elements() - 1 - junk] == std::numeric_limits<IT_>::max() && junk < _used_elements())
            ++junk;
          _used_elements() -= junk;
        }
      }

      /**
       * \brief Retrieve the absolute maximum value of this vector.
       *
       * \return The largest absolute value.
       */
      DT_ max_abs_element() const
      {
        TimeStamp ts_start;

        Index max_abs_index = Arch::MaxAbsIndex::value(this->template elements<Perspective::pod>(), this->template  used_elements<Perspective::pod>());
        ASSERT(max_abs_index < this->template  used_elements<Perspective::pod>());
        DT_ result(this->template elements<Perspective::pod>()[max_abs_index]);
        result = Math::abs(result);

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Retrieve the absolute minimum value of this vector.
       *
       * \return The smallest absolute value.
       */
      DT_ min_abs_element() const
      {
        TimeStamp ts_start;

        Index min_abs_index = Arch::MinAbsIndex::value(this->template elements<Perspective::pod>(), this->template  used_elements<Perspective::pod>());
        ASSERT(min_abs_index < this->template  used_elements<Perspective::pod>());
        DT_ result(this->template elements<Perspective::pod>()[min_abs_index]);
        result = Math::abs(result);

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Retrieve the maximum value of this vector.
       *
       * \return The largest value.
       */
      DT_ max_element() const
      {
        TimeStamp ts_start;

        Index max_index = Arch::MaxIndex::value(this->template elements<Perspective::pod>(), this->template  used_elements<Perspective::pod>());
        ASSERT(max_index < this->template  used_elements<Perspective::pod>());
        DT_ result(this->template elements<Perspective::pod>()[max_index]);

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Retrieve the minimum value of this vector.
       *
       * \return The smallest value.
       */
      DT_ min_element() const
      {
        TimeStamp ts_start;

        Index min_index = Arch::MinIndex::value(this->template elements<Perspective::pod>(), this->template  used_elements<Perspective::pod>());
        ASSERT(min_index < this->template used_elements<Perspective::pod>());
        DT_ result(this->template elements<Perspective::pod>()[min_index]);

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Retrieve the maximum relative difference of this vector and another one
       * y.max_rel_diff(x) returns  \f$ \max_{0\leq i < n}\frac{|x_i-y_i|}{\max{|x_i|+|y_i|, eps}} \f$
       *
       * \return The largest relative difference.
       */
      DT_ max_rel_diff(const SparseVector & x) const
      {
        XASSERTM(x.used_elements() == this->used_elements(), "Nonzero count does not match!");
        TimeStamp ts_start;

        DataType max_rel_diff = Arch::MaxRelDiff::value(this->template elements<Perspective::pod>(), x.template elements<Perspective::pod>(), this->template size<Perspective::pod>());
        ASSERT(max_rel_diff < this->template size<Perspective::pod>());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return max_rel_diff;
      }

      /**
       * \brief Deserialization of complete container entity.
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Recreate a complete container entity by a single binary array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void deserialize(std::vector<char> input)
      {
        this->template _deserialize<DT2_, IT2_>(FileMode::fm_sv, input);
      }

      /**
       * \brief Serialization of complete container entity.
       *
       * \param[in] config LAFEM::SerialConfig, a struct describing the serialize configuration.
       * \note the corresponding configure flags 'zlib' and/or 'zfp' need to be added in the build-id at the configure call.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAT::LAFEM::Container::_serialize for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialize(const LAFEM::SerialConfig& config = SerialConfig()) const
      {
        return this->template _serialize<DT2_, IT2_>(FileMode::fm_sv, config);
      }

      /**
       * \brief Read in vector from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, const String& filename)
      {
        std::ios_base::openmode bin = std::ifstream::in | std::ifstream::binary;
        if(mode == FileMode::fm_mtx)
          bin = std::ifstream::in;
        std::ifstream file(filename.c_str(), bin);
        if (! file.is_open())
          XABORTM("Unable to open Vector file " + filename);
        read_from(mode, file);
        file.close();
      }

      /**
       * \brief Read in vector from stream.
       *
       * \param[in] mode The used file format.
       * \param[in] file The stream that shall be read in.
       */
      void read_from(FileMode mode, std::istream& file)
      {
        switch(mode)
        {
          case FileMode::fm_binary:
          case FileMode::fm_sv:
            this->template _deserialize<double, std::uint64_t>(FileMode::fm_sv, file);
            break;
          case FileMode::fm_mtx:
          {
            this->clear();
            this->_scalar_index.push_back(0);
            this->_scalar_index.push_back(0);
            this->_scalar_index.push_back(Math::min<Index>(0, 1000));
            this->_scalar_index.push_back(1);

            Index rows;
            Index nnz;
            String line;
            std::getline(file, line);
            if (line.find("%%MatrixMarket matrix coordinate real general") == String::npos)
            {
              XABORTM("Input-file is not a compatible mtx-file");
            }
            while(!file.eof())
            {
              std::getline(file,line);
              if (file.eof())
                XABORTM("Input-file is empty");

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
                XABORTM("Input-file is no sparse-vector-file");

              begin = line.find_first_not_of(" ");
              line.erase(0, begin);
              end = line.find_first_of(" ");
              String snnz(line, 0, end);
              nnz = (Index)atol(snnz.c_str());
              line.erase(0, end);
            }

            DenseVector<IT_, IT_> ind(nnz);
            DenseVector<DT_, IT_> val(nnz);

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
            this->move(SparseVector<DT_, IT_>(rows, val, ind, false));
            break;
          }
          default:
            XABORTM("Filemode not supported!");
        }
      }

      /**
       * \brief Write out vector to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the matrix shall be stored.
       */
      void write_out(FileMode mode, const String& filename) const
      {
        std::ios_base::openmode bin = std::ofstream::out | std::ofstream::binary;
        if(mode == FileMode::fm_mtx)
          bin = std::ofstream::out;
        std::ofstream file;
        char* buff = nullptr;
        if(mode == FileMode::fm_mtx)
        {
          buff = new char[LAFEM::FileOutStreamBufferSize];
          file.rdbuf()->pubsetbuf(buff, LAFEM::FileOutStreamBufferSize);
        }
        file.open(filename.c_str(), bin);
        if(! file.is_open())
          XABORTM("Unable to open Matrix file " + filename);
        write_out(mode, file);
        file.close();
        delete[] buff;
      }

      /**
       * \brief Write out vector to file.
       *
       * \param[in] mode The used file format.
       * \param[in] file The stream that shall be written to.
       */
      void write_out(FileMode mode, std::ostream& file) const
      {
        switch(mode)
        {
          case FileMode::fm_binary:
          case FileMode::fm_sv:
            this->template _serialize<double, std::uint64_t>(FileMode::fm_sv, file);
            break;
          case FileMode::fm_mtx:
          {

            file << "%%MatrixMarket matrix coordinate real general\n";
            file << this->size() << " " << 1 << " " << this->used_elements() << "\n";

            const Index u_elem(this->used_elements());
            const IT_ * pind(this->indices());
            const DT_ * pval(this->elements());
            for (Index i(0) ; i < u_elem ; ++i, ++pind, ++pval)
            {
              file << *pind+1 << " " << 1 << " " << stringify_fp_sci(*pval) << "\n";
            }
            break;
          }
          default:
            XABORTM("Filemode not supported!");
        }
      }

      /**
       * \brief Retrieve non zero element count.
       *
       * \returns Non zero element count.
       */
      template <Perspective = Perspective::native>
      Index used_elements() const
      {
        if(this->_scalar_index.empty())
          return Index(0);
        if (sorted() == 0)
          const_cast<SparseVector *>(this)->sort();
        return this->_scalar_index.at(1);
      }

      /**
       * \brief Retrieve amount of allocated elements.
       *
       * \return Allocated element count.
       */
      Index allocated_elements() const
      {
        if(this->_scalar_index.empty())
          return Index(0);
        return this->_scalar_index.at(2);
      }

      /**
       * \brief Retrieve allocation incrementation value.
       *
       * \return Allocation increment.
       */
      Index alloc_increment() const
      {
        if(this->_scalar_index.empty())
          return Index(0);
        return this->_scalar_index.at(3);
      }

      /**
       * \brief Retrieve status of element sorting.
       *
       * \return Sorting status.
       */
      Index sorted() const
      {
        if(this->_scalar_index.empty())
          return Index(0);
        return this->_scalar_index.at(4);
      }

      /// Permutate vector according to the given Permutation
      void permute(Adjacency::Permutation & perm)
      {
        if (perm.size() == 0)
          return;

        XASSERTM(perm.size() == this->size(), "Container size does not match permutation size");

        SparseVector<DT_, IT_> target(this->size());

        auto inv = perm.inverse();
        const Index * const inv_pos(inv.get_perm_pos());
        for (Index i(0) ; i < this->used_elements() ; ++i)
        {
          const Index col = this->indices()[i];
          target(inv_pos[col], (*this)(col));
        }

        target.sort();
        this->move(std::move(target));
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
       * \brief Checks whether the layout of this vector is identical to another sparse vector.
       *
       * \param[in] other
       * A \transient reference to the other vector to compare to.
       *
       * \returns
       * \c true, if both vectors have the same layout, otherwise \c false.
       */
      bool compare_layout(const SparseVector& other) const
      {
        if(this->size() != other.size())
          return false;
        if(this->used_elements() != other.used_elements())
          return false;

        Index n = this->used_elements();
        if(n == Index(0))
          return true; // both vectors are empty

        const IT_* a = this->indices();
        const IT_* b = this->indices();

        // shallow copy? (aka same array)
        if(a == b)
          return true;

        // check all array entries
        for(Index i(0); i < n; ++i)
        {
          if(a[i] != b[i])
            return false;
        }

        // okay, arrays are identical
        return true;
      }

      /**
       * \brief SparseVector comparison operator
       *
       * \param[in] a A vector to compare with.
       * \param[in] b A vector to compare with.
       */
      friend bool operator== (const SparseVector & a, const SparseVector & b)
      {
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
          lhs << "  " << stringify(b(i));
        }
        lhs << "]";

        return lhs;
      }
    }; // class SparseVector<...>

#ifdef FEAT_EICKT
    extern template class SparseVector<float, std::uint32_t>;
    extern template class SparseVector<double, std::uint32_t>;
    extern template class SparseVector<float, std::uint64_t>;
    extern template class SparseVector<double, std::uint64_t>;
#endif


  } // namespace LAFEM
} // namespace FEAT
