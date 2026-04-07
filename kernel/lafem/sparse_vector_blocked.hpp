// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/adjacency/permutation.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/arch/min_max_value_dense.hpp>
#include <kernel/lafem/arch/max_rel_diff_dense.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Sparse vector class template.
     *
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indexing type to be used.
     * \tparam block_size_ The size of the represented blocks
     *
     * This class represents a vector with non zero element blocks in a sparse layout. \n
     * Logical, the data are organized in small blocks of block_size_ length.\n\n
     * Data survey: \n
     * _elements[0]: raw number values \n
     * _indices[0]: non zero indices \n
     * _scalar_index[0]: container size - aka block count \n
     * _scalar_index[1]: non zero element count (used elements) \n
     * _scalar_index[2]: allocated elements \n
     * _scalar_index[3]: allocation size increment \n
     * _scalar_index[4]: boolean flag, if container is sorted \n
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename DT_, typename IT_, int block_size_>
    class SparseVectorBlocked : public Container<DT_, IT_>
    {
    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our size of a single block
      static constexpr int block_size = block_size_;
      /// Our value type
      typedef Tiny::Vector<DT_, block_size_> ValueType;

    private:
      Index & _size()
      {
        return this->_scalar_index.at(0);
      }

      Index & _nonzeros()
      {
        return this->_scalar_index.at(1);
      }

    public:
      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional vector.
       */
      explicit SparseVectorBlocked()
      {
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);
      }

      /**
       * \brief Constructor
       *
       * \param[in] size_in The size of the created vector.
       *
       * Creates a vector with a given size.
       */
      explicit SparseVectorBlocked(Index size_in, Index nnze_in)
      {
        this->_scalar_index.push_back(size_in);
        this->_scalar_index.push_back(nnze_in);
        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * nnze_in * Index(block_size_)));
        this->_elements_size.push_back(nnze_in * Index(block_size_));
        this->_indices.push_back(Memory::Arbiter(sizeof(IT_) * nnze_in));
        this->_indices_size.push_back(nnze_in);
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
      explicit SparseVectorBlocked(Index size_in, DenseVectorBlocked<DT_, IT_, block_size_> & elements_in,
        DenseVector<IT_, IT_> & indices_in)
      {
        XASSERT(size_in != Index(0));
        XASSERTM(indices_in.size() == elements_in.size(), "Vector size mismatch!");

        this->_scalar_index.push_back(size_in);
        this->_scalar_index.push_back(elements_in.size());

        this->_elements.push_back(elements_in.elements_arbiter().attach());
        this->_elements_size.push_back(elements_in.size() * Index(block_size_));
        this->_indices.push_back(indices_in.elements_arbiter().attach());
        this->_indices_size.push_back(indices_in.size());
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit SparseVectorBlocked(std::vector<char> input)
      {
        this->deserialize<DT2_,IT2_>(input);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a vector based on the source file.
       */
      explicit SparseVectorBlocked(FileMode mode, const String& filename)
      {
        this->read_from(mode, filename);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The source filestream.
       *
       * Creates a vector based on the source filestream.
       */
      explicit SparseVectorBlocked(FileMode mode, std::istream& file)
      {
        this->read_from(mode, file);
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source vector.
       *
       * Moves another vector to this vector.
       */
      SparseVectorBlocked(SparseVectorBlocked && other) :
        Container<DT_, IT_>(std::forward<SparseVectorBlocked>(other))
      {
      }

      /**
       * \brief Assignment move operator
       *
       * \param[in] other The source vector.
       *
       * Moves another vector to the target vector.
       */
      SparseVectorBlocked & operator= (SparseVectorBlocked && other)
      {
        this->move(std::forward<SparseVectorBlocked>(other));

        return *this;
      }

      /** \brief Clone operation
       *
       * Create a deep clone of this container.
       *
       * \param[in] clone_mode The actual cloning procedure.
       * \returns The created clone.
       *
       */
      SparseVectorBlocked clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        SparseVectorBlocked t;
        t.clone(*this, clone_mode);
        return t;
      }

      /** \brief Clone operation
       *
       * Create a deep clone of this container.
       *
       * \param[in] other The source container to create the clone from.
       * \param[in] clone_mode The actual cloning procedure.
       *
       */
      template<typename DT2_, typename IT2_>
      void clone(const SparseVectorBlocked<DT2_, IT2_, block_size_> & other, CloneMode clone_mode = CloneMode::Deep)
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
      void convert(const SparseVectorBlocked<DT2_, IT2_, block_size_> & other)
      {
        this->clone(other);
      }

      /// Returns the native size of the sparse vector
      Index size() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(0);
      }

      /// Returns the raw size of the sparse vector
      Index size_raw() const
      {
        return this->size() * Index(block_size_);
      }

      /// Returns the number of non-zero elements in the sparse vector
      Index num_nzes() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(1);
      }

      /// Returns the raw number of non-zero elements in the sparse vector
      Index num_nzes_raw() const
      {
        return this->num_nzes() * Index(block_size_);
      }

      /// Checks whether the size of the sparse vector is zero
      bool empty() const
      {
        return this->_scalar_index.empty() || (this->_scalar_index.at(0) == Index(0));
      }

      /// Checks whether the sparse vector does not contain any non-zero elements
      bool hollow() const
      {
        return this->_scalar_index.empty() || (this->_scalar_index.at(1) == Index(0));
      }

      Memory::Arbiter& elements_arbiter()
      {
        return this->_elements.front();
      }

      const Memory::Arbiter& elements_arbiter() const
      {
        return this->_elements.front();
      }

      Memory::TypedView<ValueType> elements_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_elements.empty())
          return Memory::TypedView<ValueType>();
        return Memory::TypedView<ValueType>(this->_elements.at(0).view(loc, Memory::Access::read));
      }

      Memory::TypedView<ValueType> elements_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<ValueType>();
        return Memory::TypedView<ValueType>(this->_elements.at(0).view(loc, Memory::Access::write));
      }

      Memory::TypedView<ValueType> elements_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<ValueType>();
        return Memory::TypedView<ValueType>(this->_elements.at(0).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<ValueType> elements_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_elements.empty())
          return Memory::TypedView<ValueType>();
        return Memory::TypedView<ValueType>(this->_elements.at(0).view(loc, acc));
      }

      Memory::Arbiter& indices_arbiter()
      {
        return this->_indices.front();
      }

      const Memory::Arbiter& indices_arbiter() const
      {
        return this->_indices.front();
      }

      Memory::TypedView<IT_> indices_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::read));
      }

      Memory::TypedView<IT_> indices_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::write));
      }

      Memory::TypedView<IT_> indices_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<IT_> indices_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, acc));
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
        this->template _deserialize<DT2_, IT2_>(FileMode::fm_svb, input);
      }

      bool element_exists(IndexType ind) const
      {
        const Memory::TypedView<IT_> idx_view = this->indices_view_r();
        const IT_* pindices = idx_view.get_r();
        return std::find(pindices, pindices + this->_scalar_index.at(1), ind) != pindices + this->_scalar_index.at(1);
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
        return this->template _serialize<DT2_, IT2_>(FileMode::fm_svb, config);
      }

      /**
       * \brief Read in vector from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, const String& filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
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
          case FileMode::fm_svb:
            this->template _deserialize<double, std::uint64_t>(FileMode::fm_svb, file);
            break;
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
          buff = new char[LAFEM::file_out_stream_buffer_size];
          file.rdbuf()->pubsetbuf(buff, LAFEM::file_out_stream_buffer_size);
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
          [[fallthrough]];

        case FileMode::fm_svb:
          this->template _serialize<double, std::uint64_t>(FileMode::fm_svb, file);
          break;

        case FileMode::fm_mtx:
          {

            file << "%%MatrixMarket matrix coordinate real general\n";
            file << this->size()*block_size << " " << 1 << " " << this->num_nzes()*block_size << "\n";
            file << std::scientific << std::setprecision(Type::Traits<DT_>::format_precision);

            const Index u_elem(this->num_nzes());
            const Memory::TypedView<IT_> idx_view = this->indices_view_r();
            const Memory::TypedView<ValueType> val_view = this->elements_view_r();
            for (Index i(0) ; i < u_elem ; ++i)
            {
              for(int k(0); k < block_size; ++k)
                file << (idx_view[i]*Index(block_size+k+1)) << " " << 1 << " " << val_view[i][k] << "\n";
            }
            break;
          }

          default:
            XABORTM("Filemode not supported!");
        }
      }

      ///@name Linear algebra operations
      ///@{

      /**
       * \brief Retrieve the absolute maximum value of this vector.
       *
       * \return The largest absolute value.
       */
      DT_ max_abs_element() const
      {
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->elements_arbiter(), this->num_nzes_raw(), false, true);

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
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->elements_arbiter(), this->num_nzes_raw(), true, true);

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
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->elements_arbiter(), this->num_nzes_raw(), false, false);

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
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->elements_arbiter(), this->num_nzes_raw(), true, false);

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
      DT_ max_rel_diff(const SparseVectorBlocked & x) const
      {
        XASSERTM(x.num_nzes() == this->num_nzes(), "Nonzero count does not match!");
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DataType result = Arch::MaxRelDiffDense::template exec<DT_>(this->elements_arbiter(), x.elements_arbiter(), this->num_nzes_raw());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      ///@}

      /// Permute vector according to the given Permutation
      void permute(Adjacency::Permutation & perm)
      {
        if (perm.size() == 0)
          return;

        std::map<IT_, ValueType> xm;

        auto ip = perm.inverse();
        const Index * const inv_pos(ip.get_perm_pos());

        Memory::TypedView<ValueType> xv(this->elements_view_rw());
        Memory::TypedView<IT_> xi(this->indices_view_rw());
        for (Index i(0) ; i < this->num_nzes() ; ++i)
        {
          xm.emplace(IT_(inv_pos[xi[i]]), xv[i]);
        }

        Index k(0);
        for(const auto& x : xm)
        {
          xi[k] = x.first;
          xv[k] = x.second;
          ++k;
        }
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "SparseVectorBlocked";
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
      bool same_layout(const SparseVectorBlocked& other) const
      {
        if(this->size() != other.size())
          return false;
        if(this->num_nzes() != other.num_nzes())
          return false;

        if(this->num_nzes() == Index(0))
          return true;

        // check if the indices arbiters are the same
        if(this->indices_arbiter() == other.indices_arbiter())
          return true;

        const Memory::TypedView<IT_> idx_a = this->indices_view_r();
        const Memory::TypedView<IT_> idx_b = other.indices_view_r();

        // check all array entries
        const Index n = this->num_nzes();
        for(Index i(0); i < n; ++i)
        {
          if(idx_a[i] != idx_b[i])
            return false;
        }

        // okay, arrays are identical
        return true;
      }

      /**
       * \brief SparseVectorBlocked streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The vector to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & lhs, const SparseVectorBlocked & b)
      {
        const Memory::TypedView<IT_> idx_view(b.indices_view_r());
        const Memory::TypedView<ValueType> val_view(b.elements_view_r());
        const Index n = b.size();
        const Index nze = b.num_nzes();
        lhs << "[";
        Index k = 0;
        for (Index i(0) ; i < nze ; ++i)
        {
          for(; k < idx_view[i]; ++k)
          {
            for(int d = 0; d < block_size_; ++d)
              lhs << "  " << DT_(0);
          }

          for(int d = 0; d < block_size_; ++d)
            lhs << "  " << val_view[i][d];
          ++k;
        }
        for(; k < n; ++k)
          for(int d = 0; d < block_size_; ++d)
            lhs << "  " << DT_(0);
        lhs << "]";

        return lhs;
      }
    }; // class SparseVectorBlocked<...>
  } // namespace LAFEM
} // namespace FEAT
