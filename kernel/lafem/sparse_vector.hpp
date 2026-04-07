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
     *
     * This class represents a vector with non zero elements in a sparse layout. \n \n
     * Data survey: \n
     * _elements[0]: raw number values \n
     * _indices[0]: non zero indices \n
     * _scalar_index[0]: container size \n
     * _scalar_index[1]: non zero element count (used elements) \n
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock, Peter Zajac
     */
    template <typename DT_, typename IT_ = Index>
    class SparseVector : public Container<DT_, IT_>
    {
    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;

    protected:
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
      SparseVector()
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
      explicit SparseVector(Index size_in, Index nnze_in)
      {
        this->_scalar_index.push_back(size_in);
        this->_scalar_index.push_back(nnze_in);
        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * nnze_in));
        this->_elements_size.push_back(nnze_in);
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
      explicit SparseVector(Index size_in, DenseVector<DT_, IT_> & elements_in,
        DenseVector<IT_, IT_> & indices_in)
      {
        XASSERT(size_in != Index(0));
        XASSERTM(indices_in.size() == elements_in.size(), "Vector size mismatch!");

        this->_scalar_index.push_back(size_in);
        this->_scalar_index.push_back(elements_in.size());

        this->_elements.push_back(elements_in.elements_arbiter().attach());
        this->_elements_size.push_back(elements_in.size());
        this->_indices.push_back(indices_in.elements_arbiter().attach());
        this->_indices_size.push_back(indices_in.size());
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Creates a vector from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit SparseVector(std::vector<char> input)
      {
        this->deserialize<DT2_, IT2_>(input);
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
      explicit SparseVector(FileMode mode, const String& filename)
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
      explicit SparseVector(FileMode mode, std::istream& file)
      {
        this->read_from(mode, file);
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
        return this->size();
      }

      /// Returns the number of non-zero elements in the sparse vector
      Index num_nzes() const
      {
        return this->_scalar_index.empty() ? Index(0) : this->_scalar_index.at(1);
      }

      /// Returns the raw number of non-zero elements in the sparse vector
      Index num_nzes_raw() const
      {
        return this->num_nzes();
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

      /// Returns a reference to the elements memory arbiter
      Memory::Arbiter& elements_arbiter()
      {
        XASSERT(!this->_elements.empty());
        return this->_elements.front();
      }

      /// Returns a const reference to the elements memory arbiter
      const Memory::Arbiter& elements_arbiter() const
      {
        XASSERT(!this->_elements.empty());
        return this->_elements.front();
      }

      /// Creates a read-only memory view for the elements array for a given memory location
      Memory::TypedView<DT_> elements_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::read));
      }

      /// Creates a write-only memory view for the elements array for a given memory location
      Memory::TypedView<DT_> elements_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::write));
      }

      /// Creates a read-write memory view for the elements array for a given memory location
      Memory::TypedView<DT_> elements_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::read_write));
      }

      /**
       * \brief Creates a memory view for the elements array
       *
       * \param[in] loc
       * The memory location for which the view is to be created
       *
       * \param[in] acc
       * A combination of access rights to grant for the view
       */
      Memory::TypedView<DT_> elements_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, acc));
      }

      /// Returns a reference to the indices memory arbiter
      Memory::Arbiter& indices_arbiter()
      {
        return this->_indices.front();
      }

      /// Returns a const reference to the indices memory arbiter
      const Memory::Arbiter& indices_arbiter() const
      {
        return this->_indices.front();
      }

      /// Creates a read-only memory view for the indices array for a given memory location
      Memory::TypedView<IT_> indices_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::read));
      }

      /// Creates a write-only memory view for the indices array for a given memory location
      Memory::TypedView<IT_> indices_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::write));
      }

      /// Creates a read-write memory view for the indices array for a given memory location
      Memory::TypedView<IT_> indices_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, Memory::Access::read_write));
      }

      /**
       * \brief Creates a memory view for the indices array
       *
       * \param[in] loc
       * The memory location for which the view is to be created
       *
       * \param[in] acc
       * A combination of access rights to grant for the view
       */
      Memory::TypedView<IT_> indices_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_indices.empty())
          return Memory::TypedView<IT_>();
        return Memory::TypedView<IT_>(this->_indices.at(0).view(loc, acc));
      }

      bool element_exists(IndexType ind) const
      {
        const Memory::TypedView<IT_> idx_view = this->indices_view_r();
        const IT_* pindices = idx_view.get_r();

        return std::find(pindices, pindices + this->_scalar_index.at(1), ind) != pindices + this->_scalar_index.at(1);
      }

      /**
       * \brief Retrieve the absolute maximum value of this vector.
       *
       * \return The largest absolute value.
       */
      DT_ max_abs_element() const
      {
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->elements_arbiter(), this->num_nzes(), false, true);

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

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->elements_arbiter(), this->num_nzes(), true, true);

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

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->elements_arbiter(), this->num_nzes(), false, false);

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

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->elements_arbiter(), this->num_nzes(), true, false);

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
        XASSERTM(x.num_nzes() == this->num_nzes(), "Nonzero count does not match!");
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DataType result = Arch::MaxRelDiffDense::template exec<DT_>(this->elements_arbiter(), x.elements_arbiter(), this->num_nzes());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
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
          [[fallthrough]];

        case FileMode::fm_sv:
          this->template _deserialize<double, std::uint64_t>(FileMode::fm_sv, file);
          break;

        case FileMode::fm_mtx:
          {
            this->clear();

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

            this->move(SparseVector<DT_, IT_>(rows, nnz));

            Memory::TypedView<IT_> idx_view(this->indices_view_w());
            Memory::TypedView<DT_> val_view(this->elements_view_w());

            IT_ * pind(idx_view.get_w());
            DT_ * pval(val_view.get_w());

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
          }
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

        case FileMode::fm_sv:
          this->template _serialize<double, std::uint64_t>(FileMode::fm_sv, file);
          break;

        case FileMode::fm_mtx:
          {

            file << "%%MatrixMarket matrix coordinate real general\n";
            file << this->size() << " " << 1 << " " << this->num_nzes() << "\n";
            file << std::scientific << std::setprecision(Type::Traits<DT_>::format_precision);

            const Index u_elem(this->num_nzes());
            const Memory::TypedView<IT_> idx_view = this->indices_view_r();
            const Memory::TypedView<DT_> val_view = this->elements_view_r();
            for (Index i(0) ; i < u_elem ; ++i)
            {
              file << (idx_view[i]+1) << " " << 1 << " " << val_view[i] << "\n";
            }
          }
          break;

        default:
          XABORTM("Filemode not supported!");
        }
      }

      /// Permute vector according to the given Permutation
      void permute(Adjacency::Permutation & perm)
      {
        if (perm.size() == 0)
          return;

        std::map<IT_, DT_> xm;

        auto ip = perm.inverse();
        const Index * const inv_pos(ip.get_perm_pos());

        Memory::TypedView<DT_> xv(this->elements_view_rw());
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
      bool same_layout(const SparseVector& other) const
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
       * \brief Prints this vector in human readable format to an output stream
       *
       * \param[in] os
       * The target stream.
       *
       * \param[in] print_dense
       * Specifies whether to print the vector as a dense vector or in a compressed format
       */
      void print(std::ostream & os, bool print_dense) const
      {
        if(this->empty())
        {
          os << "[]";
          return;
        }
        const Memory::TypedView<IT_> idx_view(this->indices_view_r());
        const Memory::TypedView<DT_> val_view(this->elements_view_r());
        const Index n = this->size();
        const Index nze = this->num_nzes();
        os << "[";
        Index k = 0;
        for (Index i(0) ; i < nze ; ++i)
        {
          if(print_dense)
          {
            for(; k < idx_view[i]; ++k)
              os << "  " << DT_(0);
          }
          else
            os << "  " << idx_view[i] << ':';

          os << "  " << val_view[i];
          ++k;
        }
        if(print_dense)
        {
          for(; k < n; ++k)
            os << "  " << DT_(0);
        }
        os << "]";
      }

      /**
       * \brief SparseVector streaming operator
       *
       * \param[in] os The target stream.
       * \param[in] b The vector to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & os, const SparseVector & b)
      {
        b.print(os, true);
        return os;
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
