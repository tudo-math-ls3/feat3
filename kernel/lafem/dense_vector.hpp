// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/likwid_marker.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/random.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/adjacency/permutation.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/arch/axpy_dense.hpp>
#include <kernel/lafem/arch/component_invert_dense.hpp>
#include <kernel/lafem/arch/component_product_dense.hpp>
#include <kernel/lafem/arch/dot_product_dense.hpp>
#include <kernel/lafem/arch/max_rel_diff_dense.hpp>
#include <kernel/lafem/arch/min_max_index_dense.hpp>
#include <kernel/lafem/arch/min_max_value_dense.hpp>
#include <kernel/lafem/arch/norm2_dense.hpp>
#include <kernel/lafem/arch/scale_dense.hpp>
#include <kernel/lafem/arch/triple_dot_product_dense.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Dense data vector class template.
     *
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indextype to be used.
     *
     * This class represents a vector of continuous data in memory. \n \n
     * Data survey: \n
     * _elements[0]: raw number values \n
     * _scalar_index[0]: container size
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename DT_, typename IT_ = Index>
    class DenseVector : public Container<DT_, IT_>
    {
    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our value type
      typedef DT_ ValueType;
      /// Our 'base' class type
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = DenseVector<DT2_, IT2_>;

      /// this typedef lets you create a vector container with different Data and Index types
      template <typename DataType2_, typename IndexType2_>
      using ContainerTypeByDI = ContainerType<DataType2_, IndexType2_>;

    protected:
      Index& _size()
      {
        return this->_elements_size.at(0);
      }

    public:
      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional vector.
       */
      DenseVector() = default;

      /**
       * \brief Constructor
       *
       * \param[in] size_in The size of the created vector.
       *
       * Creates a vector with a given size.
       */
      explicit DenseVector(Index size_in)
      {
        if (size_in == Index(0))
          return;

        this->_elements.push_back(Memory::Arbiter(sizeof(DT_) * size_in));
        this->_elements_size.push_back(size_in);
      }

      /**
       * \brief Constructor
       *
       * \param[in] size_in The size of the created vector.
       * \param[in] value The value, each element will be set to.
       *
       * Creates a vector with given size and value.
       */
      explicit DenseVector(Index size_in, DT_ value) :
        DenseVector(size_in)
      {
        this->format(value);
      }

      /**
       * \brief Constructor
       *
       * \param[in] size_in The size of the created vector.
       * \param[in] data An array containing the value data.
       *
       * Creates a vector with given size and given data.
       *
       * \note The array must be allocated by FEAT's own memory pool
       *
       * \note Obviously, the pointer must point to data in the same memory.
       *
       * \warning If you use a data pointer from another container, both containers will be able to modify the array and will effect the other container,
       * as they share one and the same array.
       */
      explicit DenseVector(Index size_in, Memory::Arbiter&& data)
      {
        if (size_in == Index(0))
          return;

        this->_elements.push_back(std::forward<Memory::Arbiter>(data));
        this->_elements_size.push_back(size_in);
      }

      /**
       * \brief Constructor
       *
       * \param[in] dv_in The source DenseVector
       * \param[in] size_in The size of the created vector range.
       * \param[in] offset_in The starting element of the created vector range in relation to the source vector.
       *
       * Creates a vector range from a given DenseVector
       *
       * \note The created DenseVector has no own memory management nor own allocated memory and should be used carefully!
       */
      explicit DenseVector(const DenseVector & dv_in, Index size_in, Index offset_in)
      {
        XASSERT(size_in > Index(0));
        XASSERTM(size_in + offset_in <= dv_in.size(), "Ranged vector part exceeds original vector size!");

        this->_elements.push_back(dv_in._elements.front().attach(offset_in * sizeof(DT_), size_in * sizeof(DT_)));
        this->_elements_size.push_back(size_in);
      }

      /**
       * \brief Constructor
       *
       * \param[in] other The source blocked vector
       *
       * Creates a vector from a given source blocked vector
       */
      template <int block_size_>
      explicit DenseVector(const DenseVectorBlocked<DT_, IT_, block_size_> & other)
      {
        this->convert(other);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a vector from the given source file.
       */
      explicit DenseVector(FileMode mode, String filename)
      {
        this->read_from(mode, filename);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The stream that is to be read from.
       *
       * Creates a vector from the given source file.
       */
      explicit DenseVector(FileMode mode, std::istream& file)
      {
        this->read_from(mode, file);
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte.
       *
       * Creates a vector from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit DenseVector(std::vector<char> input)
      {
        this->deserialize<DT2_, IT2_>(input);
      }

      /**
       * \brief Constructor
       *
       * \param[in] rng The random number generator.
       * \param[in] size_in The vector size.
       * \param[in] min Lower rng bound.
       * \param[in] max Upper rng bound.
       *
       * Creates a vector from the given random number generator.
       */
      explicit DenseVector(Random & rng, Index size_in, DataType min, DataType max) :
        DenseVector(size_in)
      {
        this->format(rng, min, max);
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source vector.
       *
       * Moves another vector to this vector.
       */
      DenseVector(DenseVector && other) :
        Container<DT_, IT_>(std::forward<DenseVector>(other))
      {
      }

      /**
       * \brief Destructor
       */
      virtual ~DenseVector() = default;

      /**
       * \brief Assignment move operator
       *
       * \param[in] other The source vector.
       *
       * Moves another vector to the target vector.
       */
      DenseVector & operator=(DenseVector && other)
      {
        this->move(std::forward<DenseVector>(other));
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
      DenseVector clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        DenseVector t;
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
      template <typename DT2_, typename IT2_>
      void clone(const DenseVector<DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Deep)
      {
        Container<DT_, IT_>::clone(other, clone_mode);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source vector.
       *
       * Use source vector content as content of current vector
       */
      template <typename DT2_, typename IT2_>
      void convert(const DenseVector<DT2_, IT2_> & other)
      {
        this->assign(other);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source vector.
       *
       * Use source vector content as content of current vector
       */
      template <typename DT2_, typename IT2_, int BS2_>
      void convert(const DenseVectorBlocked< DT2_, IT2_, BS2_> & other)
      {
        this->assign(other);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] a The input vector.
       *
       * Converts any vector to DenseVector-format
       */
      template <typename VT_>
      void convert(const VT_ & a)
      {
        DenseVector vec(a.size_raw());
        vec.copy(a);
        this->move(std::move(vec));
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] a The vector to be copied (could be of any format; must have same size).
       */
      template <typename VT_>
      void copy(const VT_ & a)
      {
        XASSERTM(this->size_raw() == a.size_raw(), "Vectors have not the same size!");
        Memory::TypedView<DT_> x(this->elements_view_w());
        a.get_values(x.get_w());
      }

      /**
       * \brief Performs \f$x \leftarrow this\f$.
       *
       * \param[in] a The target-vector to be copied to (could be of any format; must have same size).
       */
      template <typename VT_>
      void copy_to(VT_ & a) const
      {
        XASSERTM(this->size_raw() == a.size_raw(), "Vectors have not the same size!");
        Memory::TypedView<DT_> x(this->elements_view_r());
        a.set_values(x.get_r());
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
        this->template _deserialize<DT2_, IT2_>(FileMode::fm_dv, input);
      }

      /**
       * \brief Serialization of complete container entity.
       *
       * \param[in] config LAFEM::SerialConfig, a struct describing the serialize configuration.
       * \note the corresponding configure flags 'zlib' and/or 'zfp' need to be added in the build-id at the configure call.
       *
       * \returns A std::vector, containing the byte array.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAT::LAFEM::Container::_serialize for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialize(const LAFEM::SerialConfig& config = LAFEM::SerialConfig()) const
      {
        return this->template _serialize<DT2_, IT2_>(FileMode::fm_dv, config);
      }

      /**
       * \brief Expand DenseVector to DenseVectorBlocked.
       *
       * Inflate the DenseVector to DenseVectorBlocked by filling each complete block with the input scalar vector entry.
       *
       * \note The resulting DenseVectorBlocked will have size == DenseVector.size() * block_size_.
       */
      template <int block_size_>
      DenseVectorBlocked<DT_, IT_, block_size_> inflate_to_blocks()
      {
        const Index n = this->size();
        DenseVectorBlocked<DT_, IT_, block_size_> result(this->size());

        auto vb = result.elements_view_w();
        const auto vx = this->elements_view_r();
        for(Index i(0); i < n; ++i)
          vb[i] = vx(i);

        return result;
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
        if (mode == FileMode::fm_mtx)
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
      void read_from(FileMode mode, std::istream & file)
      {
        switch (mode)
        {
        case FileMode::fm_mtx:
        {
          Index rows;
          String line;
          std::getline(file, line);
          if (line.find("%%MatrixMarket matrix array real general") == String::npos)
          {
            XABORTM("Input-file is not a compatible mtx-vector-file");
          }
          while (! file.eof())
          {
            std::getline(file, line);
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
            XASSERTM(cols == 1, "Input-file is no dense-vector-file");
          }

          DenseVector<DT_, IT_> tmp(rows);
          Memory::TypedView<DT_> tmp_view(tmp.elements_view_w());
          DT_* pval = tmp_view.get_w();

          while (! file.eof())
          {
            std::getline(file, line);
            if (file.eof())
              break;

            String::size_type begin(line.find_first_not_of(" "));
            line.erase(0, begin);
            String::size_type end(line.find_first_of(" "));
            String sval(line, 0, end);
            DT_ tval((DT_)atof(sval.c_str()));

            *pval = tval;
            ++pval;
          }
          tmp_view.release();
          this->move(std::move(tmp));
          break;
        }
        case FileMode::fm_exp:
        {
          std::vector<DT_> data;

          while (! file.eof())
          {
            std::string line;
            std::getline(file, line);
            if (line.find("#", 0) < line.npos)
              continue;
            if (file.eof())
              break;

            std::string n_z_s;

            std::string::size_type first_digit(line.find_first_not_of(" "));
            line.erase(0, first_digit);
            std::string::size_type eol(line.length());
            for (std::string::size_type i(0); i < eol; ++i)
            {
              n_z_s.append(1, line[i]);
            }

            DT_ n_z((DT_)atof(n_z_s.c_str()));

            data.push_back(n_z);
          }

          DenseVector<DT_, IT_> tmp(data.size());
          {
            Memory::TypedView<DT_> tmp_view(tmp.elements_view_w());
            tmp_view.convert_from(data.data());
          }

          this->move(std::move(tmp));
          break;
        }

        case FileMode::fm_dv:
          [[fallthrough]];

        case FileMode::fm_binary:
          this->template _deserialize<double, std::uint64_t>(FileMode::fm_dv, file);
          break;

        default:
          XABORTM("Filemode not supported!");
        }
      }

      /**
       * \brief Write out vector to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the vector shall be stored.
       */
      void write_out(FileMode mode, const String& filename) const
      {
        std::ios_base::openmode bin = std::ofstream::out | std::ofstream::binary;
        if (mode == FileMode::fm_mtx || mode == FileMode::fm_exp)
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
      void write_out(FileMode mode, std::ostream & file) const
      {
        switch (mode)
        {
        case FileMode::fm_mtx:
        {
          const Index tsize(this->size_raw());
          file << "%%MatrixMarket matrix array real general\n";
          file << tsize << " " << 1 << "\n";
          file << std::scientific << std::setprecision(Type::Traits<DT_>::format_precision);
          const Memory::TypedView<DT_> elem_view(this->elements_view_r());
          for (Index i(0); i < tsize; ++i)
          {
            file << elem_view(i) << "\n";
          }
          break;
        }

        case FileMode::fm_exp:
        {
          const Index tsize(this->size_raw());
          file << std::scientific << std::setprecision(Type::Traits<DT_>::format_precision);
          const Memory::TypedView<DT_> elem_view(this->elements_view_r());
          for (Index i(0); i < tsize; ++i)
          {
            file << elem_view(i) << "\n";
          }
          break;
        }

        case FileMode::fm_dv:
          [[fallthrough]];

        case FileMode::fm_binary:
          this->template _serialize<double, std::uint64_t>(FileMode::fm_dv, file);
          break;

        default:
          XABORTM("Filemode not supported!");
        }
      }

      /// Returns the size of this vector, i.e. its number of entries
      Index size() const
      {
        return this->_elements_size.empty() ? Index(0) : this->_elements_size.at(0);
      }

      /// Returns the size of this vector in scalar data type entries
      Index size_raw() const
      {
        return this->size();
      }

      /// Checks whether the vector is empty, i.e. if it has size 0
      bool empty() const
      {
        return this->_elements_size.empty() || (this->_elements_size.at(0) <= Index(0));
      }

      /// Returns a reference to the element array arbiter
      Memory::Arbiter& elements_arbiter()
      {
        return this->_elements.front();
      }

      /// Returns a reference to the element array arbiter
      const Memory::Arbiter& elements_arbiter() const
      {
        return this->_elements.front();
      }

      Memory::TypedView<DT_> elements_view_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::read));
      }

      Memory::TypedView<DT_> elements_view_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::write));
      }

      Memory::TypedView<DT_> elements_view_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<DT_> elements_view(Memory::Location loc, Memory::Access acc)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, acc));
      }

      Memory::TypedView<DT_> elements_view_raw_r(Memory::Location loc = Memory::Location::main) const
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::read));
      }

      Memory::TypedView<DT_> elements_view_raw_w(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::write));
      }

      Memory::TypedView<DT_> elements_view_raw_rw(Memory::Location loc = Memory::Location::main)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, Memory::Access::read_write));
      }

      Memory::TypedView<DT_> elements_view_raw(Memory::Location loc, Memory::Access acc)
      {
        if(this->_elements.empty())
          return Memory::TypedView<DT_>();
        return Memory::TypedView<DT_>(this->_elements.at(0).view(loc, acc));
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "DenseVector";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The vector to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      void copy(const DenseVector & x, bool full = false)
      {
        this->_copy_content(x, full);
      }

      ///@name Linear algebra operations
      ///@{
      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x + this\f$
       *
       * \param[in] x The first summand vector to be scaled.
       * \param[in] y The second summand vector
       * \param[in] alpha A scalar to multiply x with.
       */
      void axpy(
        const DenseVector & x,
        const DT_ alpha = DT_(1))
      {
        XASSERTM(x.size() == this->size(), "Vector size does not match!");

        if(this->empty())
          return;

        TimeStamp ts_start;
        FEAT_KERNEL_MARKER_START("DV_axpy");

        Arch::AxpyDense::template exec<DT_>(this->elements_arbiter(), alpha, x.elements_arbiter(), this->size());

        FEAT_KERNEL_MARKER_STOP("DV_axpy");
        TimeStamp ts_stop;

        Statistics::add_flops(this->size() * 2);
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this_i \leftarrow x_i \cdot y_i\f$
       *
       * \param[in] x The first factor.
       * \param[in] y The second factor.
       */
      void component_product(const DenseVector & x, const DenseVector & y)
      {
        XASSERTM(this->size() == x.size(), "Vector size does not match!");
        XASSERTM(this->size() == y.size(), "Vector size does not match!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::ComponentProductDense::template exec<DT_>(this->elements_arbiter(), x.elements_arbiter(), y.elements_arbiter(), this->size());

        TimeStamp ts_stop;

        Statistics::add_flops(this->size());
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Performs \f$ this_i \leftarrow \alpha / x_i \f$
       *
       * \param[in] x
       * The vector whose components serve as denominators.
       *
       * \param[in] alpha
       * The nominator.
       */
      void component_invert(const DenseVector & x, const DT_ alpha = DT_(1))
      {
        XASSERTM(this->size() == x.size(), "Vector size does not match!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::ComponentInvertDense::template exec<DT_>(this->elements_arbiter(), x.elements_arbiter(), alpha, this->size());

        TimeStamp ts_stop;

        Statistics::add_flops(this->size());
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x \f$
       *
       * \param[in] x The vector to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const DenseVector & x, const DT_ alpha)
      {
        XASSERTM(x.size() == this->size(), "Vector size does not match!");

        if(this->empty())
          return;

        TimeStamp ts_start;

        Arch::ScaleDense::template exec<DT_>(this->elements_arbiter(), x.elements_arbiter(), alpha, this->size());

        TimeStamp ts_stop;

        Statistics::add_flops(this->size());
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$result \leftarrow x^T \mathrm{diag}(this) y \f$
       *
       * \param[in] x The first vector.
       *
       * \param[in] y The second vector.
       *
       * \return The computed triple dot product.
       */
      DataType triple_dot(const DenseVector & x, const DenseVector & y) const
      {
        XASSERTM(x.size() == this->size(), "Vector size does not match!");
        XASSERTM(y.size() == this->size(), "Vector size does not match!");
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DataType result = Arch::TripleDotProductDense::template exec<DT_>(this->elements_arbiter(), x.elements_arbiter(), y.elements_arbiter(), this->size());

        TimeStamp ts_stop;

        Statistics::add_flops(this->size() * 3);
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Calculate \f$result \leftarrow this \cdot x\f$
       *
       * \param[in] x The other vector.
       *
       * \return The computed dot product.
       */
      DataType dot(const DenseVector & x) const
      {
        XASSERTM(x.size() == this->size(), "Vector size does not match!");
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DataType result = Arch::DotProductDense::template exec<DT_>(this->elements_arbiter(), x.elements_arbiter(), this->size());

        TimeStamp ts_stop;

        Statistics::add_flops(this->size() * 2);
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Calculates and returns the euclid norm of this vector.
       *
       */
      DT_ norm2() const
      {
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DataType result = Arch::Norm2Dense::template exec<DT_>(this->elements_arbiter(), this->size());

        TimeStamp ts_stop;

        Statistics::add_flops(this->size() * 2);
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Calculates and returns the squared euclid norm of this vector.
       *
       * \return The computed norm.
       *
       */
      DT_ norm2sqr() const
      {
        // fallback
        return Math::sqr(this->norm2());
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

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->elements_arbiter(), this->size(), false, true);

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

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->elements_arbiter(), this->size(), true, true);

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

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->elements_arbiter(), this->size(), false, false);

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

        DT_ result = Arch::MinMaxValueDense::template exec<DT_>(this->elements_arbiter(), this->size(), true, false);

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
      DT_ max_rel_diff(const DenseVector & x) const
      {
        XASSERTM(x.size() == this->size(), "vector size mismatch!");
        XASSERT(!this->empty());

        TimeStamp ts_start;

        DT_ result = Arch::MaxRelDiffDense::template exec<DT_>(this->elements_arbiter(), x.elements_arbiter(), this->size());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Checks if the structural layout of this vector matches that of another vector.
       * This excludes comparison of the actual data values.
       *
       * \param[in] x The vector to compare this vector to
       *
       * \returns true if the layouts match, false otherwise.
       */
      bool same_layout(const DenseVector& x) const
      {
        return this->size() == x.size();
      }

      ///@}

      /// Permute vector according to the given Permutation
      void permute(Adjacency::Permutation & perm)
      {
        if(perm.empty())
          return;

        XASSERTM(perm.size() == this->size(), "Container size does not match permutation size");
        XASSERT(!this->empty());

        Memory::TypedView<DT_> view(this->elements_view(Memory::Location::main, Memory::Access::read_write));

        perm.apply(view.get_w());
      }

      /// \cond internal
      /**
       * \brief Extracts the values of this vector
       *
       * \param[out] pvals
       * A \transient array that receives the values
       *
       * \returns The number of values extracted
       */
      Index get_values(DT_ * const pval_set) const
      {
        this->_elements.front().copy_to(pval_set, Memory::Location::main);
        return this->size_raw();
      }

      /**
       * \brief Overwrites the values of this vector
       *
       * \param[out] pvals
       * A \transient array containing the values to write to the vector
       *
       * \returns The number of values written
       */
      Index set_values(const DT_ * const pval_set)
      {
        this->_elements.front().copy(pval_set, Memory::Location::main);
        return this->size_raw();
      }
      /// \endcond

      /**
       * \brief DenseVector streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The vector to be streamed.
       */
      friend std::ostream & operator<<(std::ostream & lhs, const DenseVector & b)
      {
        const Memory::TypedView<DT_> view(b.elements_view_r());
        const Index n = b.size();
        lhs << "[";
        for (Index i(0); i < n; ++i)
        {
          lhs << "  " << view[i];
        }
        lhs << "]";

        return lhs;
      }

      /**
       * \brief Scatter-Axpy operation for DenseVector
       *
       * \author Peter Zajac
       */
      class ScatterAxpy
      {
      public:
        typedef LAFEM::DenseVector<DT_, IT_> VectorType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

      private:
        const Index _num_entries;
        Memory::TypedView<DT_> _data_view;

      public:
        explicit ScatterAxpy(VectorType & vector) :
          _num_entries(vector.size()),
          _data_view(vector.elements_view(Memory::Location::main, Memory::Access::read_write | Memory::Access::overlap))
        {
        }

        template <typename LocalVector_, typename Mapping_>
        void operator()(const LocalVector_ & loc_vec, const Mapping_ & mapping, DT_ alpha = DT_(1))
        {
          // loop over all local entries
          for (int i(0); i < mapping.get_num_local_dofs(); ++i)
          {
            // get dof index
            Index dof_idx = mapping.get_index(i);
            ASSERT(dof_idx < _num_entries);

            // update vector entry
            _data_view[dof_idx] += alpha * loc_vec[i];
          }
        }
      }; // class ScatterAxpy

      /**
       * \brief Gather-Axpy operation for DenseVector
       *
       * \author Peter Zajac
       */
      class GatherAxpy
      {
      public:
        typedef LAFEM::DenseVector<DT_, IT_> VectorType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

      private:
        const Index _num_entries;
        const Memory::TypedView<DT_> _data_view;

      public:
        explicit GatherAxpy(const VectorType & vector) :
          _num_entries(vector.size()),
          _data_view(vector.elements_view_r())
        {
        }

        template <typename LocalVector_, typename Mapping_>
        void operator()(LocalVector_ & loc_vec, const Mapping_ & mapping, DT_ alpha = DT_(1))
        {
          // loop over all local entries
          for (int i(0); i < mapping.get_num_local_dofs(); ++i)
          {
            // get dof index
            Index dof_idx = mapping.get_index(i);
            ASSERT(dof_idx < _num_entries);

            // update local vector data
            loc_vec[i] += alpha * _data_view(dof_idx);
          }
        }
      }; // class GatherAxpy
    }; // class DenseVector<...>

  #ifdef FEAT_EICKT
    extern template class DenseVector<float, std::uint32_t>;
    extern template class DenseVector<double, std::uint32_t>;
    extern template class DenseVector<float, std::uint64_t>;
    extern template class DenseVector<double, std::uint64_t>;
  #endif

  } // namespace LAFEM
} // namespace FEAT
