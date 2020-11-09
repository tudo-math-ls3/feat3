// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_DENSE_VECTOR_HPP
  #define KERNEL_LAFEM_DENSE_VECTOR_HPP 1

  // includes, FEAT
  #include <kernel/base_header.hpp>
  #include <kernel/lafem/forward.hpp>
  #include <kernel/util/assertion.hpp>
  #include <kernel/util/type_traits.hpp>
  #include <kernel/util/math.hpp>
  #include <kernel/util/random.hpp>
  #include <kernel/lafem/container.hpp>
  #include <kernel/lafem/dense_vector_blocked.hpp>
  #include <kernel/lafem/edi.hpp>
  #include <kernel/lafem/arch/component_invert.hpp>
  #include <kernel/lafem/arch/dot_product.hpp>
  #include <kernel/lafem/arch/norm.hpp>
  #include <kernel/lafem/arch/max_abs_index.hpp>
  #include <kernel/lafem/arch/min_abs_index.hpp>
  #include <kernel/lafem/arch/max_index.hpp>
  #include <kernel/lafem/arch/min_index.hpp>
  #include <kernel/lafem/arch/scale.hpp>
  #include <kernel/lafem/arch/axpy.hpp>
  #include <kernel/lafem/arch/component_product.hpp>
  #include <kernel/adjacency/permutation.hpp>
  #include <kernel/util/statistics.hpp>
  #include <kernel/util/time_stamp.hpp>

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
     * \tparam Mem_ The \ref FEAT::Mem "memory architecture" to be used.
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
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class DenseVector : public Container<Mem_, DT_, IT_>
    {
    public:
      /**
       * \brief Scatter-Axpy operation for DenseVector
       *
       * \author Peter Zajac
       */
      class ScatterAxpy
      {
      public:
        typedef LAFEM::DenseVector<Mem::Main, DT_, IT_> VectorType;
        typedef Mem::Main MemType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

      private:
        Index _num_entries;
        DT_ * _data;

      public:
        explicit ScatterAxpy(VectorType & vector) :
          _num_entries(vector.size()),
          _data(vector.elements())
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
            _data[dof_idx] += alpha * loc_vec[i];
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
        typedef LAFEM::DenseVector<Mem::Main, DT_, IT_> VectorType;
        typedef Mem::Main MemType;
        typedef DT_ DataType;
        typedef IT_ IndexType;

      private:
        Index _num_entries;
        const DT_ * _data;

      public:
        explicit GatherAxpy(const VectorType & vector) :
          _num_entries(vector.size()),
          _data(vector.elements())
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
            loc_vec[i] += alpha * _data[dof_idx];
          }
        }
      }; // class GatherAxpy

    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our memory architecture type
      typedef Mem_ MemType;
      /// Our value type
      typedef DT_ ValueType;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = DenseVector<Mem2_, DT2_, IT2_>;

      /// this typedef lets you create a vector container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional vector.
       */
      explicit DenseVector() :
        Container<Mem_, DT_, IT_>(0)
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] size_in The size of the created vector.
       * \param[in] pinning True if the memory should be allocated in a pinned manner.
       *
       * \warning Pinned memory allocation is only possible in main memory and needs cuda support.
       *
       * Creates a vector with a given size.
       */
      explicit DenseVector(Index size_in, Pinning pinning = Pinning::disabled) :
        Container<Mem_, DT_, IT_>(size_in)
      {
        if (size_in == Index(0))
          return;

        XASSERTM(! (pinning == Pinning::enabled && (typeid(Mem_) != typeid(Mem::Main))), "Pinned memory allocation only possible in main memory!");

        if (pinning == Pinning::enabled)
        {
  #ifdef FEAT_HAVE_CUDA
          this->_elements.push_back(MemoryPool<Mem::Main>::template allocate_pinned_memory<DT_>(size_in));
  #else
          // no cuda support enabled - we cannot serve and do not need pinned memory support
          this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(size_in));
  #endif
        }
        else
        {
          this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(size_in));
        }

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
        Container<Mem_, DT_, IT_>(size_in)
      {
        if (size_in == Index(0))
          return;

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(size_in));
        this->_elements_size.push_back(size_in);

        MemoryPool<Mem_>::set_memory(this->_elements.at(0), value, size_in);
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
      explicit DenseVector(Index size_in, DT_ * data) :
        Container<Mem_, DT_, IT_>(size_in)
      {
        if (size_in == Index(0))
          return;

        this->_elements.push_back(data);
        this->_elements_size.push_back(size_in);

        for (Index i(0); i < this->_elements.size(); ++i)
          MemoryPool<Mem_>::increase_memory(this->_elements.at(i));
        for (Index i(0); i < this->_indices.size(); ++i)
          MemoryPool<Mem_>::increase_memory(this->_indices.at(i));
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
      explicit DenseVector(const DenseVector & dv_in, Index size_in, Index offset_in) :
        Container<Mem_, DT_, IT_>(size_in)
      {
        XASSERT(size_in > Index(0));
        XASSERTM(size_in + offset_in <= dv_in.size(), "Ranged vector part exceeds original vector size!");

        this->_foreign_memory = true;

        DT_ * te(const_cast<DT_ *>(dv_in.elements()));
        this->_elements.push_back(te + offset_in);
        this->_elements_size.push_back(size_in);
      }

      /**
       * \brief Constructor
       *
       * \param[in] other The source blocked vector
       *
       * Creates a vector from a given source blocked vector
       */
      template <int BS_>
      explicit DenseVector(const DenseVectorBlocked<Mem_, DT_, IT_, BS_> & other) :
        Container<Mem_, DT_, IT_>(other.template size<Perspective::pod>())
      {
        convert(other);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a vector from the given source file.
       */
      explicit DenseVector(FileMode mode, String filename) :
        Container<Mem_, DT_, IT_>(0)
      {
        read_from(mode, filename);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] file The stream that is to be read from.
       *
       * Creates a vector from the given source file.
       */
      explicit DenseVector(FileMode mode, std::istream & file) :
        Container<Mem_, DT_, IT_>(0)
      {
        read_from(mode, file);
      }

      /**
       * \brief Constructor
       *
       * \param[in] input A std::vector, containing the byte.
       *
       * Creates a vector from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit DenseVector(std::vector<char> input) :
        Container<Mem_, DT_, IT_>(0)
      {
        deserialize<DT2_, IT2_>(input);
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
        Container<Mem_, DT_, IT_>(size_in)
      {
        if (size_in == Index(0))
          return;

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(size_in));
        this->_elements_size.push_back(size_in);

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
        Container<Mem_, DT_, IT_>(std::forward<DenseVector>(other))
      {
      }

      /**
       * \brief Destructor
       */
      virtual ~DenseVector()
      {
      }

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
      template <typename Mem2_, typename DT2_, typename IT2_>
      void clone(const DenseVector<Mem2_, DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Deep)
      {
        Container<Mem_, DT_, IT_>::clone(other, clone_mode);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source vector.
       *
       * Use source vector content as content of current vector
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const DenseVector<Mem2_, DT2_, IT2_> & other)
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
      template <typename Mem2_, typename DT2_, typename IT2_, int BS2_>
      void convert(const DenseVectorBlocked<Mem2_, DT2_, IT2_, BS2_> & other)
      {
        this->clear();

        this->_scalar_index.push_back(other.template size<Perspective::pod>());
        this->_elements.push_back(other.get_elements().at(0));
        this->_elements_size.push_back(this->size());

        for (Index i(0); i < this->_elements.size(); ++i)
          MemoryPool<Mem_>::increase_memory(this->_elements.at(i));
        for (Index i(0); i < this->_indices.size(); ++i)
          MemoryPool<Mem_>::increase_memory(this->_indices.at(i));
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
        typedef typename VT_::MemType Mem2_;

        if (std::is_same<Mem_, Mem2_>::value)
        {
          DenseVector vec(a.template size<Perspective::pod>());
          a.set_vec(vec.elements());
          this->assign(vec);
        }
        else
        {
          typename VT_::template ContainerType<Mem_, DT_, IT_> ta;
          ta.convert(a);
          DenseVector vec(ta.template size<Perspective::pod>());
          ta.set_vec(vec.elements());
          this->assign(vec);
        }
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
      std::vector<char> serialize(const LAFEM::SerialConfig& config = LAFEM::SerialConfig())
      {
        return this->template _serialize<DT2_, IT2_>(FileMode::fm_dv, config);
      }

      /**
       * \brief Expand DenseVector to DenseVectorBlocked.
       *
       * Inflate the DenseVector to DenseVectorBlocked by filling each complete block with the input scalar vector entry.
       *
       * \note The resulting DenseVectorBlocked will have size == DenseVector.size() * BlockSize_.
       */
      template <int BlockSize_>
      DenseVectorBlocked<Mem_, DT_, IT_, BlockSize_> inflate_to_blocks()
      {
        DenseVector<Mem::Main, DT_, IT_> main;
        main.convert(*this);
        DenseVectorBlocked<Mem::Main, DT_, IT_, BlockSize_> result_main(main.size());

        const auto * mp = main.elements();
        auto * rmp = result_main.template elements<Perspective::native>();
        for (Index i(0); i < main.size(); ++i)
        {
          rmp[i] = typename decltype(result_main)::ValueType(mp[i]);
        }
        DenseVectorBlocked<Mem_, DT_, IT_, BlockSize_> result;
        result.convert(result_main);
        return result;
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] a The vector to be copied (could be of any format; must have same size).
       */
      template <typename VT_>
      void copy(const VT_ & a)
      {
        XASSERTM(this->template size<Perspective::pod>() == a.template size<Perspective::pod>(), "Vectors have not the same size!");

        typedef typename VT_::MemType Mem2_;
        typedef typename VT_::IndexType IT2_;

        if (std::is_same<Mem_, Mem2_>::value)
        {
          a.set_vec(this->elements());
        }
        else
        {
          typename VT_::template ContainerType<Mem_, DT_, IT2_> ta;
          ta.convert(a);
          ta.set_vec(this->elements());
        }
      }

      /**
       * \brief Performs \f$x \leftarrow this\f$.
       *
       * \param[in] a The target-vector to be copied to (could be of any format; must have same size).
       */
      template <typename VT_>
      void copy_inv(VT_ & a) const
      {
        XASSERTM(this->template size<Perspective::pod>() == a.template size<Perspective::pod>(), "Vectors have not the same size!");

        typedef typename VT_::MemType Mem2_;
        if (std::is_same<Mem_, Mem2_>::value)
        {
          a.set_vec_inv(this->elements());
        }
        else
        {
          DenseVector<Mem2_, DT_, IT_> t_this;
          t_this.convert(*this);
          a.set_vec_inv(t_this.elements());
        }
      }

      /**
       * \brief Read in vector from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */

      void read_from(FileMode mode, String filename)
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
          this->clear();
          this->_scalar_index.push_back(0);

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

          DenseVector<Mem::Main, DT_, IT_> tmp(rows);
          DT_ * pval(tmp.elements());

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
          this->assign(tmp);
          break;
        }
        case FileMode::fm_exp:
        {
          this->clear();
          this->_scalar_index.push_back(0);

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
            for (unsigned long i(0); i < eol; ++i)
            {
              n_z_s.append(1, line[i]);
            }

            DT_ n_z((DT_)atof(n_z_s.c_str()));

            data.push_back(n_z);
          }

          this->_scalar_index.at(0) = Index(data.size());
          this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(Index(data.size())));
          this->_elements_size.push_back(Index(data.size()));
          MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), &data[0], Index(data.size()));
          break;
        }
        case FileMode::fm_dv:
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
      void write_out(FileMode mode, String filename) const
      {
        std::ios_base::openmode bin = std::ofstream::out | std::ofstream::binary;
        if (mode == FileMode::fm_mtx || mode == FileMode::fm_exp)
          bin = std::ofstream::out;
        std::ofstream file(filename.c_str(), bin);
        if (! file.is_open())
          XABORTM("Unable to open Vector file " + filename);

        write_out(mode, file);
        file.close();
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
          DenseVector<Mem::Main, DT_, IT_> temp;
          temp.convert(*this);

          const Index tsize(temp.size());
          file << "%%MatrixMarket matrix array real general" << std::endl;
          file << tsize << " " << 1 << std::endl;

          const DT_ * pval(temp.elements());
          for (Index i(0); i < tsize; ++i, ++pval)
          {
            file << std::scientific << *pval << std::endl;
          }
          break;
        }
        case FileMode::fm_exp:
        {
          DT_ * temp = MemoryPool<Mem::Main>::template allocate_memory<DT_>((this->size()));
          MemoryPool<Mem_>::template download<DT_>(temp, this->_elements.at(0), this->size());

          for (Index i(0); i < this->size(); ++i)
          {
            file << std::scientific << temp[i] << std::endl;
          }
          MemoryPool<Mem::Main>::release_memory(temp);
          break;
        }
        case FileMode::fm_dv:
        case FileMode::fm_binary:
          this->template _serialize<double, std::uint64_t>(FileMode::fm_dv, file);
          break;
        default:
          XABORTM("Filemode not supported!");
        }
      }

      /**
       * \brief Get a pointer to the data array.
       *
       * \returns Pointer to the data array.
       */
      template <Perspective = Perspective::native>
      DT_ * elements()
      {
        if (this->_elements.size() == 0)
          return nullptr;

        return this->_elements.at(0);
      }

      template <Perspective = Perspective::native>
      DT_ const * elements() const
      {
        if (this->_elements.size() == 0)
          return nullptr;

        return this->_elements.at(0);
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
        ASSERT(index < this->size());
        return MemoryPool<Mem_>::get_element(this->_elements.at(0), index);
      }

      /**
       * \brief Set specific vector element.
       *
       * \param[in] index The index of the vector element.
       * \param[in] value The value to be set.
       */
      void operator()(Index index, DT_ value)
      {
        ASSERT(index < this->size());
        MemoryPool<Mem_>::set_memory(this->_elements.at(0) + index, value);
      }

      /**
       * \brief Create temporary object for direct data manipulation.
       * \warning Be aware, that any synchronization only takes place, when the object is destroyed!
       *
       * \param[in] index The index of the vector element.
       */
      EDI<Mem_, DT_> edi(Index index)
      {
        EDI<Mem_, DT_> t(MemoryPool<Mem_>::get_element(this->_elements.at(0), index), this->_elements.at(0) + index);
        return t;
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

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The vector to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      template <typename Mem2_>
      void copy(const DenseVector<Mem2_, DT_, IT_> & x, bool full = false)
      {
        this->_copy_content(x, full);
      }

      ///@name Linear algebra operations
      ///@{
      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x + y\f$
       *
       * \param[in] x The first summand vector to be scaled.
       * \param[in] y The second summand vector
       * \param[in] alpha A scalar to multiply x with.
       */
      void axpy(
        const DenseVector & x,
        const DenseVector & y,
        const DT_ alpha = DT_(1))
      {
        XASSERTM(x.size() == y.size(), "Vector size does not match!");
        XASSERTM(x.size() == this->size(), "Vector size does not match!");

        if (Math::abs(alpha) < Math::eps<DT_>())
        {
          this->copy(y);
          //y.scale(beta);
          return;
        }

        TimeStamp ts_start;

        Statistics::add_flops(this->size() * 2);
        Arch::Axpy<Mem_>::dv(this->elements(), alpha, x.elements(), y.elements(), this->size());

        TimeStamp ts_stop;
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

        TimeStamp ts_start;

        Statistics::add_flops(this->size());
        Arch::ComponentProduct<Mem_>::value(this->elements(), x.elements(), y.elements(), this->size());

        TimeStamp ts_stop;
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

        TimeStamp ts_start;

        Statistics::add_flops(this->size());
        Arch::ComponentInvert<Mem_>::value(this->elements(), x.elements(), alpha, this->size());

        TimeStamp ts_stop;
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

        TimeStamp ts_start;

        Statistics::add_flops(this->size());
        Arch::Scale<Mem_>::value(this->elements(), x.elements(), alpha, this->size());

        TimeStamp ts_stop;
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

        TimeStamp ts_start;

        Statistics::add_flops(this->size() * 3);
        DataType result = Arch::TripleDotProduct<Mem_>::value(this->elements(), x.elements(), y.elements(), this->size());

        TimeStamp ts_stop;
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

        TimeStamp ts_start;

        Statistics::add_flops(this->size() * 2);
        DataType result = Arch::DotProduct<Mem_>::value(this->elements(), x.elements(), this->size());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Calculates and returns the euclid norm of this vector.
       *
       */
      DT_ norm2() const
      {
        TimeStamp ts_start;
        Statistics::add_flops(this->size() * 2);

        DT_ result = Arch::Norm2<Mem_>::value(this->elements(), this->size());

        TimeStamp ts_stop;
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
        TimeStamp ts_start;

        Index max_abs_index = Arch::MaxAbsIndex<Mem_>::value(this->template elements<Perspective::pod>(), this->template size<Perspective::pod>());
        ASSERT(max_abs_index < this->template size<Perspective::pod>());
        DT_ result;
        MemoryPool<Mem_>::template download<DT_>(&result, this->template elements<Perspective::pod>() + max_abs_index, 1);
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

        Index min_abs_index = Arch::MinAbsIndex<Mem_>::value(this->template elements<Perspective::pod>(), this->template size<Perspective::pod>());
        ASSERT(min_abs_index < this->template size<Perspective::pod>());
        DT_ result;
        MemoryPool<Mem_>::template download<DT_>(&result, this->template elements<Perspective::pod>() + min_abs_index, 1);
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

        Index max_index = Arch::MaxIndex<Mem_>::value(this->template elements<Perspective::pod>(), this->template size<Perspective::pod>());
        ASSERT(max_index < this->template size<Perspective::pod>());
        DT_ result;
        MemoryPool<Mem_>::template download<DT_>(&result, this->template elements<Perspective::pod>() + max_index, 1);

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

        Index min_index = Arch::MinIndex<Mem_>::value(this->template elements<Perspective::pod>(), this->template size<Perspective::pod>());
        ASSERT(min_index < this->template size<Perspective::pod>());
        DT_ result;
        MemoryPool<Mem_>::template download<DT_>(&result, this->template elements<Perspective::pod>() + min_index, 1);

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      ///@}

      /// Permutate vector according to the given Permutation
      void permute(Adjacency::Permutation & perm)
      {
        if (perm.size() == 0)
          return;

        XASSERTM(perm.size() == this->size(), "Container size does not match permutation size");

        DenseVector<Mem::Main, DT_, IT_> local;
        local.convert(*this);
        perm(local.elements());
        this->assign(local);
      }

      /// \cond internal
      /// Writes the vector-entries in an allocated array
      void set_vec(DT_ * const pval_set) const
      {
        MemoryPool<Mem_>::copy(pval_set, this->elements(), this->size());
      }

      /// Writes data of an array in the vector
      void set_vec_inv(const DT_ * const pval_set)
      {
        MemoryPool<Mem_>::copy(this->elements(), pval_set, this->size());
      }
      /// \endcond

      /**
       * \brief DenseVector comparison operator
       *
       * \param[in] a A vector to compare with.
       * \param[in] b A vector to compare with.
       */
      template <typename Mem2_>
      friend bool operator==(const DenseVector & a, const DenseVector<Mem2_, DT_, IT_> & b)
      {
        if (a.size() != b.size())
          return false;
        if (a.get_elements().size() != b.get_elements().size())
          return false;
        if (a.get_indices().size() != b.get_indices().size())
          return false;

        if (a.size() == 0 && b.size() == 0 && a.get_elements().size() == 0 && b.get_elements().size() == 0)
          return true;

        bool ret(true);

        DT_ * ta;
        DT_ * tb;

        if (std::is_same<Mem::Main, Mem_>::value)
        {
          ta = const_cast<DT_ *>(a.elements());
        }
        else
        {
          ta = new DT_[a.size()];
          MemoryPool<Mem_>::template download<DT_>(ta, a.elements(), a.size());
        }

        if (std::is_same<Mem::Main, Mem2_>::value)
        {
          tb = const_cast<DT_ *>(b.elements());
        }
        else
        {
          tb = new DT_[b.size()];
          MemoryPool<Mem2_>::template download<DT_>(tb, b.elements(), b.size());
        }

        for (Index i(0); i < a.size(); ++i)
        {
          if (ta[i] != tb[i])
          {
            ret = false;
            break;
          }
        }

        if (! std::is_same<Mem::Main, Mem_>::value)
          delete[] ta;
        if (! std::is_same<Mem::Main, Mem2_>::value)
          delete[] tb;

        return ret;
      }

      /**
       * \brief DenseVector streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The vector to be streamed.
       */
      friend std::ostream & operator<<(std::ostream & lhs, const DenseVector & b)
      {
        lhs << "[";
        for (Index i(0); i < b.size(); ++i)
        {
          lhs << "  " << b(i);
        }
        lhs << "]";

        return lhs;
      }
    }; // class DenseVector<...>

  #ifdef FEAT_EICKT
    extern template class DenseVector<Mem::Main, float, unsigned int>;
    extern template class DenseVector<Mem::Main, double, unsigned int>;
    extern template class DenseVector<Mem::Main, float, unsigned long>;
    extern template class DenseVector<Mem::Main, double, unsigned long>;
    #ifdef FEAT_HAVE_CUDA
    extern template class DenseVector<Mem::CUDA, float, unsigned int>;
    extern template class DenseVector<Mem::CUDA, double, unsigned int>;
    extern template class DenseVector<Mem::CUDA, float, unsigned long>;
    extern template class DenseVector<Mem::CUDA, double, unsigned long>;
    #endif
  #endif

  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_DENSE_VECTOR_HPP
