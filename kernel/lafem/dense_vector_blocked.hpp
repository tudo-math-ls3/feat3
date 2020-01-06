// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_DENSE_VECTOR_BLOCKED_HPP
  #define KERNEL_LAFEM_DENSE_VECTOR_BLOCKED_HPP 1

  // includes, FEAT
  #include <kernel/base_header.hpp>
  #include <kernel/lafem/forward.hpp>
  #include <kernel/util/assertion.hpp>
  #include <kernel/util/type_traits.hpp>
  #include <kernel/util/math.hpp>
  #include <kernel/util/random.hpp>
  #include <kernel/lafem/container.hpp>
  #include <kernel/lafem/dense_vector.hpp>
  #include <kernel/lafem/arch/dot_product.hpp>
  #include <kernel/lafem/arch/norm.hpp>
  #include <kernel/lafem/arch/scale.hpp>
  #include <kernel/lafem/arch/axpy.hpp>
  #include <kernel/lafem/arch/component_product.hpp>
  #include <kernel/lafem/arch/component_invert.hpp>
  #include <kernel/lafem/arch/max_element.hpp>
  #include <kernel/util/tiny_algebra.hpp>
  #include <kernel/util/statistics.hpp>
  #include <kernel/util/time_stamp.hpp>
  #include <kernel/adjacency/permutation.hpp>

  #include <iostream>
  #include <fstream>
  #include <string>
  #include <stdint.h>

namespace FEAT
{
  namespace LAFEM
  {
    /// \cond internal
    namespace Intern
    {
      template <typename DT_, int BlockSize_, Perspective perspective_>
      struct DenseVectorBlockedPerspectiveHelper
      {
        typedef Tiny::Vector<DT_, BlockSize_> Type;
      };

      template <typename DT_, int BlockSize_>
      struct DenseVectorBlockedPerspectiveHelper<DT_, BlockSize_, Perspective::pod>
      {
        typedef DT_ Type;
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief Blocked Dense data vector class template.
     *
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indextype to be used.
     * \tparam BlockSize_ The size of the represented blocks
     *
     * This class represents a vector of continuous data in memory.\n
     * Logical, the data are organized in small blocks of BlockSize_ length.\n\n
     * Data survey: \n
     * _elements[0]: raw number values \n
     * _scalar_index[0]: container size - aka block count
     * _scalar_index[1]: boolean flag, signaling DenseVectorRange usage
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename IT_, int BlockSize_>
    class DenseVectorBlocked : public Container<Mem_, DT_, IT_>
    {
    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our memory architecture type
      typedef Mem_ MemType;
      /// Our size of a single block
      static constexpr int BlockSize = BlockSize_;
      /// Our value type
      typedef Tiny::Vector<DT_, BlockSize_> ValueType;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = DenseVectorBlocked<Mem2_, DT2_, IT2_, BlockSize_>;

      /// this typedef lets you create a vector container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

      /**
       * \brief Scatter-Axpy operation for DenseVectorBlocked
       *
       * \author Peter Zajac
       */
      class ScatterAxpy
      {
      public:
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, BlockSize_> VectorType;
        typedef Mem::Main MemType;
        typedef DT_ DataType;
        typedef IT_ IndexType;
        typedef Tiny::Vector<DT_, BlockSize_> ValueType;

      private:
        IT_ _num_entries;
        ValueType * _data;

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

            // update vector data
            _data[dof_idx] += alpha * loc_vec[i];
          }
        }
      }; // class ScatterAxpy

      /**
       * \brief Gather-Axpy operation for DenseVectorBlocked
       *
       * \author Peter Zajac
       */
      class GatherAxpy
      {
      public:
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, BlockSize_> VectorType;
        typedef Mem::Main MemType;
        typedef DT_ DataType;
        typedef IT_ IndexType;
        typedef Tiny::Vector<DT_, BlockSize_> ValueType;

      private:
        IT_ _num_entries;
        const ValueType * _data;

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
            loc_vec[i].axpy(alpha, _data[dof_idx]);
          }
        }
      }; // class GatherAxpy

    private:
      Index & _size()
      {
        return this->_scalar_index.at(0);
      }

    public:
      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional vector.
       */
      explicit DenseVectorBlocked() :
        Container<Mem_, DT_, IT_>(0)
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] size_in
       * The size of the created vector. aka block count
       * \param[in] pinning True if the memory should be allocated in a pinned manner.
       *
       * \warning Pinned memory allocation is only possible in main memory and needs cuda support.
       *
       * Creates a vector with a given block count.
       */
      explicit DenseVectorBlocked(Index size_in, Pinning pinning = Pinning::disabled) :
        Container<Mem_, DT_, IT_>(size_in)
      {
        if (size_in == Index(0))
          return;

        XASSERTM(! (pinning == Pinning::enabled && (typeid(Mem_) != typeid(Mem::Main))), "Pinned memory allocation only possible in main memory!");

        if (pinning == Pinning::enabled)
        {
  #ifdef FEAT_HAVE_CUDA
          this->_elements.push_back(MemoryPool<Mem::Main>::template allocate_pinned_memory<DT_>(size<Perspective::pod>()));
  #else
          // no cuda support enabled - we cannot serve and do not need pinned memory support
          this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(size<Perspective::pod>()));
  #endif
        }
        else
        {
          this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(size<Perspective::pod>()));
        }

        this->_elements_size.push_back(size<Perspective::pod>());
      }

      /**
       * \brief Constructor
       *
       * \param[in] size_in
       * The size of the created vector.
       *
       * \param[in] value
       * The value, each element will be set to.
       *
       * Creates a vector with given size and value.
       */
      explicit DenseVectorBlocked(Index size_in, DT_ value) :
        Container<Mem_, DT_, IT_>(size_in)
      {
        if (size_in == Index(0))
          return;

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(size<Perspective::pod>()));
        this->_elements_size.push_back(size<Perspective::pod>());

        MemoryPool<Mem_>::set_memory(this->_elements.at(0), value, size<Perspective::pod>());
      }

      /**
       * \brief Constructor
       *
       * \param[in] size_in The block count of the created vector.
       * \param[in] data An array containing the value data.
       *
       * Creates a vector with given size and given data.
       */
      explicit DenseVectorBlocked(Index size_in, DT_ * data) :
        Container<Mem_, DT_, IT_>(size_in)
      {
        if (size_in == Index(0))
          return;

        this->_elements.push_back(data);
        this->_elements_size.push_back(size<Perspective::pod>());

        for (Index i(0); i < this->_elements.size(); ++i)
          MemoryPool<Mem_>::increase_memory(this->_elements.at(i));
        for (Index i(0); i < this->_indices.size(); ++i)
          MemoryPool<Mem_>::increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Constructor
       *
       * \param[in] other The source DenseVector.
       *
       * Creates a vector from a DenseVector source
       */
      explicit DenseVectorBlocked(const DenseVector<Mem_, DT_, IT_> & other) :
        Container<Mem_, DT_, IT_>(other.size() / Index(BlockSize_))
      {
        convert(other);
      }

      /**
       * \brief Constructor
       *
       * \param[in] dv_in The source DenseVectorBlocked
       * \param[in] size_in The (native) size of the created vector range.
       * \param[in] offset_in The (native) starting element of the created vector range in relation to the source vector.
       *
       * Creates a vector range from a given DenseVectorBlocked
       *
       * \note The created DenseVectorBlocked has no own memory management nor own allocated memory and should be used carefully!
       */
      explicit DenseVectorBlocked(const DenseVectorBlocked & dv_in, Index size_in, Index offset_in) :
        Container<Mem_, DT_, IT_>(size_in)
      {
        this->_foreign_memory = true;

        DT_ * te(const_cast<DT_ *>(dv_in.template elements<Perspective::pod>()));
        this->_elements.push_back(te + offset_in * Index(BlockSize_));
        this->_elements_size.push_back(size<Perspective::pod>());
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file.
       *
       * Creates a vector from the given source file.
       */
      explicit DenseVectorBlocked(FileMode mode, String filename) :
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
      explicit DenseVectorBlocked(FileMode mode, std::istream & file) :
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
      explicit DenseVectorBlocked(std::vector<char> input) :
        Container<Mem_, DT_, IT_>(0)
      {
        deserialise<DT2_, IT2_>(input);
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
      explicit DenseVectorBlocked(Random & rng, Index size_in, DataType min, DataType max) :
        Container<Mem_, DT_, IT_>(size_in)
      {
        if (size_in == Index(0))
          return;

        this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(size<Perspective::pod>()));
        this->_elements_size.push_back(size<Perspective::pod>());

        this->format(rng, min, max);
      }

      /**
       * \brief Move Constructor
       *
       * \param[in] other The source vector.
       *
       * Moves another vector to this vector.
       */
      DenseVectorBlocked(DenseVectorBlocked && other) :
        Container<Mem_, DT_, IT_>(std::forward<DenseVectorBlocked>(other))
      {
      }

      /**
       * \brief Destructor
       */
      virtual ~DenseVectorBlocked()
      {
      }

      /**
       * \brief Assignment move operator
       *
       * \param[in] other The source vector.
       *
       * Moves another vector to the target vector.
       *
       * \returns The vector moved to.
       */
      DenseVectorBlocked & operator=(DenseVectorBlocked && other)
      {
        this->move(std::forward<DenseVectorBlocked>(other));

        return *this;
      }

      /** \brief Clone operation
       *
       * Create a clone of this container.
       *
       * \param[in] clone_mode The actual cloning procedure.
       *
       */
      DenseVectorBlocked clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        DenseVectorBlocked t;
        t.clone(*this, clone_mode);
        return t;
      }

      /** \brief Clone operation
       *
       * Create a clone of this container.
       *
       * \param[in] other The source container to create the clone from.
       * \param[in] clone_mode The actual cloning procedure.
       *
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void clone(const DenseVectorBlocked<Mem2_, DT2_, IT2_, BlockSize_> & other, CloneMode clone_mode = CloneMode::Deep)
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
      void convert(const DenseVectorBlocked<Mem2_, DT2_, IT2_, BlockSize_> & other)
      {
        this->assign(other);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source scalar vector.
       *
       * Use source scalar vector content as content of current vector
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const DenseVector<Mem2_, DT2_, IT2_> & other)
      {
        XASSERTM(other.size() % Index(BlockSize_) == 0, "DenseVector cannot be converted to given blocksize!");

        this->clear();

        this->_scalar_index.push_back(other.size() / Index(BlockSize_));

        this->_elements.push_back(other.get_elements().at(0));
        this->_elements_size.push_back(size<Perspective::pod>());

        for (Index i(0); i < this->_elements.size(); ++i)
          MemoryPool<Mem_>::increase_memory(this->_elements.at(i));
        for (Index i(0); i < this->_indices.size(); ++i)
          MemoryPool<Mem_>::increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Deserialisation of complete container entity.
       *
       * \param[in] input A std::vector, containing the byte array.
       *
       * Recreate a complete container entity by a single binary array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void deserialise(std::vector<char> input)
      {
        this->template _deserialise<DT2_, IT2_>(FileMode::fm_dvb, input);
      }

      /**
       * \brief Serialisation of complete container entity.
       *
       * \returns A std::vector, containing the byte array.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * See \ref FEAT::LAFEM::Container::_serialise for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> serialise()
      {
        return this->template _serialise<DT2_, IT2_>(FileMode::fm_dvb);
      }

      /**
       * \brief Retrieve a pointer to the data array.
       *
       * \tparam perspective_ template parameter to choose the return value type
       *
       * \returns Non zero element array if perspective_ = Perspective::native, e.g. treat every block as one block.
       * \returns Raw non zero element array if perspective_ = Perspective::pod, e.g. treat every entry of a block separated.
       */
      template <Perspective perspective_ = Perspective::native>
      auto elements() const -> const typename Intern::DenseVectorBlockedPerspectiveHelper<DT_, BlockSize_, perspective_>::Type *
      {
        if (this->size() == 0)
          return nullptr;

        return (const typename Intern::DenseVectorBlockedPerspectiveHelper<DT_, BlockSize_, perspective_>::Type *)(this->_elements.at(0));
      }

      /// \copydoc val()
      /// const version.
      template <Perspective perspective_ = Perspective::native>
      auto elements() -> typename Intern::DenseVectorBlockedPerspectiveHelper<DT_, BlockSize_, perspective_>::Type *
      {
        if (this->size() == 0)
          return nullptr;

        return (typename Intern::DenseVectorBlockedPerspectiveHelper<DT_, BlockSize_, perspective_>::Type *)(this->_elements.at(0));
      }

      /**
       * \brief The number of elements
       *
       * \returns number of elements of type Tiny::Vector<DT_, Blocksize_> if perspective_ = false, e.g. count every block as one entry.
       * \returns Raw number of elements of type DT_ if perspective_ = true, e.g. size * BlockSize_
       */
      template <Perspective perspective_ = Perspective::native>
      Index size() const
      {
        if (perspective_ == Perspective::pod)
          return static_cast<const Container<Mem_, DT_, IT_> *>(this)->size() * Index(BlockSize_);
        else
          return static_cast<const Container<Mem_, DT_, IT_> *>(this)->size();
      }

      /**
       * \brief Retrieve specific vector element.
       *
       * \param[in] index The index of the vector element.
       *
       * \returns Specific Tiny::Vector element.
       */
      const Tiny::Vector<DT_, BlockSize_> operator()(Index index) const
      {
        ASSERT(index < this->size());
        Tiny::Vector<DT_, BlockSize_> t;
        MemoryPool<Mem_>::download(t.v, this->_elements.at(0) + index * Index(BlockSize_), Index(BlockSize_));
        return t;
      }

      /**
       * \brief Set specific vector element.
       *
       * \param[in] index The index of the vector element.
       * \param[in] value The value to be set.
       */
      void operator()(Index index, const Tiny::Vector<DT_, BlockSize_> & value)
      {
        ASSERT(index < this->size());
        MemoryPool<Mem_>::upload(this->_elements.at(0) + index * Index(BlockSize_), value.v, Index(BlockSize_));
      }

      /**
       * \brief Returns a descriptive string.
       *
       * \returns A string describing the container.
       */
      static String name()
      {
        return "DenseVectorBlocked";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The vector to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      void copy(const DenseVectorBlocked & x, bool full = false)
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
      void copy(const DenseVectorBlocked<Mem2_, DT_, IT_, BlockSize_> & x, bool full = false)
      {
        this->_copy_content(x, full);
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
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
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
            throw InternalError(__func__, __FILE__, __LINE__, "Input-file is not a compatible mtx-vector-file");
          }
          while (! file.eof())
          {
            std::getline(file, line);
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
              throw InternalError(__func__, __FILE__, __LINE__, "Input-file is no dense-vector-file");
          }

          DenseVectorBlocked<Mem::Main, DT_, IT_, BlockSize_> tmp(rows / BlockSize_);
          DT_ * pval(tmp.template elements<Perspective::pod>());

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

          _size() = Index(data.size()) / Index(BlockSize_);
          this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(Index(data.size())));
          this->_elements_size.push_back(Index(data.size()));
          MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(0), &data[0], Index(data.size()));
          break;
        }
        case FileMode::fm_dvb:
        case FileMode::fm_binary:
        {
          this->clear();

          this->template _deserialise<double, uint64_t>(FileMode::fm_dvb, file);
          break;
        }
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
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
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
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
          DenseVectorBlocked<Mem::Main, DT_, IT_, BlockSize_> temp;
          temp.convert(*this);

          const Index tsize(temp.template size<Perspective::pod>());
          file << "%%MatrixMarket matrix array real general" << std::endl;
          file << tsize << " " << 1 << std::endl;

          const DT_ * pval(temp.template elements<Perspective::pod>());
          for (Index i(0); i < tsize; ++i, ++pval)
          {
            file << std::scientific << *pval << std::endl;
          }
          break;
        }
        case FileMode::fm_exp:
        {
          DT_ * temp = MemoryPool<Mem::Main>::template allocate_memory<DT_>((this->size<Perspective::pod>()));
          MemoryPool<Mem_>::template download<DT_>(temp, this->template elements<Perspective::pod>(), this->size<Perspective::pod>());

          for (Index i(0); i < this->size<Perspective::pod>(); ++i)
          {
            file << std::scientific << temp[i] << std::endl;
          }

          MemoryPool<Mem::Main>::release_memory(temp);
          break;
        }
        case FileMode::fm_dvb:
        case FileMode::fm_binary:
          this->template _serialise<double, uint64_t>(FileMode::fm_dvb, file);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
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
        const DenseVectorBlocked & x,
        const DenseVectorBlocked & y,
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

        Statistics::add_flops(this->size<Perspective::pod>() * 2);
        Arch::Axpy<Mem_>::dv(elements<Perspective::pod>(), alpha, x.template elements<Perspective::pod>(), y.template elements<Perspective::pod>(), this->size<Perspective::pod>());

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this_i \leftarrow x_i \cdot y_i\f$
       *
       * \param[in] x The first factor.
       * \param[in] y The second factor.
       */
      void component_product(const DenseVectorBlocked & x, const DenseVectorBlocked & y)
      {
        XASSERTM(this->size() == x.size(), "Vector size does not match!");
        XASSERTM(this->size() == y.size(), "Vector size does not match!");

        TimeStamp ts_start;

        Arch::ComponentProduct<Mem_>::value(elements<Perspective::pod>(), x.template elements<Perspective::pod>(), y.template elements<Perspective::pod>(), this->size<Perspective::pod>());
        Statistics::add_flops(this->size<Perspective::pod>());

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
      void component_invert(const DenseVectorBlocked & x, const DT_ alpha = DT_(1))
      {
        XASSERTM(this->size() == x.size(), "Vector size does not match!");

        TimeStamp ts_start;

        Arch::ComponentInvert<Mem_>::value(this->template elements<Perspective::pod>(), x.template elements<Perspective::pod>(), alpha, this->size<Perspective::pod>());
        Statistics::add_flops(this->size<Perspective::pod>());

        TimeStamp ts_stop;
        Statistics::add_time_axpy(ts_stop.elapsed(ts_start));
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha~ x \f$
       *
       * \param[in] x The vector to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const DenseVectorBlocked & x, const DT_ alpha)
      {
        XASSERTM(x.size() == this->size(), "Vector size does not match!");

        TimeStamp ts_start;

        Arch::Scale<Mem_>::value(elements<Perspective::pod>(), x.template elements<Perspective::pod>(), alpha, this->size<Perspective::pod>());
        Statistics::add_flops(this->size<Perspective::pod>());

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
      DataType triple_dot(const DenseVectorBlocked & x, const DenseVectorBlocked & y) const
      {
        XASSERTM(x.template size<Perspective::pod>() != this->template size<Perspective::pod>() || y.template size<Perspective::pod>() == this->template size<Perspective::pod>(), "Vector sizes does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->template size<Perspective::pod>() * 3);
        DataType result = Arch::TripleDotProduct<Mem_>::value(this->template elements<Perspective::pod>(), x.template elements<Perspective::pod>(), y.template elements<Perspective::pod>(), this->template size<Perspective::pod>());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Calculate \f$this \leftarrow this \cdot x\f$
       *
       * \param[in] x The other vector.
       *
       * \return The calculated dot product.
       */
      DataType dot(const DenseVectorBlocked & x) const
      {
        XASSERTM(x.size() == this->size(), "Vector size does not match!");

        TimeStamp ts_start;

        Statistics::add_flops(this->size<Perspective::pod>() * 2);
        DataType result = Arch::DotProduct<Mem_>::value(elements<Perspective::pod>(), x.template elements<Perspective::pod>(), this->size<Perspective::pod>());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Calculates and returns the euclid norm of this vector.
       *
       * \return The calculated norm.
       */
      DT_ norm2() const
      {
        TimeStamp ts_start;

        Statistics::add_flops(this->size<Perspective::pod>() * 2);
        DataType result = Arch::Norm2<Mem_>::value(elements<Perspective::pod>(), this->size<Perspective::pod>());

        TimeStamp ts_stop;
        Statistics::add_time_reduction(ts_stop.elapsed(ts_start));

        return result;
      }

      /**
       * \brief Calculates and returns the squared euclid norm of this vector.
       *
       * \return The calculated norm.
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
      DT_ max_element() const
      {
        TimeStamp ts_start;

        Index max_index = Arch::MaxElement<Mem_>::value(this->elements<Perspective::pod>(), this->size<Perspective::pod>());
        ASSERT(max_index < this->size<Perspective::pod>());
        DT_ result;
        MemoryPool<Mem_>::template download<DT_>(&result, this->elements<Perspective::pod>() + max_index, 1);
        result = Math::abs(result);

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

        DenseVectorBlocked<Mem::Main, DT_, IT_, BlockSize_> local;
        local.convert(*this);
        perm(local.template elements<Perspective::native>());
        this->assign(local);
      }

      /// \cond internal
      /// Writes the vector-entries in an allocated array
      void set_vec(DT_ * const pval_set) const
      {
        MemoryPool<Mem_>::copy(pval_set, this->template elements<Perspective::pod>(), this->size<Perspective::pod>());
      }

      /// Writes data of an array in the vector
      void set_vec_inv(const DT_ * const pval_set)
      {
        MemoryPool<Mem_>::copy(this->template elements<Perspective::pod>(), pval_set, this->size<Perspective::pod>());
      }
      /// \endcond

      /**
       * \brief DenseVectorBlocked comparison operator
       *
       * \param[in] a A vector to compare with.
       * \param[in] b A vector to compare with.
       */
      template <typename Mem2_>
      friend bool operator==(const DenseVectorBlocked & a, const DenseVectorBlocked<Mem2_, DT_, IT_, BlockSize_> & b)
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
          ta = const_cast<DT_ *>(a.template elements<Perspective::pod>());
        }
        else
        {
          ta = new DT_[a.template size<Perspective::pod>()];
          MemoryPool<Mem_>::template download<DT_>(ta, a.template elements<Perspective::pod>(), a.template size<Perspective::pod>());
        }

        if (std::is_same<Mem::Main, Mem2_>::value)
        {
          tb = const_cast<DT_ *>(b.template elements<Perspective::pod>());
        }
        else
        {
          tb = new DT_[b.template size<Perspective::pod>()];
          MemoryPool<Mem2_>::template download<DT_>(tb, b.template elements<Perspective::pod>(), b.template size<Perspective::pod>());
        }

        for (Index i(0); i < a.template size<Perspective::pod>(); ++i)
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
       * \brief DenseVectorBlocked streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The vector to be streamed.
       */
      friend std::ostream & operator<<(std::ostream & lhs, const DenseVectorBlocked & b)
      {
        lhs << "[";
        for (Index i(0); i < b.size(); ++i)
        {
          Tiny::Vector<DT_, BlockSize_> t = b(i);
          for (int j(0); j < BlockSize_; ++j)
            lhs << "  " << t[j];
        }
        lhs << "]";

        return lhs;
      }
    }; // class DenseVectorBlocked<...>

  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_DENSE_VECTOR_BLOCKED_HPP
