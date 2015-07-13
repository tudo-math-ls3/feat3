#pragma once
#ifndef KERNEL_LAFEM_DENSE_VECTOR_HPP
#define KERNEL_LAFEM_DENSE_VECTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/lafem/forward.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/vector_base.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/edi.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/lafem/arch/component_invert.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/dot_product.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/component_product.hpp>
#include <kernel/lafem/arch/synch.hpp>
#include <kernel/adjacency/permutation.hpp>

#include <iostream>
#include <fstream>
#include <string>
#include <stdint.h>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Dense data vector class template.
     *
     * \tparam Mem_ The \ref FEAST::Mem "memory architecture" to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indextype to be used.
     *
     * This class represents a vector of continuous data in memory. \n \n
     * Data survey: \n
     * _elements[0]: raw number values \n
     * _scalar_index[0]: container size
     * _scalar_index[1]: boolean flag, signaling DenseVectorRange usage
     *
     * Refer to \ref lafem_design for general usage informations.
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class DenseVector : public Container<Mem_, DT_, IT_>, public VectorBase
    {
    private:
      Index & _size()
      {
        return this->_scalar_index.at(0);
      }

      /// \cond internal
      template <typename VT_, typename IT2_>
      void _convert(const typename VT_::template ContainerType<Mem_, DT_, IT2_> & a)
      {
        DenseVector<Mem_, DT_, IT_> vec(a.size());
        a.set_vec(vec.elements());

        this->assign(vec);
      }

      template <typename VT_>
      void _convert(const VT_ & a)
      {
        typename VT_::template ContainerType<Mem_, DT_, IT_> ta;
        ta.convert(a);

        this->convert(ta);
      }

      template <typename VT_, typename Mem2_, typename IT2_>
      void _copy(const typename VT_::template ContainerType<Mem2_, DT_, IT2_> & a)
      {
        if (std::is_same<Mem_, Mem2_>::value)
        {
          a.set_vec(this->elements());
        }
        else
        {
          typename VT_::template ContainerType<Mem_, DT_, IT2_> ta;
          ta.convert(a);

          this->copy(ta);
        }
      }

      template <typename VT_, typename Mem2_, typename IT2_>
      void _copy_inv(typename VT_::template ContainerType<Mem2_, DT_, IT2_> & a) const
      {
        if (std::is_same<Mem_, Mem2_>::value)
        {
          a.set_vec_inv(this->elements());
        }
        else
        {
          DenseVector<Mem2_, DT_, IT_> t_this;
          t_this.convert(*this);

          t_this.copy_inv(a);
        }
      }
      /// \endcond

    public:
      /// Our datatype
      typedef DT_ DataType;
      /// Our indextype
      typedef IT_ IndexType;
      /// Our memory architecture type
      typedef Mem_ MemType;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using ContainerType = class DenseVector<Mem2_, DT2_, IT2_>;

      /**
       * \brief Constructor
       *
       * Creates an empty non dimensional vector.
       */
      explicit DenseVector() :
        Container<Mem_, DT_, IT_> (0)
      {
        CONTEXT("When creating DenseVector");
        this->_scalar_index.push_back(0);
        this->_elements.push_back(nullptr);
      }

      /**
       * \brief Constructor
       *
       * \param[in] size The size of the created vector.
       *
       * Creates a vector with a given size.
       */
      explicit DenseVector(Index size_in) :
        Container<Mem_, DT_, IT_>(size_in)
      {
        CONTEXT("When creating DenseVector");

        this->_scalar_index.push_back(0);
        this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(size_in));
        this->_elements_size.push_back(size_in);
      }

      /**
       * \brief Constructor
       *
       * \param[in] size The size of the created vector.
       * \param[in] value The value, each element will be set to.
       *
       * Creates a vector with given size and value.
       */
      explicit DenseVector(Index size_in, DT_ value) :
        Container<Mem_, DT_, IT_>(size_in)
      {
        CONTEXT("When creating DenseVector");

        this->_scalar_index.push_back(0);
        this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(size_in));
        this->_elements_size.push_back(size_in);

        Util::MemoryPool<Mem_>::instance()->set_memory(this->_elements.at(0), value, size_in);
      }

      /**
       * \brief Constructor
       *
       * \param[in] size The size of the created vector.
       * \param[in] data An array containing the value data.
       *
       * Creates a vector with given size and given data.
       *
       * \note The array must be allocated by FEAST's own memory pool
       */
      explicit DenseVector(Index size_in, DT_ * data) :
        Container<Mem_, DT_, IT_>(size_in)
      {
        CONTEXT("When creating DenseVector");

        this->_scalar_index.push_back(0);
        this->_elements.push_back(data);
        this->_elements_size.push_back(size_in);

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
      }

      /**
       * \brief Constructor
       *
       * \param[in] dv_ind The source DenseVector
       * \param[in] size_in The size of the created vector range.
       * \param[in] offset The starting element of the created vector range in relation to the source vector.
       *
       * Creates a vector range from a given DenseVector
       *
       * \note The created DenseVector has no own memory management nor own allocated memory and should be used carefully!
       */
      explicit DenseVector(const DenseVector & dv_in, Index size_in, Index offset_in) :
        Container<Mem_, DT_, IT_>(size_in)
      {
        CONTEXT("When creating DenseVector");

        this->_scalar_index.push_back(1);
        DT_ * te(const_cast<DT_*>(dv_in.elements()));
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
        Container<Mem_, DT_, IT_>(other.raw_size())
      {
        CONTEXT("When creating DenseVector");

        this->_scalar_index.push_back(0);
        convert(other);
      }

      /**
       * \brief Constructor
       *
       * \param[in] mode The used file format.
       * \param[in] filename The source file in EXP format.
       *
       * Creates a vector from the given source file.
       */
      explicit DenseVector(FileMode mode, String filename) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating DenseVector");

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
      explicit DenseVector(FileMode mode, std::istream& file) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating DenseVector");

        read_from(mode, file);
      }

      /**
       * \brief Constructor
       *
       * \param[in] std::vector<char> A std::vector, containing the byte.
       *
       * Creates a vector from the given byte array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      explicit DenseVector(std::vector<char> input) :
        Container<Mem_, DT_, IT_>(0)
      {
        CONTEXT("When creating DenseVector");
        deserialise<DT2_, IT2_>(input);
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
        CONTEXT("When moving DenseVector");
      }

      /**
       * \brief Destructor
       *
       * Destroys the DenseVector and releases all of its used arrays if its not marked a range vector.
       */
      virtual ~DenseVector()
      {
        CONTEXT("When destroying DenseVector");

        // avoid releasing memory by base class destructore, because we do not own the referenced memory
        if (this->_scalar_index.size() > 0 && this->_scalar_index.at(1) == 1)
        {
          for (Index i(0) ; i < this->_elements.size() ; ++i)
            this->_elements.at(i) = nullptr;
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            this->_indices.at(i) = nullptr;
        }
      }

      /**
       * \brief Assignment move operator
       *
       * \param[in] other The source vector.
       *
       * Moves another vector to the target vector.
       */
      DenseVector & operator= (DenseVector && other)
      {
        CONTEXT("When moving DenseVector");

        this->move(std::forward<DenseVector>(other));

        return *this;
      }

      InsertWeakClone( DenseVector );

      /** \brief Shallow copy operation
       *
       * Create a shallow copy of itself.
       *
       */
      DenseVector shared() const
      {
        CONTEXT("When sharing DenseVector");
        DenseVector r;
        r.assign(*this);
        return r;
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
        CONTEXT("When converting DenseVector");
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
        CONTEXT("When converting DenseVector");

        this->clear();

        this->_scalar_index.push_back(other.raw_size());
        this->_scalar_index.push_back(0);
        this->_elements.push_back(other.get_elements().at(0));
        this->_elements_size.push_back(this->size());

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          Util::MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
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
        CONTEXT("When converting DenseVector");

        this->template _convert<VT_>(a);
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
        this->template _deserialise<DT2_, IT2_>(FileMode::fm_dv, input);
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
        return this->template _serialise<DT2_, IT2_>(FileMode::fm_dv);
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The vector to be copied (could be of any format; must have same size).
       */
      template<typename VT_>
      void copy(const VT_ & a)
      {
        if (this->size() != a.size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vectors have not the same size!");

        this->template _copy<VT_>(a);
      }

      /**
       * \brief Performs \f$x \leftarrow this\f$.
       *
       * \param[in] x The target-vector to be copied to (could be of any format; must have same size).
       */
      template<typename VT_>
      void copy_inv(VT_ & a) const
      {
        if (this->size() != a.size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vectors have not the same size!");

        this->template _copy_inv<VT_>(a);
      }


      /**
       * \brief Read in vector from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, String filename)
      {
        CONTEXT("When reading in DenseVector");

        switch(mode)
        {
        case FileMode::fm_mtx:
          read_from_mtx(filename);
          break;
        case FileMode::fm_exp:
          read_from_exp(filename);
          break;
        case FileMode::fm_dv:
          read_from_dv(filename);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Read in vector from stream.
       *
       * \param[in] mode The used file format.
       * \param[in] file The stream that shall be read in.
       */
      void read_from(FileMode mode, std::istream& file)
      {
        CONTEXT("When reading in DenseVector");

        switch(mode)
        {
        case FileMode::fm_mtx:
          read_from_mtx(file);
          break;
        case FileMode::fm_exp:
          read_from_exp(file);
          break;
        case FileMode::fm_dv:
          read_from_dv(file);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Read in vector from MatrixMarket mtx file.
       *
       * \param[in] filename The file that shall be read in.
       */
      void read_from_mtx(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
        read_from_mtx(file);
        file.close();
      }

      /**
       * \brief Read in vector from MatrixMarket mtx stream.
       *
       * \param[in] file The stream that shall be read in.
       */
      void read_from_mtx(std::istream& file)
      {
        this->clear();
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);

        Index rows;
        String line;
        std::getline(file, line);
        if (line.find("%%MatrixMarket matrix array real general") == String::npos)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Input-file is not a compatible mtx-vector-file");
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
            throw InternalError(__func__, __FILE__, __LINE__, "Input-file is no dense-vector-file");
        }

        DenseVector<Mem::Main, DT_, IT_> tmp(rows);
        DT_ * pval(tmp.elements());

        while(!file.eof())
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
      }

      /**
       * \brief Read in vector from ASCII file.
       *
       * \param[in] filename The file that shall be read in.
       */
      void read_from_exp(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
        read_from_exp(file);
        file.close();
      }

      /**
       * \brief Read in vector from ASCII stream.
       *
       * \param[in] file The stream that shall be read in.
       */
      void read_from_exp(std::istream& file)
      {
        this->clear();
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);

        std::vector<DT_> data;

        while(!file.eof())
        {
          std::string line;
          std::getline(file, line);
          if(line.find("#", 0) < line.npos)
            continue;
          if(file.eof())
            break;

          std::string n_z_s;

          std::string::size_type first_digit(line.find_first_not_of(" "));
          line.erase(0, first_digit);
          std::string::size_type eol(line.length());
          for(unsigned long i(0) ; i < eol ; ++i)
          {
            n_z_s.append(1, line[i]);
          }

          DT_ n_z((DT_)atof(n_z_s.c_str()));

          data.push_back(n_z);

        }

        _size() = Index(data.size());
        this->_elements.push_back(Util::MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(Index(data.size())));
        this->_elements_size.push_back(Index(data.size()));
        Util::MemoryPool<Mem_>::instance()->template upload<DT_>(this->_elements.at(0), &data[0], Index(data.size()));
      }

      /**
       * \brief Read in vector from binary file.
       *
       * \param[in] filename The file that shall be read in.
       */
      void read_from_dv(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
        read_from_dv(file);
        file.close();
      }

      /**
       * \brief Read in vector from binary stream.
       *
       * \param[in] file The stream that shall be read in.
       */
      void read_from_dv(std::istream& file)
      {
        this->clear();
        this->_scalar_index.push_back(0);
        this->_scalar_index.push_back(0);

        this->template _deserialise<double, uint64_t>(FileMode::fm_dv, file);
      }


      /**
       * \brief Write out vector to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the vector shall be stored.
       */
      void write_out(FileMode mode, String filename) const
      {
        CONTEXT("When writing out DenseVector");

        switch(mode)
        {
        case FileMode::fm_mtx:
          write_out_mtx(filename);
          break;
        case FileMode::fm_exp:
          write_out_exp(filename);
          break;
        case FileMode::fm_dv:
          write_out_dv(filename);
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
        CONTEXT("When writing out DenseVector");

        switch(mode)
        {
        case FileMode::fm_mtx:
          write_out_mtx(file);
          break;
        case FileMode::fm_exp:
          write_out_exp(file);
          break;
        case FileMode::fm_dv:
          write_out_dv(file);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Write out vector to MatrixMarket mtx file.
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
       * \brief Write out vector to MatrixMarket mtx file.
       *
       * \param[in] file The stream that shall be written to.
       */
      void write_out_mtx(std::ostream& file) const
      {
        DenseVector<Mem::Main, DT_, IT_> temp;
        temp.convert(*this);

        const Index tsize(temp.size());
        file << "%%MatrixMarket matrix array real general" << std::endl;
        file << tsize << " " << 1 << std::endl;

        const DT_ * pval(temp.elements());
        for (Index i(0) ; i < tsize ; ++i, ++pval)
        {
          file << std::scientific << *pval << std::endl;
        }
      }

      /**
       * \brief Write out vector to file.
       *
       * \param[in] filename The file where the vector shall be stored.
       */
      void write_out_exp(String filename) const
      {
        std::ofstream file(filename.c_str(), std::ofstream::out);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
        write_out_exp(file);
        file.close();
      }

      /**
       * \brief Write out vector to file.
       *
       * \param[in] file The stream that shall be written to.
       */
      void write_out_exp(std::ostream& file) const
      {
        DT_ * temp = Util::MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>((this->size()));
        Util::MemoryPool<Mem_>::template download<DT_>(temp, this->_elements.at(0), this->size());

        for (Index i(0) ; i < this->size() ; ++i)
        {
          file << std::scientific << temp[i] << std::endl;
        }

        Util::MemoryPool<Mem::Main>::instance()->release_memory(temp);
      }

      /**
       * \brief Write out vector to file.
       *
       * \param[in] filename The file where the vector shall be stored.
       */
      void write_out_dv(String filename) const
      {
        std::ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
        write_out_dv(file);
        file.close();
      }

      /**
       * \brief Write out vector to file.
       *
       * \param[in] file The stream that shall be written to.
       */
      void write_out_dv(std::ostream& file) const
      {
        if (! std::is_same<DT_, double>::value)
          std::cout<<"Warning: You are writing out a dense vector that is not double precision!"<<std::endl;

        this->template _serialise<double, uint64_t>(FileMode::fm_dv, file);
      }

      /**
       * \brief Get a pointer to the data array.
       *
       * \returns Pointer to the data array.
       */
      DT_ * elements()
      {
        return this->_elements.at(0);
      }

      DT_ const * elements() const
      {
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
        CONTEXT("When retrieving DenseVector element");

        ASSERT(index < this->size(), "Error: " + stringify(index) + " exceeds dense vector size " + stringify(this->size()) + " !");
        return Util::MemoryPool<Mem_>::get_element(this->_elements.at(0), index);
      }

      /**
       * \brief Set specific vector element.
       *
       * \param[in] index The index of the vector element.
       * \param[in] value The value to be set.
       */
      void operator()(Index index, DT_ value)
      {
        CONTEXT("When setting DenseVector element");

        ASSERT(index < this->size(), "Error: " + stringify(index) + " exceeds dense vector size " + stringify(this->size()) + " !");
        Util::MemoryPool<Mem_>::set_memory(this->_elements.at(0) + index, value);
      }

      /**
       * \brief Create temporary object for direct data manipulation.
       * \warning Be aware, that any synchronisation only takes place, when the object is destroyed!
       *
       * \param[in] index The index of the vector element.
       */
      EDI<Mem_, DT_> edi(Index index)
      {
        EDI<Mem_, DT_> t(Util::MemoryPool<Mem_>::get_element(this->_elements.at(0), index), this->_elements.at(0) + index);
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
       */
      void copy(const DenseVector & x)
      {
        this->_copy_content(x);
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The vector to be copied.
       */
      template <typename Mem2_>
      void copy(const DenseVector<Mem2_, DT_, IT_> & x)
      {
        this->_copy_content(x);
      }

      ///@name Linear algebra operations
      ///@{
      /**
       * \brief Calculate \f$this \leftarrow \alpha x + y\f$
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
        if (x.size() != y.size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");
        if (x.size() != this->size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

        // check for special cases
        // r <- x + y
        if(Math::abs(alpha - DT_(1)) < Math::eps<DT_>())
          Arch::Sum<Mem_>::value(this->elements(), x.elements(), y.elements(), this->size());
        // r <- y - x
        else if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
          Arch::Difference<Mem_>::value(this->elements(), y.elements(), x.elements(), this->size());
        // r <- y
        else if(Math::abs(alpha) < Math::eps<DT_>())
          this->copy(y);
        // r <- y + alpha*x
        else
          Arch::Axpy<Mem_>::dv(this->elements(), alpha, x.elements(), y.elements(), this->size());
      }

      /**
       * \brief Calculate \f$this_i \leftarrow x_i \cdot y_i\f$
       *
       * \param[in] x The first factor.
       * \param[in] y The second factor.
       */
      void component_product(const DenseVector & x, const DenseVector & y)
      {
        if (this->size() != x.size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");
        if (this->size() != y.size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

        Arch::ComponentProduct<Mem_>::value(this->elements(), x.elements(), y.elements(), this->size());
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
        if (this->size() != x.size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

        Arch::ComponentInvert<Mem_>::value(this->elements(), x.elements(), alpha, this->size());
      }

      /**
       * \brief Calculate \f$this \leftarrow \alpha x \f$
       *
       * \param[in] x The vector to be scaled.
       * \param[in] alpha A scalar to scale x with.
       */
      void scale(const DenseVector & x, const DT_ alpha)
      {
        if (x.size() != this->size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

        Arch::Scale<Mem_>::value(this->elements(), x.elements(), alpha, this->size());
      }

      /**
       * \brief Calculate \f$result \leftarrow x^T \mathrm{diag}(this) y \f$
       *
       * \param[in] x The first vector.
       *
       * \param[in] x The second vector.
       *
       * \return The computed triple dot product.
       */
      DataType triple_dot(const DenseVector & x, const DenseVector & y) const
      {
        if (x.size() != this->size() || y.size() != this->size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector sizes does not match!");

        return Arch::TripleDotProduct<Mem_>::value(this->elements(), x.elements(), y.elements(), this->size());
      }

      /**
       * \brief Calculate \f$result \leftarrow this \cdot this\f$
       *
       * \param[in] x The other vector.
       *
       * \return The computed dot product.
       */
      DataType dot(const DenseVector & x) const
      {
        if (x.size() != this->size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

        return Arch::DotProduct<Mem_>::value(this->elements(), x.elements(), this->size());
      }

      /**
       * \brief Calculate \f$r \leftarrow this \cdot x\f$
       *
       * \param[out] r The dot product result.
       * \param[in] x The other vector.
       * \param[in] gate The gateway to be used for the global operation.
       *
       */
      DataType dot(const DenseVector & x, const Arch::DotGatewayBase<Mem_, DenseVector>* gate) const
      {
        if (x.size() != this->size())
          throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

        return gate->value(*this, x);
      }

      /**
       * \brief Calculates and returns the euclid norm of this vector, global version.
       *
       * \return The computed norm.
       *
       */
      DT_ norm2(const Arch::Norm2GatewayBase<Mem_, DenseVector>* gate) const
      {
        return gate->value(*this);
      }

      /**
       * \brief Calculates and returns the euclid norm of this vector.
       *
       */
      DT_ norm2() const
      {
        return Arch::Norm2<Mem_>::value(this->elements(), this->size());
      }

      /**
       * \brief Calculates and returns the euclid norm of this vector, global version.
       *
       */
      DT_ norm2sqr(const Arch::Norm2SquaredGatewayBase<Mem_, DenseVector>* gate) const
      {
        return gate->value(*this);
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
       * \brief Synchronise type-0 vector
       *
       */
      DenseVector& synch0(const Arch::SynchVec0GatewayBase<Mem_, DenseVector>* gate)
      {
        return gate->value(*this);
      }

      /**
       * \brief Synchronise type-1 vector
       *
       */
      DenseVector& synch1(const Arch::SynchVec1GatewayBase<Mem_, DenseVector>* gate)
      {
        return gate->value(*this);
      }
      ///@}

      /// Permutate vector according to the given Permutation
      void permute(Adjacency::Permutation & perm)
      {
        ASSERT(perm.size() == this->size(), "Error: Container size " + stringify(this->size()) + " does not match permutation size " + stringify(perm.size()) + " !");

        DenseVector<Mem::Main, DT_, IT_> local;
        local.convert(*this);
        perm(local.elements());
        this->assign(local);
      }

      /// Writes the vector-entries in an allocated array
      void set_vec(DT_ * const pval_set) const
      {
        Util::MemoryPool<Mem_>::copy(pval_set, this->elements(), this->size());
      }

      /// Writes data of an array in the vector
      void set_vec_inv(const DT_ * const pval_set)
      {
        Util::MemoryPool<Mem_>::copy(this->elements(), pval_set, this->size());
      }
      /// \endcond

      /**
       * \brief DenseVector comparison operator
       *
       * \param[in] a A vector to compare with.
       * \param[in] b A vector to compare with.
       */
      template <typename Mem2_> friend bool operator== (const DenseVector & a, const DenseVector<Mem2_, DT_, IT_> & b)
      {
        CONTEXT("When comparing DenseVectors");

        if (a.size() != b.size())
          return false;
        if (a.get_elements().size() != b.get_elements().size())
          return false;
        if (a.get_indices().size() != b.get_indices().size())
          return false;

        if(a.size() == 0 && b.size() == 0)
          return true;

        bool ret(true);

        DT_ * ta;
        DT_ * tb;

        if(std::is_same<Mem::Main, Mem_>::value)
          ta = (DT_*)a.elements();
        else
        {
          ta = new DT_[a.size()];
          Util::MemoryPool<Mem_>::instance()->template download<DT_>(ta, a.elements(), a.size());
        }
        if(std::is_same<Mem::Main, Mem2_>::value)
          tb = (DT_*)b.elements();
        else
        {
          tb = new DT_[b.size()];
          Util::MemoryPool<Mem2_>::instance()->template download<DT_>(tb, b.elements(), b.size());
        }

        for (Index i(0) ; i < a.size() ; ++i)
          if (ta[i] != tb[i])
          {
            ret = false;
            break;
          }

        if(! std::is_same<Mem::Main, Mem_>::value)
          delete[] ta;
        if(! std::is_same<Mem::Main, Mem2_>::value)
          delete[] tb;

        return ret;
      }

      /**
       * \brief DenseVector streaming operator
       *
       * \param[in] lhs The target stream.
       * \param[in] b The vector to be streamed.
       */
      friend std::ostream & operator<< (std::ostream & lhs, const DenseVector & b)
      {
        lhs << "[";
        for (Index i(0) ; i < b.size() ; ++i)
        {
          lhs << "  " << b(i);
        }
        lhs << "]";

        return lhs;
      }
    }; // class DenseVector<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_DENSE_VECTOR_HPP
