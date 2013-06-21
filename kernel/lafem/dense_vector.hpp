#pragma once
#ifndef KERNEL_LAFEM_DENSE_VECTOR_HPP
#define KERNEL_LAFEM_DENSE_VECTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>
#include <kernel/lafem/vector_base.hpp>

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
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a vector of continuous data in memory. \n \n
     * Data survey: \n
     * _elements[0]: raw number values \n
     * _scalar_index[0]: container size
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_>
    class DenseVector : public Container<Mem_, DT_>, public VectorBase
    {
      private:
        void _read_from_exp(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in);
          if (! file.is_open())
            throw InternalError("Unable to open Vector file " + filename);
          _read_from_exp(file);
          file.close();
        }

        void _read_from_exp(std::istream& file)
        {
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

            DT_ n_z = (DT_)atof(n_z_s.c_str());

            data.push_back(n_z);

          }

          this->_scalar_index.at(0) = Index(data.size());
          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(Index(data.size() * sizeof(DT_))));
          this->_elements_size.push_back(Index(data.size()));
          MemoryPool<Mem_>::instance()->upload(this->_elements.at(0), &data[0], Index(data.size() * sizeof(DT_)));
        }

        void _read_from_dv(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
          if (! file.is_open())
            throw InternalError("Unable to open Vector file " + filename);
          _read_from_dv(file);
          file.close();
        }

        void _read_from_dv(std::istream& file)
        {
          uint64_t size;
          file.read((char *)&size, sizeof(uint64_t));
          this->_scalar_index.at(0) = (Index)size;
          this->_elements_size.push_back(this->_scalar_index.at(0));

          double * ctemp = new double[size];
          file.read((char *)ctemp, size * sizeof(double));

          DT_ * temp = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(0)) * sizeof(DT_));
          for (Index i(0) ; i < size ; ++i)
          {
            temp[i] = (DT_)ctemp[i];
          }
          delete[] ctemp;
          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(this->_scalar_index.at(0) * sizeof(DT_)));
          MemoryPool<Mem_>::instance()->upload(this->_elements.at(0), temp, this->_scalar_index.at(0) * sizeof(DT_));
          MemoryPool<Mem::Main>::instance()->release_memory(temp);
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
        explicit DenseVector() :
          Container<Mem_, DT_> (0)
        {
          CONTEXT("When creating DenseVector");
        }

        /**
         * \brief Constructor
         *
         * \param[in] size The size of the created vector.
         *
         * Creates a vector with a given size.
         */
        explicit DenseVector(Index size) :
          Container<Mem_, DT_>(size)
        {
          CONTEXT("When creating DenseVector");

          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(size * sizeof(DT_)));
          this->_elements_size.push_back(size);
        }

        /**
         * \brief Constructor
         *
         * \param[in] size The size of the created vector.
         * \param[in] value The value, each element will be set to.
         *
         * Creates a vector with given size and value.
         */
        explicit DenseVector(Index size, DT_ value) :
          Container<Mem_, DT_>(size)
        {
          CONTEXT("When creating DenseVector");

          this->_elements.push_back((DT_*)MemoryPool<Mem_>::instance()->allocate_memory(size * sizeof(DT_)));
          this->_elements_size.push_back(size);

          MemoryPool<Mem_>::instance()->set_memory(this->_elements.at(0), value, size);
        }

        /**
         * \brief Constructor
         *
         * \param[in] size The size of the created vector.
         * \param[in] data An array containing the value data.
         *
         * Creates a vector with given size and given data.
         */
        explicit DenseVector(Index size, DT_ * data) :
          Container<Mem_, DT_>(size)
        {
          CONTEXT("When creating DenseVector");

          this->_elements.push_back(data);
          this->_elements_size.push_back(size);

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
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
          Container<Mem_, DT_>(0)
        {
          CONTEXT("When creating DenseVector");

          switch(mode)
          {
            case fm_exp:
              _read_from_exp(filename);
              break;
            case fm_dv:
              _read_from_dv(filename);
              break;
            default:
              throw InternalError("Filemode not supported!");
          }
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
          Container<Mem_, DT_>(0)
        {
          CONTEXT("When creating DenseVector");

          switch(mode)
          {
            case fm_exp:
              _read_from_exp(file);
              break;
            case fm_dv:
              _read_from_dv(file);
              break;
            default:
              throw InternalError("Filemode not supported!");
          }
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source vector.
         *
         * Creates a shallow copy of a given vector.
         */
        DenseVector(const DenseVector<Mem_, DT_> & other) :
          Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying DenseVector");
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source vector.
         *
         * Creates a copy of a given vector from another memory architecture.
         */
        template <typename Mem2_, typename DT2_>
        DenseVector(const DenseVector<Mem2_, DT2_> & other) :
            Container<Mem_, DT_>(other)
        {
          CONTEXT("When copying DenseVector");
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this vector.
         */
        DenseVector<Mem_, DT_> clone()
        {
          CONTEXT("When cloning DenseVector");

          DenseVector<Mem_, DT_> t;
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
        DenseVector<Mem_, DT_> & operator= (const DenseVector<Mem_, DT_> & other)
        {
          CONTEXT("When assigning DenseVector");

          assign(other);

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
        DenseVector<Mem_, DT_> & operator= (const DenseVector<Mem2_, DT2_> & other)
        {
          CONTEXT("When assigning DenseVector");

          assign(other);

          return *this;
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
            case fm_exp:
              write_out_exp(filename);
              break;
            case fm_dv:
              write_out_dv(filename);
              break;
            default:
              throw InternalError("Filemode not supported!");
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
            case fm_exp:
              write_out_exp(file);
              break;
            case fm_dv:
              write_out_dv(file);
              break;
            default:
              throw InternalError("Filemode not supported!");
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
            throw InternalError("Unable to open Vector file " + filename);
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
          DT_ * temp = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(0)) * sizeof(DT_));
          MemoryPool<Mem_>::download(temp, this->_elements.at(0), this->_scalar_index.at(0) * sizeof(DT_));

          for (Index i(0) ; i < this->_scalar_index.at(0) ; ++i)
          {
            file << std::scientific << (double)temp[i] << std::endl;
          }

          MemoryPool<Mem::Main>::instance()->release_memory(temp);
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
            throw InternalError("Unable to open Vector file " + filename);
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
          DT_ * temp = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_scalar_index.at(0)) * sizeof(DT_));
          MemoryPool<Mem_>::download(temp, this->_elements.at(0), this->_scalar_index.at(0) * sizeof(DT_));
          double * ctemp = new double[this->_scalar_index.at(0)];
          for (Index i(0) ; i < this->_scalar_index.at(0) ; ++i)
          {
            ctemp[i] = temp[i];
          }
          MemoryPool<Mem::Main>::instance()->release_memory(temp);

          uint64_t size(this->_scalar_index.at(0));
          file.write((const char *)&size, sizeof(uint64_t));
          file.write((const char *)ctemp, size * sizeof(double));

          delete[] ctemp;
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

        const DT_ * elements() const
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

          ASSERT(index < this->_scalar_index.at(0), "Error: " + stringify(index) + " exceeds dense vector size " + stringify(this->_scalar_index.at(0)) + " !");
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
          CONTEXT("When setting DenseVector element");

          ASSERT(index < this->_scalar_index.at(0), "Error: " + stringify(index) + " exceeds dense vector size " + stringify(this->_scalar_index.at(0)) + " !");
          MemoryPool<Mem_>::set_memory(this->_elements.at(0) + index, value);
        }

        /**
         * \brief Returns a descriptive string.
         *
         * \returns A string describing the container.
         */
        static String type_name()
        {
          return "DenseVector";
        }
    };

    /**
     * \brief DenseVector comparison operator
     *
     * \param[in] a A vector to compare with.
     * \param[in] b A vector to compare with.
     */
    template <typename Mem_, typename Mem2_, typename DT_> bool operator== (const DenseVector<Mem_, DT_> & a, const DenseVector<Mem2_, DT_> & b)
    {
      CONTEXT("When comparing DenseVectors");

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
     * \brief DenseVector streaming operator
     *
     * \param[in] lhs The target stream.
     * \param[in] b The vector to be streamed.
     */
    template <typename Mem_, typename DT_>
    std::ostream &
    operator<< (std::ostream & lhs, const DenseVector<Mem_, DT_> & b)
    {
      lhs << "[";
      for (Index i(0) ; i < b.size() ; ++i)
      {
        lhs << "  " << b(i);
      }
      lhs << "]" << std::endl;

      return lhs;
    }

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_DENSE_VECTOR_HPP
