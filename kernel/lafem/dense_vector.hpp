#pragma once
#ifndef KERNEL_LAFEM_DENSE_VECTOR_HPP
#define KERNEL_LAFEM_DENSE_VECTOR_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/lafem/container.hpp>

#include <iostream>
#include <fstream>
#include <string>


namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Dense data vector class template.
     *
     * \tparam Arch_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This class represents a vector of continuous data in memory.
     *
     * \author Dirk Ribbrock
     */
    template <typename Arch_, typename DT_>
    class DenseVector : public Container<Arch_, DT_>
    {
      private:
        /// Pointer to our elements.
        DT_ * _pelements;

      public:
        /// Our datatype
        typedef DT_ DataType;
        /// Our memory architecture type
        typedef Arch_ MemType;

        /**
         * \brief Constructor
         *
         * Creates an empty non dimensional vector.
         */
        explicit DenseVector() :
          Container<Arch_, DT_> (0)
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
          Container<Arch_, DT_>(size)
        {
          CONTEXT("When creating DenseVector");

          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(size * sizeof(DT_)));
          this->_elements_size.push_back(size);
          this->_pelements = this->_elements.at(0);
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
          Container<Arch_, DT_>(size)
        {
          CONTEXT("When creating DenseVector");

          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(size * sizeof(DT_)));
          this->_elements_size.push_back(size);
          this->_pelements = this->_elements.at(0);

          MemoryPool<Arch_>::instance()->set_memory(_pelements, value, size);
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
          Container<Arch_, DT_>(size)
        {
          CONTEXT("When creating DenseVector");

          this->_elements.push_back(data);
          this->_elements_size.push_back(size);
          this->_pelements = this->_elements.at(0);

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_indices.at(i));
        }

        /**
         * \brief Constructor
         *
         * \param[in] filename The source file in EXP format.
         *
         * Creates a vector from the given source file.
         */
        explicit DenseVector(String filename) :
          Container<Arch_, DT_>(0)
        {
          CONTEXT("When creating DenseVector");

          std::vector<DT_> data;

          std::ifstream file(filename.c_str(), std::ifstream::in);
          if (! file.is_open())
            throw InternalError("Unable to open Vector file " + filename);

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
          file.close();

          this->_size = Index(data.size());
          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(Index(data.size() * sizeof(DT_))));
          this->_elements_size.push_back(Index(data.size()));
          this->_pelements = this->_elements.at(0);
          MemoryPool<Arch_>::instance()->upload(this->_pelements, &data[0], Index(data.size() * sizeof(DT_)));
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source vector.
         *
         * Creates a shallow copy of a given vector.
         */
        DenseVector(const DenseVector<Arch_, DT_> & other) :
          Container<Arch_, DT_>(other)
        {
          CONTEXT("When copying DenseVector");

          this->_pelements = this->_elements.at(0);
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source vector.
         *
         * Creates a copy of a given vector from another memory architecture.
         */
        template <typename Arch2_, typename DT2_>
        DenseVector(const DenseVector<Arch2_, DT2_> & other) :
            Container<Arch_, DT_>(other)
        {
          CONTEXT("When copying DenseVector");

          this->_pelements = this->_elements.at(0);
        }

        /** \brief Clone operation
         *
         * Creates a deep copy of this vector.
         */
        DenseVector<Arch_, DT_> clone()
        {
          CONTEXT("When cloning DenseVector");

          DenseVector<Arch_, DT_> t(this->_size);

          void * pdest(t.elements());
          const void * psrc(this->elements());
          MemoryPool<Arch_>::copy(pdest, psrc, this->_size * sizeof(DT_));

          return t;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source vector.
         *
         * Assigns another vector to the target vector.
         */
        DenseVector<Arch_, DT_> & operator= (const DenseVector<Arch_, DT_> & other)
        {
          CONTEXT("When assigning DenseVector");

          if (this == &other)
            return *this;

          this->_size = other.size();

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_indices.at(i));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();

          std::vector<DT_ *> new_elements = other.get_elements();
          std::vector<unsigned long *> new_indices = other.get_indices();

          this->_elements.assign(new_elements.begin(), new_elements.end());
          this->_indices.assign(new_indices.begin(), new_indices.end());
          this->_elements_size.assign(other.get_elements_size().begin(), other.get_elements_size().end());
          this->_indices_size.assign(other.get_indices_size().begin(), other.get_indices_size().end());

          _pelements = this->_elements.at(0);

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(this->_indices.at(i));

          return *this;
        }

        /**
         * \brief Assignment operator
         *
         * \param[in] other The source vector.
         *
         * Assigns a vector from another memory architecture to the target vector.
         */
        template <typename Arch2_, typename DT2_>
        DenseVector<Arch_, DT_> & operator= (const DenseVector<Arch2_, DT2_> & other)
        {
          CONTEXT("When assigning DenseVector");

          this->_size = other.size();

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(this->_indices.at(i));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();


          this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(other.size() * sizeof(DT_)));
          this->_elements_size.push_back(this->_size);
          this->_pelements = this->_elements.at(0);

          Index src_size(other.get_elements_size().at(0) * sizeof(DT2_));
          Index dest_size(other.get_elements_size().at(0) * sizeof(DT_));
          void * temp(::malloc(src_size));
          MemoryPool<Arch2_>::download(temp, other.get_elements().at(0), src_size);
          MemoryPool<Arch_>::upload(this->get_elements().at(0), temp, dest_size);
          ::free(temp);

          return *this;
        }

        /**
         * \brief Write out vector to file.
         *
         * \param[in] filename The file where the vector shall be stored.
         */
        void write_out(String filename) const
        {
          CONTEXT("When writing out DenseVector");

          DT_ * temp = (DT_*)MemoryPool<Mem::Main>::instance()->allocate_memory((this->_size) * sizeof(DT_));
          MemoryPool<Arch_>::download(temp, _pelements, this->_size * sizeof(DT_));

          std::ofstream file(filename.c_str(), std::ofstream::out);
          if (! file.is_open())
            throw InternalError("Unable to open Vector file " + filename);

          for (Index i(0) ; i < this->_size ; ++i)
          {
            file << std::scientific << (double)temp[i] << std::endl;
          }

          file.close();
          MemoryPool<Mem::Main>::instance()->release_memory(temp);
        }

        /**
         * \brief Get a pointer to the data array.
         *
         * \returns Pointer to the data array.
         */
        DT_ * elements()
        {
          return _pelements;
        }

        const DT_ * elements() const
        {
          return _pelements;
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

          ASSERT(index < this->_size, "Error: " + stringify(index) + " exceeds dense vector size " + stringify(this->_size) + " !");
          return MemoryPool<Arch_>::get_element(_pelements, index);
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

          ASSERT(index < this->_size, "Error: " + stringify(index) + " exceeds dense vector size " + stringify(this->_size) + " !");
          MemoryPool<Arch_>::modify_element(_pelements, index, value);
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
    template <typename Arch_, typename Arch2_, typename DT_> bool operator== (const DenseVector<Arch_, DT_> & a, const DenseVector<Arch2_, DT_> & b)
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
    template <typename Arch_, typename DT_>
    std::ostream &
    operator<< (std::ostream & lhs, const DenseVector<Arch_, DT_> & b)
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
