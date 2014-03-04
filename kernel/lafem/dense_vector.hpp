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
#include <kernel/lafem/edi.hpp>
#include <kernel/lafem/arch/sum.hpp>
#include <kernel/lafem/arch/difference.hpp>
#include <kernel/lafem/arch/dot_product.hpp>
#include <kernel/lafem/arch/norm.hpp>
#include <kernel/lafem/arch/scale.hpp>
#include <kernel/lafem/arch/axpy.hpp>
#include <kernel/lafem/arch/component_product.hpp>

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
     * \tparam IT_ The indextype to be used.
     *
     * This class represents a vector of continuous data in memory. \n \n
     * Data survey: \n
     * _elements[0]: raw number values \n
     * _scalar_index[0]: container size
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename IT_ = Index>
    class DenseVector : public Container<Mem_, DT_, IT_>, public VectorBase
    {
      private:
        void _read_from_exp(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in);
          if (! file.is_open())
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
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

            DT_ n_z((DT_)atof(n_z_s.c_str()));

            data.push_back(n_z);

          }

          this->_scalar_index.at(0) = Index(data.size());
          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(Index(data.size())));
          this->_elements_size.push_back(Index(data.size()));
          MemoryPool<Mem_>::instance()->template upload<DT_>(this->_elements.at(0), &data[0], Index(data.size()));
        }

        void _read_from_dv(String filename)
        {
          std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
          if (! file.is_open())
            throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
          _read_from_dv(file);
          file.close();
        }

        void _read_from_dv(std::istream& file)
        {
          uint64_t size;
          file.read((char *)&size, (long)(sizeof(uint64_t)));
          this->_scalar_index.at(0) = (Index)size;
          this->_elements_size.push_back(this->_scalar_index.at(0));

          double * ctemp = new double[std::size_t(size)];
          file.read((char *)ctemp, (long)(size * sizeof(double)));

          DT_ * temp = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>((this->_scalar_index.at(0)));
          for (Index i(0) ; i < size ; ++i)
          {
            temp[i] = (DT_)ctemp[i];
          }
          delete[] ctemp;
          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_scalar_index.at(0)));
          MemoryPool<Mem_>::instance()->template upload<DT_>(this->_elements.at(0), temp, this->_scalar_index.at(0));
          MemoryPool<Mem::Main>::instance()->release_memory(temp);
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
        explicit DenseVector() :
          Container<Mem_, DT_, IT_> (0)
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
          Container<Mem_, DT_, IT_>(size)
        {
          CONTEXT("When creating DenseVector");

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(size));
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
          Container<Mem_, DT_, IT_>(size)
        {
          CONTEXT("When creating DenseVector");

          this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(size));
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
          Container<Mem_, DT_, IT_>(size)
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
          Container<Mem_, DT_, IT_>(0)
        {
          CONTEXT("When creating DenseVector");

          switch(mode)
          {
            case FileMode::fm_exp:
              _read_from_exp(filename);
              break;
            case FileMode::fm_dv:
              _read_from_dv(filename);
              break;
            default:
              throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
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
          Container<Mem_, DT_, IT_>(0)
        {
          CONTEXT("When creating DenseVector");

          switch(mode)
          {
            case FileMode::fm_exp:
              _read_from_exp(file);
              break;
            case FileMode::fm_dv:
              _read_from_dv(file);
              break;
            default:
              throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
          }
        }

        /**
         * \brief Move Constructor
         *
         * \param[in] other The source vector.
         *
         * Moves another vector to this vector.
         */
        DenseVector(DenseVector && other) :
          Container<Mem_, DT_, IT_>(std::move(other))
        {
          CONTEXT("When moving DenseVector");
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

          this->move(std::move(other));

          return *this;
        }

        /**
         * \brief Convertion method
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
          DT_ * temp = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>((this->_scalar_index.at(0)));
          MemoryPool<Mem_>::template download<DT_>(temp, this->_elements.at(0), this->_scalar_index.at(0));

          for (Index i(0) ; i < this->_scalar_index.at(0) ; ++i)
          {
            file << std::scientific << temp[i] << std::endl;
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
            std::cout<<"Warning: You are writing out an dense vector with less than double precission!"<<std::endl;

          const Index csize(this->_scalar_index.at(0));
          DT_ * temp = MemoryPool<Mem::Main>::instance()->template allocate_memory<DT_>(csize);
          MemoryPool<Mem_>::template download<DT_>(temp, this->_elements.at(0), csize);
          double * ctemp = new double[csize];
          for (Index i(0) ; i < csize ; ++i)
          {
            ctemp[i] = (double)temp[i];
          }
          MemoryPool<Mem::Main>::instance()->release_memory(temp);

          uint64_t size(csize);
          file.write((const char *)&size, sizeof(uint64_t));
          file.write((const char *)ctemp, (long)(size * sizeof(double)));

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
         * \brief Create temporary object for direct data manipulation.
         * \warning Be aware, that any synchronisation only takes place, when the object is destroyed!
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

        /**
         * \brief Calculate \f$this \leftarrow \alpha x + y\f$
         *
         * \param[in] x The first summand vector to be scaled.
         * \param[in] y The second summand vector
         * \param[in] alpha A scalar to multiply x with.
         */
        template <typename Algo_>
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
            Arch::Sum<Mem_, Algo_>::value(this->elements(), x.elements(), y.elements(), this->size());
          // r <- y - x
          else if(Math::abs(alpha + DT_(1)) < Math::eps<DT_>())
            Arch::Difference<Mem_, Algo_>::value(this->elements(), y.elements(), x.elements(), this->size());
          // r <- y
          else if(Math::abs(alpha) < Math::eps<DT_>())
            this->copy(y);
          // r <- y + alpha*x
          else
            Arch::Axpy<Mem_, Algo_>::dv(this->elements(), alpha, x.elements(), y.elements(), this->size());
        }

        /**
         * \brief Calculate \f$this_i \leftarrow x_i \cdot y_i\f$
         *
         * \param[in] x The first factor.
         * \param[in] y The second factor.
         */
        template <typename Algo_>
        void component_product(const DenseVector & x, const DenseVector & y)
        {
          if (this->size() != x.size())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");
          if (this->size() != y.size())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

          Arch::ComponentProduct<Mem_, Algo_>::value(this->elements(), x.elements(), y.elements(), this->size());
        }

        /**
         * \brief Calculate \f$this_i \leftarrow x_i \cdot y_i + z_i\f$
         *
         * \param[in] x The first factor.
         * \param[in] y The second factor.
         * \param[in] z The second summand.
         */
        template <typename Algo_>
        void component_product(
          const DenseVector & x,
          const DenseVector & y,
          const DenseVector & z)
        {
          if (this->size() != x.size())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");
          if (this->size() != y.size())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");
          if (this->size() != z.size())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

          Arch::Axpy<Mem_, Algo_>::dv(this->elements(), x.elements(), y.elements(), z.elements(), this->size());
        }

        /**
         * \brief Calculate \f$this \leftarrow \alpha x \f$
         *
         * \param[in] x The vector to be scaled.
         * \param[in] alpha A scalar to scale x with.
         */
        template <typename Algo_>
        void scale(const DenseVector & x, const DT_ alpha)
        {
          if (x.size() != this->size())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

          Arch::Scale<Mem_, Algo_>::value(this->elements(), x.elements(), alpha, this->size());
        }

        /**
         * \brief Calculate \f$r \leftarrow this \cdot x\f$
         *
         * \param[out] r The dot product result.
         * \param[in] x The other vector.
         */
        template <typename Algo_>
        DataType dot(const DenseVector & x) const
        {
          if (x.size() != this->size())
            throw InternalError(__func__, __FILE__, __LINE__, "Vector size does not match!");

          return Arch::DotProduct<Mem_, Algo_>::value(this->elements(), x.elements(), this->size());
        }

        /**
         * \brief Calculates and returns the euclid norm of this vector.
         */
        template <typename Algo_>
        DT_ norm2() const
        {
          return Arch::Norm2<Mem_, Algo_>::value(this->elements(), this->size());
        }

        /**
         * \brief Calculates and returns the squared euclid norm of this vector.
         */
        template <typename Algo_>
        DT_ norm2sqr() const
        {
          // fallback
          return Math::sqr(this->norm2<Algo_>());
        }
    }; // class DenseVector<...>


    /**
     * \brief DenseVector comparison operator
     *
     * \param[in] a A vector to compare with.
     * \param[in] b A vector to compare with.
     */
    template <typename Mem_, typename Mem2_, typename DT_, typename IT_> bool operator== (const DenseVector<Mem_, DT_, IT_> & a, const DenseVector<Mem2_, DT_, IT_> & b)
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
    template <typename Mem_, typename DT_, typename IT_>
    std::ostream &
    operator<< (std::ostream & lhs, const DenseVector<Mem_, DT_, IT_> & b)
    {
      lhs << "[";
      for (Index i(0) ; i < b.size() ; ++i)
      {
        lhs << "  " << b(i);
      }
      lhs << "]";

      return lhs;
    }

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_DENSE_VECTOR_HPP
