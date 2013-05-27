#pragma once
#ifndef KERNEL_LAFEM_CONTAINER_HPP
#define KERNEL_LAFEM_CONTAINER_HPP 1

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/lafem/memory_pool.hpp>
#include <kernel/archs.hpp>

#include <vector>
#include <limits>
#include <cmath>
#include <typeinfo>


namespace FEAST
{
  /**
   * \brief LAFEM namespace
   */
  namespace LAFEM
  {
      /**
       * Supported File modes.
       */
      enum FileMode
      {
        fm_exp = 0,
        fm_dv,
        fm_m,
        fm_mtx,
        fm_ell,
        fm_csr,
        fm_coo
      };

    /**
     * \brief Container base class.
     *
     * \tparam Arch_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     *
     * This is the base class of all inheritated containers.
     *
     * \author Dirk Ribbrock
     */
    template <typename Arch_, typename DT_>
    class Container
    {
      protected:
        Index _size;
        /// List of pointers to all datatype dependent arrays.
        std::vector<DT_*> _elements;
        /// List of pointers to all Index dependent arrays.
        std::vector<Index*> _indices;
        /// List of corresponding datatype array sizes.
        std::vector<Index> _elements_size;
        /// List of corresponding Index array sizes.
        std::vector<Index> _indices_size;
        /// List of scalars with datatype index.
        std::vector<Index> _scalar_index;
        /// List of scalars with datatype DT_
        std::vector<DT_> _scalar_dt;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] size The size of the created container.
         *
         * Creates a container with a given size.
         */
        Container(Index size) :
          _size(size)
        {
          CONTEXT("When creating Container");
        }

        /**
         * \brief Destructor
         *
         * Destroys a container and releases all of its used arrays.
         */
        ~Container()
        {
          CONTEXT("When destroying Container");

          for (Index i(0) ; i < _elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(_elements.at(i));
          for (Index i(0) ; i < _indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->release_memory(_indices.at(i));
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source container.
         *
         * Creates a shallow copy of a given container.
         */
        Container(const Container<Arch_, DT_> & other) :
          _size(other._size),
          _elements(other._elements),
          _indices(other._indices),
          _elements_size(other._elements_size),
          _indices_size(other._indices_size),
          _scalar_index(other._scalar_index),
          _scalar_dt(other._scalar_dt)
        {
          CONTEXT("When copying Container");

          for (Index i(0) ; i < _elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(_elements.at(i));
          for (Index i(0) ; i < _indices.size() ; ++i)
            MemoryPool<Arch_>::instance()->increase_memory(_indices.at(i));
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source container.
         *
         * Creates a copy of a given container from another memory architecture.
         */
        template <typename Arch2_, typename DT2_>
        Container(const Container<Arch2_, DT2_> & other) :
          _size(other.size()),
          _scalar_index(other.get_scalar_index()),
          _scalar_dt(other.get_scalar_dt())
        {
          CONTEXT("When copying Container");

          if (typeid(DT_) != typeid(DT2_))
            throw InternalError("type conversion not supported yet!");


          for (Index i(0) ; i < other.get_elements().size() ; ++i)
            this->_elements.push_back((DT_*)MemoryPool<Arch_>::instance()->allocate_memory(other.get_elements_size().at(i) * sizeof(DT_)));
          for (Index i(0) ; i < other.get_indices().size() ; ++i)
            this->_indices.push_back((Index*)MemoryPool<Arch_>::instance()->allocate_memory(other.get_indices_size().at(i) * sizeof(Index)));

          for (Index i(0) ; i < other.get_elements_size().size() ; ++i)
            this->_elements_size.push_back(other.get_elements_size().at(i));
          for (Index i(0) ; i < other.get_indices_size().size() ; ++i)
            this->_indices_size.push_back(other.get_indices_size().at(i));

          for (Index i(0) ; i < (other.get_elements()).size() ; ++i)
          {
            Index src_size(other.get_elements_size().at(i) * sizeof(DT2_));
            Index dest_size(other.get_elements_size().at(i) * sizeof(DT_));
            void * temp(::malloc(src_size));
            MemoryPool<Arch2_>::download(temp, other.get_elements().at(i), src_size);
            MemoryPool<Arch_>::upload(this->get_elements().at(i), temp, dest_size);
            ::free(temp);
          }
          for (Index i(0) ; i < other.get_indices().size() ; ++i)
          {
            Index src_size(other.get_indices_size().at(i) * sizeof(Index));
            Index dest_size(other.get_indices_size().at(i) * sizeof(Index));
            void * temp(::malloc(src_size));
            MemoryPool<Arch2_>::download(temp, other.get_indices().at(i), src_size);
            MemoryPool<Arch_>::upload(this->get_indices().at(i), temp, dest_size);
            ::free(temp);
          }
        }

        /**
         * \brief Reset all elements of the container to a given value or zero if missing.
         */
        void clear(DT_ value = 0)
        {
          CONTEXT("When clearing Container");

          for (Index i(0) ; i < _elements.size() ; ++i)
            MemoryPool<Arch_>::instance()->set_memory(_elements.at(i), value, _elements_size.at(i));
        }

        /**
         * \brief Returns a list of all data arrays.
         *
         * \returns A list of all data arrays.
         */
        const std::vector<DT_*> & get_elements() const
        {
          return _elements;
        }

        /**
         * \brief Returns a list of all Index arrays.
         *
         * \returns A list of all Index arrays.
         */
        const std::vector<Index*> & get_indices() const
        {
          return _indices;
        }

        /**
         * \brief Returns a list of all data array sizes.
         *
         * \returns A list of all data array sizes.
         */
        const std::vector<Index> & get_elements_size() const
        {
          return _elements_size;
        }

        /**
         * \brief Returns a list of all Index array sizes.
         *
         * \returns A list of all Index array sizes.
         */
        const std::vector<Index> & get_indices_size() const
        {
          return _indices_size;
        }

        /**
         * \brief Returns a list of all scalar values with datatype index.
         *
         * \returns A list of all scalars with datatype index.
         */
        const std::vector<Index> & get_scalar_index() const
        {
          return _scalar_index;
        }

        /**
         * \brief Returns a list of all scalar values with datatype dt.
         *
         * \returns A list of all scalars with datatype dt.
         */
        const std::vector<DT_> & get_scalar_dt() const
        {
          return _scalar_dt;
        }

        /**
         * \brief Returns the containers size.
         *
         * \returns The containers size.
         */
        Index size() const
        {
          return _size;
        }

        /**
         * \brief Returns a descriptive string.
         *
         * \returns A string describing the container.
         */
        static String type_name()
        {
          return "Container";
        }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_CONTAINER_HPP
