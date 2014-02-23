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
      enum class FileMode
      {
        fm_exp = 0, /**< Exponential ascii */
        fm_dv, /**< Binary data */
        fm_m, /**< Matlab ascii */
        fm_mtx, /**< Matrix market ascii */
        fm_ell, /**< Binary ell data */
        fm_csr, /**< Binary csr data */
        fm_coo /**< Binary coo data */
      };

    /**
     * \brief Container base class.
     *
     * \tparam Mem_ The memory architecture to be used.
     * \tparam DT_ The datatype to be used.
     * \tparam IT_ The indextype to be used.
     *
     * This is the base class of all inheritated containers. \n\n
     * Data survey: \n
     * _scalar_index[0]: container size
     *
     * \author Dirk Ribbrock
     */
    template <typename Mem_, typename DT_, typename ID_>
    class Container
    {
      protected:
        /// List of pointers to all datatype dependent arrays.
        std::vector<DT_*> _elements;
        /// List of pointers to all ID_ dependent arrays.
        std::vector<ID_*> _indices;
        /// List of corresponding datatype array sizes.
        std::vector<ID_> _elements_size;
        /// List of corresponding ID_ array sizes.
        std::vector<ID_> _indices_size;
        /// List of scalars with datatype index.
        std::vector<ID_> _scalar_index;
        /// List of scalars with datatype DT_
        std::vector<DT_> _scalar_dt;

        void _copy_content(const Container & other)
        {
          // avoid self-copy
          if(this == &other)
            return;

          if (_elements.size() != other.get_elements().size())
            throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");
          if (_indices.size() != other.get_indices().size())
            throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");
          if (_scalar_index.size() != other.get_scalar_index().size())
            throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");
          if (_scalar_dt.size() != other.get_scalar_dt().size())
            throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");

          for (ID_ i(0) ; i < _elements.size() ; ++i)
          {
            if (_elements_size.at(i) != other.get_elements_size().at(i))
              throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");
            MemoryPool<Mem_>::template copy<DT_>(_elements.at(i), other.get_elements().at(i), _elements_size.at(i));
          }

          for (ID_ i(0) ; i < _indices.size() ; ++i)
          {
            if (_indices_size.at(i) != other.get_indices_size().at(i))
              throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");
            MemoryPool<Mem_>::template copy<ID_>(_indices.at(i), other.get_indices().at(i), _indices_size.at(i));
          }

          this->_scalar_index.assign(other._scalar_index.begin(), other._scalar_index.end());
          this->_scalar_dt.assign(other._scalar_dt.begin(), other._scalar_dt.end());
        }

        template <typename Mem2_>
        void _copy_content(const Container<Mem2_, DT_, ID_> & other)
        {
          Container temp(other);
          this->_copy_content(temp);
        }

      public:
        /**
         * \brief Constructor
         *
         * \param[in] size The size of the created container.
         *
         * Creates a container with a given size.
         */
        explicit Container(ID_ size)
        {
          CONTEXT("When creating Container");
          _scalar_index.push_back(size);
        }

        /**
         * \brief Destructor
         *
         * Destroys a container and releases all of its used arrays.
         */
        virtual ~Container()
        {
          CONTEXT("When destroying Container");

          for (ID_ i(0) ; i < _elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(_elements.at(i));
          for (ID_ i(0) ; i < _indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(_indices.at(i));
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source container.
         *
         * Creates a shallow copy of a given container.
         */
        Container(const Container & other) :
          _elements(other._elements),
          _indices(other._indices),
          _elements_size(other._elements_size),
          _indices_size(other._indices_size),
          _scalar_index(other._scalar_index),
          _scalar_dt(other._scalar_dt)
        {
          CONTEXT("When copying Container");

          for (ID_ i(0) ; i < _elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(_elements.at(i));
          for (ID_ i(0) ; i < _indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(_indices.at(i));
        }

        /**
         * \brief Move Constructor
         *
         * \param[in] other The source container.
         *
         * Moves another container to this container.
         */
        Container(Container && other) :
          _elements(std::move(other._elements)),
          _indices(std::move(other._indices)),
          _elements_size(std::move(other._elements_size)),
          _indices_size(std::move(other._indices_size)),
          _scalar_index(std::move(other._scalar_index)),
          _scalar_dt(std::move(other._scalar_dt))
        {
          CONTEXT("When moving Container");
          other._elements.clear();
          other._indices.clear();
          other._elements_size.clear();
          other._indices_size.clear();
          other._scalar_index.clear();
          other._scalar_dt.clear();
        }

        /**
         * \brief Copy Constructor
         *
         * \param[in] other The source container.
         *
         * Creates a copy of a given container from another memory architecture.
         */
        template <typename Mem2_, typename DT2_>
        explicit Container(const Container<Mem2_, DT2_, ID_> & other) :
          _scalar_index(other.get_scalar_index()),
          _scalar_dt(other.get_scalar_dt())
        {
          CONTEXT("When copying Container");

          if (! std::is_same<DT_, DT2_>::value)
            throw InternalError(__func__, __FILE__, __LINE__, "type conversion not supported yet!");


          for (ID_ i(0) ; i < other.get_elements().size() ; ++i)
            this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(other.get_elements_size().at(i)));
          for (ID_ i(0) ; i < other.get_indices().size() ; ++i)
            this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<ID_>(other.get_indices_size().at(i)));

          for (ID_ i(0) ; i < other.get_elements_size().size() ; ++i)
            this->_elements_size.push_back(other.get_elements_size().at(i));
          for (ID_ i(0) ; i < other.get_indices_size().size() ; ++i)
            this->_indices_size.push_back(other.get_indices_size().at(i));

          for (ID_ i(0) ; i < (other.get_elements()).size() ; ++i)
          {
            const unsigned long size(other.get_elements_size().at(i));
            if (std::is_same<Mem_, Mem2_>::value)
            {
              MemoryPool<Mem_>::template copy<DT_>(this->_elements.at(i), other.get_elements().at(i), size);
            }
            else
            {
              DT2_ * temp((DT2_*)::malloc(size * sizeof(DT2_)));
              MemoryPool<Mem2_>::template download<DT2_>(temp, other.get_elements().at(i), size);
              MemoryPool<Mem_>::template upload<DT_>(this->get_elements().at(i), temp, size);
              ::free(temp);
            }
          }
          for (ID_ i(0) ; i < other.get_indices().size() ; ++i)
          {
            const unsigned long size(other.get_indices_size().at(i));
            if (std::is_same<Mem_, Mem2_>::value)
            {
              MemoryPool<Mem_>::template copy<ID_>(this->_indices.at(i), other.get_indices().at(i), size);
            }
            else
            {
              ID_ * temp((ID_*)::malloc(size * sizeof(ID_)));
              MemoryPool<Mem2_>::template download<ID_>(temp, other.get_indices().at(i), size);
              MemoryPool<Mem_>::template upload<ID_>(this->get_indices().at(i), temp, size);
              ::free(temp);
            }
          }
        }

        /**
         * \brief Reset all elements of the container to a given value or zero if missing.
         *
         * \param[in] value The value to be set (defaults to 0)
         *
         */
        void format(DT_ value = 0)
        {
          CONTEXT("When formating Container");

          for (ID_ i(0) ; i < _elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->set_memory(_elements.at(i), value, _elements_size.at(i));
        }

        /**
         * \brief Free all allocated arrays
         *
         */
        virtual void clear()
        {
          CONTEXT("When clearing Container");

          for (Index i(0) ; i < _elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < _indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(i));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();
        }

        /** \brief Clone operation
         *
         * Become a deep copy of a given container.
         *
         * \param[in] other The source container.
         *
         */
        void clone(const Container & other)
        {
          CONTEXT("When cloning Container");

          this->_scalar_index.assign(other._scalar_index.begin(), other._scalar_index.end());
          this->_scalar_dt.assign(other._scalar_dt.begin(), other._scalar_dt.end());
          this->_elements_size.assign(other._elements_size.begin(), other._elements_size.end());
          this->_indices_size.assign(other._indices_size.begin(), other._indices_size.end());

          for (ID_ i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(i));
          for (ID_ i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(i));
          this->_elements.clear();
          this->_indices.clear();

          for (ID_ i(0) ; i < other._elements.size() ; ++i)
          {
            this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_elements_size.at(i)));
            MemoryPool<Mem_>::template copy<DT_>(this->_elements.at(i), other._elements.at(i), this->_elements_size.at(i));
          }

          for (ID_ i(0) ; i < other._indices.size() ; ++i)
          {
            this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<ID_>(this->_indices_size.at(i)));
            MemoryPool<Mem_>::template copy<ID_>(this->_indices.at(i), other._indices.at(i), this->_indices_size.at(i));
          }
        }

        /** \brief Assignment operation
         *
         * Assign another container to the current one.
         *
         * \param[in] other The source container.
         *
         */
        void assign(const Container & other)
        {
          CONTEXT("When assigning Container");

          if (this == &other)
            return;

          for (ID_ i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(i));
          for (ID_ i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(i));

          std::vector<DT_ *> new_elements = other.get_elements();
          std::vector<ID_*> new_indices = other.get_indices();

          this->_elements.assign(new_elements.begin(), new_elements.end());
          this->_indices.assign(new_indices.begin(), new_indices.end());
          this->_elements_size.assign(other.get_elements_size().begin(), other.get_elements_size().end());
          this->_indices_size.assign(other.get_indices_size().begin(), other.get_indices_size().end());
          this->_scalar_index.assign(other.get_scalar_index().begin(), other.get_scalar_index().end());
          this->_scalar_dt.assign(other.get_scalar_dt().begin(), other.get_scalar_dt().end());

          for (ID_ i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_elements.at(i));
          for (ID_ i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->increase_memory(this->_indices.at(i));
        }

        /** \brief Assignment move operation
         *
         * Move another container to the current one.
         *
         * \param[in] other The source container.
         *
         */
        void move(Container && other)
        {
          CONTEXT("When moving Container");

          if (this == &other)
            return;

          for (ID_ i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(i));
          for (ID_ i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(i));

          this->_elements = std::move(other._elements);
          this->_indices = std::move(other._indices);
          this->_elements_size = std::move(other._elements_size);
          this->_indices_size = std::move(other._indices_size);
          this->_scalar_index = std::move(other._scalar_index);
          this->_scalar_dt = std::move(other._scalar_dt);
        }

        /** \brief Assignment operation
         *
         * Assigns a container from another memory architecture to the current one.
         *
         * \param[in] other The source container.
         *
         */
        template <typename Mem2_, typename DT2_>
        void assign(const Container<Mem2_, DT2_, ID_> & other)
        {
          CONTEXT("When assigning Container");

          if (! std::is_same<DT_, DT2_>::value)
            throw InternalError(__func__, __FILE__, __LINE__, "type conversion not supported yet!");

          for (ID_ i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(i));
          for (ID_ i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(i));

          this->_elements.clear();
          this->_indices.clear();
          this->_elements_size.clear();
          this->_indices_size.clear();
          this->_scalar_index.clear();
          this->_scalar_dt.clear();

          this->_elements_size.assign(other.get_elements_size().begin(), other.get_elements_size().end());
          this->_indices_size.assign(other.get_indices_size().begin(), other.get_indices_size().end());
          this->_scalar_index.assign(other.get_scalar_index().begin(), other.get_scalar_index().end());
          this->_scalar_dt.assign(other.get_scalar_dt().begin(), other.get_scalar_dt().end());


          for (ID_ i(0) ; i < this->_elements_size.size() ; ++i)
          {
            const unsigned long size(this->_elements_size.at(i));
            this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(size));
            if (std::is_same<Mem_, Mem2_>::value)
            {
              MemoryPool<Mem_>::template copy<DT_>(this->_elements.at(i), other.get_elements().at(i), size);
            }
            else
            {
              DT2_ * temp((DT2_*)::malloc(size * sizeof(DT2_)));
              MemoryPool<Mem2_>::template download<DT2_>(temp, other.get_elements().at(i), size);
              MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(i), temp, size);
              ::free(temp);
            }
          }

          for (ID_ i(0) ; i < this->_indices_size.size() ; ++i)
          {
            const unsigned long size(this->_indices_size.at(i));
            this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<ID_>(size));
            if (std::is_same<Mem_, Mem2_>::value)
            {
              MemoryPool<Mem_>::template copy<ID_>(this->_indices.at(i), other.get_indices().at(i), size);
            }
            else
            {
              ID_ * temp((ID_*)::malloc(size * sizeof(ID_)));
              MemoryPool<Mem2_>::template download<ID_>(temp, other.get_indices().at(i), size);
              MemoryPool<Mem_>::template upload<ID_>(this->_indices.at(i), temp, size);
              ::free(temp);
            }
          }

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
         * \brief Returns a list of all ID_ arrays.
         *
         * \returns A list of all ID_ arrays.
         */
        const std::vector<ID_*> & get_indices() const
        {
          return _indices;
        }

        /**
         * \brief Returns a list of all data array sizes.
         *
         * \returns A list of all data array sizes.
         */
        const std::vector<ID_> & get_elements_size() const
        {
          return _elements_size;
        }

        /**
         * \brief Returns a list of all ID_ array sizes.
         *
         * \returns A list of all ID_ array sizes.
         */
        const std::vector<ID_> & get_indices_size() const
        {
          return _indices_size;
        }

        /**
         * \brief Returns a list of all scalar values with datatype index.
         *
         * \returns A list of all scalars with datatype index.
         */
        const std::vector<ID_> & get_scalar_index() const
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
        ID_ size() const
        {
          if (_scalar_index.size() > 0)
            return _scalar_index.at(0);
          else
            return ID_(0);
        }

        /**
         * \brief Returns the number of effective stored elements.
         *
         * \returns The number of data values.
         */
        virtual ID_ used_elements() const
        {
          return this->size();
        }

        /**
         * \brief Returns a descriptive string.
         *
         * \returns A string describing the container.
         */
        static String name()
        {
          return "Container";
        }
    };

  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_CONTAINER_HPP
