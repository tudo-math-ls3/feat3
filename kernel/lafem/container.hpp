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
    namespace Intern
    {
      struct AssignStruct
      {
        template <typename MT_, typename T_, typename T2_>
          static void assign(std::vector<T_ *> & own, const std::vector<T_ *> & other)
          {
            own.assign(other.begin(), other.end());

            for (Index i(0) ; i < own.size() ; ++i)
              MemoryPool<MT_>::instance()->increase_memory(own.at(i));
          }

        template <typename MT_, typename T_, typename T2_>
          static void assign(std::vector<T_ *> &, const std::vector<T2_ *> &)
          {
            throw InternalError(__func__, __FILE__, __LINE__, "Should never be reached!");
          }
      };
    } //namespace Intern

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
    template <typename Mem_, typename DT_, typename IT_>
    class Container
    {
      private:
        /**
         * \brief Constructor
         *
         * Creates an empty container.
         *
         * \note Internal use only
         */
        Container() {}

      protected:
        /// List of pointers to all datatype dependent arrays.
        std::vector<DT_*> _elements;
        /// List of pointers to all IT_ dependent arrays.
        std::vector<IT_*> _indices;
        /// List of corresponding datatype array sizes.
        std::vector<Index> _elements_size;
        /// List of corresponding IT_ array sizes.
        std::vector<Index> _indices_size;
        /// List of scalars with datatype index.
        std::vector<Index> _scalar_index;
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

          this->_scalar_index.assign(other._scalar_index.begin(), other._scalar_index.end());
          this->_scalar_dt.assign(other._scalar_dt.begin(), other._scalar_dt.end());

          for (Index i(0) ; i < _elements.size() ; ++i)
          {
            if (_elements_size.at(i) != other.get_elements_size().at(i))
              throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");
            MemoryPool<Mem_>::template copy<DT_>(_elements.at(i), other.get_elements().at(i), _elements_size.at(i));
          }

          for (Index i(0) ; i < _indices.size() ; ++i)
          {
            if (_indices_size.at(i) != other.get_indices_size().at(i))
              throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");
            MemoryPool<Mem_>::template copy<IT_>(_indices.at(i), other.get_indices().at(i), _indices_size.at(i));
          }
        }

        template <typename Mem2_>
        void _copy_content(const Container<Mem2_, DT_, IT_> & other)
        {
          if (_elements.size() != other.get_elements().size())
            throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");
          if (_indices.size() != other.get_indices().size())
            throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");
          if (_scalar_index.size() != other.get_scalar_index().size())
            throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");
          if (_scalar_dt.size() != other.get_scalar_dt().size())
            throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");

          this->_scalar_index.assign(other.get_scalar_index().begin(), other.get_scalar_index().end());
          this->_scalar_dt.assign(other.get_scalar_dt().begin(), other.get_scalar_dt().end());

          for (Index i(0) ; i < _elements.size() ; ++i)
          {
            if (_elements_size.at(i) != other.get_elements_size().at(i))
              throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");
            if (std::is_same<Mem_, Mem::Main>::value && std::is_same<Mem2_, Mem::CUDA>::value)
              MemoryPool<Mem2_>::template download<DT_>(_elements.at(i), other.get_elements().at(i), _elements_size.at(i));
            else if (std::is_same<Mem_, Mem::CUDA>::value && std::is_same<Mem2_, Mem::Main>::value)
              MemoryPool<Mem_>::template upload<DT_>(_elements.at(i), other.get_elements().at(i), _elements_size.at(i));
            else
              throw InternalError(__func__, __FILE__, __LINE__, "Memory Backend not known!");
          }

          for (Index i(0) ; i < _indices.size() ; ++i)
          {
            if (_indices_size.at(i) != other.get_indices_size().at(i))
              throw InternalError(__func__, __FILE__, __LINE__, "Container size missmatch!");
            if (std::is_same<Mem_, Mem::Main>::value && std::is_same<Mem2_, Mem::CUDA>::value)
              MemoryPool<Mem2_>::template download<IT_>(_indices.at(i), other.get_indices().at(i), _indices_size.at(i));
            else if (std::is_same<Mem_, Mem::CUDA>::value && std::is_same<Mem2_, Mem::Main>::value)
              MemoryPool<Mem_>::template upload<IT_>(_indices.at(i), other.get_indices().at(i), _indices_size.at(i));
            else
              throw InternalError(__func__, __FILE__, __LINE__, "Memory Backend not known!");
          }
        }

        /** \brief Assignment operation
         *
         * Assigns contents of another container
         *
         * \param[in] other The source container.
         *
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void assign(const Container<Mem2_, DT2_, IT2_> & other)
        {
          CONTEXT("When assigning Container");

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
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


          if (std::is_same<Mem_, Mem2_>::value && std::is_same<DT_, DT2_>::value)
          {
            Intern::AssignStruct::template assign<Mem_, DT_, DT2_>(this->_elements, other.get_elements());
          }
          else
          {
            for (Index i(0) ; i < this->_elements_size.size() ; ++i)
            {
              const Index tsize(this->_elements_size.at(i));
              this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(tsize));

              DT_ * pthis(nullptr);
              DT2_ * pother(nullptr);
              if (std::is_same<Mem_, Mem::Main>::value)
              {
                pthis = this->_elements.at(i);
              }
              else
              {
                pthis = new DT_[tsize];
              }
              if (std::is_same<Mem2_, Mem::Main>::value)
              {
                pother = other.get_elements().at(i);
              }
              else
              {
                pother = new DT2_[tsize];
                MemoryPool<Mem2_>::template download<DT2_>(pother, other.get_elements().at(i), tsize);
              }

              for (Index j(0) ; j < tsize ; ++j)
                pthis[j] = DT_(pother[j]);

              if (! std::is_same<Mem_, Mem::Main>::value)
              {
                MemoryPool<Mem_>::template upload<DT_>(this->_elements.at(i), pthis, tsize);
                delete[] pthis;
              }
              if (!std::is_same<Mem2_, Mem::Main>::value)
                delete[] pother;
            }
          }

          if (std::is_same<Mem_, Mem2_>::value && std::is_same<IT_, IT2_>::value)
          {
            Intern::AssignStruct::template assign<Mem_, IT_, IT2_>(this->_indices, other.get_indices());
          }
          else
          {
            for (Index i(0) ; i < this->_indices_size.size() ; ++i)
            {
              const Index tsize(this->_indices_size.at(i));
              this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(tsize));

              IT_ * pthis(nullptr);
              IT2_ * pother(nullptr);

              if (std::is_same<Mem_, Mem::Main>::value)
              {
                pthis = this->_indices.at(i);
              }
              else
              {
                pthis = new IT_[tsize];
              }

              if (std::is_same<Mem2_, Mem::Main>::value)
              {
                pother = other.get_indices().at(i);
              }
              else
              {
                pother = new IT2_[tsize];
                MemoryPool<Mem2_>::template download<IT2_>(pother, other.get_indices().at(i), tsize);
              }

              for (Index j(0) ; j < tsize ; ++j)
                pthis[j] = IT_(pother[j]);

              if (! std::is_same<Mem_, Mem::Main>::value)
              {
                MemoryPool<Mem_>::template upload<IT_>(this->_indices.at(i), pthis, tsize);
                delete[] pthis;
              }
              if (!std::is_same<Mem2_, Mem::Main>::value)
                delete[] pother;
            }
          }
        }

      public:
        /**
         * \brief Constructor
         *
         * \param[in] size The size of the created container.
         *
         * Creates a container with a given size.
         */
        explicit Container(Index size_in)
        {
          CONTEXT("When creating Container");
          _scalar_index.push_back(size_in);
        }

        /**
         * \brief Destructor
         *
         * Destroys a container and releases all of its used arrays.
         */
        virtual ~Container()
        {
          CONTEXT("When destroying Container");

          for (Index i(0) ; i < _elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(_elements.at(i));
          for (Index i(0) ; i < _indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(_indices.at(i));
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
         * \brief Reset all elements of the container to a given value or zero if missing.
         *
         * \param[in] value The value to be set (defaults to 0)
         *
         */
        void format(DT_ value = DT_(0))
        {
          CONTEXT("When formating Container");

          for (Index i(0) ; i < _elements.size() ; ++i)
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
          this->_scalar_index.clear();
          this->_scalar_dt.clear();
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

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(i));
          this->_elements.clear();
          this->_indices.clear();

          for (Index i(0) ; i < other._elements.size() ; ++i)
          {
            this->_elements.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<DT_>(this->_elements_size.at(i)));
            MemoryPool<Mem_>::template copy<DT_>(this->_elements.at(i), other._elements.at(i), this->_elements_size.at(i));
          }

          for (Index i(0) ; i < other._indices.size() ; ++i)
          {
            this->_indices.push_back(MemoryPool<Mem_>::instance()->template allocate_memory<IT_>(this->_indices_size.at(i)));
            MemoryPool<Mem_>::template copy<IT_>(this->_indices.at(i), other._indices.at(i), this->_indices_size.at(i));
          }
        }

        /** \brief Clone operation
         *
         * Become a deep copy of a given container.
         *
         * \param[in] other The source container.
         *
         */
        template <typename Mem2_, typename DT2_, typename IT2_>
        void clone(const Container<Mem2_, DT2_, IT2_> & other)
        {
          CONTEXT("When cloning Container");
          Container t;
          t.assign(other);
          clone(t);
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

          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_elements.at(i));
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::instance()->release_memory(this->_indices.at(i));

          this->_elements = std::move(other._elements);
          this->_indices = std::move(other._indices);
          this->_elements_size = std::move(other._elements_size);
          this->_indices_size = std::move(other._indices_size);
          this->_scalar_index = std::move(other._scalar_index);
          this->_scalar_dt = std::move(other._scalar_dt);
        }

        /**
         * \brief Returns the total amount of bytes allocated.
         *
         * \returns The amount of bytes allocated in all arrays
         */
        Index bytes_allocated() const
        {
          Index bytes(0);

          for (Index i(0) ; i < _elements_size.size() ; ++i)
          {
            bytes += Index(_elements_size.at(i) * sizeof(DT_));
          }

          for (Index i(0) ; i < _indices_size.size() ; ++i)
          {
            bytes += Index(_indices_size.at(i) * sizeof(IT_));
          }

          bytes += Index(_scalar_index.size() * sizeof(IT_));
          bytes += Index(_scalar_dt.size() * sizeof(DT_));

          return bytes;
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
        const std::vector<IT_*> & get_indices() const
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
          if (_scalar_index.size() > 0)
            return _scalar_index.at(0);
          else
            return Index(0);
        }

        /**
         * \brief Returns the number of effective stored elements.
         *
         * \returns The number of data values.
         */
        virtual const Index & used_elements() const
        {
          return std::forward<Index>(this->size());
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
