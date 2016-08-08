#pragma once
#ifndef KERNEL_LAFEM_CONTAINER_HPP
#define KERNEL_LAFEM_CONTAINER_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/base.hpp>

#include <vector>
#include <limits>
#include <cmath>
#include <typeinfo>
#include <string>
#include <type_traits>
#include <cstdlib>
#include <stdint.h>

#define InsertWeakClone( TContainer )                             \
  /** \brief Clone operation                                      \
   *                                                              \
   * Create a weak clone of this container.                       \
   *                                                              \
   * \param[in] clone_mode The actual cloning procedure.          \
   *                                                              \
   */                                                             \
  TContainer clone(CloneMode clone_mode = CloneMode::Weak) const  \
  {                                                               \
    TContainer t;                                                 \
    t.clone(*this, clone_mode);                                   \
    return t;                                                     \
  }                                                               \
  using Container<Mem_, DT_, IT_>::clone;

#define InsertDeepClone( TContainer )                                   \
  /** \brief Clone operation                                            \
   *                                                                    \
   * Create a deep clone of this container.                             \
   *                                                                    \
   * \param[in] clone_mode The actual cloning procedure.                \
   *                                                                    \
   */                                                                   \
  TContainer clone(CloneMode clone_mode = CloneMode::Deep) const        \
  {                                                                     \
    TContainer t;                                                       \
    t.clone(*this, clone_mode);                                         \
    return t;                                                           \
  }                                                                     \
  /** \brief Clone operation                                            \
   *                                                                    \
   * Create a deep clone of this container.                             \
   *                                                                    \
   * \param[in] other The source container to create the clone from.    \
   * \param[in] clone_mode The actual cloning procedure.                \
   *                                                                    \
   */                                                                   \
  template<typename Mem2_, typename DT2_, typename IT2_>                \
  void clone(const TContainer<Mem2_, DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Deep) \
  {                                                                     \
    Container<Mem_, DT_, IT_>::clone(other, clone_mode);                \
  }                                                                     \


namespace FEAT
{
  /**
   * \brief LAFEM namespace
   */
  namespace LAFEM
  {
    /// \cond internal
    namespace Intern
    {
      struct AssignStruct
      {
        template <typename MT_, typename T1_>
        static void assign(std::vector<T1_ *> & own, const std::vector<T1_ *> & other)
        {
          own.assign(other.begin(), other.end());

          for (Index i(0) ; i < own.size() ; ++i)
            MemoryPool<MT_>::increase_memory(own.at(i));
        }

        template <typename MT_, typename T1_, typename T2_>
        static void assign(std::vector<T1_ *> &, const std::vector<T2_ *> &)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Should never be reached!");
        }

        template <typename T_>
        static void assign_scalar(std::vector<T_> & own, const std::vector<T_> & other)
        {
          own(other.begin(), other.end());
        }

        template <typename T1_, typename T2_>
        static void assign_scalar(std::vector<T1_> & own, const std::vector<T2_> & other)
        {
          for(auto t : other)
          {
            own.push_back(T1_(t));
          }
        }
      };
    } //namespace Intern
    /// \endcond

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
      template <typename Mem2_, typename DT2_, typename IT2_>
      friend class Container;

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
      /// do we use memory that we did not allocate, nor are we allowed to free it - this mostly holds true, if the container is a ranged slice of another one
      bool _foreign_memory;

      void _copy_content(const Container & other, bool full)
      {
        // avoid self-copy
        if(this == &other)
          return;

        if (_elements.size() != other.get_elements().size())
          throw InternalError(__func__, __FILE__, __LINE__, "Container size mismatch!");
        if (_indices.size() != other.get_indices().size())
          throw InternalError(__func__, __FILE__, __LINE__, "Container size mismatch!");
        if (_scalar_index.size() != other.get_scalar_index().size())
          throw InternalError(__func__, __FILE__, __LINE__, "Container size mismatch!");
        if (_scalar_dt.size() != other.get_scalar_dt().size())
          throw InternalError(__func__, __FILE__, __LINE__, "Container size mismatch!");

        if (full)
        {
          this->_scalar_index.assign(other._scalar_index.begin(), other._scalar_index.end());
          this->_scalar_dt.assign(other._scalar_dt.begin(), other._scalar_dt.end());

          for (Index i(0) ; i < _indices.size() ; ++i)
          {
            if (_indices_size.at(i) != other.get_indices_size().at(i))
              throw InternalError(__func__, __FILE__, __LINE__, "Container size mismatch!");
            MemoryPool<Mem_>::template copy<IT_>(_indices.at(i), other.get_indices().at(i), _indices_size.at(i));
          }
        }

        for (Index i(0) ; i < _elements.size() ; ++i)
        {
          if (_elements_size.at(i) != other.get_elements_size().at(i))
            throw InternalError(__func__, __FILE__, __LINE__, "Container size mismatch!");
          MemoryPool<Mem_>::template copy<DT_>(_elements.at(i), other.get_elements().at(i), _elements_size.at(i));
        }

      }

      template <typename Mem2_>
      void _copy_content(const Container<Mem2_, DT_, IT_> & other, bool full)
      {
        if (_elements.size() != other.get_elements().size())
          throw InternalError(__func__, __FILE__, __LINE__, "Container size mismatch!");
        if (_indices.size() != other.get_indices().size())
          throw InternalError(__func__, __FILE__, __LINE__, "Container size mismatch!");
        if (_scalar_index.size() != other.get_scalar_index().size())
          throw InternalError(__func__, __FILE__, __LINE__, "Container size mismatch!");
        if (_scalar_dt.size() != other.get_scalar_dt().size())
          throw InternalError(__func__, __FILE__, __LINE__, "Container size mismatch!");

        if (full)
        {
          this->_scalar_index.assign(other.get_scalar_index().begin(), other.get_scalar_index().end());
          this->_scalar_dt.assign(other.get_scalar_dt().begin(), other.get_scalar_dt().end());

          for (Index i(0) ; i < _indices.size() ; ++i)
          {
            if (_indices_size.at(i) != other.get_indices_size().at(i))
              throw InternalError(__func__, __FILE__, __LINE__, "Container size mismatch!");
            if (std::is_same<Mem_, Mem::Main>::value && std::is_same<Mem2_, Mem::CUDA>::value)
              MemoryPool<Mem2_>::template download<IT_>(_indices.at(i), other.get_indices().at(i), _indices_size.at(i));
            else if (std::is_same<Mem_, Mem::CUDA>::value && std::is_same<Mem2_, Mem::Main>::value)
              MemoryPool<Mem_>::template upload<IT_>(_indices.at(i), other.get_indices().at(i), _indices_size.at(i));
            else
              throw InternalError(__func__, __FILE__, __LINE__, "Memory Backend not known!");
          }
        }

        for (Index i(0) ; i < _elements.size() ; ++i)
        {
          if (_elements_size.at(i) != other.get_elements_size().at(i))
            throw InternalError(__func__, __FILE__, __LINE__, "Container size mismatch!");
          if (std::is_same<Mem_, Mem::Main>::value && std::is_same<Mem2_, Mem::CUDA>::value)
            MemoryPool<Mem2_>::template download<DT_>(_elements.at(i), other.get_elements().at(i), _elements_size.at(i));
          else if (std::is_same<Mem_, Mem::CUDA>::value && std::is_same<Mem2_, Mem::Main>::value)
            MemoryPool<Mem_>::template upload<DT_>(_elements.at(i), other.get_elements().at(i), _elements_size.at(i));
          else
            throw InternalError(__func__, __FILE__, __LINE__, "Memory Backend not known!");
        }
      }

      /**
       * \brief Assignment operation
       *
       * Assigns contents of another container
       *
       * \param[in] other The source container.
       *
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void assign(const Container<Mem2_, DT2_, IT2_> & other)
      {
        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Mem_>::release_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Mem_>::release_memory(this->_indices.at(i));

        this->_elements.clear();
        this->_indices.clear();
        this->_elements_size.clear();
        this->_indices_size.clear();
        this->_scalar_index.clear();
        this->_scalar_dt.clear();

        this->_elements_size.assign(other.get_elements_size().begin(), other.get_elements_size().end());
        this->_indices_size.assign(other.get_indices_size().begin(), other.get_indices_size().end());
        this->_scalar_index.assign(other.get_scalar_index().begin(), other.get_scalar_index().end());
        Intern::AssignStruct::template assign_scalar<DT_, DT2_>(this->_scalar_dt, other.get_scalar_dt());


        if (std::is_same<Mem_, Mem2_>::value && std::is_same<DT_, DT2_>::value)
        {
          Intern::AssignStruct::template assign<Mem_, DT_>(this->_elements, other.get_elements());
        }
        else if (std::is_same<Mem_, Mem2_>::value)
        {
          for (Index i(0) ; i < this->_elements_size.size() ; ++i)
          {
            const Index tsize(this->_elements_size.at(i));
            this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(tsize));
            MemoryPool<Mem_>::convert(this->_elements.at(i), other.get_elements().at(i), tsize);
          }
        }
        else
        {
          for (Index i(0) ; i < this->_elements_size.size() ; ++i)
          {
            const Index tsize(this->_elements_size.at(i));
            this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(tsize));

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

            MemoryPool<Mem::Main>::convert(pthis, pother, tsize);

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
          Intern::AssignStruct::template assign<Mem_, IT_>(this->_indices, other.get_indices());
        }
        else if (std::is_same<Mem_, Mem2_>::value)
        {
          for (Index i(0) ; i < this->_indices_size.size() ; ++i)
          {
            const Index tsize(this->_indices_size.at(i));
            this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(tsize));
            MemoryPool<Mem_>::convert(this->_indices.at(i), other.get_indices().at(i), tsize);
          }
        }
        else
        {
          for (Index i(0) ; i < this->_indices_size.size() ; ++i)
          {
            const Index tsize(this->_indices_size.at(i));
            this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(tsize));

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

            MemoryPool<Mem::Main>::convert(pthis, pother, tsize);

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

      /**
       * \brief Serialisation of complete container entity.
       *
       * \param[in] mode FileMode enum, describing the actual container specialisation.
       * \param[out] std::vector<char> A std::vector, containing the byte array.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * Array data layout:
       * \code
       * raw array size in bytes (uint64_t)
       * magic number (uint64_t)
       * _elements.size() (uint64_t)
       * _indices.size() (uint64_t)
       * _elements_size.size() (uint64_t)
       * _indices_size.size() (uint64_t)
       * _scalar_index.size() (uint64_t)
       * _scalar_dt.size() (uint64_t)
       *
       * _elements_size array (uint64_t)
       * _indices_size array (uint64_t)
       * _scalar_index array (uint64_t)
       * _scalar_dt array (DT2_)
       * _elements arrays (DT2_)
       * _indices arrays (IT2_)
       * \endcode
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> _serialise(FileMode mode) const
      {
        Container<Mem::Main, DT2_, IT2_> tc(0);
        tc.assign(*this);

        uint64_t gsize(2 * sizeof(uint64_t)); //raw array size + magic number
        gsize += 6 * sizeof(uint64_t); // size of all six stl containers
        gsize += tc._elements_size.size() * sizeof(uint64_t); // _elements_size contents
        gsize += tc._indices_size.size() * sizeof(uint64_t); // _indices_size contents
        gsize += tc._scalar_index.size() * sizeof(uint64_t); // _scalar_index contents
        gsize += tc._scalar_dt.size() * sizeof(DT2_); // _scalar_dt contents

        for (Index i(0) ; i < tc._elements_size.size() ; ++i)
        {
          gsize += tc._elements_size.at(i) * sizeof(DT2_); // actual array sizes
        }
        for (Index i(0) ; i < tc._indices_size.size() ; ++i)
        {
          gsize += tc._indices_size.at(i) * sizeof(IT2_); // actual array sizes
        }
        gsize +=16; //padding for datatype alignment mismatch

        std::vector<char> result((size_t(gsize)));
        char * array(result.data());
        uint64_t * uiarray(reinterpret_cast<uint64_t *>(array));
        DT2_ * dtarray(reinterpret_cast<DT2_ *>(array));
        IT2_ * itarray(reinterpret_cast<IT2_ *>(array));

        /// \compilerhack clang/icc seems to use older libc++ without std::underlying_type
#if defined(FEAT_COMPILER_CLANG) || defined(FEAT_COMPILER_INTEL)
        uint64_t magic = (uint64_t)static_cast<__underlying_type(FileMode)>(mode);
#else
        uint64_t magic = (uint64_t)static_cast<typename std::underlying_type<FileMode>::type>(mode);
#endif
        uiarray[0] = gsize;
        uiarray[1] = magic;
        uiarray[2] = tc._elements.size();
        uiarray[3] = tc._indices.size();
        uiarray[4] = tc._elements_size.size();
        uiarray[5] = tc._indices_size.size();
        uiarray[6] = tc._scalar_index.size();
        uiarray[7] = tc._scalar_dt.size();

        Index global_i(8); // count how many elements of uint64_t have been inserted so far

        for (Index i(0) ; i < tc._elements_size.size() ; ++i)
        {
          uiarray[i + global_i] = tc._elements_size.at(i);
        }
        global_i += Index(tc._elements_size.size());

        for (Index i(0) ; i < tc._indices_size.size() ; ++i)
        {
          uiarray[i + global_i] = tc._indices_size.at(i);
        }
        global_i += Index(tc._indices.size());

        for (Index i(0) ; i < tc._scalar_index.size() ; ++i)
        {
          uiarray[i + global_i] = tc._scalar_index.at(i);
        }
        global_i += Index(tc._scalar_index.size());

        global_i = (Index)std::ceil(((float)global_i * (float)sizeof(uint64_t)) / (float)sizeof(DT2_)); // now counting how many DT2 have been inserted so far

        for (Index i(0) ; i < tc._scalar_dt.size() ; ++i)
        {
          dtarray[i + global_i] = tc._scalar_dt.at(i);
        }
        global_i += Index(tc._scalar_dt.size());

        for (Index i(0) ; i < tc._elements.size() ; ++i)
        {
          std::memcpy(&dtarray[global_i], tc._elements.at(i), tc._elements_size.at(i) * sizeof(DT2_));
          global_i += tc._elements_size.at(i);
        }

        global_i = (Index)std::ceil(((float)global_i * (float)sizeof(DT2_)) / (float)sizeof(IT2_)); // now counting IT2_ elements
        for (Index i(0) ; i < tc._indices.size() ; ++i)
        {
          std::memcpy(&itarray[global_i], tc._indices.at(i), tc._indices_size.at(i) * sizeof(IT2_));
          global_i += tc._indices_size.at(i);
        }

        return result;
      }

      /**
       * \brief Serialisation of complete container entity.
       *
       * \param[in] mode FileMode enum, describing the actual container specialisation.
       * \param[in] file The output stream to write data into.
       *
       * Serialize a complete container entity into a single binary file.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void _serialise(FileMode mode, std::ostream & file) const
      {
        auto temp(this->template _serialise<DT2_, IT2_>(mode));
        file.write(temp.data(), long(temp.size()));
        if (!file.good())
          throw InternalError(__func__, __FILE__, __LINE__, "Error in _serialise - file ostream is not good anymore!");
      }

      /**
       * \brief Deserialisation of complete container entity.
       *
       * \param[in] mode FileMode enum, describing the actual container specialisation.
       * \param[in] std::vector<char> A std::vector containing the byte array.
       *
       * Recreate a complete container entity by a single binary array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void _deserialise(FileMode mode, std::vector<char> & input)
      {
        this->clear();
        Container<Mem::Main, DT2_, IT2_> tc(0);
        tc.clear();

        char * array(input.data());
        uint64_t * uiarray(reinterpret_cast<uint64_t *>(array));
        DT2_ * dtarray(reinterpret_cast<DT2_ *>(array));
        IT2_ * itarray(reinterpret_cast<IT2_ *>(array));

        // clang/icc seems to use older libc++ without std::underlying_type
#if defined(FEAT_COMPILER_CLANG) || defined(FEAT_COMPILER_INTEL)
        uint64_t magic = (uint64_t)static_cast<__underlying_type(FileMode)>(mode);
#else
        uint64_t magic = (uint64_t)static_cast<typename std::underlying_type<FileMode>::type>(mode);
#endif
        if (magic != uiarray[1])
          throw InternalError(__func__, __FILE__, __LINE__, "_deserialise: given FileMode incompatible with given array!");


        Index global_i(8);
        for (uint64_t i(0) ; i < uiarray[4] ; ++i)
        {
          tc._elements_size.push_back(Index(uiarray[i + global_i]));
        }
        global_i += Index(uiarray[4]);

        for (uint64_t i(0) ; i < uiarray[5] ; ++i)
        {
          tc._indices_size.push_back(Index(uiarray[i + global_i]));
        }
        global_i += Index(uiarray[5]);

        for (uint64_t i(0) ; i < uiarray[6] ; ++i)
        {
          tc._scalar_index.push_back(Index(uiarray[i + global_i]));
        }
        global_i += Index(uiarray[6]);

        global_i = (Index)std::ceil(((float)global_i * (float)sizeof(uint64_t)) / (float)sizeof(DT2_));
        for (uint64_t i(0) ; i < uiarray[7] ; ++i)
        {
          tc._scalar_dt.push_back(dtarray[i + global_i]);
        }
        global_i += Index(uiarray[7]);

        for (Index i(0) ; i < Index(uiarray[2]) ; ++i)
        {
          tc._elements.push_back(MemoryPool<Mem::Main>::template allocate_memory<DT2_>(tc._elements_size.at(i)));
          MemoryPool<Mem::Main>::template upload<DT2_>(tc._elements.at(i), &dtarray[global_i], tc._elements_size.at(i));
          global_i += tc._elements_size.at(i);
        }

        global_i = (Index)std::ceil(((float)global_i * (float)sizeof(DT2_)) / (float)sizeof(IT2_));
        for (Index i(0) ; i < Index(uiarray[3]) ; ++i)
        {
          tc._indices.push_back(MemoryPool<Mem::Main>::template allocate_memory<IT2_>(tc._indices_size.at(i)));
          MemoryPool<Mem::Main>::template upload<IT2_>(tc._indices.at(i), &itarray[global_i], tc._indices_size.at(i));
          global_i += tc._indices_size.at(i);
        }

        this->assign(tc);
      }

      /**
       * \brief Deserialisation of complete container entity.
       *
       * \param[in] mode FileMode enum, describing the actual container specialisation.
       * \param[in] file std::istream, pointing to the input data.
       *
       * Recreate a complete container entity by a single binary file.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void _deserialise(FileMode mode, std::istream & file)
      {
        uint64_t tsize;
        file.read((char *)&tsize, (long)(sizeof(uint64_t)));
        std::vector<char> temp((size_t(tsize)));
        file.seekg(-(long)sizeof(uint64_t), file.cur);
        file.read(temp.data(), (long)(tsize));
        if (!file.good())
          throw InternalError(__func__, __FILE__, __LINE__, "Error in _deserialise - file istream is not good anymore!");
        this->template _deserialise<DT2_, IT2_>(mode, temp);
      }

    public:
      /**
       * \brief Constructor
       *
       * \param[in] size The size of the created container.
       *
       * Creates a container with a given size.
       */
      explicit Container(Index size_in) :
        _foreign_memory(false)
      {
        _scalar_index.push_back(size_in);
      }

      /**
       * \brief Destructor
       *
       * Destroys a container and releases all of its used arrays.
       */
      virtual ~Container()
      {
        if(! _foreign_memory)
        {
          for (Index i(0) ; i < _elements.size() ; ++i)
            MemoryPool<Mem_>::release_memory(_elements.at(i));
          for (Index i(0) ; i < _indices.size() ; ++i)
            MemoryPool<Mem_>::release_memory(_indices.at(i));
        }
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
        _scalar_dt(std::move(other._scalar_dt)),
        _foreign_memory(other._foreign_memory)
      {
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
        for (Index i(0) ; i < _elements.size() ; ++i)
          MemoryPool<Mem_>::set_memory(_elements.at(i), value, _elements_size.at(i));
      }

      /**
       * \brief Free all allocated arrays
       *
       */
      virtual void clear()
      {
        if (! _foreign_memory)
        {
          for (Index i(0) ; i < _elements.size() ; ++i)
            MemoryPool<Mem_>::release_memory(this->_elements.at(i));
          for (Index i(0) ; i < _indices.size() ; ++i)
          MemoryPool<Mem_>::release_memory(this->_indices.at(i));
        }

        this->_elements.clear();
        this->_indices.clear();
        this->_elements_size.clear();
        this->_indices_size.clear();
        this->_scalar_index.clear();
        this->_scalar_dt.clear();
        this->_foreign_memory = false;
      }

      /** \brief Clone operation
       *
       * Become a copy of a given container.
       *
       * \param[in] other The source container.
       * \param[in] clone_mode The actual cloning procedure
       *
       */
      void clone(const Container & other, CloneMode clone_mode = CloneMode::Weak)
      {
        if (this == &other)
        {
          throw InternalError(__func__, __FILE__, __LINE__, "Trying to self-clone a lafem container!");
        }

        this->clear();

        this->_scalar_index.assign(other._scalar_index.begin(), other._scalar_index.end());
        this->_scalar_dt.assign(other._scalar_dt.begin(), other._scalar_dt.end());
        this->_elements_size.assign(other._elements_size.begin(), other._elements_size.end());
        this->_indices_size.assign(other._indices_size.begin(), other._indices_size.end());


        if (clone_mode == CloneMode::Deep)
        {
          for (Index i(0) ; i < other._indices.size() ; ++i)
          {
            this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(this->_indices_size.at(i)));
            MemoryPool<Mem_>::template copy<IT_>(this->_indices.at(i), other._indices.at(i), this->_indices_size.at(i));
          }

          for (Index i(0) ; i < other._elements.size() ; ++i)
          {
            this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(this->_elements_size.at(i)));
            MemoryPool<Mem_>::template copy<DT_>(this->_elements.at(i), other._elements.at(i), this->_elements_size.at(i));
          }

          return;
        }
        else
        {
          this->_indices.assign(other._indices.begin(), other._indices.end());
          for (Index i(0) ; i < this->_indices.size() ; ++i)
            MemoryPool<Mem_>::increase_memory(this->_indices.at(i));
        }

        if(clone_mode == CloneMode::Shallow)
        {
          this->_elements.assign(other._elements.begin(), other._elements.end());
          for (Index i(0) ; i < this->_elements.size() ; ++i)
            MemoryPool<Mem_>::increase_memory(this->_elements.at(i));

          return;
        }
        else if(clone_mode == CloneMode::Layout)
        {
          for (Index i(0) ; i < other._elements.size() ; ++i)
          {
            this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(this->_elements_size.at(i)));
          }

          return;
        }
        else if(clone_mode == CloneMode::Weak)
        {
          for (Index i(0) ; i < other._elements.size() ; ++i)
          {
            this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(this->_elements_size.at(i)));
            MemoryPool<Mem_>::template copy<DT_>(this->_elements.at(i), other._elements.at(i), this->_elements_size.at(i));
          }

          return;
        }
      }

      /// \copydoc clone(const Container&,CloneMode)
      template <typename Mem2_, typename DT2_, typename IT2_>
      void clone(const Container<Mem2_, DT2_, IT2_> & other, CloneMode clone_mode = CloneMode::Weak)
      {
        Container t(other.size());
        t.assign(other);
        clone(t, clone_mode);
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
        if (this == &other)
          return;

        for (Index i(0) ; i < this->_elements.size() ; ++i)
          MemoryPool<Mem_>::release_memory(this->_elements.at(i));
        for (Index i(0) ; i < this->_indices.size() ; ++i)
          MemoryPool<Mem_>::release_memory(this->_indices.at(i));

        this->_elements = std::move(other._elements);
        this->_indices = std::move(other._indices);
        this->_elements_size = std::move(other._elements_size);
        this->_indices_size = std::move(other._indices_size);
        this->_scalar_index = std::move(other._scalar_index);
        this->_scalar_dt = std::move(other._scalar_dt);

        this->_foreign_memory = other._foreign_memory;
      }

      /**
       * \brief Returns the total amount of bytes allocated.
       *
       * \returns The amount of bytes allocated in all arrays
       */
      std::size_t bytes() const
      {
        std::size_t tbytes(0);

        for (Index i(0) ; i < _elements_size.size() ; ++i)
        {
          tbytes += std::size_t(_elements_size.at(i) * sizeof(DT_));
        }

        for (Index i(0) ; i < _indices_size.size() ; ++i)
        {
          tbytes += std::size_t(_indices_size.at(i) * sizeof(IT_));
        }

        tbytes += std::size_t(_scalar_index.size() * sizeof(Index));
        tbytes += std::size_t(_scalar_dt.size() * sizeof(DT_));

        return tbytes;
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
       *
       * The template parameter of type Perspective is used in
       * blocked containers to switch between
       * the raw size and
       * the size, when every block is treated as one entry.
       *
       * Most containers like SparseMatrixCSR ignore this template parameter.
       */
      template <Perspective = Perspective::native>
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
       *
       * The template parameter of type Perspective is used in
       * block containers to switch between
       * the raw number of used elements and
       * the number of used elements, when every block is counted as one entry.
       *
       * Most containers like SparseMatrixCSR ignore this template parameter.
       */
      template <Perspective = Perspective::native>
      Index used_elements() const
      {
        return this->size();
      }

      /**
       * \brief Checks whether the container is empty.
       *
       * \returns \c true if the container is empty, otherwise \c false.
       */
      bool empty() const
      {
        return (this->size() == Index(0));
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
} // namespace FEAT

#endif // KERNEL_LAFEM_CONTAINER_HPP
