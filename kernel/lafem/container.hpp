// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_CONTAINER_HPP
#define KERNEL_LAFEM_CONTAINER_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/memory_pool.hpp>
#include <kernel/archs.hpp>
#include <kernel/lafem/base.hpp>
#include <kernel/util/type_traits.hpp>
#include <kernel/util/random.hpp>
#include <kernel/util/pack.hpp>


#include <vector>
#include <limits>
#include <cmath>
#include <typeinfo>
#include <string>
#include <type_traits>
#include <cstdlib>
#include <stdint.h>

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
          XABORTM("Should never be reached!");
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
     * \brief Config class for serialize parameter
     *
     * Data survey:
     * elements_compression LAFEM::SerializeMode, The compression mode for the elements array. Refer to CompressionModes for details.
     * indices_compression LAFEM::SerializeMode, The compression mode for the indices array
     * tolerance FEAT::Real, The max error that should be occur while compressing with zfp.
     * \note tolerance should be unset(meaning -1.) if elements_compression != CompressionModes::elements_zfp
     *          and strictly greater 0 if elements_compression == CompressionModes::elements_zfp.
     *
     * \note the corresponding configure flags 'zlib' and/or 'zfp' need to be added in the build-id at the configure call.
     *
     */
    class SerialConfig
    {
    private:
      CompressionModes elements_compression, indices_compression;
      FEAT::Real tolerance;
    public:
      /**
       * \brief Default Constructor
       *
       * Sets default values for compression, depending whether zlib is loaded or not
       */
#ifdef FEAT_HAVE_ZLIB
      explicit SerialConfig() :
      elements_compression(CompressionModes::elements_zlib),
      indices_compression(CompressionModes::indices_zlib),
      tolerance(FEAT::Real(-1.))
#else
      explicit SerialConfig() :
      elements_compression(CompressionModes::elements_off),
      indices_compression(CompressionModes::indices_off),
      tolerance(FEAT::Real(-1.))
#endif
      {
      }
      /**
       * \brief Constructor
       *
       * \param[in] zlib_compression bool, if true zlib will be used to compress indice arrays and if zfp is not loaded,
       *                                     will also be used to compress element arrays.
       * \param[in] zfp_compression bool, if true zfp will be used to compress element arrays.
       * \param[in] tolerance FEAT::Real, the maximum error done compressing with zfp. If zfp is enabled, has to be greater 0.
       *
       * \warning if zfp_compression is true and tolerance <= 0, an error will be raised.
       *
       * Set values depending on the configuration of the user
       */
      explicit SerialConfig(bool zlib_compression, bool zfp_compression, FEAT::Real tol = FEAT::Real(-1.))
      {
        set_elements_compression(CompressionModes::elements_off);
        set_indices_compression(CompressionModes::indices_off);
        if(zlib_compression)
        {
          set_elements_compression(CompressionModes::elements_zlib);
          set_indices_compression(CompressionModes::indices_zlib);
        }
        if(zfp_compression)
        {
          set_tolerance(tol);
          set_elements_compression(CompressionModes::elements_zfp);
        }
      }

      /**
       * \brief method for setting elements_compression
       *
       * \param[in] comp CompressionModes, describing the compression mode for _elements. See CompressionModes for details.
       *
       */
      void set_elements_compression(CompressionModes comp)
      {
        CompressionModes temp = comp & CompressionModes::elements_mask;
        switch(temp)
        {
          case CompressionModes::elements_off:
            break;
          case CompressionModes::elements_zlib:
#ifdef FEAT_HAVE_ZLIB
            break;
#else
            XABORTM("Zlib compression was enabled for elements, but zlib bibliography was not loaded. Please configure FEAT with -zlib flag!");
#endif
          case CompressionModes::elements_zfp:
#ifdef FEAT_HAVE_ZFP
            XASSERTM(tolerance > FEAT::Real(0.), "tolerance ist smaller or equall 0 while zfp is enabled!");
            break;
#else
            XABORTM("Zfp compression was enabled for elements, but zfp bibliography was not loaded. Please configure FEAT with -zfp flag!");
#endif
          default:
            XABORTM("No legal CompressionModes was given");
        }
        elements_compression = temp;
      }
      /**
       * \brief method for setting indices_compression
       *
       * \param[in] comp CompressionModes, describing the compression mode for _indices. See CompressionModes for details.
       *
       */
      void set_indices_compression(CompressionModes comp)
      {
        CompressionModes temp = comp & CompressionModes::indices_mask;
        switch(temp)
        {
          case CompressionModes::indices_off:
            break;
          case CompressionModes::indices_zlib:
#ifdef FEAT_HAVE_ZLIB
            break;
#else
            XABORTM("Zlib compression was enabled for indices, but zlib bibliography was not loaded. Please configure FEAT with -zlib flag!");
#endif
          default:
            XABORTM("No legal CompressionModes was given");
        }
        indices_compression = temp;
      }

      /**
       * \brief method for setting tolerance
       *
       * \param[in] tol FEAT::Real, the maximum error to be made by zfp compression. If zfp is enabled this has to be greater 0.
       *
       */
     void set_tolerance(FEAT::Real tol)
     {
       if(elements_compression == CompressionModes::elements_zfp)
       {
         XASSERTM(tol > FEAT::Real(0.), "Zfp enabled, but tolerance is smaller or equal to 0!");
       }
       tolerance = tol;
     }
     /**
       * \brief method for getting elements_compression
       *
       * \return CompressionModes, The compression mode for elements.
       *
       */
     CompressionModes get_elements_compression() const
     {
       return elements_compression;
     }

     /**
       * \brief method for getting indices_compression
       *
       * \return CompressionModes, The compression mode for indices.
       *
       */
     CompressionModes get_indices_compression() const
     {
       return indices_compression;
     }

     /**
       * \brief method for getting tolerance
       *
       * \return FEAT::Real, the maximum error to be made by zfp compression. If zfp is enabled this has to be greater 0.
       *
       */
     FEAT::Real get_tolerance() const
     {
       return tolerance;
     }

    }; //class SerialConfig


    /**
     * \brief Container base class.
     *
     * All LAFEM-Container are derived from this base class. \sa \ref lafem_design
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

        XASSERTM(_elements.size() == other.get_elements().size(), "Container size mismatch!");
        XASSERTM(_indices.size() == other.get_indices().size(), "Container size mismatch!");
        XASSERTM(_scalar_index.size() == other.get_scalar_index().size(), "Container size mismatch!");
        XASSERTM(_scalar_dt.size() == other.get_scalar_dt().size(), "Container size mismatch!");

        if (full)
        {
          this->_scalar_index.assign(other._scalar_index.begin(), other._scalar_index.end());
          this->_scalar_dt.assign(other._scalar_dt.begin(), other._scalar_dt.end());

          for (Index i(0) ; i < _indices.size() ; ++i)
          {
            XASSERTM(_indices_size.at(i) == other.get_indices_size().at(i), "Container size mismatch!");
            MemoryPool<Mem_>::template copy<IT_>(_indices.at(i), other.get_indices().at(i), _indices_size.at(i));
          }
        }

        for (Index i(0) ; i < _elements.size() ; ++i)
        {
          XASSERTM(_elements_size.at(i) == other.get_elements_size().at(i), "Container size mismatch!");
          MemoryPool<Mem_>::template copy<DT_>(_elements.at(i), other.get_elements().at(i), _elements_size.at(i));
        }

      }

      template <typename Mem2_>
      void _copy_content(const Container<Mem2_, DT_, IT_> & other, bool full)
      {
        XASSERTM(_elements.size() == other.get_elements().size(), "Container size mismatch!");
        XASSERTM(_indices.size() == other.get_indices().size(), "Container size mismatch!");
        XASSERTM(_scalar_index.size() == other.get_scalar_index().size(), "Container size mismatch!");
        XASSERTM(_scalar_dt.size() == other.get_scalar_dt().size(), "Container size mismatch!");

        if (full)
        {
          this->_scalar_index.assign(other.get_scalar_index().begin(), other.get_scalar_index().end());
          this->_scalar_dt.assign(other.get_scalar_dt().begin(), other.get_scalar_dt().end());

          for (Index i(0) ; i < _indices.size() ; ++i)
          {
            XASSERTM(_indices_size.at(i) == other.get_indices_size().at(i), "Container size mismatch!");
            if (std::is_same<Mem_, Mem::Main>::value && std::is_same<Mem2_, Mem::CUDA>::value)
              MemoryPool<Mem2_>::template download<IT_>(_indices.at(i), other.get_indices().at(i), _indices_size.at(i));
            else if (std::is_same<Mem_, Mem::CUDA>::value && std::is_same<Mem2_, Mem::Main>::value)
              MemoryPool<Mem_>::template upload<IT_>(_indices.at(i), other.get_indices().at(i), _indices_size.at(i));
            else
              XABORTM("Memory Backend not known!");
          }
        }

        for (Index i(0) ; i < _elements.size() ; ++i)
        {
          XASSERTM(_elements_size.at(i) == other.get_elements_size().at(i), "Container size mismatch!");
          if (std::is_same<Mem_, Mem::Main>::value && std::is_same<Mem2_, Mem::CUDA>::value)
            MemoryPool<Mem2_>::template download<DT_>(_elements.at(i), other.get_elements().at(i), _elements_size.at(i));
          else if (std::is_same<Mem_, Mem::CUDA>::value && std::is_same<Mem2_, Mem::Main>::value)
            MemoryPool<Mem_>::template upload<DT_>(_elements.at(i), other.get_elements().at(i), _elements_size.at(i));
          else
            XABORTM("Memory Backend not known!");
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
       * \brief Calculation of the serialized size with optional compression of complete container entity.
       *
       * \param[in] config LAFEM::SerialConfig, a struct describing the serialize configuration.
       *
       * \return An unsigned int containing the size of serialized container.
       *
       * Calculate the size of the serialized container entity.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::uint64_t _serialized_size(const LAFEM::SerialConfig& config = LAFEM::SerialConfig()) const
      {
        Container<Mem::Main, DT2_, IT2_> tc(0);
        tc.assign(*this);

        std::uint64_t gsize(4 * sizeof(std::uint64_t)); //raw array size + magic number + type_index DT_ + type_index IT_
        gsize += 7 * sizeof(std::uint64_t); // size of all seven stl containers
        gsize += 2 * tc._elements_size.size() * sizeof(std::uint64_t); // _elements_size contents + _elements_arrays bytesize
        gsize += 2 * tc._indices_size.size() * sizeof(std::uint64_t); // _indices_size contents + _indizes_arrays bytesize
        gsize += tc._scalar_index.size() * sizeof(std::uint64_t); // _scalar_index contents
        gsize += tc._scalar_dt.size() * sizeof(DT2_); // _scalar_dt contents

        CompressionModes compress = config.get_elements_compression() | config.get_indices_compression();
        Pack::Type compression_type_elements = Pack::deduct_type<DT2_>();
        Pack::Type compression_type_index = Pack::deduct_type<IT2_>();
        if((compress & CompressionModes::elements_mask) == CompressionModes::elements_zfp) // If zfp is used for elements
        {
          compression_type_elements = compression_type_elements | Pack::Type::Mask_P;
        }
        else if((compress & CompressionModes::elements_mask) == CompressionModes::elements_zlib) //If zlib is used and zfp is not used elements
        {
          compression_type_elements = compression_type_elements | Pack::Type::Mask_Z;
        }
        if((compress & CompressionModes::indices_mask) == CompressionModes::indices_zlib) // If zlib is used for indices
            compression_type_index = compression_type_index | Pack::Type::Mask_Z;


        FEAT::Real tolerance = config.get_tolerance();

        for (Index i(0) ; i < tc._elements_size.size() ; ++i)
        {
          gsize += Pack::estimate_size(tc._elements_size.at(i), compression_type_elements, (double)tolerance); // upper_bound for (compressed) elements
        }

        for (Index i(0) ; i < tc._indices_size.size() ; ++i)
        {
          gsize += Pack::estimate_size(tc._indices_size.at(i), compression_type_index); // upper_bound for (compressed) indizes
        }
        gsize += 16; //padding for datatype alignment mismatch

        return gsize;
      }

      /**
       *
       * \brief Serialization of complete container entity(with options of lossless and lossy compression of data arrays).
       *
       * \param[in] mode FileMode enum, describing the actual container specialization.
       * \param[in] config LAFEM::SerialConfig, a struct describing the serialize configuration.
       * \note the corresponding configure flags 'zlib' and/or 'zfp' need to be added in the build-id at the configure call.
       *
       * \returns A std::vector, containing the byte array.
       *
       * Serialize a complete container entity into a single binary array.
       *
       * Array data layout:
       * \code
       * raw array size in bytes (std::uint64_t)
       * magic number (std::uint64_t)
       * DT_ 'hash' value (sizeof(DT_) and floating/integral distinction (packed std::uint64_t)
       * IT_ 'hash' value (sizeof(IT_) and floating/integral distinction (packed std::uint64_t)
       * _elements.size() (std::uint64_t)
       * _indices.size() (std::uint64_t)
       * _elements_size.size() (std::uint64_t)
       * _indices_size.size() (std::uint64_t)
       * _scalar_index.size() (std::uint64_t)
       * _scalar_dt.size() (std::uint64_t)
       * _compressed (std::uint64_t)
       *
       * _elements_size array (std::uint64_t)
       * _elements arrays byte_size (std::uint64_t) just temporary to read in data
       * _indices_size array (std::uint64_t)
       * _indices arrays byte_size (std::uint64_t)
       * _scalar_index array (std::uint64_t)
       * _scalar_dt array (DT2_)
       * _elements arrays (DT2_)
       * _indices arrays (IT2_)
       * \endcode
       *
       * See \ref FEAT::LAFEM::SerialConfig for details.
       */

      template <typename DT2_ = DT_, typename IT2_ = IT_>
      std::vector<char> _serialize(FileMode mode, const SerialConfig& config = SerialConfig()) const
      {
        std::uint64_t raw_size = 0u;
        FEAT::Real tolerance = config.get_tolerance();

        CompressionModes compress = config.get_elements_compression() | config.get_indices_compression();
        Pack::Type compression_type_elements = Pack::deduct_type<DT2_>();
        Pack::Type compression_type_index = Pack::deduct_type<IT2_>();
        if((compress & CompressionModes::elements_mask) == CompressionModes::elements_zfp) // If zfp is used for elements
        {
          compression_type_elements = compression_type_elements | Pack::Type::Mask_P;
        }
        else if((compress & CompressionModes::elements_mask) == CompressionModes::elements_zlib) //If zlib is used and zfp is not used elements
        {
          compression_type_elements = compression_type_elements | Pack::Type::Mask_Z;
        }
        if((compress & CompressionModes::indices_mask) == CompressionModes::indices_zlib) // If zlib is used for indices
            compression_type_index = compression_type_index | Pack::Type::Mask_Z;

        Container<Mem::Main, DT2_, IT2_> tc(0);
        tc.assign(*this);

        std::uint64_t gsize = this->template _serialized_size<DT2_, IT2_>(config);

        std::vector<char> result((size_t(gsize)));
        char * array(result.data());
        std::uint64_t * uiarray(reinterpret_cast<std::uint64_t *>(array));
        DT2_ * dtarray(reinterpret_cast<DT2_ *>(array));
        IT2_ * itarray(reinterpret_cast<IT2_ *>(array));

        /// \compilerhack clang seems to use older libc++ without std::underlying_type
        #if defined(FEAT_COMPILER_CLANG)
        std::uint64_t magic = (std::uint64_t)static_cast<__underlying_type(FileMode)>(mode);
        #else
        std::uint64_t magic = (std::uint64_t)static_cast<typename std::underlying_type<FileMode>::type>(mode);
        #endif
        uiarray[0] = gsize; //we will overwrite this at the end with the finalized data...
        uiarray[1] = magic;
        uiarray[2] = Type::Traits<DT_>::feature_hash();
        uiarray[3] = Type::Traits<IT_>::feature_hash();
        uiarray[4] = tc._elements.size();
        uiarray[5] = tc._indices.size();
        uiarray[6] = tc._elements_size.size();
        uiarray[7] = tc._indices_size.size();
        uiarray[8] = tc._scalar_index.size();
        uiarray[9] = tc._scalar_dt.size();
        uiarray[10] = (std::uint64_t)compress;
        raw_size += 11 * (std::uint64_t)sizeof(std::uint64_t);

        Index global_i(11); // count how many elements of std::uint64_t have been inserted so far
        Index local_elements(0); // count where elements array bytesize starts
        Index local_index(0); // same for index...

        for (Index i(0) ; i < tc._elements_size.size() ; ++i)
        {
          uiarray[i + global_i] = tc._elements_size.at(i);
        }
        global_i += Index(tc._elements_size.size());
        local_elements = global_i;
        for (Index i(0) ; i < tc._elements_size.size() ; ++i)
        {
          uiarray[i + global_i] = (std::uint64_t)0u; //placeholder, until we calculate the real data...
        }
        global_i += Index(tc._elements_size.size());
        raw_size += 2 * std::uint64_t(tc._elements_size.size()) * std::uint64_t(sizeof(std::uint64_t));

        for (Index i(0) ; i < tc._indices_size.size() ; ++i)
        {
          uiarray[i + global_i] = tc._indices_size.at(i);
        }
        global_i += Index(tc._indices_size.size());
        local_index = global_i;
        for (Index i(0) ; i < tc._indices_size.size() ; ++i)
        {
          uiarray[i + global_i] = 0;
        }
        global_i += Index(tc._indices_size.size());
        raw_size += 2 * std::uint64_t(tc._indices_size.size()) * std::uint64_t(sizeof(std::uint64_t));

        for (Index i(0) ; i < tc._scalar_index.size() ; ++i)
        {
          uiarray[i + global_i] = tc._scalar_index.at(i);
        }
        global_i += Index(tc._scalar_index.size());
        raw_size += std::uint64_t(tc._scalar_index.size()) * std::uint64_t(sizeof(std::uint64_t));

        global_i = Index((global_i * sizeof(std::uint64_t) + sizeof(DT2_) - 1u) / sizeof(DT2_)); // now counting how many DT2 have been inserted so far

        for (Index i(0) ; i < tc._scalar_dt.size() ; ++i)
        {
          dtarray[i + global_i] = tc._scalar_dt.at(i);
        }
        global_i += Index(tc._scalar_dt.size());
        raw_size += std::uint64_t(tc._scalar_dt.size()) * std::uint64_t(sizeof(DT2_));

        if((compress & CompressionModes::elements_mask) != CompressionModes::elements_off)
        {
          global_i = Index((global_i * sizeof(DT2_) + sizeof(char) - 1u) / sizeof(char)); // now counting how many bytes/chars have been inserted so far
          for (Index i(0) ; i < tc._elements.size() ; ++i)
          {
            char* dest = &result[global_i];
            std::size_t est_buff_size = Pack::estimate_size(tc._elements_size.at(i), compression_type_elements, (double)tolerance);
            std::size_t real_size = Pack::encode<DT2_>(dest, tc._elements.at(i), est_buff_size, tc._elements_size.at(i), compression_type_elements, false, (double)tolerance);
            XASSERTM(real_size != 0u, "could not encode elements");
            uiarray[i + local_elements] = (std::uint64_t)real_size; //save used memory for compressed array
            global_i += (Index)real_size; //go forward the number of bytes(chars)
            raw_size += (std::uint64_t)real_size * (std::uint64_t)sizeof(char);
          }
          global_i = Index((global_i * sizeof(char) + sizeof(DT2_) - 1u) / sizeof(DT2_)); // go back to DT2_ indexing
        }
        else
        {
          for (Index i(0) ; i < tc._elements.size() ; ++i)
          {
            std::memcpy(&dtarray[global_i], tc._elements.at(i), tc._elements_size.at(i) * sizeof(DT2_));
            uiarray[i + local_elements] = (std::uint64_t) tc._elements_size.at(i) * sizeof(DT2_);
            global_i += tc._elements_size.at(i);
            raw_size += (std::uint64_t) tc._elements_size.at(i) * sizeof(DT2_);
          }
        }
        if((compress & CompressionModes::indices_mask) != CompressionModes::indices_off)
        {
          global_i = Index((global_i * sizeof(DT2_) + sizeof(char) - 1u) / sizeof(char)); // now counting how many bytes/chars have been inserted so far
          for (Index i(0) ; i < tc._indices.size() ; ++i)
          {
            char* dest = &result[global_i];
            std::size_t est_buff_size = Pack::estimate_size(tc._indices_size.at(i), compression_type_index);
            std::size_t real_size = Pack::encode<IT2_>(dest, tc._indices.at(i), est_buff_size, tc._indices_size.at(i), compression_type_index, false);
            XASSERTM(real_size != 0u, "could not encode indices");
            uiarray[i + local_index] = (std::uint64_t) real_size;
            global_i += Index(real_size);
            raw_size += std::uint64_t(real_size) * std::uint64_t(sizeof(char));
          }
        }
        else
        {
          global_i = Index((global_i * sizeof(DT2_) + sizeof(IT2_) - 1u) / sizeof(IT2_)); // now counting IT2_ elements
          for (Index i(0) ; i < tc._indices.size() ; ++i)
          {
            std::memcpy(&itarray[global_i], tc._indices.at(i), tc._indices_size.at(i) * sizeof(IT2_));
            uiarray[i + local_index] = (std::uint64_t) tc._indices_size.at(i) * sizeof(IT2_);
            global_i += tc._indices_size.at(i);
            raw_size += (std::uint64_t) tc._indices_size.at(i) * sizeof(IT2_);
          }
        }
        uiarray[0] = raw_size + 16u; //padding
        //std::cout << "Compressed size is: " << uiarray[0] << std::endl;

        result.resize(std::size_t(uiarray[0]));
        return result;
      }

      /**
       *
       * \brief Serialization of complete container entity(with options of lossless and lossy compression of data arrays).
       *
       * \param[in] mode FileMode enum, describing the actual container specialization.
       * \param[in] file The output stream to write data into.
       * \param[in] config LAFEM::SerialConfig, a struct describing the serialize configuration.
       * \note the corresponding configure flags 'zlib' and/or 'zfp' need to be added in the build-id at the configure call.
       *
       * Serialize a complete container entity into a single binary file with compression enabled.
       * See \ref FEAT::LAFEM::SerialConfig for details.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void _serialize(FileMode mode, std::ostream & file, const SerialConfig& config = SerialConfig()) const
      {
        auto temp(this->template _serialize<DT2_, IT2_>(mode, config));
        file.write(temp.data(), long(temp.size()));
        if (!file.good())
          XABORTM("Error in _serialize - file ostream is not good anymore!");
      }

       /**
       * \brief Deserialization of complete container entity with possible decompression.
       *
       * \param[in] mode FileMode enum, describing the actual container specialization.
       * \param[in] input A std::vector containing the byte array.
       *
       * Recreate a complete container entity by a single binary array.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void _deserialize(FileMode mode, std::vector<char> & input)
      {
        this->clear();
        Container<Mem::Main, DT2_, IT2_> tc(0);
        tc.clear();

        char * array(input.data());
        std::uint64_t * uiarray(reinterpret_cast<std::uint64_t *>(array));
        DT2_ * dtarray(reinterpret_cast<DT2_ *>(array));
        IT2_ * itarray(reinterpret_cast<IT2_ *>(array));

        /// \compilerhack clang seems to use older libc++ without std::underlying_type
#if defined(FEAT_COMPILER_CLANG)
        std::uint64_t magic = (std::uint64_t)static_cast<__underlying_type(FileMode)>(mode);
#else
        std::uint64_t magic = (std::uint64_t)static_cast<typename std::underlying_type<FileMode>::type>(mode);
#endif
        XASSERTM(magic == uiarray[1], "_deserialize: given FileMode incompatible with given array!");

        //ensure that we have the same integral/floating type configuration, that was used when storing the serialized data
        XASSERT(Type::Traits<DT_>::is_int == Type::Helper::extract_intness(uiarray[2]));
        XASSERT(Type::Traits<DT_>::is_float == Type::Helper::extract_floatness(uiarray[2]));
        XASSERT(Type::Traits<DT_>::is_signed == Type::Helper::extract_signedness(uiarray[2]));
        XASSERT(Type::Traits<IT_>::is_int == Type::Helper::extract_intness(uiarray[3]));
        XASSERT(Type::Traits<IT_>::is_float == Type::Helper::extract_floatness(uiarray[3]));
        XASSERT(Type::Traits<IT_>::is_signed == Type::Helper::extract_signedness(uiarray[3]));

        if (sizeof(DT_) > Type::Helper::extract_type_size(uiarray[2]))
          std::cerr<<"Warning: You are reading a container floating point in higher precision then it was saved before!"<<std::endl;

        if (sizeof(IT_) > Type::Helper::extract_type_size(uiarray[3]))
          std::cerr<<"Warning: You are reading a container integral type in higher precision then it was saved before!"<<std::endl;

        CompressionModes compress = (CompressionModes)uiarray[10];
        //test whether system has right third_party software loaded:
#ifndef FEAT_HAVE_ZLIB
        XASSERTM((compress & CompressionModes::elements_mask) != CompressionModes::elements_zlib,
                 "Data was compressed with ZLIB! To read in data, please configure FEAT with zlib enabled.\n For more information see FEAT documentation.");
        XASSERTM((compress & CompressionModes::indices_mask) != CompressionModes::indices_zlib,
                 "Data was compressed with ZLIB! To read in data, please configure FEAT with zlib enabled.\n For more information see FEAT documentation.");
#endif
#ifndef FEAT_HAVE_ZFP
        XASSERTM((compress & CompressionModes::elements_mask) != CompressionModes::elements_zfp,
                 "Data was compressed with zfp! To read in data, please configure FEAT with zfp enabled.\n For more information see FEAT documentation.");
#endif
        std::vector<std::uint64_t> elements_bytes(uiarray[6]); //temp arrays to keep track of bytessize of compressed arrays
        std::vector<std::uint64_t> indices_bytes(uiarray[7]);
        Index global_i(11);
        for (std::uint64_t i(0) ; i < uiarray[6] ; ++i)
        {
          tc._elements_size.push_back(Index(uiarray[i + global_i]));
        }
        global_i += Index(uiarray[6]);
        for (std::uint64_t i(0) ; i < uiarray[6] ; ++i)
        {
          elements_bytes[i] = uiarray[i + global_i];
        }
        global_i += Index(uiarray[6]);

        for (std::uint64_t i(0) ; i < uiarray[7] ; ++i)
        {
          tc._indices_size.push_back(Index(uiarray[i + global_i]));
        }
        global_i += Index(uiarray[7]);
        for (std::uint64_t i(0) ; i < uiarray[7] ; ++i)
        {
          indices_bytes[i] = uiarray[i + global_i];
        }
        global_i += Index(uiarray[7]);

        for (std::uint64_t i(0) ; i < uiarray[8] ; ++i)
        {
          tc._scalar_index.push_back(Index(uiarray[i + global_i]));
        }
        global_i += Index(uiarray[8]);

        global_i = Index((global_i * sizeof(std::uint64_t) + sizeof(DT2_) - 1u) / sizeof(DT2_));
        for (std::uint64_t i(0) ; i < uiarray[9] ; ++i)
        {
          tc._scalar_dt.push_back(dtarray[i + global_i]);
        }
        global_i += Index(uiarray[9]);

        Pack::Type compression_type_elements = Pack::deduct_type<DT2_>();
        Pack::Type compression_type_index = Pack::deduct_type<IT2_>();
        if((compress & CompressionModes::elements_mask) == CompressionModes::elements_zfp)// If zfp is used
        {
          compression_type_elements = compression_type_elements | Pack::Type::Mask_P;
        }
        else if((compress & CompressionModes::elements_mask) == CompressionModes::elements_zlib) //If zlib is loaded and zfp is not used
        {
          compression_type_elements = compression_type_elements | Pack::Type::Mask_Z;
        }
        if((compress & CompressionModes::indices_mask) == CompressionModes::indices_zlib)
        {
          compression_type_index = compression_type_index | Pack::Type::Mask_Z;
        }

        if((compress & CompressionModes::elements_mask) != CompressionModes::elements_off)
        {
          global_i = Index((global_i * sizeof(DT2_) + sizeof(char) - 1u) / sizeof(char));
          for(Index i(0); i < Index(uiarray[4]); ++i)
          {
            tc._elements.push_back(MemoryPool<Mem::Main>::template allocate_memory<DT2_>(tc._elements_size.at(i)));
            if(0u == Pack::decode(tc._elements.at(i), &array[global_i], (std::size_t) tc._elements_size.at(i), (std::size_t) elements_bytes[i], compression_type_elements, false))
              XABORTM("Cannot decode compressed elements array");
            global_i += (Index)elements_bytes[i];
          }
          global_i = Index((global_i * sizeof(char) + sizeof(DT2_) - 1u) / sizeof(DT2_));
        }
        else
        {
          for (Index i(0) ; i < Index(uiarray[4]) ; ++i)
          {
            tc._elements.push_back(MemoryPool<Mem::Main>::template allocate_memory<DT2_>(tc._elements_size.at(i)));
            MemoryPool<Mem::Main>::template upload<DT2_>(tc._elements.at(i), &dtarray[global_i], tc._elements_size.at(i));
            global_i += tc._elements_size.at(i);
          }
        }

        if((compress & CompressionModes::indices_mask) != CompressionModes::indices_off)
        {
          global_i = Index((global_i * sizeof(DT2_) + sizeof(char) - 1u) / sizeof(char));
          for (Index i(0) ; i < Index(uiarray[5]) ; ++i)
          {
            tc._indices.push_back(MemoryPool<Mem::Main>::template allocate_memory<IT2_>(tc._indices_size.at(i)));
            if(0u == Pack::decode(tc._indices.at(i), &array[global_i], (std::size_t) tc._indices_size.at(i), (std::size_t) indices_bytes[i], compression_type_index, false))
              XABORTM("Cannot decode compressed indice array");
            global_i += (Index)indices_bytes[i];
          }
        }
        else
        {
          global_i = Index((global_i * sizeof(DT2_) + sizeof(IT2_) - 1u) / sizeof(IT2_));
          for (Index i(0) ; i < Index(uiarray[5]) ; ++i)
          {
            tc._indices.push_back(MemoryPool<Mem::Main>::template allocate_memory<IT2_>(tc._indices_size.at(i)));
            MemoryPool<Mem::Main>::template upload<IT2_>(tc._indices.at(i), &itarray[global_i], tc._indices_size.at(i));
            global_i += tc._indices_size.at(i);
          }
        }

        this->assign(tc);
      }

      /**
       * \brief Deserialization of complete container entity with possible decompression.
       *
       * \param[in] mode FileMode enum, describing the actual container specialization.
       * \param[in] file std::istream, pointing to the input data.
       *
       * Recreate a complete container entity by a single binary file.
       */
      template <typename DT2_ = DT_, typename IT2_ = IT_>
      void _deserialize(FileMode mode, std::istream & file)
      {
        std::uint64_t tsize;
        file.read((char *)&tsize, (long)(sizeof(std::uint64_t)));
        std::vector<char> temp((size_t(tsize)));
        file.seekg(-(long)sizeof(std::uint64_t), file.cur);
        file.read(temp.data(), (long)(tsize));
        if (!file.good())
          XABORTM("Error in _deserialize - file istream is not good anymore!");
        this->template _deserialize<DT2_, IT2_>(mode, temp);
      }

    public:
      /**
       * \brief Constructor
       *
       * \param[in] size_in The size of the created container.
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
       */
      void format(DT_ value = DT_(0))
      {
        for (Index i(0) ; i < _elements.size() ; ++i)
          MemoryPool<Mem_>::set_memory(_elements.at(i), value, _elements_size.at(i));
      }

      /**
       * \brief Reset all elements of the container to random values.
       *
       * \param[in] rng The random number generator.
       * \param[in] min Lower rng bound.
       * \param[in] max Upper rng bound.
       */
      void format(Random & rng, DT_ min, DT_ max)
      {
        for (Index e(0) ; e < _elements.size() ; ++e)
        {
          for (Index i(0) ; i < _elements_size.at(e); ++i)
          {
            MemoryPool<Mem_>::set_memory(this->_elements.at(e) + i, rng(min, max));
          }
        }
      }

      /**
       * \brief Free all allocated arrays
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
          XABORTM("Trying to self-clone a lafem container!");
        }

        XASSERTM(other._foreign_memory == false || clone_mode == CloneMode::Deep, "Must use deep cloning with ranged based source containers");

        this->clear();

        this->_scalar_index.assign(other._scalar_index.begin(), other._scalar_index.end());
        this->_scalar_dt.assign(other._scalar_dt.begin(), other._scalar_dt.end());
        this->_elements_size.assign(other._elements_size.begin(), other._elements_size.end());
        this->_indices_size.assign(other._indices_size.begin(), other._indices_size.end());


        if (clone_mode == CloneMode::Deep || clone_mode == CloneMode::Allocate)
        {
          for (Index i(0) ; i < other._indices.size() ; ++i)
          {
            this->_indices.push_back(MemoryPool<Mem_>::template allocate_memory<IT_>(this->_indices_size.at(i)));
            if (clone_mode == CloneMode::Deep)
              MemoryPool<Mem_>::template copy<IT_>(this->_indices.at(i), other._indices.at(i), this->_indices_size.at(i));
          }

          for (Index i(0) ; i < other._elements.size() ; ++i)
          {
            this->_elements.push_back(MemoryPool<Mem_>::template allocate_memory<DT_>(this->_elements_size.at(i)));
            if (clone_mode == CloneMode::Deep)
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

      /// \copydoc FEAT::Control::Checkpointable::get_checkpoint_size()
      std::uint64_t get_checkpoint_size(LAFEM::SerialConfig& config)
      {
        return this->template _serialized_size<>(config);
      }

      /// \copydoc FEAT::Control::Checkpointable::restore_from_checkpoint_data(std::vector<char>&)
      void restore_from_checkpoint_data(std::vector<char> & data)
      {
        this->template _deserialize<>(FileMode::fm_binary, data);
      }

      /// \copydoc FEAT::Control::Checkpointable::set_checkpoint_data(std::vector<char>&)
      std::uint64_t set_checkpoint_data(std::vector<char>& data, LAFEM::SerialConfig& config)
      {
        auto buffer = this->template _serialize<>(FileMode::fm_binary, config);
        data.insert(std::end(data), std::begin(buffer), std::end(buffer));
        return std::uint64_t(buffer.size());
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
