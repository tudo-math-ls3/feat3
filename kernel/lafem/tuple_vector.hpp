#pragma once
#ifndef KERNEL_LAFEM_TUPLE_VECTOR_HPP
#define KERNEL_LAFEM_TUPLE_VECTOR_HPP 1

// includes, FEAT
#include <kernel/lafem/meta_element.hpp>

// includes, system
#include <iostream>
#include <type_traits>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Variadic TupleVector class template
     *
     * This class template implements a composition of sub-vectors of arbitrary classes.
     *
     * \tparam First_, ...Rest_
     * A sequence of (meta) vector classes which are to be composed.
     *
     * \author Peter Zajac
     */
    template<
      typename First_,
      typename... Rest_>
    class TupleVector
    {
    private:
      template<typename,typename...>
      friend class TupleVector;

      typedef TupleVector<Rest_...> RestClass;

      /// read binary data of _first and _rest to file.
      void _read_from_binary(std::istream& file)
      {
        //process first
        _first.read_from(FileMode::fm_binary, file);

        // append rest
        _rest._read_from_binary(file);
      }

      /// write binary data of _first and _rest to file.
      void _write_out_binary(std::ostream& file) const
      {
        //process first
        _first.write_out(FileMode::fm_binary, file);

        // append rest
        _rest._write_out_binary(file);
      }

    public:
      /// number of vector blocks
      static constexpr int num_blocks = TupleVector<Rest_...>::num_blocks + 1;

      /// sub-vector mem-type
      typedef typename First_::MemType MemType;
      /// sub-vector data-type
      typedef typename First_::DataType DataType;
      /// sub-vector index-type
      typedef typename First_::IndexType IndexType;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = TupleVector<
        typename First_::template ContainerType<Mem2_, DT2_, IT2_>,
        typename Rest_::template ContainerType<Mem2_, DT2_, IT2_>...>;

      /// this typedef lets you create a vector container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

      // ensure that all sub-vector have the same mem- and data-type
      static_assert(std::is_same<MemType, typename RestClass::MemType>::value,
                    "sub-vectors have different mem-types");
      static_assert(std::is_same<DataType, typename RestClass::DataType>::value,
                    "sub-vectors have different data-types");
      static_assert(std::is_same<IndexType, typename RestClass::IndexType>::value,
                    "sub-vectors have different index-types");

    protected:
      /// the first sub-vector
      First_ _first;
      /// all remaining sub-vectors
      RestClass _rest;

      /// Returns a list of all sub-vector type names
      static String sub_name_list()
      {
        return First_::name() + "," + RestClass::sub_name_list();
      }

    public:
      /// default CTOR
      TupleVector()
      {
      }

      /// rest-class emplacement ctor; for internal use only
      explicit TupleVector(First_&& the_first, RestClass&& the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest))
      {
      }

      /// Sub-Vector emplacement constructor
      explicit TupleVector(First_&& the_first, Rest_&&... the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest...))
      {
      }

      /// move ctor
      TupleVector(TupleVector&& other) :
        _first(std::move(other._first)),
        _rest(std::move(other._rest))
      {
      }

      /// move-assign operator
      TupleVector& operator=(TupleVector&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
          _rest = std::move(other._rest);
        }
        return *this;
      }

      /// deleted copy-ctor
      TupleVector(const TupleVector&) = delete;
      /// deleted copy-assign operator
      TupleVector& operator=(const TupleVector&) = delete;

      /**
       * \brief Creates and returns a copy of this vector
       *
       * \param[in] mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       *
       */
      TupleVector clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return TupleVector(_first.clone(mode), _rest.clone(mode));
      }

      /**
       * \brief Turns this vector into a clone of other
       *
       * \param[in] clone_mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       *
       */
      void clone(const TupleVector& other, CloneMode clone_mode)
      {
        _first.clone(other._first, clone_mode);
        _rest.clone(other._rest, clone_mode);
      }

      /**
       * \brief Turns this vector into a clone of other
       *
       * As the default CloneMode of the underlying containers is unknown, this has to be separate.
       *
       */
      void clone(const TupleVector& other)
      {
        _first.clone(other._first);
        _rest.clone(other._rest);
      }

      /// \cond internal
      First_& first()
      {
        return _first;
      }

      const First_& first() const
      {
        return _first;
      }

      TupleVector<Rest_...>& rest()
      {
        return _rest;
      }

      const TupleVector<Rest_...>& rest() const
      {
        return _rest;
      }
      /// \endcond

      /**
       * \brief Returns a sub-vector block.
       *
       * \tparam i_
       * The index of the sub-vector block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-vector at position \p i_.
       */
      template<int i_>
      typename TupleElement<i_, First_, Rest_...>::Type& at()
      {
        static_assert((0 <= i_) && (i_ < num_blocks), "invalid sub-vector index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc at() */
      template<int i_>
      typename TupleElement<i_, First_, Rest_...>::Type const& at() const
      {
        static_assert((0 <= i_) && (i_ < num_blocks), "invalid sub-vector index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /// Returns the total size of this tuple-vector.
      template <Perspective perspective_ = Perspective::native>
      Index size() const
      {
        return _first.template size<perspective_>() + _rest.template size<perspective_>();
      }

      /// Returns the number of blocks in this tuple-vector.
      int blocks() const
      {
        return num_blocks;
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("TupleVector<") + sub_name_list() + ">";
      }

      /**
       * \brief Reset all elements of the container to a given value or zero if missing.
       *
       * \param[in] value
       * The value to which the vector's entries are to be set to.
       */
      void format(DataType value = DataType(0))
      {
        first().format(value);
        rest().format(value);
      }

      /// Free all allocated arrays
      void clear()
      {
        first().clear();
        rest().clear();
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The vector to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      //template<typename First2_, typename... Rest2_>
      void copy(const TupleVector/*<First2_, Rest2_...>*/& x, bool full = false)
      {
        first().copy(x.first(), full);
        rest().copy(x.rest(), full);
      }

      void axpy(const TupleVector& x, const TupleVector& y, DataType alpha = DataType(1))
      {
        first().axpy(x.first(), y.first(), alpha);
        rest().axpy(x.rest(), y.rest(), alpha);
      }

      void component_product(const TupleVector & x, const TupleVector & y)
      {
        first().component_product(x.first(), y.first());
        rest().component_product(x.rest(), y.rest());
      }

      void component_invert(const TupleVector& x, DataType alpha = DataType(1))
      {
        first().component_invert(x.first(), alpha);
        rest().component_invert(x.rest(), alpha);
      }

      void scale(const TupleVector& x, DataType alpha)
      {
        first().scale(x.first(), alpha);
        rest().scale(x.rest(), alpha);
      }

      DataType dot(const TupleVector& x) const
      {
        return first().dot(x.first()) + rest().dot(x.rest());
      }

      /**
       * \copydoc LAFEM::DenseVector::triple_dot()
       **/
      DataType triple_dot(const TupleVector& x, const TupleVector& y) const
      {
        return first().triple_dot(x.first(), y.first())
          + rest().triple_dot(x.rest(), y.rest());
      }

      DataType norm2sqr() const
      {
        return first().norm2sqr() + rest().norm2sqr();
      }

      DataType norm2() const
      {
        return Math::sqrt(norm2sqr());
      }

      /**
       * \brief Retrieve the absolute maximum value of this vector.
       *
       * \return The largest absolute value.
       */
      DataType max_element() const
      {
        return Math::max(first().max_element(), rest().max_element());
      }

      /// \cond internal
      /// Writes the vector-entries in an allocated array
      void set_vec(DataType * const pval_set) const
      {
        this->first().set_vec(pval_set);
        this->rest().set_vec(pval_set + this->first().template size<Perspective::pod>());
      }

      /// Writes data of an array in the vector
      void set_vec_inv(const DataType * const pval_set)
      {
        this->first().set_vec_inv(pval_set);
        this->rest().set_vec_inv(pval_set + this->first().template size<Perspective::pod>());
      }
      /// \endcond

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Vector.
       *
       * Use source vector content as content of current vector
       */
      template <typename First2_, typename... Rest2_>
      void convert(const TupleVector<First2_, Rest2_...>& other)
      {
        this->first().convert(other.first());
        this->rest().convert(other.rest());
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes() + _rest.bytes();
      }

      /**
       * \brief Read in vector from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, String filename)
      {
        switch(mode)
        {
        case FileMode::fm_binary:
          read_from_binary(filename);
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
        switch(mode)
        {
        case FileMode::fm_binary:
          read_from_binary(file);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Read in vector from binary file.
       *
       * \param[in] filename The file that shall be read in.
       */
      void read_from_binary(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
        read_from_binary(file);
        file.close();
      }

      /**
       * \brief Read in vector from binary stream.
       *
       * \param[in] file The stream that shall be read in.
       */
      void read_from_binary(std::istream& file)
      {
        uint64_t magic; // magic_number
        file.read((char *)&magic, (long)(sizeof(uint64_t)));
        if (magic != 101)
          throw InternalError(__func__, __FILE__, __LINE__, "Given file or file component is no TupleVector!");
        uint64_t count; // subvector count
        file.read((char *)&count, (long)(sizeof(uint64_t)));
        //if (count != num_blocks)
        //  throw InternalError(__func__, __FILE__, __LINE__, "TupleVector file read in component count mismatch: class has " + stringify(num_blocks) + "- " + stringify(count) + " read in!");

        _read_from_binary(file);
      }

      /**
       * \brief Write out vector to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the vector shall be stored.
       */
      void write_out(FileMode mode, String filename) const
      {
        switch(mode)
        {
        case FileMode::fm_binary:
          write_out_binary(filename);
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
        switch(mode)
        {
        case FileMode::fm_binary:
          write_out_binary(file);
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
      void write_out_binary(String filename) const
      {
        std::ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
        write_out_binary(file);
        file.close();
      }

      /**
       * \brief Write out vector to file.
       *
       * \param[in] file The stream that shall be written to.
       *
       * Creates a binary file, that consists of a small header describing the type and all subvector's binary dumps.
       */
      void write_out_binary(std::ostream& file) const
      {
        size_t gsize(2 * sizeof(uint64_t)); // magic_number and subvector count
        std::vector<char> result(gsize);
        char * array(result.data());
        uint64_t * uiarray(reinterpret_cast<uint64_t *>(array));
        uiarray[0] = 101; /// \todo globale liste anlegen
        uiarray[1] = num_blocks;

        file.write(result.data(), long(result.size()));

        _write_out_binary(file);
      }
    }; // class TupleVector<...>

    /// \cond internal
    template<typename First_>
    class TupleVector<First_>
    {
    private:
      template<typename,typename...>
      friend class TupleVector;

      /// read binary data of _first to file.
      void _read_from_binary(std::istream& file)
      {
        //process first
        _first.read_from(FileMode::fm_binary, file);
      }

      /// write binary data of _first to file.
      void _write_out_binary(std::ostream& file) const
      {
        //process first
        _first.write_out(FileMode::fm_binary, file);
      }

    public:
      static constexpr int num_blocks = 1;

      typedef typename First_::MemType MemType;
      typedef typename First_::DataType DataType;
      typedef typename First_::IndexType IndexType;

      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = class TupleVector<typename First_::template ContainerType<Mem2_, DT2_, IT2_> >;

      /// this typedef lets you create a vector container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

    protected:
      First_ _first;

      static String sub_name_list()
      {
        return First_::name();
      }

    public:
      TupleVector()
      {
      }

      explicit TupleVector(First_&& the_first) :
        _first(std::move(the_first))
      {
      }

      /// move-ctor
      TupleVector(TupleVector&& other) :
        _first(std::move(other._first))
      {
      }

      /// move-assign operator
      TupleVector& operator=(TupleVector&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
        }
        return *this;
      }

      /// deleted copy-ctor
      TupleVector(const TupleVector&) = delete;
      /// deleted copy-assign operator
      TupleVector& operator=(const TupleVector&) = delete;

      /**
       * \brief Creates and returns a copy of this vector
       *
       * \param[in] mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       *
       */
      TupleVector clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return TupleVector(_first.clone(mode));
      }

      /**
       * \brief Turns this vector into a clone of other
       *
       * \param[in] clone_mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       *
       */
      void clone(const TupleVector& other, CloneMode clone_mode)
      {
        _first.clone(other._first, clone_mode);
      }

      /**
       * \brief Turns this vector into a clone of other
       *
       * As the default CloneMode of the underlying containers is unknown, this has to be separate.
       *
       */
      void clone(const TupleVector& other)
      {
        _first.clone(other._first);
      }

      First_& first()
      {
        return _first;
      }

      const First_& first() const
      {
        return _first;
      }
      template <Perspective perspective_ = Perspective::native>
      Index size() const
      {
        return _first.template size<perspective_>();
      }

      /// Returns the number of blocks in this tuple-vector.
      int blocks() const
      {
        return num_blocks;
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("TupleVector<") + sub_name_list() + ">";
      }

      template<int i_>
      typename TupleElement<i_, First_>::Type& at()
      {
        static_assert(i_ == 0, "invalid sub-vector index");
        return first();
      }

      template<int i_>
      typename TupleElement<i_, First_>::Type const& at() const
      {
        static_assert(i_ == 0, "invalid sub-vector index");
        return first();
      }

      void format(DataType value = DataType(0))
      {
        _first.format(value);
      }

      void clear()
      {
        _first.clear();
      }

      //template<typename First2_>
      void copy(const TupleVector/*<First2_>*/& x, bool full = false)
      {
        first().copy(x.first(), full);
      }

      void axpy(const TupleVector& x, const TupleVector& y, DataType alpha = DataType(1))
      {
        first().axpy(x.first(), y.first(), alpha);
      }

      void component_product(const TupleVector & x, const TupleVector & y)
      {
        first().component_product(x.first(), y.first());
      }

      void component_invert(const TupleVector& x, DataType alpha = DataType(1))
      {
        first().component_invert(x.first(), alpha);
      }

      void scale(const TupleVector& x, DataType alpha)
      {
        first().scale(x.first(), alpha);
      }

      DataType dot(const TupleVector& x) const
      {
        return first().dot(x.first());
      }

      DataType triple_dot(const TupleVector& x, const TupleVector& y) const
      {
        return first().triple_dot(x.first(), y.first());
      }

      DataType triple_dot_i(const TupleVector& x, const TupleVector& y) const
      {
        return first().triple_dot_i(x.first(), y.first());
      }

      DataType norm2sqr() const
      {
        return first().norm2sqr();
      }

      DataType norm2() const
      {
        return first().norm2();
      }

      DataType max_element() const
      {
        return first().max_element();
      }

      const typename First_::DataType operator()(Index index) const
      {
        ASSERT(index < size());
        return first()(index);
      }

      void operator()(Index index, typename First_::DataType value)
      {
        ASSERT(index < size());
        first()(index, value);
      }

      void set_vec(DataType * const pval_set) const
      {
        this->first().set_vec(pval_set);
      }

      void set_vec_inv(const DataType * const pval_set)
      {
        this->first().set_vec_inv(pval_set);
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other The source Vector.
       *
       * Use source vector content as content of current vector
       */
      template <typename First2_>
      void convert(const TupleVector<First2_>& other)
      {
        this->first().convert(other.first());
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _first.bytes();
      }

      /**
       * \brief Read in vector from file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file that shall be read in.
       */
      void read_from(FileMode mode, String filename)
      {
        switch(mode)
        {
        case FileMode::fm_binary:
          read_from_binary(filename);
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
        switch(mode)
        {
        case FileMode::fm_binary:
          read_from_binary(file);
          break;
        default:
          throw InternalError(__func__, __FILE__, __LINE__, "Filemode not supported!");
        }
      }

      /**
       * \brief Read in vector from binary file.
       *
       * \param[in] filename The file that shall be read in.
       */
      void read_from_binary(String filename)
      {
        std::ifstream file(filename.c_str(), std::ifstream::in | std::ifstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
        read_from_binary(file);
        file.close();
      }

      /**
       * \brief Read in vector from binary stream.
       *
       * \param[in] file The stream that shall be read in.
       */
      void read_from_binary(std::istream& file)
      {
        uint64_t magic; // magic_number
        file.read((char *)&magic, (long)(sizeof(uint64_t)));
        if (magic != 101)
          throw InternalError(__func__, __FILE__, __LINE__, "Given file or file component is no TupleVector!");
        uint64_t count; // subvector count
        file.read((char *)&count, (long)(sizeof(uint64_t)));
        if (count != 1)
          throw InternalError(__func__, __FILE__, __LINE__, "PowerVector file read in component count mismatch: class has 1 - " + stringify(count) + " read in!");

        _read_from_binary(file);
      }

      /**
       * \brief Write out vector to file.
       *
       * \param[in] mode The used file format.
       * \param[in] filename The file where the vector shall be stored.
       */
      void write_out(FileMode mode, String filename) const
      {
        switch(mode)
        {
        case FileMode::fm_binary:
          write_out_binary(filename);
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
        switch(mode)
        {
        case FileMode::fm_binary:
          write_out_binary(file);
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
      void write_out_binary(String filename) const
      {
        std::ofstream file(filename.c_str(), std::ofstream::out | std::ofstream::binary);
        if (! file.is_open())
          throw InternalError(__func__, __FILE__, __LINE__, "Unable to open Vector file " + filename);
        write_out_binary(file);
        file.close();
      }

      /**
       * \brief Write out vector to file.
       *
       * \param[in] file The stream that shall be written to.
       *
       * Creates a binary file, that consists of a small header describing the type and all subvector's binary dumps.
       */
      void write_out_binary(std::ostream& file) const
      {
        size_t gsize(2 * sizeof(uint64_t)); // magic_number and subvector count
        std::vector<char> result(gsize);
        char * array(result.data());
        uint64_t * uiarray(reinterpret_cast<uint64_t *>(array));
        uiarray[0] = 101; /// \todo globale liste anlegen
        uiarray[1] = 1; //fixed num_blocks for tuple vector specialisation

        file.write(result.data(), long(result.size()));

        _write_out_binary(file);
      }
    };
    /// \endcond

    /// \cond internal
    template <typename First_>
    inline void dump_tuple_vector(std::ostream & os, const TupleVector<First_>& x)
    {
      os << x.first();
    }

    template <typename First_, typename... Rest_>
    inline void dump_tuple_vector(std::ostream & os, const TupleVector<First_, Rest_...>& x)
    {
      os << x.first() << ",";
      dump_tuple_vector<Rest_...>(os, x.rest());
    }

    template <typename First_, typename... Rest_>
    inline std::ostream& operator<< (std::ostream & os, const TupleVector<First_, Rest_...>& x)
    {
      os << "[";
      dump_tuple_vector(os, x);
      os << "]";
      return os;
    }
    /// \endcond
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_TUPLE_VECTOR_HPP
