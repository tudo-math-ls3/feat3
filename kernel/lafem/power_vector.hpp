#pragma once
#ifndef KERNEL_LAFEM_POWER_VECTOR_HPP
#define KERNEL_LAFEM_POWER_VECTOR_HPP 1

// includes, FEAST
#include <kernel/lafem/meta_element.hpp>
#include <kernel/lafem/container.hpp>

// includes, system
#include <iostream>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Power-Vector meta class template
     *
     * This class template implements a composition of \e n sub-vectors of the same class.
     *
     * \note For a composed meta-vector of different sub-vector classes, see the TupleVector class template.
     *
     * \tparam SubType_
     * The type of the sub-vector.
     *
     * \tparam count_
     * The number of sub-vector blocks.
     *
     * \author Peter Zajac
     */
    template<
      typename SubType_,
      int count_>
    class PowerVector
    {
    private:
      // Note: the case = 1 is specialised below
      static_assert(count_ > 1, "invalid block size");

      // declare this class template as a friend for recursive inheritance
      template<typename, int>
      friend class PowerVector;

      /// base-class typedef
      typedef PowerVector<SubType_, count_ - 1> RestClass;

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
      /// sub-vector type
      typedef SubType_ SubVectorType;
      /// sub-vector memory type
      typedef typename SubVectorType::MemType MemType;
      /// sub-vector data type
      typedef typename SubVectorType::DataType DataType;
      /// sub-vector index type
      typedef typename SubVectorType::IndexType IndexType;
      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = class PowerVector<typename SubType_::template ContainerType<Mem2_, DT2_, IT2_>, count_>;

      /// number of vector blocks
      static constexpr int num_blocks = count_;

    protected:
      /// the first sub-vector
      SubVectorType _first;
      /// the remaining part
      RestClass _rest;

    public:
      /// default CTOR
      PowerVector()
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] sub_size
       * The size of any sub-vector of this power-vector.
       */
      explicit PowerVector(Index sub_size) :
        _first(sub_size),
        _rest(sub_size)
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] sub_size
       * The size of any sub-vector of this power-vector.
       *
       * \param[in] value
       * The value that the sub-vectors are to be initialised with.
       */
      explicit PowerVector(Index sub_size, DataType value) :
        _first(sub_size, value),
        _rest(sub_size, value)
      {
      }

      /// base-class constructor; for internal use only
      explicit PowerVector(SubVectorType&& the_first, RestClass&& the_rest) :
        _first(std::move(the_first)),
        _rest(std::move(the_rest))
      {
      }

      /// move CTOR
      PowerVector(PowerVector&& other) :
        _first(std::move(other._first)),
        _rest(std::move(other._rest))
      {
      }

      /// move-assign operator
      PowerVector& operator=(PowerVector&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
          _rest = std::move(other._rest);
        }
        return *this;
      }

      /// deleted copy-ctor
      PowerVector(const PowerVector&) = delete;
      /// deleted copy-assign operator
      PowerVector& operator=(const PowerVector&) = delete;

      /// virtual destructor
      virtual ~PowerVector()
      {
      }

      /**
       * \brief Creates and returns a copy of this vector
       *
       * \param[in] mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       *
       */
      PowerVector clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return PowerVector(_first.clone(mode), _rest.clone(mode));
      }

      /**
       * \brief Turns this vector into a clone of other
       *
       * \param[in] clone_mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       *
       */
      void clone(const PowerVector& other, CloneMode clone_mode)
      {
        _first.clone(other._first, clone_mode);
        _rest.clone(other._rest, clone_mode);
      }

      /**
       * \brief Turns this vector into a clone of other
       *
       * As the default CloneMode of the underlying containers is unknown, this has to be seperate.
       *
       */
      void clone(const PowerVector& other)
      {
        _first.clone(other._first);
        _rest.clone(other._rest);
      }

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
      SubVectorType& at()
      {
        static_assert((0 <= i_) && (i_ < count_), "invalid sub-vector index");
        return PowerElement<i_, SubVectorType>::get(*this);
      }

      /** \copydoc at() */
      template<int i_>
      const SubVectorType& at() const
      {
        static_assert((0 <= i_) && (i_ < count_), "invalid sub-vector index");
        return PowerElement<i_, SubVectorType>::get(*this);
      }

      /**
       * \brief Returns a sub-vector block.
       *
       * \param[in] i
       * The index of the sub-vector block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-vector at position \p i.
       */
      SubVectorType& get(int i)
      {
        ASSERT_((0 <= i) && (i < count_));
        return (i == 0) ? _first : _rest.get(i-1);
      }

      /** \copydoc get() */
      const SubVectorType& get(int i) const
      {
        ASSERT_((0 <= i) && (i < count_));
        return (i == 0) ? _first : _rest.get(i-1);
      }

      /// \cond internal
      SubVectorType& first()
      {
        return _first;
      }

      const SubVectorType& first() const
      {
        return _first;
      }

      RestClass& rest()
      {
        return _rest;
      }

      const RestClass& rest() const
      {
        return _rest;
      }
      /// \endcond

      /// Returns the number of blocks in this power-vector.
      int blocks() const
      {
        return num_blocks;
      }

      /// Returns the total number of scalar entries of this power-vector.
      template <Perspective perspective_ = Perspective::native>
      Index size() const
      {
        return first().template size<perspective_>() + rest().template size<perspective_>();
      }

      /**
       * \brief Clears this vector.
       *
       * \param[in] value
       * The value to which the vector's entries are to be set to.
       */
      void format(DataType value = DataType(0))
      {
        first().format(value);
        rest().format(value);
      }

      /// Clears the vector.
      void clear()
      {
        first().clear();
        rest().clear();
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("PowerVector<") + SubVectorType::name() + "," + stringify(count_) + ">";
      }

      /**
       * \brief Performs \f$this \leftarrow x\f$.
       *
       * \param[in] x The vector to be copied.
       * \param[in] full Shall we create a full copy, including scalars and index arrays?
       */
      template<typename SubType2_>
      void copy(const PowerVector<SubType2_, count_>& x, bool full = false)
      {
        first().copy(x.first(), full);
        rest().copy(x.rest(), full);
      }

      /**
       * \brief Performs \f$ this \leftarrow \alpha\cdot x + y \f$
       *
       * \param[in] x
       * The first summand vector.
       *
       * \param[in] y
       * The second summand vector.
       *
       * \param[in] alpha
       * The scaling factor for \p x.
       */
      void axpy(const PowerVector& x, const PowerVector& y, DataType alpha = DataType(1))
      {
        first().axpy(x.first(), y.first(), alpha);
        rest().axpy(x.rest(), y.rest(), alpha);
      }

      /**
       * \brief Performs \f$this_i \leftarrow x_i \cdot y_i\f$
       *
       * \param[in] x The first factor.
       * \param[in] y The second factor.
       */
      void component_product(const PowerVector & x, const PowerVector & y)
      {
        first().component_product(x.first(), y.first());
        rest().component_product(x.rest(), y.rest());
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
      void component_invert(const PowerVector& x, DataType alpha = DataType(1))
      {
        first().component_invert(x.first(), alpha);
        rest().component_invert(x.rest(), alpha);
      }

      /**
       * \brief Performs \f$ this \leftarrow \alpha\cdot x \f$
       *
       * \param[in] x
       * The vector to be scaled.
       *
       * \param[in] alpha
       * The scaling factor for \p x.
       */
      void scale(const PowerVector& x, DataType alpha)
      {
        first().scale(x.first(), alpha);
        rest().scale(x.rest(), alpha);
      }

      /**
       * \brief Returns the dot-product of this vector and another vector.
       *
       * \param[in] x
       * The second vector for the dot-product.
       */
      DataType dot(const PowerVector& x) const
      {
        return first().dot(x.first()) + rest().dot(x.rest());
      }

      /**
       * \copydoc LAFEM::DenseVector::triple_dot()
       **/
      DataType triple_dot(const PowerVector& x, const PowerVector& y) const
      {
        return first().triple_dot(x.first(), y.first())
          + rest().triple_dot(x.rest(), y.rest());
      }

      /**
       * \brief Returns the squared euclid norm of this vector.
       */
      DataType norm2sqr() const
      {
        return first().norm2sqr() + rest().norm2sqr();
      }

      /**
       * \brief Returns the euclid norm of this vector.
       */
      DataType norm2() const
      {
        return Math::sqrt(norm2sqr());
      }

      /**
       * \brief Retrieve specific PowerVector element.
       *
       * \param[in] index The index of the vector element.
       *
       * \returns Specific vector element.
       */
      const DataType operator()(Index index) const
      {
        CONTEXT("When retrieving PowerVector element");

        ASSERT(index < size(), "Error: " + stringify(index) + " exceeds power vector size " + stringify(size()) + " !");

        if (index < first().size())
        {
          return first()(index);
        }
        else
        {
          return rest()(index - first().size());
        }
      }

      /**
       * \brief Set specific PowerVector element.
       *
       * \param[in] index The index of the vector element.
       * \param[in] value The value to be set.
       */
      void operator()(Index index, DataType value)
      {
        CONTEXT("When retrieving PowerVector element");

        ASSERT(index < size(), "Error: " + stringify(index) + " exceeds power vector size " + stringify(size()) + " !");

        if (index < first().size())
        {
          first()(index, value);
        }
        else
        {
          rest()(index - first().size(), value);
        }
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
       * \brief Convertion method
       *
       * \param[in] other The source Vector.
       *
       * Use source vector content as content of current vector
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const ContainerType<Mem2_, DT2_, IT2_> & other)
      {
        CONTEXT("When converting PowerVector");

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
        CONTEXT("When reading in PowerVector");

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
        CONTEXT("When reading in PowerVector");

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
        if (magic != 100)
          throw InternalError(__func__, __FILE__, __LINE__, "Given file or file component is no PowerVector!");
        uint64_t count; // subvector count
        file.read((char *)&count, (long)(sizeof(uint64_t)));
        if (count != count_)
          throw InternalError(__func__, __FILE__, __LINE__, "PowerVector file read in component count missmatch: class has " + stringify(count_) + "- " + stringify(count) + " read in!");

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
        CONTEXT("When writing out PowerVector");

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
        CONTEXT("When writing out PowerVector");

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
        uiarray[0] = 100; /// \todo globale liste anlegen
        uiarray[1] = count_;

        file.write(result.data(), long(result.size()));

        _write_out_binary(file);
      }
    }; // class PowerVector<...>

    /// \cond internal
    /**
     * \brief Specialisation of PowerVector for only one block
     *
     * \author Peter Zajac
     */
    template<typename SubType_>
    class PowerVector<SubType_, 1>
    {
    private:

      template<typename, int>
      friend class PowerVector;

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
      typedef SubType_ SubVectorType;

      typedef typename SubVectorType::MemType MemType;
      typedef typename SubVectorType::DataType DataType;
      typedef typename SubVectorType::IndexType IndexType;
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = class PowerVector<typename SubType_::template ContainerType<Mem2_, DT2_, IT2_>, Index(1)>;

      static constexpr int num_blocks = 1;

    protected:
      SubVectorType _first;

    public:
      PowerVector()
      {
      }

      explicit PowerVector(Index sub_size)
        : _first(sub_size)
      {
      }

      explicit PowerVector(Index sub_size, DataType value)
        : _first(sub_size, value)
      {
      }

      explicit PowerVector(SubVectorType&& the_first) :
        _first(std::move(the_first))
      {
      }

      PowerVector(PowerVector&& other) :
        _first(std::move(other._first))
      {
      }

      PowerVector& operator=(PowerVector&& other)
      {
        if(this != &other)
        {
          _first = std::move(other._first);
        }
        return *this;
      }

      /// deleted copy-ctor
      PowerVector(const PowerVector&) = delete;
      /// deleted copy-assign operator
      PowerVector& operator=(const PowerVector&) = delete;

      virtual ~PowerVector()
      {
      }

      /**
       * \brief Creates and returns a copy of this vector
       *
       * \param[in] mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       *
       */
      PowerVector clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return PowerVector(_first.clone(mode));
      }

      /**
       * \brief Turns this vector into a clone of other
       *
       * \param[in] clone_mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       *
       */
      void clone(const PowerVector& other, CloneMode clone_mode)
      {
        _first.clone(other._first, clone_mode);
      }

      /**
       * \brief Turns this vector into a clone of other
       *
       * As the default CloneMode of the underlying containers is unknown, this has to be seperate.
       *
       */
      void clone(const PowerVector& other)
      {
        _first.clone(other._first);
      }

      template<int i>
      SubVectorType& at()
      {
        static_assert(i == 0, "invalid sub-vector index");
        return _first;
      }

      template<int i>
      const SubVectorType& at() const
      {
        static_assert(i == 0, "invalid sub-vector index");
        return _first;
      }

      SubVectorType& get(int i)
      {
#ifndef DEBUG
        (void)i;
#endif
        ASSERT_(i == 0);
        return _first;
      }

      const SubVectorType& get(int i) const
      {
#ifndef DEBUG
        (void)i;
#endif
        ASSERT_(i == 0);
        return _first;
      }

      SubVectorType& first()
      {
        return _first;
      }

      const SubVectorType& first() const
      {
        return _first;
      }

      int blocks() const
      {
        return 1;
      }

      template <Perspective perspective_ = Perspective::native>
      Index size() const
      {
        return _first.template size<perspective_>();
      }

      void format(DataType value = DataType(0))
      {
        _first.format(value);
      }

      void clear()
      {
        first().clear();
      }

      static String name()
      {
        return String("PowerVector<") + SubVectorType::name() + ",1>";
      }

      template<typename SubType2_>
      void copy(const PowerVector<SubType2_, 1>& x, bool full = false)
      {
        first().copy(x.first(), full);
      }

      void axpy(const PowerVector& x, const PowerVector& y, DataType alpha = DataType(1))
      {
        first().axpy(x.first(), y.first(), alpha);
      }

      void component_product(const PowerVector & x, const PowerVector & y)
      {
        first().component_product(x.first(), y.first());
      }

      void component_invert(const PowerVector& x, DataType alpha = DataType(1))
      {
        first().component_invert(x.first(), alpha);
      }

      void scale(const PowerVector& x, DataType alpha)
      {
        first().scale(x.first(), alpha);
      }

      DataType dot(const PowerVector& x) const
      {
        return first().dot(x.first());
      }

      DataType triple_dot(const PowerVector& x, const PowerVector& y) const
      {
        return first().triple_dot(x.first(), y.first());
      }

      DataType triple_dot_i(const PowerVector& x, const PowerVector& y) const
      {
        return first().triple_dot_i(x.first(), y.first());
      }

      DataType norm2sqr() const
      {
        return first().norm2sqr();
      }

      DataType norm2() const
      {
        return Math::sqrt(norm2sqr());
      }

      const DataType operator()(Index index) const
      {
        CONTEXT("When retrieving PowerVector element");

        ASSERT(index < size(), "Error: " + stringify(index) + " exceeds power vector size " + stringify(size()) + " !");

        return first()(index);
      }

      void operator()(Index index, DataType value)
      {
        CONTEXT("When retrieving PowerVector element");

        ASSERT(index < size(), "Error: " + stringify(index) + " exceeds power vector size " + stringify(size()) + " !");

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
       * \brief Convertion method
       *
       * \param[in] other The source Vector.
       *
       * Use source vector content as content of current vector
       */
      template <typename Mem2_, typename DT2_, typename IT2_>
      void convert(const ContainerType<Mem2_, DT2_, IT2_> & other)
      {
        CONTEXT("When converting PowerVector");

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
        CONTEXT("When reading in PowerVector");

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
        CONTEXT("When reading in PowerVector");

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
        if (magic != 100)
          throw InternalError(__func__, __FILE__, __LINE__, "Given file or file component is no PowerVector!");
        uint64_t count; // subvector count
        file.read((char *)&count, (long)(sizeof(uint64_t)));
        if (count != 1)
          throw InternalError(__func__, __FILE__, __LINE__, "PowerVector file read in component count missmatch: class has 1 - " + stringify(count) + " read in!");

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
        CONTEXT("When writing out PowerVector");

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
        CONTEXT("When writing out PowerVector");

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
        uiarray[0] = 100; /// \todo globale liste anlegen
        uiarray[1] = 1; //fixed count_ for power vector specialisation

        file.write(result.data(), long(result.size()));

        _write_out_binary(file);
      }
    }; // class PowerVector<...,1>
    /// \cond internal

    /// \cond internal
    template <typename SubType_, Index count_>
    inline void dump_power_vector(std::ostream & os, const PowerVector<SubType_, count_>& x)
    {
      os << "," << x.first();
      dump_power_vector(os, x.rest());
    }

    template <typename SubType_>
    inline void dump_power_vector(std::ostream & os, const PowerVector<SubType_, 1>& x)
    {
      os << x.first();
    }

    template <typename SubType_, Index count_>
    inline std::ostream& operator<< (std::ostream & os, const PowerVector<SubType_, count_>& x)
    {
      os << "[";
      dump_power_vector(os, x);
      os << "]";
      return os;
    }
    /// \endcond
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_POWER_VECTOR_HPP
