// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_TUPLE_MATRIX_HPP
#define KERNEL_LAFEM_TUPLE_MATRIX_HPP 1

#include <kernel/lafem/tuple_vector.hpp>


namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief TupleMatrix row helper class template
     *
     * This class template is a helper class for the TupleMatrix class template,
     * which is responsible for the management of a single "meta-row" of the
     * outer tuple matrix.
     *
     * \note
     * The implementation of this class is quite similar to that of TupleVector.
     *
     * \tparam First_, ...Rest_
     * A sequence of (meta) matrix classes which are to be composed.
     *
     * \author Peter Zajac
     */
    template<typename First_, typename... Rest_>
    class TupleMatrixRow
    {
    private:
      template<typename,typename...>
      friend class TupleMatrixRow;

      typedef TupleMatrixRow<Rest_...> RestClass;

    public:
      /// number of matrix columns
      static constexpr int num_blocks = TupleMatrixRow<Rest_...>::num_blocks + 1;

      /// sub-matrix mem-type
      typedef typename First_::MemType MemType;
      /// sub-matrix data-type
      typedef typename First_::DataType DataType;
      /// sub-matrix index-type
      typedef typename First_::IndexType IndexType;

      // ensure that all sub-matrix have the same mem- and data-type
      static_assert(std::is_same<MemType, typename RestClass::MemType>::value,
                    "sub-matrices have different mem-types");
      static_assert(std::is_same<DataType, typename RestClass::DataType>::value,
                    "sub-matrices have different data-types");
      static_assert(std::is_same<IndexType, typename RestClass::IndexType>::value,
                    "sub-matrices have different index-types");

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = TupleMatrixRow<
        typename First_::template ContainerType<Mem2_, DT2_, IT2_>,
        typename Rest_::template ContainerType<Mem2_, DT2_, IT2_>...>;

      /// this typedef lets you create a vector container with new Memory, Datatype and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

      /// Compatible L-vector type
      typedef typename First_::VectorTypeL VectorTypeL;
      /// Compatible R-vector type
      typedef TupleVector<typename First_::VectorTypeR, typename Rest_::VectorTypeR...> VectorTypeR;

    protected:
      /// the first sub-matrix
      First_ _first;
      /// all remaining sub-matrices
      RestClass _rest;

      /// Returns a list of all sub-matrix type names
      static String sub_name_list()
      {
        return First_::name() + "," + RestClass::sub_name_list();
      }

    public:
      /// default CTOR
      TupleMatrixRow()
      {
      }

      /// rest-class emplacement ctor; for internal use only
      explicit TupleMatrixRow(First_&& the_first, RestClass&& the_rest) :
        _first(std::forward<First_>(the_first)),
        _rest(std::forward<RestClass>(the_rest))
      {
      }

      /// Sub-Matrix emplacement constructor
      explicit TupleMatrixRow(First_&& the_first, Rest_&&... the_rest) :
        _first(std::forward<First_>(the_first)),
        _rest(std::forward<Rest_>(the_rest)...)
      {
      }

      /// move ctor
      TupleMatrixRow(TupleMatrixRow&& other) :
        _first(std::forward<First_>(other._first)),
        _rest(std::forward<RestClass>(other._rest))
      {
      }

      /// move-assign operator
      TupleMatrixRow& operator=(TupleMatrixRow&& other)
      {
        if(this != &other)
        {
          _first = std::forward<First_>(other._first);
          _rest = std::forward<RestClass>(other._rest);
        }
        return *this;
      }

      /// deleted copy-ctor
      TupleMatrixRow(const TupleMatrixRow&) = delete;
      /// deleted copy-assign operator
      TupleMatrixRow& operator=(const TupleMatrixRow&) = delete;

      /**
       * \brief Creates and returns a copy of this matrix
       *
       * \param[in] mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       */
      TupleMatrixRow clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return TupleMatrixRow(_first.clone(mode), _rest.clone(mode));
      }

      /**
       * \brief Turns this matrix into a clone of other
       *
       * \param[in] clone_mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       */
      void clone(const TupleMatrixRow& other, CloneMode clone_mode)
      {
        _first.clone(other._first, clone_mode);
        _rest.clone(other._rest, clone_mode);
      }

      /**
       * \brief Turns this matrix into a clone of other
       *
       * As the default CloneMode of the underlying containers is unknown, this has to be separate.
       */
      void clone(const TupleMatrixRow& other)
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

      TupleMatrixRow<Rest_...>& rest()
      {
        return _rest;
      }

      const TupleMatrixRow<Rest_...>& rest() const
      {
        return _rest;
      }
      /// \endcond

      /**
       * \brief Returns a sub-matrix block.
       *
       * \tparam i_
       * The index of the sub-matrix block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-matrix at position \p i_.
       */
      template<int i_>
      typename TupleElement<i_, First_, Rest_...>::Type& at()
      {
        static_assert((0 <= i_) && (i_ < num_blocks), "invalid sub-matrix index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /** \copydoc at() */
      template<int i_>
      typename TupleElement<i_, First_, Rest_...>::Type const& at() const
      {
        static_assert((0 <= i_) && (i_ < num_blocks), "invalid sub-matrix index");
        return TupleElement<i_, First_, Rest_...>::get(*this);
      }

      /**
       * \brief Returns the total number of rows in this matrix.
       *
       * \returns Matrix row count if perspective_ = false.
       * \returns Raw matrix row count if perspective_ = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index rows() const
      {
        const Index rows_f = _first.template rows<perspective_>();
        const Index rows_r = _rest.template rows<perspective_>();
        XASSERTM(rows_f == rows_r, "row count mismatch");
        return rows_f;
      }

      /**
       * \brief Returns the total number of columns in this matrix.
       *
       * \returns Matrix column count if perspective_ = false.
       * \returns Raw matrix column count if perspective_ = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index columns() const
      {
        return _first.template columns<perspective_>() + _rest.template columns<perspective_>();
      }

      /**
       * \brief Returns the total number of non-zeros in this matrix.
       *
       * \returns Matrix non zero element count if perspective_ = false.
       * \returns Raw matrix non zero element count if perspective_ = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index used_elements() const
      {
        return _first.template used_elements<perspective_>() + _rest.template used_elements<perspective_>();
      }

      /**
       * \brief Returns the total amount of bytes allocated.
       *
       * \returns The amount of bytes allocated in all arrays
       */
      std::size_t bytes() const
      {
        return _first.bytes() + _rest.bytes();
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("TupleMatrixRow<") + sub_name_list() + ">";
      }

      /**
       * \brief Reset all elements of the container to a given value or zero if missing.
       *
       * \param[in] value
       * The value to which the matrix's entries are to be set to.
       */
      void format(DataType value = DataType(0))
      {
        first().format(value);
        rest().format(value);
      }

      /**
       * \brief Reset all elements of the container to random values.
       *
       * \param[in] rng The random number generator.
       * \param[in] min Lower rng bound.
       * \param[in] max Upper rng bound.
       */
      void format(Random & rng, DataType min, DataType max)
      {
        first().format(rng, min, max);
        rest().format(rng, min, max);
      }

      /// Free all allocated arrays
      void clear()
      {
        first().clear();
        rest().clear();
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return first().create_vector_l();
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(first().create_vector_r(), rest().create_vector_r());
      }

      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        first().apply(r, x.first());
        rest().apply(r, x.rest(), r, DataType(1));
      }

      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        first().apply(r, x.first(), y, alpha);
        rest().apply(r, x.rest(), r, alpha);
      }

      Index get_length_of_line(const Index row) const
      {
        return first().get_length_of_line(row) + rest().get_length_of_line(row);
      }

      void set_line(const Index row, DataType * const pval_set, IndexType * const pcol_set, const Index col_start, const Index stride = 1) const
      {
        const Index first_length(first().get_length_of_line(row));

        first().set_line(row, pval_set, pcol_set, col_start, stride);
        rest().set_line(row, pval_set + stride * first_length, pcol_set + stride * first_length, col_start + first().template columns<Perspective::pod>(), stride);
      }

      void set_line_reverse(const Index row, DataType * const pval_set, const Index stride = 1)
      {
        const Index first_length(first().get_length_of_line(row));

        first().set_line_reverse(row, pval_set, stride);
        rest().set_line_reverse(row, pval_set + stride * first_length, stride);
      }

      /// \copydoc FEAT::Control::Checkpointable::get_checkpoint_size()
      std::uint64_t get_checkpoint_size(SerialConfig& config)
      {
        return sizeof(std::uint64_t) + _first.get_checkpoint_size(config) + _rest.get_checkpoint_size(config); //sizeof(std::uint64_t) bits needed to store lenght of checkpointed _first
      }

      /// \copydoc FEAT::Control::Checkpointable::restore_from_checkpoint_data(std::vector<char>&)
      void restore_from_checkpoint_data(std::vector<char> & data)
      {
        std::uint64_t isize = *(std::uint64_t*) data.data(); //get size of checkpointed _first
        std::vector<char>::iterator start = std::begin(data) + sizeof(std::uint64_t);
        std::vector<char>::iterator last_of_first = std::begin(data) + sizeof(std::uint64_t) + (int) isize;
        std::vector<char> buffer_first(start, last_of_first);
        _first.restore_from_checkpoint_data(buffer_first);

        data.erase(std::begin(data), last_of_first);
        _rest.restore_from_checkpoint_data(data);
      }

      /// \copydoc FEAT::Control::Checkpointable::set_checkpoint_data(std::vector<char>&)
      std::uint64_t set_checkpoint_data(std::vector<char>& data, SerialConfig& config)
      {
        std::size_t old_size = data.size();
        data.insert(std::end(data), sizeof(std::uint64_t), 0); //add placeholder
        std::uint64_t ireal_size = _first.set_checkpoint_data(data, config); //add data of _first to the overall checkpoint and save its size
        char* csize = reinterpret_cast<char*>(&ireal_size);
        for(std::size_t i(0) ; i < sizeof(std::uint64_t) ; ++i)  //overwrite the guessed datalength
        {
          data[old_size + i] = csize[i];
        }

        return sizeof(std::uint64_t) + ireal_size + _rest.set_checkpoint_data(data, config); //generate and add checkpoint data for the _rest
      }

      template<typename OtherFirst_, typename... OtherRest_>
      void convert(const TupleMatrixRow<OtherFirst_, OtherRest_...>& other)
      {
        first().convert(other.first());
        rest().convert(other.rest());
      }

      template<typename OtherFirst_, typename... OtherRest_>
      void convert_reverse(TupleMatrixRow<OtherFirst_, OtherRest_...>& other) const
      {
        first().convert_reverse(other.first());
        rest().convert_reverse(other.rest());
      }
    }; // template class TupleMatrixRow

    /// \cond internal
    // specialization for a single sub-matrix
    template<typename First_>
    class TupleMatrixRow<First_>
    {
    private:
      template<typename,typename...>
      friend class TupleMatrixRow;

    public:
      /// number of matrix columns
      static constexpr int num_blocks = 1;

      /// sub-matrix mem-type
      typedef typename First_::MemType MemType;
      /// sub-matrix data-type
      typedef typename First_::DataType DataType;
      /// sub-matrix index-type
      typedef typename First_::IndexType IndexType;

      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = TupleMatrixRow<typename First_::template ContainerType<Mem2_, DT2_, IT2_> >;

      /// this typedef lets you create a vector container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

      /// Compatible L-vector type
      typedef typename First_::VectorTypeL VectorTypeL;
      /// Compatible R-vector type
      typedef TupleVector<typename First_::VectorTypeR> VectorTypeR;

    protected:
      /// the first sub-matrix
      First_ _first;

      /// Returns a list of all sub-matrix type names
      static String sub_name_list()
      {
        return First_::name();
      }

    public:
      /// default CTOR
      TupleMatrixRow()
      {
      }

      /// Sub-Vector emplacement constructor
      explicit TupleMatrixRow(First_&& the_first) :
        _first(std::forward<First_>(the_first))
      {
      }

      /// move ctor
      TupleMatrixRow(TupleMatrixRow&& other) :
        _first(std::forward<First_>(other._first))
      {
      }

      /// move-assign operator
      TupleMatrixRow& operator=(TupleMatrixRow&& other)
      {
        if(this != &other)
        {
          _first = std::forward<First_>(other._first);
        }
        return *this;
      }

      /// deleted copy-ctor
      TupleMatrixRow(const TupleMatrixRow&) = delete;
      /// deleted copy-assign operator
      TupleMatrixRow& operator=(const TupleMatrixRow&) = delete;

      TupleMatrixRow clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return TupleMatrixRow(_first.clone(mode));
      }

      void clone(const TupleMatrixRow& other, CloneMode clone_mode)
      {
        _first.clone(other._first, clone_mode);
      }

      void clone(const TupleMatrixRow& other)
      {
        _first.clone(other._first);
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
      /// \endcond

      template<int i_>
      First_& at()
      {
        static_assert(i_ == 0, "invalid sub-matrix index");
        return first();
      }

      /** \copydoc at() */
      template<int i_>
      const First_& at() const
      {
        static_assert(i_ == 0, "invalid sub-matrix index");
        return first();
      }

      template <Perspective perspective_ = Perspective::native>
      Index rows() const
      {
        return _first.template rows<perspective_>();
      }

      template <Perspective perspective_ = Perspective::native>
      Index columns() const
      {
        return _first.template columns<perspective_>();
      }

      template <Perspective perspective_ = Perspective::native>
      Index used_elements() const
      {
        return _first.template used_elements<perspective_>();
      }

      std::size_t bytes() const
      {
        return _first.bytes();
      }

      static String name()
      {
        return String("TupleMatrixRow<") + sub_name_list() + ">";
      }

      void format(DataType value = DataType(0))
      {
        first().format(value);
      }

      void format(Random & rng, DataType min, DataType max)
      {
        first().format(rng, min, max);
      }

      void clear()
      {
        first().clear();
      }

      VectorTypeL create_vector_l() const
      {
        return first().create_vector_l();
      }

      VectorTypeR create_vector_r() const
      {
        return VectorTypeR(first().create_vector_r());
      }

      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        first().apply(r, x.first());
      }

      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        first().apply(r, x.first(), y, alpha);
      }

      Index get_length_of_line(const Index row) const
      {
        return first().get_length_of_line(row);
      }

      void set_line(const Index row, DataType * const pval_set, IndexType * const pcol_set, const Index col_start, const Index stride = 1) const
      {
        first().set_line(row, pval_set, pcol_set, col_start, stride);
      }

      void set_line_reverse(const Index row, DataType * const pval_set, const Index stride = 1)
      {
        first().set_line_reverse(row, pval_set, stride);
      }

      /// \copydoc FEAT::Control::Checkpointable::get_checkpoint_size()
      std::uint64_t get_checkpoint_size(SerialConfig& config)
      {
        return _first.get_checkpoint_size(config);
      }

      /// \copydoc FEAT::Control::Checkpointable::restore_from_checkpoint_data(std::vector<char>&)
      void restore_from_checkpoint_data(std::vector<char> & data)
      {
        _first.restore_from_checkpoint_data(data);
      }

      /// \copydoc FEAT::Control::Checkpointable::set_checkpoint_data(std::vector<char>&)
      std::uint64_t set_checkpoint_data(std::vector<char>& data, SerialConfig& config)
      {
        return _first.set_checkpoint_data(data, config);
      }

      template<typename OtherFirst_>
      void convert(const TupleMatrixRow<OtherFirst_>& other)
      {
        first().convert(other.first());
      }

      template<typename OtherFirst_>
      void convert_reverse(TupleMatrixRow<OtherFirst_>& other) const
      {
        first().convert_reverse(other.first());
      }
    }; // template class TupleMatrixRow
    /// \endcond

    /// \cond internal
    namespace Intern
    {
      // helper class: checks whether 'T_' is an instance of 'TupleMatrixRow<...>'
      template<typename T_>
      struct IsTMR_
      {
        static constexpr bool value = false;
      };

      template<typename... T_>
      struct IsTMR_<TupleMatrixRow<T_...>>
      {
        static constexpr bool value = true;
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief TupleMatrix element helper class template
     *
     * This class template is a helper which is used for the implementation
     * of the 'at' function template of the TupleMatrix class.
     *
     * \author Peter Zajac
     */
    template<int i_, int j_, typename FirstRow_, typename... RestRows_>
    class TupleMatrixElement
    {
    public:
      static_assert(i_ > 0, "invalid row index"); // i_=0 is specialized below
      static_assert(j_ >= 0, "invalid column index");

      typedef typename TupleMatrixElement<i_-1, j_, RestRows_...>::Type Type;

      template<typename Meta_>
      static Type& get(Meta_& meta)
      {
        return TupleMatrixElement<i_-1, j_, RestRows_...>::get(meta.rest());
      }

      template<typename Meta_>
      static const Type& get(const Meta_& meta)
      {
        return TupleMatrixElement<i_-1, j_, RestRows_...>::get(meta.rest());
      }
    };

    /// \cond internal
    // specialization for i_ = 0;
    // Note: This specialization explicitly assumes that 'FirstRow_' is an instance of
    // 'TupleMatrixRow', which is ensured by a static_assert in the TupleMatrix class.
    template<int j_, typename First_, typename... Rest_, typename... RestRows_>
    class TupleMatrixElement<0, j_, TupleMatrixRow<First_, Rest_...>, RestRows_... >
    {
    public:
      static_assert(j_ >= 0, "invalid column index");

      typedef typename TupleElement<j_, First_, Rest_...>::Type Type;

      template<typename Meta_>
      static Type& get(Meta_& meta)
      {
        return TupleElement<j_, First_, Rest_...>::get(meta.first());
      }

      template<typename Meta_>
      static const Type& get(const Meta_& meta)
      {
        return TupleElement<j_, First_, Rest_...>::get(meta.first());
      }
    };
    /// \endcond

    /**
     * \brief Variadic TupleMatrix class template
     *
     * \tparam FirstRow_, ...RestRows_
     * A sequence of TupleMatrixRow template class instances which are to be composed
     *
     * \attention
     * All template parameters supplied to this class \b must be instances of the
     * TupleMatrixRow class template, even if the TupleMatrix only has 1 meta-column!
     *
     * \author Peter Zajac
     */
    template<typename FirstRow_, typename... RestRows_>
    class TupleMatrix
    {
    private:
      template<typename,typename...>
      friend class TupleMatrix;

      typedef TupleMatrix<RestRows_...> RestClass;

      // sanity check: ensure that 'FirstRow_' is an instance of TupleMatrixRow
      // note: this is essential for the implementation of TupleMatrixElement!
      static_assert(Intern::IsTMR_<FirstRow_>::value,
        "template parameter of TupleMatrix must be TupleMatrixRow<...>");

    public:
      /// number of row blocks (vertical size)
      static constexpr int num_row_blocks = TupleMatrix<RestRows_...>::num_row_blocks + 1;
      /// number of column blocks (horizontal size)
      static constexpr int num_col_blocks = FirstRow_::num_blocks;

      // sanity check: number of column blocks must match for all rows
      static_assert(num_col_blocks == RestClass::num_col_blocks, "column blocks mismatch");

      /// sub-matrix mem-type
      typedef typename FirstRow_::MemType MemType;
      /// sub-matrix data-type
      typedef typename FirstRow_::DataType DataType;
      /// sub-matrix index-type
      typedef typename FirstRow_::IndexType IndexType;

      // ensure that all sub-matrix have the same mem- and data-type
      static_assert(std::is_same<MemType, typename RestClass::MemType>::value,
                    "sub-matrices have different mem-types");
      static_assert(std::is_same<DataType, typename RestClass::DataType>::value,
                    "sub-matrices have different data-types");
      static_assert(std::is_same<IndexType, typename RestClass::IndexType>::value,
                    "sub-matrices have different index-types");

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = TupleMatrix<
        typename FirstRow_::template ContainerType<Mem2_, DT2_, IT2_>,
        typename RestRows_::template ContainerType<Mem2_, DT2_, IT2_>...>;

      /// this typedef lets you create a vector container with new Memory, Datatype and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

      /// Compatible L-vector type
      typedef TupleVector<typename FirstRow_::VectorTypeL, typename RestRows_::VectorTypeL...> VectorTypeL;
      /// Compatible R-vector type
      typedef typename FirstRow_::VectorTypeR VectorTypeR;

      /// sanity check: vector types must match for all rows
      static_assert(std::is_same<VectorTypeR, typename RestClass::VectorTypeR>::value, "vector type mismatch");

    protected:
      /// the first sub-matrix row
      FirstRow_ _first;
      /// all remaining sub-matrix rows
      RestClass _rest;

      /// Returns a list of all sub-matrix type names
      static String sub_name_list()
      {
        return FirstRow_::name() + "," + RestClass::sub_name_list();
      }

    public:
      /// default CTOR
      TupleMatrix()
      {
      }

      /// rest-class emplacement ctor; for internal use only
      explicit TupleMatrix(FirstRow_&& the_first, RestClass&& the_rest) :
        _first(std::forward<FirstRow_>(the_first)),
        _rest(std::forward<RestClass>(the_rest))
      {
      }

      /// Sub-Matrix emplacement constructor
      explicit TupleMatrix(FirstRow_&& the_first, RestRows_&&... the_rest) :
        _first(std::forward<FirstRow_>(the_first)),
        _rest(std::forward<RestRows_>(the_rest)...)
      {
      }

      /// move ctor
      TupleMatrix(TupleMatrix&& other) :
        _first(std::forward<FirstRow_>(other._first)),
        _rest(std::forward<RestClass>(other._rest))
      {
      }

      /// move-assign operator
      TupleMatrix& operator=(TupleMatrix&& other)
      {
        if(this != &other)
        {
          _first = std::forward<FirstRow_>(other._first);
          _rest = std::forward<RestClass>(other._rest);
        }
        return *this;
      }

      /// deleted copy-ctor
      TupleMatrix(const TupleMatrix&) = delete;
      /// deleted copy-assign operator
      TupleMatrix& operator=(const TupleMatrix&) = delete;

      /**
       * \brief Creates and returns a copy of this matrix
       *
       * \param[in] mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       */
      TupleMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return TupleMatrix(_first.clone(mode), _rest.clone(mode));
      }

      /**
       * \brief Turns this matrix into a clone of other
       *
       * \param[in] clone_mode
       * Determines the type of clone returned (shallow, weak, layout, deep)
       */
      void clone(const TupleMatrix& other, CloneMode clone_mode)
      {
        _first.clone(other._first, clone_mode);
        _rest.clone(other._rest, clone_mode);
      }

      /**
       * \brief Turns this matrix into a clone of other
       *
       * As the default CloneMode of the underlying containers is unknown, this has to be separate.
       */
      void clone(const TupleMatrix& other)
      {
        _first.clone(other._first);
        _rest.clone(other._rest);
      }

      /// \cond internal
      FirstRow_& first()
      {
        return _first;
      }

      const FirstRow_& first() const
      {
        return _first;
      }

      TupleMatrix<RestRows_...>& rest()
      {
        return _rest;
      }

      const TupleMatrix<RestRows_...>& rest() const
      {
        return _rest;
      }
      /// \endcond

      /**
       * \brief Returns a sub-matrix block.
       *
       * \tparam i_, j_
       * The indices of the sub-matrix block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-matrix at position (i_, j_).
       */
      template<int i_, int j_>
      typename TupleMatrixElement<i_, j_, FirstRow_, RestRows_...>::Type& at()
      {
        static_assert((0 <= i_) && (i_ < num_row_blocks), "invalid sub-matrix row index");
        static_assert((0 <= j_) && (j_ < num_col_blocks), "invalid sub-matrix column index");
        return TupleMatrixElement<i_, j_, FirstRow_, RestRows_...>::get(*this);
      }

      /** \copydoc at() */
      template<int i_, int j_>
      const typename TupleMatrixElement<i_, j_, FirstRow_, RestRows_...>::Type& at() const
      {
        static_assert((0 <= i_) && (i_ < num_row_blocks), "invalid sub-matrix row index");
        static_assert((0 <= j_) && (j_ < num_col_blocks), "invalid sub-matrix column index");
        return TupleMatrixElement<i_, j_, FirstRow_, RestRows_...>::get(*this);
      }

      /**
       * \brief Returns the total number of rows in this matrix.
       *
       * \returns Matrix row count if perspective_ = false.
       * \returns Raw matrix row count if perspective_ = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index rows() const
      {
        return _first.template rows<perspective_>() + _rest.template rows<perspective_>();
      }

      /**
       * \brief Returns the total number of columns in this matrix.
       *
       * \returns Matrix column count if perspective_ = false.
       * \returns Raw matrix column count if perspective_ = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index columns() const
      {
        const Index cols_f = _first.template columns<perspective_>();
        const Index cols_r = _rest.template columns<perspective_>();
        XASSERTM(cols_f == cols_r, "column count mismatch");
        return cols_f;
      }

      /**
       * \brief Returns the total number of non-zeros in this matrix.
       *
       * \returns Matrix non zero element count if perspective_ = false.
       * \returns Raw matrix non zero element count if perspective_ = true.
       */
      template <Perspective perspective_ = Perspective::native>
      Index used_elements() const
      {
        return _first.template used_elements<perspective_>() + _rest.template used_elements<perspective_>();
      }

      /**
       * \brief Returns the total amount of bytes allocated.
       *
       * \returns The amount of bytes allocated in all arrays
       */
      std::size_t bytes() const
      {
        return _first.bytes() + _rest.bytes();
      }

      /// Returns a descriptive string for this container.
      static String name()
      {
        return String("TupleMatrix<") + sub_name_list() + ">";
      }

      /**
       * \brief Reset all elements of the container to a given value or zero if missing.
       *
       * \param[in] value
       * The value to which the matrix's entries are to be set to.
       */
      void format(DataType value = DataType(0))
      {
        first().format(value);
        rest().format(value);
      }

      /**
       * \brief Reset all elements of the container to random values.
       *
       * \param[in] rng The random number generator.
       * \param[in] min Lower rng bound.
       * \param[in] max Upper rng bound.
       */
      void format(Random & rng, DataType min, DataType max)
      {
        first().format(rng, min, max);
        rest().format(rng, min, max);
      }

      /// Free all allocated arrays
      void clear()
      {
        first().clear();
        rest().clear();
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(first().create_vector_l(), rest().create_vector_l());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return first().create_vector_r();
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        first().apply(r.first(), x);
        rest().apply(r.rest(), x);
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        first().apply(r.first(), x, y.first(), alpha);
        rest().apply(r.rest(), x, y.rest(), alpha);
      }

      Index get_length_of_line(const Index row) const
      {
        const Index first_rows = first().template rows<Perspective::pod>();
        if(row < first_rows)
          return first().get_length_of_line(row);
        else
          return rest().get_length_of_line(row - first_rows);
      }

      void set_line(const Index row, DataType * const pval_set, IndexType * const pcol_set, const Index col_start, const Index stride = 1) const
      {
        const Index first_rows = first().template rows<Perspective::pod>();
        if(row < first_rows)
          first().set_line(row, pval_set, pcol_set, col_start, stride);
        else
          rest().set_line(row - first_rows, pval_set, pcol_set, col_start, stride);
      }

      void set_line_reverse(const Index row, DataType * const pval_set, const Index stride = 1)
      {
        const Index first_rows = first().template rows<Perspective::pod>();
        if(row < first_rows)
          first().set_line_reverse(row, pval_set, stride);
        else
          rest().set_line_reverse(row - first_rows, pval_set, stride);
      }

      /// \copydoc FEAT::Control::Checkpointable::get_checkpoint_size()
      std::uint64_t get_checkpoint_size(SerialConfig& config)
      {
        return sizeof(std::uint64_t) + _first.get_checkpoint_size(config) + _rest.get_checkpoint_size(config); //sizeof(std::uint64_t) bits needed to store lenght of checkpointed _first
      }

      /// \copydoc FEAT::Control::Checkpointable::restore_from_checkpoint_data(std::vector<char>&)
      void restore_from_checkpoint_data(std::vector<char> & data)
      {
        std::uint64_t isize = *(std::uint64_t*) data.data(); //get size of checkpointed _first
        std::vector<char>::iterator start = std::begin(data) + sizeof(std::uint64_t);
        std::vector<char>::iterator last_of_first = std::begin(data) + sizeof(std::uint64_t) + (int) isize;
        std::vector<char> buffer_first(start, last_of_first);
        _first.restore_from_checkpoint_data(buffer_first);

        data.erase(std::begin(data), last_of_first);
        _rest.restore_from_checkpoint_data(data);
      }

      /// \copydoc FEAT::Control::Checkpointable::set_checkpoint_data(std::vector<char>&)
      std::uint64_t set_checkpoint_data(std::vector<char>& data, SerialConfig& config)
      {
        std::size_t old_size = data.size();
        data.insert(std::end(data), sizeof(std::uint64_t), 0); //add placeholder
        std::uint64_t ireal_size = _first.set_checkpoint_data(data, config); //add data of _first to the overall checkpoint and save its size
        char* csize = reinterpret_cast<char*>(&ireal_size);
        for(std::size_t i(0) ; i < sizeof(std::uint64_t) ; ++i)  //overwrite the guessed datalength
        {
          data[old_size +i] = csize[i];
        }
        return sizeof(std::uint64_t) + ireal_size + _rest.set_checkpoint_data(data, config); //generate and add checkpoint data for the _rest
      }

      template<typename OtherFirstRow_, typename... OtherRestRows_>
      void convert(const TupleMatrix<OtherFirstRow_, OtherRestRows_...>& other)
      {
        first().convert(other.first());
        rest().convert(other.rest());
      }

      template<typename OtherFirstRow_, typename... OtherRestRows_>
      void convert_reverse(TupleMatrix<OtherFirstRow_, OtherRestRows_...>& other) const
      {
        first().convert_reverse(other.first());
        rest().convert_reverse(other.rest());
      }
    }; // template class TupleMatrix

    /// \cond internal
    // specialization for 1 row
    template<typename FirstRow_>
    class TupleMatrix<FirstRow_>
    {
    private:
      template<typename,typename...>
      friend class TupleMatrix;

      // sanity check: ensure that 'FirstRow_' is an instance of TupleMatrixRow
      static_assert(Intern::IsTMR_<FirstRow_>::value,
        "template parameter of TupleMatrix must be TupleMatrixRow<...>");

    public:
      /// number of matrix rows
      static constexpr int num_blocks = 1;

      /// number of row blocks (vertical size)
      static constexpr int num_row_blocks = num_blocks;
      /// number of column blocks (horizontal size)
      static constexpr int num_col_blocks = FirstRow_::num_blocks;


      /// sub-matrix mem-type
      typedef typename FirstRow_::MemType MemType;
      /// sub-matrix data-type
      typedef typename FirstRow_::DataType DataType;
      /// sub-matrix index-type
      typedef typename FirstRow_::IndexType IndexType;

      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = TupleMatrix<typename FirstRow_::template ContainerType<Mem2_, DT2_, IT2_> >;

      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

      /// Compatible L-vector type
      typedef TupleVector<typename FirstRow_::VectorTypeL> VectorTypeL;
      /// Compatible R-vector type
      typedef typename FirstRow_::VectorTypeR VectorTypeR;

    protected:
      /// the first sub-matrix
      FirstRow_ _first;

      /// Returns a list of all sub-matrix type names
      static String sub_name_list()
      {
        return FirstRow_::name();
      }

    public:
      /// default CTOR
      TupleMatrix()
      {
      }

      /// Sub-Vector emplacement constructor
      explicit TupleMatrix(FirstRow_&& the_first) :
        _first(std::forward<FirstRow_>(the_first))
      {
      }

      /// move ctor
      TupleMatrix(TupleMatrix&& other) :
        _first(std::forward<FirstRow_>(other._first))
      {
      }

      /// move-assign operator
      TupleMatrix& operator=(TupleMatrix&& other)
      {
        if(this != &other)
        {
          _first = std::forward<FirstRow_>(other._first);
        }
        return *this;
      }

      /// deleted copy-ctor
      TupleMatrix(const TupleMatrix&) = delete;
      /// deleted copy-assign operator
      TupleMatrix& operator=(const TupleMatrix&) = delete;

      TupleMatrix clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return TupleMatrix(_first.clone(mode));
      }

      void clone(const TupleMatrix& other, CloneMode clone_mode)
      {
        _first.clone(other._first, clone_mode);
      }

      void clone(const TupleMatrix& other)
      {
        _first.clone(other._first);
      }

      FirstRow_& first()
      {
        return _first;
      }

      const FirstRow_& first() const
      {
        return _first;
      }

      template<int i_, int j_>
      typename TupleMatrixElement<i_, j_, FirstRow_>::Type& at()
      {
        static_assert(0 == i_, "invalid sub-matrix row index");
        static_assert((0 <= j_) && (j_ < num_col_blocks), "invalid sub-matrix column index");
        return TupleMatrixElement<i_, j_, FirstRow_>::get(*this);
      }

      template<int i_, int j_>
      const typename TupleMatrixElement<i_, j_, FirstRow_>::Type& at() const
      {
        static_assert(0 == i_, "invalid sub-matrix row index");
        static_assert((0 <= j_) && (j_ < num_col_blocks), "invalid sub-matrix column index");
        return TupleMatrixElement<i_, j_, FirstRow_>::get(*this);
      }

      template <Perspective perspective_ = Perspective::native>
      Index rows() const
      {
        return _first.template rows<perspective_>();
      }

      template <Perspective perspective_ = Perspective::native>
      Index columns() const
      {
        return _first.template columns<perspective_>();
      }

      template <Perspective perspective_ = Perspective::native>
      Index used_elements() const
      {
        return _first.template used_elements<perspective_>();
      }

      std::size_t bytes() const
      {
        return _first.bytes();
      }

      static String name()
      {
        return String("TupleMatrix<") + sub_name_list() + ">";
      }

      void format(DataType value = DataType(0))
      {
        first().format(value);
      }

      void format(Random & rng, DataType min, DataType max)
      {
        first().format(rng, min, max);
      }

      void clear()
      {
        first().clear();
      }

      /// Returns a new compatible L-Vector.
      VectorTypeL create_vector_l() const
      {
        return VectorTypeL(first().create_vector_l());
      }

      /// Returns a new compatible R-Vector.
      VectorTypeR create_vector_r() const
      {
        return first().create_vector_r();
      }

      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        first().apply(r.first(), x);
      }

      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, DataType alpha = DataType(1)) const
      {
        first().apply(r.first(), x, y.first(), alpha);
      }

      Index get_length_of_line(const Index row) const
      {
        return first().get_length_of_line(row);
      }

      void set_line(const Index row, DataType * const pval_set, IndexType * const pcol_set, const Index col_start, const Index stride = 1) const
      {
        first().set_line(row, pval_set, pcol_set, col_start, stride);
      }

      void set_line_reverse(const Index row, DataType * const pval_set, const Index stride = 1)
      {
        first().set_line_reverse(row, pval_set, stride);
      }

      /// \copydoc FEAT::Control::Checkpointable::get_checkpoint_size()
      std::uint64_t get_checkpoint_size(SerialConfig& config)
      {
        return _first.get_checkpoint_size(config);
      }

      /// \copydoc FEAT::Control::Checkpointable::restore_from_checkpoint_data(std::vector<char>&)
      void restore_from_checkpoint_data(std::vector<char> & data)
      {
        _first.restore_from_checkpoint_data(data);
      }

      /// \copydoc FEAT::Control::Checkpointable::set_checkpoint_data(std::vector<char>&)
      std::uint64_t set_checkpoint_data(std::vector<char>& data, SerialConfig& config)
      {
        return _first.set_checkpoint_data(data, config);
      }

      template<typename OtherFirstRow_>
      void convert(const TupleMatrix<OtherFirstRow_>& other)
      {
        first().convert(other.first());
      }

      template<typename OtherFirstRow_>
      void convert_reverse(TupleMatrix<OtherFirstRow_>& other) const
      {
        first().convert_reverse(other.first());
      }
    }; // template class TupleMatrix
    /// \endcond
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_TUPLE_MATRIX_HPP
