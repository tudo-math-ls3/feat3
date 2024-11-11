// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <deque>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Sequence of filters of the same type
     *
     * \tparam Filter_
     * Type of the filter that gets applied in sequence
     *
     * This is a thin wrapper around an std::deque<std::pair<String, Filter_>> that just adds filter operations to
     * it. We use an std::pair<String, Filter_> so that every member of the sequence has an identifier (usually the
     * name of the MeshPart it acts on) so one can find the correct assembler for it.
     * We cannot use an std::map because the order of the filters in the sequence might be important (i.e. if one has
     * singularities in the essential boundary conditions enforced by a UnitFilter, as it is the case for the
     * Stokes problem with driven cavity).
     *
     * \author Jordi Paul
     */
    template<typename Filter_>
    class FilterSequence : public std::deque<std::pair<String, Filter_>>
    {
    public:
      /// Our baseclass
      typedef std::deque<std::pair<String, Filter_>> BaseClass;

      /// The type of the filter that gets applied in sequence
      typedef Filter_ InnerFilterType;
      /// The filter's floating point type
      typedef typename Filter_::DataType DataType;
      /// The filter's index type
      typedef typename Filter_::IndexType IndexType;

      /// The type of vector the filter can be applied to
      typedef typename Filter_::VectorType VectorType;

      /// Use the baseclass' iterator
      typedef typename BaseClass::iterator iterator;
      /// Use the baseclass' const_iterator
      typedef typename BaseClass::const_iterator const_iterator;

      /// Our 'base' class type
      template <typename DT2_, typename IT2_>
      using FilterType = FilterSequence<typename Filter_::template FilterType<DT2_, IT2_> >;

      /// this typedef lets you create a filter sequence with new Data and Index types
      template <typename DT2_, typename IT2_>
      using FilterTypeByDI = FilterType<DT2_, IT2_>;

    public:
      /// Empty default constructor
      FilterSequence()
      {
      }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~FilterSequence()
      {
      }

      /**
       * \brief Constructs a FilterSequence of empty filters
       *
       * \param[in] ids_
       * Identifiers for the empty filters.
       *
       */
      explicit FilterSequence(const std::deque<String>& ids_)
      {
        for(auto it = ids_.begin(); it != ids_.end(); ++it)
        {
          String id(*it);
          InnerFilterType new_filter;
          BaseClass::push_back(std::make_pair<String, Filter_>(std::move(id), std::move(new_filter)));
        }
      }

      FilterSequence(FilterSequence&& other) :
        BaseClass(std::forward<BaseClass>(other))
      {
      }

      /// move-assign operator
      FilterSequence& operator=(FilterSequence&& other)
      {
        if(this != &other)
        {
          static_cast<BaseClass&>(*this) = std::forward<BaseClass>(other);
        }

        return *this;
      }

      /// \brief Creates a clone of itself
      FilterSequence clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        FilterSequence result;
        result.resize(this->size());

        for(size_t i(0); i < result.size(); ++i)
        {
          result.at(i).first = this->at(i).first;
          result.at(i).second.clone(this->at(i).second, clone_mode);
        }

        return result;
      }

      /// \brief Clones data from another FilterSequence
      void clone(const FilterSequence& other, CloneMode clone_mode = CloneMode::Deep)
      {
        BaseClass::clear();

        for(const_iterator it(other.begin()); it != other.end(); ++it)
        {
          BaseClass::push_back(std::make_pair<String, InnerFilterType>(it->first, it->second->clone(clone_mode)));
        }
      }

      /**
       * \brief Converts another FilterSequence to this type
       */
      template<typename OtherFilter_>
      void convert(const FilterSequence<OtherFilter_>& other)
      {
        BaseClass::clear();

        BaseClass::resize(other.size());

        for(size_t i(0); i < this->size(); ++i)
        {
          this->at(i).first = other.at(i).first;
          this->at(i).second.convert(other.at(i).second);
        }
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        std::size_t my_bytes(0);

        for(size_t i(0); i < this->size(); ++i)
        {
          my_bytes += this->at(i).second.bytes();
        }

        return my_bytes;
      }

      Filter_& find_or_add(const String& name)
      {
        for(iterator it(BaseClass::begin()); it != BaseClass::end(); ++it)
        {
          if(it->first == name)
            return it->second;
        }
        // if we come out here, we have to add a new sub-filter
        this->push_back(std::make_pair(name, Filter_()));
        return this->back().second;
      }


      /** \copydoc UnitFilter::filter_rhs() */
      template<typename Vector_>
      void filter_rhs(Vector_& vector) const
      {
        for(const_iterator it(BaseClass::begin()); it != BaseClass::end(); ++it)
          it->second.filter_rhs(vector);
      }

      /** \copydoc UnitFilter::filter_sol() */
      template<typename Vector_>
      void filter_sol(Vector_& vector) const
      {
        for(const_iterator it(BaseClass::begin()); it != BaseClass::end(); ++it)
          it->second.filter_sol(vector);
      }

      /** \copydoc UnitFilter::filter_def() */
      template<typename Vector_>
      void filter_def(Vector_& vector) const
      {
        for(const_iterator it(BaseClass::begin()); it != BaseClass::end(); ++it)
          it->second.filter_def(vector);
      }

      /** \copydoc UnitFilter::filter_cor() */
      template<typename Vector_>
      void filter_cor(Vector_& vector) const
      {
        for(const_iterator it(BaseClass::begin()); it != BaseClass::end(); ++it)
          it->second.filter_cor(vector);
      }

      /** \copydoc UnitFilter::filter_mat() */
      template<typename Matrix_>
      void filter_mat(Matrix_& matrix) const
      {
        for(const_iterator it(BaseClass::begin()); it != BaseClass::end(); ++it)
          it->second.filter_mat(matrix);
      }
    }; // class FilterSequence

  } // namespace LAFEM
} // namespace FEAT
