// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /// \cond internal
    namespace Intern
    {
      // PowerMirror helper class: this implements the actual recursion for the gather/scatter methods
      template<int i_>
      struct PowerMirrorHelper
      {
        template<typename SM_, typename BV_, typename PV_>
        static void gather(const SM_& sm, BV_& bv, const PV_& pv, const Index bo)
        {
          sm.gather(bv, pv.first(), bo);
          PowerMirrorHelper<i_-1>::gather(sm, bv, pv.rest(), bo + sm.buffer_size(pv.first()));
        }
        template<typename SM_, typename PV_, typename BV_, typename AT_>
        static void scatter_axpy(const SM_& sm, PV_& pv, const BV_& bv, const AT_ a, const Index bo)
        {
          sm.scatter_axpy(pv.first(), bv, a, bo);
          PowerMirrorHelper<i_-1>::scatter_axpy(sm, pv.rest(), bv, a, bo + sm.buffer_size(pv.first()));
        }
        template<typename SM_, typename PV_>
        static Index buffer_size(const SM_& sm, PV_& pv)
        {
          return sm.buffer_size(pv.first()) + PowerMirrorHelper<i_-1>::buffer_size(sm, pv.rest());
        }
        template<Perspective perspective_, typename SM_, typename Pv_>
        static Index mask_scatter(const SM_& sm, const Pv_& pv, std::vector<int>& mask, const int value, const Index offset)
        {
          Index nf = sm.template mask_scatter<perspective_>(pv.first(), mask, value, offset);
          Index nr = PowerMirrorHelper<i_-1>::template mask_scatter<perspective_>(sm, pv.rest(), mask, value, offset + nf);
          return nf + nr;
        }
      };

      template<>
      struct PowerMirrorHelper<1>
      {
        template<typename SM_, typename BV_, typename PV_>
        static void gather(const SM_& sm, BV_& bv, const PV_& pv, const Index bo)
        {
          sm.gather(bv, pv.first(), bo);
        }
        template<typename SM_, typename PV_, typename BV_, typename AT_>
        static void scatter_axpy(const SM_& sm, PV_& pv, const BV_& bv, const AT_ a, const Index bo)
        {
          sm.scatter_axpy(pv.first(), bv, a, bo);
        }
        template<typename SM_, typename PV_>
        static Index buffer_size(const SM_& sm, PV_& pv)
        {
          return sm.buffer_size(pv.first());
        }
        template<Perspective perspective_, typename SM_, typename Pv_>
        static Index mask_scatter(const SM_& sm, const Pv_& pv, std::vector<int>& mask, const int value, const Index offset)
        {
          return sm.template mask_scatter<perspective_>(pv.first(), mask, value, offset);
        }
      };
    } // namespace Intern
    /// \endcond

    /**
     * \brief PowerVector meta-mirror class template
     *
     * This class template implements a vector-mirror which is capable of mirroring
     * PowerVector objects.
     *
     * \attention
     * This class template uses the same mirror for all blocks of a PowerVector object!
     *
     * \tparam SubMirror_
     * The type of the sub-mirror.
     *
     * \tparam count_
     * The number of sub-mirror blocks.
     *
     * \author Peter Zajac
     */
    template<typename SubMirror_, int count_>
    class PowerMirror
    {
    public:
      static_assert(count_ > 0, "invalid block size");

      /// sub-mirror type
      typedef SubMirror_ SubMirrorType;
      /// sub-mirror data-type
      typedef typename SubMirrorType::DataType DataType;
      /// sub-mirror index-type
      typedef typename SubMirrorType::IndexType IndexType;

      /// total number of sub-mirror blocks
      static constexpr int num_blocks = count_;

      /// Our 'base' class type
      template <typename DT2_ = DataType, typename IT2_ = IndexType>
      using MirrorType = PowerMirror<typename SubMirrorType::template MirrorType<DT2_, IT2_>, count_>;

      /// this typedef lets you create a mirror with new Data and Index types
      template <typename DataType2_, typename IndexType2_>
      using MirrorTypeByDI = MirrorType<DataType2_, IndexType2_>;

      /// the one and only sub-mirror object
      SubMirrorType _sub_mirror;

      /// default ctor
      PowerMirror()
      {
      }

      /// sub-mirror move-ctor
      explicit PowerMirror(SubMirrorType&& sub_mirror) :
        _sub_mirror(std::move(sub_mirror))
      {
      }

      /// move-ctor
      PowerMirror(PowerMirror&& other) :
        _sub_mirror(std::move(other._sub_mirror))
      {
      }

      /// move-assign operator
      PowerMirror& operator=(PowerMirror&& other)
      {
        if(this != &other)
        {
          _sub_mirror = std::move(other._sub_mirror);
        }
        return *this;
      }

      /**
       * \brief Clone operation
       *
       * \returns A deep copy of this power mirror.
       */
      PowerMirror clone(CloneMode clone_mode = CloneMode::Weak) const
      {
        return PowerMirror(_sub_mirror.clone(clone_mode));
      }

      /**
       * \brief Conversion method
       *
       * \param[in] other
       * The source mirror.
       */
      template<typename SubMirror2_>
      void convert(const PowerMirror<SubMirror2_, count_>& other)
      {
        this->_sub_mirror.convert(other._sub_mirror);
      }

      /**
       * \brief Creates and returns an empty mirror
       *
       * \param[in] tmpl_vec
       * A \transient reference to a template vector to get the correct size from
       *
       * \returns The new empty mirror
       */
      template<typename SubVector_>
      static PowerMirror make_empty(const PowerVector<SubVector_, count_>& tmpl_vec)
      {
        return PowerMirror(SubMirrorType::make_empty(tmpl_vec.first()));
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _sub_mirror.bytes();
      }

      /**
       * \brief Checks whether the mirror is empty.
       *
       * \returns \c true, if there are no indices in the mirror, otherwise \c false.
       */
      bool empty() const
      {
        return this->_sub_mirror.empty();
      }

      /**
       * \brief Computes the required buffer size for a PowerVector.
       *
       * \param[in] vector
       * The vector whose buffer size is to be computed.
       */
      template<typename SubVector_>
      Index buffer_size(const PowerVector<SubVector_, count_>& vector) const
      {
        return Intern::PowerMirrorHelper<count_>::buffer_size(_sub_mirror, vector);
      }

      /**
       * \brief Creates a new buffer vector for a PowerVector.
       *
       * \param[in] vector
       * The vector for which the buffer is to be created.
       */
      template<typename SubVector_>
      DenseVector<DataType, IndexType> create_buffer(const PowerVector<SubVector_, count_>& vector) const
      {
        return DenseVector<DataType, IndexType>(buffer_size(vector));
      }

      /**
       * \brief Returns a sub-mirror block.
       *
       * \tparam i_
       * The index of the sub-mirror block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-mirror at position \p i_.
       */
      template<int i_>
      SubMirrorType& at()
      {
        static_assert((0 <= i_) && (i_ < count_), "invalid sub-mirror index");
        return _sub_mirror;
      }

      /** \copydoc at() */
      template<int i_>
      const SubMirrorType& at() const
      {
        static_assert((0 <= i_) && (i_ < count_), "invalid sub-mirror index");
        return _sub_mirror;
      }

      /**
       * \brief Returns a sub-mirror block.
       *
       * \param[in] i
       * The index of the sub-mirror block that is to be returned.
       *
       * \returns
       * A (const) reference to the sub-mirror at position \p i.
       */
      SubMirrorType& get(int i)
      {
        XASSERTM((0 <= i) && (i < count_), "invalid sub-mirror index");
        return _sub_mirror;
      }

      /** \copydoc get() */
      const SubMirrorType& get(int i) const
      {
        XASSERTM((0 <= i) && (i < count_), "invalid sub-mirror index");
        return _sub_mirror;
      }

      /** \copydoc VectorMirror::gather() */
      template<typename Tv_>
      void gather(
        LAFEM::DenseVector<DataType, IndexType>& buffer,
        const LAFEM::PowerVector<Tv_, count_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::gather(_sub_mirror, buffer, vector, buffer_offset);
      }

      /** \copydoc VectorMirror::scatter_axpy() */
      template<typename Tv_>
      void scatter_axpy(
        LAFEM::PowerVector<Tv_, count_>& vector,
        const LAFEM::DenseVector<DataType, IndexType>& buffer,
        const DataType alpha = DataType(1),
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::scatter_axpy(_sub_mirror, vector, buffer, alpha, buffer_offset);
      }

      /** \copydoc VectorMirror::mask_scatter() */
      template<Perspective perspective_, typename Tv_>
      Index mask_scatter(const LAFEM::PowerVector<Tv_, count_>& vector, std::vector<int>& mask,
        const int value, const Index offset = Index(0)) const
      {
        return Intern::PowerMirrorHelper<count_>::template mask_scatter<perspective_>(_sub_mirror, vector, mask, value, offset);
      }

    }; // class PowerMirror<...>
  } // namespace LAFEM
} // namespace FEAT
