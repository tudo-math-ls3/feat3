#pragma once
#ifndef KERNEL_LAFEM_POWER_MIRROR_HPP
#define KERNEL_LAFEM_POWER_MIRROR_HPP 1

#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /// \cond internal
    namespace Intern
    {
      // PowerMirror helper class: this implements the actual recursion for the gather/scatter methods
      template<Index i_>
      struct PowerMirrorHelper
      {
        template<typename SM_, typename BV_, typename PV_>
        static void gather_prim(const SM_& sm, BV_& bv, const PV_& pv, const Index bo)
        {
          sm.gather_prim(bv, pv.first(), bo);
          PowerMirrorHelper<i_-1>::gather_prim(sm, bv, pv.rest(), bo + sm.size());
        }
        template<typename SM_, typename BV_, typename PV_>
        static void gather_dual(const SM_& sm, BV_& bv, const PV_& pv, const Index bo)
        {
          sm.gather_dual(bv, pv.first(), bo);
          PowerMirrorHelper<i_-1>::gather_dual(sm, bv, pv.rest(), bo + sm.size());
        }
        template<typename SM_, typename PV_, typename BV_>
        static void scatter_prim(const SM_& sm, PV_& pv, const BV_& bv, const Index bo)
        {
          sm.scatter_prim(pv.first(), bv, bo);
          PowerMirrorHelper<i_-1>::scatter_prim(sm, pv.rest(), bv, bo + sm.size());
        }
        template<typename SM_, typename PV_, typename BV_>
        static void scatter_dual(const SM_& sm, PV_& pv, const BV_& bv, const Index bo)
        {
          sm.scatter_dual(pv.first(), bv, bo);
          PowerMirrorHelper<i_-1>::scatter_dual(sm, pv.rest(), bv, bo + sm.size());
        }
      };
      template<>
      struct PowerMirrorHelper<Index(1)>
      {
        template<typename SM_, typename BV_, typename PV_>
        static void gather_prim(const SM_& sm, BV_& bv, const PV_& pv, const Index bo)
        {
          sm.gather_prim(bv, pv.first(), bo);
        }
        template<typename SM_, typename BV_, typename PV_>
        static void gather_dual(const SM_& sm, BV_& bv, const PV_& pv, const Index bo)
        {
          sm.gather_dual(bv, pv.first(), bo);
        }
        template<typename SM_, typename PV_, typename BV_>
        static void scatter_prim(const SM_& sm, PV_& pv, const BV_& bv, const Index bo)
        {
          sm.scatter_prim(pv.first(), bv, bo);
        }
        template<typename SM_, typename PV_, typename BV_>
        static void scatter_dual(const SM_& sm, PV_& pv, const BV_& bv, const Index bo)
        {
          sm.scatter_dual(pv.first(), bv, bo);
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
    template<
      typename SubMirror_,
      Index count_>
    class PowerMirror
    {
    public:
      /// sub-mirror type
      typedef SubMirror_ SubMirrorType;
      /// sub-mirror mem-type
      typedef typename SubMirrorType::MemType MemType;
      /// sub-mirror data-type
      typedef typename SubMirrorType::DataType DataType;

      /// dummy enum
      enum
      {
        /// total number of sub-mirror blocks
        num_blocks = count_
      };

    protected:
      /// the one and only sub-mirror object
      SubMirrorType _sub_mirror;

    public:
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
          _sub_mirror(std::move(other._sub_mirror));
        }
        return *this;
      }

      /// Returns the total size of the mirror.
      Index size() const
      {
        return count_ * _sub_mirror.size();
      }

      /** \copydoc VectorMirror::gather_prim() */
      template<typename Tx_, typename Tv_>
      void gather_prim(
        LAFEM::DenseVector<MemType, Tx_>& buffer,
        const LAFEM::PowerVector<Tv_, count_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::gather_prim(_sub_mirror, buffer, vector, buffer_offset);
      }

      /** \copydoc VectorMirror::scatter_prim() */
      template<typename Tv_, typename Tx_>
      void scatter_prim(
        LAFEM::PowerVector<Tv_, count_>& vector,
        const LAFEM::DenseVector<MemType, Tx_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::scatter_prim(_sub_mirror, vector, buffer, buffer_offset);
      }

      /** \copydoc VectorMirror::gather_dual() */
      template<typename Tx_, typename Tv_>
      void gather_dual(
        LAFEM::DenseVector<MemType, Tx_>& buffer,
        const LAFEM::PowerVector<Tv_, count_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::gather_dual(_sub_mirror, buffer, vector, buffer_offset);
      }

      /** \copydoc VectorMirror::scatter_dual() */
      template<typename Tv_, typename Tx_>
      void scatter_dual(
        LAFEM::PowerVector<Tv_, count_>& vector,
        const LAFEM::DenseVector<MemType, Tx_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::scatter_dual(_sub_mirror, vector, buffer, buffer_offset);
      }
    }; // class PowerMirror<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_POWER_MIRROR_HPP
