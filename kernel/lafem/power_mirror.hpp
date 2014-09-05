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
        template<typename Algo_, typename SM_, typename BV_, typename PV_>
        static void gather_prim(const SM_& sm, BV_& bv, const PV_& pv, const Index bo)
        {
          sm.template gather_prim<Algo_>(bv, pv.first(), bo);
          PowerMirrorHelper<i_-1>::template gather_prim<Algo_>(sm, bv, pv.rest(), bo + sm.size());
        }
        template<typename Algo_, typename SM_, typename BV_, typename PV_>
        static void gather_dual(const SM_& sm, BV_& bv, const PV_& pv, const Index bo)
        {
          sm.template gather_dual<Algo_>(bv, pv.first(), bo);
          PowerMirrorHelper<i_-1>::template gather_dual<Algo_>(sm, bv, pv.rest(), bo + sm.size());
        }
        template<typename Algo_, typename SM_, typename BV_, typename PV_, typename AT_>
        static void gather_axpy_prim(const SM_& sm, BV_& bv, const PV_& pv, const AT_ a, const Index bo)
        {
          sm.template gather_axpy_prim<Algo_>(bv, pv.first(), a, bo);
          PowerMirrorHelper<i_-1>::template gather_axpy_prim<Algo_>(sm, bv, pv.rest(), a, bo + sm.size());
        }
        template<typename Algo_, typename SM_, typename BV_, typename PV_, typename AT_>
        static void gather_axpy_dual(const SM_& sm, BV_& bv, const PV_& pv, const AT_ a, const Index bo)
        {
          sm.template gather_axpy_dual<Algo_>(bv, pv.first(), a, bo);
          PowerMirrorHelper<i_-1>::template gather_axpy_dual<Algo_>(sm, bv, pv.rest(), a, bo + sm.size());
        }
        template<typename Algo_, typename SM_, typename PV_, typename BV_>
        static void scatter_prim(const SM_& sm, PV_& pv, const BV_& bv, const Index bo)
        {
          sm.template scatter_prim<Algo_>(pv.first(), bv, bo);
          PowerMirrorHelper<i_-1>::template scatter_prim<Algo_>(sm, pv.rest(), bv, bo + sm.size());
        }
        template<typename Algo_, typename SM_, typename PV_, typename BV_>
        static void scatter_dual(const SM_& sm, PV_& pv, const BV_& bv, const Index bo)
        {
          sm.template scatter_dual<Algo_>(pv.first(), bv, bo);
          PowerMirrorHelper<i_-1>::template scatter_dual<Algo_>(sm, pv.rest(), bv, bo + sm.size());
        }
        template<typename Algo_, typename SM_, typename PV_, typename BV_, typename AT_>
        static void scatter_axpy_prim(const SM_& sm, PV_& pv, const BV_& bv, const AT_ a, const Index bo)
        {
          sm.template scatter_axpy_prim<Algo_>(pv.first(), bv, a, bo);
          PowerMirrorHelper<i_-1>::template scatter_axpy_prim<Algo_>(sm, pv.rest(), bv, a, bo + sm.size());
        }
        template<typename Algo_, typename SM_, typename PV_, typename BV_, typename AT_>
        static void scatter_axpy_dual(const SM_& sm, PV_& pv, const BV_& bv, const AT_ a, const Index bo)
        {
          sm.template scatter_axpy_dual<Algo_>(pv.first(), bv, a, bo);
          PowerMirrorHelper<i_-1>::template scatter_axpy_dual<Algo_>(sm, pv.rest(), bv, a, bo + sm.size());
        }
      };
      template<>
      struct PowerMirrorHelper<Index(1)>
      {
        template<typename Algo_, typename SM_, typename BV_, typename PV_>
        static void gather_prim(const SM_& sm, BV_& bv, const PV_& pv, const Index bo)
        {
          sm.template gather_prim<Algo_>(bv, pv.first(), bo);
        }
        template<typename Algo_, typename SM_, typename BV_, typename PV_>
        static void gather_dual(const SM_& sm, BV_& bv, const PV_& pv, const Index bo)
        {
          sm.template gather_dual<Algo_>(bv, pv.first(), bo);
        }
        template<typename Algo_, typename SM_, typename BV_, typename PV_, typename AT_>
        static void gather_axpy_prim(const SM_& sm, BV_& bv, const PV_& pv, const AT_ a, const Index bo)
        {
          sm.template gather_axpy_prim<Algo_>(bv, pv.first(), a, bo);
        }
        template<typename Algo_, typename SM_, typename BV_, typename PV_, typename AT_>
        static void gather_axpy_dual(const SM_& sm, BV_& bv, const PV_& pv, const AT_ a, const Index bo)
        {
          sm.template gather_axpy_dual<Algo_>(bv, pv.first(), a, bo);
        }
        template<typename Algo_, typename SM_, typename PV_, typename BV_>
        static void scatter_prim(const SM_& sm, PV_& pv, const BV_& bv, const Index bo)
        {
          sm.template scatter_prim<Algo_>(pv.first(), bv, bo);
        }
        template<typename Algo_, typename SM_, typename PV_, typename BV_>
        static void scatter_dual(const SM_& sm, PV_& pv, const BV_& bv, const Index bo)
        {
          sm.template scatter_dual<Algo_>(pv.first(), bv, bo);
        }
        template<typename Algo_, typename SM_, typename PV_, typename BV_, typename AT_>
        static void scatter_axpy_prim(const SM_& sm, PV_& pv, const BV_& bv, const AT_ a, const Index bo)
        {
          sm.template scatter_axpy_prim<Algo_>(pv.first(), bv, a, bo);
        }
        template<typename Algo_, typename SM_, typename PV_, typename BV_, typename AT_>
        static void scatter_axpy_dual(const SM_& sm, PV_& pv, const BV_& bv, const AT_ a, const Index bo)
        {
          sm.template scatter_axpy_dual<Algo_>(pv.first(), bv, a, bo);
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
      /// sub-mirror index-type
      typedef typename SubMirrorType::IndexType IndexType;

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

      /**
       * \brief Clone operation
       *
       * \returns A deep copy of this power mirror.
       */
      PowerMirror clone() const
      {
        return PowerMirror(_sub_mirror.clone());
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

      /// Returns the total size of the mirror.
      Index size() const
      {
        return count_ * _sub_mirror.size();
      }

      /** \copydoc VectorMirror::template gather_prim() */
      template<typename Algo_, typename Tx_, typename Ix_, typename Tv_>
      void gather_prim(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::PowerVector<Tv_, count_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::template gather_prim<Algo_>(_sub_mirror, buffer, vector, buffer_offset);
      }

      /** \copydoc VectorMirror::template gather_axpy_prim() */
      template<typename Algo_, typename Tx_, typename Ix_, typename Tv_>
      void gather_axpy_prim(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::PowerVector<Tv_, count_>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::template gather_axpy_prim<Algo_>(_sub_mirror, buffer, vector, alpha, buffer_offset);
      }

      /** \copydoc VectorMirror::template scatter_prim() */
      template<typename Algo_, typename Tv_, typename Tx_, typename Ix_>
      void scatter_prim(
        LAFEM::PowerVector<Tv_, count_>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::template scatter_prim<Algo_>(_sub_mirror, vector, buffer, buffer_offset);
      }

      /** \copydoc VectorMirror::template scatter_axpy_prim() */
      template<typename Algo_, typename Tv_, typename Tx_, typename Ix_>
      void scatter_axpy_prim(
        LAFEM::PowerVector<Tv_, count_>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::template scatter_axpy_prim<Algo_>(_sub_mirror, vector, buffer, alpha, buffer_offset);
      }

      /** \copydoc VectorMirror::template gather_dual() */
      template<typename Algo_, typename Tx_, typename Ix_, typename Tv_>
      void gather_dual(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::PowerVector<Tv_, count_>& vector,
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::template gather_dual<Algo_>(_sub_mirror, buffer, vector, buffer_offset);
      }

      /** \copydoc VectorMirror::template gather_axpy_dual() */
      template<typename Algo_, typename Tx_, typename Ix_, typename Tv_>
      void gather_axpy_dual(
        LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const LAFEM::PowerVector<Tv_, count_>& vector,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::template gather_axpy_dual<Algo_>(_sub_mirror, buffer, vector, alpha, buffer_offset);
      }

      /** \copydoc VectorMirror::template scatter_dual() */
      template<typename Algo_, typename Tv_, typename Tx_, typename Ix_>
      void scatter_dual(
        LAFEM::PowerVector<Tv_, count_>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::template scatter_dual<Algo_>(_sub_mirror, vector, buffer, buffer_offset);
      }

      /** \copydoc VectorMirror::template scatter_axpy_dual() */
      template<typename Algo_, typename Tv_, typename Tx_, typename Ix_>
      void scatter_axpy_dual(
        LAFEM::PowerVector<Tv_, count_>& vector,
        const LAFEM::DenseVector<MemType, Tx_, Ix_>& buffer,
        const Tx_ alpha = Tx_(1),
        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::template scatter_axpy_dual<Algo_>(_sub_mirror, vector, buffer, alpha, buffer_offset);
      }
    }; // class PowerMirror<...>
  } // namespace LAFEM
} // namespace FEAST

#endif // KERNEL_LAFEM_POWER_MIRROR_HPP
