#pragma once
#ifndef KERNEL_LAFEM_POWER_MIRROR_HPP
#define KERNEL_LAFEM_POWER_MIRROR_HPP 1

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
        template<typename SM_, typename BV_, typename PV_, typename AT_>
        static void gather_axpy_prim(const SM_& sm, BV_& bv, const PV_& pv, const AT_ a, const Index bo)
        {
          sm.gather_axpy_prim(bv, pv.first(), a, bo);
          PowerMirrorHelper<i_-1>::gather_axpy_prim(sm, bv, pv.rest(), a, bo + sm.size());
        }
        template<typename SM_, typename BV_, typename PV_, typename AT_>
        static void gather_axpy_dual(const SM_& sm, BV_& bv, const PV_& pv, const AT_ a, const Index bo)
        {
          sm.gather_axpy_dual(bv, pv.first(), a, bo);
          PowerMirrorHelper<i_-1>::gather_axpy_dual(sm, bv, pv.rest(), a, bo + sm.size());
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
        template<typename SM_, typename PV_, typename BV_, typename AT_>
        static void scatter_axpy_prim(const SM_& sm, PV_& pv, const BV_& bv, const AT_ a, const Index bo)
        {
          sm.scatter_axpy_prim(pv.first(), bv, a, bo);
          PowerMirrorHelper<i_-1>::scatter_axpy_prim(sm, pv.rest(), bv, a, bo + sm.size());
        }
        template<typename SM_, typename PV_, typename BV_, typename AT_>
        static void scatter_axpy_dual(const SM_& sm, PV_& pv, const BV_& bv, const AT_ a, const Index bo)
        {
          sm.scatter_axpy_dual(pv.first(), bv, a, bo);
          PowerMirrorHelper<i_-1>::scatter_axpy_dual(sm, pv.rest(), bv, a, bo + sm.size());
        }
        template<typename SM_, typename PV_>
        static void create_vector(const SM_& sm, PV_& pv)
        {
          pv.first() = sm.create_vector();
          PowerMirrorHelper<i_-1>::create_vector(sm, pv.rest());
        }
      };
      template<>
      struct PowerMirrorHelper<1>
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
        template<typename SM_, typename BV_, typename PV_, typename AT_>
        static void gather_axpy_prim(const SM_& sm, BV_& bv, const PV_& pv, const AT_ a, const Index bo)
        {
          sm.gather_axpy_prim(bv, pv.first(), a, bo);
        }
        template<typename SM_, typename BV_, typename PV_, typename AT_>
        static void gather_axpy_dual(const SM_& sm, BV_& bv, const PV_& pv, const AT_ a, const Index bo)
        {
          sm.gather_axpy_dual(bv, pv.first(), a, bo);
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
        template<typename SM_, typename PV_, typename BV_, typename AT_>
        static void scatter_axpy_prim(const SM_& sm, PV_& pv, const BV_& bv, const AT_ a, const Index bo)
        {
          sm.scatter_axpy_prim(pv.first(), bv, a, bo);
        }
        template<typename SM_, typename PV_, typename BV_, typename AT_>
        static void scatter_axpy_dual(const SM_& sm, PV_& pv, const BV_& bv, const AT_ a, const Index bo)
        {
          sm.scatter_axpy_dual(pv.first(), bv, a, bo);
        }
        template<typename SM_, typename PV_>
        static void create_vector(const SM_& sm, PV_& pv)
        {
          pv.first() = sm.create_vector();
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
      int count_>
    class PowerMirror
    {
    public:
      // Note: the case = 1 is specialised below
      static_assert(count_ > 1, "invalid block size");

      /// sub-mirror type
      typedef SubMirror_ SubMirrorType;
      /// sub-mirror mem-type
      typedef typename SubMirrorType::MemType MemType;
      /// sub-mirror data-type
      typedef typename SubMirrorType::DataType DataType;
      /// sub-mirror index-type
      typedef typename SubMirrorType::IndexType IndexType;

      /// total number of sub-mirror blocks
      static constexpr int num_blocks = count_;

      /// corresponding vector type
      typedef PowerVector<typename SubMirrorType::VectorType, num_blocks> VectorType;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType, typename IT2_ = IndexType>
      using ContainerType = class PowerMirror<typename SubMirrorType::template ContainerType<Mem2_, DT2_, IT2_>, count_>;

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

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return _sub_mirror.bytes();
      }

      /// Returns the total size of the mirror.
      Index size() const
      {
        return Index(count_) * _sub_mirror.size();
      }

      /**
       * \brief Creates a new buffer vector.
       */
      DenseVector<Mem::Main, DataType, IndexType> create_buffer_vector() const
      {
        return DenseVector<Mem::Main, DataType, IndexType>(size(), true);
      }

      /**
       * \brief Creates a new (local) vector.
       */
      VectorType create_vector() const
      {
        VectorType vector;
        Intern::PowerMirrorHelper<count_>::create_vector(_sub_mirror, vector) ;
        return vector;
      }

      template<int i_>
      SubMirrorType& at()
      {
        static_assert((0 <= i_) && (i_ < count_), "invalid sub-mirror index");
        return _sub_mirror;
      }

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

      /** \copydoc VectorMirror::gather_prim() */
      template<typename Tx_, typename Ix_, typename Tv_>
      void gather_prim(
                       LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                       const LAFEM::PowerVector<Tv_, count_>& vector,
                       const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::gather_prim(_sub_mirror, buffer, vector, buffer_offset);
      }

      /** \copydoc VectorMirror::gather_axpy_prim() */
      template<typename Tx_, typename Ix_, typename Tv_>
      void gather_axpy_prim(
                            LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                            const LAFEM::PowerVector<Tv_, count_>& vector,
                            const Tx_ alpha = Tx_(1),
                            const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::gather_axpy_prim(_sub_mirror, buffer, vector, alpha, buffer_offset);
      }

      /** \copydoc VectorMirror::scatter_prim() */
      template<typename Tv_, typename Tx_, typename Ix_>
      void scatter_prim(
                        LAFEM::PowerVector<Tv_, count_>& vector,
                        const LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::scatter_prim(_sub_mirror, vector, buffer, buffer_offset);
      }

      /** \copydoc VectorMirror::scatter_axpy_prim() */
      template<typename Tv_, typename Tx_, typename Ix_>
      void scatter_axpy_prim(
                             LAFEM::PowerVector<Tv_, count_>& vector,
                             const LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                             const Tx_ alpha = Tx_(1),
                             const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::scatter_axpy_prim(_sub_mirror, vector, buffer, alpha, buffer_offset);
      }

      /** \copydoc VectorMirror::gather_dual() */
      template<typename Tx_, typename Ix_, typename Tv_>
      void gather_dual(
                       LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                       const LAFEM::PowerVector<Tv_, count_>& vector,
                       const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::gather_dual(_sub_mirror, buffer, vector, buffer_offset);
      }

      /** \copydoc VectorMirror::gather_axpy_dual() */
      template<typename Tx_, typename Ix_, typename Tv_>
      void gather_axpy_dual(
                            LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                            const LAFEM::PowerVector<Tv_, count_>& vector,
                            const Tx_ alpha = Tx_(1),
                            const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::gather_axpy_dual(_sub_mirror, buffer, vector, alpha, buffer_offset);
      }

      /** \copydoc VectorMirror::scatter_dual() */
      template<typename Tv_, typename Tx_, typename Ix_>
      void scatter_dual(
                        LAFEM::PowerVector<Tv_, count_>& vector,
                        const LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                        const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::scatter_dual(_sub_mirror, vector, buffer, buffer_offset);
      }

      /** \copydoc VectorMirror::scatter_axpy_dual() */
      template<typename Tv_, typename Tx_, typename Ix_>
      void scatter_axpy_dual(
                             LAFEM::PowerVector<Tv_, count_>& vector,
                             const LAFEM::DenseVector<Mem::Main, Tx_, Ix_>& buffer,
                             const Tx_ alpha = Tx_(1),
                             const Index buffer_offset = Index(0)) const
      {
        Intern::PowerMirrorHelper<count_>::scatter_axpy_dual(_sub_mirror, vector, buffer, alpha, buffer_offset);
      }
    }; // class PowerMirror<...>
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_POWER_MIRROR_HPP
