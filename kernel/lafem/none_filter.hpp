#pragma once
#ifndef KERNEL_LAFEM_NONE_FILTER_HPP
#define KERNEL_LAFEM_NONE_FILTER_HPP 1

// includes, FEAT
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief None Filter class template.
     *
     * This class implements a filter that does nothing.
     *
     * \note Not to be confused with <em>do-nothing boundary conditions</em>.
     *
     * \author Peter Zajac
     */
    template<
      typename MemType_,
      typename DataType_,
      typename IndexType_ = Index>
    class NoneFilter
    {
    public:
      /// mem-type typedef
      typedef MemType_ MemType;
      /// data-type typedef
      typedef DataType_ DataType;
      /// index-type typedef
      typedef IndexType_ IndexType;

      /// our supported vector type
      typedef DenseVector<MemType, DataType, IndexType> VectorType;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType_, typename IT2_ = IndexType_>
      using FilterType = class NoneFilter<Mem2_, DT2_, IT2_>;

      /// this typedef lets you create a filter with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using FilterTypeByMDI = FilterType<Mem2_, DataType2_, IndexType2_>;

      static constexpr bool is_global = false;
      static constexpr bool is_local = true;

      /// \brief Creates a (empty) clone of itself
      NoneFilter clone(CloneMode /*clone_mode*/ = CloneMode::Deep) const
      {
        return NoneFilter();
      }

      /// \brief Clones data from another NoneFilter
      void clone(const NoneFilter & /*other*/, CloneMode /*clone_mode*/ = CloneMode::Deep)
      {
        // do nothing
      }

      template<typename MT2_, typename DT2_, typename IT2_>
      void convert(const NoneFilter<MT2_, DT2_, IT2_>&)
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return 0;
      }

      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      void filter_rhs(VectorType& DOXY(vector)) const
      {
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      void filter_sol(VectorType& DOXY(vector)) const
      {
      }

      /**
       * \brief Applies the filter onto a defect vector.
       *
       * \param[in,out] vector
       * A reference to the defect vector to be filtered.
       */
      void filter_def(VectorType& DOXY(vector)) const
      {
      }

      /**
       * \brief Applies the filter onto a correction vector.
       *
       * \param[in,out] vector
       * A reference to the correction vector to be filtered.
       */
      void filter_cor(VectorType& DOXY(vector)) const
      {
      }
    }; // class NoneFilter<...>

    /**
     * \brief Blocked None Filter class template.
     *
     * This class implements a filter that does nothing.
     *
     * \note Not to be confused with <em>do-nothing boundary conditions</em>.
     *
     * \author Peter Zajac
     */
    template<
      typename MemType_,
      typename DataType_,
      typename IndexType_,
      int BlockSize_
    >
    class NoneFilterBlocked
    {
    public:
      /// mem-type typedef
      typedef MemType_ MemType;
      /// data-type typedef
      typedef DataType_ DataType;
      /// index-type typedef
      typedef IndexType_ IndexType;

      static constexpr int BlockSize = BlockSize_;

      /// our supported vector type
      typedef DenseVectorBlocked<MemType, DataType, IndexType, BlockSize> VectorType;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DataType_, typename IT2_ = IndexType_, int BS_ = BlockSize_>
      using FilterType = class NoneFilter<Mem2_, DT2_, IT2_>;

      /// this typedef lets you create a filter with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_, int BlockSize2_>
      using FilterTypeByMDI = FilterType<Mem2_, DataType2_, IndexType2_>;

      /// \brief Creates a (empty) clone of itself
      NoneFilterBlocked clone(CloneMode /*clone_mode*/ = CloneMode::Deep) const
      {
        return NoneFilterBlocked();
      }

      /// \brief Clones data from another NoneFilterBlocked
      void clone(const NoneFilterBlocked & /*other*/, CloneMode /*clone_mode*/ = CloneMode::Deep)
      {
        // do nothing
      }

      template<typename MT2_, typename DT2_, typename IT2_>
      void convert(const NoneFilterBlocked<MT2_, DT2_, IT2_, BlockSize>&)
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return 0;
      }

      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      void filter_rhs(VectorType& DOXY(vector)) const
      {
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      void filter_sol(VectorType& DOXY(vector)) const
      {
      }

      /**
       * \brief Applies the filter onto a defect vector.
       *
       * \param[in,out] vector
       * A reference to the defect vector to be filtered.
       */
      void filter_def(VectorType& DOXY(vector)) const
      {
      }

      /**
       * \brief Applies the filter onto a correction vector.
       *
       * \param[in,out] vector
       * A reference to the correction vector to be filtered.
       */
      void filter_cor(VectorType& DOXY(vector)) const
      {
      }
    }; // class NoneFilterBlocked<...>
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_NONE_FILTER_HPP
