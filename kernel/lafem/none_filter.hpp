// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
      typename DataType_,
      typename IndexType_ = Index>
    class NoneFilter
    {
    public:
      /// data-type typedef
      typedef DataType_ DataType;
      /// index-type typedef
      typedef IndexType_ IndexType;

      /// our supported vector type
      typedef DenseVector<DataType, IndexType> VectorType;

      /// Our 'base' class type
      template <typename DT2_ = DataType_, typename IT2_ = IndexType_>
      using FilterType = NoneFilter<DT2_, IT2_>;

      /// this typedef lets you create a filter with different Data and Index types
      template <typename DataType2_, typename IndexType2_>
      using FilterTypeByDI = FilterType<DataType2_, IndexType2_>;

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

      template<typename DT2_, typename IT2_>
      void convert(const NoneFilter<DT2_, IT2_>&)
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

      /**
       * \brief Applies the filter onto a system matrix.
       *
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       */
      template<typename MT_>
      void filter_mat(MT_& DOXY(matrix)) const
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
      typename DataType_,
      typename IndexType_,
      int BlockSize_
    >
    class NoneFilterBlocked
    {
    public:
      /// data-type typedef
      typedef DataType_ DataType;
      /// index-type typedef
      typedef IndexType_ IndexType;

      static constexpr int BlockSize = BlockSize_;

      /// our supported vector type
      typedef DenseVectorBlocked<DataType, IndexType, BlockSize> VectorType;

      /// Our 'base' class type
      template <typename DT2_ = DataType_, typename IT2_ = IndexType_, int BS_ = BlockSize_>
      using FilterType = NoneFilterBlocked<DT2_, IT2_, BS_>;

      /// this typedef lets you create a filter with different Data and Index types
      template <typename DataType2_, typename IndexType2_, int BlockSize2_>
      using FilterTypeByDI = FilterType<DataType2_, IndexType2_, BlockSize2_>;

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

      template<typename DT2_, typename IT2_>
      void convert(const NoneFilterBlocked<DT2_, IT2_, BlockSize>&)
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

      /**
       * \brief Applies the filter onto a system matrix.
       *
       * \param[in,out] matrix
       * A reference to the matrix to be filtered.
       */
      template<typename MT_>
      void filter_mat(MT_& DOXY(matrix)) const
      {
      }
    }; // class NoneFilterBlocked<...>
  } // namespace LAFEM
} // namespace FEAT

#endif // KERNEL_LAFEM_NONE_FILTER_HPP
