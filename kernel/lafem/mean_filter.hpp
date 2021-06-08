// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_LAFEM_MEAN_FILTER_HPP
#define KERNEL_LAFEM_MEAN_FILTER_HPP 1

// includes, FEAT
#include <kernel/lafem/dense_vector.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Mean Filter class template.
     *
     * \author Peter Zajac
     */
    template<typename Mem_, typename DT_, typename IT_>
    class MeanFilter
    {
    public:
      /// vector-type typedef
      typedef DenseVector<Mem_, DT_, IT_> VectorType;
      /// mem-type typedef
      typedef typename VectorType::MemType MemType;
      /// data-type typedef
      typedef typename VectorType::DataType DataType;
      /// index-type typedef
      typedef typename VectorType::IndexType IndexType;

      /// Our 'base' class type
      template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
      using FilterType = MeanFilter<Mem2_, DT2_, IT2_>;

      /// this typedef lets you create a filter with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using FilterTypeByMDI = FilterType<Mem2_, DataType2_, IndexType2_>;

    protected:
      /// primal weighting vector
      VectorType _vec_prim;
      /// dual weighting vector
      VectorType _vec_dual;
      /// weight volume
      DataType _volume;
      /// desired solution vector mean
      DataType _sol_mean;

    public:
      // default CTOR
      MeanFilter() :
        _vec_prim(),
        _vec_dual(),
        _volume(),
        _sol_mean()
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] vec_prim, vec_dual
       * The primal-dual weighting vector pair for the mean filter.
       *
       * \param[in] sol_mean
       * The desired integral mean of the solution vector. Defaults to 0.
       */
      explicit MeanFilter(VectorType && vec_prim, VectorType && vec_dual,
        DataType sol_mean = DataType(0)) :
        _vec_prim(std::forward<VectorType>(vec_prim)),
        _vec_dual(std::forward<VectorType>(vec_dual)),
        _volume(_vec_prim.dot(_vec_dual)),
        _sol_mean(sol_mean)
      {
        XASSERTM(_volume > Math::eps<DataType>(), "domain volume must not be zero");
      }

      /**
       * \brief Constructor
       *
       * \param[in] vec_prim, vec_dual
       * The primal-dual weighting vector pair for the mean filter.
       *
       * \param[in] sol_mean
       * The desired integral mean of the solution vector.
       *
       * \param[in] volume
       * The domain volume. This is simply the dot-product of the two vectors.
       */
      explicit MeanFilter(VectorType && vec_prim, VectorType && vec_dual, DataType sol_mean, DataType volume) :
        _vec_prim(std::forward<VectorType>(vec_prim)),
        _vec_dual(std::forward<VectorType>(vec_dual)),
        _volume(volume),
        _sol_mean(sol_mean)
      {
        XASSERTM(_volume > Math::eps<DataType>(), "domain volume must not be zero");
      }

      /// move ctor
      MeanFilter(MeanFilter && other) :
        _vec_prim(std::move(other._vec_prim)),
        _vec_dual(std::move(other._vec_dual)),
        _volume(other._volume),
        _sol_mean(other._sol_mean)
      {
      }

      /// move-assign operator
      MeanFilter & operator=(MeanFilter && other)
      {
        if(this != &other)
        {
          _vec_prim = std::forward<VectorType>(other._vec_prim);
          _vec_dual = std::forward<VectorType>(other._vec_dual);
          _volume = other._volume;
          _sol_mean = other._sol_mean;
        }
        return *this;
      }

      /// virtual destructor
      virtual ~MeanFilter()
      {
      }

      /**
       * \brief Creates a clone of itself
       * \warning _volume and _sol_mean will always be a deep copy
       */
      MeanFilter clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        return MeanFilter(_vec_prim.clone(clone_mode), _vec_dual.clone(clone_mode), _sol_mean, _volume);
      }

      /**
       * \brief Clones data from another MeanFilter
       * \warning _volume will always be a deep copy
       */
      void clone(const MeanFilter & other, CloneMode clone_mode = CloneMode::Deep)
      {
        _vec_prim.clone(other.get_vec_prim(), clone_mode);
        _vec_dual.clone(other.get_vec_dual(), clone_mode);
        _volume = other._volume;
        _sol_mean = other._sol_mean;
      }

      /// Conversion method
      template<typename Mem2_, typename DT2_, typename IT2_>
      void convert(const MeanFilter<Mem2_, DT2_, IT2_>& other)
      {
        _vec_prim.convert(other.get_vec_prim());
        _vec_dual.convert(other.get_vec_dual());
        _volume = DataType(other.get_volume());
        _sol_mean = DataType(other.get_sol_mean());
      }

      /// \cond internal
      VectorType & get_vec_prim()
      {
        return _vec_prim;
      }

      const VectorType & get_vec_prim() const
      {
        return _vec_prim;
      }

      VectorType & get_vec_dual()
      {
        return _vec_dual;
      }

      const VectorType & get_vec_dual() const
      {
        return _vec_dual;
      }

      DataType get_volume() const
      {
        return _volume;
      }

      DataType get_sol_mean() const
      {
        return _sol_mean;
      }
      /// \endcond

      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      void filter_rhs(VectorType& vector) const
      {
        // subtract to dual integral mean zero
        if(!_vec_prim.empty())
          vector.axpy(_vec_dual, vector, -vector.dot(_vec_prim) / _volume);
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      void filter_sol(VectorType& vector) const
      {
        // given: v=vec_prim, w=vec_dual, vol=v*w, input: x
        // note that int(v) = v*w and int(x) = x*w
        // we want:
        //      1/vol * int (x + alpha*v) = sol_mean
        // <==>     int(x) + alpha*int(v) = sol_mean*vol
        // <==>           x*w + alpha*v*w = sol_mean*v*w
        // <==>                 alpha*v*w = sol_mean*v*w - x*w
        // <==>                     alpha = sol_mean - (x*w/v*w)
        if(!_vec_prim.empty())
          vector.axpy(_vec_prim, vector, _sol_mean - vector.dot(_vec_dual) / _volume);
      }

      /**
       * \brief Applies the filter onto a defect vector.
       *
       * \param[in,out] vector
       * A reference to the defect vector to be filtered.
       */
      void filter_def(VectorType& vector) const
      {
        // same as rhs
        filter_rhs(vector);
      }

      /**
       * \brief Applies the filter onto a correction vector.
       *
       * \param[in,out] vector
       * A reference to the correction vector to be filtered.
       */
      void filter_cor(VectorType& vector) const
      {
        // subtract to primal integral mean zero
        if(!_vec_prim.empty())
          vector.axpy(_vec_prim, vector, -vector.dot(_vec_dual) / _volume);
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
        // nothing to do here
      }
    }; // class MeanFilter<...>
  } // namespace LAFEM
} // namespace FEAT


#endif // KERNEL_LAFEM_MEAN_FILTER_HPP
