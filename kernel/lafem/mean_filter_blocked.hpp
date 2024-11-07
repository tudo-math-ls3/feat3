// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

// includes, FEAT
#include <kernel/lafem/dense_vector_blocked.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
    * \brief Mean Filter Blocked class template.
    *
    * \author Pia Ritter
    */
    template<typename DT_, typename IT_, int BlockSize_>
    class MeanFilterBlocked
    {
    public:
      /// vector-type typedef
      typedef DenseVectorBlocked<DT_, IT_, BlockSize_> VectorType;
      /// data-type typedef
      typedef typename VectorType::DataType DataType;
      /// index-type typedef
      typedef typename VectorType::IndexType IndexType;
      /// value-type typedef
      typedef typename VectorType::ValueType ValueType;

      /// Our 'base' class type
      template <typename DT2_ = DT_, typename IT2_ = IT_, int BlockSize2_ = BlockSize_>
      using FilterType = MeanFilterBlocked<DT2_, IT2_, BlockSize2_>;

      /// this typedef lets you create a filter with different Data and Index types
      template <typename DataType2_, typename IndexType2_, int BS2_>
      using FilterTypeByDI = FilterType<DataType2_, IndexType2_, BS2_>;

    protected:
      /// primal weighting vector
      VectorType _vec_prim;
      /// dual weighting vector
      VectorType _vec_dual;
      /// weight volume
      ValueType _volume;
      /// desired solution vector mean
      ValueType _sol_mean;

    public:
      // default CTOR
      MeanFilterBlocked() :
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
      explicit MeanFilterBlocked(VectorType && vec_prim, VectorType && vec_dual, ValueType sol_mean = ValueType(0)) :
        _vec_prim(std::forward<VectorType>(vec_prim)),
        _vec_dual(std::forward<VectorType>(vec_dual)),
        _volume(_vec_prim.dot_blocked(_vec_dual)),
        _sol_mean(sol_mean)
      {
        if(!_vec_prim.empty())
        {
          XASSERTM(_volume.norm_euclid_sqr()  > Math::eps<DataType>(), "domain volume must not be zero");
        }
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
      explicit MeanFilterBlocked(VectorType && vec_prim, VectorType && vec_dual, ValueType sol_mean, ValueType volume) :
        _vec_prim(std::forward<VectorType>(vec_prim)),
        _vec_dual(std::forward<VectorType>(vec_dual)),
        _volume(volume),
        _sol_mean(sol_mean)
      {
        if(!_vec_prim.empty())
        {
          XASSERTM(_volume.norm_euclid_sqr() > Math::eps<DataType>(), "domain volume must not be zero");
        }
      }

      /// move ctor
      MeanFilterBlocked(MeanFilterBlocked && other) :
        _vec_prim(std::move(other._vec_prim)),
        _vec_dual(std::move(other._vec_dual)),
        _volume(other._volume),
        _sol_mean(other._sol_mean)
      {
      }

      /// move-assign operator
      MeanFilterBlocked & operator=(MeanFilterBlocked && other)
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
      virtual ~MeanFilterBlocked()
      {
      }

      /**
      * \brief Creates a clone of itself
      * \warning _volume and _sol_mean will always be a deep copy
      */
      MeanFilterBlocked clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        return MeanFilterBlocked(_vec_prim.clone(clone_mode), _vec_dual.clone(clone_mode), _sol_mean, _volume);
      }

      /**
      * \brief Clones data from another MeanFilterBlocked
      * \warning _volume and _sol_mean will always be a deep copy
      */
      void clone(const MeanFilterBlocked & other, CloneMode clone_mode = CloneMode::Deep)
      {
        _vec_prim.clone(other.get_vec_prim(), clone_mode);
        _vec_dual.clone(other.get_vec_dual(), clone_mode);
        _volume = other._volume;
        _sol_mean = other._sol_mean;
      }

      /// Conversion method
      template<typename DT2_, typename IT2_, int BlockSize2_>
      void convert(const MeanFilterBlocked<DT2_, IT2_, BlockSize2_>& other)
      {
        _vec_prim.convert(other.get_vec_prim());
        _vec_dual.convert(other.get_vec_dual());
        _volume = DataType(other.get_volume());
        _sol_mean = DataType(other.get_sol_mean());
        if(!_vec_prim.empty())
        {
          XASSERTM(_volume.norm_euclid_sqr() > Math::eps<DataType>(), "domain volume must not be zero");
        }
      }

      /// Clears this filter
      void clear()
      {
        _vec_prim.clear();
        _vec_dual.clear();
        _volume.clear();
        _sol_mean.clear();
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

      ValueType get_volume() const
      {
        return _volume;
      }

      ValueType get_sol_mean() const
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
        {
          ValueType tmp(vector.dot_blocked(_vec_prim));
          for(int i(0); i<ValueType::n; ++i)
          {
            tmp(i) /= -_volume(i);
          }
          vector.axpy_blocked(_vec_dual, tmp);
        }
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
        {
          ValueType tmp(vector.dot_blocked(_vec_dual));
          for(int i(0); i<ValueType::n; ++i)
          {
            tmp(i)= _sol_mean(i) - tmp(i) * DT_(DT_(1)/_volume(i));
          }
          vector.axpy_blocked(_vec_prim, tmp);
        }
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
        {
          ValueType tmp(vector.dot_blocked(_vec_dual));
          for(int i(0); i<ValueType::n; ++i)
          {
            tmp(i) /= -_volume(i);
          }
          vector.axpy_blocked(_vec_prim, tmp);
        }
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
    }; // class MeanFilterBlocked<...>
  } // namespace LAFEM
} // namespace FEAT
