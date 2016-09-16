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
      using FilterType = class MeanFilter<Mem2_, DT2_, IT2_>;

    protected:
      /// primal weighting vector
      VectorType _vec_prim;
      /// dual weighting vector
      VectorType _vec_dual;
      /// weight volume
      DataType _volume;

    public:
      // default CTOR
      MeanFilter()
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] vec_prim, vec_dual
       * The primal-dual weighting vector pair for the mean filter.
       */
      explicit MeanFilter(VectorType && vec_prim, VectorType && vec_dual) :
        _vec_prim(std::forward<VectorType>(vec_prim)),
        _vec_dual(std::forward<VectorType>(vec_dual)),
        _volume(_vec_prim.dot(_vec_dual))
      {
        XASSERTM(_volume > Math::eps<DataType>(), "domain volume must not be zero");
      }

      /**
       * \brief Constructor
       *
       * \param[in] vec_prim, vec_dual
       * The primal-dual weighting vector pair for the mean filter.
       *
       * \param[in] volume
       * The domain volume. This is simply the dot-product of the two vectors.
       */
      explicit MeanFilter(VectorType && vec_prim, VectorType && vec_dual, DataType volume) :
        _vec_prim(std::forward<VectorType>(vec_prim)),
        _vec_dual(std::forward<VectorType>(vec_dual)),
        _volume(volume)
      {
        XASSERTM(_volume > Math::eps<DataType>(), "domain volume must not be zero");
      }

      /// move ctor
      MeanFilter(MeanFilter && other) :
        _vec_prim(std::move(other._vec_prim)),
        _vec_dual(std::move(other._vec_dual)),
        _volume(other._volume)
      {
      }

      /// move-assign operator
      MeanFilter & operator=(MeanFilter && other)
      {
        if(this != &other)
        {
          _vec_prim = std::move(other._vec_prim);
          _vec_dual = std::move(other._vec_dual);
          _volume = other._volume;
        }
        return *this;
      }

      /// virtual destructor
      virtual ~MeanFilter()
      {
      }

      /**
       * \brief Creates a clone of itself
       * \warning _volume will always be a deep copy
       */
      MeanFilter clone(CloneMode clone_mode = CloneMode::Deep) const
      {
        return MeanFilter(_vec_prim.clone(clone_mode), _vec_dual.clone(clone_mode), _volume);
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
      }


      /// Conversion method
      template<typename Mem2_, typename DT2_, typename IT2_>
      void convert(const MeanFilter<Mem2_, DT2_, IT2_>& other)
      {
        _vec_prim.convert(other.get_vec_prim());
        _vec_dual.convert(other.get_vec_dual());
        _volume = DataType(other.get_volume());
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
      /// \endcond

      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      void filter_rhs(VectorType& vector) const
      {
        // compute dual integral
        DataType integ = vector.dot(_vec_prim);
        // subtract mean
        vector.axpy(_vec_dual, vector, -integ / _volume);
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      void filter_sol(VectorType& vector) const
      {
        // compute primal integral
        DataType integ = vector.dot(_vec_dual);
        // subtract mean
        vector.axpy(_vec_prim, vector, -integ / _volume);
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
        // same as sol
        filter_sol(vector);
      }
    }; // class MeanFilter<...>
  } // namespace LAFEM
} // namespace FEAT


#endif // KERNEL_LAFEM_MEAN_FILTER_HPP
