#pragma once
#ifndef KERNEL_GLOBAL_MEAN_FILTER_HPP
#define KERNEL_GLOBAL_MEAN_FILTER_HPP 1

// includes, FEAT
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/global/synch_scal.hpp>

namespace FEAT
{
  namespace Global
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
      typedef LAFEM::DenseVector<Mem_, DT_, IT_> VectorType;
      /// mem-type typedef
      typedef typename VectorType::MemType MemType;
      /// data-type typedef
      typedef typename VectorType::DataType DataType;
      /// index-type typedef
      typedef typename VectorType::IndexType IndexType;

    protected:
      /// primal weighting vector
      VectorType _vec_prim;
      /// dual weighting vector
      VectorType _vec_dual;
      /// frequency vector
      VectorType _vec_freq;
      /// weight volume
      DataType _volume;

    public:
      // default CTOR
      MeanFilter() :
        _vec_prim(),
        _vec_dual(),
        _vec_freq(),
        _volume(DataType(0))
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] vec_prim, vec_dual
       * The primal-dual weighting vector pair for the mean filter.
       *
       * \param[in] vec_freq
       * The frequency vector for the dot-product.
       */
      explicit MeanFilter(VectorType&& vec_prim, VectorType&& vec_dual, VectorType&& vec_freq) :
        _vec_prim(std::forward<VectorType>(vec_prim)),
        _vec_dual(std::forward<VectorType>(vec_dual)),
        _vec_freq(std::forward<VectorType>(vec_freq)),
        _volume()
      {
        if(!_vec_freq.empty())
        {
          _volume = _vec_freq.triple_dot(_vec_prim, _vec_dual);
          _volume = Global::SynchScal0::value(_volume);
        }
        else
          _volume = _vec_prim.dot(_vec_dual);

        XASSERTM(_volume > Math::eps<DataType>(), "domain volume must not be zero");
      }

      /// move ctor
      MeanFilter(MeanFilter && other) :
        _vec_prim(std::move(other._vec_prim)),
        _vec_dual(std::move(other._vec_dual)),
        _vec_freq(std::move(other._vec_freq)),
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
          _vec_freq = std::move(other._vec_freq);
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
      MeanFilter clone(LAFEM::CloneMode clone_mode = LAFEM::CloneMode::Deep) const
      {
        return MeanFilter(_vec_prim.clone(clone_mode), _vec_dual.clone(clone_mode), _vec_freq.clone(clone_mode));
      }

      /**
       * \brief Clones data from another MeanFilter
       * \warning _volume will always be a deep copy
       */
      void clone(const MeanFilter & other, LAFEM::CloneMode clone_mode = LAFEM::CloneMode::Deep)
      {
        _vec_prim.clone(other.get_vec_prim(), clone_mode);
        _vec_dual.clone(other.get_vec_dual(), clone_mode);
        _vec_freq.clone(other.get_vec_freq(), clone_mode);
        _volume = other.get_volume();
      }

      /// Conversion method
      template<typename Mem2_, typename DT2_, typename IT2_>
      void convert(const MeanFilter<Mem2_, DT2_, IT2_>& other)
      {
        _vec_prim.convert(other.get_vec_prim());
        _vec_dual.convert(other.get_vec_dual());
        _vec_freq.convert(other.get_vec_freq());
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

      VectorType & get_vec_freq()
      {
        return _vec_freq;
      }

      const VectorType & get_vec_freq() const
      {
        return _vec_freq;
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
        DataType integ = DataType(0);
        if(_vec_freq.empty())
          integ = vector.dot(_vec_prim);
        else
        {
          integ = _vec_freq.triple_dot(vector, _vec_prim);
          integ = Global::SynchScal0::value(integ);
        }
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
        DataType integ = DataType(0);
        if(_vec_freq.empty())
          integ = vector.dot(_vec_dual);
        else
        {
          integ = _vec_freq.triple_dot(vector, _vec_dual);
          integ = Global::SynchScal0::value(integ);
        }
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
  } // namespace Global
} // namespace FEAT


#endif // KERNEL_GLOBAL_MEAN_FILTER_HPP
