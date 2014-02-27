#pragma once
#ifndef KERNEL_LAFEM_MEAN_FILTER_HPP
#define KERNEL_LAFEM_MEAN_FILTER_HPP 1

// includes, FEAST
#include <kernel/lafem/dense_vector.hpp>

namespace FEAST
{
  namespace LAFEM
  {
    /**
     * \brief Mean Filter class template.
     *
     * \author Peter Zajac
     */
    template<
      typename Mem_,
      typename DT_>
    class MeanFilter
    {
    public:
      /// mem-type typedef
      typedef Mem_ MemType;
      /// data-type typedef
      typedef DT_ DataType;
      /// vector-type typedef
      typedef DenseVector<Mem_, DT_> VectorType;

    protected:
      /// primal weighting vector
      VectorType _vec_prim;
      /// dual weighting vector
      VectorType _vec_dual;
      /// weight volume
      DT_ _volume;

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
       *
       * \param[in] volume
       * The weight volume. This is simply the dot-product of the primal and dual weighting vector.
       */
      explicit MeanFilter(VectorType && vec_prim, VectorType && vec_dual, DT_ volume) :
        _vec_prim(vec_prim),
        _vec_dual(vec_dual),
        _volume(volume)
      {
        ASSERT(_volume > Math::Limits<DT_>::epsilon(), "domain volume must not be zero");
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

      /// Creates and returns a deep copy of this filter.
      MeanFilter clone() const
      {
        return MeanFilter(_vec_prim.clone(), _vec_dual.clone(), _volume);
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
      /// \endcond

      template<typename Algo_, typename Matrix_>
      void filter_mat(Matrix_ &) const
      {
        // nothing to do here
      }

      template<typename Algo_, typename Matrix_>
      void filter_offdiag_row_mat(Matrix_ &) const
      {
        // nothing to do here
      }

      template<typename Algo_, typename Matrix_>
      void filter_offdiag_col_mat(Matrix_ &) const
      {
        // nothing to do here
      }

      /**
       * \brief Applies the filter onto the right-hand-side vector.
       *
       * \param[in,out] vector
       * A reference to the right-hand-side vector to be filtered.
       */
      template<typename Algo_>
      void filter_rhs(DenseVector<Mem_,DT_> & vector) const
      {
        // compute dual integral
        DT_ integ = vector.template dot<Algo_>(_vec_prim);
        // subtract mean
        vector.template axpy<Algo_>(_vec_dual, vector, -integ / _volume);
      }

      /**
       * \brief Applies the filter onto the solution vector.
       *
       * \param[in,out] vector
       * A reference to the solution vector to be filtered.
       */
      template<typename Algo_>
      void filter_sol(DenseVector<Mem_,DT_> & vector) const
      {
        // compute primal integral
        DT_ integ = vector.template dot<Algo_>(_vec_dual);
        // subtract mean
        vector.template axpy<Algo_>(_vec_prim, vector, -integ / _volume);
      }

      /**
       * \brief Applies the filter onto a defect vector.
       *
       * \param[in,out] vector
       * A reference to the defect vector to be filtered.
       */
      template<typename Algo_>
      void filter_def(DenseVector<Mem_,DT_> & vector) const
      {
        // same as rhs
        filter_rhs<Algo_>(vector);
      }

      /**
       * \brief Applies the filter onto a correction vector.
       *
       * \param[in,out] vector
       * A reference to the correction vector to be filtered.
       */
      template<typename Algo_>
      void filter_cor(DenseVector<Mem_,DT_> & vector) const
      {
        // same as sol
        filter_sol<Algo_>(vector);
      }
    }; // class MeanFilter<...>
  } // namespace LAFEM
} // namespace FEAST


#endif // KERNEL_LAFEM_MEAN_FILTER_HPP
