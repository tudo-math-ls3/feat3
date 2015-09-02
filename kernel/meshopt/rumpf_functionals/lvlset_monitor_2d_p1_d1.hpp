#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_FUNCTIONALS_LVLSET_MONITOR_2D_P1_D1_HPP
#define KERNEL_MESHOPT_RUMPF_FUNCTIONALS_LVLSET_MONITOR_2D_P1_D1_HPP 1

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/meshopt/rumpf_functional.hpp>
// For HeavisideReg
#include <kernel/assembly/common_functions.hpp>

namespace FEAST
{
  namespace Meshopt
  {
    /// \cond internal

    /**
     * @brief Class template for Rumpf functionals with levelset term and monitor function , 2d Q1
     **/
    template<typename MemoryType_, typename DataType_>
    class RumpfFunctionalLevelsetMonitor < MemoryType_, DataType_, Shape::Simplex<2>> :
    public RumpfFunctionalBase<MemoryType_, DataType_>
    {
      public:
        /// Our data type
        typedef DataType_ DataType;
        /// Shape type of the underlying transformation
        typedef Shape::Simplex<2> ShapeType;
        /// Our baseclass
        typedef RumpfFunctionalBase<MemoryType_, DataType> BaseClass;
        typedef FEAST::Tiny::Matrix<DataType, 3, 2> MonitorGradType;
        /// Factor for the levelset penalty term
        DataType fac_lvlset;
        /// Factor for making the regularised Heaviside function steeper
        static DataType heaviside_reg_fac(){ return DataType(10); };
        //static constexpr DataType heaviside_reg_fac = DataType(10);

      private:
        static constexpr DataType _mu_1 = DataType(1e3);
        static constexpr DataType _mu_2 = DataType(1e8);

        /// Value of the monitor function on the current cell
        DataType _monitor;
        /// Derivative of the monitor function on the current cell
        MonitorGradType _monitor_grad;

        /**
         * \brief Evaluates the monitor function
         *
         * \param[in] dist
         * Distance to evaluate the monitor function at.
         *
         * \returns
         * The value of the monitor function.
         */
        DataType compute_monitor_val(DataType dist)
        {
          return DataType(1)/( DataType(1) + _mu_1 /
              Math::sqrt( DataType(1) + Math::sqr(_mu_2 * dist)));
        }

        /**
         * \brief Evaluates the monitor function's derivative
         *
         * \param[in] dist
         * Distance to evaluate the monitor function's derivative at.
         *
         * \returns
         * The value of the monitor function's derivative wrt. the distance
         */
        DataType compute_monitor_der(DataType dist)
        {
          DataType tmp = (Math::sqrt(DataType(1) + Math::sqr(_mu_2 * dist) ));
          return _mu_1*Math::sqr(_mu_2)*dist / (tmp * Math::sqr(_mu_1 + tmp));
        }

      public:
        template<typename Tl_>
        void prepare_lvlset(const Tl_& lvlset_vals)
        {
          DataType tmp(DataType(0));
          for(Index i(0); i < Tl_::n; ++i)
            tmp += lvlset_vals(i);
          tmp /= DataType(Tl_::n);

          _monitor = compute_monitor_val(tmp);
        }

        template<typename Tl_, typename Tlg_>
        void prepare_lvlset(const Tl_& lvlset_vals, const Tlg_& lvlset_grad_vals)
        {
          DataType tmp(DataType(0));
          for(int i(0); i < Tl_::n; ++i)
            tmp += lvlset_vals(i);
          tmp /= DataType(Tl_::n);

          _monitor = compute_monitor_val(tmp);

          MonitorGradType tmp_grad(DataType(0));
          for(int i(0); i < Tlg_::n; ++i)
            _monitor_grad[i] = compute_monitor_der(tmp) / DataType(Tlg_::m) * lvlset_grad_vals[i] ;

        }
        /**
         * \brief Constructor
         **/
        RumpfFunctionalLevelsetMonitor(
          const DataType fac_norm_,
          const DataType fac_det_,
          const DataType fac_cof_,
          const DataType fac_reg_):
          BaseClass( fac_norm_,
          fac_det_/( Math::sqrt( Math::sqr(fac_reg_) + DataType(1) ) + Math::sqr(fac_reg_) + DataType(1)),
          fac_det_,
          fac_cof_,
          fac_reg_),
          fac_lvlset(DataType(1)),
          _monitor(DataType(1)),
          _monitor_grad(DataType(0))
          {
          }

        /**
         * \brief Computes the additional levelset penalty term for the Rumpf functional
         **/
        template<typename Tl_>
        DataType compute_lvlset_penalty(const Tl_& lvlset_vals)
        {
          DataType penalty(0);
          for(int i(0); i < 3; ++i)
          {
            for(int j(0); j < i; ++j)
              penalty += FEAST::Assembly::Common::template HeavisideRegStatic<DataType>::eval(-heaviside_reg_fac()*lvlset_vals(i)*lvlset_vals(j));
          }

          return penalty;
        } // compute_lvlset_penalty

        /**
         * \brief Adds the gradient of the additional levelset penalty term
         **/
        template<typename Tl_, typename Tlg_, typename Tgrad_>
        void add_lvlset_penalty_grad(const Tl_& lvlset_vals, const Tlg_& lvlset_grad_vals, Tgrad_& grad, DataType& lvlset_constraint_last)
        {
          // Compute local gradient
          for(int i(0); i < 3; ++i)
          {
            for(int j(0); j < i; ++j)
            {
              auto lvlset_prod = -heaviside_reg_fac()*lvlset_vals(i)*lvlset_vals(j);
              // Derivative of the heaviside function
              auto heaviside_der = FEAST::Assembly::Common::template HeavisideRegStatic<DataType>::der_x(lvlset_prod);
                grad[i] -= heaviside_reg_fac() * fac_lvlset * lvlset_constraint_last * heaviside_der
                  * lvlset_vals(j) lvlset_grad_vals[i];
            }
          }

        } // add_levelset_penalty_grad

        /**
         *  \brief Computes value the Rumpf functional on one element.
         **/
        template<typename Tx_, typename Th_>
        DataType compute_local_functional(const Tx_& x, const Th_& h)
        {
          DataType norm_A = compute_norm_A(x,h);
          DataType det_A = compute_det_A(x,h);
          DataType rec_det_A = compute_rec_det_A(x,h);

          return this->_fac_norm*norm_A + this->_fac_det*det_A + this->_fac_rec_det*rec_det_A;

        }

        /**
         * \copydoc compute_local_functional
         * Debug variant that also returns the contributions of the different terms seperately.
         **/
        template<typename Tx_, typename Th_>
        DataType compute_local_functional(const Tx_& x, const Th_& h,
        DataType& func_norm,
        DataType& func_det,
        DataType& func_rec_det)
        {
          func_norm = this->_fac_norm*compute_norm_A(x,h);
          func_det = this->_fac_det*compute_det_A(x,h);
          func_rec_det = this->_fac_rec_det*compute_rec_det_A(x,h);

          return func_norm + func_det + func_rec_det;
        }

        /**
         * \brief Computes the det term on one element
         **/
        template<typename Tx_, typename Th_ >
        DataType compute_det_A( const Tx_& x, const Th_& h)
        {
          DataType det_;
          det_ =
          return det_;
        }

        /**
         * \brief Computes the 1/det term on one element
         **/
        template<typename Tx_, typename Th_ >
        DataType compute_rec_det_A( const Tx_& x, const Th_& h)
        {
          DataType rec_det_;
          rec_det_ =
          return rec_det_;
        }

        /**
         * \brief Computes the Frobenius norm term for one cell
         **/
        template<typename Tx_, typename Th_ >
        DataType compute_norm_A( const Tx_& x, const Th_& h)
        {
          DataType norm_;
          norm_ =
          return norm_;

        }

        /**
         * \brief Computes the functional gradient for one cell
         **/
        template<typename Tx_, typename Th_, typename Tgrad_>
        void compute_local_grad( const Tx_& x, const Th_& h, Tgrad_& grad)
        {

          grad(0,0) =
          grad(0,1) =
          grad(1,0) =
          grad(1,1) =
          grad(2,0) =
          grad(2,1) =

          return;
        }

    }; // class RumpfFunctionalLevelsetMonitor
    /// \endcond
  } // namespace Meshopt
} // namespace FEAST

#endif // KERNEL_MESHOPT_RUMPF_FUNCTIONALS_LVLSET_MONITOR_2D_P1_D1_HPP
