#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_LVLSET_MONITOR_2D_P1_HPP
#define KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_LVLSET_MONITOR_2D_P1_HPP 1

#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_p1.hpp>
// For HeavisideReg
#include <kernel/assembly/common_functions.hpp>

namespace FEAST
{
  namespace Geometry
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
        typedef FEAST::Tiny::Matrix<DataType, 2, 3> MonitorGradType;
        /// Factor for the levelset penalty term
        DataType fac_lvlset;
        /// Factor for making the regularised Heaviside function steeper
        static DataType heaviside_reg_fac(){ return DataType(10); };
        //static constexpr DataType heaviside_reg_fac = DataType(10);

      private:
        static constexpr DataType _mu_1 = DataType(1e3);
        static constexpr DataType _mu_2 = DataType(1e8);

        DataType _monitor;
        MonitorGradType _monitor_grad;

        DataType compute_monitor_val(DataType dist)
        {
          return DataType(1)/( DataType(1) + _mu_1 / Math::sqrt( DataType(1) + Math::sqr(_mu_2 * dist)));
          //return DataType(1);
          //return DataType(1)/(_mu_1 + _mu_2*Math::sqr(dist));
        }

        DataType compute_monitor_der(DataType dist)
        {
          DataType tmp = (Math::sqrt(DataType(1) + Math::sqr(_mu_2 * dist) ));
          return _mu_1*Math::sqr(_mu_2)*dist / (tmp * Math::sqr(_mu_1 + tmp));
          //return DataType(0);
          //return DataType(2)*_mu_2*dist/Math::sqr(_mu_1 + _mu_2*Math::sqr(dist));
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
          for(Index i(0); i < Tl_::n; ++i)
            tmp += lvlset_vals(i);
          tmp /= DataType(Tl_::n);

          _monitor = compute_monitor_val(tmp);

          MonitorGradType tmp_grad(DataType(0));
          for(Index d(0); d < Tlg_::m; ++d)
          {
            for(Index i(0); i < Tlg_::n; ++i)
              _monitor_grad(d,i) = compute_monitor_der(tmp) * lvlset_grad_vals(d,i) / DataType(Tlg_::m);
          }

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
          for(Index i(0); i < Index(3); ++i)
          {
            for(Index j(0); j < i; ++j)
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
          for(Index i(0); i < Index(3); ++i)
          {
            for(Index j(0); j < i; ++j)
            {
              auto lvlset_prod = -heaviside_reg_fac()*lvlset_vals(i)*lvlset_vals(j);
              // Derivative of the heaviside function
              auto heaviside_der = FEAST::Assembly::Common::template HeavisideRegStatic<DataType>::der_x(lvlset_prod);
              for(Index d(0); d < Index(2); ++d)
                grad(d,i) -= heaviside_reg_fac() * fac_lvlset * lvlset_constraint_last * heaviside_der *
                  (lvlset_grad_vals(d,i) * lvlset_vals(j));
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
          det_ = DataType( DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)));
          return det_;
        }

        /**
         * \brief Computes the 1/det term on one element
         **/
        template<typename Tx_, typename Th_ >
        DataType compute_rec_det_A( const Tx_& x, const Th_& h)
        {
          DataType rec_det_;
          rec_det_ = DataType( DataType(0.1e1) / (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1)));
          return rec_det_;
        }

        /**
         * \brief Computes the Frobenius norm term for one cell
         **/
        template<typename Tx_, typename Th_ >
        DataType compute_norm_A( const Tx_& x, const Th_& h)
        {
          DataType norm_;
          norm_ = DataType( DataType(0.4e1) + DataType(0.4e1) / DataType(0.9e1) * Math::pow(-DataType(0.2e1) * x(0,0) * x(0,1) + DataType(0.2e1) * Math::pow(x(0,1), DataType(0.2e1)) - DataType(0.2e1) * x(0,0) * x(0,2) + DataType(0.2e1) * Math::pow(x(0,0), DataType(0.2e1)) + DataType(0.2e1) * Math::pow(x(0,2), DataType(0.2e1)) + DataType(0.2e1) * Math::pow(x(1,0), DataType(0.2e1)) - DataType(0.2e1) * x(1,0) * x(1,1) - DataType(0.2e1) * x(0,1) * x(0,2) - DataType(0.2e1) * x(1,0) * x(1,2) - DataType(0.2e1) * x(1,1) * x(1,2) + DataType(0.2e1) * Math::pow(x(1,2), DataType(0.2e1)) + DataType(0.2e1) * Math::pow(x(1,1), DataType(0.2e1)), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.4e1)) * _monitor * _monitor - DataType(0.8e1) / DataType(0.3e1) * Math::pow(h(0), -DataType(0.2e1)) * (-DataType(0.2e1) * x(0,0) * x(0,1) + DataType(0.2e1) * Math::pow(x(0,1), DataType(0.2e1)) - DataType(0.2e1) * x(0,0) * x(0,2) + DataType(0.2e1) * Math::pow(x(0,0), DataType(0.2e1)) + DataType(0.2e1) * Math::pow(x(0,2), DataType(0.2e1)) + DataType(0.2e1) * Math::pow(x(1,0), DataType(0.2e1)) - DataType(0.2e1) * x(1,0) * x(1,1) - DataType(0.2e1) * x(0,1) * x(0,2) - DataType(0.2e1) * x(1,0) * x(1,2) - DataType(0.2e1) * x(1,1) * x(1,2) + DataType(0.2e1) * Math::pow(x(1,2), DataType(0.2e1)) + DataType(0.2e1) * Math::pow(x(1,1), DataType(0.2e1))) * _monitor);
          return norm_;

        }


        /**
         * \brief Computes the functional gradient for one cell
         **/
        template<typename Tx_, typename Th_, typename Tgrad_>
        void compute_local_grad( const Tx_& x, const Th_& h, Tgrad_& grad)
        {

          grad(0,0) = DataType( DataType(0.2e1) * this->_fac_norm * (_monitor * (Math::pow(x(0,0) - x(0,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1) + Math::pow(x(1,0) - x(1,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1)) - DataType(0.2e1)) * (_monitor_grad(0, 0) * _monitor * (Math::pow(x(0,0) - x(0,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1) + Math::pow(x(1,0) - x(1,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1)) + _monitor * _monitor * (DataType(0.2e1) * (x(0,0) - x(0,1)) * Math::pow(h(0), -DataType(0.2e1)) + DataType(0.2e1) / DataType(0.3e1) * (x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2)) * Math::pow(h(0), -DataType(0.2e1)))) + DataType(0.2e1) / DataType(0.3e1) * this->_fac_det * Math::sqrt(DataType(0.3e1)) * (x(1,1) - x(1,2)) * Math::pow(h(1), -DataType(0.2e1)) - DataType(0.1e1) * this->_fac_rec_det * Math::pow(DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1), -DataType(0.2e1)) * (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(1,1) - x(1,2)) * Math::pow(h(1), -DataType(0.2e1)) + DataType(0.4e1) * Math::pow(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1)), -DataType(0.1e1) / DataType(0.2e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.4e1)) * (x(1,1) - x(1,2))));
          grad(0,1) = DataType( DataType(0.2e1) * this->_fac_norm * (_monitor * (Math::pow(x(0,0) - x(0,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1) + Math::pow(x(1,0) - x(1,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1)) - DataType(0.2e1)) * (_monitor_grad(0, 1) * _monitor * (Math::pow(x(0,0) - x(0,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1) + Math::pow(x(1,0) - x(1,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1)) + _monitor * _monitor * (-DataType(0.2e1) * (x(0,0) - x(0,1)) * Math::pow(h(0), -DataType(0.2e1)) + DataType(0.2e1) / DataType(0.3e1) * (x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2)) * Math::pow(h(0), -DataType(0.2e1)))) + DataType(0.2e1) / DataType(0.3e1) * this->_fac_det * Math::sqrt(DataType(0.3e1)) * (-x(1,0) + x(1,2)) * Math::pow(h(1), -DataType(0.2e1)) - DataType(0.1e1) * this->_fac_rec_det * Math::pow(DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1), -DataType(0.2e1)) * (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (-x(1,0) + x(1,2)) * Math::pow(h(1), -DataType(0.2e1)) + DataType(0.4e1) * Math::pow(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1)), -DataType(0.1e1) / DataType(0.2e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.4e1)) * (-x(1,0) + x(1,2))));
          grad(1,0) = DataType( DataType(0.2e1) * this->_fac_norm * (_monitor * (Math::pow(x(0,0) - x(0,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1) + Math::pow(x(1,0) - x(1,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1)) - DataType(0.2e1)) * (_monitor_grad(1, 0) * _monitor * (Math::pow(x(0,0) - x(0,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1) + Math::pow(x(1,0) - x(1,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1)) + _monitor * _monitor * (DataType(0.2e1) * (x(1,0) - x(1,1)) * Math::pow(h(0), -DataType(0.2e1)) + DataType(0.2e1) / DataType(0.3e1) * (x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2)) * Math::pow(h(0), -DataType(0.2e1)))) + DataType(0.2e1) / DataType(0.3e1) * this->_fac_det * Math::sqrt(DataType(0.3e1)) * (-x(0,1) + x(0,2)) * Math::pow(h(1), -DataType(0.2e1)) - DataType(0.1e1) * this->_fac_rec_det * Math::pow(DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1), -DataType(0.2e1)) * (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (-x(0,1) + x(0,2)) * Math::pow(h(1), -DataType(0.2e1)) + DataType(0.4e1) * Math::pow(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1)), -DataType(0.1e1) / DataType(0.2e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.4e1)) * (-x(0,1) + x(0,2))));
          grad(1,1) = DataType( DataType(0.2e1) * this->_fac_norm * (_monitor * (Math::pow(x(0,0) - x(0,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1) + Math::pow(x(1,0) - x(1,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1)) - DataType(0.2e1)) * (_monitor_grad(1, 1) * _monitor * (Math::pow(x(0,0) - x(0,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1) + Math::pow(x(1,0) - x(1,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1)) + _monitor * _monitor * (-DataType(0.2e1) * (x(1,0) - x(1,1)) * Math::pow(h(0), -DataType(0.2e1)) + DataType(0.2e1) / DataType(0.3e1) * (x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2)) * Math::pow(h(0), -DataType(0.2e1)))) + DataType(0.2e1) / DataType(0.3e1) * this->_fac_det * Math::sqrt(DataType(0.3e1)) * (x(0,0) - x(0,2)) * Math::pow(h(1), -DataType(0.2e1)) - DataType(0.1e1) * this->_fac_rec_det * Math::pow(DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1), -DataType(0.2e1)) * (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) - x(0,2)) * Math::pow(h(1), -DataType(0.2e1)) + DataType(0.4e1) * Math::pow(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1)), -DataType(0.1e1) / DataType(0.2e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.4e1)) * (x(0,0) - x(0,2))));
          grad(0,2) = DataType( DataType(0.2e1) * this->_fac_norm * (_monitor * (Math::pow(x(0,0) - x(0,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1) + Math::pow(x(1,0) - x(1,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1)) - DataType(0.2e1)) * (_monitor_grad(0, 2) * _monitor * (Math::pow(x(0,0) - x(0,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1) + Math::pow(x(1,0) - x(1,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1)) - DataType(0.4e1) / DataType(0.3e1) * _monitor * _monitor * (x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2)) * Math::pow(h(0), -DataType(0.2e1))) + DataType(0.2e1) / DataType(0.3e1) * this->_fac_det * Math::sqrt(DataType(0.3e1)) * (x(1,0) - x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) - DataType(0.1e1) * this->_fac_rec_det * Math::pow(DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1), -DataType(0.2e1)) * (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(1,0) - x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + DataType(0.4e1) * Math::pow(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1)), -DataType(0.1e1) / DataType(0.2e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.4e1)) * (x(1,0) - x(1,1))));
          grad(1,2) = DataType( DataType(0.2e1) * this->_fac_norm * (_monitor * (Math::pow(x(0,0) - x(0,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1) + Math::pow(x(1,0) - x(1,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1)) - DataType(0.2e1)) * (_monitor_grad(1, 2) * _monitor * (Math::pow(x(0,0) - x(0,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(0,0) + x(0,1) - DataType(0.2e1) * x(0,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1) + Math::pow(x(1,0) - x(1,1), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) + Math::pow(x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.2e1)) / DataType(0.3e1)) - DataType(0.4e1) / DataType(0.3e1) * _monitor * _monitor * (x(1,0) + x(1,1) - DataType(0.2e1) * x(1,2)) * Math::pow(h(0), -DataType(0.2e1))) + DataType(0.2e1) / DataType(0.3e1) * this->_fac_det * Math::sqrt(DataType(0.3e1)) * (-x(0,0) + x(0,1)) * Math::pow(h(1), -DataType(0.2e1)) - DataType(0.1e1) * this->_fac_rec_det * Math::pow(DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1), -DataType(0.2e1)) * (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (-x(0,0) + x(0,1)) * Math::pow(h(1), -DataType(0.2e1)) + DataType(0.4e1) * Math::pow(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1)), -DataType(0.1e1) / DataType(0.2e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.4e1)) * (-x(0,0) + x(0,1))));


          return;
        }

    }; // class RumpfFunctionalLevelsetMonitor
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_LVLSET_MONITOR_2D_P1_HPP
