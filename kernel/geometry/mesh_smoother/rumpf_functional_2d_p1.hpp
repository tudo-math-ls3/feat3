#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_2D_P1_HPP
#define KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_2D_P1_HPP 1

#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal

    /**
     * @brief Class template for Rumpf functionals, 2d P1
     **/
    template<typename DataType_>
    class RumpfFunctional<DataType_, Shape::Simplex<2> > :
    public RumpfFunctionalBase<DataType_>
    {

      public:
        /// Our data type
        typedef DataType_ DataType;
        /// Shape type of the underlying transformation
        typedef Shape::Simplex<2> ShapeType;
        /// Our baseclass
        typedef RumpfFunctionalBase<DataType> BaseClass;

        /**
         * \brief Constructor
         **/
        RumpfFunctional(
          const DataType fac_norm_,
          const DataType fac_det_,
          const DataType fac_cof_,
          const DataType fac_reg_) :
          BaseClass( fac_norm_,
          fac_det_,
          fac_det_*( Math::sqrt( Math::sqr(fac_reg_) + DataType(1) ) + Math::sqr(fac_reg_) + DataType(1)),
          fac_cof_,
          fac_reg_)
          {
          }

        /**
         * \brief Computes value the Rumpf functional on one element.
         **/
        template<typename Tx_, typename Th_>
        DataType compute_local_functional(const Tx_& x, const Th_& h)
        {
          DataType norm_A = compute_norm_A(x,h);
          DataType det_A = compute_det_A(x,h);
          DataType rec_det_A = compute_rec_det_A(x,h);

          return this->_fac_norm*norm_A
            + this->_fac_det*det_A
            + this->_fac_rec_det*rec_det_A;

        }

        /**
         * \copydoc compute_local_functional()
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
        DataType compute_norm_A(const Tx_& x, const Th_& h)
        {
          DataType norm_;
          norm_ = DataType( DataType(0.4e1) / DataType(0.9e1) * Math::pow(-DataType(0.2e1) * Math::pow(x(0,0), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,1) - DataType(0.2e1) * Math::pow(x(0,1), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,2) + DataType(0.2e1) * x(0,1) * x(0,2) - DataType(0.2e1) * Math::pow(x(0,2), DataType(0.2e1)) - DataType(0.2e1) * Math::pow(x(1,0), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,1) - DataType(0.2e1) * Math::pow(x(1,1), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,2) + DataType(0.2e1) * x(1,1) * x(1,2) - DataType(0.2e1) * Math::pow(x(1,2), DataType(0.2e1)) + DataType(0.3e1) * Math::pow(h(0), DataType(0.2e1)), DataType(0.2e1)) * Math::pow(h(0), -DataType(0.4e1)));
          return norm_;
        }

        /**
         * \brief Computes the functional gradient for one cell
         **/
        template<typename Tx_, typename Th_, typename Tgrad_>
        void compute_local_grad( const Tx_& x, const Th_& h, Tgrad_& grad)
        {
          grad(0,0) = DataType( DataType(0.8e1) / DataType(0.9e1) * this->_fac_norm * (-DataType(0.2e1) * Math::pow(x(0,0), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,1) - DataType(0.2e1) * Math::pow(x(0,1), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,2) + DataType(0.2e1) * x(0,1) * x(0,2) - DataType(0.2e1) * Math::pow(x(0,2), DataType(0.2e1)) - DataType(0.2e1) * Math::pow(x(1,0), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,1) - DataType(0.2e1) * Math::pow(x(1,1), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,2) + DataType(0.2e1) * x(1,1) * x(1,2) - DataType(0.2e1) * Math::pow(x(1,2), DataType(0.2e1)) + DataType(0.3e1) * Math::pow(h(0), DataType(0.2e1))) * Math::pow(h(0), -DataType(0.4e1)) * (-DataType(0.4e1) * x(0,0) + DataType(0.2e1) * x(0,1) + DataType(0.2e1) * x(0,2)) + DataType(0.2e1) / DataType(0.3e1) * this->_fac_det * Math::sqrt(DataType(0.3e1)) * (x(1,1) - x(1,2)) * Math::pow(h(1), -DataType(0.2e1)) - this->_fac_rec_det * Math::pow(DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1), -DataType(0.2e1)) * (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(1,1) - x(1,2)) * Math::pow(h(1), -DataType(0.2e1)) + DataType(0.4e1) * Math::pow(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1)), -DataType(0.1e1) / DataType(0.2e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.4e1)) * (x(1,1) - x(1,2))));
          grad(0,1) = DataType( DataType(0.8e1) / DataType(0.9e1) * this->_fac_norm * (-DataType(0.2e1) * Math::pow(x(0,0), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,1) - DataType(0.2e1) * Math::pow(x(0,1), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,2) + DataType(0.2e1) * x(0,1) * x(0,2) - DataType(0.2e1) * Math::pow(x(0,2), DataType(0.2e1)) - DataType(0.2e1) * Math::pow(x(1,0), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,1) - DataType(0.2e1) * Math::pow(x(1,1), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,2) + DataType(0.2e1) * x(1,1) * x(1,2) - DataType(0.2e1) * Math::pow(x(1,2), DataType(0.2e1)) + DataType(0.3e1) * Math::pow(h(0), DataType(0.2e1))) * Math::pow(h(0), -DataType(0.4e1)) * (DataType(0.2e1) * x(0,0) - DataType(0.4e1) * x(0,1) + DataType(0.2e1) * x(0,2)) + DataType(0.2e1) / DataType(0.3e1) * this->_fac_det * Math::sqrt(DataType(0.3e1)) * (-x(1,0) + x(1,2)) * Math::pow(h(1), -DataType(0.2e1)) - this->_fac_rec_det * Math::pow(DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1), -DataType(0.2e1)) * (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (-x(1,0) + x(1,2)) * Math::pow(h(1), -DataType(0.2e1)) + DataType(0.4e1) * Math::pow(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1)), -DataType(0.1e1) / DataType(0.2e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.4e1)) * (-x(1,0) + x(1,2))));
          grad(1,0) = DataType( DataType(0.8e1) / DataType(0.9e1) * this->_fac_norm * (-DataType(0.2e1) * Math::pow(x(0,0), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,1) - DataType(0.2e1) * Math::pow(x(0,1), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,2) + DataType(0.2e1) * x(0,1) * x(0,2) - DataType(0.2e1) * Math::pow(x(0,2), DataType(0.2e1)) - DataType(0.2e1) * Math::pow(x(1,0), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,1) - DataType(0.2e1) * Math::pow(x(1,1), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,2) + DataType(0.2e1) * x(1,1) * x(1,2) - DataType(0.2e1) * Math::pow(x(1,2), DataType(0.2e1)) + DataType(0.3e1) * Math::pow(h(0), DataType(0.2e1))) * Math::pow(h(0), -DataType(0.4e1)) * (-DataType(0.4e1) * x(1,0) + DataType(0.2e1) * x(1,1) + DataType(0.2e1) * x(1,2)) + DataType(0.2e1) / DataType(0.3e1) * this->_fac_det * Math::sqrt(DataType(0.3e1)) * (-x(0,1) + x(0,2)) * Math::pow(h(1), -DataType(0.2e1)) - this->_fac_rec_det * Math::pow(DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1), -DataType(0.2e1)) * (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (-x(0,1) + x(0,2)) * Math::pow(h(1), -DataType(0.2e1)) + DataType(0.4e1) * Math::pow(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1)), -DataType(0.1e1) / DataType(0.2e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.4e1)) * (-x(0,1) + x(0,2))));
          grad(1,1) = DataType( DataType(0.8e1) / DataType(0.9e1) * this->_fac_norm * (-DataType(0.2e1) * Math::pow(x(0,0), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,1) - DataType(0.2e1) * Math::pow(x(0,1), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,2) + DataType(0.2e1) * x(0,1) * x(0,2) - DataType(0.2e1) * Math::pow(x(0,2), DataType(0.2e1)) - DataType(0.2e1) * Math::pow(x(1,0), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,1) - DataType(0.2e1) * Math::pow(x(1,1), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,2) + DataType(0.2e1) * x(1,1) * x(1,2) - DataType(0.2e1) * Math::pow(x(1,2), DataType(0.2e1)) + DataType(0.3e1) * Math::pow(h(0), DataType(0.2e1))) * Math::pow(h(0), -DataType(0.4e1)) * (DataType(0.2e1) * x(1,0) - DataType(0.4e1) * x(1,1) + DataType(0.2e1) * x(1,2)) + DataType(0.2e1) / DataType(0.3e1) * this->_fac_det * Math::sqrt(DataType(0.3e1)) * (x(0,0) - x(0,2)) * Math::pow(h(1), -DataType(0.2e1)) - this->_fac_rec_det * Math::pow(DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1), -DataType(0.2e1)) * (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) - x(0,2)) * Math::pow(h(1), -DataType(0.2e1)) + DataType(0.4e1) * Math::pow(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1)), -DataType(0.1e1) / DataType(0.2e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.4e1)) * (x(0,0) - x(0,2))));
          grad(0,2) = DataType( DataType(0.8e1) / DataType(0.9e1) * this->_fac_norm * (-DataType(0.2e1) * Math::pow(x(0,0), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,1) - DataType(0.2e1) * Math::pow(x(0,1), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,2) + DataType(0.2e1) * x(0,1) * x(0,2) - DataType(0.2e1) * Math::pow(x(0,2), DataType(0.2e1)) - DataType(0.2e1) * Math::pow(x(1,0), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,1) - DataType(0.2e1) * Math::pow(x(1,1), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,2) + DataType(0.2e1) * x(1,1) * x(1,2) - DataType(0.2e1) * Math::pow(x(1,2), DataType(0.2e1)) + DataType(0.3e1) * Math::pow(h(0), DataType(0.2e1))) * Math::pow(h(0), -DataType(0.4e1)) * (DataType(0.2e1) * x(0,0) + DataType(0.2e1) * x(0,1) - DataType(0.4e1) * x(0,2)) + DataType(0.2e1) / DataType(0.3e1) * this->_fac_det * Math::sqrt(DataType(0.3e1)) * (x(1,0) - x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) - this->_fac_rec_det * Math::pow(DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1), -DataType(0.2e1)) * (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(1,0) - x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + DataType(0.4e1) * Math::pow(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1)), -DataType(0.1e1) / DataType(0.2e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.4e1)) * (x(1,0) - x(1,1))));
          grad(1,2) = DataType( DataType(0.8e1) / DataType(0.9e1) * this->_fac_norm * (-DataType(0.2e1) * Math::pow(x(0,0), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,1) - DataType(0.2e1) * Math::pow(x(0,1), DataType(0.2e1)) + DataType(0.2e1) * x(0,0) * x(0,2) + DataType(0.2e1) * x(0,1) * x(0,2) - DataType(0.2e1) * Math::pow(x(0,2), DataType(0.2e1)) - DataType(0.2e1) * Math::pow(x(1,0), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,1) - DataType(0.2e1) * Math::pow(x(1,1), DataType(0.2e1)) + DataType(0.2e1) * x(1,0) * x(1,2) + DataType(0.2e1) * x(1,1) * x(1,2) - DataType(0.2e1) * Math::pow(x(1,2), DataType(0.2e1)) + DataType(0.3e1) * Math::pow(h(0), DataType(0.2e1))) * Math::pow(h(0), -DataType(0.4e1)) * (DataType(0.2e1) * x(1,0) + DataType(0.2e1) * x(1,1) - DataType(0.4e1) * x(1,2)) + DataType(0.2e1) / DataType(0.3e1) * this->_fac_det * Math::sqrt(DataType(0.3e1)) * (-x(0,0) + x(0,1)) * Math::pow(h(1), -DataType(0.2e1)) - this->_fac_rec_det * Math::pow(DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.2e1)) + Math::sqrt(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1))) / DataType(0.3e1), -DataType(0.2e1)) * (DataType(0.2e1) / DataType(0.3e1) * Math::sqrt(DataType(0.3e1)) * (-x(0,0) + x(0,1)) * Math::pow(h(1), -DataType(0.2e1)) + DataType(0.4e1) * Math::pow(DataType(0.9e1) * this->_fac_reg * this->_fac_reg + DataType(0.12e2) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(0.2e1)) * Math::pow(h(1), -DataType(0.4e1)), -DataType(0.1e1) / DataType(0.2e1)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(0.4e1)) * (-x(0,0) + x(0,1))));

          return;
        }

    }; // class RumpfFunctional
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_2D_P1_HPP
