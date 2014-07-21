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
    template<typename MemoryType_, typename DataType_>
    class RumpfFunctional<MemoryType_, DataType_, Shape::Simplex<2> > :
    public RumpfFunctionalBase<MemoryType_, DataType_>
    {

      public:
        /// Shape type of the underlying transformation
        typedef Shape::Simplex<2> ShapeType;
        /// Our baseclass
        typedef RumpfFunctionalBase<MemoryType_, DataType_> BaseClass;

        /**
         * \brief Constructor
         **/
        RumpfFunctional(
          const DataType_ fac_norm_,
          const DataType_ fac_det_,
          const DataType_ fac_cof_,
          const DataType_ fac_reg_) :
          BaseClass( fac_norm_,
          fac_det_/( Math::sqrt( Math::sqr(fac_reg_) + DataType_(1) ) + Math::sqr(fac_reg_) + DataType_(1)),
          fac_det_,
          fac_cof_,
          fac_reg_)
          {
          }

        /**
         * \brief Computes value the Rumpf functional on one element.
         **/
        template<typename Tx_, typename Th_>
        DataType_ compute_local_functional(const Tx_& x, const Th_& h)
        {
          DataType_ norm_A = compute_norm_A(x,h);
          DataType_ det_A = compute_det_A(x,h);
          DataType_ det2_A = compute_det2_A(x,h);

          return this->_fac_norm*norm_A
            + this->_fac_det*det_A
            + this->_fac_det2*det2_A;

        }

        /**
         * \copydoc compute_local_functional()
         **/
        template<typename Tx_, typename Th_>
        DataType_ compute_local_functional(const Tx_& x, const Th_& h,
        DataType_& func_norm,
        DataType_& func_det,
        DataType_& func_det2)
        {
          func_norm = this->_fac_norm*compute_norm_A(x,h);
          func_det = this->_fac_det*compute_det_A(x,h);
          func_det2 = this->_fac_det2*compute_det2_A(x,h);

          return func_norm + func_det + func_det2;

        }

        /**
         * \brief Computes the det term on one element
         **/
        template<typename Tx_, typename Th_ >
        DataType_ compute_det_A( const Tx_& x, const Th_& h)
        {
          return DataType_(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.2e1));
        }

        /**
         * \brief Computes the 1/det term on one element
         **/
        template<typename Tx_, typename Th_ >
        DataType_ compute_det2_A( const Tx_& x, const Th_& h)
        {
          return DataType_(0.1e1 / (0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1)) / 0.3e1));
        }

        /**
         * \brief Computes the Frobenius norm term for one cell
         **/
        template<typename Tx_, typename Th_ >
        DataType_ compute_norm_A(const Tx_& x, const Th_& h)
        {

          return DataType_(pow(pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1, 0.2e1));

        }

        /**
         * \brief Computes the functional gradient for one cell
         **/
        template<typename Tx_, typename Th_, typename Tgrad_>
        void compute_local_grad( const Tx_& x, const Th_& h, Tgrad_& grad)
        {
          grad(0,0) = DataType_(0.2e1 * this->_fac_norm * (pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1) * (-0.2e1 * pow(h(0), -0.2e1) * (-x(0,0) + x(0,1)) - 0.2e1 / 0.3e1 * pow(h(0), -0.2e1) * (-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1)) * sqrt(0.3e1)) + 0.2e1 / 0.3e1 * this->_fac_det * sqrt(0.3e1) * (x(1,1) - x(1,2)) * pow(h(1), -0.2e1) - 0.1e1 * this->_fac_det2 * pow(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1)) / 0.3e1, -0.2e1) * (0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(1,1) - x(1,2)) * pow(h(1), -0.2e1) + 0.4e1 * pow((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1), -0.1e1 / 0.2e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.4e1) * (x(1,1) - x(1,2))));
          grad(0,1) = DataType_(0.2e1 * this->_fac_norm * (pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1) * (0.2e1 * pow(h(0), -0.2e1) * (-x(0,0) + x(0,1)) - 0.2e1 / 0.3e1 * pow(h(0), -0.2e1) * (-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1)) * sqrt(0.3e1)) + 0.2e1 / 0.3e1 * this->_fac_det * sqrt(0.3e1) * (-x(1,0) + x(1,2)) * pow(h(1), -0.2e1) - 0.1e1 * this->_fac_det2 * pow(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1)) / 0.3e1, -0.2e1) * (0.2e1 / 0.3e1 * sqrt(0.3e1) * (-x(1,0) + x(1,2)) * pow(h(1), -0.2e1) + 0.4e1 * pow((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1), -0.1e1 / 0.2e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.4e1) * (-x(1,0) + x(1,2))));
          grad(1,0) = DataType_(0.2e1 * this->_fac_norm * (pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1) * (-0.2e1 * pow(h(0), -0.2e1) * (-x(1,0) + x(1,1)) - 0.2e1 / 0.3e1 * pow(h(0), -0.2e1) * (-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1)) * sqrt(0.3e1)) + 0.2e1 / 0.3e1 * this->_fac_det * sqrt(0.3e1) * (-x(0,1) + x(0,2)) * pow(h(1), -0.2e1) - 0.1e1 * this->_fac_det2 * pow(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1)) / 0.3e1, -0.2e1) * (0.2e1 / 0.3e1 * sqrt(0.3e1) * (-x(0,1) + x(0,2)) * pow(h(1), -0.2e1) + 0.4e1 * pow((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1), -0.1e1 / 0.2e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.4e1) * (-x(0,1) + x(0,2))));
          grad(1,1) = DataType_(0.2e1 * this->_fac_norm * (pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1) * (0.2e1 * pow(h(0), -0.2e1) * (-x(1,0) + x(1,1)) - 0.2e1 / 0.3e1 * pow(h(0), -0.2e1) * (-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1)) * sqrt(0.3e1)) + 0.2e1 / 0.3e1 * this->_fac_det * sqrt(0.3e1) * (x(0,0) - x(0,2)) * pow(h(1), -0.2e1) - 0.1e1 * this->_fac_det2 * pow(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1)) / 0.3e1, -0.2e1) * (0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) - x(0,2)) * pow(h(1), -0.2e1) + 0.4e1 * pow((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1), -0.1e1 / 0.2e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.4e1) * (x(0,0) - x(0,2))));
          grad(0,2) = DataType_(0.8e1 / 0.3e1 * this->_fac_norm * (pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1) * pow(h(0), -0.2e1) * (-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1)) * sqrt(0.3e1) + 0.2e1 / 0.3e1 * this->_fac_det * sqrt(0.3e1) * (x(1,0) - x(1,1)) * pow(h(1), -0.2e1) - 0.1e1 * this->_fac_det2 * pow(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1)) / 0.3e1, -0.2e1) * (0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(1,0) - x(1,1)) * pow(h(1), -0.2e1) + 0.4e1 * pow((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1), -0.1e1 / 0.2e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.4e1) * (x(1,0) - x(1,1))));
          grad(1,2) = DataType_(0.8e1 / 0.3e1 * this->_fac_norm * (pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1) * pow(h(0), -0.2e1) * (-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1)) * sqrt(0.3e1) + 0.2e1 / 0.3e1 * this->_fac_det * sqrt(0.3e1) * (-x(0,0) + x(0,1)) * pow(h(1), -0.2e1) - 0.1e1 * this->_fac_det2 * pow(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1)) / 0.3e1, -0.2e1) * (0.2e1 / 0.3e1 * sqrt(0.3e1) * (-x(0,0) + x(0,1)) * pow(h(1), -0.2e1) + 0.4e1 * pow((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(1), -0.4e1), -0.1e1 / 0.2e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(1), -0.4e1) * (-x(0,0) + x(0,1))));

          return;
        }

    }; // class RumpfFunctional
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_2D_P1_HPP
