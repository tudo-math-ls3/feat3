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
    /**
     * @brief Class template for Rumpf functionals, 2d P1
     *
     * \tparam MemoryType_
     * Memory architecture
     *
     * \tparam DataType_
     * Our data type
     *
     **/
    template<typename MemoryType_, typename DataType_>
    class RumpfFunctional
    <
      MemoryType_,
      DataType_,
      FEAST::Trafo::Standard::Mapping
      <
        Geometry::ConformalMesh
        <
          Shape::Simplex<2>,
          Shape::Simplex<2>::dimension,
          Shape::Simplex<2>::dimension,
          DataType_
        >
      >
    > : public RumpfFunctionalBase<MemoryType_, DataType_>
    {

      public:
        /// Our baseclass
        typedef RumpfFunctionalBase<MemoryType_, DataType_> BaseClass;

      public:
        /**
         * \copydoc RumpfFunctionalBase()
         *
         * Because we want the functional to have a minimum if the mesh consists of correctly scaled versions of the
         * optimal element (that is, every local transformation is the identity modulo translation), fac_det and
         * fac_det2 are now coupled.
         *
         * \f$ f(d) = d + \frac{c_s}{d^2 + \sqrt{ \mathrm{fac_reg}^2 + d^2}} \f$
         * We want to have \f$ f'(d) = 0 \Leftrightarrow c_s = \frac{1}{d \sqrt{ \mathrm{fac_reg}^2 + d^2} +
         * \mathrm{fac_reg}^2 + d^2) \f$ for \f$ d = 1 \f$.
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

      public:

        /**
         * \brief Computes value the Rumpf functional on one element.
         *
         * \param[in] x
         * (world_dim \f$ \times \f$ (number of vertices per cell)) - matrix that holds the local coordinates that
         * define the cell \f$ K \f$.
         *
         * \param[in] h
         * (world_dim) - vector that holds the dimensions of the target reference cell
         *
         * \returns
         * Local contribution to functional value
         *
         * Computes \f$ det(A), \| A \|^2_F \f$ but not \f$ \| Cof(A) \|^2_F \f$ for one cell defined by the coordinates x with regard
         * to the reference cell which depends on the shape type. The functional value is a function of these quantities,
         * in this case
         * \f$ f( det(A), \| A \|^2_F) = fac_norm (\| A \|^2_F - 2)^2 + fac_det det(A)^2 + fac_det2 \frac{1}{det(A)} \f$
         *
         **/
        template<typename Tx_, typename Th_>
        DataType_ compute_local_functional(const Tx_& x, const Th_& h)
        {
          DataType_ norm_A = compute_norm_A(x,h);
          DataType_ det_A = compute_det_A(x,h);
          DataType_ det2_A = compute_det2_A(x,h);

          return this->_fac_norm*norm_A
            + this->_fac_det*Math::sqr(det_A)
            + this->_fac_det2*det2_A;

        }

        /**
         * \copydoc compute_local_functional()
         *
         * \param[in] func_norm
         * Contribution from the Frobenius norm
         *
         * \param[in] func_det
         * Contribution from the det term
         *
         * \param[in] func_det2
         * Contribution from the 1/det term
         *
         * Debug variant that also returns the contributions of the different terms seperately.
         *
         **/
        template<typename Tx_, typename Th_>
        DataType_ compute_local_functional(const Tx_& x, const Th_& h,
        DataType_& func_norm,
        DataType_& func_det,
        DataType_& func_det2)
        {
          func_norm = this->_fac_norm*compute_norm_A(x,h);
          func_det = this->_fac_det*Math::sqr(compute_det_A(x,h));
          func_det2 = this->_fac_det2*compute_det2_A(x,h);

          return func_norm + func_det + func_det2;

        }

        /**
         * \brief Computes the det term on one element
         *
         * \param[in] x
         * (world_dim \f$ \times \f$ (number of vertices per cell)) - matrix that holds the local coordinates that
         * define the cell \f$ K \f$.
         *
         * \param[in] h
         * (world_dim) - vector that holds the dimensions of the target reference cell
         *
         * \returns
         * The local contribution to the 1/det term
         *
         * Computes \f$ \int_K \sqrt(det(A^T A)) dx \f$ for the cell defined by the coordinates x with regard to the reference
         * cell \f$ K^* \f$.
         *
         **/
        template<typename Tx_, typename Th_ >
        DataType_ compute_det_A( const Tx_& x, const Th_& h)
        {
          return DataType_(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.2e1));
        }

        /**
         * \brief Computes the 1/det term on one element
         *
         * \param[in] x
         * (world_dim \f$ \times \f$ (number of vertices per cell)) - matrix that holds the local coordinates that
         * define the cell \f$ K \f$.
         *
         * \param[in] h
         * (world_dim) - vector that holds the dimensions of the target reference cell
         *
         * \returns
         * The local contribution to the 1/det term
         *
         * Computes \f$ \int_K \sqrt{\frac{1}{det(A^T A)}} dx \f$ for the cell defined by the coordinates x with regard to the reference
         * cell \f$ K^* \f$.
         *
         **/
        template<typename Tx_, typename Th_ >
        DataType_ compute_det2_A( const Tx_& x, const Th_& h)
        {
          return DataType_(0.1e1 / (0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1)) / 0.3e1));
        }

        /**
         * \brief Computes \f$ \int_K \| A \|^2_F  dx \f$ for the cell defined by the coordinates x with regard to the reference
         * cell \f$ K^* \f$.
         *
         * \param[in] x
         * (world_dim \f$ \times \f$ (number of vertices per cell)) - matrix that holds the local coordinates that
         * define the cell \f$ K \f$.
         *
         * \param[in] h
         * (world_dim) - vector that holds the dimensions of the target reference cell
         *
         * \returns
         * The local contribution to the Frobenius norm term
         *
         **/
        template<typename Tx_, typename Th_ >
        DataType_ compute_norm_A(const Tx_& x, const Th_& h)
        {

          return DataType_(pow(pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1, 0.2e1));

        }

        /**
         * \brief Computes the functional gradient for the cell defined by the coordinates x with regard to the reference
         * cell \f$ K^* \f$.
         *
         * \param[in] x
         * (world_dim \f$ \times \f$ (number of vertices per cell)) - matrix that holds the local coordinates that
         * define the cell \f$ K \f$.
         *
         * \param[in] h
         * (world_dim) - vector that holds the dimensions of the target reference cell
         *
         * \param[in] grad
         * (world_dim \f$ \times \f$ (number of vertices per cell)) - matrix that holds the local contribution to the
         * global functional gradient
         *
         * \returns
         * The local contribution grad to global functional gradient
         *
         **/
        template<typename Tx_, typename Th_, typename Tgrad_>
        void compute_local_grad( const Tx_& x, const Th_& h, Tgrad_& grad)
        {

          grad(0,0) = DataType(0.2e1 * this->_fac_norm * (pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1) * (-0.2e1 * pow(h(0), -0.2e1) * (-x(0,0) + x(0,1)) - 0.2e1 / 0.3e1 * pow(h(0), -0.2e1) * (-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1)) * sqrt(0.3e1)) + 0.2e1 / 0.3e1 * this->_fac_det * sqrt(0.3e1) * (x(1,1) - x(1,2)) * pow(h(0), -0.2e1) - 0.1e1 * this->_fac_det2 * pow(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1)) / 0.3e1, -0.2e1) * (0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(1,1) - x(1,2)) * pow(h(0), -0.2e1) + 0.4e1 * pow((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1), -0.1e1 / 0.2e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.4e1) * (x(1,1) - x(1,2))));

          grad(0,1) = DataType(0.2e1 * this->_fac_norm * (pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1) * (0.2e1 * pow(h(0), -0.2e1) * (-x(0,0) + x(0,1)) - 0.2e1 / 0.3e1 * pow(h(0), -0.2e1) * (-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1)) * sqrt(0.3e1)) + 0.2e1 / 0.3e1 * this->_fac_det * sqrt(0.3e1) * (-x(1,0) + x(1,2)) * pow(h(0), -0.2e1) - 0.1e1 * this->_fac_det2 * pow(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1)) / 0.3e1, -0.2e1) * (0.2e1 / 0.3e1 * sqrt(0.3e1) * (-x(1,0) + x(1,2)) * pow(h(0), -0.2e1) + 0.4e1 * pow((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1), -0.1e1 / 0.2e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.4e1) * (-x(1,0) + x(1,2))));

          grad(1,0) = DataType(0.2e1 * this->_fac_norm * (pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1) * (-0.2e1 * pow(h(0), -0.2e1) * (-x(1,0) + x(1,1)) - 0.2e1 / 0.3e1 * pow(h(0), -0.2e1) * (-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1)) * sqrt(0.3e1)) + 0.2e1 / 0.3e1 * this->_fac_det * sqrt(0.3e1) * (-x(0,1) + x(0,2)) * pow(h(0), -0.2e1) - 0.1e1 * this->_fac_det2 * pow(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1)) / 0.3e1, -0.2e1) * (0.2e1 / 0.3e1 * sqrt(0.3e1) * (-x(0,1) + x(0,2)) * pow(h(0), -0.2e1) + 0.4e1 * pow((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1), -0.1e1 / 0.2e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.4e1) * (-x(0,1) + x(0,2))));

          grad(1,1) = DataType(0.2e1 * this->_fac_norm * (pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1) * (0.2e1 * pow(h(0), -0.2e1) * (-x(1,0) + x(1,1)) - 0.2e1 / 0.3e1 * pow(h(0), -0.2e1) * (-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1)) * sqrt(0.3e1)) + 0.2e1 / 0.3e1 * this->_fac_det * sqrt(0.3e1) * (x(0,0) - x(0,2)) * pow(h(0), -0.2e1) - 0.1e1 * this->_fac_det2 * pow(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1)) / 0.3e1, -0.2e1) * (0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) - x(0,2)) * pow(h(0), -0.2e1) + 0.4e1 * pow((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1), -0.1e1 / 0.2e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.4e1) * (x(0,0) - x(0,2))));

          grad(0,2) = DataType(0.8e1 / 0.3e1 * this->_fac_norm * (pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1) * pow(h(0), -0.2e1) * (-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1)) * sqrt(0.3e1) + 0.2e1 / 0.3e1 * this->_fac_det * sqrt(0.3e1) * (x(1,0) - x(1,1)) * pow(h(0), -0.2e1) - 0.1e1 * this->_fac_det2 * pow(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1)) / 0.3e1, -0.2e1) * (0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(1,0) - x(1,1)) * pow(h(0), -0.2e1) + 0.4e1 * pow((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1), -0.1e1 / 0.2e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.4e1) * (x(1,0) - x(1,1))));

          grad(1,2) = DataType(0.8e1 / 0.3e1 * this->_fac_norm * (pow(h(0), -0.2e1) * pow(-x(0,0) + x(0,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(0,0) + x(0,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(0,0) + x(0,2)) * sqrt(0.3e1), 0.2e1) + pow(h(0), -0.2e1) * pow(-x(1,0) + x(1,1), 0.2e1) + pow(h(0), -0.2e1) * pow(-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1), 0.2e1) - 0.2e1) * pow(h(0), -0.2e1) * (-(-x(1,0) + x(1,1)) * sqrt(0.3e1) / 0.3e1 + 0.2e1 / 0.3e1 * (-x(1,0) + x(1,2)) * sqrt(0.3e1)) * sqrt(0.3e1) + 0.2e1 / 0.3e1 * this->_fac_det * sqrt(0.3e1) * (-x(0,0) + x(0,1)) * pow(h(0), -0.2e1) - 0.1e1 * this->_fac_det2 * pow(0.2e1 / 0.3e1 * sqrt(0.3e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.2e1) + sqrt((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1)) / 0.3e1, -0.2e1) * (0.2e1 / 0.3e1 * sqrt(0.3e1) * (-x(0,0) + x(0,1)) * pow(h(0), -0.2e1) + 0.4e1 * pow((double) (9 * this->_fac_reg * this->_fac_reg) + 0.12e2 * pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), 0.2e1) * pow(h(0), -0.4e1), -0.1e1 / 0.2e1) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * pow(h(0), -0.4e1) * (-x(0,0) + x(0,1))));

          return;
        }

    }; // class RumpfFunctional< DataType_, FEAST::Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>>> >
    } // namespace Geometry
  } // namespace FEAST

#endif // KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_2D_P1_HPP
