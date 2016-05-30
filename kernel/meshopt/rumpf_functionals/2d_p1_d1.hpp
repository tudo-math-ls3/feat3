#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_FUNCTIONALS_2D_P1_D1_HPP
#define KERNEL_MESHOPT_RUMPF_FUNCTIONALS_2D_P1_D1_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/shape.hpp>
#include <kernel/meshopt//rumpf_functional.hpp>

namespace FEAT
{
  namespace Meshopt
  {
    /// \cond internal

    /**
     * @brief Rumpf functional, 2d P1, linear det term
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
         */
        explicit RumpfFunctional(
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
         * \brief The class name
         *
         * \returns String with the class name
         */
        static String name()
        {
          return "RumpfFunctional<"+ShapeType::name()+">";
        }

        /**
         * \brief Prints object parameters
         */
        void print()
        {
          std::cout << name() << std::endl;
          BaseClass::print();
        }

        /**
         * \brief Computes value the Rumpf functional on one element.
         */
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
        DataType NOINLINE compute_det_A( const Tx_& x, const Th_& h)
        {
          DataType det_;
          det_ = DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(2));
          return det_;
        }

        /**
         * \brief Computes the 1/det term on one element
         **/
        template<typename Tx_, typename Th_ >
        DataType NOINLINE compute_rec_det_A( const Tx_& x, const Th_& h)
        {
          DataType rec_det_;
          rec_det_ = Math::pow(DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(2)) + Math::sqrt(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4))) / DataType(3), -DataType(2));
          return rec_det_;
        }

        /**
         * \brief Computes the Frobenius norm term for one cell
         **/
        template<typename Tx_, typename Th_ >
        DataType NOINLINE compute_norm_A(const Tx_& x, const Th_& h)
        {
          DataType norm_;
          norm_ = DataType(4) / DataType(9) * Math::pow(DataType(3) * Math::pow(h(0), DataType(2)) - DataType(2) * Math::pow(x(0,0), DataType(2)) + DataType(2) * x(0,0) * x(1,0) + DataType(2) * x(0,0) * x(2,0) - DataType(2) * Math::pow(x(0,1), DataType(2)) + DataType(2) * x(0,1) * x(1,1) + DataType(2) * x(0,1) * x(2,1) - DataType(2) * Math::pow(x(1,0), DataType(2)) + DataType(2) * x(1,0) * x(2,0) - DataType(2) * Math::pow(x(1,1), DataType(2)) + DataType(2) * x(1,1) * x(2,1) - DataType(2) * Math::pow(x(2,0), DataType(2)) - DataType(2) * Math::pow(x(2,1), DataType(2)), DataType(2)) * Math::pow(h(0), -DataType(4));
          return norm_;
        }

        /**
         * \brief Computes the functional gradient for one cell
         **/
        template<typename Tx_, typename Th_, typename Tgrad_>
        void NOINLINE compute_local_grad( const Tx_& x, const Th_& h, Tgrad_& grad)
        {
          grad(0,0) = DataType(2) / DataType(3) * this->_fac_det * Math::sqrt(DataType(3)) * (x(1,1) - x(2,1)) * Math::pow(h(1), -DataType(2)) + DataType(8) / DataType(9) * this->_fac_norm * (DataType(3) * Math::pow(h(0), DataType(2)) - DataType(2) * Math::pow(x(0,0), DataType(2)) + DataType(2) * x(0,0) * x(1,0) + DataType(2) * x(0,0) * x(2,0) - DataType(2) * Math::pow(x(0,1), DataType(2)) + DataType(2) * x(0,1) * x(1,1) + DataType(2) * x(0,1) * x(2,1) - DataType(2) * Math::pow(x(1,0), DataType(2)) + DataType(2) * x(1,0) * x(2,0) - DataType(2) * Math::pow(x(1,1), DataType(2)) + DataType(2) * x(1,1) * x(2,1) - DataType(2) * Math::pow(x(2,0), DataType(2)) - DataType(2) * Math::pow(x(2,1), DataType(2))) * Math::pow(h(0), -DataType(4)) * (-DataType(4) * x(0,0) + DataType(2) * x(1,0) + DataType(2) * x(2,0)) - DataType(2) * this->_fac_rec_det * Math::pow(DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(2)) + Math::sqrt(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4))) / DataType(3), -DataType(3)) * (DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (x(1,1) - x(2,1)) * Math::pow(h(1), -DataType(2)) + DataType(4) * Math::pow(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4)), -DataType(1) / DataType(2)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(4)) * (x(1,1) - x(2,1)));
          grad(0,1) = DataType(2) / DataType(3) * this->_fac_det * Math::sqrt(DataType(3)) * (-x(1,0) + x(2,0)) * Math::pow(h(1), -DataType(2)) + DataType(8) / DataType(9) * this->_fac_norm * (DataType(3) * Math::pow(h(0), DataType(2)) - DataType(2) * Math::pow(x(0,0), DataType(2)) + DataType(2) * x(0,0) * x(1,0) + DataType(2) * x(0,0) * x(2,0) - DataType(2) * Math::pow(x(0,1), DataType(2)) + DataType(2) * x(0,1) * x(1,1) + DataType(2) * x(0,1) * x(2,1) - DataType(2) * Math::pow(x(1,0), DataType(2)) + DataType(2) * x(1,0) * x(2,0) - DataType(2) * Math::pow(x(1,1), DataType(2)) + DataType(2) * x(1,1) * x(2,1) - DataType(2) * Math::pow(x(2,0), DataType(2)) - DataType(2) * Math::pow(x(2,1), DataType(2))) * Math::pow(h(0), -DataType(4)) * (-DataType(4) * x(0,1) + DataType(2) * x(1,1) + DataType(2) * x(2,1)) - DataType(2) * this->_fac_rec_det * Math::pow(DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(2)) + Math::sqrt(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4))) / DataType(3), -DataType(3)) * (DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (-x(1,0) + x(2,0)) * Math::pow(h(1), -DataType(2)) + DataType(4) * Math::pow(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4)), -DataType(1) / DataType(2)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(4)) * (-x(1,0) + x(2,0)));
          grad(1,0) = DataType(2) / DataType(3) * this->_fac_det * Math::sqrt(DataType(3)) * (-x(0,1) + x(2,1)) * Math::pow(h(1), -DataType(2)) + DataType(8) / DataType(9) * this->_fac_norm * (DataType(3) * Math::pow(h(0), DataType(2)) - DataType(2) * Math::pow(x(0,0), DataType(2)) + DataType(2) * x(0,0) * x(1,0) + DataType(2) * x(0,0) * x(2,0) - DataType(2) * Math::pow(x(0,1), DataType(2)) + DataType(2) * x(0,1) * x(1,1) + DataType(2) * x(0,1) * x(2,1) - DataType(2) * Math::pow(x(1,0), DataType(2)) + DataType(2) * x(1,0) * x(2,0) - DataType(2) * Math::pow(x(1,1), DataType(2)) + DataType(2) * x(1,1) * x(2,1) - DataType(2) * Math::pow(x(2,0), DataType(2)) - DataType(2) * Math::pow(x(2,1), DataType(2))) * Math::pow(h(0), -DataType(4)) * (DataType(2) * x(0,0) - DataType(4) * x(1,0) + DataType(2) * x(2,0)) - DataType(2) * this->_fac_rec_det * Math::pow(DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(2)) + Math::sqrt(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4))) / DataType(3), -DataType(3)) * (DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (-x(0,1) + x(2,1)) * Math::pow(h(1), -DataType(2)) + DataType(4) * Math::pow(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4)), -DataType(1) / DataType(2)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(4)) * (-x(0,1) + x(2,1)));
          grad(1,1) = DataType(2) / DataType(3) * this->_fac_det * Math::sqrt(DataType(3)) * (x(0,0) - x(2,0)) * Math::pow(h(1), -DataType(2)) + DataType(8) / DataType(9) * this->_fac_norm * (DataType(3) * Math::pow(h(0), DataType(2)) - DataType(2) * Math::pow(x(0,0), DataType(2)) + DataType(2) * x(0,0) * x(1,0) + DataType(2) * x(0,0) * x(2,0) - DataType(2) * Math::pow(x(0,1), DataType(2)) + DataType(2) * x(0,1) * x(1,1) + DataType(2) * x(0,1) * x(2,1) - DataType(2) * Math::pow(x(1,0), DataType(2)) + DataType(2) * x(1,0) * x(2,0) - DataType(2) * Math::pow(x(1,1), DataType(2)) + DataType(2) * x(1,1) * x(2,1) - DataType(2) * Math::pow(x(2,0), DataType(2)) - DataType(2) * Math::pow(x(2,1), DataType(2))) * Math::pow(h(0), -DataType(4)) * (DataType(2) * x(0,1) - DataType(4) * x(1,1) + DataType(2) * x(2,1)) - DataType(2) * this->_fac_rec_det * Math::pow(DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(2)) + Math::sqrt(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4))) / DataType(3), -DataType(3)) * (DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (x(0,0) - x(2,0)) * Math::pow(h(1), -DataType(2)) + DataType(4) * Math::pow(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4)), -DataType(1) / DataType(2)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(4)) * (x(0,0) - x(2,0)));
          grad(2,0) = DataType(2) / DataType(3) * this->_fac_det * Math::sqrt(DataType(3)) * (x(0,1) - x(1,1)) * Math::pow(h(1), -DataType(2)) + DataType(8) / DataType(9) * this->_fac_norm * (DataType(3) * Math::pow(h(0), DataType(2)) - DataType(2) * Math::pow(x(0,0), DataType(2)) + DataType(2) * x(0,0) * x(1,0) + DataType(2) * x(0,0) * x(2,0) - DataType(2) * Math::pow(x(0,1), DataType(2)) + DataType(2) * x(0,1) * x(1,1) + DataType(2) * x(0,1) * x(2,1) - DataType(2) * Math::pow(x(1,0), DataType(2)) + DataType(2) * x(1,0) * x(2,0) - DataType(2) * Math::pow(x(1,1), DataType(2)) + DataType(2) * x(1,1) * x(2,1) - DataType(2) * Math::pow(x(2,0), DataType(2)) - DataType(2) * Math::pow(x(2,1), DataType(2))) * Math::pow(h(0), -DataType(4)) * (DataType(2) * x(0,0) + DataType(2) * x(1,0) - DataType(4) * x(2,0)) - DataType(2) * this->_fac_rec_det * Math::pow(DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(2)) + Math::sqrt(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4))) / DataType(3), -DataType(3)) * (DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (x(0,1) - x(1,1)) * Math::pow(h(1), -DataType(2)) + DataType(4) * Math::pow(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4)), -DataType(1) / DataType(2)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(4)) * (x(0,1) - x(1,1)));
          grad(2,1) = DataType(2) / DataType(3) * this->_fac_det * Math::sqrt(DataType(3)) * (-x(0,0) + x(1,0)) * Math::pow(h(1), -DataType(2)) + DataType(8) / DataType(9) * this->_fac_norm * (DataType(3) * Math::pow(h(0), DataType(2)) - DataType(2) * Math::pow(x(0,0), DataType(2)) + DataType(2) * x(0,0) * x(1,0) + DataType(2) * x(0,0) * x(2,0) - DataType(2) * Math::pow(x(0,1), DataType(2)) + DataType(2) * x(0,1) * x(1,1) + DataType(2) * x(0,1) * x(2,1) - DataType(2) * Math::pow(x(1,0), DataType(2)) + DataType(2) * x(1,0) * x(2,0) - DataType(2) * Math::pow(x(1,1), DataType(2)) + DataType(2) * x(1,1) * x(2,1) - DataType(2) * Math::pow(x(2,0), DataType(2)) - DataType(2) * Math::pow(x(2,1), DataType(2))) * Math::pow(h(0), -DataType(4)) * (DataType(2) * x(0,1) + DataType(2) * x(1,1) - DataType(4) * x(2,1)) - DataType(2) * this->_fac_rec_det * Math::pow(DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(2)) + Math::sqrt(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4))) / DataType(3), -DataType(3)) * (DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (-x(0,0) + x(1,0)) * Math::pow(h(1), -DataType(2)) + DataType(4) * Math::pow(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0), DataType(2)) * Math::pow(h(1), -DataType(4)), -DataType(1) / DataType(2)) * (x(0,0) * x(1,1) - x(0,0) * x(2,1) - x(1,0) * x(0,1) + x(0,1) * x(2,0) + x(1,0) * x(2,1) - x(1,1) * x(2,0)) * Math::pow(h(1), -DataType(4)) * (-x(0,0) + x(1,0)));

          return;
        }

        // This is for the RumpfSmootherSplit variants
        ///**
        // * \brief Computes the gradient of the Frobenius norm term for one cell
        // **/
        //template<typename Tx_, typename Th_, typename Tgrad_>
        //void compute_grad_norm( const Tx_& x, const Th_& h, Tgrad_& grad_norm)
        //{
        //  grad_norm(0,0) =
        //    grad_norm(0,1) =
        //    grad_norm(1,0) =
        //    grad_norm(1,1) =
        //    grad_norm(2,0) =
        //    grad_norm(2,1) =
        //}

    }; // class RumpfFunctional
    /// \endcond
  } // namespace Meshopt
} // namespace FEAT

#endif // KERNEL_MESHOPT_RUMPF_FUNCTIONALS_2D_P1_D1_HPP
