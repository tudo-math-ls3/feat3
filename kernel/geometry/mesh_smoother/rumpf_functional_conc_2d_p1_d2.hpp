#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_CONC_2D_P1_D2_HPP
#define KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_CONC_2D_P1_D2_HPP 1

#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_p1_d2.hpp>
// For HeavisideReg
#include <kernel/assembly/common_functions.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal

    /**
     * \brief Class template for Rumpf functionals with levelset term, 2d P1
     **/
    template<typename DataType_>
    class RumpfFunctionalConc_D2<DataType_, Shape::Simplex<2> > :
    public RumpfFunctional_D2<DataType_, Shape::Simplex<2> >
    {
      public:
        /// Our data Type
        typedef DataType_ DataType;
        /// Shape type of the underlying transformation
        typedef Shape::Simplex<2> ShapeType;
        /// Our base class
        typedef RumpfFunctional_D2<DataType, ShapeType> BaseClass;

      public:
        /**
         * \brief Constructor
         **/
        RumpfFunctionalConc_D2(
          const DataType fac_norm_,
          const DataType fac_det_,
          const DataType fac_cof_,
          const DataType fac_reg_) :
          BaseClass( fac_norm_, fac_det_, fac_cof_, fac_reg_)
          {
          }

        /**
         * \brief
         **/
        template<typename Tgrad_, typename Tx_, typename Th_, typename Tgradh_>
        void add_grad_h_part(Tgrad_& grad, const Tx_& x, const Th_& h, const Tgradh_& grad_h)
        {
          DataType der_h_(0);
          der_h_ = this->_fac_norm * (DataType(16) / DataType(3) * (-DataType(2) * Math::pow(x(0,0), DataType(2)) + DataType(2) * x(0,0) * x(0,1) - DataType(2) * Math::pow(x(0,1), DataType(2)) + DataType(2) * x(0,0) * x(0,2) + DataType(2) * x(0,1) * x(0,2) - DataType(2) * Math::pow(x(0,2), DataType(2)) - DataType(2) * Math::pow(x(1,0), DataType(2)) + DataType(2) * x(1,0) * x(1,1) - DataType(2) * Math::pow(x(1,1), DataType(2)) + DataType(2) * x(1,0) * x(1,2) + DataType(2) * x(1,1) * x(1,2) - DataType(2) * Math::pow(x(1,2), DataType(2)) + DataType(3) * Math::pow(h(0), DataType(2))) * Math::pow(h(0), -DataType(3)) - DataType(16) / DataType(9) * Math::pow(-DataType(2) * Math::pow(x(0,0), DataType(2)) + DataType(2) * x(0,0) * x(0,1) - DataType(2) * Math::pow(x(0,1), DataType(2)) + DataType(2) * x(0,0) * x(0,2) + DataType(2) * x(0,1) * x(0,2) - DataType(2) * Math::pow(x(0,2), DataType(2)) - DataType(2) * Math::pow(x(1,0), DataType(2)) + DataType(2) * x(1,0) * x(1,1) - DataType(2) * Math::pow(x(1,1), DataType(2)) + DataType(2) * x(1,0) * x(1,2) + DataType(2) * x(1,1) * x(1,2) - DataType(2) * Math::pow(x(1,2), DataType(2)) + DataType(3) * Math::pow(h(0), DataType(2)), DataType(2)) * Math::pow(h(0), -DataType(5))) - DataType(16) / DataType(3) * this->_fac_det * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(2)) * Math::pow(h(1), -DataType(5)) - DataType(2) * this->_fac_rec_det * Math::pow(DataType(2) / DataType(3) * Math::sqrt(DataType(3)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(2)) + Math::sqrt(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(2)) * Math::pow(h(1), -DataType(4))) / DataType(3), -DataType(3)) * (-DataType(4) / DataType(3) * Math::sqrt(DataType(3)) * (x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1)) * Math::pow(h(1), -DataType(3)) - DataType(8) * Math::pow(DataType(9) * this->_fac_reg * this->_fac_reg + DataType(12) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(2)) * Math::pow(h(1), -DataType(4)), -DataType(1) / DataType(2)) * Math::pow(x(0,0) * x(1,1) - x(0,0) * x(1,2) - x(0,1) * x(1,0) + x(0,1) * x(1,2) + x(0,2) * x(1,0) - x(0,2) * x(1,1), DataType(2)) * Math::pow(h(1), -DataType(5)));

          for(Index d(0); d < Tgrad_::m; ++d)
          {
            for(Index i(0); i < Tgrad_::n; ++i)
            {
              grad(d,i) += der_h_*grad_h(d*Tgrad_::n+i);
            }
          }
        } // add_grad_h_part

    }; // class RumpfFunctionalConc_D2
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_CONC_2D_P1_D2_HPP
