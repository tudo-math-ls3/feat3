#pragma once
#ifndef KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_2D_Q1HACK_HPP
#define KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_2D_Q1HACK_HPP 1


#include <kernel/geometry/mesh_smoother/rumpf_functional.hpp>
// For HeavisideReg
#include <kernel/assembly/common_functions.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal

    /**
     * \brief Class template for Rumpf functionals in 2d using the Q1 hack
     **/
    template<typename DataType_, template<typename, typename> class BaseClass_>
    class RumpfFunctionalQ1Hack<DataType_, Shape::Hypercube<2>, BaseClass_> :
    public BaseClass_<DataType_, Shape::Simplex<2>>
    {
      public:
        /// Our data type
        typedef DataType_ DataType;
        /// Our shape type
        typedef Shape::Hypercube<2> ShapeType;
        /// Our base class
        typedef BaseClass_<DataType_, Shape::Simplex<2>> BaseClass;

        /**
         * \brief Constructor
         **/
        RumpfFunctionalQ1Hack( DataType fac_norm_, DataType fac_det_, DataType fac_cof_, DataType fac_reg_) :
          BaseClass( fac_norm_,
          fac_det_,
          fac_cof_,
          fac_reg_)
          {
          }

    }; // class RumpfFunctionalQ1Hack
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_RUMPF_FUNCTIONAL_LVLSET_2D_Q1HACK_HPP
