#pragma once
#ifndef KERNEL_MESHOPT_RUMPF_FUNCTIONALS_Q1HACK_HPP
#define KERNEL_MESHOPT_RUMPF_FUNCTIONALS_Q1HACK_HPP 1

#include <kernel/base_header.hpp>

namespace FEAST
{
  namespace Meshopt
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
         */
        explicit RumpfFunctionalQ1Hack(DataType fac_norm_, DataType fac_det_, DataType fac_cof_, DataType fac_reg_) :
          BaseClass(fac_norm_,
          fac_det_,
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
          return "RumpfFunctionalQ1Hack<"+ShapeType::name()+", "+BaseClass::name()+">";
        }

        /**
         * \brief Prints object parameters
         */
        void print()
        {
          std::cout << name() << std::endl;
          BaseClass::print();
        }


    }; // class RumpfFunctionalQ1Hack
    /// \endcond
  } // namespace Meshopt
} // namespace FEAST

#endif // KERNEL_MESHOPT_RUMPF_FUNCTIONALS_Q1HACK_HPP
