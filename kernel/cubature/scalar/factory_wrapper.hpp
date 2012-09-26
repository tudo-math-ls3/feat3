#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_FACTORY_WRAPPER_HPP
#define KERNEL_CUBATURE_SCALAR_FACTORY_WRAPPER_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/driver_factory.hpp>
#include <kernel/cubature/scalar/gauss_legendre_driver.hpp>
#include <kernel/cubature/scalar/pulcherrima_driver.hpp>
#include <kernel/cubature/scalar/trapezoidal_driver.hpp>

namespace FEAST
{
  namespace Cubature
  {
    namespace Scalar
    {
      /**
       * \brief Scalar Cubature Factory Wrapper class template
       *
       * \author Peter Zajac
       */
      template<
        typename Weight_,
        typename Coord_>
      class FactoryWrapper
      {
      public:
        template<typename Functor_>
        static void driver(Functor_& functor)
        {
          // >>> CUBATURE DRIVER LIST >>>
          // TODO: add you new scalar cubature driver at the end of the list below, e.g.
          // functor.template driver<YourDriverName>();
          functor.template driver<GaussLegendreDriver>();
          functor.template driver<PulcherrimaDriver>();
          functor.template driver<TrapezoidalDriver>();

          // <<< END OF CUBATURE DRIVER LIST <<<
        }

        template<typename Functor_>
        static void factory(Functor_& functor)
        {
          // >>> CUBATURE FACTORY LIST >>>
          // TODO: add you new scalar cubature factory at the end of the list below, e.g.
          // functor.template factory<YourFactoryName>();

          // <<< END OF CUBATURE FACTORY LIST <<<

          // call driver factory functor
          DriverFactoryFunctor<Functor_> driver_functor(functor);
          driver(driver_functor);
        }

        /// \cond internal
      private:
        template<typename Functor_>
        class DriverFactoryFunctor
        {
        protected:
          Functor_& _functor;

        public:
          explicit DriverFactoryFunctor(Functor_& functor) :
            _functor(functor)
          {
          }

          template<template<typename,typename> class Driver_>
          void driver()
          {
            _functor.template factory< DriverFactory<Driver_, Weight_, Coord_> >();
          }
        };
        /// \endcond
      };
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_DRIVER_WRAPPER_HPP
