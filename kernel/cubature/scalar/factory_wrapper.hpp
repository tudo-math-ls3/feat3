#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_FACTORY_WRAPPER_HPP
#define KERNEL_CUBATURE_SCALAR_FACTORY_WRAPPER_HPP 1

// includes, FEAST
#include <kernel/cubature/scalar/driver_factory.hpp>
#include <kernel/cubature/scalar/gauss_legendre_driver.hpp>
#include <kernel/cubature/scalar/pulcherima_driver.hpp>
#include <kernel/cubature/scalar/trapezoidal_driver.hpp>

namespace FEAST
{
  namespace Cubature
  {
    namespace Scalar
    {
      template<
        typename Weight_,
        typename Coord_>
      class FactoryWrapper
      {
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

      public:
        template<typename Functor_>
        static void driver(Functor_& functor)
        {
          // add scalar cubature drivers
          functor.template driver<GaussLegendreDriver>();
          functor.template driver<PulcherimaDriver>();
          functor.template driver<TrapezoidalDriver>();

          // TODO: add you new scalar cubature driver in the list above, e.g.
          // functor.template driver<YourDriverName>();
        }

        template<typename Functor_>
        static void factory(Functor_& functor)
        {
          // TODO: add you new scalar cubature driver in the list above, e.g.
          // functor.template factory<YourFactoryName>();

          // call driver factory functor
          DriverFactoryFunctor<Functor_> driver_functor(functor);
          driver(driver_functor);
        }
      };
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_SCALAR_DRIVER_WRAPPER_HPP
