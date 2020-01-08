// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_CUBATURE_SCALAR_FACTORY_WRAPPER_HPP
#define KERNEL_CUBATURE_SCALAR_FACTORY_WRAPPER_HPP 1

// includes, FEAT
#include <kernel/cubature/scalar/driver_factory.hpp>
#include <kernel/cubature/scalar/gauss_legendre_driver.hpp>
#include <kernel/cubature/scalar/gauss_lobatto_driver.hpp>
#include <kernel/cubature/scalar/maclaurin_driver.hpp>
#include <kernel/cubature/scalar/midpoint_driver.hpp>
#include <kernel/cubature/scalar/newton_cotes_closed_driver.hpp>
#include <kernel/cubature/scalar/newton_cotes_open_driver.hpp>
#include <kernel/cubature/scalar/trapezoidal_driver.hpp>

namespace FEAT
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
      class FactoryWrapper
      {
      protected:
        template<typename Functor_>
        static void _driver_list(Functor_& functor)
        {
          // >>> CUBATURE DRIVER LIST >>>
          // TODO: add you new scalar cubature driver at the end of the list below, e.g.
          // functor.template driver<YourDriverName>();
          functor.template driver<GaussLegendreDriver>();
          functor.template driver<GaussLobattoDriver>();
          functor.template driver<MaclaurinDriver>();
          functor.template driver<MidpointDriver>();
          functor.template driver<NewtonCotesClosedDriver>();
          functor.template driver<NewtonCotesOpenDriver>();
          functor.template driver<TrapezoidalDriver>();

          // <<< END OF CUBATURE DRIVER LIST <<<
        }

        template<typename Functor_>
        static void _factory_list(Functor_& /*functor*/)
        {
          // >>> CUBATURE FACTORY LIST >>>
          // TODO: add you new scalar cubature factory at the end of the list below, e.g.
          // functor.template factory<YourFactoryName>();

          // <<< END OF CUBATURE FACTORY LIST <<<
        }

      public:
        template<typename Functor_>
        static void driver(Functor_& functor)
        {
          // call driver list
          _driver_list(functor);
        }

        template<typename Functor_>
        static void factory(Functor_& functor)
        {
          // call factory list
          _factory_list(functor);

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

          template<typename Driver_>
          void driver()
          {
            _functor.template factory< DriverFactory<Driver_> >();
          }
        };
        /// \endcond
      };
    } // namespace Scalar
  } // namespace Cubature
} // namespace FEAT

#endif // KERNEL_CUBATURE_SCALAR_DRIVER_WRAPPER_HPP
