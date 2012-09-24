#pragma once
#ifndef KERNEL_CUBATURE_FACTORY_WRAPPER_HPP
#define KERNEL_CUBATURE_FACTORY_WRAPPER_HPP 1

// includes, FEAST
#include <kernel/cubature/driver_factory.hpp>
#include <kernel/cubature/barycentre_driver.hpp>
#include <kernel/cubature/trapezoidal_driver.hpp>
#include <kernel/cubature/tensor_product_factory.hpp>
#include <kernel/cubature/scalar/factory_wrapper.hpp>

namespace FEAST
{
  namespace Cubature
  {
    /// \cond internal
    namespace Intern
    {
      template<
        template<typename,typename> class ScalarDriver_,
        typename Shape_,
        typename Weight_,
        typename Coord_,
        typename Point_,
        bool tensorise_ = (ScalarDriver_<Weight_, Coord_>::tensorise != 0)>
      class TensorProductFunctorHelper
      {
      public:
        template<typename Functor_>
        static void scalar_driver(Functor_& functor)
        {
          TensorProductDriver<Shape_, Weight_, Coord_, Point_>::template scalar_driver<ScalarDriver_>(functor);
        }
      };

      template<
        template<typename,typename> class ScalarDriver_,
        typename Shape_,
        typename Weight_,
        typename Coord_,
        typename Point_>
      class TensorProductFunctorHelper<ScalarDriver_, Shape_, Weight_, Coord_, Point_, false>
      {
      public:
        template<typename Functor_>
        static void scalar_driver(Functor_&)
        {
          // do nothing
        }
      };
    } // namespace Intern
    /// \endcond

    template<
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
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

        template<template<typename,typename,typename,typename> class Driver_>
        void driver()
        {
          _functor.template factory< DriverFactory<Driver_, Shape_, Weight_, Coord_, Point_> >();
        }
      };

      template<typename Functor_>
      class TensorProductFunctor
      {
      protected:
        Functor_& _functor;

      public:
        explicit TensorProductFunctor(Functor_& functor) :
          _functor(functor)
        {
        }

        template<template<typename,typename> class ScalarDriver_>
        void driver()
        {
          Intern::TensorProductFunctorHelper<ScalarDriver_, Shape_, Weight_, Coord_, Point_>::scalar_driver(_functor);
        }
      };

    public:
      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // add cubature drivers
        functor.template driver<BarycentreDriver>();
        functor.template driver<TrapezoidalDriver>();

        // TODO: add you new cubature driver in the list above, e.g.
        // functor.template driver<YourDriverName>();
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // TODO: add you new cubature driver in the list above, e.g.
        // functor.template factory<YourFactoryName>();

        // call driver factory functor
        DriverFactoryFunctor<Functor_> driver_functor(functor);
        driver(driver_functor);

        // call tensor-product functor
        TensorProductFunctor<Functor_> tensor_functor(functor);
        tensor(tensor_functor);
      }

    protected:
      template<typename Functor_>
      static void tensor(TensorProductFunctor<Functor_>& functor)
      {
        // call the scalar factory wrapper
        Scalar::FactoryWrapper<Weight_, Coord_>::driver(functor);
      }
    }; // class FactoryWrapper<...>
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_FACTORY_WRAPPER_HPP
