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
        typename Policy_,
        typename ScalarFactory_,
        bool tensorise_ = (ScalarFactory_::tensorise != 0)>
      class TensorProductFunctorHelper
      {
      public:
        template<typename Functor_>
        static void factory(Functor_& functor)
        {
          TensorProductDriver<Policy_>::template factory<ScalarFactory_>(functor);
        }
      };

      template<
        typename Policy_,
        typename ScalarFactory_>
      class TensorProductFunctorHelper<Policy_, ScalarFactory_, false>
      {
      public:
        template<typename Functor_>
        static void factory(Functor_&)
        {
          // do nothing
        }
      };
    } // namespace Intern
    /// \endcond

    template<typename Policy_>
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
          _functor.template factory< DriverFactory<Driver_, Policy_> >();
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

        template<typename ScalarFactory_>
        void factory()
        {
          Intern::TensorProductFunctorHelper<Policy_, ScalarFactory_>::factory(_functor);
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
        typedef Rule<Policy_> RuleType;
        typedef typename RuleType::WeightType WeightType;
        typedef typename RuleType::CoordType CoordType;
        typedef Scalar::FactoryWrapper<WeightType, CoordType> ScalarFactoryWrapper;
        // call the scalar factory wrapper
        ScalarFactoryWrapper::factory(functor);
      }
    }; // class FactoryWrapper<...>
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_FACTORY_WRAPPER_HPP
