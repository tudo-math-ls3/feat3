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
        typename Shape_,
        typename Weight_,
        typename Coord_,
        typename Point_,
        typename Functor_>
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

      template<
        typename Shape_,
        typename Weight_,
        typename Coord_,
        typename Point_,
        typename Functor_>
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
    } // namespace Intern
    /// \endcond

    /**
     * \brief Explicitly Specialised Cubature Factory Wrapper class template
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class FactoryExplicitWrapper;

    /**
     * \brief Explicit specialisation for Simplex<1> shape
     *
     * \author Peter Zajac
     */
    template<
      typename Weight_,
      typename Coord_,
      typename Point_>
    class FactoryExplicitWrapper<Shape::Simplex<1>, Weight_, Coord_, Point_>
    {
    public:
      typedef Shape::Simplex<1> ShapeType;

      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
        // functor.template factory<YourFactoryName>();

        // <<< END OF CUBATURE FACTORY LIST <<<

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        driver(driver_functor);
      }
    }; // class FactoryExplicitWrapper<Simplex<1>,...>

    /**
     * \brief Explicit specialisation for Simplex<2> shape
     *
     * \author Peter Zajac
     */
    template<
      typename Weight_,
      typename Coord_,
      typename Point_>
    class FactoryExplicitWrapper<Shape::Simplex<2>, Weight_, Coord_, Point_>
    {
    public:
      typedef Shape::Simplex<2> ShapeType;

      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
        // functor.template factory<YourFactoryName>();

        // <<< END OF CUBATURE FACTORY LIST <<<

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        driver(driver_functor);
      }
    }; // class FactoryExplicitWrapper<Simplex<2>,...>

    /**
     * \brief Explicit specialisation for Simplex<3> shape
     *
     * \author Peter Zajac
     */
    template<
      typename Weight_,
      typename Coord_,
      typename Point_>
    class FactoryExplicitWrapper<Shape::Simplex<3>, Weight_, Coord_, Point_>
    {
    public:
      typedef Shape::Simplex<3> ShapeType;

      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
        // functor.template factory<YourFactoryName>();

        // <<< END OF CUBATURE FACTORY LIST <<<

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        driver(driver_functor);
      }
    }; // class FactoryExplicitWrapper<Simplex<3>,...>

    /**
     * \brief Explicit specialisation for Hypercube<1> shape
     *
     * \author Peter Zajac
     */
    template<
      typename Weight_,
      typename Coord_,
      typename Point_>
    class FactoryExplicitWrapper<Shape::Hypercube<1>, Weight_, Coord_, Point_>
    {
    public:
      typedef Shape::Hypercube<1> ShapeType;

      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
        // functor.template factory<YourFactoryName>();

        // <<< END OF CUBATURE FACTORY LIST <<<

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        driver(driver_functor);
      }
    }; // class FactoryExplicitWrapper<Hypercube<1>,...>

    /**
     * \brief Explicit specialisation for Hypercube<2> shape
     *
     * \author Peter Zajac
     */
    template<
      typename Weight_,
      typename Coord_,
      typename Point_>
    class FactoryExplicitWrapper<Shape::Hypercube<2>, Weight_, Coord_, Point_>
    {
    public:
      typedef Shape::Hypercube<2> ShapeType;

      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
        // functor.template factory<YourFactoryName>();

        // <<< END OF CUBATURE FACTORY LIST <<<

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        driver(driver_functor);
      }
    }; // class FactoryExplicitWrapper<Hypercube<2>,...>

    /**
     * \brief Explicit specialisation for Hypercube<3> shape
     *
     * \author Peter Zajac
     */
    template<
      typename Weight_,
      typename Coord_,
      typename Point_>
    class FactoryExplicitWrapper<Shape::Hypercube<3>, Weight_, Coord_, Point_>
    {
    public:
      typedef Shape::Hypercube<3> ShapeType;

      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
        // functor.template factory<YourFactoryName>();

        // <<< END OF CUBATURE FACTORY LIST <<<

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        driver(driver_functor);
      }
    }; // class FactoryExplicitWrapper<Hypercube<3>,...>

    /**
     * \brief Partially Specialised Cubature Factory Wrapper class template
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class FactoryPartialWrapper;

    /**
     * \brief Partial specialisation for Simplex shapes
     *
     * \author Peter Zajac
     */
    template<
      int dim_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class FactoryPartialWrapper<Shape::Simplex<dim_>, Weight_, Coord_, Point_> :
      public FactoryExplicitWrapper<Shape::Simplex<dim_>, Weight_, Coord_, Point_>
    {
    public:
      typedef Shape::Simplex<dim_> ShapeType;

      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // first: call the base class driver function template
        FactoryExplicitWrapper<ShapeType, Weight_, Coord_, Point_>::driver(functor);

        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // first: call the base class factory function template
        FactoryExplicitWrapper<ShapeType, Weight_, Coord_, Point_>::factory(functor);

        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
        // functor.template factory<YourFactoryName>();

        // <<< END OF CUBATURE FACTORY LIST <<<

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        driver(driver_functor);
      }
    }; // class FactoryPartialWrapper<Simplex<...>,...>

    /**
     * \brief Partial specialisation for Hypercube shapes
     *
     * \author Peter Zajac
     */
    template<
      int dim_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class FactoryPartialWrapper<Shape::Hypercube<dim_>, Weight_, Coord_, Point_> :
      public FactoryExplicitWrapper<Shape::Hypercube<dim_>, Weight_, Coord_, Point_>
    {
    public:
      typedef Shape::Hypercube<dim_> ShapeType;

      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // first: call the base class driver function template
        FactoryExplicitWrapper<ShapeType, Weight_, Coord_, Point_>::driver(functor);

        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // first: call the base class factory function template
        FactoryExplicitWrapper<ShapeType, Weight_, Coord_, Point_>::factory(functor);

        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
        // functor.template factory<YourFactoryName>();

        // <<< END OF CUBATURE FACTORY LIST <<<

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        driver(driver_functor);
      }
    }; // class FactoryPartialWrapper<Hypercube<...>,...>

    /**
     * \brief Generic Cubature Factory Wrapper class template
     *
     * \author Peter Zajac
     */
    template<
      typename Shape_,
      typename Weight_,
      typename Coord_,
      typename Point_>
    class FactoryWrapper :
      public FactoryPartialWrapper<Shape_, Weight_, Coord_, Point_>
    {
    public:
      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // first: call the base class driver function template
        FactoryPartialWrapper<Shape_, Weight_, Coord_, Point_>::driver(functor);

        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();
        functor.template driver<BarycentreDriver>();
        functor.template driver<TrapezoidalDriver>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // first: call the base class factory function template
        FactoryPartialWrapper<Shape_, Weight_, Coord_, Point_>::factory(functor);

        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
        // functor.template factory<YourFactoryName>();

        // <<< END OF CUBATURE FACTORY LIST <<<

        // call driver factory functor
        Intern::DriverFactoryFunctor<Shape_, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        driver(driver_functor);

        // call tensor-product functor
        Intern::TensorProductFunctor<Shape_, Weight_, Coord_, Point_, Functor_> tensor_functor(functor);
        Scalar::FactoryWrapper<Weight_, Coord_>::driver(tensor_functor);
      }
    }; // class FactoryWrapper<...>
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_FACTORY_WRAPPER_HPP
