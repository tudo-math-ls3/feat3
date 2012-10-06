#pragma once
#ifndef KERNEL_CUBATURE_FACTORY_WRAPPER_HPP
#define KERNEL_CUBATURE_FACTORY_WRAPPER_HPP 1

// includes, FEAST
#include <kernel/cubature/driver_factory.hpp>
#include <kernel/cubature/refine_factory.hpp>
#include <kernel/cubature/barycentre_driver.hpp>
#include <kernel/cubature/trapezoidal_driver.hpp>
#include <kernel/cubature/hammer_stroud_d2_driver.hpp>
#include <kernel/cubature/hammer_stroud_d3_driver.hpp>
#include <kernel/cubature/hammer_stroud_d5_driver.hpp>
#include <kernel/cubature/lauffer_d2_driver.hpp>
#include <kernel/cubature/lauffer_d4_driver.hpp>
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
      class TensorProductFunctorHelper;
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
    protected:
      typedef Shape::Simplex<1> ShapeType;

      template<typename Functor_>
      static void _driver_list(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void _factory_list(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
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

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        _driver_list(driver_functor);
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
    protected:
      typedef Shape::Simplex<2> ShapeType;

      template<typename Functor_>
      static void _driver_list(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();
        functor.template driver<HammerStroudD2Driver>();
        functor.template driver<HammerStroudD3Driver>();
        functor.template driver<LaufferD2Driver>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void _factory_list(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
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

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        _driver_list(driver_functor);
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
    protected:
      typedef Shape::Simplex<3> ShapeType;

      template<typename Functor_>
      static void _driver_list(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();
        functor.template driver<HammerStroudD2Driver>();
        functor.template driver<HammerStroudD3Driver>();
        functor.template driver<HammerStroudD5Driver>();
        functor.template driver<LaufferD2Driver>();
        functor.template driver<LaufferD4Driver>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void _factory_list(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
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

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        _driver_list(driver_functor);
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
    protected:
      typedef Shape::Hypercube<1> ShapeType;

      template<typename Functor_>
      static void _driver_list(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void _factory_list(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
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

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        _driver_list(driver_functor);
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
    protected:
      typedef Shape::Hypercube<2> ShapeType;

      template<typename Functor_>
      static void _driver_list(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void _factory_list(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
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

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        _driver_list(driver_functor);
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
    protected:
      typedef Shape::Hypercube<3> ShapeType;

      template<typename Functor_>
      static void _driver_list(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void _factory_list(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
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

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        _driver_list(driver_functor);
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
    protected:
      typedef Shape::Simplex<dim_> ShapeType;

      template<typename Functor_>
      static void _driver_list(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();


        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void _factory_list(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
        // functor.template factory<YourFactoryName>();

        // <<< END OF CUBATURE FACTORY LIST <<<
      }

    public:
      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // first: call the base class driver function template
        FactoryExplicitWrapper<ShapeType, Weight_, Coord_, Point_>::driver(functor);

        // call driver list
        _driver_list(functor);
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // first: call the base class factory function template
        FactoryExplicitWrapper<ShapeType, Weight_, Coord_, Point_>::factory(functor);

        // call factory list
        _factory_list(functor);

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        _driver_list(driver_functor);
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
    protected:
      typedef Shape::Hypercube<dim_> ShapeType;

      template<typename Functor_>
      static void _driver_list(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void _factory_list(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
        // functor.template factory<YourFactoryName>();

        // <<< END OF CUBATURE FACTORY LIST <<<
      }

    public:
      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // first: call the base class driver function template
        FactoryExplicitWrapper<ShapeType, Weight_, Coord_, Point_>::driver(functor);

        // call driver list
        _driver_list(functor);
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // first: call the base class factory function template
        FactoryExplicitWrapper<ShapeType, Weight_, Coord_, Point_>::factory(functor);

        // call factory list
        _factory_list(functor);

        // call tensor-product functor
        TensorProductFunctor<Functor_> tensor_functor(functor);
        Scalar::FactoryWrapper<Weight_, Coord_>::driver(tensor_functor);

        // last: call driver factory functor
        Intern::DriverFactoryFunctor<ShapeType, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        _driver_list(driver_functor);
      }
      /// \cond internal
    private:
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
          Intern::TensorProductFunctorHelper<ScalarDriver_, Shape::Hypercube<dim_>, Weight_, Coord_, Point_>::scalar_driver(_functor);
        }
      };
      /// \endcond
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
    protected:
      template<typename Functor_>
      static void _driver_list(Functor_& functor)
      {
        // >>> CUBATURE DRIVER LIST >>>
        // TODO: add you new cubature driver at the end of the list below, e.g.
        // functor.template driver<YourDriverName>();
        functor.template driver<BarycentreDriver>();
        functor.template driver<TrapezoidalDriver>();

        // <<< END OF CUBATURE DRIVER LIST <<<
      }

      template<typename Functor_>
      static void _factory_list(Functor_& functor)
      {
        // >>> CUBATURE FACTORY LIST >>>
        // TODO: add you new cubature factory at the end of the list below, e.g.
        // functor.template factory<YourFactoryName>();

        // <<< END OF CUBATURE FACTORY LIST <<<
      }

    public:
      template<typename Functor_>
      static void factory_no_refine(Functor_& functor)
      {
        // first: call the base class factory function template
        FactoryPartialWrapper<Shape_, Weight_, Coord_, Point_>::factory(functor);

        // call factory list
        _factory_list(functor);

        // call driver factory functor
        Intern::DriverFactoryFunctor<Shape_, Weight_, Coord_, Point_, Functor_> driver_functor(functor);
        _driver_list(driver_functor);
      }

      template<typename Functor_>
      static void driver(Functor_& functor)
      {
        // first: call the base class driver function template
        FactoryPartialWrapper<Shape_, Weight_, Coord_, Point_>::driver(functor);

        // call driver list
        _driver_list(functor);
      }

      template<typename Functor_>
      static void factory(Functor_& functor)
      {
        // call non-refine factory list
        factory_no_refine(functor);

        // call refinement factory functor
        RefineFactoryFunctor<Functor_> refine_functor(functor);
        factory_no_refine(refine_functor);
      }

      /// \cond internal
    private:
      template<typename Functor_>
      class RefineFactoryFunctor
      {
      protected:
        Functor_& _functor;

      public:
        explicit RefineFactoryFunctor(Functor_& functor) :
          _functor(functor)
        {
        }

        template<typename Factory_>
        void factory()
        {
          _functor.template factory< RefineFactory<Factory_> >();
        }
      };
      /// \endcond
    }; // class FactoryWrapper<...>

    /// \cond internal
    namespace Intern
    {
      template<
        template<typename,typename> class ScalarDriver_,
        typename Shape_,
        typename Weight_,
        typename Coord_,
        typename Point_,
        bool tensorise_>
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
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_FACTORY_WRAPPER_HPP