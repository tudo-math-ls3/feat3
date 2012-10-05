#pragma once
#ifndef KERNEL_CUBATURE_AUTO_ALIAS_HPP
#define KERNEL_CUBATURE_AUTO_ALIAS_HPP 1

// includes, FEAST
#include <kernel/shape.hpp>
#include <kernel/cubature/scalar/gauss_legendre_driver.hpp>

// includes, STL
#include <algorithm>

namespace FEAST
{
  namespace Cubature
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      class AutoDegree;
    };
    /// \endcond

    template<typename Shape_>
    class AutoAlias
    {
    public:
      static String map(const String& name)
      {
        // all auto-aliases have the format "auto-<mode>:<param>"

        // try to find a colon within the string
        String::size_type l = name.find_first_of('-');
        String::size_type k = name.find_first_of(':');
        if((l == name.npos) || (k == name.npos))
          return name;

        // the '-' must preceed the ':'
        if(l >= k)
          return name;

        // separate the strings
        String head(name.substr(0, l));
        String mode(name.substr(l+1, k-l-1));
        String param(name.substr(k+1));

        // check head
        if(head.trim().compare_no_case("auto") != 0)
          return name;

        // check mode
        if(mode.trim().compare_no_case("degree") == 0)
        {
          // auto-degree
          Index degree = 0;
          if(!param.parse(degree))
            return name;

          // choose auto-degree alias
          return Intern::AutoDegree<Shape_>::choose(degree);
        }

        // unknown auto alias
        return name;
      }
    };

    // \cond internal
    namespace Intern
    {
      template<int dim_>
      class AutoDegree< Shape::Simplex<dim_> >
      {
      public:
        static String choose(Index degree)
        {
          if(dim_ == 1) // do the same as for hypercubes
          {
            // k-point Gauss-Legendre cubature is exact up to a degree of 2*k-1,
            // so for a degree of n we need k := (n+2)/2 = n/2 + 1
            Index k = degree/2 + 1;

            // Also, ensure that we respect the currently implemented borders...
            k = std::max(k, Index(Scalar::GaussLegendreDriverBase::min_points));
            k = std::min(k, Index(Scalar::GaussLegendreDriverBase::max_points));

            // and build the name
            return
#ifdef FEAST_CUBATURE_TENSOR_PREFIX
            "tensor:" +
#endif
            Scalar::GaussLegendreDriverBase::name() + ":" + stringify(k);
          }
          else if(degree <= 1) // barycentre
          {
            return BarycentreDriverBase::name();
          }
          else if(degree == 2) // Hammer-Stroud of degree 2
          {
            return HammerStroudD2DriverBase::name();
          }
          else if(degree == 3) // Hammer-Stroud of degree 3
          {
            return HammerStroudD3DriverBase::name();
          }
          else if(dim_ == 3)
          {
            if(degree <= 4) // Lauffer formula of degree 4 only for tetrahedron
            {
              return LaufferD4DriverBase::name();
            }
            else // Hammer-Stroud of degree 5 only for tetrahedra
            {
              return HammerStroudD5DriverBase::name();
            }
          }
          else
            return HammerStroudD3DriverBase::name();
        }
      };

      template<int dim_>
      class AutoDegree< Shape::Hypercube<dim_> >
      {
      public:
        static String choose(Index degree)
        {
          // k-point Gauss-Legendre cubature is exact up to a degree of 2*k-1,
          // so for a degree of n we need k := (n+2)/2 = n/2 + 1
          Index k = degree/2 + 1;

          // Also, ensure that we respect the currently implemented borders...
          k = std::max(k, Index(Scalar::GaussLegendreDriverBase::min_points));
          k = std::min(k, Index(Scalar::GaussLegendreDriverBase::max_points));

          // and build the name
          return
#ifdef FEAST_CUBATURE_TENSOR_PREFIX
            "tensor:" +
#endif
            Scalar::GaussLegendreDriverBase::name() + ":" + stringify(k);
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Cubature
} // namespace FEAST

#endif // KERNEL_CUBATURE_AUTO_ALIAS_HPP
