#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_CONGRUENCY_TRAFO_HPP
#define KERNEL_GEOMETRY_INTERN_CONGRUENCY_TRAFO_HPP 1

#include <kernel/shape.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Congruency trafo class template.
       *
       * \todo detailed description
       *
       * \author Peter Zajac
       */
      template<typename Shape_>
#ifndef DOXYGEN
      class CongruencyTrafo;
#else
      class CongruencyTrafo
      {
      public:
        /**
         * \brief Computes the transformation matrix for an orientation code.
         *
         * \param[out] a
         * The transformation matrix for the orientation code.
         *
         * \param[out] b
         * The transformation vector for the orientation code.
         *
         * \param[in] orient
         * An orientation code as returned by the CongruencySampler<Shape_>::compare() function.
         */
        template<typename DataType_, int sm_, int sn_, int sl_>
        static void compute(Tiny::Matrix<DataType_, m, n, sm_, sn_>& a, Tiny::Vector<DataType_, n, sl_>& b, int orient);
      };
#endif // DOXYGEN

      /**
       * \brief Congruency trafo implementation for Simplex<1> shape
       */
      template<>
      class CongruencyTrafo<Shape::Simplex<1>>
      {
      public:
        template<typename DataType_, int sm_, int sn_, int sl_>
        static void compute(Tiny::Matrix<DataType_, 1, 1, sm_, sn_>& a, Tiny::Vector<DataType_, 1, sl_>& b, int orient)
        {
          switch(orient)
          {
          case 0:
            a(Index(0),Index(0)) = +DataType_(1);
            b(Index(0)) =  DataType_(0);
            break;
          case 1:
            a(Index(0),Index(0)) = -DataType_(1);
            b(Index(0)) = +DataType_(1);
            break;
          }
        }
      };

      /**
       * \brief Congruency trafo implementation for Simplex<2> shape
       */
      template<>
      class CongruencyTrafo<Shape::Simplex<2>>
      {
      public:
        template<typename DataType_, int sm_, int sn_, int sl_>
        static void compute(Tiny::Matrix<DataType_, 2, 2, sm_, sn_>& a, Tiny::Vector<DataType_, 2, sl_>& b, int orient)
        {
          switch(orient)
          {
          case 0:
            a(Index(0),Index(0)) = +DataType_(1);
            a(Index(1),Index(0)) =  DataType_(0);
            a(Index(0),Index(1)) =  DataType_(0);
            a(Index(1),Index(1)) = +DataType_(1);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(0);
            break;
          case 1:
            a(Index(0),Index(0)) = -DataType_(1);
            a(Index(1),Index(0)) = +DataType_(1);
            a(Index(0),Index(1)) = -DataType_(1);
            a(Index(1),Index(1)) =  DataType_(0);
            b(Index(0)) = +DataType_(1);
            b(Index(1)) =  DataType_(0);
            break;
         case 2:
            a(Index(0),Index(0)) =  DataType_(0);
            a(Index(1),Index(0)) = -DataType_(1);
            a(Index(0),Index(1)) =  DataType_(1);
            a(Index(1),Index(1)) = -DataType_(1);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) = +DataType_(1);
            break;
         case 4:
            a(Index(0),Index(0)) =  DataType_(0);
            a(Index(1),Index(0)) = +DataType_(1);
            a(Index(0),Index(1)) = +DataType_(1);
            a(Index(1),Index(1)) =  DataType_(0);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) =  DataType_(0);
            break;
         case 5:
            a(Index(0),Index(0)) = -DataType_(1);
            a(Index(1),Index(0)) =  DataType_(0);
            a(Index(0),Index(1)) = -DataType_(1);
            a(Index(1),Index(1)) = +DataType_(1);
            b(Index(0)) = +DataType_(1);
            b(Index(1)) =  DataType_(0);
            break;
         case 6:
            a(Index(0),Index(0)) = +DataType_(1);
            a(Index(1),Index(0)) = -DataType_(1);
            a(Index(0),Index(1)) =  DataType_(0);
            a(Index(1),Index(1)) = -DataType_(1);
            b(Index(0)) =  DataType_(0);
            b(Index(1)) = +DataType_(1);
            break;
          }
        }
      };

      /**
       * \brief Congruency trafo implementation for Hypercube<1> shape
       */
      template<>
      class CongruencyTrafo<Shape::Hypercube<1>>
      {
      public:
        template<typename DataType_, int sm_, int sn_, int sl_>
        static void compute(Tiny::Matrix<DataType_, 1, 1, sm_, sn_>& a, Tiny::Vector<DataType_, 1, sl_>& b, int orient)
        {
          switch(orient)
          {
          case 0:
            a(Index(0),Index(0)) = +DataType_(1);
            break;
          case 1:
            a(Index(0),Index(0)) = -DataType_(1);
            break;
          }
          b(Index(0)) = DataType_(0);
        }
      };

      /**
       * \brief Congruency trafo implementation for Hypercube<2> shape
       */
      template<>
      class CongruencyTrafo<Shape::Hypercube<2>>
      {
      public:
        template<typename DataType_, int sm_, int sn_, int sl_>
        static void compute(Tiny::Matrix<DataType_, 2, 2, sm_, sn_>& a, Tiny::Vector<DataType_, 2, sl_>& b, int orient)
        {
          switch(orient)
          {
          case 0:
            a(Index(0),Index(0)) = +DataType_(1);
            a(Index(1),Index(0)) =  DataType_(0);
            a(Index(0),Index(1)) =  DataType_(0);
            a(Index(1),Index(1)) = +DataType_(1);
            break;
          case 1:
            a(Index(0),Index(0)) =  DataType_(0);
            a(Index(1),Index(0)) = +DataType_(1);
            a(Index(0),Index(1)) = -DataType_(1);
            a(Index(1),Index(1)) =  DataType_(0);
            break;
          case 2:
            a(Index(0),Index(0)) =  DataType_(0);
            a(Index(1),Index(0)) = -DataType_(1);
            a(Index(0),Index(1)) = +DataType_(1);
            a(Index(1),Index(1)) =  DataType_(0);
            break;
          case 3:
            a(Index(0),Index(0)) = -DataType_(1);
            a(Index(1),Index(0)) =  DataType_(0);
            a(Index(0),Index(1)) =  DataType_(0);
            a(Index(1),Index(1)) = -DataType_(1);
            break;
          case 4:
            a(Index(0),Index(0)) =  DataType_(0);
            a(Index(1),Index(0)) = +DataType_(1);
            a(Index(0),Index(1)) = +DataType_(1);
            a(Index(1),Index(1)) =  DataType_(0);
            break;
          case 5:
            a(Index(0),Index(0)) = -DataType_(1);
            a(Index(1),Index(0)) =  DataType_(0);
            a(Index(0),Index(1)) =  DataType_(0);
            a(Index(1),Index(1)) = +DataType_(1);
            break;
          case 6:
            a(Index(0),Index(0)) = +DataType_(1);
            a(Index(1),Index(0)) =  DataType_(0);
            a(Index(0),Index(1)) =  DataType_(0);
            a(Index(1),Index(1)) = -DataType_(1);
            break;
          case 7:
            a(Index(0),Index(0)) =  DataType_(0);
            a(Index(1),Index(0)) = -DataType_(1);
            a(Index(0),Index(1)) = -DataType_(1);
            a(Index(1),Index(1)) =  DataType_(0);
            break;
          }
          b(Index(0)) = DataType_(0);
          b(Index(0)) = DataType_(0);
        }
      };

    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_INTERN_CONGRUENCY_TRAFO_HPP
