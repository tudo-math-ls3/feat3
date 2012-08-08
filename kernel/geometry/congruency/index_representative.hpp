#pragma once
#ifndef KERNEL_GEOMETRY_CONGRUENCY_INDEX_REPRESENTATIVE_HPP
#define KERNEL_GEOMETRY_CONGRUENCY_INDEX_REPRESENTATIVE_HPP 1

// includes, FEAST
#include <kernel/geometry/shape.hpp>

// includes, system
#include <algorithm>

namespace FEAST
{
  namespace Geometry
  {
    namespace Congruency
    {

      /**
       * \brief Index representative calculator.
       *
       * \author Constantin Christof
       */
      template<typename Shape_>
#ifndef DOXYGEN
      struct IndexRepresentative;
#else
      struct IndexRepresentative
      {
      public:
        /**
         * \brief Calculates the index vector representative.
         *
         * \param[out] ivo
         * The representative of \p ivi.
         *
         * \param[in] ivi
         * The index vector whose representative is to be computed.
         */
        template<
          typename IndexVectorIn,
          typename IndexVectorOut>
        static void compute(IndexVectorOut& ivo, const IndexVectorIn& ivi);
      };
#endif // DOXYGEN

      /// \cond internal
      // specialisation for Hypercube<1>
      template<>
      struct IndexRepresentative< Shape::Hypercube<1> >
      {
        template<
          typename IndexVectorIn,
          typename IndexVectorOut>
        static void compute(IndexVectorOut& ivo, const IndexVectorIn& ivi)
        {
          ivo[0] = std::min(ivi[0], ivi[1]);
          ivo[1] = std::max(ivi[0], ivi[1]);
        }
      }; // IndexRepresentative< Shape::Hypercube<1> >

      // specialisation for Hypercube<2>
      template<>
      struct IndexRepresentative< Shape::Hypercube<2> >
      {
        template<
          typename IndexVectorIn,
          typename IndexVectorOut>
        static void compute(IndexVectorOut& ivo, const IndexVectorIn& ivi)
        {
          ivo[0] = std::min(std::min(ivi[0], ivi[1]), std::min(ivi[2], ivi[3]));

          int i = 0;

          while(ivi[i] != ivo[0])
          {
            ++i;
          }

          if((i == 0) || (i == 3))
          {
            ivo[1] = std::min(ivi[1], ivi[2]);
            ivo[2] = std::max(ivi[1], ivi[2]);
          }
          else
          {
            ivo[1] = std::min(ivi[0], ivi[3]);
            ivo[2] = std::max(ivi[0], ivi[3]);
          }

          ivo[3] = ivi[0] + ivi[1] + ivi[2] + ivi[3] - ivo[0] - ivo[1] - ivo[2];
        }
      }; // IndexRepresentative< Shape::Hypercube<2> >

      // specialisation for Simplex<1>
      template<>
      struct IndexRepresentative< Shape::Simplex<1> >
      {
        template<
          typename IndexVectorIn,
          typename IndexVectorOut>
        static void compute(IndexVectorOut& ivo, const IndexVectorIn& ivi)
        {
          ivo[0] = std::min(ivi[0], ivi[1]);
          ivo[1] = std::max(ivi[0], ivi[1]);
        }
      }; // IndexRepresentative< Shape::Simplex<1> >

      // specialisation for Simplex<2>
      template<>
      struct IndexRepresentative< Shape::Simplex<2> >
      {
        template<
          typename IndexVectorIn,
          typename IndexVectorOut>
        static void compute(IndexVectorOut& ivo, const IndexVectorIn& ivi)
        {
          ivo[0] = std::min(std::min(ivi[0], ivi[1]), ivi[2]);
          ivo[2] = std::max(std::max(ivi[0], ivi[1]), ivi[2]);
          ivo[1] = ivi[0] + ivi[1] + ivi[2] - ivo[0] - ivo[2];
        }
      }; // IndexRepresentative< Shape::Simplex<2> >
      /// \endcond
    } // namespace Congruency
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONGRUENCY_INDEX_REPRESENTATIVE_HPP
