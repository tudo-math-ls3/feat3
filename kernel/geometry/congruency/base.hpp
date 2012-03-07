#pragma once
#ifndef KERNEL_GEOMETRY_CONGRUENCY_BASE_HPP
#define KERNEL_GEOMETRY_CONGRUENCY_BASE_HPP 1

// includes, FEAST
#include <kernel/geometry/shape.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Congruency namespace
     */
    namespace Congruency
    {
      /**
       * \brief Congruency sampler class template.
       *
       * \todo detailed description
       *
       * \author Peter Zajac
       */
      template<typename Shape_>
#ifndef DOXYGEN
      class Sampler;
#else
      class Sampler
      {
      public:
        /**
         * \brief Compares two vertex index vectors and returns the orientation.
         *
         * \param[in] src
         * A reference to the source index vector that is to be compared.
         *
         * \param[in] trg
         * A reference to the target index vector that is to be compared against.
         *
         * \returns
         * An orientation code describing the relation between the source and target index vectors.
         */
        template<
          typename Source_,
          typename Target_>
        static int compare(const Source_& src, const Target_& trg);
      };
#endif // DOXYGEN

      /**
       * \brief Congruency index mapping class template.
       *
       * \todo detailed description
       *
       * \tparam Shape_
       * The shape tag class for the index mapping.
       *
       * \tparam face_dim_
       * The dimension of the faces for the index mapping.
       *
       * \author Peter Zajac
       */
      template<
        typename Shape_,
        int face_dim_>
#ifndef DOXYGEN
      class IndexMapping;
#else
      class IndexMapping
      {
      public:
        /**
         * \brief Maps a local index based on an orientation code.
         *
         * \param[in] orient
         * An orientation code as returned by the Sampler<Shape_>::compare() function.
         *
         * \param[in] idx
         * The local face index that is to be mapped.
         *
         * \returns
         * The mapped local face index.
         */
        static int map(int orient, int idx);
      };
#endif // DOXYGEN
    } // namespace Congruency
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_CONGRUENCY_BASE_HPP
