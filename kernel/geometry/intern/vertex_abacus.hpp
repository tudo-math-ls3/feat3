#pragma once
#ifndef KERNEL_GEOMETRY_INTERN_VERTEX_ABACUS_HPP
#define KERNEL_GEOMETRY_INTERN_VERTEX_ABACUS_HPP 1

// includes, FEAST
#include <kernel/geometry/vertex_set.hpp>

namespace FEAST
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
    {
      /**
       * \brief Vertex Abacus class template.
       *
       * This class template is a helper class that provides functions for commonly used operations on vertices
       * of a vertex set.
       *
       * \author Peter Zajac
       */
      template<typename VertexSet_>
      class VertexAbacus
      {
      public:
        /// vertex set type
        typedef VertexSet_ VertexSetType;

        /// vertex coordinate type
        typedef typename VertexSetType::CoordType CoordType;

        /// vertex reference type
        typedef typename VertexSetType::VertexReference VertexReference;

        /// const vertex reference type
        typedef typename VertexSetType::ConstVertexReference ConstVertexReference;

      protected:
        /// number of coordinates per vertex.
        int _num_coords;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] vertex_set
         * A reference to the vertex set whose vertices are to be operated on.
         */
        explicit VertexAbacus(const VertexSetType& vertex_set) :
          _num_coords(vertex_set.get_num_coords())
        {
        }

        /**
         * \brief Clears a vertex.
         *
         * This function sets all coordinates of a vertex to one value.
         *
         * \param[in] v
         * A reference to the vertex to be cleared.
         *
         * \param[in] alpha
         * The value that the coordinates of \p v are to be set to.
         */
        void clear(VertexReference v, CoordType alpha = CoordType(0)) const
        {
          for(int i(0); i < _num_coords; ++i)
          {
            v[i] = alpha;
          }
        }

        /**
         * \brief Copies a vertex.
         *
         * \param[out] v
         * A reference to the destination vertex.
         *
         * \param[in] w
         * A const reference to the source vertex.
         */
        void copy(VertexReference v, ConstVertexReference w) const
        {
          for(int i(0); i < _num_coords; ++i)
          {
            v[i] = w[i];
          }
        }

        /**
         * \brief Adds one vertex onto another.
         *
         * \param[in,out] v
         * A reference to the destination vertex.
         *
         * \param[in] w
         * A const reference to the summand vertex.
         */
        void add(VertexReference v, ConstVertexReference w) const
        {
          for(int i(0); i < _num_coords; ++i)
          {
            v[i] += w[i];
          }
        }

        /**
         * \brief Scales a vertex.
         *
         * This function multiples all coordinates of a vertex by a scalar value.
         *
         * \param[in,out] v
         * A reference to the vertex to be scaled.
         *
         * \param[in] alpha
         * The value that \p v is to be scaled by.
         */
        void scale(VertexReference v, CoordType alpha) const
        {
          for(int i(0); i < _num_coords; ++i)
          {
            v[i] *= alpha;
          }
        }

        /**
         * \brief Calculates the mean of two vertices.
         *
         * \param[out] v
         * A reference to the destination vertex.
         *
         * \param[in] w0, w1
         * Const references to the vertices whose mean is to be calculated.
         */
        void mean(VertexReference v, ConstVertexReference w0, ConstVertexReference w1) const
        {
          for(int i(0); i < _num_coords; ++i)
          {
            v[i] = CoordType(0.5) * (w0[i] + w1[i]);
          }
        }

        /**
         * \brief Calculates the mean of four vertices.
         *
         * \param[out] v
         * A reference to the destination vertex.
         *
         * \param[in] w0, w1, w2, w3
         * Const references to the vertices whose mean is to be calculated.
         */
        void mean(VertexReference v,
          ConstVertexReference w0, ConstVertexReference w1, ConstVertexReference w2, ConstVertexReference w3) const
        {
          for(int i(0); i < _num_coords; ++i)
          {
            v[i] = CoordType(0.25) * (w0[i] + w1[i] + w2[i] + w3[i]);
          }
        }

        /**
         * \brief Calculates the mean of eight vertices.
         *
         * \param[out] v
         * A reference to the destination vertex.
         *
         * \param[in] w0, w1, w2, w3, w4, w5, w6, w7
         * Const references to the vertices whose mean is to be calculated.
         */
        void mean(VertexReference v,
          ConstVertexReference w0, ConstVertexReference w1, ConstVertexReference w2, ConstVertexReference w3,
          ConstVertexReference w4, ConstVertexReference w5, ConstVertexReference w6, ConstVertexReference w7) const
        {
          for(int i(0); i < _num_coords; ++i)
          {
            v[i] = CoordType(0.125) * (w0[i] + w1[i] + w2[i] + w3[i] + w4[i] + w5[i] + w6[i] + w7[i]);
          }
        }
      }; // class VertexAbacus<...>
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_INTERN_VERTEX_ABACUS_HPP
