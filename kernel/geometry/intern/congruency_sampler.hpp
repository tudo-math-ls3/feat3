// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/shape.hpp>

namespace FEAT
{
  namespace Geometry
  {
    /// \cond internal
    namespace Intern
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
      class CongruencySampler;
#else
      class CongruencySampler
      {
      public:
        /**
         * \brief Compares two vertex index vectors and returns the orientation code.
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

        /**
         * \brief Maps an orientation code to positive or negative
         *
         * From the orientation code of two representatives of the same vertex index vectors we can tell if they
         * describe the same entity in the same or different orientation.
         *
         * \returns 1 if the orientation code belongs to two orientations of the same sign, -1 if they are of
         * different sign or 0 if the orientation code was -1.
         */
        static int orientation(int orientation_code);
      };
#endif // DOXYGEN

      /**
       * \brief Congruency sampler implementation for Simplex<1> shape
       *
         \verbatim
               0        1
             0---1    1---0
         \endverbatim
       *
       * \author Peter Zajac
       */
      template<>
      class CongruencySampler< Shape::Simplex<1> >
      {
      public:
        template<
          typename Source_,
          typename Target_>
        static int compare(const Source_& src, const Target_& trg)
        {
          return compare_with(src, trg, [](const auto& a, const auto& b) { return a == b; });
        }

        template<
          typename Source_,
          typename Target_,
          typename Comparator_>
        static int compare_with(const Source_& src, const Target_& trg, const Comparator_& comp)
        {
          auto sv0 = src[0];
          if(comp(sv0, trg[0]))
          {
            return 0;
          }
          else if(comp(sv0, trg[1]))
          {
            return 1;
          }
          // invalid
          return -1;
        }

        static int orientation(int orientation_code)
        {
          if(orientation_code == 0)
            return 1;
          if(orientation_code == 1)
            return -1;
          return 0;
        }
      }; // class CongruencySampler<Simplex<1>>

      /**
       * \brief Congruency sampler implementation for Simplex<2> shape
       *
         \verbatim
           0        1       2
             2        1       0
            / \      / \     / \
           0---1    2---0   1---2

           4        5       6
             1        2       0
            / \      / \     / \
           0---2    1---0   2---1
         \endverbatim
       *
       * \note The lack of the number \c 3 is intentional.
       *
       * \author Peter Zajac
       */
      template<>
      class CongruencySampler< Shape::Simplex<2> >
      {
      public:
        template<
          typename Source_,
          typename Target_>
        static int compare(const Source_& src, const Target_& trg)
        {
          return compare_with(src, trg, [](const auto& a, const auto& b) { return a == b; });
        }

        template<
          typename Source_,
          typename Target_,
          typename Comparator_>
        static int compare_with(const Source_& src, const Target_& trg, const Comparator_& comp)
        {
          auto sv0 = src[0];
          auto sv1 = src[1];
          if(comp(sv0, trg[0]))
          {
            if(comp(sv1, trg[1]))
            {
              return 0;
            }
            else if(comp(sv1, trg[2]))
            {
              return 4;
            }
          }
          else if(comp(sv0, trg[1]))
          {
            if(comp(sv1, trg[2]))
            {
              return 1;
            }
            else if(comp(sv1, trg[0]))
            {
              return 5;
            }
          }
          else if(comp(sv0, trg[2]))
          {
            if(comp(sv1, trg[0]))
            {
              return 2;
            }
            else if(comp(sv1, trg[1]))
            {
              return 6;
            }
          }
          // invalid
          return -1;
        }

        static int orientation(int orientation_code)
        {
          if(orientation_code == -1)
            return 0;
          if(orientation_code < 3)
            return 1;
          else
            return -1;
        }
      }; // class CongruencySampler<Simplex<2>>

      /**
       * \brief Congruency sampler implementation for Hypercube<1> shape
       *
         \verbatim
               0        1
             0---1    1---0
         \endverbatim
       *
       * \author Peter Zajac
       */
      template<>
      class CongruencySampler< Shape::Hypercube<1> >
      {
      public:
        template<
          typename Source_,
          typename Target_>
        static int compare(const Source_& src, const Target_& trg)
        {
          return compare_with(src, trg, [](const auto& a, const auto& b) { return a == b; });
        }

        template<
          typename Source_,
          typename Target_,
          typename Comparator_>
        static int compare_with(const Source_& src, const Target_& trg, const Comparator_& comp)
        {
          auto sv0 = src[0];
          if(comp(sv0, trg[0]))
          {
            return 0;
          }
          else if(comp(sv0, trg[1]))
          {
            return 1;
          }
          // invalid
          return -1;
        }

        static int orientation(int orientation_code)
        {
          if(orientation_code == 0)
            return 1;
          if(orientation_code == 1)
            return -1;
          return 0;
        }
      }; // class CongruencySampler<Hypercube<1>>

      /**
       * \brief Congruency sampler implementation for Hypercube<2> shape
       *
         \verbatim
               0        1        2        3
             2---3    3---1    0---2    1---0
             |   |    |   |    |   |    |   |
             0---1    2---0    1---3    3---2

               4        5        6        7
             1---3    3---2    0---1    2---0
             |   |    |   |    |   |    |   |
             0---2    1---0    2---3    3---1
         \endverbatim
       *
       * \author Peter Zajac
       */
      template<>
      class CongruencySampler< Shape::Hypercube<2> >
      {
      public:
        template<
          typename Source_,
          typename Target_>
        static int compare(const Source_& src, const Target_& trg)
        {
          return compare_with(src, trg, [](const auto& a, const auto& b) { return a == b; });
        }

        template<
          typename Source_,
          typename Target_,
          typename Comparator_>
        static int compare_with(const Source_& src, const Target_& trg, const Comparator_& comp)
        {
          auto sv0 = src[0];
          auto sv1 = src[1];
          if(comp(sv0, trg[0]))
          {
            if(comp(sv1, trg[1]))
            {
              return 0;
            }
            else if(comp(sv1, trg[2]))
            {
              return 4;
            }
          }
          else if(comp(sv0, trg[1]))
          {
            if(comp(sv1, trg[3]))
            {
              return 1;
            }
            else if(comp(sv1, trg[0]))
            {
              return 5;
            }
          }
          else if(comp(sv0, trg[2]))
          {
            if(comp(sv1, trg[0]))
            {
              return 2;
            }
            else if(comp(sv1, trg[3]))
            {
              return 6;
            }
          }
          else if(comp(sv0, trg[3]))
          {
            if(comp(sv1, trg[2]))
            {
              return 3;
            }
            else if(comp(sv1, trg[1]))
            {
              return 7;
            }
          }
          // invalid
          return -1;
        }

        static int orientation(int orientation_code)
        {
          if(orientation_code == -1)
            return 0;
          if(orientation_code < 4)
            return 1;
          else
            return -1;
        }
      }; // class CongruencySampler<Hypercube<2>>
    } // namespace Intern
    /// \endcond
  } // namespace Geometry
} // namespace FEAT
