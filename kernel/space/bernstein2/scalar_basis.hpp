// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SPACE_BERNSTEIN2_SCALAR_BASIS_HPP
#define KERNEL_SPACE_BERNSTEIN2_SCALAR_BASIS_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/space/node_functional_base.hpp>
#include <kernel/cubature/dynamic_factory.hpp>


namespace FEAT
{
  namespace Space
  {
    namespace Bernstein2
    {
      /// \cond internal
      namespace Intern
      {
        // p0, p1 and p2 are the 1D basis functions on the reference interval [-1,+1].
        // These are used for the tensor-product approach in the Hypercube evaluators.

        // basis function for left vertex
        template<typename T_>
        inline T_ p0(T_ x)
        {
          return T_(0.25) * (T_(1) - T_(2) * x + x * x);
        }

        // basis function for right vertex
        template<typename T_>
        inline T_ p1(T_ x)
        {
          return T_(0.25) * (T_(1) + T_(2) * x + x * x);
        }

        // basis function for edge midpoint
        template<typename T_>
        inline T_ p2(T_ x)
        {

          return T_(0.5) * (T_(1) - x * x);
        }

        // first order derivatives

        template<typename T_>
        inline T_ d1p0(T_ x)
        {
          return T_(0.5) * (x - T_(1));
        }

        template<typename T_>
        inline T_ d1p1(T_ x)
        {
          return  T_(0.5) * (x + T_(1));
        }

        template<typename T_>
        inline T_ d1p2(T_ x)
        {

          return -x;
        }

        // second order derivatives

        template<typename T_>
        inline T_ d2p0(T_)
        {
          return T_(0.5);
        }

        template<typename T_>
        inline T_ d2p1(T_)
        {
          return T_(0.5);
        }

        template<typename T_>
        inline T_ d2p2(T_)
        {

          return T_(-1);
        }


        // q0, q1 and q2 are the 1D basis functions on the reference interval [-1,+1].
        // These are used for the tensor-product approach in the Hypercube evaluators.

        // basis function for left vertex
        template<typename T_>
        inline T_ q0(T_ x)
        {
          return T_(-0.75) - T_(1.5) * x + T_(3.75) * x * x;
        }

        // basis function for right vertex
        template<typename T_>
        inline T_ q1(T_ x)
        {
          return T_(-0.75) + T_(1.5) * x + T_(3.75) * x * x;
        }

        // basis function for edge midpoint
        template<typename T_>
        inline T_ q2(T_ x)
        {
          return T_(3) - T_(7.5) * x * x;
        }

        // first order derivatives

        template<typename T_>
        inline T_ d1q0(T_ x)
        {
          return T_(7.5) * x - T_(1.5);
        }

        template<typename T_>
        inline T_ d1q1(T_ x)
        {
          return  T_(7.5) * x + T_(1.5);
        }

        template<typename T_>
        inline T_ d1q2(T_ x)
        {

          return T_(-15) * x;
        }

        // second order derivatives

        template<typename T_>
        inline T_ d2q0(T_)
        {
          return T_(7.5);
        }

        template<typename T_>
        inline T_ d2q1(T_)
        {
          return T_(7.5);
        }

        template<typename T_>
        inline T_ d2q2(T_)
        {

          return T_(-15);
        }
      } // namespace Intern
      /// \endcond

    } // namespace Bernstein2
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_BERNSTEIN2_SCALAR_BASIS_HPP
