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
      /**
       * \brief Bernstein-2 Basisfunctions
       *
       * \author Gesa Pottbrock
       */
      namespace Intern
      {
        // p0, p1 and p2 are the 1D Bernstein basis functions on the reference interval [-1,+1].
        // These are used for the tensor-product approach in the Hypercube evaluators.

        // Bernstein basis function for left vertex.
        template<typename T_>
        inline T_ p0(T_ x)
        {
          return T_(0.25) * (T_(1) - T_(2) * x + x * x);
        }

        // Bernstein basis function for right vertex.
        template<typename T_>
        inline T_ p1(T_ x)
        {
          return T_(0.25) * (T_(1) + T_(2) * x + x * x);
        }

        // Bernstein basis function for edge midpoint.
        template<typename T_>
        inline T_ p2(T_ x)
        {

          return T_(0.5) * (T_(1) - x * x);
        }

        // The first order derivatives of the basis functions.

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

        // The second order derivatives of the basis functions.

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

        // Templatefunction to evaluate the basisfunctions without having to call the function manually.
        // This is especially great, if you need to evaluate the basisfunctions inside a for-loop.
        template<typename T_>
        T_ p(int i,T_ x)
        {
          switch (i)
          {
          case 0:
            return p0(x);
            break;
          case 1:
            return p1(x);
            break;
          case 2:
            return p2(x);
          default:
            XABORTM("You can only use p0,p1,p2. The basisfunction p"+stringify( i)+" is not available.");
            std::abort();
          }
        }

        // q0, q1 and q2 are the 1D basis functions on the reference interval [-1,+1].
        // These are used for the tensor-product approach in the Hypercube evaluators.
        // The basisfunctions satisfy  int_(-1,1) p_i(x)*q_j(x) dx = kron_i,j

        // Orthogonal basis function for left vertex.
        template<typename T_>
        inline T_ q0(T_ x)
        {
          return T_(-0.75) - T_(1.5) * x + T_(3.75) * x * x;
        }

        // Orthogonal basis function for right vertex.
        template<typename T_>
        inline T_ q1(T_ x)
        {
          return T_(-0.75) + T_(1.5) * x + T_(3.75) * x * x;
        }

        // Orthogonal basis function for edge midpoint.
        template<typename T_>
        inline T_ q2(T_ x)
        {
          return T_(3) - T_(7.5) * x * x;
        }

        // The first order derivatives of the orthogonal basis functions.

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

        // The second order derivatives of the orthogonal basis functions.

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

        // Templatefunction to evaluate the orthogonal functions without having to call the function manually.
        // This is especially great, if you need to evaluate the orthogonal functions inside a for-loop.
        template<typename T_>
        T_ q(int i, T_ x)
        {
          switch (i)
          {
          case 0:
            return q0(x);
            break;
          case 1:
            return q1(x);
            break;
          case 2:
            return q2(x);
          default:
            XABORTM("You can only use q0,q1,q2. The basisfunction q"+ stringify(i)+" is not available.");
            std::abort();
          }
        }

        // This function evaluates sum_{i=0}^2 koef_i *q_i(x) and is used in node_functional.
        template<typename T_>
        T_ base_q(T_ x, Tiny::Vector<T_, 3> koef)
        {
          T_ rsl = T_(0);
          for (int qx = 0; qx < 3; ++qx)
          {
              rsl += koef(qx ) * q(qx, x);
          } // end for qx
          return rsl;
        }

        // This function evaluates sum_{i,j=0}^2 koef_i,j *q_i(x)*q_j(y) and is used in node_functional.
        template<typename T_>
        T_ base_q(T_ x, T_ y, Tiny::Vector<T_, 9> koef)
        {
          T_ rsl = T_(0);
          for (int qx = 0; qx < 3; ++qx)
          {
            for (int qy = 0; qy < 3; ++qy)
            {
              rsl+=koef( qx * 3 + qy)* q(qx,x)*q(qy,y);
            } // end for qy
          } // end for qx
          return rsl;
        }

        // This function evaluates sum_{i,j,k=0}^2 koef_i,j,k *q_i(x)*q_j(y)*q_k(z) and is used in node_functional.
        template<typename T_>
        T_ base_q(T_ x, T_ y, T_ z, Tiny::Vector<T_, 27> koef)
        {
          T_ rsl = T_(0);
          for (int qz = 0; qz < 3; ++qz)
          {
            for (int qx = 0; qx < 3; ++qx)
            {
              for (int qy = 0; qy < 3; ++qy)
              {
                rsl += koef(qz*9+qx * 3 + qy) * q(qx, x) * q(qy, y)*q(qz,z);
              } // end for qy
            } // end for qx
          } // end for qz
          return rsl;
        }
      } // namespace Intern
    } // namespace Bernstein2
  } // namespace Space
} // namespace FEAT

#endif // KERNEL_SPACE_BERNSTEIN2_SCALAR_BASIS_HPP
