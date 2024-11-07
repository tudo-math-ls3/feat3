// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// This file defines detail implementations as templated functions.
// This is mainly done due to avoid code duplication for the voxel assembly.
// \author Maximilian Esser, Peter Zajac

#include <kernel/base_header.hpp>
#include <kernel/shape.hpp>

namespace FEAT
{
  namespace Space
  {
    namespace Lagrange2
    {
      /// \cond internal
      namespace Intern
      {
        // p0, p1 and p2 are the 1D basis functions on the reference interval [-1,+1].
        // These are used for the tensor-product approach in the Hypercube evaluators.

        // basis function for left vertex
        template<typename T_>
        CUDA_HOST_DEVICE inline T_ p0(T_ x)
        {
          return T_(0.5) * x * (x - T_(1));
        }

        // basis function for right vertex
        template<typename T_>
        CUDA_HOST_DEVICE inline T_ p1(T_ x)
        {
          return T_(0.5) * x * (x + T_(1));
        }

        // basis function for edge midpoint
        template<typename T_>
        CUDA_HOST_DEVICE inline T_ p2(T_ x)
        {
          return (T_(1) - x) * (T_(1) + x);
        }

        // first order derivatives

        template<typename T_>
        CUDA_HOST_DEVICE inline T_ d1p0(T_ x)
        {
          return x - T_(0.5);
        }

        template<typename T_>
        CUDA_HOST_DEVICE inline T_ d1p1(T_ x)
        {
          return x + T_(0.5);
        }

        template<typename T_>
        CUDA_HOST_DEVICE inline T_ d1p2(T_ x)
        {
          return -T_(2) * x;
        }

        // second order derivatives

        template<typename T_>
        CUDA_HOST_DEVICE inline T_ d2p0(T_)
        {
          return T_(1);
        }

        template<typename T_>
        CUDA_HOST_DEVICE inline T_ d2p1(T_)
        {
          return T_(1);
        }

        template<typename T_>
        CUDA_HOST_DEVICE inline T_ d2p2(T_)
        {
          return -T_(2);
        }
      } // namespace Intern
      /// \endcond

      /// Helper class for evaluation of basis function, gradients and hessians
      /// Specialized for each shapetype
      template<typename PointType_, typename DataType_, typename Shape_>
      struct EvalHelper;

      template<typename PointType_, typename DataType_>
      struct EvalHelper<PointType_, DataType_, Shape::Simplex<2>>
      {
        typedef DataType_ DataType;

        CUDA_HOST_DEVICE static constexpr int get_num_local_dofs()
        {
          return 6;
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_values(EvalData_& data, const PointType_& point)
        {
            // vertex dofs
            data.phi[0].ref_value = DataType(2) * (point[0] + point[1] - DataType(0.5)) * (point[0] + point[1] - DataType(1));
            data.phi[1].ref_value = DataType(2) * point[0] * (point[0] - DataType(0.5));
            data.phi[2].ref_value = DataType(2) * point[1] * (point[1] - DataType(0.5));
            // edge dofs
            data.phi[3].ref_value = DataType(4) * point[0] * point[1];
            data.phi[4].ref_value = DataType(4) * point[1] * (DataType(1) - point[0] - point[1]);
            data.phi[5].ref_value = DataType(4) * point[0] * (DataType(1) - point[0] - point[1]);
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_gradients(
          EvalData_& data,
          const PointType_& point)
        {
          // vertex dofs
          data.phi[0].ref_grad[0] = DataType(4) * (point[0] + point[1]) - DataType(3);
          data.phi[0].ref_grad[1] = DataType(4) * (point[0] + point[1]) - DataType(3);
          data.phi[1].ref_grad[0] = DataType(4) * point[0] - DataType(1);
          data.phi[1].ref_grad[1] = DataType(0);
          data.phi[2].ref_grad[0] = DataType(0);
          data.phi[2].ref_grad[1] = DataType(4) * point[1] - DataType(1);
          // edge dofs
          data.phi[3].ref_grad[0] = DataType(4) * point[1];
          data.phi[3].ref_grad[1] = DataType(4) * point[0];
          data.phi[4].ref_grad[0] = -DataType(4) * point[1];
          data.phi[4].ref_grad[1] = -DataType(4) * (DataType(2)*point[1] + point[0] - DataType(1));
          data.phi[5].ref_grad[0] = -DataType(4) * (DataType(2)*point[0] + point[1] - DataType(1));
          data.phi[5].ref_grad[1] = -DataType(4) * point[0];
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_hessians(
          EvalData_& data,
          const PointType_& DOXY(point))
        {
          // vertex dofs
          data.phi[0].ref_hess[0][0] = data.phi[0].ref_hess[1][1] =
          data.phi[0].ref_hess[0][1] = data.phi[0].ref_hess[1][0] = DataType(4);
          data.phi[1].ref_hess[0][0] = DataType(4);
          data.phi[1].ref_hess[1][1] = data.phi[1].ref_hess[0][1] = data.phi[1].ref_hess[1][0] = DataType(0);
          data.phi[2].ref_hess[1][1] = DataType(4);
          data.phi[2].ref_hess[0][0] = data.phi[2].ref_hess[0][1] = data.phi[2].ref_hess[1][0] = DataType(0);
          // edge dofs
          data.phi[3].ref_hess[0][0] = data.phi[3].ref_hess[1][1] = DataType(0);
          data.phi[3].ref_hess[0][1] = data.phi[3].ref_hess[1][0] = DataType(4);
          data.phi[4].ref_hess[0][0] = DataType(0);
          data.phi[4].ref_hess[1][1] = -DataType(8);
          data.phi[4].ref_hess[0][1] = data.phi[4].ref_hess[1][0] = -DataType(4);
          data.phi[5].ref_hess[0][0] = -DataType(8);
          data.phi[5].ref_hess[1][1] = DataType(0);
          data.phi[5].ref_hess[0][1] = data.phi[5].ref_hess[1][0] = -DataType(4);
        }
      };

      template<typename PointType_, typename DataType_>
      struct EvalHelper<PointType_, DataType_, Shape::Simplex<3>>
      {
        typedef DataType_ DataType;

        CUDA_HOST_DEVICE static constexpr int get_num_local_dofs()
        {
          return 10;
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_values(EvalData_& data, const PointType_& point)
        {
          // vertex dofs
          data.phi[0].ref_value = DataType(2) * (point[0] + point[1] + point[2] - DataType(0.5)) * (point[0] + point[1] + point[2] - DataType(1));
          data.phi[1].ref_value = DataType(2) * point[0] * (point[0] - DataType(0.5));
          data.phi[2].ref_value = DataType(2) * point[1] * (point[1] - DataType(0.5));
          data.phi[3].ref_value = DataType(2) * point[2] * (point[2] - DataType(0.5));
          // edge dofs
          data.phi[4].ref_value = DataType(4) * point[0] * (DataType(1) - point[0] - point[1] - point[2]);
          data.phi[5].ref_value = DataType(4) * point[1] * (DataType(1) - point[0] - point[1] - point[2]);
          data.phi[6].ref_value = DataType(4) * point[2] * (DataType(1) - point[0] - point[1] - point[2]);
          data.phi[7].ref_value = DataType(4) * point[0] * point[1];
          data.phi[8].ref_value = DataType(4) * point[0] * point[2];
          data.phi[9].ref_value = DataType(4) * point[1] * point[2];
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_gradients(
          EvalData_& data,
          const PointType_& point)
        {
          // vertex dofs
          data.phi[0].ref_grad[0] = DataType(4) * (point[0] + point[1] + point[2]) - DataType(3);
          data.phi[0].ref_grad[1] = DataType(4) * (point[0] + point[1] + point[2]) - DataType(3);
          data.phi[0].ref_grad[2] = DataType(4) * (point[0] + point[1] + point[2]) - DataType(3);
          data.phi[1].ref_grad[0] = DataType(4) * point[0] - DataType(1);
          data.phi[1].ref_grad[1] = DataType(0);
          data.phi[1].ref_grad[2] = DataType(0);
          data.phi[2].ref_grad[0] = DataType(0);
          data.phi[2].ref_grad[1] = DataType(4) * point[1] - DataType(1);
          data.phi[2].ref_grad[2] = DataType(0);
          data.phi[3].ref_grad[0] = DataType(0);
          data.phi[3].ref_grad[1] = DataType(0);
          data.phi[3].ref_grad[2] = DataType(4) * point[2] - DataType(1);
          // edge dofs
          data.phi[4].ref_grad[0] = -DataType(4) * (DataType(2)*point[0] + point[1] + point[2] - DataType(1));
          data.phi[4].ref_grad[1] = -DataType(4) * point[0];
          data.phi[4].ref_grad[2] = -DataType(4) * point[0];
          data.phi[5].ref_grad[0] = -DataType(4) * point[1];
          data.phi[5].ref_grad[1] = -DataType(4) * (DataType(2)*point[1] + point[0] + point[2] - DataType(1));
          data.phi[5].ref_grad[2] = -DataType(4) * point[1];
          data.phi[6].ref_grad[0] = -DataType(4) * point[2];
          data.phi[6].ref_grad[1] = -DataType(4) * point[2];
          data.phi[6].ref_grad[2] = -DataType(4) * (DataType(2)*point[2] + point[0] + point[1] - DataType(1));
          data.phi[7].ref_grad[0] = DataType(4) * point[1];
          data.phi[7].ref_grad[1] = DataType(4) * point[0];
          data.phi[7].ref_grad[2] = DataType(0);
          data.phi[8].ref_grad[0] = DataType(4) * point[2];
          data.phi[8].ref_grad[1] = DataType(0);
          data.phi[8].ref_grad[2] = DataType(4) * point[0];
          data.phi[9].ref_grad[0] = DataType(0);
          data.phi[9].ref_grad[1] = DataType(4) * point[2];
          data.phi[9].ref_grad[2] = DataType(4) * point[1];
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_hessians(
          EvalData_& data,
          const PointType_& DOXY(point))
        {
          // vertex dofs
          data.phi[0].ref_hess[0][0] = DataType(4);
          data.phi[0].ref_hess[1][1] = DataType(4);
          data.phi[0].ref_hess[2][2] = DataType(4);
          data.phi[0].ref_hess[0][1] = data.phi[0].ref_hess[1][0] = DataType(4);
          data.phi[0].ref_hess[0][2] = data.phi[0].ref_hess[2][0] = DataType(4);
          data.phi[0].ref_hess[1][2] = data.phi[0].ref_hess[2][1] = DataType(4);

          data.phi[1].ref_hess[0][0] = DataType(4);
          data.phi[1].ref_hess[1][1] = DataType(0);
          data.phi[1].ref_hess[2][2] = DataType(0);
          data.phi[1].ref_hess[0][1] = data.phi[1].ref_hess[1][0] = DataType(0);
          data.phi[1].ref_hess[0][2] = data.phi[1].ref_hess[2][0] = DataType(0);
          data.phi[1].ref_hess[1][2] = data.phi[1].ref_hess[2][1] = DataType(0);

          data.phi[2].ref_hess[0][0] = DataType(0);
          data.phi[2].ref_hess[1][1] = DataType(4);
          data.phi[2].ref_hess[2][2] = DataType(0);
          data.phi[2].ref_hess[0][1] = data.phi[2].ref_hess[1][0] = DataType(0);
          data.phi[2].ref_hess[0][2] = data.phi[2].ref_hess[2][0] = DataType(0);
          data.phi[2].ref_hess[1][2] = data.phi[2].ref_hess[2][1] = DataType(0);

          data.phi[3].ref_hess[0][0] = DataType(0);
          data.phi[3].ref_hess[1][1] = DataType(0);
          data.phi[3].ref_hess[2][2] = DataType(4);
          data.phi[3].ref_hess[0][1] = data.phi[3].ref_hess[1][0] = DataType(0);
          data.phi[3].ref_hess[0][2] = data.phi[3].ref_hess[2][0] = DataType(0);
          data.phi[3].ref_hess[1][2] = data.phi[3].ref_hess[2][1] = DataType(0);

          // edge dofs
          data.phi[4].ref_hess[0][0] = -DataType(8);
          data.phi[4].ref_hess[1][1] = DataType(0);
          data.phi[4].ref_hess[2][2] = DataType(0);
          data.phi[4].ref_hess[0][1] = data.phi[4].ref_hess[1][0] = -DataType(4);
          data.phi[4].ref_hess[0][2] = data.phi[4].ref_hess[2][0] = -DataType(4);
          data.phi[4].ref_hess[1][2] = data.phi[4].ref_hess[2][1] = DataType(0);

          data.phi[5].ref_hess[0][0] = DataType(0);
          data.phi[5].ref_hess[1][1] = -DataType(8);
          data.phi[5].ref_hess[2][2] = DataType(0);
          data.phi[5].ref_hess[0][1] = data.phi[5].ref_hess[1][0] = -DataType(4);
          data.phi[5].ref_hess[0][2] = data.phi[5].ref_hess[2][0] = DataType(0);
          data.phi[5].ref_hess[1][2] = data.phi[5].ref_hess[2][1] = -DataType(4);

          data.phi[6].ref_hess[0][0] = DataType(0);
          data.phi[6].ref_hess[1][1] = DataType(0);
          data.phi[6].ref_hess[2][2] = -DataType(8);
          data.phi[6].ref_hess[0][1] = data.phi[6].ref_hess[1][0] = DataType(0);
          data.phi[6].ref_hess[0][2] = data.phi[6].ref_hess[2][0] = -DataType(4);
          data.phi[6].ref_hess[1][2] = data.phi[6].ref_hess[2][1] = -DataType(4);

          data.phi[7].ref_hess[0][0] = DataType(0);
          data.phi[7].ref_hess[1][1] = DataType(0);
          data.phi[7].ref_hess[2][2] = DataType(0);
          data.phi[7].ref_hess[0][1] = data.phi[7].ref_hess[1][0] = DataType(4);
          data.phi[7].ref_hess[0][2] = data.phi[7].ref_hess[2][0] = DataType(0);
          data.phi[7].ref_hess[1][2] = data.phi[7].ref_hess[2][1] = DataType(0);

          data.phi[8].ref_hess[0][0] = DataType(0);
          data.phi[8].ref_hess[1][1] = DataType(0);
          data.phi[8].ref_hess[2][2] = DataType(0);
          data.phi[8].ref_hess[0][1] = data.phi[8].ref_hess[1][0] = DataType(0);
          data.phi[8].ref_hess[0][2] = data.phi[8].ref_hess[2][0] = DataType(4);
          data.phi[8].ref_hess[1][2] = data.phi[8].ref_hess[2][1] = DataType(0);

          data.phi[9].ref_hess[0][0] = DataType(0);
          data.phi[9].ref_hess[1][1] = DataType(0);
          data.phi[9].ref_hess[2][2] = DataType(0);
          data.phi[9].ref_hess[0][1] = data.phi[9].ref_hess[1][0] = DataType(0);
          data.phi[9].ref_hess[0][2] = data.phi[9].ref_hess[2][0] = DataType(0);
          data.phi[9].ref_hess[1][2] = data.phi[9].ref_hess[2][1] = DataType(4);
        }
      };

      template<typename PointType_, typename DataType_>
      struct EvalHelper<PointType_, DataType_, Shape::Hypercube<1>>
      {
        typedef DataType_ DataType;

        CUDA_HOST_DEVICE static constexpr int get_num_local_dofs()
        {
          return 3;
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_values(EvalData_& data, const PointType_& point)
        {
          data.phi[0].ref_value = Intern::p0(point[0]);
          data.phi[1].ref_value = Intern::p1(point[0]);
          data.phi[2].ref_value = Intern::p2(point[0]);
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_gradients(
          EvalData_& data,
          const PointType_& point)
        {
          data.phi[0].ref_grad[0] = Intern::d1p0(point[0]);
          data.phi[1].ref_grad[0] = Intern::d1p1(point[0]);
          data.phi[2].ref_grad[0] = Intern::d1p2(point[0]);
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_hessians(
          EvalData_& data,
          const PointType_& point)
        {
          data.phi[0].ref_hess[0][0] = Intern::d2p0(point[0]);
          data.phi[1].ref_hess[0][0] = Intern::d2p1(point[0]);
          data.phi[2].ref_hess[0][0] = Intern::d2p2(point[0]);
        }
      };

      template<typename PointType_, typename DataType_>
      struct EvalHelper<PointType_, DataType_, Shape::Hypercube<2>>
      {
        typedef DataType_ DataType;

        CUDA_HOST_DEVICE static constexpr int get_num_local_dofs()
        {
          return 9;
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_values(EvalData_& data, const PointType_& point)
        {
          using namespace Lagrange2::Intern;

          // vertex dofs
          data.phi[0].ref_value = p0(point[0]) * p0(point[1]);
          data.phi[1].ref_value = p1(point[0]) * p0(point[1]);
          data.phi[2].ref_value = p0(point[0]) * p1(point[1]);
          data.phi[3].ref_value = p1(point[0]) * p1(point[1]);
          // edge dofs
          data.phi[4].ref_value = p2(point[0]) * p0(point[1]);
          data.phi[5].ref_value = p2(point[0]) * p1(point[1]);
          data.phi[6].ref_value = p0(point[0]) * p2(point[1]);
          data.phi[7].ref_value = p1(point[0]) * p2(point[1]);
          // center dof
          data.phi[8].ref_value = p2(point[0]) * p2(point[1]);
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_gradients(
          EvalData_& data,
          const PointType_& point)
        {
          using namespace Lagrange2::Intern;

          // vertex dofs
          data.phi[0].ref_grad[0] = d1p0(point[0]) * p0(point[1]);
          data.phi[0].ref_grad[1] = p0(point[0]) * d1p0(point[1]);
          data.phi[1].ref_grad[0] = d1p1(point[0]) * p0(point[1]);
          data.phi[1].ref_grad[1] = p1(point[0]) * d1p0(point[1]);
          data.phi[2].ref_grad[0] = d1p0(point[0]) * p1(point[1]);
          data.phi[2].ref_grad[1] = p0(point[0]) * d1p1(point[1]);
          data.phi[3].ref_grad[0] = d1p1(point[0]) * p1(point[1]);
          data.phi[3].ref_grad[1] = p1(point[0]) * d1p1(point[1]);
          // edge dofs
          data.phi[4].ref_grad[0] = d1p2(point[0]) * p0(point[1]);
          data.phi[4].ref_grad[1] = p2(point[0]) * d1p0(point[1]);
          data.phi[5].ref_grad[0] = d1p2(point[0]) * p1(point[1]);
          data.phi[5].ref_grad[1] = p2(point[0]) * d1p1(point[1]);
          data.phi[6].ref_grad[0] = d1p0(point[0]) * p2(point[1]);
          data.phi[6].ref_grad[1] = p0(point[0]) * d1p2(point[1]);
          data.phi[7].ref_grad[0] = d1p1(point[0]) * p2(point[1]);
          data.phi[7].ref_grad[1] = p1(point[0]) * d1p2(point[1]);
          // center dof
          data.phi[8].ref_grad[0] = d1p2(point[0]) * p2(point[1]);
          data.phi[8].ref_grad[1] = p2(point[0]) * d1p2(point[1]);
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_hessians(
          EvalData_& data,
          const PointType_& point)
        {
          using namespace Lagrange2::Intern;

          // vertex dofs
          data.phi[0].ref_hess[0][0] = d2p0(point[0]) * p0(point[1]);
          data.phi[0].ref_hess[1][1] = p0(point[0]) * d2p0(point[1]);
          data.phi[0].ref_hess[1][0] =
          data.phi[0].ref_hess[0][1] = d1p0(point[0]) * d1p0(point[1]);
          data.phi[1].ref_hess[0][0] = d2p1(point[0]) * p0(point[1]);
          data.phi[1].ref_hess[1][1] = p1(point[0]) * d2p0(point[1]);
          data.phi[1].ref_hess[1][0] =
          data.phi[1].ref_hess[0][1] = d1p1(point[0]) * d1p0(point[1]);
          data.phi[2].ref_hess[0][0] = d2p0(point[0]) * p1(point[1]);
          data.phi[2].ref_hess[1][1] = p0(point[0]) * d2p1(point[1]);
          data.phi[2].ref_hess[1][0] =
          data.phi[2].ref_hess[0][1] = d1p0(point[0]) * d1p1(point[1]);
          data.phi[3].ref_hess[0][0] = d2p1(point[0]) * p1(point[1]);
          data.phi[3].ref_hess[1][1] = p1(point[0]) * d2p1(point[1]);
          data.phi[3].ref_hess[1][0] =
          data.phi[3].ref_hess[0][1] = d1p1(point[0]) * d1p1(point[1]);
          // edge dofs
          data.phi[4].ref_hess[0][0] = d2p2(point[0]) * p0(point[1]);
          data.phi[4].ref_hess[1][1] = p2(point[0]) * d2p0(point[1]);
          data.phi[4].ref_hess[1][0] =
          data.phi[4].ref_hess[0][1] = d1p2(point[0]) * d1p0(point[1]);
          data.phi[5].ref_hess[0][0] = d2p2(point[0]) * p1(point[1]);
          data.phi[5].ref_hess[1][1] = p2(point[0]) * d2p1(point[1]);
          data.phi[5].ref_hess[1][0] =
          data.phi[5].ref_hess[0][1] = d1p2(point[0]) * d1p1(point[1]);
          data.phi[6].ref_hess[0][0] = d2p0(point[0]) * p2(point[1]);
          data.phi[6].ref_hess[1][1] = p0(point[0]) * d2p2(point[1]);
          data.phi[6].ref_hess[1][0] =
          data.phi[6].ref_hess[0][1] = d1p0(point[0]) * d1p2(point[1]);
          data.phi[7].ref_hess[0][0] = d2p1(point[0]) * p2(point[1]);
          data.phi[7].ref_hess[1][1] = p1(point[0]) * d2p2(point[1]);
          data.phi[7].ref_hess[1][0] =
          data.phi[7].ref_hess[0][1] = d1p1(point[0]) * d1p2(point[1]);
          // center dof
          data.phi[8].ref_hess[0][0] = d2p2(point[0]) * p2(point[1]);
          data.phi[8].ref_hess[1][1] = p2(point[0]) * d2p2(point[1]);
          data.phi[8].ref_hess[1][0] =
          data.phi[8].ref_hess[0][1] = d1p2(point[0]) * d1p2(point[1]);
        }
      };

      template<typename PointType_, typename DataType_>
      struct EvalHelper<PointType_, DataType_, Shape::Hypercube<3>>
      {
        typedef DataType_ DataType;

        CUDA_HOST_DEVICE static constexpr int get_num_local_dofs()
        {
          return 27;
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_values(EvalData_& data, const PointType_& point)
        {
          using namespace Lagrange2::Intern;

          // vertex dofs
          data.phi[0].ref_value = p0(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[1].ref_value = p1(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[2].ref_value = p0(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[3].ref_value = p1(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[4].ref_value = p0(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[5].ref_value = p1(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[6].ref_value = p0(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[7].ref_value = p1(point[0]) * p1(point[1]) * p1(point[2]);
          // edge dofs
          data.phi[ 8].ref_value = p2(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ 9].ref_value = p2(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[10].ref_value = p2(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[11].ref_value = p2(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[12].ref_value = p0(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[13].ref_value = p1(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[14].ref_value = p0(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[15].ref_value = p1(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[16].ref_value = p0(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[17].ref_value = p1(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[18].ref_value = p0(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[19].ref_value = p1(point[0]) * p1(point[1]) * p2(point[2]);
          // face dofs
          data.phi[20].ref_value = p2(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[21].ref_value = p2(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[22].ref_value = p2(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[23].ref_value = p2(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[24].ref_value = p0(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[25].ref_value = p1(point[0]) * p2(point[1]) * p2(point[2]);
          // center dof
          data.phi[26].ref_value = p2(point[0]) * p2(point[1]) * p2(point[2]);
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_gradients(
          EvalData_& data,
          const PointType_& point)
        {
          using namespace Lagrange2::Intern;
          // vertex dofs
          data.phi[0].ref_grad[0] = d1p0(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[0].ref_grad[1] = p0(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[0].ref_grad[2] = p0(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[1].ref_grad[0] = d1p1(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[1].ref_grad[1] = p1(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[1].ref_grad[2] = p1(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[2].ref_grad[0] = d1p0(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[2].ref_grad[1] = p0(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[2].ref_grad[2] = p0(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[3].ref_grad[0] = d1p1(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[3].ref_grad[1] = p1(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[3].ref_grad[2] = p1(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[4].ref_grad[0] = d1p0(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[4].ref_grad[1] = p0(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[4].ref_grad[2] = p0(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[5].ref_grad[0] = d1p1(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[5].ref_grad[1] = p1(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[5].ref_grad[2] = p1(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[6].ref_grad[0] = d1p0(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[6].ref_grad[1] = p0(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[6].ref_grad[2] = p0(point[0]) * p1(point[1]) * d1p1(point[2]);
          data.phi[7].ref_grad[0] = d1p1(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[7].ref_grad[1] = p1(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[7].ref_grad[2] = p1(point[0]) * p1(point[1]) * d1p1(point[2]);
          // edge dofs
          data.phi[ 8].ref_grad[0] = d1p2(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ 8].ref_grad[1] = p2(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[ 8].ref_grad[2] = p2(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[ 9].ref_grad[0] = d1p2(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[ 9].ref_grad[1] = p2(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[ 9].ref_grad[2] = p2(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[10].ref_grad[0] = d1p2(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[10].ref_grad[1] = p2(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[10].ref_grad[2] = p2(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[11].ref_grad[0] = d1p2(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[11].ref_grad[1] = p2(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[11].ref_grad[2] = p2(point[0]) * p1(point[1]) * d1p1(point[2]);
          data.phi[12].ref_grad[0] = d1p0(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[12].ref_grad[1] = p0(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[12].ref_grad[2] = p0(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[13].ref_grad[0] = d1p1(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[13].ref_grad[1] = p1(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[13].ref_grad[2] = p1(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[14].ref_grad[0] = d1p0(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[14].ref_grad[1] = p0(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[14].ref_grad[2] = p0(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[15].ref_grad[0] = d1p1(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[15].ref_grad[1] = p1(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[15].ref_grad[2] = p1(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[16].ref_grad[0] = d1p0(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[16].ref_grad[1] = p0(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[16].ref_grad[2] = p0(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[17].ref_grad[0] = d1p1(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[17].ref_grad[1] = p1(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[17].ref_grad[2] = p1(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[18].ref_grad[0] = d1p0(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[18].ref_grad[1] = p0(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[18].ref_grad[2] = p0(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[19].ref_grad[0] = d1p1(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[19].ref_grad[1] = p1(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[19].ref_grad[2] = p1(point[0]) * p1(point[1]) * d1p2(point[2]);
          // face dofs
          data.phi[20].ref_grad[0] = d1p2(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[20].ref_grad[1] = p2(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[20].ref_grad[2] = p2(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[21].ref_grad[0] = d1p2(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[21].ref_grad[1] = p2(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[21].ref_grad[2] = p2(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[22].ref_grad[0] = d1p2(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[22].ref_grad[1] = p2(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[22].ref_grad[2] = p2(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[23].ref_grad[0] = d1p2(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[23].ref_grad[1] = p2(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[23].ref_grad[2] = p2(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[24].ref_grad[0] = d1p0(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[24].ref_grad[1] = p0(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[24].ref_grad[2] = p0(point[0]) * p2(point[1]) * d1p2(point[2]);
          data.phi[25].ref_grad[0] = d1p1(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[25].ref_grad[1] = p1(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[25].ref_grad[2] = p1(point[0]) * p2(point[1]) * d1p2(point[2]);
          // center dof
          data.phi[26].ref_grad[0] = d1p2(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[26].ref_grad[1] = p2(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[26].ref_grad[2] = p2(point[0]) * p2(point[1]) * d1p2(point[2]);
        }

        template<typename EvalData_>
        CUDA_HOST_DEVICE static inline void eval_ref_hessians(
          EvalData_& data,
          const PointType_& point)
        {
          using namespace Lagrange2::Intern;
          // vertex dofs
          data.phi[ 0].ref_hess[0][0] = d2p0(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ 0].ref_hess[1][1] = p0(point[0]) * d2p0(point[1]) * p0(point[2]);
          data.phi[ 0].ref_hess[2][2] = p0(point[0]) * p0(point[1]) * d2p0(point[2]);
          data.phi[ 0].ref_hess[1][0] =
          data.phi[ 0].ref_hess[0][1] = d1p0(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[ 0].ref_hess[2][0] =
          data.phi[ 0].ref_hess[0][2] = d1p0(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[ 0].ref_hess[2][1] =
          data.phi[ 0].ref_hess[1][2] = p0(point[0]) * d1p0(point[1]) * d1p0(point[2]);
          data.phi[ 1].ref_hess[0][0] = d2p1(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ 1].ref_hess[1][1] = p1(point[0]) * d2p0(point[1]) * p0(point[2]);
          data.phi[ 1].ref_hess[2][2] = p1(point[0]) * p0(point[1]) * d2p0(point[2]);
          data.phi[ 1].ref_hess[1][0] =
          data.phi[ 1].ref_hess[0][1] = d1p1(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[ 1].ref_hess[2][0] =
          data.phi[ 1].ref_hess[0][2] = d1p1(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[ 1].ref_hess[2][1] =
          data.phi[ 1].ref_hess[1][2] = p1(point[0]) * d1p0(point[1]) * d1p0(point[2]);
          data.phi[ 2].ref_hess[0][0] = d2p0(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[ 2].ref_hess[1][1] = p0(point[0]) * d2p1(point[1]) * p0(point[2]);
          data.phi[ 2].ref_hess[2][2] = p0(point[0]) * p1(point[1]) * d2p0(point[2]);
          data.phi[ 2].ref_hess[1][0] =
          data.phi[ 2].ref_hess[0][1] = d1p0(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[ 2].ref_hess[2][0] =
          data.phi[ 2].ref_hess[0][2] = d1p0(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[ 2].ref_hess[2][1] =
          data.phi[ 2].ref_hess[1][2] = p0(point[0]) * d1p1(point[1]) * d1p0(point[2]);
          data.phi[ 3].ref_hess[0][0] = d2p1(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[ 3].ref_hess[1][1] = p1(point[0]) * d2p1(point[1]) * p0(point[2]);
          data.phi[ 3].ref_hess[2][2] = p1(point[0]) * p1(point[1]) * d2p0(point[2]);
          data.phi[ 3].ref_hess[1][0] =
          data.phi[ 3].ref_hess[0][1] = d1p1(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[ 3].ref_hess[2][0] =
          data.phi[ 3].ref_hess[0][2] = d1p1(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[ 3].ref_hess[2][1] =
          data.phi[ 3].ref_hess[1][2] = p1(point[0]) * d1p1(point[1]) * d1p0(point[2]);
          data.phi[ 4].ref_hess[0][0] = d2p0(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[ 4].ref_hess[1][1] = p0(point[0]) * d2p0(point[1]) * p1(point[2]);
          data.phi[ 4].ref_hess[2][2] = p0(point[0]) * p0(point[1]) * d2p1(point[2]);
          data.phi[ 4].ref_hess[1][0] =
          data.phi[ 4].ref_hess[0][1] = d1p0(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[ 4].ref_hess[2][0] =
          data.phi[ 4].ref_hess[0][2] = d1p0(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[ 4].ref_hess[2][1] =
          data.phi[ 4].ref_hess[1][2] = p0(point[0]) * d1p0(point[1]) * d1p1(point[2]);
          data.phi[ 5].ref_hess[0][0] = d2p1(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[ 5].ref_hess[1][1] = p1(point[0]) * d2p0(point[1]) * p1(point[2]);
          data.phi[ 5].ref_hess[2][2] = p1(point[0]) * p0(point[1]) * d2p1(point[2]);
          data.phi[ 5].ref_hess[1][0] =
          data.phi[ 5].ref_hess[0][1] = d1p1(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[ 5].ref_hess[2][0] =
          data.phi[ 5].ref_hess[0][2] = d1p1(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[ 5].ref_hess[2][1] =
          data.phi[ 5].ref_hess[1][2] = p1(point[0]) * d1p0(point[1]) * d1p1(point[2]);
          data.phi[ 6].ref_hess[0][0] = d2p0(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[ 6].ref_hess[1][1] = p0(point[0]) * d2p1(point[1]) * p1(point[2]);
          data.phi[ 6].ref_hess[2][2] = p0(point[0]) * p1(point[1]) * d2p1(point[2]);
          data.phi[ 6].ref_hess[1][0] =
          data.phi[ 6].ref_hess[0][1] = d1p0(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[ 6].ref_hess[2][0] =
          data.phi[ 6].ref_hess[0][2] = d1p0(point[0]) * p1(point[1]) * d1p1(point[2]);
          data.phi[ 6].ref_hess[2][1] =
          data.phi[ 6].ref_hess[1][2] = p0(point[0]) * d1p1(point[1]) * d1p1(point[2]);
          data.phi[ 7].ref_hess[0][0] = d2p1(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[ 7].ref_hess[1][1] = p1(point[0]) * d2p1(point[1]) * p1(point[2]);
          data.phi[ 7].ref_hess[2][2] = p1(point[0]) * p1(point[1]) * d2p1(point[2]);
          data.phi[ 7].ref_hess[1][0] =
          data.phi[ 7].ref_hess[0][1] = d1p1(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[ 7].ref_hess[2][0] =
          data.phi[ 7].ref_hess[0][2] = d1p1(point[0]) * p1(point[1]) * d1p1(point[2]);
          data.phi[ 7].ref_hess[2][1] =
          data.phi[ 7].ref_hess[1][2] = p1(point[0]) * d1p1(point[1]) * d1p1(point[2]);

          // edge dofs
          data.phi[ 8].ref_hess[0][0] = d2p2(point[0]) * p0(point[1]) * p0(point[2]);
          data.phi[ 8].ref_hess[1][1] = p2(point[0]) * d2p0(point[1]) * p0(point[2]);
          data.phi[ 8].ref_hess[2][2] = p2(point[0]) * p0(point[1]) * d2p0(point[2]);
          data.phi[ 8].ref_hess[1][0] =
          data.phi[ 8].ref_hess[0][1] = d1p2(point[0]) * d1p0(point[1]) * p0(point[2]);
          data.phi[ 8].ref_hess[2][0] =
          data.phi[ 8].ref_hess[0][2] = d1p2(point[0]) * p0(point[1]) * d1p0(point[2]);
          data.phi[ 8].ref_hess[2][1] =
          data.phi[ 8].ref_hess[1][2] = p2(point[0]) * d1p0(point[1]) * d1p0(point[2]);
          data.phi[ 9].ref_hess[0][0] = d2p2(point[0]) * p1(point[1]) * p0(point[2]);
          data.phi[ 9].ref_hess[1][1] = p2(point[0]) * d2p1(point[1]) * p0(point[2]);
          data.phi[ 9].ref_hess[2][2] = p2(point[0]) * p1(point[1]) * d2p0(point[2]);
          data.phi[ 9].ref_hess[1][0] =
          data.phi[ 9].ref_hess[0][1] = d1p2(point[0]) * d1p1(point[1]) * p0(point[2]);
          data.phi[ 9].ref_hess[2][0] =
          data.phi[ 9].ref_hess[0][2] = d1p2(point[0]) * p1(point[1]) * d1p0(point[2]);
          data.phi[ 9].ref_hess[2][1] =
          data.phi[ 9].ref_hess[1][2] = p2(point[0]) * d1p1(point[1]) * d1p0(point[2]);
          data.phi[10].ref_hess[0][0] = d2p2(point[0]) * p0(point[1]) * p1(point[2]);
          data.phi[10].ref_hess[1][1] = p2(point[0]) * d2p0(point[1]) * p1(point[2]);
          data.phi[10].ref_hess[2][2] = p2(point[0]) * p0(point[1]) * d2p1(point[2]);
          data.phi[10].ref_hess[1][0] =
          data.phi[10].ref_hess[0][1] = d1p2(point[0]) * d1p0(point[1]) * p1(point[2]);
          data.phi[10].ref_hess[2][0] =
          data.phi[10].ref_hess[0][2] = d1p2(point[0]) * p0(point[1]) * d1p1(point[2]);
          data.phi[10].ref_hess[2][1] =
          data.phi[10].ref_hess[1][2] = p2(point[0]) * d1p0(point[1]) * d1p1(point[2]);
          data.phi[11].ref_hess[0][0] = d2p2(point[0]) * p1(point[1]) * p1(point[2]);
          data.phi[11].ref_hess[1][1] = p2(point[0]) * d2p1(point[1]) * p1(point[2]);
          data.phi[11].ref_hess[2][2] = p2(point[0]) * p1(point[1]) * d2p1(point[2]);
          data.phi[11].ref_hess[1][0] =
          data.phi[11].ref_hess[0][1] = d1p2(point[0]) * d1p1(point[1]) * p1(point[2]);
          data.phi[11].ref_hess[2][0] =
          data.phi[11].ref_hess[0][2] = d1p2(point[0]) * p1(point[1]) * d1p1(point[2]);
          data.phi[11].ref_hess[2][1] =
          data.phi[11].ref_hess[1][2] = p2(point[0]) * d1p1(point[1]) * d1p1(point[2]);
          data.phi[12].ref_hess[0][0] = d2p0(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[12].ref_hess[1][1] = p0(point[0]) * d2p2(point[1]) * p0(point[2]);
          data.phi[12].ref_hess[2][2] = p0(point[0]) * p2(point[1]) * d2p0(point[2]);
          data.phi[12].ref_hess[1][0] =
          data.phi[12].ref_hess[0][1] = d1p0(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[12].ref_hess[2][0] =
          data.phi[12].ref_hess[0][2] = d1p0(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[12].ref_hess[2][1] =
          data.phi[12].ref_hess[1][2] = p0(point[0]) * d1p2(point[1]) * d1p0(point[2]);
          data.phi[13].ref_hess[0][0] = d2p1(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[13].ref_hess[1][1] = p1(point[0]) * d2p2(point[1]) * p0(point[2]);
          data.phi[13].ref_hess[2][2] = p1(point[0]) * p2(point[1]) * d2p0(point[2]);
          data.phi[13].ref_hess[1][0] =
          data.phi[13].ref_hess[0][1] = d1p1(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[13].ref_hess[2][0] =
          data.phi[13].ref_hess[0][2] = d1p1(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[13].ref_hess[2][1] =
          data.phi[13].ref_hess[1][2] = p1(point[0]) * d1p2(point[1]) * d1p0(point[2]);
          data.phi[14].ref_hess[0][0] = d2p0(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[14].ref_hess[1][1] = p0(point[0]) * d2p2(point[1]) * p1(point[2]);
          data.phi[14].ref_hess[2][2] = p0(point[0]) * p2(point[1]) * d2p1(point[2]);
          data.phi[14].ref_hess[1][0] =
          data.phi[14].ref_hess[0][1] = d1p0(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[14].ref_hess[2][0] =
          data.phi[14].ref_hess[0][2] = d1p0(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[14].ref_hess[2][1] =
          data.phi[14].ref_hess[1][2] = p0(point[0]) * d1p2(point[1]) * d1p1(point[2]);
          data.phi[15].ref_hess[0][0] = d2p1(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[15].ref_hess[1][1] = p1(point[0]) * d2p2(point[1]) * p1(point[2]);
          data.phi[15].ref_hess[2][2] = p1(point[0]) * p2(point[1]) * d2p1(point[2]);
          data.phi[15].ref_hess[1][0] =
          data.phi[15].ref_hess[0][1] = d1p1(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[15].ref_hess[2][0] =
          data.phi[15].ref_hess[0][2] = d1p1(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[15].ref_hess[2][1] =
          data.phi[15].ref_hess[1][2] = p1(point[0]) * d1p2(point[1]) * d1p1(point[2]);
          data.phi[16].ref_hess[0][0] = d2p0(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[16].ref_hess[1][1] = p0(point[0]) * d2p0(point[1]) * p2(point[2]);
          data.phi[16].ref_hess[2][2] = p0(point[0]) * p0(point[1]) * d2p2(point[2]);
          data.phi[16].ref_hess[1][0] =
          data.phi[16].ref_hess[0][1] = d1p0(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[16].ref_hess[2][0] =
          data.phi[16].ref_hess[0][2] = d1p0(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[16].ref_hess[2][1] =
          data.phi[16].ref_hess[1][2] = p0(point[0]) * d1p0(point[1]) * d1p2(point[2]);
          data.phi[17].ref_hess[0][0] = d2p1(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[17].ref_hess[1][1] = p1(point[0]) * d2p0(point[1]) * p2(point[2]);
          data.phi[17].ref_hess[2][2] = p1(point[0]) * p0(point[1]) * d2p2(point[2]);
          data.phi[17].ref_hess[1][0] =
          data.phi[17].ref_hess[0][1] = d1p1(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[17].ref_hess[2][0] =
          data.phi[17].ref_hess[0][2] = d1p1(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[17].ref_hess[2][1] =
          data.phi[17].ref_hess[1][2] = p1(point[0]) * d1p0(point[1]) * d1p2(point[2]);
          data.phi[18].ref_hess[0][0] = d2p0(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[18].ref_hess[1][1] = p0(point[0]) * d2p1(point[1]) * p2(point[2]);
          data.phi[18].ref_hess[2][2] = p0(point[0]) * p1(point[1]) * d2p2(point[2]);
          data.phi[18].ref_hess[1][0] =
          data.phi[18].ref_hess[0][1] = d1p0(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[18].ref_hess[2][0] =
          data.phi[18].ref_hess[0][2] = d1p0(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[18].ref_hess[2][1] =
          data.phi[18].ref_hess[1][2] = p0(point[0]) * d1p1(point[1]) * d1p2(point[2]);
          data.phi[19].ref_hess[0][0] = d2p1(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[19].ref_hess[1][1] = p1(point[0]) * d2p1(point[1]) * p2(point[2]);
          data.phi[19].ref_hess[2][2] = p1(point[0]) * p1(point[1]) * d2p2(point[2]);
          data.phi[19].ref_hess[1][0] =
          data.phi[19].ref_hess[0][1] = d1p1(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[19].ref_hess[2][0] =
          data.phi[19].ref_hess[0][2] = d1p1(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[19].ref_hess[2][1] =
          data.phi[19].ref_hess[1][2] = p1(point[0]) * d1p1(point[1]) * d1p2(point[2]);

          // face dofs
          data.phi[20].ref_hess[0][0] = d2p2(point[0]) * p2(point[1]) * p0(point[2]);
          data.phi[20].ref_hess[1][1] = p2(point[0]) * d2p2(point[1]) * p0(point[2]);
          data.phi[20].ref_hess[2][2] = p2(point[0]) * p2(point[1]) * d2p0(point[2]);
          data.phi[20].ref_hess[1][0] =
          data.phi[20].ref_hess[0][1] = d1p2(point[0]) * d1p2(point[1]) * p0(point[2]);
          data.phi[20].ref_hess[2][0] =
          data.phi[20].ref_hess[0][2] = d1p2(point[0]) * p2(point[1]) * d1p0(point[2]);
          data.phi[20].ref_hess[2][1] =
          data.phi[20].ref_hess[1][2] = p2(point[0]) * d1p2(point[1]) * d1p0(point[2]);
          data.phi[21].ref_hess[0][0] = d2p2(point[0]) * p2(point[1]) * p1(point[2]);
          data.phi[21].ref_hess[1][1] = p2(point[0]) * d2p2(point[1]) * p1(point[2]);
          data.phi[21].ref_hess[2][2] = p2(point[0]) * p2(point[1]) * d2p1(point[2]);
          data.phi[21].ref_hess[1][0] =
          data.phi[21].ref_hess[0][1] = d1p2(point[0]) * d1p2(point[1]) * p1(point[2]);
          data.phi[21].ref_hess[2][0] =
          data.phi[21].ref_hess[0][2] = d1p2(point[0]) * p2(point[1]) * d1p1(point[2]);
          data.phi[21].ref_hess[2][1] =
          data.phi[21].ref_hess[1][2] = p2(point[0]) * d1p2(point[1]) * d1p1(point[2]);
          data.phi[22].ref_hess[0][0] = d2p2(point[0]) * p0(point[1]) * p2(point[2]);
          data.phi[22].ref_hess[1][1] = p2(point[0]) * d2p0(point[1]) * p2(point[2]);
          data.phi[22].ref_hess[2][2] = p2(point[0]) * p0(point[1]) * d2p2(point[2]);
          data.phi[22].ref_hess[1][0] =
          data.phi[22].ref_hess[0][1] = d1p2(point[0]) * d1p0(point[1]) * p2(point[2]);
          data.phi[22].ref_hess[2][0] =
          data.phi[22].ref_hess[0][2] = d1p2(point[0]) * p0(point[1]) * d1p2(point[2]);
          data.phi[22].ref_hess[2][1] =
          data.phi[22].ref_hess[1][2] = p2(point[0]) * d1p0(point[1]) * d1p2(point[2]);
          data.phi[23].ref_hess[0][0] = d2p2(point[0]) * p1(point[1]) * p2(point[2]);
          data.phi[23].ref_hess[1][1] = p2(point[0]) * d2p1(point[1]) * p2(point[2]);
          data.phi[23].ref_hess[2][2] = p2(point[0]) * p1(point[1]) * d2p2(point[2]);
          data.phi[23].ref_hess[1][0] =
          data.phi[23].ref_hess[0][1] = d1p2(point[0]) * d1p1(point[1]) * p2(point[2]);
          data.phi[23].ref_hess[2][0] =
          data.phi[23].ref_hess[0][2] = d1p2(point[0]) * p1(point[1]) * d1p2(point[2]);
          data.phi[23].ref_hess[2][1] =
          data.phi[23].ref_hess[1][2] = p2(point[0]) * d1p1(point[1]) * d1p2(point[2]);
          data.phi[24].ref_hess[0][0] = d2p0(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[24].ref_hess[1][1] = p0(point[0]) * d2p2(point[1]) * p2(point[2]);
          data.phi[24].ref_hess[2][2] = p0(point[0]) * p2(point[1]) * d2p2(point[2]);
          data.phi[24].ref_hess[1][0] =
          data.phi[24].ref_hess[0][1] = d1p0(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[24].ref_hess[2][0] =
          data.phi[24].ref_hess[0][2] = d1p0(point[0]) * p2(point[1]) * d1p2(point[2]);
          data.phi[24].ref_hess[2][1] =
          data.phi[24].ref_hess[1][2] = p0(point[0]) * d1p2(point[1]) * d1p2(point[2]);
          data.phi[25].ref_hess[0][0] = d2p1(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[25].ref_hess[1][1] = p1(point[0]) * d2p2(point[1]) * p2(point[2]);
          data.phi[25].ref_hess[2][2] = p1(point[0]) * p2(point[1]) * d2p2(point[2]);
          data.phi[25].ref_hess[1][0] =
          data.phi[25].ref_hess[0][1] = d1p1(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[25].ref_hess[2][0] =
          data.phi[25].ref_hess[0][2] = d1p1(point[0]) * p2(point[1]) * d1p2(point[2]);
          data.phi[25].ref_hess[2][1] =
          data.phi[25].ref_hess[1][2] = p1(point[0]) * d1p2(point[1]) * d1p2(point[2]);

          // center dof
          data.phi[26].ref_hess[0][0] = d2p2(point[0]) * p2(point[1]) * p2(point[2]);
          data.phi[26].ref_hess[1][1] = p2(point[0]) * d2p2(point[1]) * p2(point[2]);
          data.phi[26].ref_hess[2][2] = p2(point[0]) * p2(point[1]) * d2p2(point[2]);
          data.phi[26].ref_hess[1][0] =
          data.phi[26].ref_hess[0][1] = d1p2(point[0]) * d1p2(point[1]) * p2(point[2]);
          data.phi[26].ref_hess[2][0] =
          data.phi[26].ref_hess[0][2] = d1p2(point[0]) * p2(point[1]) * d1p2(point[2]);
          data.phi[26].ref_hess[2][1] =
          data.phi[26].ref_hess[1][2] = p2(point[0]) * d1p2(point[1]) * d1p2(point[2]);
        }
      };
    }
  }
}
