#pragma once
#ifndef KERNEL_SPACE_DETAILS_HPP
#define KERNEL_SPACE_DETAILS_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/tiny_algebra.hpp>

namespace FEAT
{
  namespace Space
  {
    struct ParametricEvalHelper
    {
      template<typename SpaceData_>
      CUDA_FUNC static void trans_values(SpaceData_& space_data)
      {
        // loop over all basis functions
        for(int i(0); i < SpaceData_::max_local_dofs; ++i)
        {
          // and copy the basis function value
          space_data.phi[i].value = space_data.phi[i].ref_value;
        }
      }

      template<typename SpaceData_, typename JacInvType_>
      CUDA_FUNC static void trans_gradients(SpaceData_& space_data, const JacInvType_& jac_inv)
      {
        // loop over all basis functions
        for(int i(0); i < SpaceData_::max_local_dofs; ++i)
        {
          // and apply the first-order chain rule
          space_data.phi[i].grad.set_vec_mat_mult(space_data.phi[i].ref_grad, jac_inv);
        }
      }

      template<typename SpaceData_, typename JacInvType_, typename HessInvType_>
      CUDA_FUNC static void trans_hessians(SpaceData_& space_data, const JacInvType_& jac_inv, const HessInvType_& hess_inv)
      {
        // loop over all basis functions
        for(int i(0); i < SpaceData_::max_local_dofs; ++i)
        {
          // and apply the second-order chain rule
          space_data.phi[i].hess.set_double_mat_mult(space_data.phi[i].ref_hess, jac_inv, jac_inv);
          space_data.phi[i].hess.add_vec_tensor_mult(space_data.phi[i].ref_grad, hess_inv);
        }
      }
    };
  }
}

#endif