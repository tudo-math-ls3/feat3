// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include "base.hpp"
#include <kernel/analytic/lambda_function.hpp>

namespace Turb
{
  // ================================================================================================
  // Class BackWardConstant
  // ------------------------------------------------------------------------------------------------
  template<typename MeshNodeType>
  class BackWardConstant
  {
  private:
    const DataType NaN = Math::nan<DataType>();
    DataType tol = 1e-10;
    typedef Geometry::MeshPart<MeshType> MeshPartType;

    // Mesh dependent parameters
    DataType L = 2.0; // characteristic length of domain
    DataType wall_dist = 0.1;
    DataType l_max = 3.0;

  public:
    // Constants
    DataType u_mean = 1.0;
    DataType C_mu = 0.09;
    DataType c_bc = 0.007; // c_bc in [0.003,0.01]
    DataType l_0 = DataType(0.08);

    // Mesh Parts for Stokes Solver
    String no_flow_parts = "";
    String slip_parts = "";
    String flow_parts = "bnd:l1 | bnd:l bnd:b bnd:b1 bnd:t";
    std::deque<MeshPartType*> trace_mesh_parts;

    BackWardConstant(MeshNodeType& mesh_node)
    {
      trace_mesh_parts.push_back(mesh_node.find_mesh_part("bnd:t"));
      trace_mesh_parts.push_back(mesh_node.find_mesh_part("bnd:b"));
      trace_mesh_parts.push_back(mesh_node.find_mesh_part("bnd:b1"));
      trace_mesh_parts.push_back(mesh_node.find_mesh_part("bnd:l1"));
    }

    auto flow_func_2d ()
    {
      return Analytic::create_lambda_function_vector_2d(
        // x-component
        [this](DataType x, DataType y) -> DataType // bnd:l
        {
          if (std::abs(x - 0.0) < this->tol && y >= 1.0 && y <= 3.0)
          {
            return this->u_mean;
          }
          if (std::abs(x - 5.0) < this->tol && y > 0.0 && y < 1.0) return 0.0;                    // bnd:l1
          if (std::abs(y - 1.0) < this->tol && x > 0.0 && x < 5.0) return this->NaN;              // bnd:b1
          if (std::abs(y - 0.0) < this->tol && x > 5.0 && x <= 25.0) return this->NaN;            // bnd:b
          if (std::abs(y - 3.0) < this->tol && x > 0.0 && x <= 25.0) return this->NaN;            // bnd:t
          if (std::abs(x - 5.0) < this->tol && std::abs(y - 1.0) < this->tol) return 0.0;         //e1
          if (std::abs(x - 5.0) < this->tol && std::abs(y - 0.0) < this->tol) return 0.0;         //e2
          return 0.0;
        },
        // y-component
        [this](DataType x, DataType y) -> DataType
        {
          if (std::abs(x - 0.0) < this->tol && y >= 1.0 && y <= 3.0) return 0.0;                  // bnd:l
          if (std::abs(x - 5.0) < this->tol && y > 0.0 && y < 1.0) return this->NaN;              // bnd:l1
          if (std::abs(y - 1.0) < this->tol && x > 0.0 && x < 5.0) return 0.0;                    // bnd:b1
          if (std::abs(y - 0.0) < this->tol && x > 5.0 && x <= 25.0) return 0.0;                  // bnd:b
          if (std::abs(y - 3.0) < this->tol && x > 0.0 && x <= 25.0) return 0.0;                  // bnd:t
          if (std::abs(x - 5.0) < this->tol && std::abs(y - 1.0) < this->tol) return 0.0;         //e1
          if (std::abs(x - 5.0) < this->tol && std::abs(y - 0.0) < this->tol) return 0.0;         //e2
          return 0.0;
        }
      );
    }

    auto k_sol_func()
    {
      return FEAT::Analytic::create_lambda_function_scalar_2d(
        [this](DataType, DataType)
        {
          return this->c_bc * this->u_mean * this->u_mean;
        }
      );
    }

    auto e_sol_func()
    {
      return FEAT::Analytic::create_lambda_function_scalar_2d(
        [this](DataType, DataType)
        {
          return this->C_mu * pow((this->c_bc * this->u_mean*this->u_mean), DataType(3.0 / 2.0)) / this->l_0;
        }
      );
    }

    DataType get_l()
    {
      return L;
    }

    DataType get_wall_dist()
    {
      return wall_dist;
    }

    DataType get_l_max()
    {
      return l_max;
    }
  };
} // namespace Turb
