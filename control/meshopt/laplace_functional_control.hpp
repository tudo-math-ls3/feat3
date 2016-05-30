#pragma once
#ifndef FEAT_CONTROL_MESHOPT_LAPLACE_FUNCTIONAL_CONTROL_HPP
#define FEAT_CONTROL_MESHOPT_LAPLACE_FUNCTIONAL_CONTROL_HPP 1
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>

#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/foundation_gate.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/nonlinear_functional.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/meshopt/laplace_smoother.hpp>

#include <control/domain/domain_control.hpp>
#include <control/meshopt/meshopt_control.hpp>
#include <control/meshopt/meshopt_solver_factory.hpp>

namespace FEAT
{
  namespace Control
  {
    namespace Meshopt
    {
      template<typename Space_>
      class LaplaceSmootherAssemblerLevel : public MeshoptAssemblerLevel<Space_>
      {
        public:
          typedef Control::Meshopt::MeshoptAssemblerLevel<Space_> BaseClass;
          static constexpr int shape_dim = BaseClass::MeshType::ShapeType::dimension;

          // Copy baseclass constructors
          using BaseClass::BaseClass;

          template<typename SystemLevel_>
          void assemble_system_matrix(SystemLevel_& sys_level)
          {
            // get the global matrix
            typename SystemLevel_::GlobalOperator& mat_glob = sys_level.op_sys;

            // get the local matrix
            typename SystemLevel_::LocalOperator::BaseClass& mat_loc = *mat_glob;

            // assemble matrix structure?
            if (mat_loc.empty())
            {
              Assembly::SymbolicAssembler::assemble_matrix_std1(mat_loc, this->trafo_space);
            }

            // assemble velocity laplace matrix
            {
              mat_loc.format();
              Assembly::Common::LaplaceOperator dudv_op;
              Assembly::BilinearOperatorAssembler::assemble_matrix1(
                mat_loc, dudv_op, this->trafo_space, this->cubature);
            }
          }
      }; // class LaplaceSmootherAssemblerLevel
    } // namespace Meshopt
  } // namespace Control
} // namespace FEAT

#endif// FEAT_CONTROL_MESHOPT_DUDV_FUNCTIONAL_CONTROL_HPP
