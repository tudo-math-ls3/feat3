// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "stokes_level.hpp"

#include <kernel/analytic/parsed_function.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>


namespace CCNDSimple
{
  // helper function: add meshparts to assembler
  template<typename Asm_>
  String aux_add_mesh_parts(Asm_& assm, DomainLevel& domain_level, const String& mesh_parts)
  {
    // split mesh part names and sort them
    std::deque<String> parts_in = mesh_parts.split_by_whitespaces();
    std::sort(parts_in.begin(), parts_in.end());

    // loop over all mesh parts
    for(String& s : parts_in)
    {
      // try to find the mesh part node
      auto* part_node = domain_level.get_mesh_node()->find_mesh_part_node(s);
      XASSERT(part_node != nullptr);

      // get the mesh part and add it to the assembler
      auto* mpart = part_node->get_mesh();
      if(mpart != nullptr)
        assm.add_mesh_part(*mpart);
    }

    // join all part names
    return stringify_join(parts_in, "##");
  }

  void StokesLevel::assemble_velocity_laplace_matrix(const String& cubature, const DataType nu, bool defo)
  {
    auto& loc_a = this->matrix_a.local();
    loc_a.format();
    Assembly::Common::LaplaceOperatorBlocked<dim> lapl_op;
    Assembly::Common::DuDvOperatorBlocked<dim> dudv_op;
    if(defo)
      Assembly::assemble_bilinear_operator_matrix_1(this->domain_level.domain_asm, loc_a, dudv_op, this->domain_level.space_velo, cubature, nu);
    else
      Assembly::assemble_bilinear_operator_matrix_1(this->domain_level.domain_asm, loc_a, lapl_op, this->domain_level.space_velo, cubature, nu);
  }

  void StokesLevel::assemble_velocity_mass_matrix(const String& cubature)
  {
    local_velo_mass_matrix = this->matrix_a.local().clone(LAFEM::CloneMode::Weak);
    local_velo_mass_matrix.format();

    Assembly::Common::IdentityOperatorBlocked<dim> id_op;
    Assembly::assemble_bilinear_operator_matrix_1(this->domain_level.domain_asm, local_velo_mass_matrix, id_op, this->domain_level.space_velo, cubature);
  }

  void StokesLevel::assemble_grad_div_matrices(const String& cubature)
  {
    BaseClass::assemble_grad_div_matrices(this->domain_level.domain_asm, this->domain_level.space_velo, this->domain_level.space_pres, cubature);
  }

  void StokesLevel::compile_local_matrix()
  {
    // convert local matrices
    this->local_matrix_sys.block_a() = this->matrix_a.convert_to_1();
    this->local_matrix_sys.block_b() = this->matrix_b.local().clone(LAFEM::CloneMode::Weak);
    this->local_matrix_sys.block_d() = this->matrix_d.local().clone(LAFEM::CloneMode::Weak);

    // apply velocity unit filters to matrices A and B
    const typename BaseClass::LocalVeloUnitFilterSeq& fil_seq = this->get_local_velo_unit_filter_seq();
    for(const auto& filter : fil_seq)
    {
      filter.second.filter_mat(this->local_matrix_sys.block_a());
      filter.second.filter_offdiag_row_mat(this->local_matrix_sys.block_b());
    }
  }

  void StokesLevel::assemble_inflow_bc(const String& mesh_parts, const String& formula)
  {
    Assembly::UnitFilterAssembler<MeshType> unit_filter_asm;
    Analytic::ParsedVectorFunction<dim> parsed_function(formula);

    // get the filter name
    String filter_name = aux_add_mesh_parts(unit_filter_asm, this->domain_level, mesh_parts);

    // get the unit filter and clear it
    BaseClass::LocalVeloUnitFilter& filter = this->get_local_velo_unit_filter_seq().find_or_add(filter_name);

    // clear and assemble
    filter.clear();
    unit_filter_asm.assemble(filter, this->domain_level.space_velo, parsed_function);
  }

  void StokesLevel::assemble_noflow_bc(const String& mesh_parts)
  {
    Assembly::UnitFilterAssembler<MeshType> unit_filter_asm;

    // get the filter name
    String filter_name = aux_add_mesh_parts(unit_filter_asm, this->domain_level, mesh_parts);

    // get the unit filter and clear it
    BaseClass::LocalVeloUnitFilter& filter = this->get_local_velo_unit_filter_seq().find_or_add(filter_name);

    // clear and assemble
    filter.clear();
    unit_filter_asm.assemble(filter, this->domain_level.space_velo);
  }

  void StokesLevel::assemble_slip_bc(const String& mesh_parts)
  {
    Assembly::SlipFilterAssembler<TrafoType> slip_filter_asm(this->domain_level.trafo);

    // get the filter name
    String filter_name = aux_add_mesh_parts(slip_filter_asm, this->domain_level, mesh_parts);

    // get the slip filter
    typename BaseClass::LocalVeloSlipFilter& filter = this->get_local_velo_slip_filter_seq().find_or_add(filter_name);

    // clear and assemble the filter
    filter.clear();
    slip_filter_asm.assemble(filter, this->domain_level.space_velo);

    // synchronize the slip filters
    this->sync_velocity_slip_filters();
  }

  void StokesLevel::assemble_pressure_mean_filter(const String& cubature)
  {
    // check if the filter has already been assembled
    if(this->get_local_pres_mean_filter().empty())
      BaseClass::assemble_pressure_mean_filter(this->domain_level.space_pres, cubature);
  }
} // namespace CCNDSimple
