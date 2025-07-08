// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This header/source file pair defines the StokesLevel class, which is derived from an instance of the
// Control::StokesBlockedCombinedSystemLevel class template. This class contains all linear algebra containers which
// are required on each level of the multigrid hierarchy, i.e. all gates, muxers, matrices, filters and grid transfers.
// This class also provides a set of member functions for the assembly of matrices and filters.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "base.hpp"
#include "domain_control.hpp"

#include <control/stokes_blocked.hpp>

namespace CCNDSimple
{
  // helper function: add meshparts to assembler
  template<typename Asm_>
  String aux_add_mesh_parts(Asm_& assembler, DomainLevel& domain_level, const String& mesh_parts);

  /**
   * \brief (Navier-)Stokes System Level class
   *
   * This class contains all linear algebra structures defining the Navier-Stokes system on a single level.
   *
   * See the base class (and its base class) in <control/stokes_blocked.hpp> for more details
   *
   * \author Peter Zajac
   */
  class StokesLevel :
    public Control::StokesBlockedCombinedSystemLevel<dim, DataType, IndexType>
  {
  public:
    /// our base-class
    typedef Control::StokesBlockedCombinedSystemLevel<dim, DataType, IndexType> BaseClass;

    /// a reference to the corresponding domain level
    DomainLevel& domain_level;

    /// the filtered local system matrix for Vanka/UMFPACK
    typename BaseClass::LocalSystemMatrix local_matrix_sys;

    /// the local velocity mass matrix (used for the RHS assembly in unsteady Stokes)
    typename BaseClass::LocalMatrixBlockA local_velo_mass_matrix;

    /// the global Stokes vector type
    typedef BaseClass::GlobalSystemVector GlobalStokesVector;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  public:
    /**
     * \brief Constructor
     *
     * \param[in] domain_lvl
     * The domain level that this Stokes level corresponds to
     */
    explicit StokesLevel(DomainLevel& domain_lvl) :
      domain_level(domain_lvl)
    {
    }

    /**
     * \brief Assembles the Laplace matrix into the matrix block A
     *
     * \param[in] cubature
     * The cubature rule to use for assembly
     *
     * \param[in] nu
     * The viscosity parameter nu
     *
     * \param[in] deformation
     * Specifies whether to use the deformation tensor formulation or the gradient tensor formulation.
     *
     * \note
     * This function only updates the matrix_a member object, i.e. you have to call the compile_system_matrix()
     * function afterwards to update the system matrix object.
     */
    void assemble_velocity_laplace_matrix(const String& cubature, const DataType nu, bool deformation);

    // assembles the local_velo_mass_matrix object
    void assemble_velocity_mass_matrix(const String& cubature);

    /**
     * \brief Assembles the pressure gradient / velocity divergence matrices B and D
     *
     * \param[in] cubature
     * The cubature rule to use for assembly
     *
     * \note
     * This function only updates the matrix_b and matrix_d member objects, i.e. you have to call the
     * compile_system_matrix() function afterwards to update the system matrix object.
     */
    void assemble_grad_div_matrices(const String& cubature);

    /**
     * \brief Compiles the local system matrix
     *
     * This function convert the type-0 global system matrix to the type-1 local system matrix and applies the
     * local unit filters onto it, so that it can be used as a matrix for a local AmaVanka smoother or UMFPACK
     * solver.
     */
    void compile_local_matrix();

    /**
     * \brief Assembles a flow boundary condition for the velocity.
     *
     * This function will create a single unit filter in the filter sequence for the union of all
     * the mesh-parts given.
     *
     * \param[in] mesh_parts
     * The names of the mesh-parts on which the boundary conditions are to be assembled,
     * separated by whitespaces.
     *
     * \param[in] flow_function
     * The analytic function to be used for the flow, given as a string.
     */
    template<typename FlowFunc_>
    void assemble_flow_bc(const String& mesh_parts, const FlowFunc_& flow_function)
    {
      Assembly::UnitFilterAssembler<MeshType> unit_filter_asm;

      // get the filter name
      String filter_name = aux_add_mesh_parts(unit_filter_asm, this->domain_level, mesh_parts);

      // get the unit filter and clear it
      BaseClass::LocalVeloUnitFilter& filter = this->get_local_velo_unit_filter_seq().find_or_add(filter_name);

      // clear and assemble
      filter.clear();
      unit_filter_asm.assemble(filter, this->domain_level.space_velo, flow_function);
    }

    /**
     * \brief Assembles a no-flow boundary condition for the velocity.
     *
     * This function will create a single unit filter in the filter sequence for the union of all
     * the mesh-parts given.
     *
     * \param[in] mesh_parts
     * The names of the mesh-parts on which the boundary conditions are to be assembled,
     * separated by whitespaces.
     */
    void assemble_noflow_bc(const String& mesh_parts);

    /**
     * \brief Assembles a slip boundary condition for the velocity.
     *
     * This function will create a single slip filter in the filter sequence for the union of all
     * the mesh-parts given.
     *
     * \param[in] mesh_parts
     * The names of the mesh-parts on which the boundary conditions are to be assembled,
     * separated by whitespaces.
     */
    void assemble_slip_bc(const String& mesh_parts);

    /**
     * \brief Assembles an integral-mean filter for the pressure.
     *
     * This function assembles a filter to make the pressure definite in the case where the domain
     * has no outflow boundary, i.e. where the flow domain is closed.
     *
     * \param[in] cubature
     * The cubature rule to use for integration.
     */
    void assemble_pressure_mean_filter(const String& cubature);

    /**
     * \brief Updates the RHS vector by adding the terms that depend on previous time step solutions
     *
     * \param[inout] vec_rhs
     * The RHS vector that should be updated
     *
     * \param[in] theta
     * The theta factors from the time stepping scheme.
     *
     * \param[in] vec_sol_1, vec_sol_2
     * The solution vectors of the last two time steps.
     */
    void update_unsteady_rhs(GlobalStokesVector& vec_rhs, const GlobalStokesVector& vec_tmp)
    {
      local_velo_mass_matrix.apply(vec_rhs.local().first(), vec_tmp.local().first(), vec_rhs.local().first());
    }
  }; // class StokesLevel

  /// the global Stokes matrix type
  typedef StokesLevel::GlobalSystemMatrix GlobalStokesMatrix;

  /// the global Stokes vector type
  typedef StokesLevel::GlobalSystemVector GlobalStokesVector;

  /// the global Stokes filter type
  typedef StokesLevel::GlobalSystemFilter GlobalStokesFilter;

  /// the global velocity vector type
  typedef StokesLevel::GlobalVeloVector GlobalVeloVector;

  /// the local velocity vector type
  typedef StokesLevel::LocalVeloVector LocalVeloVector;

  /// the local matrix A block
  typedef StokesLevel::LocalMatrixBlockA LocalMatrixBlockA;

  // helper function: add meshparts to assembler
  template<typename Asm_>
  String aux_add_mesh_parts(Asm_& assembler, DomainLevel& domain_level, const String& mesh_parts)
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
        assembler.add_mesh_part(*mpart);
    }

    // join all part names
    return stringify_join(parts_in, "##");
  }
} // namespace CCNDSimple
