// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef FEAT_APPLICATIONS_CCND_SIMPLE_STOKES_LEVEL_HPP
#define FEAT_APPLICATIONS_CCND_SIMPLE_STOKES_LEVEL_HPP 1

#include "domain.hpp"

#include <control/stokes_blocked.hpp>

namespace CCNDSimple
{
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
     * \brief Assembles an inflow boundary condition for the velocity.
     *
     * This function will create a single unit filter in the filter sequence for the union of all
     * the mesh-parts given.
     *
     * \param[in] mesh_parts
     * The names of the mesh-parts on which the boundary conditions are to be assembled,
     * separated by whitespaces.
     *
     * \param[in] formula
     * The formula to be used for the inflow, given as a string.
     */
    void assemble_inflow_bc(const String& mesh_parts, const String& formula);

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
} // namespace CCNDSimple

#endif // FEAT_APPLICATIONS_CCND_SIMPLE_STOKES_LEVEL_HPP
