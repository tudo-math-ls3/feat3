// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef FEAT_APPLICATIONS_CCND_SIMPLE_VTK_WRITER_HPP
#define FEAT_APPLICATIONS_CCND_SIMPLE_VTK_WRITER_HPP 1

#include "domain.hpp"
#include "stokes_level.hpp"

#include <kernel/geometry/export_vtk.hpp>

namespace CCNDSimple
{
  /**
   * \brief VTK Writer class
   *
   * This class provides a slightly more sophisticated VTK writer, which enables to write a VTK on a
   * once refined mesh to exploit the higher precision of the Lagrange-2/P1dc element pair.
   *
   * \author Peter Zajac
   */
  class VtkWriter
  {
  public:
    /// the domain controller
    const DomainControl& domain;
    /// the once refined mesh (nullptr if we're exporting on non-refined mesh)
    std::unique_ptr<MeshType> refined_mesh;
    /// the actual VTK exporter that does the dirty work
    std::unique_ptr<Geometry::ExportVTK<MeshType>> exporter;

    // did the user want a refined export?
    bool want_refined = true;

    /// VTK name prefix
    String name_prefix;
    /// VTK filename
    String vtk_name;
    /// VTK stepping
    Index stepping = Index(0);

  public:
    /**
     * \brief Constructor
     *
     * \param[in] domain_
     * A \resident reference to the domain controller
     *
     * \param[in] name_prefix_
     * The default filename prefix to be used if the caller did not specify one via --vtk <filename>
     *
     * \param[in] refined_
     * Specifies whether to export on a once refined mesh
     */
    explicit VtkWriter(const DomainControl& domain_, const String& name_prefix_, bool refined_ = true);

    /// destructor
    virtual ~VtkWriter();

    /// adds all supported arguments to the argument parser
    static void add_supported_args(SimpleArgParser& args);

    /// parse arguments from command line
    virtual bool parse_args(SimpleArgParser& args);

    // print configuration to console
    virtual void print_config();

    /**
     * \brief Prepares the writer for the output in a steady simulation
     *
     * \returns
     * \c true, if a VTK file should be written, or \c false, if the caller did not want a VTK file to be written.
     */
    virtual bool prepare_write();

    /**
     * \brief Prepares the writer for the output of a time step in an unsteady simulation
     *
     * \returns
     * \c true, if a VTK file should be written in this time step, or \c false, if the caller did
     * not want a VTK file to be written in this time step (or at all).
     */
    virtual bool prepare_write(Index step);

    /**
     * \brief Adds a Stokes vector to the exporter
     *
     * \param[in] vector
     * A \transient reference to the Stokes vector to be exported
     *
     * \param[in] v_name, p_name
     * The names of the velocity and pressure components in the VTK file
     */
    void add_stokes_vector(const GlobalStokesVector& vector, const String& v_name = "v", const String& p_name = "p");

    /**
     * \brief Adds a scalar/blocked Lagrange-1/2 vector to the exported
     *
     * \param[in] vector
     * A \transient reference to the vector to be exported
     *
     * \param[in] name
     * The name of the vector in the VTK file
     */
    void add_lagrange1_vector(const LocalScalarVectorType& vector, const String& name);
    void add_lagrange2_vector(const LocalScalarVectorType& vector, const String& name);
    void add_lagrange1_vector(const LocalFieldVectorType& vector, const String& name);
    void add_lagrange2_vector(const LocalFieldVectorType& vector, const String& name);

    /**
    * \brief Adds a scalar P0dc/P1dc vector to the exported
    *
    * \param[in] vector
    * A \transient reference to the vector to be exported
    *
    * \param[in] name
    * The name of the vector in the VTK file
    */
    void add_p0_vector(const LocalScalarVectorType& vector, const String& name);
    void add_p1dc_vector(const LocalScalarVectorType& vector, const String& name);

    /**
     * \brief Writes the actual VTK file (for the current time step)
     *
     * This function also clears the exporter from all previously added vectors, so that the
     * same VtkWriter object can be used to write several VTKs, e.g. in an unsteady simulation.
     */
    void write();
  }; // class VtkWriter
} // namespace CCNDSimple

#endif // FEAT_APPLICATIONS_CCND_SIMPLE_VTK_WRITER_HPP
