// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This header/source file pair defines the VtkWriter class, which is an extension of the simple Geometry::ExportVTK
// class template. The VtkWriter class allows to export of the solutions on a once refined mesh, which is generally
// useful since the spatial discretization uses a Q2/P1dc FE space pair, which contains more DOFs than could be
// exported into a VTK file on the same level.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "base.hpp"
#include "domain_control.hpp"
#include "stokes_level.hpp"
#include "time_stepping.hpp"

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

    // ----------------
    // input attributes
    // ----------------

    /// write refined VTK files?
    bool want_refined = true;

    /// VTK name prefix
    String name_prefix;

    /// VTK filename
    String vtk_name;

    /// VTK stepping
    Index stepping = Index(0);

    // -----------------
    // output attributes
    // -----------------

    /// plot line for short plot
    String plot_line;

    /// watch for VTK output
    StopWatch watch_vtk_write;

    // ----------------
    // state attributes
    // ----------------

    /// the once refined mesh (nullptr if we're exporting on non-refined mesh)
    std::unique_ptr<MeshType> refined_mesh;

    /// the actual VTK exporter that does the dirty work
    std::unique_ptr<Geometry::ExportVTK<MeshType>> exporter;

    /// a map of Stokes vectors
    std::map<String, GlobalStokesVector*> stokes_vectors;

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  public:
    /**
     * \brief Constructor
     *
     * \param[in] domain_
     * A \resident reference to the domain controller
     *
     * \param[in] refined_
     * Specifies whether to export on a once refined mesh
     */
    explicit VtkWriter(const DomainControl& domain_, bool refined_ = true);

    /// destructor
    virtual ~VtkWriter();

    /// adds all supported arguments to the argument parser
    static void add_supported_args(SimpleArgParser& args);

    /// parse arguments from command line
    virtual bool parse_args(SimpleArgParser& args);

    // print configuration to console
    virtual void print_config();

    /// adds a Stokes vector reference to the writer
    virtual bool register_stokes_vector(const String& name, GlobalStokesVector& vector);

    /// removes all previously registered Stokes vectors from the writer
    virtual void unregister_stokes_vectors();

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
    void add_p0dc_vector(const LocalScalarVectorType& vector, const String& name);
    void add_p1dc_vector(const LocalScalarVectorType& vector, const String& name);

    /**
     * \brief Writes the actual VTK file (for the current time step)
     *
     * This function also clears the exporter from all previously added vectors, so that the
     * same VtkWriter object can be used to write several VTKs, e.g. in an unsteady simulation.
     */
    virtual void write();

    /**
     * \brief Writes all registered vectors to a VTK file
     *
     * This function checks whether a VTK file is to be written for the current time step and,
     * if so, writes all vectors that have been registered to this writer to a VTK file.
     *
     * \returns
     * \c true, if a VTK file was written or \c false, if no VTK file is to be written in this
     * time step.
     */
    virtual bool write_registered(const TimeStepping& time_stepping);

    /// prints runtime information for this object
    virtual void print_runtime(double total_time);

  }; // class VtkWriter
} // namespace CCNDSimple
