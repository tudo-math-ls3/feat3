// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//
// This header/source file pair defines the DomainLevel class, which is a required class for the use of the
// Control::Domain::PartiDomainControl class, as well as a custom DomainControl class, which derives from the
// PartiDomainControl class. The DomainLevel class merely defines the basic mandatory contents of a domain level, i.e.
// the trafo and the FE spaces a well as a DomainAssembler object. The DomainControl class extends the
// PartiDomainControl with some additional output. The DomainControl class is used by the application to read in the
// mesh from a mesh file and to partition and distribute it to the individual processes in an MPI-parallel simulation.
//
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#pragma once

#include "base.hpp"

#include <kernel/assembly/domain_assembler.hpp>
#include <control/domain/parti_domain_control.hpp>

namespace CCNDSimple
{
  /**
   * \brief Domain level class for CCND apps
   *
   * This class stores the mesh node, the trafo, the FE spaces as well as the domain assembler for
   * a single refinement level.
   *
   * See the base class in <control/domain/domain_level.hpp> for details.
   *
   * \author Peter Zajac
   */
  class DomainLevel :
    public Control::Domain::DomainLevel<MeshType>
  {
  public:
    /// our base class type
    typedef Control::Domain::DomainLevel<MeshType> BaseClass;

    /// our trafo
    TrafoType trafo;

    /// our velocity space: Lagrange-2
    SpaceVeloType space_velo;

    /// our pressure space: discontinuous P1
    SpacePresType space_pres;

    /// auxiliary Q1 space
    SpaceTypeQ1 space_q1;

    /// our domain assembler object
    Assembly::DomainAssembler<TrafoType> domain_asm;

  public:
    /**
     * \brief Constructor
     *
     * \param[in] lvl_idx
     * The level index of this domain level.
     *
     * \param[in] node
     * The mesh node for this domain level
     */
    explicit DomainLevel(int lvl_idx, std::unique_ptr<Geometry::RootMeshNode<MeshType>> node);

    /// destructor
    virtual ~DomainLevel();
  }; // class DomainLevel<...>

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * \brief Domain controller class for simple CCND apps
   *
   * This class creates and manages the distributed domain and mesh hierarchy
   *
   * See the base class (and its base class) in <control/domain/parti_domain_control.hpp> for details.
   *
   * \author Peter Zajac
   */
  class DomainControl :
    public Control::Domain::PartiDomainControl<DomainLevel>
  {
  public:
    /// our base class
    typedef Control::Domain::PartiDomainControl<DomainLevel> BaseClass;

    /// the names of the mesh files
    std::deque<String> mesh_file_names;

    /// watch for domain creation
    StopWatch watch_domain_create;

  public:
    /**
     * \brief Constructor
     *
     * \param[in] comm_
     * The communicator for this domain control
     */
    explicit DomainControl(const Dist::Comm& comm_);

    /// destructor
    virtual ~DomainControl();

    /// adds all supported arguments to the argument parser
    static void add_supported_args(SimpleArgParser& args);

    /// parses the arguments
    virtual bool parse_args(SimpleArgParser& args) override;

    /// creates the domain controller based on the arguments
    virtual void create_domain(SimpleArgParser& args);

    /// prints info on chosen partitioning and levels after creation
    virtual void print_info();

    /// prints runtime information for this object
    virtual void print_runtime(double total_time);

  }; // class DomainControl

  /// get a typedef for our domain layer type
  typedef typename DomainControl::LayerType DomainLayer;
} // namespace CCNDSimple
