// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "domain_control.hpp"

namespace CCNDSimple
{
  DomainLevel::DomainLevel(int lvl_idx, std::unique_ptr<Geometry::RootMeshNode<MeshType>> node) :
    BaseClass(lvl_idx, std::move(node)),
    trafo(BaseClass::get_mesh()),
    space_velo(trafo),
    space_pres(trafo),
    space_q1(trafo),
    domain_asm(trafo)
  {
  }

  DomainLevel::~DomainLevel()
  {
  }

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  DomainControl::DomainControl(const Dist::Comm& comm_) :
    BaseClass(comm_, true)
  {
  }

  DomainControl::~DomainControl()
  {
  }

  void DomainControl::add_supported_args(SimpleArgParser& args)
  {
    //          1         2         3         4         5         6         7         8         9         0
    // 123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
    Control::Domain::add_supported_pdc_args(args);
    args.support("level", "<level-max> [[<levels-med...>] <level-min>]\n"
      "Specifies the refinement levels for the multigrid hierarchy.\n"
      "The mandatory first parameter <level-max> specifies the maximum refinement level on which the system\n"
      "is actually solved.\n"
      "The optional final parameter <level-min> specifies the minimum refinement level, which represents\n"
      "the coarse mesh level for the multigrid solver.\n"
      "The optional intermediary parameters <level-mid...> specify a sequence of partition change levels\n"
      "of the form 'level:ranks' where 'level' is the refinement level index and 'ranks' specifies the\n"
      "number of MPI ranks that should be used from this level to the next coarser partitioning change\n"
      "level or the coarse level, if no further partition change levels were given."
    );
    args.support("mesh", "<meshfiles...>\nSpecifies the meshfile(s) to be read in for the coarse mesh.");
  }

  bool DomainControl::parse_args(SimpleArgParser& args)
  {
    if(!BaseClass::parse_args(args))
      return false;

    return true;
  }

  void DomainControl::create_domain(SimpleArgParser& args)
  {
    watch_domain_create.start();

    if(args.check("mesh") < 1)
    {
      _comm.print(std::cerr, "ERROR: Mandatory option '--mesh <mesh-file>' is missing!");
      FEAT::Runtime::abort();
    }
    if(args.check("level") < 1)
    {
      _comm.print(std::cerr, "ERROR: Mandatory option '--level <levels>' is missing!");
      FEAT::Runtime::abort();
    }

    this->parse_args(args);
    this->set_desired_levels(args.query("level")->second);

    // if we want to write/read a joined solution, we also need to tell the domain control to keep the base levels
    if((args.check("save-joined-sol") >= 0) || (args.check("load-joined-sol") >= 0))
      this->keep_base_levels();

    // Our mesh file reader
    Geometry::MeshFileReader mesh_reader;

    // read in the mesh files
    mesh_file_names = args.query("mesh")->second;
    mesh_reader.add_mesh_files(_comm, mesh_file_names);

    // read the mesh file root markups
    mesh_reader.read_root_markup();
    String mesh_type = mesh_reader.get_meshtype_string();

    // run 2D or 3D ?
    if constexpr(dim == 2)
    {
      if(mesh_type != "conformal:hypercube:2:2")
      {
        _comm.print("ERROR: invalid mesh type '" + mesh_type + "', expected 'conformal:hypercube:2:2'");
        FEAT::Runtime::abort();
      }
    }
    if constexpr(dim == 3)
    {
      if(mesh_type != "conformal:hypercube:3:3")
      {
        _comm.print("ERROR: invalid mesh type '" + mesh_type + "', expected 'conformal:hypercube:3:3'");
        FEAT::Runtime::abort();
      }
    }

    // create the domain control
    BaseClass::create(mesh_reader);

    // compile all domain assemblers
    for(Index i(0); i < this->size_physical(); ++i)
      this->at(i)->domain_asm.compile_all_elements();

    // okay, we're done here
    watch_domain_create.stop();
  }

  void DomainControl::print_info()
  {
    print_pad(_comm, "Mesh Files", stringify_join(this->mesh_file_names, " "));
    auto parti_info = this->get_chosen_parti_info().split_by_charset("\n");
    for(auto& pi : parti_info)
      print_pad(_comm, "Partitioning Info", pi);
    print_pad(_comm, "Desired Levels", this->format_desired_levels());
    print_pad(_comm, "Chosen Levels", this->format_chosen_levels());
  }

  void DomainControl::print_runtime(double total_time)
  {
    print_time(_comm, "Domain Create Time", watch_domain_create.elapsed(), total_time);
  }
} // namespace CCNDSimple
