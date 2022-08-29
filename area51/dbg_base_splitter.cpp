// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/runtime.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/scalar_basic.hpp>

namespace DbgBaseSplitter
{
  using namespace FEAT;

  void main(int argc, char* argv [])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    Control::Domain::add_supported_pdc_args(args);
    args.support("mesh");
    args.support("level");
    args.support("in");
    args.support("in-b");
    args.support("out");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
        comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

      comm.print(std::cerr, "Supported Options are:");
      comm.print(std::cerr, args.get_supported_help());

      // abort
      FEAT::Runtime::abort();
    }

    if(args.check("mesh") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--mesh <mesh-file>' is missing!");
      FEAT::Runtime::abort();
    }
    if(args.check("level") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--level <levels>' is missing!");
      FEAT::Runtime::abort();
    }

    // define our mesh type
    typedef Shape::Hypercube<2> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceType;

    // create our domain control
    typedef Control::Domain::SimpleDomainLevel<MeshType, TrafoType, SpaceType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

    domain.keep_base_levels();
    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);
    domain.create(args.query("mesh")->second);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // plot our levels
    comm.print("Desired Levels: " + domain.format_desired_levels());
    comm.print("Chosen  Levels: " + domain.format_chosen_levels());

    // define our arch types
    typedef Mem::Main MemType;
    typedef double DataType;
    typedef Index IndexType;

    // choose our desired analytical solution
    Analytic::Common::SineBubbleFunction<ShapeType::dimension> sol_func;

    // define our system level
    typedef Control::ScalarUnitFilterSystemLevel<MemType, DataType, IndexType> SystemLevelType;

    // get our assembled vector type
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;

    typedef typename SystemLevelType::SystemMirror MirrorType;

    typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, 2> VectorTypeBlocked;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType the_system_level;
    the_system_level.assemble_gate(domain.front());

    //__debugbreak();
    the_system_level.assemble_base_splitter(domain.front());

    // create new vector
    GlobalSystemVector vector(&the_system_level.gate_sys);
    Assembly::Interpolator::project(vector.local(), sol_func, the_domain_level.space);


    GlobalSystemVector vector2 = vector.clone(LAFEM::CloneMode::Layout);
    vector2.format();

    LAFEM::FileMode file_mode = LAFEM::FileMode::fm_binary;

    // write out a vector?
    if(args.check("out") > 0)
    {
      String filename;
      args.parse("out", filename);

      // join vector for base level
      comm.print("Writing output file '" + filename + "'");
      the_system_level.base_splitter_sys.join_write_out(vector, filename, file_mode);
    }

    // read in a vector?
    if(args.check("in") > 0)
    {
      String filename;
      args.parse("in", filename);

      comm.print("Reading input file '" + filename + "'");
      the_system_level.base_splitter_sys.split_read_from(vector2, filename, file_mode);

      // compute errors
      vector2.axpy(vector, vector2, -1.0);
      double ad = vector2.max_abs_element();
      comm.print("Max absolute difference: " + stringify_fp_sci(ad));
    }

    // read in a blocked vector?
    if(args.check("in-b") > 0)
    {
      //if(comm.rank() == 0) __debugbreak();

      VectorTypeBlocked vector_b(the_domain_level.space.get_num_dofs());

      // convert gate (just for kicks)
      Global::Gate<VectorTypeBlocked, MirrorType> gate_blocked;
      gate_blocked.convert(the_system_level.gate_sys, vector_b.clone());

      // convert splitter
      Global::Splitter<VectorTypeBlocked, MirrorType> splitter_blocked;
      splitter_blocked.convert(the_system_level.base_splitter_sys, vector_b.clone());

      String filename;
      args.parse("in-b", filename);

      comm.print("Reading input file '" + filename + "'");
      splitter_blocked.split_read_from(vector_b, filename, file_mode);
    }
  }
} // namespace DbgBaseSplitter

int main(int argc, char* argv [])
{
  FEAT::Runtime::initialize(argc, argv);
  //try
  {
    DbgBaseSplitter::main(argc, argv);
  }
  /*catch (const std::exception& exc)
  {
    std::cerr << "ERROR: unhandled exception: " << exc.what() << std::endl;
    FEAT::Runtime::abort();
  }
  catch (...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    FEAT::Runtime::abort();
  }*/
  return FEAT::Runtime::finalize();
}
