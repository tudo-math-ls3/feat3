// includes, system
#include <iostream> // for std::ostream
#include <cstdlib> // for exit()

// includes, FEAST
#include <kernel/base_header.hpp>
// COMMENT_HILMAR: wenn Fusion mit MPI vollendet, dann nutze das hier:
//#include <kernel/error_handler.hpp>
#include <kernel/base_mesh/file_parser.hpp>
#include <kernel/base_mesh/bm.hpp>

using namespace FEAST;
using namespace BaseMesh;

#define WDIM 3
#define SDIM 3

int main (int argc, char **argv)
{
  BM<SDIM, WDIM> bm;
  FileParser<SDIM, WDIM> parser;
  try
  {
    parser.parse("dummy", &bm);
  }
  catch(InternalError* e)
  {
    std::cerr << e->message() << std::endl;
    exit(1);
  }
  // set cell numbers (currently equal to indices since all cells are active)
  bm.set_cell_numbers();
  // create base mesh's graph structure
  bm.create_graph();
  // print base mesh
  bm.print(std::cout);
  // validate base mesh
  bm.validate();

  // subdivide cell 0
  std::cout << "******************" << std::endl;
  std::cout << "Subdividing cell 0" << std::endl;
  std::cout << "******************" << std::endl;
  SubdivisionData<3, SDIM, WDIM>* subdiv_data = new SubdivisionData<3, SDIM, WDIM>(NONCONFORM_SAME_TYPE);
  bm.cell(0)->subdivide(subdiv_data);

  // add created cells and subcells to the corresponding base mesh vectors
  bm.add_created_items(bm.cell(0)->subdiv_data());
  // set cell numbers (now they differ from indices)
  bm.set_cell_numbers();
  // print base mesh
  bm.print(std::cout);
  // validate base mesh
  bm.validate();

  // TODO: neighbourhood update

  std::cout << "!!! Neighbourhood update after subdivision not implemented yet!!!" << std::endl;
  std::cout << "!!! Neighbourhood update after subdivision not implemented yet!!!" << std::endl;
  std::cout << "!!! DTORS not checked yet! Possible memory holes! Not 'valgrinded' yet !!!" << std::endl;
  std::cout << "!!! DTORS not checked yet! Possible memory holes! Not 'valgrinded' yet !!!" << std::endl;
}
