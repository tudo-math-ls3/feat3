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

#define WDIM 2
#define SDIM 2
int main (int argc, char **argv)
{

  if (argc < 2)
  {
    std::cerr << "Call the program with \"" << argv[0] << " <relative_path_to_mesh_file>\"" << std::endl;
    return 1;
  }

  // COMMENT_HILMAR: FEAST1 file format, which we currently use, cannot distinguish space and world dimension. Hence,
  // the case WDIM = 3 and SDIM = 2 is currently not supported.
  assert(WDIM == SDIM);

  std::string name_mesh_file(argv[1]);
  BM<SDIM, WDIM> bm;
  FileParser<SDIM, WDIM> parser;
  try
  {
    parser.parse(name_mesh_file, &bm);
  }
  catch(InternalError* e)
  {
    std::cerr << e->message() << std::endl;
    exit(1);
  }
  // set cell numbers (equal to indices since all cells are active)
  bm.set_cell_numbers();
  // create base mesh's graph structure
  bm.create_graph();
  // print base mesh
  bm.print(std::cout);

  // subdivide cell 0
  std::cout << "Subdividing cell 0" << std::endl;
  bm.cell(0)->init_subdiv_data(NONCONFORM_SAME_TYPE);
  bm.cell(0)->subdivide();
  bm.add_created_items(bm.cell(0)->subdiv_data());
  // set cell numbers (now they differ from indices)
  bm.set_cell_numbers();
  // print base mesh
  bm.print(std::cout);

 // validate base mesh
  bm.validate();

  std::cout << "!!! Neighbourhood update after subdivision not implemented yet!!!" << std::endl;
  std::cout << "!!! Neighbourhood update after subdivision not implemented yet!!!" << std::endl;
  std::cout << "!!! DTORS not checked yet! Possible memory holes! Not 'valgrinded' yet !!!" << std::endl;
  std::cout << "!!! DTORS not checked yet! Possible memory holes! Not 'valgrinded' yet !!!" << std::endl;


// old test code
//  BaseMesh2D<WDIM> bm;
//  std::cout << "Base Mesh in ASCII art (v3 is at (0,0) and all edges are unit length):" << std::endl;
//  std::cout << "v0---e0---v1---e1---v2 \\         " << std::endl;  // stupid '//' is necessary to escape the space after /
//  std::cout << " |         |         |   \\       " << std::endl;
//  std::cout << "e2   c0   e3   c1   e4  c2 \e5    " << std::endl;
//  std::cout << " |         |         |       \\   " << std::endl;
//  std::cout << "v3---e6---v4---e7---v5----e8---v6 " << std::endl;
//  std::cout << "                 /   |         |  " << std::endl;
//  std::cout << "            e9 / c3 e10   c4  e11 " << std::endl;
//  std::cout << "            /        |         |  " << std::endl;
//  std::cout << "          v7---e12--v8---e13---v9 " << std::endl;
//
//  // subdivide cell 0
//  std::cout << "Subdividing cell 0" << std::endl;
//  bm.cell(0)->init_subdiv_data(NONCONFORM_SAME_TYPE);
//  bm.cell(0)->subdivide();
//  bm.add_created_items(bm.cell(0)->subdiv_data());
//  // set cell numbers (equal to indices since all cells are active)
//  bm.set_cell_numbers();
//
//  // TODO: neighbourhood update
//
//  // create base mesh's graph structure
//  bm.create_graph();
//
//  std::cout << "Cell 0 should now look like this:" << std::endl;
//  std::cout << "v00---e16---v11---e17---v01" << std::endl;
//  std::cout << " |           |           | " << std::endl;
//  std::cout << "e19   c7    e23   c8    e21" << std::endl;
//  std::cout << " |           |           | " << std::endl;
//  std::cout << "v12---e24---v14---e25---v13" << std::endl;
//  std::cout << " |           |           | " << std::endl;
//  std::cout << "e18   c5    e22   c6    e20" << std::endl;
//  std::cout << " |           |           | " << std::endl;
//  std::cout << "v03---e14---v10---e15---v04" << std::endl;
//  bm.print(std::cout);
}
