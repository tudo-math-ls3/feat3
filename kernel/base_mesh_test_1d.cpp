// includes, system
#include <iostream> // for std::ostream
#include <cstdlib> // for exit()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh_1d.hpp>
#include <kernel/base_mesh_cell.hpp>

using namespace FEAST;
using namespace BaseMesh;

#define WDIM 1
#define SDIM 1

int main (int argc, char **argv)
{
  BaseMesh1D<WDIM> bm;
  std::cout << "Base mesh in ASCII art (v0 is at (0,0) and all edges are unit length):" << std::endl;
  std::cout << "v0--------e0--------v1--------e1--------v2--------e2--------v3" << std::endl;
  bm.print(std::cout);

  // create a couple of actions and apply them
  std::cout << "Subdividing cell 0" << std::endl;

  SubdivisionData<1, WDIM, SDIM> subdiv_data;
  bm.cell(1)->subdivide(subdiv_data);
  bm.add_created_items(subdiv_data);
  // TODO: neighbourhood update

  std::cout << "Base mesh should now look like this:" << std::endl;
  std::cout << "v0--------e0--------v1---e3---v4---e4---v2--------e2--------v3" << std::endl;
  bm.print(std::cout);

  // validate base mesh
  bm.validate();

  std::cout << "!!! Neighbourhood update after subdivision not implemented yet!!!" << std::endl;
  std::cout << "!!! Neighbourhood update after subdivision not implemented yet!!!" << std::endl;
  std::cout << "!!! DTORS not checked yet! Possible memory holes! Not 'valgrinded' yet !!!" << std::endl;
  std::cout << "!!! DTORS not checked yet! Possible memory holes! Not 'valgrinded' yet !!!" << std::endl;

}
