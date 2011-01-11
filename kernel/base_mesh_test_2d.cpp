// includes, system
#include <iostream> // for std::ostream
#include <cstdlib> // for exit()

// includes, FEAST
#include <kernel/base_header.hpp>
#include <kernel/base_mesh_2d.hpp>
#include <kernel/base_mesh_cell.hpp>

using namespace FEAST;
using namespace BaseMesh;

#define WDIM 2
#define SDIM 2

int main (int argc, char **argv)
{
  BaseMesh2D<WDIM> bm;
  std::cout << "Base Mesh in ASCII art (v3 is at (0,0) and all edges are unit length):" << std::endl;
  std::cout << "v0---e0---v1---e1---v2 \\         " << std::endl;  // stupid '//' is necessary to escape the space after /
  std::cout << " |         |         |   \\       " << std::endl;
  std::cout << "e2   c0   e3   c1   e4  c2 \e5    " << std::endl;
  std::cout << " |         |         |       \\   " << std::endl;
  std::cout << "v3---e6---v4---e7---v5----e8---v6 " << std::endl;
  std::cout << "                 /   |         |  " << std::endl;
  std::cout << "            e9 / c3 e10   c4  e11 " << std::endl;
  std::cout << "            /        |         |  " << std::endl;
  std::cout << "          v7---e12--v8---e13---v9 " << std::endl;
  bm.print(std::cout);

  // subdivide cell 0
  std::cout << "Subdividing cell 0" << std::endl;
  SubdivisionData<2, WDIM, SDIM> subdiv_data;
  bm.cell(0)->subdivide(subdiv_data);
  bm.add_created_items(subdiv_data);
  // TODO: neighbourhood update

  std::cout << "Cell 0 should now look like this:" << std::endl;
  std::cout << "v00---e16---v11---e17---v01" << std::endl;
  std::cout << " |           |           | " << std::endl;
  std::cout << "e19   c7    e23   c8    e21" << std::endl;
  std::cout << " |           |           | " << std::endl;
  std::cout << "v12---e24---v14---e25---v13" << std::endl;
  std::cout << " |           |           | " << std::endl;
  std::cout << "e18   c5    e22   c6    e20" << std::endl;
  std::cout << " |           |           | " << std::endl;
  std::cout << "v03---e14---v10---e15---v04" << std::endl;
  bm.print(std::cout);

  // validate base mesh
  bm.validate();

  std::cout << "!!! Neighbourhood update after subdivision not implemented yet!!!" << std::endl;
  std::cout << "!!! Neighbourhood update after subdivision not implemented yet!!!" << std::endl;
  std::cout << "!!! DTORS not checked yet! Possible memory holes! Not 'valgrinded' yet !!!" << std::endl;
  std::cout << "!!! DTORS not checked yet! Possible memory holes! Not 'valgrinded' yet !!!" << std::endl;
}
