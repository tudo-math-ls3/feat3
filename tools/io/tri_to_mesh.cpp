#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_file_writer.hpp>

#include <iostream>
#include <fstream>
#include <list>
#include <algorithm>
#include <cstdlib>

using namespace FEAT;
using CMESH = Geometry::ConformalMesh<Shape::Hexahedron>;
using RMESH = Geometry::RootMeshNode<CMESH>;

void add_mesh_part(RMESH * rmesh, String & file)
{
    std::ifstream par_file(file);
    if (!par_file)
    {
      std::cerr<<"file " << file << " not found!"<<std::endl;
    }

    String header;
    getline(par_file, header);
    auto words = header.split_by_charset(String::whitespaces());
    Index nodes = std::strtoul(words.at(0).c_str(), nullptr, 0);
    getline(par_file, header); // skip bulk line

    Index num_entities[4];
    num_entities[0] = nodes;
    num_entities[1] = 0;
    num_entities[2] = 0;
    num_entities[3] = 0;

    auto mesh_part = new typename RMESH::MeshPartType(num_entities);
    Index counter(0);
    for (String line ; getline(par_file, line); )
    {
      mesh_part->get_target_set<0>()[counter] = std::strtoul(line.c_str(), nullptr, 0) - 1;
      ++counter;
    }

    mesh_part->deduct_target_sets_from_bottom<0>(rmesh->get_mesh()->get_index_set_holder());
    String part_name = file.lower();
    part_name.replace_all(".par", "");
    rmesh->add_mesh_part(part_name, mesh_part);

    par_file.close();
}

int main(int argc, char ** argv)
{
    if (argc < 3)
    {
        std::cout<<"Usage 'otto2mesh file.prj mesh.xml [coord-scaling-factor]'"<<std::endl;
        std::cout<<"Must be executed in the file.prj directory, as no folder vodoo takes place internally."<<std::endl;
        exit(EXIT_FAILURE);
    }

    String input(argv[1]);
    String output(argv[2]);

    double scaling(1.);
    if (argc == 4)
    {
      scaling = std::strtod(argv[3], nullptr);
    }

    std::list<String> file_list;
    std::ifstream file_prj(input);
    if (!file_prj)
    {
      std::cerr<<"file " << input << " not found!"<<std::endl;
      exit(1);
    }
    for (String line ; getline(file_prj, line); )
    {
      file_list.push_back(line);
    }
    file_prj.close();

    auto it = file_list.begin();
    while (it != file_list.end() && !it->ends_with(".tri"))
    {
      ++it;
    }
    if (it == file_list.end())
    {
      std::cerr<<"no .tri file found in " << input << " !"<<std::endl;
      exit(1);
    }

    std::ifstream mesh_tri(*it);
    if (!mesh_tri)
    {
      std::cerr<<"file " << *it << " not found!"<<std::endl;
    }

    file_list.erase(it);

    Index num_entities[4];
    {
      String header;
      getline(mesh_tri, header); // skip bulk line
      getline(mesh_tri, header); // skip bulk line
      getline(mesh_tri, header);
      header.trim_me();
      auto words = header.split_by_charset(String::whitespaces());
      num_entities[0] = std::strtoul(words.at(1).c_str(), nullptr, 0);
      num_entities[1] = Index(0);
      num_entities[2] = Index(0);
      num_entities[3] = std::strtoul(words.at(0).c_str(), nullptr, 0);
    }

    auto cmesh = new CMESH(num_entities);

    Index counter(0);
    for (String line ; getline(mesh_tri, line); )
    {
      line.trim_me();
      if (line == "DCORVG")
        continue;
      if (line == "KVERT")
        break;
      auto words = line.split_by_charset(String::whitespaces());
      Tiny::Vector<CMESH::CoordType, 3> v;
      v[0] = std::strtod(words.at(0).c_str(), nullptr) * scaling;
      v[1] = std::strtod(words.at(1).c_str(), nullptr) * scaling;
      v[2] = std::strtod(words.at(2).c_str(), nullptr) * scaling;
      cmesh->get_vertex_set()[counter] = v;
      ++counter;
    }

    counter = 0;
    for (String line ; getline(mesh_tri, line); )
    {
      line.trim_me();
      if (line == "KNPR")
        break;

      auto words = line.split_by_charset(String::whitespaces());
      cmesh->get_index_set<3, 0>()(counter, 0) = std::strtoul(words.at(0).c_str(), nullptr, 0) - 1;
      cmesh->get_index_set<3, 0>()(counter, 1) = std::strtoul(words.at(1).c_str(), nullptr, 0) - 1;
      cmesh->get_index_set<3, 0>()(counter, 2) = std::strtoul(words.at(3).c_str(), nullptr, 0) - 1;
      cmesh->get_index_set<3, 0>()(counter, 3) = std::strtoul(words.at(2).c_str(), nullptr, 0) - 1;
      cmesh->get_index_set<3, 0>()(counter, 4) = std::strtoul(words.at(4).c_str(), nullptr, 0) - 1;
      cmesh->get_index_set<3, 0>()(counter, 5) = std::strtoul(words.at(5).c_str(), nullptr, 0) - 1;
      cmesh->get_index_set<3, 0>()(counter, 6) = std::strtoul(words.at(7).c_str(), nullptr, 0) - 1;
      cmesh->get_index_set<3, 0>()(counter, 7) = std::strtoul(words.at(6).c_str(), nullptr, 0) - 1;
      ++counter;
    }
    mesh_tri.close();

    cmesh->deduct_topology_from_top();
    auto rmesh = new Geometry::RootMeshNode<CMESH>(cmesh);

    for (auto file : file_list)
    {
      if (!file.ends_with(".par"))
      {
        std::cerr<<file<<" has unsupported type in file.prj!"<<std::endl;
        exit(1);
      }
      add_mesh_part(rmesh, file);
    }

    //Geometry::ExportVTK<CMESH> exporter(*(rmesh->get_mesh()));
    //exporter.write("test.vtk");
    std::ofstream mesh_out(output);
    Geometry::MeshFileWriter mesh_writer(mesh_out);
    mesh_writer.write(rmesh);

    delete rmesh;

    return 0;
}
