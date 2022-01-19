// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

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
    Index nodes(0);
    words.at(0).parse(nodes);
    getline(par_file, header); // skip bulk line

    Index num_entities[4];
    num_entities[0] = nodes;
    num_entities[1] = 0;
    num_entities[2] = 0;
    num_entities[3] = 0;

    std::unique_ptr<typename RMESH::MeshPartType> mesh_part(new typename RMESH::MeshPartType(num_entities));
    Index counter(0);
    for (String line ; getline(par_file, line); )
    {
      line.parse(mesh_part->get_target_set<0>()[counter]);
      mesh_part->get_target_set<0>()[counter] -= 1;

      ++counter;
    }

    mesh_part->deduct_target_sets_from_bottom<0>(rmesh->get_mesh()->get_index_set_holder());
    String part_name = file.lower();
    part_name.replace_all(".par", "");
    rmesh->add_mesh_part(part_name, std::move(mesh_part));

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
      String sscaling(argv[3]);
      sscaling.parse(scaling);
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
      words.at(1).parse(num_entities[0]);
      num_entities[1] = Index(0);
      num_entities[2] = Index(0);
      words.at(0).parse(num_entities[3]);
    }

    std::unique_ptr<CMESH> cmesh(new CMESH(num_entities));

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
      words.at(0).parse(v[0]);
      v[0] *= scaling;
      words.at(1).parse(v[1]);
      v[1] *= scaling;
      words.at(2).parse(v[2]);
      v[2] *= scaling;
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
      words.at(0).parse(cmesh->get_index_set<3, 0>()(counter, 0));
      cmesh->get_index_set<3, 0>()(counter, 0) -= 1;
      words.at(1).parse(cmesh->get_index_set<3, 0>()(counter, 1));
      cmesh->get_index_set<3, 0>()(counter, 1) -= 1;
      words.at(3).parse(cmesh->get_index_set<3, 0>()(counter, 2));
      cmesh->get_index_set<3, 0>()(counter, 2) -= 1;
      words.at(2).parse(cmesh->get_index_set<3, 0>()(counter, 3));
      cmesh->get_index_set<3, 0>()(counter, 3) -= 1;
      words.at(4).parse(cmesh->get_index_set<3, 0>()(counter, 4));
      cmesh->get_index_set<3, 0>()(counter, 4) -= 1;
      words.at(5).parse(cmesh->get_index_set<3, 0>()(counter, 5));
      cmesh->get_index_set<3, 0>()(counter, 5) -= 1;
      words.at(7).parse(cmesh->get_index_set<3, 0>()(counter, 6));
      cmesh->get_index_set<3, 0>()(counter, 6) -= 1;
      words.at(6).parse(cmesh->get_index_set<3, 0>()(counter, 7));
      cmesh->get_index_set<3, 0>()(counter, 7) -= 1;
      ++counter;
    }
    mesh_tri.close();

    cmesh->deduct_topology_from_top();
    std::unique_ptr<Geometry::RootMeshNode<CMESH>> rmesh(new Geometry::RootMeshNode<CMESH>(std::move(cmesh)));

    for (auto file : file_list)
    {
      if (!file.ends_with(".par"))
      {
        std::cerr<<file<<" has unsupported type in file.prj!"<<std::endl;
        exit(1);
      }
      add_mesh_part(rmesh.get(), file);
    }

    //Geometry::ExportVTK<CMESH> exporter(*(rmesh->get_mesh()));
    //exporter.write("test.vtk");
    std::ofstream mesh_out(output);
    Geometry::MeshFileWriter mesh_writer(mesh_out);
    mesh_writer.write(rmesh.get());

    return 0;
}
