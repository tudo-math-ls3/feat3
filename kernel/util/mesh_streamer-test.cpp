#include <kernel/util/mesh_streamer.hpp>
#include <test_system/test_system.hpp>
#include <sstream>

using namespace FEAST;
using namespace FEAST::TestSystem;

/**
 * \brief Test class for the MeshStreamer class.
 *
 * \test Tests the MeshStreamer class.
 *
 * \author Constantin Christof
 */

class MeshStreamerTest
  : public TaggedTest<Archs::None, Archs::None>
{
public:
  MeshStreamerTest() :
    TaggedTest<Archs::None, Archs::None>("mesh_streamer_test")
  {
  }

  void test_0() const
  {
    using namespace std;
    stringstream ioss;

    // let's write an (awfully ugly) mesh file into the stream
    ioss << "<feast_mesh_file>" << endl;
    ioss << "<header>  " << endl;
    ioss << " version 1" << endl;
    ioss << " chart_file unit_quad_chart.txt" << endl;
    ioss << " meshparts 2" << endl;
    ioss << "</header>" << endl;
    ioss << "   <info>    " << endl;
    ioss << "This file contains a simple unit-square..." << endl;
    ioss << "</info>" << endl;
    ioss << "<mesh>" << endl;
    ioss << " <header>" << endl;
    ioss << "  type conformal" << endl;
    ioss << "  shape quad" << endl;
    ioss << "  coords 2" << endl;

    ioss << " </header>" << endl;
    ioss << " <info>" << endl;
    ioss << "test" << endl;
    ioss << " </info>" << endl;
    ioss << " <counts>" << endl;
    ioss << "  verts 4" << endl;
    ioss << "  quads 1" << endl;
    ioss << "  edges 4" << endl;
    ioss << " </counts>" << endl;
    ioss << " <coords>" << endl;
    ioss << " 0.0 0.0" << endl;
    ioss << "  1.0 0.0" << endl;
    ioss << "0.0 1.0" << endl;
    ioss << "  1.0 1.0       " << endl;
    ioss << " </coords>" << endl;
    ioss << " <vert@edge>" << endl;
    ioss << "  0 1" << endl;
    ioss << "2   3" << endl;
    ioss << "  0   2   " << endl;
    ioss << "1                    3" << endl;
    ioss << " </vert@edge>" << endl;
    ioss << " <vert@quad>" << endl;
    ioss << "  0 1 2 3" << endl;
    ioss << " </vert@quad>" << endl;
    ioss << "</mesh>" << endl;
    // cellset
    ioss << "<meshpart>" << endl;
    ioss << " <header>" << endl;
    ioss << "  name cellset_1" << endl;
    ioss << "  parent root" << endl;
    ioss << "  type conformal" << endl;
    ioss << "  shape vertex" << endl;
    ioss << " </header>" << endl;
    ioss << " <info>" << endl;
    ioss << " I am a cellset!" << endl;
    ioss << " </info>" << endl;
    ioss << " <counts>" << endl;
    ioss << "  verts 4" << endl;
    ioss << " </counts>" << endl;
    ioss << " <vert_idx>" << endl;
    ioss << "1" << endl;
    ioss << "2" << endl;
    ioss << "3" << endl;
    ioss << "0" << endl;
    ioss << " </vert_idx>" << endl;
    ioss << "</meshpart>" << endl;
    // submesh
    ioss << "<meshpart>   " << endl;
    ioss << " <header>" << endl;
    ioss << "  name outer" << endl;
    ioss << "  parent root" << endl;
    ioss << "  type conformal  " << endl;
    ioss << "  shape edge " << endl;
    ioss << "  attribute_sets 1" << endl;
    ioss << "</header>" << endl;
    ioss << " <info>  " << endl;
    ioss << " This is a submesh that..." << endl;
    ioss << " </info> " << endl;
    ioss << " <counts>" << endl;
    ioss << "   verts    5   " << endl;
    ioss << "  edges 4" << endl;
    ioss << "  </counts>" << endl;
    ioss << " <attribute> " << endl;
    ioss << "  <header> " << endl;
    ioss << "   dimension 0" << endl;
    ioss << "   name Parametrisation" << endl;
    ioss << "   value_dim 1" << endl;
    ioss << "   value_count 5" << endl;
    ioss << "  </header> " << endl;
    ioss << "   <values> " << endl;
    ioss << "    0.0" << endl;
    ioss << "    1.0" << endl;
    ioss << "    2.0" << endl;
    ioss << "3.0" << endl;
    ioss << "    4.0 " << endl;
    ioss << " </values> " << endl;
    ioss << " </attribute> " << endl;
    ioss << " <vert@edge>" << endl;
    ioss << " 0  1" << endl;
    ioss << "1  2   " << endl;
    ioss << "  2 3" << endl;
    ioss << "3 4" << endl;
    ioss << "</vert@edge>" << endl;
    ioss << " <vert_idx>" << endl;
    ioss << "0" << endl;
    ioss << "1" << endl;
    ioss << "2" << endl;
    ioss << "3" << endl;
    ioss << "0" << endl;
    ioss << " </vert_idx>" << endl;
    ioss << " <edge_idx>" << endl;
    ioss << "0" << endl;
    ioss << "3" << endl;
    ioss << "1" << endl;
    ioss << "2" << endl;
    ioss << " </edge_idx>" << endl;
    ioss << "</meshpart>" << endl;
    ioss << "</feast_mesh_file>";

    // two mesh reader objects for reading and writing the mesh data
    MeshStreamer reader, writer;

    // parse the stream
    writer.parse_mesh_file(ioss);

    // temporary MeshNodes to test the _insert_sub_mesh function
    MeshStreamer::MeshNode* tempMN1 = new MeshStreamer::MeshNode();
    MeshStreamer::MeshDataContainer& tempMDC1 = tempMN1->mesh_data;
    tempMDC1.name = "tempSubmesh1";
    tempMDC1.parent = "root";
    tempMDC1.info = "super info1";

    MeshStreamer::MeshNode* tempMN2 = new MeshStreamer::MeshNode();
    MeshStreamer::MeshDataContainer& tempMDC2 = tempMN2->mesh_data;
    tempMDC2.name = "tempSubmesh2";
    tempMDC2.parent = "tempSubmesh1";
    tempMDC2.info = "super info2";

    // insert the temporary MeshNodes
    writer._insert_sub_mesh(tempMN1);
    TEST_CHECK_EQUAL(writer.get_num_submeshes(), Index(3));
    TEST_CHECK_EQUAL((writer.get_root_mesh_node())->find_sub_mesh("tempSubmesh1"), tempMN1);
    TEST_CHECK_EQUAL((writer.get_root_mesh_node())->find_sub_mesh("tempSubmesh1")->mesh_data.name, "tempSubmesh1");
    TEST_CHECK_EQUAL((writer.get_root_mesh_node())->find_sub_mesh("tempSubmesh1")->mesh_data.info, "super info1");

    writer._insert_sub_mesh(tempMN2);
    TEST_CHECK_EQUAL(writer.get_num_submeshes(), Index(4));
    TEST_CHECK_EQUAL((writer.get_root_mesh_node())->find_sub_mesh("tempSubmesh2"), tempMN2);
    TEST_CHECK_EQUAL((writer.get_root_mesh_node())->find_sub_mesh("tempSubmesh2")->mesh_data.name, "tempSubmesh2");
    TEST_CHECK_EQUAL((writer.get_root_mesh_node())->find_sub_mesh("tempSubmesh2")->mesh_data.info, "super info2");

    // remove the temporary submeshes
    writer._delete_sub_mesh("tempSubmesh1");

    // drop the data into an auxiliary file
    stringstream ioss2;
    writer.write_mesh_file(ioss2);

    // parse the data with the other mesh reader
    reader.parse_mesh_file(ioss2);

    // now check if everything is ok

    // check members
    TEST_CHECK_EQUAL(reader.get_chart_path(), "unit_quad_chart.txt");
    TEST_CHECK_EQUAL(reader.get_num_submeshes(), Index(2));
    TEST_CHECK_EQUAL(reader.get_info(), "This file contains a simple unit-square..." );

    // test the get_mesh function
    MeshStreamer::MeshDataContainer* root_mesh_ptr = reader.get_mesh("rooot");
    TEST_CHECK_EQUAL(root_mesh_ptr, ((MeshStreamer::MeshDataContainer*)nullptr));

    root_mesh_ptr = reader.get_mesh("root");
    TEST_CHECK_NOT_EQUAL(root_mesh_ptr, ((MeshStreamer::MeshDataContainer*)nullptr));

    MeshStreamer::MeshDataContainer& root_mesh(*root_mesh_ptr);

    // check the root mesh data
    TEST_CHECK_EQUAL(root_mesh.name , "root");
    TEST_CHECK_EQUAL(root_mesh.parent , "none");
    TEST_CHECK_EQUAL(root_mesh.chart , "none");
    TEST_CHECK_EQUAL(root_mesh.convert_mesh_type(root_mesh.mesh_type) , "conformal");
    TEST_CHECK_EQUAL(root_mesh.convert_shape_type(root_mesh.shape_type) , "quad");
    TEST_CHECK_EQUAL(root_mesh.coord_per_vertex , Index(2));
    TEST_CHECK_EQUAL(root_mesh.vertex_count , Index(4));
    TEST_CHECK_EQUAL(root_mesh.edge_count , Index(4));
    TEST_CHECK_EQUAL(root_mesh.tria_count , Index(0));
    TEST_CHECK_EQUAL(root_mesh.quad_count , Index(1));
    TEST_CHECK_EQUAL(root_mesh.tetra_count , Index(0));
    TEST_CHECK_EQUAL(root_mesh.hexa_count, Index(0));

    TEST_CHECK_EQUAL(root_mesh.info, "test");

    // check the root mesh coordinates
    std::vector<std::vector<double> >& coords = root_mesh.coords;

    // reference coordinates
    double coords_ref[] =
    {
      0.0, 0.0,
      1.0, 0.0,
      0.0, 1.0,
      1.0, 1.0
    };

    // loop through the coordinates
    Index count(0);
    bool error = false;
    for(Index i(0); i < coords.size(); ++i)
    {
      for(Index j(0); j < root_mesh.coord_per_vertex; ++j)
      {
        if((coords.at(i)).at(j) != coords_ref[count])
        {
          error = true;
          break;
        }
        ++count;
      }
      if(error)
      {
        break;
      }
    }

    // check if an error occured
    TEST_CHECK_EQUAL(error, false);

    // check the adjacencies
    std::vector<std::vector<Index> > adj_stack;

    // reference adjacencies
    Index adj_ref_01[] =
    {
      0, 1,
      2, 3,
      0, 2,
      1, 3
    };

    Index adj_ref_02[] =
    {
      0, 1, 2, 3
    };

    // check vertex at edge adjacencies
    adj_stack = root_mesh.adjacencies[0][1];
    count = 0;
    for(Index i(0); i < 4; ++i)
    {
      for(Index j(0); j < 2; ++j)
      {
        if((adj_stack.at(i)).at(j) != adj_ref_01[count])
        {
          error = true;
          break;
        }
        ++count;
      }
      if(error)
      {
        break;
      }
    }
    TEST_CHECK_EQUAL(error, false);

    // check vertex at quad adjacencies
    adj_stack = root_mesh.adjacencies[0][2];
    count = 0;
    for(Index j(0); j < 4; ++j)
    {
      if((adj_stack.at(0)).at(j) != adj_ref_02[count])
      {
        error = true;
        break;
      }
      ++count;
    }

    // check if an error occured
    TEST_CHECK_EQUAL(error, false);

    // check if the rest is emtpy
    TEST_CHECK_EQUAL((root_mesh.adjacencies[0][3]).empty(), true);

    // check parent indices
    TEST_CHECK_EQUAL((root_mesh.parent_indices[0]).empty(), true);
    TEST_CHECK_EQUAL((root_mesh.parent_indices[1]).empty(), true);
    TEST_CHECK_EQUAL((root_mesh.parent_indices[2]).empty(), true);
    TEST_CHECK_EQUAL((root_mesh.parent_indices[3]).empty(), true);

    //
    // check the cellset
    //

    String cellset_name = "cellset_1";
    MeshStreamer::MeshDataContainer* cellset = reader.get_mesh(cellset_name);

    TEST_CHECK_EQUAL(cellset->name, "cellset_1");
    TEST_CHECK_EQUAL(cellset->parent, "root");
    TEST_CHECK_EQUAL(cellset->vertex_count, Index(4));
    TEST_CHECK_EQUAL(cellset->edge_count, Index(0));
    TEST_CHECK_EQUAL(cellset->quad_count, Index(0));
    TEST_CHECK_EQUAL(cellset->tria_count, Index(0));
    TEST_CHECK_EQUAL(cellset->hexa_count, Index(0));
    TEST_CHECK_EQUAL(cellset->tetra_count, Index(0));
    TEST_CHECK_EQUAL(cellset->info, "I am a cellset!");

    // check parent indices of the cellset
    std::vector<Index> par_vert = cellset->parent_indices[0];

    // reference
    Index par_vert_ref[] =
    {
      1, 2, 3, 0
    };

    // check vertex parents
    error = false;
    for(Index j(0); j < 4; ++j)
    {
      if(par_vert_ref[j] != par_vert.at(j))
      {
        error = true;
        break;
      }
    }
    TEST_CHECK_EQUAL(error, false);

    // check if the rest is empty
    TEST_CHECK_EQUAL((cellset->parent_indices[1]).empty(), true);
    TEST_CHECK_EQUAL((cellset->parent_indices[2]).empty(), true);
    TEST_CHECK_EQUAL((cellset->parent_indices[3]).empty(), true);

    //
    // check the sub mesh data
    //

    MeshStreamer::MeshDataContainer* sub_mesh_ptr = reader.get_mesh("outer");
    TEST_CHECK_NOT_EQUAL(sub_mesh_ptr, ((MeshStreamer::MeshDataContainer*)nullptr));

    MeshStreamer::MeshDataContainer& sub_mesh(*sub_mesh_ptr);
    TEST_CHECK_EQUAL(sub_mesh.name , "outer");
    TEST_CHECK_EQUAL(sub_mesh.parent , "root");
    TEST_CHECK_EQUAL(sub_mesh.chart , "");
    TEST_CHECK_EQUAL(sub_mesh.convert_mesh_type(sub_mesh.mesh_type) , "conformal");
    TEST_CHECK_EQUAL(sub_mesh.convert_shape_type(sub_mesh.shape_type) , "edge");
    TEST_CHECK_EQUAL(sub_mesh.vertex_count , Index(5));
    TEST_CHECK_EQUAL(sub_mesh.edge_count , Index(4));
    TEST_CHECK_EQUAL(sub_mesh.tria_count , Index(0));
    TEST_CHECK_EQUAL(sub_mesh.quad_count , Index(0));
    TEST_CHECK_EQUAL(sub_mesh.tetra_count , Index(0));
    TEST_CHECK_EQUAL(sub_mesh.hexa_count, Index(0));
    TEST_CHECK_EQUAL(sub_mesh.attributes[0][0].value_dim, Index(1));

    TEST_CHECK_EQUAL(sub_mesh.info, "This is a submesh that...");

    // check the sub mesh coordinates
    std::vector<std::vector<double> >& coords_sub = sub_mesh.coords;

    // reference coordinates
    double coords_sub_ref[] =
    {
      0.0, 1.0, 2.0, 3.0, 4.0
    };

    // loop through the coordinates
    count = 0;
    error = false;
    for(Index i(0); i < 5; ++i)
    {
      for(Index j(0); j < sub_mesh.coord_per_vertex; ++j)
      {
        if((coords_sub.at(i)).at(j) != coords_sub_ref[count])
        {
          error = true;
          break;
        }
        ++count;
      }
      if(error)
      {
        break;
      }
    }

    // check if an error occured
    TEST_CHECK_EQUAL(error, false);

    // check the adjacencies
    std::vector<std::vector<Index> > adj_stack_sub;

    // reference adjacencies
    Index adj_sub_ref_01[] =
    {
      0, 1,
      1, 2,
      2, 3,
      3, 4
    };

    // check vertex at edge adjacencies
    adj_stack_sub = sub_mesh.adjacencies[0][1];
    count = 0;
    for(Index i(0); i < 4; ++i)
    {
      for(Index j(0); j < 2; ++j)
      {
        if((adj_stack_sub.at(i)).at(j) != adj_sub_ref_01[count])
        {
          error = true;
          break;
        }
        ++count;
      }
      if(error)
      {
        break;
      }
    }
    TEST_CHECK_EQUAL(error, false);

    // check if the rest is empty
    TEST_CHECK_EQUAL((sub_mesh.adjacencies[0][2]).empty(), true);
    TEST_CHECK_EQUAL((sub_mesh.adjacencies[0][3]).empty(), true);

    // check parent indices of the sub mesh
    std::vector<Index> par0 = sub_mesh.parent_indices[0];
    std::vector<Index> par1 = sub_mesh.parent_indices[1];

    // reference
    Index par0_ref[] =
    {
      0, 1, 2, 3, 0
    };

    Index par1_ref[] =
    {
      0, 3, 1, 2
    };

    // check vertex parents
    error = false;
    for(Index j(0); j < 5; ++j)
    {
      if(par0_ref[j] != par0.at(j))
      {
        error = true;
        break;
      }
    }
    TEST_CHECK_EQUAL(error, false);

    // check edge parents
    for(Index j(0); j < 4; ++j)
    {
      if(par1_ref[j] != par1.at(j))
      {
        error = true;
        break;
      }
    }
    TEST_CHECK_EQUAL(error, false);

    // check if the rest is empty
    TEST_CHECK_EQUAL((root_mesh.parent_indices[2]).empty(), true);
    TEST_CHECK_EQUAL((root_mesh.parent_indices[3]).empty(), true);

    // ok, everything is right
  } // test_0

  void test_1() const
  {
    using namespace std;
    stringstream ioss1, ioss2, ioss3, ioss4;

    // ioss1: missing "</mesh>" flag
    ioss1 << "<feast_mesh_file>" << endl;
    ioss1 << "<header>  " << endl;
    ioss1 << "  version 1" << endl;
    ioss1 << "  meshparts 0" << endl;
    ioss1 << "</header>" << endl;
    ioss1 << "<mesh>" << endl;
    ioss1 << "  <header>" << endl;
    ioss1 << "    type conformal" << endl;
    ioss1 << "    coords 2" << endl;
    ioss1 << "    shape quad" << endl;
    ioss1 << "  </header>" << endl;
    ioss1 << "  <counts>" << endl;
    ioss1 << "    verts 4" << endl;
    ioss1 << "    quads 1" << endl;
    ioss1 << "    edges 4" << endl;
    ioss1 << "  </counts>" << endl;
    ioss1 << "  <coords>" << endl;
    ioss1 << "    0.0 0.0" << endl;
    ioss1 << "    1.0 0.0" << endl;
    ioss1 << "    0.0 1.0" << endl;
    ioss1 << "    1.0 1.0" << endl;
    ioss1 << "  </coords>" << endl;
    ioss1 << "</feast_mesh_file>" << endl;

    // parse the stream and check if an exception is thrown
    MeshStreamer reader1;
    TEST_CHECK_THROWS(reader1.parse_mesh_file(ioss1), SyntaxError);

    // ioss2: wrong number of coordinates
    ioss2 << "<feast_mesh_file>" << endl;
    ioss2 << "<header>  " << endl;
    ioss2 << "  version 1" << endl;
    ioss2 << "  meshparts 0" << endl;
    ioss2 << "</header>" << endl;
    ioss2 << "<mesh>" << endl;
    ioss2 << "  <header>" << endl;
    ioss2 << "    type conformal" << endl;
    ioss2 << "    coords 2" << endl;
    ioss2 << "    shape quad" << endl;
    ioss2 << "  </header>" << endl;
    ioss2 << "  <counts>" << endl;
    ioss2 << "    verts 4" << endl;
    ioss2 << "    quads 1" << endl;
    ioss2 << "    edges 4" << endl;
    ioss2 << "  </counts>" << endl;
    ioss2 << "  <coords>" << endl;
    ioss2 << "    0.0 0.0" << endl;
    ioss2 << "    1.0 0.0" << endl;
    ioss2 << "    0.0 1.0" << endl;
    ioss2 << "    1.0 1.0 42.23" << endl;
    ioss2 << "  </coords>" << endl;
    ioss2 << "</mesh>" << endl;
    ioss2 << "</feast_mesh_file>" << endl;

    // parse the stream and check if an exception is thrown
    MeshStreamer reader2;
    TEST_CHECK_THROWS(reader2.parse_mesh_file(ioss2), SyntaxError);

    // ioss3: missing version entry
    ioss3 << "<feast_mesh_file>" << endl;
    ioss3 << "<header>  " << endl;
    ioss3 << "  version   " << endl;
    ioss3 << "  meshparts 0" << endl;
    ioss3 << "</header>" << endl;
    ioss3 << "<mesh>" << endl;
    ioss3 << "  <header>" << endl;
    ioss3 << "    type conformal" << endl;
    ioss3 << "    coords 2" << endl;
    ioss3 << "    shape quad" << endl;
    ioss3 << "  </header>" << endl;
    ioss3 << "  <counts>" << endl;
    ioss3 << "    verts 4" << endl;
    ioss3 << "    quads 1" << endl;
    ioss3 << "    edges 4" << endl;
    ioss3 << "  </counts>" << endl;
    ioss3 << "  <coords>" << endl;
    ioss3 << "    0.0 0.0" << endl;
    ioss3 << "    1.0 0.0" << endl;
    ioss3 << "    0.0 1.0" << endl;
    ioss3 << "    1.0 1.0" << endl;
    ioss3 << "  </coords>" << endl;
    ioss3 << "</mesh>" << endl;
    ioss3 << "</feast_mesh_file>" << endl;

    // parse the stream and check if an exception is thrown
    MeshStreamer reader3;
    TEST_CHECK_THROWS(reader3.parse_mesh_file(ioss3), SyntaxError);

    // ioss4: nonsense
    ioss4 << "<feast_mesh_file>" << endl;
    ioss4 << "<header>  " << endl;
    ioss4 << "  version 1" << endl;
    ioss4 << "  meshparts 0" << endl;
    ioss4 << "</header>" << endl;
    ioss4 << "<mesh>" << endl;
    ioss4 << "  <header>" << endl;
    ioss4 << "    type conformal" << endl;
    ioss4 << "    coords 2" << endl;
    ioss4 << "    shape quad" << endl;
    ioss4 << "  </header>" << endl;
    ioss4 << "  <counts>" << endl;
    ioss4 << "    verts 4" << endl;
    ioss4 << "    quads 1" << endl;
    ioss4 << "    edges 4" << endl;
    ioss4 << "  </counts>" << endl;
    ioss4 << "  <coords>" << endl;
    ioss4 << "    0.0 0.0" << endl;
    ioss4 << "    1.0 0.0" << endl;
    ioss4 << "    0.0 1.0" << endl;
    ioss4 << "    1.0 1.0" << endl;
    ioss4 << "  </coords> blubb" << endl;
    ioss4 << "</mesh>" << endl;
    ioss4 << "</feast_mesh_file>" << endl;

    // parse the stream and check if an exception is thrown
    MeshStreamer reader4;
    TEST_CHECK_THROWS(reader4.parse_mesh_file(ioss4), SyntaxError);

    // okay
  } // test_1

  // test for structured mesh
  void test_2() const
  {
    using namespace std;
    stringstream ioss;

    // let's write an (awfully ugly) mesh file into the stream
    ioss << "<feast_mesh_file>" << endl;
    ioss << "<header>  " << endl;
    ioss << " version 1" << endl;
    ioss << " chart_file unit_quad_chart.txt" << endl;
    ioss << " meshparts 0" << endl;
    ioss << "</header>" << endl;
    ioss << "   <info>    " << endl;
    ioss << "This file contains a simple unit-square..." << endl;
    ioss << "</info>" << endl;
    ioss << "<mesh>" << endl;
    ioss << " <header>" << endl;
    ioss << "  type structured" << endl;
    ioss << "  shape quad" << endl;

    //ioss << "  coord_file /home/doncamillo/swahlers/Desktop/MeshReader_Beispieldateien/coorddatei.txt" << endl;

    ioss << "  coords 2" << endl;
    ioss << " </header>" << endl;
    ioss << " <info>" << endl;
    ioss << "test" << endl;
    ioss << " </info>" << endl;
    ioss << " <counts>" << endl;
    ioss << "  slices 1 1" << endl;
    ioss << " </counts>" << endl;
    ioss << " <coords>" << endl;
    ioss << " 0.0 0.0" << endl;
    ioss << "  1.0 0.0" << endl;
    ioss << "0.0 1.0" << endl;
    ioss << "  1.0 1.0       " << endl;
    ioss << " </coords>" << endl;
    ioss << "</mesh>" << endl;
    ioss << "</feast_mesh_file>";


    // two mesh reader objects for reading and writing the mesh data
    MeshStreamer reader, writer;

    // pares the stream
    writer.parse_mesh_file(ioss);

    // drop the data into an auxiliary file
    stringstream ioss2;
    writer.write_mesh_file(ioss2);

    // parse the data with the other mesh reader
    reader.parse_mesh_file(ioss2);

    // now check if everything is ok

    // check members
    TEST_CHECK_EQUAL(reader.get_chart_path(), "unit_quad_chart.txt");
    TEST_CHECK_EQUAL(reader.get_num_submeshes(), Index(0));
    TEST_CHECK_EQUAL(reader.get_info(), "This file contains a simple unit-square..." );

    // test the get_mesh function
    MeshStreamer::MeshDataContainer* root_mesh_ptr = reader.get_mesh("rooot");
    TEST_CHECK_EQUAL(root_mesh_ptr, ((MeshStreamer::MeshDataContainer*)nullptr));

    root_mesh_ptr = reader.get_mesh("root");
    TEST_CHECK_NOT_EQUAL(root_mesh_ptr, ((MeshStreamer::MeshDataContainer*)nullptr));

    MeshStreamer::MeshDataContainer& root_mesh(*root_mesh_ptr);

    // check the root mesh data
    TEST_CHECK_EQUAL(root_mesh.name , "root");
    TEST_CHECK_EQUAL(root_mesh.parent , "none");
    TEST_CHECK_EQUAL(root_mesh.chart , "none");
    TEST_CHECK_EQUAL(root_mesh.convert_mesh_type(root_mesh.mesh_type) , "structured");
    TEST_CHECK_EQUAL(root_mesh.convert_shape_type(root_mesh.shape_type) , "quad");
    TEST_CHECK_EQUAL(root_mesh.coord_per_vertex , Index(2));
    TEST_CHECK_EQUAL(root_mesh.vertex_count , Index(4));
    TEST_CHECK_EQUAL(root_mesh.edge_count , Index(0));
    TEST_CHECK_EQUAL(root_mesh.tria_count , Index(0));
    TEST_CHECK_EQUAL(root_mesh.quad_count , Index(0));
    TEST_CHECK_EQUAL(root_mesh.tetra_count , Index(0));
    TEST_CHECK_EQUAL(root_mesh.hexa_count, Index(0));

    TEST_CHECK_EQUAL(root_mesh.slices.size(), Index(2));
    TEST_CHECK_EQUAL((root_mesh.slices).at(0), Index(1));
    TEST_CHECK_EQUAL((root_mesh.slices).at(1), Index(1));

    TEST_CHECK_EQUAL(root_mesh.info, "test");

    // check the root mesh coordinates
    std::vector<std::vector<double> >& coords = root_mesh.coords;

    // reference coordinates
    double coords_ref[] =
    {
      0.0, 0.0,
      1.0, 0.0,
      0.0, 1.0,
      1.0, 1.0
    };

    // loop through the coordinates
    Index count(0);
    bool error = false;
    for(Index i(0); i < 4; ++i)
    {
      for(Index j(0); j < root_mesh.coord_per_vertex; ++j)
      {
        if((coords.at(i)).at(j) != coords_ref[count])
        {
          error = true;
          break;
        }
        ++count;
      }
      if(error)
      {
        break;
      }
    }

    // check if an error occured
    TEST_CHECK_EQUAL(error, false);


    // check if the rest is emtpy
    TEST_CHECK_EQUAL((root_mesh.adjacencies[0][1]).empty(), true);
    TEST_CHECK_EQUAL((root_mesh.adjacencies[0][2]).empty(), true);
    TEST_CHECK_EQUAL((root_mesh.adjacencies[0][3]).empty(), true);

    // check parent indices
    TEST_CHECK_EQUAL((root_mesh.parent_indices[0]).empty(), true);
    TEST_CHECK_EQUAL((root_mesh.parent_indices[1]).empty(), true);
    TEST_CHECK_EQUAL((root_mesh.parent_indices[2]).empty(), true);
    TEST_CHECK_EQUAL((root_mesh.parent_indices[3]).empty(), true);

    // ok, everything is right
  } // test_2


  virtual void run() const
  {
    // run test #0 (checks parsing)
    test_0();
    // run test #1 (checks if errors are found)
    test_1();
    // run test #2 (checks structured mesh)
    test_2();
  }
} mesh_reader_test;
