// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/intern/coarse_fine_cell_mapping.hpp>
#include <kernel/shape.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/global/inverse_mapping.hpp>
#include <kernel/util/dist.hpp>
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

class InverseMappingTest : public UnitTest
{
  using ShapeType = Shape::Quadrilateral;

  using MeshType = Geometry::ConformalMesh<ShapeType>;

  using DataType = typename MeshType::CoordType;

  using TrafoType = Trafo::Standard::Mapping<MeshType>;

  using VertexType = typename MeshType::VertexType;

  using LocalInverseMappingType = Trafo::InverseMapping<TrafoType, DataType>;

  using GlobalInverseMappingType = Global::InverseMapping<LocalInverseMappingType>;

  using GlobalInverseMappingDataType = typename GlobalInverseMappingType::InvMapDataType;
public:

  InverseMappingTest() : UnitTest("GlobalInverseMappingTest")
  {

  }

  void run() const override
  {
    Dist::Comm comm(Dist::Comm::world());

    if(comm.size() > 1)
    {
      test_dist(comm);
    }
  }

  void test_dist(const Dist::Comm& comm) const
  {
    Geometry::RefinedUnitCubeFactory<MeshType> factory(1);
    MeshType mesh(factory);

    // Translate mesh
    mesh.transform(VertexType(0.0), VertexType(0.0), VertexType{DataType(comm.rank()), 0.0});

    // Globally, we now have a row of subdivided unit squares, one per rank, with cells 0 to 3 on each rank
    //   Rank 0  Rank 1
    // *---*---*---*---*
    // | 2 | 3 | 2 | 3 |
    // *---*---*---*---*
    // | 0 | 1 | 0 | 1 |
    // *---*---*---*---*

    TrafoType trafo(mesh);
    LocalInverseMappingType inv_mapping(trafo);
    GlobalInverseMappingType global_mapping(comm, inv_mapping);

    // Rank 0 searches for points
    {
      if(comm.rank() == 0)
      {
        std::vector<VertexType> query{{0.75, 0.75}};
        auto result = global_mapping.unmap_points(query);

        TEST_CHECK(result.size() > 0);

        GlobalInverseMappingDataType& data = result[0];

        // Cell should have been found in cell 3 of rank 0
        TEST_CHECK_EQUAL(data.cells.size(), 1);
        TEST_CHECK_EQUAL(data.cells.front(), 3);
        TEST_CHECK_EQUAL(data.ranks.front(), 0);
      }
      else
      {
        std::vector<VertexType> query;
        auto data = global_mapping.unmap_points(query);
      }

      if(comm.rank() == 0)
      {
        std::vector<VertexType> query{{1.75, 0.25}};
        auto result = global_mapping.unmap_points(query);

        TEST_CHECK(result.size() > 0);

        GlobalInverseMappingDataType& data = result[0];

        // Cell should have been found in cell 1 of rank 1
        TEST_CHECK_EQUAL(data.cells.size(), 1);
        TEST_CHECK_EQUAL(data.cells.front(), 1);
        TEST_CHECK_EQUAL(data.ranks.front(), 1);
      }
      else
      {
        std::vector<VertexType> query;
        auto data = global_mapping.unmap_points(query);
      }

      if(comm.rank() == 0)
      {
        std::vector<VertexType> query{{1.5, 0.25}};
        auto result = global_mapping.unmap_points(query);


        GlobalInverseMappingDataType& data = result[0];

        // Cell should have been found in cell 0 of rank 1 and cell 1 of rank 1
        TEST_CHECK_EQUAL(data.cells.size(), 2);
        TEST_CHECK_EQUAL(data.cells.front(), 0);
        TEST_CHECK_EQUAL(data.ranks.front(), 1);
        TEST_CHECK_EQUAL(data.cells.back(), 1);
        TEST_CHECK_EQUAL(data.ranks.back(), 1);
      }
      else
      {
        std::vector<VertexType> query;
        auto data = global_mapping.unmap_points(query);
      }
    }

    // All ranks search for points
    {
      std::vector<VertexType> query
      {
        VertexType{0.75, 0.75},
        VertexType{1.75, 0.25}
      };

      auto result = global_mapping.unmap_points(query);

      TEST_CHECK_EQUAL(result.size(), 2);

      GlobalInverseMappingDataType& data_a = result[0];

      // Cell should have been found in cell 3 of rank 0
      TEST_CHECK_EQUAL(data_a.cells.size(), 1);
      TEST_CHECK_EQUAL(data_a.cells.front(), 3);
      TEST_CHECK_EQUAL(data_a.ranks.front(), 0);

      GlobalInverseMappingDataType& data_b = result[1];

      // Cell should have been found in cell 1 of rank 1
      TEST_CHECK_EQUAL(data_b.cells.size(), 1);
      TEST_CHECK_EQUAL(data_b.cells.front(), 1);
      TEST_CHECK_EQUAL(data_b.ranks.front(), 1);
    }
  }
};

static const InverseMappingTest inverse_mapping_test;
