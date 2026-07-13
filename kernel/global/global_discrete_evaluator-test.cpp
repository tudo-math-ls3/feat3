// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/base_header.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/intern/coarse_fine_cell_mapping.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/shape.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/global/discrete_evaluator.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

class DiscreteEvaluatorTest : public UnitTest
{
  using ShapeType = Shape::Quadrilateral;

  using MeshType = Geometry::ConformalMesh<ShapeType>;

  using DataType = typename MeshType::CoordType;

  using TrafoType = Trafo::Standard::Mapping<MeshType>;

  using SpaceType = Space::Lagrange1::Element<TrafoType>;

  using VertexType = typename MeshType::VertexType;

  using LocalEvaluatorType = Assembly::DiscreteEvaluator;

  using GlobalEvaluatorType = Global::DiscreteEvaluator;

  using TinyVec = Tiny::Vector<Real, 2>;

public:

  DiscreteEvaluatorTest() : UnitTest("GlobalDiscreteEvaluatorTest")
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
    SpaceType space(trafo);

    LAFEM::DenseVector<Real, Index> scalar_values(space.get_num_dofs());
    LAFEM::DenseVectorBlocked<Real, Index, 2> vector_values(space.get_num_dofs());

    {
      const auto& vertices = mesh.get_vertex_set();
      auto scalar_view = scalar_values.elements_view_w();
      auto vector_view = vector_values.elements_view_w();
      for(Index i(0); i < mesh.get_num_vertices(); i++)
      {
        scalar_view[i] = vertices[i][0];
        vector_view[i] = vertices[i];
      }
    }

    // Rank 0 evaluates points
    {
      if(comm.rank() == 0)
      {
        std::vector<VertexType> query{{0.75, 0.75}};
        {
          auto result = Global::DiscreteEvaluator::eval_fe_function(comm, query, scalar_values, space);
          TEST_CHECK(result.size() > 0);
          TEST_CHECK_EQUAL(result[0].mean_value(), DataType(0.75));
        }
        {
          auto result = Global::DiscreteEvaluator::eval_fe_function(comm, query, vector_values, space);
          TEST_CHECK(result.size() > 0);
          TEST_CHECK_EQUAL(result[0].mean_value()[0], DataType(0.75));
          TEST_CHECK_EQUAL(result[0].mean_value()[1], DataType(0.75));
        }
        {
          auto result = Global::DiscreteEvaluator::eval_fe_gradient(comm, query, scalar_values, space);
          TEST_CHECK(result.size() > 0);
          TEST_CHECK_EQUAL(result[0].mean_value()[0], DataType(1));
          TEST_CHECK_EQUAL(result[0].mean_value()[1], DataType(0));
        }
        {
          auto result = Global::DiscreteEvaluator::eval_fe_gradient(comm, query, vector_values, space);
          TEST_CHECK(result.size() > 0);
          TEST_CHECK_EQUAL(result[0].mean_value()[0][0], DataType(1));
          TEST_CHECK_EQUAL(result[0].mean_value()[1][0], DataType(0));
          TEST_CHECK_EQUAL(result[0].mean_value()[0][1], DataType(0));
          TEST_CHECK_EQUAL(result[0].mean_value()[1][1], DataType(1));
        }
      }
      else
      {
        std::vector<VertexType> query;
        {
          auto data = Global::DiscreteEvaluator::eval_fe_function(comm, query, scalar_values, space);
        }
        {
          auto data = Global::DiscreteEvaluator::eval_fe_function(comm, query, vector_values, space);
        }
        {
          auto data = Global::DiscreteEvaluator::eval_fe_gradient(comm, query, scalar_values, space);
        }
        {
          auto data = Global::DiscreteEvaluator::eval_fe_gradient(comm, query, vector_values, space);
        }

      }

      if(comm.rank() == 0)
      {
        std::vector<VertexType> query{{1.75, 0.25}};
        {
          auto result = Global::DiscreteEvaluator::eval_fe_function(comm, query, scalar_values, space);
          TEST_CHECK(result.size() > 0);
          TEST_CHECK_EQUAL(result[0].mean_value(), DataType(1.75));
        }
        {
          auto result = Global::DiscreteEvaluator::eval_fe_function(comm, query, vector_values, space);
          TEST_CHECK(result.size() > 0);
          TEST_CHECK_EQUAL(result[0].mean_value()[0], DataType(1.75));
          TEST_CHECK_EQUAL(result[0].mean_value()[1], DataType(0.25));
        }
        {
          auto result = Global::DiscreteEvaluator::eval_fe_gradient(comm, query, scalar_values, space);
          TEST_CHECK(result.size() > 0);
          TEST_CHECK_EQUAL(result[0].mean_value()[0], DataType(1));
          TEST_CHECK_EQUAL(result[0].mean_value()[1], DataType(0));
        }
        {
          auto result = Global::DiscreteEvaluator::eval_fe_gradient(comm, query, vector_values, space);
          TEST_CHECK(result.size() > 0);
          TEST_CHECK_EQUAL(result[0].mean_value()[0][0], DataType(1));
          TEST_CHECK_EQUAL(result[0].mean_value()[1][0], DataType(0));
          TEST_CHECK_EQUAL(result[0].mean_value()[0][1], DataType(0));
          TEST_CHECK_EQUAL(result[0].mean_value()[1][1], DataType(1));
        }
      }
      else
      {
        std::vector<VertexType> query;
        {
          auto data = Global::DiscreteEvaluator::eval_fe_function(comm, query, scalar_values, space);
        }
        {
          auto data = Global::DiscreteEvaluator::eval_fe_function(comm, query, vector_values, space);
        }
        {
          auto data = Global::DiscreteEvaluator::eval_fe_gradient(comm, query, scalar_values, space);
        }
        {
          auto data = Global::DiscreteEvaluator::eval_fe_gradient(comm, query, vector_values, space);
        }
      }
    }

    // All ranks search for points
    {
      std::vector<VertexType> query
      {
        VertexType{0.75, 0.75},
        VertexType{1.75, 0.25}
      };

      {
        auto result = Global::DiscreteEvaluator::eval_fe_function(comm, query, scalar_values, space);
        TEST_CHECK_EQUAL(result.size(), 2);
        TEST_CHECK_EQUAL(result[0].mean_value(), DataType(0.75));
        TEST_CHECK_EQUAL(result[1].mean_value(), DataType(1.75));
      }
      {
        auto result = Global::DiscreteEvaluator::eval_fe_function(comm, query, vector_values, space);
        TEST_CHECK_EQUAL(result.size(), 2);
        TEST_CHECK_EQUAL(result[0].mean_value()[0], DataType(0.75));
        TEST_CHECK_EQUAL(result[0].mean_value()[1], DataType(0.75));
        TEST_CHECK_EQUAL(result[1].mean_value()[0], DataType(1.75));
        TEST_CHECK_EQUAL(result[1].mean_value()[1], DataType(0.25));
      }
      {
        auto result = Global::DiscreteEvaluator::eval_fe_gradient(comm, query, scalar_values, space);
        TEST_CHECK_EQUAL(result.size(), 2);
        TEST_CHECK_EQUAL(result[0].mean_value()[0], DataType(1));
        TEST_CHECK_EQUAL(result[0].mean_value()[1], DataType(0));
        TEST_CHECK_EQUAL(result[1].mean_value()[0], DataType(1));
        TEST_CHECK_EQUAL(result[1].mean_value()[1], DataType(0));
      }
      {
        auto result = Global::DiscreteEvaluator::eval_fe_gradient(comm, query, vector_values, space);
        TEST_CHECK_EQUAL(result.size(), 2);
        TEST_CHECK_EQUAL(result[0].mean_value()[0][0], DataType(1));
        TEST_CHECK_EQUAL(result[0].mean_value()[1][0], DataType(0));
        TEST_CHECK_EQUAL(result[0].mean_value()[0][1], DataType(0));
        TEST_CHECK_EQUAL(result[0].mean_value()[1][1], DataType(1));
        TEST_CHECK_EQUAL(result[1].mean_value()[0][0], DataType(1));
        TEST_CHECK_EQUAL(result[1].mean_value()[1][0], DataType(0));
        TEST_CHECK_EQUAL(result[1].mean_value()[0][1], DataType(0));
        TEST_CHECK_EQUAL(result[1].mean_value()[1][1], DataType(1));
      }
    }
  }
};

static const DiscreteEvaluatorTest discrete_evaluator_test;
