// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/factory.hpp>
#include <kernel/geometry/intern/coarse_fine_cell_mapping.hpp>
#include <kernel/shape.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/random.hpp>
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

template<typename Shape_>
class RandomPointHelper;

template<int dim_>
class RandomPointHelper<Shape::Simplex<dim_>>
{
public:
  template<typename DT_>
  static void inner(Tiny::Vector<DT_, dim_, dim_>& p, Random& rng)
  {
    DT_ v[dim_+1];
    DT_ s = DT_(0);
    for(int i(0); i <= dim_; ++i)
      s += (v[i] = rng(DT_(0.05), DT_(0.95)));
    for(int i(0); i < dim_; ++i)
      p[i] = v[i] / s;
  }
};

template<int dim_>
class RandomPointHelper<Shape::Hypercube<dim_>>
{
public:
  template<typename DT_>
  static void inner(Tiny::Vector<DT_, dim_, dim_>& p, Random& rng)
  {
    for(int i(0); i < dim_; ++i)
      p[i] = rng(-DT_(0.95), DT_(0.95));
  }
};

template<typename DataType_>
class InverseMappingTest :
  public TestSystem::UnitTest
{
public:
  InverseMappingTest() :
    TestSystem::UnitTest("InverseMappingTest", Type::Traits<DataType_>::name())
  {
  }

  virtual void run() const override
  {
    test_unmap<Shape::Simplex<2>>(2, 0.05);
    test_unmap<Shape::Simplex<3>>(1, 0.03);
    test_unmap<Shape::Hypercube<2>>(3, 0.03);
    test_unmap<Shape::Hypercube<3>>(2, 0.06);
  }

  template<typename Mesh_>
  static void distort_mesh(Mesh_& mesh, Random& rng, DataType_ dsh)
  {
    const DataType_ pi2 = DataType_(2) * Math::pi<DataType_>();

    auto& vtx = mesh.get_vertex_set();
    if(vtx.get_num_coords() == 2)
    {
      // 2D
      for(Index i(0); i < vtx.get_num_vertices(); ++i)
      {
        DataType_ t = rng(DataType_(0), pi2);
        if((vtx[i][0] > DataType_(0.001)) && (vtx[i][0] < DataType_(0.999)))
          vtx[i][0] += dsh * Math::cos(t);
        if((vtx[i][1] > DataType_(0.001)) && (vtx[i][1] < DataType_(0.999)))
          vtx[i][1] += dsh * Math::sin(t);
      }
    }
    else
    {
      // 3D
      for(Index i(0); i < vtx.get_num_vertices(); ++i)
      {
        DataType_ t = rng(DataType_(0), pi2);
        DataType_ r = rng(DataType_(0), DataType_(1));
        DataType_ b = DataType_(int(rng.next() & 0x1) != 0 ? +1 : -1);

        if((vtx[i][0] > DataType_(0.001)) && (vtx[i][0] < DataType_(0.999)))
          vtx[i][0] += dsh * r * Math::cos(t);
        if((vtx[i][1] > DataType_(0.001)) && (vtx[i][1] < DataType_(0.999)))
          vtx[i][1] += dsh * r * Math::sin(t);
        if((vtx[i][2] > DataType_(0.001)) && (vtx[i][2] < DataType_(0.999)))
          vtx[i][2] += dsh * b * Math::sqrt(DataType_(1) - r*r);
      }
    }
  }

  template<typename Shape_>
  void test_unmap(const Index level, const DataType_ disto) const
  {
    typedef Shape_ ShapeType;
    static constexpr int dim = ShapeType::dimension;
    typedef Geometry::ConformalMesh<ShapeType, dim, DataType_> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // compute tolerance
    const DataType_ tol = TestSystem::tol<DataType_>();

    // create a quad mesh
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);

    const Index num_elems = mesh.get_num_elements();

    // distort mesh
    Random rng;
    std::cout << "RNG Seed: " << rng.get_seed() << "\n";
    distort_mesh(mesh, rng, disto);

    // create a quad-trafo
    TrafoType trafo(mesh);

    // create a trafo evaluator for forward mapping
    typedef typename TrafoType::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    typedef typename TrafoEvaluator::template ConfigTraits<TrafoTags::img_point>::EvalDataType TrafoEvalData;
    TrafoEvaluator trafo_eval(trafo);
    TrafoEvalData trafo_data;

    typename TrafoEvaluator::DomainPointType dom_point;

    // create an inverse trafo mapping
    Trafo::InverseMapping<TrafoType, DataType_> inv_mapping(trafo);

    // loop over all elements
    for(Index elem(0); elem < num_elems; ++elem)
    {
      // prepare evaluator
      trafo_eval.prepare(elem);

      // test a few random inner points
      for(int j(0); j < 7; ++j)
      {
        // create random inner point
        RandomPointHelper<ShapeType>::inner(dom_point, rng);

        // map to world coords
        trafo_eval(trafo_data, dom_point);

        // unmap point (ignore failures)
        auto inv_data = inv_mapping.unmap_point(trafo_data.img_point, true);

        // make sure that we have exactly one element
        TEST_CHECK_EQUAL(inv_data.size(), std::size_t(1));

        // make sure that this is our element
        TEST_CHECK_EQUAL(inv_data.cells.front(), elem);

        // compute domain point distance
        auto dist_point = dom_point - inv_data.dom_points.front();
        DataType_ dist_norm = dist_point.norm_euclid();

        // check domain point
        TEST_CHECK(dist_norm < tol);
      }

      trafo_eval.finish();
    }
  }
}; // class InverseMappingTest

InverseMappingTest<double> inverse_mapping_test_double;

class HierarchicalInverseMappingTest : public UnitTest
{
  using ShapeType = Shape::Quadrilateral;

  using MeshType = Geometry::ConformalMesh<ShapeType>;

  using DataType = typename MeshType::CoordType;

  using TrafoType = Trafo::Standard::Mapping<MeshType>;

  using VertexType = typename MeshType::VertexType;

  using InverseMappingType = Trafo::InverseMapping<TrafoType, DataType>;

  using HierarchicalInverseMappingType = Trafo::HierarchicalInverseMapping<InverseMappingType, Geometry::Intern::CoarseFineCellMapping<MeshType>>;

  using InverseMappingDataType = typename HierarchicalInverseMappingType::InvMapDataType;
public:

  HierarchicalInverseMappingTest() : UnitTest("HierarchicalInverseMappingTest")
  {

  }

  void run() const override
  {
    test_single_mesh();
    test_mg();
  }

  void test_single_mesh() const
  {
    Geometry::RefinedUnitCubeFactory<MeshType> factory(1);

    MeshType mesh(factory);
    TrafoType trafo(mesh);
    Trafo::InverseMapping<TrafoType, DataType> inverse_mapping(trafo);

    HierarchicalInverseMappingType hi_mapping;
    hi_mapping.push_level(inverse_mapping);

    {
      VertexType point{0.25, 0.25};
      InverseMappingDataType locator_result = hi_mapping.unmap_point(point);
      InverseMappingDataType im_result = inverse_mapping.unmap_point(point);

      TEST_CHECK(compare(locator_result, im_result));
    }

    {
      VertexType point{0.5, 0.25};
      InverseMappingDataType locator_result = hi_mapping.unmap_point(point);
      InverseMappingDataType im_result = inverse_mapping.unmap_point(point);

      TEST_CHECK(compare(locator_result, im_result));
    }

    {
      VertexType point{0.5, 0.25};
      InverseMappingDataType locator_result = hi_mapping.unmap_point(point);
      InverseMappingDataType im_result = inverse_mapping.unmap_point(point);

      TEST_CHECK(compare(locator_result, im_result));
    }

    {
      VertexType point{0.5, 0.5};
      InverseMappingDataType locator_result = hi_mapping.unmap_point(point);
      InverseMappingDataType im_result = inverse_mapping.unmap_point(point);

      TEST_CHECK(compare(locator_result, im_result));
    }
  }

  void test_mg() const
  {
    std::vector<MeshType> meshes;
    meshes.reserve(2);
    std::vector<TrafoType> trafos;
    trafos.reserve(2);
    std::vector<InverseMappingType> mappings;
    mappings.reserve(2);

    Geometry::UnitCubeFactory<MeshType> factory;
    meshes.emplace_back(factory);
    meshes.emplace_back(Geometry::StandardRefinery<MeshType>(meshes.back()).make());

    trafos.emplace_back(meshes.front());
    trafos.emplace_back(meshes.back());

    mappings.emplace_back(trafos.front());
    mappings.emplace_back(trafos.back());

    Geometry::Intern::CoarseFineCellMapping<MeshType> cf(meshes.back(), meshes.front());

    HierarchicalInverseMappingType hi_mapping;
    // Push coarse level
    hi_mapping.push_level(mappings.front());
    // Push fine level
    hi_mapping.push_level(mappings.back(), cf);

    InverseMappingType& inverse_mapping = mappings.back();

    {
      VertexType point{0.25, 0.25};
      InverseMappingDataType locator_result = hi_mapping.unmap_point(point);
      InverseMappingDataType im_result = inverse_mapping.unmap_point(point);

      TEST_CHECK(compare(locator_result, im_result));
    }

    {
      VertexType point{0.5, 0.25};
      InverseMappingDataType locator_result = hi_mapping.unmap_point(point);
      InverseMappingDataType im_result = inverse_mapping.unmap_point(point);

      TEST_CHECK(compare(locator_result, im_result));
    }

    {
      VertexType point{0.5, 0.25};
      InverseMappingDataType locator_result = hi_mapping.unmap_point(point);
      InverseMappingDataType im_result = inverse_mapping.unmap_point(point);

      TEST_CHECK(compare(locator_result, im_result));
    }

    {
      VertexType point{0.5, 0.5};
      InverseMappingDataType locator_result = hi_mapping.unmap_point(point);
      InverseMappingDataType im_result = inverse_mapping.unmap_point(point);

      TEST_CHECK(compare(locator_result, im_result));
    }
  }

  static bool compare(const InverseMappingDataType& a, const InverseMappingDataType& b)
  {
    std::set<Index> cells_a;
    std::set<Index> cells_b;

    for(Index cell : a.cells)
    {
      cells_a.insert(cell);
    }

    for(Index cell : b.cells)
    {
      cells_b.insert(cell);
    }

    if(cells_a != cells_b)
    {
      return false;
    }

    if(a.dom_points.size() != b.dom_points.size())
    {
      return false;
    }

    for(const VertexType& p : a.dom_points)
    {
      bool found_match = false;
      for(const VertexType& other : b.dom_points)
      {
        found_match |= (p - other).norm_euclid() < 1E-12;
      }

      if(!found_match)
      {
        return false;
      }
    }

    return true;
  }
};

static const HierarchicalInverseMappingTest hierarchical_inverse_mapping_test;
