// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/geometry/shape_convert_factory.hpp>
#include <kernel/geometry/test_aux/copy_comp_set.hpp>
#include <kernel/geometry/atlas/bezier.hpp>
#include <kernel/trafo/isoparam/mapping.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

template<typename DataType_>
class IsoparamTrafoVolumeTest
  : public TestSystem::TaggedTest<Archs::None, DataType_>
{
public:
  IsoparamTrafoVolumeTest() :
    TestSystem::TaggedTest<Archs::None, DataType_>("IsoparamTrafoVolumeTest")
  {
  }

  virtual ~IsoparamTrafoVolumeTest()
  {
  }

  virtual void run() const override
  {
    test_2d_quad();
  }

  void test_2d_quad() const
  {
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.4));

    typedef Shape::Hypercube<2> ShapeType;
    typedef Shape::Hypercube<1> FacetType;
    typedef Geometry::ConformalMesh<ShapeType, 2, DataType_> MeshType;
    typedef Geometry::MeshPart<MeshType> MeshPartType;
    typedef Trafo::Isoparam::Mapping<MeshType, 2> TrafoType;
    typedef Tiny::Vector<DataType_, 2> PointType;

    Geometry::ReferenceCellFactory<ShapeType> factory;

    MeshType mesh(factory);

    // define the corner vertex coordinates
    static const DataType_ vtx[4*2] =
    {
      DataType_(-4), DataType_(-2),
      DataType_( 5), DataType_(-3),
      DataType_(-4), DataType_( 4),
      DataType_( 2), DataType_( 6)
    };
    Geometry::TestAux::copy_vtx((&mesh)->get_vertex_set(), vtx);

    // create boundary mesh part
    const Index num_mp_ents[] = {5,4,0};
    MeshPartType mesh_part(num_mp_ents, true);
    {
      auto& trg_v = mesh_part.template get_target_set<0>();
      trg_v[0] = Index(0);
      trg_v[1] = Index(1);
      trg_v[2] = Index(3);
      trg_v[3] = Index(2);
      trg_v[4] = Index(0);
      auto& trg_e = mesh_part.template get_target_set<1>();
      trg_e[0] = Index(0);
      trg_e[1] = Index(3);
      trg_e[2] = Index(1);
      trg_e[3] = Index(2);
      auto& idx_e = mesh_part.template get_index_set<1,0>();
      idx_e(0u, 0) = Index(0);
      idx_e(0u, 1) = Index(1);
      idx_e(1u, 0) = Index(1);
      idx_e(1u, 1) = Index(2);
      idx_e(2u, 0) = Index(2);
      idx_e(2u, 1) = Index(3);
      idx_e(3u, 0) = Index(3);
      idx_e(3u, 1) = Index(4);
    }

    // create Bezier chart
    Geometry::Atlas::Bezier<MeshType> chart(true);
    chart.push_vertex( PointType{DataType_(-4), DataType_(-2)});
    chart.push_vertex( PointType{DataType_( 5), DataType_(-3)});
    chart.push_control(PointType{DataType_(6.5), DataType_(2.5)});
    chart.push_vertex( PointType{DataType_( 2), DataType_( 6)});
    chart.push_control(PointType{DataType_(-2), DataType_(8)});
    chart.push_vertex( PointType{DataType_(-4), DataType_( 4)});
    chart.push_control(PointType{DataType_(-2), DataType_(1)});
    chart.push_vertex( PointType{DataType_(-4), DataType_(-2)});

    TrafoType trafo(mesh);

    // compute quad volume before adding edge midpoints to trafo
    {
      typename TrafoType::template Evaluator<ShapeType, DataType_>::Type trafo_eval(trafo);
      trafo_eval.prepare(Index(0));
      const DataType_ vol = trafo_eval.volume();
      trafo_eval.finish();
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, DataType_(57), eps);
    }

    // compute edge lengths before adding edge midpoints to trafo
    {
      typename TrafoType::template Evaluator<FacetType, DataType_>::Type trafo_eval(trafo);
      trafo_eval.prepare(Index(0));
      const DataType_ len_0 = trafo_eval.volume();
      trafo_eval.finish();
      trafo_eval.prepare(Index(1));
      const DataType_ len_1 = trafo_eval.volume();
      trafo_eval.finish();
      trafo_eval.prepare(Index(2));
      const DataType_ len_2 = trafo_eval.volume();
      trafo_eval.finish();
      trafo_eval.prepare(Index(3));
      const DataType_ len_3 = trafo_eval.volume();
      trafo_eval.finish();
      TEST_CHECK_EQUAL_WITHIN_EPS(len_0, Math::sqrt(DataType_(82)), eps);
      TEST_CHECK_EQUAL_WITHIN_EPS(len_1, Math::sqrt(DataType_(40)), eps);
      TEST_CHECK_EQUAL_WITHIN_EPS(len_2, Math::sqrt(DataType_(36)), eps);
      TEST_CHECK_EQUAL_WITHIN_EPS(len_3, Math::sqrt(DataType_(90)), eps);
    }

    // add chart to trafo
    trafo.add_meshpart_chart(mesh_part, chart);

    // compute quad volume of isoparametrically deformed quad
    {
      typename TrafoType::template Evaluator<ShapeType, DataType_>::Type trafo_eval(trafo);
      trafo_eval.prepare(Index(0));
      const DataType_ vol = trafo_eval.volume();
      trafo_eval.finish();
      TEST_CHECK_EQUAL_WITHIN_EPS(vol, DataType_(209)/DataType_(3), eps);
    }

    // compute edge lengths of isoparametrically deformed quad
    {
      typename TrafoType::template Evaluator<FacetType, DataType_>::Type trafo_eval(trafo);
      trafo_eval.prepare(Index(0));
      const DataType_ len_0 = trafo_eval.volume();
      trafo_eval.finish();
      trafo_eval.prepare(Index(1));
      const DataType_ len_1 = trafo_eval.volume();
      trafo_eval.finish();
      trafo_eval.prepare(Index(2));
      const DataType_ len_2 = trafo_eval.volume();
      trafo_eval.finish();
      trafo_eval.prepare(Index(3));
      const DataType_ len_3 = trafo_eval.volume();
      trafo_eval.finish();
      //std::cout << std::setprecision(15) << len_0 << " / " << len_1 << " / " << len_2 << " / " << len_3 << std::endl;
      TEST_CHECK_EQUAL_WITHIN_EPS(len_0, DataType_(9.0553851381374166266), eps);
      TEST_CHECK_EQUAL_WITHIN_EPS(len_1, DataType_(7.2592839594939513807), eps);
      TEST_CHECK_EQUAL_WITHIN_EPS(len_2, DataType_(6.4187043030908643875), eps);
      TEST_CHECK_EQUAL_WITHIN_EPS(len_3, DataType_(10.148862612445443273), eps);
    }
  }
};

IsoparamTrafoVolumeTest<double> isoparam_trafo_volume_test_double;
