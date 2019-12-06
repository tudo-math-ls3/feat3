// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// This debug tool tests the order of the isoparametric trafo for a
// 3D sphere domain by computing its volume and surface area
// and by computing the mean distance error of the isoparameteric
// surface to the analytic sphere boundary.

#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_extruder.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/isoparam/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/lagrange3/element.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/runtime.hpp>

using namespace FEAT;
using namespace FEAT::Geometry;

namespace DbgIsoParam
{
#ifdef FEAT_HAVE_QUADMATH
  typedef __float128 DT;
#else
  typedef double DT;
#endif
  typedef Shape::Hexahedron ShapeType3D;
  typedef ConformalMesh<ShapeType3D, 3, DT> MeshType3D;
  typedef MeshPart<MeshType3D> MeshPartType3D;
  typedef MeshAtlas<MeshType3D> AtlasType3D;
  typedef RootMeshNode<MeshType3D> MeshNodeType3D;

  std::shared_ptr<AtlasType3D> make_atlas_3d()
  {
    auto atlas = std::make_shared<Geometry::MeshAtlas<MeshType3D>>();
    auto sphere = new Geometry::Atlas::Sphere<MeshType3D>(0.0, 0.0, 0.0, 1.0);
    atlas->add_mesh_chart("sphere", sphere);
    return atlas;
  }

  std::shared_ptr<MeshNodeType3D> make_mesh_node(std::shared_ptr<AtlasType3D> atlas)
  {
    Geometry::UnitCubeFactory<MeshType3D> factory_3d;
    auto mesh = new MeshType3D(factory_3d);

    // move coords from [0,1] to [-0.7,+0.7]
    {
      auto& vtx = mesh->get_vertex_set();
      for(Index i(0); i < vtx.get_num_vertices(); ++i)
      {
        (vtx[i][0] *= 1.4) -= 0.7;
        (vtx[i][1] *= 1.4) -= 0.7;
        (vtx[i][2] *= 1.4) -= 0.7;
      }
    }

    auto node = std::make_shared<Geometry::RootMeshNode<MeshType3D>>(mesh, atlas.get());

    Geometry::BoundaryFactory<MeshType3D> bnd_factory(*mesh);
    auto part = new MeshPartType3D(bnd_factory);

    auto chart = atlas->find_mesh_chart("sphere");
    node->add_mesh_part("bnd", part, "sphere", chart);

    node->adapt();

    return node;
  }

  // calculate the volume/area of the domain; should converge to pi
  template<typename Trafo_>
  DT calc_vol(Trafo_& trafo, int degree)
  {
    typedef DT DataType;
    typedef typename Trafo_::ShapeType ShapeType;
    typedef typename Trafo_::template Evaluator<ShapeType, DataType>::Type TrafoEvaluator;
    typename TrafoEvaluator::template ConfigTraits<TrafoTags::jac_det>::EvalDataType trafo_data;
    TrafoEvaluator trafo_eval(trafo);

    Cubature::DynamicFactory cubature_factory("gauss-legendre:" + stringify(degree));
    typename Assembly::Intern::CubatureTraits<TrafoEvaluator>::RuleType
      cubature_rule(Cubature::ctor_factory, cubature_factory);

    DataType vol = DataType(0);

    for(Index cell = trafo_eval.begin(); cell < trafo_eval.end(); ++cell)
    {
      trafo_eval.prepare(cell);
      for(int k(0); k < cubature_rule.get_num_points(); ++k)
      {
        trafo_eval(trafo_data, cubature_rule.get_point(k));
        vol += trafo_data.jac_det * cubature_rule.get_weight(k);
      }
      trafo_eval.finish();
    }

    return vol;
  }

  // calculate the circumference of the domain; should converge to 2*pi
  template<typename Trafo_, typename MeshPartType_>
  DT calc_circum(Trafo_& trafo, const MeshPartType_& mpart, int degree)
  {
    typedef DT DataType;

    typedef typename Trafo_::ShapeType ShapeType;
    static constexpr int facet_dim = ShapeType::dimension - 1;

    typedef typename Trafo_::template Evaluator<Shape::Hypercube<facet_dim>, DataType>::Type TrafoEvaluator;
    typename TrafoEvaluator::template ConfigTraits<TrafoTags::jac_det>::EvalDataType trafo_data;
    TrafoEvaluator trafo_eval(trafo);

    Cubature::DynamicFactory cubature_factory("gauss-legendre:" + stringify(degree));
    typename Assembly::Intern::CubatureTraits<TrafoEvaluator>::RuleType
      cubature_rule(Cubature::ctor_factory, cubature_factory);

    DataType vol = DataType(0);

    const auto& trg = mpart.template get_target_set<facet_dim>();

    //for(Index cell = trafo_eval.begin(); cell < trafo_eval.end(); ++cell)
    for(Index i(0); i < trg.get_num_entities(); ++i)
    {
      trafo_eval.prepare(trg[i]);
      for(int k(0); k < cubature_rule.get_num_points(); ++k)
      {
        trafo_eval(trafo_data, cubature_rule.get_point(k));
        vol += trafo_data.jac_det * cubature_rule.get_weight(k);
      }
      trafo_eval.finish();
    }

    return vol;
  }

  // calculate the mean distance to the boundary; should converge to 0
  template<typename Trafo_, typename MeshPartType_>
  DT calc_mean_dist(Trafo_& trafo, const MeshPartType_& mpart, int degree)
  {
    typedef DT DataType;

    typedef typename Trafo_::ShapeType ShapeType;
    static constexpr int facet_dim = ShapeType::dimension - 1;

    typedef typename Trafo_::template Evaluator<Shape::Hypercube<facet_dim>, DataType>::Type TrafoEvaluator;
    typename TrafoEvaluator::template ConfigTraits<TrafoTags::img_point>::EvalDataType trafo_data;
    TrafoEvaluator trafo_eval(trafo);

    Cubature::DynamicFactory cubature_factory("gauss-legendre:" + stringify(degree));
    typename Assembly::Intern::CubatureTraits<TrafoEvaluator>::RuleType
      cubature_rule(Cubature::ctor_factory, cubature_factory);

    typename TrafoEvaluator::ImagePointType img_point;

    DataType err = DataType(0);

    const auto& trg = mpart.template get_target_set<facet_dim>();

    for(Index i(0); i < trg.get_num_entities(); ++i)
    {
      trafo_eval.prepare(trg[i]);
      for(int k(0); k < cubature_rule.get_num_points(); ++k)
      {
        trafo_eval(trafo_data, cubature_rule.get_point(k));
        img_point = trafo_data.img_point;
        img_point.normalise();
        img_point -= trafo_data.img_point;
        err += img_point.norm_euclid_sqr();
        //err += img_point.norm_euclid();
      }
      trafo_eval.finish();
    }

    // compute mean error
    return Math::sqrt(err / DataType(trg.get_num_entities()));
    //return err / DataType(trg.get_num_entities());
  }

  void run(int argc, char* argv[])
  {
    SimpleArgParser args(argc, argv);

    auto atlas = make_atlas_3d();
    typedef typename std::remove_reference<decltype(*atlas)>::type AtlasType;
    typedef typename AtlasType::MeshType MeshType;

    typedef Trafo::Standard::Mapping<MeshType> TrafoDeg1;
    typedef Trafo::Isoparam::Mapping<MeshType, 2> TrafoDeg2;
    typedef Trafo::Isoparam::Mapping<MeshType, 3> TrafoDeg3;

    // create an empty atlas and a root mesh node
    auto node = make_mesh_node(atlas);

    const auto& chart = *atlas->find_mesh_chart("sphere");

    Index lvl_min = 0;
    Index lvl_max = 5;

    args.parse("level", lvl_max);

    static constexpr std::size_t ns = 3;
    std::vector<std::array<DT,ns>> vols(lvl_max+1), verrs(lvl_max+1);
    std::vector<std::array<DT,ns>> cirs(lvl_max+1), cerrs(lvl_max+1), merrs(lvl_max+1);

    static const DT pi = Math::pi<DT>();

    // refine
    for(Index lvl(0); lvl <= lvl_max; ++lvl)
    {
      if(lvl > 0)
      {
        node = std::shared_ptr<Geometry::RootMeshNode<MeshType>>(node->refine());
      }

      // get our mesh
      MeshType& mesh = *node->get_mesh();
      const auto& mpart = *node->find_mesh_part("bnd");

      // create trafos
      TrafoDeg1 trafo_1(mesh);
      TrafoDeg2 trafo_2(mesh);
      TrafoDeg3 trafo_3(mesh);

      // add chart and boundary mesh part
      trafo_2.add_meshpart_chart(mpart, chart);
      trafo_3.add_meshpart_chart(mpart, chart);

      // rerference volume and circumference
      const DT ref_vol = (pi * DT(4)) / DT(3);
      const DT ref_cir = pi * DT(4);

      // volume errors
      verrs.at(lvl)[0] = Math::abs(ref_vol - calc_vol(trafo_1, 2));
      verrs.at(lvl)[1] = Math::abs(ref_vol - calc_vol(trafo_2, 3));
      verrs.at(lvl)[2] = Math::abs(ref_vol - calc_vol(trafo_3, 4));

      // circumference (surface area) errors
      cerrs.at(lvl)[0] = Math::abs(ref_cir - calc_circum(trafo_1, mpart, 2));
      cerrs.at(lvl)[1] = Math::abs(ref_cir - calc_circum(trafo_2, mpart, 3));
      cerrs.at(lvl)[2] = Math::abs(ref_cir - calc_circum(trafo_3, mpart, 4));

      // mean distance errors
      merrs.at(lvl)[0] = calc_mean_dist(trafo_1, mpart, 2);
      merrs.at(lvl)[1] = calc_mean_dist(trafo_2, mpart, 3);
      merrs.at(lvl)[2] = calc_mean_dist(trafo_3, mpart, 4);
    }

    std::cout << "Volume Errors" << std::endl;
    for(Index lvl(0); lvl <= lvl_max; ++lvl)
    {
      std::cout << stringify(lvl).pad_front(2) << ":";

      for(std::size_t i(0); i < ns; ++i)
        std::cout << stringify_fp_sci(verrs.at(lvl)[i], 4, 12);
      if(lvl > lvl_min)
      {
        std::cout << " ||";
        for(std::size_t i(0); i < ns; ++i)
          std::cout << stringify_fp_fix(verrs.at(lvl-1)[i] / verrs.at(lvl)[i], 3, 8);
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Circumference Errors" << std::endl;
    for(Index lvl(0); lvl <= lvl_max; ++lvl)
    {
      std::cout << stringify(lvl).pad_front(2) << ":";

      for(std::size_t i(0); i < ns; ++i)
        std::cout << stringify_fp_sci(cerrs.at(lvl)[i], 4, 12);
      if(lvl > lvl_min)
      {
        std::cout << " ||";
        for(std::size_t i(0); i < ns; ++i)
          std::cout << stringify_fp_fix(cerrs.at(lvl-1)[i] / cerrs.at(lvl)[i], 3, 8);
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Mean Distance" << std::endl;
    for(Index lvl(0); lvl <= lvl_max; ++lvl)
    {
      std::cout << stringify(lvl).pad_front(2) << ":";

      for(std::size_t i(0); i < ns; ++i)
        std::cout << stringify_fp_sci(merrs.at(lvl)[i], 4, 12);
      if(lvl > lvl_min)
      {
        std::cout << " ||";
        for(std::size_t i(0); i < ns; ++i)
          std::cout << stringify_fp_fix(merrs.at(lvl-1)[i] / merrs.at(lvl)[i], 3, 8);
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

  }
} // namespace DbgIsoParam

int main(int argc, char* argv[])
{
  Runtime::initialise(argc, argv);
  DbgIsoParam::run(argc, argv);
  return Runtime::finalise();
}
