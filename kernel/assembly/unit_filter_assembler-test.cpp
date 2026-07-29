// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/analytic/lambda_function.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/geometry/export_vtk.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for UnitFilter assembly
 *
 * \author Peter Zajac
 */
template<typename DT_, typename IT_>
class UnitFilterAssemblyTest
  : public UnitTest
{
public:
  explicit UnitFilterAssemblyTest(PreferredBackend backend)
    : UnitTest("UnitFilterAssemblyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));
    const DT_ pi = Math::pi<DT_>();

    static constexpr int world_dim = 2;
    typedef Shape::Simplex<world_dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, world_dim, DT_> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceType;

    typedef LAFEM::DenseVector<DT_, IT_> VectorType;
    typedef LAFEM::UnitFilter<DT_, IT_> FilterType;

    std::stringstream ioss;

    ioss << "<FeatMeshFile version=\"1\" mesh=\"conformal:simplex:2:2\">\n";
    ioss << "  <Mesh type=\"conformal:simplex:2:2\" size=\"5 8 4\">\n";
    ioss << "    <Vertices>\n";
    ioss << "      0 0\n";
    ioss << "      1 0\n";
    ioss << "      1 1\n";
    ioss << "      0 1\n";
    ioss << "      0.5 0.5\n";
    ioss << "    </Vertices>\n";
    ioss << "    <Topology dim=\"1\">\n";
    ioss << "      0 1\n";
    ioss << "      1 2\n";
    ioss << "      2 3\n";
    ioss << "      3 0\n";
    ioss << "      0 4\n";
    ioss << "      1 4\n";
    ioss << "      2 4\n";
    ioss << "      3 4\n";
    ioss << "    </Topology>\n";
    ioss << "    <Topology dim=\"2\">\n";
    ioss << "      0 1 4\n";
    ioss << "      1 2 4\n";
    ioss << "      2 3 4\n";
    ioss << "      3 0 4\n";
    ioss << "    </Topology>\n";
    ioss << "  </Mesh>\n";
    ioss << "  <MeshPart name=\"bnd:b\" parent=\"root\" topology=\"none\" size=\"2 1\">\n";
    ioss << "    <Mapping dim=\"0\">\n";
    ioss << "      0\n";
    ioss << "      1\n";
    ioss << "    </Mapping>\n";
    ioss << "    <Mapping dim=\"1\">\n";
    ioss << "      0\n";
    ioss << "    </Mapping>\n";
    ioss << "  </MeshPart>\n";
    ioss << "  <MeshPart name=\"bnd:r\" parent=\"root\" topology=\"none\" size=\"2 1\">\n";
    ioss << "    <Mapping dim=\"0\">\n";
    ioss << "      1\n";
    ioss << "      2\n";
    ioss << "    </Mapping>\n";
    ioss << "    <Mapping dim=\"1\">\n";
    ioss << "      1\n";
    ioss << "    </Mapping>\n";
    ioss << "  </MeshPart>\n";
    ioss << "  <MeshPart name=\"bnd:t\" parent=\"root\" topology=\"none\" size=\"2 1\">\n";
    ioss << "    <Mapping dim=\"0\">\n";
    ioss << "      2\n";
    ioss << "      3\n";
    ioss << "    </Mapping>\n";
    ioss << "    <Mapping dim=\"1\">\n";
    ioss << "      2\n";
    ioss << "    </Mapping>\n";
    ioss << "  </MeshPart>\n";
    ioss << "  <MeshPart name=\"bnd:l\" parent=\"root\" topology=\"none\" size=\"2 1\">\n";
    ioss << "    <Mapping dim=\"0\">\n";
    ioss << "      3\n";
    ioss << "      0\n";
    ioss << "    </Mapping>\n";
    ioss << "    <Mapping dim=\"1\">\n";
    ioss << "      3\n";
    ioss << "    </Mapping>\n";
    ioss << "  </MeshPart>\n";
    ioss << "</FeatMeshFile>\n";


    // create a reader and read the root markup
    Geometry::MeshFileReader reader(ioss);
    reader.read_root_markup();

    // create an empty atlas and a root mesh node
    std::unique_ptr<Geometry::MeshAtlas<MeshType>> atlas = Geometry::MeshAtlas<MeshType>::make_unique();
    std::unique_ptr<Geometry::RootMeshNode<MeshType>> node = Geometry::RootMeshNode<MeshType>::make_unique(nullptr, atlas.get());

    reader.parse(*node, *atlas);

    // Refine the MeshNode so the MeshParts get refined, too.
    Index lvl_max(2);
    for(Index lvl(0); lvl <= lvl_max; ++lvl)
    {
      node = node->refine_unique();
    }

    // Trafo and space
    MeshType& mesh = *(node->get_mesh());
    TrafoType trafo(mesh);
    SpaceType space(trafo);

    // create test and reference vectors
    VectorType vec_ref1(space.get_num_dofs());
    VectorType vec_ref2(space.get_num_dofs());
    vec_ref1.format();
    vec_ref2.format();
    {
      const Index num_v = mesh.get_num_vertices();
      const Index num_e = mesh.get_num_entities(1);
      const auto& vtx = mesh.get_vertex_set();
      const auto& idx = mesh.template get_index_set<1,0>();
      Memory::TypedView<DT_> vref1 = vec_ref1.elements_view_rw();
      Memory::TypedView<DT_> vref2 = vec_ref2.elements_view_rw();

      for(Index i = 0; i < num_v; ++i)
      {
        // vertical boundary edge
        if(vtx[i][0] < DT_(0.0001))
        {
          vref1[i] = DT_(0);
          vref2[i] = DT_(4) * vtx[i][1] * (DT_(1) - vtx[i][1]); // = 4*y*(1-y)
        }
        if(vtx[i][0] > DT_(0.9999))
        {
          vref1[i] = DT_(0);
          vref2[i] = DT_(4) * vtx[i][1] * (DT_(1) - vtx[i][1]); // = 4*y*(1-y)
        }
        // horizontal boundary edge
        if(vtx[i][1] < DT_(0.0001))
        {
          vref1[i] = Math::sin(pi*vtx[i][0]); // = sin(pi*x)*cos(pi*y)
          vref2[i] = vref1[i];
        }
        if(vtx[i][1] > DT_(0.9999))
        {
          vref1[i] = -Math::sin(pi*vtx[i][0]); // = sin(pi*x)*cos(pi*y)
          vref2[i] = vref1[i];
        }
      }

      for(Index i = 0; i < num_e; ++i)
      {
        // compute edge midpoint
        Tiny::Vector<DT_, 2> pt = DT_(0.5) * (vtx[idx(i,0)] + vtx[idx(i,1)]);
        // vertical boundary edge
        if(pt[0] < DT_(0.0001))
        {
          vref1[num_v+i] = DT_(0);
          vref2[num_v+i] = DT_(4) * pt[1] * (DT_(1) - pt[1]); // = 4*y*(1-y)
        }
        if(pt[0] > DT_(0.9999))
        {
          vref1[num_v+i] = DT_(0);
          vref2[num_v+i] = DT_(4) * pt[1] * (DT_(1) - pt[1]); // = 4*y*(1-y)
        }
        // horizontal boundary edge
        if(pt[1] < DT_(0.0001))
        {
          vref1[num_v+i] = Math::sin(pi*pt[0]); // = sin(pi*x)*cos(pi*y)
          vref2[num_v+i] = vref1[num_v+i];
        }
        if(pt[1] > DT_(0.9999))
        {
          vref1[num_v+i] = -Math::sin(pi*pt[0]); // = sin(pi*x)*cos(pi*y)
          vref2[num_v+i] = vref1[num_v+i];
        }
      }
    }

    Assembly::UnitFilterAssembler<MeshType> unit_asm_1, unit_asm_2;
    unit_asm_1.add_mesh_part(*(node->find_mesh_part("bnd:b")));
    unit_asm_1.add_mesh_part(*(node->find_mesh_part("bnd:t")));
    unit_asm_1.add_mesh_part(*(node->find_mesh_part("bnd:l")));
    unit_asm_1.add_mesh_part(*(node->find_mesh_part("bnd:r")));
    unit_asm_2.add_mesh_part(*(node->find_mesh_part("bnd:l")));
    unit_asm_2.add_mesh_part(*(node->find_mesh_part("bnd:r")));

    auto func_sincos = Analytic::create_lambda_function_scalar_2d([pi](DT_ x, DT_ y){return Math::sin(pi*x)*Math::cos(pi*y);});
    auto func_parpro = Analytic::create_lambda_function_scalar_2d([](DT_, DT_ y){return DT_(4)*y*(DT_(1)-y);});

    VectorType vec_sincos, vec_parpro;
    Assembly::Interpolator::project(vec_sincos, func_sincos, space);
    Assembly::Interpolator::project(vec_parpro, func_parpro, space);

    // assemble filter, part 1
    FilterType filter_a, filter_b;
    unit_asm_1.assemble(filter_a, space, func_sincos);
    unit_asm_1.assemble(filter_b, space, vec_sincos);

    // apply filter, part 1.a
    VectorType vector1a(space.get_num_dofs());
    vector1a.format();
    filter_a.filter_sol(vector1a);

    // check filter apply, part 1.a
    TEST_CHECK_LESS_THAN(vector1a.max_rel_diff(vec_ref1), tol);

    // apply filter, part 1.b
    VectorType vector1b(space.get_num_dofs());
    vector1b.format();
    filter_b.filter_sol(vector1b);

    // check filter apply, part 1.b
    TEST_CHECK_LESS_THAN(vector1b.max_rel_diff(vec_ref1), tol);

    // use same filter and only overwrite values on vertical edges
    // the values on the horizontal edges shall not change

    // assemble filter, part 2
    unit_asm_2.assemble(filter_a, space, func_parpro);
    unit_asm_2.assemble(filter_b, space, vec_parpro);

    // apply filter, part 2.a
    VectorType vector2a(space.get_num_dofs());
    vector2a.format();
    filter_a.filter_sol(vector2a);

    // check filter apply, part 2.a
    TEST_CHECK_LESS_THAN(vector2a.max_rel_diff(vec_ref2), tol);

    // apply filter, part 2.b
    VectorType vector2b(space.get_num_dofs());
    vector2b.format();
    filter_b.filter_sol(vector2b);

    // check filter apply, part 2.b
    TEST_CHECK_LESS_THAN(vector2b.max_rel_diff(vec_ref2), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
//SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, Half, std::uint32_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterAssemblyTest, double, std::uint32_t, PreferredBackend::cuda);
#endif


/**
* \brief Test class for UnitFilterBlocked assembly
*
* \author Peter Zajac
*/
template<typename DT_, typename IT_>
class UnitFilterBlockedAssemblyTest
  : public UnitTest
{
public:
  explicit UnitFilterBlockedAssemblyTest(PreferredBackend backend)
    : UnitTest("UnitFilterBlockedAssemblyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));
    const DT_ pi = Math::pi<DT_>();

    static constexpr int dim = 2;
    typedef Shape::Simplex<dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, dim, DT_> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceType;

    typedef LAFEM::DenseVectorBlocked<DT_, IT_, dim> VectorType;
    typedef LAFEM::UnitFilterBlocked<DT_, IT_, dim> FilterType;

    std::stringstream ioss;

    ioss << "<FeatMeshFile version=\"1\" mesh=\"conformal:simplex:2:2\">\n";
    ioss << "  <Mesh type=\"conformal:simplex:2:2\" size=\"5 8 4\">\n";
    ioss << "    <Vertices>\n";
    ioss << "      0 0\n";
    ioss << "      1 0\n";
    ioss << "      1 1\n";
    ioss << "      0 1\n";
    ioss << "      0.5 0.5\n";
    ioss << "    </Vertices>\n";
    ioss << "    <Topology dim=\"1\">\n";
    ioss << "      0 1\n";
    ioss << "      1 2\n";
    ioss << "      2 3\n";
    ioss << "      3 0\n";
    ioss << "      0 4\n";
    ioss << "      1 4\n";
    ioss << "      2 4\n";
    ioss << "      3 4\n";
    ioss << "    </Topology>\n";
    ioss << "    <Topology dim=\"2\">\n";
    ioss << "      0 1 4\n";
    ioss << "      1 2 4\n";
    ioss << "      2 3 4\n";
    ioss << "      3 0 4\n";
    ioss << "    </Topology>\n";
    ioss << "  </Mesh>\n";
    ioss << "  <MeshPart name=\"bnd:b\" parent=\"root\" topology=\"none\" size=\"2 1\">\n";
    ioss << "    <Mapping dim=\"0\">\n";
    ioss << "      0\n";
    ioss << "      1\n";
    ioss << "    </Mapping>\n";
    ioss << "    <Mapping dim=\"1\">\n";
    ioss << "      0\n";
    ioss << "    </Mapping>\n";
    ioss << "  </MeshPart>\n";
    ioss << "  <MeshPart name=\"bnd:r\" parent=\"root\" topology=\"none\" size=\"2 1\">\n";
    ioss << "    <Mapping dim=\"0\">\n";
    ioss << "      1\n";
    ioss << "      2\n";
    ioss << "    </Mapping>\n";
    ioss << "    <Mapping dim=\"1\">\n";
    ioss << "      1\n";
    ioss << "    </Mapping>\n";
    ioss << "  </MeshPart>\n";
    ioss << "  <MeshPart name=\"bnd:t\" parent=\"root\" topology=\"none\" size=\"2 1\">\n";
    ioss << "    <Mapping dim=\"0\">\n";
    ioss << "      2\n";
    ioss << "      3\n";
    ioss << "    </Mapping>\n";
    ioss << "    <Mapping dim=\"1\">\n";
    ioss << "      2\n";
    ioss << "    </Mapping>\n";
    ioss << "  </MeshPart>\n";
    ioss << "  <MeshPart name=\"bnd:l\" parent=\"root\" topology=\"none\" size=\"2 1\">\n";
    ioss << "    <Mapping dim=\"0\">\n";
    ioss << "      3\n";
    ioss << "      0\n";
    ioss << "    </Mapping>\n";
    ioss << "    <Mapping dim=\"1\">\n";
    ioss << "      3\n";
    ioss << "    </Mapping>\n";
    ioss << "  </MeshPart>\n";
    ioss << "</FeatMeshFile>\n";


    // create a reader and read the root markup
    Geometry::MeshFileReader reader(ioss);
    reader.read_root_markup();

    // create an empty atlas and a root mesh node
    std::unique_ptr<Geometry::MeshAtlas<MeshType>> atlas = Geometry::MeshAtlas<MeshType>::make_unique();
    std::unique_ptr<Geometry::RootMeshNode<MeshType>> node = Geometry::RootMeshNode<MeshType>::make_unique(nullptr, atlas.get());

    reader.parse(*node, *atlas);

    // Refine the MeshNode so the MeshParts get refined, too.
    Index lvl_max(2);
    for(Index lvl(0); lvl <= lvl_max; ++lvl)
    {
      node = node->refine_unique();
    }

    // Trafo and space
    MeshType& mesh = *(node->get_mesh());
    TrafoType trafo(mesh);
    SpaceType space(trafo);

    // create test and reference vectors
    VectorType vec_ref1(space.get_num_dofs());
    VectorType vec_ref2(space.get_num_dofs());
    vec_ref1.format();
    vec_ref2.format();
    {
      const Index num_v = mesh.get_num_vertices();
      const Index num_e = mesh.get_num_entities(1);
      const auto& vtx = mesh.get_vertex_set();
      const auto& idx = mesh.template get_index_set<1,0>();
      auto vref1 = vec_ref1.elements_view_rw();
      auto vref2 = vec_ref2.elements_view_rw();

      for(Index i = 0; i < num_v; ++i)
      {
        // vertical boundary edge
        if(vtx[i][0] < DT_(0.0001))
        {
          vref1[i][1] = -Math::sin(pi*vtx[i][1]); // = -cos(pi*x)*sin(pi*y)
          vref2[i][0] = DT_(4) * vtx[i][1] * (DT_(1) - vtx[i][1]); // = 4*y*(1-y)*(1-x)
        }
        if(vtx[i][0] > DT_(0.9999))
        {
          vref1[i][1] =  Math::sin(pi*vtx[i][1]); // = -cos(pi*x)*sin(pi*y)
          vref2[i][1] = DT_(4) * vtx[i][1] * (DT_(1) - vtx[i][1]); // = 4*y*(1-y)*x
        }
        // horizontal boundary edge
        if(vtx[i][1] < DT_(0.0001))
        {
          vref1[i][0] =  Math::sin(pi*vtx[i][0]); // =  sin(pi*x)*cos(pi*y)
          vref2[i] = vref1[i];
        }
        if(vtx[i][1] > DT_(0.9999))
        {
          vref1[i][0] = -Math::sin(pi*vtx[i][0]); // =  sin(pi*x)*cos(pi*y)
          vref2[i] = vref1[i];
        }
      }

      for(Index i = 0; i < num_e; ++i)
      {
        // compute edge midpoint
        const Tiny::Vector<DT_, 2> pt = DT_(0.5) * (vtx[idx(i,0)] + vtx[idx(i,1)]);

        // vertical boundary edge
        if(pt[0] < DT_(0.0001))
        {
          vref1[num_v+i][1] = -Math::sin(pi*pt[1]); // = -cos(pi*x)*sin(pi*y)
          vref2[num_v+i][0] = DT_(4) * pt[1] * (DT_(1) - pt[1]); // = 4*y*(1-y)*(1-x)
        }
        if(pt[0] > DT_(0.9999))
        {
          vref1[num_v+i][1] =  Math::sin(pi*pt[1]); // = -cos(pi*x)*sin(pi*y)
          vref2[num_v+i][1] = DT_(4) * pt[1] * (DT_(1) - pt[1]); // = 4*y*(1-y)*x
        }
        // horizontal boundary edge
        if(pt[1] < DT_(0.0001))
        {
          vref1[num_v+i][0] =  Math::sin(pi*pt[0]); // =  sin(pi*x)*cos(pi*y)
          vref2[num_v+i] = vref1[num_v+i];
        }
        if(pt[1] > DT_(0.9999))
        {
          vref1[num_v+i][0] = -Math::sin(pi*pt[0]); // =  sin(pi*x)*cos(pi*y)
          vref2[num_v+i] = vref1[num_v+i];
        }
      }
    }

    Assembly::UnitFilterAssembler<MeshType> unit_asm_1, unit_asm_2;
    unit_asm_1.add_mesh_part(*(node->find_mesh_part("bnd:b")));
    unit_asm_1.add_mesh_part(*(node->find_mesh_part("bnd:t")));
    unit_asm_1.add_mesh_part(*(node->find_mesh_part("bnd:l")));
    unit_asm_1.add_mesh_part(*(node->find_mesh_part("bnd:r")));
    unit_asm_2.add_mesh_part(*(node->find_mesh_part("bnd:l")));
    unit_asm_2.add_mesh_part(*(node->find_mesh_part("bnd:r")));

    auto func_sincos = Analytic::create_lambda_function_vector_2d(
      [pi](DT_ x, DT_ y){return  Math::sin(pi*x)*Math::cos(pi*y);},
      [pi](DT_ x, DT_ y){return -Math::cos(pi*x)*Math::sin(pi*y);}
    );
    auto func_parpro = Analytic::create_lambda_function_vector_2d(
      [](DT_ x, DT_ y){return DT_(4)*y*(DT_(1)-y)*(DT_(1)-x);},
      [](DT_ x, DT_ y){return DT_(4)*y*(DT_(1)-y)*x;}
    );

    //auto refnode = node->refine_unique();
    //Geometry::ExportVTK<MeshType> vtk(*refnode->get_mesh());

    VectorType vec_sincos, vec_parpro;
    Assembly::Interpolator::project(vec_sincos, func_sincos, space);
    Assembly::Interpolator::project(vec_parpro, func_parpro, space);

    // assemble filter, part 1
    FilterType filter_a, filter_b;
    unit_asm_1.assemble(filter_a, space, func_sincos);
    unit_asm_1.assemble(filter_b, space, vec_sincos);

    // apply filter, part 1.a
    VectorType vector1a(space.get_num_dofs());
    vector1a.format();
    filter_a.filter_sol(vector1a);

    // check filter apply, part 1.a
    TEST_CHECK_LESS_THAN(vector1a.max_rel_diff(vec_ref1), tol);
    //vtk.add_vertex_vector("vec1a", vector1a);
    //vtk.add_vertex_vector("ref1", vec_ref1);

    // apply filter, part 1.b
    VectorType vector1b(space.get_num_dofs());
    vector1b.format();
    filter_b.filter_sol(vector1b);

    // check filter apply, part 1.b
    TEST_CHECK_LESS_THAN(vector1b.max_rel_diff(vec_ref1), tol);
    //vtk.add_vertex_vector("vec1b", vector1b);

    // use same filter and only overwrite values on vertical edges
    // the values on the horizontal edges shall not change

    // assemble filter, part 2
    unit_asm_2.assemble(filter_a, space, func_parpro);
    unit_asm_2.assemble(filter_b, space, vec_parpro);

    // apply filter, part 2.a
    VectorType vector2a(space.get_num_dofs());
    vector2a.format();
    filter_a.filter_sol(vector2a);

    // check filter apply, part 2.a
    TEST_CHECK_LESS_THAN(vector2a.max_rel_diff(vec_ref2), tol);
    //vtk.add_vertex_vector("vec2a", vector2a);
    //vtk.add_vertex_vector("ref2", vec_ref2);

    // apply filter, part 2.b
    VectorType vector2b(space.get_num_dofs());
    vector2b.format();
    filter_b.filter_sol(vector2b);

    // check filter apply, part 2.b
    TEST_CHECK_LESS_THAN(vector2b.max_rel_diff(vec_ref2), tol);
    //vtk.add_vertex_vector("vec2b", vector2b);
    //vtk.write("unit_asm_blocked");
  }
};

SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
//SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, Half, std::uint32_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(UnitFilterBlockedAssemblyTest, double, std::uint32_t, PreferredBackend::cuda);
#endif
