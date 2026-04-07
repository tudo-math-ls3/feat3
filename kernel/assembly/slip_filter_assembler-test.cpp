// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/unit_cube_patch_generator.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/trafo/standard/mapping.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for SlipFilter assembly
 *
 * Create a mesh, some MeshParts, assemble the filter and apply it to a bogus function.
 *
 * \author Jordi Paul
 *
 */
template
<
  typename DT_,
  typename IT_
>
class SlipFilterAssemblyTest
  : public UnitTest
{
public:
  explicit SlipFilterAssemblyTest(PreferredBackend backend)
    : UnitTest("SlipFilterAssemblyTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  /**
   * Runs a test in 2d.
   */
  template<template<typename> class SpaceType_>
  void run_2d() const
  {
    static constexpr int world_dim = 2;
    typedef Shape::Simplex<world_dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, world_dim, DT_> MeshType;

    typedef Tiny::Vector<DT_, world_dim> ValueType;
    typedef DenseVectorBlocked<DT_, IT_, world_dim> VectorType;

    typedef SlipFilter<DT_, IT_, world_dim> FilterType;

    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef SpaceType_<TrafoType> SpaceType;

    std::stringstream ioss;

    // Dump the content of ./data/meshes/unit-circle-tria.txt into the stream
    ioss << "<FeatMeshFile version=\"1\" mesh=\"conformal:simplex:2:2\">\n";
    ioss << "  <Info>\n";
    ioss << "   This is the unit-circle mesh consisting of four triangles.\n";
    ioss << "  </Info>\n";
    ioss << "  <Chart name=\"outer\">\n";
    ioss << "    <Circle radius=\"1\" midpoint=\"0 0\" domain=\"0 4\" />\n";
    ioss << "  </Chart>\n";
    ioss << "  <Mesh type=\"conformal:simplex:2:2\" size=\"5 8 4\">\n";
    ioss << "    <Vertices>\n";
    ioss << "      1 0\n";
    ioss << "      0 1\n";
    ioss << "      -1 0\n";
    ioss << "      0 -1\n";
    ioss << "      0 0\n";
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
    ioss << "  <MeshPart name=\"outer\" parent=\"root\" chart=\"outer\" topology=\"full\" size=\"5 4\">\n";
    ioss << "    <Mapping dim=\"0\">\n";
    ioss << "      0\n";
    ioss << "      1\n";
    ioss << "      2\n";
    ioss << "      3\n";
    ioss << "      0\n";
    ioss << "    </Mapping>\n";
    ioss << "    <Mapping dim=\"1\">\n";
    ioss << "      0\n";
    ioss << "      1\n";
    ioss << "      2\n";
    ioss << "      3\n";
    ioss << "    </Mapping>\n";
    ioss << "    <Topology dim=\"1\">\n";
    ioss << "      0 1\n";
    ioss << "      1 2\n";
    ioss << "      2 3\n";
    ioss << "      3 4\n";
    ioss << "    </Topology>\n";
    ioss << "    <Attribute name=\"param\" dim=\"1\">\n";
    ioss << "      0\n";
    ioss << "      1\n";
    ioss << "      2\n";
    ioss << "      3\n";
    ioss << "      4\n";
    ioss << "    </Attribute>\n";
    ioss << "  </MeshPart>\n";
    ioss << "</FeatMeshFile>\n";

    // create a reader and read the root markup
    Geometry::MeshFileReader reader(ioss);
    reader.read_root_markup();

    // create an empty atlas and a root mesh node
    std::unique_ptr<Geometry::MeshAtlas<MeshType>> atlas = Geometry::MeshAtlas<MeshType>::make_unique();
    std::unique_ptr<Geometry::RootMeshNode<MeshType>> node = Geometry::RootMeshNode<MeshType>::make_unique(nullptr, atlas.get());

    reader.parse(*node, *atlas);

    node->adapt();

    // Refine the MeshNode so the MeshParts get refined, too.
    Index lvl_max(0);
    for(Index lvl(0); lvl <= lvl_max; ++lvl)
    {
      node = node->refine_unique();
    }

    // Trafo and space
    TrafoType my_trafo(*(node->get_mesh()));
    SpaceType my_space(my_trafo);

    // Create the analytic function component wise
    LAFEM::DenseVector<DT_, IT_>comp0(my_space.get_num_dofs());
    Analytic::Common::CosineWaveFunction<2> func0;
    Assembly::Interpolator::project(comp0, func0, my_space);

    LAFEM::DenseVector<DT_, IT_>comp1(my_space.get_num_dofs());
    Analytic::Common::SineBubbleFunction<2> func1;
    Assembly::Interpolator::project(comp1, func1, my_space);

    // Paste the components into a blocked vector and keep a copy for checking
    VectorType vec_org(my_space.get_num_dofs(), DT_(0));
    {
      Memory::TypedView<ValueType> vorg = vec_org.elements_view_w();
      const Memory::TypedView<DT_> v0 = comp0.elements_view_r();
      const Memory::TypedView<DT_> v1 = comp1.elements_view_r();
      for(Index i(0); i < my_space.get_num_dofs(); ++i)
      {
        vorg[i][0] = v0[i];
        vorg[i][1] = v1[i];
      }
    }
    VectorType vec(my_space.get_num_dofs(), DT_(0));
    vec.clone(vec_org);

    // The assembler
    Assembly::SlipFilterAssembler<TrafoType> slip_filter_assembler(my_trafo);

    FilterType my_filter;
    slip_filter_assembler.add_mesh_part(*(node->find_mesh_part("outer")));
    slip_filter_assembler.assemble(my_filter, my_space);

    // Apply the filter
    my_filter.filter_sol(vec);

    // Check results
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));

    // First check all filtered entries if they are really orthogonal to the normal vector saved in the filter
    {
      const Memory::TypedView<IT_> idx = my_filter.get_filter_vector().indices_view_r();
      const Memory::TypedView<ValueType> nu = my_filter.get_filter_vector().elements_view_r();
      const Memory::TypedView<ValueType> vorg = vec_org.elements_view_r();
      Memory::TypedView<ValueType> vv = vec.elements_view_rw();
      for(Index i(0); i < my_filter.num_nzes(); ++i)
      {
        Index j(idx[i]);
        TEST_CHECK_EQUAL_WITHIN_EPS(Tiny::dot(vv[j], nu(i)), DT_(0), tol);
        // If this was ok, replace with the original value so we can check the whole vector without bothering with
        // identifying the filtered values below
        vv[j] = vorg[j];
      }
    }

    // Now check all values in the vector to make sure the filter did not touch the rest
    TEST_CHECK_LESS_THAN(vec.max_rel_diff(vec_org), tol);
  }

  /**
   * Runs the test in 3d. Creates a unit cube Hypercube<3> mesh, adds 3 of its faces to the filter and filters
   * an interpolation of an analytic function.
   */
  template<template<typename> class SpaceType_>
  void run_3d() const
  {
    static constexpr int world_dim = 3;
    typedef Shape::Hypercube<world_dim> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType, world_dim, DT_> MeshType;

    typedef Tiny::Vector<DT_, world_dim> ValueType;
    typedef DenseVectorBlocked<DT_, IT_, world_dim> VectorType;

    typedef SlipFilter<DT_, IT_, world_dim> FilterType;

    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    // The SlipFilter is implemented for Lagrange 1/2 only
    typedef SpaceType_<TrafoType> SpaceType;

    // This is for creating the mesh and its MeshParts
    std::unique_ptr<Geometry::RootMeshNode<MeshType>> node;
    std::vector<int> ranks;
    Geometry::UnitCubePatchGenerator<MeshType>::create_unique(0, 1, node, ranks);

    // Refine the MeshNode so the MeshParts get refined, too.
    Index lvl_max(2);
    for(Index lvl(0); lvl <= lvl_max; ++lvl)
    {
      node = node->refine_unique();
    }

    // Trafo and space
    TrafoType my_trafo(*(node->get_mesh()));
    SpaceType my_space(my_trafo);

    // Create the analytic function component wise
    LAFEM::DenseVector<DT_, IT_>comp0(my_space.get_num_dofs());
    Analytic::Common::ConstantFunction<3> func0(-DT_(0.5));
    Assembly::Interpolator::project(comp0, func0, my_space);

    LAFEM::DenseVector<DT_, IT_>comp1(my_space.get_num_dofs());
    Analytic::Common::SineBubbleFunction<3> func1;
    Assembly::Interpolator::project(comp1, func1, my_space);

    LAFEM::DenseVector<DT_, IT_>comp2(my_space.get_num_dofs());
    Analytic::Common::CosineWaveFunction<3> func2;
    Assembly::Interpolator::project(comp2, func2, my_space);

    // Paste the components into a blocked vector and keep a copy for checking
    VectorType vec_org(my_space.get_num_dofs(), DT_(0));
    {
      Memory::TypedView<ValueType> vorg = vec_org.elements_view_w();
      const Memory::TypedView<DT_> v0 = comp0.elements_view_r();
      const Memory::TypedView<DT_> v1 = comp1.elements_view_r();
      const Memory::TypedView<DT_> v2 = comp2.elements_view_r();
      for(Index i(0); i < my_space.get_num_dofs(); ++i)
      {
        vorg[i][0] = v0[i];
        vorg[i][1] = v1[i];
        vorg[i][2] = v2[i];
      }
    }
    VectorType vec(my_space.get_num_dofs(), DT_(0));
    vec.clone(vec_org);

    // The assembler
    Assembly::SlipFilterAssembler<TrafoType>slip_filter_assembler(my_trafo) ;

    // Create the filter, add 3 of the 6 MeshParts and assemble the filter
    FilterType my_filter;
    slip_filter_assembler.add_mesh_part(*(node->find_mesh_part("bnd:0")));
    slip_filter_assembler.add_mesh_part(*(node->find_mesh_part("bnd:2")));
    slip_filter_assembler.add_mesh_part(*(node->find_mesh_part("bnd:4")));
    slip_filter_assembler.assemble(my_filter, my_space);

    // Apply the filter
    my_filter.filter_sol(vec);

    // Check results
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));

    // First check all filtered entries if they are really orthogonal to the normal vector saved in the filter
    {
      const Memory::TypedView<IT_> idx = my_filter.get_filter_vector().indices_view_r();
      const Memory::TypedView<ValueType> nu = my_filter.get_filter_vector().elements_view_r();
      const Memory::TypedView<ValueType> vorg = vec_org.elements_view_r();
      Memory::TypedView<ValueType> vv = vec.elements_view_rw();
      for(Index i(0); i < my_filter.num_nzes(); ++i)
      {
        Index j(idx[i]);
        TEST_CHECK_EQUAL_WITHIN_EPS(Tiny::dot(vv[j], nu(i)), DT_(0), tol);
        // If this was ok, replace with the original value so we can check the whole vector without bothering with
        // identifying the filtered values below
        vv[j] = vorg[j];
      }
    }

    // Now check all values in the vector to make sure the filter did not touch the rest
    TEST_CHECK_LESS_THAN(vec.max_rel_diff(vec_org), tol);
  }

  void run() const override
  {
    // The SlipFilter is implemented for Lagrange 1/2 only
    run_2d<Space::Lagrange1::Element>();
    run_2d<Space::Lagrange2::Element>();
    run_3d<Space::Lagrange1::Element>();
    run_3d<Space::Lagrange2::Element>();
  }
};

SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, double, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, double, std::uint64_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
//SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, Half, std::uint32_t, PreferredBackend::generic);
//SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(SlipFilterAssemblyTest, double, std::uint32_t, PreferredBackend::cuda);
#endif
