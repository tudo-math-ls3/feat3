// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
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
 * \brief Test class for SlipFilter class template
 *
 * Primitive tests based on add()-ing values to the filter
 *
 * \author Jordi Paul
 *
 */
template
<
  typename MemType_,
  typename DT_,
  typename IT_,
  int BlockSize_
>
class SlipFilterVectorTest
: public FullTaggedTest<MemType_, DT_, IT_>
{
  typedef Tiny::Vector<DT_, BlockSize_> ValueType;
  typedef DenseVectorBlocked<MemType_, DT_, IT_, BlockSize_> VectorType;
  typedef DenseVectorBlocked<MemType_, IT_, IT_, BlockSize_> IVectorType;
  typedef SlipFilter<MemType_, DT_, IT_, BlockSize_> FilterType;

  public:
  SlipFilterVectorTest()
    : FullTaggedTest<MemType_, DT_, IT_>("SlipFilterVectorTest")
    {
    }

  virtual ~SlipFilterVectorTest()
  {
  }

  virtual void run() const override
  {
    const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.75));

    const IT_ nn(100);
    FilterType my_filter(nn, nn);

    IT_ jj[8];

    for(IT_ i(0); i < IT_(8); ++i)
      jj[i] = i*(IT_(3) + i);

    for(IT_ i(0); i < IT_(8); ++i)
    {
      IT_ j = jj[i];

      ValueType tmp(Math::sqrt(DT_(j+1)));
      tmp(0) = tmp(0)*DT_(Math::pow(-DT_(0.5), DT_(i)));
      tmp(BlockSize_-1) = -DT_(i);

      //tmp.normalise();

      my_filter.add(j,tmp);
    }

    VectorType my_vector(nn, DT_(-2));

    IT_ j(0);
    // 0: Set element to 0
    my_vector(j, ValueType(0));
    // 1: Set element to value of the filter
    j = jj[1];
    ValueType tmp1(my_filter.get_filter_vector()(j));
    my_vector(j, tmp1);
    // 2: Set element to negative value of the filter
    j = jj[2];
    ValueType tmp2(my_filter.get_filter_vector()(j));
    for(int d(0); d < BlockSize_; ++d)
    {
      tmp2(d) = -tmp2(d);
    }
    my_vector(j, tmp2);
    // 3: Set element to first unit vector
    j = jj[3];
    ValueType tmp3(DT_(0));
    tmp3(0) = DT_(1);
    my_vector(j, tmp3);
    // 4: Set element to -(last unit vector)
    j = jj[4];
    ValueType tmp4(DT_(0));
    tmp4(BlockSize_-1) = -DT_(1);
    my_vector(j, tmp4);
    // 5: Set element to twice the negative value of the filter
    j = jj[5];
    ValueType tmp5(my_filter.get_filter_vector()(j));
    for(int d(0); d < BlockSize_; ++d)
    {
      tmp5(d) = -DT_(2)*tmp5(d);
    }
    my_vector(j, tmp5);

    // Filter vector
    my_filter.filter_def(my_vector);

    // Check results
    for(IT_ i(0); i < nn; ++i)
    {
      // Check if we have a filtered entry
      bool filtered(false);
      for(int k(0); k < 8; ++k)
      {
        if(i == jj[k])
        {
          filtered = true;
        }
      }

      if(filtered)
      {
        ValueType nu(my_filter.get_filter_vector()(i));
        DT_ res(Tiny::dot(nu, my_vector(i)));
        TEST_CHECK_EQUAL_WITHIN_EPS(res, DT_(0), tol);
      }
      else
      {
        for(int d(0); d < BlockSize_; ++d)
          TEST_CHECK_EQUAL_WITHIN_EPS(my_vector(i)(d), -DT_(2), tol);
      }
    }
  }
};

SlipFilterVectorTest<Mem::Main, float, Index, 2> component_filter_test_generic_fi_2;
SlipFilterVectorTest<Mem::Main, double, unsigned int, 2> component_filter_test_generic_di_2;
SlipFilterVectorTest<Mem::Main, float, unsigned int, 3> component_filter_test_generic_fi_3;
SlipFilterVectorTest<Mem::Main, double, Index, 3> component_filter_test_generic_di_3;
#ifdef FEAT_HAVE_QUADMATH
SlipFilterVectorTest<Mem::Main, __float128, unsigned int, 2> component_filter_test_generic_q_2;
SlipFilterVectorTest<Mem::Main, __float128, Index, 3> component_filter_test_generic_q_3;
#endif
#ifdef FEAT_HAVE_CUDA
SlipFilterVectorTest<Mem::CUDA, float, unsigned int, 2> component_filter_test_cuda_fi_2;
SlipFilterVectorTest<Mem::CUDA, float, Index, 3> component_filter_test_cuda_fi_3;
SlipFilterVectorTest<Mem::CUDA, double, Index, 2> component_filter_test_cuda_di_2;
SlipFilterVectorTest<Mem::CUDA, double, unsigned int, 3> component_filter_test_cuda_di_3;
#endif

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
  typename MemType_,
  typename DT_,
  typename IT_
>
class SlipFilterAssemblyTest
: public FullTaggedTest<MemType_, DT_, IT_>
{

  public:
    SlipFilterAssemblyTest()
      : FullTaggedTest<MemType_, DT_, IT_>("SlipFilterAssemblyTest")
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
      typedef DenseVectorBlocked<MemType_, DT_, IT_, world_dim> VectorType;
      typedef DenseVectorBlocked<Mem::Main, DT_, IT_, world_dim> CheckVectorType;

      typedef SlipFilter<MemType_, DT_, IT_, world_dim> FilterType;
      typedef SlipFilter<Mem::Main, DT_, IT_, world_dim> CheckFilterType;

      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      typedef SpaceType_<TrafoType> SpaceType;

      std::stringstream ioss;

      // Dump the content of ./data/meshes/unit-circle-tria.txt into the stream
      ioss << "<FeatMeshFile version=\"1\" mesh=\"conformal:simplex:2:2\">" << std::endl;
      ioss << "  <Info>" << std::endl;
      ioss << "   This is the unit-circle mesh consisting of four triangles." << std::endl;
      ioss << "  </Info>" << std::endl;
      ioss << "  <Chart name=\"outer\">" << std::endl;
      ioss << "    <Circle radius=\"1\" midpoint=\"0 0\" domain=\"0 4\" />" << std::endl;
      ioss << "  </Chart>" << std::endl;
      ioss << "  <Mesh type=\"conformal:simplex:2:2\" size=\"5 8 4\">" << std::endl;
      ioss << "    <Vertices>" << std::endl;
      ioss << "      1 0" << std::endl;
      ioss << "      0 1" << std::endl;
      ioss << "      -1 0" << std::endl;
      ioss << "      0 -1" << std::endl;
      ioss << "      0 0" << std::endl;
      ioss << "    </Vertices>" << std::endl;
      ioss << "    <Topology dim=\"1\">" << std::endl;
      ioss << "      0 1" << std::endl;
      ioss << "      1 2" << std::endl;
      ioss << "      2 3" << std::endl;
      ioss << "      3 0" << std::endl;
      ioss << "      0 4" << std::endl;
      ioss << "      1 4" << std::endl;
      ioss << "      2 4" << std::endl;
      ioss << "      3 4" << std::endl;
      ioss << "    </Topology>" << std::endl;
      ioss << "    <Topology dim=\"2\">" << std::endl;
      ioss << "      0 1 4" << std::endl;
      ioss << "      1 2 4" << std::endl;
      ioss << "      2 3 4" << std::endl;
      ioss << "      3 0 4" << std::endl;
      ioss << "    </Topology>" << std::endl;
      ioss << "  </Mesh>" << std::endl;
      ioss << "  <MeshPart name=\"outer\" parent=\"root\" chart=\"outer\" topology=\"full\" size=\"5 4\">" << std::endl;
      ioss << "    <Mapping dim=\"0\">" << std::endl;
      ioss << "      0" << std::endl;
      ioss << "      1" << std::endl;
      ioss << "      2" << std::endl;
      ioss << "      3" << std::endl;
      ioss << "      0" << std::endl;
      ioss << "    </Mapping>" << std::endl;
      ioss << "    <Mapping dim=\"1\">" << std::endl;
      ioss << "      0" << std::endl;
      ioss << "      1" << std::endl;
      ioss << "      2" << std::endl;
      ioss << "      3" << std::endl;
      ioss << "    </Mapping>" << std::endl;
      ioss << "    <Topology dim=\"1\">" << std::endl;
      ioss << "      0 1" << std::endl;
      ioss << "      1 2" << std::endl;
      ioss << "      2 3" << std::endl;
      ioss << "      3 4" << std::endl;
      ioss << "    </Topology>" << std::endl;
      ioss << "    <Attribute name=\"param\" dim=\"1\">" << std::endl;
      ioss << "      0" << std::endl;
      ioss << "      1" << std::endl;
      ioss << "      2" << std::endl;
      ioss << "      3" << std::endl;
      ioss << "      4" << std::endl;
      ioss << "    </Attribute>" << std::endl;
      ioss << "  </MeshPart>" << std::endl;
      ioss << "</FeatMeshFile>" << std::endl;

      // create a reader and read the root markup
      Geometry::MeshFileReader reader(ioss);
      reader.read_root_markup();

      // create an empty atlas and a root mesh node
      Geometry::MeshAtlas<MeshType>* atlas = new Geometry::MeshAtlas<MeshType>();
      Geometry::RootMeshNode<MeshType>* node = new Geometry::RootMeshNode<MeshType>(nullptr, atlas);

      reader.parse(*node, *atlas);

      node->adapt();

      // Refine the MeshNode so the MeshParts get refined, too.
      Index lvl_max(0);
      for(Index lvl(0); lvl <= lvl_max; ++lvl)
      {
        auto* old = node;
        node = old->refine();
        delete old;
      }

      // Trafo and space
      TrafoType my_trafo(*(node->get_mesh()));
      SpaceType my_space(my_trafo);

      // Create the analytic function component wise
      LAFEM::DenseVector<MemType_, DT_, IT_>comp0(my_space.get_num_dofs());
      Analytic::Common::CosineWaveFunction<2> func0;
      Assembly::Interpolator::project(comp0, func0, my_space);

      LAFEM::DenseVector<MemType_, DT_, IT_>comp1(my_space.get_num_dofs());
      Analytic::Common::SineBubbleFunction<2> func1;
      Assembly::Interpolator::project(comp1, func1, my_space);

      // Paste the components into a blocked vector and keep a copy for checking
      CheckVectorType vec_org(my_space.get_num_dofs(), DT_(0));
      for(Index i(0); i < my_space.get_num_dofs(); ++i)
      {
        ValueType tmp(DT_(0));
        tmp(0) = comp0(i);
        tmp(1) = comp1(i);
        vec_org(i, tmp);
      }
      VectorType vec(my_space.get_num_dofs(), DT_(0));
      vec.clone(vec_org);

      // The assembler
      Assembly::SlipFilterAssembler<MeshType>slip_filter_assembler(my_trafo.get_mesh());

      FilterType my_filter;
      slip_filter_assembler.add_mesh_part(*(node->find_mesh_part("outer")));
      slip_filter_assembler.assemble(my_filter, my_space);

      // Apply the filter
      my_filter.filter_sol(vec);

      // Check results
      const DT_ tol = Math::pow(Math::eps<DT_>(), DT_(0.9));

      // Download the filtered vector if necessary
      CheckVectorType check_vec;
      check_vec.convert(vec);

      // Download the filter if necessary
      CheckFilterType check_filter;
      check_filter.convert(my_filter);

      // First check all filtered entries if they are really orthogonal to the normal vector saved in the filter
      for(Index i(0); i < check_filter.used_elements(); ++i)
      {
        Index j(check_filter.get_indices()[i]);
        //TEST_CHECK_EQUAL_WITHIN_EPS(check_filter.get_filter_vector()(j).norm_euclid(), DT_(1), tol);
        TEST_CHECK_EQUAL_WITHIN_EPS(Tiny::dot(check_vec(j),check_filter.get_filter_vector()(j)), DT_(0), tol);
        // If this was ok, replace with the original value so we can check the whole vector without bothering with
        // identifying the filtered values below
        check_vec(j, vec_org(j));
      }

      // Now check all values in the vector to make sure the filter did not touch the rest
      for(Index i(0); i < my_space.get_num_dofs(); ++i)
      {
        for(int d(0); d < world_dim; ++d)
          TEST_CHECK_EQUAL(check_vec(i)(d), vec_org(i)(d));
      }

      // Clean up
      delete node;
      delete atlas;
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
      typedef DenseVectorBlocked<MemType_, DT_, IT_, world_dim> VectorType;

      typedef DenseVectorBlocked<Mem::Main, DT_, IT_, world_dim> CheckVectorType;

      typedef SlipFilter<MemType_, DT_, IT_, world_dim> FilterType;
      typedef SlipFilter<Mem::Main, DT_, IT_, world_dim> CheckFilterType;

      typedef Trafo::Standard::Mapping<MeshType> TrafoType;
      // The SlipFilter is implemented for Lagrange 1/2 only
      typedef SpaceType_<TrafoType> SpaceType;

      // This is for creating the mesh and its MeshParts
      Geometry::RootMeshNode<MeshType>* node = nullptr;
      std::vector<Index> ranks;
      std::vector<Index> ctags;
      Geometry::UnitCubePatchGenerator<MeshType>::create(0, 1, node, ranks, ctags);

      // Refine the MeshNode so the MeshParts get refined, too.
      Index lvl_max(2);
      for(Index lvl(0); lvl <= lvl_max; ++lvl)
      {
        auto* old = node;
        node = old->refine();
        delete old;
      }

      // Trafo and space
      TrafoType my_trafo(*(node->get_mesh()));
      SpaceType my_space(my_trafo);

      // Create the analytic function component wise
      LAFEM::DenseVector<MemType_, DT_, IT_>comp0(my_space.get_num_dofs());
      Analytic::Common::ConstantFunction<3> func0(-DT_(0.5));
      Assembly::Interpolator::project(comp0, func0, my_space);

      LAFEM::DenseVector<MemType_, DT_, IT_>comp1(my_space.get_num_dofs());
      Analytic::Common::SineBubbleFunction<3> func1;
      Assembly::Interpolator::project(comp1, func1, my_space);

      LAFEM::DenseVector<MemType_, DT_, IT_>comp2(my_space.get_num_dofs());
      Analytic::Common::CosineWaveFunction<3> func2;
      Assembly::Interpolator::project(comp2, func2, my_space);

      // Paste the components into a blocked vector and keep a copy for checking
      VectorType vec(my_space.get_num_dofs(), DT_(0));
      CheckVectorType vec_org(my_space.get_num_dofs(), DT_(0));
      for(Index i(0); i < my_space.get_num_dofs(); ++i)
      {
        ValueType tmp(DT_(0));
        tmp(0) = comp0(i);
        tmp(1) = comp1(i);
        tmp(2) = comp2(i);
        vec(i, tmp);
      }
      vec_org.clone(vec);


      // The assembler
      Assembly::SlipFilterAssembler<MeshType>slip_filter_assembler(my_trafo.get_mesh()) ;

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

      // Download the filter if necessary
      CheckFilterType check_filter;
      check_filter.convert(my_filter);

      // Download the filtered vector if necessary
      CheckVectorType check_vec;
      check_vec.convert(vec);

      // First check all filtered entries if they are really orthogonal to the normal vector saved in the filter
      for(Index i(0); i < check_filter.used_elements(); ++i)
      {
        Index j(check_filter.get_indices()[i]);
        //TEST_CHECK_EQUAL_WITHIN_EPS(check_filter.get_filter_vector()(j).norm_euclid(), DT_(1), tol);
        TEST_CHECK_EQUAL_WITHIN_EPS(Tiny::dot(check_vec(j),check_filter.get_filter_vector()(j)), DT_(0), tol);
        // If this was ok, replace with the original value so we can check the whole vector without bothering with
        // identifying the filtered values below
        check_vec(j, vec_org(j));
      }

      // Now check all values in the vector to make sure the filter did not touch the rest
      for(Index i(0); i < my_space.get_num_dofs(); ++i)
      {
        for(int d(0); d < world_dim; ++d)
          TEST_CHECK_EQUAL(check_vec(i)(d), vec_org(i)(d));
      }

      // Clean up
      delete node;
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

SlipFilterAssemblyTest<Mem::Main, float, Index> sfat_f;
SlipFilterAssemblyTest<Mem::Main, double, unsigned int> sfat_d;
#ifdef FEAT_HAVE_QUADMATH
SlipFilterAssemblyTest<Mem::Main, __float128, Index> sfat_q;
#endif
#ifdef FEAT_HAVE_CUDA
SlipFilterAssemblyTest<Mem::CUDA, float, unsigned int> sfat_f_cuda;
SlipFilterAssemblyTest<Mem::CUDA, double, Index> sfat_d_cuda;
#endif
