// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/isoparam/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/lagrange3/element.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/assembly/common_functionals.hpp>
#include <kernel/assembly/linear_functional_assembler.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/solver/umfpack.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/jacobi_precond.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/runtime.hpp>

using namespace FEAT;
using namespace FEAT::Geometry;

namespace DbgIsoParam
{
#ifdef FEAT_HAVE_QUADMATH
  typedef __float128 DT;
#else
  typedef double DT;
#endif

  typedef Shape::Quadrilateral ShapeType2D;
  typedef ConformalMesh<ShapeType2D, 2, DT> MeshType2D;
  typedef MeshPart<MeshType2D> MeshPartType2D;
  typedef MeshAtlas<MeshType2D> AtlasType2D;
  typedef RootMeshNode<MeshType2D> MeshNodeType2D;

  typedef Shape::Hexahedron ShapeType3D;
  typedef ConformalMesh<ShapeType3D, 3, DT> MeshType3D;
  typedef MeshPart<MeshType3D> MeshPartType3D;
  typedef MeshAtlas<MeshType3D> AtlasType3D;
  typedef RootMeshNode<MeshType3D> MeshNodeType3D;

  typedef LAFEM::SparseMatrixCSR<DT, Index> MatrixType;
  typedef LAFEM::DenseVector<DT, Index> VectorType;
  typedef LAFEM::UnitFilter<DT, Index> FilterType;

  static constexpr int nt = 4;
  static constexpr int ns = 3;

  typedef Tiny::Tensor3<DT, nt, ns, 3> ErrorTensor;
  typedef Tiny::Matrix<DT, ns, 3> ErrorMatrix;
  typedef Tiny::Vector<DT, 3> ErrorVector;

  /*template<typename T_, int ai = 2>
  struct StaticSolFunction
  {
    static_assert(ai >= 2, "ai must be >= 2");
    static constexpr T_ a = T_(ai);

    static T_ eval(T_ x, T_ y)
    {
      return Math::sqrt(T_(2) - x*x - y*y) - Math::sqrt(a - T_(1));
    }

    static T_ der_x(T_ x, T_ y)
    {
      return  -x / Math::sqrt(a - x*x - y*y);
    }

    static T_ der_y(T_ x, T_ y)
    {
      return  -y / Math::sqrt(a - x*x - y*y);
    }

    static T_ der_xx(T_ x, T_ y)
    {
      return  (y*y - a) / Math::pow(a - x*x - y*y, T_(1.5));
    }

    static T_ der_yy(T_ x, T_ y)
    {
      return  (x*x - a) / Math::pow(a - x*x - y*y, T_(1.5));
    }

    static T_ der_xy(T_ x, T_ y)
    {
      return  (-x*y) / Math::pow(a - x*x - y*y, T_(1.5));
    }

    static T_ der_yx(T_ x, T_ y)
    {
      return  (-x*y) / Math::pow(a - x*x - y*y, T_(1.5));
    }

    ///////////////////////////////////////////////////////////////////////////////////////////////

    static T_ eval(T_ x, T_ y, T_ z)
    {
      return Math::sqrt(a - x*x - y*y - z*z) - Math::sqrt(a - T_(1));
    }

    static T_ der_x(T_ x, T_ y, T_ z)
    {
      return  -x / Math::sqrt(a - x*x - y*y - z*z);
    }

    static T_ der_y(T_ x, T_ y, T_ z)
    {
      return  -y / Math::sqrt(a - x*x - y*y - z*z);
    }

    static T_ der_z(T_ x, T_ y, T_ z)
    {
      return  -z / Math::sqrt(a - x*x - y*y - z*z);
    }

    static T_ der_xx(T_ x, T_ y, T_ z)
    {
      return  (y*y + z*z - a) / Math::pow(a - x*x - y*y - z*z, T_(1.5));
    }

    static T_ der_yy(T_ x, T_ y, T_ z)
    {
      return  (x*x + z*z - a) / Math::pow(a - x*x - y*y - z*z, T_(1.5));
    }

    static T_ der_zz(T_ x, T_ y, T_ z)
    {
      return  (x*x + y*y - a) / Math::pow(a - x*x - y*y - z*z, T_(1.5));
    }

    static T_ der_xy(T_ x, T_ y, T_ z)
    {
      return  (-x*y) / Math::pow(a - x*x - y*y - z*z, T_(1.5));
    }

    static T_ der_xz(T_ x, T_ y, T_ z)
    {
      return  (-x*z) / Math::pow(a - x*x - y*y - z*z, T_(1.5));
    }

    static T_ der_yx(T_ x, T_ y, T_ z)
    {
      return  (-x*y) / Math::pow(a - x*x - y*y - z*z, T_(1.5));
    }

    static T_ der_yz(T_ x, T_ y, T_ z)
    {
      return  (-y*z) / Math::pow(a - x*x - y*y - z*z, T_(1.5));
    }

    static T_ der_zx(T_ x, T_ y, T_ z)
    {
      return  (-x*z) / Math::pow(a - x*x - y*y - z*z, T_(1.5));
    }

    static T_ der_zy(T_ x, T_ y, T_ z)
    {
      return  (-y*z) / Math::pow(a - x*x - y*y - z*z, T_(1.5));
    }
  };*/

  std::unique_ptr<AtlasType2D> make_atlas_2d()
  {
    auto atlas = Geometry::MeshAtlas<MeshType2D>::make_unique();
    atlas->add_mesh_chart("chart",
      std::unique_ptr<Geometry::Atlas::Circle<MeshType2D>>(new Geometry::Atlas::Circle<MeshType2D>(0.0, 0.0, 1.0)));
    return atlas;
  }

  std::unique_ptr<AtlasType3D> make_atlas_3d()
  {
    auto atlas = Geometry::MeshAtlas<MeshType3D>::make_unique();
    atlas->add_mesh_chart("chart",
      std::unique_ptr<Geometry::Atlas::Sphere<MeshType3D>>(new Geometry::Atlas::Sphere<MeshType3D>(0.0, 0.0, 0.0, 1.0)));
    return atlas;
  }


  std::unique_ptr<MeshNodeType2D> make_mesh_node(AtlasType2D* atlas)
  {
    Geometry::UnitStarCubeFactory<MeshType2D> factory;
    auto mesh = factory.make_unique();

    // move coords from [0,1] to [-0.7,+0.7]
    {
      auto& vtx = mesh->get_vertex_set();
      for(Index i(0); i < vtx.get_num_vertices(); ++i)
      {
        (vtx[i][0] *= 1.4) -= 0.7;
        (vtx[i][1] *= 1.4) -= 0.7;
      }
    }

    auto node = Geometry::RootMeshNode<MeshType2D>::make_unique(std::move(mesh), atlas);

    Geometry::BoundaryFactory<MeshType2D> bnd_factory(*node->get_mesh());
    node->add_mesh_part("bnd", bnd_factory.make_unique(), "chart", atlas->find_mesh_chart("chart"));

    return node;
  }

  std::unique_ptr<MeshNodeType3D> make_mesh_node(AtlasType3D* atlas)
  {
    Geometry::UnitCubeFactory<MeshType3D> factory_3d;
    auto mesh = factory_3d.make_unique();

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

    auto node = Geometry::RootMeshNode<MeshType3D>::make_unique(std::move(mesh), atlas);

    Geometry::BoundaryFactory<MeshType3D> bnd_factory(*node->get_mesh());
    node->add_mesh_part("bnd", bnd_factory.make_unique(), "chart", atlas->find_mesh_chart("chart"));

    node->adapt();

    return node;
  }

  //template<typename T_>
  //using FooFunction = StaticSolFunction<T_, 2>;

  template<typename Space_, typename MeshPart_>
  void run_space(ErrorVector& err, Space_& space, const MeshPart_& bnd, const String& name)
  {
    typedef typename Space_::MeshType MeshType;
    Analytic::Common::SphereCapFunction<MeshType::shape_dim, 2> sol_function;
    //Analytic::StaticWrapperFunction<MeshType::shape_dim, FooFunction, true, true, true> sol_function;

    Cubature::DynamicFactory cubature("gauss-legendre:5");

    std::cout << String(80, '#') << "\n";
    std::cout << name << "\n";
    //VectorType vec_sol;
    //Assembly::Interpolator::project(vec_sol, sol_function, space);

    MatrixType matrix;
    Assembly::SymbolicAssembler::assemble_matrix_std1(matrix, space);
    matrix.format();

    Assembly::Common::LaplaceOperator oper;
    Assembly::BilinearOperatorAssembler::assemble_matrix1(matrix, oper, space, cubature);

    VectorType vec_sol = matrix.create_vector_l();
    VectorType vec_rhs = matrix.create_vector_l();

    vec_sol.format();
    vec_rhs.format();

    Assembly::Common::LaplaceFunctional<decltype(sol_function)> rhs_force(sol_function);
    Assembly::LinearFunctionalAssembler::assemble_vector(vec_rhs, rhs_force, space, cubature);

    FilterType filter;

    Assembly::UnitFilterAssembler<typename Space_::MeshType> unit_asm;
    unit_asm.add_mesh_part(bnd);
    unit_asm.assemble(filter, space);

    filter.filter_mat(matrix);
    filter.filter_rhs(vec_rhs);
    filter.filter_sol(vec_sol);

    auto jacobi = Solver::new_jacobi_precond(matrix, filter);
    auto solver = Solver::new_pcg(matrix, filter, jacobi);
    solver->set_tol_rel(Math::sqrt(Math::eps<DT>()));
    solver->set_max_iter(50000);
    solver->set_plot_mode(Solver::PlotMode::summary);

    solver->init();
    solver->correct(vec_sol, vec_rhs);
    solver->done();

    // compute errors
    static constexpr int max_err = (Space_::local_degree >= 2) ? 2 : 1;
    auto e = Assembly::ScalarErrorComputer<max_err>::compute(vec_sol, sol_function, space, cubature);

    err[0] = e.norm_h0;
    err[1] = e.norm_h1;
    err[2] = e.norm_h2;
  }

  template<typename Trafo_, typename MeshPart_>
  void run_trafo(ErrorMatrix& err, Trafo_& trafo, const MeshPart_& bnd, const String& name)
  {
    typedef Space::Lagrange1::Element<Trafo_> SpaceQ1;
    typedef Space::Lagrange2::Element<Trafo_> SpaceQ2;
    typedef Space::Lagrange3::Element<Trafo_> SpaceQ3;

    {
      SpaceQ1 space_q1(trafo);
      run_space(err[0], space_q1, bnd, String("Q1:") + name);
    }
    {
      SpaceQ2 space_q2(trafo);
      run_space(err[1], space_q2, bnd, String("Q2:") + name);
    }
    {
      SpaceQ3 space_q3(trafo);
      run_space(err[2], space_q3, bnd, String("Q3:") + name);
    }
  }

  void run(int argc, char* argv[])
  {
    SimpleArgParser args(argc, argv);

    //auto atlas = make_atlas_2d();
    auto atlas = make_atlas_3d();

    typedef typename std::remove_reference<decltype(*atlas)>::type AtlasType;
    typedef typename AtlasType::MeshType MeshType;

    static constexpr int dim = MeshType::shape_dim;

    typedef Trafo::Standard::Mapping<MeshType> TrafoStd;
    typedef Trafo::Isoparam::Mapping<MeshType, 1> TrafoDeg1;
    typedef Trafo::Isoparam::Mapping<MeshType, 2> TrafoDeg2;
    typedef Trafo::Isoparam::Mapping<MeshType, 3> TrafoDeg3;

    // create an empty atlas and a root mesh node
    auto node = make_mesh_node(atlas.get());

    const auto& chart = *atlas->find_mesh_chart("chart");

    Index lvl_min = 0;
    Index lvl_max = 7 - dim;

    std::vector<ErrorTensor> errs(lvl_max+1u);

    // refine
    for(Index lvl(0); lvl <= lvl_max; ++lvl)
    {
      if(lvl > 0)
      {
        node = node->refine_unique();
      }

      ErrorTensor& err = errs.at(lvl);

      // get our mesh
      MeshType& mesh = *node->get_mesh();
      const auto& mpart = *node->find_mesh_part("bnd");

      // create trafos
      TrafoStd trafo_0(mesh);
      TrafoDeg1 trafo_1(mesh);
      TrafoDeg2 trafo_2(mesh);
      TrafoDeg3 trafo_3(mesh);

      // add chart and boundary mesh part
      trafo_1.add_meshpart_chart(mpart, chart);
      trafo_2.add_meshpart_chart(mpart, chart);
      trafo_3.add_meshpart_chart(mpart, chart);

      run_trafo(err[0], trafo_0, mpart, String("Trafo0 on Level ") + stringify(lvl));
      run_trafo(err[1], trafo_1, mpart, String("Trafo1 on Level ") + stringify(lvl));
      run_trafo(err[2], trafo_2, mpart, String("Trafo2 on Level ") + stringify(lvl));
      run_trafo(err[3], trafo_3, mpart, String("Trafo3 on Level ") + stringify(lvl));
    }

    std::cout << String(80, '#') << "\n";
    std::cout << String(80, '#') << "\n";
    std::cout << String(80, '#') << "\n";

    std::cout << "H0 Errors" << "\n" << "     ";
    for(int k(0); k < 12; ++k)
      std::cout << "T" << (k/3) << ":Q" << (1+k%3) << "       ";
    std::cout << "   ";
    for(int k(0); k < 12; ++k)
      std::cout << "T" << (k/3) << ":Q" << (1+k%3) << "  ";
    std::cout << "\n";
    for(Index lvl(0); lvl <= lvl_max; ++lvl)
    {
      std::cout << stringify(lvl).pad_front(2) << ":";
      for(int i(0); i < nt; ++i)
        for(int  j(0); j < ns; ++j)
          std::cout << stringify_fp_sci(errs.at(lvl)(i,j,0), 4, 12);

      if(lvl > lvl_min)
      {
        std::cout << " ||";
        for(int i(0); i < nt; ++i)
          for(int  j(0); j < ns; ++j)
            std::cout << stringify_fp_fix(errs.at(lvl-1)(i,j,0) / errs.at(lvl)(i,j,0), 3, 7);
      }
      std::cout << "\n";
    }
    std::cout << "\n";


    std::cout << "H1 Errors" << "\n" << "     ";
    for(int k(0); k < 12; ++k)
      std::cout << "T" << (k/3) << ":Q" << (1+k%3) << "       ";
    std::cout << "   ";
    for(int k(0); k < 12; ++k)
      std::cout << "T" << (k/3) << ":Q" << (1+k%3) << "  ";
    std::cout << "\n";
    for(Index lvl(0); lvl <= lvl_max; ++lvl)
    {
      std::cout << stringify(lvl).pad_front(2) << ":";
      for(int i(0); i < nt; ++i)
        for(int  j(0); j < ns; ++j)
          std::cout << stringify_fp_sci(errs.at(lvl)(i,j,1), 4, 12);

      if(lvl > lvl_min)
      {
        std::cout << " ||";
        for(int i(0); i < nt; ++i)
          for(int  j(0); j < ns; ++j)
            std::cout << stringify_fp_fix(errs.at(lvl-1)(i,j,1) / errs.at(lvl)(i,j,1), 3, 7);
      }
      std::cout << "\n";
    }
    std::cout << "\n";

    std::cout << "H2 Errors" << "\n" << "     ";
    for(int k(0); k < 12; ++k)
      std::cout << "T" << (k/3) << ":Q" << (1+k%3) << "       ";
    std::cout << "   ";
    for(int k(0); k < 12; ++k)
      std::cout << "T" << (k/3) << ":Q" << (1+k%3) << "  ";
    std::cout << "\n";
    for(Index lvl(0); lvl <= lvl_max; ++lvl)
    {
      std::cout << stringify(lvl).pad_front(2) << ":";
      for(int i(0); i < nt; ++i)
        for(int  j(0); j < ns; ++j)
          std::cout << stringify_fp_sci(errs.at(lvl)(i,j,2), 4, 12);

      if(lvl > lvl_min)
      {
        std::cout << " ||";
        for(int i(0); i < nt; ++i)
          for(int  j(0); j < ns; ++j)
            std::cout << stringify_fp_fix(errs.at(lvl)(i,j,2) > 0.0 ? errs.at(lvl-1)(i,j,2) / errs.at(lvl)(i,j,2) : 0.0, 3, 7);
      }
      std::cout << "\n";
    }
  }
} // namespace DbgIsoParam

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  DbgIsoParam::run(argc, argv);
  return 0;
}
