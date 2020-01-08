// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/runtime.hpp>
#include <kernel/util/random.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/trafo/standard/mapping.hpp>

using namespace FEAT;

namespace DbgCubature
{
  // the data type used in this code
#ifdef FEAT_HAVE_QUADMATH
  typedef __float128 DataType;
#else
  typedef double DataType;
#endif

  template<typename Mesh_>
  void disturb_mesh(Mesh_& mesh, DataType disto)
  {
    auto& vtx = mesh.get_vertex_set();

    const DataType pi2 = DataType(2) * Math::pi<DataType>();
    Random rng;

    for(Index i(0); i < vtx.get_num_vertices(); ++i)
    {
      // generate rotation matrix
      DataType t = rng(DataType(0), pi2);
      if((vtx[i][0] > DataType(0.0001)) && (vtx[i][0] < DataType(0.9999)))
        vtx[i][0] += disto * Math::cos(t);
      if((vtx[i][1] > DataType(0.0001)) && (vtx[i][1] < DataType(0.9999)))
        vtx[i][1] += disto * Math::sin(t);
    }
  }

  template<typename Shape_>
  std::vector<std::array<int, Shape_::dimension>> build_powers(const int degree);

  template<>
  std::vector<std::array<int, 2>> build_powers<Shape::Simplex<2>>(const int degree)
  {
    std::vector<std::array<int, 2>> p;
    std::array<int, 2> q;

    for(int i = 0; i <= degree; ++i)
    {
      for(int j = 0; j <= i; ++j)
      {
        q[0] = i-j;
        q[1] = j;
        p.push_back(q);
      }
    }

    return p;
  }

  template<>
  std::vector<std::array<int, 3>> build_powers<Shape::Simplex<3>>(const int degree)
  {
    std::vector<std::array<int, 3>> p;
    std::array<int, 3> q;

    for(int i = 0; i <= degree; ++i)
    {
      for(int j = 0; j <= i; ++j)
      {
        for(int k = 0; j+k <= i; ++k)
        {
          q[0] = i-j-k;
          q[1] = k;
          q[2] = j;
          p.push_back(q);
        }
      }
    }

    return p;
  }

  template<>
  std::vector<std::array<int, 2>> build_powers<Shape::Hypercube<2>>(const int degree)
  {
    std::vector<std::array<int, 2>> p;
    std::array<int, 2> q;

    for(int i = 0; i <= degree; ++i)
    {
      for(int j = 0; j <= i; ++j)
      {
        q[0] = i;
        q[1] = j;
        p.push_back(q);
      }
      for(int j = 1; j <= i; ++j)
      {
        q[0] = i-j;
        q[1] = i;
        p.push_back(q);
      }
    }

    return p;
  }

  template<>
  std::vector<std::array<int, 3>> build_powers<Shape::Hypercube<3>>(const int degree)
  {
    std::vector<std::array<int, 3>> p;
    std::array<int, 3> q;

    for(int i = 0; i <= degree; ++i)
    {
      for(int j = 0; j <= i; ++j)
      {
        for(int k = 0; k <= i; ++k)
        {
          q[0] = i;
          q[1] = j;
          q[2] = k;
          p.push_back(q);
        }
      }
      for(int j = 1; j <= i; ++j)
      {
        for(int k = 0; k <= i; ++k)
        {
          q[0] = i-j;
          q[1] = i;
          q[2] = k;
          p.push_back(q);
        }
      }
      for(int j = 0; j < i; ++j)
      {
        for(int k = 0; k < i; ++k)
        {
          q[0] = j;
          q[1] = k;
          q[2] = i;
          p.push_back(q);
        }
      }
    }

    return p;
  }

  template<typename Shape_>
  int mono_degree(const std::array<int, Shape_::dimension>& p);

  template<>
  int mono_degree<Shape::Simplex<2>>(const std::array<int, 2>& p)
  {
    return p[0] + p[1];
  }

  template<>
  int mono_degree<Shape::Simplex<3>>(const std::array<int, 3>& p)
  {
    return p[0] + p[1] + p[2];
  }

  template<>
  int mono_degree<Shape::Hypercube<2>>(const std::array<int, 2>& p)
  {
    return Math::max(p[0], p[1]);
  }

  template<>
  int mono_degree<Shape::Hypercube<3>>(const std::array<int, 3>& p)
  {
    return Math::max(p[0], Math::max(p[1], p[2]));
  }

  template<std::size_t dim_>
  void print_monomial(const std::array<int, dim_>& p)
  {
    static const char vs[3] = {'x', 'y', 'z'};
    int deg = 0;
    for(std::size_t d = 0u; d < dim_; ++d)
    {
      std::cout << ((p[d] > 0) && (deg > 0) ? "* " : "  ");
      if(p[d] > 1)
        std::cout << vs[d] << '^' << stringify(p[d]).pad_back(2);
      else if(p[d] > 0)
        std::cout << vs[d] << "   ";
      else
        std::cout << "    ";
      deg += p[d];
    }
  }

  template<std::size_t dim_, int n_>
  DataType eval_monomial(const std::array<int, dim_>& powers, const Tiny::Vector<DataType, n_, n_>& x)
  {
    DataType v = DataType(1);
    for(std::size_t i(0u); i < dim_; ++i)
      for(int j(0); j < powers[i]; ++j)
        v *= x[int(i)];
    return v;
  }

  template<std::size_t dim_>
  DataType int_monomial(const std::array<int, dim_>& powers)
  {
    DataType v = DataType(1);
    for(std::size_t i(0u); i < dim_; ++i)
    {
      v /= DataType(powers[i]+1);
    }
    return v;
  }

  template<typename Shape_>
  void run(SimpleArgParser& args)
  {
    static constexpr int dim = Shape_::dimension;

    typedef Geometry::ConformalMesh<Shape_, Shape_::dimension, DataType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

    // define trafo evaluator and eval data types
    typedef typename TrafoType::template Evaluator<Shape_, DataType>::Type TrafoEvaluator;
    static constexpr TrafoTags trafo_config = TrafoTags::img_point | TrafoTags::jac_det;
    typedef typename TrafoEvaluator::template ConfigTraits<trafo_config>::EvalDataType TrafoEvalData;


    Index level(0u);
    args.parse("level", level);

    DataType disto(0.0);
    args.parse("disto", disto);
    disto /= DataType(1 << level);

    int degree(0);
    args.parse("degree", degree);

    DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.7));
    args.parse("tol", tol);


    std::cout << "Shape...: " << Shape_::name() << std::endl;
    std::cout << "DataType: " << Type::Traits<DataType>::name() << std::endl;
    std::cout << "Level...: " << level << std::endl;
    std::cout << "Disto...: " << disto << std::endl;
    std::cout << "Degree..: " << degree << std::endl;
    std::cout << "Tol.....: " << tol << std::endl;

    // create mesh
    Geometry::RefinedUnitCubeFactory<MeshType> factory(level);
    MeshType mesh(factory);
    if(disto > DataType(1E-12))
      disturb_mesh(mesh, disto);
    TrafoType trafo(mesh);

    const Index num_elems = mesh.get_num_elements();
    std::cout << "Elements: " << num_elems << std::endl;
    std::cout << std::endl;

    // create cubature rule
    String cub_rule;
    args.parse("rule", cub_rule);
    std::cout << "Rule Name Given.: " << cub_rule << std::endl;

    // create cubature factory
    Cubature::DynamicFactory cubature_factory(cub_rule);

    // create actual rule
    Cubature::Rule<Shape_, DataType, DataType> cubature;
    try
    {
      cubature = Cubature::Rule<Shape_, DataType, DataType>(Cubature::ctor_factory, cubature_factory);
    }
    catch(std::exception& e)
    {
      std::cout << "ERROR: Failed to create cubature rule" << std::endl;
      std::cout << e.what() << std::endl;
      return;
    }
    std::cout << "Rule Name Mapped: " << cubature.get_name() << std::endl;
    std::cout << "Number of Points: " << cubature.get_num_points() << std::endl;
    std::cout << std::endl;


    // build powers vector
    std::vector<std::array<int,dim>> powers = build_powers<Shape_>(degree);

    // monomial integrals
    std::vector<DataType> monomial_ints(powers.size(), DataType(0));

    // create the trafo evaluator and its data
    TrafoEvaluator trafo_eval(trafo);
    TrafoEvalData trafo_data;

    // loop over all elements of the mesh
    for(Index ielem = 0; ielem < num_elems; ++ielem)
    {
      trafo_eval.prepare(ielem);

      // loop over all cubature points
      for(int ipoint = 0; ipoint < cubature.get_num_points(); ++ipoint)
      {
        // evaluate trafo
        trafo_eval(trafo_data, cubature.get_point(ipoint));

        // get integration weight
        const DataType weight = cubature.get_weight(ipoint) * trafo_data.jac_det;

        // loop over all monomials
        for(std::size_t imono(0u); imono < powers.size(); ++imono)
        {
          monomial_ints[imono] += weight * eval_monomial(powers[imono], trafo_data.img_point);
        }
      }

      trafo_eval.finish();
    }

    int max_degree = degree;

    std::cout << "  Monomial" << String(6*dim-8, ' ');
    std::cout << "Deg   Cubature       Exact          Relative Error" << std::endl;

    // loop over all monomials
    for(std::size_t imono(0u); imono < powers.size(); ++imono)
    {
      // print monomial to console
      print_monomial(powers[imono]);
      const int deg = mono_degree<Shape_>(powers[imono]);
      std::cout << " : " << stringify(deg).pad_front(2) << " > ";

      // compute exact integral of monomial
      const DataType exact_int = int_monomial<dim>(powers[imono]);

      // compute relative integration error
      const DataType int_err = Math::abs(exact_int - monomial_ints[imono]) / exact_int;

      // print integral value
      std::cout << stringify_fp_fix(monomial_ints[imono], 10);
      std::cout << " - " << stringify_fp_fix(exact_int, 10);
      std::cout << " = " << stringify_fp_sci(int_err);

      // significant integration error?
      if(int_err > tol)
      {
        std::cout << " < !!!";
        max_degree = Math::min(max_degree, deg-1);
      }

      std::cout << std::endl;
    }

    std::cout << std::endl;
    std::cout << "Maximum Exact Degree: " << max_degree << std::endl;

    // number of cubature points
  }

  void main(int argc, char** argv)
  {
    SimpleArgParser args(argc, argv);

    args.support("shape");
    args.support("level");
    args.support("tol");
    args.support("degree");
    args.support("disto");
    args.support("rule");

    std::deque<std::pair<int,String>> unsupported = args.query_unsupported();
    if(!unsupported.empty())
    {
      std::cerr << std::endl;
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unsupported option #" << (*it).first << " '--" << (*it).second << "'" << std::endl;
      return;
    }

    if(args.check("shape") <= 0)
    {
      std::cout << "ERROR: mandatory shape parameter '--shape <shape>' is missing!" << std::endl;
      std::cout << "Valid shape parameters are:" << std::endl;
      std::cout << "s2         Simplex<2>" << std::endl;
      std::cout << "s3         Simplex<3>" << std::endl;
      std::cout << "h2         Hypercube<2>" << std::endl;
      std::cout << "h3         Hypercube<3>" << std::endl;
      return;
    }

    String shape;
    args.parse("shape", shape);
    if(shape.compare_no_case("s2") == 0) run<Shape::Simplex<2>>(args); else
    if(shape.compare_no_case("s3") == 0) run<Shape::Simplex<3>>(args); else
    if(shape.compare_no_case("h2") == 0) run<Shape::Hypercube<2>>(args); else
    if(shape.compare_no_case("h3") == 0) run<Shape::Hypercube<3>>(args); else
    {
      std::cout << "ERROR: Failed to parse '" << shape << "' as shape parameter" << std::endl;
      std::cout << "Valid shape parameters are:" << std::endl;
      std::cout << "s2         Simplex<2>" << std::endl;
      std::cout << "s3         Simplex<3>" << std::endl;
      std::cout << "h2         Hypercube<2>" << std::endl;
      std::cout << "h3         Hypercube<3>" << std::endl;
      return;
    }
  }
} // namespace DbgCubature

int main(int argc, char** argv)
{
  Runtime::initialise(argc, argv);
  DbgCubature::main(argc, argv);
  return Runtime::finalise();
}
