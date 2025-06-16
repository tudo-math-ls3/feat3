// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/runtime.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/trafo/standard/mapping.hpp>

using namespace FEAT;
using namespace FEAT::Geometry;

template<typename Shape_, int shape_dim_, int face_dim_>
static void validate_index_set_bounds(const Geometry::IndexSetHolder<Shape_>& ish)
{
  std::cout << "IndexSet<"  << shape_dim_ << "," << face_dim_ << ">: validating bounds..." << std::endl;

  const auto& idx = ish.template get_index_set<shape_dim_, face_dim_>();

  const Index n = idx.get_num_entities();
  const Index bound = idx.get_index_bound();

  for(Index i(0); i < n; ++i)
  {
    for(int j(0); j < idx.num_indices; ++j)
    {
      if(idx(i,j) >= bound)
      {
        std::cout << "ERROR: idx(" << i << "," << "j" << ") = " << idx(i,j) << " >= " << bound << " (bound)" << std::endl;
      }
    }
  }
}

template<typename Shape_, int shape_dim_, int face_dim_, int test_dim_>
static void validate_index_set_indices(const Geometry::IndexSetHolder<Shape_>& ish)
{
  std::cout << "IndexSet<"  << shape_dim_ << "," << face_dim_ << ">: validating " << test_dim_ << "-face indices..." << std::endl;

  const auto& idx = ish.template get_index_set<shape_dim_, face_dim_>();
  const auto& idx_s = ish.template get_index_set<shape_dim_, test_dim_>();
  const auto& idx_f = ish.template get_index_set<face_dim_, test_dim_>();

  const Index n = idx.get_num_entities();

  for(Index cell(0); cell < n; ++cell)
  {
    // loop over all faces
    for(int j(0); j < idx.num_indices; ++j)
    {
      const Index face = idx(cell,j);
      // loop over all vertices of the face
      for(int k(0); k < idx_f.num_indices; ++k)
      {
        const Index vtx = idx_f(face, k);
        bool found = false;
        // loop over all vertices of the cell
        for(int i(0); i < idx_s.num_indices; ++i)
        {
          found = found || (idx_s(cell, i) == vtx);
        }
        if(!found)
        {
          std::cout << "ERROR: vertex " << vtx << " present" << std::endl;
        }
      }
    }
  }
}

// checks whether all facet are adjacent to exactly 1 or 2 elements
template<typename Shape_, int shape_dim_ = Shape_::dimension>
static void validate_index_set_facet_conformity(const Geometry::IndexSetHolder<Shape_>& ish)
{
  const auto& idx = ish.template get_index_set<shape_dim_, shape_dim_-1>();

  std::vector<int> counts(idx.get_index_bound(), 0);

  for(Index i(0); i < idx.get_num_entities(); ++i)
  {
    for(int j(0); j < idx.num_indices; ++j)
      ++counts.at(idx(i,j));
  }

  for(std::size_t i(0); i < counts.size(); ++i)
  {
    // note: count = 0 violates surjectivity, so it is checked in another function
    /*if(counts[i] <= 0)
    {
      std::cout << "ERROR: facet " << i << " is not adjacent to any element!" << std::endl;
    }
    else*/
    if(counts[i] > 2)
    {
      std::cout << "ERROR: facet " << i << " is adjacent to more than 2 elements!" << std::endl;
    }
  }
}

template<typename Shape_, int shape_dim_, int face_dim_>
static void validate_index_set_surjectivity(const Geometry::IndexSetHolder<Shape_>& ish)
{
  Adjacency::Graph graph(Adjacency::RenderType::transpose, ish.template get_index_set<shape_dim_, face_dim_>());
  for(Index i(0); i < graph.get_num_nodes_domain(); ++i)
  {
    if(graph.degree(i) <= Index(0))
    {
      std::cout << "ERROR: " << face_dim_ << "-face " << i << " has no adjacent " << shape_dim_ << "-shape!" << std::endl;
    }
  }
}

template<typename Shape_>
class IndexSetValidator;

template<>
class IndexSetValidator<Shape::Hypercube<1>>
{
public:
  typedef Shape::Hypercube<1> ShapeType;
  static void validate(const Geometry::IndexSetHolder<ShapeType>& ish)
  {
    validate_index_set_bounds<ShapeType, 1, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 1, 0>(ish);
    validate_index_set_facet_conformity<ShapeType>(ish);
  }
};

template<>
class IndexSetValidator<Shape::Hypercube<2>>
{
public:
  typedef Shape::Hypercube<2> ShapeType;
  static void validate(const Geometry::IndexSetHolder<ShapeType>& ish)
  {
    validate_index_set_bounds<ShapeType, 1, 0>(ish);
    validate_index_set_bounds<ShapeType, 2, 0>(ish);
    validate_index_set_bounds<ShapeType, 2, 1>(ish);
    validate_index_set_indices<ShapeType, 2, 1, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 1, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 2, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 2, 1>(ish);
    validate_index_set_facet_conformity<ShapeType>(ish);
  }
};

template<>
class IndexSetValidator<Shape::Hypercube<3>>
{
public:
  typedef Shape::Hypercube<3> ShapeType;
  static void validate(const Geometry::IndexSetHolder<ShapeType>& ish)
  {
    validate_index_set_bounds<ShapeType, 1, 0>(ish);
    validate_index_set_bounds<ShapeType, 2, 0>(ish);
    validate_index_set_bounds<ShapeType, 3, 0>(ish);
    validate_index_set_bounds<ShapeType, 2, 1>(ish);
    validate_index_set_bounds<ShapeType, 3, 1>(ish);
    validate_index_set_bounds<ShapeType, 3, 2>(ish);
    validate_index_set_indices<ShapeType, 2, 1, 0>(ish);
    validate_index_set_indices<ShapeType, 3, 1, 0>(ish);
    validate_index_set_indices<ShapeType, 3, 2, 0>(ish);
    validate_index_set_indices<ShapeType, 3, 2, 1>(ish);
    validate_index_set_surjectivity<ShapeType, 1, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 2, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 3, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 2, 1>(ish);
    validate_index_set_surjectivity<ShapeType, 3, 1>(ish);
    validate_index_set_surjectivity<ShapeType, 3, 2>(ish);
    validate_index_set_facet_conformity<ShapeType>(ish);
  }
};

template<>
class IndexSetValidator<Shape::Simplex<2>>
{
public:
  typedef Shape::Simplex<2> ShapeType;
  static void validate(const Geometry::IndexSetHolder<ShapeType>& ish)
  {
    validate_index_set_bounds<ShapeType, 1, 0>(ish);
    validate_index_set_bounds<ShapeType, 2, 0>(ish);
    validate_index_set_bounds<ShapeType, 2, 1>(ish);
    validate_index_set_indices<ShapeType, 2, 1, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 1, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 2, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 2, 1>(ish);
    validate_index_set_facet_conformity<ShapeType>(ish);
  }
};

template<>
class IndexSetValidator<Shape::Simplex<3>>
{
public:
  typedef Shape::Simplex<3> ShapeType;
  static void validate(const Geometry::IndexSetHolder<ShapeType>& ish)
  {
    validate_index_set_bounds<ShapeType, 1, 0>(ish);
    validate_index_set_bounds<ShapeType, 2, 0>(ish);
    validate_index_set_bounds<ShapeType, 3, 0>(ish);
    validate_index_set_bounds<ShapeType, 2, 1>(ish);
    validate_index_set_bounds<ShapeType, 3, 1>(ish);
    validate_index_set_bounds<ShapeType, 3, 2>(ish);
    validate_index_set_indices<ShapeType, 2, 1, 0>(ish);
    validate_index_set_indices<ShapeType, 3, 1, 0>(ish);
    validate_index_set_indices<ShapeType, 3, 2, 0>(ish);
    validate_index_set_indices<ShapeType, 3, 2, 1>(ish);
    validate_index_set_surjectivity<ShapeType, 1, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 2, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 3, 0>(ish);
    validate_index_set_surjectivity<ShapeType, 2, 1>(ish);
    validate_index_set_surjectivity<ShapeType, 3, 1>(ish);
    validate_index_set_surjectivity<ShapeType, 3, 2>(ish);
    validate_index_set_facet_conformity<ShapeType>(ish);
  }
};

template<typename Mesh_>
void validate_jacobian_determinants(Mesh_& mesh)
{
  std::cout << "Validating Jacobian determinants...\n";
  typedef Trafo::Standard::Mapping<Mesh_> TrafoType;
  TrafoType trafo(mesh);

  typedef typename TrafoType::template Evaluator<>::Type TrafoEvaluator;
  TrafoEvaluator trafo_eval(trafo);

  typename TrafoEvaluator::JacobianMatrixType jac_mat;
  typename TrafoEvaluator::DomainPointType dom_point(0.0);

  const Index n = trafo_eval.get_num_cells();
  for(Index i = 0; i < n; ++i)
  {
    trafo_eval.prepare(i);
    trafo_eval.calc_jac_mat(jac_mat, dom_point);
    double det = jac_mat.det();
    if(det < 1E-8)
    {
      std::cout << "ERROR: element " << i << " has non-positive determinant with value " << det << "\n";
    }
    trafo_eval.finish();
  }
}

template<typename Mesh_>
int run_xml(Geometry::MeshFileReader& mesh_reader)
{
  typedef typename Mesh_::ShapeType ShapeType;

  // create an empty atlas and a root mesh node
  Geometry::MeshAtlas<Mesh_> atlas;
  auto node = Geometry::RootMeshNode<Mesh_>::make_unique(nullptr, &atlas);
  Geometry::PartitionSet part_set;

  // try to parse the mesh file
  try
  {
    std::cout << "Parsing mesh files..." << std::endl;
    // Now parse the mesh file
    mesh_reader.parse(*node, atlas, &part_set);
  }
  catch(std::exception& exc)
  {
    std::cerr << "ERROR: " << exc.what() << std::endl;
    return 1;
  }
  catch(...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    return 1;
  }
  std::cout << "Mesh parsed successfully!" << std::endl << std::endl;

  // get the mesh
  const Mesh_& mesh = *node->get_mesh();

  static constexpr int shape_dim = Mesh_::shape_dim;
  static constexpr int world_dim = Mesh_::world_dim;

  std::cout << "Number of mesh entities:" << std::endl;
  for(int i(0); i <= ShapeType::dimension; ++i)
  {
    std::cout << i << ":" << stringify(mesh.get_num_entities(i)).pad_front(12) << std::endl;
  }
  std::cout << std::endl;

  // validate index sets
  IndexSetValidator<ShapeType>::validate(mesh.get_index_set_holder());

  // validate jacobian determinants
  if constexpr (shape_dim == world_dim)
  {
    validate_jacobian_determinants(mesh);
  }

  return 0;
}

int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);

  // This is the list of all supported meshes that could appear in the mesh file
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 2, Real> S2M2D;
  typedef Geometry::ConformalMesh<Shape::Simplex<2>, 3, Real> S2M3D;
  typedef Geometry::ConformalMesh<Shape::Simplex<3>, 3, Real> S3M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 1, Real> H1M1D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 2, Real> H1M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<1>, 3, Real> H1M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 3, Real> H2M3D;
  typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Real> H3M3D;

  // create an empty mesh file reader
  Geometry::MeshFileReader mesh_reader;

  // add all arguments as mesh files
  for(int i(1); i < argc; ++i)
  {
    std::cout << "Adding mesh file '" << argv[i] << "'..." << std::endl;
    mesh_reader.add_mesh_file(argv[i]);
  }

  std::cout << "Reading mesh file root markup..." << std::endl;

  // read root markup
  mesh_reader.read_root_markup();

  // get mesh type
  const String mtype = mesh_reader.get_meshtype_string();

  std::cout << "Mesh Type: " << mtype << std::endl;

  if(mtype == "conformal:hypercube:1:1") return run_xml<H1M1D>(mesh_reader); else
  if(mtype == "conformal:hypercube:1:2") return run_xml<H1M2D>(mesh_reader); else
  if(mtype == "conformal:hypercube:1:3") return run_xml<H1M3D>(mesh_reader); else
  if(mtype == "conformal:hypercube:2:2") return run_xml<H2M2D>(mesh_reader); else
  if(mtype == "conformal:hypercube:2:3") return run_xml<H2M3D>(mesh_reader); else
  if(mtype == "conformal:hypercube:3:3") return run_xml<H3M3D>(mesh_reader); else
  if(mtype == "conformal:simplex:2:2")   return run_xml<S2M2D>(mesh_reader); else
  if(mtype == "conformal:simplex:2:3")   return run_xml<S2M3D>(mesh_reader); else
  if(mtype == "conformal:simplex:3:3")   return run_xml<S3M3D>(mesh_reader); else
  {
    std::cout << "ERROR: unsupported mesh type!" << std::endl;
    return 1;
  }
}
