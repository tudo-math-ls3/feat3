#include <iostream>
#include <string>
#include <kernel/archs.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_smoother.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_smoother_q1hack.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_q1.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_q1_d2.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_q1hack.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_p1.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_p1_d2.hpp>
#include <kernel/geometry/export_vtk.hpp>

using namespace FEAST;

/// \cond internal
// Helper class for resizing tests
template<typename ShapeType_>
struct helperclass;

// Specialisation for hypercubes
template<int shape_dim_>
struct helperclass< FEAST::Shape::Hypercube<shape_dim_> >
{
  // Sets coordinates so we deal the the reference element
  template<typename VectorType_, typename DataType_>
  static void set_coords(VectorType_& coords_, const DataType_& scaling)
  {
    for(Index i(0); i < Index(1 << shape_dim_); ++i)
    {
      for(Index d(0); d < Index(shape_dim_); ++d)
        coords_[d](i, (DataType_(((i >> d) & 1) << 1) - DataType_(1)) * scaling );
    }
  }

  // Prints the ShapeType
  static std::string print_typename()
  {
    return "Hypercube<" + stringify(shape_dim_) +">";
  }
};

// Specialisation for 2d simplices
template<>
struct helperclass< FEAST::Shape::Simplex<2> >
{
  // \brief Sets coordinates so we deal the the Rumpf reference element
  template<typename VectorType_, typename DataType_>
  static void set_coords(VectorType_& coords_, const DataType_& scaling)
  {
    coords_[0](0, DataType_(0));
    coords_[1](0, DataType_(0));

    coords_[0](1, DataType_(1) * scaling);
    coords_[1](1, DataType_(0));

    coords_[0](2, DataType_(0.5) * scaling);
    coords_[1](2, DataType_(0.5)*Math::sqrt(DataType_(3))*scaling);
  }
  // Prints the ShapeType
  static std::string print_typename()
  {
    return "Simplex<2>";
  }
};

/// \endcond

/**
 * \brief Wrapper struct as functions do not seem to agree with template template parameters
 **/
template
<
  typename DataType_,
  typename MemType_,
  typename ShapeType_,
  template<typename, typename> class FunctionalType_,
  template<typename ... > class RumpfSmootherType_
> struct ResizeApp
{
/**
 * @brief Runs mesh smoother stuff
 *
 **/
static void run()
{
  typedef DataType_ DataType;
  typedef MemType_ MemType;
  typedef ShapeType_ ShapeType;

  typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension,ShapeType::dimension, DataType> MeshType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;

  typedef FunctionalType_<DataType, ShapeType> FunctionalType;
  typedef RumpfSmootherType_<DataType_, MemType, TrafoType, FunctionalType> RumpfSmootherType;
  typedef DataType_ DataType;


  // Mesh and trafo
  Geometry::ReferenceCellFactory<ShapeType, DataType> mesh_factory;
  MeshType mesh(mesh_factory);
  TrafoType trafo(mesh);

  DataType fac_norm = DataType(1e-2),fac_det = DataType(1e0),fac_cof = DataType(0), fac_reg(DataType(0e0));

  FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);

  // The smoother in all its template glory
  RumpfSmootherType rumpflpumpfl(trafo, my_functional);

  // Set bdry_id to 0 again so the mesh smoother can resize the single element
  for(Index i(0); i < mesh.get_num_entities(0); ++i)
    rumpflpumpfl._bdry_id[i] = 0;

  // Set initial coordinates by scaling the original Rumpf reference cell by ...
  DataType scaling(5);
  helperclass<ShapeType>::set_coords(rumpflpumpfl._coords, scaling);

  // Since we changed the internal _coords, they have to be copied back to the mesh
  rumpflpumpfl.set_coords();
  rumpflpumpfl.init();

  rumpflpumpfl.print();

  // Set target edge lengths
  rumpflpumpfl._h[0](0, DataType(4));
  rumpflpumpfl._h[1](0, DataType(4));
  //rumpflpumpfl._h[0](0, DataType(4)/Math::sqrt(DataType(3)));
  //rumpflpumpfl._h[1](0, DataType(4)/(Math::sqrt(DataType(2))*Math::sqrt(Math::sqrt(DataType(3)))));

  DataType fval(0);
  DataType* func_norm(new DataType[mesh.get_num_entities(2)]);
  DataType* func_det(new DataType[mesh.get_num_entities(2)]);
  DataType* func_det2(new DataType[mesh.get_num_entities(2)]);

  // Evaluates the levelset function and its gradient
  rumpflpumpfl.prepare();
  // Compute initial functional value
  fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_det2);
  std::cout << "fval pre optimisation = " << scientify(fval) << std::endl;
  // Compute initial functional gradient
  rumpflpumpfl.compute_gradient();

  std::string filename;
  filename = "pre_" + helperclass<ShapeType>::print_typename() + ".vtk";

  Geometry::ExportVTK<MeshType> writer_pre(mesh);
  writer_pre.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
  writer_pre.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh.get_num_entities(0)]);
  writer_pre.add_scalar_cell("norm_A", func_norm);
  writer_pre.add_scalar_cell("det_A", func_det);
  writer_pre.add_scalar_cell("det2_A", func_det2);
  writer_pre.add_scalar_cell("h_0", rumpflpumpfl._h[0].elements() );
  writer_pre.add_scalar_cell("h_1", rumpflpumpfl._h[1].elements() );
  writer_pre.write(filename);

  rumpflpumpfl.optimise();
  fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_det2);
  std::cout << "fval post optimisation = " << scientify(fval) << std::endl;

  rumpflpumpfl.prepare();
  rumpflpumpfl.compute_gradient();

  filename = "post_" + helperclass<ShapeType>::print_typename() + ".vtk";

  Geometry::ExportVTK<MeshType> writer_post(mesh);
  writer_post.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
  writer_post.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh.get_num_entities(0)]);
  writer_post.add_scalar_cell("norm_A", func_norm);
  writer_post.add_scalar_cell("det_A", func_det);
  writer_post.add_scalar_cell("det2_A", func_det2);
  writer_post.add_scalar_cell("h_0", rumpflpumpfl._h[0].elements() );
  writer_post.add_scalar_cell("h_1", rumpflpumpfl._h[1].elements() );
  writer_post.write(filename);
}
};

  template<typename A, typename B, typename C, typename D>
  using MySmoother = Geometry::RumpfSmoother<A, B, C, D>;

  template<typename A, typename B, typename C, typename D>
  using MySmootherQ1Hack = Geometry::RumpfSmootherQ1Hack<A, B, C, D>;
int main()
{
  typedef Mem::Main MemType;

  ResizeApp<double, MemType, Shape::Hypercube<2>, Geometry::RumpfFunctional_D2, MySmoother>::run();

  return 0;
}
