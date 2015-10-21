#include <kernel/base_header.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/meshopt/rumpf_smoother.hpp>
#include <kernel/meshopt/rumpf_smoother_q1hack.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1hack.hpp>
#include <kernel/geometry/export_vtk.hpp>

using namespace FEAST;
/**
 * \brief This application demonstrates the usage of some of the RumpfSmoother classes for boundary deformations
 *
 * \note Because an application of the (nonlinear) Rumpf smoother requires operations similar to a matrix assembly,
 * Rumpf smoothers are implemented for Mem::Main only.
 *
 * \author Jordi Paul
 *
 * \tparam DT_
 * The precision of the mesh etc.
 *
 * \tparam ShapeType
 * The shape of the mesh's cells
 *
 * \tparam FunctionalType
 * The Rumpf functional variant to use
 *
 * \tparam RumpfSmootherType_
 * The Rumpf smoother variant to use
 *
 **/
template
<
  typename DT_,
  typename ShapeType_,
  template<typename, typename> class FunctionalType_,
  template<typename ... > class RumpfSmootherType_
> struct BdryDeformApp
{
  /**
   * @brief Runs mesh smoother stuff
   *
   **/
  static void run(Index level)
  {
    /// Precision for meshes etc, everything else uses the same data type
    typedef DT_ DataType;
    /// Rumpf Smoothers are implemented for Mem::Main only
    typedef Mem::Main MemType;
    /// So we use Index
    typedef Index IndexType;
    /// Shape of the mesh cells
    typedef ShapeType_ ShapeType;
    /// The complete mesh type
    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension,ShapeType::dimension, DataType> MeshType;
    /// The corresponding transformation
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    /// Our functional type
    typedef FunctionalType_<DataType, ShapeType> FunctionalType;
    /// The Rumpf smoother
    typedef RumpfSmootherType_<TrafoType, FunctionalType> RumpfSmootherType;

    // Create a unit cube, translate and deform it so it becomes a unit circle
    // For hypercubes, we need the UnitStarCubeFactory so elements do not degenerate to triangles at the boundary
    Geometry::RefineFactory<MeshType,Geometry::UnitStarCubeFactory> mesh_factory(level);
    // Create the mesh
    MeshType* mesh(new MeshType(mesh_factory));

    // Create a basic RootMeshNode from this mesh
    Geometry::RootMeshNode<MeshType>* rmn(new Geometry::RootMeshNode<MeshType>(mesh, nullptr));

    // Create a MeshPart for all of the outer boundary
    typedef typename Geometry::MeshPart<MeshType> BoundaryType;
    typedef typename Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;
    BoundaryFactoryType boundary_factory(*mesh);
    BoundaryType* boundary(new BoundaryType(boundary_factory));

    Geometry::TargetSet& boundary_set = boundary->template get_target_set<0>();

    // Add the boundary to the RootMeshNode
    rmn->add_mesh_part("boundary", boundary, nullptr);
    // Boundary stuff: Dirichlet boundary conditions on the outer boundary
    std::deque<String> dirichlet_list;
    dirichlet_list.push_back("boundary");
    std::deque<String> slip_list;

    // Parameters for the Rumpf functional
    DataType fac_norm = DataType(1e-0),fac_det = DataType(1e0), fac_cof = DataType(0), fac_reg(DataType(1e-8));
    // Create the functional
    FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);

    // Create the smoother
    RumpfSmootherType rumpflpumpfl(rmn, dirichlet_list, slip_list, my_functional);
    rumpflpumpfl.init();
    rumpflpumpfl.print();

    // The domain is [0,1] x [0,1], so translate every vertex by (0.5, 0.5)
    // Note that we modify the _coords member of the smoother, which means that we modify boundary conditions but
    // not the acutal mesh. This is done by calling set_coords
    Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> translation(DataType(0.5));

    for(Index i(0); i < mesh->get_num_entities(0); ++i)
      rumpflpumpfl._coords(i, rumpflpumpfl._coords(i) - translation);

    // Now normalise all boundary vertices
    for(Index i(0); i < boundary->get_num_entities(0); ++i)
    {
      Index j = boundary_set[i];
      Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> tmp(rumpflpumpfl._coords(j));
      tmp.normalise();
      rumpflpumpfl._coords(j, tmp);
    }

    // Now write the changes we made to _coords to the mesh
    rumpflpumpfl.set_coords();

    // Arrays for saving the contributions of the different Rumpf functional parts
    DataType* func_norm(new DataType[mesh->get_num_entities(2)]);
    DataType* func_det(new DataType[mesh->get_num_entities(2)]);
    DataType* func_rec_det(new DataType[mesh->get_num_entities(2)]);

    // Compute initial functional value
    DataType fval(0);
    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
    std::cout << "fval pre optimisation = " << stringify_fp_sci(fval) << std::endl;

    // Compute initial functional gradient
    rumpflpumpfl.compute_gradient();

    // Write initial state to file
    Geometry::ExportVTK<MeshType> writer_initial_pre(*mesh);
    writer_initial_pre.add_field_cell_blocked_vector("h", rumpflpumpfl._h);
    writer_initial_pre.add_field_cell("fval", func_norm, func_det, func_rec_det);
    writer_initial_pre.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
    writer_initial_pre.write("pre_initial");

    // Smooth the mesh
    rumpflpumpfl.optimise();

    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
    std::cout << "fval post optimisation = " << stringify_fp_sci(fval) << std::endl;

    // Call prepare() again because the mesh changed due to the optimisation and it was not called again after the
    // last iteration
    rumpflpumpfl.prepare();
    rumpflpumpfl.compute_gradient();

    // Write optimised initial mesh
    Geometry::ExportVTK<MeshType> writer_initial_post(*mesh);
    writer_initial_post.add_field_cell_blocked_vector("h", rumpflpumpfl._h);
    writer_initial_post.add_field_cell("fval", func_norm, func_det, func_rec_det);
    writer_initial_post.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
    writer_initial_post.write("post_initial");

    // For saving the old coordinates
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim>
      coords_old(mesh->get_num_entities(0), DataType(0));
    // For computing the mesh velocity
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim>
      mesh_velocity(mesh->get_num_entities(0), DataType(0));

    // Initial time
    DataType time(0);
    // Timestep size
    DataType deltat(DataType(1e-3));
    std::cout << "deltat = " << stringify_fp_sci(deltat) << std::endl;

    // Counter for timesteps
    Index n(0);
    // Filename for writing .vtu output
    std::string filename;

    while(time < DataType(1e-1))
    {
      std::cout << "timestep " << n << std::endl;
      time+= deltat;

      // Save old vertex coordinates
      coords_old.clone(rumpflpumpfl._coords);

      // Update the boundary
      for(Index i(0); i < boundary->get_num_entities(0); ++i)
      {
        Index j = boundary_set[i];
        Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> tmp0(rumpflpumpfl._coords(j));
        Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> tmp1(rumpflpumpfl._coords(j));

        tmp1(0) += deltat * DataType(0.25)*((DataType(4)*tmp0(0) - DataType(2))
            + Math::pow((DataType(4)*tmp0(1) - DataType(2)),DataType(3) ) );
        tmp1(1) -= deltat*DataType(0.25)*((DataType(4)*tmp1(0) - DataType(2))
            + Math::pow((DataType(4)*tmp0(0) - DataType(2)),DataType(3) ) );

        rumpflpumpfl._coords(j, tmp1);
      }

      // Write new boundary to mesh
      rumpflpumpfl.set_coords();

      // Compute pre-optimisation functional value and gradient
      rumpflpumpfl.prepare();
      fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval pre optimisation = " << stringify_fp_sci(fval) << std::endl;

      // Write pre-optimisation mesh
      filename = "pre_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_pre(*mesh);
      writer_pre.add_field_cell_blocked_vector("h", rumpflpumpfl._h);
      writer_pre.add_field_cell("fval", func_norm, func_det, func_rec_det);
      writer_pre.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
      writer_pre.add_field_vertex_blocked_vector("mesh_velocity", mesh_velocity);
      std::cout << "Writing " << filename << std::endl;
      writer_pre.write(filename);

      // Optimise the mesh
      rumpflpumpfl.optimise();

      // Compute max. mesh velocity
      DataType max_mesh_velocity(-1e10);
      DataType ideltat = DataType(1)/deltat;

      for(Index i(0); i < mesh->get_num_entities(0); ++i)
      {
        mesh_velocity(i, ideltat*(rumpflpumpfl._coords(i) - coords_old(i)));

        DataType my_mesh_velocity(mesh_velocity(i).norm_euclid());

        if(my_mesh_velocity > max_mesh_velocity)
          max_mesh_velocity = my_mesh_velocity;
      }
      std::cout << "max mesh velocity = " << stringify_fp_sci(max_mesh_velocity) << std::endl;

      // Compute post optimisation functional value and gradient
      rumpflpumpfl.prepare();
      fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
      rumpflpumpfl.compute_gradient();
      std::cout << "fval post optimisation = " << stringify_fp_sci(fval) << std::endl;

      // Write post-optimisation mesh
      filename = "post_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_post(*mesh);
      writer_post.add_field_cell_blocked_vector("h", rumpflpumpfl._h);
      writer_post.add_field_cell("fval", func_norm, func_det, func_rec_det);
      writer_post.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
      writer_post.add_field_vertex_blocked_vector("mesh_velocity", mesh_velocity);
      std::cout << "Writing " << filename << std::endl;
      writer_post.write(filename);

      n++;
    }

    // Clean up
    delete rmn;
    delete[] func_norm;
    delete[] func_det;
    delete[] func_rec_det;

  }
}; // struct BdryDeformApp

// Template aliases to easier switch between variants

// Vanilla Rumpf smoother
template<typename A, typename B>
using MySmoother = Meshopt::RumpfSmoother<A, B>;

// Using the Q1 hack
template<typename A, typename B>
using MySmootherQ1Hack = Meshopt::RumpfSmootherQ1Hack<A, B>;

// For the Q1 hack, the functional is a bit more complicated
template<typename A, typename B>
using MyFunctionalQ1Hack = Meshopt::RumpfFunctionalQ1Hack<A, B, Meshopt::RumpfFunctional>;

int main()
{
  Index level(3);

  BdryDeformApp<double, Shape::Hypercube<2>, MyFunctionalQ1Hack, MySmootherQ1Hack>::run(level);
  return 0;
}
