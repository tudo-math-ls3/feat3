#include <kernel/base_header.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/meshopt/dudv_smoother.hpp>
#include <kernel/meshopt/laplace_smoother.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/runtime.hpp>

using namespace FEAST;
/**
 * \brief This application demonstrates the usage of linear variational mesh smoothers boundary deformations
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
 * \tparam SmootherType_
 * The smoother variant to use
 *
 * Possible smoother variants are:
 * LaplaceSmoother
 * DuDvSmoother
 *
 **/
template
<
  typename Mem_,
  typename DT_,
  typename IT_,
  typename ShapeType_,
  template <typename, typename, typename, typename > class SmootherType_
> struct BdryDeformApp
{
  /**
   * @brief Runs mesh smoother stuff
   *
   **/
  static void run(Index level)
  {
    /// Memory architecture for the solver
    typedef Mem_ MemType;
    /// Precision for the solver
    typedef DT_ DataType;
    /// Indexing type
    typedef IT_ IndexType;
    /// Shape of the mesh cells
    typedef ShapeType_ ShapeType;
    /// The complete mesh type
    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension,ShapeType::dimension, DataType> MeshType;
    /// The corresponding transformation
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    /// The full mesh smoother class
    typedef SmootherType_<Mem_, DT_, IT_, TrafoType> SmootherType;

    // Create a unit cube, translate and deform it so it becomes a unit circle
    // For hypercubes, we need the UnitStarCubeFactory so elements do not degenerate to triangles at the boundary
    Geometry::RefineFactory<MeshType,Geometry::UnitStarCubeFactory> mesh_factory(level);
    // Create the mesh
    MeshType* mesh(new MeshType(mesh_factory));

    // Create simple RootMeshNode from this mesh
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

    // Creat the mesh smoother
    SmootherType mr_laplace(rmn, dirichlet_list, slip_list);
    mr_laplace.init();

    // The domain is [0,1] x [0,1], so translate every vertex by (0.5, 0.5)
    // Note that we modify the _coords member of the smoother, which means that we modify boundary conditions but
    // not the acutal mesh. This is done by calling set_coords
    Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> translation(DataType(0.5));

    for(Index i(0); i < mesh->get_num_entities(0); ++i)
      mr_laplace._coords(i, mr_laplace._coords(i) - translation);

    // Now normalise all boundary vertices
    for(Index i(0); i < boundary->get_num_entities(0); ++i)
    {
      Index j = boundary_set[i];
      Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> tmp(mr_laplace._coords(j));
      tmp.normalise();
      mr_laplace._coords(j, tmp);
    }

    // Important: Do not call set_coords, as we need to solve on the old mesh using the values in _coords as
    // boundary values

    // Write initial state to file
    Geometry::ExportVTK<MeshType> writer_initial_pre(*mesh);
    writer_initial_pre.write("pre_initial");

    // Smooth the mesh
    mr_laplace.optimise();

    // Write optimised initial mesh
    Geometry::ExportVTK<MeshType> writer_initial_post(*mesh);
    writer_initial_post.write("post_initial");

    // For saving the old coordinates
    LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, MeshType::world_dim>
      coords_old(mesh->get_num_entities(0),DataType(0));
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
      coords_old.clone(mr_laplace._coords);

      // Update the boundary
      for(Index i(0); i < boundary->get_num_entities(0); ++i)
      {
        Index j = boundary_set[i];
        Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> tmp0(mr_laplace._coords(j));
        Tiny::Vector<DataType, MeshType::world_dim, MeshType::world_dim> tmp1(mr_laplace._coords(j));

        tmp1(0) += deltat * DataType(0.25)*((DataType(4)*tmp0(0) - DataType(2))
            + Math::pow((DataType(4)*tmp0(1) - DataType(2)),DataType(3) ) );
        tmp1(1) -= deltat*DataType(0.25)*((DataType(4)*tmp1(0) - DataType(2))
            + Math::pow((DataType(4)*tmp0(0) - DataType(2)),DataType(3) ) );

        mr_laplace._coords(j, tmp1);
      }

      // Important: Do not call set_coords, as we need to solve on the old mesh using the values in _coords as
      // boundary values

      // Write pre-optimisation mesh
      filename = "pre_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_pre(*mesh);
      writer_pre.add_field_vertex_blocked_vector("mesh_velocity", mesh_velocity);
      std::cout << "Writing " << filename << std::endl;
      writer_pre.write(filename);

      // Optimise the mesh
      mr_laplace.optimise();

      // Compute max. mesh velocity
      DataType max_mesh_velocity(-1e10);
      DataType ideltat = DataType(1)/deltat;

      for(Index i(0); i < mesh->get_num_entities(0); ++i)
      {
        mesh_velocity(i, ideltat*(mr_laplace._coords(i) - coords_old(i)));

        DataType my_mesh_velocity(mesh_velocity(i).norm_euclid());

        if(my_mesh_velocity > max_mesh_velocity)
          max_mesh_velocity = my_mesh_velocity;
      }
      std::cout << "max mesh velocity = " << stringify_fp_sci(max_mesh_velocity) << std::endl;

      // Write post-optimisation mesh
      filename = "post_" + stringify(n);
      Geometry::ExportVTK<MeshType> writer_post(*mesh);
      writer_post.add_field_vertex_blocked_vector("mesh_velocity", mesh_velocity);
      std::cout << "Writing " << filename << std::endl;
      writer_post.write(filename);

      n++;

    } // time loop

    // Clean up
    delete rmn;
  } // void run()
}; // struct BdryDeformApp

int main(int argc, char** argv)
{
  FEAST::Runtime::initialise(argc, argv);
  Index level(3);

  BdryDeformApp<Mem::Main, double, Index, Shape::Hypercube<2>, Meshopt::DuDvSmoother>::run(level);
  return FEAST::Runtime::finalise();
}
