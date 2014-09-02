#include <iostream>
#include <kernel/archs.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/mesh_smoother/dudv_smoother.hpp>
#include <kernel/geometry/export_vtk.hpp>

using namespace FEAST;
/**
 * \brief Wrapper struct as functions do not seem to agree with template template parameters
 **/
template
<
  typename ShapeType_,
  typename DataType_,
  typename MemType_
> struct BdryDeformApp
{
  /**
   * @brief Runs mesh smoother stuff
   *
   **/
  static void run()
  {
    typedef MemType_ MemType;
    typedef DataType_ DataType;
    typedef ShapeType_ ShapeType;

    typedef Geometry::ConformalMesh<ShapeType, ShapeType::dimension, ShapeType::dimension, DataType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;

//    typedef DuDvSmoother<TrafoType, DataType, MemType> SmootherType;
    typedef DataType_ DataType;

    // Mesh and trafo
    Index level(3);
    Geometry::RefineFactory<MeshType,Geometry::UnitStarCubeFactory> mesh_factory(level);
    MeshType mesh(mesh_factory);
    TrafoType trafo(mesh);

    //DataType pi(Math::pi<DataType>());
    DataType deltat(DataType(1e-3));

    // The smoother in all its template glory
    Geometry::DuDvSmoother<TrafoType, DataType, MemType> mr_dudv(trafo);
    //SmootherType mr_dudv(trafo);

    // Call init before tinkering with the boundary coordinates
    mr_dudv.init();

    // Boundary stuff
    typedef typename Geometry::CellSubSet<ShapeType> BoundaryType;
    typedef typename Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;
    BoundaryFactoryType boundary_factory(mesh);
    BoundaryType boundary(boundary_factory);

    Geometry::TargetSet boundary_set = boundary.template get_target_set<0>();

    //Initial boundary deformation
    for(Index i(0); i < boundary.get_num_entities(0); ++i)
    {
      Index j = boundary_set[i];
      DataType tmp0 = mr_dudv._coords[0](j);
      DataType tmp1 = mr_dudv._coords[1](j);
      DataType norm = Math::sqrt(DataType(0.5)) / Math::sqrt( Math::sqr(tmp0 - DataType(0.5)) + Math::sqr(tmp1 - DataType(0.5)) );
      mr_dudv._coords[0](j, (tmp0 - DataType(0.5)) * norm );
      mr_dudv._coords[1](j, (tmp1 - DataType(0.5)) * norm );
      //mr_dudv._coords[0](j, (tmp0 - DataType(0.5)) * norm - Math::sqrt(DataType(0.5)));
      //mr_dudv._coords[1](j, (tmp1 - DataType(0.5)) * norm - Math::sqrt(DataType(0.5)));
      //mr_dudv._coords[0](j, tmp0 - ( Math::sin(DataType(2)*pi*tmp1) )/DataType(1 << (level+2)));
      //mr_dudv._coords[1](j, tmp1 + ( Math::sin(DataType(2)*pi*tmp0) )/DataType(1 << (level+2)));
    }

    // Set coords since we used the boundary deformation
    // mr_dudv.set_coords();
    Geometry::ExportVTK<MeshType> writer_initial_pre(mesh);
    writer_initial_pre.write("pre_initial.vtk");

    mr_dudv.optimise();

    Geometry::ExportVTK<MeshType> writer_initial_post(mesh);
    writer_initial_post.write("post_initial.vtk");

    std::string filename;

    DataType time(0);
    Index n(0);

    DataType* mesh_velocity(new DataType[mesh.get_num_entities(0)]);

    LAFEM::DenseVector<MemType, DataType_> coords_old[MeshType::world_dim];
    for(Index d = 0; d < MeshType::world_dim; ++d)
      coords_old[d]= std::move(LAFEM::DenseVector<MemType, DataType>(mesh.get_num_entities(0)));


    Index outputstep(1);
    deltat /= DataType(outputstep);
    std::cout << "deltat = " << scientify(deltat) << ", outputstep = " << outputstep << std::endl;

    while(time < DataType(0.2))
    {

      std::cout << "timestep " << n << std::endl;
      time+= deltat;

      // Save old vertex coordinates
      for(Index d(0); d < MeshType::world_dim; ++d)
      {
        for(Index i(0); i < mesh.get_num_entities(0); ++i)
          coords_old[d](i, mr_dudv._coords[d](i));
      }

      // Boundary update case
      for(Index i(0); i < boundary.get_num_entities(0); ++i)
      {
        Index j = boundary_set[i];
        DataType tmp0 = mr_dudv._coords[0](j);
        DataType tmp1 = mr_dudv._coords[1](j);

        //DataType norm = Math::sqrt( Math::sqr(tmp0 - 0.5) + Math::sqr(tmp1 - 0.5));
        //mr_dudv._coords[0](j, tmp0 + deltat*tmp1*norm );
        //mr_dudv._coords[1](j, tmp1 - deltat*tmp0*norm );
        mr_dudv._coords[0](j, tmp0 + deltat*DataType(0.25)*((DataType(4)*tmp0 - DataType(2)) + Math::pow((DataType(4)*tmp1 - DataType(2)),DataType(3) ) )) ;
        mr_dudv._coords[1](j, tmp1 - deltat*DataType(0.25)*((DataType(4)*tmp1 - DataType(2)) + Math::pow((DataType(4)*tmp0 - DataType(2)),DataType(3) ) )) ;
        //mr_dudv.set_coords();

      }

      if( n%outputstep || outputstep==1)
      {
        filename = "pre_" + stringify(n) + ".vtk";
        Geometry::ExportVTK<MeshType> writer_pre(mesh);
//        mr_dudv.prepare();
        writer_pre.add_scalar_vertex("mesh_velocity", mesh_velocity);
        std::cout << "Writing " << filename << std::endl;
        writer_pre.write(filename);
      }

      mr_dudv.optimise();

      DataType max_mesh_velocity(-1e10);
      DataType ideltat = DataType(1)/deltat;
      // Compute grid velocity
      for(Index i(0); i < mesh.get_num_entities(0); ++i)
      {
        mesh_velocity[i] = DataType(0);
        for(Index d(0); d < MeshType::world_dim; ++d)
          mesh_velocity[i] += Math::sqr(ideltat*(coords_old[d](i) - mr_dudv._coords[d](i)));

        mesh_velocity[i] = Math::sqrt(mesh_velocity[i]);
        if(mesh_velocity[i] > max_mesh_velocity)
          max_mesh_velocity = mesh_velocity[i];
      }
      std::cout << "max mesh velocity = " << scientify(max_mesh_velocity) << std::endl;

      if( n%outputstep || outputstep==1)
      {

        filename = "post_" + stringify(n) + ".vtk";
        Geometry::ExportVTK<MeshType> writer_post(mesh);
        writer_post.add_scalar_vertex("mesh_velocity", mesh_velocity);
        std::cout << "Writing " << filename << std::endl;
        writer_post.write(filename);

      }

      n++;
    }

    delete mesh_velocity;

  }
}; // struct BdryDeformApp

int main()
{
  typedef Mem::Main MemType;

  BdryDeformApp<Shape::Hypercube<2>, double, MemType>::run();
  return 0;
}
