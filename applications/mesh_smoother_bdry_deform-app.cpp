#include <iostream>
#include <kernel/archs.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_smoother.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_q1.hpp>
#include <kernel/geometry/mesh_smoother/rumpf_functional_2d_p1.hpp>
//#include <kernel/geometry/mesh_smoother/rumpf_functional_levelset_2d_p1.hpp>
//#include <kernel/geometry/mesh_smoother/rumpf_functional_levelset_2d_q1.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/assembly/interpolator.hpp>
#include <kernel/assembly/common_functions.hpp>
#include <kernel/assembly/discrete_projector.hpp>

using namespace FEAST;

  /**
   * @brief Runs mesh smoother stuff
   *
   **/
  template<typename DataType_, typename ShapeType_ >
  void run()
  {
    typedef DataType_ DataType;
    typedef Mem::Main MemType;
    typedef ShapeType_ ShapeType;

    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef typename Geometry::RumpfFunctional<MemType, DataType, TrafoType> FunctionalType;

    typedef LAFEM::DenseVector<MemType, DataType> VectorType;


    // Mesh and trafo
    Index level(3);
    Geometry::RefinedUnitCubeFactory<MeshType> mesh_factory(level);
    MeshType mesh(mesh_factory);
    TrafoType trafo(mesh);

    DataType pi(Math::pi<DataType>());
    DataType deltat(DataType(0.0125));

    DataType fac_norm = DataType(1e0),fac_det = DataType(1e0), fac_cof = DataType(0), fac_reg(DataType(1e-8));
    FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);

    // The smoother in all its template glory
    Geometry::RumpfSmoother
    <
      FunctionalType,
      TrafoType,
      DataType,
      MemType
    > rumpflpumpfl(trafo, my_functional);

    // Call init before tinkering with the boundary coordinates
    rumpflpumpfl.init();

    // Boundary stuff
    typedef typename Geometry::CellSubSet<ShapeType> BoundaryType;
    typedef typename Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;
    BoundaryFactoryType boundary_factory(mesh);
    BoundaryType boundary(boundary_factory);

    Geometry::TargetSet boundary_set = boundary.template get_target_set<0>();

    //Initial boundary deformation
    //for(Index i(0); i < boundary.get_num_entities(0); ++i)
    //{
    //  Index j = boundary_set[i];
    //  DataType tmp0 = rumpflpumpfl._coords[0](j);
    //  DataType tmp1 = rumpflpumpfl._coords[1](j);
    //  rumpflpumpfl._coords[0](j, tmp0 - ( Math::sin(DataType(2)*pi*tmp1) )/DataType(1 << (level+2)));
    //  rumpflpumpfl._coords[1](j, tmp1 + ( Math::sin(DataType(2)*pi*tmp0) )/DataType(1 << (level+2)));
    //}

    DataType* func_norm(new DataType[mesh.get_num_entities(2)]);
    DataType* func_det(new DataType[mesh.get_num_entities(2)]);
    DataType* func_det2(new DataType[mesh.get_num_entities(2)]);

    // Set coords since we used the boundary deformation
    rumpflpumpfl.set_coords();
    // Compute initial functional value
    DataType fval(0);
    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_det2);
    std::cout << "fval pre optimisation = " << scientify(fval) << std::endl;
    // Compute initial functional gradient
    rumpflpumpfl.compute_gradient();

    Geometry::ExportVTK<MeshType> writer_pre(mesh);
    writer_pre.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
    writer_pre.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh.get_num_entities(0)]);
    writer_pre.add_scalar_cell("norm_A", func_norm);
    writer_pre.add_scalar_cell("det_A", func_det);
    writer_pre.add_scalar_cell("det2_A", func_det2);
    writer_pre.add_scalar_cell("h", rumpflpumpfl._h[0].elements() );
    writer_pre.write("pre_initial.vtk");

    rumpflpumpfl.optimise();
    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_det2);
    std::cout << "fval post optimisation = " << scientify(fval) << std::endl;

    rumpflpumpfl.prepare();
    rumpflpumpfl.compute_gradient();

    Geometry::ExportVTK<MeshType> writer_post(mesh);
    writer_post.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
    writer_post.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh.get_num_entities(0)]);
    writer_post.add_scalar_cell("norm_A", func_norm);
    writer_post.add_scalar_cell("det_A", func_det);
    writer_post.add_scalar_cell("det2_A", func_det2);
    writer_post.add_scalar_cell("h", rumpflpumpfl._h[0].elements() );
    writer_post.write("post_initial.vtk");


    /*std::string filename;

    DataType time(0);
    Index n(0);

    DataType* mesh_velocity(new DataType[mesh.get_num_entities(0)]);

    LAFEM::DenseVector<MemType, DataType_> coords_old[MeshType::world_dim];
    for(Index d = 0; d < MeshType::world_dim; ++d)
      coords_old[d]= std::move(LAFEM::DenseVector<MemType, DataType>(mesh.get_num_entities(0)));


    Index outputstep(4);
    deltat /= DataType(outputstep);
    std::cout << "deltat = " << scientify(deltat) << ", outputstep = " << outputstep << std::endl;
    std::cout << "fac_norm = " << scientify(fac_norm) << ", fac_det= " << scientify(fac_det)
    << ", fac_reg = " << scientify(fac_reg) << std::endl;

    while(time < DataType(6))
    {

      std::cout << "timestep " << n << std::endl;
      time+= deltat;

      // Save old vertex coordinates
      for(Index d(0); d < MeshType::world_dim; ++d)
      {
        for(Index i(0); i < mesh.get_num_entities(0); ++i)
          coords_old[d](i, rumpflpumpfl._coords[d](i));
      }

      // Boundary update case
      for(Index i(0); i < boundary.get_num_entities(0); ++i)
      {
        Index j = boundary_set[i];
        DataType tmp0 = rumpflpumpfl._coords[0](j);
        DataType tmp1 = rumpflpumpfl._coords[1](j);

        //DataType norm = Math::sqrt( Math::sqr(tmp0 - 0.5) + Math::sqr(tmp1 - 0.5));
        //rumpflpumpfl._coords[0](j, tmp0 + deltat*tmp1*norm );
        //rumpflpumpfl._coords[1](j, tmp1 - deltat*tmp0*norm );
        rumpflpumpfl._coords[0](j, tmp0 + deltat*DataType(0.25)*((DataType(4)*tmp0 - DataType(2)) + Math::pow((DataType(4)*tmp1 - DataType(2)),DataType(3) ) )) ;
        rumpflpumpfl._coords[1](j, tmp1 - deltat*DataType(0.25)*((DataType(4)*tmp1 - DataType(2)) + Math::pow((DataType(4)*tmp0 - DataType(2)),DataType(3) ) )) ;
        rumpflpumpfl.set_coords();

      }

      if( n%outputstep || outputstep==1)
      {
        filename = "pre_" + stringify(n) + ".vtk";
        Geometry::ExportVTK<MeshType> writer_pre(mesh);
        rumpflpumpfl.prepare();
        fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_det2);
        std::cout << "fval pre optimisation = " << scientify(fval) << std::endl;
        writer_pre.add_scalar_cell("norm", func_norm);
        writer_pre.add_scalar_cell("det", func_det);
        writer_pre.add_scalar_cell("det2", func_det2);
        writer_pre.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
        writer_pre.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh.get_num_entities(0)]);
        writer_pre.add_scalar_cell("h", rumpflpumpfl._h[0].elements() );
        writer_pre.add_scalar_vertex("mesh_velocity", mesh_velocity);
        std::cout << "Writing " << filename << std::endl;
        writer_pre.write(filename);
      }

      rumpflpumpfl.optimise();

      DataType max_mesh_velocity(-1e10);
      DataType ideltat = DataType(1)/deltat;
      // Compute grid velocity
      for(Index i(0); i < mesh.get_num_entities(0); ++i)
      {
        mesh_velocity[i] = DataType(0);
        for(Index d(0); d < MeshType::world_dim; ++d)
          mesh_velocity[i] += Math::sqr(ideltat*(coords_old[d](i) - rumpflpumpfl._coords[d](i)));

        mesh_velocity[i] = Math::sqrt(mesh_velocity[i]);
        if(mesh_velocity[i] > max_mesh_velocity)
          max_mesh_velocity = mesh_velocity[i];
      }
      std::cout << "max mesh velocity = " << scientify(max_mesh_velocity) << std::endl;

      if( n%outputstep || outputstep==1)
      {

        filename = "post_" + stringify(n) + ".vtk";
        Geometry::ExportVTK<MeshType> writer_post(mesh);
        fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_det2);
        std::cout << "fval post optimisation = " << scientify(fval) << std::endl;
        rumpflpumpfl.compute_gradient();
        writer_post.add_scalar_cell("norm", func_norm);
        writer_post.add_scalar_cell("det", func_det);
        writer_post.add_scalar_cell("det2", func_det2);
        writer_post.add_scalar_vertex("grad_0", &rumpflpumpfl._grad[0]);
        writer_post.add_scalar_vertex("grad_1", &rumpflpumpfl._grad[mesh.get_num_entities(0)]);
        writer_post.add_scalar_vertex("mesh_velocity", mesh_velocity);
        writer_post.add_scalar_cell("h", rumpflpumpfl._h[0].elements() );
        std::cout << "Writing " << filename << std::endl;
        writer_post.write(filename);

      }

      n++;
    }

    delete func_norm;
    delete func_det;
    delete func_det2;

    delete mesh_velocity;
  */
  }

int main()
{
  run<double,Shape::Hypercube<2>>();
  return 0;
}
