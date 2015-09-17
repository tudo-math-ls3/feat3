#include <kernel/base_header.hpp>
#include <kernel/util/math.hpp>
#include <kernel/geometry/boundary_factory.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/reference_cell_factory.hpp>
#include <kernel/meshopt/rumpf_smoother.hpp>
#include <kernel/meshopt/rumpf_smoother_q1hack.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1_d2.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_q1hack.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d1.hpp>
#include <kernel/meshopt/rumpf_functionals/2d_p1_d2.hpp>

using namespace FEAST;

/// \cond internal
// Helper class for resizing tests
template<typename ShapeType_>
struct helperclass;

/// \endcond

/**
 * \brief This application demonstrates the usage of some of the RumpfSmoother classes to resize reference cells
 *
 * This is excellent for error checking in smoothers and functionals, as a correct implementation should resize
 * a Rumpf reference cell to the desired volume and achieve a global minimum of the functional value.
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
  > struct ResizeApp
{
  /**
   * @brief Runs mesh smoother stuff
   *
   **/
  static void run()
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

    // Create a single reference cell of the shape type
    Geometry::ReferenceCellFactory<ShapeType, DataType> mesh_factory;
    // Create the mesh
    MeshType mesh(mesh_factory);
    // Create the transformation
    TrafoType trafo(mesh);

    // Parameters for the Rumpf functional
    DataType fac_norm = DataType(1e0),fac_det = DataType(1e0),fac_cof = DataType(0), fac_reg(DataType(0e0));
    // Create the functional with these parameters
    FunctionalType my_functional(fac_norm, fac_det, fac_cof, fac_reg);
    // This is the correct scaling for the Q1Hack with _D2 functional
    //my_functional._fac_rec_det = DataType(1);
    // Create empty unit filter of the correct size so we can create a smoother with no boundary conditions
    LAFEM::UnitFilterBlocked<MemType, DataType, IndexType, MeshType::world_dim> filter(mesh.get_num_entities(0));

    // Create the smoother
    RumpfSmootherType rumpflpumpfl(trafo, my_functional, filter);
    // Print information
    rumpflpumpfl.print();

    // Set initial coordinates by scaling the original Rumpf reference cell by ...
    DataType target_scaling(6);
    helperclass<ShapeType>::set_coords(rumpflpumpfl._coords, target_scaling);

    // Since we changed the internal _coords, they have to be copied back to the mesh
    rumpflpumpfl.set_coords();
    rumpflpumpfl.init();
    // This is the correct scaling for the Q1Hack with _D2 functional
    //Tiny::Vector<DataType, MeshType::world_dim> tmp;
    //tmp(0) = target_scaling*DataType(4)/Math::sqrt(DataType(3));
    //tmp(1) = target_scaling*DataType(2)*Math::sqrt(DataType(2))*Math::pow(DataType(2)/DataType(3), DataType(1)/DataType(4));
    //rumpflpumpfl._h(0, tmp);

    DataType scaling(2.5);
    // This transforms the unit element to the Rumpf reference element
    helperclass<ShapeType>::set_coords(rumpflpumpfl._coords, scaling);

    // Since we changed the internal _coords, they have to be copied back to the mesh
    rumpflpumpfl.set_coords();

    // Arrays for saving the contributions of the different Rumpf functional parts
    DataType* func_norm(new DataType[mesh.get_num_entities(MeshType::shape_dim)]);
    DataType* func_det(new DataType[mesh.get_num_entities(MeshType::shape_dim)]);
    DataType* func_rec_det(new DataType[mesh.get_num_entities(MeshType::shape_dim)]);

    // Compute initial functional value
    DataType fval(0);
    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
    std::cout << "fval pre optimisation = " << scientify(fval) << std::endl;

    // Compute initial functional gradient
    rumpflpumpfl.compute_gradient();

    std::string filename;
    // Write initial state to file
    filename = "pre_" + helperclass<ShapeType>::print_typename();
    Geometry::ExportVTK<MeshType> writer_initial_pre(mesh);
    writer_initial_pre.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
    writer_initial_pre.add_scalar_cell("norm_A", func_norm);
    writer_initial_pre.add_scalar_cell("det_A", func_det);
    writer_initial_pre.add_scalar_cell("rec_det_A", func_rec_det);
    writer_initial_pre.add_field_cell_blocked_vector("h", rumpflpumpfl._h);
    writer_initial_pre.write(filename);

    // Smooth the mesh
    rumpflpumpfl.optimise();

    fval = rumpflpumpfl.compute_functional(func_norm, func_det, func_rec_det);
    std::cout << "fval post optimisation = " << scientify(fval) << std::endl;

    // Call prepare() again because the mesh changed due to the optimisation and it was not called again after the
    // last iteration
    rumpflpumpfl.prepare();
    rumpflpumpfl.compute_gradient();

    // Write optimised initial mesh
    filename = "post_" + helperclass<ShapeType>::print_typename();
    Geometry::ExportVTK<MeshType> writer_initial_post(mesh);
    writer_initial_post.add_field_vertex_blocked_vector("grad", rumpflpumpfl._grad);
    writer_initial_post.add_scalar_cell("norm_A", func_norm);
    writer_initial_post.add_scalar_cell("det_A", func_det);
    writer_initial_post.add_scalar_cell("rec_det_A", func_rec_det);
    writer_initial_post.add_field_cell_blocked_vector("h", rumpflpumpfl._h);
    writer_initial_post.write(filename);

  }
};

// Template aliases to easier switch between variants

// Vanilla Rumpf smoother
template<typename A, typename B>
using MySmoother = Meshopt::RumpfSmoother<A, B>;

template<typename A, typename B>
using MyFunctional = Meshopt::RumpfFunctional_D2<A, B>;

// Using the Q1 hack
template<typename A, typename B>
using MySmootherQ1Hack = Meshopt::RumpfSmootherQ1Hack<A, B>;

// For the Q1 hack, the functional is a bit more complicated
template<typename A, typename B>
using MyFunctionalQ1Hack = Meshopt::RumpfFunctionalQ1Hack<A, B, MyFunctional>;

int main()
{
  ResizeApp<double, Shape::Hypercube<2>, MyFunctional, MySmoother>::run();
  return 0;
}

/// \brief Specialisation for hypercubes
template<int shape_dim_>
struct helperclass< FEAST::Shape::Hypercube<shape_dim_> >
{
  /// \brief Sets coordinates so we deal the the reference element
  template<typename VectorType_, typename DataType_>
  static void set_coords(VectorType_& coords_, const DataType_& scaling)
  {
    for(Index i(0); i < Index(1 << shape_dim_); ++i)
    {
      Tiny::Vector<DataType_, VectorType_::BlockSize, VectorType_::BlockSize> tmp;
      for(int d(0); d < shape_dim_; ++d)
        tmp(d) = (DataType_(((i >> d) & 1) << 1) - DataType_(1)) * scaling ;

      coords_(i, tmp);
    }
  }

  // Prints the ShapeType
  static std::string print_typename()
  {
    return "Hypercube<" + stringify(shape_dim_) +">";
  }
};

/// \brief Specialisation for 2d simplices
template<>
struct helperclass< FEAST::Shape::Simplex<2> >
{
  /// \brief Sets coordinates so we deal the the Rumpf reference element
  template<typename VectorType_, typename DataType_>
  static void set_coords(VectorType_& coords_, const DataType_& scaling)
  {
    Tiny::Vector<DataType_, 2, 2> tmp(0);
    coords_(0, tmp);

    tmp(0) = scaling;
    coords_(1, tmp);

    tmp(0) = DataType_(0.5) * scaling;
    tmp(1) = DataType_(0.5)*Math::sqrt(DataType_(3))*scaling;
    coords_(2, tmp);
  }

  // Prints the ShapeType
  static std::string print_typename()
  {
    return "Simplex<2>";
  }
};
