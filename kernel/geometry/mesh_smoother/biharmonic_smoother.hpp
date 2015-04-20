#pragma once
#ifndef KERNEL_GEOMETRY_BIHARMONIC_SMOOTHER_HPP
#define KERNEL_GEOMETRY_BIHARMONIC_SMOOTHER_HPP 1

#include <kernel/archs.hpp>

#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicMatrixAssembler
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/lafem/proto_solver.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/preconditioner.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/tuple_filter.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/geometry/mesh_smoother/mesh_smoother.hpp>
#include <kernel/space/lagrange1/element.hpp>

namespace FEAST
{
  namespace Geometry
  {

    /**
     * \brief Mesh optimiser based on minimisation of biharmonic energy
     *
     * \tparam TrafoType_
     * Type of the underlying transformation.
     *
     * \tparam DataType_
     * Our datatype.
     *
     * \tparam MemType_
     * Memory architecture.
     *
     * \author Jordi Paul
     *
     */
    template<typename DataType_, typename MemType_, typename TrafoType_>
    class BiharmonicSmoother :
      public MeshSmoother<DataType_, MemType_, TrafoType_>
    {
      public:
        /// Type for the transformation
        typedef TrafoType_ TrafoType;
        /// Our datatype
        typedef DataType_ DataType;
        /// Memory architecture
        typedef MemType_ MemType;
        /// Index type, hardwired for now
        typedef Index IndexType;

        /// Our base class
        typedef MeshSmoother<DataType, MemType, TrafoType> BaseClass;

        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;

        /// FE space the transformation lives in
        typedef Space::Lagrange1::Element<TrafoType> SpaceType;

        /// Type for the vectors to hold coordinates etc.
        typedef LAFEM::DenseVector<MemType, DataType> SubVectorType;
        typedef LAFEM::SparseMatrixCSR<MemType, DataType> SubMatrixType;

        // Chose your poison: UnitFilter means bogus UnitFilter values for the Laplacian, if they are zero then the
        // whole method is equivalent to just using the Laplace smoother.
        // Or use a filter for Neuman BVs for the Laplacian. Unfortunately, only homogeneous Neumann BVs are
        // implemented now, and this does not make any sense for the mesh coordinate distribution point of view.
        typedef LAFEM::NoneFilter<MemType, DataType, IndexType> SubFilterType0;
        typedef LAFEM::UnitFilter<MemType, DataType, IndexType> SubFilterType1;

        typedef LAFEM::TupleFilter<SubFilterType0, SubFilterType1> FilterType;
        typedef LAFEM::TupleVector<SubVectorType, SubVectorType> VectorType;
        typedef LAFEM::SaddlePointMatrix<SubMatrixType> MatrixType;

      protected:
        SpaceType _trafo_space;
        MatrixType _sys_matrix;
        VectorType _vec_rhs;
        Cubature::DynamicFactory _cubature_factory;
        FilterType _filter;
        Assembly::UnitFilterAssembler<MeshType> _dirichlet_asm;

      public:
        /// Constructor
        explicit BiharmonicSmoother(TrafoType& trafo_) :
          BaseClass(trafo_),
          _trafo_space(trafo_),
          _sys_matrix(),
          _vec_rhs(),
          _cubature_factory("auto-degree:5"),
          _filter(),
          _dirichlet_asm()
          {
            /// Type for the boundary mesh
            typedef typename Geometry::MeshPart<MeshType> BoundaryType;
            /// Factory for the boundary mesh
            typedef typename Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;

            // Total number of vertices in the mesh
            Index nvertices(this->_mesh.get_num_entities(0));

            // Get the boundary set
            BoundaryFactoryType boundary_factory(this->_mesh);
            BoundaryType boundary(boundary_factory);

            // This does not have to be set if we use a NoneFilter
            //_filter.template at<0>() = std::move(SubFilterType0(nvertices));
            _filter.template at<1>() = std::move(SubFilterType1(nvertices));

            _vec_rhs.template at<0>() = std::move(SubVectorType(nvertices));
            _vec_rhs.template at<1>() = std::move(SubVectorType(nvertices));

            _dirichlet_asm.add_mesh_part(boundary);

          }

        virtual ~BiharmonicSmoother()
        {
        }

        /// \brief Initialises parts of the MeshSmoother not set in in the constructor
        virtual void init() override
        {
          BaseClass::init();
          Assembly::SymbolicMatrixAssembler<>::assemble1(_sys_matrix.block_a(), _trafo_space);
          _sys_matrix.block_b() = _sys_matrix.block_a().clone(LAFEM::CloneMode::Weak);
          _sys_matrix.block_d() = _sys_matrix.block_a().clone(LAFEM::CloneMode::Weak);
        }

        /// \brief Optimises the mesh according to the criteria implemented in the mesh smoother.
        virtual void optimise() override
        {
          // Total number of vertices in the mesh
          Index nvertices(this->_mesh.get_num_entities(0));

          prepare();
          // Create a SSOR preconditioner
          // No usable preconditioner for SaddlePointMatrix yet, so pass a nullptr to use NO preconditioner
          //LAFEM::PreconWrapper<MatrixType, LAFEM::SSORPreconditioner> precond(_sys_matrix);
          // Create a BiCGStab solver
          LAFEM::BiCGStabSolver<MatrixType, FilterType> solver(_sys_matrix, _filter, nullptr);// &precond);
          // Enable convergence plot
          solver.set_plot(false);
          solver.set_max_iter(5000);
          // Initialise the solver
          solver.init();

          prepare();

          for(int d(0); d < MeshType::world_dim; ++d)
          {
            SubVectorType tmp(nvertices,DataType(0));
            VectorType v(SubVectorType(nvertices,DataType(0)), this->_coords[d].clone());
            // This does not need to be assembled if we use a NoneFilter
            //_dirichlet_asm.assemble(_filter.template at<0>(), _trafo_space,tmp);
            _dirichlet_asm.assemble(_filter.template at<1>(), _trafo_space, this->_coords[d]);

            _vec_rhs.template at<0>().format();
            _vec_rhs.template at<1>().format();

            _filter.filter_rhs(_vec_rhs);
            _filter.filter_sol(v);

            // Correct our initial solution vector
            solver.correct(v, _vec_rhs);

            this->_coords[d] = std::move(v.template at<1>());
          }
          this->set_coords();

          return;
        }

        /// \copydoc MeshSmoother::prepare()
        virtual void prepare() override
        {
          //_sys_matrix.format();
          //Assembly::Common::LaplaceOperator laplace_operator;
          //Assembly::BilinearOperatorAssembler::assemble_matrix1(
          //  _sys_matrix,           // the matrix that receives the assembled operator
          //  laplace_operator, // the operator that is to be assembled
          //  _trafo_space,            // the finite element space in use
          //  _cubature_factory  // the cubature factory to be used for integration
          //  );

          _sys_matrix.block_a().format();
          Assembly::Common::IdentityOperator identity_operator;
          Assembly::BilinearOperatorAssembler::assemble_matrix1(
            _sys_matrix.block_a(),           // the matrix that receives the assembled operator
            identity_operator, // the operator that is to be assembled
            _trafo_space,            // the finite element space in use
            _cubature_factory,  // the cubature factory to be used for integration
            DataType(-1)
            );

          _sys_matrix.block_b().format();
          Assembly::Common::LaplaceOperator laplace_operator;
          Assembly::BilinearOperatorAssembler::assemble_matrix1(
            _sys_matrix.block_b(),           // the matrix that receives the assembled operator
            laplace_operator, // the operator that is to be assembled
            _trafo_space,            // the finite element space in use
            _cubature_factory  // the cubature factory to be used for integration
            );

          _sys_matrix.block_d() = _sys_matrix.block_b().shared();

          // DEBUG: Dump matrices
          //_sys_matrix.block_a().write_out_mtx("A.mtx");
          //_sys_matrix.block_b().write_out_mtx("B.mtx");
          //_sys_matrix.block_d().write_out_mtx("D.mtx");

          return;
        }

    }; // class BiharmonicSmoother

  } // namespace Geometry
} // namespace FEAST
#endif // KERNEL_GEOMETRY_BIHARMONIC_SMOOTHER_HPP
