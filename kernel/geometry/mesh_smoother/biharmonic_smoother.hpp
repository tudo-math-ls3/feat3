#pragma once
#ifndef KERNEL_GEOMETRY_BIHARMONIC_SMOOTHER_HPP
#define KERNEL_GEOMETRY_BIHARMONIC_SMOOTHER_HPP 1

#include <kernel/archs.hpp>

#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/common_operators.hpp>            // for LaplaceOperator
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicMatrixAssembler
#include <kernel/assembly/dirichlet_assembler.hpp>         // for DirichletAssembler
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/lafem/bicgstab.hpp>
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
    template<typename TrafoType_, typename DataType_, typename MemType_>
    class BiharmonicSmoother :
      public MeshSmoother<TrafoType_, DataType_, MemType_>
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
        typedef MeshSmoother<TrafoType, DataType, MemType> BaseClass;

        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;

        /// FE space the transformation lives in
        typedef Space::Lagrange1::Element<TrafoType> SpaceType;

        /// Type for the vectors to hold coordinates etc.
        typedef LAFEM::DenseVector<MemType, DataType> SubVectorType;
        typedef LAFEM::SparseMatrixCSR<MemType, DataType> SubMatrixType;

        typedef LAFEM::NoneFilter<MemType, DataType, IndexType> SubFilterType0;
        typedef LAFEM::UnitFilter<MemType, DataType, IndexType> SubFilterType1;

        typedef LAFEM::TupleFilter<SubFilterType0, SubFilterType1> FilterType;
        typedef LAFEM::TupleVector<SubVectorType, SubVectorType> VectorType;
        typedef LAFEM::SaddlePointMatrix<SubMatrixType> MatrixType;

        typedef Algo::Generic AlgoType;

      protected:
        SpaceType _trafo_space;
        MatrixType _sys_matrix;
        VectorType _vec_rhs;
        Cubature::DynamicFactory _cubature_factory;
        FilterType _filter;
        Assembly::DirichletAssembler<SpaceType> _dirichlet_asm;

      public:
        /// Constructor
        explicit BiharmonicSmoother(TrafoType& trafo_) :
          BaseClass(trafo_),
          _trafo_space(trafo_),
          _sys_matrix(),
          _vec_rhs(),
          _cubature_factory("auto-degree:5"),
          _filter(),
          _dirichlet_asm(_trafo_space)
          {
            /// Type for the boundary mesh
            typedef typename Geometry::CellSubSet<ShapeType> BoundaryType;
            /// Factory for the boundary mesh
            typedef typename Geometry::BoundaryFactory<MeshType> BoundaryFactoryType;

            // Get the boundary set
            BoundaryFactoryType boundary_factory(this->_mesh);
            BoundaryType boundary(boundary_factory);

            //_filter.template at<0>() = std::move(SubFilterType0(this->_nk));
            _filter.template at<1>() = std::move(SubFilterType1(this->_nk));

            _vec_rhs.template at<0>() = std::move(SubVectorType(this->_nk));
            _vec_rhs.template at<1>() = std::move(SubVectorType(this->_nk));

            _dirichlet_asm.add_cell_set(boundary);

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
          prepare();
          // Create a dummy preconditioner
          LAFEM::NonePreconditioner<AlgoType, MatrixType, VectorType> precond;

          for(Index d(0); d < this->_world_dim; ++d)
          {
            VectorType v(SubVectorType(this->_nk,DataType(0)), this->_coords[d].clone());
            //_dirichlet_asm.assemble(_filter.template at<0>(),SubVectorType(this->_nk,DataType(0)));
            _dirichlet_asm.assemble(_filter.template at<1>(),this->_coords[d]);

            // DEBUG: Dump filters
            // _filter.template at<0>()._sv.write_out_mtx("filter0.mtx");
            // _filter.template at<1>()._sv.write_out_mtx("filter1.mtx");

            _vec_rhs.template at<0>().format();
            _vec_rhs.template at<1>().format();

            _filter.template filter_rhs<AlgoType>(_vec_rhs);
            _filter.template filter_sol<AlgoType>(v);

            // Fire up the BiCGStab solver
            LAFEM::BiCGStab<AlgoType>::value(
              v,    // the initial solution vector
              _sys_matrix,     // the system matrix
              _vec_rhs,    // the right-hand-side vector
              _filter,
              precond,    // the dummy preconditioner
              1000,        // maximum number of iterations
              1E-8        // relative tolerance to achieve
              );

            // DEBUG: Dump solution
            // v.template at<0>().write_out_mtx("w.mtx");
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
