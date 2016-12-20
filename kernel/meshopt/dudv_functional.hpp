#pragma once
#ifndef KERNEL_MESHOPT_DUDV_FUNCTIONAL_HPP
#define KERNEL_MESHOPT_DUDV_FUNCTIONAL_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/archs.hpp>
#include <kernel/assembly/bilinear_operator_assembler.hpp> // for BilinearOperatorAssembler
#include <kernel/assembly/common_operators.hpp>            // for DuDvOperator
#include <kernel/assembly/symbolic_assembler.hpp>          // for SymbolicMatrixAssembler
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/cubature/dynamic_factory.hpp>             // for DynamicFactory
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/filter_chain.hpp>
#include <kernel/lafem/filter_sequence.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/meshopt/mesh_quality_functional.hpp>
#include <kernel/space/lagrange1/element.hpp>

namespace FEAT
{
  namespace Meshopt
  {
    /**
     * \brief Mesh optimiser based on minimisation of harmonic energy
     *
     * \tparam Mem_
     * Memory architecture for the solver (not the mesh)
     *
     * \tparam DT_
     * Data type for the solver (not the mesh)
     *
     * \tparam IT_
     * Index type for the solver (not the mesh)
     *
     * \tparam TrafoType_
     * Type of the underlying transformation.
     *
     * \author Jordi Paul
     *
     */
    template
    <
      typename Mem_, typename DT_, typename IT_, typename TrafoType_,
      template<typename, typename, typename, int, int> class MatrixType_ = LAFEM::SparseMatrixBCSR
    >
    class DuDvFunctional:
      public MeshQualityFunctional<typename TrafoType_::MeshType>
    {
      public:
        /// Memory architecture
        typedef Mem_ MemType;
        /// Our datatype
        typedef DT_ DataType;
        /// Our index type
        typedef IT_ IndexType;
        /// Type for the transformation
        typedef TrafoType_ TrafoType;
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// The shape type of the mesh
        typedef typename MeshType::ShapeType ShapeType;

        /// Type for the system matrix
        typedef MatrixType_<Mem_, DT_, IT_, MeshType::world_dim, MeshType::world_dim> MatrixType;
        /// Blockheight of the system matrix
        static constexpr int BlockHeight = MatrixType::BlockHeight;
        /// Blockwidth of the system matrix
        static constexpr int BlockWidth = MatrixType::BlockWidth;

        /// Our 'base' class type
        template <typename Mem2_, typename DT2_ = DT_, typename IT2_ = IT_>
        using ContainerType = DuDvFunctional<Mem2_, DT2_, IT2_, TrafoType_, MatrixType_>;

        /// this typedef lets you create a matrix container with new Memory, Datatape and Index types
        template <typename Mem2_, typename DataType2_, typename IndexType2_>
        using ContainerTypeByMDI = ContainerType<Mem2_, DataType2_, IndexType2_>;

        /// Our base class
        typedef MeshQualityFunctional<MeshType> BaseClass;
        /// Type for exchanging information between state variable and mesh
        typedef typename BaseClass::CoordsBufferType CoordsBufferType;

        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;

        //template<typename A, typename B, typename C>
        //using MatrixTemplate = MatrixType_<A, B, C, MeshType::world_dim, MeshType::world_dim>;

        /// Type for vectors from the dual space
        typedef typename MatrixType::VectorTypeL VectorTypeL;
        /// Type for vectors from the primal space
        typedef typename MatrixType::VectorTypeR VectorTypeR;
        /// Type for i.e. cell vectors
        typedef LAFEM::DenseVector<Mem_, DT_, IT_> ScalarVectorType;

        /// Filter for Dirichlet boundary conditions
        typedef LAFEM::UnitFilterBlocked<Mem_, DT_, IT_, MeshType::world_dim> DirichletFilterType;
        /// Sequence of Dirichlet filters for several different boundary parts
        typedef LAFEM::FilterSequence<DirichletFilterType> DirichletFilterSequence;
        /// Filter for slip boundary conditions
        typedef LAFEM::SlipFilter<Mem_, DT_, IT_, MeshType::world_dim> SlipFilterType;
        /// Sequence of Slip filters for several different boundary parts
        typedef LAFEM::FilterSequence<SlipFilterType> SlipFilterSequence;
        /// Combined filter
        typedef LAFEM::FilterChain<SlipFilterSequence, DirichletFilterSequence> FilterType;

        /// Finite Element space for the transformation
        typedef typename Intern::TrafoFE<TrafoType>::Space SpaceType;
        /// Maximum polynomial degree
        // 2 * (degree of trafo) for both trial and test spaces, +1 for safety reasons
        // This could be decreased by the degree of the operator, i.e. 2 for Du:Dv
        static constexpr int _local_degree = 4*SpaceType::local_degree + 1;

        /// The system matrix
        MatrixType matrix_sys;

      protected:
        /// The transformation defining the physical mesh
        TrafoType& _trafo;
        /// Vector saving the cell sizes on the reference mesh
        ScalarVectorType _lambda;
        /// Assembler for Dirichlet boundary conditions
        std::map<String, std::shared_ptr<Assembly::UnitFilterAssembler<MeshType>>> _dirichlet_asm;
        /// Assembler for slip boundary conditions
        std::map<String, std::shared_ptr<Assembly::SlipFilterAssembler<MeshType>>> _slip_asm;
        /// Cubature factory, for P1/Q1 transformations in 2d degree 5 is enough
        Cubature::DynamicFactory _cubature_factory;

      public:
        /// The FE space for the transformation, needed for filtering
        SpaceType trafo_space;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] rmn_
         * The RootMeshNode representing the tree of root mesh, all of its MeshParts and Charts
         *
         * \param[in] trafo_space_
         * The FE space the transformation lives in.
         *
         * \param[in] dirichlet_asm_
         * The set of unit filter assemblers.
         *
         * \param[in] slip_asm_
         * The set of slip filter assemblers.
         *
         */
        explicit DuDvFunctional(
          Geometry::RootMeshNode<MeshType>* rmn_,
          TrafoType& trafo,
          const std::deque<String>& dirichlet_list, const std::deque<String>& slip_list):
          BaseClass(rmn_),
          matrix_sys(),
          _trafo(trafo),
          _lambda(rmn_->get_mesh()->get_num_entities(MeshType::shape_dim)),
          _dirichlet_asm(),
          _slip_asm(),
          _cubature_factory("auto-degree:"+stringify(int(_local_degree))),
          trafo_space(trafo)
          {
            // The symbolic assembly has to happen in the constructor because the matrix is shallow copied immediately
            // after construction in the Control::QuadraticSystemLevel
            Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_sys, trafo_space);

            // For every MeshPart specified in the list of Dirichlet boundaries, create a UnitFilterAssembler and
            // insert it into the std::map
            for(auto& it : dirichlet_list)
            {
              // Create empty assembler
              auto new_asm = std::make_shared<Assembly::UnitFilterAssembler<MeshType>>();

              // Add the MeshPart to the assembler if it is there. There are legimate reasons for it NOT to be
              // there, i.e. we are in parallel and our patch is not adjacent to that MeshPart
              auto* mpp = rmn_->find_mesh_part(it);
              if(mpp != nullptr)
              {
                new_asm->add_mesh_part(*mpp);
              }

              // Insert into the map
              String identifier(it);
              _dirichlet_asm.emplace(identifier, new_asm);
            }

            // For every MeshPart specified in the list of slip boundaries, create a SlipFilterAssembler and
            // insert it into the std::map
            for(auto& it : slip_list)
            {
              // Create empty assembler
              auto new_asm = std::make_shared<Assembly::SlipFilterAssembler<MeshType>>(*this->get_mesh());

              // Add the MeshPart to the assembler if it is there. There are legimate reasons for it NOT to be
              // there, i.e. we are in parallel and our patch is not adjacent to that MeshPart
              auto* mpp = rmn_->find_mesh_part(it);
              if(mpp != nullptr)
              {
                new_asm->add_mesh_part(*mpp);
              }

              // Insert into the map
              String identifier(it);
              _slip_asm.emplace(identifier, new_asm);
            }

          }

        /// Explicitly delete copy constructor
        DuDvFunctional(const DuDvFunctional&) = delete;

        /// \brief Virtual destructor
        virtual ~DuDvFunctional()
        {
        }

        /**
         * \brief Performs one-time initialisations
         *
         * This is not done in the constructor for the case that the system matrix gets overwritten by a derived
         * class, so the unused system matrix of THIS class is not assembled symbolically
         */
        virtual void init() override
        {
          //Assembly::SymbolicAssembler::assemble_matrix_std1(matrix_sys, trafo_space);
          assemble_system_matrix();
        }

        /// \copydoc BaseClass::name()
        static String name()
        {
          return "DuDvFunctional<"+MeshType::name()+">";
        }

        /**
         * \brief Adds relevant quantities of this object to a VTK exporter
         *
         * \param[in, out] exporter
         * The exporter to add our data to.
         *
         */
        virtual void add_to_vtk_exporter(Geometry::ExportVTK<MeshType>& exporter) const override
        {
          exporter.add_cell_scalar("lambda", this->_lambda.elements());
        }

        /**
         * \brief Assembles the system matrix
         */
        virtual void assemble_system_matrix()
        {
          matrix_sys.format();

          Assembly::Common::DuDvOperatorBlocked<MeshType::world_dim> my_operator;

          Assembly::BilinearOperatorAssembler::assemble_block_matrix1(
            matrix_sys,           // the matrix that receives the assembled operator
            my_operator, // the operator that is to be assembled
            trafo_space,            // the finite element space in use
            _cubature_factory  // the cubature factory to be used for integration
            );
        }

        /**
         * \brief Prepares the functional for evaluation.
         *
         * \param[in, out] vec_state
         * The state vector. This is not const because it might be filtered using a slip filter.
         *
         * \param[in, out] filter
         * The filter. This is not const because the slip filter might need to be reassembled.
         *
         * Needs to be called whenever any data like the mesh etc. changed.
         *
         */
        virtual void prepare(VectorTypeR& vec_state, FilterType& filter)
        {
          for(Index cell(0); cell < this->get_mesh()->get_num_entities(MeshType::shape_dim); ++cell)
          {
            _lambda(cell, DataType(_trafo.template compute_vol<typename MeshType::ShapeType>(cell)));
          }

          auto& dirichlet_filters = filter.template at<1>();

          for(auto& it : dirichlet_filters)
          {
            const auto& assembler = _dirichlet_asm.find(it.first);
            if(assembler == _dirichlet_asm.end())
            {
              throw InternalError(__func__,__FILE__,__LINE__,
              "Could not find dirichlet assembler for filter with key "+it.first);
            }

            assembler->second->assemble(it.second, trafo_space, vec_state);
          }

          // The slip filter contains the outer unit normal, so reassemble it
          auto& slip_filters = filter.template at<0>();

          for(auto& it : slip_filters)
          {
            const auto& assembler = _slip_asm.find(it.first);
            if(assembler == _slip_asm.end())
            {
              throw InternalError(__func__,__FILE__,__LINE__,
              "Could not find slip filter assembler for filter with key "+it.first);
            }

            assembler->second->assemble(it.second, trafo_space);
          }

          for(const auto& it:slip_filters)
          {
            this->_mesh_node->adapt_by_name(it.first);
          }

        }

        /**
         * \brief Computes a quality indicator concerning the cell sizes
         *
         * \param[out] lambda_min
         * Minimum of the optimal cell size lambda over all cells
         *
         * \param[out] lambda_max
         * Maximum of the optimal cell size lambda over all cells
         *
         * \param[out] vol_min
         * Minimum cell volume
         *
         * \param[out] vol_max
         * Maximum cell volume
         *
         * \param[out] vol
         * Total volume of the domain given by the mesh
         *
         * In a truly optimal mesh (consisting ONLY of reference cells of the right size), every cell's volume is
         * exactly lambda(cell). This is especially the goal for r-adaptivity.
         * So in an optimal mesh,
         * \f[
         *   \forall K \in \mathcal{T}_h: |K|/|\Omega| = \lambda(K)
         * \f]
         * so we compute the 1-norm of the vector
         * \f$(v)_i = \left| \frac{|K_i|}{\sum_j |K_j|} - \lambda(K_i) \right| \f$.
         *
         * \returns The relative cell size quality indicator.
         *
         * \note lambda_min, lambda_max, vol_min, and vol_max are all volume fractions.
         *
         * \note As these quantities are global, this function has a pre_sync and post_sync part.
         *
         * \see Control::DuDvFunctionalControl::compute_cell_size_defect()
         *
         */
        virtual CoordType compute_cell_size_defect(CoordType& lambda_min, CoordType& lambda_max, CoordType& vol_min, CoordType& vol_max, CoordType& vol) const override
        {

          compute_cell_size_defect_pre_sync(vol_min, vol_max, vol);
          return compute_cell_size_defect_post_sync(lambda_min, lambda_max, vol_min, vol_max, vol);

        }

        /**
         * \brief Computes a quality indicator concerning the cell sizes, pre synchronisation phase
         *
         * \param[out] vol_min
         * Minimum cell volume
         *
         * \param[out] vol_max
         * Maximum cell volume
         *
         * \param[out] vol
         * Total volume of the domain given by the mesh
         *
         * \returns The relative cell size quality indicator.
         *
         * \note lambda_min, lambda_max, vol_min, and vol_max are all volume fractions.
         *
         * \note As these quantities are global, this function has a pre_sync and post_sync part.
         *
         * \see Control::DuDvFunctionalControl::compute_cell_size_defect()
         *
         */
        void compute_cell_size_defect_pre_sync(CoordType& vol_min, CoordType& vol_max, CoordType& vol) const
        {
          vol = CoordType(0);

          vol_min = Math::huge<CoordType>();
          vol_max = CoordType(0);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            CoordType my_vol = this->_trafo.template compute_vol<ShapeType, CoordType>(cell);
            vol_min = Math::min(vol_min, my_vol);
            vol_max = Math::max(vol_min, my_vol);
            vol += my_vol;
          }
        }

        /**
         * \brief Computes a quality indicator concerning the cell sizes, pre synchronisation phase
         *
         * \param[out] lambda_min
         * Minimum of the optimal cell size lambda over all cells
         *
         * \param[out] lambda_max
         * Maximum of the optimal cell size lambda over all cells
         *
         * \param[out] vol_min
         * Minimum cell volume
         *
         * \param[out] vol_max
         * Maximum cell volume
         *
         * \param[out] vol
         * Total volume of the domain given by the mesh
         *
         * \returns The relative cell size quality indicator.
         *
         * \note lambda_min, lambda_max, vol_min, and vol_max are all volume fractions.
         *
         * \note As these quantities are global, this function has a pre_sync and post_sync part.
         *
         * \see Control::DuDvFunctionalControl::compute_cell_size_defect()
         *
         */
        virtual CoordType compute_cell_size_defect_post_sync(CoordType& lambda_min, CoordType& lambda_max, CoordType& vol_min, CoordType& vol_max, const CoordType& vol) const
        {
          CoordType size_defect(0);
          lambda_min = Math::huge<CoordType>();
          lambda_max = CoordType(0);

          for(Index cell(0); cell < this->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            size_defect += Math::abs(this->_trafo.template compute_vol<ShapeType, CoordType>(cell)/vol - this->_lambda(cell));
            lambda_min = Math::min(lambda_min, CoordType(this->_lambda(cell)));
            lambda_max = Math::max(lambda_max, CoordType(this->_lambda(cell)));
          }

          vol_min /= vol;
          vol_max/= vol;

          return size_defect;
        }

        /**
         * \brief Checks if the functional is empty (= the null functional
         *
         * \returns True if the number of DoFs is zero.)
         */
        bool empty() const
        {
          return matrix_sys.empty();
        }

        /**
         * \brief Creates an L-vector for the functional's gradient
         */
        VectorTypeL create_vector_l() const
        {
          return VectorTypeL(trafo_space.get_num_dofs());
        }

        /**
         * \brief Creates an R-vector for the functional and its gradient
         */
        VectorTypeR create_vector_r() const
        {
          return VectorTypeR(trafo_space.get_num_dofs());
        }


    }; // class DuDvFunctional

#ifdef FEAT_EICKT
    extern template class DuDvFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Simplex<2>, 2, 2, double>>,
      LAFEM::SparseMatrixBCSR
    >;

    extern template class DuDvFunctional
    <
      Mem::Main, double, Index,
      Trafo::Standard::Mapping<Geometry::ConformalMesh<Shape::Hypercube<2>, 2, 2, double>>,
      LAFEM::SparseMatrixBCSR
    >;
#endif // FEAT_EICKT
  } // namespace Meshopt
} // namespace FEAT
#endif // KERNEL_MESHOPT_DUDV_FUNCTIONAL_HPP
