#pragma once
#ifndef CONTROL_SCALAR_BASIC_HPP
#define CONTROL_SCALAR_BASIC_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/comm_base.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/mean_filter.hpp>
#include <kernel/lafem/none_filter.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/transfer.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/muxer.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/transfer.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/mean_filter.hpp>
#include <kernel/assembly/symbolic_assembler.hpp>
#include <kernel/assembly/grid_transfer.hpp>
#include <kernel/assembly/mirror_assembler.hpp>
#include <kernel/assembly/mean_filter_assembler.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/iterative.hpp>
#include <kernel/util/property_map.hpp>

#include <control/domain/domain_control.hpp>

#include <deque>

namespace FEAT
{
  namespace Control
  {
    template<
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_>,
      typename TransferMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    class ScalarBasicSystemLevel
    {
    public:
      static_assert(std::is_same<MemType_, typename ScalarMatrix_::MemType>::value, "MemType mismatch!");
      static_assert(std::is_same<DataType_, typename ScalarMatrix_::DataType>::value, "DataType mismatch!");
      static_assert(std::is_same<IndexType_, typename ScalarMatrix_::IndexType>::value, "IndexType mismatch!");

      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define local scalar matrix type
      typedef ScalarMatrix_ LocalScalarMatrix;

      /// define local system matrix type
      typedef ScalarMatrix_ LocalSystemMatrix;

      /// define local transfer matrix type
      typedef TransferMatrix_ LocalSystemTransferMatrix;

      /// define local system vector type
      typedef typename LocalSystemMatrix::VectorTypeR LocalSystemVector;

      /// define local system transfer operator type
      typedef LAFEM::Transfer<LocalSystemTransferMatrix> LocalSystemTransfer;

      /// define system mirror type
      typedef LAFEM::VectorMirror<MemType_, DataType, IndexType> SystemMirror;

      /// define system gate
      typedef Global::Gate<LocalSystemVector, SystemMirror> SystemGate;

      /// define system muxer
      typedef Global::Muxer<LocalSystemVector, SystemMirror> SystemMuxer;

      /// define global system vector type
      typedef Global::Vector<LocalSystemVector, SystemMirror> GlobalSystemVector;

      /// define global system matrix type
      typedef Global::Matrix<LocalSystemMatrix, SystemMirror, SystemMirror> GlobalSystemMatrix;

      /// define global system transfer operator type
      typedef Global::Transfer<LocalSystemTransfer, SystemMirror> GlobalSystemTransfer;

      /* ***************************************************************************************** */

      /// our system gate
      SystemGate gate_sys;

      /// our coarse-level system muxer
      SystemMuxer coarse_muxer_sys;

      /// our global system matrix
      GlobalSystemMatrix matrix_sys;

      /// our global transfer operator
      GlobalSystemTransfer transfer_sys;

      /// CTOR
      ScalarBasicSystemLevel() :
        matrix_sys(&gate_sys, &gate_sys),
        transfer_sys(&coarse_muxer_sys)
      {
      }

      virtual ~ScalarBasicSystemLevel()
      {
      }

      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const ScalarBasicSystemLevel<M_, D_, I_, SM_> & other)
      {
        gate_sys.convert(other.gate_sys);
        coarse_muxer_sys.convert(other.coarse_muxer_sys);
        matrix_sys.convert(&gate_sys, &gate_sys, other.matrix_sys);
        transfer_sys.convert(&coarse_muxer_sys, other.transfer_sys);
      }

      template<typename DomainLayer_, typename DomainLevel_, typename Space_>
      void assemble_gate(const DomainLayer_& dom_layer, const DomainLevel_& dom_level, const Space_& space)
      {
        // set the gate comm
        this->gate_sys.set_comm(dom_layer.get_comm());

        // loop over all ranks
        for(Index i(0); i < dom_layer.size(); ++i)
        {
          Index rank = dom_layer.get_rank(i);

          // try to find our halo
          auto* halo = dom_level.find_halo_part(rank);
          XASSERT(halo != nullptr);

          // assemble the mirror
          SystemMirror mirror_sys;
          Assembly::MirrorAssembler::assemble_mirror(mirror_sys, space, *halo);

          // push mirror into gate
          this->gate_sys.push(int(rank), std::move(mirror_sys));
        }

        // create local template vector
        LocalSystemVector tmpl_s(space.get_num_dofs());

        // compile gate
        this->gate_sys.compile(std::move(tmpl_s));
      }

      template<typename Space_, typename Cubature_>
      void assemble_transfer(const Space_& space_fine, const Space_& space_coarse, const Cubature_& cubature)
      {
        // get local transfer operator
        LocalSystemTransfer& loc_trans = this->transfer_sys.local();

        // get local transfer matrices
        LocalSystemTransferMatrix& loc_prol = loc_trans.get_mat_prol();
        LocalSystemTransferMatrix& loc_rest = loc_trans.get_mat_rest();

        // assemble structure?
        if (loc_prol.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_2lvl(loc_prol, space_fine, space_coarse);
        }

        // create local weight vector
        LocalSystemVector loc_vec_weight = loc_prol.create_vector_l();

        // assemble prolongation matrix
        {
          loc_prol.format();
          loc_vec_weight.format();

          // assemble prolongation matrix
          Assembly::GridTransfer::assemble_prolongation(loc_prol, loc_vec_weight,
            space_fine, space_coarse, cubature);

          // synchronise weight vector using the gate
          this->gate_sys.sync_0(loc_vec_weight);

          // invert components
          loc_vec_weight.component_invert(loc_vec_weight);

          // scale prolongation matrix
          loc_prol.scale_rows(loc_prol, loc_vec_weight);

          // copy and transpose
          loc_rest = loc_prol.transpose();
        }
      }

      template<typename Space_, typename Cubature_>
      void assemble_laplace_matrix(const Space_& space, const Cubature_& cubature, const DataType nu = DataType(1))
      {
        // get local matrix
        auto& loc_matrix = this->matrix_sys.local();

        // assemble structure?
        if(loc_matrix.empty())
        {
          Assembly::SymbolicAssembler::assemble_matrix_std1(loc_matrix, space);
        }

        // format and assemble laplace
        loc_matrix.format();
        Assembly::Common::LaplaceOperator laplace_op;
        Assembly::BilinearOperatorAssembler::assemble_matrix1(loc_matrix, laplace_op, space, cubature, nu);
      }
    }; // class ScalarBasicSystemLevel<...>

    template<
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    class ScalarUnitFilterSystemLevel :
      public ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_>
    {
    public:
      typedef ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_> BaseClass;

      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define system mirror type
      typedef LAFEM::VectorMirror<MemType_, DataType, IndexType> SystemMirror;

      /// define local filter type
      typedef LAFEM::UnitFilter<MemType_, DataType_, IndexType_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter, SystemMirror> GlobalSystemFilter;

      /// Our class base type
      template <typename Mem2_, typename DT2_, typename IT2_, typename ScalarMatrix2_>
      using BaseType = class ScalarUnitFilterSystemLevel<Mem2_, DT2_, IT2_, ScalarMatrix2_>;

      /// our global system filter
      GlobalSystemFilter filter_sys;

      /// CTOR
      ScalarUnitFilterSystemLevel() :
        BaseClass()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return (*this->matrix_sys).bytes () + this->coarse_muxer_sys.bytes() + (*this->filter_sys).bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source ScalarUnitFilterSystemLevel content as content of current ScalarUnitFilterSystemLevel.
       *
       */
      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const ScalarUnitFilterSystemLevel<M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_sys.convert(other.filter_sys);
      }

      template<typename DomainLevel_, typename Space_>
      void assemble_homogeneous_unit_filter(const DomainLevel_& dom_level, const Space_& space)
      {
        auto& loc_filter = this->filter_sys.local();

        // create unit-filter assembler
        Assembly::UnitFilterAssembler<typename DomainLevel_::MeshType> unit_asm;

        std::deque<String> part_names = dom_level.get_mesh_node()->get_mesh_part_names(true);
        for(const auto& name : part_names)
        {
          auto* mesh_part_node = dom_level.get_mesh_node()->find_mesh_part_node(name);
          XASSERT(mesh_part_node != nullptr);

          // let's see if we have that mesh part
          // if it is nullptr, then our patch is not adjacent to that boundary part
          auto* mesh_part = mesh_part_node->get_mesh();
          if (mesh_part != nullptr)
          {
            // add to boundary assembler
            unit_asm.add_mesh_part(*mesh_part);
          }
        }

        // finally, assemble the filter
        unit_asm.assemble(loc_filter, space);
      }

    }; // class ScalarUnitFilterSystemLevel<...>

    template<
      typename MemType_ = Mem::Main,
      typename DataType_ = Real,
      typename IndexType_ = Index,
      typename ScalarMatrix_ = LAFEM::SparseMatrixCSR<MemType_, DataType_, IndexType_> >
    class ScalarMeanFilterSystemLevel :
      public ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_>
    {
    public:
      typedef ScalarBasicSystemLevel<MemType_, DataType_, IndexType_, ScalarMatrix_> BaseClass;

      // basic types
      typedef MemType_ MemType;
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      /// define system mirror type
      typedef LAFEM::VectorMirror<MemType_, DataType, IndexType> SystemMirror;

      /// define local filter type
      typedef Global::MeanFilter<MemType_, DataType_, IndexType_> LocalSystemFilter;

      /// define global filter type
      typedef Global::Filter<LocalSystemFilter, SystemMirror> GlobalSystemFilter;

      /// Our class base type
      template <typename Mem2_, typename DT2_, typename IT2_, typename ScalarMatrix2_>
      using BaseType = class ScalarMeanFilterSystemLevel<Mem2_, DT2_, IT2_, ScalarMatrix2_>;

      /// our global system filter
      GlobalSystemFilter filter_sys;

      /// CTOR
      ScalarMeanFilterSystemLevel() :
        BaseClass()
      {
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        return (*this->matrix_sys).bytes () + this->coarse_muxer_sys.bytes() + (*this->filter_sys).bytes();
      }

      /**
       *
       * \brief Conversion method
       *
       * Use source ScalarMeanFilterSystemLevel content as content of current ScalarMeanFilterSystemLevel.
       *
       */
      template<typename M_, typename D_, typename I_, typename SM_>
      void convert(const ScalarMeanFilterSystemLevel<M_, D_, I_, SM_> & other)
      {
        BaseClass::convert(other);
        filter_sys.convert(other.filter_sys);
      }

      template<typename Space_, typename Cubature_>
      void assemble_mean_filter(const Space_& space, const Cubature_& cubature)
      {
        // get our local system filter
        LocalSystemFilter& loc_filter = this->filter_sys.local();

        // create two global vectors
        typename BaseClass::GlobalSystemVector vec_glob_v(&this->gate_sys), vec_glob_w(&this->gate_sys);

        // fetch the local vectors
        typename BaseClass::LocalSystemVector& vec_loc_v = vec_glob_v.local();
        typename BaseClass::LocalSystemVector& vec_loc_w = vec_glob_w.local();

        // fetch the frequency vector of the system gate
        typename BaseClass::LocalSystemVector& vec_loc_f = this->gate_sys._freqs;

        // assemble the mean filter
        Assembly::MeanFilterAssembler::assemble(vec_loc_v, vec_loc_w, space, cubature);

        // synchronise the vectors
        vec_glob_v.sync_1();
        vec_glob_w.sync_0();

        // build the mean filter
        loc_filter = LocalSystemFilter(vec_loc_v.clone(), vec_loc_w.clone(), vec_loc_f.clone(), this->gate_sys.get_comm());
      }
    }; // class ScalarMeanFilterSystemLevel<...>

    template<typename Space_>
    class ScalarBasicAssemblerLevel
    {
    public:
      typedef Space_ SpaceType;
      typedef typename SpaceType::TrafoType TrafoType;
      typedef typename TrafoType::MeshType MeshType;
      typedef Control::Domain::DomainLevel<MeshType> DomainLevelType;
      typedef Control::Domain::DomainLayer<MeshType> DomainLayerType;

      DomainLevelType& domain_level;
      MeshType& mesh;
      TrafoType trafo;
      SpaceType space;

      explicit ScalarBasicAssemblerLevel(DomainLevelType& dom_lvl) :
        domain_level(dom_lvl),
        mesh(domain_level.get_mesh()),
        trafo(mesh),
        space(trafo)
      {
      }

      virtual ~ScalarBasicAssemblerLevel()
      {
      }
    }; // class ScalarBasicAssemblerLevel<...>

  } // namespace Control
} // namespace FEAT

#endif // CONTROL_SCALAR_BASIC_HPP
