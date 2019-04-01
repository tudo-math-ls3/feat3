// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_AMG_HPP
#define KERNEL_SOLVER_AMG_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/transfer.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief AMG Factory, creating algebraic grid hierarchies
     *
     * This factory is able to create coarser level based on a fine input grid
     * and to update the non zero values of an existing coarse grid, while
     * reusing the transfer operator and matrix layout.
     *
     * The principal implemention is based on \cite Metsch:2013
     *
     * The page \ref amg elaborates on how to integrate this into your own application.
     *
     * \todo Truncation of Interpolation (2.8.9)
     * \todo Aggressive coarsening (2.7)&& Multipass Interpolation (2.8.7)
     * \todo Extended Interpolation (2.8.6
     * \todo Direct +- Interpolation (2.8.3) eq (2.45)
     * \todo (2.12.3) eq (2.109) block interpolation && block smoothing
     * \todo SaddlePoint AMG (4)
     */
    template <typename SystemMatrix_, typename SystemFilter_, typename TransferOperator_>
    class AMGFactory
    {
      typedef typename SystemMatrix_::MemType MemType;
      typedef typename SystemMatrix_::DataType DataType;
      typedef typename SystemMatrix_::IndexType IndexType;

    private:
      template <typename SM_, typename = typename std::enable_if<SM_::is_global>::type >
      static typename SM_::LocalMatrix get_local_matrix_type()
      {
        typename SM_::LocalMatrix dummy;
        return dummy;
      }

      template <typename SM_, typename = typename std::enable_if<SM_::is_local>::type >
      static SM_ get_local_matrix_type()
      {
        SM_ dummy;
        return dummy;
      }

      template <typename SF_, typename = typename std::enable_if<SF_::is_global>::type >
      static typename SF_::LocalFilter get_local_filter_type()
      {
        typename SF_::LocalFilter dummy;
        return dummy;
      }

      template <typename SF_, typename = typename std::enable_if<SF_::is_local>::type >
      static SF_ get_local_filter_type()
      {
        SF_ dummy;
        return dummy;
      }


      template <typename TO_, typename = typename std::enable_if<TO_::is_local>::type >
      static TO_ get_local_transfer_operator_type()
      {
        TO_ dummy;
        return dummy;
      }

      template <typename TO_, typename = typename std::enable_if<TO_::is_global>::type >
      static typename TO_::LocalTransfer get_local_transfer_operator_type()
      {
        typename TO_::LocalFilter dummy;
        return dummy;
      }

      typedef LAFEM::VectorMirror<typename SystemMatrix_::MemType, DataType, IndexType> MirrorType;

      template <typename SM_, typename = typename std::enable_if<SM_::is_global>::type >
      static SM_ _get_global_matrix_type()
      {
        SM_ dummy;
        return dummy;
      }

      template <typename SM_, typename = typename std::enable_if<SM_::is_local>::type >
      static Global::Matrix<SM_,
        LAFEM::VectorMirror<typename SystemMatrix_::MemType, DataType, IndexType>,
        LAFEM::VectorMirror<typename SystemMatrix_::MemType, DataType, IndexType>
        > _get_global_matrix_type()
      {
        Global::Matrix<SM_,
        LAFEM::VectorMirror<typename SystemMatrix_::MemType, DataType, IndexType>,
        LAFEM::VectorMirror<typename SystemMatrix_::MemType, DataType, IndexType>
        > dummy;
        return dummy;
      }

      template <typename SF_, typename = typename std::enable_if<SF_::is_global>::type >
      static SF_ _get_global_filter_type()
      {
        SF_ dummy;
        return dummy;
      }

      template <typename SF_, typename = typename std::enable_if<SF_::is_local>::type >
      static Global::Filter<SF_, MirrorType> _get_global_filter_type()
      {
        Global::Filter<SF_, MirrorType> dummy;
        return dummy;
      }

      template <typename ST_, typename = typename std::enable_if<ST_::is_global>::type >
      static ST_ _get_global_transfer_operator_type()
      {
        ST_ dummy;
        return dummy;
      }

      template <typename ST_, typename = typename std::enable_if<ST_::is_local>::type >
      static Global::Transfer<ST_, MirrorType> _get_global_transfer_operator_type()
      {
        Global::Transfer<ST_, MirrorType> dummy;
        return dummy;
      }

      using LocalMatrixType = decltype(get_local_matrix_type<SystemMatrix_>());
      using LocalFilterType = decltype(get_local_filter_type<SystemFilter_>());
      using LocalTransferOperatorType = decltype(get_local_transfer_operator_type<TransferOperator_>());
      typedef Global::Gate<typename LocalMatrixType::VectorTypeL, MirrorType> GateRowType;
      using GlobalMatrixType = decltype(_get_global_matrix_type<SystemMatrix_>());
      using GlobalFilterType = decltype(_get_global_filter_type<SystemFilter_>());
      using GlobalTransferOperatorType = decltype(_get_global_transfer_operator_type<TransferOperator_>());


      template <typename MT_, typename = typename std::enable_if<MT_::is_global>::type >
      static MT_ _make_matrix_global(const MT_ & input)
      {
        return input.clone(LAFEM::CloneMode::Shallow);
      }

      template <typename MT_, typename = typename std::enable_if<MT_::is_local>::type >
      static Global::Matrix<MT_, MirrorType, MirrorType> _make_matrix_global(const MT_ & input)
      {
        Global::Matrix<MT_, MirrorType, MirrorType> matrix(nullptr, nullptr, input.clone(LAFEM::CloneMode::Shallow));
        return matrix;
      }

      template <typename FT_, typename = typename std::enable_if<FT_::is_global>::type >
      static FT_ _make_filter_global(const FT_ & input)
      {
        return input.clone(LAFEM::CloneMode::Shallow);
      }

      template <typename FT_, typename = typename std::enable_if<FT_::is_local>::type >
      static Global::Filter<FT_, MirrorType> _make_filter_global(const FT_ & input)
      {
        Global::Filter<FT_, MirrorType> filter(input.clone(LAFEM::CloneMode::Shallow));
        return filter;
      }

      template <typename TT_, typename = typename std::enable_if<TT_::is_global>::type >
      static TT_ _make_transfer_global(const TT_ & input)
      {
        return input.clone(LAFEM::CloneMode::Shallow);
      }

      template <typename TT_, typename = typename std::enable_if<TT_::is_local>::type >
      static Global::Transfer<TT_, MirrorType> _make_transfer_global(const TT_ & input)
      {
        Global::Transfer<TT_, MirrorType> transfer(nullptr, input.clone(LAFEM::CloneMode::Shallow));
        return transfer;
      }

      template <typename TargetType_, typename InputType_>
      static void _adaptive_clone(TargetType_ & output, const InputType_ & input,
        typename std::enable_if<InputType_::is_global>::type* = nullptr,
        typename std::enable_if<TargetType_::is_global>::type* = nullptr)
      {
        output = input.clone(LAFEM::CloneMode::Shallow);
      }

      template <typename TargetType_, typename InputType_>
      static void _adaptive_clone(TargetType_ & output, const InputType_ & input,
        typename std::enable_if<InputType_::is_local>::type* = nullptr,
        typename std::enable_if<TargetType_::is_local>::type* = nullptr)
      {
        output = input.clone(LAFEM::CloneMode::Shallow);
      }

      template <typename TargetType_, typename InputType_>
      static void _adaptive_clone(TargetType_ & output, const InputType_ & input,
        typename std::enable_if<InputType_::is_global>::type* = nullptr,
        typename std::enable_if<TargetType_::is_local>::type* = nullptr)
      {
        output = input.local().clone(LAFEM::CloneMode::Shallow);
      }

      public:
      /**
       * \brief Creates a new coarse grid
       *
       * This function is used to create a new coarse grid level, based on a provided fine grid level.
       *
       * \param[in] fine_grid_input The fine grid matrix
       * \param[in] fine_filter_input The fine grid filter
       * \param[in] comm The communicator to use
       * \param[in] theta Dropping value for coupling computation
       * \param[out] coarse_grid_out The new coarse grid matrix
       * \param[out] coarse_filter_out The new coarse grid filter
       * \param[out] transfer_out The transfer operator between the old fine and the new coarse grid level
       * \param[out] gate_out The new coarse grid gate
       */
        static void new_coarse_level(const SystemMatrix_ & fine_grid_input, const SystemFilter_ & fine_filter_input,  double theta,
            SystemMatrix_ & coarse_grid_out, SystemFilter_ & coarse_filter_out, TransferOperator_ & transfer_out,
            const Dist::Comm * comm, typename GlobalMatrixType::GateRowType * gate_out = nullptr)
        {
          auto gate_coarse = std::make_shared<typename GlobalMatrixType::GateRowType>();
          gate_coarse->set_comm(comm);
          const GlobalMatrixType matrix_fine(_make_matrix_global<SystemMatrix_>(fine_grid_input));
          const GlobalFilterType filter_fine(_make_filter_global<SystemFilter_>(fine_filter_input));

          typename LocalMatrixType::template ContainerTypeByMDI<Mem::Main, DataType, IndexType> matrix_fine_local;
          matrix_fine_local.convert(matrix_fine.local());

          // list of all points, that point i strongly depends on, aka S_i
          std::map<Index, std::set<Index>> depends_on;
          // which points are strongly influenced by point i?, aka S_iT
          std::map<Index, std::set<Index>> influences;

          //calculate strong couplings
          _get_couplings(_interpolate_scalar_matrix(matrix_fine_local), depends_on, influences, theta);

          //lambda_i is the number of points strongly influenced by point i
          std::vector<Index> lambdas(matrix_fine_local.rows(), Index(0));
          for (Index i(0) ; i < matrix_fine_local.rows() ; ++i)
          {
            lambdas.at(i) = Index(influences[i].size());
          }

          //designated coarse and fine points and yet undecided
          std::set<Index> coarse;
          std::set<Index> fine;
          std::set<Index> undecided;
          for (Index i(0) ; i < matrix_fine_local.rows() ; ++i)
          {
            undecided.insert(i);
          }


          /// \todo use grid coloring for parallel processing
          Index rank = (Index)comm->rank();
          for (Index nextrank = 0 ; nextrank < (Index)comm->size() ; ++nextrank)
          {
            //comm->print(stringify(nextrank));
            if (nextrank == rank)
            {
              //first sweep

              //iterator to point with maximum remaining lambda value
              auto max_it = max_element(lambdas.begin(), lambdas.end());
              //loop until all lambdas are zero (all points have been designated)
              while (*max_it != Index(0))
              {
                //next max point selected, designate it as coarse point
                Index max_pos = (Index)std::distance(lambdas.begin(), max_it);
                coarse.insert(max_pos);
                undecided.erase(max_pos);
                *max_it = Index(0);

                //loop over all points that strongly depend on the currently choosen coarse point
                //these points should be fine points
                for (auto influencer : influences[max_pos])
                {
                  fine.insert(influencer);
                  undecided.erase(influencer);
                  lambdas.at(influencer) = Index(0);
                  //for each new fine point, increase the lambda of all remaining unassigned points that strongly influence this fine point
                  for (auto& depends : depends_on[influencer])
                  {
                    if (lambdas.at(depends) != Index(0))
                    {
                      lambdas.at(depends) += Index(1);
                    }
                  }
                }
                //decrease weihgts of points that influence max_pos
                for (auto depends : depends_on[max_pos])
                {
                  if (lambdas.at(depends) != Index(0))
                    lambdas.at(depends) -= Index(1);
                }
                //choose new unassigned point with maxium labmda
                max_it = max_element(lambdas.begin(), lambdas.end());
              }
              for (auto& spare : undecided)
              {
                fine.insert(spare);
              }
              undecided.clear();

              //second sweep, Algorithm 2.6 AmgPhaseII
              //loop over all fine points
              for (auto fine_it(fine.begin()) ; fine_it != fine.end() ; )
              {
                bool skip = false;
                Index i(*fine_it);
                std::set<Index> c_new;
                std::set<Index> Si_vs_F;

                //setup Sj_vs_F
                //loop over all strong influencing points
                for (auto& depend : depends_on[*fine_it])
                {
                  //check if influencing point is fine point
                  if (fine.count(depend) > 0)
                  {
                    Si_vs_F.insert(depend);
                  }
                }

                for (auto& j : Si_vs_F)
                {
                  //setup Sj_vs_Si_vs_C
                  std::set<Index> Sj_vs_Si_vs_C;
                  for (auto& k : depends_on[j])
                  {
                    if (depends_on[i].count(k) > 0 && coarse.count(k) > 0)
                    {
                      Sj_vs_Si_vs_C.insert(k);
                    }
                  }

                  if (Sj_vs_Si_vs_C.size() == 0)
                  {
                    if (c_new.size() != 0)
                    {
                      coarse.insert(i);
                      fine_it = fine.erase(fine_it);
                      skip = true;
                      break; //exit j loop
                    }
                    else
                    {
                      c_new.insert(j);
                    }
                  }
                }
                if(c_new.size() == 0)
                {
                  ++fine_it;
                  continue;
                }
                if (!skip)
                {
                  XASSERT(c_new.size() == 1);
                  coarse.insert(*c_new.begin());
                  if (i != *c_new.begin())
                  {
                    fine.erase(*c_new.begin());
                    ++fine_it;
                  }
                  else
                    fine_it = fine.erase(fine_it);
                }
              }

              //finished local c/f splitting
              //now set c/f points in neighbouring patches
              // 1 designates coarse and 0 designates fine
              /// \todo store buffers and use asynch mpi ops
              auto gate = matrix_fine.get_row_gate();
              if (gate != nullptr)
              {
                for (Index ni(0) ; ni < gate->_mirrors.size() ; ++ni)
                {
                  auto& mirror = gate->_mirrors.at(ni);
                  int* buffer = new int[mirror.num_indices()];
                  for (Index i(0) ; i < mirror.num_indices() ; ++i)
                  {
                    if (coarse.count(mirror.indices()[i]) > 0)
                      buffer[i] = 1;
                    else
                      buffer[i] = 0;
                  }
                  comm->send(buffer, mirror.num_indices(), gate->_ranks.at(ni));
                  delete[] buffer;
                }
              }
            }
            else
            {
              //add coarse/fine ; remove undecided ; set lambda to zero
              auto gate = matrix_fine.get_row_gate();
              for (Index ni(0) ; ni < gate->_mirrors.size() ; ++ni)
              {
                if (gate->_ranks.at(ni) == (int)nextrank)
                {
                  auto& mirror = gate->_mirrors.at(ni);
                  int* buffer = new int[mirror.num_indices()];
                  comm->recv(buffer, mirror.num_indices(), gate->_ranks.at(ni));
                  for (Index i(0) ; i < mirror.num_indices() ; ++i)
                  {
                    Index current_element(mirror.indices()[i]);
                    undecided.erase(current_element);
                    lambdas.at(current_element) = Index(0);
                    if (buffer[i] == 1)
                    {
                      fine.erase(current_element);
                      coarse.insert(current_element);
                      //loop over all points that strongly depend on the currently choosen coarse point
                      //these points should be fine points
                      for (auto influencer : influences[current_element])
                      {
                        fine.insert(influencer);
                        undecided.erase(influencer);
                        lambdas.at(influencer) = Index(0);
                        //for each new fine point, increase the lambda of all remaining unassigned points that strongly influence this fine point
                        for (auto& depends : depends_on[influencer])
                        {
                          if (lambdas.at(depends) != Index(0))
                          {
                            lambdas.at(depends) += Index(1);
                          }
                        }
                      }
                      //decrease weihgts of points that influence current_element
                      for (auto depends : depends_on[current_element])
                      {
                        if (lambdas.at(depends) != Index(0))
                          lambdas.at(depends) -= Index(1);
                      }
                    }
                    else
                    {
                      fine.insert(current_element);
                    }
                  }
                  delete[] buffer;
                }
              }
            }
          }

          // create interpolation matrix
          std::vector<IndexType> vrow_ptr;
          std::vector<IndexType> vcol_ind;
          //respect the following cases:
          // system bcsr && transfer bcsr -> no interpolate
          // system bcsr && transfer csr -> interpolate
          // system csr && transfer csr -> no interpolate
          std::vector<typename LocalTransferOperatorType::MatrixType::ValueType> vval; //exploit that BWrappedCSR has ValueType==DataType
          _create_interpolation_direct(
              _interpolate_scalar_matrix_if<!std::is_same<typename LocalTransferOperatorType::MatrixType::ValueType, typename LocalMatrixType::ValueType>::value>
              (matrix_fine_local), coarse, depends_on, vrow_ptr, vcol_ind, vval);

          typename LocalTransferOperatorType::MatrixType::template ContainerTypeByMDI<Mem::Main, DataType, IndexType> prolongation_main(matrix_fine_local.rows(), Index(coarse.size()), Index(vval.size()));
          for (Index i(0) ; i < vrow_ptr.size() ; ++i)
          {
            prolongation_main.row_ptr()[i] = vrow_ptr.at(i);
          }
          for (Index i(0) ; i < vval.size() ; ++i)
          {
            prolongation_main.val()[i] = vval.at(i);
            prolongation_main.col_ind()[i] = vcol_ind.at(i);
          }
          auto restriction_main = prolongation_main.transpose();

          // compute coarse matrix sparsity pattern
          Adjacency::Graph graph_tmp(Adjacency::RenderType::injectify, matrix_fine_local, prolongation_main);
          Adjacency::Graph graph_crs(Adjacency::RenderType::injectify, restriction_main, graph_tmp);
          typename LocalMatrixType::template ContainerTypeByMDI<Mem::Main, DataType, IndexType> matrix_coarse_main(graph_crs);

          // compute coarse matrix entries
          matrix_coarse_main.format(DataType(0));
          matrix_coarse_main.add_double_mat_mult(restriction_main, matrix_fine_local, prolongation_main);

          //create coarse unit filter
          std::vector<Index> new_filter_entries;
          for (Index i(0) ; i < filter_fine.local().used_elements() ; ++i)
          {
            auto it = coarse.find(filter_fine.local().get_indices()[i]);
            if (it != coarse.end())
            {
              new_filter_entries.push_back((Index)std::distance(coarse.begin(), it));
            }
          }

          GlobalFilterType gnf;
          if (new_filter_entries.size() > 0)
          {
            typename LocalFilterType::VectorType::template ContainerTypeByMDI<Mem::Main, DataType, Index> uf_values(Index(new_filter_entries.size()));
            uf_values.format(); //can be set to zero, as we are only working on defects and not any actual solution vector
            LAFEM::DenseVector<Mem::Main, Index, Index> uf_indices(Index(new_filter_entries.size()));
            for (Index i(0) ; i < uf_indices.size() ; ++i)
            {
              uf_indices(i, new_filter_entries.at(i));
            }
            //this call would not be valid with a zero sized (but allocated) uf_values vector
            LocalFilterType uf(Index(coarse.size()), uf_values, uf_indices);
            GlobalFilterType gnf2(uf.clone(LAFEM::CloneMode::Shallow));
            gnf.convert(gnf2);
          }
          else
          {
            LocalFilterType uf(Index(coarse.size()));
            GlobalFilterType gnf2(uf.clone(LAFEM::CloneMode::Shallow));
            gnf.convert(gnf2);
          }

          GlobalTransferOperatorType trans;
          typename LocalTransferOperatorType::MatrixType local_prolongation;
          local_prolongation.convert(prolongation_main);
          typename LocalTransferOperatorType::MatrixType local_restriction;
          local_restriction.convert(restriction_main);
          trans.local().get_mat_prol().convert(local_prolongation);
          trans.local().get_mat_rest().convert(local_restriction);

          //setup gate/mirrors
          //typename GlobalMatrixType::GateRowType * coarse_gate = nullptr;
          if (matrix_fine.get_row_gate() != nullptr)
          {
            //coarse_gate = new typename GlobalMatrixType::GateRowType(*comm);
            auto fine_gate = matrix_fine.get_row_gate();
            for (Index ni(0) ; ni < fine_gate->_mirrors.size() ; ++ni)
            {
              auto& fine_mirror = fine_gate->_mirrors.at(ni);
              //check if any coarse entry lies in this mirror
              bool use(false);
              for (Index i(0) ; i < fine_mirror.num_indices() ; ++i)
              {
                if (coarse.count(fine_mirror.indices()[i]) > 0)
                {
                  use = true;
                  break;
                }
              }
              if (!use)
                continue;

              //collect new mirror entries
              std::vector<Index> mirror_entries;
              for (Index i(0) ; i < fine_mirror.num_indices() ; ++i)
              {
                auto find_it = coarse.find(fine_mirror.indices()[i]);
                if (find_it != coarse.end())
                {
                  mirror_entries.push_back((Index)std::distance(coarse.begin(), find_it)); //index into coarse vector is index in coarse elements only
                }
              }

              //fill mirror indices
              gate_coarse->_mirrors.emplace_back(matrix_coarse_main.rows(), Index(mirror_entries.size()));
              auto& coarse_mirror = gate_coarse->_mirrors.back();
              for (Index i(0) ; i < coarse_mirror.num_indices() ; ++i)
              {
                coarse_mirror.indices()[i] = mirror_entries.at(i);
              }

              gate_coarse->_ranks.push_back(fine_gate->_ranks.at(ni));
            }
            typename LocalMatrixType::VectorTypeR freqs(matrix_coarse_main.rows());
            gate_coarse->compile(freqs.clone(LAFEM::CloneMode::Shallow));
          }
          LocalMatrixType local_matrix_coarse;
          local_matrix_coarse.convert(matrix_coarse_main);

          if (gate_out != nullptr)
          {
            gate_out->convert(*gate_coarse);
          }
            GlobalMatrixType matrix_coarse(gate_out, gate_out, local_matrix_coarse.clone(LAFEM::CloneMode::Shallow));
            _adaptive_clone(coarse_grid_out, matrix_coarse);
            _adaptive_clone(coarse_filter_out, gnf);
            _adaptive_clone(transfer_out, trans);
        }

      /**
       * \brief Update coarse grid level
       *
       * This function is used to update numerical values of a given coarse grid, whilst preserving
       * the old transfer operator and matrix layout.
       *
       * \note While this method saves some effort by not calculating a complete new coarsening scheme,
       * it may happen, that the old transfer operator and coarsening scheme is not optimal suited for the current
       * matrix entries.
       *
       * \param[in] fine_grid_input The fine grid matrix
       * \param[in] fine_transfer_input The transfer operator between the old fine and the old coarse grid level
       * \param[out] coarse_grid_out The updated coarse grid matrix
       */
        static void update_coarse_level(const SystemMatrix_ & fine_grid_input, TransferOperator_ & fine_transfer_input, SystemMatrix_ & coarse_grid_out)
        {
          const GlobalMatrixType matrix_fine(_make_matrix_global<SystemMatrix_>(fine_grid_input));
          typename LocalMatrixType::template ContainerTypeByMDI<Mem::Main, DataType, IndexType> matrix_fine_local;
          matrix_fine_local.convert(matrix_fine.local());
          GlobalMatrixType matrix_coarse(_make_matrix_global<SystemMatrix_>(coarse_grid_out));

          const GlobalTransferOperatorType transfer_fine(_make_transfer_global<TransferOperator_>(fine_transfer_input));
          typename LocalTransferOperatorType::template TransferTypeByMDI<Mem::Main, DataType, IndexType> transfer_fine_local;
          transfer_fine_local.convert(transfer_fine.local());

          // compute coarse matrix sparsity pattern
          Adjacency::Graph graph_tmp(Adjacency::RenderType::injectify, matrix_fine_local, transfer_fine_local.get_mat_prol());
          Adjacency::Graph graph_crs(Adjacency::RenderType::injectify, transfer_fine_local.get_mat_rest(), graph_tmp);
          typename LocalMatrixType::template ContainerTypeByMDI<Mem::Main, DataType, IndexType> matrix_coarse_main(graph_crs);

          // compute coarse matrix entries
          matrix_coarse_main.format(DataType(0));
          matrix_coarse_main.add_double_mat_mult(transfer_fine_local.get_mat_rest(), matrix_fine_local, transfer_fine_local.get_mat_prol());

          LocalMatrixType local_matrix_coarse;
          local_matrix_coarse.convert(matrix_coarse_main);

          matrix_coarse.local().clone(local_matrix_coarse, LAFEM::CloneMode::Shallow);
          _adaptive_clone(coarse_grid_out, matrix_coarse);
        }

      private:

        template <bool B_, typename MT_, typename = typename std::enable_if<B_>::type>
        static LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> _interpolate_scalar_matrix_if(const MT_ & matrix)
        {
          return _interpolate_scalar_matrix(matrix);
        }

        template <bool B_, typename MT_, typename = typename std::enable_if<!B_>::type>
        static MT_ _interpolate_scalar_matrix_if(const MT_ & matrix)
        {
          return matrix.clone(LAFEM::CloneMode::Shallow);
        }

        static LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> _interpolate_scalar_matrix(const LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> & matrix_fine_local)
        {
          return matrix_fine_local.clone(LAFEM::CloneMode::Shallow);
        }

        // (2.97) in diss metsch
        template <int BlockDim_>
        static LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> _interpolate_scalar_matrix(const LAFEM::SparseMatrixBCSR<Mem::Main, DataType, IndexType, BlockDim_, BlockDim_> & block_matrix)
        {
          LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> new_scalar_matrix(block_matrix.layout());

          for (Index row(0) ; row < new_scalar_matrix.rows() ; ++row)
          {
            for (Index col(new_scalar_matrix.row_ptr()[row]) ; col < new_scalar_matrix.row_ptr()[row+1] ; ++col)
            {
              Index j(new_scalar_matrix.col_ind()[col]);
              if (row != j)
              {
                new_scalar_matrix.val()[col] = -block_matrix.val()[col].norm_frobenius();
              }
              else
              {
                new_scalar_matrix.val()[col] = block_matrix.val()[col].norm_frobenius();
              }
            }
          }

          return new_scalar_matrix.clone(LAFEM::CloneMode::Shallow);
        }

        // (2.97) in diss metsch with alternative main diagonal entry
        /*template <int BlockDim_>
        static LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> _interpolate_scalar_matrix(const LAFEM::SparseMatrixBCSR<Mem::Main, DataType, IndexType, BlockDim_, BlockDim_> & block_matrix)
        {
          LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> new_scalar_matrix(block_matrix.layout());

          for (Index row(0) ; row < new_scalar_matrix.rows() ; ++row)
          {
            for (Index col(new_scalar_matrix.row_ptr()[row]) ; col < new_scalar_matrix.row_ptr()[row+1] ; ++col)
            {
              Index j(new_scalar_matrix.col_ind()[col]);
              if (row != j)
              {
                new_scalar_matrix.val()[col] = -block_matrix.val()[col].norm_frobenius();
              }
            }
          }

          for (Index row(0) ; row < new_scalar_matrix.rows() ; ++row)
          {
            for (Index col(new_scalar_matrix.row_ptr()[row]) ; col < new_scalar_matrix.row_ptr()[row+1] ; ++col)
            {
              Index j(new_scalar_matrix.col_ind()[col]);
              if (row == j)
              {
                DataType temp(0);
                // compute negative sum of all off-diagonal entries
                for (Index col2(new_scalar_matrix.row_ptr()[row]) ; col2 < new_scalar_matrix.row_ptr()[row+1] ; ++col2)
                {
                  Index j2(new_scalar_matrix.col_ind()[col2]);
                  if (row != j2)
                  {
                    temp += new_scalar_matrix.val()[col];
                  }
                  temp = -temp;
                }
                new_scalar_matrix.val()[col] = (temp != DataType(0)) ? temp : DataType(1);
                break;
              }
            }
          }

          return new_scalar_matrix.clone(LAFEM::CloneMode::Shallow);
        }*/

        // diss mesh page 77, strategy 2: only take value [0,0] from each block
        /*template <int BlockDim_>
        static LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> _interpolate_scalar_matrix(const LAFEM::SparseMatrixBCSR<Mem::Main, DataType, IndexType, BlockDim_, BlockDim_> & block_matrix)
        {
          LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> new_scalar_matrix(block_matrix.layout());

          for (Index row(0) ; row < new_scalar_matrix.rows() ; ++row)
          {
            for (Index col(new_scalar_matrix.row_ptr()[row]) ; col < new_scalar_matrix.row_ptr()[row+1] ; ++col)
            {
              new_scalar_matrix.val()[col] = block_matrix.val()[col][0][0];
            }
          }

          return new_scalar_matrix.clone(LAFEM::CloneMode::Shallow);
        }*/

        static void _get_couplings(const LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> & matrix_fine_local, std::map<Index, std::set<Index>> & depends_on,
            std::map<Index, std::set<Index>> & influences, DataType theta)
        {
          //set up dependency/influence maps
          for (Index row(0) ; row < matrix_fine_local.rows() ; ++row)
          {
            DataType max(0);
            for (Index col(matrix_fine_local.row_ptr()[row]) ; col < matrix_fine_local.row_ptr()[row+1] ; ++col)
            {
              if (matrix_fine_local.col_ind()[col] != row)
                max = Math::max(max, Math::abs(matrix_fine_local.val()[col]));
            }
            for (Index col(matrix_fine_local.row_ptr()[row]); col < matrix_fine_local.row_ptr()[row+1] ; ++col)
            {
              if (matrix_fine_local.col_ind()[col] != row)
              {
                if (Math::abs(matrix_fine_local.val()[col]) >= theta * max && Math::abs(matrix_fine_local.val()[col]) != DataType(0))
                {
                  depends_on[row].insert(matrix_fine_local.col_ind()[col]);
                  influences[matrix_fine_local.col_ind()[col]].insert(row);
                }
              }
            }
          }
        }

        //classical (diss metsch) interpolation - this is the interpolation from "A multigrid turoial" by Bricks et al
        static void _create_interpolation_classical(const LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> & matrix_fine_local, const std::set<Index> & coarse,
            const std::map<Index, std::set<Index>> & depends_on, std::vector<Index> & vrow_ptr, std::vector<Index> & vcol_ind, std::vector<DataType> & vval)
        {
          vrow_ptr.push_back(IndexType(0));
          for (Index row(0) ; row < matrix_fine_local.rows() ; ++row)
          {
            //check if row is present in coarse set
            auto find_it = coarse.find(row);
            if(find_it != coarse.end())
            {
              //row is in coarse, simply transfer error to fine grid
              vrow_ptr.push_back(vrow_ptr.back() + Index(1));
              vcol_ind.push_back((Index)std::distance(coarse.begin(), find_it)); //index into coarse vector is index in coarse elements only
              vval.push_back(DataType(1));
            }
            else
            {
              //row is not in coarse, interpolation via omega weights is necessary
              //setup neighbourhood sets
              std::set<Index> c_i;
              for (const auto& depends : depends_on.at(row))
              {
                if (coarse.count(depends) > 0)
                {
                  c_i.insert(depends);
                }
              }
              std::set<Index> d_is;
              for (const auto& depends : depends_on.at(row))
              {
                if (coarse.count(depends) == 0)
                {
                  d_is.insert(depends);
                }
              }
              std::set<Index> d_iw;
              for (Index col(matrix_fine_local.row_ptr()[row]); col < matrix_fine_local.row_ptr()[row+1] ; ++col)
              {
                Index col_idx(matrix_fine_local.col_ind()[col]);
                if (col_idx != row && matrix_fine_local.val()[col] != DataType(0) && c_i.count(col_idx) == 0 && d_is.count(col_idx) == 0)
                {
                  d_iw.insert(col_idx);
                }
              }

              //calculate omega for every member in c_i
              for (auto c_it(c_i.begin()) ; c_it != c_i.end() ; ++c_it)
              {
                Index i(row);
                Index j(*c_it);
                DataType omega(matrix_fine_local(i, j));

                DataType temp;
                for (auto& m : d_is)
                {
                  temp = matrix_fine_local(i, m) * matrix_fine_local(m, j);
                  DataType temp2(DataType(0));
                  for (auto& k : c_i)
                  {
                    temp2 += matrix_fine_local(m, k);
                  }
                  temp /= temp2;
                  omega += temp;
                }

                temp = matrix_fine_local(i, i);
                for (auto& n : d_iw)
                {
                  temp += matrix_fine_local(i, n);
                }

                omega /= temp;
                omega *= DataType(-1);

                vval.push_back(omega);
                vcol_ind.push_back((Index)std::distance(coarse.begin(), coarse.find(*c_it))); //index into coarse vector is index in coarse elements only
              }
              vrow_ptr.push_back(vrow_ptr.back() + c_i.size());
            }
          }
        }

        //direct (diss metsch) interpolation
        static void _create_interpolation_direct(const LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> & matrix_fine_local, const std::set<Index> & coarse,
            const std::map<Index, std::set<Index>> & depends_on, std::vector<Index> & vrow_ptr, std::vector<Index> & vcol_ind, std::vector<DataType> & vval)
        {
          vrow_ptr.push_back(IndexType(0));
          for (Index row(0) ; row < matrix_fine_local.rows() ; ++row)
          {
            //check if row is present in coarse set
            auto find_it = coarse.find(row);
            if(find_it != coarse.end())
            {
              //row is in coarse, simply transfer error to fine grid
              vrow_ptr.push_back(vrow_ptr.back() + Index(1));
              vcol_ind.push_back((Index)std::distance(coarse.begin(), find_it)); //index into coarse vector is index in coarse elements only
              vval.push_back(DataType(1));
            }
            else
            {
              //row is not in coarse, interpolation via omega weights is necessary
              //setup neighbourhood sets
              std::set<Index> c_i;
              for (const auto& depends : depends_on.at(row))
              {
                if (coarse.count(depends) > 0)
                {
                  c_i.insert(depends);
                }
              }
              std::set<Index> e_i;
              for (Index col(matrix_fine_local.row_ptr()[row]); col < matrix_fine_local.row_ptr()[row+1] ; ++col)
              {
                Index col_idx(matrix_fine_local.col_ind()[col]);
                if (col_idx != row && matrix_fine_local.val()[col] != DataType(0))
                {
                  e_i.insert(col_idx);
                }
              }

              //calculate omega for every member in c_i
              for (auto c_it(c_i.begin()) ; c_it != c_i.end() ; ++c_it)
              {
                Index i(row);
                Index j(*c_it);
                DataType omega(DataType(1) / matrix_fine_local(i, i));

                DataType temp(0);
                for (auto& k : e_i)
                {
                  temp += matrix_fine_local(i, k);
                }
                omega *= temp;

                temp = DataType(0);
                for (auto& k : c_i)
                {
                  temp += matrix_fine_local(i, k);
                }
                omega /= temp;

                omega *= matrix_fine_local(i, j);
                omega *= DataType(-1);

                vval.push_back(omega);
                vcol_ind.push_back((Index)std::distance(coarse.begin(), coarse.find(*c_it))); //index into coarse vector is index in coarse elements only
              }
              vrow_ptr.push_back(vrow_ptr.back() + Index(c_i.size()));
            }
          }
        }

        //block direct (diss metsch) interpolation (2.103)
        template <int BlockDim_>
        static void _create_interpolation_direct(const LAFEM::SparseMatrixBCSR<Mem::Main, DataType, IndexType, BlockDim_, BlockDim_> & matrix_fine_local, const std::set<Index> & coarse,
            const std::map<Index, std::set<Index>> & depends_on, std::vector<Index> & vrow_ptr, std::vector<Index> & vcol_ind,
            std::vector<Tiny::Matrix<DataType, BlockDim_, BlockDim_>> & vval)
        {
          using ValueType = Tiny::Matrix<DataType, BlockDim_, BlockDim_>;

          vrow_ptr.push_back(IndexType(0));
          for (Index row(0) ; row < matrix_fine_local.rows() ; ++row)
          {
            //check if row is present in coarse set
            auto find_it = coarse.find(row);
            if(find_it != coarse.end())
            {
              //row is in coarse, simply transfer error to fine grid
              vrow_ptr.push_back(vrow_ptr.back() + Index(1));
              vcol_ind.push_back((Index)std::distance(coarse.begin(), find_it)); //index into coarse vector is index in coarse elements only

              /*vval.push_back(ValueType(1));
              vval.back()[0][1] = 0;
              vval.back()[1][0] = 0;*/

              vval.emplace_back();
              for (Index col(matrix_fine_local.row_ptr()[row]) ; col < matrix_fine_local.row_ptr()[row+1] ; ++col)
              {
                if (matrix_fine_local.col_ind()[col] == row)
                {
                  for (int i(0) ; i < BlockDim_ ; ++i)
                  {
                    for (int j(0) ; j < BlockDim_ ; ++j)
                    {
                      vval.back()[i][j] = Math::abs(matrix_fine_local.val()[col][i][j]) < Math::eps<DataType>() ? DataType(0) : DataType(1);
                    }
                  }
                  break;
                }
              }

            }
            else
            {
              //row is not in coarse, interpolation via omega weights is necessary
              //setup neighbourhood sets
              std::set<Index> c_i;
              for (const auto& depends : depends_on.at(row))
              {
                if (coarse.count(depends) > 0)
                {
                  c_i.insert(depends);
                }
              }
              std::set<Index> e_i;
              for (Index col(matrix_fine_local.row_ptr()[row]) ; col < matrix_fine_local.row_ptr()[row+1] ; ++col)
              {
                Index col_idx(matrix_fine_local.col_ind()[col]);
                if (col_idx != row /*&& matrix_fine_local.val()[col].norm_frobenius() == DataType(0)*/)
                {
                  e_i.insert(col_idx);
                }
              }

              //calculate omega for every member in c_i
              for (auto c_it(c_i.begin()) ; c_it != c_i.end() ; ++c_it)
              {
                Index i(row);
                Index j(*c_it);
                ValueType omega;
                ValueType omega2;
                omega.set_inverse(matrix_fine_local(i, i));

                ValueType temp(0);
                for (auto& k : e_i)
                {
                  temp += matrix_fine_local(i, k);
                }
                omega2.set_mat_mat_mult(omega, temp);

                temp.format();
                for (auto& k : c_i)
                {
                  temp += matrix_fine_local(i, k);
                }
                ValueType temp2;
                temp2.set_inverse(temp);
                omega.set_mat_mat_mult(omega2, temp2);

                omega2.set_mat_mat_mult(omega, matrix_fine_local(i, j));
                omega2 *= DataType(-1);

                vval.push_back(omega2);
                vcol_ind.push_back((Index)std::distance(coarse.begin(), coarse.find(*c_it))); //index into coarse vector is index in coarse elements only
              }
              vrow_ptr.push_back(vrow_ptr.back() + c_i.size());
            }
          }
        }
    }; // class AMGFactory<...>
  }
}

#endif // KERNEL_SOLVER_AMG_HPP
