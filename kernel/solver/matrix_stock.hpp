// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_MATRIX_STOCK_HPP
#define KERNEL_SOLVER_MATRIX_STOCK_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/global/synch_mat.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/muxer.hpp>
#include <kernel/global/transfer.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/transfer.hpp>

namespace FEAT
{
  namespace Solver
  {
    /// Stock of matrices, based on global inpupt contains
    template <typename MatrixType_, typename FilterType_, typename TransferType_>
    class MatrixStock
    {
      public:
        using MatrixType = MatrixType_;
        using FilterType = FilterType_;
        using TransferType = TransferType_;
        using VectorType = typename MatrixType_::VectorTypeL;

        using GateRowType_ = typename MatrixType_::GateRowType;
        using GateColType_ = typename MatrixType_::GateColType;

        using MuxerType = typename TransferType_::MuxerType;

        /// virtual hierarchy sizes
        std::size_t size_virtual;

        // input deques with original containers
        std::deque<MatrixType_> systems;
        std::deque<GateRowType_*> gates_row;
        std::deque<GateColType_*> gates_col;
        std::deque<MuxerType*> muxers;
        std::deque<FilterType_> filters;
        std::deque<TransferType_> transfers;

        using MT_main_float_ulong = typename MatrixType_::template ContainerTypeByMDI<Mem::Main, float, unsigned long>;
        using MT_main_double_ulong = typename MatrixType_::template ContainerTypeByMDI<Mem::Main, double, unsigned long>;
#ifdef FEAT_HAVE_CUDA
        using MT_cuda_float_ulong = typename MatrixType_::template ContainerTypeByMDI<Mem::CUDA, float, unsigned long>;
        using MT_cuda_double_ulong = typename MatrixType_::template ContainerTypeByMDI<Mem::CUDA, double, unsigned long>;
#endif

        using FT_main_float_ulong = typename FilterType_::template FilterTypeByMDI<Mem::Main, float, unsigned long>;
        using FT_main_double_ulong = typename FilterType_::template FilterTypeByMDI<Mem::Main, double, unsigned long>;
#ifdef FEAT_HAVE_CUDA
        using FT_cuda_float_ulong = typename FilterType_::template FilterTypeByMDI<Mem::CUDA, float, unsigned long>;
        using FT_cuda_double_ulong = typename FilterType_::template FilterTypeByMDI<Mem::CUDA, double, unsigned long>;
#endif

        using TT_main_float_ulong = typename TransferType_::template TransferTypeByMDI<Mem::Main, float, unsigned long>;
        using TT_main_double_ulong = typename TransferType_::template TransferTypeByMDI<Mem::Main, double, unsigned long>;
#ifdef FEAT_HAVE_CUDA
        using TT_cuda_float_ulong = typename TransferType_::template TransferTypeByMDI<Mem::CUDA, float, unsigned long>;
        using TT_cuda_double_ulong = typename TransferType_::template TransferTypeByMDI<Mem::CUDA, double, unsigned long>;
#endif

        using GT_row_main_float_ulong = typename GateRowType_::template GateTypeByMDI<Mem::Main, float, unsigned long>;
        using GT_col_main_float_ulong = typename GateColType_::template GateTypeByMDI<Mem::Main, float, unsigned long>;
        using GT_row_main_double_ulong = typename GateRowType_::template GateTypeByMDI<Mem::Main, double, unsigned long>;
        using GT_col_main_double_ulong = typename GateColType_::template GateTypeByMDI<Mem::Main, double, unsigned long>;
#ifdef FEAT_HAVE_CUDA
        using GT_row_cuda_float_ulong = typename GateRowType_::template GateTypeByMDI<Mem::CUDA, float, unsigned long>;
        using GT_col_cuda_float_ulong = typename GateColType_::template GateTypeByMDI<Mem::CUDA, float, unsigned long>;
        using GT_row_cuda_double_ulong = typename GateRowType_::template GateTypeByMDI<Mem::CUDA, double, unsigned long>;
        using GT_col_cuda_double_ulong = typename GateColType_::template GateTypeByMDI<Mem::CUDA, double, unsigned long>;
#endif

        using MXT_main_float_ulong = typename MuxerType::template MuxerTypeByMDI<Mem::Main, float, unsigned long>;
        using MXT_main_double_ulong = typename MuxerType::template MuxerTypeByMDI<Mem::Main, double, unsigned long>;
#ifdef FEAT_HAVE_CUDA
        using MXT_cuda_float_ulong = typename MuxerType::template MuxerTypeByMDI<Mem::CUDA, float, unsigned long>;
        using MXT_cuda_double_ulong = typename MuxerType::template MuxerTypeByMDI<Mem::CUDA, double, unsigned long>;
#endif

        using MT_main_float_uint = typename MatrixType_::template ContainerTypeByMDI<Mem::Main, float, unsigned int>;
        using MT_main_double_uint = typename MatrixType_::template ContainerTypeByMDI<Mem::Main, double, unsigned int>;
#ifdef FEAT_HAVE_CUDA
        using MT_cuda_float_uint = typename MatrixType_::template ContainerTypeByMDI<Mem::CUDA, float, unsigned int>;
        using MT_cuda_double_uint = typename MatrixType_::template ContainerTypeByMDI<Mem::CUDA, double, unsigned int>;
#endif

        using FT_main_float_uint = typename FilterType_::template FilterTypeByMDI<Mem::Main, float, unsigned int>;
        using FT_main_double_uint = typename FilterType_::template FilterTypeByMDI<Mem::Main, double, unsigned int>;
#ifdef FEAT_HAVE_CUDA
        using FT_cuda_float_uint = typename FilterType_::template FilterTypeByMDI<Mem::CUDA, float, unsigned int>;
        using FT_cuda_double_uint = typename FilterType_::template FilterTypeByMDI<Mem::CUDA, double, unsigned int>;
#endif

        using TT_main_float_uint = typename TransferType_::template TransferTypeByMDI<Mem::Main, float, unsigned int>;
        using TT_main_double_uint = typename TransferType_::template TransferTypeByMDI<Mem::Main, double, unsigned int>;
#ifdef FEAT_HAVE_CUDA
        using TT_cuda_float_uint = typename TransferType_::template TransferTypeByMDI<Mem::CUDA, float, unsigned int>;
        using TT_cuda_double_uint = typename TransferType_::template TransferTypeByMDI<Mem::CUDA, double, unsigned int>;
#endif

        using GT_row_main_float_uint = typename GateRowType_::template GateTypeByMDI<Mem::Main, float, unsigned int>;
        using GT_col_main_float_uint = typename GateColType_::template GateTypeByMDI<Mem::Main, float, unsigned int>;
        using GT_row_main_double_uint = typename GateRowType_::template GateTypeByMDI<Mem::Main, double, unsigned int>;
        using GT_col_main_double_uint = typename GateColType_::template GateTypeByMDI<Mem::Main, double, unsigned int>;
#ifdef FEAT_HAVE_CUDA
        using GT_row_cuda_float_uint = typename GateRowType_::template GateTypeByMDI<Mem::CUDA, float, unsigned int>;
        using GT_col_cuda_float_uint = typename GateColType_::template GateTypeByMDI<Mem::CUDA, float, unsigned int>;
        using GT_row_cuda_double_uint = typename GateRowType_::template GateTypeByMDI<Mem::CUDA, double, unsigned int>;
        using GT_col_cuda_double_uint = typename GateColType_::template GateTypeByMDI<Mem::CUDA, double, unsigned int>;
#endif

        using MXT_main_float_uint = typename MuxerType::template MuxerTypeByMDI<Mem::Main, float, unsigned int>;
        using MXT_main_double_uint = typename MuxerType::template MuxerTypeByMDI<Mem::Main, double, unsigned int>;
#ifdef FEAT_HAVE_CUDA
        using MXT_cuda_float_uint = typename MuxerType::template MuxerTypeByMDI<Mem::CUDA, float, unsigned int>;
        using MXT_cuda_double_uint = typename MuxerType::template MuxerTypeByMDI<Mem::CUDA, double, unsigned int>;
#endif

        std::deque<MT_main_float_ulong> systems_main_float_ulong;
        std::deque<MT_main_double_ulong> systems_main_double_ulong;
#ifdef FEAT_HAVE_CUDA
        std::deque<MT_cuda_float_ulong> systems_cuda_float_ulong;
        std::deque<MT_cuda_double_ulong> systems_cuda_double_ulong;
#endif

        std::deque<FT_main_float_ulong> filters_main_float_ulong;
        std::deque<FT_main_double_ulong> filters_main_double_ulong;
#ifdef FEAT_HAVE_CUDA
        std::deque<FT_cuda_float_ulong> filters_cuda_float_ulong;
        std::deque<FT_cuda_double_ulong> filters_cuda_double_ulong;
#endif

        std::deque<TT_main_float_ulong> transfers_main_float_ulong;
        std::deque<TT_main_double_ulong> transfers_main_double_ulong;
#ifdef FEAT_HAVE_CUDA
        std::deque<TT_cuda_float_ulong> transfers_cuda_float_ulong;
        std::deque<TT_cuda_double_ulong> transfers_cuda_double_ulong;
#endif

        std::deque<std::shared_ptr<GT_row_main_float_ulong>> gates_row_main_float_ulong;
        std::deque<std::shared_ptr<GT_col_main_float_ulong>> gates_col_main_float_ulong;
        std::deque<std::shared_ptr<GT_row_main_double_ulong>> gates_row_main_double_ulong;
        std::deque<std::shared_ptr<GT_col_main_double_ulong>> gates_col_main_double_ulong;
#ifdef FEAT_HAVE_CUDA
        std::deque<std::shared_ptr<GT_row_cuda_float_ulong>> gates_row_cuda_float_ulong;
        std::deque<std::shared_ptr<GT_col_cuda_float_ulong>> gates_col_cuda_float_ulong;
        std::deque<std::shared_ptr<GT_row_cuda_double_ulong>> gates_row_cuda_double_ulong;
        std::deque<std::shared_ptr<GT_col_cuda_double_ulong>> gates_col_cuda_double_ulong;
#endif

        std::deque<std::shared_ptr<MXT_main_float_ulong>> muxers_main_float_ulong;
        std::deque<std::shared_ptr<MXT_main_double_ulong>> muxers_main_double_ulong;
#ifdef FEAT_HAVE_CUDA
        std::deque<std::shared_ptr<MXT_cuda_float_ulong>> muxers_cuda_float_ulong;
        std::deque<std::shared_ptr<MXT_cuda_double_ulong>> muxers_cuda_double_ulong;
#endif

        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_main_float_ulong, FT_main_float_ulong, TT_main_float_ulong> > > hierarchy_map_main_float_ulong;
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_main_double_ulong, FT_main_double_ulong, TT_main_double_ulong> > > hierarchy_map_main_double_ulong;
#ifdef FEAT_HAVE_CUDA
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_cuda_float_ulong, FT_cuda_float_ulong, TT_cuda_float_ulong> > > hierarchy_map_cuda_float_ulong;
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_cuda_double_ulong, FT_cuda_double_ulong, TT_cuda_double_ulong> > > hierarchy_map_cuda_double_ulong;
#endif

        std::deque<typename MT_main_float_ulong::LocalMatrixType> local_systems_main_float_ulong;
        std::deque<typename MT_main_double_ulong::LocalMatrixType> local_systems_main_double_ulong;
#ifdef FEAT_HAVE_CUDA
        std::deque<typename MT_cuda_float_ulong::LocalMatrixType> local_systems_cuda_float_ulong;
        std::deque<typename MT_cuda_double_ulong::LocalMatrixType> local_systems_cuda_double_ulong;
#endif

        std::deque<typename FT_main_float_ulong::LocalFilterType> local_filters_main_float_ulong;
        std::deque<typename FT_main_double_ulong::LocalFilterType> local_filters_main_double_ulong;
#ifdef FEAT_HAVE_CUDA
        std::deque<typename FT_cuda_float_ulong::LocalFilterType> local_filters_cuda_float_ulong;
        std::deque<typename FT_cuda_double_ulong::LocalFilterType> local_filters_cuda_double_ulong;
#endif

        std::deque<typename TT_main_float_ulong::LocalTransferType> local_transfers_main_float_ulong;
        std::deque<typename TT_main_double_ulong::LocalTransferType> local_transfers_main_double_ulong;
#ifdef FEAT_HAVE_CUDA
        std::deque<typename TT_cuda_float_ulong::LocalTransferType> local_transfers_cuda_float_ulong;
        std::deque<typename TT_cuda_double_ulong::LocalTransferType> local_transfers_cuda_double_ulong;
#endif

        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_main_float_ulong::LocalMatrixType, typename FT_main_float_ulong::LocalFilterType, typename TT_main_float_ulong::LocalTransferType>
          > > local_hierarchy_map_main_float_ulong;
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_main_double_ulong::LocalMatrixType, typename FT_main_double_ulong::LocalFilterType, typename TT_main_double_ulong::LocalTransferType>
          > > local_hierarchy_map_main_double_ulong;
#ifdef FEAT_HAVE_CUDA
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_cuda_float_ulong::LocalMatrixType, typename FT_cuda_float_ulong::LocalFilterType, typename TT_cuda_float_ulong::LocalTransferType>
          > > local_hierarchy_map_cuda_float_ulong;
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_cuda_double_ulong::LocalMatrixType, typename FT_cuda_double_ulong::LocalFilterType, typename TT_cuda_double_ulong::LocalTransferType>
          > > local_hierarchy_map_cuda_double_ulong;
#endif

        std::deque<MT_main_float_uint> systems_main_float_uint;
        std::deque<MT_main_double_uint> systems_main_double_uint;
#ifdef FEAT_HAVE_CUDA
        std::deque<MT_cuda_float_uint> systems_cuda_float_uint;
        std::deque<MT_cuda_double_uint> systems_cuda_double_uint;
#endif

        std::deque<FT_main_float_uint> filters_main_float_uint;
        std::deque<FT_main_double_uint> filters_main_double_uint;
#ifdef FEAT_HAVE_CUDA
        std::deque<FT_cuda_float_uint> filters_cuda_float_uint;
        std::deque<FT_cuda_double_uint> filters_cuda_double_uint;
#endif

        std::deque<TT_main_float_uint> transfers_main_float_uint;
        std::deque<TT_main_double_uint> transfers_main_double_uint;
#ifdef FEAT_HAVE_CUDA
        std::deque<TT_cuda_float_uint> transfers_cuda_float_uint;
        std::deque<TT_cuda_double_uint> transfers_cuda_double_uint;
#endif

        std::deque<std::shared_ptr<GT_row_main_float_uint>> gates_row_main_float_uint;
        std::deque<std::shared_ptr<GT_col_main_float_uint>> gates_col_main_float_uint;
        std::deque<std::shared_ptr<GT_row_main_double_uint>> gates_row_main_double_uint;
        std::deque<std::shared_ptr<GT_col_main_double_uint>> gates_col_main_double_uint;
#ifdef FEAT_HAVE_CUDA
        std::deque<std::shared_ptr<GT_row_cuda_float_uint>> gates_row_cuda_float_uint;
        std::deque<std::shared_ptr<GT_col_cuda_float_uint>> gates_col_cuda_float_uint;
        std::deque<std::shared_ptr<GT_row_cuda_double_uint>> gates_row_cuda_double_uint;
        std::deque<std::shared_ptr<GT_col_cuda_double_uint>> gates_col_cuda_double_uint;
#endif

        std::deque<std::shared_ptr<MXT_main_float_uint>> muxers_main_float_uint;
        std::deque<std::shared_ptr<MXT_main_double_uint>> muxers_main_double_uint;
#ifdef FEAT_HAVE_CUDA
        std::deque<std::shared_ptr<MXT_cuda_float_uint>> muxers_cuda_float_uint;
        std::deque<std::shared_ptr<MXT_cuda_double_uint>> muxers_cuda_double_uint;
#endif

        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_main_float_uint, FT_main_float_uint, TT_main_float_uint> > > hierarchy_map_main_float_uint;
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_main_double_uint, FT_main_double_uint, TT_main_double_uint> > > hierarchy_map_main_double_uint;
#ifdef FEAT_HAVE_CUDA
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_cuda_float_uint, FT_cuda_float_uint, TT_cuda_float_uint> > > hierarchy_map_cuda_float_uint;
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_cuda_double_uint, FT_cuda_double_uint, TT_cuda_double_uint> > > hierarchy_map_cuda_double_uint;
#endif

        std::deque<typename MT_main_float_uint::LocalMatrixType> local_systems_main_float_uint;
        std::deque<typename MT_main_double_uint::LocalMatrixType> local_systems_main_double_uint;
#ifdef FEAT_HAVE_CUDA
        std::deque<typename MT_cuda_float_uint::LocalMatrixType> local_systems_cuda_float_uint;
        std::deque<typename MT_cuda_double_uint::LocalMatrixType> local_systems_cuda_double_uint;
#endif

        std::deque<typename FT_main_float_uint::LocalFilterType> local_filters_main_float_uint;
        std::deque<typename FT_main_double_uint::LocalFilterType> local_filters_main_double_uint;
#ifdef FEAT_HAVE_CUDA
        std::deque<typename FT_cuda_float_uint::LocalFilterType> local_filters_cuda_float_uint;
        std::deque<typename FT_cuda_double_uint::LocalFilterType> local_filters_cuda_double_uint;
#endif

        std::deque<typename TT_main_float_uint::LocalTransferType> local_transfers_main_float_uint;
        std::deque<typename TT_main_double_uint::LocalTransferType> local_transfers_main_double_uint;
#ifdef FEAT_HAVE_CUDA
        std::deque<typename TT_cuda_float_uint::LocalTransferType> local_transfers_cuda_float_uint;
        std::deque<typename TT_cuda_double_uint::LocalTransferType> local_transfers_cuda_double_uint;
#endif

        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_main_float_uint::LocalMatrixType, typename FT_main_float_uint::LocalFilterType, typename TT_main_float_uint::LocalTransferType>
          > > local_hierarchy_map_main_float_uint;
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_main_double_uint::LocalMatrixType, typename FT_main_double_uint::LocalFilterType, typename TT_main_double_uint::LocalTransferType>
          > > local_hierarchy_map_main_double_uint;
#ifdef FEAT_HAVE_CUDA
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_cuda_float_uint::LocalMatrixType, typename FT_cuda_float_uint::LocalFilterType, typename TT_cuda_float_uint::LocalTransferType>
          > > local_hierarchy_map_cuda_float_uint;
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_cuda_double_uint::LocalMatrixType, typename FT_cuda_double_uint::LocalFilterType, typename TT_cuda_double_uint::LocalTransferType>
          > > local_hierarchy_map_cuda_double_uint;
#endif

        explicit MatrixStock(std::size_t size_virt) :
          size_virtual(size_virt)
        {
        }

        /// compile global systems
        template <typename SolverVectorType_>
        void compile_systems(typename SolverVectorType_::GateType *)
        {
          using Mem_ = typename SolverVectorType_::MemType;
          using DT_ = typename SolverVectorType_::DataType;
          using IT_ = typename SolverVectorType_::IndexType;

          if (typeid(Mem_) == typeid(Mem::Main) && typeid(DT_) == typeid(float) && typeid(IT_) == typeid(unsigned long))
          {
            if (! systems_main_float_ulong.empty())
              return;

            for (auto& gate_row : gates_row)
            {
              gates_row_main_float_ulong.push_back(std::make_shared<GT_row_main_float_ulong>());
              gates_row_main_float_ulong.back()->convert(*gate_row);
            }

            for (auto& gate_col : gates_col)
            {
              gates_col_main_float_ulong.push_back(std::make_shared<GT_col_main_float_ulong>());
              gates_col_main_float_ulong.back()->convert(*gate_col);
            }

            auto it_row = gates_row_main_float_ulong.begin();
            auto it_col = gates_col_main_float_ulong.begin();
            for (auto system_it = systems.begin() ; system_it != systems.end() ; ++system_it, ++it_row, ++it_col)
            {
              systems_main_float_ulong.emplace_back();
              systems_main_float_ulong.back().convert((*it_row).get(), (*it_col).get(), *system_it);
            }

            for (auto & filter : filters)
            {
              filters_main_float_ulong.emplace_back();
              filters_main_float_ulong.back().convert(filter);
            }

            for(auto & muxer : muxers)
            {
              muxers_main_float_ulong.emplace_back(std::make_shared<MXT_main_float_ulong>());
              muxers_main_float_ulong.back()->convert(*muxer);
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_main_float_ulong.emplace_back();
              transfers_main_float_ulong.back().convert(muxers_main_float_ulong.at(i).get(), transfers.at(i));
            }
          }

          else if (typeid(Mem_) == typeid(Mem::Main) && typeid(DT_) == typeid(double) && typeid(IT_) == typeid(unsigned long))
          {
            if (! systems_main_double_ulong.empty())
              return;

            for (auto& gate_row : gates_row)
            {
              gates_row_main_double_ulong.push_back(std::make_shared<GT_row_main_double_ulong>());
              gates_row_main_double_ulong.back()->convert(*gate_row);
            }

            for (auto& gate_col : gates_col)
            {
              gates_col_main_double_ulong.push_back(std::make_shared<GT_col_main_double_ulong>());
              gates_col_main_double_ulong.back()->convert(*gate_col);
            }

            auto it_row = gates_row_main_double_ulong.begin();
            auto it_col = gates_col_main_double_ulong.begin();
            for (auto system_it = systems.begin() ; system_it != systems.end() ; ++system_it, ++it_row, ++it_col)
            {
              systems_main_double_ulong.emplace_back();
              systems_main_double_ulong.back().convert((*it_row).get(), (*it_col).get(), *system_it);
            }

            for (auto & filter : filters)
            {
              filters_main_double_ulong.emplace_back();
              filters_main_double_ulong.back().convert(filter);
            }

            for(auto & muxer : muxers)
            {
              muxers_main_double_ulong.emplace_back(std::make_shared<MXT_main_double_ulong>());
              muxers_main_double_ulong.back()->convert(*muxer);
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_main_double_ulong.emplace_back();
              transfers_main_double_ulong.back().convert(muxers_main_double_ulong.at(i).get(), transfers.at(i));
            }
          }

#ifdef FEAT_HAVE_CUDA
          else if (typeid(Mem_) == typeid(Mem::CUDA) && typeid(DT_) == typeid(float) && typeid(IT_) == typeid(unsigned long))
          {
            if (! systems_cuda_float_ulong.empty())
              return;

            for (auto& gate_row : gates_row)
            {
              gates_row_cuda_float_ulong.push_back(std::make_shared<GT_row_cuda_float_ulong>());
              gates_row_cuda_float_ulong.back()->convert(*gate_row);
            }

            for (auto& gate_col : gates_col)
            {
              gates_col_cuda_float_ulong.push_back(std::make_shared<GT_col_cuda_float_ulong>());
              gates_col_cuda_float_ulong.back()->convert(*gate_col);
            }

            auto it_row = gates_row_cuda_float_ulong.begin();
            auto it_col = gates_col_cuda_float_ulong.begin();
            for (auto system_it = systems.begin() ; system_it != systems.end() ; ++system_it, ++it_row, ++it_col)
            {
              systems_cuda_float_ulong.emplace_back();
              systems_cuda_float_ulong.back().convert((*it_row).get(), (*it_col).get(), *system_it);
            }

            for (auto & filter : filters)
            {
              filters_cuda_float_ulong.emplace_back();
              filters_cuda_float_ulong.back().convert(filter);
            }

            for(auto & muxer : muxers)
            {
              muxers_cuda_float_ulong.emplace_back(std::make_shared<MXT_cuda_float_ulong>());
              muxers_cuda_float_ulong.back()->convert(*muxer);
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_cuda_float_ulong.emplace_back();
              transfers_cuda_float_ulong.back().convert(muxers_cuda_float_ulong.at(i).get(), transfers.at(i));
            }
          }

          else if (typeid(Mem_) == typeid(Mem::CUDA) && typeid(DT_) == typeid(double) && typeid(IT_) == typeid(unsigned long))
          {
            if (! systems_cuda_double_ulong.empty())
              return;

            for (auto& gate_row : gates_row)
            {
              gates_row_cuda_double_ulong.push_back(std::make_shared<GT_row_cuda_double_ulong>());
              gates_row_cuda_double_ulong.back()->convert(*gate_row);
            }

            for (auto& gate_col : gates_col)
            {
              gates_col_cuda_double_ulong.push_back(std::make_shared<GT_col_cuda_double_ulong>());
              gates_col_cuda_double_ulong.back()->convert(*gate_col);
            }

            auto it_row = gates_row_cuda_double_ulong.begin();
            auto it_col = gates_col_cuda_double_ulong.begin();
            for (auto system_it = systems.begin() ; system_it != systems.end() ; ++system_it, ++it_row, ++it_col)
            {
              systems_cuda_double_ulong.emplace_back();
              systems_cuda_double_ulong.back().convert((*it_row).get(), (*it_col).get(), *system_it);
            }

            for (auto & filter : filters)
            {
              filters_cuda_double_ulong.emplace_back();
              filters_cuda_double_ulong.back().convert(filter);
            }

            for(auto & muxer : muxers)
            {
              muxers_cuda_double_ulong.emplace_back(std::make_shared<MXT_cuda_double_ulong>());
              muxers_cuda_double_ulong.back()->convert(*muxer);
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_cuda_double_ulong.emplace_back();
              transfers_cuda_double_ulong.back().convert(muxers_cuda_double_ulong.at(i).get(), transfers.at(i));
            }
          }
#endif

          else if (typeid(Mem_) == typeid(Mem::Main) && typeid(DT_) == typeid(float) && typeid(IT_) == typeid(unsigned int))
          {
            if (! systems_main_float_uint.empty())
              return;

            for (auto& gate_row : gates_row)
            {
              gates_row_main_float_uint.push_back(std::make_shared<GT_row_main_float_uint>());
              gates_row_main_float_uint.back()->convert(*gate_row);
            }

            for (auto& gate_col : gates_col)
            {
              gates_col_main_float_uint.push_back(std::make_shared<GT_col_main_float_uint>());
              gates_col_main_float_uint.back()->convert(*gate_col);
            }

            auto it_row = gates_row_main_float_uint.begin();
            auto it_col = gates_col_main_float_uint.begin();
            for (auto system_it = systems.begin() ; system_it != systems.end() ; ++system_it, ++it_row, ++it_col)
            {
              systems_main_float_uint.emplace_back();
              systems_main_float_uint.back().convert((*it_row).get(), (*it_col).get(), *system_it);
            }

            for (auto & filter : filters)
            {
              filters_main_float_uint.emplace_back();
              filters_main_float_uint.back().convert(filter);
            }

            for(auto & muxer : muxers)
            {
              muxers_main_float_uint.emplace_back(std::make_shared<MXT_main_float_uint>());
              muxers_main_float_uint.back()->convert(*muxer);
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_main_float_uint.emplace_back();
              transfers_main_float_uint.back().convert(muxers_main_float_uint.at(i).get(), transfers.at(i));
            }
          }

          else if (typeid(Mem_) == typeid(Mem::Main) && typeid(DT_) == typeid(double) && typeid(IT_) == typeid(unsigned int))
          {
            if (! systems_main_double_uint.empty())
              return;

            for (auto& gate_row : gates_row)
            {
              gates_row_main_double_uint.push_back(std::make_shared<GT_row_main_double_uint>());
              gates_row_main_double_uint.back()->convert(*gate_row);
            }

            for (auto& gate_col : gates_col)
            {
              gates_col_main_double_uint.push_back(std::make_shared<GT_col_main_double_uint>());
              gates_col_main_double_uint.back()->convert(*gate_col);
            }

            auto it_row = gates_row_main_double_uint.begin();
            auto it_col = gates_col_main_double_uint.begin();
            for (auto system_it = systems.begin() ; system_it != systems.end() ; ++system_it, ++it_row, ++it_col)
            {
              systems_main_double_uint.emplace_back();
              systems_main_double_uint.back().convert((*it_row).get(), (*it_col).get(), *system_it);
            }

            for (auto & filter : filters)
            {
              filters_main_double_uint.emplace_back();
              filters_main_double_uint.back().convert(filter);
            }

            for(auto & muxer : muxers)
            {
              muxers_main_double_uint.emplace_back(std::make_shared<MXT_main_double_uint>());
              muxers_main_double_uint.back()->convert(*muxer);
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_main_double_uint.emplace_back();
              transfers_main_double_uint.back().convert(muxers_main_double_uint.at(i).get(), transfers.at(i));
            }
          }

#ifdef FEAT_HAVE_CUDA
          else if (typeid(Mem_) == typeid(Mem::CUDA) && typeid(DT_) == typeid(float) && typeid(IT_) == typeid(unsigned int))
          {
            if (! systems_cuda_float_uint.empty())
              return;

            for (auto& gate_row : gates_row)
            {
              gates_row_cuda_float_uint.push_back(std::make_shared<GT_row_cuda_float_uint>());
              gates_row_cuda_float_uint.back()->convert(*gate_row);
            }

            for (auto& gate_col : gates_col)
            {
              gates_col_cuda_float_uint.push_back(std::make_shared<GT_col_cuda_float_uint>());
              gates_col_cuda_float_uint.back()->convert(*gate_col);
            }

            auto it_row = gates_row_cuda_float_uint.begin();
            auto it_col = gates_col_cuda_float_uint.begin();
            for (auto system_it = systems.begin() ; system_it != systems.end() ; ++system_it, ++it_row, ++it_col)
            {
              systems_cuda_float_uint.emplace_back();
              systems_cuda_float_uint.back().convert((*it_row).get(), (*it_col).get(), *system_it);
            }

            for (auto & filter : filters)
            {
              filters_cuda_float_uint.emplace_back();
              filters_cuda_float_uint.back().convert(filter);
            }

            for(auto & muxer : muxers)
            {
              muxers_cuda_float_uint.emplace_back(std::make_shared<MXT_cuda_float_uint>());
              muxers_cuda_float_uint.back()->convert(*muxer);
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_cuda_float_uint.emplace_back();
              transfers_cuda_float_uint.back().convert(muxers_cuda_float_uint.at(i).get(), transfers.at(i));
            }
          }

          else if (typeid(Mem_) == typeid(Mem::CUDA) && typeid(DT_) == typeid(double) && typeid(IT_) == typeid(unsigned int))
          {
            if (! systems_cuda_double_uint.empty())
              return;

            for (auto& gate_row : gates_row)
            {
              gates_row_cuda_double_uint.push_back(std::make_shared<GT_row_cuda_double_uint>());
              gates_row_cuda_double_uint.back()->convert(*gate_row);
            }

            for (auto& gate_col : gates_col)
            {
              gates_col_cuda_double_uint.push_back(std::make_shared<GT_col_cuda_double_uint>());
              gates_col_cuda_double_uint.back()->convert(*gate_col);
            }

            auto it_row = gates_row_cuda_double_uint.begin();
            auto it_col = gates_col_cuda_double_uint.begin();
            for (auto system_it = systems.begin() ; system_it != systems.end() ; ++system_it, ++it_row, ++it_col)
            {
              systems_cuda_double_uint.emplace_back();
              systems_cuda_double_uint.back().convert((*it_row).get(), (*it_col).get(), *system_it);
            }

            for (auto & filter : filters)
            {
              filters_cuda_double_uint.emplace_back();
              filters_cuda_double_uint.back().convert(filter);
            }

            for(auto & muxer : muxers)
            {
              muxers_cuda_double_uint.emplace_back(std::make_shared<MXT_cuda_double_uint>());
              muxers_cuda_double_uint.back()->convert(*muxer);
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_cuda_double_uint.emplace_back();
              transfers_cuda_double_uint.back().convert(muxers_cuda_double_uint.at(i).get(), transfers.at(i));
            }
          }
#endif
        }

        /// compile local systems
        template <typename SolverVectorType_>
        void compile_systems(...)
        {
          using Mem_ = typename SolverVectorType_::MemType;
          using DT_ = typename SolverVectorType_::DataType;
          using IT_ = typename SolverVectorType_::IndexType;

          if (typeid(Mem_) == typeid(Mem::Main) && typeid(DT_) == typeid(float) && typeid(IT_) == typeid(unsigned long))
          {
            compile_systems<typename decltype(systems_main_float_ulong)::value_type::VectorTypeL>(nullptr);

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_main_float_ulong.emplace_back();
              local_systems_main_float_ulong.back().convert(local_matrix);
            }

            for (auto& filter : filters_main_float_ulong)
            {
              local_filters_main_float_ulong.push_back(filter.local().clone(LAFEM::CloneMode::Shallow));
            }

            for (auto& transfer : transfers_main_float_ulong)
            {
              local_transfers_main_float_ulong.push_back(transfer.local().clone(LAFEM::CloneMode::Shallow));
            }

          }

          else if (typeid(Mem_) == typeid(Mem::Main) && typeid(DT_) == typeid(double) && typeid(IT_) == typeid(unsigned long))
          {
            compile_systems<typename decltype(systems_main_double_ulong)::value_type::VectorTypeL>(nullptr);

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_main_double_ulong.emplace_back();
              local_systems_main_double_ulong.back().convert(local_matrix);
            }

            for (auto& filter : filters_main_double_ulong)
            {
              local_filters_main_double_ulong.push_back(filter.local().clone(LAFEM::CloneMode::Shallow));
            }

            for (auto& transfer : transfers_main_double_ulong)
            {
              local_transfers_main_double_ulong.push_back(transfer.local().clone(LAFEM::CloneMode::Shallow));
            }

          }

#ifdef FEAT_HAVE_CUDA
          else if (typeid(Mem_) == typeid(Mem::CUDA) && typeid(DT_) == typeid(float) && typeid(IT_) == typeid(unsigned long))
          {
            compile_systems<typename decltype(systems_cuda_float_ulong)::value_type::VectorTypeL>(nullptr);

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_cuda_float_ulong.emplace_back();
              local_systems_cuda_float_ulong.back().convert(local_matrix);
            }

            for (auto& filter : filters_cuda_float_ulong)
            {
              local_filters_cuda_float_ulong.push_back(filter.local().clone(LAFEM::CloneMode::Shallow));
            }

            for (auto& transfer : transfers_cuda_float_ulong)
            {
              local_transfers_cuda_float_ulong.push_back(transfer.local().clone(LAFEM::CloneMode::Shallow));
            }

          }

          else if (typeid(Mem_) == typeid(Mem::CUDA) && typeid(DT_) == typeid(double) && typeid(IT_) == typeid(unsigned long))
          {
            compile_systems<typename decltype(systems_cuda_double_ulong)::value_type::VectorTypeL>(nullptr);

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_cuda_double_ulong.emplace_back();
              local_systems_cuda_double_ulong.back().convert(local_matrix);
            }

            for (auto& filter : filters_cuda_double_ulong)
            {
              local_filters_cuda_double_ulong.push_back(filter.local().clone(LAFEM::CloneMode::Shallow));
            }

            for (auto& transfer : transfers_cuda_double_ulong)
            {
              local_transfers_cuda_double_ulong.push_back(transfer.local().clone(LAFEM::CloneMode::Shallow));
            }

          }
#endif

          else if (typeid(Mem_) == typeid(Mem::Main) && typeid(DT_) == typeid(float) && typeid(IT_) == typeid(unsigned int))
          {
            compile_systems<typename decltype(systems_main_float_uint)::value_type::VectorTypeL>(nullptr);

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_main_float_uint.emplace_back();
              local_systems_main_float_uint.back().convert(local_matrix);
            }

            for (auto& filter : filters_main_float_uint)
            {
              local_filters_main_float_uint.push_back(filter.local().clone(LAFEM::CloneMode::Shallow));
            }

            for (auto& transfer : transfers_main_float_uint)
            {
              local_transfers_main_float_uint.push_back(transfer.local().clone(LAFEM::CloneMode::Shallow));
            }

          }

          else if (typeid(Mem_) == typeid(Mem::Main) && typeid(DT_) == typeid(double) && typeid(IT_) == typeid(unsigned int))
          {
            compile_systems<typename decltype(systems_main_double_uint)::value_type::VectorTypeL>(nullptr);

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_main_double_uint.emplace_back();
              local_systems_main_double_uint.back().convert(local_matrix);
            }

            for (auto& filter : filters_main_double_uint)
            {
              local_filters_main_double_uint.push_back(filter.local().clone(LAFEM::CloneMode::Shallow));
            }

            for (auto& transfer : transfers_main_double_uint)
            {
              local_transfers_main_double_uint.push_back(transfer.local().clone(LAFEM::CloneMode::Shallow));
            }

          }

#ifdef FEAT_HAVE_CUDA
          else if (typeid(Mem_) == typeid(Mem::CUDA) && typeid(DT_) == typeid(float) && typeid(IT_) == typeid(unsigned int))
          {
            compile_systems<typename decltype(systems_cuda_float_uint)::value_type::VectorTypeL>(nullptr);

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_cuda_float_uint.emplace_back();
              local_systems_cuda_float_uint.back().convert(local_matrix);
            }

            for (auto& filter : filters_cuda_float_uint)
            {
              local_filters_cuda_float_uint.push_back(filter.local().clone(LAFEM::CloneMode::Shallow));
            }

            for (auto& transfer : transfers_cuda_float_uint)
            {
              local_transfers_cuda_float_uint.push_back(transfer.local().clone(LAFEM::CloneMode::Shallow));
            }

          }

          else if (typeid(Mem_) == typeid(Mem::CUDA) && typeid(DT_) == typeid(double) && typeid(IT_) == typeid(unsigned int))
          {
            compile_systems<typename decltype(systems_cuda_double_uint)::value_type::VectorTypeL>(nullptr);

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_cuda_double_uint.emplace_back();
              local_systems_cuda_double_uint.back().convert(local_matrix);
            }

            for (auto& filter : filters_cuda_double_uint)
            {
              local_filters_cuda_double_uint.push_back(filter.local().clone(LAFEM::CloneMode::Shallow));
            }

            for (auto& transfer : transfers_cuda_double_uint)
            {
              local_transfers_cuda_double_uint.push_back(transfer.local().clone(LAFEM::CloneMode::Shallow));
            }

          }
#endif
        }

        /// rebuilt all data in use from the original sources (that may have been changed)
        void refresh()
        {
          if (! systems_main_float_ulong.empty())
          {
            for (Index i(0) ; i < gates_row.size() ; ++i)
            {
              gates_row_main_float_ulong.at(i)->convert(*(gates_row.at(i)));
            }

            for (Index i(0) ; i < gates_col.size() ; ++i)
            {
              gates_col_main_float_ulong.at(i)->convert(*(gates_col.at(i)));
            }

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              systems_main_float_ulong.at(i).convert((gates_row_main_float_ulong.at(i)).get(), (gates_col_main_float_ulong.at(i)).get(), systems.at(i));
            }

            for (Index i(0) ; i < filters.size() ; ++i)
            {
              filters_main_float_ulong.at(i).convert(filters.at(i));
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_main_float_ulong.at(i).convert(muxers_main_float_ulong.at(i).get(), transfers.at(i));
            }
          }

          if (! systems_main_double_ulong.empty())
          {
            for (Index i(0) ; i < gates_row.size() ; ++i)
            {
              gates_row_main_double_ulong.at(i)->convert(*(gates_row.at(i)));
            }

            for (Index i(0) ; i < gates_col.size() ; ++i)
            {
              gates_col_main_double_ulong.at(i)->convert(*(gates_col.at(i)));
            }

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              systems_main_double_ulong.at(i).convert((gates_row_main_double_ulong.at(i)).get(), (gates_col_main_double_ulong.at(i)).get(), systems.at(i));
            }

            for (Index i(0) ; i < filters.size() ; ++i)
            {
              filters_main_double_ulong.at(i).convert(filters.at(i));
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_main_double_ulong.at(i).convert(muxers_main_double_ulong.at(i).get(), transfers.at(i));
            }
          }

#ifdef FEAT_HAVE_CUDA
          if (! systems_cuda_float_ulong.empty())
          {
            for (Index i(0) ; i < gates_row.size() ; ++i)
            {
              gates_row_cuda_float_ulong.at(i)->convert(*(gates_row.at(i)));
            }

            for (Index i(0) ; i < gates_col.size() ; ++i)
            {
              gates_col_cuda_float_ulong.at(i)->convert(*(gates_col.at(i)));
            }

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              systems_cuda_float_ulong.at(i).convert((gates_row_cuda_float_ulong.at(i)).get(), (gates_col_cuda_float_ulong.at(i)).get(), systems.at(i));
            }

            for (Index i(0) ; i < filters.size() ; ++i)
            {
              filters_cuda_float_ulong.at(i).convert(filters.at(i));
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_cuda_float_ulong.at(i).convert(muxers_cuda_float_ulong.at(i).get(), transfers.at(i));
            }
          }

          if (! systems_cuda_double_ulong.empty())
          {
            for (Index i(0) ; i < gates_row.size() ; ++i)
            {
              gates_row_cuda_double_ulong.at(i)->convert(*(gates_row.at(i)));
            }

            for (Index i(0) ; i < gates_col.size() ; ++i)
            {
              gates_col_cuda_double_ulong.at(i)->convert(*(gates_col.at(i)));
            }

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              systems_cuda_double_ulong.at(i).convert((gates_row_cuda_double_ulong.at(i)).get(), (gates_col_cuda_double_ulong.at(i)).get(), systems.at(i));
            }

            for (Index i(0) ; i < filters.size() ; ++i)
            {
              filters_cuda_double_ulong.at(i).convert(filters.at(i));
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_cuda_double_ulong.at(i).convert(muxers_cuda_double_ulong.at(i).get(), transfers.at(i));
            }
          }
#endif

          if (! systems_main_float_uint.empty())
          {
            for (Index i(0) ; i < gates_row.size() ; ++i)
            {
              gates_row_main_float_uint.at(i)->convert(*(gates_row.at(i)));
            }

            for (Index i(0) ; i < gates_col.size() ; ++i)
            {
              gates_col_main_float_uint.at(i)->convert(*(gates_col.at(i)));
            }

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              systems_main_float_uint.at(i).convert((gates_row_main_float_uint.at(i)).get(), (gates_col_main_float_uint.at(i)).get(), systems.at(i));
            }

            for (Index i(0) ; i < filters.size() ; ++i)
            {
              filters_main_float_uint.at(i).convert(filters.at(i));
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_main_float_uint.at(i).convert(muxers_main_float_uint.at(i).get(), transfers.at(i));
            }
          }

          if (! systems_main_double_uint.empty())
          {
            for (Index i(0) ; i < gates_row.size() ; ++i)
            {
              gates_row_main_double_uint.at(i)->convert(*(gates_row.at(i)));
            }

            for (Index i(0) ; i < gates_col.size() ; ++i)
            {
              gates_col_main_double_uint.at(i)->convert(*(gates_col.at(i)));
            }

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              systems_main_double_uint.at(i).convert((gates_row_main_double_uint.at(i)).get(), (gates_col_main_double_uint.at(i)).get(), systems.at(i));
            }

            for (Index i(0) ; i < filters.size() ; ++i)
            {
              filters_main_double_uint.at(i).convert(filters.at(i));
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_main_double_uint.at(i).convert(muxers_main_double_uint.at(i).get(), transfers.at(i));
            }
          }

#ifdef FEAT_HAVE_CUDA
          if (! systems_cuda_float_uint.empty())
          {
            for (Index i(0) ; i < gates_row.size() ; ++i)
            {
              gates_row_cuda_float_uint.at(i)->convert(*(gates_row.at(i)));
            }

            for (Index i(0) ; i < gates_col.size() ; ++i)
            {
              gates_col_cuda_float_uint.at(i)->convert(*(gates_col.at(i)));
            }

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              systems_cuda_float_uint.at(i).convert((gates_row_cuda_float_uint.at(i)).get(), (gates_col_cuda_float_uint.at(i)).get(), systems.at(i));
            }

            for (Index i(0) ; i < filters.size() ; ++i)
            {
              filters_cuda_float_uint.at(i).convert(filters.at(i));
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_cuda_float_uint.at(i).convert(muxers_cuda_float_uint.at(i).get(), transfers.at(i));
            }
          }

          if (! systems_cuda_double_uint.empty())
          {
            for (Index i(0) ; i < gates_row.size() ; ++i)
            {
              gates_row_cuda_double_uint.at(i)->convert(*(gates_row.at(i)));
            }

            for (Index i(0) ; i < gates_col.size() ; ++i)
            {
              gates_col_cuda_double_uint.at(i)->convert(*(gates_col.at(i)));
            }

            for (Index i(0) ; i < systems.size() ; ++i)
            {
              systems_cuda_double_uint.at(i).convert((gates_row_cuda_double_uint.at(i)).get(), (gates_col_cuda_double_uint.at(i)).get(), systems.at(i));
            }

            for (Index i(0) ; i < filters.size() ; ++i)
            {
              filters_cuda_double_uint.at(i).convert(filters.at(i));
            }

            for (Index i(0) ; i < transfers.size() ; ++i)
            {
              transfers_cuda_double_uint.at(i).convert(muxers_cuda_double_uint.at(i).get(), transfers.at(i));
            }
          }
#endif

          if (! local_systems_main_float_ulong.empty())
          {
            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_main_float_ulong.at(i).convert(local_matrix);
            }

            for (Index i(0) ; i < filters_main_float_ulong.size() ; ++i)
            {
              local_filters_main_float_ulong.at(i).convert(filters_main_float_ulong.at(i).local());
            }

            for (Index i(0) ; i < transfers_main_float_ulong.size() ; ++i)
            {
              local_transfers_main_float_ulong.at(i).convert(transfers_main_float_ulong.at(i).local());
            }

          }

          if (! local_systems_main_double_ulong.empty())
          {
            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_main_double_ulong.at(i).convert(local_matrix);
            }

            for (Index i(0) ; i < filters_main_double_ulong.size() ; ++i)
            {
              local_filters_main_double_ulong.at(i).convert(filters_main_double_ulong.at(i).local());
            }

            for (Index i(0) ; i < transfers_main_double_ulong.size() ; ++i)
            {
              local_transfers_main_double_ulong.at(i).convert(transfers_main_double_ulong.at(i).local());
            }

          }

#ifdef FEAT_HAVE_CUDA
          if (! local_systems_cuda_float_ulong.empty())
          {
            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_cuda_float_ulong.at(i).convert(local_matrix);
            }

            for (Index i(0) ; i < filters_cuda_float_ulong.size() ; ++i)
            {
              local_filters_cuda_float_ulong.at(i).convert(filters_cuda_float_ulong.at(i).local());
            }

            for (Index i(0) ; i < transfers_cuda_float_ulong.size() ; ++i)
            {
              local_transfers_cuda_float_ulong.at(i).convert(transfers_cuda_float_ulong.at(i).local());
            }

          }

          if (! local_systems_cuda_double_ulong.empty())
          {
            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_cuda_double_ulong.at(i).convert(local_matrix);
            }

            for (Index i(0) ; i < filters_cuda_double_ulong.size() ; ++i)
            {
              local_filters_cuda_double_ulong.at(i).convert(filters_cuda_double_ulong.at(i).local());
            }

            for (Index i(0) ; i < transfers_cuda_double_ulong.size() ; ++i)
            {
              local_transfers_cuda_double_ulong.at(i).convert(transfers_cuda_double_ulong.at(i).local());
            }

          }
#endif

          if (! local_systems_main_float_uint.empty())
          {
            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_main_float_uint.at(i).convert(local_matrix);
            }

            for (Index i(0) ; i < filters_main_float_uint.size() ; ++i)
            {
              local_filters_main_float_uint.at(i).convert(filters_main_float_uint.at(i).local());
            }

            for (Index i(0) ; i < transfers_main_float_uint.size() ; ++i)
            {
              local_transfers_main_float_uint.at(i).convert(transfers_main_float_uint.at(i).local());
            }

          }

          if (! local_systems_main_double_uint.empty())
          {
            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_main_double_uint.at(i).convert(local_matrix);
            }

            for (Index i(0) ; i < filters_main_double_uint.size() ; ++i)
            {
              local_filters_main_double_uint.at(i).convert(filters_main_double_uint.at(i).local());
            }

            for (Index i(0) ; i < transfers_main_double_uint.size() ; ++i)
            {
              local_transfers_main_double_uint.at(i).convert(transfers_main_double_uint.at(i).local());
            }

          }

#ifdef FEAT_HAVE_CUDA
          if (! local_systems_cuda_float_uint.empty())
          {
            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_cuda_float_uint.at(i).convert(local_matrix);
            }

            for (Index i(0) ; i < filters_cuda_float_uint.size() ; ++i)
            {
              local_filters_cuda_float_uint.at(i).convert(filters_cuda_float_uint.at(i).local());
            }

            for (Index i(0) ; i < transfers_cuda_float_uint.size() ; ++i)
            {
              local_transfers_cuda_float_uint.at(i).convert(transfers_cuda_float_uint.at(i).local());
            }

          }

          if (! local_systems_cuda_double_uint.empty())
          {
            for (Index i(0) ; i < systems.size() ; ++i)
            {
              auto local_matrix = (systems.at(i)).local().clone();
              Global::synch_matrix(local_matrix, *gates_row.at(i)->_comm, gates_row.at(i)->_ranks, gates_row.at(i)->_mirrors, gates_col.at(i)->_mirrors);
              local_systems_cuda_double_uint.at(i).convert(local_matrix);
            }

            for (Index i(0) ; i < filters_cuda_double_uint.size() ; ++i)
            {
              local_filters_cuda_double_uint.at(i).convert(filters_cuda_double_uint.at(i).local());
            }

            for (Index i(0) ; i < transfers_cuda_double_uint.size() ; ++i)
            {
              local_transfers_cuda_double_uint.at(i).convert(transfers_cuda_double_uint.at(i).local());
            }

          }
#endif
        }

        /// initialise any multigrid hierarchy in use
        void hierarchy_init()
        {
          for (auto& h : local_hierarchy_map_main_float_ulong)
            h.second->init();
          for (auto& h : local_hierarchy_map_main_double_ulong)
            h.second->init();
          for (auto& h : local_hierarchy_map_main_float_uint)
            h.second->init();
          for (auto& h : local_hierarchy_map_main_double_uint)
            h.second->init();
#ifdef FEAT_HAVE_CUDA
          for (auto& h : local_hierarchy_map_cuda_float_ulong)
            h.second->init();
          for (auto& h : local_hierarchy_map_cuda_double_ulong)
            h.second->init();
          for (auto& h : local_hierarchy_map_cuda_float_uint)
            h.second->init();
          for (auto& h : local_hierarchy_map_cuda_double_uint)
            h.second->init();
#endif

          for (auto& h : hierarchy_map_main_float_ulong)
            h.second->init();
          for (auto& h : hierarchy_map_main_double_ulong)
            h.second->init();
          for (auto& h : hierarchy_map_main_float_uint)
            h.second->init();
          for (auto& h : hierarchy_map_main_double_uint)
            h.second->init();
#ifdef FEAT_HAVE_CUDA
          for (auto& h : hierarchy_map_cuda_float_ulong)
            h.second->init();
          for (auto& h : hierarchy_map_cuda_double_ulong)
            h.second->init();
          for (auto& h : hierarchy_map_cuda_float_uint)
            h.second->init();
          for (auto& h : hierarchy_map_cuda_double_uint)
            h.second->init();
#endif
        }

        /// initialise any multigrid hierarchy in use
        void hierarchy_init_symbolic()
        {
          for (auto& h : local_hierarchy_map_main_float_ulong)
            h.second->init_symbolic();
          for (auto& h : local_hierarchy_map_main_double_ulong)
            h.second->init_symbolic();
          for (auto& h : local_hierarchy_map_main_float_uint)
            h.second->init_symbolic();
          for (auto& h : local_hierarchy_map_main_double_uint)
            h.second->init_symbolic();
#ifdef FEAT_HAVE_CUDA
          for (auto& h : local_hierarchy_map_cuda_float_ulong)
            h.second->init_symbolic();
          for (auto& h : local_hierarchy_map_cuda_double_ulong)
            h.second->init_symbolic();
          for (auto& h : local_hierarchy_map_cuda_float_uint)
            h.second->init_symbolic();
          for (auto& h : local_hierarchy_map_cuda_double_uint)
            h.second->init_symbolic();
#endif

          for (auto& h : hierarchy_map_main_float_ulong)
            h.second->init_symbolic();
          for (auto& h : hierarchy_map_main_double_ulong)
            h.second->init_symbolic();
          for (auto& h : hierarchy_map_main_float_uint)
            h.second->init_symbolic();
          for (auto& h : hierarchy_map_main_double_uint)
            h.second->init_symbolic();
#ifdef FEAT_HAVE_CUDA
          for (auto& h : hierarchy_map_cuda_float_ulong)
            h.second->init_symbolic();
          for (auto& h : hierarchy_map_cuda_double_ulong)
            h.second->init_symbolic();
          for (auto& h : hierarchy_map_cuda_float_uint)
            h.second->init_symbolic();
          for (auto& h : hierarchy_map_cuda_double_uint)
            h.second->init_symbolic();
#endif
        }

        /// initialise any multigrid hierarchy in use
        /// \note this includes a call to refresh() to update the actual matrices in use
        void hierarchy_init_numeric()
        {
          this->refresh();

          for (auto& h : local_hierarchy_map_main_float_ulong)
            h.second->init_numeric();
          for (auto& h : local_hierarchy_map_main_double_ulong)
            h.second->init_numeric();
          for (auto& h : local_hierarchy_map_main_float_uint)
            h.second->init_numeric();
          for (auto& h : local_hierarchy_map_main_double_uint)
            h.second->init_numeric();
#ifdef FEAT_HAVE_CUDA
          for (auto& h : local_hierarchy_map_cuda_float_ulong)
            h.second->init_numeric();
          for (auto& h : local_hierarchy_map_cuda_double_ulong)
            h.second->init_numeric();
          for (auto& h : local_hierarchy_map_cuda_float_uint)
            h.second->init_numeric();
          for (auto& h : local_hierarchy_map_cuda_double_uint)
            h.second->init_numeric();
#endif

          for (auto& h : hierarchy_map_main_float_ulong)
            h.second->init_numeric();
          for (auto& h : hierarchy_map_main_double_ulong)
            h.second->init_numeric();
          for (auto& h : hierarchy_map_main_float_uint)
            h.second->init_numeric();
          for (auto& h : hierarchy_map_main_double_uint)
            h.second->init_numeric();
#ifdef FEAT_HAVE_CUDA
          for (auto& h : hierarchy_map_cuda_float_ulong)
            h.second->init_numeric();
          for (auto& h : hierarchy_map_cuda_double_ulong)
            h.second->init_numeric();
          for (auto& h : hierarchy_map_cuda_float_uint)
            h.second->init_numeric();
          for (auto& h : hierarchy_map_cuda_double_uint)
            h.second->init_numeric();
#endif
        }

        /// finish any multigrid hierarchy in use
        void hierarchy_done()
        {
          for (auto& h : hierarchy_map_main_float_ulong)
            h.second->done();
          for (auto& h : hierarchy_map_main_double_ulong)
            h.second->done();
          for (auto& h : hierarchy_map_main_float_uint)
            h.second->done();
          for (auto& h : hierarchy_map_main_double_uint)
            h.second->done();
#ifdef FEAT_HAVE_CUDA
          for (auto& h : hierarchy_map_cuda_float_ulong)
            h.second->done();
          for (auto& h : hierarchy_map_cuda_double_ulong)
            h.second->done();
          for (auto& h : hierarchy_map_cuda_float_uint)
            h.second->done();
          for (auto& h : hierarchy_map_cuda_double_uint)
            h.second->done();
#endif

          for (auto& h : local_hierarchy_map_main_float_ulong)
            h.second->done();
          for (auto& h : local_hierarchy_map_main_double_ulong)
            h.second->done();
          for (auto& h : local_hierarchy_map_main_float_uint)
            h.second->done();
          for (auto& h : local_hierarchy_map_main_double_uint)
            h.second->done();
#ifdef FEAT_HAVE_CUDA
          for (auto& h : local_hierarchy_map_cuda_float_ulong)
            h.second->done();
          for (auto& h : local_hierarchy_map_cuda_double_ulong)
            h.second->done();
          for (auto& h : local_hierarchy_map_cuda_float_uint)
            h.second->done();
          for (auto& h : local_hierarchy_map_cuda_double_uint)
            h.second->done();
#endif
        }

        /// finish any multigrid hierarchy in use
        void hierarchy_done_symbolic()
        {
          for (auto& h : hierarchy_map_main_float_ulong)
            h.second->done_symbolic();
          for (auto& h : hierarchy_map_main_double_ulong)
            h.second->done_symbolic();
          for (auto& h : hierarchy_map_main_float_uint)
            h.second->done_symbolic();
          for (auto& h : hierarchy_map_main_double_uint)
            h.second->done_symbolic();
#ifdef FEAT_HAVE_CUDA
          for (auto& h : hierarchy_map_cuda_float_ulong)
            h.second->done_symbolic();
          for (auto& h : hierarchy_map_cuda_double_ulong)
            h.second->done_symbolic();
          for (auto& h : hierarchy_map_cuda_float_uint)
            h.second->done_symbolic();
          for (auto& h : hierarchy_map_cuda_double_uint)
            h.second->done_symbolic();
#endif

          for (auto& h : local_hierarchy_map_main_float_ulong)
            h.second->done_symbolic();
          for (auto& h : local_hierarchy_map_main_double_ulong)
            h.second->done_symbolic();
          for (auto& h : local_hierarchy_map_main_float_uint)
            h.second->done_symbolic();
          for (auto& h : local_hierarchy_map_main_double_uint)
            h.second->done_symbolic();
#ifdef FEAT_HAVE_CUDA
          for (auto& h : local_hierarchy_map_cuda_float_ulong)
            h.second->done_symbolic();
          for (auto& h : local_hierarchy_map_cuda_double_ulong)
            h.second->done_symbolic();
          for (auto& h : local_hierarchy_map_cuda_float_uint)
            h.second->done_symbolic();
          for (auto& h : local_hierarchy_map_cuda_double_uint)
            h.second->done_symbolic();
#endif
        }

        /// finish any multigrid hierarchy in use
        void hierarchy_done_numeric()
        {
          for (auto& h : hierarchy_map_main_float_ulong)
            h.second->done_numeric();
          for (auto& h : hierarchy_map_main_double_ulong)
            h.second->done_numeric();
          for (auto& h : hierarchy_map_main_float_uint)
            h.second->done_numeric();
          for (auto& h : hierarchy_map_main_double_uint)
            h.second->done_numeric();
#ifdef FEAT_HAVE_CUDA
          for (auto& h : hierarchy_map_cuda_float_ulong)
            h.second->done_numeric();
          for (auto& h : hierarchy_map_cuda_double_ulong)
            h.second->done_numeric();
          for (auto& h : hierarchy_map_cuda_float_uint)
            h.second->done_numeric();
          for (auto& h : hierarchy_map_cuda_double_uint)
            h.second->done_numeric();
#endif

          for (auto& h : local_hierarchy_map_main_float_ulong)
            h.second->done_numeric();
          for (auto& h : local_hierarchy_map_main_double_ulong)
            h.second->done_numeric();
          for (auto& h : local_hierarchy_map_main_float_uint)
            h.second->done_numeric();
          for (auto& h : local_hierarchy_map_main_double_uint)
            h.second->done_numeric();
#ifdef FEAT_HAVE_CUDA
          for (auto& h : local_hierarchy_map_cuda_float_ulong)
            h.second->done_numeric();
          for (auto& h : local_hierarchy_map_cuda_double_ulong)
            h.second->done_numeric();
          for (auto& h : local_hierarchy_map_cuda_float_uint)
            h.second->done_numeric();
          for (auto& h : local_hierarchy_map_cuda_double_uint)
            h.second->done_numeric();
#endif
        }


        template <typename SolverVectorType_>
        std::deque<typename MatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->systems_main_float_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::LocalMatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
          {
            return this->local_systems_main_float_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->systems_main_double_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::LocalMatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
          {
            return this->local_systems_main_double_ulong;
          }

#ifdef FEAT_HAVE_CUDA
        template <typename SolverVectorType_>
        std::deque<typename MatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->systems_cuda_float_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::LocalMatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
          {
            return this->local_systems_cuda_float_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->systems_cuda_double_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::LocalMatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
          {
            return this->local_systems_cuda_double_ulong;
          }

        template <typename MemType_, typename DataType_, typename IndexType_>
          std::deque<typename MatrixType::template ContainerTypeByMDI<MemType_, DataType_, IndexType_> > &
          get_systems(
              typename std::enable_if<std::is_same<MemType_, Mem::CUDA>::value>::type * = nullptr,
              typename std::enable_if<std::is_same<DataType_, double>::value>::type * = nullptr,
              typename std::enable_if<std::is_same<IndexType_, unsigned long>::value>::type * = nullptr)
          {
            return this->systems_cuda_double_ulong;
          }
#endif

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->systems_main_float_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::LocalMatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
          {
            return this->local_systems_main_float_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->systems_main_double_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::LocalMatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
          {
            return this->local_systems_main_double_uint;
          }

#ifdef FEAT_HAVE_CUDA
        template <typename SolverVectorType_>
        std::deque<typename MatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->systems_cuda_float_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::LocalMatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
          {
            return this->local_systems_cuda_float_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->systems_cuda_double_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename MatrixType::LocalMatrixType::template ContainerTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_systems(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
          {
            return this->local_systems_cuda_double_uint;
          }
#endif

        template <typename SolverVectorType_>
        std::deque<typename FilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->filters_main_float_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename FilterType::LocalFilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
          {
            return this->local_filters_main_float_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename FilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->filters_main_double_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename FilterType::LocalFilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
          {
            return this->local_filters_main_double_ulong;
          }

#ifdef FEAT_HAVE_CUDA
        template <typename SolverVectorType_>
        std::deque<typename FilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->filters_cuda_float_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename FilterType::LocalFilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
          {
            return this->local_filters_cuda_float_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename FilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->filters_cuda_double_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename FilterType::LocalFilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
          {
            return this->local_filters_cuda_double_ulong;
          }
#endif

        template <typename SolverVectorType_>
        std::deque<typename FilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->filters_main_float_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename FilterType::LocalFilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
          {
            return this->local_filters_main_float_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename FilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->filters_main_double_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename FilterType::LocalFilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
          {
            return this->local_filters_main_double_uint;
          }

#ifdef FEAT_HAVE_CUDA
        template <typename SolverVectorType_>
        std::deque<typename FilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->filters_cuda_float_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename FilterType::LocalFilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
          {
            return this->local_filters_cuda_float_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename FilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->filters_cuda_double_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename FilterType::LocalFilterType::template FilterTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_filters(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
          {
            return this->local_filters_cuda_double_uint;
          }
#endif

        template <typename SolverVectorType_>
        std::deque<typename TransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->transfers_main_float_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename TransferType::LocalTransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
          {
            return this->local_transfers_main_float_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename TransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->transfers_main_double_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename TransferType::LocalTransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
          {
            return this->local_transfers_main_double_ulong;
          }

#ifdef FEAT_HAVE_CUDA
        template <typename SolverVectorType_>
        std::deque<typename TransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->transfers_cuda_float_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename TransferType::LocalTransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
          {
            return this->local_transfers_cuda_float_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename TransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->transfers_cuda_double_ulong;
          }

        template <typename SolverVectorType_>
        std::deque<typename TransferType::LocalTransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
          {
            return this->local_transfers_cuda_double_ulong;
          }
#endif

        template <typename SolverVectorType_>
        std::deque<typename TransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->transfers_main_float_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename TransferType::LocalTransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
          {
            return this->local_transfers_main_float_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename TransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->transfers_main_double_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename TransferType::LocalTransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
          {
            return this->local_transfers_main_double_uint;
          }

#ifdef FEAT_HAVE_CUDA
        template <typename SolverVectorType_>
        std::deque<typename TransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->transfers_cuda_float_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename TransferType::LocalTransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
          {
            return this->local_transfers_cuda_float_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename TransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
          {
            return this->transfers_cuda_double_uint;
          }

        template <typename SolverVectorType_>
        std::deque<typename TransferType::LocalTransferType::template TransferTypeByMDI<typename SolverVectorType_::MemType,
                                                                    typename SolverVectorType_::DataType,
                                                                    typename SolverVectorType_::IndexType> > &
          get_transfers(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
          {
            return this->local_transfers_cuda_double_uint;
          }
#endif

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_main_float_ulong, FT_main_float_ulong, TT_main_float_ulong> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
        {
          return this->hierarchy_map_main_float_ulong;
        }

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_main_float_ulong::LocalMatrixType,
          typename FT_main_float_ulong::LocalFilterType,
          typename TT_main_float_ulong::LocalTransferType
          > > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
        {
          return this->local_hierarchy_map_main_float_ulong;
        }

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_main_double_ulong, FT_main_double_ulong, TT_main_double_ulong> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
        {
          return this->hierarchy_map_main_double_ulong;
        }

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_main_double_ulong::LocalMatrixType,
          typename FT_main_double_ulong::LocalFilterType,
          typename TT_main_double_ulong::LocalTransferType> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
        {
          return this->local_hierarchy_map_main_double_ulong;
        }

#ifdef FEAT_HAVE_CUDA
        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_cuda_float_ulong, FT_cuda_float_ulong, TT_cuda_float_ulong> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
        {
          return this->hierarchy_map_cuda_float_ulong;
        }

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_cuda_float_ulong::LocalMatrixType,
          typename FT_cuda_float_ulong::LocalFilterType,
          typename TT_cuda_float_ulong::LocalTransferType> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
        {
          return this->local_hierarchy_map_cuda_float_ulong;
        }

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_cuda_double_ulong, FT_cuda_double_ulong, TT_cuda_double_ulong> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              typename SolverVectorType_::GateType *)
        {
          return this->hierarchy_map_cuda_double_ulong;
        }

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_cuda_double_ulong::LocalMatrixType,
          typename FT_cuda_double_ulong::LocalFilterType,
          typename TT_cuda_double_ulong::LocalTransferType> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned long>::value>::type *,
              ...)
        {
          return this->local_hierarchy_map_cuda_double_ulong;
        }
#endif

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_main_float_uint, FT_main_float_uint, TT_main_float_uint> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
        {
          return this->hierarchy_map_main_float_uint;
        }

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_main_float_uint::LocalMatrixType,
          typename FT_main_float_uint::LocalFilterType,
          typename TT_main_float_uint::LocalTransferType> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
        {
          return this->local_hierarchy_map_main_float_uint;
        }

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_main_double_uint, FT_main_double_uint, TT_main_double_uint> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
        {
          return this->hierarchy_map_main_double_uint;
        }

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_main_double_uint::LocalMatrixType,
          typename FT_main_double_uint::LocalFilterType,
          typename TT_main_double_uint::LocalTransferType> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::Main>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
        {
          return this->local_hierarchy_map_main_double_uint;
        }

#ifdef FEAT_HAVE_CUDA
        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_cuda_float_uint, FT_cuda_float_uint, TT_cuda_float_uint> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
        {
          return this->hierarchy_map_cuda_float_uint;
        }

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_cuda_float_uint::LocalMatrixType,
          typename FT_cuda_float_uint::LocalFilterType,
          typename TT_cuda_float_uint::LocalTransferType> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, float>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
        {
          return this->local_hierarchy_map_cuda_float_uint;
        }

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<MT_cuda_double_uint, FT_cuda_double_uint, TT_cuda_double_uint> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              typename SolverVectorType_::GateType *)
        {
          return this->hierarchy_map_cuda_double_uint;
        }

        template <typename SolverVectorType_>
        std::map<String, std::shared_ptr<Solver::MultiGridHierarchy<
          typename MT_cuda_double_uint::LocalMatrixType,
          typename FT_cuda_double_uint::LocalFilterType,
          typename TT_cuda_double_uint::LocalTransferType> > > &
        get_hierarchy_map(
              typename std::enable_if<std::is_same<typename SolverVectorType_::MemType, Mem::CUDA>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::DataType, double>::value>::type *,
              typename std::enable_if<std::is_same<typename SolverVectorType_::IndexType, unsigned int>::value>::type *,
              ...)
        {
          return this->local_hierarchy_map_cuda_double_uint;
        }
#endif

    };
#ifdef FEAT_EICKT
    extern template class MatrixStock<
      Global::Matrix<LAFEM::SparseMatrixCSR<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      Global::Filter<LAFEM::UnitFilter<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      Global::Transfer<LAFEM::Transfer<LAFEM::SparseMatrixCSR<Mem::Main, double, Index>>, LAFEM::VectorMirror<Mem::Main, double, Index>>
        >;

    extern template class MatrixStock<
      Global::Matrix<LAFEM::SparseMatrixBCSR<Mem::Main, double, Index, 2, 2>, LAFEM::VectorMirror<Mem::Main, double, Index>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      Global::Filter<LAFEM::UnitFilterBlocked<Mem::Main, double, Index, 2>, LAFEM::VectorMirror<Mem::Main, double, Index>>,
      Global::Transfer<LAFEM::Transfer<LAFEM::SparseMatrixBCSR<Mem::Main, double, Index, 2, 2>>, LAFEM::VectorMirror<Mem::Main, double, Index>>
        >;
#endif

  }
}

#endif // KERNEL_SOLVER_MATRIX_STOCK_HPP
