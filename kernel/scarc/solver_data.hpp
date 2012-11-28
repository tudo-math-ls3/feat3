#pragma once
#ifndef SCARC_GUARD_SOLVER_DATA_HH
#define SCARC_GUARD_SOLVER_DATA_HH 1

#include<kernel/lafem/dense_vector.hpp>
#include<kernel/lafem/sparse_matrix_csr.hpp>

using namespace FEAST::LAFEM;
using namespace FEAST;

namespace FEAST
{
  namespace ScaRC
  {

    struct SolverDataBase
    {
      virtual const std::string type_name() = 0;
    };

    ///only store data that is *referenced* by solver functors, do not store constants
    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class VectorType_ = DenseVector,
             template<typename, typename> class MatrixType_ = SparseMatrixCSR,
             typename TransferContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             typename PreconContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector>
    struct SolverData : public SolverDataBase
    {
      ///type exports
      typedef VectorType_<MemTag_, DataType_> vector_type_;
      typedef MatrixType_<MemTag_, DataType_> matrix_type_;
      typedef StorageType_<VectorType_<MemTag_, DataType_>, std::allocator<VectorType_<MemTag_, DataType_> > > vector_storage_type_;
      typedef StorageType_<MatrixType_<MemTag_, DataType_>, std::allocator<VectorType_<MemTag_, DataType_> > > matrix_storage_type_;
      typedef StorageType_<TransferContType_, std::allocator<TransferContType_> > transfercont_storage_type_;
      typedef StorageType_<PreconContType_, std::allocator<PreconContType_> > preconcont_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "SolverData";
      };

      ///CTOR from system data
      SolverData(matrix_type_& A,
                 vector_type_& x,
                 vector_type_& b,
                 Index num_temp_vectors = 0) :
        stored_sys(A),
        stored_precon(),
        stored_b(b),
        stored_x(x),
        stored_temp(num_temp_vectors, vector_type_(x.size(), DataType_(0))),
        stored_norm_0(DataType_(0)),
        stored_norm(DataType_(0))
      {
      }

      SolverData(matrix_type_& A,
                 matrix_type_& P,
                 vector_type_& x,
                 vector_type_& b,
                 Index num_temp_vectors = 0) :
        stored_sys(A),
        stored_precon(P),
        stored_b(b),
        stored_x(x),
        stored_temp(num_temp_vectors, vector_type_(x.size(), DataType_(0))),
        stored_norm_0(DataType_(0)),
        stored_norm(DataType_(0))
      {
      }

      ///MemTag_ memory to store matrices and vectors
      matrix_type_ stored_sys;
      PreconContType_ stored_precon;
      vector_type_ stored_b;
      vector_type_ stored_x;
      vector_storage_type_ stored_temp;

      ///host memory to store global scalars
      DataType_ stored_norm_0;
      DataType_ stored_norm;
    };

    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class VectorType_ = DenseVector,
             template<typename, typename> class MatrixType_ = SparseMatrixCSR,
             typename TransferContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             typename PreconContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector>
    struct MultiLevelSolverData : public SolverDataBase
    {
      ///type exports
      typedef StorageType_<VectorType_<MemTag_, DataType_>, std::allocator<VectorType_<MemTag_, DataType_> > > vector_storage_type_;
      typedef StorageType_<MatrixType_<MemTag_, DataType_>, std::allocator<VectorType_<MemTag_, DataType_> > > matrix_storage_type_;
      typedef StorageType_<TransferContType_, std::allocator<TransferContType_> > transfercont_storage_type_;
      typedef StorageType_<PreconContType_, std::allocator<PreconContType_> > preconcont_storage_type_;

      typedef StorageType_<
                StorageType_<VectorType_<MemTag_, DataType_>,
                             std::allocator<VectorType_<MemTag_, DataType_> > >,
                std::allocator<StorageType_<VectorType_<MemTag_, DataType_>, std::allocator<VectorType_<MemTag_, DataType_> > > >
      > vector_storage_storage_type_;

      typedef StorageType_<Index, std::allocator<Index> > index_storage_type_;
      typedef StorageType_<StorageType_<Index, std::allocator<Index> >, std::allocator<StorageType_<Index, std::allocator<Index> > > > index_storage_storage_type_;
      typedef StorageType_<StorageType_<DataType_, std::allocator<DataType_> >, std::allocator<StorageType_<DataType_, std::allocator<DataType_> > > > datatype_storage_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "MultiLevelSolverData";
      };

      ///std CTOR overwrite
      MultiLevelSolverData() :
        stored_sys(),
        stored_res(),
        stored_pro(),
        stored_precon(),
        stored_b(),
        stored_x(),
        stored_d(),
        stored_c(),
        stored_temp(),
        max_iters_subsolver(),
        stored_used_iters_subsolver(),
        eps_relative_subsolver(),
        stored_temp_subsolver(),
        stored_norm_0(DataType_(0)),
        stored_norm(DataType_(0)),
        max_iters(Index(0)),
        num_levels(Index(0)),
        min_level(0),
        eps_relative(DataType_(0))
      {
      }

      ///MemTag_ memory to store matrices and vectors for each level
      matrix_storage_type_ stored_sys;
      transfercont_storage_type_ stored_res;
      transfercont_storage_type_ stored_pro;
      preconcont_storage_type_ stored_precon;
      vector_storage_type_ stored_b;
      vector_storage_type_ stored_x;
      vector_storage_type_ stored_d;
      vector_storage_type_ stored_c;
      vector_storage_type_ stored_temp;

      ///MemTag_ memory to store scalars for each subsolver
      index_storage_storage_type_ max_iters_subsolver; ///max_iters_subsolver[i][j] = k means subsolver number j of level i has max_iter k
      index_storage_storage_type_ stored_used_iters_subsolver;   ///analogously
      datatype_storage_storage_type_ eps_relative_subsolver;

      ///MemTag_ memory to store temp vectors for each subsolver
      vector_storage_storage_type_ stored_temp_subsolver;

      ///host memory to store global scalars
      DataType_ stored_norm_0;
      DataType_ stored_norm;

      ///global constants
      Index num_levels;
      Index min_level;
      Index max_iters;
      DataType_ eps_relative;
    };
  }
}

#endif
