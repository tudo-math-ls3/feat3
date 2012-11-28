#pragma once
#ifndef SCARC_GUARD_SOLVER_DATA_HH
#define SCARC_GUARD_SOLVER_DATA_HH 1

#include<kernel/lafem/dense_vector.hpp>
#include<kernel/lafem/sparse_matrix_csr.hpp>
#include<kernel/util/cpp11_smart_pointer.hpp>

using namespace FEAST::LAFEM;
using namespace FEAST;

namespace FEAST
{
  namespace ScaRC
  {

    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class VectorType_ = DenseVector,
             template<typename, typename> class MatrixType_ = SparseMatrixCSR,
             template<typename, typename> class StorageType_ = std::vector>
    struct SolverDataBase
    {
      ///type exports
      typedef VectorType_<MemTag_, DataType_> vector_type_;
      typedef MatrixType_<MemTag_, DataType_> matrix_type_;
      typedef StorageType_<VectorType_<MemTag_, DataType_>, std::allocator<VectorType_<MemTag_, DataType_> > > vector_storage_type_;
      typedef StorageType_<MatrixType_<MemTag_, DataType_>, std::allocator<VectorType_<MemTag_, DataType_> > > matrix_storage_type_;

      virtual const std::string type_name() = 0;

      ///MemTag_ memory to store matrices and vectors
      matrix_type_ stored_sys;
      vector_type_ stored_rhs;
      vector_type_ stored_sol;
      vector_storage_type_ stored_temp;

      ///host memory to store global scalars
      DataType_ stored_norm_0;
      DataType_ stored_norm;

      virtual ~SolverDataBase()
      {
      }
    };

    ///only store data that is *referenced* by solver functors, do not store constants
    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class VectorType_ = DenseVector,
             template<typename, typename> class MatrixType_ = SparseMatrixCSR,
             template<typename, typename> class StorageType_ = std::vector>
    struct SolverData : public SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>
    {
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::matrix_type_ matrix_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_type_ vector_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_storage_type_ vector_storage_type_;
      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "SolverData";
      }

      ///CTOR from system data
      SolverData(matrix_type_& A,
                 vector_type_& x,
                 vector_type_& b,
                 Index num_temp_vectors = 0)
      {
        this->stored_sys = A;
        this->stored_rhs = b;
        this->stored_sol = x;
        this->stored_temp = vector_storage_type_(num_temp_vectors, vector_type_(x.size()));
        this->stored_norm_0 = DataType_(0);
        this->stored_norm = DataType_(0);
      }

      ///copy CTOR
      SolverData(const SolverData& other)
      {
        this->stored_sys = other.stored_sys;
        this->stored_rhs = other.stored_rhs;
        this->stored_sol = other.stored_sol;
        this->stored_temp = other.stored_temp;
        this->stored_norm_0 = other.stored_norm_0;
        this->stored_norm = other.stored_norm;
      }

      ///assignment operator overload
      SolverData& operator=(const SolverData& other)
      {
          if(this == &other)
            return *this;

        this->stored_sys = other.stored_sys;
        this->stored_rhs = other.stored_rhs;
        this->stored_sol = other.stored_sol;
        this->stored_temp = other.stored_temp;
        this->stored_norm_0 = other.stored_norm_0;
        this->stored_norm = other.stored_norm;

        return *this;
      }

      ///CTOR from any SolverDataBase
      template<typename DT_,
               typename Tag_,
               template<typename, typename> class VT_,
               template<typename, typename> class MT_,
               template<typename, typename> class StoreT_>
      SolverData(const SolverDataBase<DT_, Tag_, VT_, MT_, StoreT_>& other)
      {
        this->stored_sys = other.stored_sys;
        this->stored_rhs = other.stored_rhs;
        this->stored_sol = other.stored_sol;
        this->stored_temp = other.stored_temp;
        this->stored_norm_0 = other.stored_norm_0;
        this->stored_norm = other.stored_norm;
      }

    };

    ///only store data that is *referenced* by solver functors, do not store constants
    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class VectorType_ = DenseVector,
             template<typename, typename> class MatrixType_ = SparseMatrixCSR,
             typename PreconContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector>
    struct PreconditionedSolverData : public SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>
    {
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::matrix_type_ matrix_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_type_ vector_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_storage_type_ vector_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "PreconditionedSolverData";
      }

      ///CTOR from system data
      PreconditionedSolverData(matrix_type_& A,
                              matrix_type_& P,
                              vector_type_& x,
                              vector_type_& b,
                              Index num_temp_vectors = 0) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(A, x, b, num_temp_vectors),
        stored_prec(P)
      {
      }

      ///copy CTOR
      PreconditionedSolverData(const PreconditionedSolverData& other) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(other),
        stored_prec(other.stored_prec)
      {
      }

      ///assignment operator overload
      PreconditionedSolverData& operator=(const PreconditionedSolverData& other)
      {
          if(this == &other)
            return *this;

        this->stored_sys = other.stored_sys;
        this->stored_rhs = other.stored_rhs;
        this->stored_sol = other.stored_sol;
        this->stored_temp = other.stored_temp;
        this->stored_norm_0 = other.stored_norm_0;
        this->stored_norm = other.stored_norm;

        this->stored_prec = other.stored_prec;

        return *this;
      }

      PreconContType_ stored_prec;
    };

    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class VectorType_ = DenseVector,
             template<typename, typename> class MatrixType_ = SparseMatrixCSR,
             typename TransferContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             typename PreconContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector>
    struct MultiLevelSolverData : public SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>
    {
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::matrix_type_ matrix_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_type_ vector_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_storage_type_ vector_storage_type_;

      typedef StorageType_<std::shared_ptr<SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_> >, std::allocator<std::shared_ptr<SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_> > > > leveldata_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "MultiLevelSolverData";
      }

      ///CTOR from system
      MultiLevelSolverData(matrix_type_& A,
                           vector_type_& x,
                           vector_type_& b,
                           Index num_temp_vectors = 0) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(A, x, b, num_temp_vectors),
        stored_level_data()
      {
        ///TODO
      }

      ///copy CTOR
      MultiLevelSolverData(const MultiLevelSolverData& other) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(other),
        stored_level_data(other.stored_level_data)
      {
      }

      ///assignment operator overload
      MultiLevelSolverData& operator=(const MultiLevelSolverData& other)
      {
          if(this == &other)
            return *this;

        this->stored_sys = other.stored_sys;
        this->stored_rhs = other.stored_rhs;
        this->stored_sol = other.stored_sol;
        this->stored_temp = other.stored_temp;
        this->stored_norm_0 = other.stored_norm_0;
        this->stored_norm = other.stored_norm;

        this->stored_level_data = other.stored_level_data;

        return *this;
      }

      leveldata_storage_type_ stored_level_data;
    };
  }
}

#endif
