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
      public:
        ///type exports
        typedef VectorType_<MemTag_, DataType_> vector_type_;
        typedef MatrixType_<MemTag_, DataType_> matrix_type_;
        typedef StorageType_<VectorType_<MemTag_, DataType_>, std::allocator<VectorType_<MemTag_, DataType_> > > vector_storage_type_;
        typedef StorageType_<MatrixType_<MemTag_, DataType_>, std::allocator<VectorType_<MemTag_, DataType_> > > matrix_storage_type_;
        typedef StorageType_<DataType_, std::allocator<DataType_> > scalar_storage_type_;

        virtual const std::string type_name() = 0;

        virtual matrix_type_& sys()
        {
          return this->_stored_sys;
        }

        virtual const matrix_type_& sys() const
        {
          return this->_stored_sys;
        }

        virtual vector_type_& rhs()
        {
          return this->_stored_rhs;
        }

        virtual const vector_type_& rhs() const
        {
          return this->_stored_rhs;
        }

        virtual vector_type_& sol()
        {
          return this->_stored_sol;
        }

        virtual const vector_type_& sol() const
        {
          return this->_stored_sol;
        }

        virtual vector_storage_type_& temp()
        {
          return _stored_temp;
        }

        virtual const vector_storage_type_& temp() const
        {
          return _stored_temp;
        }

        virtual DataType_& norm_0()
        {
          return _stored_norm_0;
        }

        virtual const DataType_& norm_0() const
        {
          return _stored_norm_0;
        }

        virtual DataType_& norm()
        {
          return _stored_norm;
        }

        virtual const DataType_& norm() const
        {
          return _stored_norm;
        }

        virtual scalar_storage_type_& scalars()
        {
          return _stored_scalars;
        }

        virtual const scalar_storage_type_& scalars() const
        {
          return _stored_scalars;
        }

        virtual DataType_& eps()
        {
          return _stored_eps;
        }

        virtual const DataType_& eps() const
        {
          return _stored_eps;
        }

        virtual Index& max_iters()
        {
          return _stored_max_iters;
        }

        virtual const Index& max_iters() const
        {
          return _stored_max_iters;
        }

        virtual Index& used_iters()
        {
          return _stored_used_iters;
        }

        virtual const Index& used_iters() const
        {
          return _stored_used_iters;
        }

        virtual ~SolverDataBase()
        {
        }

      protected:
        ///MemTag_ memory to store matrices and vectors and scalars
        matrix_type_ _stored_sys;
        vector_type_ _stored_rhs;
        vector_type_ _stored_sol;
        vector_storage_type_ _stored_temp;
        scalar_storage_type_ _stored_scalars;

        ///host memory to store global scalars
        DataType_ _stored_norm_0;
        DataType_ _stored_norm;
        DataType_ _stored_eps;
        Index _stored_max_iters;
        Index _stored_used_iters;

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
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::scalar_storage_type_ scalar_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "SolverData";
      }

      ///CTOR from system data
      SolverData(matrix_type_& A,
                 vector_type_& x,
                 vector_type_& b,
                 Index num_temp_vectors = 0,
                 Index num_temp_scalars = 0)
      {
        this->_stored_sys = A;
        this->_stored_rhs = b;
        this->_stored_sol = x;
        this->_stored_temp = vector_storage_type_(num_temp_vectors, vector_type_(x.size()));
        this->_stored_scalars = scalar_storage_type_(num_temp_scalars, DataType_(0));
        this->_stored_norm_0 = DataType_(0);
        this->_stored_norm = DataType_(0);
        this->_stored_eps = DataType_(0);
        this->_stored_max_iters = Index(0);
        this->_stored_used_iters = Index(0);
      }

      ///copy CTOR
      SolverData(const SolverData& other)
      {
        this->_stored_sys = other._stored_sys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_norm_0 = other._stored_norm_0;
        this->_stored_norm = other._stored_norm;
        this->_stored_eps = other._stored_eps;
        this->_stored_max_iters = Index(0);
        this->_stored_used_iters = Index(0);
      }

      ///assignment operator overload
      SolverData& operator=(const SolverData& other)
      {
          if(this == &other)
            return *this;

        this->_stored_sys = other._stored_sys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_norm_0 = other._stored_norm_0;
        this->_stored_norm = other._stored_norm;
        this->_stored_eps = other._stored_eps;
        this->_stored_max_iters = Index(0);
        this->_stored_used_iters = Index(0);

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
        this->_stored_sys = other._stored_sys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_norm_0 = other._stored_norm_0;
        this->_stored_norm = other._stored_norm;
        this->_stored_eps = other._stored_eps;
        this->_stored_max_iters = Index(0);
        this->_stored_used_iters = Index(0);
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
                              Index num_temp_vectors = 0,
                              Index num_temp_scalars = 0) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(A, x, b, num_temp_vectors, num_temp_scalars),
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

        this->_stored_sys = other._stored_sys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_norm_0 = other._stored_norm_0;
        this->_stored_norm = other._stored_norm;
        this->_stored_eps = other._stored_eps;
        this->_stored_max_iters = Index(0);
        this->_stored_used_iters = Index(0);

        this->stored_prec = other._stored_prec;

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
                           Index num_temp_vectors = 0,
                           Index num_temp_scalars = 0) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(A, x, b, num_temp_vectors, num_temp_scalars),
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

        this->_stored_sys = other._stored_sys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_norm_0 = other._stored_norm_0;
        this->_stored_norm = other._stored_norm;
        this->_stored_eps = other._stored_eps;
        this->_stored_max_iters = Index(0);
        this->_stored_used_iters = Index(0);

        this->stored_level_data = other._stored_level_data;

        return *this;
      }

      leveldata_storage_type_ stored_level_data;
    };
  }
}

#endif
