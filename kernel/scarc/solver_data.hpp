#pragma once
#ifndef SCARC_GUARD_SOLVER_DATA_HH
#define SCARC_GUARD_SOLVER_DATA_HH 1

#include<kernel/lafem/dense_vector.hpp>
#include<kernel/lafem/vector_mirror.hpp>
#include<kernel/lafem/sparse_matrix_csr.hpp>
#include<kernel/lafem/unit_filter.hpp>
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

    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class PreconContType_ = SparseMatrixCSR>
    struct PreconditionerDataContainer
    {
      public:
        typedef PreconContType_<MemTag_, DataType_> precon_type_;

        virtual precon_type_& precon()
        {
          return _stored_precon;
        }

        virtual const precon_type_& precon() const
        {
          return _stored_precon;
        }

      protected:
        ///CTORs to be used in subclasses
        PreconditionerDataContainer(const precon_type_& precon) :
          _stored_precon(precon)
        {
        }

        PreconditionerDataContainer(const PreconditionerDataContainer& other)
        {
          this->_stored_precon = other._stored_precon;
        }

        precon_type_ _stored_precon;
    };

    ///only store data that is *referenced* by solver functors, do not store constants
    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class VectorType_ = DenseVector,
             template<typename, typename> class MatrixType_ = SparseMatrixCSR,
             template<typename, typename> class PreconContType_ = SparseMatrixCSR,
             template<typename, typename> class StorageType_ = std::vector>
    struct PreconditionedSolverData :
      public SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>,
      public PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>
    {
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::matrix_type_ matrix_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_type_ vector_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_storage_type_ vector_storage_type_;

      typedef PreconContType_<MemTag_, DataType_> precon_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "PreconditionedSolverData";
      }

      ///CTOR from system data
      PreconditionedSolverData(matrix_type_& A,
                               precon_type_& P,
                               vector_type_& x,
                               vector_type_& b,
                               Index num_temp_vectors = 0,
                               Index num_temp_scalars = 0) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(A, x, b, num_temp_vectors, num_temp_scalars),
        PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>(P)
      {
      }

      ///copy CTOR
      PreconditionedSolverData(const PreconditionedSolverData& other) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(other),
        PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>(other)
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

        this->_stored_precon = other._stored_precon;

        return *this;
      }
    };

    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class VectorType_ = DenseVector,
             template<typename, typename> class VectorMirrorType_ = VectorMirror,
             template<typename, typename> class StorageType_ = std::vector>
    struct SynchronizationDataContainer
    {
      public:
        typedef VectorType_<MemTag_, DataType_> vector_type_;
        typedef VectorMirrorType_<MemTag_, DataType_> vector_mirror_type_;
        typedef VectorType_<MemTag_, DataType_> vector_storage_type_;
        typedef StorageType_<VectorMirrorType_<MemTag_, DataType_>, std::allocator<VectorMirrorType_<MemTag_, DataType_> > > vector_mirror_storage_type_;
        typedef StorageType_<Index, std::allocator<Index> > index_storage_type_;

        virtual vector_mirror_storage_type_& vector_mirrors()
        {
          return _stored_vector_mirrors;
        }
        virtual const vector_mirror_storage_type_& vector_mirrors() const
        {
          return _stored_vector_mirrors;
        }

        virtual vector_storage_type_& vector_mirror_sendbufs()
        {
          return _stored_vector_mirror_sendbufs;
        }
        virtual const vector_storage_type_& vector_mirror_sendbufs() const
        {
          return _stored_vector_mirror_sendbufs;
        }

        virtual vector_storage_type_& vector_mirror_recvbufs()
        {
          return _stored_vector_mirror_recvbufs;
        }
        virtual const vector_storage_type_& vector_mirror_recvbufs() const
        {
          return _stored_vector_mirror_recvbufs;
        }

        virtual index_storage_type_& dest_ranks()
        {
          return _stored_dest_ranks;
        }
        virtual const index_storage_type_& dest_ranks() const
        {
          return _stored_dest_ranks;
        }

        virtual index_storage_type_& source_ranks()
        {
          return _stored_source_ranks;
        }
        virtual const index_storage_type_& source_ranks() const
        {
          return _stored_source_ranks;
        }

      protected:
        ///CTORs to be used in subclasses
        SynchronizationDataContainer() :
          _stored_vector_mirrors(),
          _stored_vector_mirror_sendbufs(),
          _stored_vector_mirror_recvbufs(),
          _stored_dest_ranks(),
          _stored_source_ranks()
        {
        }

        SynchronizationDataContainer(const SynchronizationDataContainer& other)
        {
          this->_stored_vector_mirrors = other._stored_vector_mirrors;
          this->_stored_vector_mirror_sendbufs = other._stored_vector_mirror_sendbufs;
          this->_stored_vector_mirror_recvbufs = other._stored_vector_mirror_recvbufs;
          this->_stored_dest_ranks = other.stored_dest_ranks;
          this->_stored_source_ranks = other.stored_source_ranks;
        }

        vector_mirror_storage_type_ _stored_vector_mirrors;
        vector_storage_type_ _stored_vector_mirror_sendbufs;
        vector_storage_type_ _stored_vector_mirror_recvbufs;
        index_storage_type_ _stored_dest_ranks;
        index_storage_type_ _stored_source_ranks;
    };

    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class VectorType_ = DenseVector,
             template<typename, typename> class VectorMirrorType_ = VectorMirror,
             template<typename, typename> class MatrixType_ = SparseMatrixCSR,
             template<typename, typename> class StorageType_ = std::vector>
    struct SynchronisedSolverData :
      public SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>,
      public SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>
    {
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::matrix_type_ matrix_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_type_ vector_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_storage_type_ vector_storage_type_;
      typedef StorageType_<VectorMirrorType_<MemTag_, DataType_>, std::allocator<VectorMirrorType_<MemTag_, DataType_> > > vector_mirror_storage_type_;
      typedef StorageType_<Index, std::allocator<Index> > index_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "SynchronisedSolverData";
      }

      ///CTOR from system data
      SynchronisedSolverData(matrix_type_& A,
                             vector_type_& x,
                             vector_type_& b,
                             Index num_temp_vectors = 0,
                             Index num_temp_scalars = 0) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(A, x, b, num_temp_vectors, num_temp_scalars),
        SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>()
      {
      }

      ///copy CTOR
      SynchronisedSolverData(const SynchronisedSolverData& other) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(other),
        SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>(other)
      {
      }

      ///assignment operator overload
      SynchronisedSolverData& operator=(const SynchronisedSolverData& other)
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

        this->_stored_vector_mirrors = other._stored_vector_mirrors;
        this->_stored_vector_mirror_sendbufs = other._stored_vector_mirror_sendbufs;
        this->_stored_vector_mirror_recvbufs = other._stored_vector_mirror_recvbufs;
        this->_stored_dest_ranks = other._stored_dest_ranks;
        this->_stored_source_ranks = other._stored_source_ranks;

        return *this;
      }
    };
/*
    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class VectorType_ = DenseVector,
             template<typename, typename> class VectorMirrorType_ = VectorMirror,
             template<typename, typename> class MatrixType_ = SparseMatrixCSR,
             typename PreconContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector>
    struct SynchronisedPreconditionedSolverData : public SynchronisedSolverData<DataType_, MemTag_, VectorType_, VectorMirrorType_, MatrixType_, StorageType_>
    {
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::matrix_type_ matrix_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_type_ vector_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_storage_type_ vector_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "SynchronisedPreconditionedSolverData";
      }

      ///CTOR from system data
      SynchronisedPreconditionedSolverData(matrix_type_& A,
                             PreconContType_& P,
                             vector_type_& x,
                             vector_type_& b,
                             Index num_temp_vectors = 0,
                             Index num_temp_scalars = 0) :
        SynchronisedSolverData<DataType_, MemTag_, VectorType_, VectorMirrorType_, MatrixType_, StorageType_>(A, x, b, num_temp_vectors, num_temp_scalars),
        stored_prec(P)
      {
      }

      ///copy CTOR
      SynchronisedPreconditionedSolverData(const SynchronisedPreconditionedSolverData& other) :
        SynchronisedSolverData<DataType_, MemTag_, VectorType_, VectorMirrorType_, MatrixType_, StorageType_>(other),
        stored_prec(other.stored_prec)
      {
      }

      ///assignment operator overload
      SynchronisedPreconditionedSolverData& operator=(const SynchronisedPreconditionedSolverData& other)
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

        this->stored_mirrors = other.stored_mirrors;
        this->stored_mirror_sendbufs = other.stored_mirror_sendbufs;
        this->stored_mirror_recvbufs = other.stored_mirror_recvbufs;
        this->stored_dest_ranks = other.stored_dest_ranks;
        this->stored_source_ranks = other.stored_source_ranks;

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
    };*/
  }
}

#endif
