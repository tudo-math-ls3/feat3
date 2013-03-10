/**
 * \file
 * \brief FEAST milestone 1 ScaRC data implementations
 * \author Markus Geveler
 * \date 2012 - 2013
 *
 * See class documentation.
 *
 */

#pragma once
#ifndef SCARC_GUARD_SOLVER_DATA_HH
#define SCARC_GUARD_SOLVER_DATA_HH 1

#include<kernel/lafem/dense_vector.hpp>
#include<kernel/lafem/vector_mirror.hpp>
#include<kernel/lafem/sparse_matrix_csr.hpp>
#include<kernel/lafem/unit_filter.hpp>
#include<kernel/util/cpp11_smart_pointer.hpp>

#include<kernel/foundation/communication.hpp>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::Foundation;

namespace FEAST
{
  namespace ScaRC
  {

    /**
     * \brief Base class for all solver data containers
     *
     * Solver data containers rely on singly (per value) or multiple (per STL) stored data.
     * Solver data conatainers are referenced by solver pattern creation functions.
     *
     * \tparam DataType_
     * data type
     *
     * \tparam MemTag_
     * memory tag
     *
     * \tparam VectorType_
     * vector type template
     *
     * \tparam MatrixType_
     * matrix type template
     *
     * \tparam StorageType_
     * storage type template (STL or -conformal)
     *
     * \author Markus Geveler
     */
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
        typedef StorageType_<Index, std::allocator<Index> > index_storage_type_;

        virtual const std::string type_name() = 0;

        virtual matrix_type_& sys()
        {
          return this->_stored_sys;
        }

        virtual const matrix_type_& sys() const
        {
          return this->_stored_sys;
        }

        virtual matrix_type_& localsys()
        {
          return this->_stored_localsys;
        }

        virtual const matrix_type_& localsys() const
        {
          return this->_stored_localsys;
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

        virtual vector_type_& def()
        {
          return this->_stored_def;
        }

        virtual const vector_type_& def() const
        {
          return this->_stored_def;
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

        virtual index_storage_type_& indices()
        {
          return _stored_indices;
        }

        virtual const index_storage_type_& indices() const
        {
          return _stored_indices;
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
        matrix_type_ _stored_localsys;
        vector_type_ _stored_rhs;
        vector_type_ _stored_sol;
        vector_type_ _stored_def;
        vector_storage_type_ _stored_temp;
        scalar_storage_type_ _stored_scalars;
        index_storage_type_ _stored_indices;

        ///host memory to store global scalars
        DataType_ _stored_norm_0;
        DataType_ _stored_norm;
        DataType_ _stored_eps;
        Index _stored_max_iters;
        Index _stored_used_iters;

    };

    /**
     * \brief Simplemost data container implementation
     *
     * See bse class.
     *
     * \tparam DataType_
     * data type
     *
     * \tparam MemTag_
     * memory tag
     *
     * \tparam VectorType_
     * vector type template
     *
     * \tparam MatrixType_
     * matrix type template
     *
     * \tparam StorageType_
     * storage type template (STL or -conformal)
     *
     * \author Markus Geveler
     */
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
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::index_storage_type_ index_storage_type_;

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
                 Index num_temp_scalars = 0,
                 Index num_temp_indices = 0)
      {
        this->_stored_sys = A;
        this->_stored_rhs = b;
        this->_stored_sol = x;
        this->_stored_def = vector_type_(x.size());
        this->_stored_temp = vector_storage_type_(num_temp_vectors, vector_type_(x.size(), DataType_(0)));
        this->_stored_scalars = scalar_storage_type_(num_temp_scalars, DataType_(0));
        this->_stored_indices = index_storage_type_(num_temp_indices, Index(0));
        this->_stored_norm_0 = DataType_(1000);
        this->_stored_norm = DataType_(1000);
        this->_stored_eps = DataType_(1e-8);
        this->_stored_max_iters = Index(1000);
        this->_stored_used_iters = Index(0);
      }

      ///copy CTOR
      SolverData(const SolverData& other)
      {
        this->_stored_sys = other._stored_sys;
        this->_stored_localsys = other._stored_localsys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_def = other._stored_def;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_indices = other._stored_indices;
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
        this->_stored_localsys = other._stored_localsys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_def = other._stored_def;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_indices = other._stored_indices;
        this->_stored_norm_0 = other._stored_norm_0;
        this->_stored_norm = other._stored_norm;
        this->_stored_eps = other._stored_eps;
        this->_stored_max_iters = Index(0);
        this->_stored_used_iters = Index(0);

        return *this;
      }

      ///CTOR from any SolverDataBase
      /*template<typename DT_,
               typename Tag_,
               template<typename, typename> class VT_,
               template<typename, typename> class MT_,
               template<typename, typename> class StoreT_>
      SolverData(const SolverDataBase<DT_, Tag_, VT_, MT_, StoreT_>& other)
      {
        this->_stored_sys = other._stored_sys;
        this->_stored_localsys = other._stored_localsys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_def = other._stored_def;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_indices = other._stored_indices;
        this->_stored_norm_0 = other._stored_norm_0;
        this->_stored_norm = other._stored_norm;
        this->_stored_eps = other._stored_eps;
        this->_stored_max_iters = Index(0);
        this->_stored_used_iters = Index(0);
      }*/

    };

    /**
     * \brief Preconditioner data interface
     *
     * Can be subclassed in addition to SolverDataBase (or subclasses thereof) in order to store
     * preconditioner data.
     *
     * \tparam DataType_
     * data type
     *
     * \tparam MemTag_
     * memory tag
     *
     * \tparam PreconContType_
     * storage type template for preconditioner data
     *
     * \author Markus Geveler
     */
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

    /**
     * \brief Simple data container with preconditioner data
     *
     * PreconditionedSolverData combines SoverData and PreconditionerData.
     * See base classes.
     *
     * \tparam DataType_
     * data type
     *
     * \tparam MemTag_
     * memory tag
     *
     * \tparam VectorType_
     * vector type template
     *
     * \tparam MatrixType_
     * matrix type template
     *
     * \tparam StorageType_
     * storage type template (STL or -conformal)
     *
     * \author Markus Geveler
     */
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
                               Index num_temp_scalars = 0,
                               Index num_temp_indices = 0) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(A, x, b, num_temp_vectors, num_temp_scalars, num_temp_indices),
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
        this->_stored_localsys = other._stored_localsys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_def = other._stored_def;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_indices = other._stored_indices;
        this->_stored_norm_0 = other._stored_norm_0;
        this->_stored_norm = other._stored_norm;
        this->_stored_eps = other._stored_eps;
        this->_stored_max_iters = Index(0);
        this->_stored_used_iters = Index(0);

        this->_stored_precon = other._stored_precon;

        return *this;
      }
    };

    /**
     * \brief Synchronisation data interface
     *
     * Can be subclassed in addition to SolverDataBase (or subclasses thereof) in order to store
     * synchronisation data (like buffers and mirrors).
     *
     * \tparam DataType_
     * data type
     *
     * \tparam MemTag_
     * memory tag
     *
     * \tparam PreconContType_
     * storage type template for preconditioner data
     *
     * \author Markus Geveler
     */
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
        typedef StorageType_<VectorType_<MemTag_, DataType_>, std::allocator<VectorType_<MemTag_, DataType_> > > vector_storage_type_;
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

    /**
     * \brief Simple data container with synchronisation data
     *
     * SynchronisedSolverData combines SoverData and SynchronizationDataContainer.
     * See base classes.
     *
     * \tparam DataType_
     * data type
     *
     * \tparam MemTag_
     * memory tag
     *
     * \tparam VectorType_
     * vector type template
     *
     * \tparam VectorMirrorType_
     * vector mirror type template
     *
     * \tparam MatrixType_
     * matrix type template
     *
     * \tparam StorageType_
     * storage type template (STL or -conformal)
     *
     * \author Markus Geveler
     */
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
                             Index num_temp_scalars = 0,
                             Index num_temp_indices = 0) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(A, x, b, num_temp_vectors, num_temp_scalars, num_temp_indices),
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
        this->_stored_localsys = other._stored_localsys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_def = other._stored_def;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_indices = other._stored_indices;
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

    /**
     * \brief Data container with synchronisation and preconditioner data
     *
     * See base classes.
     *
     * \tparam DataType_
     * data type
     *
     * \tparam MemTag_
     * memory tag
     *
     * \tparam VectorType_
     * vector type template
     *
     * \tparam VectorMirrorType_
     * vector mirror type template
     *
     * \tparam MatrixType_
     * matrix type template
     *
     * \tparam PreconContType_
     * preconditioner type template
     *
     * \tparam StorageType_
     * storage type template (STL or -conformal)
     *
     * \author Markus Geveler
     */
    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class VectorType_ = DenseVector,
             template<typename, typename> class VectorMirrorType_ = VectorMirror,
             template<typename, typename> class MatrixType_ = SparseMatrixCSR,
             template<typename, typename> class PreconContType_ = SparseMatrixCSR,
             template<typename, typename> class StorageType_ = std::vector>
    struct SynchronisedPreconditionedSolverData :
      public SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>,
      public PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>,
      public SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>
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
                             PreconContType_<MemTag_, DataType_>& P,
                             vector_type_& x,
                             vector_type_& b,
                             Index num_temp_vectors = 0,
                             Index num_temp_scalars = 0,
                             Index num_temp_indices = 0) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(A, x, b, num_temp_vectors, num_temp_scalars, num_temp_indices),
        PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>(P),
        SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>()

      {
      }

      ///copy CTOR
      SynchronisedPreconditionedSolverData(const SynchronisedPreconditionedSolverData& other) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(other),
        PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>(other),
        SynchronisedSolverData<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>(other)
      {
      }

      ///assignment operator overload
      SynchronisedPreconditionedSolverData& operator=(const SynchronisedPreconditionedSolverData& other)
      {
          if(this == &other)
            return *this;

        this->_stored_sys = other._stored_sys;
        this->_stored_localsys = other._stored_localsys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_def = other._stored_def;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_indices = other._stored_indices;
        this->_stored_norm_0 = other._stored_norm_0;
        this->_stored_norm = other._stored_norm;
        this->_stored_eps = other._stored_eps;
        this->_stored_max_iters = Index(0);
        this->_stored_used_iters = Index(0);

        this->_stored_precon = other._stored_precon;

        this->_stored_vector_mirrors = other._stored_vector_mirrors;
        this->_stored_vector_mirror_sendbufs = other._stored_vector_mirror_sendbufs;
        this->_stored_vector_mirror_recvbufs = other._stored_vector_mirror_recvbufs;
        this->_stored_dest_ranks = other._stored_dest_ranks;
        this->_stored_source_ranks = other._stored_source_ranks;

        return *this;
      }
    };

    ///multilevel solver data draft - NOTE: not subject to milestone 1 (v0.1)
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
                           Index num_temp_scalars = 0,
                           Index num_temp_indices = 0) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(A, x, b, num_temp_vectors, num_temp_scalars, num_temp_indices),
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
        this->_stored_localsys = other._stored_localsys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_def = other._stored_def;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_indices = other._stored_indices;
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

    /**
     * \brief Filter data interface
     *
     * Can be subclassed in addition to SolverDataBase (subclasses thereof) in order to store
     * filter data.
     *
     * \tparam DataType_
     * data type
     *
     * \tparam MemTag_
     * memory tag
     *
     * \tparam FilterType_
     * storage type template for filter data
     *
     * \author Markus Geveler
     */
    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class FilterType_ = UnitFilter>
    struct FilterDataContainer
    {
      public:
        typedef FilterType_<MemTag_, DataType_> filter_type_;

        virtual filter_type_& filter()
        {
          return _stored_filter;
        }

        virtual const filter_type_& filter() const
        {
          return _stored_filter;
        }

      protected:
        ///CTORs to be used in subclasses
        FilterDataContainer(const filter_type_& filter) :
          _stored_filter(filter)
        {
        }

        FilterDataContainer(const FilterDataContainer& other)
        {
          this->_stored_filter = other._stored_filter;
        }

        filter_type_ _stored_filter;
    };

    /**
     * \brief Complex container with filter data
     *
     * See base classes.
     *
     * \tparam DataType_
     * data type
     *
     * \tparam MemTag_
     * memory tag
     *
     * \tparam VectorType_
     * vector type template
     *
     * \tparam VectorMirrorType_
     * vector mirror type template
     *
     * \tparam MatrixType_
     * matrix type template
     *
     * \tparam PreconContType_
     * preconditioner data type template
     *
     * \tparam FilterType_
     * filter data type template
     *
     * \tparam StorageType_
     * storage type template (STL or -conformal)
     *
     * \author Markus Geveler
     */
    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             template<typename, typename> class VectorType_ = DenseVector,
             template<typename, typename> class VectorMirrorType_ = VectorMirror,
             template<typename, typename> class MatrixType_ = SparseMatrixCSR,
             template<typename, typename> class PreconContType_ = SparseMatrixCSR,
             template<typename, typename> class FilterType_ = UnitFilter,
             template<typename, typename> class StorageType_ = std::vector>
    struct SynchronisedPreconditionedFilteredSolverData :
      public SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>,
      public FilterDataContainer<DataType_, MemTag_, FilterType_>,
      public PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>,
      public SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>
    {
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::matrix_type_ matrix_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_type_ vector_type_;
      typedef typename SolverDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>::vector_storage_type_ vector_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "SynchronisedPreconditionedFilteredSolverData";
      }

      ///CTOR from system data
      SynchronisedPreconditionedFilteredSolverData(matrix_type_& A,
                                                   PreconContType_<MemTag_, DataType_>& P,
                                                   vector_type_& x,
                                                   vector_type_& b,
                                                   FilterType_<MemTag_, DataType_>& filter,
                                                   Index num_temp_vectors = 0,
                                                   Index num_temp_scalars = 0,
                                                   Index num_temp_indices = 0) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(A, x, b, num_temp_vectors, num_temp_scalars, num_temp_indices),
        FilterDataContainer<DataType_, MemTag_, FilterType_>(filter),
        PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>(P),
        SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>()
      {
      }

      ///copy CTOR
      SynchronisedPreconditionedFilteredSolverData(const SynchronisedPreconditionedFilteredSolverData& other) :
        SolverData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(other),
        FilterDataContainer<DataType_, MemTag_, FilterType_>(other),
        PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>(other),
        SynchronisedSolverData<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>(other)
      {
      }

      ///assignment operator overload
      SynchronisedPreconditionedFilteredSolverData& operator=(const SynchronisedPreconditionedFilteredSolverData& other)
      {
          if(this == &other)
            return *this;

        this->_stored_sys = other._stored_sys;
        this->_stored_localsys = other._stored_localsys;
        this->_stored_rhs = other._stored_rhs;
        this->_stored_sol = other._stored_sol;
        this->_stored_def = other._stored_def;
        this->_stored_temp = other._stored_temp;
        this->_stored_scalars = other._stored_scalars;
        this->_stored_indices = other._stored_indices;
        this->_stored_norm_0 = other._stored_norm_0;
        this->_stored_norm = other._stored_norm;
        this->_stored_eps = other._stored_eps;
        this->_stored_max_iters = Index(0);
        this->_stored_used_iters = Index(0);

        this->_stored_precon = other._stored_precon;

        this->_stored_vector_mirrors = other._stored_vector_mirrors;
        this->_stored_vector_mirror_sendbufs = other._stored_vector_mirror_sendbufs;
        this->_stored_vector_mirror_recvbufs = other._stored_vector_mirror_recvbufs;
        this->_stored_dest_ranks = other._stored_dest_ranks;
        this->_stored_source_ranks = other._stored_source_ranks;

        this->_stored_filter = other._stored_filter;

        return *this;
      }
    };
  }
}

#endif
