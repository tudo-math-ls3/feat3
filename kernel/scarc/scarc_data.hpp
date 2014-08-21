/**
 * \file
 * \brief FEAST milestone 2 ScaRC data implementations
 * \author Markus Geveler
 * \date 2014
 *
 * See class documentation.
 *
 */

#pragma once
#ifndef SCARC_GUARD_SCARC_DATA_HPP
#define SCARC_GUARD_SCARC_DATA_HPP 1

#include<kernel/lafem/dense_vector.hpp>
#include<kernel/lafem/vector_mirror.hpp>
#include<kernel/lafem/sparse_matrix_csr.hpp>
#include<kernel/lafem/unit_filter.hpp>

#include<kernel/foundation/communication.hpp>
#include<kernel/foundation/gateway.hpp>

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
             typename VectorType_ = DenseVector<MemTag_, DataType_>,
             typename MatrixType_ = SparseMatrixCSR<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index>
    struct ScaRCDataBase
    {
      public:
        ///type exports
        typedef MemTag_ mem_;
        typedef VectorType_ vector_type_;
        typedef MatrixType_ matrix_type_;
        typedef StorageType_<VectorType_, std::allocator<VectorType_ > > vector_storage_type_;
        typedef StorageType_<MatrixType_, std::allocator<MatrixType_ > > matrix_storage_type_;
        typedef StorageType_<DataType_, std::allocator<DataType_> > scalar_storage_type_;
        typedef StorageType_<IT_, std::allocator<IT_> > index_storage_type_;
        typedef IT_ index_type_;

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

        virtual ~ScaRCDataBase()
        {
        }

      protected:
        ///MemTag_ memory to store matrices and vectors and scalars
        matrix_type_ _stored_sys;
        matrix_type_ _stored_localsys;
        vector_type_ _stored_rhs;
        vector_type_ _stored_sol;
        vector_type_ _stored_def;
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
             typename VectorType_ = DenseVector<MemTag_, DataType_>,
             typename MatrixType_ = SparseMatrixCSR<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index>
    struct ScaRCData : public ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>
    {
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::mem_ mem_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::matrix_type_ matrix_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::vector_type_ vector_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::vector_storage_type_ vector_storage_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::scalar_storage_type_ scalar_storage_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::index_storage_type_ index_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "ScaRCData";
      }

      ///CTOR from system data, moves data to solver data
      ScaRCData(matrix_type_&& A,
                 vector_type_&& x,
                 vector_type_&& b)
      {
        this->_stored_sys = std::move(A);
        this->_stored_rhs = std::move(b);
        this->_stored_sol = std::move(x);
        this->_stored_def = vector_type_(this->_stored_rhs.size());
      }

      ///move CTOR
      ScaRCData(ScaRCData&& other)
      {
        this->_stored_sys = std::move(other._stored_sys);
        this->_stored_localsys = std::move(other._stored_localsys);
        this->_stored_rhs = std::move(other._stored_rhs);
        this->_stored_sol = std::move(other._stored_sol);
        this->_stored_def = std::move(other._stored_def);
      }

      ///move assignment operator overload
      ScaRCData& operator=(ScaRCData&& other)
      {
          if(this == &other)
            return *this;

        this->_stored_sys = std::move(other._stored_sys);
        this->_stored_localsys = std::move(other._stored_localsys);
        this->_stored_rhs = std::move(other._stored_rhs);
        this->_stored_sol = std::move(other._stored_sol);
        this->_stored_def = std::move(other._stored_def);

        return *this;
      }

    };

    /**
     * \brief Preconditioner data interface
     *
     * Can be subclassed in addition to ScaRCDataBase (or subclasses thereof) in order to store
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
             typename PreconContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             typename IT_ = Index>
    struct PreconditionerDataContainer
    {
      public:
        typedef PreconContType_ precon_type_;

        virtual precon_type_& precon()
        {
          return _stored_precon;
        }

        virtual const precon_type_& precon() const
        {
          return _stored_precon;
        }

        /// virtual destructor
        virtual ~PreconditionerDataContainer()
        {
        }

      protected:
        ///CTORs to be used in subclasses, moves data into interface
        PreconditionerDataContainer(precon_type_&& p) :
          _stored_precon(std::move(p))
        {
        }

        PreconditionerDataContainer(PreconditionerDataContainer&& other)
        {
          this->_stored_precon = std::move(other._stored_precon);
        }

        precon_type_ _stored_precon;
    };

    /**
     * \brief Simple data container with preconditioner data
     *
     * PreconditionedScaRCData combines SoverData and PreconditionerData.
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
             typename VectorType_ = DenseVector<MemTag_, DataType_>,
             typename MatrixType_ = SparseMatrixCSR<MemTag_, DataType_>,
             typename PreconContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index>
    struct PreconditionedScaRCData :
      public ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>,
      public PreconditionerDataContainer<DataType_, MemTag_, PreconContType_, IT_>
    {
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::mem_ mem_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::matrix_type_ matrix_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::vector_type_ vector_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::vector_storage_type_ vector_storage_type_;
      typedef PreconContType_ precon_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "PreconditionedScaRCData";
      }

      ///CTOR from system data
      PreconditionedScaRCData(matrix_type_&& A,
                               precon_type_&& P,
                               vector_type_&& x,
                               vector_type_&& b) :
        ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(std::move(A), std::move(x), std::move(b)),
        PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>(std::move(P))
      {
      }

      ///move CTOR
      PreconditionedScaRCData(PreconditionedScaRCData&& other) :
        ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(other),
        PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>(other)
      {
      }

      ///move assignment operator overload
      PreconditionedScaRCData& operator=(PreconditionedScaRCData&& other)
      {
          if(this == &other)
            return *this;

        this->_stored_sys = std::move(other._stored_sys);
        this->_stored_localsys = std::move(other._stored_localsys);
        this->_stored_rhs = std::move(other._stored_rhs);
        this->_stored_sol = std::move(other._stored_sol);
        this->_stored_def = std::move(other._stored_def);

        this->_stored_precon = std::move(other._stored_precon);

        return *this;
      }
    };

    /**
     * \brief Synchronisation data interface
     *
     * Can be subclassed in addition to ScaRCDataBase (or subclasses thereof) in order to store
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
             typename VectorType_ = DenseVector<MemTag_, DataType_>,
             typename VectorMirrorType_ = VectorMirror<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index>
    struct SynchronizationDataContainer
    {
      public:
        typedef VectorType_ vector_type_;
        typedef VectorMirrorType_ vector_mirror_type_;
        typedef StorageType_<VectorType_, std::allocator<VectorType_ > > vector_storage_type_;
        typedef StorageType_<VectorMirrorType_, std::allocator<VectorMirrorType_ > > vector_mirror_storage_type_;
        typedef StorageType_<IT_, std::allocator<IT_> > index_storage_type_;

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

        virtual vector_type_& halo_frequencies()
        {
          return _halo_frequencies;
        }

        virtual const vector_type_& halo_frequencies() const
        {
          return _halo_frequencies;
        }

        virtual index_storage_type_& tags()
        {
          return _tags;
        }

        virtual const index_storage_type_& tags() const
        {
          return _tags;
        }


        virtual StorageType_<Communicator, std::allocator<Communicator> >& communicators()
        {
          return _communicators;
        }

        virtual const StorageType_<Communicator, std::allocator<Communicator> >& communicators() const
        {
          return _communicators;
        }


        virtual ~SynchronizationDataContainer()
        {
        }

      protected:
        ///CTORs to be used in subclasses
        SynchronizationDataContainer() :
          _stored_vector_mirrors(),
          _stored_vector_mirror_sendbufs(),
          _stored_vector_mirror_recvbufs(),
          _stored_dest_ranks(),
          _stored_source_ranks(),
          _tags(),
          _communicators()
        {
        }

        SynchronizationDataContainer(SynchronizationDataContainer&& other)
        {
          this->_stored_vector_mirrors = std::move(other._stored_vector_mirrors);
          this->_stored_vector_mirror_sendbufs = std::move(other._stored_vector_mirror_sendbufs);
          this->_stored_vector_mirror_recvbufs = std::move(other._stored_vector_mirror_recvbufs);
          this->_stored_dest_ranks = std::move(other._stored_dest_ranks);
          this->_stored_source_ranks = std::move(other._stored_source_ranks);
          this->_halo_frequencies = std::move(other._halo_frequencies);
          this->_tags = std::move(other._tags);
          this->_communicators = std::move(other._communicators);
        }

        vector_mirror_storage_type_ _stored_vector_mirrors;
        vector_storage_type_ _stored_vector_mirror_sendbufs;
        vector_storage_type_ _stored_vector_mirror_recvbufs;
        index_storage_type_ _stored_dest_ranks;
        index_storage_type_ _stored_source_ranks;

        vector_type_ _halo_frequencies;

        index_storage_type_ _tags;

        ///communicators associated with scarc layers of the same index
        StorageType_<Communicator, std::allocator<Communicator> > _communicators;
    };

    /**
     * \brief Simple data container with synchronisation data
     *
     * SynchronisedScaRCData combines SoverData and SynchronizationDataContainer.
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
             typename VectorType_ = DenseVector<MemTag_, DataType_>,
             typename VectorMirrorType_ = VectorMirror<MemTag_, DataType_>,
             typename MatrixType_ = SparseMatrixCSR<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index>
    struct SynchronisedScaRCData :
      public ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>,
      public SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_, IT_>
    {
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::mem_ mem_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::matrix_type_ matrix_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::vector_type_ vector_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::vector_storage_type_ vector_storage_type_;
      typedef StorageType_<VectorMirrorType_, std::allocator<VectorMirrorType_> > vector_mirror_storage_type_;
      typedef StorageType_<IT_, std::allocator<IT_> > index_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "SynchronisedScaRCData";
      }

      ///CTOR from system data
      SynchronisedScaRCData(matrix_type_&& A,
                             vector_type_&& x,
                             vector_type_&& b) :
        ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(std::move(A), std::move(x), std::move(b)),
        SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>()
      {
      }

      ///move CTOR
      SynchronisedScaRCData(SynchronisedScaRCData&& other) :
        ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(other),
        SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>(other)
      {
      }

      ///move assignment operator overload
      SynchronisedScaRCData& operator=(SynchronisedScaRCData&& other)
      {
          if(this == &other)
            return *this;

        this->_stored_sys = std::move(other._stored_sys);
        this->_stored_localsys = std::move(other._stored_localsys);
        this->_stored_rhs = std::move(other._stored_rhs);
        this->_stored_sol = std::move(other._stored_sol);
        this->_stored_def = std::move(other._stored_def);

        this->_stored_vector_mirrors = std::move(other._stored_vector_mirrors);
        this->_stored_vector_mirror_sendbufs = std::move(other._stored_vector_mirror_sendbufs);
        this->_stored_vector_mirror_recvbufs = std::move(other._stored_vector_mirror_recvbufs);
        this->_stored_dest_ranks = std::move(other._stored_dest_ranks);
        this->_stored_source_ranks = std::move(other._stored_source_ranks);
        this->_halo_frequencies = std::move(other._halo_frequencies);

        this->_tags = std::move(other._tags);
        this->_communicators = std::move(other._communicators);

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
             typename VectorType_ = DenseVector<MemTag_, DataType_>,
             typename VectorMirrorType_ = VectorMirror<MemTag_, DataType_>,
             typename MatrixType_ = SparseMatrixCSR<MemTag_, DataType_>,
             typename PreconContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index>
    struct SynchronisedPreconditionedScaRCData :
      public ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>,
      public PreconditionerDataContainer<DataType_, MemTag_, PreconContType_, IT_>,
      public SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_, IT_>
    {
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::mem_ mem_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::matrix_type_ matrix_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::vector_type_ vector_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::vector_storage_type_ vector_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "SynchronisedPreconditionedScaRCData";
      }

      ///CTOR from system data
      SynchronisedPreconditionedScaRCData(matrix_type_&& A,
                                           PreconContType_&& P,
                                           vector_type_&& x,
                                           vector_type_&& b) :
        ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(std::move(A), std::move(x), std::move(b)),
        PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>(std::move(P)),
        SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>()

      {
      }

      ///move CTOR
      SynchronisedPreconditionedScaRCData(SynchronisedPreconditionedScaRCData&& other) :
        ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(std::move(other)),
        PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>(std::move(other)),
        SynchronisedScaRCData<DataType_, MemTag_, VectorType_, VectorMirrorType_, MatrixType_, StorageType_>(std::move(other))
      {
      }

      ///assignment operator overload
      SynchronisedPreconditionedScaRCData& operator=(SynchronisedPreconditionedScaRCData&& other)
      {
          if(this == &other)
            return *this;

        this->_stored_sys = std::move(other._stored_sys);
        this->_stored_localsys = std::move(other._stored_localsys);
        this->_stored_rhs = std::move(other._stored_rhs);
        this->_stored_sol = std::move(other._stored_sol);
        this->_stored_def = std::move(other._stored_def);
        this->_stored_precon = std::move(other._stored_precon);
        this->_stored_vector_mirrors = std::move(other._stored_vector_mirrors);
        this->_stored_vector_mirror_sendbufs = std::move(other._stored_vector_mirror_sendbufs);
        this->_stored_vector_mirror_recvbufs = std::move(other._stored_vector_mirror_recvbufs);
        this->_stored_dest_ranks = std::move(other._stored_dest_ranks);
        this->_stored_source_ranks = std::move(other._stored_source_ranks);
        this->_halo_frequencies = std::move(other._halo_frequencies);

        this->_tags = std::move(other._tags);
        this->_communicators = std::move(other._communicators);
        return *this;
      }
    };

    /**
     * \brief Filter data interface
     *
     * Can be subclassed in addition to ScaRCDataBase (subclasses thereof) in order to store
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
             typename FilterType_ = UnitFilter<MemTag_, DataType_>,
             typename IT_ = Index>
    struct FilterDataContainer
    {
      public:
        typedef FilterType_ filter_type_;

        virtual filter_type_& filter()
        {
          return _stored_filter;
        }

        virtual const filter_type_& filter() const
        {
          return _stored_filter;
        }

        virtual ~FilterDataContainer()
        {
        }

      protected:
        ///CTORs to be used in subclasses
        FilterDataContainer(filter_type_&& f) :
          _stored_filter(std::move(f))
        {
        }

        FilterDataContainer(FilterDataContainer&& other)
        {
          this->_stored_filter = std::move(other._stored_filter);
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
             typename VectorType_ = DenseVector<MemTag_, DataType_>,
             typename VectorMirrorType_ = VectorMirror<MemTag_, DataType_>,
             typename MatrixType_ = SparseMatrixCSR<MemTag_, DataType_>,
             typename PreconContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             typename FilterType_ = UnitFilter<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index>
    struct SynchronisedPreconditionedFilteredScaRCData :
      public ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>,
      public FilterDataContainer<DataType_, MemTag_, FilterType_, IT_>,
      public PreconditionerDataContainer<DataType_, MemTag_, PreconContType_, IT_>,
      public SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_, IT_>
    {
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::mem_ mem_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::matrix_type_ matrix_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::vector_type_ vector_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::vector_storage_type_ vector_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "SynchronisedPreconditionedFilteredScaRCData";
      }

      ///CTOR from system data
      SynchronisedPreconditionedFilteredScaRCData(matrix_type_&& A,
                                                   PreconContType_&& P,
                                                   vector_type_&& x,
                                                   vector_type_&& b,
                                                   FilterType_&& f) :
        ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(std::move(A), std::move(x), std::move(b)),
        FilterDataContainer<DataType_, MemTag_, FilterType_, IT_>(std::move(f)),
        PreconditionerDataContainer<DataType_, MemTag_, PreconContType_>(std::move(P)),
        SynchronizationDataContainer<DataType_, MemTag_, VectorType_, VectorMirrorType_, StorageType_>()
      {
      }

      ///move CTOR
      SynchronisedPreconditionedFilteredScaRCData(SynchronisedPreconditionedFilteredScaRCData&& other) :
        ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>(std::move(other)),
        FilterDataContainer<DataType_, MemTag_, FilterType_, IT_>(std::move(other)),
        PreconditionerDataContainer<DataType_, MemTag_, PreconContType_, IT_>(std::move(other)),
        SynchronisedScaRCData<DataType_, MemTag_, VectorType_, VectorMirrorType_, MatrixType_, StorageType_, IT_>(std::move(other))
      {
      }

      ///move assignment operator overload
      SynchronisedPreconditionedFilteredScaRCData& operator=(SynchronisedPreconditionedFilteredScaRCData&& other)
      {
          if(this == &other)
            return *this;

        this->_stored_sys = std::move(other._stored_sys);
        this->_stored_localsys = std::move(other._stored_localsys);
        this->_stored_rhs = std::move(other._stored_rhs);
        this->_stored_sol = std::move(other._stored_sol);
        this->_stored_def = std::move(other._stored_def);
        this->_stored_precon = std::move(other._stored_precon);
        this->_stored_vector_mirrors = std::move(other._stored_vector_mirrors);
        this->_stored_vector_mirror_sendbufs = std::move(other._stored_vector_mirror_sendbufs);
        this->_stored_vector_mirror_recvbufs = std::move(other._stored_vector_mirror_recvbufs);
        this->_stored_dest_ranks = std::move(other._stored_dest_ranks);
        this->_stored_source_ranks = std::move(other._stored_source_ranks);
        this->_halo_frequencies = std::move(other._halo_frequencies);

        this->_tags = std::move(other._tags);
        this->_communicators = std::move(other._communicators);
        this->_stored_filter = std::move(other._stored_filter);

        return *this;
      }
    };

    template<typename DataType_ = double,
             typename MemTag_ = Mem::Main,
             typename VectorType_ = DenseVector<MemTag_, DataType_>,
             typename MatrixType_ = SparseMatrixCSR<MemTag_, DataType_>,
             typename TransferContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             typename PreconContType_ = SparseMatrixCSR<MemTag_, DataType_>,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index>
    struct MultiLevelScaRCData : public ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>
    {
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::mem_ mem_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::matrix_type_ matrix_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::vector_type_ vector_type_;
      typedef typename ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_>::vector_storage_type_ vector_storage_type_;

      typedef StorageType_<std::shared_ptr<ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_> >, std::allocator<std::shared_ptr<ScaRCDataBase<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_, IT_> > > > leveldata_storage_type_;

      ///fulfill pure virtual
      virtual const std::string type_name()
      {
        return "MultiLevelScaRCData";
      }

      ///CTOR from system
      MultiLevelScaRCData(matrix_type_&& A,
                           vector_type_&& x,
                           vector_type_&& b) :
        ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(std::move(A), std::move(x), std::move(b)),
        stored_level_data()
      {
        ///TODO
      }

      ///move CTOR
      MultiLevelScaRCData(MultiLevelScaRCData&& other) :
        ScaRCData<DataType_, MemTag_, VectorType_, MatrixType_, StorageType_>(std::move(other)),
        stored_level_data(std::move(other.stored_level_data))
      {
      }

      ///move assignment operator overload
      MultiLevelScaRCData& operator=(MultiLevelScaRCData&& other)
      {
          if(this == &other)
            return *this;

        this->_stored_sys = std::move(other._stored_sys);
        this->_stored_localsys = std::move(other._stored_localsys);
        this->_stored_rhs = std::move(other._stored_rhs);
        this->_stored_sol = std::move(other._stored_sol);
        this->_stored_def = std::move(other._stored_def);
        this->_stored_temp = std::move(other._stored_temp);
        this->_stored_scalars = std::move(other._stored_scalars);
        this->_stored_indices = std::move(other._stored_indices);
        this->_stored_norm_0 = other._stored_norm_0;
        this->_stored_norm = other._stored_norm;
        this->_stored_eps = other._stored_eps;
        this->_stored_max_iters = IT_(0);
        this->_stored_used_iters = IT_(0);

        this->stored_level_data = std::move(other._stored_level_data);

        return *this;
      }

      leveldata_storage_type_ stored_level_data;
    };

  }
}

#endif
