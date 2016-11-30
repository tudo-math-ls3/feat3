#pragma once
#ifndef GLOBAL_SYNCH_MAT_HPP
#define GLOBAL_SYNCH_MAT_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/lafem/matrix_mirror.hpp>
#include <kernel/lafem/power_diag_matrix.hpp>

#include <vector>
#include <array>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Ticket class for asynchronous global matrix conversion
     *
     * \todo statistics
     *
     * \tparam MT_
     * The matrix type
     *
     * \tparam VMT_
     * The vector mirror type
     *
     * \author Peter Zajac
     */
    template<typename MT_, typename VMT_>
    class SynchMatrix
    {
    public:
      /// the matrix mirror type
      typedef LAFEM::MatrixMirror<Mem::Main, typename MT_::DataType, typename MT_::IndexType> MatrixMirrorType;
      /// the buffer matrix type
      typedef LAFEM::MatrixMirrorBuffer<Mem::Main, typename MT_::DataType, typename MT_::IndexType> BufferMatrixType;

    protected:
      bool _initialised;
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      /// our communicator
      const Dist::Comm& _comm;
      /// the neighbour ranks
      std::vector<int> _ranks;
      /// the matrix mirrors
      std::vector<MatrixMirrorType> _mirrors;
      /// send and receive request vectors
      Dist::RequestVector _send_reqs, _recv_reqs;
      /// send and receive buffers
      std::vector<BufferMatrixType> _send_bufs, _recv_bufs;
#endif // FEAT_HAVE_MPI || DOXYGEN

    public:
      /**
       * \brief Constructor
       *
       * \param[in] comm
       * The communicator
       *
       * \param[in] ranks
       * The neighbour ranks within the communicator
       *
       * \param[in] mirrors_row
       * The row vector mirrors to be used for synchronisation
       *
       * \param[in] mirrors_col
       * The column vector mirrors to be used for synchronisation
       */
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      SynchMatrix(const Dist::Comm& comm, const std::vector<int>& ranks,
        const std::vector<VMT_>& mirrors_row, const std::vector<VMT_>& mirrors_col) :
        _initialised(false),
        _comm(comm),
        _ranks(ranks),
        _send_reqs(ranks.size()),
        _recv_reqs(ranks.size()),
        _send_bufs(ranks.size()),
        _recv_bufs(ranks.size())
      {
        const std::size_t n = ranks.size();

        XASSERTM(mirrors_row.size() == n, "invalid row vector mirror count");
        XASSERTM(mirrors_col.size() == n, "invalid column vector mirror count");

        _mirrors.reserve(n);

        // create matrix mirrors and buffers
        for(std::size_t i(0); i < n; ++i)
        {
          const VMT_& mir_r = mirrors_row.at(i);
          const VMT_& mir_c = mirrors_col.at(i);

          // create matrix mirror
          _mirrors.push_back(MatrixMirrorType(mir_r, mir_c));
        }
      }
#else // non-MPI version
      SynchMatrix(const Dist::Comm&, const std::vector<int>& ranks, const std::vector<VMT_>&, const std::vector<VMT_>&) :
        _initialised(false)
      {
        XASSERT(ranks.empty());
      }
#endif // FEAT_HAVE_MPI

      /// Unwanted copy constructor: Do not implement!
      SynchMatrix(const SynchMatrix &) = delete;
      /// Unwanted copy assignment operator: Do not implement!
      SynchMatrix & operator=(const SynchMatrix &) = delete;

      /**
       * \brief Initialises the internal buffers for synchronisation
       *
       * \param[in] matrix
       * The matrix to be used as a template for the buffers. The structure must be initialised,
       * but the numerical content of the matrix is ignored.
       */
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      void init(const MT_& matrix)
      {
        XASSERTM(!_initialised, "SynchMatrix object is already initialised");

        const std::size_t n = _ranks.size();

        // create matrix mirror buffers
        for(std::size_t i(0); i < n; ++i)
        {
          _send_bufs.at(i) = _mirrors.at(i).create_buffer(matrix);
        }

        // receive buffer dimensions vector
        std::vector<std::array<Index,4>> recv_dims(n);

        // post send-buffer dimension receives
        for(std::size_t i(0); i < n; ++i)
        {
          _recv_reqs[i] = _comm.irecv(recv_dims.at(i).data(), std::size_t(4), _ranks.at(i));
        }

        // send send-buffer dimensions
        for(std::size_t i(0); i < n; ++i)
        {
          const BufferMatrixType& sbuf = _send_bufs.at(i);
          Index dims[4] =
          {
            sbuf.rows(),
            sbuf.columns(),
            sbuf.entries_per_nonzero(),
            sbuf.used_elements()
          };
          _send_reqs[i] = _comm.isend(dims, std::size_t(4), _ranks.at(i));
        }

        // wait for all receives to finish
        _recv_reqs.wait_all();

        // create receive buffers and post receives
        for(std::size_t i(0); i < n; ++i)
        {
          // get the receive buffer dimensions
          Index nrows = recv_dims.at(i)[0];
          Index ncols = recv_dims.at(i)[1];
          Index nepnz = recv_dims.at(i)[2];
          Index nnze  = recv_dims.at(i)[3];

          // allocate receive buffer
          _recv_bufs.at(i) = BufferMatrixType(nrows, ncols, nnze, nepnz);
        }

        // post buffer row-pointer array receives
        for(std::size_t i(0); i < n; ++i)
        {
          _recv_reqs[i] = _comm.irecv(_recv_bufs.at(i).row_ptr(), _recv_bufs.at(i).rows()+std::size_t(1), _ranks.at(i));
        }

        // wait for all previous sends to finish
        _send_reqs.wait_all();

        // post buffer row-pointer array sends
        for(std::size_t i(0); i < n; ++i)
        {
          _send_reqs[i] = _comm.isend(_send_bufs.at(i).row_ptr(), _send_bufs.at(i).rows()+std::size_t(1), _ranks.at(i));
        }

        // wait for all previous receives to finish
        _recv_reqs.wait_all();

        // post buffer row-pointer array receives
        for(std::size_t i(0); i < n; ++i)
        {
          _recv_reqs[i] = _comm.irecv(_recv_bufs.at(i).col_ind(), _recv_bufs.at(i).used_elements(), _ranks.at(i));
        }

        // wait for all previous sends to finish
        _send_reqs.wait_all();

        // post buffer row-pointer array sends
        for(std::size_t i(0); i < n; ++i)
        {
          _send_reqs[i] = _comm.isend(_send_bufs.at(i).col_ind(), _send_bufs.at(i).used_elements(), _ranks.at(i));
        }

        // wait for all receives and sends to finish
        _recv_reqs.wait_all();
        _send_reqs.wait_all();

        _initialised = true;
      }
#else // non-MPI version
      void init(const MT_&)
      {
        XASSERTM(!_initialised, "SynchMatrix object is already initialised");
        _initialised = true;
      }
#endif // FEAT_HAVE_MPI

      /**
       * \brief Converts a type-0 matrix to a type-1 matrix.
       *
       * \param[inout] matrix
       * The type-0 matrix that is to be converted to type-1.
       */
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      void exec(MT_& matrix)
      {
        XASSERTM(_initialised, "SynchMatrix object has not been initialised");

        const std::size_t n = _ranks.size();

        // create receive buffers and post receives
        for(std::size_t i(0); i < n; ++i)
        {
          BufferMatrixType& buf = _recv_bufs.at(i);

          // post receive
          _recv_reqs.get_request(i) = _comm.irecv(buf.val(), buf.val_size(), _ranks.at(i));
        }

        // post sends
        for(std::size_t i(0); i < n; ++i)
        {
          BufferMatrixType& buf = _send_bufs.at(i);

          // gather from mirror
          _mirrors.at(i).gather(buf, matrix);

          // post send
          _send_reqs.get_request(i) = _comm.isend(buf.val(), buf.val_size(), _ranks.at(i));
        }

        // process all pending receives
        for(std::size_t idx; _recv_reqs.wait_any(idx); )
        {
          // scatter the receive buffer
          _mirrors.at(idx).scatter_axpy(matrix, _recv_bufs.at(idx));
        }

        // wait for all sends to finish
        _send_reqs.wait_all();
      }
#else // non-MPI version
      void exec(MT_&)
      {
        XASSERTM(_initialised, "SynchMatrix object has not been initialised");
      }
#endif // FEAT_HAVE_MPI
    }; // class SynchMatrix

    /**
     * \brief Synchronises a type-0 matrix
     *
     * \param[inout] target
     * The type-0 matrix to be synchronised
     *
     * \param[in] comm
     * The communicator
     *
     * \param[in] ranks
     * The neighbour ranks within the communicator
     *
     * \param[in] mirrors_row
     * The row vector mirrors to be used for synchronisation
     *
     * \param[in] mirrors_col
     * The column vector mirrors to be used for synchronisation
     */
    template<typename MT_, typename VMT_>
    void synch_matrix(MT_& target, const Dist::Comm& comm, const std::vector<int>& ranks,
        const std::vector<VMT_>& mirrors_row, const std::vector<VMT_>& mirrors_col)
    {
      SynchMatrix<MT_, VMT_> synch(comm, ranks, mirrors_row, mirrors_col);
      synch.init(target);
      synch.exec(target);
    }

    template <typename MT_, typename SVMT_, int blocks_>
    class SynchMatrix<LAFEM::PowerDiagMatrix<MT_, blocks_>, SVMT_>
    {
    public:
      using VMT_ = typename SVMT_::SubMirrorType;
      using SMT_ = LAFEM::PowerDiagMatrix<MT_, blocks_>;
      bool _initialised;

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      std::vector<std::shared_ptr<SynchMatrix<MT_, VMT_>>> synch_matrix_list;
      std::vector<std::vector<VMT_>> mirrors_row_split;
      std::vector<std::vector<VMT_>> mirrors_col_split;
#endif // FEAT_HAVE_MPI || DOXYGEN

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      SynchMatrix(const Dist::Comm& comm, const std::vector<int>& ranks,
        const std::vector<SVMT_>& mirrors_row, const std::vector<SVMT_>& mirrors_col) :
        _initialised(false),
        synch_matrix_list(blocks_),
        mirrors_row_split(blocks_),
        mirrors_col_split(blocks_)
      {
        for (int block(0) ; block < blocks_ ; ++block)
        {
          for (Index i(0) ; i < mirrors_row.size() ; ++i)
          {
            mirrors_row_split.at((size_t)block).push_back(mirrors_row.at(i).get(block).clone(LAFEM::CloneMode::Shallow));
          }
          for (Index i(0) ; i < mirrors_col.size() ; ++i)
          {
            mirrors_col_split.at((size_t)block).push_back(mirrors_col.at(i).get(block).clone(LAFEM::CloneMode::Shallow));
          }

        synch_matrix_list.at((size_t)block) = std::make_shared<SynchMatrix<MT_, VMT_>>(comm, ranks, mirrors_row_split.at((size_t)block), mirrors_col_split.at((size_t)block));

        }
      }
#else // non-MPI version
      SynchMatrix(const Dist::Comm&, const std::vector<int>& ranks, const std::vector<SVMT_>&, const std::vector<SVMT_>&) :
        _initialised(false)
      {
        XASSERT(ranks.empty());
      }
#endif // FEAT_HAVE_MPI

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      void init(const SMT_& matrix)
      {
        XASSERTM(!_initialised, "SynchMatrix object is already initialised");

        for (int block(0) ; block < blocks_ ; ++block)
        {
          synch_matrix_list.at((size_t)block)->init(matrix.get(block, block));
        }

        _initialised = true;
      }
#else // non-MPI version
      void init(const MT_&)
      {
        XASSERTM(!_initialised, "SynchMatrix object is already initialised");
        _initialised = true;
      }
#endif // FEAT_HAVE_MPI

#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)
      void exec(SMT_& matrix)
      {
        XASSERTM(_initialised, "SynchMatrix object has not been initialised");

        for (int block(0) ; block < blocks_ ; ++block)
        {
          synch_matrix_list.at((size_t)block)->exec(matrix.get(block, block));
        }
      }
#else // non-MPI version
      void exec(SMT_&)
      {
        XASSERTM(_initialised, "SynchMatrix object has not been initialised");
      }
#endif // FEAT_HAVE_MPI
    };

  } // namespace Global
} // namespace FEAT

#endif // GLOBAL_SYNCH_MAT_HPP
