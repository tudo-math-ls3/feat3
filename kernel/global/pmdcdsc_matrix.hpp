#pragma once
#ifndef KERNEL_GLOBAL_PMDCDSC_MATRIX_HPP
#define KERNEL_GLOBAL_PMDCDSC_MATRIX_HPP 1

// includes, FEAT
#include <kernel/global/vector.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/util/stop_watch.hpp>

#include <vector>
#include <map>
#include <set>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Pre-Multiplied Discontinuous Diagonal Schur-Complement Matrix
     */
    template<typename MatrixB_, typename MatrixD_>
    class PMDCDSCMatrix;

    /**
     * \brief Pre-Multiplied Discontinuous Diagonal Schur-Complement Matrix
     *
     * This class implements an efficient pre-multiplied distributed Schur-Complement matrix
     *
     * \f[S := D\cdot A\cdot B\f]
     *
     * where
     * - \e A is a diagonal matrix (usually an inverted lumped velocity mass matrix)
     * - \e D is a 1-by-d BCSR matrix (usually the velocity-divergence matrix of a Stokes system)
     * - \e B is a d-by-1 BCSR matrix (usually the pressure-gradient matrix of a Stokes system)
     *
     * \note
     * This class effectively performs the same task as the SymmetricLumpedSchurMatrix, however,
     * the internal implementation is different and relies on explicitly pre-multiplying the
     * Schur-complement matrices rather than performing three consecutive matrix-vector products.
     *
     * \attention
     * This class can only be used if the pressure space is discontinuous! Any attempts to use
     * this class for matrices deriving from a continuous pressure space will be aborted by a
     * failed assertion. Please note that this limitation is not just for kicks, but the design
     * of this class and its underlying algorithms simply <b>can not</b> work in this case.
     *
     * \author Peter Zajac
     */
    template<typename DT_, typename IT_, int dim_, typename MirrorV_, typename MirrorP_>
    class PMDCDSCMatrix<
      Global::Matrix<LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, dim_, 1>, MirrorV_, MirrorP_>,
      Global::Matrix<LAFEM::SparseMatrixBCSR<Mem::Main, DT_, IT_, 1, dim_>, MirrorV_, MirrorP_>>
    {
    public:
      typedef Mem::Main MemType;
      typedef DT_ DataType;
      typedef IT_ IndexType;
      static constexpr int dim = dim_;

      typedef MirrorV_ MirrorTypeV;
      typedef MirrorP_ MirrorTypeP;

      typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, dim, 1> LocalMatrixTypeB;
      typedef LAFEM::SparseMatrixBCSR<MemType, DataType, IndexType, 1, dim> LocalMatrixTypeD;
      typedef LAFEM::SparseMatrixCSR<MemType, DataType, IndexType> LocalMatrixTypeS;

      typedef Global::Matrix<LocalMatrixTypeB, MirrorTypeV, MirrorTypeP> GlobalMatrixTypeB;
      typedef Global::Matrix<LocalMatrixTypeD, MirrorTypeP, MirrorTypeV> GlobalMatrixTypeD;

      typedef LAFEM::DenseVectorBlocked<MemType, DataType, IndexType, dim> LocalVectorTypeV;
      typedef LAFEM::DenseVector<MemType, DataType, IndexType> LocalVectorTypeP;
      typedef LAFEM::DenseVector<MemType, DataType, IndexType> BufferVectorType;

      typedef Global::Vector<LocalVectorTypeV, MirrorTypeV> GlobalVectorTypeV;
      typedef Global::Vector<LocalVectorTypeP, MirrorTypeP> GlobalVectorTypeP;

      typedef Global::Gate<LocalVectorTypeV, MirrorTypeV> GateTypeV;
      typedef Global::Gate<LocalVectorTypeP, MirrorTypeP> GateTypeP;

      typedef GlobalVectorTypeP VectorTypeL;
      typedef GlobalVectorTypeP VectorTypeR;

      typedef typename LocalVectorTypeV::ValueType ValueTypeA;
      typedef typename LocalMatrixTypeB::ValueType ValueTypeB;
      typedef typename LocalMatrixTypeD::ValueType ValueTypeD;

      /// The row-gate type (used by SFINAE)
      typedef typename GlobalMatrixTypeD::GateRowType GateRowType;
      /// The column-gate type (used by SFINAE)
      typedef typename GlobalMatrixTypeB::GateColType GateColType;

      static constexpr bool is_global = true;
      static constexpr bool is_local = false;

    protected:
      /// the diagonal matrix a (stored as a vector)
      const GlobalVectorTypeV& _diagonal_a;
      /// the pressure-gradient matrix B
      const GlobalMatrixTypeB& _matrix_b;
      /// the velocity-divergence matrix D
      const GlobalMatrixTypeD& _matrix_d;

      /**
       * \brief The local Schur-complement matrix
       *
       * This matrix stores the product (D*A*B) restricted onto the patch of this process,
       * i.e. it is simply the product of those three local matrices.
       */
      LocalMatrixTypeS _matrix_s;

      /**
       * \brief Pre-Multiplied transposed matrix-product (D*A)^T
       *
       * This matrix stores the product (D*A) of the two given input matrices D and A in
       * transposed form, which corresponds to the CSC format rather than the CSR format.
       * As D is assumed to have the transposed layout of B, the layout of this transposed
       * matrix therefore coincides with the layout of B, which also explains why the
       * type of the matrix is LocalMatrixTypeB rather than LocalMatrixTypeD.
       *
       * This pre-multiplied matrix is required for the \e efficient numerical assembly of
       * the Neighbour Schur-Matrices, which is performed by the #_asm_neighbour_schur_matrix()
       * function. Unfortunately, it does not seem to be possible to implement that assembly
       * without the transpose of D in an efficient way. A previous attempt to implement that
       * numerical assembly without this explicit transpose has lead to a slow-down of factor
       * 20000 (in words: twenty thousand).
       *
       * \remark Although D^T is \e usually equal to B, we do not want to exploit this here,
       * as we cannot say if the user has modified the D/B matrices in any way, which would
       * break this relation and thus this class would not work correctly.\n
       * Just in case that you cannot imagine any case where D != B^T, here you go:
       * - The B-matrix may be filtered due to boundary conditions.
       * - The B-matrix may contain additional terms due to a Newton method applied onto
       *   non-newtonian Navier-Stokes equations with pressure-dependent viscosity.
       * - The B-matrix may contain additional stuff like stabilisation or penalty terms.
       */
      LocalMatrixTypeB _matrix_da;

      /// neighbour process ranks
      std::vector<int> _ranks;

      /**
       * \brief Pressure DOF mirrors
       *
       * For each neighbour process, this mirror stores the indices of all pressure DOFs
       * that have to be send to that particular neighbour process, so that it can perform
       * its matrix-vector product with its corresponding neighbour Schur-complement matrix.
       * This mirror contains all pressure DOFs which are defined on elements which are
       * directly adjacent to the communication halo mesh-part of the neighbour.
       */
      std::vector<MirrorTypeP> _pres_mirrors;

      /**
       * \brief B-matrix data mirrors
       *
       * For each neighbour process, this mirror stores the indices of all B-matrix entries
       * that have to be send to that particular neighbour process, so that it can pre-compute
       * its neighbour Schur-complement matrix during numeric initialisation.
       */
      std::vector<MirrorTypeP> _data_mirrors;

      /**
       * \brief Neighbour B-matrix adjacency graphs
       *
       * For each neighbour process, this mirror stores the adjacency graph of the neighbour's
       * B-matrix restricted to the pressure DOFs that this particular neighbour send to us.
       */
      std::vector<Adjacency::Graph> _neighbour_graphs;

      /**
       * \brief Pre-multiplied neighbour Schur-complement matrices
       *
       * For each neighbour process, this mirror stores the pre-multiplied Schur-complement
       * matrix (D*A*B_k), where B_k is the part of the neighbour's B-matrix, whose rows
       * correlate to the velocity DOFs shared by this process and the neighbour process.
       */
      std::vector<LocalMatrixTypeS> _neighbour_matrices;

      /// total symbolic/numeric init time
      StopWatch watch_init_symbolic;
      /// local Schur-Matrix symbolic assembly time
      StopWatch watch_init_sym_matrix_loc;
      /// pressure mirror assembly time
      StopWatch watch_init_sym_pres_mirror;
      /// reduced B-matrix symbolic assembly time
      StopWatch watch_init_sym_reduced_b;
      /// reduced B-matrix data mirror assembly time
      StopWatch watch_init_sym_data_mirror;
      /// neighbour Schur-Matrix symbolic assembly time
      StopWatch watch_init_sym_neighbour_s;

      /// total numeric init time
      StopWatch watch_init_numeric;
      /// local Schur-Matrix numeric assembly time
      StopWatch watch_init_num_matrix_loc;
      /// reduced-B matrix data gather time
      StopWatch watch_init_num_gather_b;
      /// pre-multiply D*A time
      StopWatch watch_init_num_premult_da;
      /// neighbour Schur-matrix numeric assembly time
      StopWatch watch_init_num_neighbour_s;

      /// total apply time
      mutable StopWatch watch_apply;
      /// local Schur-Matrix apply time
      mutable StopWatch watch_apply_matrix_loc;
      /// neighbour Matrix apply time
      mutable StopWatch watch_apply_neighbour_s;

    public:
      /**
       * \brief Constructor
       *
       * \note This constructor just saves the references, but does \b not perform any further initialisation.
       * This has to be done by hand by calling the #init_symbolic() and #init_numeric() functions.
       *
       * \param[in] diagonal_a
       * A reference to the vector representing the diagonal matrix A.
       *
       * \param[in] matrix_b
       * A reference to the pressure-gradient matrix B.
       *
       * \param[in] matrix_d
       * A reference to the velocity-divergence matrix D.
       */
      explicit PMDCDSCMatrix(
        const GlobalVectorTypeV& diagonal_a,
        const GlobalMatrixTypeB& matrix_b,
        const GlobalMatrixTypeD& matrix_d) :
        _diagonal_a(diagonal_a),
        _matrix_b(matrix_b),
        _matrix_d(matrix_d)
      {
      }

      /// virtual destructor
      virtual ~PMDCDSCMatrix()
      {
      }

      /// Returns a pointer to the underlying communicator.
      const Dist::Comm* get_comm() const
      {
        return _diagonal_a.get_comm();
      }

      /// \returns A new vector compatible L-vector.
      VectorTypeL create_vector_l() const
      {
        return _matrix_d.create_vector_l();
      }

      /// \returns A new vector compatible R-vector.
      VectorTypeR create_vector_r() const
      {
        return _matrix_b.create_vector_r();
      }

      /// \returns A reference to the local Schur-Complement matrix
      const LocalMatrixTypeS& get_local_schur_matrix() const
      {
        return _matrix_s;
      }

      /**
       * \brief Extracts the main diagonal
       *
       * \note
       * The main diagonal vector does not need to be synchronised, because the
       * underlying pressure space is discontinuous by definition/assumption.
       * However, nothing bad will happen if you synchronise the vector anyway.
       *
       * \param[out] vec_diag
       * A vector that receives the main diagonal of our matrix (D*A*B).
       */
      void extract_diag(GlobalVectorTypeP& vec_diag) const
      {
        _matrix_s.extract_diag(vec_diag.local());
      }

      /**
       * \brief Extracts the main diagonal
       *
       * \note
       * The main diagonal vector does not need to be synchronised, because the
       * underlying pressure space is discontinuous by definition/assumption.
       * However, nothing bad will happen if you synchronise the vector anyway.
       *
       * \returns
       * A vector that contains the main diagonal of our matrix (D*A*B).
       */
      GlobalVectorTypeP extract_diag() const
      {
        GlobalVectorTypeP vec_diag = _matrix_d.create_vector_l();
        _matrix_s.extract_diag(vec_diag.local());
        return vec_diag;
      }

      /// Resets the internal stop watches used for collecting timing statistics.
      void reset_timings()
      {
        watch_init_symbolic.reset();
        watch_init_sym_matrix_loc.reset();
        watch_init_sym_pres_mirror.reset();
        watch_init_sym_reduced_b.reset();
        watch_init_sym_data_mirror.reset();
        watch_init_sym_neighbour_s.reset();
        watch_init_numeric.reset();
        watch_init_num_matrix_loc.reset();
        watch_init_num_gather_b.reset();
        watch_init_num_premult_da.reset();
        watch_init_num_neighbour_s.reset();
        watch_apply.reset();
        watch_apply_matrix_loc.reset();
        watch_apply_neighbour_s.reset();
      }

      /**
       * \brief Returns a string that contains the formatted timing statistics.
       *
       * \attention This function must be called by all processes, because it contains a collective call.
       *
       * \returns
       * A string containing the timing statistics.
       */
      String format_timings() const
      {
        static constexpr std::size_t nt = 14;
        double tsum[nt], tmax[nt], tloc[nt] =
        {
          watch_init_symbolic.elapsed(),
          watch_init_sym_matrix_loc.elapsed(),
          watch_init_sym_pres_mirror.elapsed(),
          watch_init_sym_reduced_b.elapsed(),
          watch_init_sym_data_mirror.elapsed(),
          watch_init_sym_neighbour_s.elapsed(),
          watch_init_numeric.elapsed(),
          watch_init_num_matrix_loc.elapsed(),
          watch_init_num_gather_b.elapsed(),
          watch_init_num_premult_da.elapsed(),
          watch_init_num_neighbour_s.elapsed(),
          watch_apply.elapsed(),
          watch_apply_matrix_loc.elapsed(),
          watch_apply_neighbour_s.elapsed()
        };

        this->get_comm()->allreduce(tloc, tsum, nt, Dist::op_sum);
        this->get_comm()->allreduce(tloc, tmax, nt, Dist::op_max);

        // divide sum by number of ranks to obtain mean
        {
          const double ds = 1.0 / double(this->get_comm()->size());
          for(std::size_t i(0); i < nt; ++i)
            tsum[i] *= ds;
        }

        String s;
        s += String(34, ' ') + "Mean Time      Max Time\n";
        s += _fmt_time(tsum[0], tmax[0], "Total Symbolic Initialisation");
        s += _fmt_time(tsum[0], tmax[0], tsum[1], tmax[1], "Local Schur Matrix Structure");
        s += _fmt_time(tsum[0], tmax[0], tsum[2], tmax[2], "Pressure Mirror");
        s += _fmt_time(tsum[0], tmax[0], tsum[3], tmax[3], "Reduced-B Matrix Structure");
        s += _fmt_time(tsum[0], tmax[0], tsum[4], tmax[4], "Reduced-B Data Mirror");
        s += _fmt_time(tsum[0], tmax[0], tsum[5], tmax[5], "Neighbour Matrix Structure");
        double tsym_other_sum = tsum[0] - tsum[1] - tsum[2] - tsum[3] - tsum[4] - tsum[5];
        double tsym_other_max = tmax[0] - tmax[1] - tmax[2] - tmax[3] - tmax[4] - tmax[5];
        s += _fmt_time(tsum[0], tmax[0], tsym_other_sum, tsym_other_max, "Other Symbolic");

        s += _fmt_time(tsum[6], tmax[6], "Total Numeric Initialisation");
        s += _fmt_time(tsum[6], tmax[6], tsum[7], tmax[7], "Local Schur Matrix Values");
        s += _fmt_time(tsum[6], tmax[6], tsum[8], tmax[8], "Reduced-B Gather");
        s += _fmt_time(tsum[6], tmax[6], tsum[9], tmax[9], "Pre-Multiply D*A");
        s += _fmt_time(tsum[6], tmax[6], tsum[10], tmax[10], "Neighbour Matrix Values");
        double tnum_other_sum = tsum[6] - tsum[7] - tsum[8] - tsum[9] - tsum[10];
        double tnum_other_max = tmax[6] - tmax[7] - tmax[8] - tmax[9] - tmax[10];
        s += _fmt_time(tsum[6], tmax[6], tnum_other_sum, tnum_other_max, "Other Numeric");

        s += _fmt_time(tsum[11], tmax[11], "Total Matrix Apply");
        s += _fmt_time(tsum[11], tmax[11], tsum[12], tmax[12], "Local Schur Matrix");
        s += _fmt_time(tsum[11], tmax[11], tsum[13], tmax[13], "Neighbour Schur Matrix");
        double tapp_other_sum = tsum[11] - tsum[12] - tsum[13];
        double tapp_other_max = tmax[11] - tmax[12] - tmax[13];
        s += _fmt_time(tsum[11], tmax[11], tapp_other_sum, tapp_other_max, "Other Apply");
        return s;
      }

      /**
       * \brief Performs both symbolic and numeric initialisation of the matrix.
       */
      void init()
      {
        init_symbolic();
        init_numeric();
      }

      /**
       * \brief Performs the symbolic initialisation of the matrix.
       *
       * This function performs the following internal tasks:
       * -# Pre-compute the layout of the local Schur-complement matrix (D*A*B)
       * -# For each neighbour, compute a mirror of all pressure DOFs, which couple
       *    to at least one velocity DOF shared with that neighbour via our B-matrix.
       *    These are the pressure DOFs that have to be send to that neighbour during
       *    the matrix-vector multiplication in the #apply() function.
       * -# For each neighbour, compute the layout of the reduced B-matrix, which
       *    contains all B-matrix entries, whose column indices correspond to the
       *    pressure DOFs in the previously assembled pressure mirror. Translate the
       *    column indices from "local pressure DOF indices" to pressure mirror DOF
       *    indices during this step.
       * -# For each neighbour, compute the data-mirror of the reduced B-matrix.
       * -# Send the reduced B-matrix layouts to the corresponding neighbour processes.
       * -# For each neighbour, pre-compute the layout of the neighbour Schur-complement
       *    matrix S_k := (D*A*B_k), where B_k is the reduced B-matrix received from that
       *    particular neighbour process.
       */
      void init_symbolic()
      {
        watch_init_symbolic.start();

        // get the velocity and pressure gates
        const GateTypeV* gate_v = this->_matrix_b.get_row_gate();
        const GateTypeP* gate_p = this->_matrix_d.get_row_gate();
        XASSERT(gate_v != nullptr);
        XASSERT(gate_p != nullptr);

        // the pressure gate must be empty, otherwise the pressure space is not discontinuous
        XASSERTM(gate_p->_ranks.empty(), "pressure space is not discontinuous");

        // compute the local matrix structure of S by D * B
        {
          watch_init_sym_matrix_loc.start();

          // compose structures of D and B
          Adjacency::Graph graph_s(Adjacency::RenderType::injectify, _matrix_d.local(), _matrix_b.local());
          // sort column indices
          graph_s.sort_indices();
          // create the matrix layout of S
          _matrix_s = LocalMatrixTypeS(graph_s);

          watch_init_sym_matrix_loc.stop();
        }

        // get our communicator
        const Dist::Comm& comm = *gate_v->get_comm();

        // get neighbour ranks
        this->_ranks = gate_v->_ranks;

        // get the number of our neighbours
        const std::size_t num_neighs = this->_ranks.size();
        if(num_neighs <= std::size_t(0))
        {
          watch_init_symbolic.stop();
          return; // no neighbours, no problems :)
        }

        // copy the layout of B into (D*A)^T
        this->_matrix_da = this->_matrix_b.local().clone(LAFEM::CloneMode::Layout);

        // resize our member arrays
        this->_pres_mirrors.resize(num_neighs);
        this->_data_mirrors.resize(num_neighs);
        this->_neighbour_graphs.resize(num_neighs);
        this->_neighbour_matrices.resize(num_neighs);

        // allocate a vector of graphs for B
        std::vector<Adjacency::Graph> my_graphs(num_neighs);

        // loop over all neighbour processes
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          // get the velocity mirror
          const MirrorTypeV& mirror_v = gate_v->_mirrors.at(i);

          // assemble the pressure mirror for this neighbour
          watch_init_sym_pres_mirror.start();
          MirrorTypeP& mirror_p = this->_pres_mirrors.at(i);
          this->_asm_pres_mirror(mirror_p, mirror_v, this->_matrix_b.local());
          watch_init_sym_pres_mirror.stop();

          // assemble reduced B-matrix graph
          watch_init_sym_reduced_b.start();
          my_graphs.at(i) = this->_asm_reduced_b( mirror_p, mirror_v, this->_matrix_b.local());
          watch_init_sym_reduced_b.stop();

          // assemble B' data mirror
          watch_init_sym_data_mirror.start();
          this->_data_mirrors.at(i) = this->_asm_data_mirror(mirror_p, mirror_v, this->_matrix_b.local(), my_graphs.at(i));
          watch_init_sym_data_mirror.stop();
        }

        // dimension send/receive buffers and requests
        std::vector<std::array<Index,3>> recv_dims(num_neighs), send_dims(num_neighs);
        Dist::RequestVector recv_reqs(num_neighs), send_reqs(num_neighs);

        // post receive requests for dimensions
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          recv_reqs[i] = comm.irecv(recv_dims.at(i).data(), std::size_t(3), this->_ranks.at(i));
        }

        // send dimensions
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          const Adjacency::Graph& g = my_graphs.at(i);
          auto& sdim = send_dims.at(i);
          sdim[0] = g.get_num_nodes_domain(); // corresponds to velocity mirror index size
          sdim[1] = g.get_num_nodes_image();
          sdim[2] = g.get_num_indices();
          send_reqs[i] = comm.isend(sdim.data(), std::size_t(3), this->_ranks.at(i));
        }

        // process all pending receives
        for(std::size_t i; recv_reqs.wait_any(i); )
        {
          // get received dimension
          auto& rdim = recv_dims.at(i);

          // the first dimension must match our velocity mirror index set size
          XASSERT(rdim[0] == gate_v->_mirrors.at(i).num_indices());

          // allocate graph of corresponding dimensions
          this->_neighbour_graphs.at(i) = Adjacency::Graph(rdim[0], rdim[1], rdim[2]);
        }

        // post domain-pointer array receives
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          recv_reqs[i] = comm.irecv(this->_neighbour_graphs.at(i).get_domain_ptr(),
            recv_dims.at(i)[0] + std::size_t(1),  this->_ranks.at(i));
        }

        // wait for all previous sends to finish
        send_reqs.wait_all();

        // post domain-pointer array sends
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          send_reqs[i] = comm.isend(my_graphs.at(i).get_domain_ptr(),
            send_dims.at(i)[0] + std::size_t(1), this->_ranks.at(i));
        }

        // wait for all pending receives to finish
        recv_reqs.wait_all();

        // post image-index array receives
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          recv_reqs[i] = comm.irecv(this->_neighbour_graphs.at(i).get_image_idx(),
            recv_dims.at(i)[2],  this->_ranks.at(i));
        }

        // wait for all previous sends to finish
        send_reqs.wait_all();

        // post image-index array sends
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          send_reqs[i] = comm.isend(my_graphs.at(i).get_image_idx(), send_dims.at(i)[2], this->_ranks.at(i));
        }

        // wait for all pending receives to finish
        recv_reqs.wait_all();

        // wait for all previous sends to finish
        send_reqs.wait_all();

        // compute Schur-matrix structures for neighbours
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          watch_init_sym_neighbour_s.start();

          // D*M^T = (M*B)^T
          Adjacency::Graph graph_dm(Adjacency::RenderType::injectify_transpose, gate_v->_mirrors.at(i), this->_matrix_b.local());

          // S = (D*M^T) * B'
          Adjacency::Graph graph_s(Adjacency::RenderType::injectify, graph_dm, this->_neighbour_graphs.at(i));
          graph_s.sort_indices();

          // allocate Schur-matrix
          this->_neighbour_matrices.at(i) = LocalMatrixTypeS(graph_s);

          watch_init_sym_neighbour_s.stop();
        }

        // that's it
        watch_init_symbolic.stop();
      }

      /**
       * \brief Performs the numeric initialisation of the matrix.
       *
       * This function performs the following internal tasks:
       * -# Compute the numerical values of the local Schur-complement matrix S = (D*A*B).
       * -# For each neighbour, extract all matrix entries of our local matrix B, which
       *    coincide with entries of the reduced B-matrix for that particular neighbour,
       *    and send these matrix entries over to that neighbour process.
       * -# For each neighbour, receive the entries of the reduced B-matrix and pre-compute
       *    the numerical values of the neighbour Schur-complement matrix S_k := (D*A*B_k).
       */
      void init_numeric()
      {
        watch_init_numeric.start();

        // pre-multiply local matrix product
        watch_init_num_matrix_loc.start();
        this->_matrix_s.format();
        _asm_local_schur_matrix(this->_matrix_s, this->_matrix_d.local(), this->_diagonal_a.local(), this->_matrix_b.local());
        watch_init_num_matrix_loc.stop();

        // get the number of our neighbours
        const std::size_t num_neighs = this->_ranks.size();
        if(num_neighs <= std::size_t(0))
        {
          watch_init_numeric.stop();
          return; // no neighbours, no problems :)
        }

        // get our communicator
        const Dist::Comm& comm = *this->_matrix_b.get_comm();

        // send/receive buffers and requests
        std::vector<BufferVectorType> recv_bufs(num_neighs), send_bufs(num_neighs);
        Dist::RequestVector recv_reqs(num_neighs), send_reqs(num_neighs);

        // allocate receive buffer matrices B' and post receives
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          recv_bufs.at(i) = BufferVectorType(Index(dim) * this->_neighbour_graphs.at(i).get_num_indices());
          recv_reqs[i] = comm.irecv(recv_bufs.at(i).elements(), recv_bufs.at(i).size(), this->_ranks.at(i));
        }

        // extract reduced matrix data and post send
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          watch_init_num_gather_b.start();
          send_bufs.at(i) = _gather_b(this->_data_mirrors.at(i), this->_matrix_b.local());
          watch_init_num_gather_b.stop();
          send_reqs[i] = comm.isend(send_bufs.at(i).elements(), send_bufs.at(i).size(), this->_ranks.at(i));
        }

        // pre-multiply D*A and store in transposed form, i.e. CSC rather than CSR
        watch_init_num_premult_da.start();
        _premult_da(this->_matrix_da, this->_matrix_d.local(), this->_diagonal_a.local());
        watch_init_num_premult_da.stop();

        // process receives and compute neighbour schur matrices
        for(std::size_t i; recv_reqs.wait_any(i); )
        {
          watch_init_num_neighbour_s.start();
          this->_neighbour_matrices.at(i).format();
          _asm_neighbour_schur_matrix(this->_neighbour_matrices.at(i), this->_matrix_da,
            this->_matrix_b.get_row_gate()->_mirrors.at(i), this->_neighbour_graphs.at(i), recv_bufs.at(i));
          watch_init_num_neighbour_s.stop();
        }

        // wait for all previous sends to finish
        send_reqs.wait_all();
        watch_init_numeric.stop();
      }

      /**
       * \brief Calculate \f$ r \leftarrow this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       */
      void apply(VectorTypeL& r, const VectorTypeR& x) const
      {
        watch_apply.start();
        this->_apply(r.local(), x.local(), r.local(), DataType(1), true);
        watch_apply.stop();
      }

      /**
       * \brief Calculate \f$ r \leftarrow y + \alpha~ this\cdot x \f$
       *
       * \param[out] r The vector that receives the result.
       * \param[in] x The vector to be multiplied by this matrix.
       * \param[in] y The summand vector.
       * \param[in] alpha A scalar to scale the product with.
       */
      void apply(VectorTypeL& r, const VectorTypeR& x, const VectorTypeL& y, const DataType alpha = DataType(1)) const
      {
        watch_apply.start();
        this->_apply(r.local(), x.local(), y.local(), alpha, false);
        watch_apply.stop();
      }

      /**
       * \brief Assembles the matrix structure of an algebraic DOF partitioned matrix.
       *
       * This function computes the layout of the ADP matrix which corresponds to this
       * PMDCDSC matrix, i.e. the matrix, which contains all rows of the virtual global
       * matrix, which correspond to the global DOFs owned by this process.
       * This function is used by the Solver::ADPSolverBase class specialisation for
       * this class.
       *
       * \param[out] glob_dof_offset
       * Receives the global DOF offset of this process.
       *
       * \param[out] glob_dof_count
       * Receives the global DOF count over all processes.
       *
       * \returns
       * A new local matrix object containing the rows of all DOFs that belong to this process.
       */
      LocalMatrixTypeS asm_adp_symbolic(Index& glob_dof_offset, Index& glob_dof_count) const
      {
        // no neighbours?
        if(_ranks.empty())
        {
          glob_dof_offset = Index(0);
          glob_dof_count = _matrix_s.rows();
          return _matrix_s.clone(LAFEM::CloneMode::Layout);
        }

        // get our communicator
        const Dist::Comm& comm = *this->get_comm();

        const std::size_t num_neighs = this->_ranks.size();

        // get number of local DOFs
        const Index num_loc_dofs = _matrix_s.rows();

        // compute our global DOF offset and count
        comm.exscan(&num_loc_dofs, &glob_dof_offset, std::size_t(1), Dist::op_sum);
        comm.allreduce(&num_loc_dofs, &glob_dof_count, std::size_t(1), Dist::op_sum);

        // The columns of our neighbour matrices correspond to the entries in the pressure mirror.
        // However, for the desired ADP matrix, we have to translate these into global DOF indices.
        // For this, each process has to map the DOF in its pressure mirrors to global DOFs and then
        // send these DOF indices to the corresponding neighbour, so that it can map the column
        // indices of its neighbour matrix to global DOF indices.

        // send/receive mirrors and requests
        std::vector<std::vector<IndexType>> recv_dofs(num_neighs), send_dofs(num_neighs);
        Dist::RequestVector recv_reqs(num_neighs), send_reqs(num_neighs);

        // allocate receive vectors and post receives
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          recv_dofs.at(i).resize(_neighbour_graphs.at(i).get_num_nodes_image());
          recv_reqs[i] = comm.irecv(recv_dofs.at(i).data(), recv_dofs.at(i).size(), this->_ranks.at(i));
        }

        // translate our pressure mirrors to global DOF vectors and post sends
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          const Index num_idx = _pres_mirrors.at(i).num_indices();
          const IndexType* pidx = _pres_mirrors.at(i).indices();
          send_dofs.at(i).resize(num_idx);
          IndexType* sidx = send_dofs.at(i).data();
          for(Index k(0); k < num_idx; ++k)
            sidx[k] = glob_dof_offset + pidx[k];
          send_reqs[i] = comm.isend(send_dofs.at(i).data(), send_dofs.at(i).size(), this->_ranks.at(i));
        }

        // get local matrix stuff
        const Index num_rows = num_loc_dofs;
        const IndexType* row_ptr_s = _matrix_s.row_ptr();
        const IndexType* col_idx_s = _matrix_s.col_ind();

        // compute number of non-zeros per row and total
        Index num_nzes = _matrix_s.used_elements();
        std::vector<IndexType> row_aux(num_rows, IndexType(0));
        for(Index i(0); i  < num_rows; ++i)
          row_aux[i] = (row_ptr_s[i+1] - row_ptr_s[i]);

        for(const auto& x : _neighbour_matrices)
        {
          num_nzes += x.used_elements();
          const IndexType* row_ptr_x = x.row_ptr();
          for(Index i(0); i  < num_rows; ++i)
            row_aux[i] += (row_ptr_x[i+1] - row_ptr_x[i]);
        }

        // allocate our matrix
        LocalMatrixTypeS matrix_g(num_rows, glob_dof_count, num_nzes);

        // get our index arrays
        IndexType* row_ptr_g = matrix_g.row_ptr();
        IndexType* col_idx_g = matrix_g.col_ind();

        // compute row pointer array and store backup in aux
        row_ptr_g[0] = IndexType(0);
        for(Index i(0); i  < num_rows; ++i)
        {
          row_ptr_g[i+1] = row_ptr_g[i] + row_aux[i];
          row_aux[i] = row_ptr_g[i];
        }

        // Note: For the sake of compatibility with picky third-party libraries, we want to
        // ensure that the column indices of the output matrix are in ascending order.
        // For this, we have to combine the matrix layout from our own local matrix and
        // the matrices of our neighbours in rank-ascending order. So we first create two
        // rank maps for our neighbours with lower and higher ranks, so that we can easily
        // loop over all matrices in rank-order.

        // create two neighbours maps: one of all lower ranks and one of all higher ranks
        std::map<int, std::size_t> neigh_map_l, neigh_map_h;
        for(std::size_t ineigh(0); ineigh < num_neighs; ++ineigh)
        {
          if(_ranks.at(ineigh) < comm.rank())
            neigh_map_l.emplace(_ranks.at(ineigh), ineigh);
          else
            neigh_map_h.emplace(_ranks.at(ineigh), ineigh);
        }

        // wait for all receive requests to finish
        recv_reqs.wait_all();

        // first, insert all neighbour matrices with a lower rank in rank-ascending order
        for(auto it = neigh_map_l.begin(); it != neigh_map_l.end(); ++it)
        {
          std::size_t ineigh = it->second;
          const IndexType* row_ptr_x = _neighbour_matrices.at(ineigh).row_ptr();
          const IndexType* col_idx_x = _neighbour_matrices.at(ineigh).col_ind();
          const IndexType* dof_idx_x = recv_dofs.at(ineigh).data();
          for(Index i(0); i < num_rows; ++i)
          {
            IndexType k = row_aux[i];
            for(IndexType j(row_ptr_x[i]); j < row_ptr_x[i + 1]; ++j, ++k)
              col_idx_g[k] = dof_idx_x[col_idx_x[j]];
            row_aux[i] = k;
          }
        }

        // now insert our local matrix S
        for(Index i(0); i < num_rows; ++i)
        {
          IndexType k = row_aux[i];
          for(IndexType j(row_ptr_s[i]); j < row_ptr_s[i + 1]; ++j, ++k)
            col_idx_g[k] = glob_dof_offset + col_idx_s[j];
          row_aux[i] = k;
        }

        // finally, insert all neighbour matrices with a higher rank in rank-ascending order
        for(auto it = neigh_map_h.begin(); it != neigh_map_h.end(); ++it)
        {
          std::size_t ineigh = it->second;
          const IndexType* row_ptr_x = _neighbour_matrices.at(ineigh).row_ptr();
          const IndexType* col_idx_x = _neighbour_matrices.at(ineigh).col_ind();
          const IndexType* dof_idx_x = recv_dofs.at(ineigh).data();
          for(Index i(0); i < num_rows; ++i)
          {
            IndexType k = row_aux[i];
            for(IndexType j(row_ptr_x[i]); j < row_ptr_x[i + 1]; ++j, ++k)
              col_idx_g[k] = dof_idx_x[col_idx_x[j]];
            row_aux[i] = k;
          }
        }

        // wait for all send requests to finish
        send_reqs.wait_all();

#ifdef DEBUG
        // sanity check: ensure that the column indices are sorted correctly
        for(Index i = 0; i < num_rows; ++i)
        {
          for(IndexType j(row_ptr_g[i]); j + 1 < row_ptr_g[i + 1]; ++j)
          {
            ASSERT(col_idx_g[j] < col_idx_g[j+1]);
          }
        }
#endif // DEBUG

        // that's it
        return matrix_g;
      }

      /**
       * \brief Copies the numeric values of the matrix into an ADP matrix.
       *
       * \param[inout] matrix
       * A reference to the ADP matrix which receives the copy of the values,
       * as returned by the #asm_adp_symoblic() function.
       */
      void asm_adp_numeric(LocalMatrixTypeS& matrix) const
      {
        // no neighbours?
        if(_ranks.empty())
        {
          // copy values
          matrix.copy(_matrix_s);
          return;
        }

        // get number of local DOFs
        const Index num_rows = _matrix_s.rows();
        XASSERT(matrix.rows() == _matrix_s.rows());

        // get row pointer arrays
        const IndexType* row_ptr_s = _matrix_s.row_ptr();
        const IndexType* row_ptr_g = matrix.row_ptr();
        const DataType* val_s = _matrix_s.val();
        DataType* val_g = matrix.val();

        // make a copy of the row pointer
        std::vector<IndexType> row_aux(num_rows);

        // copy the row-pointer array
        for(Index i(0); i < num_rows; ++i)
          row_aux[i] = row_ptr_g[i];

        // get this process's rank
        const int my_rank = this->get_comm()->rank();

        // create two neighbours maps: one of all lower ranks and one of all higher ranks
        std::map<int, std::size_t> neigh_map_l, neigh_map_h;
        for(std::size_t ineigh(0); ineigh < _ranks.size(); ++ineigh)
        {
          if(_ranks.at(ineigh) < my_rank)
            neigh_map_l.emplace(_ranks.at(ineigh), ineigh);
          else
            neigh_map_h.emplace(_ranks.at(ineigh), ineigh);
        }

        // first, copy all neighbour matrices with a lower rank in rank-ascending order
        for(auto it = neigh_map_l.begin(); it != neigh_map_l.end(); ++it)
        {
          std::size_t ineigh = it->second;
          const IndexType* row_ptr_x = _neighbour_matrices.at(ineigh).row_ptr();
          const DataType* val_x = _neighbour_matrices.at(ineigh).val();
          for(Index i(0); i < num_rows; ++i)
          {
            IndexType k = row_aux[i];
            for(IndexType j(row_ptr_x[i]); j < row_ptr_x[i + 1]; ++j, ++k)
              val_g[k] = val_x[j];
            row_aux[i] = k;
          }
        }

        // now copy our own matrix
        for(Index i(0); i < num_rows; ++i)
        {
          Index k = row_aux[i];
          for(Index j(row_ptr_s[i]); j < row_ptr_s[i+1]; ++j, ++k)
            val_g[k] = val_s[j];
          row_aux[i] = k;
        }

        // finally, copy all neighbour matrices with a higher rank in rank-ascending order
        for(auto it = neigh_map_h.begin(); it != neigh_map_h.end(); ++it)
        {
          std::size_t ineigh = it->second;
          const IndexType* row_ptr_x = _neighbour_matrices.at(ineigh).row_ptr();
          const DataType* val_x = _neighbour_matrices.at(ineigh).val();
          for(Index i(0); i < num_rows; ++i)
          {
            IndexType k = row_aux[i];
            for(IndexType j(row_ptr_x[i]); j < row_ptr_x[i + 1]; ++j, ++k)
              val_g[k] = val_x[j];
            row_aux[i] = k;
          }
        }

#ifdef DEBUG
        // sanity check: ensure that all entries have been processed
        for(Index i(0); i < num_rows; ++i)
        {
          ASSERT(row_aux[i] == row_ptr_g[i+1]);
        }
#endif // DEBUG
      }

    protected:
      /// auxiliary function: format a time line
      static String _fmt_time(double tsum_total, double tmax_total, String st)
      {
        String s = st.pad_back(30, '.') + ":";
        s += stringify_fp_fix(tsum_total, 6, 12) + " :";
        s += stringify_fp_fix(tmax_total, 6, 12) + "\n";
        return s;
      }
      /// auxiliary function: format a time line
      static String _fmt_time(double tsum_total, double tmax_total, double tsum, double tmax, String st)
      {
        String s = st.pad_back(30, '.') + ":";
        s += stringify_fp_fix(tsum, 6, 12) + " :";
        s += stringify_fp_fix(tmax, 6, 12) + " [";
        if(tsum_total > 1E-8 * Math::abs(tsum))
          s += stringify_fp_fix(100.0*tsum/tsum_total, 2, 6) + "% :";
        else
          s += "    --- :";
        if(tmax_total > 1E-8 * Math::abs(tmax))
          s += stringify_fp_fix(100.0*tmax/tmax_total, 2, 6) + "% ]\n";
        else
          s += "    --- ]\n";
        return s;
      }

      /**
       * \brief Auxiliary function: assembles a pressure mirror from a velocity mirror and the B-matrix
       *
       * \param[out] mirror_p
       * The mirror of all pressure DOFs that are coupled by at least one non-zero entry (each) of the
       * given B-matrix to the velocity DOFs stored in the velocity mirror.
       *
       * \param[in] mirror_v
       * The velocity mirror that contains the indices of the rows of B that are to be extracted.
       *
       * \param[in] matrix_b
       * The B-matrix which couples velocity and pressure DOFs.
       */
      static void _asm_pres_mirror(MirrorTypeP& mirror_p, const MirrorTypeV& mirror_v, const LocalMatrixTypeB& matrix_b)
      {
        const Index num_dof_mir_v = mirror_v.num_indices();
        const IndexType* velo_idx = mirror_v.indices();
        const IndexType* row_ptr = matrix_b.row_ptr();
        const IndexType* col_idx = matrix_b.col_ind();

        // loop over all rows of B, which are indexed in the velocity mirror,
        // and add all column indices (pressure DOFs) into the pressure DOF set
        std::set<IndexType> dof_set;
        for(Index i(0); i < num_dof_mir_v; ++i)
        {
          const Index irow = velo_idx[i];
          for(IndexType j(row_ptr[irow]); j < row_ptr[irow+1]; ++j)
            dof_set.insert(col_idx[j]);
        }

        // convert DOF set into a mirror
        mirror_p = MirrorTypeP(matrix_b.columns(), Index(dof_set.size()));
        IndexType* pidx = mirror_p.indices();
        for(auto it = dof_set.begin(); it != dof_set.end(); ++it, ++pidx)
          *pidx = *it;
      }

      /**
       * \brief Auxiliary Function: reduces the matrix B to B'
       *
       * \param[in] mirror_p
       * The pressure mirror as returned by the _asm_pres_mirror() function.
       *
       * \param[in] mirror_v
       * The velocity mirror that contains the indices of the rows of B that are to be extracted to B'.
       *
       * \param[in] matrix_b
       * The B-matrix which couples velocity and pressure DOFs.
       *
       * \returns
       * An adjacency graph storing the structure of the reduced matrix B'
       */
      static Adjacency::Graph _asm_reduced_b(const MirrorTypeP& mirror_p, const MirrorTypeV& mirror_v,
        const LocalMatrixTypeB& matrix_b)
      {
        const Index num_dof_mir_v = mirror_v.num_indices();
        const Index num_dof_mir_p = mirror_p.num_indices();
        const IndexType* velo_idx = mirror_v.indices();
        const IndexType* pres_idx = mirror_p.indices();
        const IndexType* row_ptr = matrix_b.row_ptr();
        const IndexType* col_idx = matrix_b.col_ind();

        // count number of non-zeros in indexed rows of B = non-zeros in B'
        Index num_red_nzes = Index(0);
        for(Index i(0); i < num_dof_mir_v; ++i)
        {
          const Index irow = velo_idx[i];
          num_red_nzes += Index(row_ptr[irow+1] - row_ptr[irow]);
        }

        // allocate matrix graph for reduced part B' of B
        Adjacency::Graph graph(num_dof_mir_v, num_dof_mir_p, num_red_nzes);
        Index* dom_ptr = graph.get_domain_ptr();
        Index* img_idx = graph.get_image_idx();

        // loop over all rows of reduced part B'
        dom_ptr[0] = Index(0);
        for(Index i(0); i < num_dof_mir_v; ++i)
        {
          // get the B-row index of the i-th row of B'
          const IndexType irow = velo_idx[i];
          // loop over all non-zeroes in the row of B / B'
          for(IndexType j(row_ptr[irow]), k(dom_ptr[i]); j < row_ptr[irow+1]; ++j, ++k)
          {
            // initialise invalid index for assertion below
            img_idx[k] = ~IndexType(0);

            // try to find this column index (=pressure DOF) in our pressure DOF mirror
            for(Index l(0); l < num_dof_mir_p; ++l)
            {
              if(col_idx[j] == pres_idx[l])
              {
                // that's our pressure DOF, so store its index as column index of B'
                img_idx[k] = l;
                break;
              }
            }
            ASSERT(img_idx[k] != ~IndexType(0));
          }
          // set next row pointer of B'
          dom_ptr[i+1] = dom_ptr[i] + (row_ptr[irow+1] - row_ptr[irow]);
        }

        // sort indices of B' and return graph
        graph.sort_indices();
        return graph;
      }

      /**
       * \brief Auxiliary Function: computes a data mirror from B to B'
       *
       * \param[in] mirror_p
       * The pressure mirror as returned by the _asm_pres_mirror() function.
       *
       * \param[in] mirror_v
       * The velocity mirror that contains the indices of the rows of B that are to be extracted to B'.
       *
       * \param[in] matrix_b
       * The B-matrix which couples velocity and pressure DOFs.
       *
       * \param[in] graph
       * The adjacency graph storing the structure of the reduced matrix B'
       *
       * \returns
       * A mirror for the data array of B to B'
       */
      static MirrorTypeP _asm_data_mirror(const MirrorTypeP& mirror_p, const MirrorTypeV& mirror_v,
        const LocalMatrixTypeB& matrix_b, const Adjacency::Graph& graph)
      {
        const Index num_dof_mir_v = mirror_v.num_indices();
        const IndexType* velo_idx = mirror_v.indices();
        const IndexType* pres_idx = mirror_p.indices();
        const IndexType* row_ptr = matrix_b.row_ptr();
        const IndexType* col_idx = matrix_b.col_ind();
        const Index* dom_ptr = graph.get_domain_ptr();
        const Index* img_idx = graph.get_image_idx();

        // allocate mirror for data array of B'
        MirrorTypeP data_mirror(matrix_b.used_elements(), graph.get_num_indices());
        IndexType* dat_idx = data_mirror.indices();

        // loop over all rows of B' again
        for(Index i(0); i < num_dof_mir_v; ++i)
        {
          // get the B-row index (=velocity DOF index)
          const IndexType irow = velo_idx[i];

          // loop over all non-zeroes in the row of B / B'
          for(Index k(dom_ptr[i]); k < dom_ptr[i+1]; ++k)
          {
            // initialise invalid index for assertion below
            dat_idx[k] = ~IndexType(0);

            // get the B-column index (=pressure DOF index)
            const Index jcol = pres_idx[img_idx[k]];

            // try to find that column in our matrix
            for(IndexType j(row_ptr[irow]); j < row_ptr[irow + 1]; ++j)
            {
              if(jcol == col_idx[j])
              {
                // that's the entry we were looking for
                dat_idx[k] = j;
                break;
              }
            }
            ASSERT(dat_idx[k] != ~IndexType(0));
          }
        }

        // that's it
        return data_mirror;
      }

      /**
       * \brief Auxiliary Function: Gathers the data array from B to B' using the data mirror
       *
       * \param[in] data_mirror
       * The data mirror from B to B' as returned by the #_asm_data_mirror() function.
       *
       * \param[in] matrix_b
       * The input matrix B from which the values for B' are to be gathered.
       *
       * \returns
       * A buffer vector containing the gathered data array of B'.
       */
      static BufferVectorType _gather_b(const MirrorTypeP& data_mirror, const LocalMatrixTypeB& matrix_b)
      {
        const Index num_idx = data_mirror.num_indices();
        const IndexType* idx = data_mirror.indices();
        const ValueTypeB* mat_val = matrix_b.val();

        BufferVectorType buf(Index(dim)*num_idx);
        ValueTypeB* buf_val = reinterpret_cast<ValueTypeB*>(buf.elements());

        for(Index i(0); i < num_idx; ++i)
        {
          buf_val[i] = mat_val[idx[i]];
        }

        return buf;
      }

      /// auxiliary function for _asm_local_schur_matrix: multiply two values D*A
      inline static ValueTypeD _mat_mult_d_a(const ValueTypeD& val_d, const ValueTypeA& val_a)
      {
        ValueTypeD da;
        for(int i(0); i < dim; ++i)
          da(0, i) = val_d(0, i) * val_a(i);
        return da;
      }

      /// auxiliary function for _asm_local_schur_matrix: multiply two values DA*B
      inline static DataType _mat_mult_da_b(const ValueTypeD& val_da, const ValueTypeB& val_b)
      {
        DataType s = DataType(0);
        for(int i(0); i < dim; ++i)
          s += val_da(0, i) * val_b(i, 0);
        return s;
      }

      /// auxiliary function for _asm_neighbour_schur_matrix: multiply two values D*A and store in transposed form
      inline static void _mat_mult_d_a(ValueTypeB& val_da, const ValueTypeD& val_d, const ValueTypeA& val_a)
      {
        for(int i(0); i < dim; ++i)
          val_da(i, 0) = val_d(0, i) * val_a(i);
      }

      /// auxiliary function for _asm_neighbour_schur_matrix: multiply two values DA*B with DA in transposed form
      inline static DataType _mat_mult_da_b(const ValueTypeB& val_da, const ValueTypeB& val_b)
      {
        DataType s = DataType(0);
        for(int i(0); i < dim; ++i)
          s += val_da(i, 0) * val_b(i, 0);
        return s;
      }

      /**
       * \brief Auxiliary Function: Computes the product (D*A) and stores the result in a CSC matrix.
       *
       * \param[inout] mat_da
       * The matrix that receives the product (D*A).
       *
       * \param[in] mat_d
       * The multiplicand matrix D.
       *
       * \param[in] mat_a
       * A vector representing the diagonal matrix A.
       */
      static void _premult_da(LocalMatrixTypeB& mat_da, const LocalMatrixTypeD& mat_d, const LocalVectorTypeV& mat_a)
      {
        const Index num_rows = mat_d.rows();
        const Index num_cols = mat_d.columns();

        const IndexType* row_ptr = mat_d.row_ptr();
        const IndexType* col_idx = mat_d.col_ind();
        const IndexType* row_ptr_da = mat_da.row_ptr();

        ValueTypeB* val_da = mat_da.val();
        const ValueTypeD* val_d = mat_d.val();
        const ValueTypeA* val_a = mat_a.elements();

        // create a temporary copy of the row-pointer of D^T
        std::vector<Index> ptr(num_cols);
        for(Index i(0); i < num_cols; ++i)
          ptr[i] = row_ptr_da[i];

        // transpose matrix
        for(Index i(0); i < num_rows; ++i)
        {
          for(Index j(row_ptr[i]); j < row_ptr[i + 1]; ++j)
          {
            const Index col = col_idx[j];
            _mat_mult_d_a(val_da[ptr[col]++], val_d[j], val_a[col]);
          }
        }
      }

      /**
       * \brief Assembles the local Schur-Matrix S = (D*A*B)
       *
       * \param[inout] s
       * The Schur-Matrix S that is to be assembled.
       *
       * \param[in] d, a, b
       * The three multiplicand matrices.
       */
      static void _asm_local_schur_matrix(LocalMatrixTypeS& s, const LocalMatrixTypeD& d,
        const LocalVectorTypeV& a, const LocalMatrixTypeB& b)
      {
        // Note: this is a modified version of SparseMatrixCSR::add_double_mat_mult for BCSR multiplicands

        // validate matrix dimensions
        XASSERT(s.rows() == d.rows());
        XASSERT(d.columns() == a.size());
        XASSERT(a.size() == b.rows());
        XASSERT(b.columns() == s.columns());

        // fetch matrix arrays:
        DataType* data_s = s.val();
        const ValueTypeD* data_d = d.val();
        const ValueTypeA* data_a = a.elements();
        const ValueTypeB* data_b = b.val();
        const IndexType* row_ptr_s = s.row_ptr();
        const IndexType* col_idx_s = s.col_ind();
        const IndexType* row_ptr_d = d.row_ptr();
        const IndexType* col_idx_d = d.col_ind();
        const IndexType* row_ptr_b = b.row_ptr();
        const IndexType* col_idx_b = b.col_ind();

        // loop over all rows of D and S, resp.
        for(IndexType i(0); i < IndexType(s.rows()); ++i)
        {
          // loop over all non-zeros D_ik in row i of D
          for(IndexType ik(row_ptr_d[i]); ik  < row_ptr_d[i+1]; ++ik)
          {
            // get column index k
            const IndexType k = col_idx_d[ik];

            // pre-compute (D_ik * A_kk)
            const ValueTypeD val_da = _mat_mult_d_a(data_d[ik], data_a[k]);

            //   S_i. += (D_ik * A_kk) * B_k.
            for(IndexType ij(row_ptr_s[i]), kj(row_ptr_b[k]); kj < row_ptr_b[k+1]; ++ij)
            {
              ASSERT(ij < row_ptr_s[i+1]);
              ASSERT(col_idx_s[ij] <= col_idx_b[kj]);
              if(col_idx_s[ij] == col_idx_b[kj])
              {
                data_s[ij] += _mat_mult_da_b(val_da, data_b[kj]);
                ++kj;
              }
            }
          }
        }
      }

      /**
       * \brief Assembles a neighbour Schur-Matrix S_k = (D*A*M_k^T*B_k)
       *
       * \param[inout] s
       * The neighbour Schur-Matrix S_k that is to be assembled.
       *
       * \param[in] da
       * The pre-multiplied matrix (D*A) stored in CSC format (transposed CSR)
       *
       * \param[in] mirror_v
       * The velocity mirror representing the mirror matrix M_k.
       *
       * \param[in] graph_b
       * The adjacency graph representing the layout of the matrix B_k.
       *
       * \param[in] buffer_b
       * A buffer vector containing the data array of the matrix B_k.
       */
      static void _asm_neighbour_schur_matrix(LocalMatrixTypeS& s, const LocalMatrixTypeB& da,
        const MirrorTypeV& mirror_v, const Adjacency::Graph& graph_b, const BufferVectorType& buffer_b)
      {
        // validate matrix dimensions
        XASSERT(s.rows() == da.columns());
        XASSERT(da.rows() == mirror_v.size());
        XASSERT(graph_b.get_num_nodes_domain() == mirror_v.num_indices());
        XASSERT(graph_b.get_num_nodes_image() == s.columns());
        XASSERT(mirror_v.size() == da.rows());

        // fetch matrix arrays:
        DataType* data_s = s.val();
        const IndexType* row_ptr_s = s.row_ptr();
        const IndexType* col_idx_s = s.col_ind();

        // we use CSR for (D*A)^T here, which is effectively a CSC storage of (D*A),
        // thus rows and columns swap their meaning here
        const ValueTypeB* data_da = da.val();
        const IndexType* col_ptr_da = da.row_ptr();
        const IndexType* row_idx_da = da.col_ind();

        const Index num_mir_idx = mirror_v.num_indices();
        const IndexType* mir_idx = mirror_v.indices();

        const ValueTypeB* data_b = reinterpret_cast<const ValueTypeB*>(buffer_b.elements());
        const Index* dom_ptr_b = graph_b.get_domain_ptr();
        const Index* img_idx_b = graph_b.get_image_idx();

        // loop over all velocity mirror entries, which correspond to the rows of B and the columns of D
        for(Index l(0); l < num_mir_idx; ++l)
        {
          // get velocity DOF index k = column index of D = row index of B
          const Index k = mir_idx[l];

          // loop over all columns k of D
          for(Index ik(col_ptr_da[k]); ik < col_ptr_da[k + 1]; ++ik)
          {
            // get the row index i of D_ik
            const Index i = row_idx_da[ik];

            //   S_i. += (D_ik * A_kk) * B_l.
            for(IndexType ij(row_ptr_s[i]), kj(dom_ptr_b[l]); kj < dom_ptr_b[l+1]; ++ij)
            {
              ASSERT(ij < row_ptr_s[i+1]);
              ASSERT(col_idx_s[ij] <= img_idx_b[kj]);
              if(col_idx_s[ij] == img_idx_b[kj])
              {
                data_s[ij] += _mat_mult_da_b(data_da[ik], data_b[kj]);
                ++kj;
              }
            }
          }
        }
      }

      /**
       * \brief Applies this matrix onto a vector.
       *
       * \param[inout] r
       * The vector that receives the result.
       *
       * \param[in] x
       * The multiplicand vector.
       *
       * \param[in] x
       * The summand vector.
       *
       * \param[in] alpha
       * The scaling factor for the product.
       *
       * \param[in] only_Ax
       * If \c true, then compute r := A*x, otherwise compute r := y + alpha*A*x.
       */
      void _apply(LocalVectorTypeP& r, const LocalVectorTypeP& x, const LocalVectorTypeP& y, const DataType alpha, bool only_Ax) const
      {
        // get the number of our neighbours
        const std::size_t num_neighs = this->_ranks.size();

        // no neighbours?
        if(num_neighs <= std::size_t(0))
        {
          // multiply by local schur matrix and return
          watch_apply_matrix_loc.start();
          if(only_Ax)
            this->_matrix_s.apply(r, x);
          else
            this->_matrix_s.apply(r, x, y, alpha);
          watch_apply_matrix_loc.stop();
          return;
        }

        // get our communicator
        const Dist::Comm& comm = *this->_matrix_b.get_comm();

        // send/receive buffers and requests
        std::vector<BufferVectorType> recv_bufs(num_neighs), send_bufs(num_neighs);
        Dist::RequestVector recv_reqs(num_neighs), send_reqs(num_neighs);

        // allocate receive buffer vectors and post receives
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          recv_bufs.at(i) = BufferVectorType(this->_neighbour_matrices.at(i).columns());
          recv_reqs[i] = comm.irecv(recv_bufs.at(i).elements(), recv_bufs.at(i).size(), this->_ranks.at(i));
        }

        // extract pressure dofs and post sends
        for(std::size_t i(0); i < num_neighs; ++i)
        {
          send_bufs.at(i) = BufferVectorType(this->_pres_mirrors.at(i).num_indices());
          this->_pres_mirrors.at(i).gather(send_bufs.at(i), x);
          send_reqs[i] = comm.isend(send_bufs.at(i).elements(), send_bufs.at(i).size(), this->_ranks.at(i));
        }

        // multiply by local schur matrix
        watch_apply_matrix_loc.start();
        if(only_Ax)
          this->_matrix_s.apply(r, x);
        else
          this->_matrix_s.apply(r, x, y, alpha);
        watch_apply_matrix_loc.stop();

        // process receives and mulitply by neighbour schur matrices
        for(std::size_t i; recv_reqs.wait_any(i); )
        {
          watch_apply_neighbour_s.start();
          this->_neighbour_matrices.at(i).apply(r, recv_bufs.at(i), r, (only_Ax ? DataType(1) : alpha));
          watch_apply_neighbour_s.stop();
        }

        // wait for all previous sends to finish
        send_reqs.wait_all();
      }
    }; // class PMDCDSCMatrix<...>
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_PMDCDSC_MATRIX_HPP
