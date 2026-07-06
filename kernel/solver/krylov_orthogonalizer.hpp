// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/lafem/arch/batched_vector_apply.hpp>
#include <kernel/lafem/base.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>

#include <vector>

namespace FEAT::Solver
{
  namespace Intern
  {

    template<typename VectorType_>
    class BatchedHelper;
  }

  /**
   * \brief Enumeration of the available Krylov subspace orthogonalization strategies
   */
  enum class KrylovOrthogonalizerType
  {
    /// modified Gram-Schmidt: most stable, but fully sequentializes the inner dot products, i.e. it
    /// requires one global reduction per basis vector (see ModifiedGramSchmidtOrthogonalizer)
    modified_gram_schmidt,
    /// classical Gram-Schmidt: all inner dot products of one step are computed in a single batched
    /// reduction, trading stability for a single global reduction (see ClassicGramSchmidtOrthogonalizer)
    classical_gram_schmidt,
    /// classical Gram-Schmidt with a (conditional) reorthogonalization sweep, recovering a stability
    /// close to modified Gram-Schmidt while keeping the reductions batched (see CGS2Orthogonalizer)
    cgs2
  };

  /**
   * \brief Base class for orthogonlizing a vector of a Krylov subspace
   */
  template<typename Vector_>
  class KrylovSubspaceOrthogonalizer
  {
  public:
   typedef Vector_ VectorType;
   typedef typename Vector_::DataType DataType;

  protected:
    /// The component holds the krylov subspace (and potentially additonal helper vectors)
    std::vector<VectorType> _subspace;
    /// The Krylov dimension
    Index _krylov_dimension = 0u;


  public:

    KrylovSubspaceOrthogonalizer(const VectorType& vec_ref, Index krylov_dimension) :
     _krylov_dimension(krylov_dimension)
    {
      _subspace.reserve(_krylov_dimension+1);
      for(std::size_t n = 0; n <= krylov_dimension; ++n)
      {
        _subspace.emplace_back(vec_ref.clone(LAFEM::CloneMode::Layout));
      }
    }

    explicit KrylovSubspaceOrthogonalizer() = default;
    KrylovSubspaceOrthogonalizer(const KrylovSubspaceOrthogonalizer&) = delete;
    KrylovSubspaceOrthogonalizer& operator=(const KrylovSubspaceOrthogonalizer&) = delete;
    KrylovSubspaceOrthogonalizer(KrylovSubspaceOrthogonalizer&&) = default;
    KrylovSubspaceOrthogonalizer& operator=(KrylovSubspaceOrthogonalizer&&) = default;

    virtual ~KrylovSubspaceOrthogonalizer()=default;

    /**
     * \brief Orthogonalizes v_[ti] against v_0,...v_{ti-1}
     *
     * \param[out] dot_product A \transisitive array of size target_index in which the dot products should be written
     * \param[in] target_index The index of target inside the underlying vector array.
     *
     * \returns The euclid norm of the orthogonalized target
     *
     * \note target has to be in the vector array provided to this component.
     * \note You should always call this interface in a sequential manner regarding target_index, else
     *       the subspace is not guaranteed to be orthogonal
     */
    virtual DataType orthogonalize(DataType* dot_products, const Index target_index) = 0;

    /**
     * \brief Performs the commulative sum on vector target with the internal vector array
     */
    virtual void correct_vector(VectorType& target, const DataType* local_scales, const Index num_axpy, const DataType alpha) const
    {
      XASSERT(num_axpy <= _krylov_dimension);
      for(Index k = 0; k < num_axpy; ++k)
      {
        target.axpy(_subspace[k], alpha * local_scales[k]);
      }
    }

    /**
     * \brief The krylov dimension
     */
    Index get_krylov_dim() const
    {
      return _krylov_dimension;
    }

    /**
     * \brief Get a reference to the subspace vectors
     */
    std::vector<VectorType>& get_subspace()
    {
      return _subspace;
    }

    const std::vector<VectorType>& get_subspace() const
    {
      return _subspace;
    }

    VectorType& at(const Index idx)
    {
      return _subspace.at(idx);
    }

    const VectorType& at(const Index idx) const
    {
      return _subspace.at(idx);
    }
  }; // class KrylovSubspaceOrthogonalizer<...>

  /**
   * \brief Simple abstraction of the modified Gram Schmidt procedure
   */
  template<typename Vector_>
  class ModifiedGramSchmidtOrthogonalizer : public KrylovSubspaceOrthogonalizer<Vector_>
  {
  public:
    typedef KrylovSubspaceOrthogonalizer<Vector_> BaseClass;
    typedef typename BaseClass::VectorType VectorType;
    typedef typename BaseClass::DataType DataType;

  public:
    explicit ModifiedGramSchmidtOrthogonalizer() = default;

    ModifiedGramSchmidtOrthogonalizer(const Vector_& vec_ref, Index krylov_dimension) :
     BaseClass(vec_ref, krylov_dimension)
    {
    }

    ModifiedGramSchmidtOrthogonalizer(const ModifiedGramSchmidtOrthogonalizer&) = delete;
    ModifiedGramSchmidtOrthogonalizer& operator=(const ModifiedGramSchmidtOrthogonalizer&) = delete;
    ModifiedGramSchmidtOrthogonalizer(ModifiedGramSchmidtOrthogonalizer&&) = default;
    ModifiedGramSchmidtOrthogonalizer& operator=(ModifiedGramSchmidtOrthogonalizer&&) = default;

    virtual ~ModifiedGramSchmidtOrthogonalizer()=default;

    DataType orthogonalize(DataType* dot_products, const Index target_index) override
    {
      XASSERT(target_index <= this->_krylov_dimension);
      // Gram-Schmidt process
      for(Index k(0); k < target_index; ++k)
      {
        const auto dot = this->_subspace[target_index].dot(this->_subspace[k]);
        this->_subspace[target_index].axpy(this->_subspace[k], -dot);
        dot_products[k] = dot;
      }
      return this->_subspace[target_index].norm2();
    }
  };

  /**
   * \brief Abstraction of the classical Gram-Schmidt procedure with a batched reduction
   *
   * In contrast to the modified Gram-Schmidt procedure, the classical variant computes all inner
   * dot products of the target vector against the existing basis vectors independently, which allows
   * them to be evaluated in a \e single batched reduction (and thus a single global all-reduce in the
   * distributed case) instead of one reduction per basis vector.
   *
   * In addition, the self-dot \f$\langle w, w\rangle\f$ is appended to the same batched reduction, so
   * the euclidean norm of the orthogonalized target can be recovered algebraically as
   * \f$\sqrt{\langle w,w\rangle - \sum_k h_k^2}\f$ without a separate norm2 reduction. The algebraic
   * norm is only trusted while it stays safely positive; on heavy cancellation we fall back to an
   * explicit norm2, which guarantees a positive result so the subsequent normalization cannot produce
   * a NaN.
   */
  template<typename Vector_>
  class ClassicGramSchmidtOrthogonalizer : public KrylovSubspaceOrthogonalizer<Vector_>
  {
  public:
    typedef KrylovSubspaceOrthogonalizer<Vector_> BaseClass;
    typedef typename BaseClass::VectorType VectorType;
    typedef typename BaseClass::DataType DataType;

  protected:
    /// flattened element arbiters of the persistent basis vectors, consumed by the batched kernels
    std::vector<Memory::Arbiter const*> _vec_arbiter;
    /// device-resident accumulation buffer the batched dot kernels write into; the arbiter mediates
    /// the host<->device transfers and lets the kernels accumulate across (meta-)components
    Memory::Arbiter _dots;

    /**
     * \brief Computes the batched dot products of \p target against the leading basis vectors
     *
     * A single batched reduction fills <c>dots_out[0..num_basis)</c> with
     * \f$\langle \text{target}, v_k\rangle\f$ and returns the self-dot
     * \f$\langle \text{target}, \text{target}\rangle\f$ (basis slot \p num_basis) in \p self_dot.
     *
     * \param[in] target The vector all dot products are computed with.
     * \param[in] num_basis The number of leading basis vectors to project against.
     * \param[out] dot_out A host array of size \p num_basis receiving the inner dot products.
     * \param[out] self_dot Receives the self-dot of \p target.
     */
    void _compute_dots(const VectorType& target, const Index num_basis, DataType* dot_out, DataType& self_dot)
    {
      // zero the accumulation arbiter on all locations, then accumulate the batched dots onto it
      _dots.format();
      // the batched dot also picks up the self-dot: target == _subspace[num_basis], so slot num_basis
      // of the reduction holds <target,target>
      Intern::BatchedHelper<VectorType>::template batched_dot<DataType>(
        _dots, target, _vec_arbiter.data(), num_basis + 1, this->_subspace.size());

      // read the results back into the host output (triggers a device->host copy if needed)
      Memory::TypedView<DataType> dv(_dots.view(Memory::Location::main, Memory::Access::read));
      const DataType* dp = dv.get_r();
      for(Index k(0); k < num_basis; ++k)
        dot_out[k] = dp[k];
      self_dot = dp[num_basis];
    }

    /**
     * \brief Recovers the orthogonalized norm from the algebraic identity with a norm2 fallback
     */
    DataType _finalize_norm(const DataType nrm_sqr, const VectorType& target) const
    {
      // trust the algebraic norm only when it is safely positive; otherwise fall back to an explicit
      // norm2 (accurate, and guarantees a positive result so the normalization cannot produce a NaN)
      if(nrm_sqr <= Math::Limits<DataType>::epsilon() * DataType(100))
        return target.norm2();
      return Math::sqrt(nrm_sqr);
    }

  public:
    explicit ClassicGramSchmidtOrthogonalizer() = default;

    ClassicGramSchmidtOrthogonalizer(const Vector_& vec_ref, Index krylov_dimension) :
     BaseClass(vec_ref, krylov_dimension)
    {
      // one slot per basis vector plus one for the appended self-dot <w,w>; the largest reduction
      // projects against all krylov_dimension basis vectors, so krylov_dimension+1 slots suffice
      _dots = Memory::Arbiter((this->_krylov_dimension + 1) * sizeof(DataType));
      // pass pointers to the persistent krylov basis vectors; gather_arbiter records the addresses
      // of their (nested) element arbiters, which must remain valid for the solver's lifetime
      std::vector<const VectorType*> vec_v_ptrs;
      vec_v_ptrs.reserve(this->_subspace.size());
      for(const auto& v : this->_subspace)
        vec_v_ptrs.push_back(&v);
      Intern::BatchedHelper<VectorType>::gather_arbiter(_vec_arbiter, vec_v_ptrs.data(), this->_subspace.size());
    }

    ClassicGramSchmidtOrthogonalizer(const ClassicGramSchmidtOrthogonalizer&) = delete;
    ClassicGramSchmidtOrthogonalizer& operator=(const ClassicGramSchmidtOrthogonalizer&) = delete;
    ClassicGramSchmidtOrthogonalizer(ClassicGramSchmidtOrthogonalizer&&) = default;
    ClassicGramSchmidtOrthogonalizer& operator=(ClassicGramSchmidtOrthogonalizer&&) = default;

    virtual ~ClassicGramSchmidtOrthogonalizer()=default;

    DataType orthogonalize(DataType* dot_products, const Index target_index) override
    {
      XASSERT(target_index <= this->_krylov_dimension);
      // orthogonalize against the target_index leading basis vectors v_0,...,v_{target_index-1}
      const Index num_basis = target_index;
      VectorType& target = this->_subspace[target_index];

      // a single batched reduction yields the orthogonalization dots h[0..num_basis) and the self-dot
      DataType self_dot(0);
      _compute_dots(target, num_basis, dot_products, self_dot);

      DataType sum_sqr(0);
      for(Index k(0); k < num_basis; ++k)
        sum_sqr += Math::sqr(dot_products[k]);

      // subtract the projection in a single batched cumulative axpy: target -= sum_k h[k] * v[k]
      this->correct_vector(target, dot_products, num_basis, DataType(-1));

      // recover the post-orthogonalization norm algebraically (with a norm2 fallback)
      return _finalize_norm(self_dot - sum_sqr, target);
    }

    /**
     * \brief Batched cumulative axpy: target += alpha * sum_k local_scales[k] * v[k]
     */
    void correct_vector(VectorType& target, const DataType* local_scales, const Index num_axpy, const DataType alpha) const override
    {
      XASSERT(num_axpy <= this->_krylov_dimension);
      Intern::BatchedHelper<VectorType>::template cumulative_axpy<DataType>(
        target, _vec_arbiter.data(), local_scales, num_axpy, alpha, this->_subspace.size());
    }
  };

  /**
   * \brief Classical Gram-Schmidt with a (conditional) reorthogonalization sweep (CGS2)
   *
   * This extends the batched classical Gram-Schmidt procedure of ClassicGramSchmidtOrthogonalizer by
   * a corrective second orthogonalization sweep, recovering a numerical stability close to modified
   * Gram-Schmidt while keeping the dot products batched (at most two reductions per step instead of
   * one reduction per basis vector).
   *
   * The corrective sweep is either always performed (see set_force_reortho()) or triggered
   * automatically by the classic DGKS "twice is enough" criterion: the sweep is redone whenever the
   * orthogonalized vector lost too much of its length, i.e. \f$\|v_\text{new}\| < \tau \|w\|\f$
   * (equivalently \f$\text{nrm\_sqr} < \tau^2 \langle w,w\rangle\f$). With the DGKS threshold the
   * algebraic norm is only trusted when the vector kept a healthy fraction of its length, and any
   * heavier cancellation pays exactly one corrective reduction.
   */
  template<typename Vector_>
  class CGS2Orthogonalizer : public ClassicGramSchmidtOrthogonalizer<Vector_>
  {
  public:
    typedef ClassicGramSchmidtOrthogonalizer<Vector_> BaseClass;
    typedef typename BaseClass::VectorType VectorType;
    typedef typename BaseClass::DataType DataType;

  protected:
    /// scratch for the corrective sweep's dot products, kept separate from the accumulated h[k]
    std::vector<DataType> _corr;
    /// norm-ratio threshold below which a corrective sweep is triggered; the classic DGKS value is
    /// 1/sqrt(2), but in practice a smaller value already suffices (matching the POC default)
    DataType _reortho_tol = DataType(0.3);
    /// if set, always perform the corrective reorthogonalization sweep regardless of the criterion
    bool _force_reortho = false;

  public:
    explicit CGS2Orthogonalizer() = default;

    CGS2Orthogonalizer(const Vector_& vec_ref, Index krylov_dimension) :
     BaseClass(vec_ref, krylov_dimension)
    {
      _corr.resize(krylov_dimension);
    }

    CGS2Orthogonalizer(const CGS2Orthogonalizer&) = delete;
    CGS2Orthogonalizer& operator=(const CGS2Orthogonalizer&) = delete;
    CGS2Orthogonalizer(CGS2Orthogonalizer&&) = default;
    CGS2Orthogonalizer& operator=(CGS2Orthogonalizer&&) = default;

    virtual ~CGS2Orthogonalizer()=default;

    /// Sets the DGKS norm-ratio threshold that triggers a corrective reorthogonalization sweep
    void set_reortho_tol(DataType tol)
    {
      XASSERT(tol >= DataType(0));
      _reortho_tol = tol;
    }

    /// Forces (or unforces) a corrective reorthogonalization sweep in every step
    void set_force_reortho(bool force)
    {
      _force_reortho = force;
    }

    DataType orthogonalize(DataType* dot_products, const Index target_index) override
    {
      XASSERT(target_index <= this->_krylov_dimension);
      const Index num_basis = target_index;
      VectorType& target = this->_subspace[target_index];

      // first classical Gram-Schmidt sweep (batched dots + self-dot, then batched subtraction)
      DataType self_dot(0);
      this->_compute_dots(target, num_basis, dot_products, self_dot);
      DataType sum_sqr(0);
      for(Index k(0); k < num_basis; ++k)
        sum_sqr += Math::sqr(dot_products[k]);
      this->correct_vector(target, dot_products, num_basis, DataType(-1));
      DataType nrm_sqr = self_dot - sum_sqr;

      // corrective reorthogonalization (CGS2), forced or when heavy cancellation is detected; the
      // second batched reduction supplies both the corrections and the updated self-dot, so the
      // corrected norm is again obtained without an extra reduction
      if(_force_reortho || (nrm_sqr < Math::sqr(_reortho_tol) * self_dot))
      {
        DataType self_dot2(0);
        this->_compute_dots(target, num_basis, _corr.data(), self_dot2);
        DataType sum_sqr2(0);
        for(Index k(0); k < num_basis; ++k)
        {
          dot_products[k] += _corr[k];
          sum_sqr2 += Math::sqr(_corr[k]);
        }
        this->correct_vector(target, _corr.data(), num_basis, DataType(-1));
        nrm_sqr = self_dot2 - sum_sqr2;
      }

      return this->_finalize_norm(nrm_sqr, target);
    }
  };




  namespace Intern
  {

    #ifdef DOXY
    template<typename VectorType_>
    class BatchedHelper
    {
    public:
      /**
      * \brief Calculates the dot products of target with all members of vectors in a batched manner
      *
      * \param[in,out] local_dots The Memory::Arbiter holding the dot-product results (values are added upon); the
      *            batched kernels acquire the host/device view they need from it, so it must already contain the
      *            zeroed or partially accumulated values.
      * \param[in] target The vector all dot products shall be calculated with
      * \param[in] vectors An array of memory arbiter objects for which the dot products shall be calculated
      * \param[in] size The number of meta elements in vectors
      *
      * \param[in] stride The distance (in arbiters) between two consecutive meta-components of the same vector inside
      *            the flat \p vectors array, i.e. the per-component block size used by gather_arbiter (the total number
      *            of gathered vectors). Used to advance to the next tuple component; it is \b not equal to \p size.
      *
      * \note It is assumed that vectors hold the unpacked memory arbiter of the underlying classes, i.e. filled by the
      *       get_arbiter function.
      *
      * \note It is assumed that local dots is allocated to the correct size
      *
      */
      template<typename DT_>
      static void batched_dot(Memory::Arbiter& local_dots, const VectorType_& target, Memory::Arbiter const * const * vectors, Index size, Index stride);

      /**
      * \brief Calculates the commulative axpy onto target with all members of vectors in a batched manner
      *
      * \param[in,out] target The vector all axpy shall be added on
      * \param[in] vectors An array of memory arbiter objects for which the axpys shall be executed
      * \param[in] local_scales An array of scaling factors to be used.
      * \param[in] size The number of elements in vectors
      * \param[in] alpha An additional global scaling factor
      *
      * \note It is assumed that vectors hold the unpacked memory arbiter of the underlying classes, i.e. filled by the
      *       appropriate function
      *
      * \note It is assumed that local dots is allocated to the correct size
      *
      */
      template<typename DT_>
      static void cumulative_axpy(VectorType_& target, Memory::Arbiter const * const * vectors, DT_ const * const local_scales, Index size, DT_ alpha, Index stride);

      template<typename DT_>
      static void gather_arbiter(std::vector<Memory::Arbiter const*>& arbiter, VectorType_ const * const * vectors, Index size);
    };
    #endif

    template<typename LocalVector_, typename Mirror_>
    class BatchedHelper<Global::Vector<LocalVector_, Mirror_>>
    {
    public:
      typedef Global::Vector<LocalVector_, Mirror_> VectorType;

      template<typename DT_>
      static void batched_dot(Memory::Arbiter& local_dots, const VectorType& target, Memory::Arbiter const * const * vectors, Index size, Index stride)
      {
        // create vector of local vectors
        if(!target.get_comm() || target.get_comm()->size() == 1)
        {
          // call local batched dot
          BatchedHelper<LocalVector_>::template batched_dot<DT_>(local_dots, target.local(), vectors, size, stride);
        }
        else
        {
          BatchedHelper<LocalVector_>::template batched_triple_dot<DT_>(local_dots, target.local(), vectors, target.get_gate()->get_freqs(), size, stride);

          {
            // and now synchronize dots on the host: acquire the main-memory view (triggers a device->host
            // copy if the local dots were just computed on the device) and all-reduce in place
            Memory::TypedView<DT_> dots_view(local_dots.view(Memory::Location::main, Memory::Access::read_write));
            DT_* dots_ptr = dots_view.get_w();
            target.get_comm()->allreduce(dots_ptr, dots_ptr, size, Dist::op_sum);
          }
          // TODO: allreduce should decide itself where the communication should take place... can the view simply be on the last most update point?
          //       Simply take the backend?

        }
      }

      template<typename DT_>
      static void cumulative_axpy(VectorType& target, Memory::Arbiter const * const * vectors, DT_ const * const local_scales, Index size, DT_ alpha, Index stride)
      {
        // call local batched dot
        BatchedHelper<LocalVector_>::cumulative_axpy(target.local(), vectors, local_scales, size, alpha, stride);
      }

      static void gather_arbiter(std::vector<Memory::Arbiter const*>& arbiter, VectorType const * const * vectors, Index size)
      {
        // collect pointers into the persistent local vectors; no cloning, so the arbiter addresses
        // we gather remain valid for the lifetime of the originals
        std::vector<const LocalVector_*> locals;
        locals.reserve(size);
        for(Index k = 0; k < size; ++k)
        {
          locals.push_back(&vectors[k]->local());
        }
        BatchedHelper<LocalVector_>::gather_arbiter(arbiter, locals.data(), size);
      }
    };

    template<typename First_, typename... Rest_>
    class BatchedHelper<LAFEM::TupleVector<First_, Rest_...>>
    {
    public:
      typedef LAFEM::TupleVector<First_, Rest_...> VectorType;

      template<typename DT_>
      static void batched_dot(Memory::Arbiter& local_dots, const VectorType& target, Memory::Arbiter const * const * vectors, Index size, Index stride)
      {
        {
          BatchedHelper<First_>::template batched_dot<DT_>(local_dots, target.first(), vectors, size, stride);
        }
        {
          BatchedHelper<LAFEM::TupleVector<Rest_...>>::template batched_dot<DT_>(local_dots, target.rest(), vectors+stride, size, stride);
        }
      }

      template<typename DT_>
      static void batched_triple_dot(Memory::Arbiter& local_dots, const VectorType& target, Memory::Arbiter const * const * vectors, const VectorType& freqs, Index size, Index stride)
      {
        {
          BatchedHelper<First_>::template batched_triple_dot<DT_>(local_dots, target.first(), vectors, freqs.first(), size, stride);
        }
        {
          BatchedHelper<LAFEM::TupleVector<Rest_...>>::template batched_triple_dot<DT_>(local_dots, target.rest(), vectors+stride, freqs.rest(), size, stride);
        }
      }

      template<typename DT_>
      static void cumulative_axpy(VectorType& target, Memory::Arbiter const * const * vectors, DT_ const * const local_scales, Index size, DT_ alpha, Index stride)
      {
        BatchedHelper<First_>::cumulative_axpy(target.first(), vectors, local_scales, size, alpha, stride);
        BatchedHelper<LAFEM::TupleVector<Rest_...>>::cumulative_axpy(target.rest(), vectors+stride, local_scales, size, alpha, stride);
      }

      static void gather_arbiter(std::vector<Memory::Arbiter const*>& arbiter, VectorType const * const * vectors, Index size)
      {
        {
          std::vector<const First_*> firsts;
          firsts.reserve(size);
          for(Index k = 0; k < size; ++k)
          {
            firsts.push_back(&vectors[k]->first());
          }
          BatchedHelper<First_>::gather_arbiter(arbiter, firsts.data(), size);
        }
        {
          std::vector<const LAFEM::TupleVector<Rest_...>*> rests;
          rests.reserve(size);
          for(Index k = 0; k < size; ++k)
          {
            rests.push_back(&vectors[k]->rest());
          }
          BatchedHelper<LAFEM::TupleVector<Rest_...>>::gather_arbiter(arbiter, rests.data(), size);
        }
      }
    };

    template<typename First_>
    class BatchedHelper<LAFEM::TupleVector<First_>>
    {
    public:
      typedef LAFEM::TupleVector<First_> VectorType;

      template<typename DT_>
      static void batched_dot(Memory::Arbiter& local_dots, const VectorType& target, Memory::Arbiter const * const * vectors, Index size, Index stride)
      {
        {
          BatchedHelper<First_>::template batched_dot<DT_>(local_dots, target.first(), vectors, size, stride);
        }
      }

      template<typename DT_>
      static void batched_triple_dot(Memory::Arbiter& local_dots, const VectorType& target, Memory::Arbiter const * const * vectors, const VectorType& freqs, Index size, Index stride)
      {
        {
          BatchedHelper<First_>::template batched_triple_dot<DT_>(local_dots, target.first(), vectors, freqs.first(), size, stride);
        }
      }

      template<typename DT_>
      static void cumulative_axpy(VectorType& target, Memory::Arbiter const * const * vectors, DT_ const * const local_scales, Index size, DT_ alpha, Index stride)
      {
        BatchedHelper<First_>::cumulative_axpy(target.first(), vectors, local_scales, size, alpha, stride);
      }

      static void gather_arbiter(std::vector<Memory::Arbiter const*>& arbiter, VectorType const * const * vectors, Index size)
      {
        {
          std::vector<const First_*> firsts;
          firsts.reserve(size);
          for(Index k = 0; k < size; ++k)
          {
            firsts.push_back(&vectors[k]->first());
          }
          BatchedHelper<First_>::gather_arbiter(arbiter, firsts.data(), size);
        }
      }
    };

    template<typename DTT_, typename ITT_>
    class BatchedHelper<LAFEM::DenseVector<DTT_, ITT_>>
    {
    public:
      typedef LAFEM::DenseVector<DTT_, ITT_> VectorType;

      template<typename DT_>
      static void batched_dot(Memory::Arbiter& local_dots, const VectorType& target, Memory::Arbiter const * const * vectors, Index size, Index /*stride*/)
      {
        if(size == 0u)
          return;
        LAFEM::Arch::BatchedDotProduct::exec<DT_, DTT_>(local_dots, target.elements_arbiter(), target.size(), vectors, size);
      }

      template<typename DT_>
      static void batched_triple_dot(Memory::Arbiter& local_dots, const VectorType& target, Memory::Arbiter const * const * vectors, const VectorType& freqs, Index size, Index /*stride*/)
      {
        if(size == 0u)
          return;
        LAFEM::Arch::BatchedTripleDotProduct::exec<DT_, DTT_>(local_dots, target.elements_arbiter(), target.size(), vectors, freqs.elements_arbiter(), size);
      }

      template<typename DT_>
      static void cumulative_axpy(VectorType& target, Memory::Arbiter const * const * vectors, DT_ const * const local_scales, Index size, DT_ alpha, Index /*stride*/)
      {
        if(size == 0u)
          return;
        LAFEM::Arch::CumulativeAxpy::exec<DT_>(target.elements_arbiter(), target.size(), vectors, local_scales, size, alpha);
      }

      static void gather_arbiter(std::vector<Memory::Arbiter const*>& arbiter, VectorType const * const * vectors, Index size)
      {
        for(Index k = 0; k < size; ++k)
        {
          arbiter.emplace_back(&vectors[k]->elements_arbiter());
        }
      }
    };

    template<typename DTT_, typename ITT_, int dim_>
    class BatchedHelper<LAFEM::DenseVectorBlocked<DTT_, ITT_, dim_>>
    {
    public:
      typedef LAFEM::DenseVectorBlocked<DTT_, ITT_, dim_> VectorType;

      template<typename DT_>
      static void batched_dot(Memory::Arbiter& local_dots, const VectorType& target, Memory::Arbiter const * const * vectors, Index size, Index /*stride*/)
      {
        if(size == 0u)
          return;
        // pass the element arbiters down; the kernel acquires the (pod-typed) views for its backend itself
        LAFEM::Arch::BatchedDotProduct::exec<DT_, DTT_>(local_dots, target.elements_arbiter(), target.size_raw(), vectors, size);
      }

      template<typename DT_>
      static void batched_triple_dot(Memory::Arbiter& local_dots, const VectorType& target, Memory::Arbiter const * const * vectors, const VectorType& freqs, Index size, Index /*stride*/)
      {
        if(size == 0u)
          return;
        LAFEM::Arch::BatchedTripleDotProduct::exec<DT_, DTT_>(local_dots, target.elements_arbiter(), target.size_raw(),
                                          vectors, freqs.elements_arbiter(), size);
      }

      template<typename DT_>
      static void cumulative_axpy(VectorType& target, Memory::Arbiter const * const * vectors, DT_ const * const local_scales, Index size, DT_ alpha, Index /*stride*/)
      {
        if(size == 0u)
          return;
        LAFEM::Arch::CumulativeAxpy::exec<DT_>(target.elements_arbiter(), target.size_raw(), vectors, local_scales, size, alpha);
      }

      static void gather_arbiter(std::vector<Memory::Arbiter const*>& arbiter, VectorType const * const * vectors, Index size)
      {
        for(Index k = 0; k < size; ++k)
        {
          arbiter.emplace_back(&vectors[k]->elements_arbiter());
        }
      }
    };

  }
}
