#pragma once
#ifndef KERNEL_SOLVER_IDRS_HPP
#define KERNEL_SOLVER_IDRS_HPP 1

// includes, FEAT
#include <kernel/solver/iterative.hpp>
#include <vector>
#include <kernel/global/vector.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/random.hpp>
#include <kernel/util/dist.hpp>
namespace FEAT
{
  namespace Solver
  {

    /**
     * \brief IDR(s) solver implementation
     *
     * This class implements the induced dimension reduction solver \cite Sonneveld09.
     *
     * \see
     * P. Sonneveld and M.B. van Gijzen: IDR(s): A family of simple and fast algorithms for solving large nonsymmetric linear equations;
     * SIAM Journal on Scientific Computing, Volume 31 Issue 2, pp. 10351062-469, 2008
     *
     * \tparam Matrix_
     * The matrix class to be used by the solver.
     *
     * \tparam Filter_
     * The filter class to be used by the solver.
     *
     * \author Jonas Duennebacke
     */
    template<
      typename Matrix_,
      typename Filter_>
    class IDRS :
      public PreconditionedIterativeSolver<typename Matrix_::VectorTypeR>
    {
    public:
      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;
      typedef typename MatrixType::MemType MemType;
      typedef typename MatrixType::IndexType IndexType;
      typedef PreconditionedIterativeSolver<VectorType> BaseClass;

      typedef FEAT::LAFEM::DenseMatrix<MemType, DataType, IndexType> DMatrix;
      typedef FEAT::LAFEM::DenseVector<MemType, DataType, IndexType> DVector;
      typedef SolverBase<VectorType> PrecondType;

    protected:
      /// the matrix for the solver
      const MatrixType& _system_matrix;
      /// the filter for the solver
      const FilterType& _system_filter;
      /// krylov dimension (s)
      Index _krylov_dim;
      /// auxillary vectors from the algorithm
      VectorType _vec_q, _vec_r, _vec_v, _vec_t;
      /// real (unpreconditioned) residual
      VectorType _vec_res;
      /// n x s "matrices"
      std::vector<VectorType> _vec_P, _vec_dR, _vec_dX;
      /// dense s x s - matrices
      DMatrix _dmat_M, _dmat_Minv;
      /// local ("small") dense vectors
      DVector _dvec_m, _dvec_dm, _dvec_c;
      bool _shadow_space_setup, _random;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * A reference to the system matrix.
       *
       * \param[in] filter
       * A reference to the system filter.
       *
       * \param[in] krylov_dim
       * The maximum Krylov subspace dimension (s). Must be > 0.
       *
       * \param[in] precond
       * A pointer to the preconditioner. May be \c nullptr.
       */
      explicit IDRS(const MatrixType& matrix, const FilterType& filter, Index krylov_dim,
        std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("IDRS(" + stringify(krylov_dim) + ")", precond),
        _system_matrix(matrix),
        _system_filter(filter),
        _krylov_dim(krylov_dim),
        _shadow_space_setup(false),
        _random(true)
      {
      }

      explicit IDRS(const String& section_name, PropertyMap* section,
      const MatrixType& matrix, const FilterType& filter, std::shared_ptr<PrecondType> precond = nullptr) :
        BaseClass("IDRS", section_name, section, precond),
        _system_matrix(matrix),
        _system_filter(filter),
        _shadow_space_setup(false),
        _random(true)
      {
        // Check if we have set _krylov_dim
        auto krylov_dim_p = section->query("krylov_dim");
        if(krylov_dim_p.second)
        {
          set_krylov_dim(Index(std::stoul(krylov_dim_p.first)));
          if  (this->_plot_name == "IDRS")
          {
            this->set_plot_name("IDR("+stringify(_krylov_dim)+")");
          }
        }
        else
        {
          throw InternalError(__func__,__FILE__,__LINE__,
          "IDRS config section is missing the mandatory krylov_dim!");
        }
      }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~IDRS()
      {
      }

      /// \copydoc BaseClass::name()
      virtual String name() const override
      {
        return "IDRS";
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        _dmat_M    = DMatrix(_krylov_dim, _krylov_dim);
        _dvec_m    = DVector(_krylov_dim);
        _dvec_dm   = DVector(_krylov_dim);
        _dvec_c    = DVector(_krylov_dim);

        _vec_q   = this->_system_matrix.create_vector_r();
        _vec_r   = _vec_q.clone(LAFEM::CloneMode::Layout);
        _vec_res = _vec_q.clone(LAFEM::CloneMode::Layout);
        _vec_t   = _vec_q.clone(LAFEM::CloneMode::Layout);
        _vec_v   = _vec_q.clone(LAFEM::CloneMode::Layout);
        _vec_dX.reserve(_krylov_dim);
        _vec_dR.reserve(_krylov_dim);
        _vec_P.reserve(_krylov_dim);

        for(Index i(0); i < _krylov_dim; ++i)
        {
          _vec_dX.push_back(this->_vec_q.clone(LAFEM::CloneMode::Layout));
          _vec_dR.push_back(this->_vec_q.clone(LAFEM::CloneMode::Layout));
          _vec_P.push_back( this->_vec_q.clone(LAFEM::CloneMode::Layout));
        }
      }

      virtual void done_symbolic() override
      {
        _vec_q.clear();
        _vec_r.clear();
        _vec_res.clear();
        _vec_t.clear();
        _vec_v.clear();

        _vec_dX.clear();
        _vec_dR.clear();
        _vec_P.clear();

        BaseClass::done_symbolic();
      }

      /**
       * \brief Sets the inner Krylov space dimension
       *
       * \param[in] krylov_dim
       * The s in IDR(s)
       *
       */
      virtual void set_krylov_dim(Index krylov_dim)
      {
        XASSERT(krylov_dim > Index(0));
        _krylov_dim = krylov_dim;
      }

      /// \copydoc IterativeSolver::apply()
      virtual Status apply(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // save input rhs vector as initial defect
        this->_vec_res.copy(vec_rhs);
        //this->_system_filter.filter_def(this->_vec_v.at(0));

        // clear solution vector
        vec_sol.format();

        // apply
        Status st(_apply_intern(vec_sol, vec_rhs));
        this->plot_summary(st);
        return st;
      }

      /// \copydoc SolverBase::correct()
      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute initial defect
        this->_system_matrix.apply(this->_vec_res, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.filter_def(this->_vec_res);

        // apply
        Status st(_apply_intern(vec_sol, vec_rhs));
        this->plot_summary(st);
        return st;
      }

      /**
      * \brief Reset the shadow space before the next system will be solved
      *
      * \param[in] random
      * boolean whether the new shadow space should be random or not
      */
      void reset_shadow_space( bool random = true)
      {
        _random = random;
        _shadow_space_setup = false;
      }



    protected:

      void _set_shadow_space()
      {
        //select set of random vectors

        Random::SeedType seed;
        if (_random)
        {
          seed = Random::get_seed();
        }
        else
        {
          auto comm = Dist::Comm::world();
          seed = Random::SeedType(comm.rank());
        }
        Random rng(seed);

        for (Index l(0); l<_krylov_dim; ++l)
          _vec_P.at(l).format(rng, DataType(0), DataType(1));

        //orthonormalize vector set P (modified Gram-Schmidt)
        for (Index j(0); j<_krylov_dim; ++j)
        {
          for (Index i(0); i<j; ++i)
          {
            _vec_P.at(j).axpy(_vec_P.at(i), _vec_P.at(j), -_vec_P.at(j).dot(_vec_P.at(i)));
          }
          _vec_P.at(j).scale(_vec_P.at(j), DataType(1)/_vec_P.at(j).norm2());
        }
        _shadow_space_setup = true;
      }




      virtual Status _apply_intern(VectorType& vec_sol, const VectorType& DOXY(vec_rhs))
      {
        IterationStats pre_iter(*this);
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));
        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);

        // compute initial defect
        Status status = this->_set_initial_defect(this->_vec_res, vec_sol);
        _vec_t.copy(_vec_res);
        if(!this->_apply_precond(this->_vec_r, this->_vec_t, filter))
        {
          pre_iter.destroy();
          Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
          return Status::aborted;
        }
        //select random vector set (shadow space)
        if (!_shadow_space_setup)
          _set_shadow_space();

        pre_iter.destroy();

        //produce start vectors
        IterationStats first_iter(*this);
        DataType om(0);
        for (Index k(0); k < _krylov_dim; ++k)
        {
          matrix.apply(_vec_v, _vec_r);
          filter.filter_def(_vec_v);
          //save unpreconditioned residual update
          _vec_t.copy(_vec_v);
          if(!this->_apply_precond(this->_vec_v, this->_vec_t, filter))
          {
            first_iter.destroy();
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
            return Status::aborted;
          }
          om = _vec_v.dot(_vec_r) / _vec_v.dot(_vec_v);
          _vec_dX.at(k).scale(_vec_r,  om);
          _vec_dR.at(k).scale(_vec_v, -om);
          vec_sol.axpy(vec_sol, _vec_dX.at(k));
          _vec_r.axpy(_vec_r, _vec_dR.at(k));
          _vec_res.axpy(_vec_t, _vec_res, -om);

          // push our new defect
          status = this->_set_new_defect(_vec_res, vec_sol);

          //check for early convergence
          if (status != Status::progress)
          {
            Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
            return status;
          }
          // update k-th column of M
          for (Index l(0); l <  _krylov_dim; ++l)
            _dmat_M(l,k, _vec_P.at(l).dot(_vec_dR.at(k)) );

        }

        Index oldest(0);
        //calculate m
        for (Index l(0); l<  _krylov_dim; ++l)
          _dvec_m(l, _vec_P.at(l).dot(_vec_r) );

        first_iter.destroy();
        // outer loop
        while(status == Status::progress)
        {
          IterationStats stat(*this);

          //inner loop
          for (Index k(0); k <= _krylov_dim; ++k)
          {
            //solve c=M\m
            _dmat_Minv = _dmat_M.inverse();
            _dmat_Minv.apply(_dvec_c, _dvec_m);

            _vec_q.scale(_vec_dR.at(0), -_dvec_c(0) );
            for (Index l(1); l < _krylov_dim; ++l)
              _vec_q.axpy(_vec_dR.at(l), _vec_q, -_dvec_c(l));

            _vec_v.axpy(_vec_r, _vec_q);

            if (k == 0)
            {
              matrix.apply(_vec_dR.at(oldest), _vec_v);
              filter.filter_def(_vec_dR.at(oldest));
              if(!this->_apply_precond(this->_vec_t, this->_vec_dR.at(oldest), filter))
              {
                stat.destroy();
                Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
                return Status::aborted;
              }
              om = _vec_v.dot(_vec_t) / _vec_t.dot(_vec_t);
              _vec_dR.at(oldest).axpy(_vec_t, _vec_q, -om);

              //dX(oldest)=-dX*c+om*v
              _vec_dX.at(oldest).scale(_vec_dX.at(oldest), -_dvec_c(oldest));
              for (Index l(0); l < oldest; ++l)
                _vec_dX.at(oldest).axpy(_vec_dX.at(l), _vec_dX.at(oldest), -_dvec_c(l));
              for (Index l(oldest+1); l < _krylov_dim; ++l)
                _vec_dX.at(oldest).axpy(_vec_dX.at(l), _vec_dX.at(oldest), -_dvec_c(l));
              _vec_dX.at(oldest).axpy(_vec_v, _vec_dX.at(oldest), om);

              matrix.apply(_vec_t, _vec_dX.at(oldest));
              filter.filter_def(_vec_t);
            }
            else
            {
              //dX(oldest)=-dX*c+om*v
              _vec_dX.at(oldest).scale(_vec_dX.at(oldest), -_dvec_c(oldest));
              for (Index l(0); l < oldest; ++l)
                _vec_dX.at(oldest).axpy(_vec_dX.at(l), _vec_dX.at(oldest), -_dvec_c(l));
              for (Index l(oldest+1); l < _krylov_dim; ++l)
                _vec_dX.at(oldest).axpy(_vec_dX.at(l), _vec_dX.at(oldest), -_dvec_c(l));
              _vec_dX.at(oldest).axpy(_vec_v, _vec_dX.at(oldest), om);

              matrix.apply( _vec_dR.at(oldest), _vec_dX.at(oldest));
              filter.filter_def(_vec_dR.at(oldest));
              _vec_t.copy(_vec_dR.at(oldest));
              if(!this->_apply_precond(this->_vec_dR.at(oldest), this->_vec_t, filter))
              {
                stat.destroy();
                Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), Status::aborted, this->get_num_iter()));
                return Status::aborted;
              }
              _vec_dR.at(oldest).scale(_vec_dR.at(oldest), DataType(-1));
            }


            _vec_r.axpy( _vec_dR.at(oldest), _vec_r );
            vec_sol.axpy(_vec_dX.at(oldest), vec_sol);
            _vec_res.axpy(_vec_t, _vec_res, -1);

            status = this->_set_new_defect(_vec_res, vec_sol);
            if (status != Status::progress)
              break;

            for (Index l(0); l <  _krylov_dim; ++l)
            {
              _dvec_dm(l, _vec_P.at(l).dot(_vec_dR.at(oldest)) );
              _dmat_M( l, oldest, _dvec_dm(l));
            }
            _dvec_m.axpy(_dvec_dm, _dvec_m);
            oldest = (oldest+1) % _krylov_dim;
          } //end inner loop

        } //end outer loop
        // finished
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
        return status;
      } //_apply_intern(...)
    }; // class IDRS<...>

    /**
     * \brief Creates a new IDR(s) solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] krylov_dim
     * The maximum Krylov subspace dimension. Must be > 0.
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new IDRS object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<IDRS<Matrix_, Filter_>> new_idrs(
      const Matrix_& matrix, const Filter_& filter, Index krylov_dim)
    {
      return std::make_shared<IDRS<Matrix_, Filter_>>(matrix, filter, krylov_dim);
    }
    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<IDRS<Matrix_, Filter_>> new_idrs(
      const Matrix_& matrix, const Filter_& filter, Index krylov_dim,
      std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<IDRS<Matrix_, Filter_>>(matrix, filter, krylov_dim, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<IDRS<Matrix_, Filter_>> new_idrs(
      const Matrix_& matrix, const Filter_& filter, Index krylov_dim,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<IDRS<Matrix_, Filter_>>(matrix, filter, krylov_dim, precond);
    }
#endif

    /**
     * \brief Creates a new IDR(s) solver object using a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] precond
     * The preconditioner. May be \c nullptr.
     *
     * \returns
     * A shared pointer to a new IDRS object.
     */
     /// \compilerhack GCC < 4.9 fails to deduct shared_ptr
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 40900)
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<IDRS<Matrix_, Filter_>> new_idrs(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<IDRS<Matrix_, Filter_>>(section_name, section, matrix, filter);
    }

    template<typename Matrix_, typename Filter_, typename Precond_>
    inline std::shared_ptr<IDRS<Matrix_, Filter_>> new_idrs(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter, std::shared_ptr<Precond_> precond)
    {
      return std::make_shared<IDRS<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#else
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<IDRS<Matrix_, Filter_>> new_idrs(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter,
      std::shared_ptr<SolverBase<typename Matrix_::VectorTypeL>> precond = nullptr)
    {
      return std::make_shared<IDRS<Matrix_, Filter_>>(section_name, section, matrix, filter, precond);
    }
#endif
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_IDRS_HPP
