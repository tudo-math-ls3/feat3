#pragma once
#ifndef FEAST_KERNEL_SOLVER_LINESEARCH
#define FEAST_KERNEL_SOLVER_LINESEARCH 1
#include <kernel/base_header.hpp>
#include <kernel/solver/iterative.hpp>

#include <deque>

namespace FEAST
{
  namespace Solver
  {
    template<typename Operator_, typename Filter_>
    class NRLinesearch : public IterativeSolver<typename Operator_::VectorTypeR>
    {
      public:
        typedef typename Operator_::VectorTypeR VectorType;
        typedef typename Operator_::DataType DataType;
        typedef Filter_ FilterType;
        typedef IterativeSolver<typename Operator_::VectorTypeR> BaseClass;

      protected:
        Operator_& _op;
        const Filter_& _filter;

        /// Gradient vector
        VectorType _vec_grad;
        /// temporary vector
        VectorType _vec_tmp;

        DataType _alpha;
        DataType _sum_alpha;
        DataType _norm_dir;
        DataType _norm_sol;

      public:
        std::deque<VectorType>* iterates;

      public:
        explicit NRLinesearch(Operator_& op_, const Filter_& filter_,
        bool keep_iterates = false) :
          BaseClass("NR-LS"),
          _op(op_),
          _filter(filter_),
          _alpha(DataType(0)),
          _sum_alpha(DataType(0)),
          _norm_dir(DataType(0)),
          _norm_sol(DataType(0)),
          iterates(nullptr)
          {
            if(keep_iterates)
              iterates = new std::deque<VectorType>;

            this->_min_stag_iter = Index(2);
            this->_min_iter = Index(3);

            //this->set_plot(true);
          }

        virtual ~NRLinesearch()
        {
          if(iterates != nullptr)
            delete iterates;
        }

        virtual void init_symbolic() override
        {
          // Create temporary vector
          _vec_tmp = this->_op.create_vector_r();
          _vec_grad = this->_op.create_vector_r();
          BaseClass::init_symbolic();
        }

        virtual void done_symbolic() override
        {
          // Clear temporary vector
          _vec_tmp.clear();
          _vec_grad.clear();
          BaseClass::done_symbolic();
        }

        virtual String name() const override
        {
          return "Newton-Raphson-Linesearch";
        }

        virtual Status apply(VectorType& vec_cor, const VectorType& vec_dir) override
        {
          // clear solution vector
          vec_cor.format();
          this->_op.prepare(vec_cor);

          // apply
          return _apply_intern(vec_cor, vec_dir);
        }

        virtual Status correct(VectorType& vec_sol, const VectorType& vec_dir) override
        {
          this->_op.prepare(vec_sol);
          // apply
          Status st =_apply_intern(vec_sol, vec_dir);

          return st;
        }

        DataType get_rel_update()
        {
          return Math::abs(_sum_alpha)*_norm_dir/_norm_sol;
        }

      protected:
        virtual Status _apply_intern(VectorType& sol, const VectorType& dir)
        {
          Status status(Status::progress);
          this->_num_iter = Index(0);

          _norm_dir = dir.norm2();
          _norm_sol = sol.norm2();
          _sum_alpha = DataType(0);

          _op.prepare(sol);
          _op.compute_gradient(this->_vec_grad);
          _filter.filter_sol(this->_vec_grad);

          this->_def_init = Math::abs(dir.dot(this->_vec_grad));

          // start iterating
          while(status == Status::progress)
          {
            TimeStamp at;

            this->_op.prepare(sol);
            _op.compute_gradient(this->_vec_grad);
            _filter.filter_sol(this->_vec_grad);

            status = this->_set_new_defect(dir, this->_vec_grad);

            // Update solution: sol <- sol + _alpha*dir
            _sum_alpha += _alpha;
            sol.axpy(dir, sol, _alpha);

            if(status != Status::progress)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              //std::cout << "last alpha " << stringify_fp_sci(_alpha) << ", sum_alpha = " << stringify_fp_sci(_sum_alpha) << ", rel update " << stringify_fp_sci(get_rel_update()) << std::endl;
              return status;
            }

          }

          return Status::undefined;

        }

        Status _set_new_defect(const VectorType& dir, const VectorType& grad) override
        {
          // Increase iteration count
          ++this->_num_iter;

          // First let's see if we have to compute the defect at all
          bool calc_def = false;
          calc_def = calc_def || (this->_min_iter < this->_max_iter);
          calc_def = calc_def || this->_plot;
          calc_def = calc_def || (this->_min_stag_iter > Index(0));

          // Update defect
          this->_def_cur = Math::abs(grad.dot(dir));

          // Compute new _alpha = - grad.dot(dir) / dir.dot(Hess*dir)
          _op.apply_hess(this->_vec_tmp, dir);
          _filter.filter_sol(this->_vec_tmp);
          _alpha = -grad.dot(dir)/dir.dot(this->_vec_tmp);

          //std::cout << " alpha^2 * delta_0 = " << stringify_fp_sci(Math::sqr(_alpha)*this->_def_init) << std::endl;

          Statistics::add_solver_defect(this->_branch, this->_def_cur);

          // plot?
          if(this->_plot)
          {
            std::cout << this->_plot_name
            <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
            << " : " << stringify_fp_sci(this->_def_cur)
            << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
            << std::endl;
          }

          // ensure that the defect is neither NaN nor infinity
          if(!Math::isfinite(this->_def_cur))
            return Status::aborted;

          // is diverged?
          if(this->is_diverged())
            return Status::diverged;

          // minimum number of iterations performed?
          if(this->_num_iter < this->_min_iter)
            return Status::progress;

          // is converged?
          if(this->is_converged())
            return Status::success;

          // maximum number of iterations performed?
          if(this->_num_iter >= this->_max_iter)
            return Status::max_iter;

          // check for stagnation?
          if(this->_min_stag_iter > Index(0))
          {
            // did this iteration stagnate?
            if(Math::sqr(_alpha)*this->_def_init <= Math::eps<DataType>())
            {
              // increment stagnation count
              if(++this->_num_stag_iter >= this->_min_stag_iter)
              {
                //std::cout << this->_plot_name <<" stagnated with " <<
                //  stringify_fp_sci(Math::sqr(_alpha)*this->_def_init) << " <= " <<
                //  stringify_fp_sci(Math::sqr(this->_tol_rel)) << std::endl;
                return Status::stagnated;
              }
            }
            else
            {
              // this iteration did not stagnate
              this->_num_stag_iter = Index(0);
            }
          }

          // continue iterating
          return Status::progress;
        }

    };
    template<typename Operator_, typename Filter_>
    class SecantLinesearch : public IterativeSolver<typename Operator_::VectorTypeR>
    {
      public:
        typedef typename Operator_::VectorTypeR VectorType;
        typedef typename Operator_::DataType DataType;
        typedef Filter_ FilterType;
        typedef IterativeSolver<typename Operator_::VectorTypeR> BaseClass;

      protected:
        Operator_& _op;
        const Filter_& _filter;

        /// temporary vector
        VectorType _vec_tmp;

        DataType _sigma_0;
        DataType _alpha;
        DataType _sum_alpha;
        DataType _eta;
        DataType _norm_dir;
        DataType _norm_sol;

      public:
        std::deque<VectorType>* iterates;

      public:
        explicit SecantLinesearch(Operator_& op_, const Filter_& filter_, DataType initial_step = DataType(1e-2),
        bool keep_iterates = false) :
          BaseClass("S-LS"),
          _op(op_),
          _filter(filter_),
          _sigma_0(initial_step),
          _alpha(DataType(0)),
          _sum_alpha(DataType(0)),
          _eta(DataType(0)),
          _norm_dir(DataType(0)),
          _norm_sol(DataType(0)),
          iterates(nullptr)
          {
            if(keep_iterates)
              iterates = new std::deque<VectorType>;

            this->_min_stag_iter = Index(2);
            this->_min_iter = Index(5);

            this->set_plot(true);
          }

        virtual ~SecantLinesearch()
        {
          if(iterates != nullptr)
            delete iterates;
        }

        virtual void init_symbolic() override
        {
          // Create temporary vector
          _vec_tmp = this->_op.create_vector_r();
          BaseClass::init_symbolic();
        }

        virtual void done_symbolic() override
        {
          // Clear temporary vector
          _vec_tmp.clear();
          BaseClass::done_symbolic();
        }

        virtual String name() const override
        {
          return "SecantLinesearch";
        }

        virtual Status apply(VectorType& vec_cor, const VectorType& vec_dir) override
        {
          // clear solution vector
          vec_cor.format();
          this->_op.prepare(vec_cor);

          // apply
          return _apply_intern(vec_cor, vec_dir);
        }

        virtual Status correct(VectorType& vec_sol, const VectorType& vec_dir) override
        {
          this->_op.prepare(vec_sol);
          // apply
          Status st =_apply_intern(vec_sol, vec_dir);

          return st;
        }

        DataType get_rel_update()
        {
          return Math::abs(_sum_alpha)*_norm_dir/_norm_sol;
        }

      protected:
        virtual Status _apply_intern(VectorType& sol, const VectorType& dir)
        {
          // compute initial defect
          Status status(Status::progress);
          this->_num_iter = Index(0);

          //std::cout << "initial sol " << sol << std::endl;
          //std::cout << "        dir " << dir << std::endl;

          _alpha = -_sigma_0;
          _sum_alpha = DataType(0);
          _norm_dir = dir.norm2();
          _norm_sol = sol.norm2();

          // The first "other" point for the secant
          this->_vec_tmp.axpy(dir, sol, _sigma_0);
          //std::cout << "first secant point " << _vec_tmp << std::endl;
          _op.prepare(this->_vec_tmp);
          _op.compute_gradient(this->_vec_tmp);

          _eta = dir.dot(this->_vec_tmp);
          this->_def_init = Math::abs(_eta);

          // start iterating
          while(status == Status::progress)
          {
            TimeStamp at;

            this->_op.prepare(sol);
            _op.compute_gradient(this->_vec_tmp);
            _filter.filter_sol(this->_vec_tmp);
            //std::cout << "sol  = " << sol << std::endl;
            //std::cout << "grad = " << _vec_tmp << std::endl;

            status = this->_set_new_defect(dir, this->_vec_tmp);

            // Compute new alpha
            sol.axpy(dir, sol, _alpha);

            if(status != Status::progress)
            {
              TimeStamp bt;
              Statistics::add_solver_toe(this->_branch, bt.elapsed(at));
              //std::cout << "last alpha " << stringify_fp_sci(_alpha) << ", sum_alpha = " << stringify_fp_sci(_sum_alpha) << ", rel update " << stringify_fp_sci(get_rel_update()) << std::endl;
              return status;
            }

          }

          return Status::undefined;

        }

        Status _set_new_defect(const VectorType& dir, const VectorType& sol) override
        {
          // Increase iteration count
          ++this->_num_iter;

          DataType eta_prev = _eta;

          // Compute eta and alpha
          _eta = sol.dot(dir);
          this->_def_cur = Math::abs(_eta);

          if(Math::abs(_eta - eta_prev) < Math::eps<DataType>())///(Math::abs(_alpha*_eta)*_norm_dir))
          {
            //std::cout << "Update stagnation " << stringify_fp_sci(Math::abs(_eta - eta_prev)) << " < "
            //<< stringify_fp_sci(Math::eps<DataType>())<< std::endl;///(Math::abs(_alpha*_eta)*_norm_dir)) << std::endl;
            return Status::stagnated;
          }

          _alpha *= _eta/(eta_prev - _eta);
          _sum_alpha += _alpha;
          //std::cout << "eta_prev = " << stringify_fp_sci(eta_prev) << ", eta = " << stringify_fp_sci(_eta)
          //<< ",  alpha = " << stringify_fp_sci(_alpha);// << std::endl;
          //std::cout << " alpha^2 * delta_0 = " << stringify_fp_sci(Math::sqr(_alpha)*this->_def_init) << std::endl;

          Statistics::add_solver_defect(this->_branch, this->_def_cur);

          // plot?
          if(this->_plot)
          {
            std::cout << this->_plot_name
            <<  ": " << stringify(this->_num_iter).pad_front(this->_iter_digits)
            << " : " << stringify_fp_sci(this->_def_cur)
            << " / " << stringify_fp_sci(this->_def_cur / this->_def_init)
            << std::endl;
          }

          // ensure that the defect is neither NaN nor infinity
          if(!Math::isfinite(this->_def_cur))
            return Status::aborted;

          // is diverged?
          if(this->is_diverged())
            return Status::diverged;

          // minimum number of iterations performed?
          if(this->_num_iter < this->_min_iter)
            return Status::progress;

          // is converged?
          if(this->is_converged())
            return Status::success;

          // maximum number of iterations performed?
          if(this->_num_iter >= this->_max_iter)
            return Status::max_iter;

          // check for stagnation?
          if(this->_min_stag_iter > Index(0))
          {
            // did this iteration stagnate?
            if(Math::sqr(_alpha)*this->_def_init <= Math::eps<DataType>())
            {
              // increment stagnation count
              if(++this->_num_stag_iter >= this->_min_stag_iter)
              {
                //std::cout << this->_plot_name <<" stagnated with " <<
                //  stringify_fp_sci(Math::sqr(_alpha)*this->_def_init) << " <= " <<
                //  stringify_fp_sci(Math::sqr(this->_tol_rel)) << std::endl;
                return Status::stagnated;
              }
            }
            else
            {
              // this iteration did not stagnate
              this->_num_stag_iter = Index(0);
            }
          }

          // continue iterating
          return Status::progress;
        }

    };
  } // namespace Solver
} // namespace FEAST
#endif // FEAST_KERNEL_SOLVER_LINESEARCH
