#pragma once
#ifndef KERNEL_SOLVER_BASIC_VCYCLE_HPP
#define KERNEL_SOLVER_BASIC_VCYCLE_HPP 1

// includes, FEAST
#include <kernel/solver/base.hpp>

namespace FEAST
{
  namespace Solver
  {
    /**
     * \brief Basic Multigrid V-Cycle preconditioner class template
     *
     * This class implements a V-Cycle multigrid preconditioner.
     *
     * \attention This class does not implement a solver, but merely a preconditioner representing
     * a single cycle. If you want a "whole" multigrid solver, put an instance of this class as a
     * preconditioner into a FixPoint solver.
     *
     * \tparam Matrix_
     * The type of the system matrix.
     *
     * \tparam Filter_
     * The type of the system filter.
     *
     * \tparam TransferMatrix_
     * The type of the prolongation/restriction matrices.
     *
     * \author Peter Zajac
     */
    template<
      typename Matrix_,
      typename Filter_,
      typename TransferMatrix_>
    class BasicVCycle :
      public SolverBase<typename Matrix_::VectorTypeR>
    {
    public:
      typedef SolverBase<typename Matrix_::VectorTypeR> BaseClass;

      typedef Matrix_ MatrixType;
      typedef Filter_ FilterType;
      typedef TransferMatrix_ TransferMatrixType;
      typedef typename MatrixType::VectorTypeR VectorType;
      typedef typename MatrixType::DataType DataType;

      typedef SolverBase<VectorType> SmootherType;

    protected:
      /**
       * \brief System level structure
       */
      struct SystemLevel
      {
        /// the system matrix on this level
        const MatrixType& matrix;
        /// the filter on this level
        const FilterType& filter;
        /// temporary vectors
        VectorType vec_rhs, vec_sol, vec_def, vec_cor;
        /// the pre-smoother on this level
        std::shared_ptr<SmootherType> pre_smoother;
        /// the post-smoother on this level
        std::shared_ptr<SmootherType> post_smoother;
        /// the coarse-grid solver if this is the coarse level
        std::shared_ptr<SmootherType> coarse_solver;

        SystemLevel(const MatrixType& my_matrix, const FilterType& my_filter,
          std::shared_ptr<SmootherType> my_coarse_solver) :
          matrix(my_matrix),
          filter(my_filter),
          pre_smoother(nullptr),
          post_smoother(nullptr),
          coarse_solver(my_coarse_solver)
        {
        }

        SystemLevel(const MatrixType& my_matrix, const FilterType& my_filter,
          std::shared_ptr<SmootherType> my_pre_smoother, std::shared_ptr<SmootherType> my_post_smoother) :
          matrix(my_matrix),
          filter(my_filter),
          pre_smoother(my_pre_smoother),
          post_smoother(my_post_smoother),
          coarse_solver(nullptr)
        {
        }
      };

      /**
       * \brief Transfer level structure.
       */
      struct TransferLevel
      {
        /// the prolongation matrix on this level
        const TransferMatrixType& prol;
        /// the restriction matrix on this level
        const TransferMatrixType& rest;

        TransferLevel(const TransferMatrixType& prol_mat, const TransferMatrixType& rest_mat) :
          prol(prol_mat),
          rest(rest_mat)
        {
        }
      };

      /// our system levels
      std::deque<SystemLevel*> system_levels;
      /// our transfer levels
      std::deque<TransferLevel*> transfer_levels;

    public:
      /// constructor
      BasicVCycle()
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "VCycle";
      }

      /// destructor
      virtual ~BasicVCycle()
      {
        while(!system_levels.empty())
        {
          delete system_levels.back();
          system_levels.pop_back();
        }
        while(!transfer_levels.empty())
        {
          delete transfer_levels.back();
          transfer_levels.pop_back();
        }
      }

      /**
       * \brief Sets the coarse level.
       *
       * \param[in] matrix
       * The coarse level matrix.
       *
       * \param[in] filter
       * The coarse level filter.
       *
       * \param[in] coarse_solver
       * The coarse level solver.
       */
      void set_coarse_level(
        const MatrixType& matrix,
        const FilterType& filter,
        std::shared_ptr<SmootherType> coarse_solver)
      {
        if(!system_levels.empty())
          throw InternalError("Coarse level already exists");

        // push the coarse level
        system_levels.push_back(new SystemLevel(matrix, filter, coarse_solver));
      }

      /**
       * \brief Pushes a new level into the hierarchy.
       *
       * \param[in] matrix
       * The system matrix for the new level.
       *
       * \param[in] filter
       * The system filter for the new level.
       *
       * \param[in] prol_mat
       * The prolongation matrix from the previous level to this new level.
       *
       * \param[in] rest_mat
       * The restriction matrix from this new level to the previous level.
       *
       * \param[in] pre_smoother
       * The pre-smoother. May be \p nullptr.
       *
       * \param[in] post_smoother
       * The post-smoother. May be \p nullptr.
       *
       * \note
       * Both \p pre_smoother and \p post_smoother may point to the same object. This class takes
       * care of this scenario and will ensure that the \p init and \p done methods of those objects
       * are not called twice in this case.
       */
      void push_level(
        const MatrixType& matrix,
        const FilterType& filter,
        const TransferMatrixType& prol_mat,
        const TransferMatrixType& rest_mat,
        std::shared_ptr<SmootherType> pre_smoother,
        std::shared_ptr<SmootherType> post_smoother)
      {
        if(system_levels.empty())
          throw InternalError("Coarse level is missing");

        // push the new level
        system_levels.push_back(new SystemLevel(matrix, filter, pre_smoother, post_smoother));
        transfer_levels.push_back(new TransferLevel(prol_mat, rest_mat));
      }

      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();

        // loop over all levels
        for(auto it = system_levels.begin(); it != system_levels.end(); ++it)
        {
          (*it)->vec_rhs = (*it)->matrix.create_vector_r();
          (*it)->vec_sol = (*it)->matrix.create_vector_r();
          if(it == system_levels.begin())
          {
            if((*it)->coarse_solver)
            {
              (*it)->coarse_solver->init_symbolic();
            }
          }
          else
          {
            (*it)->vec_cor = (*it)->matrix.create_vector_r();
            (*it)->vec_def = (*it)->matrix.create_vector_r();
            if((*it)->pre_smoother)
            {
              (*it)->pre_smoother->init_symbolic();
            }
            if((*it)->post_smoother)
            {
              if((*it)->post_smoother != (*it)->pre_smoother)
              {
                (*it)->post_smoother->init_symbolic();
              }
            }
          }
        }
      }

      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        // loop over all level
        for(auto it = system_levels.begin(); it != system_levels.end(); ++it)
        {
          if(it == system_levels.begin())
          {
            if((*it)->coarse_solver)
            {
              (*it)->coarse_solver->init_numeric();
            }
          }
          else
          {
            if((*it)->pre_smoother)
            {
              (*it)->pre_smoother->init_numeric();
            }
            if((*it)->post_smoother)
            {
              if((*it)->post_smoother != (*it)->pre_smoother)
              {
                (*it)->post_smoother->init_numeric();
              }
            }
          }
        }
      }

      virtual void init_branch(String root = "") override
      {
        BaseClass::init_branch(root);

        // loop over all level
        for(auto it = system_levels.begin(); it != system_levels.end(); ++it)
        {
          if(it == system_levels.begin())
          {
            if((*it)->coarse_solver)
            {
              (*it)->coarse_solver->init_branch(root + "::" + this->name());
            }
          }
          else
          {
            if((*it)->pre_smoother)
            {
              (*it)->pre_smoother->init_branch(root + "::" + this->name());
            }
            if((*it)->post_smoother)
            {
              if((*it)->post_smoother != (*it)->pre_smoother)
              {
                (*it)->post_smoother->init_branch(root + "::" + this->name());
              }
            }
          }
        }
      }

      virtual void done_numeric() override
      {
        // loop over all level
        for(auto it = system_levels.begin(); it != system_levels.end(); ++it)
        {
          if(it == system_levels.begin())
          {
            if((*it)->coarse_solver)
            {
              (*it)->coarse_solver->done_numeric();
            }
          }
          else
          {
            if((*it)->pre_smoother)
            {
              (*it)->pre_smoother->done_numeric();
            }
            if((*it)->post_smoother)
            {
              if((*it)->post_smoother != (*it)->pre_smoother)
                (*it)->post_smoother->done_numeric();
            }
          }
        }
        BaseClass::done_numeric();
      }

      virtual void done_symbolic() override
      {
        // loop over all level
        for(auto it = system_levels.begin(); it != system_levels.end(); ++it)
        {
          (*it)->vec_rhs.clear();
          (*it)->vec_sol.clear();
          if(it == system_levels.begin())
          {
            if((*it)->coarse_solver)
            {
              (*it)->coarse_solver->done_symbolic();
            }
          }
          else
          {
            (*it)->vec_cor.clear();
            (*it)->vec_def.clear();
            if((*it)->pre_smoother)
            {
              (*it)->pre_smoother->done_symbolic();
            }
            if((*it)->post_smoother)
            {
              if((*it)->post_smoother != (*it)->pre_smoother)
                (*it)->post_smoother->done_symbolic();
            }
          }
        }
        BaseClass::done_symbolic();
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // get the number of levels
        const int nl = int(system_levels.size());

        // insert -1 to signal new starting v cycle
        Statistics::add_solver_toe(this->_branch, double(-1));

        // array containing toe for each processed level
        std::vector<double> toes((size_t)nl, double(0));

        // copy RHS vector
        system_levels.back()->vec_rhs.copy(vec_def);

        // restriction loop
        for(int i(nl-1); i > 0; --i)
        {
          TimeStamp at;

          // get our level and the coarse level
          SystemLevel& lvl = *system_levels.at(std::size_t(i));
          SystemLevel& lvl_c = *system_levels.at(std::size_t(i-1));

          // get our transfer level
          TransferLevel& trs = *transfer_levels.at(std::size_t(i-1));

          // apply pre-smoother if we have one
          if(lvl.pre_smoother)
          {
            // apply pre-smoother
            if(!status_success(lvl.pre_smoother->apply(lvl.vec_sol, lvl.vec_rhs)))
              return Status::aborted;

            // compute defect
            lvl.matrix.apply(lvl.vec_def, lvl.vec_sol, lvl.vec_rhs, -DataType(1));
          }
          else
          {
            // format solution
            lvl.vec_sol.format();

            // set our rhs as defect
            lvl.vec_def.copy(lvl.vec_rhs);
          }

          // filter our defect
          lvl.filter.filter_def(lvl.vec_def);

          // restrict onto coarse level
          trs.rest.apply(lvl_c.vec_rhs, lvl.vec_def);

          // filter coarse fefect
          lvl_c.filter.filter_def(lvl_c.vec_rhs);

          TimeStamp bt;
          toes.at((size_t)i)= bt.elapsed(at);

          // descent to prior level
        }

        // process the coarse grid level
        {
          TimeStamp at;

          SystemLevel& lvl = *system_levels.front();

          // if the have a coarse grid solver, apply it
          if(lvl.coarse_solver)
          {
            if(!status_success(lvl.coarse_solver->apply(lvl.vec_sol, lvl.vec_rhs)))
              return Status::aborted;
          }
          else
          {
            // simply copy the RHS thus emulating an identity solver
            lvl.vec_sol.copy(lvl.vec_rhs);

            // apply the correction filter
            lvl.filter.filter_cor(lvl.vec_sol);
          }

          TimeStamp bt;
          toes.at(0) = bt.elapsed(at);
        }

        // prolongation loop
        for(int i(1); i < nl; ++i)
        {
          TimeStamp at;

          // get our level and the coarse level
          SystemLevel& lvl = *system_levels.at(std::size_t(i));
          SystemLevel& lvl_c = *system_levels.at(std::size_t(i-1));

          // get our transfer level
          TransferLevel& trs = *transfer_levels.at(std::size_t(i-1));

          // prolongate coarse grid solution
          trs.prol.apply(lvl.vec_cor, lvl_c.vec_sol);

          // apply correction filter
          lvl.filter.filter_cor(lvl.vec_cor);

          // update our solution vector
          lvl.vec_sol.axpy(lvl.vec_cor, lvl.vec_sol);

          // apply post-smoother if we have one
          if(lvl.post_smoother)
          {
            // compute new defect
            lvl.matrix.apply(lvl.vec_def, lvl.vec_sol, lvl.vec_rhs, -DataType(1));

            // apply defect filter
            lvl.filter.filter_def(lvl.vec_def);

            // apply post-smoother
            if(!status_success(lvl.post_smoother->apply(lvl.vec_cor, lvl.vec_def)))
              return Status::aborted;

            // apply correction filter
            //lvl.filter.filter_cor(lvl.vec_cor);

            // update solution vector
            lvl.vec_sol.axpy(lvl.vec_cor, lvl.vec_sol);
          }

          TimeStamp bt;
          toes.at((size_t)i) += bt.elapsed(at);

          // ascend to next level
        }

        // copy sol vector
        vec_cor.copy(system_levels.back()->vec_sol);

        for (int i(0) ; i < nl ; ++i)
        {
          Statistics::add_solver_toe(this->_branch, toes.at((size_t)i));
        }

        // okay
        return Status::success;
      }
    }; // class BasicVCycle<...>
  } // namespace Solver
} // namespace FEAST

#endif // KERNEL_SOLVER_BASIC_VCYCLE_HPP
