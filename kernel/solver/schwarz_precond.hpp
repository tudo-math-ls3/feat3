// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/solver/base.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/filter.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Schwarz preconditioner class template declaration.
     *
     * \note This class is only specialized for Global::Vector types.
     * There exists no generic implementation.
     *
     * \tparam Vector_
     * The type of the vector. Must be an instance of Global::Vector.
     *
     * \tparam Filter_
     * The type of the filter. Must be an instance of Global::Filter.
     */
    template<typename Vector_, typename Filter_>
    class SchwarzPrecond;

    /**
     * \brief Schwarz preconditioner specialization for Global::Vector
     *
     * This class template implements an additive Schwarz preconditioner,
     * i.e. this class implements a global preconditioner by synchronously
     * adding (and averaging) a set of local solutions.
     *
     * \note
     * By default, this class synchronizes the returned status codes of all local solvers to ensure that all of them
     * have returned a successful status code. This synchronization requires a global reduction operation, which can
     * be very expensive when the local solver is a very cheap smoother/preconditioner, which is applied very frequently.
     * Therefore, it is possible (and recommended) to ignore the status codes of the local solvers in such a scenario,
     * since the failure of any local solver will lead to a failure of the outer solver anyways, which can then take
     * care of the error handling.
     *
     * \tparam LocalVector_
     * The type of the local vector nested in a global vector.
     *
     * \author Peter Zajac
     */
    template<typename LocalVector_, typename LocalFilter_, typename Mirror_>
    class SchwarzPrecond<Global::Vector<LocalVector_, Mirror_>, Global::Filter<LocalFilter_, Mirror_>> :
      public SolverBase<Global::Vector<LocalVector_, Mirror_>>
    {
    public:
      /// our local vector type
      typedef LocalVector_ LocalVectorType;
      /// our global vector type
      typedef Global::Vector<LocalVector_, Mirror_> GlobalVectorType;
      /// our global filter type
      typedef Global::Filter<LocalFilter_, Mirror_> GlobalFilterType;
      /// base-class typedef
      typedef SolverBase<Global::Vector<LocalVector_, Mirror_>> BaseClass;

      /// the local solver interface
      typedef SolverBase<LocalVector_> LocalSolverType;

    protected:
      /// our local solver object
      std::shared_ptr<LocalSolverType> _local_solver;
      /// our global filter
      const GlobalFilterType& _filter;
      /// whether to ignore status codes of local solvers
      bool _ignore_status;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] local_solver
       * The local solver that is to be used by the Schwarz preconditioner.
       *
       * \param[in] filter
       * The global filter that is to be applied.
       *
       * \param[in] ignore_status
       * Specifies whether the status codes of the local solvers are to be ignored.
       */
      explicit SchwarzPrecond(std::shared_ptr<LocalSolverType> local_solver, const GlobalFilterType& filter, bool ignore_status = false) :
        _local_solver(local_solver),
        _filter(filter),
        _ignore_status(ignore_status)
      {
      }

      /**
       * \brief Constructor using a PropertyMap
       *
       * \param[in] section_name
       * The name of the config section, which it does not know by itself
       *
       * \param[in] section
       * A pointer to the PropertyMap section configuring this solver
       *
       * \param[in] local_solver
       * The local solver that is to be used by the Schwarz preconditioner.
       *
       * \param[in] filter
       * The global filter that is to be applied.
       */
      explicit SchwarzPrecond(const String& section_name, PropertyMap* section,
        std::shared_ptr<LocalSolverType> local_solver, const GlobalFilterType& filter) :
        BaseClass(section_name, section),
        _local_solver(local_solver),
        _filter(filter),
        _ignore_status(false)
      {
        auto ignore_status_p = section->query("ignore_status");
        if(ignore_status_p.second)
        {
          if(ignore_status_p.first.is_one_of("yes true"))
            this->_ignore_status = true;
          else if(ignore_status_p.first.is_one_of("no false"))
            this->_ignore_status = false;
          else
            throw ParseError(section_name + ".ignore_status", ignore_status_p.first, "yes/no or true/false");
        }
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "Schwarz";
      }

      void set_ignore_status(bool ignore_status)
      {
        this->_ignore_status = ignore_status;
      }

      virtual void init_symbolic() override
      {
        _local_solver->init_symbolic();
      }

      virtual void init_numeric() override
      {
        _local_solver->init_numeric();
      }

      virtual void done_numeric() override
      {
        _local_solver->done_numeric();
      }

      virtual void done_symbolic() override
      {
        _local_solver->done_symbolic();
      }

      virtual Status apply(GlobalVectorType& vec_cor, const GlobalVectorType& vec_def) override
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        // apply local solver
        Statistics::add_solver_expression(std::make_shared<ExpressionCallPrecond>(this->name(), this->_local_solver->name()));
        Status status = _local_solver->apply(vec_cor.local(), vec_def.local());

        // ignore or synchronize status codes?
        if(this->_ignore_status)
          status = Status::success;
        else
        {
          // synchronize local status over communicator to obtain
          // a consistent status code on all processes
          const Dist::Comm* comm = vec_cor.get_comm();
          if(comm != nullptr)
          {
            // local status: 0: success, 1: failure
            int lstat = (status_success(status) ? 0 : 1);
            int gstat = -1;

            // sum up over all processes:
            comm->allreduce(&lstat, &gstat, std::size_t(1), Dist::op_sum);
            XASSERT(gstat >= 0);

            // if all processes succeeded, the sum will also be 0
            // if the sum is > 0, then at least one process failed,
            // so the global status will also be failure
            status = (gstat == 0 ? Status::success : Status::aborted);
          }
        }

        // did all processes succeed?
        if(status == Status::success)
        {
          // synchronize correction vector
          vec_cor.sync_1();

          // apply filter
          _filter.filter_cor(vec_cor);
        }

        // okay
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, 0));
        return status;
      }
    }; // class SchwarzPrecond<...>

    /**
     * \brief Creates a new SchwarzPrecond solver object
     *
     * \param[in] local_solver
     * The local solver object.
     *
     * \param[in] filter
     * The global system filter.
     *
     * \returns
     * A shared pointer to a new SchwarzPrecond object.
     */
    template<typename LocalFilter_, typename Mirror_>
    inline std::shared_ptr<SchwarzPrecond<Global::Vector<typename LocalFilter_::VectorType, Mirror_>, Global::Filter<LocalFilter_, Mirror_>>> new_schwarz_precond(
      std::shared_ptr<SolverBase<typename LocalFilter_::VectorType>> local_solver,
      const Global::Filter<LocalFilter_, Mirror_>& filter)
    {
      return std::make_shared<SchwarzPrecond<Global::Vector<typename LocalFilter_::VectorType, Mirror_>, Global::Filter<LocalFilter_, Mirror_>>>
        (local_solver, filter);
    }

    /**
     * \brief Creates a new SchwarzPrecond solver object using a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] local_solver
     * The local solver object.
     *
     * \param[in] filter
     * The global system filter.
     *
     * \returns
     * A shared pointer to a new SchwarzPrecond object.
     */
    template<typename LocalFilter_, typename Mirror_>
    inline std::shared_ptr
    <
      SchwarzPrecond
      <
        Global::Vector<typename LocalFilter_::VectorType, Mirror_>,
        Global::Filter<LocalFilter_, Mirror_>
      >
    >
    new_schwarz_precond(
      const String& section_name, PropertyMap* section,
      std::shared_ptr<SolverBase<typename LocalFilter_::VectorType>> local_solver,
      const Global::Filter<LocalFilter_, Mirror_>& filter)
    {
      return std::make_shared<SchwarzPrecond<Global::Vector<typename LocalFilter_::VectorType, Mirror_>, Global::Filter<LocalFilter_, Mirror_>>>
        (section_name, section, local_solver, filter);
    }
  } // namespace Solver
} // namespace FEAT
