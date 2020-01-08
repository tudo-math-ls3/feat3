// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_SCALE_PRECOND_HPP
#define KERNEL_SOLVER_SCALE_PRECOND_HPP 1

#include <kernel/solver/base.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Scaling preconditioner class template
     *
     * This simple preconditioner represents a scaled identity matrix.
     *
     * \author Peter Zajac
     */
    template<typename Vector_, typename Filter_>
    class ScalePrecond :
      public SolverBase<Vector_>
    {
    protected:
      /// our data type
      typedef typename Vector_::DataType DataType;
      /// the filter
      const Filter_& _filter;
      /// the scaling factor
      DataType _omega;
      /// Our base class
      typedef SolverBase<Vector_> BaseClass;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] omega
       * The scaling factor.
       */
      explicit ScalePrecond(const Filter_& filter, DataType omega) :
        _filter(filter),
        _omega(omega)
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
       * \param[in] filter
       * The system filter.
       *
       */
      explicit ScalePrecond(const String& section_name, PropertyMap* section,
        const Filter_& filter) :
        BaseClass(section_name, section),
        _filter(filter),
        _omega(0)
      {
        // Check if we have set _krylov_vim
        auto omega_p = section->query("omega");
        if(omega_p.second && !omega_p.first.parse(this->_omega))
          throw ParseError(section_name + ".omega", omega_p.first, "a positive float");
      }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~ScalePrecond()
      {
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "Scale";
      }

      /**
       * \brief Sets the damping parameter
       *
       * \param[in] omega
       * The new damping parameter.
       *
       */
      void set_omega(DataType omega)
      {
        XASSERT(omega > DataType(0));
        _omega = omega;
      }

      /** \copydoc SolverBase::apply() */
      virtual Status apply(Vector_& vec_cor, const Vector_& vec_def) override
      {
        vec_cor.scale(vec_def, this->_omega);
        _filter.filter_cor(vec_cor);
        return Status::success;
      }
    }; // class ScalePrecond<...>


    /**
     * \brief Creates a new ScalePrecond solver object
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] omega
     * The scaling parameter for the preconditioner.
     *
     * \returns
     * A shared pointer to a new ScalePrecond object.
     */
    template<typename Filter_, typename DataType_>
    inline std::shared_ptr<ScalePrecond<typename Filter_::VectorType, Filter_>> new_scale_precond(
      const Filter_& filter, DataType_ omega)
    {
      return std::make_shared<ScalePrecond<typename Filter_::VectorType, Filter_>>(filter, omega);
    }

    /**
     * \brief Creates a new ScalePrecond solver object using a PropertyMap
     *
     * \param[in] section_name
     * The name of the config section, which it does not know by itself
     *
     * \param[in] section
     * A pointer to the PropertyMap section configuring this solver
     *
     * \param[in] filter
     * The system filter.
     *
     * \returns
     * A shared pointer to a new ScalePrecond object.
     */
    template<typename Filter_>
    inline std::shared_ptr<ScalePrecond<typename Filter_::VectorType, Filter_>> new_scale_precond(
      const String& section_name, PropertyMap* section,
      const Filter_& filter)
    {
      return std::make_shared<ScalePrecond<typename Filter_::VectorType, Filter_>>(
        section_name, section, filter);
    }
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_SCALE_PRECOND_HPP
