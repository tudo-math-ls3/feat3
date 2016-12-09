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
       * \brief Reads a solver configuration from a PropertyMap
       */
      virtual void read_config(PropertyMap* section) override
      {
        BaseClass::read_config(section);

        // Check if we have set _krylov_vim
        auto omega_p = section->query("omega");
        if(omega_p.second)
        {
          set_omega(DataType(std::stod(omega_p.first)));
        }
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
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_SCALE_PRECOND_HPP
