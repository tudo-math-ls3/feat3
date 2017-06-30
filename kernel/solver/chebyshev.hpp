#pragma once
#ifndef KERNEL_SOLVER_CHEBYSHEV_PRECOND_HPP
#define KERNEL_SOLVER_CHEBYSHEV_PRECOND_HPP 1

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/solver/iterative.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Chebyshev Polynomial implementation
     *
     * This class represents the Chebyshev-Polynomial "Solver"
     *
     */
    template<typename Matrix_, typename Filter_>
    class Chebyshev :
      public IterativeSolver<typename Matrix_::VectorTypeR>
    {
    public:
      /// The matrix type
      typedef Matrix_ MatrixType;
      /// The filter type
      typedef Filter_ FilterType;
      /// The type of vector this solver can be applied to
      typedef typename MatrixType::VectorTypeL VectorType;
      /// The floating point precision
      typedef typename MatrixType::DataType DataType;
      /// Our base class
      typedef IterativeSolver<VectorType> BaseClass;

    protected:
      /// the matrix for the solver
      const MatrixType& _system_matrix;
      /// the filter for the solver
      const FilterType& _system_filter;
      /// defect vector
      VectorType _vec_def;
      /// correction vector
      VectorType _vec_cor;
      /// the minium eigenvalue of the matrix
      DataType _min_ev;
      /// the maximum eigenvalue of the matrix
      DataType _max_ev;
      /// minimum eigenvalue fraction
      DataType _fraction_min_ev;
      /// maximum eigenvalue fraction
      DataType _fraction_max_ev;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The source matrix.
       *
       * \param[in] filter
       * The system filter.
       *
       * \param[in] m
       * The order of the polynom
       *
       */
      explicit Chebyshev(const MatrixType& matrix, const FilterType& filter,
          const DataType fraction_min_ev = DataType(0.5), const DataType fraction_max_ev = DataType(0.8)) :
        BaseClass("Chebyshev"),
        _system_matrix(matrix),
        _system_filter(filter),
        _min_ev(0),
        _max_ev(0),
        _fraction_min_ev(fraction_min_ev),
        _fraction_max_ev(fraction_max_ev)
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
       * \param[in] matrix
       * The system matrix.
       *
       * \param[in] filter
       * The system filter.
       *
       * \returns
       * A shared pointer to a new Chebyshev object.
       */
      explicit Chebyshev(const String& section_name, PropertyMap* section,
      const MatrixType& matrix, const FilterType& filter) :
        BaseClass("Chebyshev", section_name, section),
        _system_matrix(matrix),
        _system_filter(filter),
        _min_ev(DataType(0)),
        _max_ev(DataType(0)),
        _fraction_min_ev(DataType(0.5)),
        _fraction_max_ev(DataType(0.8))
        {
          auto fmin_p = section->query("fraction_min_ev");
          if(fmin_p.second)
          {
            set_fraction_min_ev(DataType(std::stod(fmin_p.first)));
          }
          auto fmax_p = section->query("fraction_max_ev");
          if(fmax_p.second)
          {
            set_fraction_max_ev(DataType(std::stod(fmax_p.first)));
          }
        }

      /**
       * \brief Empty virtual destructor
       */
      virtual ~Chebyshev()
      {
      }

      /// \copydoc SolverBase::name()
      virtual String name() const override
      {
        return "Chebyshev";
      }

      /**
       * \brief Sets the minimum eigenvalue fraction
       *
       * \param[in] fraction_min_ev
       * The new minimum eigenvalue fraction.
       *
       */
      void set_fraction_min_ev(DataType fraction_min_ev)
      {
        _fraction_min_ev = fraction_min_ev;
      }

      /**
       * \brief Sets the maximum eigenvalue fraction
       *
       * \param[in] fraction_max_ev
       * The new maximum eigenvalue fraction.
       *
       */
      void set_fraction_max_ev(DataType fraction_max_ev)
      {
        _fraction_max_ev = fraction_max_ev;
      }

      /// \copydoc SolverBase::init_symbolic()
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        _vec_def = this->_system_matrix.create_vector_r();
        _vec_cor = this->_system_matrix.create_vector_r();
      }

      /// \copydoc SolverBase::done_symbolic()
      virtual void done_symbolic() override
      {
        this->_vec_cor.clear();
        this->_vec_def.clear();
        BaseClass::done_symbolic();
      }

      /// \copydoc SolverBase::init_numeric()
      virtual void init_numeric() override
      {
        BaseClass::init_numeric();

        //PowerMethod for maximum eigenvalue
        //http://blogs.sas.com/content/iml/2012/05/09/the-power-method.html
        //https://en.wikipedia.org/wiki/Power_iteration
        VectorType& v(this->_vec_cor);
        VectorType& z(this->_vec_def);
        const MatrixType& matrix(this->_system_matrix);
        v.format(DataType(3));
        DataType tolerance(DataType(1e-4));
        Index max_iters(7);
        v.scale(v, DataType(1) / v.norm2());
        DataType lambda_old(0), lambda(0);
        for (Index i(0) ; i < max_iters ; ++i)
        {
          matrix.apply(z, v);
          v.scale(z, DataType(1) / z.norm2());
          lambda = v.dot(z);
          if (Math::abs((lambda - lambda_old) / lambda) < tolerance)
            break;
          lambda_old = lambda;
        }

        _min_ev = lambda * _fraction_min_ev;
        _max_ev = lambda * _fraction_max_ev;
      }

      /// \copydoc SolverBase::write_config()
      virtual PropertyMap* write_config(PropertyMap* parent, const String& new_section_name) const override
      {

        PropertyMap* my_section = BaseClass::write_config(parent, new_section_name);

        my_section->add_entry("fraction_min_ev", stringify_fp_sci(_fraction_min_ev));
        my_section->add_entry("fraction_max_ev", stringify_fp_sci(_fraction_max_ev));

        return my_section;
      }

      virtual Status apply(VectorType& vec_cor, const VectorType& vec_def) override
      {
        // save defect
        this->_vec_def.copy(vec_def);
        //this->_system_filter.filter_def(this->_vec_def);

        // clear solution vector
        vec_cor.format();

        // apply
        Status st(_apply_intern(vec_cor, vec_def));
        this->plot_summary(st);
        return st;
      }

      virtual Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
      {
        // compute defect
        this->_system_matrix.apply(this->_vec_def, vec_sol, vec_rhs, -DataType(1));
        this->_system_filter.filter_def(this->_vec_def);

        // apply
        Status st(_apply_intern(vec_sol, vec_rhs));
        this->plot_summary(st);
        return st;
      }

    protected:
      virtual Status _apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
      {
        Statistics::add_solver_expression(std::make_shared<ExpressionStartSolve>(this->name()));

        VectorType& vec_def(this->_vec_def);
        VectorType& vec_cor(this->_vec_cor);
        const MatrixType& matrix(this->_system_matrix);
        const FilterType& filter(this->_system_filter);

        // compute initial defect
        Status status = this->_set_initial_defect(vec_def, vec_sol);

        DataType d = (_max_ev + _min_ev) / DataType(2);
        DataType c = (_max_ev - _min_ev) / DataType(2);
        DataType alpha(0);
        DataType beta(0);
        vec_cor.scale(vec_def, (DataType(1) / d));

        // start iterating
        while(status == Status::progress)
        {
          switch (this->get_num_iter())
          {
          case 0:
            alpha = DataType(1) / d;
            break;
          case 1:
            alpha = DataType(2) * d * (DataType(1) / ((DataType(2) * d * d) - (c * c)));
            break;
          default:
            alpha = DataType(1) / (d - ((alpha * c * c) / DataType(4)));
          }

          beta = (alpha * d) - DataType(1);
          vec_cor.scale(vec_cor, beta);
          vec_cor.axpy(vec_def, vec_cor, alpha);

          // update solution vector
          vec_sol.axpy(vec_cor, vec_sol);

          // compute new defect vector
          matrix.apply(vec_def, vec_sol, vec_rhs, -DataType(1));
          filter.filter_def(vec_def);

          // compute new defect norm
          status = this->_set_new_defect(vec_def, vec_sol);
        }

        // return our status
        Statistics::add_solver_expression(std::make_shared<ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
        return status;
      }
    }; // class Chebyshev<...>

    /**
     * \brief Creates a new Chebyshev solver object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter.
     *
     * \param[in] m
     * The order of the polynom
     *
     * \returns
     * A shared pointer to a new Chebyshev object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<Chebyshev<Matrix_, Filter_>> new_chebyshev(
      const Matrix_& matrix, const Filter_& filter,
      const typename Matrix_::DataType fraction_min_ev = typename Matrix_::DataType(0.5),
      const typename Matrix_::DataType fraction_max_ev = typename Matrix_::DataType(0.8))
    {
      return std::make_shared<Chebyshev<Matrix_, Filter_>>(matrix, filter, fraction_min_ev, fraction_max_ev);
    }

    /**
     * \brief Creates a new Chebyshev solver object using a PropertyMap
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
     * \returns
     * A shared pointer to a new Chebyshev object.
     */
    template<typename Matrix_, typename Filter_>
    inline std::shared_ptr<Chebyshev<Matrix_, Filter_>> new_chebyshev(
      const String& section_name, PropertyMap* section,
      const Matrix_& matrix, const Filter_& filter)
    {
      return std::make_shared<Chebyshev<Matrix_, Filter_>>(section_name, section, matrix, filter);
    }
  } // namespace Solver
} // namespace FEAT

#endif // KERNEL_SOLVER_CHEBYSHEV_PRECOND_HPP
