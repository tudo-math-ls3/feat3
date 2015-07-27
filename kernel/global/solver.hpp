#pragma once
#ifndef KERNEL_GLOBAL_SOLVER_HPP
#define KERNEL_GLOBAL_SOLVER_HPP 1

#include <kernel/lafem/proto_solver.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/vector.hpp>

namespace FEAST
{
  namespace Global
  {
    template<typename Vector_>
    class SchwarzPrecond;

    template<typename LocalVector_>
    class SchwarzPrecond<Global::Vector<LocalVector_>> :
      public LAFEM::SolverInterface<Global::Vector<LocalVector_>>
    {
    public:
      typedef LocalVector_ LocalVectorType;
      typedef Global::Vector<LocalVector_> GlobalVectorType;

      typedef LAFEM::SolverInterface<LocalVector_> LocalSolverType;

    protected:
      std::shared_ptr<LocalSolverType> _local_solver;

    public:
      explicit SchwarzPrecond(std::shared_ptr<LocalSolverType> local_solver) :
        _local_solver(local_solver)
      {
      }

      virtual bool init_symbolic() override
      {
        return _local_solver->init_symbolic();
      }

      virtual bool init_numeric() override
      {
        return _local_solver->init_numeric();
      }

      virtual void done_numeric() override
      {
        _local_solver->done_numeric();
      }

      virtual void done_symbolic() override
      {
        _local_solver->done_symbolic();
      }

      virtual LAFEM::SolverStatus solve(GlobalVectorType& vec_sol, const GlobalVectorType& vec_rhs) override
      {
        // apply local solver
        LAFEM::SolverStatus status = _local_solver->solve(*vec_sol, *vec_rhs);
        if(!LAFEM::status_success(status))
          return status;

        // sync
        vec_sol.sync_1();

        // okay
        return status;
      }
    }; // class SchwarzPrecond<...>

    //template<typename Matrix_, typename Filter_>
    //class SchurPrecond;

    template<typename MatrixA_, typename MatrixB_, typename MatrixD_, typename FilterV_, typename FilterP_>
    class SchurPrecond : //<
      //Global::Matrix<LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_>,
      //Global::Filter<LAFEM::TupleFilter<FilterV_, Filter_P_>> :
      public LAFEM::SolverInterface<Global::Vector<LAFEM::TupleVector<typename MatrixB_::VectorTypeL, typename MatrixD_::VectorTypeL>>>
    {
    public:
      typedef LAFEM::SaddlePointMatrix<MatrixA_, MatrixB_, MatrixD_> LocalMatrixType;
      typedef LAFEM::TupleFilter<FilterV_, FilterP_> LocalFilterType;
      typedef LAFEM::TupleVector<typename MatrixB_::VectorTypeL, typename MatrixD_::VectorTypeL> LocalVectorType;

      typedef Global::Matrix<LocalMatrixType> GlobalMatrixType;
      typedef Global::Filter<LocalFilterType> GlobalFilterType;
      typedef Global::Vector<LocalVectorType> GlobalVectorType;

      typedef typename GlobalMatrixType::DataType DataType;

      typedef LAFEM::SolverInterface<GlobalVectorType> BaseClass;

      typedef Global::Matrix<MatrixA_> GlobalMatrixTypeA;
      typedef Global::Matrix<MatrixB_> GlobalMatrixTypeB;
      typedef Global::Matrix<MatrixD_> GlobalMatrixTypeD;

      typedef typename MatrixB_::VectorTypeL LocalVectorTypeV;
      typedef typename MatrixD_::VectorTypeL LocalVectorTypeP;

      typedef Global::Vector<LocalVectorTypeV> GlobalVectorTypeV;
      typedef Global::Vector<LocalVectorTypeP> GlobalVectorTypeP;

      typedef Global::Filter<FilterV_> GlobalFilterTypeV;
      typedef Global::Filter<FilterP_> GlobalFilterTypeP;

      typedef Global::Gate<LocalVectorTypeV> GateTypeV;
      typedef Global::Gate<LocalVectorTypeP> GateTypeP;

      typedef LAFEM::SolverInterface<GlobalVectorTypeV> SolverA;
      typedef LAFEM::SolverInterface<GlobalVectorTypeP> SolverS;

    protected:
      //const GlobalMatrixTypeA& _matrix_a;
      const GlobalMatrixTypeB& _matrix_b;
      const GlobalMatrixTypeD& _matrix_d;
      const GlobalFilterTypeV& _filter_v;
      const GlobalFilterTypeP& _filter_p;
      GlobalVectorTypeV _vec_rhs_v, _vec_sol_v, _vec_def_v;
      GlobalVectorTypeP _vec_rhs_p, _vec_sol_p, _vec_def_p;
      std::shared_ptr<SolverA> _solver_a;
      std::shared_ptr<SolverS> _solver_s;
      LAFEM::SchurType _schur_type;

    public:
      explicit SchurPrecond(
        //const GlobalMatrixTypeA& matrix_a,
        const GlobalMatrixTypeB& matrix_b,
        const GlobalMatrixTypeD& matrix_d,
        const GlobalFilterTypeV& filter_v,
        const GlobalFilterTypeP& filter_p,
        std::shared_ptr<SolverA> solver_a,
        std::shared_ptr<SolverS> solver_s,
        LAFEM::SchurType type = LAFEM::SchurType::diagonal
        ) :
        //_matrix_a(matrix_a),
        _matrix_b(matrix_b),
        _matrix_d(matrix_d),
        _filter_v(filter_v),
        _filter_p(filter_p),
        _solver_a(solver_a),
        _solver_s(solver_s),
        _schur_type(type)
      {
      }

      virtual String name() const override
      {
        return "Global::Schur";
      }

      virtual bool init_symbolic() override
      {
        if(!BaseClass::init_symbolic())
          return false;
        if(!_solver_a->init_symbolic())
          return false;
        if(!_solver_s->init_symbolic())
          return false;

        _vec_rhs_v = _matrix_b.create_vector_l();
        _vec_sol_v = _matrix_b.create_vector_l();
        _vec_rhs_p = _matrix_d.create_vector_l();
        _vec_sol_p = _matrix_d.create_vector_l();

        if(_schur_type != LAFEM::SchurType::diagonal)
        {
          _vec_def_v = _matrix_b.create_vector_l();
          _vec_def_p = _matrix_d.create_vector_l();
        }

        return true;
      }

      virtual bool init_numeric() override
      {
        if(!BaseClass::init_numeric())
          return false;
        if(!_solver_a->init_numeric())
          return false;
        if(!_solver_s->init_numeric())
          return false;
        return true;
      }

      virtual void done_numeric() override
      {
        _solver_s->done_numeric();
        _solver_a->done_numeric();
        BaseClass::done_numeric();
      }

      virtual void done_symbolic() override
      {
        if(_schur_type != LAFEM::SchurType::diagonal)
        {
          _vec_def_p.clear();
          _vec_def_v.clear();
        }
        _vec_sol_p.clear();
        _vec_rhs_p.clear();
        _vec_sol_v.clear();
        _vec_rhs_v.clear();
        _solver_s->done_symbolic();
        _solver_a->done_symbolic();
        BaseClass::done_symbolic();
      }

      virtual LAFEM::SolverStatus solve(GlobalVectorType& vec_sol, const GlobalVectorType& vec_rhs) override
      {
        // first of all, copy RHS
        (*_vec_rhs_v).copy((*vec_rhs).template at<0>());
        (*_vec_rhs_p).copy((*vec_rhs).template at<1>());

        // now let's check the preconditioner type
        switch(_schur_type)
        {
        case LAFEM::SchurType::diagonal:
          {
            // solve A*u_v = f_v
            if(!LAFEM::status_success(_solver_a->solve(_vec_sol_v, _vec_rhs_v)))
              return LAFEM::SolverStatus::aborted;

            // solve S*u_p = f_p
            if(!LAFEM::status_success(_solver_s->solve(_vec_sol_p, _vec_rhs_p)))
              return LAFEM::SolverStatus::aborted;

            // okay
            break;
          }

        case LAFEM::SchurType::lower:
          {
            // solve A*u_v = f_v
            if(!LAFEM::status_success(_solver_a->solve(_vec_sol_v, _vec_rhs_v)))
              return LAFEM::SolverStatus::aborted;

            // compute g_p := f_p - D*u_v
            _matrix_d.apply(_vec_def_p, _vec_sol_v, _vec_rhs_p, -DataType(1));

            // apply pressure defect filter
            _filter_p.filter_def(_vec_def_p);

            // solve S*u_p = g_p
            if(!LAFEM::status_success(_solver_s->solve(_vec_sol_p, _vec_def_p)))
              return LAFEM::SolverStatus::aborted;

            // okay
            break;
          }

        case LAFEM::SchurType::upper:
          {
            // solve S*u_p = f_p
            if(!LAFEM::status_success(_solver_s->solve(_vec_sol_p, _vec_rhs_p)))
              return LAFEM::SolverStatus::aborted;

            // compute g_v := f_v - B*u_p
            _matrix_b.apply(_vec_def_v, _vec_sol_p, _vec_rhs_v, -DataType(1));

            // apply velocity defect filter
            _filter_v.filter_def(_vec_def_v);

            // solve A*u_v = g_v
            if(!LAFEM::status_success(_solver_a->solve(_vec_sol_v, _vec_def_v)))
              return LAFEM::SolverStatus::aborted;

            // okay
            break;
          }

        case LAFEM::SchurType::full:
          {
            // Note: We will use the first component of the solution vector here.
            //       It will be overwritten by the third solution step below.
            // solve A*u_v = f_v
            if(!LAFEM::status_success(_solver_a->solve(_vec_sol_v, _vec_rhs_v)))
              return LAFEM::SolverStatus::aborted;

            // compute g_p := f_p - D*u_v
            _matrix_d.apply(_vec_def_p, _vec_sol_v, _vec_rhs_p, -DataType(1));

            // apply pressure defect filter
            _filter_p.filter_def(_vec_def_p);

            // solve S*u_p = g_p
            if(!LAFEM::status_success(_solver_s->solve(_vec_sol_p, _vec_def_p)))
              return LAFEM::SolverStatus::aborted;

            // compute g_v := f_v - B*u_p
            _matrix_b.apply(_vec_def_v, _vec_sol_p, _vec_rhs_v, -DataType(1));

            // apply velocity defect filter
            _filter_v.filter_def(_vec_def_v);

            // solve A*u_v = g_v
            if(!LAFEM::status_success(_solver_a->solve(_vec_sol_v, _vec_def_v)))
              return LAFEM::SolverStatus::aborted;

            // okay
            break;
          }
        }

        // finally, copy sol
        (*vec_sol).template at<0>().copy(*_vec_sol_v);
        (*vec_sol).template at<1>().copy(*_vec_sol_p);

        // okay
        return LAFEM::SolverStatus::success;
      }
    };
  } // namespace Global
} // namespace FEAST

#endif // KERNEL_GLOBAL_SOLVER_HPP
