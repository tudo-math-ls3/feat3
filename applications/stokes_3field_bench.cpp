// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

//
// Geometric Multigrid on Recursive Partitions for 3-Field Stokes Benchmarking Application
// ---------------------------------------------------------------------------------------
// This is an application is used as a testing ground and for benchmarking geometric
// multigrid solvers for 3-field Stokes formulations.
//
// WARNING:
// Do NOT use this application as a base for your own applications unless you have
// some seriously masochistic tendencies and/or you really know what you're doing.
//
// \author Peter Zajac
//
#include <kernel/util/runtime.hpp>
#include <kernel/util/memory_usage.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/util/statistics.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/lagrange3/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/assembly/discrete_projector.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/trace_assembler.hpp>
#include <kernel/assembly/velocity_analyser.hpp>
#include <kernel/assembly/discrete_evaluator.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/schwarz_precond.hpp>
#include <kernel/solver/fgmres.hpp>
#include <kernel/solver/umfpack.hpp>
#include <kernel/solver/amavanka.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/stop_watch.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_3field.hpp>
#include <control/statistics.hpp>

#include <vector>

namespace Stokes3Field
{
  using namespace FEAT;

  template<int dim_>
  class InflowFunction;

  template<>
  class InflowFunction<2> :
    public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 2;
    typedef Analytic::Image::Vector<2> ImageType;
    static constexpr bool can_value = true;

  protected:
    Real _vmax;

  public:
    explicit InflowFunction(Real vmax) :
      _vmax(vmax)
    {
    }

    template<typename Traits_>
    class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
    {
    protected:
      typedef typename Traits_::DataType DataType;
      typedef typename Traits_::PointType PointType;
      typedef typename Traits_::ValueType ValueType;

      const DataType _vmax, _d, _den;

    public:
      explicit Evaluator(const InflowFunction& function) :
        _vmax(function._vmax),
        //_d(DataType(0.41)),
        _d(DataType(1)),
        _den(_d*_d)
      {
      }

      ValueType value(const PointType& point)
      {
        const DataType y = point[1];

        ValueType val;
        val[0] = (_vmax * DataType(4) * y * (_d - y)) / _den;
        val[1] = DataType(0);
        return val;
      }
    };
  }; // class InflowFunction<2>

  template<>
  class InflowFunction<3> :
    public Analytic::Function
  {
  public:
    static constexpr int domain_dim = 3;
    typedef Analytic::Image::Vector<3> ImageType;
    static constexpr bool can_value = true;

  protected:
    Real _vmax;

  public:
    explicit InflowFunction(Real vmax) :
      _vmax(vmax)
    {
    }

    template<typename Traits_>
    class Evaluator :
      public Analytic::Function::Evaluator<Traits_>
    {
    protected:
      typedef typename Traits_::DataType DataType;
      typedef typename Traits_::PointType PointType;
      typedef typename Traits_::ValueType ValueType;

      const DataType _vmax, _d, _den;

    public:
      explicit Evaluator(const InflowFunction& function) :
        _vmax(function._vmax),
        //_d(DataType(0.41)),
        _d(DataType(1)),
        _den(_d*_d*_d*_d)
      {
      }

      ValueType value(const PointType& point)
      {
        const DataType y = point[1];
        const DataType z = point[2];

        ValueType val;
        val[0] = (_vmax * DataType(16) * y * (_d - y) * z * (_d - z)) / _den;
        val[1] = val[2] = DataType(0);
        return val;
      }
    };
  }; // class InflowFunction


  /**
   * \brief Navier-Stokes System Level class
   */
  template<
    int dim_,
    int nsc_ = (dim_*(dim_+1))/2,
    typename MemType_ = Mem::Main,
    typename DataType_ = Real,
    typename IndexType_ = Index>
  class MyStokes3FieldSystemLevel :
    public Control::Stokes3FieldSystemLevel<dim_, nsc_, MemType_, DataType_, IndexType_>
  {
  public:
    typedef Control::Stokes3FieldSystemLevel<dim_, nsc_, MemType_, DataType_, IndexType_> BaseClass;

    // the filtered local system matrix for Vanka
    typename BaseClass::LocalSystemMatrix local_matrix_sys;

    // define local filter types
    typedef LAFEM::UnitFilterBlocked<MemType_, DataType_, IndexType_, dim_> LocalVeloFilter;
    typedef LAFEM::NoneFilter<MemType_, DataType_, IndexType_> LocalPresFilter;
    typedef LAFEM::NoneFilterBlocked<MemType_, DataType_, IndexType_, nsc_> LocalStressFilter;
    typedef LAFEM::TupleFilter<LocalVeloFilter, LocalPresFilter, LocalStressFilter> LocalSystemFilter;

    // define global filter types
    typedef Global::Filter<LocalVeloFilter, typename BaseClass::VeloMirror> GlobalVeloFilter;
    typedef Global::Filter<LocalPresFilter, typename BaseClass::PresMirror> GlobalPresFilter;
    typedef Global::Filter<LocalStressFilter, typename BaseClass::StressMirror> GlobalStressFilter;
    typedef Global::Filter<LocalSystemFilter, typename BaseClass::SystemMirror> GlobalSystemFilter;

    // (global) filters
    GlobalSystemFilter filter_sys;
    GlobalVeloFilter filter_velo;
    GlobalPresFilter filter_pres;
    GlobalStressFilter filter_stress;

    /// \brief Returns the total amount of bytes allocated.
    std::size_t bytes() const
    {
      return this->filter_sys.bytes() + BaseClass::bytes();
    }

    void compile_system_filter()
    {
      filter_sys.local().template at<0>() = filter_velo.local().clone(LAFEM::CloneMode::Shallow);
      filter_sys.local().template at<1>() = filter_pres.local().clone(LAFEM::CloneMode::Shallow);
      filter_sys.local().template at<2>() = filter_stress.local().clone(LAFEM::CloneMode::Shallow);
    }

    void compile_local_matrix()
    {
      this->local_matrix_sys.template at<0,0>() = this->matrix_a.convert_to_1();
      this->local_matrix_sys.template at<0,1>() = this->matrix_b.local().clone(LAFEM::CloneMode::Weak); // for BCs
      this->local_matrix_sys.template at<0,2>() = this->matrix_l.convert_to_1();
      this->local_matrix_sys.template at<1,0>() = this->matrix_d.local().clone(LAFEM::CloneMode::Shallow);
      this->local_matrix_sys.template at<1,1>() = this->null_matrix_pp.clone(LAFEM::CloneMode::Shallow);
      this->local_matrix_sys.template at<1,2>() = this->null_matrix_ps.clone(LAFEM::CloneMode::Shallow);
      this->local_matrix_sys.template at<2,0>() = this->matrix_k.convert_to_1();
      this->local_matrix_sys.template at<2,1>() = this->null_matrix_sp.clone(LAFEM::CloneMode::Shallow);
      this->local_matrix_sys.template at<2,2>() = this->matrix_m.convert_to_1();
      this->filter_velo.local().filter_mat(this->local_matrix_sys.template at<0,0>());
      this->filter_velo.local().filter_offdiag_row_mat(this->local_matrix_sys.template at<0,1>());
      this->filter_velo.local().filter_offdiag_row_mat(this->local_matrix_sys.template at<0,2>());
    }

    template<typename SpaceVelo_, typename SpaceStress_>
    void assemble_matrices(const SpaceVelo_& space_velo, const SpaceStress_& space_stress,
      const Cubature::DynamicFactory& cubature, const DataType_ nu, const DataType_ eta, const DataType_ chi)
    {
      this->matrix_a.local().format();
      this->matrix_m.local().format();
      this->matrix_k.local().format();
      this->matrix_l.local().format();

      Assembly::Common::DuDvOperatorBlocked<dim_> dudv_op;
      Assembly::BilinearOperatorAssembler::assemble_block_matrix1(this->matrix_a.local(), dudv_op, space_velo, cubature, nu);

      Assembly::Common::IdentityOperatorBlocked<nsc_> id_op;
      Assembly::BilinearOperatorAssembler::assemble_block_matrix1(this->matrix_m.local(), id_op, space_stress, cubature);

      Assembly::Common::StrainRateTensorOperator<dim_, nsc_> du_op;
      Assembly::BilinearOperatorAssembler::assemble_block_matrix2(this->matrix_k.local(), du_op, space_stress, space_velo, cubature, -eta);

      Assembly::Common::StressDivergenceOperator<dim_, nsc_> div_op;
      Assembly::BilinearOperatorAssembler::assemble_block_matrix2(this->matrix_l.local(), div_op, space_velo, space_stress, cubature, -chi);
    }
  };


  // 2D version with 4 stress components
  template<typename Mesh_, typename DT_, typename IT_>
  void add_stress_to_vtk(Geometry::ExportVTK<Mesh_>& exp, const LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, 4>& vs)
  {
    const std::size_t n = vs.size();
    std::vector<double> s11(n), s22(n), s12(n), s21(n);
    const auto* s = vs.elements();
    for(std::size_t i(0); i < n; ++i)
    {
      s11[i] = s[i][0];
      s12[i] = s[i][1];
      s21[i] = s[i][2];
      s22[i] = s[i][3];
    }
    exp.add_vertex_scalar("sigma_11", s11.data());
    exp.add_vertex_scalar("sigma_12", s12.data());
    exp.add_vertex_scalar("sigma_21", s21.data());
    exp.add_vertex_scalar("sigma_22", s22.data());
  }

  // 2D version with 3 stress components
  template<typename Mesh_, typename DT_, typename IT_>
  void add_stress_to_vtk(Geometry::ExportVTK<Mesh_>& exp, const LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, 3>& vs)
  {
    const std::size_t n = vs.size();
    std::vector<double> s11(n), s22(n), s12(n);
    const auto* s = vs.elements();
    for(std::size_t i(0); i < n; ++i)
    {
      s11[i] = s[i][0];
      s22[i] = s[i][1];
      s12[i] = s[i][2];
    }
    exp.add_vertex_scalar("sigma_11", s11.data());
    exp.add_vertex_scalar("sigma_22", s22.data());
    exp.add_vertex_scalar("sigma_12", s12.data());
  }

  // 3D version with 6 stress components
  template<typename Mesh_, typename DT_, typename IT_>
  void add_stress_to_vtk(Geometry::ExportVTK<Mesh_>& exp, const LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, 6>& vs)
  {
    const std::size_t n = vs.size();
    std::vector<double> s11(n), s22(n), s33(n), s12(n), s23(n), s31(n);
    const auto* s = vs.elements();
    for(std::size_t i(0); i < n; ++i)
    {
      s11[i] = s[i][0];
      s22[i] = s[i][1];
      s33[i] = s[i][2];
      s12[i] = s[i][3];
      s23[i] = s[i][4];
      s31[i] = s[i][5];
    }
    exp.add_vertex_scalar("sigma_11", s11.data());
    exp.add_vertex_scalar("sigma_22", s22.data());
    exp.add_vertex_scalar("sigma_33", s33.data());
    exp.add_vertex_scalar("sigma_12", s12.data());
    exp.add_vertex_scalar("sigma_23", s23.data());
    exp.add_vertex_scalar("sigma_31", s31.data());
  }

  template<typename Matrix_, typename Filter_>
  class S3FRichardson :
    public Solver::PreconditionedIterativeSolver<typename Matrix_::VectorTypeR>
  {
  public:
    typedef Matrix_ MatrixType;
    typedef Filter_ FilterType;
    typedef typename MatrixType::VectorTypeR VectorType;
    typedef typename MatrixType::DataType DataType;
    typedef Solver::PreconditionedIterativeSolver<VectorType> BaseClass;

  protected:
    const MatrixType& _system_matrix;
    const FilterType& _system_filter;
    VectorType _vec_def;
    VectorType _vec_cor;
    DataType _def_v, _def_p, _def_s;

  public:
    explicit S3FRichardson(const MatrixType& matrix, const FilterType& filter, std::shared_ptr<Solver::SolverBase<VectorType>> precond) :
      BaseClass("S3FRichardson", precond),
      _system_matrix(matrix),
      _system_filter(filter)
    {
      this->_set_comm_by_matrix(matrix);
    }

    virtual String name() const override
    {
      return "S3FRichardson";
    }

    virtual void init_symbolic() override
    {
      BaseClass::init_symbolic();
      // create two temporary vectors
      _vec_def = this->_system_matrix.create_vector_r();
      _vec_cor = this->_system_matrix.create_vector_r();
    }

    virtual void done_symbolic() override
    {
      this->_vec_cor.clear();
      this->_vec_def.clear();
      BaseClass::done_symbolic();
    }

    virtual Solver::Status apply(VectorType& vec_cor, const VectorType& vec_def) override
    {
      // save defect
      this->_vec_def.copy(vec_def);

      // clear solution vector
      vec_cor.format();

      // apply
      this->_status = _apply_intern(vec_cor, vec_def);
      this->plot_summary();
      return this->_status;
    }

    virtual Solver::Status correct(VectorType& vec_sol, const VectorType& vec_rhs) override
    {
      // compute defect
      this->_system_matrix.apply(this->_vec_def, vec_sol, vec_rhs, -DataType(1));
      this->_system_filter.filter_def(this->_vec_def);

      // apply
      this->_status = _apply_intern(vec_sol, vec_rhs);
      this->plot_summary();
      return this->_status;
    }

  protected:
    virtual DataType _calc_def_norm(const VectorType& vec_def, const VectorType& DOXY(vec_sol)) override
    {
      // compute defect component norms individually
      const auto& gate = *vec_def.get_gate();
      _def_v = gate.norm2(vec_def.local().template at<0>().norm2());
      _def_p = gate.norm2(vec_def.local().template at<1>().norm2());
      _def_s = gate.norm2(vec_def.local().template at<2>().norm2());
      return Math::sqrt(_def_v*_def_v + _def_p*_def_p + _def_s*_def_s);
    }

    virtual void _plot_iter_line(Index num_iter, DataType def_cur, DataType def_prev) override
    {
      // compose message line
      String msg = this->_plot_name
        +  ": " + stringify(num_iter).pad_front(this->_iter_digits)
        + " : " + stringify_fp_sci(def_cur);

      // print defect component norms
      msg += " [ ";
      msg += stringify_fp_sci(_def_v);
      msg += "  ";
      msg += stringify_fp_sci(_def_p);
      msg += "  ";
      msg += stringify_fp_sci(_def_s);
      msg += " ]";

      // not first iteration?
      if(num_iter > Index(0))
      {
        msg += " / " + stringify_fp_sci(def_cur / this->_def_init);
        msg += " / " + stringify_fp_fix(def_cur / def_prev);
      }

      this->_print_line(msg);
    }

    virtual Solver::Status _apply_intern(VectorType& vec_sol, const VectorType& vec_rhs)
    {
      Solver::IterationStats pre_iter(*this);
      Statistics::add_solver_expression(std::make_shared<Solver::ExpressionStartSolve>(this->name()));

      VectorType& vec_def(this->_vec_def);
      VectorType& vec_cor(this->_vec_cor);
      const MatrixType& matrix(this->_system_matrix);
      const FilterType& filter(this->_system_filter);

      // compute initial defect
      Solver::Status status = this->_set_initial_defect(vec_def, vec_sol);

      pre_iter.destroy();

      // start iterating
      while(status == Solver::Status::progress)
      {
        Solver::IterationStats stat(*this);

        // apply preconditioner
        if(!this->_apply_precond(vec_cor, vec_def, filter))
        {
          Statistics::add_solver_expression(std::make_shared<Solver::ExpressionEndSolve>(this->name(), Solver::Status::aborted, this->get_num_iter()));
          return Solver::Status::aborted;
        }
        //filter.filter_cor(vec_cor);

        // update solution vector
        vec_sol.axpy(vec_cor, vec_sol);

        // compute new defect vector
        //matrix.apply(vec_def, vec_sol, vec_rhs, -DataType(1));
        {
          vec_def.copy(vec_rhs);
          auto& d = vec_def.local();
          auto& x = vec_sol.local();
          auto& a = matrix.local();

        //a.template at<0,0>().apply(d.template at<0>(), x.template at<0>(), d.template at<0>(), -DataType(1)); // A*u
          a.template at<0,1>().apply(d.template at<0>(), x.template at<1>(), d.template at<0>(), -DataType(1)); // B*p
          a.template at<0,2>().apply(d.template at<0>(), x.template at<2>(), d.template at<0>(), -DataType(1)); // L*s
          a.template at<1,0>().apply(d.template at<1>(), x.template at<0>(), d.template at<1>(), -DataType(1)); // D*u
          a.template at<2,0>().apply(d.template at<2>(), x.template at<0>(), d.template at<2>(), -DataType(1)); // K*u
          a.template at<2,2>().apply(d.template at<2>(), x.template at<2>(), d.template at<2>(), -DataType(1)); // M*s

          vec_def.sync_0();
        }
        filter.filter_def(vec_def);

        // compute new defect norm
        status = this->_set_new_defect(vec_def, vec_sol);
      }

      // return our status
      Statistics::add_solver_expression(std::make_shared<Solver::ExpressionEndSolve>(this->name(), status, this->get_num_iter()));
      return status;
    }
  };


  struct Counts
  {
    static constexpr std::size_t num_ranks = 0u;
    static constexpr std::size_t num_elems = 1u;
    static constexpr std::size_t num_dofs_l_v = 2u;
    static constexpr std::size_t num_dofs_l_p = 3u;
    static constexpr std::size_t num_dofs_l_s = 4u;
    static constexpr std::size_t num_dofs_l = 5u;
    static constexpr std::size_t num_dofs_g_v = 6u;
    static constexpr std::size_t num_dofs_g_p = 7u;
    static constexpr std::size_t num_dofs_g_s = 8u;
    static constexpr std::size_t num_dofs_g = 9u;
    static constexpr std::size_t num_nze_a = 10u;
    static constexpr std::size_t num_nze_b = 11u;
    static constexpr std::size_t num_nze_d = 12u;
    static constexpr std::size_t num_nze_k = 13u;
    static constexpr std::size_t num_nze_l = 14u;
    static constexpr std::size_t num_nze_m = 15u;
    static constexpr std::size_t num_nze = 16u;
    static constexpr std::size_t bytes_domain = 17u;
    static constexpr std::size_t bytes_system = 18u;
    static constexpr std::size_t elems_mirror = 19u;
    static constexpr std::size_t elems_mirror_v = 20u;
    static constexpr std::size_t elems_mirror_p = 21u;
    static constexpr std::size_t elems_mirror_s = 22u;
    static constexpr std::size_t bytes_matrix = 23u;
    static constexpr std::size_t bytes_vanka = 24u;
    static constexpr std::size_t num_nze_vanka = 25u;
    static constexpr std::size_t count = 26u;
  };

  struct Times
  {
    static constexpr std::size_t asm_total = 0u;
    static constexpr std::size_t asm_gate = 1u;
    static constexpr std::size_t asm_muxer = 2u;
    static constexpr std::size_t asm_transfer = 3u;
    static constexpr std::size_t asm_matrix = 4u;
    static constexpr std::size_t vanka_init = 5u;
    static constexpr std::size_t vanka_apply = 6u;
    static constexpr std::size_t count = 7u;
  };

  template<typename T_>
  inline T_ sum(const std::vector<T_>& dv)
  {
    T_ t = T_(0);
    for(auto x : dv) t += x;
    return t;
  }

  template<typename T_, std::size_t n_>
  inline T_ sum(const std::vector<std::array<T_, n_>>& dv, std::size_t i)
  {
    XASSERT(i < n_);
    T_ t = T_(0);
    for(const auto& x : dv) t += x[i];
    return t;
  }

  struct BenchStats
  {
    // per [level][index]
    std::vector<std::array<unsigned long long, Counts::count>> counts, counts_sum, counts_max;

    // per [level][index]
    std::vector<std::array<double, Times::count>> times;

    // (physical, virtual)
    std::array<unsigned long long, 2> mem_use, mem_use_sum, mem_use_max;

    double toe_asm_rhs;

    explicit BenchStats(std::size_t virt_size) :
      counts(virt_size),
      counts_sum(virt_size),
      counts_max(virt_size),
      times(virt_size),
      toe_asm_rhs(0.0)
    {
      for(std::size_t i(0u); i < virt_size; ++i)
      {
        for(std::size_t j(0u); j < Counts::count; ++j)
          counts[i][j] = 0ull;
        for(std::size_t j(0u); j < Times::count; ++j)
          times[i][j] = 0.0;
      }
    }

    void sync(const Dist::Comm& comm)
    {
      comm.allreduce(counts.data(), counts_sum.data(), counts.size()*Counts::count, Dist::dt_unsigned_long_long, Dist::op_sum);
      comm.allreduce(counts.data(), counts_max.data(), counts.size()*Counts::count, Dist::dt_unsigned_long_long, Dist::op_max);

      comm.allreduce(mem_use.data(), mem_use_sum.data(), 2u, Dist::op_sum);
      comm.allreduce(mem_use.data(), mem_use_max.data(), 2u, Dist::op_max);
    }

    String format() const
    {
      String s;

      s += "\nVanka Timinigs:\n";
      s += "                Init        Apply\n";
      s += "Overall : " +
        stringify_fp_fix(sum(times, Times::vanka_init), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::vanka_apply), 6, 10) + "\n";
      for(std::size_t i(0); i < times.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(times[i][Times::vanka_init], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::vanka_apply], 6, 10) + "\n";
      }

      s += "\nAssembly Timings:\n";
      s += "                Gate        Muxer     Transfer       Matrix        Total\n";
      s += "Overall : " +
        stringify_fp_fix(sum(times, Times::asm_gate), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::asm_muxer), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::asm_transfer), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::asm_matrix), 6, 10) + " / " +
        stringify_fp_fix(sum(times, Times::asm_total), 6, 10) + "\n";
      for(std::size_t i(0); i < times.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(times[i][Times::asm_gate], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::asm_muxer], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::asm_transfer], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::asm_matrix], 6, 10) + " / " +
          stringify_fp_fix(times[i][Times::asm_total], 6, 10) + "\n";
      }

      s += "\nMesh Statistics:\n";
      s += "           Ranks       Elements [    per Patch ]           Dofs [    per Patch ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify(counts_max[i][Counts::num_ranks]).pad_front(6) + " / " +
          stringify(counts_sum[i][Counts::num_elems]).pad_front(12) + " [ " +
          stringify(counts_max[i][Counts::num_elems]).pad_front(12) + " ] / " +
          stringify(counts    [i][Counts::num_dofs_g]).pad_front(12) + " [ " +
          stringify(counts_max[i][Counts::num_dofs_l]).pad_front(12) + " ]\n";
      }

      s += "\nDegree of Freedom Statistics:\n";
      s += "                 Total [  per Patch ]       Velocity [  per Patch ]       Pressure [  per Patch ]         Stress [  per Patch ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify(counts    [i][Counts::num_dofs_g]).pad_front(12) + " [ " +
          stringify(counts_max[i][Counts::num_dofs_l]).pad_front(10) + " ] / " +
          stringify(counts    [i][Counts::num_dofs_g_v]).pad_front(12) + " [ " +
          stringify(counts_max[i][Counts::num_dofs_l_v]).pad_front(10) + " ] / " +
          stringify(counts    [i][Counts::num_dofs_g_p]).pad_front(12) + " [ " +
          stringify(counts_max[i][Counts::num_dofs_l_p]).pad_front(10) + " ] / " +
          stringify(counts    [i][Counts::num_dofs_g_s]).pad_front(12) + " [ " +
          stringify(counts_max[i][Counts::num_dofs_l_s]).pad_front(10) + " ]\n";
      }

      s += "\nMirror #elems Statistics:\n";
      s += "               Total [  per Patch ]     Velocity [  per Patch ]     Pressure [  per Patch ]       Stress [  per Patch ]\n";
      s += "Overall : " +
        stringify(sum(counts_sum, Counts::elems_mirror)).pad_front(10) + " [ " +
        stringify(sum(counts_max, Counts::elems_mirror)).pad_front(10) + " ] / " +
        stringify(sum(counts_sum, Counts::elems_mirror_v)).pad_front(10) + " [ " +
        stringify(sum(counts_max, Counts::elems_mirror_v)).pad_front(10) + " ] / " +
        stringify(sum(counts_sum, Counts::elems_mirror_p)).pad_front(10) + " [ " +
        stringify(sum(counts_max, Counts::elems_mirror_p)).pad_front(10) + " ] / " +
        stringify(sum(counts_sum, Counts::elems_mirror_s)).pad_front(10) + " [ " +
        stringify(sum(counts_max, Counts::elems_mirror_s)).pad_front(10) + " ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify(counts_sum[i][Counts::elems_mirror]).pad_front(10) + " [ " +
          stringify(counts_max[i][Counts::elems_mirror]).pad_front(10) + " ] / " +
          stringify(counts_sum[i][Counts::elems_mirror_v]).pad_front(10) + " [ " +
          stringify(counts_max[i][Counts::elems_mirror_v]).pad_front(10) + " ] / " +
          stringify(counts_sum[i][Counts::elems_mirror_p]).pad_front(10) + " [ " +
          stringify(counts_max[i][Counts::elems_mirror_p]).pad_front(10) + " ] / " +
          stringify(counts_sum[i][Counts::elems_mirror_s]).pad_front(10) + " [ " +
          stringify(counts_max[i][Counts::elems_mirror_s]).pad_front(10) + " ]\n";
      }

      s += "\nMemory Usage Statistics:\n";
      s += String("Peak Physical") +
        ": " + stringify_fp_fix(double(mem_use_sum[0])/(1024.0*1024.0*1024.0), 6, 15) + " GiB " +
        "[ " + stringify_fp_fix(double(mem_use_max[0])/(1024.0*1024.0*1024.0), 6, 15) + " GiB ]\n";
        //": " + stringify(mem_use_sum[0]).pad_front(15) +
        //" [" + stringify(mem_use_max[0]).pad_front(15) + " ]\n";
      s += String("Peak Virtual.") +
        ": " + stringify_fp_fix(double(mem_use_sum[1])/(1024.0*1024.0*1024.0), 6, 15) + " GiB " +
        "[ " + stringify_fp_fix(double(mem_use_max[1])/(1024.0*1024.0*1024.0), 6, 15) + " GiB ]\n";
        //": " + stringify(mem_use_sum[1]).pad_front(15) +
        //" [" + stringify(mem_use_max[1]).pad_front(15) + " ]\n";

      s += "\nDomain Level Bytes Statistics:\n";
      s += "Overall : " + stringify(sum(counts_sum, Counts::bytes_domain)).pad_front(15) +
        " [ " + stringify(sum(counts_max, Counts::bytes_domain)).pad_front(15) + " ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) +
          ": " + stringify(counts_sum[i][Counts::bytes_domain]).pad_front(15) +
          " [ " + stringify(counts_max[i][Counts::bytes_domain]).pad_front(15) + " ]\n";
      }

      s += "\nSystem Level Bytes Statistics:\n";
      s += "                   System [    per Patch ]          Matrix [    per Patch ]           Vanka [    per Patch ]\n";
      s += "Overall : " +
        stringify(sum(counts_sum, Counts::bytes_system)).pad_front(15) + " [ " +
        stringify(sum(counts_max, Counts::bytes_system)).pad_front(12) + " ] " +
        stringify(sum(counts_sum, Counts::bytes_matrix)).pad_front(15) + " [ " +
        stringify(sum(counts_max, Counts::bytes_matrix)).pad_front(12) + " ] " +
        stringify(sum(counts_sum, Counts::bytes_vanka)).pad_front(15) + " [ " +
        stringify(sum(counts_max, Counts::bytes_vanka)).pad_front(12) + " ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify(counts_sum[i][Counts::bytes_system]).pad_front(15) + " [ " +
          stringify(counts_max[i][Counts::bytes_system]).pad_front(12) + " ] " +
          stringify(counts_sum[i][Counts::bytes_matrix]).pad_front(15) + " [ " +
          stringify(counts_max[i][Counts::bytes_matrix]).pad_front(12) + " ] " +
          stringify(counts_sum[i][Counts::bytes_vanka]).pad_front(15) + " [ " +
          stringify(counts_max[i][Counts::bytes_vanka]).pad_front(12) + " ]\n";
      }

      s += "\nSystem Level Non-Zero Entries Statistics:\n";
      s += "                 Matrix [  per Patch ]         Vanka [  per Patch ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify(counts_sum[i][Counts::num_nze]).pad_front(13) + " [ " +
          stringify(counts_max[i][Counts::num_nze]).pad_front(10) + " ] " +
          stringify(counts_sum[i][Counts::num_nze_vanka]).pad_front(13) + " [ " +
          stringify(counts_max[i][Counts::num_nze_vanka]).pad_front(10) + " ]\n";
      }

      s += "\nMatrix Block Non-Zero Statistics:\n";
      s += "                      A [  per Patch ]           B/D [  per Patch ]             M [  per Patch ]           K/L [  per Patch ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify(counts_sum[i][Counts::num_nze_a]).pad_front(13) + " [ " +
          stringify(counts_max[i][Counts::num_nze_a]).pad_front(10) + " ] " +
          stringify(counts_sum[i][Counts::num_nze_b]).pad_front(13) + " [ " +
          stringify(counts_max[i][Counts::num_nze_b]).pad_front(10) + " ] " +
          stringify(counts_sum[i][Counts::num_nze_m]).pad_front(13) + " [ " +
          stringify(counts_max[i][Counts::num_nze_m]).pad_front(10) + " ] " +
          stringify(counts_sum[i][Counts::num_nze_k]).pad_front(13) + " [ " +
          stringify(counts_max[i][Counts::num_nze_k]).pad_front(10) + " ]\n";
      }


      return s;
    }
  };


  template<typename DomainLevel_>
  void run(SimpleArgParser& args, Control::Domain::DomainControl<DomainLevel_>& domain)
  {
    // get our main communicator
    const Dist::Comm& comm = domain.comm();

    // define our arch types
    typedef Mem::Main MemType;
    typedef double DataType;
    typedef Index IndexType;

    // define our domain type
    typedef Control::Domain::DomainControl<DomainLevel_> DomainControlType;
    typedef typename DomainControlType::LevelType DomainLevelType;

    // fetch our mesh type
    typedef typename DomainControlType::MeshType MeshType;
    typedef typename MeshType::ShapeType ShapeType;
    static constexpr int dim = ShapeType::dimension;
    static constexpr int nsc = (dim*(dim+1))/2;

    // define our system level
    typedef MyStokes3FieldSystemLevel<dim, nsc, MemType, DataType, IndexType> SystemLevelType;


    BenchStats stats(domain.size_virtual());


    {
      std::size_t n = domain.size_virtual() > domain.size_physical() ? domain.size_physical()+1 : domain.size_physical();
      for(std::size_t i(0); i < n; ++i)
      {
        const auto& vdl = domain.at(i);
        if(vdl.is_parent())
          stats.counts[i][Counts::bytes_domain] = vdl.level_c().bytes() + vdl.level_p().bytes();
        else if(vdl.is_child())
          stats.counts[i][Counts::bytes_domain] = vdl.level_c().bytes();
        else
          stats.counts[i][Counts::bytes_domain] = vdl.level().bytes();
      }
    }


    /* ****************************************************************************************** */

    DataType nu = 1.0;
    DataType eta = 1.0;
    DataType chi = 1.0;
    DataType v_max = 1.0;
    Index max_iter = 5;    // max. multigrid iterations
    Index smooth_steps = 4;
    DataType smooth_damp = 0.7;
    DataType linsol_tol = 1E-8; // rel. tolerance for linear solver
    Index gmres_k = 50;         // k of GMRES(k) coarse grid solver

    args.parse("nu", nu);
    args.parse("eta", eta);
    args.parse("chi", chi);
    args.parse("v-max", v_max);
    args.parse("max-iter", max_iter);
    args.parse("smooth-steps", smooth_steps);
    args.parse("smooth-damp", smooth_damp);
    args.parse("linsol-tol", linsol_tol);
    args.parse("gmresk", gmres_k);

#ifdef FEAT_HAVE_UMFPACK
    const bool umf_cgs = (domain.back_layer().comm().size() == 1);
#else
    const bool umf_cgs = false;
#endif

    {
      static constexpr std::size_t pl = 20u;
      static constexpr char pc = '.';
      comm.print("\nProblem Parameters:");
      comm.print(String("DIM").pad_back(pl, pc) + ": " + stringify(dim));
      comm.print(String("NSC").pad_back(pl, pc) + ": " + stringify(nsc));
      comm.print(String("Nu").pad_back(pl, pc) + ": " + stringify(nu));
      comm.print(String("Eta").pad_back(pl, pc) + ": " + stringify(eta));
      comm.print(String("Chi").pad_back(pl, pc) + ": " + stringify(chi));
      comm.print(String("V-Max").pad_back(pl, pc) + ": " + stringify(v_max));
      comm.print(String("Smoothing Steps").pad_back(pl, pc) + ": " + stringify(smooth_steps));
      comm.print(String("Vanka Damping").pad_back(pl, pc) + ": " + stringify(smooth_damp));
      comm.print(String("LinSol RelTol").pad_back(pl, pc) + ": " + stringify_fp_sci(linsol_tol));
      comm.print(String("Max MG Iter").pad_back(pl, pc) + ": " + stringify(max_iter));
      if(umf_cgs)
        comm.print(String("Coarse Solver").pad_back(pl, pc) + ": UMFPACK");
      else
        comm.print(String("Coarse Solver").pad_back(pl, pc) + ": GMRES(" + stringify(gmres_k) + ")");
    }

    /* ****************************************************************************************** */

    std::deque<std::shared_ptr<SystemLevelType>> system_levels;

    const Index num_levels = Index(domain.size_physical());

    // create system levels
    for (Index i(0); i < num_levels; ++i)
    {
      system_levels.push_back(std::make_shared<SystemLevelType>());
    }

    Cubature::DynamicFactory cubature("gauss-legendre:4");

    /* ***************************************************************************************** */

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      system_levels.at(i)->assemble_gates(domain.at(i));
      stats.times[i][Times::asm_gate] += ts.elapsed_now();

      if((i+1) < domain.size_virtual())
      {
        TimeStamp ts2;
        system_levels.at(i)->assemble_coarse_muxers(domain.at(i+1));
        stats.times[i][Times::asm_muxer] += ts2.elapsed_now();
        TimeStamp ts3;
        system_levels.at(i)->assemble_transfers(domain.at(i), domain.at(i+1), cubature);
        stats.times[i][Times::asm_transfer] += ts3.elapsed_now();
      }
      stats.times[i][Times::asm_total] += ts.elapsed_now();
    }

    /* ***************************************************************************************** */

    for(Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      system_levels.at(i)->assemble_structs(domain.at(i)->space_velo, domain.at(i)->space_pres, domain.at(i)->space_stress);
      system_levels.at(i)->assemble_grad_div_matrices(domain.at(i)->space_velo, domain.at(i)->space_pres, cubature);
      system_levels.at(i)->assemble_matrices(domain.at(i)->space_velo, domain.at(i)->space_stress, cubature, nu, eta, chi);
      system_levels.at(i)->compile_system_matrix();
      double tt = ts.elapsed_now();
      stats.times[i][Times::asm_total] += tt;
      stats.times[i][Times::asm_matrix] += tt;
    }

    /* ***************************************************************************************** */

    // our inflow BC function
    InflowFunction<dim> inflow_func(v_max);

    // the names of the mesh parts on which to assemble
    std::deque<String> part_names = domain.front()->get_mesh_node()->get_mesh_part_names(true);

    for(Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      // get our local velocity filter
      auto& fil_loc_v = system_levels.at(i)->filter_velo.local();

      // create unit-filter assembler
      Assembly::UnitFilterAssembler<MeshType> unit_asm_inflow, unit_asm_noflow;

      // loop over all boundary parts except for the right one, which is outflow
      for(const auto& name : part_names)
      {
        // skip non-boundary mesh-parts
        if(!name.starts_with("bnd:"))
          continue;

        // try to fetch the corresponding mesh part node
        auto* mesh_part_node = domain.at(i)->get_mesh_node()->find_mesh_part_node(name);
        XASSERT(mesh_part_node != nullptr);

        auto* mesh_part = mesh_part_node->get_mesh();
        if (mesh_part != nullptr)
        {
          if(name == "bnd:l")
          {
            // inflow
            unit_asm_inflow.add_mesh_part(*mesh_part);
          }
          else if(name != "bnd:r")
          {
            // outflow
            unit_asm_noflow.add_mesh_part(*mesh_part);
          }
        }
      }

      // assemble the filters
      unit_asm_inflow.assemble(fil_loc_v, domain.at(i)->space_velo, inflow_func);
      unit_asm_noflow.assemble(fil_loc_v, domain.at(i)->space_velo);

      // finally, compile the system filter
      system_levels.at(i)->compile_system_filter();
      stats.times[i][Times::asm_total] += ts.elapsed_now();
    }

    // finally, compile the local type-1 matrices
    for(Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      system_levels.at(i)->compile_local_matrix();
      stats.times[i][Times::asm_total] += ts.elapsed_now();
    }

    for(Index i(0); i < num_levels; ++i)
    {
      stats.counts[i][Counts::num_ranks] = Index(domain.at(i).layer().comm().size());
      stats.counts[i][Counts::num_elems] = domain.at(i)->get_mesh().get_num_elements();
      stats.counts[i][Counts::num_dofs_g] = system_levels.at(i)->matrix_sys.rows();
      stats.counts[i][Counts::num_dofs_l] = system_levels.at(i)->matrix_sys.local().rows();
      stats.counts[i][Counts::num_dofs_g_v] = system_levels.at(i)->matrix_a.rows();
      stats.counts[i][Counts::num_dofs_l_v] = system_levels.at(i)->matrix_a.local().rows();
      stats.counts[i][Counts::num_dofs_g_p] = system_levels.at(i)->matrix_d.rows();
      stats.counts[i][Counts::num_dofs_l_p] = system_levels.at(i)->matrix_d.local().rows();
      stats.counts[i][Counts::num_dofs_g_s] = system_levels.at(i)->matrix_m.rows();
      stats.counts[i][Counts::num_dofs_l_s] = system_levels.at(i)->matrix_m.local().rows();
      stats.counts[i][Counts::num_nze] = system_levels.at(i)->matrix_sys.local().template used_elements<LAFEM::Perspective::pod>();
      stats.counts[i][Counts::num_nze_a] = system_levels.at(i)->matrix_a.local().used_elements();
      stats.counts[i][Counts::num_nze_b] = system_levels.at(i)->matrix_b.local().used_elements();
      stats.counts[i][Counts::num_nze_d] = system_levels.at(i)->matrix_d.local().used_elements();
      stats.counts[i][Counts::num_nze_k] = system_levels.at(i)->matrix_k.local().used_elements();
      stats.counts[i][Counts::num_nze_l] = system_levels.at(i)->matrix_l.local().used_elements();
      stats.counts[i][Counts::num_nze_m] = system_levels.at(i)->matrix_m.local().used_elements();
      stats.counts[i][Counts::bytes_system] = system_levels.at(i)->bytes();
      stats.counts[i][Counts::bytes_matrix] = system_levels.at(i)->matrix_sys.local().bytes();
      stats.counts[i][Counts::elems_mirror] = 0;
      auto& sys_mirrors =  system_levels.at(i)->gate_sys._mirrors;
      for (auto& mirror : sys_mirrors)
      {
        stats.counts[i][Counts::elems_mirror_v] += mirror.template at<0>().num_indices();
        stats.counts[i][Counts::elems_mirror_p] += mirror.template at<1>().num_indices();
        stats.counts[i][Counts::elems_mirror_s] += mirror.template at<2>().num_indices();
        stats.counts[i][Counts::elems_mirror] += mirror.template at<0>().num_indices()
          + mirror.template at<1>().num_indices() + mirror.template at<2>().num_indices();
      }
    }

    /* ***************************************************************************************** */

    // get our global solver system types
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;
    typedef typename SystemLevelType::GlobalSystemMatrix GlobalSystemMatrix;
    typedef typename SystemLevelType::GlobalSystemFilter GlobalSystemFilter;
    typedef typename SystemLevelType::GlobalSystemTransfer GlobalSystemTransfer;
    //typedef typename SystemLevelType::GlobalVeloVector GlobalVeloVector;
    //typedef typename SystemLevelType::GlobalPresVector GlobalPresVector;

    // fetch our finest levels
    DomainLevelType& the_domain_level = *domain.front();
    SystemLevelType& the_system_level = *system_levels.front();

    // get our global solve matrix and filter
    GlobalSystemMatrix& matrix = the_system_level.matrix_sys;
    GlobalSystemFilter& filter = the_system_level.filter_sys;

    // create new vectors
    GlobalSystemVector vec_sol = the_system_level.matrix_sys.create_vector_r();
    GlobalSystemVector vec_rhs = the_system_level.matrix_sys.create_vector_r();

    // format the vectors
    vec_sol.format();
    vec_rhs.format();

    // and filter them
    the_system_level.filter_sys.filter_sol(vec_sol);
    the_system_level.filter_sys.filter_rhs(vec_rhs);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // create a multigrid solver
    auto multigrid_hierarchy = std::make_shared<
      Solver::MultiGridHierarchy<GlobalSystemMatrix, GlobalSystemFilter, GlobalSystemTransfer>>(domain.size_virtual());

    // array of Vankas
    std::deque<
      std::shared_ptr<
        Solver::AmaVanka<
          typename SystemLevelType::LocalSystemMatrix,
          typename SystemLevelType::LocalSystemFilter>>> ama_vankas;


    // push levels into multigrid
    for(std::size_t i(0); i < system_levels.size(); ++i)
    {
      DomainLevelType& dom_lvl = *domain.at(i);
      SystemLevelType& lvl = *system_levels.at(i);

      // Create Vanka
      auto vanka = Solver::new_amavanka(lvl.local_matrix_sys, lvl.filter_sys.local());

      vanka->push_macro_dofs(Space::DofMappingRenderer::render(dom_lvl.space_velo));
      vanka->push_macro_dofs(Space::DofMappingRenderer::render(dom_lvl.space_pres));
      vanka->push_macro_dofs(Space::DofMappingRenderer::render(dom_lvl.space_stress));

      ama_vankas.push_back(vanka);

      // create schwarz
      auto schwarz = Solver::new_schwarz_precond(vanka, lvl.filter_sys);

      if((i+1) < domain.size_virtual())
      {
        // create richardson smoother
        auto smoother = Solver::new_richardson(lvl.matrix_sys, lvl.filter_sys, smooth_damp, schwarz);
        //auto smoother = std::make_shared<MyRichardson<GlobalSystemMatrix, GlobalSystemFilter>>(lvl.matrix_sys, lvl.filter_sys, smooth_damp, schwarz);

        smoother->set_min_iter(smooth_steps);
        smoother->set_max_iter(smooth_steps);
        smoother->skip_defect_calc(true); // skip defect calculation

        //smoother->set_plot_name(String("SM:") + stringify(i));
        //smoother->set_plot_mode(Solver::PlotMode::iter);

        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);
      }
#ifdef FEAT_HAVE_UMFPACK
      else if(umf_cgs)
      {
        // create UMFPACK coarse grid solver
        auto umfpack = Solver::new_generic_umfpack(lvl.local_matrix_sys);
        auto cgsolver = Solver::new_schwarz_precond(umfpack, lvl.filter_sys);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
        ama_vankas.pop_back();
      }
#endif //  FEAT_HAVE_UMFPACK
      else
      {
        // create coarse grid solver
        auto cgsolver = Solver::new_fgmres(lvl.matrix_sys, lvl.filter_sys, gmres_k, 0.0, schwarz);
        cgsolver->set_max_iter(1000);
        cgsolver->set_tol_rel(1e-3);
        //cgsolver->set_plot_mode(Solver::PlotMode::summary);

        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);
      }
    }

    // create our multigrid solver
    auto multigrid = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V);

    // create our solver
    //auto solver = Solver::new_richardson(matrix, filter, 1.0, multigrid);
    auto solver = std::make_shared<S3FRichardson<GlobalSystemMatrix, GlobalSystemFilter>>(matrix, filter, multigrid);

    //auto solver = Solver::new_fgmres(matrix, filter, 20, 0.0, multigrid);
    //auto solver = Solver::new_bicgstab(matrix, filter, multigrid);

    solver->set_plot_name("Multigrid");

    solver->set_max_iter(max_iter);
    solver->set_tol_rel(linsol_tol);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    comm.print("");

    // initialize solver
    multigrid_hierarchy->init();
    solver->init();

    // solve Stokes system
    solver->set_plot_mode(Solver::PlotMode::iter);

    TimeStamp at;
    solver->correct(vec_sol, vec_rhs);
    const double solver_toe(at.elapsed_now());

    // store vanka times
    for(std::size_t i(0); i < ama_vankas.size(); ++i)
    {
      stats.times[i][Times::vanka_init] += ama_vankas.at(i)->time_init_symbolic();
      stats.times[i][Times::vanka_init] += ama_vankas.at(i)->time_init_numeric();
      stats.times[i][Times::vanka_apply] += ama_vankas.at(i)->time_apply();
      stats.counts[i][Counts::bytes_vanka] += ama_vankas.at(i)->bytes();
      stats.counts[i][Counts::num_nze_vanka] += ama_vankas.at(i)->data_size();
    }

    // release solvers
    solver->done_numeric();
    multigrid_hierarchy->done_numeric();

    comm.print("");
    comm.print(solver->get_summary());

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    // get memory info
    {
      MemoryUsage meminfo;
      stats.mem_use[0] = meminfo.get_peak_physical();
      stats.mem_use[1] = meminfo.get_peak_virtual();
    }

    stats.sync(comm);


    comm.print("\nMultigrid Timings:");
    comm.print("                Defect /     Smoother /     Transfer /       Coarse /        Total");
    comm.print("Overall : " +
        stringify_fp_fix(multigrid_hierarchy->get_time_defect(), 6, 12) + " / " +
        stringify_fp_fix(multigrid_hierarchy->get_time_smooth(), 6, 12) + " / " +
        stringify_fp_fix(multigrid_hierarchy->get_time_transfer(), 6, 12) + " / " +
        stringify_fp_fix(multigrid_hierarchy->get_time_coarse(), 6, 12) + " / " +
        stringify_fp_fix(multigrid_hierarchy->get_time_defect()+multigrid_hierarchy->get_time_smooth()
          +multigrid_hierarchy->get_time_transfer()+multigrid_hierarchy->get_time_coarse(), 6, 12));
    for(int i(0); i < int(multigrid_hierarchy->size_physical()); ++i)
    {
      comm.print("Level " + stringify(i).pad_front(2) + ": " +
        stringify_fp_fix(multigrid_hierarchy->get_time_defect(i), 6, 12) + " / " +
        stringify_fp_fix(multigrid_hierarchy->get_time_smooth(i), 6, 12) + " / " +
        stringify_fp_fix(multigrid_hierarchy->get_time_transfer(i), 6, 12) + " / " +
        stringify_fp_fix(multigrid_hierarchy->get_time_coarse(i), 6, 12) + " / " +
        stringify_fp_fix(multigrid_hierarchy->get_time_defect(i)+multigrid_hierarchy->get_time_smooth(i)
          +multigrid_hierarchy->get_time_transfer(i)+multigrid_hierarchy->get_time_coarse(i), 6, 12));
    }

    comm.print(stats.format());

    FEAT::Control::Statistics::report(solver_toe, 0, MeshType::ShapeType::dimension, system_levels, domain);

    comm.print(FEAT::Statistics::get_formatted_solver_internals());
    comm.print(FEAT::Statistics::get_formatted_solver_tree().trim());

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    //*
    if(args.check("vtk") >= 0)
    {
      // build VTK name
      String vtk_name;
      int npars = args.parse("vtk", vtk_name);
      if(npars < 1)
      {
        vtk_name = String("stokes-3field-") + stringify(dim) + "d";
        vtk_name += "-lvl" + stringify(the_domain_level.get_level_index());
        vtk_name += "-n" + stringify(comm.size());
      }

      comm.print("Writing '" + vtk_name + ".vtu...");

      // Create a VTK exporter for our mesh
      Geometry::ExportVTK<MeshType> exporter(the_domain_level.get_mesh());

      // project velocity
      exporter.add_vertex_vector("velocity", vec_sol.local().template at<0>());

      // project stress
      add_stress_to_vtk(exporter, vec_sol.local().template at<2>());

      // project pressure
      Cubature::DynamicFactory cub("gauss-legendre:2");
      LAFEM::DenseVector<Mem::Main, double, Index> vtx_p;
      Assembly::DiscreteCellProjector::project(vtx_p, vec_sol.local().template at<1>(), the_domain_level.space_pres, cub);

      // write pressure
      exporter.add_cell_scalar("pressure", vtx_p.elements());

      // finally, write the VTK file
      exporter.write(vtk_name, comm);
    }
  }

  template<int dim_>
  void run_dim(SimpleArgParser& args, Dist::Comm& comm, Geometry::MeshFileReader& mesh_reader)
  {
    // define our mesh type
    typedef Shape::Hypercube<dim_> ShapeType;
    typedef Geometry::ConformalMesh<ShapeType> MeshType;
    typedef Trafo::Standard::Mapping<MeshType> TrafoType;
    typedef Space::Lagrange2::Element<TrafoType> SpaceVeloType;
    typedef Space::Discontinuous::Element<TrafoType, Space::Discontinuous::Variant::StdPolyP<1>> SpacePresType;
    typedef Space::Lagrange3::Element<TrafoType> SpaceStressType;

    // create a time-stamp
    TimeStamp time_stamp;

    // create our domain control
    typedef Control::Domain::Stokes3FieldDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType, SpaceStressType> DomainLevelType;
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);

    domain.parse_args(args);
    domain.set_desired_levels(args.query("level")->second);
    domain.create(mesh_reader);

    // print partitioning info
    comm.print(domain.get_chosen_parti_info());

    // plot our levels
    comm.print("Desired Levels: " + domain.format_desired_levels());
    comm.print("Chosen  Levels: " + domain.format_chosen_levels());

    // run our application
    run(args, domain);

    // print elapsed runtime
    comm.print("Run-Time: " + time_stamp.elapsed_string_now());
  }

  void main(int argc, char* argv[])
  {
    // create world communicator
    Dist::Comm comm(Dist::Comm::world());

    // print number of processes
    comm.print("Number of Processes: " + stringify(comm.size()));

    // create arg parser
    SimpleArgParser args(argc, argv);

    // check command line arguments
    Control::Domain::add_supported_pdc_args(args);
    args.support("level");
    args.support("vtk");
    args.support("mesh");
    args.support("nu");
    args.support("eta");
    args.support("chi");
    args.support("max-iter");
    args.support("smooth-steps");
    args.support("smooth-damp");
    args.support("linsol-tol");
    args.support("v-max");
    args.support("gmresk");
    args.support("test-mode");

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
        comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

      // abort
      FEAT::Runtime::abort();
    }

    if(args.check("mesh") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--mesh <mesh-file>' is missing!");
      FEAT::Runtime::abort();
    }
    if(args.check("level") < 1)
    {
      comm.print(std::cerr, "ERROR: Mandatory option '--level <levels>' is missing!");
      FEAT::Runtime::abort();
    }

    // Our mesh file reader
    Geometry::MeshFileReader mesh_reader;

    // read in the mesh files
    mesh_reader.add_mesh_files(comm, args.query("mesh")->second);

    // read the mesh file root markups
    mesh_reader.read_root_markup();
    String mesh_type = mesh_reader.get_meshtype_string();

    // run 2D or 3D ?
    if(mesh_type == "conformal:hypercube:2:2")
      run_dim<2>(args, comm, mesh_reader);
    else if(mesh_type == "conformal:hypercube:3:3")
      run_dim<3>(args, comm, mesh_reader);
    else
    {
      comm.print(std::cerr, "ERROR: unsupported mesh type '" + mesh_type + "'");
      FEAT::Runtime::abort();
    }
  }
} // namespace Stokes3Field

int main(int argc, char* argv[])
{
  FEAT::Runtime::initialize(argc, argv);
  try
  {
    Stokes3Field::main(argc, argv);
  }
  catch(std::exception& e)
  {
    std::cerr << "ERROR: " << e.what() << std::endl;
    FEAT::Runtime::abort();
  }
  catch (...)
  {
    std::cerr << "ERROR: unknown exception" << std::endl;
    FEAT::Runtime::abort();
  }
  return FEAT::Runtime::finalize();
}
