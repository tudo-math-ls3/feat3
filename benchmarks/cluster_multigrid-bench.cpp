// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/hash.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/analytic/common.hpp>
#include <kernel/solver/richardson.hpp>
#include <kernel/solver/pcg.hpp>
#include <kernel/solver/multigrid.hpp>
#include <kernel/solver/jacobi_precond.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/scalar_basic.hpp>
#include <control/statistics.hpp>

namespace ClusterMultigridBench
{
  using namespace FEAT;

  typedef double DataType;
  typedef std::uint32_t IndexType;

  static_assert(sizeof(DataType) == 8u, "invalid data type size");
  static_assert(sizeof(IndexType) == 4u, "invalid index type size");

  static constexpr int dim = 2;

  typedef Shape::Hypercube<2> ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> MeshType;
  typedef Trafo::Standard::Mapping<MeshType> TrafoType;
  typedef Space::Lagrange1::Element<TrafoType> SpaceType;

  typedef Control::Domain::SimpleDomainLevel<MeshType, TrafoType, SpaceType> DomainLevelType;

  typedef Control::ScalarUnitFilterSystemLevel<DataType, IndexType> SystemLevelType;

  struct Counts
  {
    static constexpr std::size_t num_ranks = 0u;
    static constexpr std::size_t num_elems = 1u;
    static constexpr std::size_t num_dofs_l = 2u; // local dofs
    static constexpr std::size_t num_dofs_g = 3u; // global dofs [ != sum(num_dofs_l) ]
    static constexpr std::size_t num_nze = 4u;
    static constexpr std::size_t bytes_domain = 5u;
    static constexpr std::size_t bytes_system = 6u;
    static constexpr std::size_t bytes_solver = 7u;
    static constexpr std::size_t elems_mirror = 8u;
    static constexpr std::size_t count = 9u;
  };

  struct SystemTimes
  {
    static constexpr std::size_t asm_total = 0u;
    static constexpr std::size_t asm_gate = 1u;
    static constexpr std::size_t asm_muxer = 2u;
    static constexpr std::size_t asm_transfer = 3u;
    static constexpr std::size_t asm_matrix = 4u;
    static constexpr std::size_t asm_rhs = 5u;
    static constexpr std::size_t count = 6u;
  };

  struct SolverTimes
  {
    static constexpr std::size_t gmg_total = 0u;
    static constexpr std::size_t gmg_defect = 1u;
    static constexpr std::size_t gmg_smooth = 2u;
    static constexpr std::size_t gmg_transfer = 3u;
    static constexpr std::size_t gmg_coarse = 4u;
    static constexpr std::size_t gmg_apply = 5u;
    static constexpr std::size_t count = 6u;
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

  struct SystemStats
  {
    // per [level][index]
    std::vector<std::array<unsigned long long, Counts::count>> counts, counts_sum, counts_max;

    // per [level][index]
    std::vector<std::array<double, SystemTimes::count>> times;

    // (physical, virtual)
    std::array<unsigned long long, 2> mem_use, mem_use_sum, mem_use_max;

    explicit SystemStats(std::size_t virt_size) :
      counts(virt_size),
      counts_sum(virt_size),
      counts_max(virt_size),
      times(virt_size)
    {
      for(std::size_t i(0u); i < virt_size; ++i)
      {
        for(std::size_t j(0u); j < Counts::count; ++j)
          counts[i][j] = 0ull;
        for(std::size_t j(0u); j < SystemTimes::count; ++j)
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

      s += "\nAssembly Timings:\n";
      s += "                Gate        Muxer     Transfer       Matrix        RHS        Total\n";
      s += "Overall : " +
        stringify_fp_fix(sum(times, SystemTimes::asm_gate), 6, 10) + " / " +
        stringify_fp_fix(sum(times, SystemTimes::asm_muxer), 6, 10) + " / " +
        stringify_fp_fix(sum(times, SystemTimes::asm_transfer), 6, 10) + " / " +
        stringify_fp_fix(sum(times, SystemTimes::asm_matrix), 6, 10) + " / " +
        stringify_fp_fix(sum(times, SystemTimes::asm_rhs), 6, 10) + " / " +
        stringify_fp_fix(sum(times, SystemTimes::asm_total), 6, 10) + "\n";
      for(std::size_t i(0); i < times.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(times[i][SystemTimes::asm_gate], 6, 10) + " / " +
          stringify_fp_fix(times[i][SystemTimes::asm_muxer], 6, 10) + " / " +
          stringify_fp_fix(times[i][SystemTimes::asm_transfer], 6, 10) + " / " +
          stringify_fp_fix(times[i][SystemTimes::asm_matrix], 6, 10) + " / " +
          stringify_fp_fix(times[i][SystemTimes::asm_rhs], 6, 10) + " / " +
          stringify_fp_fix(times[i][SystemTimes::asm_total], 6, 10) + "\n";
      }

      s += "\nBasic Statistics:\n";
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

      s += "\nMirror #elems Statistics:\n";
      s += "Overall : " + stringify(sum(counts_sum, Counts::elems_mirror)).pad_front(15) +
        " [ " + stringify(sum(counts_max, Counts::elems_mirror)).pad_front(15) + " ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) +
          ": " + stringify(counts_sum[i][Counts::elems_mirror]).pad_front(15) +
          " [ " + stringify(counts_max[i][Counts::elems_mirror]).pad_front(15) + " ]\n";
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
      s += "Overall : " + stringify(sum(counts_sum, Counts::bytes_system)).pad_front(15) +
        " [ " + stringify(sum(counts_max, Counts::bytes_system)).pad_front(15) + " ]\n";
      for(std::size_t i(0); i < counts.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) +
          ": " + stringify(counts_sum[i][Counts::bytes_system]).pad_front(15) +
          " [ " + stringify(counts_max[i][Counts::bytes_system]).pad_front(15) + " ]\n";
      }

      return s;
    }

    String summary()
    {
      String s;
      const std::size_t total_elems = counts_sum[0][Counts::num_elems];
      const std::size_t total_gdofs = counts[0][Counts::num_dofs_g];
      const double total_asm_time = sum(times, SystemTimes::asm_total);
      const double asm_mdofs = (total_asm_time > 1E-5 ? 1E-6 * double(total_gdofs) / total_asm_time : 0.0);

      unsigned long long domain_bytes_sum = sum(counts_sum, Counts::bytes_domain);
      unsigned long long system_bytes_sum = sum(counts_sum, Counts::bytes_system);
      unsigned long long solver_bytes_sum = sum(counts_sum, Counts::bytes_solver);

      unsigned long long domain_bytes_max = sum(counts_max, Counts::bytes_domain);
      unsigned long long system_bytes_max = sum(counts_max, Counts::bytes_system);
      unsigned long long solver_bytes_max = sum(counts_max, Counts::bytes_solver);

      unsigned long long total_bytes_sum = domain_bytes_sum + system_bytes_sum + solver_bytes_sum + 3ull*total_gdofs*sizeof(DataType);
      unsigned long long total_bytes_max = domain_bytes_max + system_bytes_max + solver_bytes_max + 3ull*counts_max[0][Counts::num_dofs_l]*sizeof(DataType);

      s += "Overall Benchmark Summary:";
      s += "\nTotal Number of Elements........: " + stringify(total_elems).pad_front(15);
      s += "\nTotal Number of Global DOFs.....: " + stringify(total_gdofs).pad_front(15);
      s += "\nMaximum Domain Memory Usage.....: " + stringify_fp_fix(double(domain_bytes_max)/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += " [" + stringify_fp_fix(100.0*double(domain_bytes_max) / double(total_bytes_max), 2, 6) + "%]";
      s += "\nMaximum System Memory Usage.....: " + stringify_fp_fix(double(system_bytes_max)/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += " [" + stringify_fp_fix(100.0*double(system_bytes_max) / double(total_bytes_max), 2, 6) + "%]";
      s += "\nMaximum Solver Memory Usage.....: " + stringify_fp_fix(double(solver_bytes_max)/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += " [" + stringify_fp_fix(100.0*double(solver_bytes_max) / double(total_bytes_max), 2, 6) + "%]";
      s += "\nMaximum Combined Memory Usage...: " + stringify_fp_fix(double(total_bytes_max)/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += "\nMax Memory Usage Reported By OS.: " + stringify_fp_fix(double(mem_use_max[0])/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += " [" + stringify_fp_fix(100.0*(double(mem_use_max[0])/double(total_bytes_max) - 1.0), 2, 6, true) + "%]";
      s += "\nTotal Domain Memory Usage.......: " + stringify_fp_fix(double(domain_bytes_sum)/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += " [" + stringify_fp_fix(100.0*double(domain_bytes_sum) / double(total_bytes_sum), 2, 6) + "%]";
      s += "\nTotal System Memory Usage.......: " + stringify_fp_fix(double(system_bytes_sum)/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += " [" + stringify_fp_fix(100.0*double(system_bytes_sum) / double(total_bytes_sum), 2, 6) + "%]";
      s += "\nTotal Solver Memory Usage.......: " + stringify_fp_fix(double(solver_bytes_sum)/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += " [" + stringify_fp_fix(100.0*double(solver_bytes_sum) / double(total_bytes_sum), 2, 6) + "%]";
      s += "\nTotal Combined Memory Usage.....: " + stringify_fp_fix(double(total_bytes_sum)/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += "\nMemory Usage Reported By OS.....: " + stringify_fp_fix(double(mem_use_sum[0])/(1024.0*1024.0*1024.0), 6, 15) + " GiB";
      s += " [" + stringify_fp_fix(100.0*(double(mem_use_sum[0])/double(total_bytes_sum) - 1.0), 2, 6, true) + "%]";
      s += "\nTotal Assembly Runtime..........: " + stringify_fp_fix(total_asm_time, 6, 15) + " seconds";
      s += "\nAssembly Performance............: " + stringify_fp_fix(asm_mdofs, 6, 15) + " MDOF/s";

      return s;
    }
  };

  struct SolverStats
  {
    // per [level][index]
    std::vector<std::array<double, SolverTimes::count>> times;

    explicit SolverStats(std::size_t virt_size) :
      times(virt_size)
    {
      for(std::size_t i(0u); i < virt_size; ++i)
      {
        for(std::size_t j(0u); j < SolverTimes::count; ++j)
          times[i][j] = 0.0;
      }
    }

    String format() const
    {
      String s;

      s += "                Defect /     Smoother /     Transfer /       Coarse /        Total /        Apply\n";
      s += "Overall : " +
        stringify_fp_fix(sum(times, SolverTimes::gmg_defect), 6, 12) + " / " +
        stringify_fp_fix(sum(times, SolverTimes::gmg_smooth), 6, 12) + " / " +
        stringify_fp_fix(sum(times, SolverTimes::gmg_transfer), 6, 12) + " / " +
        stringify_fp_fix(sum(times, SolverTimes::gmg_coarse), 6, 12) + " / " +
        stringify_fp_fix(sum(times, SolverTimes::gmg_total), 6, 12) + " / " +
        stringify_fp_fix(sum(times, SolverTimes::gmg_apply), 6, 12) + "\n";
      for(std::size_t i(0); i < times.size(); ++i)
      {
        s += "Level " + stringify(i).pad_front(2) + ": " +
          stringify_fp_fix(times[i][SolverTimes::gmg_defect], 6, 12) + " / " +
          stringify_fp_fix(times[i][SolverTimes::gmg_smooth], 6, 12) + " / " +
          stringify_fp_fix(times[i][SolverTimes::gmg_transfer], 6, 12) + " / " +
          stringify_fp_fix(times[i][SolverTimes::gmg_coarse], 6, 12) + " / " +
          stringify_fp_fix(times[i][SolverTimes::gmg_total], 6, 12) + "\n";
      }

      return s;
    }

    String summary(String prefix, const std::size_t total_gdofs)
    {
      const double total_gmg_time = sum(times, SolverTimes::gmg_total);
      const double gmg_mdofs = (total_gmg_time > 1E-5 ? 1E-6 * double(total_gdofs) / total_gmg_time : 0.0);

      String s;
      s += (prefix + " Runtime").pad_back(32, '.') + ": " + stringify_fp_fix(total_gmg_time, 6, 15) + " seconds\n";
      s += (prefix + " Performance").pad_back(32, '.') + ": " + stringify_fp_fix(gmg_mdofs, 6, 15) + " MDOF/s";
      return s;
    }
  };

  template<typename Matrix_, typename Filter_>
  class DampedSmoother :
    public Solver::SolverBase<typename Matrix_::VectorTypeR>
  {
  public:
    typedef Matrix_ MatrixType;
    typedef Filter_ FilterType;
    typedef typename MatrixType::VectorTypeR VectorType;
    typedef typename MatrixType::DataType DataType;
    typedef Solver::SolverBase<VectorType> BaseClass;

  protected:
    /// the matrix for the solver
    const MatrixType& _system_matrix;
    /// the filter for the solver
    const FilterType& _system_filter;
    /// number of steps
    Index _steps;
    /// damping parameter
    DataType _omega;
    /// defect vector
    VectorType _vec_tmp;

  public:
    explicit DampedSmoother(const MatrixType& matrix, const FilterType& filter, Index steps, DataType omega) :
      BaseClass(),
      _system_matrix(matrix),
      _system_filter(filter),
      _steps(steps),
      _omega(omega)
    {
    }

    virtual String name() const override
    {
      return "DampedSmoother";
    }

    virtual void init_symbolic() override
    {
      BaseClass::init_symbolic();
      _vec_tmp = this->_system_matrix.create_vector_r();
    }

    virtual void done_symbolic() override
    {
      this->_vec_tmp.clear();
      BaseClass::done_symbolic();
    }

    virtual Solver::Status apply(VectorType& vec_cor, const VectorType& vec_def) override
    {
      vec_cor.scale(vec_def, _omega);
      for(Index k = 1; k < _steps; ++k)
      {
        this->_system_matrix.apply(this->_vec_tmp, vec_cor, vec_def, -DataType(1));
        this->_system_filter.filter_cor(this->_vec_tmp);
        vec_cor.axpy(this->_vec_tmp, this->_omega);
      }

      return Solver::Status::success;
    }
  }; // class DampedSmoother<...>

  void asm_transfer(SystemLevelType& sys_lvl_f, const Control::Domain::VirtualLevel<DomainLevelType>& virt_lvl_c)
  {
    const DomainLevelType& dom_lvl_c = virt_lvl_c.is_child() ? virt_lvl_c.level_c() : *virt_lvl_c;
    const MeshType& mesh_c = dom_lvl_c.get_mesh();
    const Index nvc = mesh_c.get_num_entities(0);
    const Index nve = mesh_c.get_num_entities(1);
    const Index nvf = mesh_c.get_num_entities(2);

    const Index nrows = nvc + nve + nvf;
    const Index ncols = nvc;
    const Index nnzes = nvc + 2*nve + 4*nvf;

    const auto& v_at_e = mesh_c.template get_index_set<1,0>();
    const auto& v_at_f = mesh_c.template get_index_set<2,0>();

    auto& mat_prol = sys_lvl_f.transfer_sys.get_mat_prol();
    mat_prol = LAFEM::SparseMatrixCSR<DataType, IndexType>(nrows, ncols, nnzes);

    Memory::TypedView<IndexType> row_ptr(mat_prol.row_ptr_view_w());
    Memory::TypedView<IndexType> col_idx(mat_prol.col_idx_view_w());
    Memory::TypedView<DataType> val(mat_prol.val_view_w());

    row_ptr[0] = 0u;
    Index k = 0;
    for(Index i = 0; i < nvc; ++i)
    {
      col_idx[k] = IndexType(i);
      val[k] = DataType(1);
      row_ptr[i+1] = IndexType(++k);
    }
    for(Index i = 0; i < nve; ++i)
    {
      for(int j = 0; j  < 2; ++j, ++k)
      {
        col_idx[k] = IndexType(v_at_e[i][j]);
        val[k] = DataType(0.5);
      }
      row_ptr[nvc+i+1] = IndexType(k);
    }
    for(Index i = 0; i < nvf; ++i)
    {
      for(int j = 0; j < 4; ++j, ++k)
      {
        col_idx[k] = IndexType(v_at_f[i][j]);
        val[k] = DataType(0.25);
      }
      row_ptr[nvc+nve+i+1] = IndexType(k);
    }
    val.release();
    col_idx.release();
    row_ptr.release();

    mat_prol.scale_rows(mat_prol, sys_lvl_f.gate_sys.get_freqs());
    sys_lvl_f.transfer_sys.get_mat_rest() = mat_prol.transpose();
    sys_lvl_f.transfer_sys.compile();
  }

  void asm_mat_turbo(LAFEM::SparseMatrixCSR<DataType, IndexType>& matrix, const DomainLevelType& dom_lvl)
  {
    const MeshType& mesh = dom_lvl.get_mesh();

    const auto& v_at_e = mesh.get_index_set<2,0>();

    const Index num_v = mesh.get_num_vertices();
    const Index num_e = mesh.get_num_elements();

    static constexpr Index mv = 10;
    std::vector<IndexType> count(num_v, IndexType(0)), idx(mv*num_v);
    std::vector<int> val(mv*num_v, 0);

    static const int vaa[4][4] =
    {
      { 4, -1, -1, -2},
      {-1,  4, -2, -1},
      {-1, -2,  4, -1},
      {-2, -1, -1,  4}
    };

#if 0
    for(Index i = 0; i < num_e; ++i)
    {
      for(int j = 0; j < 4; ++j)
      {
        const Index vj = v_at_e[i][j];
        for(int k = 0; k < 4; ++k)
        {
          idx[mv*vj + count[vj]++] = v_at_e[i][k];
        }
      }
    }

    Index nnzes = Index(0);
    for(Index i = 0; i < num_v; ++i)
    {
      IndexType* vi = &idx[mv*i];
      IndexType* vj = vi + count[i];
#if 1
      std::sort(vi, vj);
      Index j = 0;
      for(Index k = 1; k < count[i]; ++k)
      {
        if(vi[k] != vi[j])
          vi[++j] = vi[k];
      }
      count[i] = j+1;
#else
      std::set<IndexType> si(vi, vj);
      count[i] = IndexType(si.size());
      auto it = si.begin();
      for(IndexType j(0); j < count[i]; ++j, ++it, ++vi)
        *vi = *it;
#endif
      nnzes += Index(count[i]);
    }
#else
    for(Index i = 0; i < num_e; ++i)
    {
      for(int j = 0; j < 4; ++j)
      {
        const Index vj = IndexType(v_at_e[i][j]);
        for(int k = 0; k < 4; ++k)
        {
          const IndexType vk = IndexType(v_at_e[i][k]);
          IndexType* vb = &idx[mv*vj];
          int* vv = &val[mv*vj];

          bool found = false;
          for(IndexType l = 0; l < count[vj]; ++l)
          {
            if(vb[l] == vk)
            {
              vv[l] += vaa[j][k];
              found = true;
              break;
            }
            else if(vk < vb[l])
            {
              for(IndexType m = count[vj]; m > l; --m)
              {
                vb[m] = vb[m-1];
                vv[m] = vv[m-1];
              }
              vb[l] = vk;
              vv[l] = vaa[j][k];
              ++count[vj];
              found = true;
              break;
            }
          }
          if(!found)
          {
            vb[count[vj]] = vk;
            vv[count[vj]] = vaa[j][k];
            ++count[vj];
          }
        }
      }
    }

    Index nnzes = Index(0);
    for(Index i = 0; i < num_v; ++i)
    {
      nnzes += Index(count[i]);
    }
#endif

    matrix = LAFEM::SparseMatrixCSR<DataType, IndexType>(num_v, num_v, nnzes);
    Memory::TypedView<IndexType> row_ptr(matrix.row_ptr_view_w());
    Memory::TypedView<IndexType> col_idx(matrix.col_idx_view_w());
    Memory::TypedView<DataType> mat_val(matrix.val_view_w());
    row_ptr[0] = IndexType(0);
    for(Index i = 0, k = 0; i < num_v; ++i)
    {
      row_ptr[i+1] = row_ptr[i] + count[i];
      for(Index j = 0; j < count[i]; ++j, ++k)
      {
        col_idx[k] = idx[mv*i + j];
        mat_val[k] = DataType(val[mv*i + j]) / DataType(6);
      }
    }
    //TimeStamp st4;

    //String s = ">>>";
    //s += st2.elapsed_string(st1).pad_front(10);
    //s += st3.elapsed_string(st2).pad_front(10);
    //s += st4.elapsed_string(st3).pad_front(10);
    //s += "\n";
    //std::cout << s;
  }

  void main(int argc, char** argv)
  {
    Dist::Comm comm(Dist::Comm::world());

    TimeStamp stamp_start;

    Backend::set_preferred_backend(PreferredBackend::generic);

    // enable solver expressions
    Statistics::enable_solver_expressions = true;

    SimpleArgParser args(argc, argv);
    Control::Domain::add_supported_pdc_args(args);

    args.support("slices", "<k>\n"
      "Specifies the number of element slices of the level-0 base mesh in each dimension, so that the base mesh\n"
      "corresponds to a <k>-by-<k> rectilinear mesh; default: 1.");
    args.support("level", "<fine [coarse...]>\n"
      "Specifies the refinement level of the fine mesh on which the benchmark system is solved as well as (optionally)\n"
      "the coarse level for the multigrid hierarchy, potentially including the intermediate levels for the recursive\n"
      "partitioning layers, as required by the Control::PartiDomainControl class. If only the fine level is given,\n"
      "then the multigrid hierarchy is automatically extended down to level 0 on a single MPI process.");
    args.support("iters");
    args.support("ext-stats", "\nPrint extended solver and MPI statistics.");
    args.support("std-assembly");
    args.support("solve", "<levels>");

    // no arguments given?
    if(args.num_args() <= 1)
    {
      comm.print(args.get_supported_help());
      return;
    }

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if (!unsupported.empty())
    {
      // print all unsupported options to cerr
      for (auto it = unsupported.begin(); it != unsupported.end(); ++it)
        comm.print(std::cerr, "ERROR: unknown option '--" + (*it).second + "'");

      comm.print(std::cerr, "Supported Options are:");
      comm.print(std::cerr, args.get_supported_help());

      // abort
      FEAT::Runtime::abort();
    }

    const bool std_assembly = (args.check("std-assembly") >= 0);

    comm.print(String(100u, '#'));
    comm.print("Dist Comm World Size: " + stringify(comm.size()));
    comm.print_flush();

    int multigrid_iters = 3;
    args.parse("iters", multigrid_iters);

    int solve_levels = 3;
    args.parse("solve", solve_levels);

    // create domain control
    Control::Domain::PartiDomainControl<DomainLevelType> domain(comm, true);
    if(!domain.parse_args(args))
      Runtime::abort();

    Index base_mesh_num_slices = 1u;

    if(args.check("level") > 1)
    {
      // multiple levels given; pass them directly to the domain control
      domain.set_desired_levels(args.query("level")->second);
    }
    else if(args.check("level") == 1)
    {
      // just fine mesh level given; go down to level 0
      int level = 0;
      args.parse("level", level);
      if(comm.size() > 1)
        domain.set_desired_levels(level, 0, 0);
      else
        domain.set_desired_levels(level, 0);
    }
    else
    {
      comm.print(std::cerr, "ERROR: you must specify a fine mesh level by --level <level>");
      FEAT::Runtime::abort();
    }

    // parse slices; this may override the default setting defined by the --memory option handling above
    args.parse("slices", base_mesh_num_slices);

    // create domain from a rectilinear base mesh
    domain.create_rectilinear(base_mesh_num_slices, base_mesh_num_slices);

    // get the number of slices on the fine mesh
    Index fine_mesh_num_slices = base_mesh_num_slices << (domain.get_chosen_levels().front().first);

    // print partitioning info
    //          12345678901234567890
    comm.print("Base Mesh Dimensions: " + stringify(base_mesh_num_slices) + "^2 = " +
      stringify(base_mesh_num_slices*base_mesh_num_slices) + " elements");
    comm.print("Fine Mesh Dimensions: " + stringify(fine_mesh_num_slices) + "^2 = " +
      stringify(fine_mesh_num_slices*fine_mesh_num_slices) + " elements");
    comm.print("Desired Levels......: " + domain.format_desired_levels());
    comm.print("Chosen Levels.......: " + domain.format_chosen_levels());
    comm.print("Partitioner Info....: " + domain.get_chosen_parti_info());
    comm.print("Multigrid Iterations: " + stringify(multigrid_iters));
    comm.print("Assembly Type.......: " + String(std_assembly ? "standard" : "turbo"));
    comm.print_flush();

    // create benchmark statistics
    SystemStats sys_stats(domain.size_virtual());
    std::vector<SolverStats> vec_mg_v_stats, vec_mg_f_stats;

    {
      std::size_t n = domain.size_virtual() > domain.size_physical() ? domain.size_physical()+1 : domain.size_physical();
      for(std::size_t i(0); i < n; ++i)
      {
        const auto& vdl = domain.at(i);
        if(vdl.is_parent())
          sys_stats.counts[i][Counts::bytes_domain] = vdl.level_c().bytes() + vdl.level_p().bytes();
        else if(vdl.is_child())
          sys_stats.counts[i][Counts::bytes_domain] = vdl.level_c().bytes();
        else
          sys_stats.counts[i][Counts::bytes_domain] = vdl.level().bytes();
      }
    }

    const Index num_levels = domain.size_physical();

    // create system levels
    std::deque<std::shared_ptr<SystemLevelType>> system;

    for (Index i(0); i < num_levels; ++i)
    {
      system.push_back(std::make_shared<SystemLevelType>());
    }

    const String cubature("gauss-legendre:" + stringify(SpaceType::local_degree+1));

    /* ***************************************************************************************** */

    TimeStamp stamp_ass;

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      domain.at(i)->domain_asm.compile_all_elements();
      system.at(i)->assemble_gate(domain.at(i));
      sys_stats.times[i][SystemTimes::asm_gate] += ts.elapsed_now();
      sys_stats.times[i][SystemTimes::asm_total] += ts.elapsed_now();
    }

    for (Index i(0); (i < domain.size_physical()) && ((i+1) < domain.size_virtual()); ++i)
    {
      TimeStamp ts;
      system.at(i)->assemble_coarse_muxer(domain.at(i+1));
      sys_stats.times[i][SystemTimes::asm_muxer] += ts.elapsed_now();
      TimeStamp ts2;
      if(std_assembly)
      {
        if((i+1) < domain.size_physical())
          system.at(i)->assemble_transfer(*system.at(i+1), domain.at(i), domain.at(i+1), cubature);
        else
          system.at(i)->assemble_transfer(domain.at(i), domain.at(i+1), cubature);
        system.at(i)->transfer_sys.get_mat_prol().shrink(1E-3);
        system.at(i)->transfer_sys.get_mat_rest().shrink(1E-3);
      }
      else
        asm_transfer(*system.at(i), domain.at(i+1));

      sys_stats.times[i][SystemTimes::asm_transfer] += ts2.elapsed_now();
      sys_stats.times[i][SystemTimes::asm_total] += ts.elapsed_now();
    }

    /* ***************************************************************************************** */
    // Do symbolic and numeric assembly independently
    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      if(std_assembly)
      {
        system.at(i)->symbolic_assembly_std1(domain.at(i)->space);
        system.at(i)->assemble_laplace_matrix(domain.at(i)->domain_asm, domain.at(i)->space, cubature);
      }
      else
        asm_mat_turbo(system.at(i)->matrix_sys.local(), *domain.at(i));

      if(0)
      {
        const auto& mat = system.at(i)->matrix_sys.local();
        Memory::TypedView<IndexType> row_ptr(mat.row_ptr_view_r());
        Memory::TypedView<IndexType> col_idx(mat.col_idx_view_r());

        auto crc_rp = Hash::crc32((mat.num_rows()+1)*sizeof(IndexType), row_ptr.raw_r());
        auto crc_ci = Hash::crc32((mat.num_nzes())*sizeof(IndexType), col_idx.raw_r());
        comm.allprint(stringify(crc_rp) + " / " + stringify(crc_ci));
      }

      double tt = ts.elapsed_now();
      sys_stats.times[i][SystemTimes::asm_total] += tt;
      sys_stats.times[i][SystemTimes::asm_matrix] += tt;
      //system.at(i)->matrix_sys.local().print(std::cout, false);
    }

    /* ***************************************************************************************** */

    for (Index i(0); i < num_levels; ++i)
    {
      TimeStamp ts;
      system.at(i)->assemble_homogeneous_unit_filter(*domain.at(i), domain.at(i)->space);

      // apply filter
      system.at(i)->filter_sys.local().filter_mat(system.at(i)->matrix_sys.local());
      sys_stats.times[i][SystemTimes::asm_total] += ts.elapsed_now();

      sys_stats.counts[i][Counts::num_ranks] = Index(domain.at(i).layer().comm().size());
      sys_stats.counts[i][Counts::num_elems] = domain.at(i)->get_mesh().get_num_elements();
      sys_stats.counts[i][Counts::num_dofs_g] = system.at(i)->matrix_sys.num_rows();
      sys_stats.counts[i][Counts::num_dofs_l] = system.at(i)->matrix_sys.local().num_rows();
      sys_stats.counts[i][Counts::num_nze] = system.at(i)->matrix_sys.local().num_nzes();
      sys_stats.counts[i][Counts::bytes_system] = system.at(i)->bytes();
      sys_stats.counts[i][Counts::elems_mirror] = 0;
      auto& sys_mirrors =  system.at(i)->gate_sys._mirrors;
      for (auto& mirror : sys_mirrors)
      {
        sys_stats.counts[i][Counts::elems_mirror] += mirror.num_indices();
      }
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    auto multigrid_hierarchy = std::make_shared<
      Solver::MultiGridHierarchy<
      typename SystemLevelType::GlobalSystemMatrix,
      typename SystemLevelType::GlobalSystemFilter,
      typename SystemLevelType::GlobalSystemTransfer
      > >(domain.size_virtual());

    // push all levels except the coarse most one
    for (Index i(0); i < num_levels; ++i)
    {
      const SystemLevelType& lvl = *system.at(i);
      const std::size_t lvl_dofs = lvl.matrix_sys.local().num_rows();

      if((i+1) < domain.size_virtual())
      {
        auto smoother = std::make_shared<DampedSmoother<typename SystemLevelType::GlobalSystemMatrix,
          typename SystemLevelType::GlobalSystemFilter>>(lvl.matrix_sys, lvl.filter_sys, 4, 0.25);

        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, lvl.transfer_sys, smoother, smoother, smoother);

        // Richardson: 1 Vectors
        // Multigrid: 5 Vectors
        sys_stats.counts.at(i)[Counts::bytes_solver] = 6ull * lvl_dofs * sizeof(DataType);
      }
      else
      {
        auto cgsolver = Solver::new_jacobi_precond(lvl.matrix_sys, lvl.filter_sys, 1.0);
        multigrid_hierarchy->push_level(lvl.matrix_sys, lvl.filter_sys, cgsolver);

        // Jacobi: 1 Vector
        sys_stats.counts.at(i)[Counts::bytes_solver] = lvl_dofs * sizeof(DataType);
      }
    }

    multigrid_hierarchy->init();

    /* ***************************************************************************************** */

    // get our assembled vector type
    typedef typename SystemLevelType::GlobalSystemVector GlobalSystemVector;

    double solver_toe_v = 0.0;
    double solver_toe_f = 0.0;

    for(int ilv = 0; ilv < solve_levels; ++ilv)
    {
      vec_mg_v_stats.push_back(SolverStats(domain.size_virtual() - ilv));
      vec_mg_f_stats.push_back(SolverStats(domain.size_virtual() - ilv));
      SolverStats& mg_v_stats = vec_mg_v_stats.back();
      SolverStats& mg_f_stats = vec_mg_f_stats.back();

      // fetch our finest levels
      DomainLevelType& the_domain_level = *domain.at(ilv);
      SystemLevelType& the_system_level = *system.at(ilv);

      const int level_index = the_domain_level.get_level_index();

      comm.print(String(100, '-'));
      comm.print("Solving on Level " + stringify(level_index));

      TimeStamp ts_rhs;

      // create new vector
      GlobalSystemVector vec_sol = the_system_level.matrix_sys.create_vector_r();
      GlobalSystemVector vec_rhs = the_system_level.matrix_sys.create_vector_r();

      // format solution vector
      vec_sol.format();
      vec_rhs.format(DataType(1) / DataType(1ll << level_index));

      // and filter it
      the_system_level.filter_sys.filter_sol(vec_sol);
      the_system_level.filter_sys.filter_rhs(vec_rhs);

      sys_stats.times[ilv][SystemTimes::asm_rhs] = ts_rhs.elapsed_now();

      // create multigrid solver
      auto multigrid = Solver::new_multigrid(multigrid_hierarchy, Solver::MultiGridCycle::V, ilv);

      // initialize
      multigrid->init();

      Statistics::reset();

      comm.print("\nWarming up...");

      // warmup solve
      multigrid->apply(vec_sol, vec_rhs);

      comm.print("\nBenchmarking V-Cycle...");
      multigrid_hierarchy->reset_timings();

      multigrid->set_cycle(Solver::MultiGridCycle::V);

      TimeStamp st_v1;

      // benchmarking solve
      for(int k = 0; k < multigrid_iters; ++k)
      {
        multigrid->apply(vec_sol, vec_rhs);
      }

      TimeStamp st_v2;

      double solve_time_v = st_v2.elapsed(st_v1);
      solver_toe_v += solve_time_v;
      comm.print("Solve Apply Time: " + stringify_fp_fix(solve_time_v, 6));

      // set multigrid timings
      mg_v_stats.times[0][SolverTimes::gmg_apply] = solve_time_v;
      for(int i(0); i+ilv < int(multigrid_hierarchy->size_physical()); ++i)
      {
        mg_v_stats.times.at(i)[SolverTimes::gmg_defect] = multigrid_hierarchy->get_time_defect(ilv+int(i));
        mg_v_stats.times.at(i)[SolverTimes::gmg_smooth] = multigrid_hierarchy->get_time_smooth(ilv+int(i));
        mg_v_stats.times.at(i)[SolverTimes::gmg_transfer] = multigrid_hierarchy->get_time_transfer(ilv+int(i));
        mg_v_stats.times.at(i)[SolverTimes::gmg_coarse] = multigrid_hierarchy->get_time_coarse(ilv+int(i));
        mg_v_stats.times.at(i)[SolverTimes::gmg_total] = 
          mg_v_stats.times.at(i)[SolverTimes::gmg_defect] +
          mg_v_stats.times.at(i)[SolverTimes::gmg_smooth] +
          mg_v_stats.times.at(i)[SolverTimes::gmg_transfer] +
          mg_v_stats.times.at(i)[SolverTimes::gmg_coarse];
      }

      comm.print("\nBenchmarking F-Cycle...");
      multigrid_hierarchy->reset_timings();

      multigrid->set_cycle(Solver::MultiGridCycle::F);

      TimeStamp st_f1;

      // benchmarking solve
      for(int k = 0; k < multigrid_iters; ++k)
      {
        multigrid->apply(vec_sol, vec_rhs);
      }

      TimeStamp st_f2;
      double solve_time_f = st_f2.elapsed(st_f1);
      solver_toe_f += solve_time_f;
      comm.print("Solve Apply Time: " + stringify_fp_fix(solve_time_f, 6));

      // set multigrid timings
      mg_f_stats.times[0][SolverTimes::gmg_apply] = solve_time_f;
      for(int i(0); i+ilv < int(multigrid_hierarchy->size_physical()); ++i)
      {
        mg_f_stats.times.at(i)[SolverTimes::gmg_defect] = multigrid_hierarchy->get_time_defect(ilv+int(i));
        mg_f_stats.times.at(i)[SolverTimes::gmg_smooth] = multigrid_hierarchy->get_time_smooth(ilv+int(i));
        mg_f_stats.times.at(i)[SolverTimes::gmg_transfer] = multigrid_hierarchy->get_time_transfer(ilv+int(i));
        mg_f_stats.times.at(i)[SolverTimes::gmg_coarse] = multigrid_hierarchy->get_time_coarse(ilv+int(i));
        mg_f_stats.times.at(i)[SolverTimes::gmg_total] = 
          mg_f_stats.times.at(i)[SolverTimes::gmg_defect] +
          mg_f_stats.times.at(i)[SolverTimes::gmg_smooth] +
          mg_f_stats.times.at(i)[SolverTimes::gmg_transfer] +
          mg_f_stats.times.at(i)[SolverTimes::gmg_coarse];
      }

      // release solver
      multigrid->done();

      comm.print_flush();

      // next solver level
    }

    comm.print(String(100u, '='));

    multigrid_hierarchy->done();

    // get memory info
    {
      MemoryUsage meminfo;
      sys_stats.mem_use[0] = meminfo.get_peak_physical();
      sys_stats.mem_use[1] = meminfo.get_peak_virtual();
    }

    sys_stats.sync(comm);

    comm.print(sys_stats.format());
    for(int i = 0; i < solve_levels; ++i)
    {
      comm.print("Multigrid V-Cycle Statistics for Level " + stringify(domain.at(i)->get_level_index()));
      comm.print(vec_mg_v_stats.at(i).format());
    }
    for(int i = 0; i < solve_levels; ++i)
    {
      comm.print("Multigrid F-Cycle Statistics for Level " + stringify(domain.at(i)->get_level_index()));
      comm.print(vec_mg_f_stats.at(i).format());
    }

    comm.print("Multigrid Solve Time Analysis");
    comm.print("            Apply MG-V /  Factor MG-V /   Apply MG-F /  Factor MG-F");
    for(int i = 0; i < solve_levels; ++i)
    {
      String s = "Level" + stringify(domain.at(i)->get_level_index()).pad_front(3) + ": ";
      s += stringify_fp_fix(vec_mg_v_stats.at(i).times[0][SolverTimes::gmg_apply], 6, 12) + " / ";
      if(i > 0)
        s += stringify_fp_fix(vec_mg_v_stats.at(i-1).times[0][SolverTimes::gmg_apply] / vec_mg_v_stats.at(i).times[0][SolverTimes::gmg_apply], 3, 12) + " / ";
      else
        s += "         --- / ";
      s += stringify_fp_fix(vec_mg_f_stats.at(i).times[0][SolverTimes::gmg_apply], 6, 12) + " / ";
      if(i > 0)
        s += stringify_fp_fix(vec_mg_f_stats.at(i-1).times[0][SolverTimes::gmg_apply] / vec_mg_f_stats.at(i).times[0][SolverTimes::gmg_apply], 3, 12);
      else
        s += "         ---";
      comm.print(s);
    }
    comm.print("");
    comm.print_flush();

    comm.print(sys_stats.summary());
    comm.print(vec_mg_v_stats.front().summary("Multigrid V-Cycle", sys_stats.counts[0][Counts::num_dofs_g]));
    comm.print(vec_mg_f_stats.front().summary("Multigrid F-Cycle", sys_stats.counts[0][Counts::num_dofs_g]));
    comm.print_flush();

    if(args.check("ext-stats") >= 0)
    {
      FEAT::Control::Statistics::report(solver_toe_v+solver_toe_f, 0, MeshType::ShapeType::dimension, system, domain);
      comm.print(FEAT::Statistics::get_formatted_solver_internals());
      comm.print(FEAT::Statistics::get_formatted_solver_tree().trim());
    }

    TimeStamp stamp_end;
    comm.print("\nTotal Runtime: "  + stamp_end.elapsed_string(stamp_start));
  }
} // namespace ClusterMultigridBench

int main(int argc, char** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  ClusterMultigridBench::main(argc, argv);
  return 0;
}
