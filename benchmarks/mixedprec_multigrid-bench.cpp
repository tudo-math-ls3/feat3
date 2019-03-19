//
// Mixed-Precision Multi-Grid Benchmark
// ------------------------------------
// This code implements a hard-wired Richardson-Multigrid-Jacobi solver
// that is used to solve the 1D Poisson equation
//
//               -u_xx(x) = f(x)    for x in (0,1)
//            u(0) = u(1) = 0
//
// with the analytical solution
//
//                   u(x) = sin(pi*x)
//
// This benchmarks allows to choose different floating point types for
// 1) the assembly of the matrices and vectors as well as the H0-error computation
// 2) the outer Richardson solver iteration
// 3) the inner Multigrid preconditioner iteration
//
// This code does not use any of FEAT's sophisticated containers, but implements
// all linear algebra operations and assemblies by hand to allow for efficient
// OpenMP parallelisation, which is crucial when working with emulated quad precision.
//
// \author Peter Zajac
//
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/math.hpp>
#include <kernel/util/memory_usage.hpp>
#include <kernel/util/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/string.hpp>

#include <vector>
#include <fstream>
#include <cstring> // for memcpy

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace FEAT;

namespace MixedPrecMultiGridBench
{
  template<typename F_> struct Typo {static const char* name(); };

  // quad precision (if available)
#ifdef FEAT_HAVE_QUADMATH
  typedef __float128       f_qp;
  template<> struct Typo<f_qp> {static const char* name() {return "qp";} };
#endif

  // double precision
  typedef double           f_dp;
  template<> struct Typo<f_dp> {static const char* name() {return "dp";} };

  // single precision
  typedef float            f_sp;
  template<> struct Typo<f_sp> {static const char* name() {return "sp";} };

  // half precision (requires FloatX third-party library)
#ifdef FEAT_HAVE_FLOATX
  // Note: 'float' as backend is bugged
  typedef flx::floatx<5, 10, double> f_hp;
  template<> struct Typo<f_hp> {static const char* name() {return "hp";} };
#endif

  // use long long int for sizes
  // Note: we cannot use 'std::size_t' because OpenMP doesn't like that
  typedef long long int llint;

  // minimum vector length for OpenMP parallelisation
#define N_MIN_OMP 1000

  // assemble RHS vector by linearform evaluation
  template<typename F_>
  void assemble_rhs_sin(std::vector<F_>& vrhs)
  {
    // compute Pi and h
    const llint n = llint(vrhs.size());
    const F_ pi = Math::pi<F_>();
    const F_ h = F_(1) / F_(n+1);
    F_* vx = vrhs.data();

    // Our force functional is given by
    //
    //   f(x) := pi^2 * sin(pi*x)
    //
    // and our RHS vector entries are given by
    //
    //   b_i := \int_\Omega f(x) * \phi_i(x) dx,
    //
    // where \phi_i are the basis functions of our P1 space.
    // In our case, it can be shown that for \Omega=[0,1] we have
    //
    //   b_i = (2*n+2) * (1 - cos(pi/(n+1))) * sin(pi*(i+1)/(n+1))
    //

    // pre-compute factor
    const F_ fac = F_(2*n+2) * (F_(1) - Math::cos(pi*h));

    // initialise rhs entries
    #pragma omp parallel for if(n > N_MIN_OMP)
    for(llint i = 0; i < n; ++i)
    {
      vx[i] = fac * Math::sin(pi * F_(i+1) * h);
    }
  }

  // analyse H0 error of solution vector
  template<typename F_>
  F_ analyse_sol_sin(const std::vector<F_>& vsol)
  {
    // compute Pi and h
    const llint n = llint(vsol.size());
    const F_ pi = Math::pi<F_>();
    const F_ h = F_(1) / F_(n+1);

    // 3-point Gauss quadrature
    const F_ gp[3] =
    {
      (F_(0.5) - Math::sqrt(F_(3) / F_(20))),
      (F_(0.5) + Math::sqrt(F_(3) / F_(20))),
      (F_(0.5))
    };
    const F_ gw[3] =
    {
      h * F_(5) / F_(18),
      h * F_(5) / F_(18),
      h * F_(8) / F_(18)
    };

    F_ err_h0 = F_(0);
    const F_* vx = vsol.data();

    // loop over all intervals
    #pragma omp parallel for reduction(+:err_h0) if(n > N_MIN_OMP)
    for(llint i = 0; i <= n; ++i)
    {
      // left/right vertex coords
      const F_ xl = F_(i  ) / F_(n+1);
      const F_ xr = F_(i+1) / F_(n+1);
      // left/right unknowns
      const F_ vl = (i > 0 ? vx[i-1] : F_(0));
      const F_ vr = (i < n ? vx[i]   : F_(0));

      // integrate squared error over current sub-interval
      for(int j(0); j < 3; ++j)
      {
        const F_ gq = F_(1) - gp[j];
        err_h0 += gw[j] * Math::sqr(Math::sin(pi*(gq*xl + gp[j]*xr)) - (gq*vl + gp[j]*vr));
      }
    }

    return Math::sqrt(err_h0);
  }

  // format a vector
  template<typename F_>
  void format(std::vector<F_>& v, const F_ val = F_(0))
  {
    const llint n = llint(v.size());
    F_* x = v.data();

    #pragma omp parallel for if(n > N_MIN_OMP)
    for(llint i = 0; i < n; ++i)
      x[i] = val;
  }

  // copy/convert a vector
  template<typename Fx_, typename Fy_>
  void copy(std::vector<Fx_>& vx, const std::vector<Fy_>& vy)
  {
    XASSERT(vx.size() == vy.size());
    const llint n = llint(vx.size());
    Fx_* x = vx.data();
    const Fy_* y = vy.data();

    #pragma omp parallel for if(n > N_MIN_OMP)
    for(llint i = 0; i < n; ++i)
      x[i] = Fx_(y[i]);
  }

  // copy a vector
  template<typename F_>
  void copy(std::vector<F_>& vx, const std::vector<F_>& vy)
  {
    XASSERT(vx.size() == vy.size());
    std::memcpy(vx.data(), vy.data(), vx.size() * sizeof(F_));
  }

  // calculate vector norm
  template<typename F_>
  F_ norm2(const std::vector<F_>& v)
  {
    const llint n = llint(v.size());
    const F_* x = v.data();
    F_ r = F_(0);

    #pragma omp parallel for reduction(+:r) if(n > N_MIN_OMP)
    for(llint i = 0; i < n; ++i)
      r += x[i]*x[i];

    return Math::sqrt(r);
  }

  // perform axpy
  template<typename F_>
  void axpy(std::vector<F_>& vy, const std::vector<F_>& vx, const F_ alpha = F_(1))
  {
    XASSERT(vy.size() == vx.size());
    const llint n = llint(vy.size());

    F_* y = vy.data();
    const F_* x = vx.data();

    #pragma omp parallel for if(n > N_MIN_OMP)
    for(llint i = 0; i < n; ++i)
      y[i] += alpha * x[i];
  }

  // compute defect vector
  template<typename F_>
  void calc_def(std::vector<F_>& vd, const std::vector<F_>& vb, const std::vector<F_>& vx)
  {
    XASSERT(vd.size() == vx.size());
    XASSERT(vb.size() == vx.size());

    const llint n = llint(vx.size());
    const F_ a1 = F_(2*n+2); // main diagonal entry =  2/h =  2/(1/(n+1))
    const F_ a2 = -F_(n+1);  //  off diagonal entry = -1/h = -1/(1/(n+1))

    F_* d = vd.data();
    const F_* b = vb.data();
    const F_* x = vx.data();

    if(n == 1)
    {
      d[0] = b[0] - a1*x[0];
      return;
    }

    d[0] = b[0] - a1*x[0] - a2*x[1];
    #pragma omp parallel for if(n > N_MIN_OMP)
    for(llint i = 1; i < n-1; ++i)
    {
      d[i] = b[i] - a1*x[i] - a2*(x[i-1]+x[i+1]);
    }
    d[n-1] = b[n-1] - a1*x[n-1] - a2*x[n-2];
  }

  // restrict defect vector
  template<typename F_>
  void rest_def(std::vector<F_>& vc, const std::vector<F_>& vf)
  {
    const llint nc = llint(vc.size());
    const llint nf = llint(vf.size());
    XASSERT(nf+1 == 2*(nc+1));

    const F_* xf = vf.data();
    F_* xc = vc.data();

    //     0   1   2   3   4   5   6
    // |---|---|---|---|---|---|---|---|
    //      '\.|./' '\.|./' '\.|./'
    // |-------|-------|-------|-------|
    //         0       1       2
    #pragma omp parallel for if(nc > N_MIN_OMP)
    for(llint i = 0; i < nc; ++i)
    {
      xc[i] = xf[2*i+1] + F_(0.5) * (xf[2*i] + xf[2*i+2]);
    }
  }

  // prolongate correction vector
  template<typename F_>
  void prol_cor(std::vector<F_>& vf, const std::vector<F_>& vc)
  {
    const llint nc = llint(vc.size());
    const llint nf = llint(vf.size());
    XASSERT(nf+1 == 2*(nc+1));

    F_* xf = vf.data();
    const F_* xc = vc.data();

    //         0       1       2
    // |-------|-------|-------|-------|
    //      ./'|'\. ./'|'\. ./'|'\.
    // |---|---|---|---|---|---|---|---|
    //     0   1   2   3   4   5   6
    xf[0] = F_(0.5) * xc[0];
    #pragma omp parallel for if(nc > N_MIN_OMP)
    for(llint i = 0; i < nc-1; ++i)
    {
      xf[2*i+1] = xc[i];
      xf[2*i+2] = F_(0.5) * (xc[i] + xc[i+1]);
    }
    xf[nf-2] = xc[nc-1];
    xf[nf-1] = F_(0.5) * xc[nc-1];
  }

  // class holding four vectors: sol, rhs, def and cor
  template<typename F_>
  class LevelVectors
  {
  public:
    std::vector<F_> vec_sol, vec_rhs, vec_def, vec_cor;

    void create(int level_index)
    {
      const std::size_t n = std::size_t((1 << level_index) - 1);
      vec_sol.resize(n);
      vec_rhs.resize(n);
      vec_def.resize(n);
      vec_cor.resize(n);
    }
  }; // class LevelVectors

  // print timing
  void print_time(const String& what, double t, double total = 0.0)
  {
    std::cout << what.pad_back(25, '.') << ": " << stringify_fp_fix(t, 3, 7);
    if(total > 0.0)
      std::cout << "   [" << stringify_fp_fix(100.0 * t / total, 2, 6) << "%]";
    std::cout << std::endl;
  }

  // print memory usage
  void print_memory(const String& what, double bytes)
  {
    std::cout << what.pad_back(25, '.') << ": " << stringify_fp_fix(bytes / (1024.0*1024.0*1024.0), 3, 7) << " GB" << std::endl;
  }

  /**
   * \brief Performs a multi-precision multigrid benchmark
   *
   * \tparam AFP_
   * Assembly floating point type
   * This type is used for the initial assembly of the RHS vector
   * as well as the H0-error analysis of the solution vector.
   *
   * \tparam OFP_
   * Outer floating point type
   * This type is used for the outer Richardson solver iteration.
   *
   * \tparam IFP_
   * Inner floating point type
   * This type is used for the inner Multigrid preconditioner iteration.
   */
  template<typename AFP_, typename OFP_, typename IFP_>
  void run(const int level_max, const int max_iter, const int num_inner, const int num_smooth, const double omega)
  {
    std::cout << "Floating Types: " << Typo<AFP_>::name() << ' ' << Typo<OFP_>::name() << ' ' << Typo<IFP_>::name() << std::endl;

    // create finest level vectors
    LevelVectors<OFP_> finest;
    finest.create(level_max);

    // create inner MG levels
    std::vector<LevelVectors<IFP_>> levels(static_cast<std::size_t>(level_max));
    for(int i(0); i < level_max; ++i)
      levels.at(std::size_t(i)).create(i+1);

    // create assembly level vectors
    std::vector<AFP_> the_vec_rhs(finest.vec_sol.size());
    std::vector<AFP_> the_vec_sol(finest.vec_sol.size());
    std::vector<AFP_> tmp_vec_sol(finest.vec_sol.size());

    // assemble RHS
    assemble_rhs_sin(the_vec_rhs);

    // initialise finest level stuff
    copy(finest.vec_rhs, the_vec_rhs);
    format(finest.vec_sol);

    // compute RHS norm for relative defect computation
    const OFP_ rhs_norm = norm2(finest.vec_rhs);

    // total runtime stamp
    StopWatch watch_total, watch_inner, watch_smooth, watch_error, watch_transfer;

    // outer multigrid richardson loop
    watch_total.start();
    for(int iter(0); iter <= max_iter; /*++iter*/) // increment in inner loop below
    {
      // update solution vector and compute its error
      watch_error.start();
      copy(the_vec_sol, finest.vec_sol);
      const AFP_ sol_error = analyse_sol_sin(the_vec_sol);
      watch_error.stop();

      // compute defect and its norm
      calc_def(finest.vec_def, finest.vec_rhs, finest.vec_sol);
      const OFP_ def_norm = norm2(finest.vec_def);

      // plot
      std::cout << stringify(iter).pad_front(3) + "    : " << stringify_fp_sci(def_norm) << " / "
        << stringify_fp_sci(def_norm/rhs_norm) << " / " << stringify_fp_sci(sol_error) << std::endl;

      // maximum iterations?
      if(iter >= max_iter)
        break;

      // divergence?
      if(def_norm > OFP_(1E+10))
      {
        std::cout << std::endl << ">>> ERROR: solver diverged!" << std::endl;
        return;
      }

      // initialise multigrid vectors
      copy(levels.back().vec_rhs, finest.vec_def);
      copy(levels.back().vec_def, levels.back().vec_rhs);
      format(levels.back().vec_sol);

      // inner multigrid loop
      watch_inner.start();
      for(int inner(0); (inner < num_inner) && (iter < max_iter); ++inner, ++iter)
      {
        // update defect vector for inner iterations
        if(inner > 0)
        {
          // recompute defect and its norm
          calc_def(levels.back().vec_def, levels.back().vec_rhs, levels.back().vec_sol);
          const OFP_ def_norm_i = OFP_(norm2(levels.back().vec_def));

          // temporarily abuse finest.vec_cor for current vec_sol approximation
          copy(finest.vec_cor, levels.back().vec_sol);
          axpy(finest.vec_cor, finest.vec_sol);

          // update solution vector and compute its error
          watch_error.start();
          copy(the_vec_sol, finest.vec_cor);
          const AFP_ sol_error_i = analyse_sol_sin(the_vec_sol);
          watch_error.stop();

          // plot
          std::cout << stringify(iter).pad_front(3) + " [" << stringify(inner) << "]: "
            << stringify_fp_sci(def_norm_i) << " / " << stringify_fp_sci(def_norm_i/rhs_norm) << " / "
            << stringify_fp_sci(sol_error_i) << std::endl;
        }

        // multigrid cycle restriction
        for(std::size_t lev(levels.size()-1); lev > std::size_t(0); --lev)
        {
          LevelVectors<IFP_>& lvl = levels.at(lev);
          LevelVectors<IFP_>& lvl_c = levels.at(lev-1);
          const IFP_ smooth_damp = IFP_(0.5 * omega) / IFP_(lvl.vec_sol.size()+1);

          // apply pre-smoother
          watch_smooth.start();
          for(llint k(0); k < num_smooth; ++k)
          {
            axpy(lvl.vec_sol, lvl.vec_def, smooth_damp);
            calc_def(lvl.vec_def, lvl.vec_rhs, lvl.vec_sol);
          }
          watch_smooth.stop();

          // restrict defect
          watch_transfer.start();
          rest_def(lvl_c.vec_rhs, lvl.vec_def);
          watch_transfer.stop();
          if(lev > std::size_t(1))
          {
            copy(lvl_c.vec_def, lvl_c.vec_rhs);
            format(lvl_c.vec_sol);
          }
        }

        // "solve" coarse system
        {
          // on the coarsest level, we have a 1x1 system with matrix entry of value 4
          XASSERT(levels.front().vec_sol.size() == std::size_t(1));
          levels.front().vec_sol.front() = IFP_(0.25) * levels.front().vec_rhs.front();
        }

        // multigrid cycle prolongation
        for(std::size_t lev(1); lev < levels.size(); ++lev)
        {
          LevelVectors<IFP_>& lvl = levels.at(lev);
          const IFP_ smooth_damp = IFP_(0.5 * omega) / IFP_(lvl.vec_sol.size()+1);

          // prolongate coarse solution
          watch_transfer.start();
          prol_cor(lvl.vec_cor, levels.at(lev-1).vec_sol);
          watch_transfer.stop();

          // update solution vector
          axpy(lvl.vec_sol, lvl.vec_cor);

          // apply post-smoother
          watch_smooth.start();
          for(llint k(0); k < num_smooth; ++k)
          {
            calc_def(lvl.vec_def, lvl.vec_rhs, lvl.vec_sol);
            axpy(lvl.vec_sol, lvl.vec_def, smooth_damp);
          }
          watch_smooth.stop();
        }
      } // inner loop
      watch_inner.stop();

      // copy solution to finest correction
      copy(finest.vec_cor, levels.back().vec_sol);
      // update solution
      axpy(finest.vec_sol, finest.vec_cor);
    } // outer loop

    watch_total.stop();

    // print timing summary
    std::cout << std::endl << "Timing Statistics" << std::endl;
    print_time("Total Multigrid Time", watch_total.elapsed());
    print_time("Inner Multigrid Time", watch_inner.elapsed(), watch_total.elapsed());
    print_time("Smoothing Time", watch_smooth.elapsed(), watch_total.elapsed());
    print_time("Grid Transfer Time", watch_transfer.elapsed(), watch_total.elapsed());
    print_time("Error Computation Time", watch_error.elapsed(), watch_total.elapsed());
  }

  void main(int argc,char** argv)
  {
    SimpleArgParser args(argc, argv);

    args.support("level");
    args.support("max-iter");
    args.support("num-inner");
    args.support("num-smooth");
    args.support("omega");
    args.support("prec");

    std::deque<std::pair<int,String>> unsupported = args.query_unsupported();
    if(!unsupported.empty())
    {
      std::cerr << std::endl;
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unsupported option #" << (*it).first << " '--" << (*it).second << "'" << std::endl;
      Runtime::abort();
    }

#ifdef _OPENMP
    std::cout << "OpenMP max threads: " << omp_get_max_threads() << std::endl;
#else
    std::cout << "OpenMP not available in this build" << std::endl;
#endif


    int level = 5;
    int max_iter = 20;
    int num_inner = 1;
    int num_smooth = 4;
    double omega = 0.5;

    args.parse("level", level);
    args.parse("max-iter", max_iter);
    args.parse("num-inner", num_inner);
    args.parse("num-smooth", num_smooth);
    args.parse("omega", omega);

    // parse precision types
    String ofp("dp"), ifp("dp");
    if(args.parse("prec", ofp, ifp) < 2)
      ifp = ofp;

    // combine precision strings
    String precs = (ofp + " " + ifp).lower();

    // print argument line
    std::cout << "Args: --level " << level
      << " --max-iter " << max_iter
      << " --num-inner " << num_inner
      << " --num-smooth " << num_smooth
      << " --omega " << omega
      << " --prec " << precs << std::endl;

    // If quad-precision is available, we use that for the assembly if the solver runs in double-precision
#ifdef FEAT_HAVE_QUADMATH
    if(precs == "qp qp") run<f_qp, f_qp, f_qp>(level, max_iter, num_inner, num_smooth, omega); else
    if(precs == "qp dp") run<f_qp, f_qp, f_dp>(level, max_iter, num_inner, num_smooth, omega); else
    // use quad-prec assembly for double-prec solver
    if(precs == "dp dp") run<f_qp, f_dp, f_dp>(level, max_iter, num_inner, num_smooth, omega); else // quad-prec asm
    if(precs == "dp sp") run<f_qp, f_dp, f_sp>(level, max_iter, num_inner, num_smooth, omega); else // quad-prec asm
#else
    // use double-prec assembly for double-prec solver
    if(precs == "dp dp") run<f_dp, f_dp, f_dp>(level, max_iter, num_inner, num_smooth, omega); else // double-prec asm
    if(precs == "dp sp") run<f_dp, f_dp, f_sp>(level, max_iter, num_inner, num_smooth, omega); else // double-prec asm
#endif
    if(precs == "sp sp") run<f_dp, f_sp, f_sp>(level, max_iter, num_inner, num_smooth, omega); else
#ifdef FEAT_HAVE_FLOATX
    if(precs == "dp hp") run<f_dp, f_dp, f_hp>(level, max_iter, num_inner, num_smooth, omega); else
    if(precs == "sp hp") run<f_dp, f_sp, f_hp>(level, max_iter, num_inner, num_smooth, omega); else
    if(precs == "hp hp") run<f_dp, f_hp, f_hp>(level, max_iter, num_inner, num_smooth, omega); else // double-prec asm
#endif
    {
      std::cout << "ERROR: unsupported precision combo " << precs << std::endl;
      Runtime::abort();
    }

    // query memory usage info
    MemoryUsage mem_use;
    print_memory("Peak Physical Memory", double(mem_use.get_peak_physical()));
    print_memory("Peak Virtual  Memory", double(mem_use.get_peak_virtual()));
  }
} // namespace MixedPrecMultiGridBench


int main(int argc,char** argv)
{
  Runtime::initialise(argc, argv);
  MixedPrecMultiGridBench::main(argc, argv);
  return Runtime::finalise();
}
