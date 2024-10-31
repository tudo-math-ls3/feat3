// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

// this function implements HFEM direct solving idea 2, that is solving via
// v = Attt^-1 * (g - B^t * f - D * C^-1 * h)
// u = f - B * v
// w = C^-1 * (h - D^t * v)
// where Attt = -B^t + E - D * C^-1 * D^t

//usage: ./hfem-direct --dir-path /home/user/dribbroc/nobackup1/HFEM_FaS --backend generic --nrhs_d 2 --nrhs_f 2 --nrhs_h 2 --datatype double --basename hfem_6_1
//usage: ./hfem-direct --dir-path /home/user/dribbroc/nobackup1/HFEM_direkt --backend generic --nrhs_d 2 --nrhs_f 2 --nrhs_h 2 --datatype double  --basename hfem_256_16

#include <kernel/base_header.hpp>
#include <kernel/util/string.hpp>
#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_matrix.hpp>
#include <kernel/solver/richardson.hpp>

// -=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

using namespace FEAT;


namespace HFEM_direct
{

  template <typename DT_, typename IT_, typename DT2_, typename IT2_>
  void populate(LAFEM::DenseMatrix<DT_, IT_> & target, const LAFEM::DenseVector<DT2_, IT2_> & source)
  {
    XASSERT(target.rows() == source.size());

    for (Index row(0) ; row < target.rows() ; ++row)
    {
      MemoryPool::set_memory(target.elements() + row * target.columns() , DT_(source(row)), target.columns());
    }
  }

  template <typename DT_, typename IT_>
  void run_solver(PreferredBackend backend, Index nrhs,
      LAFEM::SparseMatrixCSR<double, Index> & D_orig,
      LAFEM::SparseMatrixCSR<double, Index> & Dt_orig,
      LAFEM::DenseMatrix<double, Index> & Atttinvd_orig,
      LAFEM::DenseMatrix<double, Index> & Bd_orig,
      LAFEM::DenseMatrix<double, Index> & Btd_orig,
      std::vector<std::shared_ptr<LAFEM::DenseMatrix<double, Index >> >& vCsubinvd_orig,
      LAFEM::DenseVector<double, Index> & f_orig,
      LAFEM::DenseVector<double, Index> & g_orig,
      LAFEM::DenseVector<double, Index> & h_orig,
      LAFEM::DenseVector<Index, Index> & distrib,
      LAFEM::DenseVector<double, Index> & vref_orig)
  {
    (void)vref_orig; //disable unused warning
    Backend::set_preferred_backend(backend);
    std::cout<<"========================================"<<"\n";
    std::cout<<"DataType: " << Type::Traits<DT_>::name()<<"\n";
    std::cout<<"nrhs: "<<nrhs<<"\n";

    LAFEM::SparseMatrixCSR<DT_, IT_> D;
    D.convert(D_orig);
    LAFEM::SparseMatrixCSR<DT_, IT_> Dt;
    Dt.convert(Dt_orig);
    LAFEM::DenseMatrix<DT_, IT_> Atttinvd;
    Atttinvd.convert(Atttinvd_orig);
    LAFEM::DenseMatrix<DT_, IT_> Bd;
    Bd.convert(Bd_orig);
    LAFEM::DenseMatrix<DT_, IT_> Btd;
    Btd.convert(Btd_orig);
    std::vector<std::shared_ptr<LAFEM::DenseMatrix<DT_, IT_> > > vCsubinvd;
    for (Index i(0) ; i < vCsubinvd_orig.size() ; ++i)
    {
      vCsubinvd.push_back(std::make_shared<LAFEM::DenseMatrix<DT_, IT_>>());
      vCsubinvd.at(i)->convert(*(vCsubinvd_orig.at(i)));
    }

    LAFEM::DenseMatrix<DT_, IT_> f(f_orig.size(), nrhs);
    populate(f, f_orig);
    LAFEM::DenseMatrix<DT_, IT_> g(g_orig.size(), nrhs);
    populate(g, g_orig);
    LAFEM::DenseMatrix<DT_, IT_> h(h_orig.size(), nrhs);
    populate(h, h_orig);

    LAFEM::DenseMatrix<DT_, IT_> u(f.rows(), nrhs);
    u.copy(f);
    LAFEM::DenseMatrix<DT_, IT_> v(g.rows(), nrhs);
    LAFEM::DenseMatrix<DT_, IT_> w(h.rows(), nrhs);
    LAFEM::DenseMatrix<DT_, IT_> tmpvec(w.rows(), w.columns());
    LAFEM::DenseMatrix<DT_, IT_> htmp(h.rows(), h.columns());
    htmp.copy(h);
    LAFEM::DenseMatrix<DT_, IT_> vtmp(v.rows(), h.columns());

    TimeStamp at;
    Index iters(25);
    for (Index iter(0) ; iter < iters ; ++iter)
    {
      ////////////// SOLVE
      // STEP 1: v = Attt^-1 * (g - B^t * f - D * C^-1 * h)

      //axpy_diag(T(1.f), Csubinv, h, T(0.f), tmpvec); // tmpvec = C^-1 * h
      //Index nblocks = h.rows() / Csubinvd.columns();
      Index coffset(0);
      Index roffset(0);
      for (Index submatrix(0); submatrix < distrib.size(); ++submatrix)
      {
        for (Index block(0) ; block < distrib(submatrix) ; ++block)
        {
          //LAFEM::Arch::ProductMatMat::dense(tmpvec.elements() + block * Csubinvd.rows(), DT_(1.0), DT_(0.0), Csubinvd.elements(), h.elements() + block * Csubinvd.columns(), tmpvec.elements() + block * Csubinvd.rows(), Csubinvd.rows(), nrhs, Csubinvd.columns());
          LAFEM::Arch::ProductMatMat::dense(tmpvec.elements() + roffset, DT_(1.0), DT_(0.0), vCsubinvd.at(submatrix)->elements(), h.elements() + coffset, tmpvec.elements() + roffset, vCsubinvd.at(submatrix)->rows(), h.columns(), vCsubinvd.at(submatrix)->columns());
          coffset += vCsubinvd.at(submatrix)->columns();
          roffset += vCsubinvd.at(submatrix)->rows();
        }
      }

      //axpy(T(1.f), D, tmpvec, T(0.f), v);            // v = D * tmpvec
      v.multiply(D, tmpvec);

      //axpy(T(-1.0), Bt, f, T(-1.0), v);              // v = -B^t * f - v
      //LAFEM::Arch::Apply::dense(v.elements(), DT_(-1.), DT_(-1.), v.elements(), Btd.elements(),
      //                           f.elements(), Btd.rows(), Btd.columns());
      //LAFEM::Arch::ProductMatMat::dense(v.elements(), DT_(-1.), DT_(-1.), Btd.elements(), f.elements(), v.elements(), Btd.rows(), Btd.columns(), Btd.columns());
      v.multiply(Btd, f, v, DT_(-1.), DT_(-1.));

      //axpy(T(1.0), g, v);                            // v = g + v
      v.axpy(g);

      //vtmp.copy_from(v);                             // vtmp = v
      vtmp.copy(v);

      //axpy(T(1.0), tttAinv, vtmp, T(0.0), v);        // v = Attt^-1 * vtmp*/
      v.multiply(Atttinvd, vtmp);

      /*MemoryPool::synchronize();
      for (Index i(0) ; i < 25 ; ++i)
      {
        std::cout<<v(i,0) << " " << vref_orig(i)<<"\n";
      }
      exit(0);*/


      // STEP 2: u = f - B * v
      // note: u is already equal to f
      //axpy(T(-1.0), B, v, T(1.0), u);                // u = -B * v + u
      //Bd.apply(u, v, u, DT_(-1.));
      //LAFEM::Arch::ProductMatMat::dense(u.elements(), DT_(-1.), DT_(1.), Bd.elements(), v.elements(), u.elements(), Btd.rows(), Btd.columns(), Bd.columns());
      u.multiply(Bd, v, u, DT_(-1));

      // STEP 3: w = C^-1 * (h - D^t * v)
      // note: htmp is equal to h
      //axpy(T(-1.0), Dt, v, T(1.0), htmp);            // htmp = -D^t * v + htmp
      //Dt.apply(htmp, v, htmp, DT_(-1.));
      htmp.multiply(Dt, v, DT_(-1.));

      //axpy_diag(T(1.0), Csubinv, htmp, T(0.0), w);   // w = C^-1 * htmp
      /*for (Index block(0); block < nblocks; ++block)
      {
        LAFEM::Arch::ProductMatMat::dense(w.elements() + block * Csubinvd.rows(), DT_(1.0), DT_(0.0), Csubinvd.elements(), htmp.elements() + block * Csubinvd.columns(), w.elements() + block * Csubinvd.rows(), Csubinvd.rows(), nrhs, Csubinvd.columns());
      }*/
      coffset = 0;
      roffset = 0;
      for (Index submatrix(0); submatrix < distrib.size(); ++submatrix)
      {
        for (Index block(0) ; block < distrib(submatrix) ; ++block)
        {
          LAFEM::Arch::ProductMatMat::dense(w.elements() + roffset, DT_(1.0), DT_(0.0), vCsubinvd.at(submatrix)->elements(), htmp.elements() + coffset, w.elements() + roffset, vCsubinvd.at(submatrix)->rows(), htmp.columns(), vCsubinvd.at(submatrix)->columns());
          coffset += vCsubinvd.at(submatrix)->columns();
          roffset += vCsubinvd.at(submatrix)->rows();
        }
      }

      if (iter == 0)
      {
        MemoryPool::synchronize();
        at.stamp();
      }
    }
    MemoryPool::synchronize();

    /*double norm(0.);
    for (Index i(0) ; i < 25 ; ++i)
    {
      norm += w(i,0);
    }
    std::cout<<"Norm: " << norm << "\n";*/

    double average(at.elapsed_now() / double(iters-1));
    std::cout<<"TOE [s]: "<<average<<"\n";
    double dofs(double((f_orig.size() + g_orig.size() + h_orig.size()) * nrhs));
    std::cout<<"MDOFs: " << dofs / 1e6 << "\n";
    double dofsps(dofs / average);
    std::cout<<"MDOF/s: "<<dofsps / 1e6 <<"\n";
  } // void run_solver

  void main(int argc, char* argv[])
  {
    SimpleArgParser args(argc, argv);
    args.support("dir-path");
    args.support("basename");
    args.support("backend");
    args.support("nrhs_d");
    args.support("nrhs_f");
    args.support("nrhs_h");
    args.support("datatype");

    std::deque<std::pair<int,String>> unsupported = args.query_unsupported();
    if(!unsupported.empty())
    {
      std::cerr << "\n";
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unsupported option #" << (*it).first << " '--" << (*it).second << "'" << "\n";
      Runtime::abort();
    }

    String dir_path("");
    if(args.check("dir-path") > 0)
    {
      args.parse("dir-path", dir_path);
    }

    String basename("");
    if(args.check("basename") > 0)
    {
      args.parse("basename", basename);
    }

    Index nrhs_d(2);
    if(args.check("nrhs_d") > 0)
    {
      args.parse("nrhs_d", nrhs_d);
    }

    Index nrhs_f(2);
    if(args.check("nrhs_f") > 0)
    {
      args.parse("nrhs_f", nrhs_f);
    }

    Index nrhs_h(2);
    if(args.check("nrhs_h") > 0)
    {
      args.parse("nrhs_h", nrhs_h);
    }

    String sdatatype("");
    if(args.check("datatype") > 0)
    {
      args.parse("datatype", sdatatype);
    }

    std::vector<PreferredBackend> backends;
    PreferredBackend backend = PreferredBackend::generic;
    if(args.check("backend") > 0)
    {
      String sel_backend;
      args.parse("backend", sel_backend);
      if (sel_backend == "generic")
      {
            backend = PreferredBackend::generic;
      }
      else if (sel_backend == "cuda")
      {
            backend = PreferredBackend::cuda;
      }
      else
      {
        throw InternalError("Unsupported backend: " + sel_backend + "!");
      }
    }

    std::cout<<"config basename: " << basename<<"\n";
    std::cout<<"preloading Data..."<<"\n";
    LAFEM::SparseMatrixCSR<double, Index> D_orig(LAFEM::FileMode::fm_csr, dir_path + "/" + basename + "_D.csr");
    LAFEM::SparseMatrixCSR<double, Index> Dt_orig(LAFEM::FileMode::fm_csr, dir_path + "/" + basename + "_Dt.csr");
    LAFEM::DenseMatrix<double, Index> Atttinvd_orig(LAFEM::FileMode::fm_dm, dir_path + "/" + basename + "_Atttinvd.dm");
    LAFEM::DenseMatrix<double, Index> Bd_orig(LAFEM::FileMode::fm_dm, dir_path + "/" + basename + "_Bd.dm");
    LAFEM::DenseMatrix<double, Index> Btd_orig(LAFEM::FileMode::fm_dm, dir_path + "/" + basename + "_Btd.dm");
    LAFEM::DenseVector<Index, Index> distrib(LAFEM::FileMode::fm_mtx, dir_path + "/" + basename + "_Cdistrib.mtx");
    std::vector<std::shared_ptr<LAFEM::DenseMatrix<double, Index> > > vCsubinvd_orig;
    for (Index i(0) ; i < distrib.size() ; ++i)
    {
      vCsubinvd_orig.push_back(std::make_shared<LAFEM::DenseMatrix<double, Index> >(LAFEM::FileMode::fm_dm, dir_path + "/" + basename + "_Csubinv" + stringify(i+1) + "d.dm"));
    }
    LAFEM::DenseVector<double, Index> f_orig(LAFEM::FileMode::fm_dv, dir_path + "/" + basename + "_f.dv");
    LAFEM::DenseVector<double, Index> g_orig(LAFEM::FileMode::fm_dv, dir_path + "/" + basename + "_g.dv");
    LAFEM::DenseVector<double, Index> h_orig(LAFEM::FileMode::fm_dv, dir_path + "/" + basename + "_h.dv");
    LAFEM::DenseVector<double, Index> vref_orig;//(LAFEM::FileMode::fm_dv, dir_path + "/" + basename + "_vref.dv");
    std::cout<<"finished"<<"\n";
    std::cout<<"datatype selection: "<<sdatatype<<"\n";
    std::cout<<"cells: "<<distrib.size()<<"\n";
    std::cout<<"blocks: "<<distrib<<"\n";
    std::cout<<"Selected backend: " << backend << "\n";

    typedef unsigned int IT_;

    if (sdatatype == "double")
      run_solver<double, IT_>(backend, nrhs_d, D_orig, Dt_orig, Atttinvd_orig, Bd_orig, Btd_orig, vCsubinvd_orig, f_orig, g_orig, h_orig, distrib, vref_orig);
    else if (sdatatype == "float")
      run_solver<float, IT_>(backend, nrhs_f, D_orig, Dt_orig, Atttinvd_orig, Bd_orig, Btd_orig, vCsubinvd_orig, f_orig, g_orig, h_orig, distrib, vref_orig);
#ifdef FEAT_HAVE_HALFMATH
    else if (sdatatype == "half")
      run_solver<Half, IT_>(backend, nrhs_h, D_orig, Dt_orig, Atttinvd_orig, Bd_orig, Btd_orig, vCsubinvd_orig, f_orig, g_orig, h_orig, distrib, vref_orig);
#endif
    else if (sdatatype == "all")
    {
      run_solver<double, IT_>(backend, 2, D_orig, Dt_orig, Atttinvd_orig, Bd_orig, Btd_orig, vCsubinvd_orig, f_orig, g_orig, h_orig, distrib, vref_orig);
      //FEAT::Util::cuda_reset_algos();
      run_solver<double, IT_>(backend, nrhs_d, D_orig, Dt_orig, Atttinvd_orig, Bd_orig, Btd_orig, vCsubinvd_orig, f_orig, g_orig, h_orig, distrib, vref_orig);
      //FEAT::Util::cuda_reset_algos();
      run_solver<float, IT_>(backend, 2, D_orig, Dt_orig, Atttinvd_orig, Bd_orig, Btd_orig, vCsubinvd_orig, f_orig, g_orig, h_orig, distrib, vref_orig);
      //FEAT::Util::cuda_reset_algos();
      run_solver<float, IT_>(backend, nrhs_f, D_orig, Dt_orig, Atttinvd_orig, Bd_orig, Btd_orig, vCsubinvd_orig, f_orig, g_orig, h_orig, distrib, vref_orig);
      //FEAT::Util::cuda_reset_algos();
#ifdef FEAT_HAVE_HALFMATH
      run_solver<Half, IT_>(backend, 2, D_orig, Dt_orig, Atttinvd_orig, Bd_orig, Btd_orig, vCsubinvd_orig, f_orig, g_orig, h_orig, distrib, vref_orig);
      //FEAT::Util::cuda_reset_algos();
      run_solver<Half, IT_>(backend, nrhs_h, D_orig, Dt_orig, Atttinvd_orig, Bd_orig, Btd_orig, vCsubinvd_orig, f_orig, g_orig, h_orig, distrib, vref_orig);
#endif
    }
    else
        throw InternalError("Unsupported datatype: " + sdatatype + "!");

  } // int main(...)
} // namespace HFEM_direct

// Here's our main function
int main(int argc, char* argv[])
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  HFEM_direct::main(argc, argv);
  return 0;
}
