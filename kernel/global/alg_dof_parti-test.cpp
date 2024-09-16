// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/lafem/pointstar_factory.hpp>
#include <kernel/global/alg_dof_parti.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/util/random.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

inline int isqrt(const int n)
{
  int m = int(Math::sqrt(double(n)));
  return (m*m == n ? m : -1);
}

template<typename DT_, typename IT_>
class AlgDofPartiTest :
  public UnitTest
{
  typedef DT_ DataType;
  typedef IT_ IndexType;

  typedef LAFEM::VectorMirror<DataType, IndexType> MirrorType;
  typedef LAFEM::DenseVector<DataType, IndexType> LocalVectorType;
  typedef LAFEM::SparseMatrixCSR<DataType, IndexType> LocalMatrixType;

  typedef Global::Gate<LocalVectorType, MirrorType> GateType;
  typedef Global::Vector<LocalVectorType, MirrorType> GlobalVectorType;
  typedef Global::Matrix<LocalMatrixType, MirrorType, MirrorType> GlobalMatrixType;

  typedef Global::AlgDofParti<LocalVectorType, MirrorType> AlgDofPartiType;
  typedef Global::AlgDofPartiVector<LocalVectorType, MirrorType> AlgDofPartiVectorType;
  typedef Global::AlgDofPartiMatrix<LocalMatrixType, MirrorType> AlgDofPartiMatrixType;


public:
  AlgDofPartiTest(PreferredBackend backend) :
    UnitTest("AlgDofPartiTest", Type::Traits<DT_>::name(), Type::Traits<IT_>::name(), backend)
  {
  }

  static MirrorType create_mirror_0(const int n, const int k)
  {
    MirrorType mirror((IT_)n, 1);
    mirror.indices()[0] = IT_(k);
    return mirror;
  }

  static MirrorType create_mirror_1(const int n, const int m, const int o, const int p)
  {
    MirrorType mirror((IT_)n, (IT_)m);
    IndexType* idx = mirror.indices();
    for(int i(0); i < m; ++i)
      idx[IT_(i)] = IT_(o + i * p);
    return mirror;
  }

  static bool create_gate(GateType& gate, const int m)
  {
    const Dist::Comm& comm = *gate.get_comm();

    // get number of processes in each direction
    const int np = int(Math::sqrt(double(comm.size())));
    if(np*np != comm.size())
      return false; // number of procs is not square

    if(np > 1)
    {
      // get our process (i,j) coords
      const int ii = comm.rank() / np;
      const int jj = comm.rank() % np;

      // add mirrors for our vertex neighbors

      // lower-left neighbor?
      if((ii > 0) && (jj > 0))
        gate.push((ii-1)*np + (jj-1), create_mirror_0(m*m, 0));

      // lower-right neighbor?
      if((ii > 0) && (jj+1 < np))
        gate.push((ii-1)*np + (jj+1), create_mirror_0(m*m, m-1));

      // upper-left neighbor?
      if((ii+1 < np) && (jj > 0))
        gate.push((ii+1)*np + (jj-1), create_mirror_0(m*m, m*(m-1)));

      // upper-right neighbor?
      if((ii+1 < np) && (jj+1 < np))
        gate.push((ii+1)*np + (jj+1), create_mirror_0(m*m, m*m-1));

      // add mirror for our edge neighbors

      // lower neighbor?
      if(ii > 0)
        gate.push((ii-1)*np + jj, create_mirror_1(m*m, m, 0, 1));

      // upper neighbor?
      if(ii+1 < np)
        gate.push((ii+1)*np + jj, create_mirror_1(m*m, m, m*(m-1), 1));

      // left neighbor?
      if(jj > 0)
        gate.push(ii*np + (jj-1), create_mirror_1(m*m, m, 0, m));

      // right neighbor?
      if(jj+1 < np)
        gate.push(ii*np + (jj+1), create_mirror_1(m*m, m, m-1, m));
    }

    // compile gate
    gate.compile(LocalVectorType(IT_(m*m)));

    return true;
  }

  virtual void run() const override
  {
    const Dist::Comm comm = Dist::Comm::world();

    // create gate
    GateType gate(comm);
    if(!create_gate(gate, 5))
      return;

    // create alg-dof-parti
    AlgDofPartiType adp;
    adp.assemble_by_gate(gate);
    adp.assemble_allgather(true);

    // test vector
    test_vector(gate, adp);

    // test matrix
    test_matrix(gate, adp, 5);
  }

  void test_vector(const GateType& gate, const AlgDofPartiType& adp) const
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.9));
    Random rng(257ull + 17ull * (unsigned long long)gate.get_comm()->rank());

    // let's create two gate vectors
    GlobalVectorType glob_vec_x(&gate, gate._freqs.clone(LAFEM::CloneMode::Layout));
    GlobalVectorType glob_vec_y(&gate, gate._freqs.clone(LAFEM::CloneMode::Layout));

    // format those vectors to random values
    glob_vec_x.local().format(rng, -1.0, +1.0);
    glob_vec_y.local().format(rng, -1.0, +1.0);

    // synchronize vectors to obtain consistency over all processes
    glob_vec_x.sync_1();
    glob_vec_y.sync_1();

    // compute dot-product of x and y
    const DataType glob_dot_x_y = glob_vec_x.dot(glob_vec_y);

    // create two adp vectors
    AlgDofPartiVectorType adp_vec_x(&adp);
    AlgDofPartiVectorType adp_vec_y(&adp);

    // upload vectors x and y
    adp_vec_x.upload(glob_vec_x);
    adp_vec_y.upload(glob_vec_y);

    // perform dot-product of owned dofs
    const DataType owned_dot_x_y = adp_vec_x.owned().dot(adp_vec_y.owned());

    // sum up owned dots over all processes
    const DataType adp_dot_x_y = gate.sum(owned_dot_x_y);

    // test dot-product error
    const DataType dot_error = Math::abs(adp_dot_x_y - glob_dot_x_y) / Math::abs(glob_dot_x_y);
    TEST_CHECK(dot_error < tol);

    // create another two global vectors
    GlobalVectorType glob_vec_x2(&gate, gate._freqs.clone(LAFEM::CloneMode::Layout));
    GlobalVectorType glob_vec_y2(&gate, gate._freqs.clone(LAFEM::CloneMode::Layout));

    // download from adp vectors
    adp_vec_x.download(glob_vec_x2);
    adp_vec_y.download(glob_vec_y2);

    // compute errors to original vectors
    glob_vec_x2.axpy(glob_vec_x, -DataType(1));
    glob_vec_y2.axpy(glob_vec_y, -DataType(1));

    // compute error
    const DataType download_error = Math::sqrt(glob_vec_x2.norm2sqr() + glob_vec_y2.norm2sqr());
    TEST_CHECK(download_error < tol);
  }

  void test_matrix(/*const*/ GateType& gate, const AlgDofPartiType& adp, const IT_ m) const
  {
    const DataType tol = Math::pow(Math::eps<DataType>(), DataType(0.9));
    Random rng(257ull + 17ull * (unsigned long long)gate.get_comm()->rank());

    // create global matrix
    GlobalMatrixType glob_mat(&gate, &gate);

    // assemble local matrix structure and initialize to random values
    LAFEM::PointstarFactoryFE<DataType, IndexType> psf(m);
    glob_mat.local() = psf.matrix_csr();
    glob_mat.local().format(rng, -1.0, +1.0);

    // scale by frequencies to obtain type-0 matrix
    glob_mat.local().scale_rows(glob_mat.local(), gate._freqs);

    // create two global vectors
    GlobalVectorType glob_vec_x = glob_mat.create_vector_r();
    GlobalVectorType glob_vec_b = glob_mat.create_vector_l();

    // fill x with random values
    glob_vec_x.format(rng, -1.0, +1.0);

    // compute b := A*x
    glob_mat.apply(glob_vec_b, glob_vec_x);

    // create ADP matrix
    AlgDofPartiMatrixType adp_mat(&adp);

    // upload matrix
    adp_mat.upload(glob_mat);

    // create two ADP vectors
    AlgDofPartiVectorType adp_vec_x(&adp);
    AlgDofPartiVectorType adp_vec_b(&adp);

    // upload vector x
    adp_vec_x.upload(glob_vec_x);

    // compute b := A*x
    adp_mat.apply(adp_vec_b, adp_vec_x);

    // create another global vector and download b
    GlobalVectorType glob_vec_b2 = glob_mat.create_vector_l();
    adp_vec_b.download(glob_vec_b2);

    // compute error vector
    glob_vec_b2.axpy(glob_vec_b, -DataType(1));
    const DataType error_norm = glob_vec_b2.norm2();
    TEST_CHECK(error_norm < tol);
  }
};

AlgDofPartiTest <float, std::uint32_t> alg_dof_parti_test_float_uint32(PreferredBackend::generic);
AlgDofPartiTest <double, std::uint32_t> alg_dof_parti_test_double_uint32(PreferredBackend::generic);
AlgDofPartiTest <float, std::uint64_t> alg_dof_parti_test_float_uint64(PreferredBackend::generic);
AlgDofPartiTest <double, std::uint64_t> alg_dof_parti_test_double_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
AlgDofPartiTest <float, std::uint64_t> mkl_alg_dof_parti_test_float_uint64(PreferredBackend::mkl);
AlgDofPartiTest <double, std::uint64_t> mkl_alg_dof_parti_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
AlgDofPartiTest <__float128, std::uint32_t> alg_dof_parti_test_float128_uint32(PreferredBackend::generic);
AlgDofPartiTest <__float128, std::uint64_t> alg_dof_parti_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
AlgDofPartiTest <Half, std::uint32_t> alg_dof_parti_test_half_uint32(PreferredBackend::generic);
AlgDofPartiTest <Half, std::uint64_t> alg_dof_parti_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
AlgDofPartiTest <float, std::uint32_t> cuda_alg_dof_parti_test_float_uint32(PreferredBackend::cuda);
AlgDofPartiTest <double, std::uint32_t> cuda_alg_dof_parti_test_double_uint32(PreferredBackend::cuda);
AlgDofPartiTest <float, std::uint64_t> cuda_alg_dof_parti_test_float_uint64(PreferredBackend::cuda);
AlgDofPartiTest <double, std::uint64_t> cuda_alg_dof_parti_test_double_uint64(PreferredBackend::cuda);
#endif
