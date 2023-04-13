// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/lafem/meta_matrix_test_base.hpp>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::TestSystem;

/**
 * \brief Meta-Matrix i/o test class
 *
 * \test The write-out and read-in operations of the following class templates:
 *  - SparseMatrixCSR
 *  - SparseMatrixCOO
 *  - SparseMatrixELL
 *  - PowerColMatrix
 *  - PowerRowMatrix
 *  - PowerDiagMatrix
 *  - PowerFullMatrix
 *  - SaddlePointMatrix
 *
 * \author Christoph Lohmann
 */
template<
  typename DataType_,
  typename IndexType_>
class MetaMatrixIOTest
  : public MetaMatrixTestBase<DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef MetaMatrixTestBase<DataType_, IndexType_> BaseClass;

   MetaMatrixIOTest(PreferredBackend backend) :
    BaseClass("MetaMatrixIOTest", Type::Traits<DataType_>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual ~MetaMatrixIOTest()
  {
  }

  virtual void run() const override
  {
    // String directory("meta_matrix-io-test.directory/");
    String directory ("./");

    // generate a test system with PowerDiagMatrix
    typename BaseClass::SystemDiagMatrix mat_diag_write;
    typename BaseClass::SystemVector vec_sol, vec_rhs;
    this->gen_system(7, mat_diag_write, vec_sol, vec_rhs);
    mat_diag_write.write_out(FileMode::fm_mtx, directory + "mat_diag.write.mtx");
    typename BaseClass::SystemDiagMatrix mat_diag_read(FileMode::fm_mtx, directory + "mat_diag.write.mtx");

    TEST_CHECK_MSG(mat_diag_write == mat_diag_read, "mat_diag_write and mat_diag_read are not the same matrices!");

    // generate a test system with PowerFullMatrix
    typename BaseClass::SystemFullMatrix mat_full_write;
    this->gen_system(7, mat_full_write, vec_sol, vec_rhs);
    mat_full_write.write_out(FileMode::fm_mtx, directory + "mat_full.write.mtx");
    typename BaseClass::SystemFullMatrix mat_full_read(FileMode::fm_mtx, directory + "mat_full.write.mtx");

    TEST_CHECK_MSG(mat_full_write == mat_full_read, "mat_full_write and mat_full_read are not the same matrices!");
  }
};

MetaMatrixIOTest<float, unsigned long> meta_matrix_io_test_generic_float_ulong(PreferredBackend::generic);
MetaMatrixIOTest<double, unsigned long> meta_matrix_io_test_generic_double_ulong(PreferredBackend::generic);
MetaMatrixIOTest<float, unsigned int> meta_matrix_io_test_generic_float_uint(PreferredBackend::generic);
MetaMatrixIOTest<double, unsigned int> meta_matrix_io_test_generic_double_uint(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
MetaMatrixIOTest<float, unsigned long> mkl_meta_matrix_io_test_float_ulong(PreferredBackend::mkl);
MetaMatrixIOTest<double, unsigned long> mkl_meta_matrix_io_test_double_ulong(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
MetaMatrixIOTest<__float128, unsigned long> meta_matrix_io_test_float128_ulong(PreferredBackend::generic);
MetaMatrixIOTest<__float128, unsigned int> meta_matrix_io_test_float128_uint(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
MetaMatrixIOTest<Half, unsigned int> meta_matrix_io_test_half_uint(PreferredBackend::generic);
MetaMatrixIOTest<Half, unsigned long> meta_matrix_io_test_half_ulong(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
MetaMatrixIOTest<float, unsigned long> cuda_meta_matrix_io_test_float_ulong(PreferredBackend::cuda);
MetaMatrixIOTest<double, unsigned long> cuda_meta_matrix_io_test_double_ulong(PreferredBackend::cuda);
MetaMatrixIOTest<float, unsigned int> cuda_meta_matrix_io_test_float_uint(PreferredBackend::cuda);
MetaMatrixIOTest<double, unsigned int> cuda_meta_matrix_io_test_double_uint(PreferredBackend::cuda);
#endif
