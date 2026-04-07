// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
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

  explicit MetaMatrixIOTest(PreferredBackend backend) :
    BaseClass("MetaMatrixIOTest", Type::Traits<DataType_>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual void run() const override
  {
    const DataType tol = TestSystem::tol<DataType>();

    // String directory("meta_matrix-io-test.directory/");
    String directory ("./");

    //std::stringstream ioss_diag, ioss_full;

    // generate a test system with PowerDiagMatrix
    typename BaseClass::SystemDiagMatrix mat_diag_write;
    typename BaseClass::SystemVector vec_sol, vec_rhs;
    this->gen_system(7, mat_diag_write, vec_sol, vec_rhs);
    mat_diag_write.write_out(FileMode::fm_mtx, directory + "mat_diag.write.mtx");
    typename BaseClass::SystemDiagMatrix mat_diag_read(FileMode::fm_mtx, directory + "mat_diag.write.mtx");
    //mat_diag_write.write_out(FileMode::fm_mtx, ioss_diag);
    //typename BaseClass::SystemDiagMatrix mat_diag_read(FileMode::fm_mtx, ioss_diag);

    TEST_CHECK_LESS_THAN(mat_diag_write.max_rel_diff(mat_diag_read), tol);

    // generate a test system with PowerFullMatrix
    typename BaseClass::SystemFullMatrix mat_full_write;
    this->gen_system(7, mat_full_write, vec_sol, vec_rhs);
    mat_full_write.write_out(FileMode::fm_mtx, directory + "mat_full.write.mtx");
    typename BaseClass::SystemFullMatrix mat_full_read(FileMode::fm_mtx, directory + "mat_full.write.mtx");
    //mat_full_write.write_out(FileMode::fm_mtx, ioss_full);
    //typename BaseClass::SystemFullMatrix mat_full_read(FileMode::fm_mtx, ioss_full);

    TEST_CHECK_LESS_THAN(mat_full_write.max_rel_diff(mat_full_read), tol);
  }
};

SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, float, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, double, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, float, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, double, std::uint32_t, PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, float, std::uint64_t, PreferredBackend::mkl);
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, double, std::uint64_t, PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, __float128, std::uint64_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, __float128, std::uint32_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, Half, std::uint32_t, PreferredBackend::generic);
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, Half, std::uint64_t, PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, float, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, double, std::uint64_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, float, std::uint32_t, PreferredBackend::cuda);
SPAWN_UNIT_TEST_2T_P(MetaMatrixIOTest, double, std::uint32_t, PreferredBackend::cuda);
#endif
