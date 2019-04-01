// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
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
template<typename DataType_, typename IndexType_>
class MetaMatrixIOTest
  : public MetaMatrixTestBase<Mem::Main, DataType_, IndexType_>
{
public:
  typedef DataType_ DataType;
  typedef MetaMatrixTestBase<Mem::Main, DataType_, IndexType_> BaseClass;

   MetaMatrixIOTest() :
    BaseClass("MetaMatrixIOTest")
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

//MetaMatrixIOTest<float, Index> meta_matrix_io_test_generic_float;
//MetaMatrixIOTest<double, Index> meta_matrix_io_test_generic_double;
