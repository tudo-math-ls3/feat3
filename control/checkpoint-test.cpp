// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>
#include <control/checkpoint_control.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/binary_stream.hpp>
#include <test_system/test_system.hpp>
#include <kernel/lafem/power_row_matrix.hpp>
#include <kernel/lafem/power_col_matrix.hpp>
#include <kernel/lafem/power_full_matrix.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/tuple_diag_matrix.hpp>
#include <kernel/lafem/tuple_matrix.hpp>
#include <kernel/lafem/tuple_vector.hpp>
#include <kernel/lafem/sparse_matrix_coo.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/pointstar_factory.hpp>



using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Test class for the checkpoint class.
 *
 * \test test description missing
 */

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class CheckpointTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  CheckpointTest()
    : FullTaggedTest<Mem_, DT_, IT_>("CheckpointTest")
  {
  }

  virtual ~CheckpointTest()
  {
  }


  virtual void run() const override
  {
    LAFEM::DenseVector<Mem_, DT_, IT_> dv1(1234);
    for (Index i(0) ; i < dv1.size() ; ++i)
      dv1(i, DT_(i) / DT_(12));

    //LAFEM::PowerRowMatrix<LAFEM::PowerRowMatrix<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>, 2>, 2> next_powrow;
    //next_powrow.get(0,1) = powrow.clone();
    //next_powrow.get(0,0) = powrow.clone();

    //simple test
    {
      auto comm = Dist::Comm::world();
      Control::CheckpointControl cp(comm);
      LAFEM::SerialConfig config(false, false);
      cp.set_config(config);
      BinaryStream bs;
      cp.add_object(String("dv1"), dv1);
      cp.save(bs);
      comm.barrier();
      bs.seekg(0);
      std::uint64_t size_1 = *(std::uint64_t *)(bs.data()); //read in the size of the output data
      XASSERTM(size_1 > 0, "Read in compressed size is zero");
      bs.seekg(0);
      cp.load(bs);
      LAFEM::DenseVector<Mem_, DT_, IT_> dv2;
      cp.restore_object(String("dv1"), dv2, false);
      TEST_CHECK_EQUAL(dv1, dv2);
      TEST_CHECK_NOT_EQUAL(cp.get_identifier_list().find("dv1"), std::string::npos);
#ifdef FEAT_HAVE_ZLIB
      LAFEM::SerialConfig config2(true, false);
      cp.clear_input();
      cp.set_config(config2);
      BinaryStream bs2;
      cp.save(bs2);
      comm.barrier();
      bs2.seekg(0);
      std::uint64_t size_2 = *(std::uint64_t *)(bs2.data()); //read in the size of the output data
      XASSERTM(size_2 < size_1, "With zlib compressed Checkpoint size is bigger than uncompressed");
      //std::cout << "Size_1 is:" << size_1 << "   size_2 is:" << size_2 << std::endl;*/
      bs2.seekg(0);
      cp.load(bs2);
      LAFEM::DenseVector<Mem_, DT_, IT_> dv3;
      cp.restore_object(String("dv1"), dv3, false);
      TEST_CHECK_EQUAL(dv1, dv3);
#ifdef FEAT_HAVE_ZFP
      LAFEM::SerialConfig config3(true, true, 1e-4);
      cp.clear_input();
      cp.set_config(config3);
      BinaryStream bs3;
      cp.save(bs3);
      comm.barrier();
      bs3.seekg(0);
      std::uint64_t size_3 = *(std::uint64_t *)(bs3.data()); //read in the size of the output data
      XASSERTM(size_3 < size_2, "With zfp compressed Checkpoint size is bigger than zlib compressed");
      //std::cout << "Size_2 is:" << size_2 << "   size_3 is:" << size_3 << std::endl;
      bs3.seekg(0);
      cp.load(bs3);
      LAFEM::DenseVector<Mem_, DT_, IT_> dv4;
      cp.restore_object(String("dv1"), dv4, false);
      TEST_CHECK_EQUAL(dv1.size(), dv4.size());
      for(Index i(0) ; i < dv1.size() ; ++i)
        TEST_CHECK_EQUAL_WITHIN_EPS(dv1(i), dv2(i), DT_(1e-3));
#endif
#endif

    }
  }
};
CheckpointTest<Mem::Main, float, unsigned int> checkpoint_test_float_uint;
CheckpointTest<Mem::Main, float, unsigned long> checkpoint_test_float_ulong;
CheckpointTest<Mem::Main, double, unsigned int> checkpoint_test_double_uint;
CheckpointTest<Mem::Main, double, unsigned long> checkpoint_test_double_ulong;

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class CheckpointPowerRowTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  CheckpointPowerRowTest()
    : FullTaggedTest<Mem_, DT_, IT_>("CheckpointPowerRowTest")
  {
  }

  virtual ~CheckpointPowerRowTest()
  {
  }


  virtual void run() const override
  {

    LAFEM::PowerRowMatrix<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>, 2> powrow;
    LAFEM::PointstarFactoryFD<DT_, IT_> five_star(40,2);
    LAFEM::PointstarFactoryFD<DT_, IT_> nine_star(40,3);
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star1 = five_star.matrix_csr();
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star2 = nine_star.matrix_csr();
    powrow.get(0,0) = star1.clone();
    powrow.get(0,1) = star2.clone();

    auto comm = Dist::Comm::world();
    Control::CheckpointControl cp(comm);
    BinaryStream bs;
    cp.add_object(String("powrow"), powrow);
    cp.save(bs);
    comm.barrier();
    bs.seekg(0);
    cp.load(bs);
    LAFEM::PowerRowMatrix<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>, 2> powrow2;
    cp.restore_object(String("powrow"), powrow2, false);
    TEST_CHECK_EQUAL(powrow.get(0,0), powrow2.get(0,0));
    TEST_CHECK_EQUAL(powrow.get(0,1), powrow2.get(0,1));
    TEST_CHECK_NOT_EQUAL(cp.get_identifier_list().find("powrow"), std::string::npos);
  }
};
CheckpointPowerRowTest<Mem::Main, float, unsigned int> checkpoint_power_row_test_float_uint;
CheckpointPowerRowTest<Mem::Main, float, unsigned long> checkpoint_power_row_test_float_ulong;
CheckpointPowerRowTest<Mem::Main, double, unsigned int> checkpoint_power_row_test_double_uint;
CheckpointPowerRowTest<Mem::Main, double, unsigned long> checkpoint_power_row_test_double_ulong;

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class CheckpointPowerColumnTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  CheckpointPowerColumnTest()
    : FullTaggedTest<Mem_, DT_, IT_>("CheckpointPowerColumnTest")
  {
  }

  virtual ~CheckpointPowerColumnTest()
  {
  }


  virtual void run() const override
  {
    LAFEM::PowerColMatrix<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>, 2> powcol;
    LAFEM::PointstarFactoryFD<DT_, IT_> five_star(40,2);
    LAFEM::PointstarFactoryFD<DT_, IT_> nine_star(40,3);
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star1 = five_star.matrix_csr();
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star2 = nine_star.matrix_csr();
    powcol.get(0,0) = star1.clone();
    powcol.get(1,0) = star2.clone();

    auto comm = Dist::Comm::world();
    Control::CheckpointControl cp(comm);
    BinaryStream bs;
    cp.add_object(String("powcol"), powcol);
    cp.save(bs);
    comm.barrier();
    bs.seekg(0);
    cp.load(bs);
    LAFEM::PowerColMatrix<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>, 2> powcol2;
    cp.restore_object(String("powcol"), powcol2, false);
    TEST_CHECK_EQUAL(powcol.get(0,0), powcol2.get(0,0));
    TEST_CHECK_EQUAL(powcol.get(1,0), powcol2.get(1,0));
    TEST_CHECK_NOT_EQUAL(cp.get_identifier_list().find("powcol"), std::string::npos);
  }
};
CheckpointPowerColumnTest<Mem::Main, float, unsigned int> checkpoint_power_column_test_float_uint;
CheckpointPowerColumnTest<Mem::Main, float, unsigned long> checkpoint_power_column_test_float_ulong;
CheckpointPowerColumnTest<Mem::Main, double, unsigned int> checkpoint_power_column_test_double_uint;
CheckpointPowerColumnTest<Mem::Main, double, unsigned long> checkpoint_power_column_test_double_ulong;

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class CheckpointPowerFullTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  CheckpointPowerFullTest()
    : FullTaggedTest<Mem_, DT_, IT_>("CheckpointPowerFullTest")
  {
  }

  virtual ~CheckpointPowerFullTest()
  {
  }


  virtual void run() const override
  {
    LAFEM::PowerFullMatrix<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>, 2,2> powful;
    LAFEM::PointstarFactoryFD<DT_, IT_> five_star(40,2);
    LAFEM::PointstarFactoryFD<DT_, IT_> nine_star(40,3);
    LAFEM::PointstarFactoryFD<DT_, IT_> nine_star2(30,3);
    LAFEM::PointstarFactoryFD<DT_, IT_> five_star2(30,2);
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star1 = five_star.matrix_csr();
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star2 = nine_star.matrix_csr();
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star3 = five_star2.matrix_csr();
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star4 = nine_star2.matrix_csr();
    powful.get(0,0) = star1.clone();
    powful.get(0,1) = star2.clone();
    powful.get(1,0) = star3.clone();
    powful.get(1,1) = star4.clone();

    auto comm = Dist::Comm::world();
    Control::CheckpointControl cp(comm);
    BinaryStream bs;
    cp.add_object(String("powful"), powful);
    cp.save(bs);
    comm.barrier();
    bs.seekg(0);
    cp.load(bs);
    LAFEM::PowerFullMatrix<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>, 2,2> powful2;
    cp.restore_object(String("powful"), powful2, false);
    TEST_CHECK_EQUAL(powful.get(0,0), powful2.get(0,0));
    TEST_CHECK_EQUAL(powful.get(1,0), powful2.get(1,0));
    TEST_CHECK_EQUAL(powful.get(0,1), powful2.get(0,1));
    TEST_CHECK_EQUAL(powful.get(1,1), powful2.get(1,1));
    TEST_CHECK_NOT_EQUAL(cp.get_identifier_list().find("powful"), std::string::npos);
  }
};
CheckpointPowerFullTest<Mem::Main, float, unsigned int> checkpoint_power_full_test_float_uint;
CheckpointPowerFullTest<Mem::Main, float, unsigned long> checkpoint_power_full_test_float_ulong;
CheckpointPowerFullTest<Mem::Main, double, unsigned int> checkpoint_power_full_test_double_uint;
CheckpointPowerFullTest<Mem::Main, double, unsigned long> checkpoint_power_full_test_double_ulong;

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class CheckpointPowerVectorTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  CheckpointPowerVectorTest()
    : FullTaggedTest<Mem_, DT_, IT_>("CheckpointPowerVectorTest")
  {
  }

  virtual ~CheckpointPowerVectorTest()
  {
  }


  virtual void run() const override
  {
    LAFEM::DenseVector<Mem_, DT_, IT_> dv1(1234);
    for (Index i(0) ; i < dv1.size() ; ++i)
      dv1(i, DT_(i) / DT_(12));

    LAFEM::DenseVector<Mem_, DT_, IT_> dv2(1289);
    for (Index i(0) ; i < dv2.size() ; ++i)
      dv2(i, DT_(i) / DT_(16));
    LAFEM::PowerVector<LAFEM::DenseVector<Mem_,DT_,IT_>, 2> powvec;
    powvec.get(0) = dv1.clone();
    powvec.get(1) = dv2.clone();

    auto comm = Dist::Comm::world();
    Control::CheckpointControl cp(comm);
    BinaryStream bs;
    cp.add_object(String("powvec"), powvec);
    cp.save(bs);
    comm.barrier();
    bs.seekg(0);
    cp.load(bs);
    LAFEM::PowerVector<LAFEM::DenseVector<Mem_, DT_, IT_>, 2> powvec2;
    cp.restore_object(String("powvec"), powvec2, false);
    TEST_CHECK_EQUAL(powvec.get(0), powvec2.get(0));
    TEST_CHECK_EQUAL(powvec.get(1), powvec2.get(1));
    TEST_CHECK_NOT_EQUAL(cp.get_identifier_list().find("powvec"), std::string::npos);
  }
};
CheckpointPowerVectorTest<Mem::Main, float, unsigned int> checkpoint_power_vector_test_float_uint;
CheckpointPowerVectorTest<Mem::Main, float, unsigned long> checkpoint_power_vector_test_float_ulong;
CheckpointPowerVectorTest<Mem::Main, double, unsigned int> checkpoint_power_vector_test_double_uint;
CheckpointPowerVectorTest<Mem::Main, double, unsigned long> checkpoint_power_vector_test_double_ulong;

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class CheckpointSaddlePointMatrixTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  CheckpointSaddlePointMatrixTest()
    : FullTaggedTest<Mem_, DT_, IT_>("CheckpointSaddlePointMatrixTest")
  {
  }

  virtual ~CheckpointSaddlePointMatrixTest()
  {
  }


  virtual void run() const override
  {
    LAFEM::SaddlePointMatrix<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>, LAFEM::SparseMatrixCOO<Mem_, DT_, IT_>> saddle;
    LAFEM::PointstarFactoryFD<DT_, IT_> five_star(40,2);
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star1 = five_star.matrix_csr();
    LAFEM::SparseMatrixCOO<Mem_, DT_, IT_> star2(star1);

    saddle.block_a() = star1.clone();
    saddle.block_b() = star2.clone();
    saddle.block_d() = star2.clone();

    auto comm = Dist::Comm::world();
    Control::CheckpointControl cp(comm);
    BinaryStream bs;
    cp.add_object(String("saddle"), saddle);
    cp.save(bs);
    comm.barrier();
    bs.seekg(0);
    cp.load(bs);
    LAFEM::SaddlePointMatrix<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>, LAFEM::SparseMatrixCOO<Mem_, DT_, IT_>> saddle2;
    cp.restore_object(String("saddle"), saddle2, false);
    TEST_CHECK_EQUAL(saddle.block_a(), saddle2.block_a());
    TEST_CHECK_EQUAL(saddle.block_b(), saddle2.block_b());
    TEST_CHECK_EQUAL(saddle.block_d(), saddle2.block_d());
    TEST_CHECK_NOT_EQUAL(cp.get_identifier_list().find("saddle"), std::string::npos);
  }
};
CheckpointSaddlePointMatrixTest<Mem::Main, float, unsigned int> checkpoint_saddle_point_matrix_test_float_uint;
CheckpointSaddlePointMatrixTest<Mem::Main, float, unsigned long> checkpoint_saddle_point_matrix_test_float_ulong;
CheckpointSaddlePointMatrixTest<Mem::Main, double, unsigned int> checkpoint_saddle_point_matrix_test_double_uint;
CheckpointSaddlePointMatrixTest<Mem::Main, double, unsigned long> checkpoint_saddle_point_matrix_test_double_ulong;

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class CheckpointTupleDiagMatrixTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  CheckpointTupleDiagMatrixTest()
    : FullTaggedTest<Mem_, DT_, IT_>("CheckpointTupleDiagMatrixTest")
  {
  }

  virtual ~CheckpointTupleDiagMatrixTest()
  {
  }


  virtual void run() const override
  {
    LAFEM::TupleDiagMatrix<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>, LAFEM::SparseMatrixCOO<Mem_, DT_, IT_>> tuple1;
    LAFEM::PointstarFactoryFD<DT_, IT_> five_star(40,2);
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star1 = five_star.matrix_csr();
    LAFEM::SparseMatrixCOO<Mem_, DT_, IT_> star2(star1);

    tuple1.template at <0,0>() = star1.clone();
    tuple1.template at <1,1>() = star2.clone();
    auto comm = Dist::Comm::world();
    Control::CheckpointControl cp(comm);
    BinaryStream bs;
    cp.add_object(String("tuple"), tuple1);
    cp.save(bs);
    comm.barrier();
    bs.seekg(0);
    cp.load(bs);
    LAFEM::TupleDiagMatrix<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>, LAFEM::SparseMatrixCOO<Mem_, DT_, IT_>> tuple2;
    cp.restore_object(String("tuple"), tuple2, false);
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star_temp = tuple2.template at<0,0>().clone();
    LAFEM::SparseMatrixCOO<Mem_, DT_, IT_> star_temp2 = tuple2.template at<1,1>().clone();
    TEST_CHECK_EQUAL(star1, star_temp);
    TEST_CHECK_EQUAL(star2, star_temp2);
    //Question: why does TEST_CHECK_EQUAL(tuple2.template at<0,0>(), tuple1.template at<0,0>()) not work?
    TEST_CHECK_NOT_EQUAL(cp.get_identifier_list().find("tuple"), std::string::npos);
  }
};
CheckpointTupleDiagMatrixTest<Mem::Main, float, unsigned int> checkpoint_tuple_diag_matrix_test_float_uint;
CheckpointTupleDiagMatrixTest<Mem::Main, float, unsigned long> checkpoint_tuple_diag_matrix_test_float_ulong;
CheckpointTupleDiagMatrixTest<Mem::Main, double, unsigned int> checkpoint_tuple_diag_matrix_test_double_uint;
CheckpointTupleDiagMatrixTest<Mem::Main, double, unsigned long> checkpoint_tuple_diag_matrix_test_double_ulong;

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class CheckpointTupleMatrixTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  CheckpointTupleMatrixTest()
    : FullTaggedTest<Mem_, DT_, IT_>("CheckpointTupleMatrixTest")
  {
  }

  virtual ~CheckpointTupleMatrixTest()
  {
  }


  virtual void run() const override
  {
    typedef LAFEM::TupleMatrixRow<LAFEM::SparseMatrixCSR<Mem_, DT_, IT_>,LAFEM::SparseMatrixCOO<Mem_, DT_, IT_>> tuprow;
    LAFEM::TupleMatrix<tuprow, tuprow> tuple1;
    LAFEM::PointstarFactoryFD<DT_, IT_> five_star(40,2);
    LAFEM::PointstarFactoryFD<DT_, IT_> nine_star(40,3);
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star1 = five_star.matrix_csr();
    LAFEM::SparseMatrixCOO<Mem_, DT_, IT_> star2(star1);
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star3 = five_star.matrix_csr();
    LAFEM::SparseMatrixCOO<Mem_, DT_, IT_> star4(star1);

    tuple1.template at <0,0>() = star1.clone();
    tuple1.template at <0,1>() = star2.clone();
    tuple1.template at <1,0>() = star3.clone();
    tuple1.template at <1,1>() = star4.clone();
    auto comm = Dist::Comm::world();
    Control::CheckpointControl cp(comm);
    BinaryStream bs;
    cp.add_object(String("tuple"), tuple1);
    cp.save(bs);
    comm.barrier();
    bs.seekg(0);
    cp.load(bs);
    LAFEM::TupleMatrix<tuprow,tuprow> tuple2;
    cp.restore_object(String("tuple"), tuple2, false);
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star_temp = tuple2.template at<0,0>().clone();
    LAFEM::SparseMatrixCOO<Mem_, DT_, IT_> star_temp2 = tuple2.template at<0,1>().clone();
    LAFEM::SparseMatrixCSR<Mem_, DT_, IT_> star_temp3 = tuple2.template at<1,0>().clone();
    LAFEM::SparseMatrixCOO<Mem_, DT_, IT_> star_temp4 = tuple2.template at<1,1>().clone();
    TEST_CHECK_EQUAL(star1, star_temp);
    TEST_CHECK_EQUAL(star2, star_temp2);
    TEST_CHECK_EQUAL(star3, star_temp3);
    TEST_CHECK_EQUAL(star4, star_temp4);
    TEST_CHECK_NOT_EQUAL(cp.get_identifier_list().find("tuple"), std::string::npos);
  }
};
CheckpointTupleMatrixTest<Mem::Main, float, unsigned int> checkpoint_tuple_matrix_test_float_uint;
CheckpointTupleMatrixTest<Mem::Main, float, unsigned long> checkpoint_tuple_matrix_test_float_ulong;
CheckpointTupleMatrixTest<Mem::Main, double, unsigned int> checkpoint_tuple_matrix_test_double_uint;
CheckpointTupleMatrixTest<Mem::Main, double, unsigned long> checkpoint_tuple_matrix_test_double_ulong;

template<
  typename Mem_,
  typename DT_,
  typename IT_>
class CheckpointTupleVectorTest
  : public FullTaggedTest<Mem_, DT_, IT_>
{
public:
  CheckpointTupleVectorTest()
    : FullTaggedTest<Mem_, DT_, IT_>("CheckpointTupleVectorTest")
  {
  }

  virtual ~CheckpointTupleVectorTest()
  {
  }


  virtual void run() const override
  {
    //Does not function with sparse_vector, as Containertype is missing <-- should we add that?
    LAFEM::TupleVector<LAFEM::DenseVector<Mem_, DT_, IT_>, LAFEM::DenseVector<Mem_, DT_, IT_>> tuple1;
    LAFEM::DenseVector<Mem_, DT_, IT_> dv1(1234);
    for (Index i(0) ; i < dv1.size() ; ++i)
      dv1(i, DT_(i) / DT_(12));

    LAFEM::DenseVector<Mem_, DT_, IT_> a(329);
    for (Index i(0) ; i < a.size() ; ++i)
      a(i, DT_(i) / DT_(7));

    tuple1.template at<0>() = dv1.clone();
    tuple1.template at<1>() = a.clone();
    auto comm = Dist::Comm::world();
    Control::CheckpointControl cp(comm);
    BinaryStream bs;
    cp.add_object(String("tuple"), tuple1);
    cp.save(bs);
    comm.barrier();
    bs.seekg(0);
    cp.load(bs);
    LAFEM::TupleVector<LAFEM::DenseVector<Mem_, DT_, IT_>, LAFEM::DenseVector<Mem_, DT_, IT_>> tuple2;
    cp.restore_object(String("tuple"), tuple2, false);
    LAFEM::DenseVector<Mem_, DT_, IT_> temp1 = tuple2.template at<0>().clone();
    LAFEM::DenseVector<Mem_, DT_, IT_> temp2 = tuple2.template at<1>().clone();
    TEST_CHECK_EQUAL(dv1, temp1);
    TEST_CHECK_EQUAL(a, temp2);
    TEST_CHECK_NOT_EQUAL(cp.get_identifier_list().find("tuple"), std::string::npos);

  }
};

CheckpointTupleVectorTest<Mem::Main, float, unsigned int> checkpoint_tuple_vector_test_float_uint;
CheckpointTupleVectorTest<Mem::Main, float, unsigned long> checkpoint_tuple_vector_test_float_ulong;
CheckpointTupleVectorTest<Mem::Main, double, unsigned int> checkpoint_tuple_vector_test_double_uint;
CheckpointTupleVectorTest<Mem::Main, double, unsigned long> checkpoint_tuple_vector_test_double_ulong;
