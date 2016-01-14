#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/scarc/scarc_data.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::TestSystem;
using namespace FEAST::LAFEM;
using namespace FEAST::ScaRC;


template<typename Tag_, typename DataType_>
class ScaRCDataTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    ScaRCDataTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("ScaRCDataTest<" + tag + ">")
    {
    }

    virtual void run() const override
    {
      DenseVector<Tag_, DataType_> x(1000, DataType_(1)), b(1000, DataType_(1));
      SparseMatrixCOO<Mem::Main, DataType_> T(1000, 1000);
      for(Index i(0) ; i < 1000 ; ++i)
        T(i, i, DataType_(1));
      SparseMatrixCSR<Tag_, DataType_> A(T);
      SparseMatrixCSR<Tag_, DataType_> P;
      P.clone(A);
      ScaRCData<> sd0(std::move(A), std::move(x), std::move(b));

      DenseVector<Tag_, DataType_> x1(1000, DataType_(1)), b1(1000, DataType_(1));
      SparseMatrixCOO<Mem::Main, DataType_> T1(1000, 1000);
      for(Index i(0) ; i < 1000 ; ++i)
        T1(i, i, DataType_(1));
      SparseMatrixCSR<Tag_, DataType_> A1(T1);
      SparseMatrixCSR<Tag_, DataType_> P1;
      P1.clone(A1);
      ScaRCData<> sd1(std::move(A1), std::move(x1), std::move(b1));

      DenseVector<Tag_, DataType_> x2(1000, DataType_(1)), b2(1000, DataType_(1));
      SparseMatrixCOO<Mem::Main, DataType_> T2(1000, 1000);
      for(Index i(0) ; i < 1000 ; ++i)
        T2(i, i, DataType_(1));
      SparseMatrixCSR<Tag_, DataType_> A2(T2);
      SparseMatrixCSR<Tag_, DataType_> P2;
      P2.clone(A2);
      PreconditionedScaRCData<> psd(std::move(A2), std::move(P2), std::move(x2), std::move(b2));

      DenseVector<Tag_, DataType_> x3(1000, DataType_(1)), b3(1000, DataType_(1));
      SparseMatrixCOO<Mem::Main, DataType_> T3(1000, 1000);
      for(Index i(0) ; i < 1000 ; ++i)
        T3(i, i, DataType_(1));
      SparseMatrixCSR<Tag_, DataType_> A3(T3);
      SparseMatrixCSR<Tag_, DataType_> P3;
      P3.clone(A3);
      SynchronisedScaRCData<> sd3(std::move(A3), std::move(x3), std::move(b3));
      TEST_CHECK_EQUAL(sd3.vector_mirrors().size(), Index(0));
      TEST_CHECK_EQUAL(sd3.vector_mirror_sendbufs().size(), Index(0));
      TEST_CHECK_EQUAL(sd3.vector_mirror_recvbufs().size(), Index(0));
      TEST_CHECK_EQUAL(sd3.dest_ranks().size(), Index(0));
      TEST_CHECK_EQUAL(sd3.source_ranks().size(), Index(0));

      DenseVector<Tag_, DataType_> x4(1000, DataType_(1)), b4(1000, DataType_(1));
      SparseMatrixCOO<Mem::Main, DataType_> T4(1000, 1000);
      for(Index i(0) ; i < 1000 ; ++i)
        T4(i, i, DataType_(1));
      SparseMatrixCSR<Tag_, DataType_> A4(T4);
      SparseMatrixCSR<Tag_, DataType_> P4;
      P4.clone(A4);
      SynchronisedPreconditionedScaRCData<> sd4(std::move(A4), std::move(P4), std::move(x4), std::move(b4));
      TEST_CHECK_EQUAL(sd4.vector_mirrors().size(), Index(0));
      TEST_CHECK_EQUAL(sd4.vector_mirror_sendbufs().size(), Index(0));
      TEST_CHECK_EQUAL(sd4.vector_mirror_recvbufs().size(), Index(0));
      TEST_CHECK_EQUAL(sd4.dest_ranks().size(), Index(0));
      TEST_CHECK_EQUAL(sd4.source_ranks().size(), Index(0));

      //-------------------------------------------------------------------------

      DenseVector<Tag_, DataType_> x5(1000, DataType_(1)), b5(1000, DataType_(1));
      SparseMatrixCOO<Mem::Main, DataType_> T5(1000, 1000);
      for(Index i(0) ; i < 1000 ; ++i)
        T4(i, i, DataType_(1));
      SparseMatrixCSR<Tag_, DataType_> A5(T5);
      SparseMatrixCSR<Tag_, DataType_> P5;
      P5.clone(A5);
      MultiLevelScaRCData<> mlsd(std::move(A5), std::move(b5), std::move(x5));
    }
};
ScaRCDataTest<Mem::Main, double> sf_cpu_double("StorageType: std::vector, DataType: double");
