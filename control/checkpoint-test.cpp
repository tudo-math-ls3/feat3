// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/util/runtime.hpp>
#include <control/checkpoint_control.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/util/binary_stream.hpp>
#include <test_system/test_system.hpp>

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

    {
      auto comm = Dist::Comm::world();
      Control::CheckpointControl cp(comm);
      BinaryStream bs;
      cp.add_object(String("dv1"), dv1);
      cp.save(bs);
      comm.barrier();
      bs.seekg(0);
      cp.load(bs);
      LAFEM::DenseVector<Mem_, DT_, IT_> dv2;
      cp.restore_object(String("dv1"), dv2, false);
      TEST_CHECK_EQUAL(dv1, dv2);
      TEST_CHECK_NOT_EQUAL(cp.get_identifier_list().find("dv1"), std::string::npos);
    }
  }
};

CheckpointTest<Mem::Main, double, unsigned int> checkpoint_test_double_uint;
