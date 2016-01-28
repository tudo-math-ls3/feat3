#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#ifdef FEAST_HAVE_PARMETIS
#include <kernel/archs.hpp>
FEAST_DISABLE_WARNINGS
#include <parmetis.h>
FEAST_RESTORE_WARNINGS

using namespace FEAST;
using namespace FEAST::TestSystem;


template<typename Tag_= Mem::Main, typename IndexType_ = Index>
class HaveParmetisTest:
  public TaggedTest<Tag_, IndexType_>
{
  public:
    HaveParmetisTest(const std::string & tag) :
      TaggedTest<Tag_, IndexType_>("HaveParmetisTest<" + tag + ">")
    {
    }

    virtual void run() const override
    {
      idx_t* vtxdist = new idx_t[2];
      idx_t* part = new idx_t[1];
      idx_t* ndims = new idx_t(1);
      real_t* xyz = new real_t[1];

      vtxdist[0] = idx_t(0);
      vtxdist[1] = idx_t(1);
      part[0] = idx_t(0);
      ndims[0] = idx_t(1);
      xyz[0] = real_t(0);

      MPI_Comm comm = MPI_COMM_WORLD;
      int result(ParMETIS_V3_PartGeom(vtxdist, ndims, xyz, part, &comm));

      TEST_CHECK_EQUAL(result, METIS_OK);

      delete[] vtxdist;
      delete[] part;
      delete ndims;
      delete[] xyz;
    }
};
HaveParmetisTest<> have_parmetis_test("None, Index");
#endif
