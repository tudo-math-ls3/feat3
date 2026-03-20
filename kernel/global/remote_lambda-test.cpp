// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/global/remote_lambda.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/string.hpp>
#include <test_system/test_system.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;
using namespace FEAT::Global;

class SelfRemoteLambdaTest : public UnitTest
{
public:
  SelfRemoteLambdaTest() : UnitTest("SelfRemoteLambdaTest")
  {
  }

  void run() const override
  {
    // Check empty call
    {
      ScalarRemoteLambda<int(int)> rl(Dist::Comm::world());

      std::vector<int> ranks{};
      std::vector<Index> indices{};
      std::vector<std::tuple<int>> args{};

      std::vector<int> results = rl.call(ranks, indices, args, [](int i) { return i; });

      TEST_CHECK(results.empty());
    }

    // Check self call
    {
      ScalarRemoteLambda<int(int)> rl(Dist::Comm::world());

      const int self = Dist::Comm::world().rank();
      std::vector<int> ranks{self, self, self, self};
      std::vector<Index> indices{0, 0, 0, 0};
      std::vector<std::tuple<int>> args{{1}};

      std::vector<int> results = rl.call(ranks, indices, args, [](int i) { return i; });

      TEST_CHECK_EQUAL(results.size(), 4);
      TEST_CHECK_EQUAL(results[0], 1);
      TEST_CHECK_EQUAL(results[1], 1);
      TEST_CHECK_EQUAL(results[2], 1);
      TEST_CHECK_EQUAL(results[3], 1);
    }
  }
};

class ScalarRemoteLambdaTest : public UnitTest
{
public:
  ScalarRemoteLambdaTest() : UnitTest("ScalarRemoteLambdaTest")
  {
  }

  void run() const override
  {
    // Tests are designed for 4 ranks
    if(Dist::Comm::world().size() < 4)
    {
      return;
    }

    // Check empty call
    {
      ScalarRemoteLambda<int(int)> rl(Dist::Comm::world());

      std::vector<int> ranks{};
      std::vector<Index> indices{};
      std::vector<std::tuple<int>> args{};

      std::vector<int> results = rl.call(ranks, indices, args, [](int i) { return i; });

      TEST_CHECK(results.empty());
    }

    // Check self call
    {
      ScalarRemoteLambda<int(int)> rl(Dist::Comm::world());

      const int self = Dist::Comm::world().rank();
      std::vector<int> ranks{self, self, self, self};
      std::vector<Index> indices{0, 0, 0, 0};
      std::vector<std::tuple<int>> args{{1}};

      std::vector<int> results = rl.call(ranks, indices, args, [](int i) { return i; });

      TEST_CHECK_EQUAL(results.size(), 4);
      TEST_CHECK_EQUAL(results[0], 1);
      TEST_CHECK_EQUAL(results[1], 1);
      TEST_CHECK_EQUAL(results[2], 1);
      TEST_CHECK_EQUAL(results[3], 1);
    }

    // Check scalar single call per rank
    {
      ScalarRemoteLambda<int(int)> rl(Dist::Comm::world());

      std::vector<int> ranks{0, 1, 2, 3};
      std::vector<Index> indices{0, 0, 1, 1};
      std::vector<std::tuple<int>> args{0, 1};

      std::vector<int> results = rl.call(ranks, indices, args, [](int i) { return i; });

      TEST_CHECK_EQUAL(results.size(), ranks.size());

      TEST_CHECK_EQUAL(results[0], 0);
      TEST_CHECK_EQUAL(results[1], 0);
      TEST_CHECK_EQUAL(results[2], 1);
      TEST_CHECK_EQUAL(results[3], 1);
    }

    // Check scalar single call per rank, in inverse order to find potential message ordering issues
    {
      ScalarRemoteLambda<int(int)> rl(Dist::Comm::world());

      std::vector<int> ranks{3, 2, 1, 0};
      std::vector<Index> indices{0, 0, 1, 1};
      std::vector<std::tuple<int>> args{0, 1};

      std::vector<int> results = rl.call(ranks, indices, args, [](int i) { return i; });

      TEST_CHECK_EQUAL(results.size(), ranks.size());

      TEST_CHECK_EQUAL(results[0], 0);
      TEST_CHECK_EQUAL(results[1], 0);
      TEST_CHECK_EQUAL(results[2], 1);
      TEST_CHECK_EQUAL(results[3], 1);
    }

    // Check scalar single call per rank with non-identity lambda and reference arguments
    {
      ScalarRemoteLambda<int(const double&, int)> rl(Dist::Comm::world());

      const double data = 10.0;
      std::vector<int> ranks{0, 1, 2, 3};
      std::vector<Index> indices{0, 0, 1, 1};
      std::vector<std::tuple<const double&, int>> args{
        {data, 5},
        {data, 10},
      };

      std::vector<int> results = rl.call(ranks, indices, args, [](const double& d, int i) { return i + int(d); });

      TEST_CHECK_EQUAL(results.size(), ranks.size());
      TEST_CHECK_EQUAL(results[0], 15);
      TEST_CHECK_EQUAL(results[1], 15);
      TEST_CHECK_EQUAL(results[2], 20);
      TEST_CHECK_EQUAL(results[3], 20);
    }

    // Check scalar multi call per rank with non-identity lambda
    {
      ScalarRemoteLambda<int(const double&, int)> rl(Dist::Comm::world());

      const double data = 10.0;
      std::vector<int> ranks{0, 1, 2, 3, 0, 1, 2, 3};
      std::vector<Index> indices{0, 0, 1, 1, 0, 0, 1, 1};
      std::vector<std::tuple<const double&, int>> args{
        {data, 5},
        {data, 10},
      };

      std::vector<int> results = rl.call(ranks, indices, args, [](const double& d, int i) { return i + int(d); });

      TEST_CHECK_EQUAL(results.size(), ranks.size());
      TEST_CHECK_EQUAL(results[0], 15);
      TEST_CHECK_EQUAL(results[1], 15);
      TEST_CHECK_EQUAL(results[2], 20);
      TEST_CHECK_EQUAL(results[3], 20);
      TEST_CHECK_EQUAL(results[4], 15);
      TEST_CHECK_EQUAL(results[5], 15);
      TEST_CHECK_EQUAL(results[6], 20);
      TEST_CHECK_EQUAL(results[7], 20);
    }
  }
};

class VectorRemoteLambdaTest : public UnitTest
{
public:
  VectorRemoteLambdaTest() : UnitTest("VectorRemoteLambdaTest")
  {
  }

  template<typename T>
  using Vector = std::vector<T>;

  void run() const override
  {
    // Tests are designed for 4 ranks
    if(Dist::Comm::world().size() < 4)
    {
      return;
    }

    // Check empty call
    {
      VectorRemoteLambda<std::vector<int>(int)> rl(Dist::Comm::world());

      std::vector<int> ranks{};
      std::vector<Index> indices{};
      std::vector<std::tuple<int>> args{};

      std::vector<std::vector<int>> results =
        rl.call(ranks, indices, args, [](int i) { return std::vector<int>(std::size_t(i), i); });

      TEST_CHECK(results.empty());
    }

    // Check self call
    {
      VectorRemoteLambda<std::vector<int>(int)> rl(Dist::Comm::world());

      const int self = Dist::Comm::world().rank();
      std::vector<int> ranks{self, self, self, self};
      std::vector<Index> indices{0, 0, 0, 0};
      std::vector<std::tuple<int>> args{{1}};

      std::vector<std::vector<int>> results =
        rl.call(ranks, indices, args, [](int i) { return std::vector(std::size_t(i), i); });

      std::vector<int> expected;
      expected.push_back(1);

      TEST_CHECK_EQUAL(results.size(), 4);
      TEST_CHECK(results[0] == expected);
      TEST_CHECK(results[1] == expected);
      TEST_CHECK(results[2] == expected);
      TEST_CHECK(results[3] == expected);
    }

    // Check empty results
    {
      VectorRemoteLambda<std::vector<int>(int)> rl(Dist::Comm::world());

      std::vector<int> ranks{0, 1, 2, 3};
      std::vector<Index> indices{0, 0, 0, 0};
      std::vector<std::tuple<int>> args{{1}};

      std::vector<std::vector<int>> results =
        rl.call(ranks, indices, args, [](int /*i*/) { return std::vector<int>{}; });

      std::vector<int> expected;

      TEST_CHECK_EQUAL(results.size(), 4);
      TEST_CHECK(results[0] == expected);
      TEST_CHECK(results[1] == expected);
      TEST_CHECK(results[2] == expected);
      TEST_CHECK(results[3] == expected);
    }

    // Check single call per rank
    {
      VectorRemoteLambda<std::vector<int>(int)> rl(Dist::Comm::world());

      std::vector<int> ranks{0, 1, 2, 3};
      std::vector<Index> indices{0, 1, 2, 3};
      std::vector<std::tuple<int>> args{0, 1, 2, 3};

      std::vector<std::vector<int>> results = rl.call(
        ranks,
        indices,
        args,
        [](int i) { return std::vector<int>(std::size_t(Dist::Comm::world().rank()), i); });

      std::vector<int> a({});
      std::vector<int> b({1});
      std::vector<int> c({2, 2});
      std::vector<int> d({3, 3, 3});

      TEST_CHECK_EQUAL(results.size(), ranks.size());
      TEST_CHECK(results[0] == a);
      TEST_CHECK(results[1] == b);
      TEST_CHECK(results[2] == c);
      TEST_CHECK(results[3] == d);
    }

    // Check single call per rank, in inverse order to find potential message ordering issues
    {
      VectorRemoteLambda<std::vector<int>(int)> rl(Dist::Comm::world());

      std::vector<int> ranks{3, 2, 1, 0};
      std::vector<Index> indices{0, 1, 2, 3};
      std::vector<std::tuple<int>> args{0, 1, 2, 3};

      std::vector<std::vector<int>> results = rl.call(
        ranks,
        indices,
        args,
        [](int i) { return std::vector<int>(std::size_t(Dist::Comm::world().rank()), i); });

      std::vector<int> a({0, 0, 0});
      std::vector<int> b({1, 1});
      std::vector<int> c({2});
      std::vector<int> d({});

      TEST_CHECK_EQUAL(results.size(), ranks.size());
      TEST_CHECK(results[0] == a);
      TEST_CHECK(results[1] == b);
      TEST_CHECK(results[2] == c);
      TEST_CHECK(results[3] == d);
    }

    // Check multi call per rank
    {
      VectorRemoteLambda<std::vector<int>(int)> rl(Dist::Comm::world());

      std::vector<int> ranks{0, 1, 2, 3, 0, 1, 2, 3};
      std::vector<Index> indices{0, 1, 2, 3, 0, 1, 2, 3};
      std::vector<std::tuple<int>> args{0, 1, 2, 3};

      std::vector<std::vector<int>> results = rl.call(
        ranks,
        indices,
        args,
        [](int i) { return std::vector<int>(std::size_t(Dist::Comm::world().rank()), i); });

      std::vector<int> a({});
      std::vector<int> b({1});
      std::vector<int> c({2, 2});
      std::vector<int> d({3, 3, 3});

      TEST_CHECK_EQUAL(results.size(), ranks.size());
      TEST_CHECK(results[0] == a);
      TEST_CHECK(results[1] == b);
      TEST_CHECK(results[2] == c);
      TEST_CHECK(results[3] == d);
      TEST_CHECK(results[4] == a);
      TEST_CHECK(results[5] == b);
      TEST_CHECK(results[6] == c);
      TEST_CHECK(results[7] == d);
    }
  }
};

static const SelfRemoteLambdaTest self_remote_lambda_test;
static const ScalarRemoteLambdaTest scalar_remote_lambda_test;
static const VectorRemoteLambdaTest vector_remote_lambda_test;
