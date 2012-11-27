#include <kernel/base_header.hpp>
#include <test_system/test_system.hpp>

#include <kernel/scarc/solver_functor.hpp>
#include <kernel/util/cpp11_smart_pointer.hpp>
#include <kernel/archs.hpp>
#include<deque>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::ScaRC;
using namespace FEAST::TestSystem;


template<typename Tag_, typename Algo_, typename DataType_>
class SolverFunctorTest:
  public TaggedTest<Tag_, DataType_>
{
  public:
    SolverFunctorTest(const std::string & tag) :
      TaggedTest<Tag_, DataType_>("SolverFunctorTest<" + tag + ">")
    {
    }

    virtual void run() const
    {
      DenseVector<Tag_, DataType_> y(1000, DataType_(0)), r(1000, DataType_(1));
      DenseVector<Tag_, DataType_> dummy;
      DenseVector<Tag_, DataType_> l(1000, DataType_(1));

      DenseVector<Tag_, DataType_> y_ref(1000, DataType_(0));

      SumFunctorProxyLeft<Algo_, DenseVector<Tag_, DataType_> > sf(y, dummy, r);
      TEST_CHECK_THROWS(sf.execute(), ScaRCError);

      sf.substitute(l);
      sf.execute();

      Sum<Algo_>::value(y_ref, l, r);
      TEST_CHECK_EQUAL(y, y_ref);

      //---------------------------------------------------------------------------------------------

      SumFunctorProxyResultLeft<Algo_, DenseVector<Tag_, DataType_> > sf1(dummy, dummy, r);
      sf1.substitute(l);
      sf1.execute();

      TEST_CHECK_EQUAL(l, y_ref);

      //---------------------------------------------------------------------------------------------

      SumFunctor<Algo_, DenseVector<Tag_, DataType_> > sf2(y, r, r);
      sf2.execute();

      TEST_CHECK_EQUAL(y, y_ref);

      //---------------------------------------------------------------------------------------------

      //testing substitution in CompoundSolverFunctor
      //a <- a + A*a := [[b <- A*a , a <- a + b]]
      DenseVector<Tag_, DataType_> a(1000, DataType_(1)), b(1000, DataType_(1));
      SparseMatrixCOO<Mem::Main, DataType_> T(1000, 1000);
      for(Index i(0) ; i < 1000 ; ++i)
        T(i, i, DataType_(1));
      SparseMatrixCSR<Tag_, DataType_> A(T);

      CompoundSolverFunctor<Algo_, DenseVector<Tag_, DataType_> > cf;
      cf.add_functor(new ProductFunctorProxyRight<Algo_, DenseVector<Tag_, DataType_>, SparseMatrixCSR<Tag_, DataType_> >(b, A, dummy));
      cf.add_functor(new SumFunctorProxyResultLeft<Algo_, DenseVector<Tag_, DataType_> >(dummy, dummy, b));
      cf.substitute(a);
      cf.execute();

      //reference
      DenseVector<Tag_, DataType_> a_ref(1000, DataType_(1));
      DenseVector<Tag_, DataType_> b_ref(1000, DataType_(1));
      Product<Algo_>::value(b_ref, A, a_ref);
      Sum<Algo_>::value(a_ref, b_ref, a_ref);

      TEST_CHECK_EQUAL(a, a_ref);


      //---------------------------------------------------------------------------------------------

      //testing preconditioner functor
      // u <- u + precon(u), precon(.) := A*(.) := [[v <- u, precon(u), u <- u + v]]
      DenseVector<Tag_, DataType_> u(1000, DataType_(47.11));
      DenseVector<Tag_, DataType_> v(1000, DataType_(0));
      CompoundSolverFunctor<Algo_, DenseVector<Tag_, DataType_> > cf1;

      cf1.add_functor(new CopyFunctor<Algo_, DenseVector<Tag_, DataType_> >(v, u));
      cf1.add_functor(new PreconFunctor<Algo_, DenseVector<Tag_, DataType_> >(u));
      cf1.add_functor(new SumFunctor<Algo_, DenseVector<Tag_, DataType_> >(u, u, v));

      std::shared_ptr<SolverFunctorBase<DenseVector<Tag_, DataType_> > > p(new CompoundSolverFunctor<Algo_, DenseVector<Tag_, DataType_> >() );
      ((CompoundSolverFunctor<Algo_, DenseVector<Tag_, DataType_> >* )(p.get()))->add_functor(new ProductFunctorProxyResultRight<Algo_, DenseVector<Tag_, DataType_> , SparseMatrixCSR<Tag_, DataType_> >(dummy, A, dummy) );

      cf1.set_preconditioner(p);
      cf1.execute();

      DenseVector<Tag_, DataType_> u_ref(1000, DataType_(47.11));
      DenseVector<Tag_, DataType_> v_ref(1000, DataType_(0));

      copy(v_ref, u_ref);
      Product<Algo_>::value(u_ref, A, u_ref);
      Sum<Algo_>::value(u_ref, v_ref, u_ref);

      TEST_CHECK_EQUAL(u, u_ref);

      //---------------------------------------------------------------------------------------------

      //testing iterate functor
      double dummy_scalar(0);
      Index used_iters(0);

      DenseVector<Tag_, DataType_> u1(1000, DataType_(1));

      std::shared_ptr<SolverFunctorBase<DenseVector<Tag_, DataType_> > > innerfunc(new SumFunctor<Algo_, DenseVector<Tag_, DataType_> >(u1, u1, u1));
      IterateFunctor<Algo_, DenseVector<Tag_, DataType_>, double> iterfunc(innerfunc, dummy_scalar, dummy_scalar, used_iters, Index(3));

      iterfunc.execute();

      TEST_CHECK_EQUAL(used_iters, 3);

      //---------------------------------------------------------------------------------------------

      //testing defect functors

      DenseVector<Tag_, DataType_> vb(1000, DataType_(1)), vy(1000), vxdummy;
      DenseVector<Tag_, DataType_> vy_ref(1000);

      DefectFunctorProxyRight<Algo_, DenseVector<Tag_, DataType_>, SparseMatrixCSR<Tag_, DataType_> > defect(vy, vb, A, vxdummy);
      DenseVector<Tag_, DataType_> vx(1000, DataType_(1));

      TEST_CHECK_THROWS(defect.execute(), ScaRCError);

      defect.substitute(vx);
      defect.execute();

      Defect<Algo_>::value(vy_ref, vb, A, vx);
      TEST_CHECK_EQUAL(vy, vy_ref);

      DefectFunctor<Algo_, DenseVector<Tag_, DataType_>, SparseMatrixCSR<Tag_, DataType_> > defect1(vy, vb, A, vx);
      TEST_CHECK_EQUAL(vy, vy_ref);

      //---------------------------------------------------------------------------------------------
    }
};
SolverFunctorTest<Mem::Main, Algo::Generic, double> sf_cpu_double("StorageType: std::vector, DataType: double");
