// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <test_system/test_system.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/cro_rav_ran_tur/element.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Ranancher-Turek Element test
 *
 * \test Tests the Ranancher-Turek Finite-Element space
 *
 * \tparam DataType_
 * The data type for the test. Shall be either double or float.
 *
 * \author Peter Zajac
 */
template<typename DataType_, typename IndexType_>
class RannacherTurekTest
  : public UnitTest
{
  typedef Shape::Quadrilateral ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::CroRavRanTur::Element<QuadTrafo> QuadSpaceStdNonPar;

  typedef Cubature::Rule<ShapeType, DataType_, DataType_, Tiny::Vector<DataType_, 2> > CubatureRule;

  static constexpr TrafoTags unit_trafo_config = TrafoTags::dom_point | TrafoTags::img_point | TrafoTags::jac_mat | TrafoTags::jac_det;

  static constexpr SpaceTags unit_space_config = SpaceTags::value | SpaceTags::grad;

public:
  RannacherTurekTest(PreferredBackend backend) :
    UnitTest("Rannacher-Turek Test", Type::Traits<DataType_>::name(), Type::Traits<IndexType_>::name(), backend)
  {
  }

  virtual ~RannacherTurekTest()
  {
  }

  virtual void run() const override
  {
    // test assembly on unit quad
    asm_unit_quad_std_non_par();
  }

  void asm_unit_quad_std_non_par() const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create a quad mesh
    Geometry::UnitCubeFactory<QuadMesh> mesh_factory;
    QuadMesh mesh(mesh_factory);

    // create a quad-trafo
    QuadTrafo trafo(mesh);

    // create a Q1~ space, standard non-parametric
    QuadSpaceStdNonPar space(trafo);

    // create a trafo evaluator
    typedef typename QuadTrafo::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);
    typename TrafoEvaluator::template ConfigTraits<unit_trafo_config>::EvalDataType trafo_data;

    // create a space evaluator
    typedef typename QuadSpaceStdNonPar::template Evaluator<TrafoEvaluator>::Type SpaceEvaluator;
    SpaceEvaluator space_eval(space);
    typename SpaceEvaluator::template ConfigTraits<unit_space_config>::EvalDataType space_data;

    // create a 3x3 Gauss-Legendre cubature formula
    CubatureRule cubature_rule(Cubature::ctor_factory, Cubature::DynamicFactory("gauss-legendre:3"));

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    // prepare space evaluator
    space_eval.prepare(trafo_eval);

    // check the number of local DOFs
    int num_loc_dofs = space_eval.get_num_local_dofs();
    TEST_CHECK_EQUAL(num_loc_dofs, 4u);

    // create local matrix assembly data
    Tiny::Matrix<DataType_, 4, 4> L, M;
    L = DataType_(0);
    M = DataType_(0);

    // loop over all 4 quadrature points and integrate
    for(int k(0); k < cubature_rule.get_num_points(); ++k)
    {
      // compute trafo data
      trafo_eval(trafo_data, cubature_rule.get_point(k));

      // compute space data
      space_eval(space_data, trafo_data);

      // test function loop
      for(int i(0); i < num_loc_dofs; ++i)
      {
        // trial function loop
        for(int j(0); j < num_loc_dofs; ++j)
        {
          // mass matrix entry
          M(i,j) += trafo_data.jac_det * cubature_rule.get_weight(k) * (
            space_data.phi[i].value * space_data.phi[j].value);

          // laplace matrix entry
          L(i,j) += trafo_data.jac_det * cubature_rule.get_weight(k) * (
            space_data.phi[i].grad[0] * space_data.phi[j].grad[0] +
            space_data.phi[i].grad[1] * space_data.phi[j].grad[1]);
          // continue with next trial function
        }
        // continue with next test function
      }
      // continue with next cubature point
    }

    // finish evaluators
    space_eval.finish();
    trafo_eval.finish();

    // test function loop
    for(int i(0); i < num_loc_dofs; ++i)
    {
      // trial function loop
      for(int j(0); j < num_loc_dofs; ++j)
      {
        // check entries
        if(i == j)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(M(i,j), DataType_(41) / DataType_(240), eps); // = 1.708p3e-1
          TEST_CHECK_EQUAL_WITHIN_EPS(L(i,j), DataType_(5) / DataType_(2), eps); // = 2.5
        }
        else if((int(i >> 1)) == int(j >> 1)) // i-j-pairs: (0-1) and (2-3)
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(M(i,j), DataType_(5) / DataType_(1200), eps); // = 4.1p6e-3
          TEST_CHECK_EQUAL_WITHIN_EPS(L(i,j), DataType_(1) / DataType_(2), eps); // = 0.5
        }
        else
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(M(i,j), DataType_(3) / DataType_(80), eps); // = 3.75e-2
          TEST_CHECK_EQUAL_WITHIN_EPS(L(i,j), DataType_(-3) / DataType_(2), eps); // = -1.5
        }
        // continue with next trial function
      }
      // continue with next test function
    }
  }
};

RannacherTurekTest <double, std::uint64_t> rannacher_turek_test_double_uint64(PreferredBackend::generic);
RannacherTurekTest <float, std::uint64_t> rannacher_turek_test_float_uint64(PreferredBackend::generic);
#ifdef FEAT_HAVE_MKL
RannacherTurekTest <float, std::uint64_t> mkl_rannacher_turek_test_float_uint64(PreferredBackend::mkl);
RannacherTurekTest <double, std::uint64_t> mkl_rannacher_turek_test_double_uint64(PreferredBackend::mkl);
#endif
#ifdef FEAT_HAVE_QUADMATH
RannacherTurekTest <__float128, std::uint64_t> rannacher_turek_test_float128_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_HALFMATH
RannacherTurekTest <Half, std::uint64_t> rannacher_turek_test_half_uint64(PreferredBackend::generic);
#endif
#ifdef FEAT_HAVE_CUDA
RannacherTurekTest <float, std::uint64_t> cuda_rannacher_turek_test_float_uint64(PreferredBackend::cuda);
RannacherTurekTest <double, std::uint64_t> cuda_rannacher_turek_test_double_uint64(PreferredBackend::cuda);
#endif
