#include <test_system/test_system.hpp>
#include <kernel/geometry/test_aux/tetris_quad.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/dof_mapping_renderer.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/adjacency/graph.hpp>
#include <kernel/geometry/conformal_factories.hpp>
#include <kernel/util/math.hpp>

using namespace FEAT;
using namespace FEAT::TestSystem;

/**
 * \brief Lagrange-1 Element test
 *
 * \test Tests the Lagrange-1 Finite-Element space
 *
 * \tparam DataType_
 * The data type for the test. Shall be either double or float.
 *
 * \author Peter Zajac
 */
template<typename DataType_>
class Lagrange1Test
  : public TestSystem::TaggedTest<Archs::None, DataType_>
{
  typedef Shape::Quadrilateral ShapeType;
  typedef Geometry::ConformalMesh<ShapeType> QuadMesh;

  typedef Trafo::Standard::Mapping<QuadMesh> QuadTrafo;

  typedef Space::Lagrange1::Element<QuadTrafo> QuadSpaceQ1;

  typedef typename QuadSpaceQ1::DofMappingType DofMapping;

  typedef Cubature::Rule<ShapeType, DataType_, DataType_, Tiny::Vector<DataType_, 2> > CubatureRule;

  static constexpr TrafoTags unit_trafo_config = TrafoTags::jac_det | TrafoTags::jac_inv;
  static constexpr SpaceTags unit_space_config = SpaceTags::value | SpaceTags::grad;

  static constexpr TrafoTags tetris_trafo_config = TrafoTags::jac_det | TrafoTags::jac_inv;
  static constexpr SpaceTags tetris_space_config = SpaceTags::grad;

public:
  Lagrange1Test() :
    TestSystem::TaggedTest<Archs::None, DataType_>("Lagrange-1 Test")
  {
  }

  virtual ~Lagrange1Test()
  {
  }

  virtual void run() const override
  {
    // test assembly on unit quad
    asm_unit_quad();

    // test assembly on tetris quad
    asm_tetris_quad();
  }

  void asm_unit_quad() const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create a quad mesh
    Geometry::UnitCubeFactory<QuadMesh> mesh_factory;
    QuadMesh mesh(mesh_factory);

    // create a quad-trafo
    QuadTrafo trafo(mesh);

    // create a Q1 space
    QuadSpaceQ1 space(trafo);

    // create a trafo evaluator
    typedef typename QuadTrafo::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);
    typename TrafoEvaluator::template ConfigTraits<unit_trafo_config>::EvalDataType trafo_data;

    // create a space evaluator
    typedef typename QuadSpaceQ1::template Evaluator<TrafoEvaluator>::Type SpaceEvaluator;
    SpaceEvaluator space_eval(space);
    typename SpaceEvaluator::template ConfigTraits<unit_space_config>::EvalDataType space_data;

    // create a 2x2 Gauss-Legendre cubature formula
    CubatureRule cubature_rule(Cubature::ctor_factory, Cubature::DynamicFactory("gauss-legendre:2"));

    // prepare trafo evaluator
    trafo_eval.prepare(0);

    // prepare space evaluator
    space_eval.prepare(trafo_eval);

    // check the number of local DOFs
    int num_loc_dofs = space_eval.get_num_local_dofs();
    TEST_CHECK_EQUAL(num_loc_dofs, 4);

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
          TEST_CHECK_EQUAL_WITHIN_EPS(M(i,j), DataType_(1) / DataType_(9), eps);
          TEST_CHECK_EQUAL_WITHIN_EPS(L(i,j), DataType_(2) / DataType_(3), eps);
        }
        else if(3 - int(i) == int(j))
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(M(i,j), DataType_(1) / DataType_(36), eps);
          TEST_CHECK_EQUAL_WITHIN_EPS(L(i,j), DataType_(-1) / DataType_(3), eps);
        }
        else
        {
          TEST_CHECK_EQUAL_WITHIN_EPS(M(i,j), DataType_(1) / DataType_(18), eps);
          TEST_CHECK_EQUAL_WITHIN_EPS(L(i,j), DataType_(-1) / DataType_(6), eps);
        }
        // continue with next trial function
      }
      // continue with next test function
    }
  }

  void asm_tetris_quad() const
  {
    // compute eps
    const DataType_ eps = Math::pow(Math::eps<DataType_>(), DataType_(0.8));

    // create a quad mesh
    QuadMesh *mesh = Geometry::TestAux::create_tetris_mesh_2d();

    // create a quad-trafo
    QuadTrafo* trafo = new QuadTrafo(*mesh);

    // create a Q1 space
    QuadSpaceQ1* space = new QuadSpaceQ1(*trafo);

    // create matrix structure
    Index neq = 0, nnze = 0;
    Index* row_ptr = nullptr;
    Index* col_idx = nullptr;

    // assemble matrix structure
    asm_mat_struct(*space, neq, nnze, row_ptr, col_idx);

    // validate sizes
    TEST_CHECK_EQUAL(neq, 10u);
    TEST_CHECK_EQUAL(nnze, 52u);

    // validate row-pointer
    Index row_ptr_ref[11] =
    {
      0,
      4,
      10,
      14,
      18,
      26,
      34,
      38,
      42,
      48,
      52
    };
    for(Index i(0); i < 11; ++i)
    {
      TEST_CHECK_EQUAL(row_ptr[i], row_ptr_ref[i]);
    }

    // validate column-indices
    Index col_idx_ref[52] =
    {
      0, 1,       4, 5,
      0, 1, 2,    4, 5, 6,
         1, 2,       5, 6,
               3, 4,       7, 8,
      0, 1,    3, 4, 5,    7, 8, 9,
      0, 1, 2,    4, 5, 6,    8, 9,
         1, 2,       5, 6,
               3, 4,       7, 8,
               3, 4, 5,    7, 8, 9,
                  4, 5,       8, 9
    };
    for(Index i(0); i < 52; ++i)
    {
      TEST_CHECK_EQUAL(col_idx[i], col_idx_ref[i]);
    }

    // create data array
    DataType_* data = new DataType_[nnze];
    for(Index i(0); i < nnze; ++i)
    {
      data[i] = DataType_(0);
    }

    // assemble laplace matrix
    asm_laplace(*space, row_ptr, col_idx, data);

    // validate data array
    static const DataType_ q16 = -DataType_(1) / DataType_(6);
    static const DataType_ q13 = -DataType_(1) / DataType_(3);
    static const DataType_ d23 = DataType_(2) / DataType_(3);
    static const DataType_ d43 = DataType_(4) / DataType_(3);
    static const DataType_ d21 = DataType_(2);
    DataType_ data_ref[52] =
    {
      d23, q16,           q16, q13,
      q16, d43, q16,      q13, q13, q13,
           q16, d23,           q13, q16,
                     d23, q16,           q16, q13,
      q16, q13,      q16, d21, q13,      q13, q13, q13,
      q13, q13, q13,      q13, d21, q16,      q13, q16,
           q13, q16,           q16, d23,
                     q16, q13,           d23, q16,
                     q13, q13, q13,      q16, d43, q16,
                          q13, q16,           q16, d23
    };

    for(Index i(0); i < 52; ++i)
    {
      TEST_CHECK_EQUAL_WITHIN_EPS(data[i], data_ref[i], eps);
    }

    // delete matrix
    delete [] data;
    delete [] col_idx;
    delete [] row_ptr;

    // delete objects
    delete space;
    delete trafo;
    delete mesh;
  }

  void asm_mat_struct(const QuadSpaceQ1& space, Index& neq, Index& nnze, Index*& row_ptr, Index*& col_idx) const
  {
    // create a dof-mapper
    Adjacency::Graph dof_mapping(Space::DofMappingRenderer::render(space));

    // fetch number of dofs
    neq = space.get_num_dofs();

    // render transposed dof-mapper to graph
    Adjacency::Graph dof_support(Adjacency::rt_transpose, dof_mapping);

    // render composite dof-support/dof-mapper graph
    Adjacency::Graph dof_adjactor(Adjacency::rt_injectify, dof_support, dof_mapping);

    // fetch graph arrays
    const Index* dom_ptr = dof_adjactor.get_domain_ptr();
    const Index* img_idx = dof_adjactor.get_image_idx();

    // fetch number of adjacencies
    nnze = dom_ptr[neq];

    // allocate arrays
    row_ptr = new Index[neq+1];
    col_idx = new Index[nnze];

    // copy arrays
    for(Index i(0); i <= neq; ++i)
    {
      row_ptr[i] = dom_ptr[i];
    }
    for(Index i(0); i < nnze; ++i)
    {
      col_idx[i] = img_idx[i];
    }

    // apply linear insertion sort on col_idx
    for(Index i(0); i < neq; ++i)
    {
      for(Index j(row_ptr[i]); j < row_ptr[i+1]; ++j)
      {
        for(Index k(j+1); k < row_ptr[i+1]; ++k)
        {
          if(col_idx[j] > col_idx[k])
          {
            Index m = col_idx[j];
            col_idx[j] = col_idx[k];
            col_idx[k] = m;
          }
        }
      }
    }
  }

  void asm_laplace(const QuadSpaceQ1& space, Index* row_ptr, Index* col_idx, DataType_* data) const
  {
    const QuadTrafo& trafo = space.get_trafo();

    // create a trafo evaluator
    typedef typename QuadTrafo::template Evaluator<ShapeType, DataType_>::Type TrafoEvaluator;
    TrafoEvaluator trafo_eval(trafo);
    typename TrafoEvaluator::template ConfigTraits<tetris_trafo_config>::EvalDataType trafo_data;

    // create a space evaluator
    typedef typename QuadSpaceQ1::template Evaluator<TrafoEvaluator>::Type SpaceEvaluator;
    SpaceEvaluator space_eval(space);
    typename SpaceEvaluator::template ConfigTraits<tetris_space_config>::EvalDataType space_data;

    // create a dof-mapper
    DofMapping dof_mapping(space);

    // create a 2x2 Gauss-Legendre cubature formula
    CubatureRule cubature_rule(Cubature::ctor_factory, Cubature::DynamicFactory("gauss-legendre:2"));

    // allocate local matrix
    Tiny::Matrix<DataType_, SpaceEvaluator::max_local_dofs, SpaceEvaluator::max_local_dofs> Lx;

    // allocate column pointer array
    Index* col_ptr = new Index[space.get_num_dofs()];

    for(Index cell(0); cell < trafo_eval.get_num_cells(); ++cell)
    {
      // prepare trafo evaluator
      trafo_eval.prepare(cell);

      // prepare space evaluator
      space_eval.prepare(trafo_eval);

      // check the number of local DOFs
      int num_loc_dofs = space_eval.get_num_local_dofs();
      TEST_CHECK_EQUAL(num_loc_dofs, 4u);

      // clear local matrix
      Lx = DataType_(0);

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
            // laplace matrix entry
            Lx(i,j) += trafo_data.jac_det * cubature_rule.get_weight(k) * (
              space_data.phi[i].grad[0] * space_data.phi[j].grad[0] +
              space_data.phi[i].grad[1] * space_data.phi[j].grad[1]);
          }
          // continue with next trial function
        }
        // continue with next test function
      }

      // finish evaluators
      space_eval.finish();
      trafo_eval.finish();

      // initialise dof-mapper
      dof_mapping.prepare(cell);

      // test function loop
      for(int i(0); i < num_loc_dofs; ++i)
      {
        // test function contribution loop
        for(int ic(0); ic < dof_mapping.get_num_contribs(i); ++ic)
        {
          // fetch test function dof index
          Index idof = dof_mapping.get_index(i, ic);

          // fetch test function dof weight
          DataType_ iw = DataType_(dof_mapping.get_weight(i, ic));

          // build column pointer for this test function contribution
          for(Index k(row_ptr[idof]); k < row_ptr[idof+1]; ++k)
          {
            col_ptr[col_idx[k]] = k;
          }

          // trial function loop
          for(int j(0); j < num_loc_dofs; ++j)
          {
            // trial function contribution loop
            for(int jc(0); jc < dof_mapping.get_num_contribs(j); ++jc)
            {
              // fetch trial function dof index
              Index jdof = dof_mapping.get_index(j, jc);

              // fetch trial function dof weight
              DataType_ jw = DataType_(dof_mapping.get_weight(j, jc));

              // incorporate data into global matrix
              data[col_ptr[jdof]] += iw * jw * Lx(i,j);

              // next trial function contribution
            }
            // continue with next trial function
          }
          // next test function contribution
        }
        // continue with next test function
      }

      // finish dof-mapper
      dof_mapping.finish();

      // continue with next cell
    }

    // delete column pointer
    delete [] col_ptr;
  }
};

Lagrange1Test<double> lagrange1_test_double;
Lagrange1Test<float> lagrange1_test_float;
