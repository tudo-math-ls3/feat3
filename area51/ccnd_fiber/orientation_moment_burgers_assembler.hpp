#pragma once

#include <area51/ccnd_fiber/tensor_operations.hpp>
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/assembly/discrete_evaluator.hpp>
#include <type_traits>

namespace FEAT
{
  namespace Assembly
  {

    template<typename Space_, typename DT_, typename IT_, int dim_>
    struct FEVelo
    {
      typedef typename LAFEM::DenseVectorBlocked<Mem::Main, DT_, IT_, dim_> LocalVectorType_;
      const Space_& _space;
      const LocalVectorType_& _local_FE_coeff;

      explicit FEVelo(const Space_& space, const LocalVectorType_& local_FE_coeff) :
      _space(space),
      _local_FE_coeff(local_FE_coeff)
      {
      }
    };

    /**
   * \brief Calculates the posible refined cell arrays a given reference point is in.
   *
   * \Attention Uses 2 level ordering of refined cells. If due to some reason, the refinement does not behave according to this rule, the point does not have to be in the returned vector.
   *
   * \Attention The input point has to be in [-1,1]^d
   *
   * \Attention Also the vector cells should have a siginifanct large buffer, bigger 4 or 8, so its inner array does not have to be reallocated...
   */
  template<typename DataType_, int dim_>
  void calc_refinement_cell(std::vector<Index>& DOXY(cells), const Tiny::Vector<DataType_, dim_>& DOXY(point), const Index DOXY(cell_ind), Index DOXY(recursion)){}

  template<typename DataType_>
  void calc_refinement_cell(std::vector<Index>& cells, const Tiny::Vector<DataType_, 2>& point, const Index cell_ind, Index recursion)
  {
    if(recursion == Index(0))
    {
      //reached our recursion depth, so nothing to refine, so just add the cell index
      cells.push_back(cell_ind);
      return;
    }
    std::vector<DataType_> new_ref_x_cord;
    std::vector<DataType_> new_ref_y_cord;
    std::vector<Index> domain_index_x;
    std::vector<Index> domain_index_y;
    //first check in which cells the next point could be
    //first check the x cordinate
    //Check if its value is near 0
    if(Math::abs(point[0]) < Math::pow(Math::eps<DataType_>(), DataType_(0.9)))
    {
      //then we push both points... whereby we shift them directly on to the edges in x-direction
      new_ref_x_cord.push_back(DataType_(1));
      new_ref_x_cord.push_back(DataType_(-1));
      domain_index_x.push_back(Index(0));
      domain_index_x.push_back(Index(1));
    }
    else if(point[0] < DataType_(0))
    {
      //shift our point by 1/2 to the right and scale it to [-1,1]
      new_ref_x_cord.push_back(DataType_(2)*(point[0] + DataType_(0.5)));
      domain_index_x.push_back(Index(0));
    }
    else if(point[0] > DataType_(0))
    {
      //shift our point by 1/2 to the right and scale it to [-1,1]
      new_ref_x_cord.push_back(DataType_(2)*(point[0] - DataType_(0.5)));
      domain_index_x.push_back(Index(1));
    }

    //first check the y cordinate
    //Check if its value is near 0
    if(Math::abs(point[1]) < Math::pow(Math::eps<DataType_>(), DataType_(0.9)))
    {
      //then we push both points... whereby we shift them directly on to the edges in x-direction
      new_ref_y_cord.push_back(DataType_(1));
      new_ref_y_cord.push_back(DataType_(-1));
      domain_index_y.push_back(Index(0));
      domain_index_y.push_back(Index(1));
    }
    else if(point[1] < DataType_(0))
    {
      //shift our point by 1/2 to the right and scale it to [-1,1]
      new_ref_y_cord.push_back(DataType_(2)*(point[1] + DataType_(0.5)));
      domain_index_y.push_back(Index(0));
    }
    else if(point[1] > DataType_(0))
    {
      //shift our point by 1/2 to the right and scale it to [-1,1]
      new_ref_y_cord.push_back(DataType_(2)*(point[1] - DataType_(0.5)));
      domain_index_y.push_back(Index(1));
    }

    Tiny::Vector<DataType_, 2> new_point;
    //now iterate over the elements
    for(Index i(0); i < domain_index_x.size(); ++i)
    {
      for(Index j(0); j < domain_index_y.size(); ++j)
      {
        //Due to the 2lvl ordering, our refinded cell of cell i is one of 4i, 4i+1, 4i+2, 4i+3
        const Index refined_cell_index = Index(4) * cell_ind + Index(2) * domain_index_y[j] + domain_index_x[i];
        new_point[0] = new_ref_x_cord[i];
        new_point[1] = new_ref_y_cord[j];
        calc_refinement_cell(cells, new_point, refined_cell_index, recursion - Index(1));
      }
    }
  }

  template<typename DataType_>
  void calc_refinement_cell(std::vector<Index>& cells, const Tiny::Vector<DataType_, 3>& point, const Index cell_ind, Index recursion)
  {
    if(recursion == Index(0))
    {
      //reach our recursion depth, so nothing to refine, so just add the cell index
      cells.push_back(cell_ind);
      return;
    }
    std::vector<DataType_> new_ref_x_cord;
    std::vector<DataType_> new_ref_y_cord;
    std::vector<DataType_> new_ref_z_cord;
    std::vector<Index> domain_index_x;
    std::vector<Index> domain_index_y;
    std::vector<Index> domain_index_z;
    //first check in which cells the next point could be
    //first check the x cordinate
    //Check if its value is near 0
    if(Math::abs(point[0]) < Math::pow(Math::eps<DataType_>(), DataType_(0.9)))
    {
      //then we push both points... whereby we shift them directly on to the edges in x-direction
      new_ref_x_cord.push_back(DataType_(1));
      new_ref_x_cord.push_back(DataType_(-1));
      domain_index_x.push_back(Index(0));
      domain_index_x.push_back(Index(1));
    }
    else if(point[0] < DataType_(0))
    {
      //shift our point by 1/2 to the right and scale it to [-1,1]
      new_ref_x_cord.push_back(DataType_(2)*(point[0] + DataType_(0.5)));
      domain_index_x.push_back(Index(0));
    }
    else if(point[0] > DataType_(0))
    {
      //shift our point by 1/2 to the right and scale it to [-1,1]
      new_ref_x_cord.push_back(DataType_(2)*(point[0] - DataType_(0.5)));
      domain_index_x.push_back(Index(1));
    }

    //first check the y cordinate
    //Check if its value is near 0
    if(Math::abs(point[1]) < Math::pow(Math::eps<DataType_>(), DataType_(0.9)))
    {
      //then we push both points... whereby we shift them directly on to the edges in x-direction
      new_ref_y_cord.push_back(DataType_(1));
      new_ref_y_cord.push_back(DataType_(-1));
      domain_index_y.push_back(Index(0));
      domain_index_y.push_back(Index(1));
    }
    else if(point[1] < DataType_(0))
    {
      //shift our point by 1/2 to the right and scale it to [-1,1]
      new_ref_y_cord.push_back(DataType_(2)*(point[1] + DataType_(0.5)));
      domain_index_y.push_back(Index(0));
    }
    else if(point[1] > DataType_(0))
    {
      //shift our point by 1/2 to the right and scale it to [-1,1]
      new_ref_y_cord.push_back(DataType_(2)*(point[1] - DataType_(0.5)));
      domain_index_y.push_back(Index(1));
    }

    //check the z cordinate
    //Check if its value is near 0
    if(Math::abs(point[2]) < Math::pow(Math::eps<DataType_>(), DataType_(0.9)))
    {
      //then we push both points... whereby we shift them directly on to the edges in x-direction
      new_ref_z_cord.push_back(DataType_(1));
      new_ref_z_cord.push_back(DataType_(-1));
      domain_index_z.push_back(Index(0));
      domain_index_z.push_back(Index(1));
    }
    else if(point[2] < DataType_(0))
    {
      //shift our point by 1/2 to the right and scale it to [-1,1]
      new_ref_z_cord.push_back(DataType_(2)*(point[2] + DataType_(0.5)));
      domain_index_z.push_back(Index(0));
    }
    else if(point[2] > DataType_(0))
    {
      //shift our point by 1/2 to the right and scale it to [-1,1]
      new_ref_z_cord.push_back(DataType_(2)*(point[2] - DataType_(0.5)));
      domain_index_z.push_back(Index(1));
    }

    Tiny::Vector<DataType_, 3> new_point;
    //now iterate over the elements
    for(Index i(0); i < domain_index_x.size(); ++i)
    {
      for(Index j(0); j < domain_index_y.size(); ++j)
      {
        for(Index k(0); k < domain_index_z.size(); ++k)
        {
        //Due to the 2lvl ordering, our refinded cell of cell i is one of 8i, 8i+1,..., 8i+7
        const Index refined_cell_index = Index(8) * cell_ind + Index(4) * domain_index_z[k] + Index(2) * domain_index_y[j] + domain_index_x[i];
        new_point[0] = new_ref_x_cord[i];
        new_point[1] = new_ref_y_cord[j];
        new_point[2] = new_ref_z_cord[k];
        calc_refinement_cell(cells, new_point, refined_cell_index, recursion - Index(1));
        }
      }
    }
  }

    /**
     * \brief This function maps the values of a symmetric 2d matrix, represented by a vector of size dim*(dim+1)/2, into a full dim_*dim_ matrix
     *
     * \param[out] orientation The orientation matrix which should be filled
     *
     * \param[in] vector The vector which holds the values
     *
     * \note Since the mapping between those is non trivial, we only implement this for 3D
     *
     */
    template<typename DataType_,int dim_>
    void fill_tensor2_matrix(Tiny::Matrix<DataType_, dim_, dim_>& DOXY(orientation), const Tiny::Vector<DataType_, dim_ * (dim_ + 1) / 2>& DOXY(vector)){
      XABORTM("Called fill_tensor for something other than dimenison 2 or 3");
    }

    /**
     * \brief Template specialisation for dim_ = 2
     *
     * \note The symmetric matrix is mapped to vector in the following way
     *
     *                      v(0) = A11
     *                      v(1) = A22
     *                      v(2) = A12
     *
     */
    template<typename DataType_>
    void fill_tensor2_matrix(Tiny::Matrix<DataType_, 2, 2>& orientation, const Tiny::Vector<DataType_, 3>& vector)
    {
      orientation[0][0] = vector[0];
      orientation[1][1] = vector[1];
      orientation[0][1] = orientation[1][0] = vector[2];
    }

    /**
     * \brief Template specialisation for dim_ = 3
     *
     * \note The symmetric matrix is mapped to vector in the following way
     *
     *                      v(0) = A11
     *                      v(1) = A22
     *                      v(2) = A33
     *                      v(3) = A12
     *                      v(4) = A13
     *                      v(5) = A23
     *
     */
    template<typename DataType_>
    void fill_tensor2_matrix(Tiny::Matrix<DataType_, 3, 3>& orientation, const Tiny::Vector<DataType_, 6>& vector)
    {
      orientation[0][0] = vector[0];
      orientation[1][1] = vector[1];
      orientation[2][2] = vector[2];
      orientation[0][1] = orientation[1][0] = vector[3];
      orientation[0][2] = orientation[2][0] = vector[4];
      orientation[1][2] = orientation[2][1] = vector[5];
    }
  /**
   * \brief Modified Burgers Operator Assembly class
   *
   * This class is responsible for assembling the scalar and vector-valued modified Burgers operators:
   *
   * \f[\mathbf{N}(v,u,\psi) := \nu \mathbf{L}(u,\psi) + \theta \mathbf{M}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \beta' \mathbf{K'(v,u,\psi)}\f]
   *
   *
   * where
   * - \b L is the diffusive operator with a modified <em>deformation tensor<em>:
   *   \f[\mathbf{L}(u,\psi) := 2 \int_\Omega (\nabla \psi) : (D + N_s(DA + AD))\f]
   *   where <b>A<b> is a 2nd order tensor, describing a 2nd momentum matrix of a probability distribution,
   *         <b>N_s<b> is an additional parameter and
   *         <b>D<b> is the symmetric gradient:
   *   \f[\mathbf{D} := \frac{1}{2}(\nabla u + \nabla u^\top)\f]
   * - \b M is the reactive operator:
   *   \f[\mathbf{M}(u,\psi) := \int_\Omega u\psi\f]
   * - \b K is the convective operator:
   *   \f[\mathbf{K}(v,u,\psi) := \int_\Omega v\cdot \nabla u \psi\f]
   * - <b>K'</b> is the Frechet derivative of the convective operator:
   *   \f[\mathbf{K'}(v,u,\psi) := \int_\Omega \nabla v u \psi\f]
   *
   * \author Maximilian Esser
   */

  template<typename DataType_, typename IndexType_, int dim_>
    class ModBurgersAssembler
    {
      public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      DataType_ nu;

      //scaling parameter for orientation part of viscosoty
      DataType_ N_s;

      /// scaling parameter for reactive operator \b M
      DataType_ theta;

      /// scaling parameter for convective operator \b K
      DataType_ beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      DataType_ frechet_beta;

      /// default constructor
      ModBurgersAssembler() :
        nu(DataType_(1)),
        N_s(DataType_(1)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix.
       *
       * \param[in,out] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] orientation
       * The orientation matrix for the deformation.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Space_, typename CubatureFactory_>
      void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Tiny::Matrix<DataType_, dim_, dim_>& orientation,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;


        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (Math::abs(frechet_beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, dim_, dim_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v;

        loc_v.format();
        mean_v.format();
        loc_grad_v.format();
        loc_A_prod_i.format();
        loc_A_prod_j.format();



        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);


          // format our local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }


            // assemble diffusion matrix?
            if(need_diff)
            {
              // assemble deformation-tensor diffusion and orientation diffusion in one go, careful, as i denotes the test-function and j denotes the trail-functions...
              // this has to be considered in my comprehensions...

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of \psi and A
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);

                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[j].grad, space_data.phi[i].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value1);

                  //add dot(grad_j, grad_i) * A to local_matrix, we can even do this smarter by adding A + Identity to our local_mat...
                  local_matrix[i][j].axpy(N_s*value1, orientation);

                  //calculate scalar product of grad(phi)^T*A*grad(psi) add this to the diag of our local matrix... and of course we can add this to operation above...:
                  const DataType value2 = nu * N_s * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);
                  local_matrix[i][j].add_scalar_main_diag(value2);

                  //calculate vector mult for grad_j
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  //add outer product of grad(phi) * A_prod_psi
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, loc_A_prod_i, N_s * nu * weight);
                  //and the other way round, but transposed
                  local_matrix[i][j].add_outer_product(loc_A_prod_j, space_data.phi[i].grad, N_s * nu * weight);
                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, space_data.phi[i].grad, nu * weight);
                }
              }

            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = frechet_beta * weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(value, loc_grad_v);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

    /**
       * \brief Assembles the Burgers operator into a vector.
       *
       * \param[in,out] vector
       * The vector to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] orientation
       * The 2 nd Order oriantation tensor.
       *
       * \param[in] primal
       * The primal vector, usually a solution vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_, typename CubatureFactory_>
      void assemble_vector(
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Tiny::Matrix<DataType_, dim_, dim_>& orientation,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& primal,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        //const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector-scatter-axpy (if needed)
        typename VectorType::ScatterAxpy scatter_vector(vector);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // create primal gather-axpy
        typename VectorType::GatherAxpy gather_prim(primal);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // local primal vector dofs
        LocalVectorType local_prim_dofs;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v;

        // our local velocity gradient
        //Tiny::Matrix<DataType, dim_, dim_> loc_grad_v;

        loc_v.format();
        //loc_grad_v.format();

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // gather our local primal dofs
          local_prim_dofs.format();
          gather_prim(local_prim_dofs, dof_mapping);

          // format our local vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            /*if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }*/

            // assemble diffusion matrix?
            if(need_diff)
            {

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of A and grad(psi)
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // compute outer product of grad(phi) and grad(psi)
                  const DataType value2 = nu * weight * Tiny::dot(local_prim_dofs[j], space_data.phi[i].grad);

                  //calculate scalar product of grad(phi)^T*A*grad(psi)
                  const DataType value3 = N_s * nu * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);

                  // compute outer product of grad(phi) and A*grad(psi)
                  const DataType value4 = N_s * nu * weight * Tiny::dot(loc_A_prod_i, local_prim_dofs[j]);

                  // compute outer product of A*grad(phi) and grad(psi)
                  const DataType value5 = N_s * value2;

                  //compute matrix vector of A and grad(phi)
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value1+value3, local_prim_dofs[j]);
                  //update dot(grad(psi), grad(phi))* A * primal onto vector... this could be performed together with above operation through adding identity...
                  local_vector[i].add_mat_vec_mult(orientation, local_prim_dofs[j], N_s*value1);
                  local_vector[i].axpy(value2 + value4, space_data.phi[j].grad);
                  local_vector[i].axpy(value5, loc_A_prod_j);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into vector
          scatter_vector(local_vector, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }
    }; //class ModBurgersAssembler

        /**
   * \brief Modified Burgers Full Tensor-formulation Space Dependent Operator Assembly class
   *
   * This class is responsible for assembling the scalar and vector-valued modified Burgers operators:
   *
   * \f[\mathbf{N}(v,u,\psi) := \nu \mathbf{L}(u,\psi) + \theta \mathbf{M}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \beta' \mathbf{K'(v,u,\psi)}\f]
   *
   *
   * where
   * - \b L is the diffusive operator with a modified <em>deformation tensor<em>:
   *   \f[\mathbf{L}(u,\psi) := 2 \int_\Omega (\nabla \psi) : (D + N_s(DA + AD) + \underline{A}:D)\f]
   *   where <b>A<b> is a 2nd order tensor, describing a 2nd momentum matrix of a probability distribution,
   *         which can be dependent on space
   *         <b>_A_<b> is a 4th order tensor, describing the fourth moments of a probability distribution,
   *         which can be space dependent
   *         <b>N_s<b> is an additional parameter and
   *         <b>D<b> is the symmetric gradient:
   *   \f[\mathbf{D} := \frac{1}{2}(\nabla u + \nabla u^\top)\f]
   * - \b M is the reactive operator:
   *   \f[\mathbf{M}(u,\psi) := \int_\Omega u\psi\f]
   * - \b K is the convective operator:
   *   \f[\mathbf{K}(v,u,\psi) := \int_\Omega v\cdot \nabla u \psi\f]
   * - <b>K'</b> is the Frechet derivative of the convective operator:
   *   \f[\mathbf{K'}(v,u,\psi) := \int_\Omega \nabla v u \psi\f]
   *
   * \author Maximilian Esser
   */

  template<typename DataType_, typename IndexType_, int dim_>
    class FullTensorBurgersSDAssembler
    {
      public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      DataType_ nu;

      /// scaling parameter for 2-nd order orientation part of viscosoty
      DataType_ N_s;

      /// scaling parameter for 4-th order orientation part of viscosity
      DataType_ N_p;

      /// scaling parameter for reactive operator \b M
      DataType_ theta;

      /// scaling parameter for convective operator \b K
      DataType_ beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      DataType_ frechet_beta;

      /// default constructor
      FullTensorBurgersSDAssembler() :
        nu(DataType_(1)),
        N_s(DataType_(1)),
        N_p(DataType_(1)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix.
       *
       * \param[in,out] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] function2
       * A dim \times dim dimensional vectorfield function which maps to a dim \times dim tensor in the
       * following way:
       * A[i][j](x) = f_{i + j*dim}(x)
       *
       * \param[in] function4
       * A dim \times dim \times dim \times dim function which maps the 4th order tensor in the following way:
       * A[g][h][i][j] = f_{g*dim^3 + h*dim^2 + i*dim + j}(x)
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename Function2_, typename Function4_>
      void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Function2_ function2,
        const Function4_ function4,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;


        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (Math::abs(frechet_beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function2_> AnalyticEvalTraits2;
        typedef Analytic::EvalTraits<DataType, Function4_> AnalyticEvalTraits4;

        // create a function evaluator
        typename Function2_::template Evaluator<AnalyticEvalTraits2> func_eval2(function2);
        typename Function4_::template Evaluator<AnalyticEvalTraits4> func_eval4(function4);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, dim_, dim_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v, orientation;

        //temp vector to save our orientation data as function output
        Tiny::Vector<DataType, dim_ * dim_> orient_vec;

        //temp vector for our tensor function output
        Tiny::Vector<DataType, dim_ * dim_ * dim_ * dim_> tensor_vec;

        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        mean_v.format();
        loc_grad_v.format();
        orientation.format();
        orient_vec.format();
        loc_A_prod_i.format();
        loc_A_prod_j.format();
        tensor_vec.format();



        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);


          // format our local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            //compute orientation matrix... we should write this as helper function...
            //we need the transformed point, as cubature points are on standarized quad...
            trafo_eval.map_point(point_dom, cubature_rule.get_point(point));
            orient_vec = func_eval2.value(point_dom);
            tensor_vec = func_eval4.value(point_dom);

            for(int i(0); i < dim_; ++i)
            {
              for(int j(0); j < dim_; ++j)
              {
                orientation[i][j] = orient_vec(i + j* dim_);
              }
            }

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }


            // assemble diffusion matrix?
            if(need_diff)
            {
              // assemble deformation-tensor diffusion and orientation diffusion in one go, careful, as i denotes the test-function and j denotes the trail-functions...
              // this has to be considered in my comprehensions...

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of \psi and A
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);

                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[j].grad, space_data.phi[i].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value1);

                  //add dot(grad_j, grad_i) * A to local_matrix, we can even do this smarter by adding A + Identity to our local_mat...
                  local_matrix[i][j].axpy(N_s*value1, orientation);

                  //calculate scalar product of grad(phi)^T*A*grad(psi) add this to the diag of our local matrix... and of course we can add this to operation above...:
                  const DataType value2 = nu * N_s * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);
                  local_matrix[i][j].add_scalar_main_diag(value2);

                  //calculate vector mult for grad_j
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  //add outer product of grad(phi) * A_prod_psi
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, loc_A_prod_i, N_s * nu * weight);
                  //and the other way round, but transposed
                  local_matrix[i][j].add_outer_product(loc_A_prod_j, space_data.phi[i].grad, N_s * nu * weight);
                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, space_data.phi[i].grad, nu * weight);


                  CCND_FIBER::add_tensor4_outer_product_contraction_24(local_matrix[i][j], space_data.phi[i].grad, space_data.phi[j].grad, tensor_vec, N_p * nu * weight);
                  CCND_FIBER::add_tensor4_outer_product_contraction_23(local_matrix[i][j], space_data.phi[i].grad, space_data.phi[j].grad, tensor_vec, N_p * nu * weight);

                }
              }

            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = frechet_beta * weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(value, loc_grad_v);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

    /**
       * \brief Assembles the Burgers operator into a vector.
       *
       * \param[in,out] vector
       * The vector to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] function2
       * Function desribing the space dependent orientation tensor.
       *
       * \param[in] function4
       *
       * \param[in] primal
       * The primal vector, usually a solution vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename Function2_, typename Function4_>
      void assemble_vector(
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Function2_& function2,
        const Function4_& function4,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& primal,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        //const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector-scatter-axpy (if needed)
        typename VectorType::ScatterAxpy scatter_vector(vector);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // create primal gather-axpy
        typename VectorType::GatherAxpy gather_prim(primal);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function2_> AnalyticEvalTraits2;
        typedef Analytic::EvalTraits<DataType, Function4_> AnalyticEvalTraits4;


        // create a function evaluator
        typename Function2_::template Evaluator<AnalyticEvalTraits2> func_eval2(function2);
        typename Function4_::template Evaluator<AnalyticEvalTraits4> func_eval4(function4);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // local primal vector dofs
        LocalVectorType local_prim_dofs;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v;

        // our 2nd order orientation tensor
        Tiny::Matrix<DataType, dim_, dim_> orientation;
        Tiny::Vector<DataType, dim_ * dim_> orient_vec;

        //temp vector for our tensor function output
        Tiny::Vector<DataType, dim_ * dim_ * dim_ * dim_> tensor_vec;
        //temp matrix so save our dot output:
        Tiny::Matrix<DataType, dim_, dim_> tensor_contraction;
        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        //loc_grad_v.format();
        orientation.format();
        orient_vec.format();
        tensor_vec.format();
        tensor_contraction.format();
        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // gather our local primal dofs
          local_prim_dofs.format();
          gather_prim(local_prim_dofs, dof_mapping);

          // format our local vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            //compute orientation matrix... we should write this as helper function...
            //we need the transformed point, as cubature points are on standarized quad...
            trafo_eval.map_point(point_dom, cubature_rule.get_point(point));
            orient_vec = func_eval2.value(point_dom);
            tensor_vec = func_eval4.value(point_dom);

            for(int i(0); i < dim_; ++i)
            {
              for(int j(0); j < dim_; ++j)
              {
                orientation[i][j] = orient_vec(i + j* dim_);
              }
            }


            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            /*if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }*/

            // assemble diffusion matrix?
            if(need_diff)
            {

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of A and grad(psi)
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // compute outer product of grad(phi) and grad(psi)
                  const DataType value2 = nu * weight * Tiny::dot(local_prim_dofs[j], space_data.phi[i].grad);

                  //calculate scalar product of grad(phi)^T*A*grad(psi)
                  const DataType value3 = N_s * nu * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);

                  // compute outer product of grad(phi) and A*grad(psi)
                  const DataType value4 = N_s * nu * weight * Tiny::dot(loc_A_prod_i, local_prim_dofs[j]);

                  // compute outer product of A*grad(phi) and grad(psi)
                  const DataType value5 = N_s * value2;

                  //compute matrix vector of A and grad(phi)
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value1+value3, local_prim_dofs[j]);
                  //update dot(grad(psi), grad(phi))* A * primal onto vector... this could be performed together with above operation through adding identity...
                  local_vector[i].add_mat_vec_mult(orientation, local_prim_dofs[j], N_s*value1);
                  local_vector[i].axpy(value2 + value4, space_data.phi[j].grad);
                  local_vector[i].axpy(value5, loc_A_prod_j);

                  //calculate the tensor dot output
                  CCND_FIBER::set_tensor4_outer_product_contraction_34(tensor_contraction, local_prim_dofs[j], space_data.phi[j].grad, tensor_vec, DataType(1.)); //right alpha?
                  //now add matrix vector prod
                  local_vector[i].add_mat_vec_mult(tensor_contraction, space_data.phi[i].grad, N_p * nu * weight);
                  local_vector[i].add_vec_mat_mult(space_data.phi[i].grad, tensor_contraction, N_p * nu * weight);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into vector
          scatter_vector(local_vector, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }
    }; //class FullTensorBurgersSDAssembler

            /**
   * \brief Modified Burgers Full Finite Element Tensor formulation Operator Assembly class
   *
   * This class is responsible for assembling the scalar and vector-valued modified Burgers operators:
   *
   * \f[\mathbf{N}(v,u,\psi) := \nu \mathbf{L}(u,\psi) + \theta \mathbf{M}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \beta' \mathbf{K'(v,u,\psi)}\f]
   *
   *
   * where
   * - \b L is the diffusive operator with a modified <em>deformation tensor<em>:
   *   \f[\mathbf{L}(u,\psi) := 2 \int_\Omega (\nabla \psi) : (D + N_s(DA + AD) + N_p\underline{A}:D)\f]
   *   where <b>A<b> is a 2nd order tensor, describing a 2nd momentum matrix of a probability distribution,
   *         which can be space dependent
   *         <b>_A_<b> is a 4th order tensor, describing the fourth moments of a probability distribution,
   *         which can be space dependent
   *         <b>N_s, N_p<b> are additional parameters and
   *         <b>D<b> is the symmetric gradient:
   *   \f[\mathbf{D} := \frac{1}{2}(\nabla u + \nabla u^\top)\f]
   * - \b M is the reactive operator:
   *   \f[\mathbf{M}(u,\psi) := \int_\Omega u\psi\f]
   * - \b K is the convective operator:
   *   \f[\mathbf{K}(v,u,\psi) := \int_\Omega v\cdot \nabla u \psi\f]
   * - <b>K'</b> is the Frechet derivative of the convective operator:
   *   \f[\mathbf{K'}(v,u,\psi) := \int_\Omega \nabla v u \psi\f]
   *
   * \author Maximilian Esser
   */

  template<typename DataType_, typename IndexType_, int dim_>
    class FullFETensorBurgersAssembler
    {
      public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      DataType_ nu;

      /// scaling parameter for 2-nd order orientation part of viscosoty
      DataType_ N_s;

      /// scaling parameter for 4-th order orientation part of viscosity
      DataType_ N_p;

      /// scaling parameter for reactive operator \b M
      DataType_ theta;

      /// scaling parameter for convective operator \b K
      DataType_ beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      DataType_ frechet_beta;

      /// default constructor
      FullFETensorBurgersAssembler() :
        nu(DataType_(1)),
        N_s(DataType_(1)),
        N_p(DataType_(1)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix.
       *
       * \param[in,out] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] velo_field
       * A coefficient vector for an input velocity(same dimension as convect, defined in velocity space).
       *
       * \param[in] function2
       * A dim \times dim dimensional vectorfield function which maps to a dim \times dim tensor in the
       * following way:
       * A[i][j](x) = f_{i + j*dim}(v_primal)
       * Thereby the function2 takes the primal velocity vector velo_field as input
       *
       * \param[in] function4
       * A dim \times dim \times dim \times dim function which maps the 4th order tensor in the following way:
       * A[g][h][i][j] = f_{g*dim^3 + h*dim^2 + i*dim + j}(x)
       * Thereby the function4 takes the primal velocity vector velo_field as input
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename Function2_, typename Function4_>
      void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& velocity_field,
        const Function2_& function2,
        const Function4_& function4,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(velocity_field.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;


        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (Math::abs(frechet_beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);
        typename VectorType::GatherAxpy gather_velocity(velocity_field);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function2_> AnalyticEvalTraits2;
        typedef Analytic::EvalTraits<DataType, Function4_> AnalyticEvalTraits4;

        // create a function evaluator
        typename Function2_::template Evaluator<AnalyticEvalTraits2> func_eval2(function2);
        typename Function4_::template Evaluator<AnalyticEvalTraits4> func_eval4(function4);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, dim_, dim_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;

        // local convection field dofs
        LocalVectorType local_conv_dofs, local_velocity_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v, loc_velocity_first;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v, orientation;

        //temp vector to save our orientation data as function output
        Tiny::Vector<DataType, dim_ * dim_> orient_vec;

        //temp vector for our tensor function output
        Tiny::Vector<DataType, dim_ * dim_ * dim_ * dim_> tensor_vec;

        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        loc_velocity_first.format();
        mean_v.format();
        loc_grad_v.format();
        orientation.format();
        orient_vec.format();
        loc_A_prod_i.format();
        loc_A_prod_j.format();
        tensor_vec.format();



        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);
          //should be the same...
          local_velocity_dofs.format();
          gather_velocity(local_velocity_dofs, dof_mapping);


          // format our local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

//             //compute orientation matrix... we should write this as helper function...
//             //we need the transformed point, as cubature points are on standarized quad...
//             trafo_eval.map_point(point_dom, cubature_rule.get_point(point));
//             orient_vec = func_eval2.value(point_dom);
//             tensor_vec = func_eval4.value(point_dom);

//             for(int i(0); i < dim_; ++i)
//             {
//               for(int j(0); j < dim_; ++j)
//               {
//                 orientation[i][j] = orient_vec(i + j* dim_);
//               }
//             }

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }
            //we always need our velocity values
            loc_velocity_first.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              //update the velocity field for our tensor input
              loc_velocity_first.axpy(space_data.phi[i].value, local_velocity_dofs[i]);
            }

            //now we can assemble the local tensor values
            orient_vec = func_eval2.value(loc_velocity_first);
            tensor_vec = func_eval4.value(loc_velocity_first);
            //write the orient_vec into a matrix
            //tensor_vec will be used on as vector
            for(int i(0); i < dim_; ++i)
            {
              for(int j(0); j < dim_; ++j)
              {
                orientation[i][j] = orient_vec(i + j* dim_);
              }
            }


            // assemble diffusion matrix?
            if(need_diff)
            {
              // assemble deformation-tensor diffusion and orientation diffusion in one go, careful, as i denotes the test-function and j denotes the trail-functions...
              // this has to be considered in my comprehensions...

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of \psi and A
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);

                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[j].grad, space_data.phi[i].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value1);

                  //add dot(grad_j, grad_i) * A to local_matrix, we can even do this smarter by adding A + Identity to our local_mat...
                  local_matrix[i][j].axpy(N_s*value1, orientation);

                  //calculate scalar product of grad(phi)^T*A*grad(psi) add this to the diag of our local matrix... and of course we can add this to operation above...:
                  const DataType value2 = nu * N_s * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);
                  local_matrix[i][j].add_scalar_main_diag(value2);

                  //calculate vector mult for grad_j
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  //add outer product of grad(phi) * A_prod_psi
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, loc_A_prod_i, N_s * nu * weight);
                  //and the other way round, but transposed
                  local_matrix[i][j].add_outer_product(loc_A_prod_j, space_data.phi[i].grad, N_s * nu * weight);
                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, space_data.phi[i].grad, nu * weight);


                  CCND_FIBER::add_tensor4_outer_product_contraction_24(local_matrix[i][j], space_data.phi[i].grad, space_data.phi[j].grad, tensor_vec, N_p * nu * weight);
                  CCND_FIBER::add_tensor4_outer_product_contraction_23(local_matrix[i][j], space_data.phi[i].grad, space_data.phi[j].grad, tensor_vec, N_p * nu * weight);

                }
              }

            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = frechet_beta * weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(value, loc_grad_v);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

    /**
       * \brief Assembles the Burgers operator into a vector.
       *
       * \param[in,out] vector
       * The vector to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] velo_field
       * A coefficient vector for an input velocity(same dimension as convect, defined in velocity space).
       *
       * \param[in] function2
       * A dim \times dim dimensional vectorfield function which maps to a dim \times dim tensor in the
       * following way:
       * A[i][j](x) = f_{i + j*dim}(v_primal)
       * Thereby the function2 takes the primal velocity vector velo_field as input
       *
       * \param[in] function4
       * A dim \times dim \times dim \times dim function which maps the 4th order tensor in the following way:
       * A[g][h][i][j] = f_{g*dim^3 + h*dim^2 + i*dim + j}(x)
       * Thereby the function4 takes the primal velocity vector velo_field as input
       *
       * \param[in] primal
       * The primal vector, usually a solution vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename Function2_, typename Function4_>
      void assemble_vector(
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& velocity_field,
        const Function2_& function2,
        const Function4_& function4,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& primal,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(velocity_field.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        //const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector-scatter-axpy (if needed)
        typename VectorType::ScatterAxpy scatter_vector(vector);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);
        typename VectorType::GatherAxpy gather_velocity(velocity_field);

        // create primal gather-axpy
        typename VectorType::GatherAxpy gather_prim(primal);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function2_> AnalyticEvalTraits2;
        typedef Analytic::EvalTraits<DataType, Function4_> AnalyticEvalTraits4;


        // create a function evaluator
        typename Function2_::template Evaluator<AnalyticEvalTraits2> func_eval2(function2);
        typename Function4_::template Evaluator<AnalyticEvalTraits4> func_eval4(function4);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs, local_velocity_dofs;

        // local primal vector dofs
        LocalVectorType local_prim_dofs;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, loc_velocity_first;

        // our 2nd order orientation tensor
        Tiny::Matrix<DataType, dim_, dim_> orientation;
        Tiny::Vector<DataType, dim_ * dim_> orient_vec;

        //temp vector for our tensor function output
        Tiny::Vector<DataType, dim_ * dim_ * dim_ * dim_> tensor_vec;
        //temp matrix so save our dot output:
        Tiny::Matrix<DataType, dim_, dim_> tensor_contraction;
        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        loc_velocity_first.format();
        //loc_grad_v.format();
        orientation.format();
        orient_vec.format();
        tensor_vec.format();
        tensor_contraction.format();
        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          //should be the same...
          local_velocity_dofs.format();
          gather_velocity(local_velocity_dofs, dof_mapping);

          // gather our local primal dofs
          local_prim_dofs.format();
          gather_prim(local_prim_dofs, dof_mapping);

          // format our local vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);



            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }

            //we always need our velocity values
            loc_velocity_first.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              //update the velocity field for our tensor input
              loc_velocity_first.axpy(space_data.phi[i].value, local_velocity_dofs[i]);
            }
            //compute orientation matrix... we should write this as helper function...
            orient_vec = func_eval2.value(loc_velocity_first);
            tensor_vec = func_eval4.value(loc_velocity_first);

            for(int i(0); i < dim_; ++i)
            {
              for(int j(0); j < dim_; ++j)
              {
                orientation[i][j] = orient_vec(i + j* dim_);
              }
            }

            /*if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }*/

            // assemble diffusion matrix?
            if(need_diff)
            {

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of A and grad(psi)
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // compute outer product of grad(phi) and grad(psi)
                  const DataType value2 = nu * weight * Tiny::dot(local_prim_dofs[j], space_data.phi[i].grad);

                  //calculate scalar product of grad(phi)^T*A*grad(psi)
                  const DataType value3 = N_s * nu * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);

                  // compute outer product of grad(phi) and A*grad(psi)
                  const DataType value4 = N_s * nu * weight * Tiny::dot(loc_A_prod_i, local_prim_dofs[j]);

                  // compute outer product of A*grad(phi) and grad(psi)
                  const DataType value5 = N_s * value2;

                  //compute matrix vector of A and grad(phi)
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value1+value3, local_prim_dofs[j]);
                  //update dot(grad(psi), grad(phi))* A * primal onto vector... this could be performed together with above operation through adding identity...
                  local_vector[i].add_mat_vec_mult(orientation, local_prim_dofs[j], N_s*value1);
                  local_vector[i].axpy(value2 + value4, space_data.phi[j].grad);
                  local_vector[i].axpy(value5, loc_A_prod_j);

                  //calculate the tensor dot output
                  CCND_FIBER::set_tensor4_outer_product_contraction_34(tensor_contraction, local_prim_dofs[j], space_data.phi[j].grad, tensor_vec, DataType(1.)); //right alpha?
                  //now add matrix vector prod
                  local_vector[i].add_mat_vec_mult(tensor_contraction, space_data.phi[i].grad, N_p * nu * weight);
                  local_vector[i].add_vec_mat_mult(space_data.phi[i].grad, tensor_contraction, N_p * nu * weight);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into vector
          scatter_vector(local_vector, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }
    }; //class FullFETensorBurgersAssembler


    /**
   * \brief Modified Burgers Full Finite Element Tensor formulation Operator Assembly class
   *
   * This class is responsible for assembling the scalar and vector-valued modified Burgers operators:
   *
   * \f[\mathbf{N}(v,u,\psi) := \nu \mathbf{L}(u,\psi) + \theta \mathbf{M}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \beta' \mathbf{K'(v,u,\psi)}\f]
   *
   *
   * where
   * - \b L is the diffusive operator with a modified <em>deformation tensor<em>:
   *   \f[\mathbf{L}(u,\psi) := 2 \int_\Omega (\nabla \psi) : (D + N_s(DA + AD) + N_p\underline{A}:D)\f]
   *   where <b>A<b> is a 2nd order tensor, describing a 2nd momentum matrix of a probability distribution,
   *         which can be space dependent
   *         <b>_A_<b> is a 4th order tensor, describing the fourth moments of a probability distribution,
   *         which can be space dependent
   *         <b>N_s, N_p<b> are additional parameters and
   *         <b>D<b> is the symmetric gradient:
   *   \f[\mathbf{D} := \frac{1}{2}(\nabla u + \nabla u^\top)\f]
   * - \b M is the reactive operator:
   *   \f[\mathbf{M}(u,\psi) := \int_\Omega u\psi\f]
   * - \b K is the convective operator:
   *   \f[\mathbf{K}(v,u,\psi) := \int_\Omega v\cdot \nabla u \psi\f]
   * - <b>K'</b> is the Frechet derivative of the convective operator:
   *   \f[\mathbf{K'}(v,u,\psi) := \int_\Omega \nabla v u \psi\f]
   *
   * \author Maximilian Esser
   */

  template<typename DataType_, typename IndexType_, int dim_>
    class FullDiscreteTensorBurgersAssembler
    {
      public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      DataType_ nu;

      /// scaling parameter for 2-nd order orientation part of viscosoty
      DataType_ N_s;

      /// scaling parameter for 4-th order orientation part of viscosity
      DataType_ N_p;

      /// scaling parameter for reactive operator \b M
      DataType_ theta;

      /// scaling parameter for convective operator \b K
      DataType_ beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      DataType_ frechet_beta;

      /// default constructor
      FullDiscreteTensorBurgersAssembler() :
        nu(DataType_(1)),
        N_s(DataType_(1)),
        N_p(DataType_(1)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix.
       *
       * \param[in,out] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] velo_obj
       * A wrapper class for the coefficient vector for an input velocity and it space
       *
       * \param[in] level_diff
       * The difference between assembly level and input FEvector level
       *
       * \param[in] function2
       * A dim \times dim dimensional vectorfield function which maps to a dim \times dim tensor in the
       * following way:
       * A[i][j](x) = f_{i + j*dim}(velo_function)
       * Thereby the function2 takes the primal velocity vector velo_field as input
       *
       * \param[in] function4
       * A dim \times dim \times dim \times dim function which maps the 4th order tensor in the following way:
       * A[g][h][i][j] = f_{g*dim^3 + h*dim^2 + i*dim + j}(x)
       * Thereby the function4 takes the primal velocity vector velo_field as input
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename VeloSpace_, typename Function2_, typename Function4_>
      void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const FEVelo<VeloSpace_, DataType_, IndexType_, dim_>& velo_func,
        const IndexType_ level_diff,
        const Function2_& function2,
        const Function4_& function4,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;


        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (Math::abs(frechet_beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function2_> AnalyticEvalTraits2;
        typedef Analytic::EvalTraits<DataType, Function4_> AnalyticEvalTraits4;

        // create a function evaluator
        typename Function2_::template Evaluator<AnalyticEvalTraits2> func_eval2(function2);
        typename Function4_::template Evaluator<AnalyticEvalTraits4> func_eval4(function4);

        //get the trafo of our input FE Vector
//         auto& FE_trafo = velo_func._space.get_trafo();
        //now create a inverse mapping object
        Trafo::InverseMapping<typename std::remove_reference<decltype(velo_func._space.get_trafo())>::type, DataType_> inverse_map(velo_func._space.get_trafo());

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, dim_, dim_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v, loc_velo_vec;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v, orientation;

        //temp vector to save our orientation data as function output
        Tiny::Vector<DataType, dim_ * dim_> orient_vec;

        //temp vector for our tensor function output
        Tiny::Vector<DataType, dim_ * dim_ * dim_ * dim_> tensor_vec;

        //Vector for possible cells for our inverse mapping
        std::vector<IndexType_> cells;
        //we will at maximum expect 8 cells in 3 dimension as candidates
        cells.reserve(8);

        //we will save our point in the domain cell in this vector...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        loc_velo_vec.format();
        mean_v.format();
        loc_grad_v.format();
        orientation.format();
        orient_vec.format();
        loc_A_prod_i.format();
        loc_A_prod_j.format();
        tensor_vec.format();



        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);


          // format our local matrix and vector
          local_matrix.format();


          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            //get the point on the refernce cell
            const auto ref_point = cubature_rule.get_point(point);
            // compute trafo data
            trafo_eval(trafo_data, ref_point);

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

//             //compute orientation matrix... we should write this as helper function...
//             //we need the transformed point, as cubature points are on standarized quad...
            point_dom.format();
            trafo_eval.map_point(point_dom, ref_point);

            //first, find candidate cells for our inverse mapping
            //pop all cells in cells vector
            cells.clear();
            calc_refinement_cell(cells, ref_point, cell, level_diff);
            XASSERTM(cells.size() > 0, "No cells returned...");
            //now we can create our inverse map point
            auto inv_data = inverse_map.unmap_point(point_dom, cells, true);
            //failsafe to much more computational heavy methode, which should guarantee to find the inverse mapping
            if(inv_data.size() == 0)
            {
              inv_data = inverse_map.unmap_point(point_dom);
            }
            //and now create an eval object from a dicrete evaluator
            auto evaldata = Assembly::DiscreteEvaluator::eval_fe_function(inv_data, velo_func._local_FE_coeff, velo_func._space);

            loc_velo_vec.format();
            //and now write the mean_value into our local velo vector
            loc_velo_vec = evaldata.mean_value();

//             orient_vec = func_eval2.value(point_dom);
//             tensor_vec = func_eval4.value(point_dom);

//             for(int i(0); i < dim_; ++i)
//             {
//               for(int j(0); j < dim_; ++j)
//               {
//                 orientation[i][j] = orient_vec(i + j* dim_);
//               }
//             }

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }

            //now we can assemble the local tensor values
            orient_vec = func_eval2.value(loc_velo_vec);
            tensor_vec = func_eval4.value(loc_velo_vec);
            //write the orient_vec into a matrix
            //tensor_vec will be used on as vector
            for(int i(0); i < dim_; ++i)
            {
              for(int j(0); j < dim_; ++j)
              {
                orientation[i][j] = orient_vec(i + j* dim_);
              }
            }


            // assemble diffusion matrix?
            if(need_diff)
            {
              // assemble deformation-tensor diffusion and orientation diffusion in one go, careful, as i denotes the test-function and j denotes the trail-functions...
              // this has to be considered in my comprehensions...

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of \psi and A
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);

                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[j].grad, space_data.phi[i].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value1);

                  //add dot(grad_j, grad_i) * A to local_matrix, we can even do this smarter by adding A + Identity to our local_mat...
                  local_matrix[i][j].axpy(N_s*value1, orientation);

                  //calculate scalar product of grad(phi)^T*A*grad(psi) add this to the diag of our local matrix... and of course we can add this to operation above...:
                  const DataType value2 = nu * N_s * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);
                  local_matrix[i][j].add_scalar_main_diag(value2);

                  //calculate vector mult for grad_j
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  //add outer product of grad(phi) * A_prod_psi
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, loc_A_prod_i, N_s * nu * weight);
                  //and the other way round, but transposed
                  local_matrix[i][j].add_outer_product(loc_A_prod_j, space_data.phi[i].grad, N_s * nu * weight);
                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, space_data.phi[i].grad, nu * weight);


                  CCND_FIBER::add_tensor4_outer_product_contraction_24(local_matrix[i][j], space_data.phi[i].grad, space_data.phi[j].grad, tensor_vec, N_p * nu * weight);
                  CCND_FIBER::add_tensor4_outer_product_contraction_23(local_matrix[i][j], space_data.phi[i].grad, space_data.phi[j].grad, tensor_vec, N_p * nu * weight);

                }
              }

            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = frechet_beta * weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(value, loc_grad_v);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

    /**
       * \brief Assembles the Burgers operator into a vector.
       *
       * \param[in,out] vector
       * The vector to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] velo_obj
       * A wrapper class for the coefficient vector for an input velocity and it space
       *
       * \param[in] level_diff
       * The difference between assembly level and input FEvector level
       *
       * \param[in] function2
       * A dim \times dim dimensional vectorfield function which maps to a dim \times dim tensor in the
       * following way:
       * A[i][j](x) = f_{i + j*dim}(v_primal)
       * Thereby the function2 takes the primal velocity vector velo_field as input
       *
       * \param[in] function4
       * A dim \times dim \times dim \times dim function which maps the 4th order tensor in the following way:
       * A[g][h][i][j] = f_{g*dim^3 + h*dim^2 + i*dim + j}(x)
       * Thereby the function4 takes the primal velocity vector velo_field as input
       *
       * \param[in] primal
       * The primal vector, usually a solution vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename VeloSpace_, typename Function2_, typename Function4_>
      void assemble_vector(
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const FEVelo<VeloSpace_, DataType_, IndexType_, dim_>& velo_func,
        const IndexType_ level_diff,
        const Function2_& function2,
        const Function4_& function4,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& primal,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        //const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector-scatter-axpy (if needed)
        typename VectorType::ScatterAxpy scatter_vector(vector);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // create primal gather-axpy
        typename VectorType::GatherAxpy gather_prim(primal);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function2_> AnalyticEvalTraits2;
        typedef Analytic::EvalTraits<DataType, Function4_> AnalyticEvalTraits4;


        // create a function evaluator
        typename Function2_::template Evaluator<AnalyticEvalTraits2> func_eval2(function2);
        typename Function4_::template Evaluator<AnalyticEvalTraits4> func_eval4(function4);

        //get the trafo of our input FE Vector
//         auto& FE_trafo = velo_func._space.get_trafo();
        //now create a inverse mapping object
        Trafo::InverseMapping<typename std::remove_reference<decltype(velo_func._space.get_trafo())>::type, DataType_> inverse_map(velo_func._space.get_trafo());
        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // local primal vector dofs
        LocalVectorType local_prim_dofs;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, loc_velo_vec;

        // our 2nd order orientation tensor
        Tiny::Matrix<DataType, dim_, dim_> orientation;
        Tiny::Vector<DataType, dim_ * dim_> orient_vec;

        //temp vector for our tensor function output
        Tiny::Vector<DataType, dim_ * dim_ * dim_ * dim_> tensor_vec;
        //temp matrix so save our dot output:
        Tiny::Matrix<DataType, dim_, dim_> tensor_contraction;

        //Vector for possible cells for our inverse mapping
        std::vector<IndexType_> cells;
        //we will at maximum expect 8 cells in 3 dimension as candidates
        cells.reserve(8);

        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        loc_velo_vec.format();
        //loc_grad_v.format();
        orientation.format();
        orient_vec.format();
        tensor_vec.format();
        tensor_contraction.format();


        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // gather our local primal dofs
          local_prim_dofs.format();
          gather_prim(local_prim_dofs, dof_mapping);

          // format our local vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            //get the point on the refernce cell
            const auto ref_point = cubature_rule.get_point(point);

            // compute trafo data
            trafo_eval(trafo_data, ref_point);

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }

            //compute orientation matrix... we should write this as helper function...
            //we need the transformed point, as cubature points are on standarized quad...
            point_dom.format();
            trafo_eval.map_point(point_dom, ref_point);

            //first, find candidate cells for our inverse mapping
            //pop all cells in cells vector
            cells.clear();
            calc_refinement_cell(cells, ref_point, cell, level_diff);
            //now we can create our inverse map point
            auto inv_data = inverse_map.unmap_point(point_dom, cells);
            //and now create an eval object from a dicrete evaluator
            auto evaldata = Assembly::DiscreteEvaluator::eval_fe_function(inv_data, velo_func._local_FE_coeff, velo_func._space);

            //we always need our velocity values
            loc_velo_vec.format();
            //and now write the mean_value into our local velo vector
            loc_velo_vec = evaldata.mean_value();

            //compute orientation matrix... we should write this as helper function...
            orient_vec = func_eval2.value(loc_velo_vec);
            tensor_vec = func_eval4.value(loc_velo_vec);

            for(int i(0); i < dim_; ++i)
            {
              for(int j(0); j < dim_; ++j)
              {
                orientation[i][j] = orient_vec(i + j* dim_);
              }
            }

            /*if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }*/

            // assemble diffusion matrix?
            if(need_diff)
            {

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of A and grad(psi)
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // compute outer product of grad(phi) and grad(psi)
                  const DataType value2 = nu * weight * Tiny::dot(local_prim_dofs[j], space_data.phi[i].grad);

                  //calculate scalar product of grad(phi)^T*A*grad(psi)
                  const DataType value3 = N_s * nu * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);

                  // compute outer product of grad(phi) and A*grad(psi)
                  const DataType value4 = N_s * nu * weight * Tiny::dot(loc_A_prod_i, local_prim_dofs[j]);

                  // compute outer product of A*grad(phi) and grad(psi)
                  const DataType value5 = N_s * value2;

                  //compute matrix vector of A and grad(phi)
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value1+value3, local_prim_dofs[j]);
                  //update dot(grad(psi), grad(phi))* A * primal onto vector... this could be performed together with above operation through adding identity...
                  local_vector[i].add_mat_vec_mult(orientation, local_prim_dofs[j], N_s*value1);
                  local_vector[i].axpy(value2 + value4, space_data.phi[j].grad);
                  local_vector[i].axpy(value5, loc_A_prod_j);

                  //calculate the tensor dot output
                  CCND_FIBER::set_tensor4_outer_product_contraction_34(tensor_contraction, local_prim_dofs[j], space_data.phi[j].grad, tensor_vec, DataType(1.)); //right alpha?
                  //now add matrix vector prod
                  local_vector[i].add_mat_vec_mult(tensor_contraction, space_data.phi[i].grad, N_p * nu * weight);
                  local_vector[i].add_vec_mat_mult(space_data.phi[i].grad, tensor_contraction, N_p * nu * weight);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into vector
          scatter_vector(local_vector, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }
    }; //class FullDiscreteTensorBurgersAssembler

   /**
   * \brief Modified Burgers Full Finite Element Tensor formulation Operator Assembly class
   *
   * This class is responsible for assembling the scalar and vector-valued modified Burgers operators:
   *
   * \f[\mathbf{N}(v,u,\psi) := \nu \mathbf{L}(u,\psi) + \theta \mathbf{M}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \beta' \mathbf{K'(v,u,\psi)}\f]
   *
   *
   * where
   * - \b L is the diffusive operator with a modified <em>deformation tensor<em>:
   *   \f[\mathbf{L}(u,\psi) := 2 \int_\Omega (\nabla \psi) : (D + N_s(DA + AD) + N_p\underline{A}:D)\f]
   *   where <b>A<b> is a 2nd order tensor, describing a 2nd momentum matrix of a probability distribution,
   *         which can be space dependent
   *         <b>_A_<b> is a 4th order tensor, describing the fourth moments of a probability distribution,
   *         which can be space dependent
   *         <b>N_s, N_p<b> are additional parameters and
   *         <b>D<b> is the symmetric gradient:
   *   \f[\mathbf{D} := \frac{1}{2}(\nabla u + \nabla u^\top)\f]
   * - \b M is the reactive operator:
   *   \f[\mathbf{M}(u,\psi) := \int_\Omega u\psi\f]
   * - \b K is the convective operator:
   *   \f[\mathbf{K}(v,u,\psi) := \int_\Omega v\cdot \nabla u \psi\f]
   * - <b>K'</b> is the Frechet derivative of the convective operator:
   *   \f[\mathbf{K'}(v,u,\psi) := \int_\Omega \nabla v u \psi\f]
   *
   * \author Maximilian Esser
   */

  template<typename DataType_, typename IndexType_, int dim_>
    class FullFiberOrientationTensorBurgersAssembler
    {
      public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      DataType_ nu;

      /// scaling parameter for 2-nd order orientation part of viscosoty
      DataType_ N_s;

      /// scaling parameter for 4-th order orientation part of viscosity
      DataType_ N_p;

      /// scaling parameter for reactive operator \b M
      DataType_ theta;

      /// scaling parameter for convective operator \b K
      DataType_ beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      DataType_ frechet_beta;

      /// default constructor
      FullFiberOrientationTensorBurgersAssembler() :
        nu(DataType_(1)),
        N_s(DataType_(0)),
        N_p(DataType_(0)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix.
       *
       * \param[in,out] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] tensor2
       * A BlockedVector of size of convect, where each block of size dim*(dim+1)/2 represents a symmetric dim\times dim matrix mapped in a specific way.
       * See the mapping function for specifics.
       *
       * \param[in] tensor4
       * A BlockedVector of size of convect, where each block of size dim*(dim+1)*(dim+2)*(dim+3)/24 rpresents a symmetric 4th order tensor of size dim * dim * dim * dim.
       * For specifcs, see mapping.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Space_, typename CubatureFactory_>
      void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_ * (dim_ +1)/2>& tensor2,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_*(dim_+1)*(dim_+2)*(dim_+3)/24>& tensor4,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(tensor2.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(tensor4.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_*(dim_+1)/2> Orient2Type;
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_*(dim_+1)*(dim_+2)*(dim_+3)/24> Orient4Type;


        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (Math::abs(frechet_beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create gather-axpy for the input coefficient vectors
        typename VectorType::GatherAxpy gather_conv(convect);
        typename Orient2Type::GatherAxpy gather_tensor2(tensor2);
        typename Orient4Type::GatherAxpy gather_tensor4(tensor4);


        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, dim_, dim_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        //local vector data for tensor2
        typedef Tiny::Vector<DataType, dim_*(dim_+1)/2> Tensor2Value;
        typedef Tiny::Vector<Tensor2Value, max_local_dofs> LocalTensor2Type;
        //local vector data for tensor4
        typedef Tiny::Vector<DataType, dim_*(dim_+1)*(dim_+2)*(dim_+3)/24> Tensor4Value;
        typedef Tiny::Vector<Tensor4Value, max_local_dofs> LocalTensor4Type;

        // local convection field dofs
        LocalVectorType local_conv_dofs;
        //local tensor2 field dofs
        LocalTensor2Type local_tensor2_dofs;
        //local tensor4 field dofs
        LocalTensor4Type local_tensor4_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v;

        //our local 2nd order matrix
        Tiny::Matrix<DataType, dim_, dim_> orientation;

        //our local vector for the 2nd tensor
        Tensor2Value loc_tensor2_vec;

        //our local vector for the 4th tensor
        Tensor4Value loc_tensor4_vec;


        loc_v.format();
        mean_v.format();
        loc_grad_v.format();
        loc_tensor2_vec.format();
        loc_tensor4_vec.format();
        orientation.format();
        loc_A_prod_i.format();
        loc_A_prod_j.format();




        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);
          //should be the same for the tensors
          local_tensor2_dofs.format();
          gather_tensor2(local_tensor2_dofs, dof_mapping);

          local_tensor4_dofs.format();
          gather_tensor4(local_tensor4_dofs, dof_mapping);



          // format our local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);



            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }
            //we always need our tensors
            loc_tensor2_vec.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              //update the velocity field for our tensor input
              loc_tensor2_vec.axpy(space_data.phi[i].value, local_tensor2_dofs[i]);
            }
            //write the entries into a matrix
            fill_tensor2_matrix(orientation, loc_tensor2_vec);

            loc_tensor4_vec.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              //update the velocity field for our tensor input
              loc_tensor4_vec.axpy(space_data.phi[i].value, local_tensor4_dofs[i]);
            }



            // assemble diffusion matrix?
            if(need_diff)
            {
              // assemble deformation-tensor diffusion and orientation diffusion in one go, careful, as i denotes the test-function and j denotes the trail-functions...

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of \psi and A
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);

                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[j].grad, space_data.phi[i].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value1);

                  //add dot(grad_j, grad_i) * A to local_matrix, we can even do this smarter by adding A + Identity to our local_mat...
                  local_matrix[i][j].axpy(N_s*value1, orientation);

                  //calculate scalar product of grad(phi)^T*A*grad(psi) add this to the diag of our local matrix... and of course we can add this to operation above...:
                  const DataType value2 = nu * N_s * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);
                  local_matrix[i][j].add_scalar_main_diag(value2);

                  //calculate vector mult for grad_j
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  //add outer product of grad(phi) * A_prod_psi
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, loc_A_prod_i, N_s * nu * weight);
                  //and the other way round, but transposed
                  local_matrix[i][j].add_outer_product(loc_A_prod_j, space_data.phi[i].grad, N_s * nu * weight);
                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, space_data.phi[i].grad, nu * weight);

                  //add the tensor4 part... since its symmetric we could reduce this to one operation...
                  CCND_FIBER::add_tensor4_outer_product_contraction_symmetric(local_matrix[i][j], space_data.phi[i].grad, space_data.phi[j].grad, loc_tensor4_vec, DataType(2) * N_p * nu * weight);
//                   CCND_FIBER::add_tensor4_outer_product_contraction_symmetric(local_matrix[i][j], space_data.phi[i].grad, space_data.phi[j].grad, loc_tensor4_vec, N_p * nu * weight);


                }
              }

            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = frechet_beta * weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(value, loc_grad_v);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

    /**
       * \brief Assembles the Burgers operator into a vector.
       *
       * \param[in,out] vector
       * The vector to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] tensor2
       * A BlockedVector of size of convect, where each block of size dim*(dim+1)/2 represents a symmetric dim\times dim matrix mapped in a specific way.
       * See the mapping function for specifics.
       *
       * \param[in] tensor4
       * A BlockedVector of size of convect, where each block of size dim*(dim+1)*(dim+2)*(dim+3)/24 rpresents a symmetric 4th order tensor of size dim * dim * dim * dim.
       * For specifcs, see mapping.
       *
       * \param[in] primal
       * The primal vector, usually a solution vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_, typename CubatureFactory_>
      void assemble_vector(
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_ * (dim_ +1)/2>& tensor2,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_*(dim_+1)*(dim_+2)*(dim_+3)/24>& tensor4,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& primal,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(tensor2.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(tensor4.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_*(dim_+1)/2> Orient2Type;
        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_*(dim_+1)*(dim_+2)*(dim_+3)/24> Orient4Type;

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        //const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector-scatter-axpy (if needed)
        typename VectorType::ScatterAxpy scatter_vector(vector);

        // create gather-axpys
        typename VectorType::GatherAxpy gather_conv(convect);
        typename Orient2Type::GatherAxpy gather_tensor2(tensor2);
        typename Orient4Type::GatherAxpy gather_tensor4(tensor4);

        // create primal gather-axpy
        typename VectorType::GatherAxpy gather_prim(primal);


        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        //local vector data for tensor2
        typedef Tiny::Vector<DataType, dim_*(dim_+1)/2> Tensor2Value;
        typedef Tiny::Vector<Tensor2Value, max_local_dofs> LocalTensor2Type;
        //local vector data for tensor4
        typedef Tiny::Vector<DataType, dim_*(dim_+1)*(dim_+2)*(dim_+3)/24> Tensor4Value;
        typedef Tiny::Vector<Tensor4Value, max_local_dofs> LocalTensor4Type;

        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs;
        //local tensor2 field dofs
        LocalTensor2Type local_tensor2_dofs;
        //local tensor4 field dofs
        LocalTensor4Type local_tensor4_dofs;

        // local primal vector dofs
        LocalVectorType local_prim_dofs;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v;

        // our 2nd order orientation tensor
        Tiny::Matrix<DataType, dim_, dim_> orientation;

        //our local vector for the 2nd tensor
        Tensor2Value loc_tensor2_vec;

        //our local vector for the 4th tensor
        Tensor4Value loc_tensor4_vec;

        //temp matrix so save our dot output:
        Tiny::Matrix<DataType, dim_, dim_> tensor_contraction;

        loc_v.format();
        //loc_grad_v.format();
        orientation.format();
        loc_tensor2_vec.format();
        loc_tensor4_vec.format();
        tensor_contraction.format();
        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);
          //should be the same for the tensors
          local_tensor2_dofs.format();
          gather_tensor2(local_tensor2_dofs, dof_mapping);

          local_tensor4_dofs.format();
          gather_tensor4(local_tensor4_dofs, dof_mapping);

          // gather our local primal dofs
          local_prim_dofs.format();
          gather_prim(local_prim_dofs, dof_mapping);

          // format our local vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);



            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            //we always need our tensors
            loc_tensor2_vec.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              //update the velocity field for our tensor input
              loc_tensor2_vec.axpy(space_data.phi[i].value, local_tensor2_dofs[i]);
            }
            //write the entries into a matrix
            fill_tensor2_matrix(orientation, loc_tensor2_vec);

            loc_tensor4_vec.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              //update the velocity field for our tensor input
              loc_tensor4_vec.axpy(space_data.phi[i].value, local_tensor4_dofs[i]);
            }


            /*if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }*/

            // assemble diffusion matrix?
            if(need_diff)
            {

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of A and grad(psi)
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // compute outer product of grad(phi) and grad(psi)
                  const DataType value2 = nu * weight * Tiny::dot(local_prim_dofs[j], space_data.phi[i].grad);

                  //calculate scalar product of grad(phi)^T*A*grad(psi)
                  const DataType value3 = N_s * nu * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);

                  // compute outer product of grad(phi) and A*grad(psi)
                  const DataType value4 = N_s * nu * weight * Tiny::dot(loc_A_prod_i, local_prim_dofs[j]);

                  // compute outer product of A*grad(phi) and grad(psi)
                  const DataType value5 = N_s * value2;

                  //compute matrix vector of A and grad(phi)
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value1+value3, local_prim_dofs[j]);
                  //update dot(grad(psi), grad(phi))* A * primal onto vector... this could be performed together with above operation through adding identity...
                  local_vector[i].add_mat_vec_mult(orientation, local_prim_dofs[j], N_s*value1);
                  local_vector[i].axpy(value2 + value4, space_data.phi[j].grad);
                  local_vector[i].axpy(value5, loc_A_prod_j);

                  //calculate the tensor dot output
                  CCND_FIBER::set_tensor4_outer_product_contraction_symmetric(tensor_contraction, local_prim_dofs[j], space_data.phi[j].grad, loc_tensor4_vec, DataType(1.)); //right alpha?
                  //now add matrix vector prod
                  local_vector[i].add_mat_vec_mult(tensor_contraction, space_data.phi[i].grad, N_p * nu * weight);
                  local_vector[i].add_vec_mat_mult(space_data.phi[i].grad, tensor_contraction, N_p * nu * weight);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into vector
          scatter_vector(local_vector, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }
    }; //class FullFiberOrientationTensorBurgersAssembler



   /**
   * \brief Modified Burgers Space Dependent Viscosity Operator Assembly class
   *
   * This class is responsible for assembling the scalar and vector-valued modified Burgers operators with space dependent Viscosity:
   *
   * \f[\mathbf{N}(v,u,\psi) := \mathbf{L}(u,\psi) + \theta \mathbf{M}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \beta' \mathbf{K'(v,u,\psi)}\f]
   *
   *
   * where
   * - \b L is the diffusive operator with a modified <em>deformation tensor<em>:
   *   \f[\mathbf{L}(u,\psi) := 2 \int_\Omega (\nabla \psi) : (\nu(x) D)\f]
   *   where <b>A<b> is a 2nd order tensor, describing a 2nd momentum matrix of a probability distribution,
   *         which can be dependent on space
   *         <b>N_s<b> is an additional parameter and
   *         <b>D<b> is the symmetric gradient:
   *   \f[\mathbf{D} := \frac{1}{2}(\nabla u + \nabla u^\top)\f]
   * - \b M is the reactive operator:
   *   \f[\mathbf{M}(u,\psi) := \int_\Omega u\psi\f]
   * - \b K is the convective operator:
   *   \f[\mathbf{K}(v,u,\psi) := \int_\Omega v\cdot \nabla u \psi\f]
   * - <b>K'</b> is the Frechet derivative of the convective operator:
   *   \f[\mathbf{K'}(v,u,\psi) := \int_\Omega \nabla v u \psi\f]
   *
   * \author Maximilian Esser
   */

  template<typename DataType_, typename IndexType_, int dim_>
    class ModBurgersSDViscosityAssembler
    {
      public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      DataType_ nu;

      /// scaling parameter for reactive operator \b M
      DataType_ theta;

      /// scaling parameter for convective operator \b K
      DataType_ beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      DataType_ frechet_beta;

      /// default constructor
      ModBurgersSDViscosityAssembler() :
        nu(DataType_(1)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix.
       *
       * \param[in,out] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] function
       * A scalar function which maps from a dim size space to a scalar
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename Function_>
      void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Function_& function,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;


        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (Math::abs(frechet_beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

        // create a function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval(function);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, dim_, dim_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v;

        //DataType to save our function output
        DataType function_output;

        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        mean_v.format();
        loc_grad_v.format();



        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);


          // format our local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            //compute orientation matrix... we should write this as helper function...
            //we need the transformed point, as cubature points are on standarized quad...
            trafo_eval.map_point(point_dom, cubature_rule.get_point(point));
            function_output = func_eval.value(point_dom);

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }


            // assemble diffusion matrix?
            if(need_diff)
            {
              // assemble deformation-tensor diffusion and orientation diffusion in one go, careful, as i denotes the test-function and j denotes the trail-functions...
              // this has to be considered in my comprehensions...

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {

                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value = function_output * weight * Tiny::dot(space_data.phi[j].grad, space_data.phi[i].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);

                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, space_data.phi[i].grad, function_output * weight);
                }
              }

            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = frechet_beta * weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(value, loc_grad_v);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

    /**
       * \brief Assembles the Burgers operator into a vector.
       *
       * \param[in,out] vector
       * The vector to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] function
       * Function desribing the space dependent orientation tensor.
       *
       * \param[in] primal
       * The primal vector, usually a solution vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename Function_>
      void assemble_vector(
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Function_& function,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& primal,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        //const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector-scatter-axpy (if needed)
        typename VectorType::ScatterAxpy scatter_vector(vector);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // create primal gather-axpy
        typename VectorType::GatherAxpy gather_prim(primal);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

        // create a function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval(function);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // local primal vector dofs
        LocalVectorType local_prim_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v;

        //DataType to save our function output
        DataType function_output;

        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        //loc_grad_v.format();

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // gather our local primal dofs
          local_prim_dofs.format();
          gather_prim(local_prim_dofs, dof_mapping);

          // format our local vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            //compute orientation matrix... we should write this as helper function...
            //we need the transformed point, as cubature points are on standarized quad...
            trafo_eval.map_point(point_dom, cubature_rule.get_point(point));
            function_output = func_eval.value(point_dom);


            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            /*if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }*/

            // assemble diffusion matrix?
            if(need_diff)
            {

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = function_output * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // compute outer product of grad(phi) and grad(psi)
                  const DataType value2 = function_output * weight * Tiny::dot(local_prim_dofs[j], space_data.phi[i].grad);

                  // update local vector
                  local_vector[i].axpy(value1, local_prim_dofs[j]);
                  local_vector[i].axpy(value2, space_data.phi[j].grad);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into vector
          scatter_vector(local_vector, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }
    }; //class ModBurgersSDViscosityAssembler




    /**
   * \brief Modified Burgers nonlinear Viscosity Operator Assembly class
   *
   * This class is responsible for assembling the scalar and vector-valued modified Burgers operators with non linear
   * viscosity:
   *
   * \f[\mathbf{N}(v,u,\psi) := \mathbf{L}(u,\psi) + \mathbf{L'}(v,u,\psi) + \theta \mathbf{M}(u,\psi)
   *    + \beta \mathbf{K}(v,u,\psi) + \beta' \mathbf{K'(v,u,\psi)}\f]
   *
   *
   * where
   * - \b L is the diffusive operator with a modified <em>deformation tensor<em>:
   *   \f[\mathbf{L}(u,\psi) := 2 \int_\Omega (\nabla \psi) : (\nu(v)D)\f]
   *   where <b>A</b> is a 2nd order tensor, describing a 2nd momentum matrix of a probability distribution,
   *         which can be dependent on space
   *         <b>N_s</b> is an additional parameter and
   *         <b>D</b> is the symmetric gradient:
   *   \f[\mathbf{D} := \frac{1}{2}(\nabla u + \nabla u^\top)\f]
   * - \f[\nu(t)\f]
   *   where is <b>\nu </b> is a (differentiable) function and
   *   <b> t </b> can be expressed as
   *   \f[t := \frac{1}{2} D(v) : D(v)\f]
   * - \b M is the reactive operator:
   *   \f[\mathbf{M}(u,\psi) := \int_\Omega u\psi\f]
   * - \b K is the convective operator:
   *   \f[\mathbf{K}(v,u,\psi) := \int_\Omega v\cdot \nabla u \psi\f]
   * - <b>K'</b> is the Frechet derivative of the convective operator:
   *   \f[\mathbf{K'}(v,u,\psi) := \int_\Omega \nabla v u \psi\f]
   * - <b>L'</b> is the Frechet derivative of the diffusion operator:
   *   \f[\mathbf{L'}(v,i,\psi) := 2 \int_\Omega (\nabla \psi) : (\nu'(v)(D(v) : D(u))D(v))\f]
   *
   * \author Maximilian Esser
   */

  template<typename DataType_, typename IndexType_, int dim_>
    class ModBurgersViscosityAssembler
    {
      public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      DataType_ nu;

      //scaling parameter for orientation part of viscosoty
      DataType_ N_s;

      /// scaling parameter for reactive operator \b M
      DataType_ theta;

      /// scaling parameter for convective operator \b K
      DataType_ beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      DataType_ frechet_beta;

      /// default constructor
      ModBurgersViscosityAssembler() :
        nu(DataType_(1)),
        N_s(DataType_(1)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix.
       *
       * \param[in,out] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] function_visco
       * A scalar to scalar function which maps  the trace of the symmetric tensor to a real, positive value
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename Function_>
      void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Function_& function,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;


        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (Math::abs(frechet_beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

        // create a function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval(function);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, dim_, dim_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v, loc_v_defo, loc_v_squared;

        //DataType to save our output of function(defo_trace)
        DataType function_output;

        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        mean_v.format();
        loc_grad_v.format();



        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);


          // format our local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            //we always need this for the viscosity
            loc_grad_v.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // update velocity gradient
              loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
            }
            //we need the strainrate tensor D(v):
            loc_v_defo.set_transpose(loc_grad_v);
            loc_v_defo.axpy(loc_grad_v, DataType(1));
            //scale by 1/2.... we should implement a scale function...
            for(IndexType_ i = 0; i < dim_; ++i)
            {
              for(IndexType_ j = 0; j < dim_; ++j)
              {
                loc_v_defo[i][j] *= DataType(0.5);
              }
            }
            //calculate our function output
            //first calculate D^T * D = D * D
            //again, this can be speed up by implementing only the necessary steps, as we are only interested in the trace
            loc_v_squared.set_mat_mat_mult(loc_v_defo, loc_v_defo);
            function_output = func_eval(DataType(0.5) * loc_v_squared.trace());

            // assemble diffusion matrix?
            if(need_diff)
            {

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {

                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value = function_output * weight * Tiny::dot(space_data.phi[j].grad, space_data.phi[i].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);

                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, space_data.phi[i].grad, function_output * weight);
                }
              }

            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = frechet_beta * weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(value, loc_grad_v);
                }
              }
            }

            //TODO add diffusive Frechet part...

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

    /**
       * \brief Assembles the Burgers operator into a vector.
       *
       * \param[in,out] vector
       * The vector to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] function
       * Function desribing the space dependent orientation tensor.
       *
       * \param[in] primal
       * The primal vector, usually a solution vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename Function_>
      void assemble_vector(
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Function_& function,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& primal,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        //const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector-scatter-axpy (if needed)
        typename VectorType::ScatterAxpy scatter_vector(vector);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // create primal gather-axpy
        typename VectorType::GatherAxpy gather_prim(primal);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

        // create a function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval(function);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // local primal vector dofs
        LocalVectorType local_prim_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v, loc_v_defo, loc_v_squared;

        //DataType to save our output of function(defo_trace)
        DataType function_output;

        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        //loc_grad_v.format();

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // gather our local primal dofs
          local_prim_dofs.format();
          gather_prim(local_prim_dofs, dof_mapping);

          // format our local vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }

            loc_grad_v.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // update velocity gradient
              loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
            }
            //we need the strainrate tensor D(v):
            loc_v_defo.set_transpose(loc_grad_v);
            loc_v_defo.axpy(loc_grad_v, DataType(1));
            //scale by 1/2.... we should implement a scale function...
            for(IndexType_ i = 0; i < dim_; ++i)
            {
              for(IndexType_ j = 0; j < dim_; ++j)
              {
                loc_v_defo[i][j] *= DataType(0.5);
              }
            }
            //calculate our function output
            //first calculate D^T * D = D * D
            //again, this can be speed up by implementing only the necessary steps, as we are only interested in the trace
            loc_v_squared.set_mat_mat_mult(loc_v_defo, loc_v_defo);
            function_output = func_eval(DataType(0.5) * loc_v_squared.trace());

            // assemble diffusion matrix?
            if(need_diff)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = function_output * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // compute outer product of grad(phi) and grad(psi)
                  const DataType value2 = function_output * weight * Tiny::dot(local_prim_dofs[j], space_data.phi[i].grad);


                  // update local vector
                  local_vector[i].axpy(value1, local_prim_dofs[j]);
                  local_vector[i].axpy(value2, space_data.phi[j].grad);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into vector
          scatter_vector(local_vector, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }
    }; //class ModBurgersViscosityAssembler

        /**
   * \brief Modified Burgers Space Dependent Operator Assembly class
   *
   * This class is responsible for assembling the scalar and vector-valued modified Burgers operators:
   *
   * \f[\mathbf{N}(v,u,\psi) := \nu \mathbf{L}(u,\psi) + \theta \mathbf{M}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \beta' \mathbf{K'(v,u,\psi)}\f]
   *
   *
   * where
   * - \b L is the diffusive operator with a modified <em>deformation tensor<em>:
   *   \f[\mathbf{L}(u,\psi) := 2 \int_\Omega (\nabla \psi) : (D + N_s(DA + AD))\f]
   *   where <b>A<b> is a 2nd order tensor, describing a 2nd momentum matrix of a probability distribution,
   *         which can be dependent on space
   *         <b>N_s<b> is an additional parameter and
   *         <b>D<b> is the symmetric gradient:
   *   \f[\mathbf{D} := \frac{1}{2}(\nabla u + \nabla u^\top)\f]
   * - \b M is the reactive operator:
   *   \f[\mathbf{M}(u,\psi) := \int_\Omega u\psi\f]
   * - \b K is the convective operator:
   *   \f[\mathbf{K}(v,u,\psi) := \int_\Omega v\cdot \nabla u \psi\f]
   * - <b>K'</b> is the Frechet derivative of the convective operator:
   *   \f[\mathbf{K'}(v,u,\psi) := \int_\Omega \nabla v u \psi\f]
   *
   * \author Maximilian Esser
   */

  template<typename DataType_, typename IndexType_, int dim_>
    class ModBurgersSDAssembler
    {
      public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      DataType_ nu;

      //scaling parameter for orientation part of viscosoty
      DataType_ N_s;

      /// scaling parameter for reactive operator \b M
      DataType_ theta;

      /// scaling parameter for convective operator \b K
      DataType_ beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      DataType_ frechet_beta;

      /// default constructor
      ModBurgersSDAssembler() :
        nu(DataType_(1)),
        N_s(DataType_(1)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix.
       *
       * \param[in,out] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] function
       * A dim \times dim dimensional vectorfield function which maps to a dim \times dim tensor in the
       * following way:
       * A[i][j](x) = f_{i + j*dim}(x)
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename Function_>
      void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Function_& function,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;


        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (Math::abs(frechet_beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

        // create a function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval(function);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, dim_, dim_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v, orientation;

        //temp vector to save our orientation data as function output
        Tiny::Vector<DataType, dim_ * dim_> orient_vec;

        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        mean_v.format();
        loc_grad_v.format();
        orientation.format();
        orient_vec.format();
        loc_A_prod_i.format();
        loc_A_prod_j.format();



        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);


          // format our local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            //compute orientation matrix... we should write this as helper function...
            //we need the transformed point, as cubature points are on standarized quad...
            trafo_eval.map_point(point_dom, cubature_rule.get_point(point));
            orient_vec = func_eval.value(point_dom);

            for(int i(0); i < dim_; ++i)
            {
              for(int j(0); j < dim_; ++j)
              {
                orientation[i][j] = orient_vec(i + j* dim_);
              }
            }

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }


            // assemble diffusion matrix?
            if(need_diff)
            {
              // assemble deformation-tensor diffusion and orientation diffusion in one go, careful, as i denotes the test-function and j denotes the trail-functions...
              // this has to be considered in my comprehensions...

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of \psi and A
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);

                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[j].grad, space_data.phi[i].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value1);

                  //add dot(grad_j, grad_i) * A to local_matrix, we can even do this smarter by adding A + Identity to our local_mat...
                  local_matrix[i][j].axpy(N_s*value1, orientation);

                  //calculate scalar product of grad(phi)^T*A*grad(psi) add this to the diag of our local matrix... and of course we can add this to operation above...:
                  const DataType value2 = nu * N_s * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);
                  local_matrix[i][j].add_scalar_main_diag(value2);

                  //calculate vector mult for grad_j
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  //add outer product of grad(phi) * A_prod_psi
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, loc_A_prod_i, N_s * nu * weight);
                  //and the other way round, but transposed
                  local_matrix[i][j].add_outer_product(loc_A_prod_j, space_data.phi[i].grad, N_s * nu * weight);
                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, space_data.phi[i].grad, nu * weight);
                }
              }

            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = frechet_beta * weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(value, loc_grad_v);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

    /**
       * \brief Assembles the Burgers operator into a vector.
       *
       * \param[in,out] vector
       * The vector to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] function
       * Function desribing the space dependent orientation tensor.
       *
       * \param[in] primal
       * The primal vector, usually a solution vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename Function_>
      void assemble_vector(
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Function_& function,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& primal,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        //const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector-scatter-axpy (if needed)
        typename VectorType::ScatterAxpy scatter_vector(vector);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // create primal gather-axpy
        typename VectorType::GatherAxpy gather_prim(primal);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

        // create a function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval(function);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // local primal vector dofs
        LocalVectorType local_prim_dofs;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> orientation;
        Tiny::Vector<DataType, dim_ * dim_> orient_vec;
        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        //loc_grad_v.format();
        orientation.format();
        orient_vec.format();
        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // gather our local primal dofs
          local_prim_dofs.format();
          gather_prim(local_prim_dofs, dof_mapping);

          // format our local vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            //compute orientation matrix... we should write this as helper function...
            //we need the transformed point, as cubature points are on standarized quad...
            trafo_eval.map_point(point_dom, cubature_rule.get_point(point));
            orient_vec = func_eval.value(point_dom);

            for(int i(0); i < dim_; ++i)
            {
              for(int j(0); j < dim_; ++j)
              {
                orientation[i][j] = orient_vec(i + j* dim_);
              }
            }


            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            /*if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }*/

            // assemble diffusion matrix?
            if(need_diff)
            {

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of A and grad(psi)
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // compute outer product of grad(phi) and grad(psi)
                  const DataType value2 = nu * weight * Tiny::dot(local_prim_dofs[j], space_data.phi[i].grad);

                  //calculate scalar product of grad(phi)^T*A*grad(psi)
                  const DataType value3 = N_s * nu * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);

                  // compute outer product of grad(phi) and A*grad(psi)
                  const DataType value4 = N_s * nu * weight * Tiny::dot(loc_A_prod_i, local_prim_dofs[j]);

                  // compute outer product of A*grad(phi) and grad(psi)
                  const DataType value5 = N_s * value2;

                  //compute matrix vector of A and grad(phi)
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value1+value3, local_prim_dofs[j]);
                  //update dot(grad(psi), grad(phi))* A * primal onto vector... this could be performed together with above operation through adding identity...
                  local_vector[i].add_mat_vec_mult(orientation, local_prim_dofs[j], N_s*value1);
                  local_vector[i].axpy(value2 + value4, space_data.phi[j].grad);
                  local_vector[i].axpy(value5, loc_A_prod_j);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into vector
          scatter_vector(local_vector, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }
    }; //class ModBurgersSDAssembler

        /**
   * \brief Modified Burgers Non Linear Defomrmation Operator Assembly class
   *
   * This class is responsible for assembling the scalar and vector-valued modified Burgers operators:
   *
   * \f[\mathbf{N}(v,u,\psi) := \nu \mathbf{L}(u,\psi) + \theta \mathbf{M}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \beta' \mathbf{K'(v,u,\psi)}\f]
   *
   *
   * where
   * - \b L is the diffusive operator with a non-linear <em>deformation tensor<em>:
   *   \f[\mathbf{L}(u, v, \psi) := 2 \int_\Omega (\nabla \psi) : (D(u) + N_s(D(u) f(v) + f(v) D(u)))\f]
   *   where <b>f<b> is a 2nd order tensor, dependent on v
   *         <b>N_s<b> is an additional parameter and
   *         <b>D<b> is the symmetric gradient:
   *   \f[\mathbf{D} := \frac{1}{2}(\nabla u + \nabla u^\top)\f]
   * - \b M is the reactive operator:
   *   \f[\mathbf{M}(u,\psi) := \int_\Omega u\psi\f]
   * - \b K is the convective operator:
   *   \f[\mathbf{K}(v,u,\psi) := \int_\Omega v\cdot \nabla u \psi\f]
   * - <b>K'</b> is the Frechet derivative of the convective operator:
   * TODO: add the frechet derivative part of the deformation
   *   \f[\mathbf{K'}(v,u,\psi) := \int_\Omega \nabla v u \psi\f]
   *
   * \author Maximilian Esser
   */

  template<typename DataType_, typename IndexType_, int dim_>
    class ModBurgersNonLinAssembler
    {
      public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      DataType_ nu;

      //scaling parameter for orientation part of viscosoty
      DataType_ N_s;

      /// scaling parameter for reactive operator \b M
      DataType_ theta;

      /// scaling parameter for convective operator \b K
      DataType_ beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      DataType_ frechet_beta;

      /// default constructor
      ModBurgersNonLinAssembler() :
        nu(DataType_(1)),
        N_s(DataType_(1)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix.
       *
       * \param[in,out] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] function
       * A dim \times dim dimensional vectorfield function which maps to a dim \times dim tensor in the
       * following way:
       * A[i][j](v) = f_{i + j*dim}(v)
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename Function_>
      void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Function_& function,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;


        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (Math::abs(frechet_beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

        // create a function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval(function);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, dim_, dim_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v, orientation;

        //temp vector to save our orientation data as function output
        Tiny::Vector<DataType, dim_ * dim_> orient_vec;

        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        mean_v.format();
        loc_grad_v.format();
        orientation.format();
        orient_vec.format();
        loc_A_prod_i.format();
        loc_A_prod_j.format();



        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);


          // format our local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);


            // evaluate convection function and its gradient (if required)

            loc_v.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // update velocity value
              loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
            }

            if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }

            //calculate the function output
            orient_vec = func_eval.value(loc_v);

            for(int i(0); i < dim_; ++i)
            {
              for(int j(0); j < dim_; ++j)
              {
                orientation[i][j] = orient_vec(i + j* dim_);
              }
            }


            // assemble diffusion matrix?
            if(need_diff)
            {
              // assemble deformation-tensor diffusion and orientation diffusion in one go, careful, as i denotes the test-function and j denotes the trail-functions...
              // this has to be considered in my comprehensions...

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of \psi and A
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);

                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[j].grad, space_data.phi[i].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value1);

                  //add dot(grad_j, grad_i) * A to local_matrix, we can even do this smarter by adding A + Identity to our local_mat...
                  local_matrix[i][j].axpy(N_s*value1, orientation);

                  //calculate scalar product of grad(phi)^T*A*grad(psi) add this to the diag of our local matrix... and of course we can add this to operation above...:
                  const DataType value2 = nu * N_s * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);
                  local_matrix[i][j].add_scalar_main_diag(value2);

                  //calculate vector mult for grad_j
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  //add outer product of grad(phi) * A_prod_psi
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, loc_A_prod_i, N_s * nu * weight);
                  //and the other way round, but transposed
                  local_matrix[i][j].add_outer_product(loc_A_prod_j, space_data.phi[i].grad, N_s * nu * weight);
                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, space_data.phi[i].grad, nu * weight);
                }
              }

            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = frechet_beta * weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(value, loc_grad_v);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

    /**
       * \brief Assembles the Burgers operator into a vector.
       *
       * \param[in,out] vector
       * The vector to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] function
       * Function desribing the space dependent orientation tensor.
       *
       * \param[in] primal
       * The primal vector, usually a solution vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_, typename CubatureFactory_, typename Function_>
      void assemble_vector(
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Function_& function,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& primal,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        //const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector-scatter-axpy (if needed)
        typename VectorType::ScatterAxpy scatter_vector(vector);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // create primal gather-axpy
        typename VectorType::GatherAxpy gather_prim(primal);

        // define our analytic evaluation traits
        typedef Analytic::EvalTraits<DataType, Function_> AnalyticEvalTraits;

        // create a function evaluator
        typename Function_::template Evaluator<AnalyticEvalTraits> func_eval(function);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // local primal vector dofs
        LocalVectorType local_prim_dofs;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> orientation;
        Tiny::Vector<DataType, dim_ * dim_> orient_vec;
        //we will save in this vector our point in the domain cell...
        //i should ask peter if this could be handled smarter...
        typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;

        loc_v.format();
        //loc_grad_v.format();
        orientation.format();
        orient_vec.format();
        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // gather our local primal dofs
          local_prim_dofs.format();
          gather_prim(local_prim_dofs, dof_mapping);

          // format our local vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);



            // evaluate convection function and its gradient (always required)

            loc_v.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // update velocity value
              loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
            }

            /*if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }*/

            //calculate the function output
            orient_vec = func_eval.value(loc_v);

            for(int i(0); i < dim_; ++i)
            {
              for(int j(0); j < dim_; ++j)
              {
                orientation[i][j] = orient_vec(i + j* dim_);
              }
            }

            // assemble diffusion matrix?
            if(need_diff)
            {

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of A and grad(psi)
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // compute outer product of grad(phi) and grad(psi)
                  const DataType value2 = nu * weight * Tiny::dot(local_prim_dofs[j], space_data.phi[i].grad);

                  //calculate scalar product of grad(phi)^T*A*grad(psi)
                  const DataType value3 = N_s * nu * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);

                  // compute outer product of grad(phi) and A*grad(psi)
                  const DataType value4 = N_s * nu * weight * Tiny::dot(loc_A_prod_i, local_prim_dofs[j]);

                  // compute outer product of A*grad(phi) and grad(psi)
                  const DataType value5 = N_s * value2;

                  //compute matrix vector of A and grad(phi)
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value1+value3, local_prim_dofs[j]);
                  //update dot(grad(psi), grad(phi))* A * primal onto vector... this could be performed together with above operation through adding identity...
                  local_vector[i].add_mat_vec_mult(orientation, local_prim_dofs[j], N_s*value1);
                  local_vector[i].axpy(value2 + value4, space_data.phi[j].grad);
                  local_vector[i].axpy(value5, loc_A_prod_j);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into vector
          scatter_vector(local_vector, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }
    }; //class ModBurgersNonLinAssembler


    /**
   * \brief SRI-Burgers Operator Assembly class (Shear-Rate increasing)
   *
   * This class is responsible for assembling the scalar and vector-valued SRI-Burgers operators:
   *
   * \f[\mathbf{N}(v,u,\psi) := \nu \mathbf{L}(u,\psi) + \theta \mathbf{M}(u,\psi) + \beta \mathbf{K}(v,u,\psi) + \beta' \mathbf{K'(v,u,\psi)}\f]
   *
   *
   * where
   * - \b L is the diffusive operator with a modified <em>deformation tensor<em>:
   *   \f[\mathbf{L}(u, v, w) := 2 \int_\Omega (\nabla w) : (D(u) + N_s(D(u)D(v) + D(v)D(u)))\f]
   *   where <b>A<b> is a 2nd order tensor, describing a 2nd momentum matrix of a probability distribution,
   *         <b>N_s<b> is an additional parameter and
   *         <b>D<b> is the symmetric gradient:
   *   \f[\mathbf{D}(u) := \frac{1}{2}(\nabla u + \nabla u^\top)\f]
   * - \b M is the reactive operator:
   *   \f[\mathbf{M}(u,\psi) := \int_\Omega u\psi\f]
   * - \b K is the convective operator:
   *   \f[\mathbf{K}(v,u,\psi) := \int_\Omega v\cdot \nabla u \psi\f]
   *TODO: Add Frecht part of L
   * - <b>K'</b> is the Frechet derivative of the convective operator:
   *   \f[\mathbf{K'}(v,u,\psi) := \int_\Omega \nabla v u \psi\f]
   *
   * \author Maximilian Esser
   */

  template<typename DataType_, typename IndexType_, int dim_>
    class SRIBurgersAssembler
    {
      public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      DataType_ nu;

      //scaling parameter for orientation part of viscosoty
      DataType_ N_s;

      /// scaling parameter for reactive operator \b M
      DataType_ theta;

      /// scaling parameter for convective operator \b K
      DataType_ beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      DataType_ frechet_beta;

      /// default constructor
      SRIBurgersAssembler() :
        nu(DataType_(1)),
        N_s(DataType_(1)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix.
       *
       * \param[in,out] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Space_, typename CubatureFactory_>
      void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;


        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (Math::abs(frechet_beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, dim_, dim_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v, orientation;

        loc_v.format();
        mean_v.format();
        loc_grad_v.format();
        orientation.format();
        loc_A_prod_i.format();
        loc_A_prod_j.format();


        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);


          // format our local matrix and vector
          local_matrix.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }

            loc_grad_v.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // update velocity gradient
              loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
            }
            //TODO: this could be speed up by implementing own operation
            orientation.set_transpose(loc_grad_v);
            orientation.axpy(DataType(1), loc_grad_v);
            for(int i = 0; i < dim_; ++i)
            {
              for(int j = 0; j < dim_; ++j)
              {
                orientation[i][j] = DataType(0.5) * orientation[i][j];
              }
            }


            // assemble diffusion matrix?
            if(need_diff)
            {
              // assemble deformation-tensor diffusion and orientation diffusion in one go, careful, as i denotes the test-function and j denotes the trail-functions...
              // this has to be considered in my comprehensions...

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of \psi and A
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);

                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[j].grad, space_data.phi[i].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value1);

                  //add dot(grad_j, grad_i) * A to local_matrix, we can even do this smarter by adding A + Identity to our local_mat...
                  local_matrix[i][j].axpy(N_s*value1, orientation);

                  //calculate scalar product of grad(phi)^T*A*grad(psi) add this to the diag of our local matrix... and of course we can add this to operation above...:
                  const DataType value2 = nu * N_s * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);
                  local_matrix[i][j].add_scalar_main_diag(value2);

                  //calculate vector mult for grad_j
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  //add outer product of grad(phi) * A_prod_psi
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, loc_A_prod_i, N_s * nu * weight);
                  //and the other way round, but transposed
                  local_matrix[i][j].add_outer_product(loc_A_prod_j, space_data.phi[i].grad, N_s * nu * weight);
                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, space_data.phi[i].grad, nu * weight);
                }
              }

            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = frechet_beta * weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(value, loc_grad_v);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

      /**
       * \brief Assembles the Burgers operator into a vector.
       *
       * \param[in,out] vector
       * The vector to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] primal
       * The primal vector, usually a solution vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] scale
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_, typename CubatureFactory_>
      void assemble_vector(
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& primal,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        //const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector-scatter-axpy (if needed)
        typename VectorType::ScatterAxpy scatter_vector(vector);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // create primal gather-axpy
        typename VectorType::GatherAxpy gather_prim(primal);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // local primal vector dofs
        LocalVectorType local_prim_dofs;

        //a local helper variable to save the matrix_vec multi of A times the basis vector grad(psi) resp. grad(phi)
        Tiny::Vector<DataType, dim_> loc_A_prod_i, loc_A_prod_j;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v, orientation;

        loc_v.format();
        loc_grad_v.format();
        orientation.format();

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // gather our local primal dofs
          local_prim_dofs.format();
          gather_prim(local_prim_dofs, dof_mapping);

          // format our local vector
          local_vector.format();

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }

            loc_grad_v.format();
            for(int i(0); i < num_loc_dofs; ++i)
            {
              // update velocity gradient
              loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
            }
            //TODO: this could be speed up by implementing own operation
            orientation.set_transpose(loc_grad_v);
            orientation.axpy(DataType(1), loc_grad_v);
            for(int i = 0; i < dim_; ++i)
            {
              for(int j = 0; j < dim_; ++j)
              {
                orientation[i][j] = DataType(0.5) * orientation[i][j];
              }
            }


            // assemble diffusion matrix?
            if(need_diff)
            {

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                //calculate the matrix vector product of A and grad(psi)
                loc_A_prod_i.set_mat_vec_mult(orientation, space_data.phi[i].grad);
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // compute outer product of grad(phi) and grad(psi)
                  const DataType value2 = nu * weight * Tiny::dot(local_prim_dofs[j], space_data.phi[i].grad);

                  //calculate scalar product of grad(phi)^T*A*grad(psi)
                  const DataType value3 = N_s * nu * weight * orientation.scalar_product(space_data.phi[j].grad, space_data.phi[i].grad);

                  // compute outer product of grad(phi) and A*grad(psi)
                  const DataType value4 = N_s * nu * weight * Tiny::dot(loc_A_prod_i, local_prim_dofs[j]);

                  // compute outer product of A*grad(phi) and grad(psi)
                  const DataType value5 = N_s * value2;

                  //compute matrix vector of A and grad(phi)
                  loc_A_prod_j.set_mat_vec_mult(orientation, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value1+value3, local_prim_dofs[j]);
                  //update dot(grad(psi), grad(phi))* A * primal onto vector... this could be performed together with above operation through adding identity...
                  local_vector[i].add_mat_vec_mult(orientation, local_prim_dofs[j], N_s*value1);
                  local_vector[i].axpy(value2 + value4, space_data.phi[j].grad);
                  local_vector[i].axpy(value5, loc_A_prod_j);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into vector
          scatter_vector(local_vector, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }
    }; //class SRIBurgersAssembler

    /**
   * \brief Modified Burgers Operator Assembly class
   *
   * This class is responsible for assembling a special variant of the burgers assembler, wich switches between the defo and gradient formulation on the domain itself
   *
   * \author Maximilian Esser
   */

  template<typename DataType_, typename IndexType_, int dim_>
    class ChangeGradToDefoBurgersAssembler
    {
      public:
      /// the datatype we use here
      typedef DataType_ DataType;

      /// scaling parameter for diffusive operator \b L (aka viscosity)
      DataType_ nu;

      /// scaling parameter for reactive operator \b M
      DataType_ theta;

      /// scaling parameter for convective operator \b K
      DataType_ beta;

      /// scaling parameter for Frechet derivative of convective operator <b>K'</b>
      DataType_ frechet_beta;

      /// default constructor
      ChangeGradToDefoBurgersAssembler() :
        nu(DataType_(1)),
        theta(DataType_(0)),
        beta(DataType_(0)),
        frechet_beta(DataType_(0))
      {
      }

      /**
       * \brief Assembles the Burgers operator into a matrix.
       *
       * \param[in,out] matrix
       * The matrix to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] x_threshold
       * The threshold for which the deformation formulation should be used.
       *
       * \param[in] scale
       * A scaling factor for the matrix to be assembled.
       */
      template<typename Space_, typename CubatureFactory_>
      void assemble_matrix(
        LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_>& matrix,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ x_threshold,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(matrix.rows() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(matrix.columns() == space.get_num_dofs(), "invalid matrix dimensions");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;
        typedef LAFEM::SparseMatrixBCSR<Mem::Main, DataType_, IndexType_, dim_, dim_> MatrixType;


        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        const bool need_conv_frechet = (Math::abs(frechet_beta) > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create matrix scatter-axpy
        typename MatrixType::ScatterAxpy scatter_matrix(matrix);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local matrix data
        typedef Tiny::Matrix<DataType, dim_, dim_> MatrixValue;
        typedef Tiny::Matrix<MatrixValue, max_local_dofs, max_local_dofs> LocalMatrixType;
        LocalMatrixType local_matrix;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v, mean_v;

        // our local velocity gradient
        Tiny::Matrix<DataType, dim_, dim_> loc_grad_v;

        loc_v.format();
        mean_v.format();
        loc_grad_v.format();



        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);


          // format our local matrix and vector
          local_matrix.format();

          bool defo = true;
          //we will save the real point_data into this vector
          typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;
          //do a small loop over all element points to check if we need to switch to the gradient formulation
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));
            //transform the reference point onto our real mesh
            trafo_eval.map_point(point_dom, cubature_rule.get_point(point));

            //check wether our x-component is bigger than our given threshold
            if(point_dom[0] >= x_threshold)
            {
               defo = false;
            }
          }

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }


            // assemble diffusion matrix?
            // assemble diffusion matrix?
            if(need_diff && !defo)
            {
              // assemble gradient-tensor diffusion

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            else if(need_diff && defo)
            {
              // assemble deformation-tensor diffusion

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value = nu * weight * Tiny::dot(space_data.phi[j].grad, space_data.phi[i].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);

                  // add outer product of grad(phi) and grad(psi)
                  local_matrix[i][j].add_outer_product(space_data.phi[j].grad, space_data.phi[i].grad, nu * weight);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // assemble convection Frechet?
            if(need_conv_frechet)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = frechet_beta * weight * space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].axpy(value, loc_grad_v);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local matrix
                  local_matrix[i][j].add_scalar_main_diag(value);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into matrix
          scatter_matrix(local_matrix, dof_mapping, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }

    /**
       * \brief Assembles the Burgers operator into a vector.
       *
       * \param[in,out] vector
       * The vector to be assembled.
       *
       * \param[in] convect
       * The transport vector for the convection.
       *
       * \param[in] primal
       * The primal vector, usually a solution vector.
       *
       * \param[in] space
       * The velocity space.
       *
       * \param[in] cubature_factory
       * The cubature factory to be used for integration.
       *
       * \param[in] x_threshold
       * The threshold for which the deformation formulation should be used.
       *
       * \param[in] scale
       * A scaling factor the the vector to be assembled.
       */
      template<typename Space_, typename CubatureFactory_>
      void assemble_vector(
        LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& vector,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& convect,
        const LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_>& primal,
        const Space_& space,
        const CubatureFactory_& cubature_factory,
        const DataType_ x_threshold,
        const DataType_ scale = DataType_(1)
        ) const
      {
        // validate matrix and vector dimensions
        XASSERTM(vector.size() == space.get_num_dofs(), "invalid vector size");
        XASSERTM(convect.size() == space.get_num_dofs(), "invalid vector size");

        typedef LAFEM::DenseVectorBlocked<Mem::Main, DataType_, IndexType_, dim_> VectorType;

        // first of all, let's see what we have to assemble
        const bool need_diff = (Math::abs(nu) > DataType(0));
        const bool need_conv = (Math::abs(beta) > DataType(0));
        //const bool need_conv_frechet = (frechet_beta > DataType(0));
        const bool need_reac = (Math::abs(theta) > DataType(0));

        // define our assembly traits
        typedef AsmTraits1<DataType_, Space_, TrafoTags::jac_det, SpaceTags::value|SpaceTags::grad> AsmTraits;

        // fetch our trafo
        const typename AsmTraits::TrafoType& trafo = space.get_trafo();

        // create a trafo evaluator
        typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

        // create a space evaluator and evaluation data
        typename AsmTraits::SpaceEvaluator space_eval(space);

        // create a dof-mapping
        typename AsmTraits::DofMapping dof_mapping(space);

        // create trafo evaluation data
        typename AsmTraits::TrafoEvalData trafo_data;

        // create space evaluation data
        typename AsmTraits::SpaceEvalData space_data;

        // create cubature rule
        typename AsmTraits::CubatureRuleType cubature_rule(Cubature::ctor_factory, cubature_factory);

        // create vector-scatter-axpy (if needed)
        typename VectorType::ScatterAxpy scatter_vector(vector);

        // create convection gather-axpy
        typename VectorType::GatherAxpy gather_conv(convect);

        // create primal gather-axpy
        typename VectorType::GatherAxpy gather_prim(primal);

        // get maximum number of local dofs
        static constexpr int max_local_dofs = AsmTraits::max_local_test_dofs;

        // create local vector data
        typedef Tiny::Vector<DataType, dim_> VectorValue;
        typedef Tiny::Vector<VectorValue, max_local_dofs> LocalVectorType;
        LocalVectorType local_vector;

        // local convection field dofs
        LocalVectorType local_conv_dofs;

        // local primal vector dofs
        LocalVectorType local_prim_dofs;

        // our local velocity value
        Tiny::Vector<DataType, dim_> loc_v;

        // our local velocity gradient
        //Tiny::Matrix<DataType, dim_, dim_> loc_grad_v;

        loc_v.format();
        //loc_grad_v.format();

        // loop over all cells of the mesh
        for(typename AsmTraits::CellIterator cell(trafo_eval.begin()); cell != trafo_eval.end(); ++cell)
        {
          // prepare trafo evaluator
          trafo_eval.prepare(cell);

          // prepare space evaluator
          space_eval.prepare(trafo_eval);

          // initialize dof-mapping
          dof_mapping.prepare(cell);

          // fetch number of local dofs
          const int num_loc_dofs = space_eval.get_num_local_dofs();

          // gather our local convection dofs
          local_conv_dofs.format();
          gather_conv(local_conv_dofs, dof_mapping);

          // gather our local primal dofs
          local_prim_dofs.format();
          gather_prim(local_prim_dofs, dof_mapping);

          // format our local vector
          local_vector.format();

          bool defo = true;
          //we will save the real point_data into this vector
          typename AsmTraits::TrafoEvaluator::ImagePointType point_dom;
          //do a small loop over all element points to check if we need to switch to the gradient formulation
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));
            //transform the reference point onto our real mesh
            trafo_eval.map_point(point_dom, cubature_rule.get_point(point));

            //check wether our x-component is bigger than our given threshold
            if(point_dom[0] >= x_threshold)
            {
               defo = false;
            }
          }

          // loop over all quadrature points and integrate
          for(int point(0); point < cubature_rule.get_num_points(); ++point)
          {
            // compute trafo data
            trafo_eval(trafo_data, cubature_rule.get_point(point));

            // compute basis function data
            space_eval(space_data, trafo_data);

            // pre-compute cubature weight
            const DataType weight = trafo_data.jac_det * cubature_rule.get_weight(point);

            // evaluate convection function and its gradient (if required)
            if(need_conv)
            {
              loc_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity value
                loc_v.axpy(space_data.phi[i].value, local_conv_dofs[i]);
              }
            }
            /*if(need_conv_frechet)
            {
              loc_grad_v.format();
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // update velocity gradient
                loc_grad_v.add_outer_product(local_conv_dofs[i], space_data.phi[i].grad);
              }
            }*/

            // assemble diffusion matrix?
            if(need_diff && !defo)
            {
              // assemble gradient-tensor diffusion

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }
            else if(need_diff && defo)
            {
              // assemble deformation-tensor diffusion

              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute inner product of grad(phi) and grad(psi)
                  const DataType value1 = nu * weight * Tiny::dot(space_data.phi[i].grad, space_data.phi[j].grad);

                  // compute outer product of grad(phi) and grad(psi)
                  const DataType value2 = nu * weight * Tiny::dot(local_prim_dofs[j], space_data.phi[i].grad);

                  // update local vector
                  local_vector[i].axpy(value1, local_prim_dofs[j]);
                  local_vector[i].axpy(value2, space_data.phi[j].grad);
                }
              }
            }

            // assemble convection?
            if(need_conv)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = beta * weight * space_data.phi[i].value * Tiny::dot(loc_v, space_data.phi[j].grad);

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // assemble reaction?
            if(need_reac)
            {
              // test function loop
              for(int i(0); i < num_loc_dofs; ++i)
              {
                // trial function loop
                for(int j(0); j < num_loc_dofs; ++j)
                {
                  // compute scalar value
                  const DataType value = theta * weight *  space_data.phi[i].value * space_data.phi[j].value;

                  // update local vector
                  local_vector[i].axpy(value, local_prim_dofs[j]);
                }
              }
            }

            // continue with next cubature point
          }

          // scatter into vector
          scatter_vector(local_vector, dof_mapping, scale);

          // finish dof mapping
          dof_mapping.finish();

          // finish evaluators
          space_eval.finish();
          trafo_eval.finish();
        }
      }
    }; //class ChangeGradToDefoBurgersAssembler
  }
}