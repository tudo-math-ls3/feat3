// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/lafem/null_matrix.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/tuple_matrix.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/solver/amavanka_base.hpp>
#include <kernel/util/stop_watch.hpp>

namespace FEAT
{
  namespace Solver
  {
    /**
     * \brief Additive Macro-wise Matrix-based Vanka preconditioner/smoother
     *
     * This class implements an additive macro-wise Vanka smoother, which stores its
     * pre-computed operator as a sparse matrix, so that each application of the Vanka
     * smoother consists of only one sparse matrix-vector multiplication.
     *
     * This class supports a whole batch of different matrix types:
     * - LAFEM::SparseMatrixCSR
     * - LAFEM::SparseMatrixBCSR
     * - LAFEM::SaddlePointMatrix with sub-blocks of type
     *   - LAFEM::SparseMatrixBCSR
     * - LAFEM::TupleMatrix with sub-blocks of type
     *   - LAFEM::SparseMatrixCSR
     *   - LAFEM::SparseMatrixBCSR
     *
     * \attention
     * In contrast to the Solver::Vanka class, one must specify the macros for the Vanka smoother for each
     * finite element space which is used to discretize the individual blocks of the solution space by
     * using the #push_macro_dofs() function. The only exception is for the LAFEM::SaddlePointMatrix,
     * as for this matrix type, there exists an automatic macro deduction algorithm, which will compute
     * the element-wise macro graphs based on the sparsity patterns of the B and D matrices, resp.
     *
     * <b>Example:</b> Assume that your solution space is the usual velocity-pressure space pair used in
     * a Stokes system, then you need to add the dof-mapping graphs of both the velocity and the pressure
     * space to the smoother in the following way:
     * \code{.cpp}
     * auto vanka = Solver::new_amavanka(system_matrix, system_filter);
     * vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_velocity));
     * vanka->push_macro_dofs(Space::DofMappingRenderer::render(space_pressure));
     * \endcode
     * If your solution space has other components (such as e.g. temperature, density or stress), then you
     * have to add the rendered dof-mappings of those spaces in the correct order as well.
     *
     * \note
     * In the case of a LAFEM::SaddlePointMatrix, the smoother implemented in this class is mathematically
     * equivalent to a Solver::Vanka smoother of type Solver::VankaType::block_full_add.
     *
     * \author Peter Zajac
     */
    template<typename Matrix_, typename Filter_>
    class AmaVanka :
      public Solver::SolverBase<typename Matrix_::VectorTypeL>
    {
    public:
      /// our base-class
      typedef Solver::SolverBase<typename Matrix_::VectorTypeL> BaseClass;

      /// our data type
      typedef typename Matrix_::DataType DataType;
      /// our index type
      typedef typename Matrix_::IndexType IndexType;
      /// our vector type
      typedef typename Matrix_::VectorTypeL VectorType;

    protected:
      /// the type of our Vanka matrix
      typedef typename Intern::AmaVankaMatrixHelper<Matrix_>::VankaMatrix VankaMatrixType;

      /// the system matrix
      const Matrix_& _matrix;
      /// the system filter
      const Filter_& _filter;
      /// the Vanka preconditioner matrix
      VankaMatrixType _vanka;
      /// deduct macro dofs automatically?
      bool _auto_macros;
      /// skip singular macros?
      bool _skip_singular;
      /// the DOF-macro graphs
      std::vector<Adjacency::Graph> _macro_dofs, _dof_macros;
      /// the macro mask
      std::vector<int> _macro_mask;
      /// number of steps
      Index _num_steps;
      /// damping parameter
      DataType _omega;
      /// temporary vectors
      VectorType _vec_c, _vec_d;

      // stop watch for symbolic factorization
      StopWatch watch_init_symbolic;
      // stop watch for numeric factorization
      StopWatch watch_init_numeric;
      // stop watch for apply time
      StopWatch watch_apply;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] matrix
       * The saddle-point system matrix.
       *
       * \param[in] filter
       * The system filter.
       *
       * \param[in] omega
       * The damping parameter.
       *
       * \param[in] num_steps
       * The number of smoothing steps to be performed.
       */
      explicit AmaVanka(const Matrix_& matrix, const Filter_& filter,
        const DataType omega = DataType(1), const Index num_steps = Index(1)) :
        _matrix(matrix),
        _filter(filter),
        _vanka(),
        _auto_macros(true),
        _skip_singular(false),
        _num_steps(num_steps),
        _omega(omega)
      {
      }

      /**
       * \brief Clears the macro dofs graphs.
       */
      void clear_macro_dofs()
      {
        this->_macro_dofs.clear();
        _auto_macros = true;
      }

      /**
       * \brief Pushes the dofs-at-macro graph of the next block.
       *
       * \param[in] dofs
       * The dofs-at-macro graph of the block.
       */
      void push_macro_dofs(Adjacency::Graph&& dofs)
      {
        _auto_macros = false;

        // make sure that we do not push more graphs than we have blocks
        if(int(_macro_dofs.size()) >= Intern::AmaVankaMatrixHelper<Matrix_>::num_blocks)
          XABORTM("all macro-dofs graphs have already been added");

        // make sure that the number of macros matches our previous graph
        if(!_macro_dofs.empty() && (_macro_dofs.back().get_num_nodes_domain() != dofs.get_num_nodes_domain()))
          XABORTM("macro count mismatch");

        // push graph into the list
        _macro_dofs.emplace_back(std::forward<Adjacency::Graph>(dofs));

        // sort the dof indices
        _macro_dofs.back().sort_indices();
      }

      /**
       * \brief Sets the number of smoothing steps.
       *
       * \param[in] num_steps
       * The number of smoothing steps to be performed.
       */
      void set_num_steps(Index num_steps)
      {
        XASSERT(num_steps > Index(0));
        this->_num_steps = num_steps;
      }

      /**
       * \brief Sets the damping parameter omega.
       *
       * \param[in] omega
       * The damping parameter.
       */
      void set_omega(DataType omega)
      {
        XASSERT(omega > DataType(0));
        this->_omega = omega;
      }

      /**
       * \brief Sets whether singular macros are to be skipped.
       *
       * \param[in] skip_sing
       * Specifies whether singular macros are to be skipped.
       */
      void set_skip_singular(bool skip_sing)
      {
        this->_skip_singular = skip_sing;
      }

      /**
       * \brief Returns the total number of bytes currently allocated in this object.
       */
      std::size_t bytes() const
      {
        std::size_t s = _vanka.bytes();
        for(const auto& g : _macro_dofs)
          s += sizeof(Index) * std::size_t(g.get_num_nodes_domain() + g.get_num_indices());
        for(const auto& g : _dof_macros)
          s += sizeof(Index) * std::size_t(g.get_num_nodes_domain() + g.get_num_indices());
        s += _vec_c.bytes();
        s += _vec_d.bytes();
        return s;
      }

      /**
       * \brief Returns the total data size used by the AmaVanka smoother.
       *
       * \returns
       * The total data size, i.e. the total number of floating point values used in the factorization.
       */
      std::size_t data_size() const
      {
        return std::size_t(_vanka.template used_elements<LAFEM::Perspective::pod>());
      }

      /**
       * \brief Resets the internal stop watches for time measurement.
       */
      void reset_timings()
      {
        watch_init_symbolic.reset();
        watch_init_numeric.reset();
        watch_apply.reset();
      }

      /**
       * \brief Returns the total accumulated time for symbolic initialization.
       */
      double time_init_symbolic() const
      {
        return watch_init_symbolic.elapsed();
      }

      /**
       * \brief Returns the total accumulated time for numeric initialization.
       */
      double time_init_numeric() const
      {
        return watch_init_numeric.elapsed();
      }

      /**
       * \brief Returns the total accumulated time for the solver application.
       */
      double time_apply() const
      {
        return watch_apply.elapsed();
      }

      /// Returns the name of the solver.
      virtual String name() const override
      {
        return "AmaVanka";
      }

      /// Performs symbolic factorization
      virtual void init_symbolic() override
      {
        watch_init_symbolic.start();

        BaseClass::init_symbolic();

        // we also need two vectors if we have to perform multiple steps
        if(this->_num_steps > Index(1))
        {
          this->_vec_c = this->_matrix.create_vector_l();
          this->_vec_d = this->_matrix.create_vector_l();
        }

        // automatically deduct macros?
        if(this->_auto_macros)
        {
          // try to deduct macros by pressure matrices in SaddlePointMatrix
          if(!Intern::AmaVankaCore::deduct_macro_dofs(this->_matrix, this->_macro_dofs))
            XABORTM("Cannot auto-deduct macros for this matrix type");
        }

        // make sure we have one macro-dof graph for each matrix block
        XASSERTM(int(_macro_dofs.size()) == Intern::AmaVankaMatrixHelper<Matrix_>::num_blocks,
           "invalid number of macro-dof graphs; did you push all of them?");

        // compute dof-macro graphs by transposing
        this->_dof_macros.resize(this->_macro_dofs.size());
        for(std::size_t i(0); i < this->_macro_dofs.size(); ++i)
        {
          // ensure that we have the same number of macros in all graphs
          XASSERT(this->_macro_dofs.at(i).get_num_nodes_domain() == this->_macro_dofs.front().get_num_nodes_domain());

          // transpose macro-dofs graph
          this->_dof_macros.at(i) = Adjacency::Graph(Adjacency::RenderType::transpose, this->_macro_dofs.at(i));
        }

        // allocate macro skip mask?
        if(this->_skip_singular)
          this->_macro_mask.resize(this->_macro_dofs.front().get_num_nodes_domain(), 0);

        Solver::Intern::AmaVankaCore::alloc(this->_vanka, this->_dof_macros, this->_macro_dofs, Index(0), Index(0));

        watch_init_symbolic.stop();
      }

      /// Releases the symbolic factorization data
      virtual void done_symbolic() override
      {
        this->_vanka.clear();
        this->_macro_mask.clear();
        this->_dof_macros.clear();
        if(this->_auto_macros)
          this->_macro_dofs.clear();
        if(this->_num_steps > Index(1))
        {
          this->_vec_d.clear();
          this->_vec_c.clear();
        }
        BaseClass::done_symbolic();
      }

      /// Performs numeric factorization
      virtual void init_numeric() override
      {
        const DataType eps = Math::eps<DataType>();

        watch_init_numeric.start();
        BaseClass::init_numeric();

        // get maximum macro size
        const Index num_macros = Index(this->_macro_dofs.front().get_num_nodes_domain());
        const Index stride = Intern::AmaVankaCore::calc_stride(this->_vanka, this->_macro_dofs);

        // allocate arrays for local matrix
        std::vector<DataType> vec_local(stride*stride, DataType(0)), vec_local_t(stride*stride, DataType(0));
        std::vector<Index> vec_pivot(stride);
        DataType* local = vec_local.data();
        DataType* local_t = vec_local_t.data();
        Index* pivot = vec_pivot.data();

        this->_vanka.format();

        // loop over all macros
        for(Index imacro(0); imacro < num_macros; ++imacro)
        {
          // gather local matrix
          const std::pair<Index,Index> nrc = Intern::AmaVankaCore::gather(this->_matrix,
            local, stride, imacro, this->_macro_dofs, Index(0), Index(0), Index(0), Index(0));

          // make sure we have gathered a square matrix
          XASSERTM(nrc.first == nrc.second, "local matrix is not square");

          // do we check for singular macros?
          if(this->_skip_singular)
          {
            // the approach used for checking the regularity of the local matrix is to check whether
            //
            //     || I - A*A^{-1} ||_F^2 < eps
            //
            // we could try to analyse the pivots returned by invert_matrix function instead, but
            // unfortunately this approach sometimes leads to false positives

            // make a backup if checking for singularity
            for(Index i(0); i < nrc.first; ++i)
              for(Index j(0); j < nrc.second; ++j)
                local_t[i*stride+j] = local[i*stride+j];

            // invert local matrix
            Math::invert_matrix(nrc.first, stride, local, pivot);

            // compute (squared) Frobenius norm of (I - A*A^{-1})
            DataType norm = DataType(0);
            for(Index i(0); i < nrc.first; ++i)
            {
              for(Index j(0); j < nrc.first; ++j)
              {
                DataType xij = DataType(i == j ? 1 : 0);
                for(Index k(0); k < nrc.first; ++k)
                  xij -= local_t[i*stride+k] * local[k*stride+j]; // A_ik * (A^{-1})_kj
                norm += xij * xij;
              }
            }

            // is the matrix block singular?
            // Note: we check for !(norm < eps) instead of (norm >= eps),
            // because the latter one evaluates to false if norm is NaN,
            // which would result in a false negative
            const bool singular = !(norm < eps);

            // set macro regularity mask
            this->_macro_mask[imacro] = (singular ? 0 : 1);

            // scatter local matrix
            if(!singular)
            {
              Intern::AmaVankaCore::scatter_add(this->_vanka, local, stride, imacro, this->_macro_dofs,
                Index(0), Index(0), Index(0), Index(0));
            }
          }
          else // no singularity check
          {
            // invert local matrix
            Math::invert_matrix(nrc.first, stride, local, pivot);

            // scatter local matrix
            Intern::AmaVankaCore::scatter_add(this->_vanka, local, stride, imacro, this->_macro_dofs,
              Index(0), Index(0), Index(0), Index(0));
          }

          // reformat local matrix
          for(Index i(0); i < nrc.first; ++i)
            for(Index j(0); j < nrc.second; ++j)
              local[i*stride+j] = DataType(0);
        }

        // scale rows of Vanka matrix
        Solver::Intern::AmaVankaCore::scale_rows(this->_vanka, this->_omega, this->_dof_macros, this->_macro_mask, Index(0), Index(0));

        watch_init_numeric.stop();
      }

      /// applies the preconditioner
      virtual Status apply(VectorType& vec_x, const VectorType& vec_b) override
      {
        watch_apply.start();

        // first step
        this->_vanka.apply(vec_x, vec_b);
        this->_filter.filter_cor(vec_x);

        // steps 2, ..., n   (if any)
        for(Index step(1); step < _num_steps; ++step)
        {
          // compute defect
          this->_matrix.apply(this->_vec_d, vec_x, vec_b, -DataType(1));
          // filter defect
          this->_filter.filter_def(this->_vec_d);
          // apply Vanka matrix
          this->_vanka.apply(this->_vec_c, this->_vec_d);
          // filter correct
          this->_filter.filter_cor(this->_vec_c);
          // update solution
          vec_x.axpy(this->_vec_c);
        }

        watch_apply.stop();

        return Status::success;
      }

      bool compare(const AmaVanka* other) const
      {
        return Intern::AmaVankaMatrixHelper<VankaMatrixType>::compare(this->_vanka, other->_vanka);
      }
    }; // class AmaVanka

    /**
     * \brief Creates a new AmaVanka smoother object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter
     *
     * \param[in] omega
     * The damping parameter.
     *
     * \param[in] num_steps
     * The number of Vanka iterations to be performed.
     *
     * \returns
     * A shared pointer to a new AmaVanka object.
     */
    template<typename Matrix_, typename Filter_>
    std::shared_ptr<AmaVanka<Matrix_, Filter_>> new_amavanka(const Matrix_& matrix, const Filter_& filter,
      typename Matrix_::DataType omega = typename Matrix_::DataType(1), Index num_steps = Index(1))
    {
      return std::make_shared<AmaVanka<Matrix_, Filter_>>(matrix, filter, omega, num_steps);
    }
  } // namespace Solver
} // namespace FEAT
