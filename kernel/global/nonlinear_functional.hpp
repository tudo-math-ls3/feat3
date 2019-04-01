// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2019 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GLOBAL_NONLINEAR_FUNCTIONAL_HPP
#define KERNEL_GLOBAL_NONLINEAR_FUNCTIONAL_HPP 1

#include <kernel/global/gate.hpp>
#include <kernel/global/vector.hpp>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Global NonlinearFunctional wrapper class template
     *
     * \tparam LocalNonlinearFunctional_
     * The class of the (patch-) local nonlinear functional
     *
     * \author Jordi Paul
     */
    template<typename LocalNonlinearFunctional_, typename RowMirrorType_, typename ColMirrorType_>
    class NonlinearFunctional
    {
      public:
        /// Memory architecture of the local functional
        typedef typename LocalNonlinearFunctional_::MemType MemType;
        /// Floating point data type of the local functional
        typedef typename LocalNonlinearFunctional_::DataType DataType;
        /// Index type of the local functional
        typedef typename LocalNonlinearFunctional_::IndexType IndexType;

        /// The local functionals left-vector type
        typedef typename LocalNonlinearFunctional_::VectorTypeL LocalVectorTypeL;
        /// The local functionals right-vector type
        typedef typename LocalNonlinearFunctional_::VectorTypeR LocalVectorTypeR;
        /// The local functionals filter type
        typedef typename LocalNonlinearFunctional_::FilterType LocalFilterType;

        /// The associated global left-vector type
        typedef Vector<LocalVectorTypeL, RowMirrorType_> VectorTypeL;
        /// The associated global right-vector type
        typedef Vector<LocalVectorTypeR, ColMirrorType_> VectorTypeR;
        /// The associated global filter type
        typedef Filter<LocalFilterType, RowMirrorType_> FilterType;

        /// Global Gate for left-vectors
        typedef Gate<LocalVectorTypeL, RowMirrorType_> GateRowType;
        /// Global Gate for right-vectors
        typedef Gate<LocalVectorTypeR, ColMirrorType_> GateColType;

        /// The a gradient vector is the output of the operators compute_grad() function and thus a left-vector
        typedef VectorTypeL GradientType;

        /// The global nonlinear functional's Blockheight ist the same as the local nonlinear functional's
        static constexpr int BlockHeight = LocalNonlinearFunctional_::BlockHeight;
        /// The global nonlinear functional's Blockwidth ist the same as the local nonlinear functional's
        static constexpr int BlockWidth = LocalNonlinearFunctional_::BlockWidth;

      protected:
        /// Gate for syncing row vectors
        GateRowType* _row_gate;
        /// Gate for syncing column vectors
        GateColType* _col_gate;
        /// The underlying local nonlinear functional
        LocalNonlinearFunctional_& _nonlinear_functional;

      public:
        /**
         * \brief Constructor
         *
         * \param[in] row_gate
         * Gate for rows
         *
         * \param[in] col_gate
         * Gate for columns
         *
         * \param[in] nonlinear_functional
         * The (patch-)local nonlinear functional
         *
         */
        explicit NonlinearFunctional(GateRowType* row_gate, GateColType* col_gate,
        LocalNonlinearFunctional_& nonlinear_functional) :
          _row_gate(row_gate),
          _col_gate(col_gate),
          _nonlinear_functional(nonlinear_functional)
          {
          }

        /// Explicitly delete default constructor
        NonlinearFunctional() = delete;
        /// Explicitly delete copy constructor
        NonlinearFunctional(const NonlinearFunctional&) = delete;

        /// \brief Empty virtual destructor
        virtual ~NonlinearFunctional()
        {
        }

        /**
         * \brief Gets the local nonlinear functional
         *
         * \returns A reference to the underlying local nonlinear functional
         */
        LocalNonlinearFunctional_& local()
        {
          return _nonlinear_functional;
        }

        /**
         * \brief Gets the local nonlinear functional, const version
         *
         * \returns A const reference to the underlying local nonlinear functional
         */
        const LocalNonlinearFunctional_& local() const
        {
          return _nonlinear_functional;
        }

        /**
         * \brief Returns the functional's communicator
         */
        const Dist::Comm* get_comm() const
        {
          return _row_gate != nullptr ? _row_gate->get_comm() : nullptr;
        }

        /**
         * \brief Get the number of times compute_func() was called
         *
         * \returns The number of functional evaluations
         */
        Index get_num_func_evals() const
        {
          return _nonlinear_functional.get_num_func_evals();
        }

        /**
         * \brief Get the number of times compute_grad() was called
         *
         * \returns The number of gradient evaluations
         */
        Index get_num_grad_evals() const
        {
          return _nonlinear_functional.get_num_grad_evals();
        }

        /**
         * \brief Get the number of times compute_hess() was called
         *
         * \returns The number of Hessian evaluations
         */
        Index get_num_hess_evals() const
        {
          return _nonlinear_functional.get_num_hess_evals();
        }

        /**
         * \brief Resets the evaluation counters
         *
         */
        void reset_num_evals()
        {
          _nonlinear_functional.reset_num_evals();
        }

        /**
         * \brief Creates an empty left-vector of the correct size
         *
         * \returns A new left-vector
         */
        VectorTypeL create_vector_l() const
        {
          return VectorTypeL(_row_gate, _nonlinear_functional.create_vector_l());
        }

        /**
         * \brief Creates an empty right-vector of the correct size
         *
         * \returns A new right-vector
         */
        VectorTypeR create_vector_r() const
        {
          return VectorTypeR(_col_gate, _nonlinear_functional.create_vector_r());
        }

        /**
         * \brief Functionality that cannot be done in the constructor
         */
        void init()
        {
          _nonlinear_functional.init_pre_sync();
          // Synchronise the necessary parts
          while(!_nonlinear_functional.sync_scalars.empty())
          {
            auto it = _nonlinear_functional.sync_scalars.begin();
            if(_col_gate != nullptr)
            {
              *((*it).second) = _col_gate->sum(*((*it).second));
            }
            _nonlinear_functional.sync_scalars.erase(it);
          }

          if(_col_gate != nullptr)
          {
            for(auto* it:_nonlinear_functional.sync_vecs)
            {
              _col_gate->sync_0(*it);
            }
          }
          _nonlinear_functional.init_post_sync();
        }

        /**
         * \brief Prepares the operator for evaluation by setting the current state
         *
         * \param[in] vec_state
         * The current state
         *
         * \param[in,out] filter
         * The filter that might be modified by the local functional's prepare
         *
         */
        void prepare(const VectorTypeR& vec_state, FilterType& filter)
        {

          // Prepare the patch local nonlinear functional, computing everything that does need synchronising
          _nonlinear_functional.prepare_pre_sync(vec_state.local(), filter.local());

          // Synchronise the necessary parts
          while(!_nonlinear_functional.sync_scalars.empty())
          {
            auto it = _nonlinear_functional.sync_scalars.begin();
            if(_col_gate != nullptr)
            {
              *((*it).second) = _col_gate->sum(*((*it).second));
            }
            _nonlinear_functional.sync_scalars.erase(it);
          }

          if(_col_gate != nullptr)
          {
            for(auto* it:_nonlinear_functional.sync_vecs)
            {
              _col_gate->sync_0(*it);
            }
          }

          // Prepare the rest of the functional
          _nonlinear_functional.prepare_post_sync(vec_state.local(), filter.local());

          // Sync the filter vector in the SlipFilter
          if(_col_gate != nullptr )
          {
            // For all slip filters...
            for(auto& it : filter.local().template at<0>())
            {
              // get the filter vector
              auto& slip_filter_vector = it.second.get_filter_vector();

              if(slip_filter_vector.used_elements() > 0)
              {
                // Temporary DenseVector for syncing
                LocalVectorTypeR tmp(slip_filter_vector.size(), DataType(0));

                auto* tmp_elements = tmp.template elements<LAFEM::Perspective::native>();
                auto* sfv_elements = slip_filter_vector.template elements<LAFEM::Perspective::native>();

                // Copy sparse filter vector contents to DenseVector
                for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
                {
                  Index idense(slip_filter_vector.indices()[isparse]);
                  tmp_elements[idense] = sfv_elements[isparse];
                }

                _col_gate->sync_0(tmp);
                // Copy sparse filter vector contents to DenseVector
                for(Index isparse(0); isparse < slip_filter_vector.used_elements(); ++isparse)
                {
                  Index idense(slip_filter_vector.indices()[isparse]);
                  tmp_elements[idense].normalise();
                  sfv_elements[isparse] = tmp_elements[idense];

                }
              }

              else
              {
                // Temporary DenseVector for syncing
                LocalVectorTypeR tmp(slip_filter_vector.size(), DataType(0));
                _col_gate->sync_0(tmp);
              }
            }
          } // col_gate
        } // prepare()

        /**
         * \brief Gets the number of columns
         *
         * \warning In parallel, this requires communication and is very expensive, so use sparingly!
         *
         * \returns The number of colums
         */
        Index columns()
        {
          // Compute total number of rows and columns
          auto vec_r = create_vector_r();
          vec_r.format(DataType(1));

          return Index(vec_r.norm2sqr());
        }

        /**
         * \brief Gets the number of rows
         *
         * \warning In parallel, this requires communication and is very expensive, so use sparingly!
         *
         * \returns The number of colums
         */
        Index rows()
        {
          // Compute total number of rows and rows
          auto vec_l = create_vector_l();
          vec_l.format(DataType(1));

          return Index(vec_l.norm2sqr());
        }

        /**
         * \brief Computes the functional's value and gradient at the current state
         *
         * \param[out] fval
         * The functional value
         *
         * \param[out] grad
         * The vector receiving the synced gradient vector
         *
         */
        void eval_fval_grad(DataType& fval, VectorTypeL& grad)
        {
          // As the penalty term is quadratic, we have to compute the local functional values separately for summing
          // them up over all ranks
          const bool add_penalty_fval(false);
          _nonlinear_functional.eval_fval_grad(fval, grad.local(), add_penalty_fval);

          // Sum up over all patches
          if(_row_gate != nullptr)
            fval = _row_gate->sum(fval);

          // Add the penalty term
          if(get_penalty_param() > DataType(0))
          {
            DataType constraint(get_constraint());
            fval += get_penalty_param()*DataType(0.5)*Math::sqr(constraint);
          }

          // Synchronise the gradient vector
          grad.sync_0();
        }

        /**
         * \brief Get the penalty parameter from the local nonlinear functional
         *
         * \returns A copy of the local nonlinear functional's penalty parameter
         */
        DataType get_penalty_param() const
        {
          return _nonlinear_functional.get_penalty_param();
        }

        /**
         * \brief Sets the local functional's penalty parameter
         */
        void set_penalty_param(const DataType fac)
        {
          _nonlinear_functional.set_penalty_param(fac);
        }

        /**
         * \brief Get the constraint from the local nonlinear functional
         *
         * \returns A copy of the global nonlinear functional's constraint
         */
        DataType get_constraint()
        {
          return _nonlinear_functional.get_constraint();
        }

        /**
         * \brief Computes the constraint and returns it
         *
         * \returns The global nonlinear functional's constraint
         */
        DataType compute_constraint()
        {
          return _nonlinear_functional.compute_constraint();
        }

    };
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_NONLINEAR_FUNCTIONAL
