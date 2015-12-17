/**
 * \file
 * \brief FEAST milestone 2 ScaRC functor implementations
 * \author Markus Geveler
 * \date 2014
 *
 * See class documentation.
 *
 */

#pragma once
#ifndef SCARC_GUARD_SCARC_FUNCTOR_HPP
#define SCARC_GUARD_SCARC_FUNCTOR_HPP 1

#include<kernel/base_header.hpp>
#include<kernel/scarc/scarc_error.hpp>
#include<kernel/scarc/scarc_data.hpp>
#include<kernel/foundation/global_defect.hpp>

namespace FEAST
{
  namespace ScaRC
  {
    /**
     * \brief Base class for all ScaRC functors
     *
     * ScaRC functors store references to variable data containers and values for constants.
     * The specific ScaRCFunctor subclass specifies the operation by overwriting its execute() member function.
     * Hence, no function pointers and bindings are required for the functor library implementation.
     *
     * \author Markus Geveler
     */
    template<typename DataType_,
             typename MemTag_,
             typename VectorType_,
             typename VectorMirrorType_,
             typename MatrixType_,
             typename PreconContType_,
             typename FilterType_,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index
             >
    class ScaRCFunctorBase
    {
      public:

        /**
         * \brief executes the functor on the stored arguments
         *
         */
        virtual void execute() = 0;

        /**
         * \brief applies the functor on specific arguments
         *
         */
        virtual void apply(VectorType_& store_to, VectorType_& apply_to, VectorType_& apply_rhs) = 0;

        ///returns a string that describes the functor
        virtual const std::string type_name() = 0;

        ///CTOR
        ScaRCFunctorBase(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                     MemTag_,
                                                                     VectorType_,
                                                                     VectorMirrorType_,
                                                                     MatrixType_,
                                                                     PreconContType_,
                                                                     FilterType_,
                                                                     StorageType_,
                                                                     IT_>& data) :
          _precon(nullptr),
          _data(&data),
          _eps(1e-8),
          _max_iters(1000),
          _conv_check(true),
          _norm_0(std::numeric_limits<DataType_>::max()),
          _norm(std::numeric_limits<DataType_>::max()),
          _used_iters(0)
        {
        }

        ///CTOR
        ScaRCFunctorBase(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                     MemTag_,
                                                                     VectorType_,
                                                                     VectorMirrorType_,
                                                                     MatrixType_,
                                                                     PreconContType_,
                                                                     FilterType_,
                                                                     StorageType_,
                                                                     IT_>& data,
                         std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                          MemTag_,
                                                          VectorType_,
                                                          VectorMirrorType_,
                                                          MatrixType_,
                                                          PreconContType_,
                                                          FilterType_,
                                                          StorageType_,
                                                          IT_> >& precon) :
          _precon(precon),
          _data(&data),
          _eps(1e-8),
          _max_iters(1000),
          _conv_check(true),
          _norm_0(std::numeric_limits<DataType_>::max()),
          _norm(std::numeric_limits<DataType_>::max()),
          _used_iters(0)
        {
        }

        ///copy-CTOR
        ScaRCFunctorBase(const ScaRCFunctorBase& other)
        {
          this->_data = other._data;
          this->_precon = other.precon;
          this->_eps = other._eps;
          this->_max_iters = other._max_iters;
          this->_conv_check = other._conv_check;
          this->_norm_0 = other._norm_0;
          this->_norm = other._norm;
          this->_used_iters = other._used_iters;
        }

        ///virtual DTOR
        virtual ~ScaRCFunctorBase()
        {
        }

        /**
         * \brief substitutes a preconditioner by a functor
         *
         * Preconditioner functors can overwrite this function. With this mechanism, solver layers can be created
         * independently and plugged into each other dynamically.
         *
         * \param[in] precon
         * the substitute
         *
         */
        virtual void reset_preconditioner(std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                                           MemTag_,
                                                                           VectorType_,
                                                                           VectorMirrorType_,
                                                                           MatrixType_,
                                                                           PreconContType_,
                                                                           FilterType_,
                                                                           StorageType_,
                                                                           IT_> >& precon)
        {
          _precon = precon;
        }

        ///use solver for different data
        virtual void reset_data(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                     MemTag_,
                                                                     VectorType_,
                                                                     VectorMirrorType_,
                                                                     MatrixType_,
                                                                     PreconContType_,
                                                                     FilterType_,
                                                                     StorageType_,
                                                                     IT_>& data)
        {
          _data = &data;
        }

        /*virtual ScaRCFunctorBase& operator=(const ScaRCFunctorBase& other)
        {
          this->_precon = other._precon;
          this->_data = other._data;
          this->_eps = other._eps;
          this->_max_iters = other._max_iters;
          this->_conv_check = other._conv_check;
          this->_norm_0 = other._norm_0;
          this->_norm = other._norm;
          this->_used_iters = other._used_iters;

          return *this;
        }*/

        DataType_& eps()
        {
          return _eps;
        }

        const DataType_& eps() const
        {
          return _eps;
        }

        IT_& max_iters()
        {
          return _max_iters;
        }

        const IT_& max_iters() const
        {
          return _max_iters;
        }

        const DataType_& norm_0() const
        {
          return _norm_0;
        }

        const DataType_& norm() const
        {
          return _norm;
        }

        const IT_& iterations() const
        {
          return _used_iters;
        }

        bool& conv_check()
        {
          return _conv_check;
        }

        const bool& conv_check() const
        {
          return _conv_check;
        }

        //virtual ScaRCFunctorBase& operator=(const ScaRCFunctorBase& other) = 0;

      protected:
        std::shared_ptr<ScaRCFunctorBase<DataType_,
                                         MemTag_,
                                         VectorType_,
                                         VectorMirrorType_,
                                         MatrixType_,
                                         PreconContType_,
                                         FilterType_,
                                         StorageType_,
                                         IT_> > _precon;
        SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                    MemTag_,
                                                    VectorType_,
                                                    VectorMirrorType_,
                                                    MatrixType_,
                                                    PreconContType_,
                                                    FilterType_,
                                                    StorageType_,
                                                    IT_>* _data;

        ///config
        DataType_ _eps;
        IT_ _max_iters;
        bool _conv_check;

        ///status
        DataType_ _norm_0;
        DataType_ _norm;
        IT_ _used_iters;
    };

    template<typename DataType_,
             typename MemTag_,
             typename VectorType_,
             typename VectorMirrorType_,
             typename MatrixType_,
             typename PreconContType_,
             typename FilterType_,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index
             >
    class CompositeScaRCFunctor : public ScaRCFunctorBase<DataType_,
                                                          MemTag_,
                                                          VectorType_,
                                                          VectorMirrorType_,
                                                          MatrixType_,
                                                          PreconContType_,
                                                          FilterType_,
                                                          StorageType_,
                                                          IT_>
    {
      public:

        virtual void execute()
        {
          for(auto& f_i : _functors)
            f_i->execute();
        }

        virtual void apply(VectorType_& store_to, VectorType_& apply_to, VectorType_& apply_rhs)
        {
          for(auto& f_i : _functors)
            f_i->apply(store_to, apply_to, apply_rhs);
        }

        ///returns a string that describes the functor
        virtual const std::string type_name()
        {
          std::string result("CompositeScaRCFunctor[");

          for(auto& f_i : _functors)
          {
            result.append((const std::string)(f_i->type_name()));
            result.append(" ");

          }
          result.append("]");
          return result;
        }

        ///CTOR
        CompositeScaRCFunctor(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                          MemTag_,
                                                                          VectorType_,
                                                                          VectorMirrorType_,
                                                                          MatrixType_,
                                                                          PreconContType_,
                                                                          FilterType_,
                                                                          StorageType_,
                                                                          IT_>& data) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data),
          _functors()
        {
        }

        ///CTOR
        CompositeScaRCFunctor(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                          MemTag_,
                                                                          VectorType_,
                                                                          VectorMirrorType_,
                                                                          MatrixType_,
                                                                          PreconContType_,
                                                                          FilterType_,
                                                                          StorageType_,
                                                                          IT_>& data,
                              std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                               MemTag_,
                                                               VectorType_,
                                                               VectorMirrorType_,
                                                               MatrixType_,
                                                               PreconContType_,
                                                               FilterType_,
                                                               StorageType_,
                                                               IT_> >& precon) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data, precon),
          _functors()
        {
        }

        ///copy-CTOR
        CompositeScaRCFunctor(const CompositeScaRCFunctor& other)
        {
          this->_data = other._data;
          this->_precon = other.precon;
          this->_functors = other._functors;
        }

        virtual CompositeScaRCFunctor& operator=(const CompositeScaRCFunctor& other)
        {
          this->_precon = other.precon;
          this->_data = other.data;
          this->_functors = other._functors;
        }

      protected:
        StorageType_<std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                      MemTag_,
                                                      VectorType_,
                                                      VectorMirrorType_,
                                                      MatrixType_,
                                                      PreconContType_,
                                                      FilterType_,
                                                      StorageType_,
                                                      IT_> >,
                     std::allocator<std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                                     MemTag_,
                                                                     VectorType_,
                                                                     VectorMirrorType_,
                                                                     MatrixType_,
                                                                     PreconContType_,
                                                                     FilterType_,
                                                                     StorageType_,
                                                                     IT_> > > > _functors;
    };

    template<typename DataType_,
             typename MemTag_,
             typename VectorType_,
             typename VectorMirrorType_,
             typename MatrixType_,
             typename PreconContType_,
             typename FilterType_,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index
             >
    class ScaRCFunctorNULL : public ScaRCFunctorBase<DataType_,
                                                     MemTag_,
                                                     VectorType_,
                                                     VectorMirrorType_,
                                                     MatrixType_,
                                                     PreconContType_,
                                                     FilterType_,
                                                     StorageType_,
                                                     IT_>
    {
      public:

        virtual void execute()
        {
        }

        virtual void apply(VectorType_& store_to, VectorType_& apply_to, VectorType_&)
        {
          store_to.convert(apply_to);
        }

        ///returns a string that describes the functor
        virtual const std::string type_name()
        {
          return "NULL\n    ";
        }

        ///CTOR
        ScaRCFunctorNULL(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                      MemTag_,
                                                                      VectorType_,
                                                                      VectorMirrorType_,
                                                                      MatrixType_,
                                                                      PreconContType_,
                                                                      FilterType_,
                                                                      StorageType_,
                                                                      IT_>& data) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data)
        {
        }

        ///CTOR
        ScaRCFunctorNULL(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                      MemTag_,
                                                                      VectorType_,
                                                                      VectorMirrorType_,
                                                                      MatrixType_,
                                                                      PreconContType_,
                                                                      FilterType_,
                                                                      StorageType_,
                                                                      IT_>& data,
                          std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                           MemTag_,
                                                           VectorType_,
                                                           VectorMirrorType_,
                                                           MatrixType_,
                                                           PreconContType_,
                                                           FilterType_,
                                                           StorageType_,
                                                           IT_> >& precon) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data, precon)
        {
        }

        ///copy CTOR
        ScaRCFunctorNULL(const ScaRCFunctorNULL& other) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(other)
        {
        }
    };

    template<typename DataType_,
             typename MemTag_,
             typename VectorType_,
             typename VectorMirrorType_,
             typename MatrixType_,
             typename PreconContType_,
             typename FilterType_,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index
             >
    class ScaRCFunctorRichardson0 : public ScaRCFunctorBase<DataType_,
                                                           MemTag_,
                                                           VectorType_,
                                                           VectorMirrorType_,
                                                           MatrixType_,
                                                           PreconContType_,
                                                           FilterType_,
                                                           StorageType_,
                                                           IT_>
    {
      public:

        virtual void execute()
        {
          ///TODO add omega to PreconData interface

          //temp0 <- FILTER(temp0);
          //this->_data->filter().filter_def(temp0);

          this->_used_iters = 0;
          _temp0.copy(this->_data->sol());
          while(true)
          {
            auto temp_x = _temp0.clone();
            ///temp0 <- b - SYNCH(Ax)
            Foundation::GlobalDefect<MemTag_>::exec(_temp0,
                                                           this->_data->rhs(),
                                                           this->_data->sys(),
                                                           temp_x,
                                                           this->_data->vector_mirrors(),
                                                           this->_data->dest_ranks(),
                                                           this->_data->vector_mirror_sendbufs(),
                                                           this->_data->vector_mirror_recvbufs(),
                                                           this->_data->tags(),
                                                           this->_data->communicators().at(0)); ///TODO ScarcLayer selection

            if(this->_conv_check)
            {
              if(this->_used_iters == 0)
              {
                Foundation::GlobalNorm2<MemTag_>::value(
                                                               this->_norm_0,
                                                               _temp0,
                                                               this->_data->halo_frequencies()
                    );
              }
              else
              {
                Foundation::GlobalNorm2<MemTag_>::value(
                    this->_norm,
                    _temp0,
                    this->_data->halo_frequencies()
                    );
              }
            }

            if(this->_conv_check)
            {
              if(this->_norm < this->_eps * this->_norm_0 || this->_used_iters >= this->_max_iters)
                break;
            }
            else
              if(this->_used_iters >= this->_max_iters)
                break;

            ++(this->_used_iters);

            this->_precon->apply(_temp1, _temp0, this->_data->rhs());
            _temp0.axpy(this->_data->sol(), _temp1);
            this->_data->sol().copy(_temp0);

          }

        }

        virtual void apply(VectorType_& store_to, VectorType_& apply_to, VectorType_& apply_rhs)
        {
          ///TODO add omega to PreconData interface

          //temp0 <- FILTER(temp0);
          //this->_data->filter().filter_def(temp0);

          this->_used_iters = 0;
          _temp0.copy(apply_to);
          store_to.copy(apply_to);
          while(true)
          {
            auto temp_x = _temp0.clone();
            ///temp0 <- b - SYNCH(Ax)
            Foundation::GlobalDefect<MemTag_>::exec(_temp0,
                                               apply_rhs,
                                               this->_data->sys(),
                                               temp_x,
                                               this->_data->vector_mirrors(),
                                               this->_data->dest_ranks(),
                                               this->_data->vector_mirror_sendbufs(),
                                               this->_data->vector_mirror_recvbufs(),
                                               this->_data->tags(),
                                               this->_data->communicators().at(0)); ///TODO ScarcLayer selection

            if(this->_conv_check)
            {
              if(this->_used_iters == 0)
              {
                Foundation::GlobalNorm2<MemTag_>::value(
                    this->_norm_0,
                    _temp0,
                    this->_data->halo_frequencies()
                    );
              }
              else
              {
                Foundation::GlobalNorm2<MemTag_>::value(
                    this->_norm,
                    _temp0,
                    this->_data->halo_frequencies()
                    );
              }
            }

            if(this->_conv_check)
            {
              if(this->_norm < this->_eps * this->_norm_0 || this->_used_iters >= this->_max_iters)
                break;
            }
            else
              if(this->_used_iters >= this->_max_iters)
                break;


            ++(this->_used_iters);

            this->_precon->apply(_temp1, _temp0, apply_rhs);
            _temp0.axpy(store_to, _temp1);
            store_to.copy(_temp0);
          }
        }

        ///returns a string that describes the functor
        virtual const std::string type_name()
        {
          return "Richardson0\n    ";
        }

        ///CTOR
        ScaRCFunctorRichardson0(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                     MemTag_,
                                                                     VectorType_,
                                                                     VectorMirrorType_,
                                                                     MatrixType_,
                                                                     PreconContType_,
                                                                     FilterType_,
                                                                     StorageType_,
                                                                     IT_>& data) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data),
          _temp0(data.sol().size()),
          _temp1(data.sol().size())
        {
          this->_precon = std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                           MemTag_,
                                                           VectorType_,
                                                           VectorMirrorType_,
                                                           MatrixType_,
                                                           PreconContType_,
                                                           FilterType_,
                                                           StorageType_,
                                                           IT_> >(new ScaRCFunctorNULL<DataType_,
                                                                                         MemTag_,
                                                                                         VectorType_,
                                                                                         VectorMirrorType_,
                                                                                         MatrixType_,
                                                                                         PreconContType_,
                                                                                         FilterType_,
                                                                                         StorageType_,
                                                                                         IT_>(data));
        }

        ///CTOR
        ScaRCFunctorRichardson0(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                           MemTag_,
                                                                           VectorType_,
                                                                           VectorMirrorType_,
                                                                           MatrixType_,
                                                                           PreconContType_,
                                                                           FilterType_,
                                                                           StorageType_,
                                                                           IT_>& data,
                               std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                                MemTag_,
                                                                VectorType_,
                                                                VectorMirrorType_,
                                                                MatrixType_,
                                                                PreconContType_,
                                                                FilterType_,
                                                                StorageType_,
                                                                IT_> >& precon) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data, precon),
          _temp0(data.sol().size()),
          _temp1(data.sol().size())
        {
        }

        ///copy CTOR
        ScaRCFunctorRichardson0(const ScaRCFunctorRichardson0& other) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(other),
          _temp0(other._data.sol().size()),
          _temp1(other._data.sol().size())
        {
        }

      private:
        VectorType_ _temp0;
        VectorType_ _temp1;
    };

    template<typename DataType_,
             typename MemTag_,
             typename VectorType_,
             typename VectorMirrorType_,
             typename MatrixType_,
             typename PreconContType_,
             typename FilterType_,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index
             >
    class ScaRCFunctorPreconSpM1V1 : public ScaRCFunctorBase<DataType_,
                                                           MemTag_,
                                                           VectorType_,
                                                           VectorMirrorType_,
                                                           MatrixType_,
                                                           PreconContType_,
                                                           FilterType_,
                                                           StorageType_,
                                                           IT_>
    {
      public:

        virtual void execute()
        {
          this->_data->precon().apply(this->_data->sol(), this->_data->sol());
        }

        virtual void apply(VectorType_& store_to, VectorType_& apply_to, VectorType_&)
        {
          this->_data->precon().apply(store_to, apply_to);
        }

        ///returns a string that describes the functor
        virtual const std::string type_name()
        {
          return "PreconSpM1V1\n    ";
        }

        ///CTOR
        ScaRCFunctorPreconSpM1V1(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                     MemTag_,
                                                                     VectorType_,
                                                                     VectorMirrorType_,
                                                                     MatrixType_,
                                                                     PreconContType_,
                                                                     FilterType_,
                                                                     StorageType_,
                                                                     IT_>& data) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data)
        {
          this->_precon = std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                           MemTag_,
                                                           VectorType_,
                                                           VectorMirrorType_,
                                                           MatrixType_,
                                                           PreconContType_,
                                                           FilterType_,
                                                           StorageType_,
                                                           IT_> >(new ScaRCFunctorNULL<DataType_,
                                                                                         MemTag_,
                                                                                         VectorType_,
                                                                                         VectorMirrorType_,
                                                                                         MatrixType_,
                                                                                         PreconContType_,
                                                                                         FilterType_,
                                                                                         StorageType_,
                                                                                         IT_>(data));
        }

        ///CTOR
        ScaRCFunctorPreconSpM1V1(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                           MemTag_,
                                                                           VectorType_,
                                                                           VectorMirrorType_,
                                                                           MatrixType_,
                                                                           PreconContType_,
                                                                           FilterType_,
                                                                           StorageType_,
                                                                           IT_>& data,
                               std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                                MemTag_,
                                                                VectorType_,
                                                                VectorMirrorType_,
                                                                MatrixType_,
                                                                PreconContType_,
                                                                FilterType_,
                                                                StorageType_,
                                                                IT_> >& precon) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data, precon)
        {
        }

        ///copy CTOR
        ScaRCFunctorPreconSpM1V1(const ScaRCFunctorPreconSpM1V1& other) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(other)
        {
        }
    };

    template<typename DataType_,
             typename MemTag_,
             typename VectorType_,
             typename VectorMirrorType_,
             typename MatrixType_,
             typename PreconContType_,
             typename FilterType_,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index
             >
    class ScaRCFunctorRichardson1 : public ScaRCFunctorBase<DataType_,
                                                           MemTag_,
                                                           VectorType_,
                                                           VectorMirrorType_,
                                                           MatrixType_,
                                                           PreconContType_,
                                                           FilterType_,
                                                           StorageType_,
                                                           IT_>
    {
      public:

        virtual void execute()
        {
          ///TODO add omega to PreconData interface

          //temp0 <- FILTER(temp0);
          //this->_data->filter().filter_def(temp0);

          this->_used_iters = 0;
          _temp0.copy(this->_data->sol());
          while(true)
          {
            auto temp_x = _temp0.clone();
            ///temp0 <- b - Ax
            this->_data->localsys().apply(_temp0, temp_x, this->_data->rhs(), -DataType_(1));
            if(this->_conv_check)
            {
              if(this->_used_iters == 0)
                this->_norm_0 = this->_temp0.norm2();
              else
                this->_norm = this->_temp0.norm2();
            }

            if(this->_conv_check)
            {
              if(this->_norm < this->_eps * this->_norm_0 || this->_used_iters >= this->_max_iters)
                break;
            }
            else
              if(this->_used_iters >= this->_max_iters)
                break;

            ++(this->_used_iters);

            this->_precon->apply(_temp1, _temp0, this->_data->rhs());
            _temp0.axpy(this->_data->sol(), _temp1);
            this->_data->sol().copy(_temp0);
          }

        }

        virtual void apply(VectorType_& store_to, VectorType_& apply_to, VectorType_& apply_rhs)
        {
          ///TODO add omega to PreconData interface

          //temp0 <- FILTER(temp0);
          //this->_data->filter().filter_def(temp0);

          this->_used_iters = 0;
          _temp0.copy(apply_to);
          store_to.copy(apply_to);
          while(true)
          {
            auto temp_x = _temp0.clone();
            ///temp0 <- b - Ax
            this->_data->localsys().apply(_temp0, temp_x, apply_rhs, -DataType_(1));

            if(this->_conv_check)
            {
              if(this->_used_iters == 0)
              {
                this->_norm_0 = _temp0.norm2();
              }
              else
                this->_norm = _temp0.norm2();
            }

            if(this->_conv_check)
            {
              if(this->_norm < this->_eps * this->_norm_0 || this->_used_iters >= this->_max_iters)
                break;
            }
            else
              if(this->_used_iters >= this->_max_iters)
                break;

            ++(this->_used_iters);

            this->_precon->apply(_temp1, _temp0, apply_rhs);
            _temp0.axpy(store_to, _temp1);
            store_to.copy(_temp0);
          }

        }

        ///returns a string that describes the functor
        virtual const std::string type_name()
        {
          return "Richardson1\n    ";
        }

        ///CTOR
        ScaRCFunctorRichardson1(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                     MemTag_,
                                                                     VectorType_,
                                                                     VectorMirrorType_,
                                                                     MatrixType_,
                                                                     PreconContType_,
                                                                     FilterType_,
                                                                     StorageType_,
                                                                     IT_>& data) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data),
          _temp0(data.sol().size()),
          _temp1(data.sol().size())
        {
          this->_precon = std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                           MemTag_,
                                                           VectorType_,
                                                           VectorMirrorType_,
                                                           MatrixType_,
                                                           PreconContType_,
                                                           FilterType_,
                                                           StorageType_,
                                                           IT_> >(new ScaRCFunctorNULL<DataType_,
                                                                                         MemTag_,
                                                                                         VectorType_,
                                                                                         VectorMirrorType_,
                                                                                         MatrixType_,
                                                                                         PreconContType_,
                                                                                         FilterType_,
                                                                                         StorageType_,
                                                                                         IT_>(data));
        }

        ///CTOR
        ScaRCFunctorRichardson1(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                           MemTag_,
                                                                           VectorType_,
                                                                           VectorMirrorType_,
                                                                           MatrixType_,
                                                                           PreconContType_,
                                                                           FilterType_,
                                                                           StorageType_,
                                                                           IT_>& data,
                               std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                                MemTag_,
                                                                VectorType_,
                                                                VectorMirrorType_,
                                                                MatrixType_,
                                                                PreconContType_,
                                                                FilterType_,
                                                                StorageType_,
                                                                IT_> >& precon) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data, precon),
          _temp0(data.sol().size()),
          _temp1(data.sol().size())
        {
        }

        ///copy CTOR
        ScaRCFunctorRichardson1(const ScaRCFunctorRichardson1& other) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(other),
          _temp0(other.data.sol().size()),
          _temp1(other.data.sol().size())
        {
        }

      private:
        VectorType_ _temp0;
        VectorType_ _temp1;
    };

    template<typename DataType_,
             typename MemTag_,
             typename VectorType_,
             typename VectorMirrorType_,
             typename MatrixType_,
             typename PreconContType_,
             typename FilterType_,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index
             >
    class ScaRCFunctorPreconBlock : public ScaRCFunctorBase<DataType_,
                                                           MemTag_,
                                                           VectorType_,
                                                           VectorMirrorType_,
                                                           MatrixType_,
                                                           PreconContType_,
                                                           FilterType_,
                                                           StorageType_,
                                                           IT_>
    {
      public:

        virtual void execute()
        {
          //global defect
          Foundation::GlobalDefect<MemTag_>::exec(this->_data->def(),
                                             this->_data->rhs(),
                                             this->_data->sys(),
                                             this->_data->sol(),
                                             this->_data->vector_mirrors(),
                                             this->_data->dest_ranks(),
                                             this->_data->vector_mirror_sendbufs(),
                                             this->_data->vector_mirror_recvbufs(),
                                             this->_data->tags(),
                                             this->_data->communicators().at(0)); ///TODO ScarcLayer selection

          _temp0.format();
          this->_precon->apply(_temp1, _temp0, this->_data->def());

          //synchronisation
          Foundation::GlobalSynchVec1<MemTag_>::exec(
                                                 _temp1,
                                                 this->_data->vector_mirrors(),
                                                 this->_data->halo_frequencies(),
                                                 this->_data->dest_ranks(),
                                                 this->_data->vector_mirror_sendbufs(),
                                                 this->_data->vector_mirror_recvbufs(),
                                                 this->_data->tags(),
                                                 this->_data->communicators().at(0)
                                                );

          //global correction
          this->_data->sol().axpy(this->_data->sol(), _temp1);

        }

        virtual void apply(VectorType_& store_to, VectorType_& apply_to, VectorType_& apply_rhs)
        {
          Foundation::GlobalDefect<MemTag_>::exec(_temp0,
                                             apply_rhs,
                                             this->_data->sys(),
                                             apply_to,
                                             this->_data->vector_mirrors(),
                                             this->_data->dest_ranks(),
                                             this->_data->vector_mirror_sendbufs(),
                                             this->_data->vector_mirror_recvbufs(),
                                             this->_data->tags(),
                                             this->_data->communicators().at(0)); ///TODO ScarcLayer selection

          //local solver, temp0 is now our rhs
          _temp1.format();
          this->_precon->apply(_temp2, _temp1, _temp0);

          //synchronisation
          Foundation::GlobalSynchVec1<MemTag_>::exec(
                                                 _temp2,
                                                 this->_data->vector_mirrors(),
                                                 this->_data->halo_frequencies(),
                                                 this->_data->dest_ranks(),
                                                 this->_data->vector_mirror_sendbufs(),
                                                 this->_data->vector_mirror_recvbufs(),
                                                 this->_data->tags(),
                                                 this->_data->communicators().at(0)
                                                );

          //global correction
          store_to.axpy(apply_to, _temp2);
        }

        ///returns a string that describes the functor
        virtual const std::string type_name()
        {
          return "PreconBlock\n    ";
        }

        ///CTOR
        ScaRCFunctorPreconBlock(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                     MemTag_,
                                                                     VectorType_,
                                                                     VectorMirrorType_,
                                                                     MatrixType_,
                                                                     PreconContType_,
                                                                     FilterType_,
                                                                     StorageType_,
                                                                     IT_>& data) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data),
          _temp0(data.sol().size()),
          _temp1(data.sol().size()),
          _temp2(data.sol().size())
        {
          this->_precon = std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                           MemTag_,
                                                           VectorType_,
                                                           VectorMirrorType_,
                                                           MatrixType_,
                                                           PreconContType_,
                                                           FilterType_,
                                                           StorageType_,
                                                           IT_> >(new ScaRCFunctorNULL<DataType_,
                                                                                         MemTag_,
                                                                                         VectorType_,
                                                                                         VectorMirrorType_,
                                                                                         MatrixType_,
                                                                                         PreconContType_,
                                                                                         FilterType_,
                                                                                         StorageType_,
                                                                                         IT_>(data));
        }

        ///CTOR
        ScaRCFunctorPreconBlock(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                           MemTag_,
                                                                           VectorType_,
                                                                           VectorMirrorType_,
                                                                           MatrixType_,
                                                                           PreconContType_,
                                                                           FilterType_,
                                                                           StorageType_,
                                                                           IT_>& data,
                               std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                                MemTag_,
                                                                VectorType_,
                                                                VectorMirrorType_,
                                                                MatrixType_,
                                                                PreconContType_,
                                                                FilterType_,
                                                                StorageType_,
                                                                IT_> >& precon) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data, precon),
          _temp0(data.sol().size()),
          _temp1(data.sol().size()),
          _temp2(data.sol().size())
        {
        }

        ///copy CTOR
        ScaRCFunctorPreconBlock(const ScaRCFunctorPreconBlock& other) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(other),
          _temp0(other.data.sol().size()),
          _temp1(other.data.sol().size()),
          _temp2(other.data.sol().size())
        {
        }

      private:
        VectorType_ _temp0;
        VectorType_ _temp1;
        VectorType_ _temp2;
    };

    template<typename DataType_,
             typename MemTag_,
             typename VectorType_,
             typename VectorMirrorType_,
             typename MatrixType_,
             typename PreconContType_,
             typename FilterType_,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index
             >
    class ScaRCFunctorPCG0 : public ScaRCFunctorBase<DataType_,
                                                      MemTag_,
                                                      VectorType_,
                                                      VectorMirrorType_,
                                                      MatrixType_,
                                                      PreconContType_,
                                                      FilterType_,
                                                      StorageType_,
                                                      IT_>
    {
      public:

        virtual void execute()
        {
          Foundation::GlobalDefect<MemTag_>::exec(_r,
                                             this->_data->rhs(),
                                             this->_data->sys(),
                                             this->_data->sol(),
                                             this->_data->vector_mirrors(),
                                             this->_data->dest_ranks(),
                                             this->_data->vector_mirror_sendbufs(),
                                             this->_data->vector_mirror_recvbufs(),
                                             this->_data->tags(),
                                             this->_data->communicators().at(0)); ///TODO ScarcLayer selection

          this->_precon->apply(_p, _r, this->_data->rhs());

          if(this->_conv_check)
          {
            Foundation::GlobalNorm2<MemTag_>::value(
                this->_norm_0,
                _r,
                this->_data->halo_frequencies()
                );
          }

          Foundation::GlobalDot<MemTag_>::value(
                                           _alpha_new,
                                           _r,
                                           _p,
                                           this->_data->halo_frequencies()
                                           );

          DataType_ temp;
          ///TODO: config objs!
          while(this->_used_iters < this->_max_iters)
          {
            Foundation::GlobalProductMat0Vec1<MemTag_>::exec(
                                                        _v,
                                                        this->_data->sys(),
                                                        _p,
                                                        this->_data->vector_mirrors(),
                                                        this->_data->dest_ranks(),
                                                        this->_data->vector_mirror_sendbufs(),
                                                        this->_data->vector_mirror_recvbufs(),
                                                        this->_data->tags(),
                                                        this->_data->communicators().at(0)
                                                       );

            Foundation::GlobalDot<MemTag_>::value(
                                             temp,
                                             _v,
                                             _p,
                                             this->_data->halo_frequencies()
                                            );

            _lambda = _alpha_new / (std::abs(temp) > std::numeric_limits<DataType_>::epsilon() ? temp : (std::numeric_limits<DataType_>::epsilon()));

            ++(this->_used_iters);

            //x <- x + lambda * p
            this->_data->sol().axpy(_p, this->_data->sol(), _lambda);
            //r <- r - lambda * v
            _r.axpy(_v, _r, -_lambda);

            if(this->_conv_check)
            {
              Foundation::GlobalNorm2<MemTag_>::value(
                  this->_norm,
                  _r,
                  this->_data->halo_frequencies()
                  );
            }

            if(this->_conv_check)
            {
              if(this->_norm < this->_eps * this->_norm_0 || this->_norm < this->_eps || this->_used_iters == this->_max_iters)
              {
                break;
              }
            }
            else
              if(this->_used_iters == this->_max_iters)
              {
                break;
              }

            this->_precon->apply(_z, _r, this->_data->rhs());

            _alpha = _alpha_new;

            Foundation::GlobalDot<MemTag_>::value(
                                             _alpha_new,
                                             _r,
                                             _z,
                                             this->_data->halo_frequencies()
                                            );

            DataType_ talpha_new(_alpha_new / (std::abs(_alpha) > std::numeric_limits<DataType_>::epsilon() ? _alpha : (std::numeric_limits<DataType_>::epsilon())));

            //p <- talpha_new * p
            _p.scale(_p, talpha_new);

            //p <- p + z
            _p.axpy(_p, _z);
          }
        }

        virtual void apply(VectorType_& store_to, VectorType_& apply_to, VectorType_& apply_rhs)
        {
          Foundation::GlobalDefect<MemTag_>::exec(_r,
                                             apply_rhs,
                                             this->_data->sys(),
                                             apply_to,
                                             this->_data->vector_mirrors(),
                                             this->_data->dest_ranks(),
                                             this->_data->vector_mirror_sendbufs(),
                                             this->_data->vector_mirror_recvbufs(),
                                             this->_data->tags(),
                                             this->_data->communicators().at(0)); ///TODO ScarcLayer selection

          this->_precon->apply(_p, _r, this->_data->rhs());

          if(this->_conv_check)
          {
            Foundation::GlobalNorm2<MemTag_>::value(
                this->_norm_0,
                _r,
                this->_data->halo_frequencies()
                );
          }

          Foundation::GlobalDot<MemTag_>::value(
                                           _alpha_new,
                                           _r,
                                           _p,
                                           this->_data->halo_frequencies()
                                           );

          DataType_ temp;
          store_to.copy(apply_to);
          while(this->_used_iters < this->_max_iters)
          {
            Foundation::GlobalProductMat0Vec1<MemTag_>::exec(
                                                        _v,
                                                        this->_data->sys(),
                                                        _p,
                                                        this->_data->vector_mirrors(),
                                                        this->_data->dest_ranks(),
                                                        this->_data->vector_mirror_sendbufs(),
                                                        this->_data->vector_mirror_recvbufs(),
                                                        this->_data->tags(),
                                                        this->_data->communicators().at(0)
                                                       );

            Foundation::GlobalDot<MemTag_>::value(
                                             temp,
                                             _v,
                                             _p,
                                             this->_data->halo_frequencies()
                                            );

            _lambda = _alpha_new / (std::abs(temp) > std::numeric_limits<DataType_>::epsilon() ? temp : (std::numeric_limits<DataType_>::epsilon()));

            ++(this->_used_iters);

            //x <- x + lambda * p
            store_to.axpy(_p, store_to, _lambda);
            //r <- r - lambda * v
            _r.axpy(_v, _r, -_lambda);

            if(this->_conv_check)
            {
              Foundation::GlobalNorm2<MemTag_>::value(
                  this->_norm,
                  _r,
                  this->_data->halo_frequencies()
                  );
            }

            if(this->_conv_check)
            {
              if(this->_norm < this->_eps * this->_norm_0 || this->_norm < this->_eps || this->_used_iters == this->_max_iters)
              {
                break;
              }
            }
            else
              if(this->_used_iters == this->_max_iters)
              {
                break;
              }

            this->_precon->apply(_z, _r, apply_rhs);

            _alpha = _alpha_new;

            Foundation::GlobalDot<MemTag_>::value(
                                             _alpha_new,
                                             _r,
                                             _z,
                                             this->_data->halo_frequencies()
                                            );

            DataType_ talpha_new(_alpha_new / (std::abs(_alpha) > std::numeric_limits<DataType_>::epsilon() ? _alpha : (std::numeric_limits<DataType_>::epsilon())));

            //p <- talpha_new * p
            _p.scale(_p, talpha_new);

            //p <- p + z
            _p.axpy(_p, _z);
          }
        }

        ///returns a string that describes the functor
        virtual const std::string type_name()
        {
          return "PCG0\n    ";
        }

        ///CTOR
        ScaRCFunctorPCG0(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                     MemTag_,
                                                                     VectorType_,
                                                                     VectorMirrorType_,
                                                                     MatrixType_,
                                                                     PreconContType_,
                                                                     FilterType_,
                                                                     StorageType_,
                                                                     IT_>& data) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data),
          _p(data.sol().size()),
          _r(data.sol().size()),
          _v(data.sol().size()),
          _z(data.sol().size()),
          _alpha(0),
          _alpha_new(0),
          _lambda(0)
        {
          this->_precon = std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                           MemTag_,
                                                           VectorType_,
                                                           VectorMirrorType_,
                                                           MatrixType_,
                                                           PreconContType_,
                                                           FilterType_,
                                                           StorageType_,
                                                           IT_> >(new ScaRCFunctorNULL<DataType_,
                                                                                         MemTag_,
                                                                                         VectorType_,
                                                                                         VectorMirrorType_,
                                                                                         MatrixType_,
                                                                                         PreconContType_,
                                                                                         FilterType_,
                                                                                         StorageType_,
                                                                                         IT_>(data));
        }

        ///CTOR
        ScaRCFunctorPCG0(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                           MemTag_,
                                                                           VectorType_,
                                                                           VectorMirrorType_,
                                                                           MatrixType_,
                                                                           PreconContType_,
                                                                           FilterType_,
                                                                           StorageType_,
                                                                           IT_>& data,
                               std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                                MemTag_,
                                                                VectorType_,
                                                                VectorMirrorType_,
                                                                MatrixType_,
                                                                PreconContType_,
                                                                FilterType_,
                                                                StorageType_,
                                                                IT_> >& precon) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data, precon),
          _p(data.sol().size()),
          _r(data.sol().size()),
          _v(data.sol().size()),
          _z(data.sol().size()),
          _alpha(0),
          _alpha_new(0),
          _lambda(0)
        {
        }

        ///copy CTOR
        ScaRCFunctorPCG0(const ScaRCFunctorPCG0& other) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(other),
          _p(other._data.sol().size()),
          _r(other._data.sol().size()),
          _v(other._data.sol().size()),
          _z(other._data.sol().size()),
          _alpha(0),
          _alpha_new(0),
          _lambda(0)
        {
        }

      private:
        VectorType_ _p;
        VectorType_ _r;
        VectorType_ _v;
        VectorType_ _z;

        DataType_ _alpha;
        DataType_ _alpha_new;
        DataType_ _lambda;

    };

    template<typename DataType_,
             typename MemTag_,
             typename VectorType_,
             typename VectorMirrorType_,
             typename MatrixType_,
             typename PreconContType_,
             typename FilterType_,
             template<typename, typename> class StorageType_ = std::vector,
             typename IT_ = Index
             >
    class ScaRCFunctorPCG1 : public ScaRCFunctorBase<DataType_,
                                                           MemTag_,
                                                           VectorType_,
                                                           VectorMirrorType_,
                                                           MatrixType_,
                                                           PreconContType_,
                                                           FilterType_,
                                                           StorageType_,
                                                           IT_>
    {
      public:

        virtual void execute()
        {
          this->_data->localsys().apply(_r, this->_data->sol(), this->_data->rhs(), -DataType_(1));

          this->_precon->apply(_p, _r, this->_data->rhs());

          if(this->_conv_check)
            this->_norm_0 = _r.norm2();

          _alpha_new = _r.dot(_p);

          DataType_ temp;
          while(this->_used_iters < this->_max_iters)
          {
            this->_data->localsys().apply(_v, _p);

            temp = _v.dot(_p);

            _lambda = _alpha_new / (std::abs(temp) > std::numeric_limits<DataType_>::epsilon() ? temp : (std::numeric_limits<DataType_>::epsilon()));

            ++(this->_used_iters);

            //x <- x + lambda * p
            this->_data->sol().axpy(_p, this->_data->sol(), _lambda);
            //r <- r - lambda * v
            _r.axpy(_v, _r, -_lambda);

            if(this->_conv_check)
            {
              this->_norm = _r.norm2();
            }

            if(this->_conv_check)
            {
              if(this->_norm < this->_eps * this->_norm_0 || this->_norm < this->_eps || this->_used_iters == this->_max_iters)
              {
                break;
              }
            }
            else
              if(this->_used_iters == this->_max_iters)
              {
                break;
              }

            this->_precon->apply(_z, _r, this->_data->rhs());

            _alpha = _alpha_new;

            _alpha_new = _r.dot(_z);

            DataType_ talpha_new(_alpha_new / (std::abs(_alpha) > std::numeric_limits<DataType_>::epsilon() ? _alpha : (std::numeric_limits<DataType_>::epsilon())));

            //p <- talpha_new * p
            _p.scale(_p, talpha_new);

            //p <- p + z
            _p.axpy(_p, _z);
          }
        }

        virtual void apply(VectorType_& store_to, VectorType_& apply_to, VectorType_& apply_rhs)
        {
          this->_data->localsys().apply(_r, apply_to, apply_rhs, -DataType_(1));

          this->_precon->apply(_p, _r, apply_rhs);

          if(this->_conv_check)
            this->_norm_0 = _r.norm2();

          _alpha_new = _r.dot(_p);

          DataType_ temp;
          store_to.copy(apply_to);
          while(this->_used_iters < this->_max_iters)
          {
            this->_data->localsys().apply(_v, _p);

            temp = _v.dot(_p);

            _lambda = _alpha_new / (std::abs(temp) > std::numeric_limits<DataType_>::epsilon() ? temp : (std::numeric_limits<DataType_>::epsilon()));

            ++(this->_used_iters);

            //x <- x + lambda * p
            store_to.axpy(_p, store_to, _lambda);
            //r <- r - lambda * v
            _r.axpy(_v, _r, -_lambda);

            if(this->_conv_check)
            {
              this->_norm = _r.norm2();
            }

            if(this->_conv_check)
            {
              if(this->_norm < this->_eps * this->_norm_0 || this->_norm < this->_eps || this->_used_iters == this->_max_iters)
              {
                break;
              }
            }
            else
              if(this->_used_iters == this->_max_iters)
              {
                break;
              }

            this->_precon->apply(_z, _r, apply_rhs);

            _alpha = _alpha_new;

            _alpha_new = _r.dot(_z);

            DataType_ talpha_new(_alpha_new / (std::abs(_alpha) > std::numeric_limits<DataType_>::epsilon() ? _alpha : (std::numeric_limits<DataType_>::epsilon())));

            //p <- talpha_new * p
            _p.scale(_p, talpha_new);

            //p <- p + z
            _p.axpy(_p, _z);
          }
        }

        ///returns a string that describes the functor
        virtual const std::string type_name()
        {
          return "PCG1\n    ";
        }

        ///CTOR
        ScaRCFunctorPCG1(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                     MemTag_,
                                                                     VectorType_,
                                                                     VectorMirrorType_,
                                                                     MatrixType_,
                                                                     PreconContType_,
                                                                     FilterType_,
                                                                     StorageType_,
                                                                     IT_>& data) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data),
          _p(data.sol().size()),
          _r(data.sol().size()),
          _v(data.sol().size()),
          _z(data.sol().size()),
          _alpha(0),
          _alpha_new(0),
          _lambda(0)
        {
          this->_precon = std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                           MemTag_,
                                                           VectorType_,
                                                           VectorMirrorType_,
                                                           MatrixType_,
                                                           PreconContType_,
                                                           FilterType_,
                                                           StorageType_,
                                                           IT_> >(new ScaRCFunctorNULL<DataType_,
                                                                                         MemTag_,
                                                                                         VectorType_,
                                                                                         VectorMirrorType_,
                                                                                         MatrixType_,
                                                                                         PreconContType_,
                                                                                         FilterType_,
                                                                                         StorageType_,
                                                                                         IT_>(data));
        }

        ///CTOR
        ScaRCFunctorPCG1(SynchronisedPreconditionedFilteredScaRCData<DataType_,
                                                                           MemTag_,
                                                                           VectorType_,
                                                                           VectorMirrorType_,
                                                                           MatrixType_,
                                                                           PreconContType_,
                                                                           FilterType_,
                                                                           StorageType_,
                                                                           IT_>& data,
                               std::shared_ptr<ScaRCFunctorBase<DataType_,
                                                                MemTag_,
                                                                VectorType_,
                                                                VectorMirrorType_,
                                                                MatrixType_,
                                                                PreconContType_,
                                                                FilterType_,
                                                                StorageType_,
                                                                IT_> >& precon) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(data, precon),
          _p(data.sol().size()),
          _r(data.sol().size()),
          _v(data.sol().size()),
          _z(data.sol().size()),
          _alpha(0),
          _alpha_new(0),
          _lambda(0)
        {
        }

        ///copy CTOR
        ScaRCFunctorPCG1(const ScaRCFunctorPCG1& other) :
          ScaRCFunctorBase<DataType_,
                           MemTag_,
                           VectorType_,
                           VectorMirrorType_,
                           MatrixType_,
                           PreconContType_,
                           FilterType_,
                           StorageType_,
                           IT_>(other),
          _p(other._data.sol().size()),
          _r(other._data.sol().size()),
          _v(other._data.sol().size()),
          _z(other._data.sol().size()),
          _alpha(0),
          _alpha_new(0),
          _lambda(0)
        {
        }

      private:
        VectorType_ _p;
        VectorType_ _r;
        VectorType_ _v;
        VectorType_ _z;

        DataType_ _alpha;
        DataType_ _alpha_new;
        DataType_ _lambda;

    };
  }
}

#endif
