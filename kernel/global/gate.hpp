#pragma once
#ifndef KERNEL_GLOBAL_GATE_HPP
#define KERNEL_GLOBAL_GATE_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/global/synch_vec.hpp>
#include <kernel/global/synch_scal.hpp>
#include <kernel/global/ticket.hpp>

namespace FEAT
{
  /**
   * \brief Global linear algebra namespace
   */
  namespace Global
  {
    /**
     * \brief Global gate implementation
     *
     * \author Peter Zajac
     */
    template<typename LocalVector_, typename Mirror_>
    class Gate
    {
    public:
      typedef typename LocalVector_::MemType MemType;
      typedef typename LocalVector_::DataType DataType;
      typedef typename LocalVector_::IndexType IndexType;
      typedef LAFEM::DenseVector<Mem::Main, DataType, IndexType> BufferVectorType;
      typedef Mirror_ MirrorType;

    public:
      /// communication ranks and tags
      mutable std::vector<Index> _ranks, _ctags;
      /// vector mirrors
      std::vector<Mirror_> _mirrors;
      /// frequency fector
      LocalVector_ _freqs;

      /// Our 'base' class type
      template <typename LocalMatrix2_, typename Mirror2_>
      using GateType = class Gate<LocalMatrix2_, Mirror2_>;

      /// this typedef lets you create a gate container with new Memory, Datatape and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using GateTypeByMDI = class Gate<typename LocalVector_::template ContainerType<Mem2_, DataType2_, IndexType2_>, typename Mirror_::template MirrorType<Mem2_, DataType2_, IndexType2_> >;

    public:
      explicit Gate()
      {
      }

      ~Gate()
      {
      }

      template<typename LVT2_, typename MT2_>
      void convert(const Gate<LVT2_, MT2_>& other)
      {
        this->_ranks.clear();
        this->_ctags.clear();
        this->_mirrors.clear();

        this->_ranks = other._ranks;
        this->_ctags = other._ctags;

        for(auto& other_mirrors_i : other._mirrors)
        {
          MirrorType mirrors_i;
          mirrors_i.convert(other_mirrors_i);
          this->_mirrors.push_back(std::move(mirrors_i));
        }

        this->_freqs.convert(other._freqs);
      }

      /// \brief Returns the total amount of bytes allocated.
      std::size_t bytes() const
      {
        size_t temp(0);
        for (auto& i : _mirrors)
        {
          temp += i.bytes();
        }
        temp += _freqs.bytes();
        temp += _ranks.size() * sizeof(Index);
        temp += _ctags.size() * sizeof(Index);

        return temp;
      }

      void push(Index rank, Index ctag, Mirror_&& mirror)
      {
        // push rank and tags
        _ranks.push_back(rank);
        _ctags.push_back(ctag);

        // push mirror
        _mirrors.push_back(std::move(mirror));
      }

      void compile(LocalVector_&& vector)
      {
        // initialise frequency vector
        _freqs = std::move(vector);
        _freqs.format(DataType(1));

        // loop over all mirrors
        for(std::size_t i(0); i < _mirrors.size(); ++i)
        {
          // sum up number of ranks per frequency entry, listed in different mirrors
          auto temp = _mirrors.at(i).create_buffer_vector();
          temp.format(DataType(1));

          // gather-axpy into frequency vector
          _mirrors.at(i).scatter_axpy_dual(_freqs, temp);
        }

        // invert frequencies
        _freqs.component_invert(_freqs);
      }

      /**
       * \brief Converts a type-1 vector into a type-0 vector.
       *
       * \param[inout] vector
       * On entry, the type-1 vector to be converted.
       * On exit, the converted type-0 vector.
       *
       * \note This function does not perform any synchronisation.
       */
      void from_1_to_0(LocalVector_& vector) const
      {
        if(!_ranks.empty())
        {
          vector.component_product(vector, _freqs);
        }
      }

      /**
       * \brief Synchronises a type-0 vector, resulting in a type-1 vector.
       *
       * \param[inout] vector
       * On entry, the type-0 vector to be synchronised.\n
       * On exit, the synchronised type-1 vector.
       */
      void sync_0(LocalVector_& vector) const
      {
        if(_ranks.empty())
          return;

        Global::SynchVec0::exec(
            vector, _mirrors, _freqs, _ranks, _ctags);
      }

      auto sync_0_async(LocalVector_& vector) const -> decltype(Global::SynchVec0Async::exec(vector, _mirrors, _freqs, _ranks, _ctags))
      {
        return Global::SynchVec0Async::exec(vector, _mirrors, _freqs, _ranks, _ctags);
      }

      /**
       * \brief Synchronises a type-1 vector, resulting in a type-1 vector.
       *
       * \param[inout] vector
       * On entry, the type-1 vector to be synchronised.\n
       * On exit, the synchronised type-1 vector.
       *
       * \note
       * This function effectively applies the from_1_to_0() and sync_0()
       * functions onto the input vector.
       */
      void sync_1(LocalVector_& vector) const
      {
        if(_ranks.empty())
          return;

        Global::SynchVec1::exec(
            vector, _mirrors, _freqs, _ranks, _ctags);
      }

      auto sync_1_async(LocalVector_& vector) const -> decltype(Global::SynchVec1Async::exec(vector, _mirrors, _freqs, _ranks, _ctags))
      {
        return Global::SynchVec1Async::exec(vector, _mirrors, _freqs, _ranks, _ctags);
      }

      /**
       * \brief Computes a synchronised dot-product of two type-1 vectors.
       *
       * \param[in] x, y
       * The two type-1 vector whose dot-product is to be computed.
       *
       * \returns
       * The dot-product of \p x and \p y.
       */
      DataType dot(const LocalVector_& x, const LocalVector_& y) const
      {
        if(_ranks.empty())
          return x.dot(y);
        return sum(_freqs.triple_dot(x, y));
      }

      std::shared_ptr<ScalTicket<DataType>> dot_async(const LocalVector_& x, const LocalVector_& y) const
      {
        return sum_async(_freqs.triple_dot(x, y));
      }

      /**
       * \brief Computes a reduced sum over all processes.
       *
       * \param[in] x
       * The value that is to be summarised over all processes.
       *
       * \returns
       * The reduced sum of all \p x.
       */
      DataType sum(DataType x) const
      {
        if(_ranks.empty())
          return x;
        else
          return Global::SynchScal0::value(x);
      }

      std::shared_ptr<ScalTicket<DataType>> sum_async(DataType x) const
      {
        return Global::SynchScal0Async::value(x);
      }

      /**
       * \brief Computes a reduced 2-norm over all processes.
       *
       * This function is equivalent to the call
       *    Math::sqrt(this->sum(x*x))
       *
       * \param[in] x
       * The value that is to be summarised over all processes.
       *
       * \returns
       * The reduced 2-norm of all \p x.
       */
      DataType norm2(DataType x) const
      {
        return Math::sqrt(sum(Math::sqr(x)));
      }

      std::shared_ptr<ScalTicket<DataType>> norm2_async(DataType x) const
      {
        auto ticket = sum_async(x);
        ticket->sqrt = true;
        return ticket;
      }

      /**
       * \brief Retrieve the absolute maximum value of this vector.
       *
       * \return The largest absolute value.
       */
      DataType max_element(const LocalVector_ & x) const
      {
        DataType t(x.max_element());
        if(_ranks.empty())
          return t;
        else
        {
          return Global::SynchScal0::value(t, Util::CommOperationMax());
        }
      }

      std::shared_ptr<ScalTicket<DataType>>  max_element_async(const LocalVector_ & x) const
      {
        DataType t(x.max_element());
        return Global::SynchScal0Async::value(t, Util::CommOperationMax());
      }
    }; // class Gate<...>
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_GATE_HPP
