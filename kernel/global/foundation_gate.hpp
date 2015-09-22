#pragma once
#ifndef KERNEL_GLOBAL_FOUNDATION_GATE_HPP
#define KERNEL_GLOBAL_FOUNDATION_GATE_HPP 1

#include <kernel/global/gate.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/foundation/global_synch_scal.hpp>
#include <kernel/foundation/global_synch_vec.hpp>

#include <vector>

namespace FEAST
{
  namespace Global
  {
    /**
     * \brief Foundation-based global Gate implementation
     *
     * \author Peter Zajac
     */
    template<typename LocalVector_, typename Mirror_>
    class FoundationGate :
      public Gate<LocalVector_>
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
      /// send/receive buffers
      mutable std::vector<BufferVectorType> _send_bufs, _recv_bufs;

    public:
      explicit FoundationGate()
      {
      }

      virtual ~FoundationGate()
      {
      }

      template<typename LVT2_, typename MT2_>
      void convert(const FoundationGate<LVT2_, MT2_>& other)
      {
        this->_ranks.clear();
        this->_ctags.clear();
        this->_mirrors.clear();
        this->_send_bufs.clear();
        this->_recv_bufs.clear();

        this->_ranks = other._ranks;
        this->_ctags = other._ctags;

        for(auto& other_mirrors_i : other._mirrors)
        {
          MirrorType mirrors_i;
          mirrors_i.convert(other_mirrors_i);
          this->_mirrors.push_back(std::move(mirrors_i));
        }

        for(auto& other_send_bufs_i : other._send_bufs)
        {
          BufferVectorType send_bufs_i;
          send_bufs_i.convert(other_send_bufs_i);
          this->_send_bufs.push_back(std::move(send_bufs_i));
        }

        for(auto& other_recv_bufs_i : other._recv_bufs)
        {
          BufferVectorType recv_bufs_i;
          recv_bufs_i.convert(other_recv_bufs_i);
          this->_recv_bufs.push_back(std::move(recv_bufs_i));
        }

        this->_freqs.convert(other._freqs);
      }

      void push(Index rank, Index ctag, Mirror_&& mirror)
      {
        // push rank and tags
        _ranks.push_back(rank);
        _ctags.push_back(ctag);

        // push mirror
        _mirrors.push_back(std::move(mirror));

        // push buffers
        _send_bufs.push_back(_mirrors.back().create_buffer_vector());
        _recv_bufs.push_back(_mirrors.back().create_buffer_vector());
      }

      void compile(LocalVector_&& vector)
      {
        // initialise frequency vector
        _freqs = std::move(vector);
        _freqs.format(DataType(1));

        // loop over all mirrors
        for(std::size_t i(0); i < _mirrors.size(); ++i)
        {
          // format receive buffer to 1
          _recv_bufs.at(i).format(DataType(1));

          // gather-axpy into frequency vector
          _mirrors.at(i).scatter_axpy_dual(_freqs, _recv_bufs.at(i));
        }

        // invert frequencies
        _freqs.component_invert(_freqs);
      }

      virtual void from_1_to_0(LocalVector_& vector) const override
      {
        if(!_ranks.empty())
        {
          vector.component_product(vector, _freqs);
        }
      }

      virtual void sync_0(LocalVector_& vector) const override
      {
        if(_ranks.empty())
          return;

        Foundation::GlobalSynchVec0<MemType>::exec(
          vector, _mirrors, _ranks, _send_bufs, _recv_bufs, _ctags);
      }

      virtual void sync_1(LocalVector_& vector) const override
      {
        if(_ranks.empty())
          return;
        Foundation::GlobalSynchVec1<MemType>::exec(
          vector, _mirrors, _freqs, _ranks, _send_bufs, _recv_bufs, _ctags);
      }

      virtual DataType dot(const LocalVector_& x, const LocalVector_& y) const override
      {
        if(_ranks.empty())
          return x.dot(y);
        return sum(_freqs.triple_dot(x, y));
      }

      virtual DataType sum(DataType x) const override
      {
        if(_ranks.empty())
          return x;
        else
          return Foundation::GlobalSynchScal0<MemType>::value(x, x);
      }

      DataType norm2(DataType x) const override
      {
        return Math::sqrt(sum(Math::sqr(x)));
      }
    }; // class FoundationGate<...>
  } // namespace Global
} // namespace FEAST

#endif // KERNEL_GLOBAL_FOUNDATION_GATE_HPP
