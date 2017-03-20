#pragma once
#ifndef KERNEL_FOUNDATION_PEXECUTOR_HPP
#define KERNEL_FOUNDATION_PEXECUTOR_HPP 1

#include<kernel/base_header.hpp>
#include<kernel/util/dist.hpp>
#include<kernel/foundation/base.hpp>
#include<kernel/foundation/pgraph.hpp>

#ifdef FEAT_HAVE_PARMETIS
FEAT_DISABLE_WARNINGS
#include<parmetis.h>
FEAT_RESTORE_WARNINGS
#endif

namespace FEAT
{
  namespace Foundation
  {

    template<typename IT_>
    class PExecutorBase
    {
      public:
        typedef IT_ IndexType;

        class PResult
        {
          public:

            typedef IT_ IndexType;

            ///Alloc CTOR
            PResult(IT_ s, const Dist::Comm& comm) :
              _size(s),
              _part(new IT_[Index(s)]),
              _vtxdist(new IT_[std::size_t(comm.size() + 1)]),
              _comm(&comm),
              _comm_ranks(),
              _comm_tags()
            {
            }

            IT_* get() const
            {
              return _part;
            }

            IT_* get_vtxdist() const
            {
              return _vtxdist;
            }

            const Dist::Comm& get_comm() const
            {
              return *_comm;
            }

            IT_ size() const
            {
              return _size;
            }

            PResult clone() const
            {
              PResult result(_size, *_comm);

              for(Index i(0) ; i < Index(_size) ; ++i)
                result.get()[i] = _part[i];

              for(Index i(0) ; i < Index(_comm->size() + 1) ; ++i)
                result.get_vtxdist()[i] = _vtxdist[i];

              result._comm_ranks = _comm_ranks;
              result._comm_tags = _comm_tags;

              return result;
            }

            ~PResult()
            {
              delete[] _part;
              delete[] _vtxdist;
            }

            PResult(PResult&& other) :
              _size(other._size),
              _part(other._part),
              _vtxdist(other._vtxdist),
              _comm(other._comm),
              _comm_ranks(other._comm_ranks),
              _comm_tags(other._comm_tags)
            {
              other._part = nullptr;
              other._vtxdist = nullptr;
            }

            PResult& operator=(PResult&& other)
            {
              if(this == &other)
                return *this;


              this->_size = other._size;
              this->_part = other._part;
              this->_vtxdist = other._vtxdist;
              this->_comm = other._comm;
              this->_comm_ranks = other._comm_ranks;
              this->_comm_tags = other._comm_tags;
              other._part = nullptr;
              other._vtxdist = nullptr;

              return *this;
            }

            Adjacency::Graph rank_at_element()
            {
              Adjacency::Graph result(Index(_size), Index(_comm->size()), Index(_size));

              Index* ptr = result.get_domain_ptr();
              for(Index i(0) ; i <= Index(_size) ; ++i)
                ptr[i] = i;

              Index* part = result.get_image_idx();
              for(Index i(0) ; i < Index(_size) ; ++i)
                part[i] = Index(_part[i]);

              return result;
            }

            IT_ size()
            {
              return _size;
            }

            void reset(IT_ s)
            {
              delete[] _part;
              _part = nullptr;
              _part = new IT_[std::size_t(s)];
              _size = s;
            }

            std::vector<Index>& get_comm_ranks()
            {
              return _comm_ranks;
            }

            std::vector<Index>& get_comm_ranks() const
            {
              return _comm_ranks;
            }

            std::vector<Index>& get_comm_tags()
            {
              return _comm_tags;
            }

            std::vector<Index>& get_comm_tags() const
            {
              return _comm_tags;
            }

          private:

            IT_ _size;
            IT_* _part;
            IT_* _vtxdist;

            // Note: this needs to be a pointer for the move-assignment op
            const Dist::Comm* _comm;

            std::vector<Index> _comm_ranks;
            std::vector<Index> _comm_tags;
        };
    };

    template<typename DT_, typename IT_>
    class PExecutorNONE : public PExecutorBase<IT_>
      {
        public:
          typedef PGraphNONE<DT_, IT_> PGraphT;
          typedef typename PExecutorBase<IT_>::PResult PResult;

          static PResult part(PGraphT& g)
          {
            PResult result(g.get_num_vtx(), g.get_comm());
            ///fill secondary data needed for synch later
            result.get_vtxdist()[0] = g.get_vtxdist()[0];
            result.get_vtxdist()[1] = g.get_vtxdist()[1];

            for(Index i(0) ; i < g.get_vtxdist()[1] - g.get_vtxdist()[0]; ++i)
              result.get()[i] = IT_(0);

            return result;
          }

          static PResult& fill_comm_structs_global(PResult& synched_part, const PGraphT& base)
          {
            if(base.ready_forexec)
              return synched_part;
            else
              throw(InternalError("EEEK!"));
          }
      };

#ifdef FEAT_HAVE_MPI
    template<typename DT_, typename IT_>
    class PExecutorFallback : public PExecutorBase<IT_>
      {
        public:
          typedef PGraphFallback<DT_, IT_> PGraphT;
          typedef typename PExecutorBase<IT_>::PResult PResult;

          static PResult part(PGraphT& g)
          {
            PResult result(g.get_num_vtx(), g.get_comm());
            ///fill secondary data needed for synch later
            for(Index i(0) ; i < Index(g.get_comm().size()) + 1 ; ++i)
              result.get_vtxdist()[i] = g.get_vtxdist()[i];

            for(Index i(0) ; i < g.get_vtxdist()[g.get_comm().rank() + 1] - g.get_vtxdist()[g.get_comm().rank()]; ++i)
              result.get()[i] = Index(g.get_comm().rank());

            return result;
          }

          static PResult& fill_comm_structs_global(PResult& synched_part, const PGraphT& base)
          {
            /// determine elems on this process (rank p)
            const Index commsize = Index(synched_part.get_comm().size());
            Index p = Index(synched_part.get_comm().rank());
            std::vector<Index> e_p;
            for(Index i(0) ; i < Index(synched_part.size()) ; ++i)
            {
              if(Index(synched_part.get()[i]) == p)
                e_p.push_back(i);
            }

            ///determine commlinks' comm ranks, i.e. (p, cr_p_i)
            std::vector<Index> cr_p;
            for(auto& e_p_i : e_p)
            {
              const Index adj_start(Index(base.get_xadj()[e_p_i]));
              const Index adj_end(Index(base.get_xadj()[e_p_i + 1]));

              for(Index i(adj_start) ; i < adj_end ; ++i)
                if(Index(synched_part.get()[ base.get_adjncy()[i] ]) != p)
                {
                  //comm link found!
                  cr_p.push_back(Index(synched_part.get()[Index(base.get_adjncy()[i])]));
                }
            }
            std::sort(cr_p.begin(), cr_p.end());
            auto last = std::unique(cr_p.begin(), cr_p.end());
            cr_p.erase(last, cr_p.end());

            ///synchronize and communicate commlinks, determine commtags
            //first, synchronize partition sizes
            Index* part_sizes_recvbuf(new Index[commsize]);
            Index part_sizes_sendbuf(Index(e_p.size()));

            //Util::Comm::allgather(&part_sizes_sendbuf, 1, part_sizes_recvbuf, 1, synched_part.get_comm());
            synched_part.get_comm().allgather(&part_sizes_sendbuf, std::size_t(1), part_sizes_recvbuf, std::size_t(1));

            //second, synchronize cr sizes
            Index* cr_sizes_recvbuf(new Index[commsize]);
            Index cr_sizes_sendbuf(Index(cr_p.size()));

            //Util::Comm::allgather(&cr_sizes_sendbuf, 1, cr_sizes_recvbuf, 1, synched_part.get_comm());
            synched_part.get_comm().allgather(&cr_sizes_sendbuf, std::size_t(1), cr_sizes_recvbuf, std::size_t(1));

            //then, synchronize cr
            Index* sendbuf = cr_p.data();

            int* sendcounts = new int[commsize];
            for(Index i(0) ; i < commsize ; ++i)
              sendcounts[i] = int(cr_p.size());

            int* sdispls = new int[commsize];
            for(Index i(0) ; i < commsize ; ++i)
              sdispls[i] = int(0);

            Index cr_sizes_sum(0);
            for(Index i(0) ; i < commsize ; ++i)
              cr_sizes_sum += cr_sizes_recvbuf[i];

            Index* recvbuf(new Index[cr_sizes_sum]);

            int* recvcounts = new int[commsize];
            for(Index i(0) ; i < commsize ; ++i)
              recvcounts[i] = int(cr_sizes_recvbuf[i]);

            int* rdispls = new int[commsize];
            rdispls[0] = 0;
            Index write_count(0);
            for(Index i(1) ; i < commsize ; ++i)
            {
              write_count += cr_sizes_recvbuf[i - 1];
              rdispls[i] = int(write_count);
            }

            //Util::Comm::alltoallv(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls, synched_part.get_comm());
            synched_part.get_comm().alltoallv(sendbuf, sendcounts, sdispls, recvbuf, recvcounts, rdispls);

            std::vector<Index> p0;
            std::vector<Index> p1;
            Index start(0);
            for(Index proc(0) ; proc < commsize ; ++proc)
            {
              for(Index i(0) ; i < cr_sizes_recvbuf[proc] ; ++i)
              {
                Index other_idx(start + i);
                Index other(recvbuf[other_idx]);
                //ct is (p, other) => search for (p, other), (other, p)
                bool found(false);
                for(Index j(0) ; j < p0.size() ; ++j)
                  if((p0.at(j) == other && p1.at(j) == proc) || (p0.at(j) == proc && p1.at(j) == other) )
                  {
                    found = true;
                  }

                if(!found)
                {
                  p0.push_back(proc);
                  p1.push_back(other);
                }
              }
              start += cr_sizes_recvbuf[proc];
            }

            for(Index i(0) ; i < p0.size() ; ++i)
            {
              if(p0.at(i) == p)
              {
                synched_part.get_comm_ranks().push_back(p1.at(i));
                synched_part.get_comm_tags().push_back(i);
              }
              else if(p1.at(i) == p)
              {
                synched_part.get_comm_ranks().push_back(p0.at(i));
                synched_part.get_comm_tags().push_back(i);
              }
            }

            delete[] part_sizes_recvbuf;
            delete[] cr_sizes_recvbuf;
            delete[] recvbuf;
            delete[] sendcounts;
            delete[] sdispls;
            delete[] recvcounts;
            delete[] rdispls;

            return synched_part;
          }
      };
#endif


#ifdef FEAT_HAVE_PARMETIS
    template<typename Mode_>
    class PExecutorParmetis : public PExecutorBase<idx_t>
    {
    };

    template<>
    class PExecutorParmetis<ParmetisModePartKway> : public PExecutorBase<idx_t>
    {
      public:
        typedef PGraphParmetis PGraphT;

        static PResult part(PGraphParmetis& g)
        {
          PResult result(g.get_num_vtx(), g.get_comm());
          ///fill secondary data needed for synch later
          for(Index i(0) ; i < Index(g.get_comm().size()+1) ; ++i)
            result.get_vtxdist()[i] = g.get_vtxdist()[i];

          ///main partitioning algorithm
          if(g.ready_forexec)
          {
            // Note: ParMETIS really wants a pointer to a MPI_Comm...
            MPI_Comm mpi_comm = result.get_comm().mpi_comm();

            auto msg = ParMETIS_V3_PartKway(
                g.get_vtxdist(),
                g.get_xadj(),
                g.get_adjncy(),
                g.get_vwgt(),
                g.get_adjwgt(),
                &g.get_wgtflag(),
                &g.get_numflag(),
                &g.get_ncon(),
                &g.get_nparts(),
                g.get_tpwgts(),
                g.get_ubvec(),
                g.get_options(),
                &g.get_edgecut(),
                result.get(),
                &mpi_comm);

            if(msg != METIS_OK)
              throw InternalError("METIS IS NOT OK WITH YOU!!!");
          }
          else
            throw InternalError("EEEK");

          ///fill comm_ranks and comm_tags
          /*for(Index i(0) ; i < Index(g.get_num_vtx()) ; ++i)
            for(Index j(Index(g.get_xadj()[i])) ; j < Index(g.get_xadj()[i + 1]) ; ++j)
            {
              Index local_edge_count(0);

              Index k(Index(g.get_adjncy()[j]));
              if(!(k >= Index(g.get_vtxdist()[Util::Comm::rank(result.get_comm())]) && k < Index(g.get_vtxdist()[Util::Comm::rank(result.get_comm()) + 1])))
              {
                // other rank
                result.get_comm_ranks().push_back(k);

                // commlink global id
                Index offset_into_global(Index(g.get_vtxdist()[Util::Comm::rank(result.get_comm())]));
                Index global_dual_edge_index(Index(base.get_adjncy()[base.get_xadj()[offset_into_global + local_edge_count]]));
                result.get_comm_tags().push_back(global_dual_edge_index);
              }
              ++local_edge_count;
            }*/

          return result;
        }

        static PResult& fill_comm_structs_global(PResult& synched_part, const PGraphParmetis& base)
        {
          /// determine elems on this process (rank p)
          const Index commsize = Index(synched_part.get_comm().size());
          Index p = Index(synched_part.get_comm().rank());
          std::vector<Index> e_p;
          for(Index i(0) ; i < Index(synched_part.size()) ; ++i)
          {
            if(Index(synched_part.get()[i]) == p)
              e_p.push_back(i);
          }

          ///determine commlinks' comm ranks, i.e. (p, cr_p_i)
          std::vector<Index> cr_p;
          for(auto& e_p_i : e_p)
          {
            const Index adj_start(Index(base.get_xadj()[e_p_i]));
            const Index adj_end(Index(base.get_xadj()[e_p_i + 1]));

            for(Index i(adj_start) ; i < adj_end ; ++i)
              if(Index(synched_part.get()[ base.get_adjncy()[i] ]) != p)
              {
                //comm link found!
                cr_p.push_back(Index(synched_part.get()[Index(base.get_adjncy()[i])]));
              }
          }
          std::sort(cr_p.begin(), cr_p.end());
          auto last = std::unique(cr_p.begin(), cr_p.end());
          cr_p.erase(last, cr_p.end());

          ///synchronize and communicate commlinks, determine commtags
          //first, synchronize partition sizes
          Index* part_sizes_recvbuf(new Index[commsize]);
          Index part_sizes_sendbuf(Index(e_p.size()));

          //Util::Comm::allgather(&part_sizes_sendbuf, 1, part_sizes_recvbuf, 1, synched_part.get_comm());
          synched_part.get_comm().allgather(&part_sizes_sendbuf, std::size_t(1), part_sizes_recvbuf, std::size_t(1));

          //second, synchronize cr sizes
          Index* cr_sizes_recvbuf(new Index[commsize]);
          Index cr_sizes_sendbuf(Index(cr_p.size()));

          //Util::Comm::allgather(&cr_sizes_sendbuf, 1, cr_sizes_recvbuf, 1, synched_part.get_comm());
          synched_part.get_comm().allgather(&cr_sizes_sendbuf, std::size_t(1), cr_sizes_recvbuf, std::size_t(1));

          //then, synchronize cr
          Index* sendbuf = cr_p.data();

          Index cr_sizes_sum(0);
          for(Index i(0) ; i < commsize ; ++i)
            cr_sizes_sum += cr_sizes_recvbuf[i];

          Index* recvbuf(new Index[cr_sizes_sum]);

          int* recvcounts = new int[commsize];
          for(Index i(0) ; i < commsize ; ++i)
            recvcounts[i] = int(cr_sizes_recvbuf[i]);

          int* rdispls = new int[commsize];
          rdispls[0] = 0;
          Index write_count(0);
          for(Index i(1) ; i < commsize ; ++i)
          {
            write_count += cr_sizes_recvbuf[i - 1];
            rdispls[i] = int(write_count);
          }

          //Util::Comm::allgatherv(sendbuf, int(cr_p.size()), recvbuf, recvcounts, rdispls, synched_part.get_comm());
          synched_part.get_comm().allgatherv(sendbuf, cr_p.size(), recvbuf, recvcounts, rdispls);

          std::vector<Index> p0;
          std::vector<Index> p1;
          Index start(0);
          for(Index proc(0) ; proc < commsize ; ++proc)
          {
            for(Index i(0) ; i < cr_sizes_recvbuf[proc] ; ++i)
            {
              Index other_idx(start + i);
              Index other(recvbuf[other_idx]);
              //ct is (p, other) => search for (p, other), (other, p)
              bool found(false);
              for(Index j(0) ; j < p0.size() ; ++j)
                if((p0.at(j) == other && p1.at(j) == proc) || (p0.at(j) == proc && p1.at(j) == other) )
                {
                  found = true;
                }

              if(!found)
              {
                p0.push_back(proc);
                p1.push_back(other);
              }
            }
            start += cr_sizes_recvbuf[proc];
          }

          for(Index i(0) ; i < p0.size() ; ++i)
          {
            if(p0.at(i) == p)
            {
              synched_part.get_comm_ranks().push_back(p1.at(i));
              synched_part.get_comm_tags().push_back(i);
            }
            else if(p1.at(i) == p)
            {
              synched_part.get_comm_ranks().push_back(p0.at(i));
              synched_part.get_comm_tags().push_back(i);
            }
          }

          delete[] part_sizes_recvbuf;
          delete[] cr_sizes_recvbuf;
          delete[] recvbuf;
          delete[] recvcounts;
          delete[] rdispls;

          return synched_part;
        }
    };
#endif
  }
}
#endif
