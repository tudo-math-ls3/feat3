#pragma once
#ifndef KERNEL_FOUNDATION_PGRAPH_HPP
#define KERNEL_FOUNDATION_PGRAPH_HPP 1

#include<memory>
#include<kernel/base_header.hpp>
#include<kernel/util/comm_base.hpp>
#include<kernel/geometry/conformal_mesh.hpp>
#include<kernel/adjacency/graph.hpp>
#include<kernel/shape.hpp>

#include<utility>

#ifdef FEAST_HAVE_PARMETIS
FEAST_DISABLE_WARNINGS
#include<parmetis.h>
FEAST_RESTORE_WARNINGS
#endif

namespace FEAST
{
  using namespace Geometry;
  using namespace Shape;

  namespace Foundation
  {
    ///base class for (distributed) graphs used for partitioning
    template<typename IT_>
    class PGraphBase
    {
      public:
        typedef IT_ IndexType;

        PGraphBase(IndexType num_local_dualvtx, IndexType num_local_dualedges, IndexType num_global_vtx, const Communicator& comm) :
          ready_forexec(false),
          ready_fordist(false),
          _num_local_dualvtx(num_local_dualvtx),
          _num_local_dualedges(num_local_dualedges),
          _vtxdist(new IndexType[Index(Comm::size(comm) + 1)]),
          _xadj(new IndexType[Index(num_local_dualvtx + 1)]),
          _adjncy(new IndexType[Index(2 * num_local_dualedges)]),
          _part(new IndexType[Index(num_local_dualvtx)]),
          _comm(comm)
        {
          IndexType k = IndexType(Comm::size(comm));
          IndexType vtxdist_end(k + 1);
          IndexType n_local_avg((num_global_vtx - (num_global_vtx % k))/k);
          for(IndexType i(0) ; i < vtxdist_end - 1 ; ++i)
            _vtxdist[i] = i * n_local_avg;

          _vtxdist[k] = (vtxdist_end - 1) * n_local_avg + (num_global_vtx % k);
        }

        PGraphBase(PGraphBase&& other) :
          ready_forexec(other.ready_forexec),
          ready_fordist(other.ready_fordist),
          _num_local_dualvtx(other._num_local_dualvtx),
          _num_local_dualedges(other._num_local_dualedges),
          _vtxdist(other._vtxdist),
          _xadj(other._xadj),
          _adjncy(other._adjncy),
          _part(other._part),
          _comm(other._comm)
        {
          other.ready_forexec = false;
          other.ready_fordist = false;
          other._num_local_dualvtx = 0;
          other._num_local_dualedges = 0;
          other._vtxdist = nullptr;
          other._xadj = nullptr;
          other._adjncy = nullptr;
          other._part = nullptr;
        }

        PGraphBase& operator=(PGraphBase&& other)
        {
          if(this == &other)
            return *this;

          this->ready_forexec = other.ready_forexec;
          this->ready_fordist = other.ready_fordist;
          this->_num_localdualvtx = other._num_localdualvtx;
          this->_num_localdualedges = other._num_localdualedges;
          this->_vtxdist = other._vtxdist;
          this->_xadj = other._xadj;
          this->_adjncy = other._adjncy;
          this->_part = other._part;
          this->_comm = other._comm;

          other.ready_forexec = false;
          other.ready_fordist = false;
          other._num_local_dualvtx = 0;
          other._num_local_dualedges = 0;
          other._vtxdist = nullptr;
          other._xadj = nullptr;
          other._adjncy = nullptr;
          other._part = nullptr;

          return *this;
        }

        virtual ~PGraphBase()
        {
          delete[] _vtxdist;
          delete[] _xadj;
          delete[] _adjncy;
          delete[] _part;
        }

        ///things all FEAST-conformal pgraphs must minimally provide

        ///number of elements
        virtual IndexType get_num_vtx() const = 0;

        ///initial distribution of elements
        virtual IndexType* get_vtxdist() = 0;
        virtual IndexType* get_vtxdist() const = 0;

        ///adjncy[ [xadj[i], xadj[i+1] ] stores the adjacency list of dual-graph vtx i, vertex indices are global
        virtual IndexType* get_xadj() = 0;
        virtual IndexType* get_xadj() const= 0;
        virtual IndexType* get_adjncy() = 0;
        virtual IndexType* get_adjncy() const = 0;

        ///part[i] contains the target process index of dual-graph vtx i
        virtual IndexType* get_part() = 0;
        virtual IndexType* get_part() const = 0;

        ///communicator used for (parallel) partitioning
        virtual Communicator get_comm() const = 0;

        bool ready_forexec;
        bool ready_fordist;

        PGraphBase() :
#ifndef SERIAL
          _comm(Communicator(MPI_COMM_WORLD))
#else
          _comm(Communicator(0))
#endif
        {
        }

        virtual void reset(IndexType* xadj, IndexType* adjncy, IndexType new_num_dualvtx, IndexType new_num_dualedges) = 0;

        virtual std::shared_ptr<PGraphBase> create_local() const = 0;

      protected:

        IndexType _num_local_dualvtx;
        IndexType _num_local_dualedges;
        IndexType* _vtxdist;
        IndexType* _xadj;
        IndexType* _adjncy;
        IndexType* _part;
        Communicator _comm;
    };

    template<typename DT_, typename IT_>
    class PGraphNONE : public PGraphBase<IT_>
    {
      public:
        typedef IT_ IndexType;
        typedef DT_ DataType;

        PGraphNONE(IndexType num_local_dualvtx, IndexType num_local_dualedges, IndexType num_globalvtx, const Communicator& comm) :
          PGraphBase<IndexType>(num_local_dualvtx, num_local_dualedges, num_globalvtx, comm)
        {
        }

        PGraphNONE(PGraphNONE&& other) :
          PGraphBase<IndexType>(std::forward<PGraphNONE>(other))
        {
        }

        PGraphNONE& operator=(PGraphNONE&& other)
        {
          if(this == &other)
            return *this;

          this->ready_forexec = other.ready_forexec;
          this->ready_fordist = other.ready_fordist;
          this->_num_local_dualvtx = other._num_local_dualvtx;
          this->_num_local_dualedges = other._num_local_dualedges;
          this->_vtxdist = other._vtxdist;
          this->_xadj = other._xadj;
          this->_adjncy = other._adjncy;
          this->_part = other._part;
          this->_comm = other._comm;

          other.ready_forexec = false;
          other.ready_fordist = false;
          other._num_local_dualvtx = 0;
          other._num_local_dualedges = 0;
          other._vtxdist = nullptr;
          other._xadj = nullptr;
          other._adjncy = nullptr;
          other._part = nullptr;

          return *this;
        }

        PGraphNONE clone() const
        {
          PGraphNONE result;

          result.ready_forexec = this->ready_forexec;
          result.ready_fordist = this->ready_fordist;
          result._num_local_dualvtx = this->_num_local_dualvtx;
          result._num_local_dualedges = this->_num_local_dualedges;

          result._vtxdist = new IndexType[Index(Comm::size(this->_comm) + 1)];
          for(Index i(0) ; i < Index(Comm::size(this->_comm) + 1) ; ++i)
            result._vtxdist[i] = this->_vtxdist[i];

          result._xadj = new IndexType[Index(this->_num_local_dualvtx + 1)];
          for(Index i(0) ; i < Index(this->_num_local_dualvtx + 1) ; ++i)
            result._xadj[i] = this->_xadj[i];

          result._adjncy = new IndexType[Index(2 * this->_num_local_dualedges)];
          for(Index i(0) ; i < Index(2 * this->_num_local_dualedges) ; ++i)
            result._adjncy[i] = this->_adjncy[i];

          result._part = new IndexType[Index(this->_num_local_dualvtx)];
          for(Index i(0) ; i < Index(this->_num_local_dualvtx) ; ++i)
            result._part[i] = this->_part[i];

          result._comm = this->_comm;

          return result;
        }

        template<typename SourceMeshT_>
        PGraphNONE(const SourceMeshT_& source, Index num_globalvtx, Communicator comm)
        {
          ///get sizes of the dual graph to create
          IndexType num_dualvtx(IndexType(source.get_num_entities(SourceMeshT_::shape_dim)));
          IndexType num_vtx(IndexType(source.get_num_entities(0)));

          //build up the graph
          auto& source_vertex_at_elem(source.template get_index_set<SourceMeshT_::shape_dim, 0>());
          Index vtxcount_at_elem(FaceTraits<typename SourceMeshT_::ShapeType, 0>::count);
          Index* xadj(new Index[Index(num_dualvtx) * vtxcount_at_elem]);
          for(IndexType i(0) ; i < num_dualvtx + 1; ++i)
          {
            xadj[i] = vtxcount_at_elem * Index(i);
          }
          Index* adjncy(new Index[Index(num_dualvtx) * vtxcount_at_elem]);
          for(IndexType i(0), j(0) ; i < num_dualvtx ; ++i, j += IndexType(vtxcount_at_elem))
          {
            for(IndexType k(0) ; k < IndexType(vtxcount_at_elem) ; ++k)
            {
              adjncy[j + k] = source_vertex_at_elem[Index(i)][Index(k)];
            }
          }

          const Index* xadjc(xadj);
          const Index* adjncyc(adjncy);

          Adjacency::Graph g((Index)num_dualvtx, Index(num_vtx), (Index)num_dualvtx * vtxcount_at_elem, xadjc, nullptr, adjncyc);
          Adjacency::Graph gt(Adjacency::rt_transpose, g);
          Adjacency::Graph n(Adjacency::rt_injectify, g, gt);

          *this = std::move(PGraphNONE<DataType, IndexType>(IndexType(n.get_num_nodes_domain()), IndexType(n.get_num_indices()), IndexType(num_globalvtx), comm));
          auto target_xadj(this->get_xadj());
          auto target_adjncy(this->get_adjncy());
          auto source_xadj(n.get_domain_ptr());
          auto source_adjncy(n.get_image_idx());

          for(IndexType i(0) ; i < num_dualvtx + 1 ; ++i)
          {
            target_xadj[i] = IndexType(source_xadj[i]);
          }
          for(IndexType i(0) ; i < IndexType(n.get_num_indices()) ; ++i)
          {
            target_adjncy[i] = IndexType(source_adjncy[i]);
          }

          delete[] xadj;
          delete[] adjncy;

          this->ready_forexec = true;
        }

        void reset(IndexType* xadj, IndexType* adjncy, IndexType new_num_dualvtx, IndexType new_num_dualedges) override
        {
          delete[] this->_xadj;
          delete[] this->_adjncy;
          delete[] this->_part;

          this->_num_local_dualvtx = new_num_dualvtx;
          this->_num_local_dualedges = new_num_dualedges;

          this->_xadj = xadj;
          this->_adjncy = adjncy;
          this->_part = new IndexType[Index(new_num_dualvtx)];
        }

        std::shared_ptr<PGraphBase<IndexType> > create_local() const override
        {
          PGraphNONE<DataType, IndexType> result;
          result = this->clone();
          result.ready_forexec = true;

          return std::shared_ptr<PGraphBase<IndexType> >(new PGraphNONE<DataType, IndexType>(std::move(result)));
        }

        virtual IndexType get_num_vtx() const override
        {
          return this->_num_local_dualvtx;
        }

        virtual IndexType* get_vtxdist() override
        {
          return this->_vtxdist;
        }

        virtual IndexType* get_vtxdist() const override
        {
          return this->_vtxdist;
        }

        virtual IndexType* get_xadj() override
        {
          return this->_xadj;
        }

        virtual IndexType* get_xadj() const override
        {
          return this->_xadj;
        }

        virtual IndexType* get_adjncy() override
        {
          return this->_adjncy;
        }

        virtual IndexType* get_adjncy() const override
        {
          return this->_adjncy;
        }

        virtual IndexType* get_part() override
        {
          return this->_part;
        }

        virtual IndexType* get_part() const override
        {
          return this->_part;
        }

        virtual Communicator get_comm()
        {
          return this->_comm;
        }

        virtual Communicator get_comm() const override
        {
          return this->_comm;
        }

        virtual ~PGraphNONE()
        {
        }

        PGraphNONE()
        {
        }
    };

    template<typename DT_, typename IT_>
    class PGraphFallback : public PGraphBase<IT_>
    {
      public:
        typedef IT_ IndexType;
        typedef DT_ DataType;

        PGraphFallback(IndexType num_local_dualvtx, IndexType num_local_dualedges, IndexType num_globalvtx, const Communicator& comm) :
          PGraphBase<IndexType>(num_local_dualvtx, num_local_dualedges, num_globalvtx, comm)
        {
        }

        PGraphFallback(PGraphFallback&& other) :
          PGraphBase<IndexType>(std::forward<PGraphFallback>(other))
        {
        }

        PGraphFallback& operator=(PGraphFallback&& other)
        {
          if(this == &other)
            return *this;

          this->ready_forexec = other.ready_forexec;
          this->ready_fordist = other.ready_fordist;
          this->_num_local_dualvtx = other._num_local_dualvtx;
          this->_num_local_dualedges = other._num_local_dualedges;
          this->_vtxdist = other._vtxdist;
          this->_xadj = other._xadj;
          this->_adjncy = other._adjncy;
          this->_part = other._part;
          this->_comm = other._comm;

          other.ready_forexec = false;
          other.ready_fordist = false;
          other._num_local_dualvtx = 0;
          other._num_local_dualedges = 0;
          other._vtxdist = nullptr;
          other._xadj = nullptr;
          other._adjncy = nullptr;
          other._part = nullptr;

          return *this;
        }

        PGraphFallback clone() const
        {
          PGraphFallback result;

          result.ready_forexec = this->ready_forexec;
          result.ready_fordist = this->ready_fordist;
          result._num_local_dualvtx = this->_num_local_dualvtx;
          result._num_local_dualedges = this->_num_local_dualedges;

          result._vtxdist = new IndexType[Index(Comm::size(this->_comm) + 1)];
          for(Index i(0) ; i < Index(Comm::size(this->_comm) + 1) ; ++i)
            result._vtxdist[i] = this->_vtxdist[i];

          result._xadj = new IndexType[Index(this->_num_local_dualvtx + 1)];
          for(Index i(0) ; i < Index(this->_num_local_dualvtx + 1) ; ++i)
            result._xadj[i] = this->_xadj[i];

          result._adjncy = new IndexType[Index(2 * this->_num_local_dualedges)];
          for(Index i(0) ; i < Index(2 * this->_num_local_dualedges) ; ++i)
            result._adjncy[i] = this->_adjncy[i];

          result._part = new IndexType[Index(this->_num_local_dualvtx)];
          for(Index i(0) ; i < Index(this->_num_local_dualvtx) ; ++i)
            result._part[i] = this->_part[i];

          result._comm = this->_comm;

          return result;
        }

        template<typename SourceMeshT_>
        PGraphFallback(const SourceMeshT_& source, Index num_globalvtx, Communicator comm)
        {
          ///get sizes of the dual graph to create
          IndexType num_dualvtx(IndexType(source.get_num_entities(SourceMeshT_::shape_dim)));
          IndexType num_vtx(IndexType(source.get_num_entities(0)));

          //build up the graph
          auto& source_vertex_at_elem(source.template get_index_set<SourceMeshT_::shape_dim, 0>());
          Index vtxcount_at_elem(FaceTraits<typename SourceMeshT_::ShapeType, 0>::count);
          Index* xadj(new Index[Index(num_dualvtx) * vtxcount_at_elem]);
          for(IndexType i(0) ; i < num_dualvtx + 1; ++i)
          {
            xadj[i] = vtxcount_at_elem * Index(i);
          }
          Index* adjncy(new Index[Index(num_dualvtx) * vtxcount_at_elem]);
          for(IndexType i(0), j(0) ; i < num_dualvtx ; ++i, j += IndexType(vtxcount_at_elem))
          {
            for(IndexType k(0) ; k < IndexType(vtxcount_at_elem) ; ++k)
            {
              adjncy[j + k] = source_vertex_at_elem[Index(i)][Index(k)];
            }
          }

          const Index* xadjc(xadj);
          const Index* adjncyc(adjncy);

          Adjacency::Graph g((Index)num_dualvtx, Index(num_vtx), (Index)num_dualvtx * vtxcount_at_elem, xadjc, nullptr, adjncyc);
          Adjacency::Graph gt(Adjacency::rt_transpose, g);
          Adjacency::Graph n(Adjacency::rt_injectify, g, gt);

          *this = std::move(PGraphFallback<DataType, IndexType>(IndexType(n.get_num_nodes_domain()), IndexType(n.get_num_indices()), IndexType(num_globalvtx), comm));
          auto target_xadj(this->get_xadj());
          auto target_adjncy(this->get_adjncy());
          auto source_xadj(n.get_domain_ptr());
          auto source_adjncy(n.get_image_idx());

          for(IndexType i(0) ; i < num_dualvtx + 1 ; ++i)
          {
            target_xadj[i] = IndexType(source_xadj[i]);
          }
          for(IndexType i(0) ; i < IndexType(n.get_num_indices()) ; ++i)
          {
            target_adjncy[i] = IndexType(source_adjncy[i]);
          }

          delete[] xadj;
          delete[] adjncy;

          this->ready_forexec = true;
        }

        void reset(IndexType* xadj, IndexType* adjncy, IndexType new_num_dualvtx, IndexType new_num_dualedges) override
        {
          delete[] this->_xadj;
          delete[] this->_adjncy;
          delete[] this->_part;

          this->_num_local_dualvtx = new_num_dualvtx;
          this->_num_local_dualedges = new_num_dualedges;

          this->_xadj = xadj;
          this->_adjncy = adjncy;
          this->_part = new IndexType[Index(new_num_dualvtx)];
        }

        std::shared_ptr<PGraphBase<IndexType> > create_local() const override
        {
          ///get sizes of the dual graph to create
          Index me(Comm::rank(this->get_comm()));
          Index start(Index(this->get_vtxdist()[me]));
          Index end(Index(this->get_vtxdist()[me + 1]));
          Index num_local_dualvtx(end - start);

          Index global_node_number(start);

          Index num_local_dualedges(0);
          std::vector<std::vector<IndexType> > g_temp((Index)(num_local_dualvtx));
          for(Index i((Index)(0)) ; i < (Index)(num_local_dualvtx) ; ++i)
          {
            Index s(Index(this->get_xadj()[start + i]));
            Index e(Index(this->get_xadj()[start + i + 1]));

            for(Index j(s) ; j < e ; ++j)
            {
              if(Index(this->get_adjncy()[j]) != global_node_number) //prevent parmetis self loops
              {
                g_temp.at(i).push_back(this->get_adjncy()[j]);
              }
            }
            ++global_node_number;
          }

          ///determine xadj and copy data
          IndexType* target_xadj(new IndexType[Index(num_local_dualvtx + 1)]);

          target_xadj[0] = 0;
          Index writecount(0);
          for(Index i(1) ; i < g_temp.size() + 1 ; ++i)
          {
            writecount += Index(g_temp.at(i - 1).size());
            target_xadj[i] = IndexType(writecount);
          }

          std::vector<IndexType> reduced_g_temp;
          for(Index i(0) ; i < g_temp.size() ; ++i)
          {
            reduced_g_temp.insert(reduced_g_temp.end(), g_temp.at(i).begin(), g_temp.at(i).end());
          }
          IndexType* target_adjncy(new IndexType[Index(reduced_g_temp.size())]); //now we have exactly m_i entries (as opposed to global, where we have 2*m_global)
          //std::copy(reduced_g_temp.begin(), reduced_g_temp.end(), target_adjncy);
          for(Index i(0) ; i < reduced_g_temp.size() ; ++i)
            target_adjncy[i] = reduced_g_temp.at(i);

          PGraphFallback<DataType, IndexType> result;
          result = this->clone();
          result.ready_forexec = true;
          result.reset(target_xadj, target_adjncy, IndexType(num_local_dualvtx), IndexType(num_local_dualedges));

          return std::shared_ptr<PGraphBase<IndexType> >(new PGraphFallback<DataType, IndexType>(std::move(result)));
        }

        virtual IndexType get_num_vtx() const override
        {
          return this->_num_local_dualvtx;
        }

        virtual IndexType* get_vtxdist() override
        {
          return this->_vtxdist;
        }

        virtual IndexType* get_vtxdist() const override
        {
          return this->_vtxdist;
        }

        virtual IndexType* get_xadj() override
        {
          return this->_xadj;
        }

        virtual IndexType* get_xadj() const override
        {
          return this->_xadj;
        }

        virtual IndexType* get_adjncy() override
        {
          return this->_adjncy;
        }

        virtual IndexType* get_adjncy() const override
        {
          return this->_adjncy;
        }

        virtual IndexType* get_part() override
        {
          return this->_part;
        }

        virtual IndexType* get_part() const override
        {
          return this->_part;
        }

        virtual Communicator get_comm()
        {
          return this->_comm;
        }

        virtual Communicator get_comm() const override
        {
          return this->_comm;
        }

        virtual ~PGraphFallback()
        {
        }

        PGraphFallback()
        {
        }
    };

#ifdef FEAST_HAVE_PARMETIS
    class PGraphParmetis : public PGraphBase<idx_t>
    {
      public:
        typedef PGraphBase<idx_t> BaseClass;
        typedef real_t DataType;

        PGraphParmetis(IndexType num_local_dualvtx, IndexType num_local_dualedges, IndexType num_globalvtx, const Communicator& comm) :
          PGraphBase(num_local_dualvtx, num_local_dualedges, num_globalvtx, comm),
          _vwgt(new IndexType[Index(num_local_dualvtx)]),
          _adjwgt(new IndexType[Index(2 * num_local_dualedges)]),
          _options(new IndexType[Index(3)]),
          _ncon(1),
          _nparts(IndexType(Comm::size(comm))),
          _wgtflag(0),
          _numflag(0),
          _tpwgts(new DataType[Index(_ncon * _nparts)]),
          _ubvec(new DataType[Index(_ncon)])
        {
          _options[0] = IndexType(0);

          for(Index i(0) ; i < Index(_ncon * _nparts) ; ++i)
            _tpwgts[i] = DataType(1. / DataType(_nparts));

          _ubvec[0] = DataType(1.05);
        }

        PGraphParmetis(PGraphParmetis&& other) :
          PGraphBase(std::forward<PGraphParmetis>(other)),
          _vwgt(other._vwgt),
          _adjwgt(other._adjwgt),
          _options(other._options),
          _ncon(other._ncon),
          _nparts(other._nparts),
          _wgtflag(other._wgtflag),
          _numflag(other._numflag),
          _tpwgts(other._tpwgts),
          _ubvec(other._ubvec)
        {
          other._vwgt = nullptr;
          other._adjwgt = nullptr;
          other._options = nullptr;
          other._tpwgts = nullptr;
          other._ubvec = nullptr;
        }

        PGraphParmetis& operator=(PGraphParmetis&& other)
        {
          if(this == &other)
            return *this;

          this->ready_forexec = other.ready_forexec;
          this->ready_fordist = other.ready_fordist;
          this->_num_local_dualvtx = other._num_local_dualvtx;
          this->_num_local_dualedges = other._num_local_dualedges;
          this->_vtxdist = other._vtxdist;
          this->_xadj = other._xadj;
          this->_adjncy = other._adjncy;
          this->_part = other._part;
          this->_comm = other._comm;

          this->_vwgt = other._vwgt;
          this->_adjwgt = other._adjwgt;
          this->_options = other._options;
          this->_ncon = other._ncon;
          this->_nparts = other._nparts;
          this->_wgtflag = other._wgtflag;
          this->_numflag = other._numflag;
          this->_tpwgts = other._tpwgts;
          this->_ubvec = other._ubvec;

          other.ready_forexec = false;
          other.ready_fordist = false;
          other._num_local_dualvtx = 0;
          other._num_local_dualedges = 0;
          other._vtxdist = nullptr;
          other._xadj = nullptr;
          other._adjncy = nullptr;
          other._part = nullptr;

          other._vwgt = nullptr;
          other._adjwgt = nullptr;
          other._options = nullptr;
          other._tpwgts = nullptr;
          other._ubvec = nullptr;

          return *this;
        }

        PGraphParmetis clone() const
        {
          PGraphParmetis result;

          result.ready_forexec = ready_forexec;
          result.ready_fordist = ready_fordist;
          result._num_local_dualvtx = _num_local_dualvtx;
          result._num_local_dualedges = _num_local_dualedges;

          result._vtxdist = new IndexType[Index(Comm::size(_comm) + 1)];
          for(Index i(0) ; i < Index(Comm::size(_comm) + 1) ; ++i)
            result._vtxdist[i] = _vtxdist[i];

          result._xadj = new IndexType[Index(_num_local_dualvtx + 1)];
          for(Index i(0) ; i < Index(_num_local_dualvtx + 1) ; ++i)
            result._xadj[i] = _xadj[i];

          result._adjncy = new IndexType[Index(2 * _num_local_dualedges)];
          for(Index i(0) ; i < Index(2 * _num_local_dualedges) ; ++i)
            result._adjncy[i] = _adjncy[i];

          result._part = new IndexType[Index(_num_local_dualvtx)];
          for(Index i(0) ; i < Index(_num_local_dualvtx) ; ++i)
            result._part[i] = _part[i];

          result._comm = _comm;

          result._vwgt = new IndexType[Index(_num_local_dualvtx)];
          for(Index i(0) ; i < Index(_num_local_dualvtx) ; ++i)
            result._vwgt[i] = _vwgt[i];

          result._adjwgt = new IndexType[Index(2 * _num_local_dualedges)];
          for(Index i(0) ; i < Index(2 * _num_local_dualedges) ; ++i)
            result._adjwgt[i] = _adjwgt[i];

          result._options = new IndexType[Index(3)];
          for(Index i(0) ; i < Index(3) ; ++i)
            result._options[i] = _options[i];

          result._ncon = _ncon;
          result._nparts = _nparts;
          result._wgtflag = _wgtflag;
          result._numflag = _numflag;

          result._tpwgts = new DataType[Index(_ncon * _nparts)];
          for(Index i(0) ; i < Index(_ncon * _nparts); ++i)
            result._tpwgts[i] = _tpwgts[i];

          result._ubvec = new DataType[Index(_ncon)];
          for(Index i(0) ; i < Index(_ncon); ++i)
            result._ubvec[i] = _ubvec[i];

          return result;
        }

        template<typename SourceMeshT_>
        PGraphParmetis(const SourceMeshT_& source, Index num_globalvtx, Communicator comm)
        {
          ///get sizes of the dual graph to create
          IndexType num_dualvtx(IndexType(source.get_num_entities(SourceMeshT_::shape_dim)));
          IndexType num_vtx(IndexType(source.get_num_entities(0)));

          //build up the graph
          auto& source_vertex_at_elem(source.template get_index_set<SourceMeshT_::shape_dim, 0>());
          Index vtxcount_at_elem(FaceTraits<typename SourceMeshT_::ShapeType, 0>::count);
          Index* xadj(new Index[Index(num_dualvtx) * vtxcount_at_elem]);
          for(IndexType i(0) ; i < num_dualvtx + 1; ++i)
          {
            xadj[i] = vtxcount_at_elem * Index(i);
          }
          Index* adjncy(new Index[Index(num_dualvtx) * vtxcount_at_elem]);
          for(IndexType i(0), j(0) ; i < num_dualvtx ; ++i, j += IndexType(vtxcount_at_elem))
          {
            for(IndexType k(0) ; k < IndexType(vtxcount_at_elem) ; ++k)
            {
              adjncy[j + k] = source_vertex_at_elem[Index(i)][Index(k)];
            }
          }

          const Index* xadjc(xadj);
          const Index* adjncyc(adjncy);

          Adjacency::Graph g((Index)num_dualvtx, Index(num_vtx), (Index)num_dualvtx * vtxcount_at_elem, xadjc, nullptr, adjncyc);
          Adjacency::Graph gt(Adjacency::rt_transpose, g);
          Adjacency::Graph n(Adjacency::rt_injectify, g, gt);

          *this = PGraphParmetis(IndexType(n.get_num_nodes_domain()), IndexType(n.get_num_indices()), IndexType(num_globalvtx), comm);
          auto target_xadj(this->get_xadj());
          auto target_adjncy(this->get_adjncy());
          auto source_xadj(n.get_domain_ptr());
          auto source_adjncy(n.get_image_idx());

          for(IndexType i(0) ; i < num_dualvtx + 1 ; ++i)
          {
            target_xadj[i] = IndexType(source_xadj[i]);
          }
          for(IndexType i(0) ; i < IndexType(n.get_num_indices()) ; ++i)
          {
            target_adjncy[i] = IndexType(source_adjncy[i]);
          }

          delete[] xadj;
          delete[] adjncy;

          this->ready_forexec = true;
        }

        void reset(IndexType* xadj, IndexType* adjncy, IndexType new_num_dualvtx, IndexType new_num_dualedges) override
        {
          delete[] _xadj;
          delete[] _adjncy;
          delete[] _part;
          delete[] _vwgt;
          delete[] _adjwgt;

          _num_local_dualvtx = new_num_dualvtx;
          _num_local_dualedges = new_num_dualedges;

          _xadj = xadj;
          _adjncy = adjncy;
          _part = new IndexType[Index(new_num_dualvtx)];
          _vwgt = new IndexType[Index(new_num_dualvtx)];
          _adjwgt = new IndexType[Index(2 * new_num_dualedges)];
        }

        std::shared_ptr<PGraphBase> create_local() const override
        {
          ///get sizes of the dual graph to create
          Index me(Comm::rank(this->get_comm()));
          Index start(Index(this->get_vtxdist()[me]));
          Index end(Index(this->get_vtxdist()[me + 1]));
          Index num_local_dualvtx(end - start);

          Index global_node_number(start);

          Index num_local_dualedges(0);
          std::vector<std::vector<IndexType> > g_temp((Index)(num_local_dualvtx));
          for(Index i((Index)(0)) ; i < (Index)(num_local_dualvtx) ; ++i)
          {
            Index s(Index(this->get_xadj()[start + i]));
            Index e(Index(this->get_xadj()[start + i + 1]));

            for(Index j(s) ; j < e ; ++j)
            {
              if(Index(this->get_adjncy()[j]) != global_node_number) //prevent parmetis self loops
              {
                g_temp.at(i).push_back(this->get_adjncy()[j]);
              }
            }
            ++global_node_number;
          }

          ///determine xadj and copy data
          IndexType* target_xadj(new IndexType[Index(num_local_dualvtx + 1)]);

          target_xadj[0] = 0;
          Index writecount(0);
          for(Index i(1) ; i < g_temp.size() + 1 ; ++i)
          {
            writecount += Index(g_temp.at(i - 1).size());
            target_xadj[i] = IndexType(writecount);
          }

          std::vector<IndexType> reduced_g_temp;
          for(Index i(0) ; i < g_temp.size() ; ++i)
          {
            reduced_g_temp.insert(reduced_g_temp.end(), g_temp.at(i).begin(), g_temp.at(i).end());
          }
          IndexType* target_adjncy(new IndexType[Index(reduced_g_temp.size())]); //now we have exactly m_i entries (as opposed to global, where we have 2*m_global)
          //std::copy(reduced_g_temp.begin(), reduced_g_temp.end(), target_adjncy);
          for(Index i(0) ; i < reduced_g_temp.size() ; ++i)
            target_adjncy[i] = reduced_g_temp.at(i);


          PGraphParmetis result;
          result = this->clone();
          result.ready_forexec = true;
          result.reset(target_xadj, target_adjncy, IndexType(num_local_dualvtx), IndexType(num_local_dualedges));

          return std::shared_ptr<PGraphBase>(new PGraphParmetis(std::move(result)));
        }

        virtual IndexType get_num_vtx() const override
        {
          return this->_num_local_dualvtx;
        }

        virtual IndexType* get_vtxdist() override
        {
          return this->_vtxdist;
        }

        virtual IndexType* get_vtxdist() const override
        {
          return this->_vtxdist;
        }

        virtual IndexType* get_xadj() override
        {
          return this->_xadj;
        }

        virtual IndexType* get_xadj() const override
        {
          return this->_xadj;
        }

        virtual IndexType* get_adjncy() override
        {
          return this->_adjncy;
        }

        virtual IndexType* get_adjncy() const override
        {
          return this->_adjncy;
        }

        virtual IndexType* get_part() override
        {
          return this->_part;
        }

        virtual IndexType* get_part() const override
        {
          return this->_part;
        }

        virtual Communicator get_comm()
        {
          return this->_comm;
        }

        virtual Communicator get_comm() const override
        {
          return this->_comm;
        }

        virtual ~PGraphParmetis()
        {
          delete[] _vwgt;
          delete[] _adjwgt;
          delete[] _options;
          delete[] _tpwgts;
          delete[] _ubvec;
        }

        IndexType* get_vwgt()
        {
          return this->_vwgt;
        }

        IndexType* get_adjwgt()
        {
          return this->_adjwgt;
        }

        IndexType& get_ncon()
        {
          return this->_ncon;
        }

        IndexType& get_nparts()
        {
          return this->_nparts;
        }

        DataType* get_tpwgts()
        {
          return this->_tpwgts;
        }

        DataType* get_ubvec()
        {
          return this->_ubvec;
        }

        IndexType& get_wgtflag()
        {
          return this->_wgtflag;
        }

        IndexType& get_numflag()
        {
          return this->_numflag;
        }

        IndexType& get_edgecut()
        {
          return this->_edgecut;
        }

        IndexType* get_options()
        {
          return this->_options;
        }

        PGraphParmetis()
        {
        }

      private:

        IndexType* _vwgt;
        IndexType* _adjwgt;

        IndexType* _options;

        IndexType _ncon;
        IndexType _nparts;

        IndexType _wgtflag;
        IndexType _numflag;
        IndexType _edgecut;

        DataType* _tpwgts; //size: ncon * nparts
        DataType* _ubvec; //size ncon
    };
#endif
  }
}
#endif
