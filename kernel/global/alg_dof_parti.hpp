// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GLOBAL_ALG_DOF_PARTI_HPP
#define KERNEL_GLOBAL_ALG_DOF_PARTI_HPP 1

#include <kernel/adjacency/dynamic_graph.hpp>
#include <kernel/global/gate.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/global/filter.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/matrix_mirror.hpp>

#include <algorithm>
#include <vector>

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Algebraic DOF Partitioning implementation class template
     *
     * See the documentation of the specialization of this class template for details.
     */
    template<typename LocalVector_, typename Mirror_>
    class AlgDofParti;

    /**
     * \brief Algebraic DOF Partitioning implementation
     *
     * This class implements the functionality to generate and manage a distributed partitioning
     * of the degrees of freedom of a (linear) system in an algebraic sense instead of the
     * usual domain-decomposition approach. With this algebraic dof partitioning approach, each
     * process in a given communicator becomes the unqiue owner of a set of consecutively enumerated
     * global dofs, whereas in the domain-decomposition approach a single global dof may be shared
     * by several processes, thus having no unique owner.
     * The primary purpose of this class and the AlgDofPartiVector and AlgDofPartiMatrix classes
     * is to offer the capability of embedding various third-party linear solver libraries as
     * solvers, preconditioners or smoothers within the solver framework used in FEAT.
     *
     * To use this class, one needs to have a fully assembled Global::Gate object, which contains
     * all the partitioning information in the usual domain-decomposition sense, which is used as
     * a "starting point" for an algebraic DOF partitioning in the assemble_by_gate() function.
     *
     * Throughout this class and the AlgDofPartiVector and AlgDofPartiMatrix classes, you will
     * often find the terms "local DOF", "shared DOF" and "owned DOF" and "donee DOF".
     * The definitions are as followed:
     * - "local" DOFs always means the DOFs as they are used in the domain decomposition world.
     * - "shared" DOFs are all local DOFs are "shared" by at least two processes in the domain
     *   decomposition world. In other words: Each local DOF, which is part of at least one
     *   neighbor mirror in the Global::Gate, is a shared DOF.
     * - "owned" DOFs are all local DOFs, which:
     *   - are not shared with any other process
     *   - are only shared with processes of higher rank
     * - "donee" DOFs are all shared DOFs, which are not owned DOFs
     *
     * Furthermore, the neighbor ranks and mirrors, which are used in the Global::Gate of the
     * domain-decomposition based approach, are converted into "owner" and "donee" neighbor
     * ranks and mirrors in this algebraic approach. The basic idea is the following: Assume
     * that three processes A, B and C share a common set of DOFs, i.e. they are neighbors in
     * the domain-decomposition sense, and assume that process A has a lower rank than process B
     * and process C. Due to the "ownership" rule (lowest rank is owner), process A becomes the
     * "owner" of the shared DOFs and thus is responsible for enumerating and managing all DOFs
     * shared by A, B and C, which  in turn makes processes B and C "donees". Now, because
     * process A owns at least one DOF that is also a local DOF of process B and a local DOF
     * of process C in the domain decomposition sense, we have that:
     * - process A is an "owner neighbor" of process B
     * - process A is an "owner neighbor" of process C
     * - process B is a "donee neighbor" of process A
     * - process C is a "donee neighbor" of process A
     *
     * \todo decide what happens if this process has no dofs at all
     * \todo make sure that this works for discontinuous elements, too
     *
     * \author Peter Zajac
     */
    template<typename DT_, typename IT_>
    class AlgDofParti<LAFEM::DenseVector<Mem::Main, DT_, IT_>, LAFEM::VectorMirror<Mem::Main, DT_, IT_>>
    {
    public:
      /// our memory tag type
      typedef Mem::Main MemType;
      /// our data type
      typedef DT_ DataType;
      /// our index type
      typedef IT_ IndexType;

      /// the local vector type
      typedef LAFEM::DenseVector<MemType, DataType, IndexType> LocalVectorType;
      /// the vector mirror type
      typedef LAFEM::VectorMirror<MemType, DataType, IndexType> MirrorType;
      /// the global vector type
      typedef Global::Vector<LocalVectorType, MirrorType> GlobalVectorType;
      /// the global gate type
      typedef Global::Gate<LocalVectorType, MirrorType> GateType;

    protected:
      /// our communicator
      const Dist::Comm* _comm;
      /// global dof offset of this process
      Index _glob_dof_offset;
      /// global dof count over all processes
      Index _glob_dof_count;
      /// global dof indices, size = number of local DOFs
      std::vector<Index> _glob_dof_idx;
      /// mirror for this process's owned DOFs
      MirrorType _owned_mirror;
      /// rank/mirror-pair of DOF-owner processes
      std::vector<std::pair<int, MirrorType>> _neighbors_owner;
      /// rank/mirror-pair of DOF-donee processes
      std::vector<std::pair<int, MirrorType>> _neighbors_donee;

      /// global DOF offsets of all processes, \see #allgather()
      std::vector<int> _all_glob_dof_offset;
      /// global DOF counts of all processes, \see #allgather()
      std::vector<int> _all_glob_dof_counts;

    public:
      AlgDofParti() :
        _comm(nullptr),
        _glob_dof_offset(0),
        _glob_dof_count(0)
      {
      }

      /// Resets the whole object
      void clear()
      {
        _comm = nullptr;
        _glob_dof_offset = Index(0);
        _glob_dof_count = Index(0);
        _glob_dof_idx.clear();
        _owned_mirror.clear();
        _neighbors_owner.clear();
        _neighbors_donee.clear();
        _all_glob_dof_offset.clear();
        _all_glob_dof_counts.clear();
      }

      /// \returns The communicator
      const Dist::Comm* get_comm() const
      {
        return _comm;
      }

      /// \returns The number of global DOFs owned by this process.
      Index get_num_owned_dofs() const
      {
        return this->_owned_mirror.num_indices();
      }

      /// \returns The number of local DOFs shared by this process.
      Index get_num_local_dofs() const
      {
        return this->_owned_mirror.size();
      }

      /// \returns The total number of global DOFs.
      Index get_num_global_dofs() const
      {
        return this->_glob_dof_count;
      }

      /// \returns The index of the first global DOF owned by this process.
      Index get_global_dof_offset() const
      {
        return this->_glob_dof_offset;
      }

      /**
       * \brief Returns the global-dof-indices array.
       */
      const std::vector<Index> get_global_dof_indices() const
      {
        return this->_glob_dof_idx;
      }

      /**
       * \brief Maps a local DOF index to the corresponding global DOF index.
       *
       * \param[in] local_dof_idx
       * The local index of the DOF whose global index is to be returned.
       *
       * \returns
       * The global index of the local DOF.
       */
      Index map_local_to_global_index(const Index local_dof_idx) const
      {
        ASSERT(local_dof_idx < this->_owned_mirror.size());
        return this->_glob_dof_idx.at(local_dof_idx);
      }

      /**
       * \brief Maps an owned DOF index to the corresponding local DOF index.
       *
       * \param[in] owned_dof_idx
       * The owned index of the DOF whose local index is to be returned.
       *
       * \returns
       * The local index of the owned DOF.
       */
      Index map_owned_to_local_index(const Index owned_dof_idx) const
      {
        ASSERT(owned_dof_idx < this->_owned_mirror.num_indices());
        return this->_owned_mirror.indices()[owned_dof_idx];
      }

      /**
       * \brief Returns the Mirror of the owned DOFs.
       */
      const MirrorType& get_owned_mirror() const
      {
        return this->_owned_mirror;
      }

      /// \returns The number of owner neighbors.
      Index get_num_owner_neighbors() const
      {
        return Index(this->_neighbors_owner.size());
      }

      /// \returns The number of donee neighbors.
      Index get_num_donee_neighbors() const
      {
        return Index(this->_neighbors_donee.size());
      }

      /**
       * \brief Returns an owner neighbor rank.
       *
       * \param[in] i
       * The index of the owner neighbor.
       *
       * \returns
       * The rank of the i-th owner neighbor.
       */
      int get_owner_rank(Index i) const
      {
        return this->_neighbors_owner.at(i).first;
      }

      /**
       * \brief Returns an donee neighbor rank.
       *
       * \param[in] i
       * The index of the donee neighbor.
       *
       * \returns
       * The rank of the i-th donee neighbor.
       */
      int get_donee_rank(Index i) const
      {
        return this->_neighbors_donee.at(i).first;
      }

      /**
       * \brief Returns an owner neighbor mirror.
       *
       * \note
       * The indices of an owner neighbor mirror are owned DOF indices.
       *
       * \param[in] i
       * The index of the owner neighbor.
       *
       * \returns
       * The mirror of the i-th owner neighbor.
       */
      const MirrorType& get_owner_mirror(Index i) const
      {
        return this->_neighbors_owner.at(i).second;
      }

      /**
       * \brief Returns a donee neighbor mirror.
       *
       * \note
       * The indices of a donee neighbor mirror are local DOF indices.
       *
       * \param[in] i
       * The index of the donee neighbor.
       *
       * \returns
       * The mirror of the i-th donee neighbor.
       */
      const MirrorType& get_donee_mirror(Index i) const
      {
        return this->_neighbors_donee.at(i).second;
      }

      /**
       * \brief Returns the all-global-dof-offsets vector
       *
       * \attention
       * This vector is only allocated after #assemble_allgather() has been called.\n
       * Don't use this unless you really know what you are doing!
       *
       * \returns The all-global-dof-offsets vector.
       */
      const std::vector<int>& get_all_global_dof_offsets() const
      {
        XASSERTM(!this->_all_glob_dof_offset.empty(), "You did not ask to assemble the global dof offsets");
        return this->_all_glob_dof_offset;
      }

      /**
       * \brief Returns the all-global-dof-counts vector
       *
       * \attention
       * This vector is only allocated after #assemble_allgather() has been called.\n
       * Don't use this unless you really know what you are doing!
       *
       * \returns The all-global-dof-counts vector.
       */
      const std::vector<int>& get_all_global_dof_counts() const
      {
        XASSERTM(!this->_all_glob_dof_counts.empty(), "You did not ask to assemble the global dof counts");
        return this->_all_glob_dof_counts;
      }

      /**
       * \brief Assembles the required data for the AlgDofPartiVector::allgather() and
       * AlgDofPartiMatrix::apply() functions.
       *
       * \warning
       * Calling this function and (more importantly) the AlgDofPartiVector::allgather() function or
       * the AlgDofPartiMatrix::apply() functions will *totally* screw up the scalability of your
       * application, so the code *will* blow up in your face on bigger clusters.\n
       * Do not use this function except for <b>small-scale debugging</b> purposes. You have been warned.
       *
       * \param[in] yes_i_really_want_this
       * An assertion will fire if this is \c false. The only purpose of this parameter is
       * to force you to read the warning in this documentation first.
       */
      void assemble_allgather(bool yes_i_really_want_this = false)
      {
        XASSERTM(yes_i_really_want_this, "You probably don't want to do this!");

        std::size_t num_ranks = std::size_t(this->_comm->size());

        int my_nod = int(this->get_num_owned_dofs());
        int my_off = int(this->get_global_dof_offset());

        // resize vectors
        this->_all_glob_dof_offset.resize(num_ranks, 0);
        this->_all_glob_dof_counts.resize(num_ranks, 0);

        // allgather offsets + counts
        this->_comm->allgather(&my_off, 1, this->_all_glob_dof_offset.data(), 1);
        this->_comm->allgather(&my_nod, 1, this->_all_glob_dof_counts.data(), 1);
      }

      /**
       * \brief Assembles the AlgDofParti object from a given Global::Gate.
       *
       * This function performs the assembly of all required internal data structures.
       */
      void assemble_by_gate(const GateType& gate)
      {
        // set our communicator
        this->_comm = gate.get_comm();
        XASSERTM(this->_comm != nullptr, "need a communicator here!");

        // get the communicator
        const Dist::Comm& comm = *this->_comm;

        // get my rank and neighbor ranks
        const int my_rank = comm.rank();
        const std::vector<int>& neighbor_ranks = gate._ranks;

        // get number of local DOFs on our patch
        const Index num_local_dofs = gate._freqs.size();

        // As a very first step, we have to decide which process will become the
        // owner of each DOF. In this implementation, each DOF is owned by by the
        // process with the lowest rank, which simplifies a lot of things.

        // So, first create a vector, which stores the rank of the owner process
        // of each of our local DOFs and format the vector to our rank, i.e. at
        // startup assume that we own all local DOFs.
        std::vector<int> dof_owners(num_local_dofs, my_rank);

        // Now loop over all our neighbor processes
        for(std::size_t i(0); i < neighbor_ranks.size(); ++i)
        {
          // Check whether the neighbors rank is less than our rank, as otherwise
          // that particular neighbor can not own any of our local DOFs.
          const int neighbor_rank = neighbor_ranks.at(i);
          if(neighbor_rank < my_rank)
          {
            // loop over all DOFs in the mirror of that neighbor
            const MirrorType& mirror = gate._mirrors.at(i);
            const Index num_indices = mirror.num_indices();
            const IndexType* idx = mirror.indices();
            for(Index j(0); j < num_indices; ++j)
            {
              // get a reference to the current owner rank of that DOF and
              // update the DOF owner if that neighbor has a lower rank
              // than the previous "owner candidate"
              int& dof_owner = dof_owners.at(idx[j]);
              if(neighbor_rank < dof_owner)
                dof_owner = neighbor_rank;
            }
          }
        }

        // Okay, at this point we (and all other processes) know the owners of
        // each of our local DOFs. Now, we have to generate a vector of all
        // local DOF indices of the DOFs that we own:
        std::vector<Index> owned_dofs;
        owned_dofs.reserve(num_local_dofs);

        // Furthermore, we need the 'inverse' information, i.e. for each of our
        // local DOFs, we have to store the "owned DOF index" if we own this
        // particular DOF. We initialize the vector with ~0, which stands for
        // "not my DOF":
        std::vector<Index> own_dof_idx(num_local_dofs, ~Index(0));

        // Loop over all local DOFs
        for(Index i(0); i < num_local_dofs; ++i)
        {
          // Am I the owner of this DOF?
          if(dof_owners[i] == my_rank)
          {
            // Yes, that's my DOF, so let's give it a new 'my owned DOF' index:
            own_dof_idx.at(i) = Index(owned_dofs.size());
            // And remember that I own this DOF:
            owned_dofs.push_back(i);
          }
        }

        /// \todo what do we do if owned_dofs==0 ?
        XASSERTM(!owned_dofs.empty(), "this process has no DOFs!");

        // get number of my owned DOFS
        const Index num_owned_dofs = Index(owned_dofs.size());

        // Next, we have to determine the global offset of our first owned DOF, which
        // is easily obtained by an exclusive-scan over each processes owned DOF count:
        comm.exscan(&num_owned_dofs, &this->_glob_dof_offset, 1, Dist::op_sum);

        // Furthermore, we also perform an allreduce to obtain the total number of
        // global DOFs, which is generally helpful and can be used for sanity checks.
        comm.allreduce(&num_owned_dofs, &this->_glob_dof_count, 1, Dist::op_sum);

        // Next, we have to set up the mirror for all of our owned DOFs, which we can
        // easily generate from our "owned_dofs" vector:
        this->_owned_mirror = _vidx2mirror(num_local_dofs, owned_dofs);

        // Now we also have to create the mirrors for each of our neighbor processes:
        for(std::size_t i(0); i < neighbor_ranks.size(); ++i)
        {
          // get the neighbor rank and its mirror
          const int neighbor_rank = neighbor_ranks.at(i);
          const MirrorType& mirror = gate._mirrors.at(i);

          // There are two cases now:
          // 1) The neighbor process has a lower rank, so it is a (potential) owner
          //    of at least one of our local DOFs.
          // 2) The neighbor process has a greater rank, so we are an owner of at least
          //    one of the neighbor process's local DOFs.
          if(neighbor_rank < my_rank)
          {
            // This is a potential "owner-neighbor", which may own one or several of
            // our local DOFs. We call a helper function, which will create a vector
            // of all *local* DOF indices, which are owned by that neighbor process.
            // Note that it is perfectly legit if the vector is empty, as this simply
            // means that all our local DOFs, which are shared with that particular
            // neighbor, are owned by some *other* neighbor process(es) and not by
            // the particular neighbor that we are currently considering.
            std::vector<Index> v = _owner_vidx(mirror, dof_owners, neighbor_rank);
            if(!v.empty())
            {
              // That neighbor owns at least one of our local DOFs, so create a mirror
              // from the index-vector and push it into our list of "owner-neighbors".
              // Note that the indices in the mirror are *local* DOF indices.
              _neighbors_owner.emplace_back(std::make_pair(neighbor_rank, _vidx2mirror(num_local_dofs, v)));
            }
          }
          else
          {
            // This is a potential "donee-neighbor", which shares one or several of
            // our local DOFs, which we might own. We call a helper function, which
            // will create a vector of all our *owned* DOF indices, which are shared
            // with that particular neighbor process.
            // Note that it is perfectly legit if the vector is empty, as this simply
            // means that all our local DOFs, which are shared with that particular
            // neighbor, are owned by some *other* neighbor process(es) and not by
            // this process (i.e. us).
            std::vector<Index> v = _donee_vidx(mirror, own_dof_idx);
            if(!v.empty())
            {
              // We own at least of one of the neighbor's local DOFs, so create a mirror
              // from the index-vector and push it into our list of "donee-neighbors".
              // Note that the indices in the mirror are *owned* DOF indices.
              _neighbors_donee.emplace_back(std::make_pair(neighbor_rank, _vidx2mirror(num_owned_dofs, v)));
            }
          }
        }

        // Finally, we have to determine the global DOF index for each of our local DOFs,
        // so that we can perform a "local-to-global" DOF index lookup. This information
        // is required for the assembly of matrices later on, so we have to store this in
        // a member-variable vector.

        // Let's create the vector of the required length, initialize all its indices to
        // ~0, which will can be used for a sanity check later on to ensure that we know
        // the global indices of all our local DOFs:
        this->_glob_dof_idx.resize(num_local_dofs, ~Index(0));

        // Next, let's loop over all DOFs that we own and assign their corresponding
        // global DOF indices, starting with our global DOF offset, which we have already
        // determined before:
        for(std::size_t i(0); i < owned_dofs.size(); ++i)
          this->_glob_dof_idx[owned_dofs[i]] = this->_glob_dof_offset + Index(i);

        // Finally, we also need to query the global DOF indices of all our local DOFs,
        // which are owned by some neighbor process (i.e. we are the donee for these DOFs).
        // Of course, we also have to send the global DOF indices of all DOFs that we own
        // to all of our "donee-neighbors".

        // Allocate index buffers and post receives for all of our owner-neighbors:
        const Index num_neigh_owner = Index(_neighbors_owner.size());
        Dist::RequestVector recv_reqs(num_neigh_owner);
        std::vector<std::vector<Index>> recv_bufs(num_neigh_owner);
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          const std::size_t num_idx = _neighbors_owner.at(i).second.num_indices();
          recv_bufs.at(i).resize(num_idx);
          recv_reqs[i] = comm.irecv(recv_bufs.at(i).data(), num_idx, _neighbors_owner.at(i).first);
        }

        // Allocate index buffers, fill them with the global DOF indices of our owned
        // DOFs and send them to the donee-neighbors:
        const Index num_neigh_donee = Index(_neighbors_donee.size());
        Dist::RequestVector send_reqs(num_neigh_donee);
        std::vector<std::vector<Index>> send_bufs(num_neigh_donee);
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          // Get the mirror for this donee-neighbor:
          const MirrorType& mirror = _neighbors_donee.at(i).second;
          const Index num_idx = mirror.num_indices();
          const Index* mir_idx = mirror.indices();
          // Allocate buffer of required size and translate the indices of our mirror,
          // which are "owned DOF indices", to global DOF indices, where:
          //   global_dof_index := global_dof_offset + owned_dof_index
          send_bufs.at(i).resize(num_idx);
          Index* jdx = send_bufs.at(i).data();
          for(Index j(0); j < num_idx; ++j)
            jdx[j] = this->_glob_dof_offset + mir_idx[j];
          send_reqs[i] = comm.isend(send_bufs.at(i).data(), num_idx, _neighbors_donee.at(i).first);
        }

        // Now process all receive requests from our owner-neighbors:
        for(std::size_t i(0u); recv_reqs.wait_any(i); )
        {
          const MirrorType& mirror = _neighbors_owner.at(i).second;
          const Index num_idx = mirror.num_indices();
          const Index* mir_idx = mirror.indices();
          const Index* jdx = recv_bufs.at(i).data();
          // Scatter the global DOF indices, which our friendly owner-neighbor has
          // provided us with, into our global DOF index vector:
          for(Index j(0); j < num_idx; ++j)
            this->_glob_dof_idx[mir_idx[j]] = jdx[j];
        }

        // wait for all sends to finish
        send_reqs.wait_all();

#ifdef DEBUG
        // Last but not least: debug-mode sanity check
        // We now should know the global DOF indices for all of our local DOFs:
        for(Index gdi : this->_glob_dof_idx)
        {
          XASSERTM(gdi != ~Index(0), "invalid global DOF index");
        }
#endif // DEBUG
      }

    private:
      // auxiliary function: create a mirror from a vector of indices
      static MirrorType _vidx2mirror(const Index size, const std::vector<Index>& vidx)
      {
        const Index nidx = Index(vidx.size());
        MirrorType mirror(size, nidx);
        IndexType* idx = mirror.indices();
        for(Index i(0); i < nidx; ++i)
          idx[i] = vidx[i];
        return mirror;
      }

      // auxiliary function: create an index vector of all *local DOFs*
      // which are owned by the process with the given neighbor rank:
      static std::vector<Index> _owner_vidx(const MirrorType& old_mir,
        const std::vector<int>& dof_owners, const int neighbor_rank)
      {
        std::vector<Index> v;
        v.reserve(old_mir.num_indices());

        const IndexType* old_idx = old_mir.indices();
        for(Index i(0); i < old_mir.num_indices(); ++i)
        {
          // is this neighbor the owner?
          if(neighbor_rank == dof_owners[old_idx[i]])
            v.push_back(old_idx[i]);
        }
        return v;
      }

      // auxiliary function: create an index vector all *owned DOFs*
      // which are local DOFs of the process whose mirror is given:
      static std::vector<Index> _donee_vidx(const MirrorType& old_mir,
        const std::vector<Index>& own_dof_idx)
      {
        std::vector<Index> v;
        v.reserve(old_mir.num_indices());

        const IndexType* old_idx = old_mir.indices();
        for(Index i(0); i < old_mir.num_indices(); ++i)
        {
          // do we own the DOF in this mirror?
          if(own_dof_idx[old_idx[i]] != ~Index(0))
            v.push_back(own_dof_idx[old_idx[i]]);
        }
        return v;
      }
    }; // class AlgDofParti

    /**
     * \brief Algebraic DOF Partitioned Vector class template
     *
     * This class implements the management of a distributed vector based on an algebraic
     * DOF partitioning. A vector of this class will contain only the owned DOFs of this
     * process. The two basic operations provided by this class are the \b upload and the
     * \b download of vectors, where:
     * - \b upload: copy the DOFs of a Global::Vector into this AlgDofPartiVector
     * - \b download: copy the DOFs of this AlgDofPartiVector into a Global::Vector
     *
     * Note that objects of this class do not require any additional initialization as
     * the constructor, which expects and AlgDofPartiObject as its parameters, will
     * automatically allocate the internal vector.
     *
     * This class also offers a function named #allgather(), which creates a vector
     * containing \e all global DOFs on each process. This function is only meant to
     * be used for small-scale debugging purposes and shall not be used for actual
     * production work, because it a) will totally screw the scalability of your code
     * on a large number of processes and b) may easily lead to fatal out-of-memory
     * events. You have been warned.
     *
     * \author Peter Zajac
     */
    template<typename LocalVector_, typename Mirror_>
    class AlgDofPartiVector
    {
    public:
      /// our memory tag type
      typedef typename LocalVector_::MemType MemType;
      /// our data type
      typedef typename LocalVector_::DataType DataType;
      /// our index type
      typedef typename LocalVector_::IndexType IndexType;

      /// the local vector type
      typedef LocalVector_ LocalVectorType;
      /// the vector mirror type
      typedef Mirror_ MirrorType;

      /// the global vector type
      typedef Global::Vector<LocalVector_, Mirror_> GlobalVectorType;

      /// the algebraic DOF partitioning type
      typedef Global::AlgDofParti<LocalVector_, Mirror_> AlgDofPartiType;

      /// the buffer vector type used for communication
      typedef LAFEM::DenseVector<MemType, DataType, IndexType> BufferVectorType;
      /// the vector type used for the internal storage of our owned DOFs
      typedef LAFEM::DenseVector<MemType, DataType, IndexType> OwnedVectorType;

    protected:
      /// the algebraic dof-partitioning
      const AlgDofPartiType* _alg_dof_parti;
      /// the internal vector object storing our owned dofs
      OwnedVectorType _vector;

    public:
      /// default constructor
      AlgDofPartiVector() :
        _alg_dof_parti(nullptr),
        _vector()
      {
      }

      /**
       * \brief Creates a ADP vector based on an algebraic dof partitioning
       *
       * Once this constructor returns, this vector does not require any further
       * initialization and is therefore ready-to-use.
       *
       * \param[in] adp_in
       * The algebraic dof partitioning to be used for this vector.
       */
      explicit AlgDofPartiVector(const AlgDofPartiType* adp_in) :
        _alg_dof_parti(adp_in),
        _vector(adp_in->get_num_owned_dofs())
      {
      }

      /// Resets the whole object
      void clear()
      {
        this->_alg_dof_parti = nullptr;
        this->_vector.clear();
      }

      /// \returns A reference to the internal owned vector object.
      OwnedVectorType& owned()
      {
        return this->_vector;
      }

      /// \returns A reference to the internal owned vector object.
      const OwnedVectorType& owned() const
      {
        return this->_vector;
      }

      /// \return A pointer to the underlying algebraic dof partitioning object.
      const AlgDofPartiType* get_alg_dof_parti() const
      {
        return this->_alg_dof_parti;
      }

      /// \return A pointer to the underlying communicator
      const Dist::Comm* get_comm() const
      {
        XASSERT(this->_alg_dof_parti != nullptr);
        return this->_alg_dof_parti->get_comm();
      }

      /**
       * \brief Gathers the full global vector on all processes
       *
       * \attention
       * This function makes use of the infamous allgatherv MPI operation, which is known
       * to completely blow the scalability of an application running on a large number of
       * processes. In addition to the scalability issues, please keep in mind that the
       * global vector may be very large -- even potentially larger that the total amount
       * of memory that this available for a single process, which would inevitably result
       * in an application crash due to an out-of-memory event.\n
       * Do not use this function except for <b>small-scale debugging</b> purposes. You have been warned.
       *
       * \param[inout] vec_full
       * A vector that receives all global dofs. Its length is assumed to be initialized to
       * the global number of degrees of freedom.
       */
      void allgather(BufferVectorType& vec_full) const
      {
        XASSERTM(!this->_alg_dof_parti->get_all_global_dof_counts().empty(),
          "You did not ask to assemble the required allgather data");

        XASSERT(vec_full.size() == this->_alg_dof_parti->get_num_global_dofs());
        this->_alg_dof_parti->get_comm()->allgatherv(this->_vector.elements(),
          this->_alg_dof_parti->get_num_owned_dofs(),
          vec_full.elements(),
          this->_alg_dof_parti->get_all_global_dof_counts().data(),
          this->_alg_dof_parti->get_all_global_dof_offsets().data());
      }

      /**
       * \brief Copies the DOF values from a global vector to this algebraic-dof-partitioned vector.
       *
       * \param[in] global_vector
       * The global vector whose values are to be copied into this vector.
       */
      void upload(const GlobalVectorType& global_vector)
      {
        this->upload(global_vector.local());
      }

      /**
       * \brief Copies the DOF values from a local vector to this algebraic-dof-partitioned vector.
       *
       * \param[in] local_vector
       * The local vector whose values are to be copied into this vector.
       */
      void upload(const LocalVectorType& local_vector)
      {
        XASSERT(this->_alg_dof_parti != nullptr);
        XASSERT(this->_vector.size() == this->_alg_dof_parti->get_num_owned_dofs());
        XASSERT(local_vector.size() == this->_alg_dof_parti->get_num_local_dofs());

        // Uploading is simple: we have a mirror that will extract all our owned DOFs
        // from our local vector, so we just have to call this mirror's gather function:
        this->_alg_dof_parti->get_owned_mirror().gather(this->_vector, local_vector);
      }

      /**
       * \brief Copies the DOF values from this algebraic-dof-partitioned vector into a global vector.
       *
       * \param[inout] global_vector
       * The global vector into which our DOF values are to be copied. The vector is assumed to be
       * allocated to the correct length, but its DOF values upon call are ignored.
       */
      void download(GlobalVectorType& global_vector) const
      {
        this->download(global_vector.local());
      }

      /**
       * \brief Copies the DOF values from this algebraic-dof-partitioned vector into a local vector.
       *
       * \param[inout] local_vector
       * The local vector into which our DOF values are to be copied. The vector is assumed to be
       * allocated to the correct length, but its DOF values upon call are ignored.
       */
      void download(LocalVectorType& local_vector) const
      {
        XASSERT(this->_alg_dof_parti != nullptr);
        XASSERT(this->_vector.size() == this->_alg_dof_parti->get_num_owned_dofs());
        XASSERT(local_vector.size() == this->_alg_dof_parti->get_num_local_dofs());

        // get our algebraic dof partitioning object
        const AlgDofPartiType& adp = *this->_alg_dof_parti;

        // get the communicator
        const Dist::Comm& comm = *this->get_comm();

        // get the number of owner- and donee-neighbors
        const Index num_neigh_owner = adp.get_num_owner_neighbors();
        const Index num_neigh_donee = adp.get_num_donee_neighbors();

        // The basic idea is simple:
        // 1) Send the shared DOFs to all of our donee-neighbors
        // 2) Receive the shared DOFs from all of our owner-neighbors
        // 3) Scatter our own owned DOFs into our local vector

        // Create receive buffers and post receives for all owner-neighbors
        Dist::RequestVector recv_reqs(num_neigh_owner);
        std::vector<BufferVectorType> recv_bufs(num_neigh_owner);
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          // create a vector buffer
          // note: owner mirrors relate to local DOF indices
          recv_bufs.at(i) = adp.get_owner_mirror(i).create_buffer(local_vector);
          // post receive from owner neighbor
          recv_reqs[i] = comm.irecv(recv_bufs.at(i).elements(), recv_bufs.at(i).size(), adp.get_owner_rank(i));
        }

        // Create send buffers, fill them with our owned DOFs and send this
        // to our donee-neighbors
        Dist::RequestVector send_reqs(num_neigh_donee);
        std::vector<BufferVectorType> send_bufs(num_neigh_donee);
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          // get mirror, create buffer and gather the shared DOFs
          // note: donee mirrors relate to owned DOF indices
          const MirrorType& mir = adp.get_donee_mirror(i);
          send_bufs.at(i) = mir.create_buffer(this->_vector);
          mir.gather(send_bufs.at(i), this->_vector);
          // post send to donee neighbor
          send_reqs[i] = comm.isend(send_bufs.at(i).elements(), send_bufs.at(i).size(), adp.get_donee_rank(i));
        }

        // format the output vector
        local_vector.format();

        // scatter our own owned DOFs into our output vector
        adp.get_owned_mirror().scatter_axpy(local_vector, this->_vector);

        // process receives from our owner-neighbors
        for(std::size_t idx(0u); recv_reqs.wait_any(idx); )
        {
          // scatter received DOFs into our output vector
          adp.get_owner_mirror(Index(idx)).scatter_axpy(local_vector, recv_bufs.at(idx));
        }

        // at this point, all receives should have finished
        XASSERT(recv_reqs.is_null());

        // wait for all sends to finish
        send_reqs.wait_all();
      }
    }; // class AlgDofPartiVector

    /**
     * \brief Algebraic DOF Partitioned CSR matrix class template
     *
     * This class implements the management of a distributed CSR matrix based on an algebraic
     * DOF partitioning. A matrix of this class will contain only the (full) rows corresponding
     * to the owned DOFs of this process, where the column indices are global DOF indices.
     *
     * In contrast to the AlgDofPartiVector, objects of this class need further initialization
     * after the call of the constructor. This additional initialization is performed by
     * uploading at least the layout of a given Global::Matrix. In analogy to the concept of
     * our (linear) solvers, the "upload" process is split into two functions:
     * - upload_symbolic(): creates all internal data structures and computes the internal
     *   matrix structure based on the given input matrix structure.
     * - upload_numeric(): copies the numerical values of the input matrix into the owned
     *   matrix and performs the required communication.
     *
     * Also in analogy to the solver concept, one may repeatedly call upload_numeric() with
     * changing matrix values without the need to call upload_symbolic() every time, as long
     * as the input matrix structure does not change. This can seriously reduce communication
     * overhead in non-linear scenarios, where the matrix values change frequently.
     *
     * Note that, in contrast to the AlgDofPartiVector class, this matrix class does not
     * offer a "download" functionality.
     *
     * \author Peter Zajac
     */
    template<typename LocalMatrix_, typename Mirror_>
    class AlgDofPartiMatrix
    {
    public:
      /// our memory tag type
      typedef typename LocalMatrix_::MemType MemType;
      /// our data type
      typedef typename LocalMatrix_::DataType DataType;
      /// our index type
      typedef typename LocalMatrix_::IndexType IndexType;

      /// the local matrix type
      typedef LocalMatrix_ LocalMatrixType;
      /// the vector mirror type
      typedef Mirror_ MirrorType;
      /// the global matrix type
      typedef Global::Matrix<LocalMatrixType, MirrorType, MirrorType> GlobalMatrixType;

      /// the local vector type
      typedef typename LocalMatrixType::VectorTypeR LocalVectorType;

      /// the buffer vector type used for communication
      typedef LAFEM::DenseVector<Mem::Main, DataType, IndexType> BufferVectorType;
      /// the buffer matrix type used for communication
      typedef LAFEM::MatrixMirrorBuffer<Mem::Main, DataType, IndexType> BufferMatrixType;

      /// the matrix type used for the internal storage of our owned matrix rows
      typedef LAFEM::SparseMatrixCSR<Mem::Main, DataType, IndexType> OwnedMatrixType;

      /// the algebraic DOF partitioning
      typedef Global::AlgDofParti<LocalVectorType, MirrorType> AlgDofPartiType;
      /// the compatible alg-dof-parti vector type
      typedef Global::AlgDofPartiVector<LocalVectorType, MirrorType> AlgDofPartiVectorType;

    protected:
      /// the algebraic dof-partitioning
      const AlgDofPartiType* _alg_dof_parti;
      /// the internal matrix storing our owned matrix rows
      OwnedMatrixType _matrix;

      /// the matrix buffers for our donee and owner neighbors
      std::vector<BufferMatrixType> _donee_bufs, _owner_bufs;
      /// the data array mirrors for our donee and owner neighbors
      std::vector<MirrorType> _data_donee_mirs, _data_owner_mirs;
      /// the data array "mirror" for this process
      std::vector<std::pair<Index,Index>> _my_data_mir;

    public:
      /// default constructor
      AlgDofPartiMatrix() :
        _alg_dof_parti(nullptr),
        _matrix()
      {
      }

      /**
       * \brief Creates a ADP matrix based on an algebraic dof partitioning
       *
       * \attention
       * In contrast to the AlgDofPartiVector class, this matrix class is \b not
       * fully initialized after this constructor returns. It is at least required
       * to set up the internal data structures by calling the #upload_symbolic()
       * or #upload() function and passing an initialized matrix to it.
       *
       * \param[in] adp_in
       * The algebraic dof partitioning to be used for this matrix.
       */
      explicit AlgDofPartiMatrix(const AlgDofPartiType* adp_in) :
        _alg_dof_parti(adp_in),
        _matrix()
      {
      }

      /// Resets the whole object
      void clear()
      {
        _alg_dof_parti = nullptr;
        _matrix.clear();
        _donee_bufs.clear();
        _owner_bufs.clear();
        _data_donee_mirs.clear();
        _data_owner_mirs.clear();
        _my_data_mir.clear();
      }

      /// \returns A reference to the internal owned matrix object.
      OwnedMatrixType& owned()
      {
        return this->_matrix;
      }

      /// \returns A reference to the internal owned matrix object.
      const OwnedMatrixType& owned() const
      {
        return this->_matrix;
      }

      /// \returns A pointer to the underlying algebraic dof partitioning object.
      const AlgDofPartiType* get_alg_dof_parti() const
      {
        return this->_alg_dof_parti;
      }

      /// \returns A pointer to the underlying communicator.
      const Dist::Comm* get_comm() const
      {
        XASSERT(this->_alg_dof_parti != nullptr);
        return this->_alg_dof_parti->get_comm();
      }

      /**
       * \brief Performs a global matrix-vector product with this matrix.
       *
       * \attention
       * As convenient and tempting as this function may look, please be aware that
       * using this function will completely destroy the scalability of our code when
       * running on a large number of processes. In addition to the scalability issues,
       * please be aware that this function may cause a fatal out-of-memory event, as
       * the underlying (naive) implementation gathers the full multiplicand vector
       * on all processes using the AlgDofPartiVector::allgather() function.\n
       * Do not use this function except for <b>small-scale debugging</b> purposes. You have been warned.
       *
       * \param[out] vec_r
       * A reference to the algebraic-dof-partitioned vector that receives
       * the result of the matrix-vector product.
       *
       * \param[in] vec_x
       * A const reference to the algebraic-dof-partitioned vector that
       * acts as a multiplicand vector.
       */
      void apply(AlgDofPartiVectorType& vec_r, const AlgDofPartiVectorType& vec_x) const
      {
        XASSERTM(!this->_alg_dof_parti->get_all_global_dof_counts().empty(),
          "You did not ask to assemble the required allgather data");

        // create full vector
        LocalVectorType vec_full(this->_alg_dof_parti->get_num_global_dofs());

        // gather full vector
        vec_x.allgather(vec_full);

        // apply our matrix
        this->_matrix.apply(vec_r.owned(), vec_full);
      }

      /**
       * \brief Initializes this algebraic dof-partitioned matrix from a global matrix.
       *
       * \param[in] mat_a
       * The global matrix whose layout and values are to be used to initialize this matrix.
       */
      void upload(const GlobalMatrixType& mat_a)
      {
        upload(mat_a.local());
      }

      /**
       * \brief Initializes this algebraic dof-partitioned matrix from a global matrix.
       *
       * \param[in] mat_a
       * The global matrix whose layout and values are to be used to initialize this matrix.
       */
      void upload(const LocalMatrixType& mat_a)
      {
        upload_symbolic(mat_a);
        upload_numeric(mat_a);
      }

      /**
       * \brief Initializes this algebraic dof-partitioned matrix from a global matrix.
       *
       * This function performs the "symbolic" initialization of this algebraic-dof-partitioned
       * matrix based on the layout of the given input matrix. After the call of ths function,
       * this matrix object is fully initialized, although the actual numerical values of
       * the matrix entries are still undefined and have to be copied from the global input
       * matrix lateron by using the #upload_numeric() function.
       *
       * \param[in] mat_a
       * The global matrix whose layout is to be used to initialize this matrix.
       */
      void upload_symbolic(const GlobalMatrixType& mat_a)
      {
        upload_symbolic(mat_a.local());
      }

      /**
       * \brief Initializes this algebraic-dof-partitioned matrix from a local matrix.
       *
       * This function performs the "symbolic" initialization of this algebraic-dof-partitioned
       * matrix based on the layout of the given input matrix. After the call of this function,
       * this matrix object is fully initialized, although the actual numerical values of
       * the matrix entries are still undefined and have to be copied from the local input
       * matrix later on by using the #upload_numeric() function.
       *
       * \param[in] mat_a
       * The local matrix whose layout is to be used to initialize this matrix.
       */
      void upload_symbolic(const LocalMatrixType& mat_a)
      {
        XASSERT(this->_alg_dof_parti != nullptr);

        this->_assemble_buffers(mat_a);
        this->_assemble_structure(mat_a);
        this->_assemble_data_mirrors(mat_a);
      }

      /**
       * \brief Copies the matrix entry values from the input matrix into this algebraic-dof-partitioned matrix.
       *
       * \param[in] mat_a
       * The input matrix whose entry values are to be copied.
       *
       * \attention
       * The input matrix is assumed to be an unfiltered type-0 matrix.
       */
      void upload_numeric(const GlobalMatrixType& mat_a)
      {
        upload_numeric(mat_a.local());
      }

      /**
       * \brief Copies the matrix entry values from the input matrix into this algebraic-dof-partitioned matrix.
       *
       * \param[in] mat_a
       * The input matrix whose entry values are to be copied.
       *
       * \attention
       * The input matrix is assumed to be an unfiltered type-0 matrix.
       */
      void upload_numeric(const LocalMatrixType& mat_a)
      {
        XASSERT(this->_alg_dof_parti != nullptr);

        // get our partitioning
        const AlgDofPartiType& adp = *this->_alg_dof_parti;

        // get our communicator
        const Dist::Comm& comm = *this->get_comm();

        // format our internal matrix
        this->_matrix.format();

        const Index num_neigh_owner = adp.get_num_owner_neighbors();
        const Index num_neigh_donee = adp.get_num_donee_neighbors();

        Dist::RequestVector recv_reqs(num_neigh_donee);
        Dist::RequestVector send_reqs(num_neigh_owner);

        std::vector<BufferVectorType> donee_vbufs(num_neigh_donee);
        std::vector<BufferVectorType> owner_vbufs(num_neigh_owner);

        // create receive buffers and post receives
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          // create buffer
          BufferVectorType& buf = donee_vbufs.at(i);
          buf = BufferVectorType(this->_data_donee_mirs.at(i).num_indices());

          // post receive
          recv_reqs[i] = comm.irecv(buf.elements(), buf.size(), adp.get_donee_rank(i));
        }

        // post sends
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          // create buffer
          BufferVectorType& buf = owner_vbufs.at(i);
          buf = BufferVectorType(this->_data_owner_mirs.at(i).num_indices());

          // gather from mirror
          this->_gather_data(buf, mat_a, this->_data_owner_mirs.at(i));

          // post send
          send_reqs[i] = comm.isend(buf.elements(), buf.size(), adp.get_owner_rank(i));
        }

        // upload our own data
        this->_upload_own(mat_a);

        // process all pending receives
        for(std::size_t idx(0u); recv_reqs.wait_any(idx); )
        {
          // scatter from buffer
          this->_scatter_data(this->_matrix, donee_vbufs.at(idx), this->_data_donee_mirs.at(idx));
        }

        // wait for all sends to finish
        send_reqs.wait_all();
      }

      /**
       * \brief Applies a global unit-filter into this matrix.
       *
       * \attention
       * You must use this function to filter the internal partitioned matrix.
       * You cannot use the filter_mat() function of the unit-filter here, because
       * the indices of the unit-filter relate to local DOF indices, whereas our
       * internal matrix object works in owned DOF indices.
       *
       * \param[in] filter
       * The unit-filter that is to be applied onto the matrix. May be empty.
       */
      template<typename LocalFilter_>
      void filter_matrix(const Global::Filter<LocalFilter_, MirrorType>& filter)
      {
        filter_matrix(filter.local());
      }

      /**
       * \brief Applies a unit-filter into this matrix.
       *
       * \attention
       * You must use this function to filter the internal partitioned matrix.
       * You cannot use the filter_mat() function of the unit-filter here, because
       * the indices of the unit-filter relate to local DOF indices, whereas our
       * internal matrix object works in owned DOF indices.
       *
       * \param[in] filter
       * The unit-filter that is to be applied onto the matrix. May be empty.
       */
      void filter_matrix(const LAFEM::UnitFilter<MemType, DataType, IndexType>& filter)
      {
        // empty filter?
        if(filter.used_elements() <= Index(0))
          return;

        // get our partitioning
        const AlgDofPartiType& adp = *this->_alg_dof_parti;

        // get global owned DOF offset and count
        const Index off_glob = adp.get_global_dof_offset();
        const Index num_glob = adp.get_num_global_dofs();

        // get local filter indices
        const Index num_idx = filter.used_elements();
        const IndexType* fil_idx = filter.get_indices();

        // get owned matrix arrays
        const IndexType* row_ptr = this->_matrix.row_ptr();
        const IndexType* col_idx = this->_matrix.col_ind();
        DataType* val = this->_matrix.val();

        // loop over all filter indices
        for(Index i(0); i < num_idx; ++i)
        {
          // translate local to global filter DOF index
          const Index gidx = adp.map_local_to_global_index(fil_idx[i]);

          // determine whether we own this DOF
          if((gidx < off_glob) || (off_glob + num_glob <= gidx))
            continue; // not my DOF

          // okay, that's my DOF, so filter the row
          const Index row = gidx - off_glob;
          for(Index j(row_ptr[row]); j < row_ptr[row+1]; ++j)
          {
            val[j] = DataType(col_idx[j] == gidx ? 1 : 0);
          }
        }
      }

    protected:
      /**
       * \brief Auxiliary function: creates a matrix buffer for a given owner mirror
       *
       * This function creates a matrix buffer for a given owner mirror, which contains
       * all matrix rows that correspond to the owned DOF indices stored in the mirror.
       * The column indices of the buffer correspond to global DOF indices.
       *
       * \param[in] local_matrix
       * The local matrix that acts as a layout template.\n
       * Its layout must be initialized, but its numerical values are ignored.
       *
       * \param[in] mirror
       * The mirror that contains the local DOF indices, which correspond to the
       * row indices of the local matrix which are to be extracted into the buffer.
       *
       * \param[in] gdi
       * An array containing the global DOF indices for all local DOFs.
       *
       * \param[in] num_glob_dofs
       * The total number of global DOFs over all processes.
       */
      static BufferMatrixType _asm_mat_buf(
        const LocalMatrixType& local_matrix, const MirrorType& mirror,
        const std::vector<Index>& gdi, const Index num_glob_dofs)
      {
        const Index num_idx = mirror.num_indices();
        const IndexType* mir_idx = mirror.indices();

        const IndexType* row_ptr = local_matrix.row_ptr();
        const IndexType* col_idx = local_matrix.col_ind();

        // count number of non-zeros in matrix buffer
        Index nnze = Index(0);
        for(Index i(0); i < num_idx; ++i)
        {
          nnze += (row_ptr[mir_idx[i]+1] - row_ptr[mir_idx[i]]);
        }

        // allocate matrix buffer of appropriate dimensions
        BufferMatrixType buf(num_idx, num_glob_dofs, nnze, 1);
        IndexType* buf_ptr = buf.row_ptr();
        IndexType* buf_idx = buf.col_ind();

        // set up buffer structure
        buf_ptr[0] = Index(0);
        for(Index i(0); i < num_idx; ++i)
        {
          const Index row = mir_idx[i];
          Index k = buf_ptr[i];
          // map local DOFs to global DOFs
          for(Index j(row_ptr[row]); j < row_ptr[row+1]; ++j, ++k)
            buf_idx[k] = gdi[col_idx[j]];
          buf_ptr[i+1] = k;
        }

        // Note: we leave the buffer indices unsorted

        // okay
        return buf;
      }

      /**
       * \brief Assembles the owner and donee matrix buffers.
       *
       * \param[in] local_matrix
       * The local matrix that acts as a layout template.\n
       * Its layout must be initialized, but its numerical values are ignored.
       */
      void _assemble_buffers(const LocalMatrixType& local_matrix)
      {
        // get our partitioning
        const AlgDofPartiType& adp = *this->_alg_dof_parti;

        // get out communicator
        const Dist::Comm& comm = *this->get_comm();

        // get number of neighbors and allocate buffer vectors
        const Index num_neigh_owner = adp.get_num_owner_neighbors();
        const Index num_neigh_donee = adp.get_num_donee_neighbors();
        _owner_bufs.resize(num_neigh_owner);
        _donee_bufs.resize(num_neigh_donee);

        // create the owner buffers by using the auxiliary helper function.
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          this->_owner_bufs.at(i) = _asm_mat_buf(local_matrix, adp.get_owner_mirror(i),
            adp.get_global_dof_indices(), adp.get_num_global_dofs());
        }

        // We have assembled our owner-neighbor matrix buffers, however, we cannot
        // assemble the donee-neighbor matrix buffers on our own. Instead, each owner
        // neighbor has to receive its donee-buffers from its donee-neighbors, because
        // only the donee knows the matrix layout of those buffers.
        // Yes, this is brain-twisting, but believe me, it cannot be realized otherwise.

        // Note:
        // The following code is effectively just copy-&-paste from the SynchMatrix::init()
        // function, as the same buffer problem also arises when converting global type-0
        // matrices to local type-1 matrices.

        Dist::RequestVector send_reqs(num_neigh_owner);
        Dist::RequestVector recv_reqs(num_neigh_donee);

        // receive buffer dimensions vector
        std::vector<std::array<Index,4>> send_dims(num_neigh_owner);
        std::vector<std::array<Index,4>> recv_dims(num_neigh_donee);

        // post send-buffer dimension receives
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          recv_reqs[i] = comm.irecv(recv_dims.at(i).data(), std::size_t(4), adp.get_donee_rank(i));
        }

        // send owner-buffer dimensions
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          const BufferMatrixType& sbuf = this->_owner_bufs.at(i);
          send_dims.at(i)[0] = sbuf.rows();
          send_dims.at(i)[1] = sbuf.columns();
          send_dims.at(i)[2] = sbuf.entries_per_nonzero();
          send_dims.at(i)[3] = sbuf.used_elements();
          send_reqs[i] = comm.isend(send_dims.at(i).data(), std::size_t(4), adp.get_owner_rank(i));
        }

        // wait for all receives to finish
        recv_reqs.wait_all();

        // allocate donee-buffers
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          // get the receive buffer dimensions
          Index nrows = recv_dims.at(i)[0];
          Index ncols = recv_dims.at(i)[1];
          Index nepnz = recv_dims.at(i)[2];
          Index nnze  = recv_dims.at(i)[3];

          // allocate receive buffer
          this->_donee_bufs.at(i) = BufferMatrixType(nrows, ncols, nnze, nepnz);
        }

        // post donee-buffer row-pointer array receives
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          recv_reqs[i] = comm.irecv(this->_donee_bufs.at(i).row_ptr(),
            this->_donee_bufs.at(i).rows()+std::size_t(1), adp.get_donee_rank(i));
        }

        // wait for all previous sends to finish
        send_reqs.wait_all();

        // post owner-buffer row-pointer array sends
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          send_reqs[i] = comm.isend(this->_owner_bufs.at(i).row_ptr(),
            this->_owner_bufs.at(i).rows()+std::size_t(1), adp.get_owner_rank(i));
        }

        // wait for all previous receives to finish
        recv_reqs.wait_all();

        // post donee-buffer column-index array receives
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          recv_reqs[i] = comm.irecv(this->_donee_bufs.at(i).col_ind(),
            this->_donee_bufs.at(i).used_elements(), adp.get_donee_rank(i));
        }

        // wait for all previous sends to finish
        send_reqs.wait_all();

        // post owner-buffer column-index array sends
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          send_reqs[i] = comm.isend(this->_owner_bufs.at(i).col_ind(),
            this->_owner_bufs.at(i).used_elements(), adp.get_owner_rank(i));
        }

        // wait for all receives and sends to finish
        recv_reqs.wait_all();
        send_reqs.wait_all();
      }

      /**
       * \brief Assembles the matrix structure of the algebraic-DOF-partitioned matrix.
       *
       * \param[in] local_matrix
       * The local matrix that acts as a layout template.\n
       * Its layout must be initialized, but its numerical values are ignored.
       */
      void _assemble_structure(const LocalMatrixType& local_matrix)
      {
        // get our partitioning
        const AlgDofPartiType& adp = *this->_alg_dof_parti;

        // get our matrix dimensions:
        // * each row of our matrix corresponds to one owned DOF
        // * each column corresponds to one global DOF
        const Index num_rows = adp.get_num_owned_dofs();
        const Index num_cols = adp.get_num_global_dofs();

        // Unfortunately, assembling the matrix structure is not that easy,
        // because it is a union of our owned matrix structure and all of
        // our donee matrix buffers. We don't have a chance to determine how
        // many duplicate (i,j) pairs we will encounter here, so we use the
        // DynamicGraph class here to keep things simple.

        // create a dynamic graph for our matrix
        Adjacency::DynamicGraph dynamic_graph(num_rows, num_cols);

        // First, let us add all local matrix rows which correspond to the
        // owned DOFs of this process, which are given by the "owned mirror"
        const IndexType* loc_row_ptr = local_matrix.row_ptr();
        const IndexType* loc_col_idx = local_matrix.col_ind();
        const IndexType* own_mir_idx = adp.get_owned_mirror().indices();
        for(Index i(0); i < num_rows; ++i)
        {
          // get the local DOF index of our i-th owned DOF
          const Index lrow = own_mir_idx[i];
          for(Index j(loc_row_ptr[lrow]); j < loc_row_ptr[lrow + 1]; ++j)
          {
            // translate the local column indices into global DOF indices
            dynamic_graph.insert(i, adp.map_local_to_global_index(loc_col_idx[j]));
          }
        }

        // process all donee-neighbor buffers
        const Index num_neigh_donee = adp.get_num_donee_neighbors();
        for(Index ineigh(0); ineigh < num_neigh_donee; ++ineigh)
        {
          // get the mirror and the matrix buffer of that donee neighbor
          // note: the donee mirror indices are owned DOF indices
          const MirrorType& mir = adp.get_donee_mirror(ineigh);
          const BufferMatrixType& buf = this->_donee_bufs.at(ineigh);
          XASSERT(mir.num_indices() == buf.rows());
          const Index num_idx = mir.num_indices();
          const IndexType* mir_idx = mir.indices();
          const IndexType* buf_ptr = buf.row_ptr();
          const IndexType* buf_idx = buf.col_ind();

          // loop over all matrix buffer rows
          for(Index i(0); i < num_idx; ++i)
          {
            // the owned DOF index of our the buffer matrix row
            const Index row = mir_idx[i];

            // loop over all buffer indices
            for(Index j(buf_ptr[i]); j < buf_ptr[i + 1]; ++j)
            {
              // insert entry into our local matrix
              dynamic_graph.insert(row, buf_idx[j]);
            }
          }
        }

        // render into a standard graph
        Adjacency::Graph graph(Adjacency::RenderType::as_is, dynamic_graph);

        // create our dof-partitoned matrix object
        this->_matrix = OwnedMatrixType(graph);
      }

      /**
       * \brief Auxiliary function: assembles a data-mirror for a donee-neighbor
       *
       * \param[in] local_matrix
       * The local matrix that acts as a layout template.\n
       * Its layout must be initialized, but its numerical values are ignored.
       *
       * \param[in] mirror
       * The mirror of the donee-neighbor.
       *
       * \param[in] buffer
       * The matrix buffer for the donee-mirror.
       */
      static MirrorType _asm_donee_data_mir(const LocalMatrixType& local_matrix,
        const MirrorType& mirror, const BufferMatrixType& buffer)
      {
        // allocate mirror
        const Index buf_rows = buffer.rows();
        const Index num_idx = buffer.used_elements();
        MirrorType dat_mir(local_matrix.used_elements(), num_idx);

        const IndexType* row_idx = mirror.indices();
        const IndexType* row_ptr = local_matrix.row_ptr();
        const IndexType* col_idx = local_matrix.col_ind();
        const IndexType* buf_ptr = buffer.row_ptr();
        const IndexType* buf_idx = buffer.col_ind();
        IndexType* dat_idx = dat_mir.indices();

        // loop over all buffer matrix rows
        for(Index i(0); i < buf_rows; ++i)
        {
          // get local matrix row index
          const Index row = row_idx[i];
          for(Index j(buf_ptr[i]); j < buf_ptr[i + 1]; ++j)
          {
            // get global column index
            const Index col = buf_idx[j];

            // loop over the corresponding row of our matrix
            for(Index k(row_ptr[row]); k < row_ptr[row + 1]; ++k)
            {
              // is this the column we are looking for?
              if(col_idx[k] == col)
              {
                // Yup, store this index in the mirror
                dat_idx[j] = k;
                break;
              }
            }
          }
        }

        return dat_mir;
      }

      /**
       * \brief Auxiliary function: assembles a data-mirror for an owner-neighbor
       *
       * \param[in] local_matrix
       * The local matrix that acts as a layout template.\n
       * Its layout must be initialized, but its numerical values are ignored.
       *
       * \param[in] mirror
       * The mirror of the owner-neighbor.
       *
       * \param[in] buffer
       * The matrix buffer for the owner-mirror.
       *
       * \param[in] glob_dof_idx
       * A vector containing the global DOF indices for all local DOFs.
       */
      static MirrorType _asm_owner_data_mir(const LocalMatrixType& local_matrix,
        const MirrorType& mirror, const BufferMatrixType& buffer,
        const std::vector<Index>& glob_dof_idx)
      {
        // allocate mirror
        const Index buf_rows = buffer.rows();
        const Index num_idx = buffer.used_elements();
        MirrorType dat_mir(local_matrix.used_elements(), num_idx);

        const IndexType* row_idx = mirror.indices();
        const IndexType* row_ptr = local_matrix.row_ptr();
        const IndexType* col_idx = local_matrix.col_ind();
        const IndexType* buf_ptr = buffer.row_ptr();
        const IndexType* buf_idx = buffer.col_ind();
        IndexType* dat_idx = dat_mir.indices();

        // loop over all buffer matrix rows
        for(Index i(0); i < buf_rows; ++i)
        {
          // get local matrix row index
          const Index row = row_idx[i];
          for(Index j(buf_ptr[i]); j < buf_ptr[i + 1]; ++j)
          {
            // get global column index
            const Index col = buf_idx[j];

            // loop over the corresponding row of our matrix
            for(Index k(row_ptr[row]); k < row_ptr[row + 1]; ++k)
            {
              // is this the column we are looking for?
              if(glob_dof_idx[col_idx[k]] == col)
              {
                // Yup, store this index in the mirror
                dat_idx[j] = k;
                break;
              }
            }
          }
        }

        return dat_mir;
      }

      /**
       * \brief Assembles the data-mirrors for all owner- and donee-neighbors.
       *
       * This functions assembles the data-mirrors, which are used for the communication
       * of the numerical values of a matrix in the #upload_numeric() function.
       * This allows that the matrix values can be uploaded and synchronized without the
       * need to re-communicate the unchanged matrix buffer layouts every time.
       *
       * \param[in] local_matrix
       * The local matrix that acts as a layout template.\n
       * Its layout must be initialized, but its numerical values are ignored.
       */
      void _assemble_data_mirrors(const LocalMatrixType& local_matrix)
      {
        // get our partitioning
        const AlgDofPartiType& adp = *this->_alg_dof_parti;

        // assemble donee data mirrors
        this->_data_donee_mirs.resize(this->_donee_bufs.size());
        for(Index i(0); i < this->_donee_bufs.size(); ++i)
          this->_data_donee_mirs.at(i) = _asm_donee_data_mir(this->_matrix,
            adp.get_donee_mirror(i), this->_donee_bufs.at(i));

        // assemble owner data mirrors
        this->_data_owner_mirs.resize(this->_owner_bufs.size());
        for(Index i(0); i < this->_owner_bufs.size(); ++i)
          this->_data_owner_mirs.at(i) = _asm_owner_data_mir(local_matrix,
            adp.get_owner_mirror(i), this->_owner_bufs.at(i),
            adp.get_global_dof_indices());

        // get matrix arrays
        const IndexType* row_ptr_x = this->_matrix.row_ptr();
        const IndexType* col_idx_x = this->_matrix.col_ind();
        const IndexType* row_ptr_a = local_matrix.row_ptr();
        const IndexType* col_idx_a = local_matrix.col_ind();

        // reset my data mirror
        _my_data_mir.clear();
        _my_data_mir.reserve(this->_matrix.used_elements());

        const Index num_owned_dofs = adp.get_num_owned_dofs();
        const IndexType* loc_dof_idx = adp.get_owned_mirror().indices();
        const std::vector<Index>& glob_dof_idx = adp.get_global_dof_indices();

        XASSERT(num_owned_dofs == this->_matrix.rows());

        // loop over all our owned DOFS
        for(Index own_dof(0); own_dof < num_owned_dofs; ++own_dof)
        {
          // get the local DOF index for this owned DOF
          const Index loc_dof = loc_dof_idx[own_dof];

          // loop over all columns of our owned DOF matrix
          for(Index j(row_ptr_x[own_dof]); j < row_ptr_x[own_dof + 1]; ++j)
          {
            // get the column index, which is a global DOF index
            const Index col_x = col_idx_x[j];

            // loop over the row of the local matrix
            for(Index k(row_ptr_a[loc_dof]); k < row_ptr_a[loc_dof + 1]; ++k)
            {
              // get the local column index and map it to a global one
              const Index col_a = glob_dof_idx[col_idx_a[k]];

              // match?
              if(col_a == col_x)
              {
                // Okay, add this matrix entry to our mirror
                _my_data_mir.push_back(std::make_pair(j, k));
                break;
              }
            }
          }
        }
      }

      /**
       * \brief Uploads the values of the local matrix into our matrix partition
       *
       * \param[in] local_matrix
       * The local matrix whose values are to be uploaded.
       */
      void _upload_own(const LocalMatrixType& local_matrix)
      {
        // get our data arrays
        const DataType* val_a = local_matrix.val();
        DataType* val_x = this->_matrix.val();

        // loop over our own data mirror
        for(auto ij : this->_my_data_mir)
          val_x[ij.first] = val_a[ij.second];
      }

      /**
       * \brief Gathers the values of a local matrix into an owner-neighbor buffer-vector
       *
       * \param[out] buffer
       * A buffer vector that receives the gathered values.
       *
       * \param[in] local_matrix
       * The local matrix whose values are to be gathered.
       *
       * \param[in] data_mirror
       * The owner-data mirror to be used.
       */
      static void _gather_data(BufferVectorType& buffer, const LocalMatrixType& local_matrix,
        const MirrorType& data_mirror)
      {
        XASSERT(buffer.size() == data_mirror.num_indices());
        XASSERT(data_mirror.size() == local_matrix.used_elements());

        DataType* buf_val = buffer.elements();
        const DataType* mat_val = local_matrix.val();
        const IndexType* mir_idx = data_mirror.indices();
        const Index n = buffer.size();
        for(Index i(0); i < n; ++i)
          buf_val[i] = mat_val[mir_idx[i]];
      }

      /**
       * \brief Scatters the values of a donee-neighbor buffer-vector into a local matrix
       *
       * \param[inout] local_matrix
       * The local matrix that receives the scattered values.
       *
       * \param[in] buffer
       * The buffer vector whose values are to be scattered.
       *
       * \param[in] data_mirror
       * The donee-data mirror to be used.
       */
      static void _scatter_data(LocalMatrixType& local_matrix, const BufferVectorType& buffer,
        const MirrorType& data_mirror)
      {
        XASSERT(buffer.size() == data_mirror.num_indices());
        XASSERT(data_mirror.size() == local_matrix.used_elements());

        DataType* mat_val = local_matrix.val();
        const DataType* buf_val = buffer.elements();
        const IndexType* mir_idx = data_mirror.indices();
        const Index n = buffer.size();
        for(Index i(0); i < n; ++i)
          mat_val[mir_idx[i]] += buf_val[i];
      }
    }; // class AlgDofPartiMatrixCSR
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_ALG_DOF_PARTI_HPP
