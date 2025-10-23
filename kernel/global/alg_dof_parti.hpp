// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/global/gate.hpp>
#include <kernel/global/vector.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/tuple_mirror.hpp>
#include <kernel/util/omp_util.hpp>

#include <vector>

namespace FEAT
{
  namespace Global
  {
    /// \cond internal
    namespace Intern
    {
      template<typename IdxVector_>
      class ADPAux;
    } // namespace Intern
    /// \endcond

    /**
     * \brief Algebraic DOF Partitioning implementation
     *
     * This class implements the functionality to generate and manage a distributed partitioning
     * of the degrees of freedom of a (linear) system in an algebraic sense instead of the usual
     * domain-decomposition approach. With this <b>algebraic dof partitioning approach (ADP)</b>,
     * each process in a given communicator becomes the unique owner of a set of consecutively
     * enumerated global dofs, whereas in the domain-decomposition approach a single global dof
     * may be shared by several processes, thus having no unique owner.
     * The primary purpose of this class and the closely related AlgDofPartiSystem class is to
     * offer the capability of embedding various third-party linear solver libraries as solvers,
     * preconditioners or smoothers within the solver framework used in FEAT.
     *
     * To use this class, one needs to have a fully assembled Global::Gate object, which contains
     * all the partitioning information in the usual domain-decomposition sense, which is used as
     * a "starting point" for an algebraic DOF partitioning in the assemble_by_gate() function.
     *
     * Throughout this class and the closely related AlgDofPartiSystem class, you will frequently
     * find the terms "local DOF", "shared DOF" and "owned DOF" and "donee DOF", which are defined
     * as follows:
     * - "local" DOFs always means the DOFs as they are used in the domain decomposition world.
     * - "shared" DOFs are all local DOFs are "shared" by at least two processes in the domain
     *   decomposition world. In other words: Each local DOF, which is part of at least one
     *   neighbor mirror in the Global::Gate, is a shared DOF.
     * - "owned" DOFs are all local DOFs, which:
     *   - are not shared with any other process
     *   - are only shared with processes of higher comm rank
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
     * <u><b>Algebraic DOF partitioning of combined discretizations:</b></u>\n
     * This class supports the algebraic dof partitioning of simple scalar FE spaces, which are
     * represented by the LAFEM::DenseVector class, as well as more complex space combinations,
     * which are represented by a hierarchy of meta containers, such as the LAFEM::TupleVector,
     * and basic containers like LAFEM::DenseVector and/or LAFEM::DenseVectorBlocked; one prominent
     * example is the velocity-pressure space pair used in Stokes systems which is typically
     * represented by a LAFEM::DenseVectorBlocked for the velocity components and a (scalar)
     * LAFEM::DenseVector for the pressure component, which are then joined together by a
     * LAFEM::TupleVector meta container to form a single class that contains the entire
     * velocity-pressure pair.
     *
     * In the case of a combined discretization (such as the velocity-pressure pair), the algebraic
     * dof partitioning distributes the DOFs of the individual components in an \b interleaved
     * fashion among the processes, i.e. each process owns a set of DOFs of the first component
     * followed by a set of DOFs of the second component, etc., so in the example of the
     * velocity-pressure pair each process owns a set of consecutively numbered velocity DOFs
     * followed by a set of consecutively numbered pressure DOFs, and so from a global perspective,
     * the global DOFs alternate between sets of velocity and pressure DOFs across all processes.
     *
     * <b>Example:</b>\n
     * Let's assume that we have 25x2=50 velocity DOFs (V0.0, V0.1, V1.0, ..., V24.0, V24.1) and 12
     * pressure DOFs (P0, ..., P11) in total partitioned among 4 processes and let's assume that each
     * process owns 3 pressure DOFs and that the four processes own 18, 12, 12 and 8 velocity dofs,
     * respectively, then the algebraic DOF partitioning of the 62 total global DOFs would look like
     * in the table below:
       \verbatim
       +------------------------------+-------------------------------+--------------------------------+----------------------------------+
       |          Process 0           |           Process 1           |            Process 2           |            Process 3             |
       +------------------------------+-------------------------------+--------------------------------+----------------------------------+
       |    0, ...,    17, 18, 19, 20 |    21, ...,    32, 33, 34, 35 |     36, ...,    47, 48, 49, 50 |    51, ...,     58, 59,  60,  61 | <- global algebraic DOF
       +------------------------------+-------------------------------+--------------------------------+----------------------------------+
       |    0, ...,    17, 18, 19, 20 |     0, ...,    11, 12, 13, 14 |      0, ...,    11, 12, 13, 14 |     0, ...,      8,  9,  10,  11 | <- owned (local) algebraic DOF
       +------------------------------+-------------------------------+--------------------------------+----------------------------------+
       | V0.0, ...,  V8.1             |  V9.0, ..., V14.1             |  V15.0, ..., V20.1             |  V21.0, ..., V24.1               | <- global velocity DOF (interpreted as blocked DOFs)
       |  V'0, ...,  V'17             |  V'18, ...,  V'29             |   V'30, ...,  V'41             |   V'42, ...,  V'49               | <- global velocity DOF (interpreted as scalar DOFs)
       |                   P0, P1, P2 |                    P3, P4, P5 |                     P6, P7, P8 |                     P9, P10, P11 | <- global pressure DOF (scalar)
       +------------------------------+-------------------------------+--------------------------------+----------------------------------+
       | v0.0, ...,  v8.1             |  v0.0, ...,  v5.1             |   v0.0, ...,  v5.1             |   v0.0, ...,  v4.1               | <- owned (local) velocity DOF (interpreted as blocked DOFs)
       |  v'0, ...,  v'17             |   v'0, ...,  v'11             |    v'0, ...,  v'11             |    v'0, ...,   v'7               | <- owned (local) velocity DOF (interpreted as scalar DOFs)
       |                   p0, p1, p2 |                    p0, p1, p2 |                     p0, p1, p2 |                     p0,  p1,  p2 | <- owned (local) pressure DOF
       +------------------------------+-------------------------------+--------------------------------+----------------------------------+
       \endverbatim
     * So in the above example, global algebraic DOF 33 corresponds to global pressure DOF 3 (P3),
     * which is owned by process 1 and which corresponds to the owned pressure DOF 1 (p1) on that
     * process, whereas global DOF 47 corresponds to the second block-component (index 1) of the
     * global velocity DOF 20 (V20.1), which is owned by process 2 and which corresponds to the
     * owned velocity DOF 5 (v5.1), etc.
     *
     * <u><b>Retrieving ADP block information for combined discretizations:</b></u>\n
     * Some third-party libraries, which are used as solver backends for ADP based solvers, do not
     * treat the system as a black box, but require information about the composition of the
     * system components within the algebraic DOF partitioning -- or in simpler terms: the library
     * needs to know which DOF is a velocity DOF and which DOF is a pressure DOF in the example above.
     * This information can be queried by calling the get_block_information() member function, which
     * returns a String object containing the relevant offsets and counts encoded in XML format.
     *
     * The returned XML string consists of one or more lines and each line contains a single XML marker,
     * which currently is either
     * - a closed <c>\<Scalar .../\></c> marker that represents a scalar component that corresponds
     *   to a LAFEM::DenseVector container
     * - a closed <c>\<Blocked .../\></c> marker that represents a blocked/vector-valued component
     *   that corresponds to a LAFEM::DenseVectorBlocked container
     * - an open <c>\<Tuple ...\></c> marker that represents a meta component that combines one or
     *   more other components and which corresponds to the LAFEM::TupleVector container
     * - a <c>\</Tuple\></c> marker that closes a previously opened <c>\<Tuple ...\></c> marker
     *
     * Each of the above XML markers contains the following five attributes in alphabetical order:
     * - <c>gc</c> ("global count"): contains the total number of global DOFs on all processes of the
     *   corresponding component; this quantity is identical on all processes
     * - <c>gf</c> ("global first"): represents the index of the first DOF of this component that is
     *   is owned by the current process
     * - <c>go</c> ("global offset"): represents the index of the first global DOF owned by this
     *   process, which corresponds to the first global DOF of the corresponding component
     * - <c>lc</c> ("local count"): contains the number of local DOFs within this component that are
     *   owned by the current process.
     * - <c>lo</c> ("local offset"): represents the index of the first local DOF owned by this
     *   process, which corresponds to the first local DOF of the corresponding component
     *
     * Additionally, the <c>\<Blocked .../\></c> marker also contains the following attribute:
     * - <c>bs</c> ("block size"): contains the total number values per DOF block
     *
     * In the example above, the four processes would contain the following block information strings:
     * - Process 0:
       \verbatim
       <Tuple gc="62" gf="0" go="0" lc="21" lo="0">
       <Blocked bs="2" gc="50" gf="0" go="0" lc="18" lo="0"/>
       <Scalar gc="12" gf="0" go="18" lc="3" lo="18"/>
       </Tuple>
       \endverbatim
     * - Process 1:
       \verbatim
       <Tuple gc="62" gf="21" go="21" lc="15" lo="0">
       <Blocked bs="2" gc="50" gf="18" go="21" lc="12" lo="0"/>
       <Scalar gc="12" gf="3" go="33" lc="3" lo="12"/>
       </Tuple>
       \endverbatim
     * - Process 2:
       \verbatim
       <Tuple gc="62" gf="36" go="36" lc="15" lo="0">
       <Blocked bs="2" gc="50" gf="30" go="36" lc="12" lo="0"/>
       <Scalar gc="12" gf="6" go="48" lc="3" lo="12"/>
       </Tuple>
       \endverbatim
     * - Process 3:
       \verbatim
       <Tuple gc="62" gf="51" go="51" lc="11" lo="0">
       <Blocked bs="2" gc="50" gf="42" go="51" lc="8" lo="0"/>
       <Scalar gc="12" gf="9" go="59" lc="3" lo="8"/>
       </Tuple>
       \endverbatim
     *
     * \todo decide what happens if this process owns no dofs at all
     * \todo make sure that this works for discontinuous elements, too
     * \todo ensure that extra local and extra shared dofs are handled appropriately
     *
     * \author Peter Zajac
     */
    template<typename LocalVector_, typename Mirror_>
    class AlgDofParti
    {
    public:
      /// the local vector type
      typedef LocalVector_ LocalVectorType;
      /// the vector mirror type
      typedef Mirror_ MirrorType;

      /// our data type
      typedef typename LocalVectorType::DataType DataType;
      /// our index type
      typedef typename LocalVectorType::IndexType IndexType;

      /// the global vector type
      typedef Global::Vector<LocalVectorType, MirrorType> GlobalVectorType;
      /// the global gate type
      typedef Global::Gate<LocalVectorType, MirrorType> GateType;

      /// the index vector type, this one is required for internal computations
      typedef typename LocalVectorType::template ContainerTypeByDI<IndexType, IndexType> IndexVectorType;

      /// auxiliary helper class type
      typedef Intern::ADPAux<IndexVectorType> ADPAuxType;

      /// the type for buffers used by our mirrors
      typedef LAFEM::DenseVector<DataType, IndexType> BufferVectorType;

      /// type of mirror for owned local dof indices
      typedef MirrorType OwnedMirrorType;
      /// type of mirror for neighbor owner local dof indices
      typedef MirrorType OwnerMirrorType;
      /// type of mirror for neighbor donee owned dof indices
      typedef LAFEM::VectorMirror<DataType, IndexType> DoneeMirrorType;

    //protected:
      /// our communicator
      const Dist::Comm* _comm;
      /// global dof offset of this process
      Index _global_dof_offset;
      /// global dof count over all processes
      Index _global_dof_count;
      /// number of local dofs of this process
      Index _local_dof_count;
      /// number of owned dofs of this process
      Index _owned_dof_count;
      /// number of extra local dofs; may be different for each process
      Index _extra_local_dof_count;
      /// number of extra shared dofs; must be equal on all processes
      Index _extra_shared_dof_count;
      /// offset of first extra DOF
      Index _extra_dof_offset;
      /// global dof indices for each local dof
      IndexVectorType _global_dof_idx;
      /// mirror for this process's owned DOFs
      MirrorType _owned_mirror;
      /// ranks of DOF-owner neighbor processes
      std::vector<int> _owner_ranks;
      /// ranks of DOF-donee neighbor processes
      std::vector<int> _donee_ranks;
      /// rank/mirror-pair of DOF-owner processes
      std::vector<OwnerMirrorType> _owner_mirrors;
      /// rank/mirror-pair of DOF-donee processes
      std::vector<DoneeMirrorType> _donee_mirrors;
      /// a string containing the block information in XML format
      String _block_information;

      // Note: the following vectors hold int objects, because these are used as recvcount in an
      // MPI_Allgatherv call, which only accepts int objects as counts; we have an XASSERT in
      // place to ensure that we do not overflow the int datatype with more than a billion DOFs

      /// global DOF offsets of all processes, \see #allgather()
      std::vector<int> _all_global_dof_offset;
      /// global DOF counts of all processes, \see #allgather()
      std::vector<int> _all_global_dof_counts;

    public:
      /// default constructor
      AlgDofParti() :
        _comm(nullptr),
        _global_dof_offset(0),
        _global_dof_count(0),
        _local_dof_count(0),
        _owned_dof_count(0),
        _extra_local_dof_count(0),
        _extra_shared_dof_count(0),
        _extra_dof_offset(0)
      {
      }

      /// no copy, no problems
      AlgDofParti(const AlgDofParti&) = delete;
      /// no copy, no problems
      AlgDofParti& operator=(const AlgDofParti&) = delete;

      /// virtual destructor
      virtual ~AlgDofParti() = default;

      /// Resets the whole object
      void clear()
      {
        _comm = nullptr;
        _global_dof_offset = IndexType(0);
        _global_dof_count = IndexType(0);
        _local_dof_count = IndexType(0);
        _owned_dof_count = IndexType(0);
        _extra_local_dof_count = IndexType(0);
        _extra_shared_dof_count = IndexType(0);
        _extra_dof_offset = Index(0);
        _global_dof_idx.clear();
        _owned_mirror.clear();
        _owner_mirrors.clear();
        _donee_mirrors.clear();
        _all_global_dof_offset.clear();
        _all_global_dof_counts.clear();
      }

      /// \returns The communicator
      const Dist::Comm* get_comm() const
      {
        return _comm;
      }

      /// \returns The size of this object in bytes
      std::size_t bytes() const
      {
        std::size_t r = 0u;
        r += _global_dof_idx.bytes();
        r += _owned_mirror.bytes();
        for(const auto& x : _owner_mirrors)
          r += x.second.bytes();
        for(const auto& x : _donee_mirrors)
          r += x.second.bytes();
        r += _all_global_dof_offset.size() * 4u;
        r += _all_global_dof_counts.size() * 4u;
        return r;
      }

      /**
       * \brief Returns the block information of the algebraic dof partitioning as an XML string
       *
       */
      String get_block_information() const
      {
        return this->_block_information;
      }

      /// \returns The number of local DOFs owned by this process.
      Index get_num_owned_dofs() const
      {
        return this->_owned_dof_count;
      }

      /// \returns The number of local DOFs shared by this process.
      Index get_num_local_dofs() const
      {
        return this->_local_dof_count;
      }

      /// \returns The total number of global DOFs.
      Index get_num_global_dofs() const
      {
        return this->_global_dof_count;
      }

      /// \returns The index of the first global DOF owned by this process.
      Index get_global_dof_offset() const
      {
        return this->_global_dof_offset;
      }

      /**
       * \brief Returns the global-dof-indices array.
       */
      const IndexVectorType& get_global_dof_indices() const
      {
        return this->_global_dof_idx;
      }

      /**
       * \brief Returns the Mirror of the owned DOFs.
       */
      const OwnedMirrorType& get_owned_mirror() const
      {
        return this->_owned_mirror;
      }

      /// \returns The number of owner neighbors.
      Index get_num_owner_neighbors() const
      {
        return Index(this->_owner_mirrors.size());
      }

      /// \returns The number of donee neighbors.
      Index get_num_donee_neighbors() const
      {
        return Index(this->_donee_mirrors.size());
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
        return this->_owner_ranks.at(i);
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
        return this->_donee_ranks.at(i);
      }

      /**
       * \brief Returns an owner neighbor mirror.
       *
       * \note
       * The indices of an owner neighbor mirror are local DOF indices.
       *
       * \param[in] i
       * The index of the owner neighbor.
       *
       * \returns
       * The mirror of the i-th owner neighbor.
       */
      const OwnerMirrorType& get_owner_mirror(Index i) const
      {
        return this->_owner_mirrors.at(i);
      }

      /**
       * \brief Returns a donee neighbor mirror.
       *
       * \note
       * The indices of a donee neighbor mirror are owned DOF indices.
       *
       * \param[in] i
       * The index of the donee neighbor.
       *
       * \returns
       * The mirror of the i-th donee neighbor.
       */
      const DoneeMirrorType& get_donee_mirror(Index i) const
      {
        return this->_donee_mirrors.at(i);
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
        XASSERTM(!this->_all_global_dof_offset.empty(), "You did not ask to assemble the global dof offsets");
        return this->_all_global_dof_offset;
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
        XASSERTM(!this->_all_global_dof_counts.empty(), "You did not ask to assemble the global dof counts");
        return this->_all_global_dof_counts;
      }

      /**
       * \brief Assembles the required data for the AlgDofPartiSystem::apply() function.
       *
       * \warning
       * Calling this function and (more importantly) the AlgDofPartiSystem::apply() functions will *totally*
       * screw up the scalability of your application, so the code *will* blow up in your face on bigger clusters.\n
       * Do not use this function except for <b>small-scale debugging</b> purposes. You have been warned.
       *
       * \param[in] yes_i_really_want_to_do_this
       * An assertion will fire if this is \c false. The only purpose of this parameter is
       * to force you to read the warning in this documentation first.
       */
      void assemble_allgather(bool yes_i_really_want_to_do_this = false)
      {
        XASSERTM(yes_i_really_want_to_do_this, "You probably don't want to do this!");

        // ensure that we don't have more than a billion DOFs
        XASSERTM(this->get_num_global_dofs() < 1'000'000'000,
          "allgather assembly not possible for more than a billion DOFs!");

        std::size_t num_ranks = std::size_t(this->_comm->size());

        int my_nod = int(this->get_num_owned_dofs());
        int my_off = int(this->get_global_dof_offset());

        // resize vectors
        this->_all_global_dof_offset.resize(num_ranks, 0);
        this->_all_global_dof_counts.resize(num_ranks, 0);

        // allgather offsets + counts
        this->_comm->allgather(&my_off, 1, this->_all_global_dof_offset.data(), 1);
        this->_comm->allgather(&my_nod, 1, this->_all_global_dof_counts.data(), 1);
      }

      /**
       * \brief Assembles the AlgDofParti object from a given Global::Gate.
       *
       * This function performs the assembly of all required internal data structures.
       *
       * \param[in] extra_local
       * The number of extra local dofs for this process.
       *
       * \param[in] extra_shared
       * The number of extra shared dofs for all processes.
       */
      void assemble_by_gate(const GateType& gate, Index extra_local = 0u, Index extra_shared = 0u)
      {
        // set our communicator
        this->_comm = gate.get_comm();
        XASSERTM(this->_comm != nullptr, "need a communicator here!");

        // get the communicator
        const Dist::Comm& comm = *this->_comm;

        // get my rank and neighbor ranks
        const int my_rank = comm.rank();
        const std::vector<int>& neighbor_ranks = gate._ranks;

        // get a template vector
        const LocalVectorType& vec_tmpl = gate.get_freqs();

        // get number of total local DOFs on our patch
        this->_local_dof_count = vec_tmpl.template size<LAFEM::Perspective::pod>();
        this->_extra_local_dof_count = extra_local;
        this->_extra_shared_dof_count = extra_shared;
        //this->_local_dof_count += extra_local + extra_shared;

        // As a very first step, we have to decide which process will become the owner of each DOF.
        // In this implementation, each DOF is owned by by the process with the lowest rank, which
        // simplifies a lot of things.

        // So, first create a vector, which stores the rank of the owner process of each of our
        // local DOFs and format the vector to our rank, i.e. at startup assume that we own all
        // local DOFs.
        IndexVectorType dof_owners;
        ADPAuxType::alloc_idx_vector(dof_owners, vec_tmpl, IndexType(my_rank));

        // Now loop over all our neighbor processes
        for(std::size_t i(0); i < neighbor_ranks.size(); ++i)
        {
          // Check whether the neighbors rank is less than our rank, as otherwise that particular
          // neighbor can not own any of our local DOFs.
          if(neighbor_ranks.at(i) < my_rank)
          {
            // this neighbor has a lower rank than our process, so update the dof owners vector and
            // update the ownership of these dofs by setting this neighbor as  their owner -- unless
            // another neighbor with an even lower rank is already a candidate for ownership
            ADPAuxType::update_owners(neighbor_ranks.at(i), dof_owners, gate.get_mirrors().at(i));
          }
        }

        // Okay, at this point we (and all other processes) know the owners of each of our local DOFs.
        // Now, we have to generate the mirror of all local DOF indices of the DOFs that we own:
        this->_owned_dof_count = ADPAuxType::build_owned_mirror(my_rank, this->_owned_mirror, dof_owners);

        // Save current owned DOF count as extra DOF offset
        this->_extra_dof_offset = this->_owned_dof_count;

        // Now add the extra local DOFs for each process
        this->_owned_dof_count += _extra_local_dof_count;

        // The shared DOFs are always owned by the last process
        const bool is_last_process = (this->_comm->rank()+1 == this->_comm->size());
        if(is_last_process)
          this->_owned_dof_count += _extra_shared_dof_count;

        /// \todo what do we do if owned_dofs==0 ?
        XASSERTM(this->_owned_dof_count > Index(0), "this process has no DOFs!");

        // Next, we have to determine the global offset of our first owned DOF, which is easily
        // obtained by an exclusive-scan over each processes owned DOF count:
        comm.exscan(&this->_owned_dof_count, &this->_global_dof_offset, 1, Dist::op_sum);

        // Furthermore, we also perform an allreduce to obtain the total number of global DOFs,
        // which is generally helpful and can be used for sanity checks.
        comm.allreduce(&this->_owned_dof_count, &this->_global_dof_count, 1, Dist::op_sum);

        // Allocate an index vector that for each local DOF contains the "owned DOF index" if we
        // own this particular DOF or ~Index(0) if this particular DOF is owned by another process:
        IndexVectorType own_dof_idx;
        ADPAuxType::alloc_idx_vector(own_dof_idx, vec_tmpl, ~IndexType(0));
        ADPAuxType::build_owned_dofs(my_rank, own_dof_idx, dof_owners, 0u);

        // Now we also have to create the mirrors for each of our neighbor processes:
        for(std::size_t i(0); i < neighbor_ranks.size(); ++i)
        {
          // get the neighbor rank and its mirror
          const int neighbor_rank = neighbor_ranks.at(i);
          const MirrorType& halo_mirror = gate.get_mirrors().at(i);

          // There are two cases now:
          // 1) The neighbor process has a lower rank, so it is a (potential) owner
          //    of at least one of our local DOFs.
          // 2) The neighbor process has a greater rank, so we are an owner of at least
          //    one of the neighbor process's local DOFs.
          if(neighbor_rank < my_rank)
          {
            // This is a potential "owner-neighbor", which may own one or several of our local DOFs.
            // We call a helper function, which will create a vector of all *local* DOF indices,
            // which are owned by that neighbor process. Note that it is perfectly legit if the
            // vector is empty, as this simply means that all our local DOFs, which are shared with
            // that particular neighbor, are owned by some *other* neighbor process(es) and not by
            // the particular neighbor that we are currently considering.
            OwnerMirrorType owner_mirror;
            if(ADPAuxType::build_owner_mirror(neighbor_rank, owner_mirror, halo_mirror, dof_owners) > Index(0))
            {
              // That neighbor owns at least one of our local DOFs, so create a mirror from the
              // index-vector and push it into our list of "owner-neighbors".
              // Note that the indices in the mirror are *local* DOF indices.
              _owner_ranks.push_back(neighbor_rank);
              _owner_mirrors.emplace_back(std::move(owner_mirror));
            }
          }
          else
          {
            // This is a potential "donee-neighbor", which shares one or several of our local DOFs,
            // which we might own. We call a helper function, which will create a vector of all our
            // *owned* DOF indices, which are shared  with that particular neighbor process.
            // Note that it is perfectly legit if the vector is empty, as this simply means that
            // all our local DOFs, which are shared with that particular neighbor, are owned by
            // some *other* neighbor process(es) and not by this process (i.e. us).
            Index num_donee_dofs = ADPAuxType::count_donee_dofs(halo_mirror, own_dof_idx);
            if(num_donee_dofs > Index(0))
            {
              // allocate and build the actual donee mirror
              DoneeMirrorType donee_mirror(this->_owned_dof_count, num_donee_dofs);
              ADPAuxType::build_donee_mirror(donee_mirror, halo_mirror, own_dof_idx, 0u);

              // We own at least of one of the neighbor's local DOFs, so create a mirror from the
              // index-vector and push it into our list of "donee-neighbors".
              // Note that the indices in the mirror are *owned* DOF indices.
              _donee_ranks.push_back(neighbor_rank);
              _donee_mirrors.emplace_back(std::move(donee_mirror));
            }
          }
        }

        // Finally, we have to determine the global DOF index for each of our local DOFs, so that
        // we can perform a "local-to-global" DOF index lookup. This information is required for
        // the assembly of matrices later on, so we have to store this in a member-variable vector.

        // Let's create the vector of the required length, initialize all its indices to ~0, which
        // can be used for a sanity check later on to ensure that we know the global indices of all
        // our local DOFs:
        ADPAuxType::alloc_idx_vector(this->_global_dof_idx, vec_tmpl, ~IndexType(0));

        // Next, let's loop over all DOFs that we own and assign their corresponding global DOF
        // indices, starting with our global DOF offset, which we have already determined before:
        ADPAuxType::init_global_dof_idx(this->_global_dof_idx, this->_owned_mirror, this->_global_dof_offset);

        // Finally, we also need to query the global DOF indices of all our local DOFs, which are
        // owned by some neighbor process (i.e. we are the donee for these DOFs).
        // Of course, we also have to send the global DOF indices of all DOFs that we own
        // to all of our "donee-neighbors".

        // Allocate index buffers and post receives for all of our owner-neighbors:
        const Index num_neigh_owner = Index(_owner_mirrors.size());
        Dist::RequestVector recv_reqs(num_neigh_owner);
        std::vector<std::vector<IndexType>> recv_bufs(num_neigh_owner);
        for(Index i(0); i < num_neigh_owner; ++i)
        {
          const std::size_t num_idx = _owner_mirrors.at(i).buffer_size(vec_tmpl);
          recv_bufs.at(i).resize(num_idx, ~IndexType(0));
          recv_reqs[i] = comm.irecv(recv_bufs.at(i).data(), num_idx, _owner_ranks.at(i));
        }

        // Allocate index buffers, fill them with the global DOF indices of our owned
        // DOFs and send them to the donee-neighbors:
        const Index num_neigh_donee = Index(_donee_mirrors.size());
        Dist::RequestVector send_reqs(num_neigh_donee);
        std::vector<std::vector<IndexType>> send_bufs(num_neigh_donee);
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          // Allocate buffer of required size and translate the indices of our mirror,
          // which are "owned DOF indices", to global DOF indices, where:
          //   global_dof_index := global_dof_offset + owned_dof_index
          const std::size_t num_idx = _donee_mirrors.at(i).num_indices();
          send_bufs.at(i).resize(num_idx, ~IndexType(0));
          _get_donee_dofs(send_bufs.at(i), _donee_mirrors.at(i), this->_global_dof_offset);
          send_reqs[i] = comm.isend(send_bufs.at(i).data(), num_idx, _donee_ranks.at(i));
        }

        // Now process all receive requests from our owner-neighbors:
        for(std::size_t i(0u); recv_reqs.wait_any(i); )
        {
          // Scatter the global DOF indices, which our friendly owner-neighbor has
          // provided us with, into our global DOF index vector:
          ADPAuxType::set_owner_dofs(this->_global_dof_idx, recv_bufs.at(i), _owner_mirrors.at(i), 0u);
        }

        // wait for all sends to finish
        send_reqs.wait_all();

        // get local counts
        std::vector<Index> local_num;
        ADPAuxType::get_local_dof_count(local_num, this->_global_dof_idx, this->_owned_mirror);

        // compute global counts
        std::vector<Index> global_num(local_num.size());
        comm.allreduce(local_num.data(), global_num.data(), local_num.size(), Dist::op_sum);

        // compute global offsets
        std::vector<Index> global_first(local_num.size());
        comm.exscan(local_num.data(), global_first.data(), local_num.size(), Dist::op_sum);

        // build deques
        std::deque<Index> loc_off, loc_num, glob_first, glob_num;
        for(Index i(0), off(0); i < Index(local_num.size()); ++i)
        {
          loc_off.push_back(off);
          loc_num.push_back(local_num[i]);
          glob_num.push_back(global_num[i]);
          glob_first.push_back(global_first[i]);
          off += local_num[i];
        }

        // finally, build the block information
        this->_block_information = ADPAuxType::build_block_info(this->_global_dof_idx, this->_global_dof_offset,
          glob_first, glob_num, loc_off, loc_num);
      }

      template<typename DT2_>
      void upload_vector(DT2_* owned_dofs, const LocalVectorType& local_vector,
        const DataType* extra_local = nullptr, const DataType* extra_shared = nullptr) const
      {
        // gather owned DOFs
        ADPAuxType::gather(owned_dofs, local_vector, this->_owned_mirror);

        // gather extra local DOFs
        if(this->_extra_local_dof_count > Index(0))
        {
          XASSERT(extra_local != nullptr);
          DT2_* dst = &owned_dofs[this->_extra_dof_offset];
          for(Index i = 0u; i < this->_extra_local_dof_count; ++i)
            dst[i] = DT2_(extra_local[i]);
        }

        // gather extra shared DOFs on last process
        if(this->_extra_shared_dof_count > Index(0))
        {
          XASSERT(extra_shared != nullptr);
          if(this->_comm->rank() + 1 == this->_comm->size())
          {
            DT2_* dst = &owned_dofs[this->_extra_dof_offset + this->_extra_local_dof_count];
            for(Index i = 0u; i < this->_extra_shared_dof_count; ++i)
              dst[i] = DT2_(extra_shared[i]);
          }
        }
      }

      template<typename DT2_>
      void download_vector(const DT2_* owned_dofs, LocalVectorType& local_vector,
        DataType* extra_local = nullptr, DataType* extra_shared = nullptr) const
      {
        // ensure we have a communicator
        XASSERT(this->_comm != nullptr);

        // get the number of owner- and donee-neighbors
        const std::size_t num_neigh_owner = this->_owner_mirrors.size();
        const std::size_t num_neigh_donee = this->_donee_mirrors.size();

        // The basic idea is simple:
        // 1) Send the shared DOFs to all of our donee-neighbors
        // 2) Receive the shared DOFs from all of our owner-neighbors
        // 3) Scatter our own owned DOFs into our local vector
        // 4) Scatter the dofs we received from our owner-neighbors

        // Create receive buffers and post receives for all owner-neighbors
        Dist::RequestVector recv_reqs(num_neigh_owner);
        std::vector<BufferVectorType> recv_bufs(num_neigh_owner);
        for(std::size_t i(0); i < num_neigh_owner; ++i)
        {
          // create a vector buffer
          // note: owner mirrors relate to local DOF indices
          recv_bufs.at(i) =  this->_owner_mirrors.at(i).create_buffer(local_vector);
          // post receive from owner neighbor
          recv_reqs[i] = this->_comm->irecv(recv_bufs.at(i).elements(), recv_bufs.at(i).size(), this->_owner_ranks.at(i));
        }

        // Create send buffers, fill them with our owned DOFs and send this
        // to our donee-neighbors
        Dist::RequestVector send_reqs(num_neigh_donee);
        std::vector<BufferVectorType> send_bufs(num_neigh_donee);
        for(Index i(0); i < num_neigh_donee; ++i)
        {
          // get mirror, create buffer and gather the shared DOFs
          // note: donee mirrors relate to owned DOF indices
          const DoneeMirrorType& mir = this->_donee_mirrors.at(i);
          send_bufs.at(i) = BufferVectorType(mir.num_indices());
          _gather(send_bufs.at(i), owned_dofs, mir);
          // post send to donee neighbor
          send_reqs[i] = this->_comm->isend(send_bufs.at(i).elements(), send_bufs.at(i).size(), this->_donee_ranks.at(i));
        }

        // format the output vector
        local_vector.format();

        // scatter our own owned DOFs into our output vector
        ADPAuxType::scatter(owned_dofs, local_vector, this->_owned_mirror);

        // process receives from our owner-neighbors
        for(std::size_t idx(0u); recv_reqs.wait_any(idx); )
        {
          // scatter received DOFs into our output vector
          ADPAuxType::scatter(recv_bufs.at(idx).elements(), local_vector, this->_owner_mirrors.at(idx));
        }

        // at this point, all receives should have finished
        XASSERT(recv_reqs.is_null());

        // wait for all sends to finish
        send_reqs.wait_all();

        // copy extra local DOFs
        if(this->_extra_local_dof_count > Index(0))
        {
          XASSERT(extra_local != nullptr);
          const DT2_* src = &owned_dofs[this->_extra_dof_offset];
          for(Index i = 0u; i < this->_extra_local_dof_count; ++i)
            extra_local[i] = DataType(src[i]);
        }

        // scatter extra shared DOFs
        if(this->_extra_shared_dof_count > Index(0))
        {
          XASSERT(extra_shared != nullptr);

          // the last process owns the shared DOFs
          if(this->_comm->rank() + 1 == this->_comm->size())
          {
            const DT2_* src = &owned_dofs[this->_extra_dof_offset + this->_extra_local_dof_count];
            for(Index i = 0u; i < this->_extra_shared_dof_count; ++i)
              extra_shared[i] = DataType(src[i]);
          }

          // broadcast shared DOFs
          this->_comm->bcast(extra_shared, std::size_t(this->_extra_shared_dof_count), this->_comm->size() - 1);
        }
      }

    protected:
      static Index _get_donee_dofs(std::vector<IndexType>& send_buf,
        const LAFEM::VectorMirror<DataType, IndexType>& donee_mir, Index offset)
      {
        const Index num_idx = donee_mir.num_indices();
        const IndexType* mir_idx = donee_mir.indices();
        for(Index j(0); j < num_idx; ++j)
          send_buf[j] = IndexType(offset) + mir_idx[j];

        return num_idx;
      }

      template<typename DT2_>
      static Index _gather(
        LAFEM::DenseVector<DataType, IndexType>& buffer,
        const DT2_* vector,
        const LAFEM::VectorMirror<DataType, IndexType>& mirror)
      {
        const Index n = mirror.num_indices();
        const IndexType* idx = mirror.indices();
        DataType* buf = buffer.elements();
        for(Index i = 0u; i < n; ++i)
          buf[i] = DT2_(vector[idx[i]]);
        return n;
      }
    }; // class AlgDofParti<...>

    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    /// \cond internal
    namespace Intern
    {
      template<typename IT_>
      class ADPAux<LAFEM::DenseVector<IT_, IT_>>
      {
      public:
        template<typename DT_>
        static Index get_local_dof_count(std::vector<Index>& v,
          const LAFEM::DenseVector<IT_, IT_>& DOXY(global_dof_idx),
          const LAFEM::VectorMirror<DT_, IT_>& owned_mirror)
        {
          v.push_back(owned_mirror.num_indices());
          return v.back();
        }

        static String build_block_info(
          const LAFEM::DenseVector<IT_, IT_>& DOXY(global_dof_idx), Index global_offset,
          std::deque<Index>& global_first, std::deque<Index>& global_num,
          std::deque<Index>& local_offset, std::deque<Index>& local_num)
        {
          String s = "<Scalar";
          s += " gc=\"" + stringify(global_num.front()) + "\"";
          s += " gf=\"" + stringify(global_first.front()) + "\"";
          s += " go=\"" + stringify(global_offset + local_offset.front()) + "\"";
          s += " lc=\"" + stringify(local_num.front()) + "\"";
          s += " lo=\"" + stringify(local_offset.front()) + "\"";
          s += "/>";
          global_first.pop_front();
          global_num.pop_front();
          local_offset.pop_front();
          local_num.pop_front();
          return s;
        }

        template<typename DT_>
        static void alloc_idx_vector(LAFEM::DenseVector<IT_, IT_>& idx_vec,
          const LAFEM::DenseVector<DT_, IT_>& tmpl_vec, IT_ value)
        {
          idx_vec = LAFEM::DenseVector<IT_, IT_>(tmpl_vec.size(), value);
        }

        template<typename DT_>
        static void update_owners(const IT_ neighbor_rank,
          LAFEM::DenseVector<IT_, IT_>& dof_owners,
          const LAFEM::VectorMirror<DT_, IT_>& mirror)
        {
          const Index num_indices = mirror.num_indices();
          const IT_* mir_idx = mirror.indices();
          IT_* dof_own = dof_owners.elements();
          for(Index j(0); j < num_indices; ++j)
          {
            IT_& dof_owner = dof_own[mir_idx[j]];
            if(neighbor_rank < dof_owner)
              dof_owner = neighbor_rank;
          }
        }

        static Index build_owned_dofs(const IT_ my_rank,
          LAFEM::DenseVector<IT_, IT_>& owned_dofs,
          const LAFEM::DenseVector<IT_, IT_>& dof_owners, Index offset)
        {
          const Index n = dof_owners.size();
          const IT_* dof_own = dof_owners.elements();
          IT_* own_idx = owned_dofs.elements();

          Index k = 0;
          for(Index i = 0; i < n; ++i)
          {
            if(dof_own[i] == my_rank)
              own_idx[i] = IT_(offset + k++);
          }

          return k;
        }

        template<typename DT_>
        static Index build_owned_mirror(const IT_ owner_rank,
          LAFEM::VectorMirror<DT_, IT_>& owner_mirror,
          const LAFEM::DenseVector<IT_, IT_>& dof_owners)
        {
          // count the number of owned DOFs
          const Index n = dof_owners.size();
          const IT_* dof_own = dof_owners.elements();
          Index num_owned = 0u;
          for(Index i = 0; i < n; ++i)
          {
            if(dof_own[i] == owner_rank)
              ++num_owned;
          }

          // allocate mirror
          owner_mirror = LAFEM::VectorMirror<DT_, IT_>(n, num_owned);
          IT_* mir_idx = owner_mirror.indices();

          // store owned DOF indices
          for(Index i = 0, k = 0; i < n; ++i)
          {
            if(dof_own[i] == owner_rank)
            {
              mir_idx[k++] = IT_(i);
            }
          }

          return num_owned;
        }

        template<typename DT_>
        static Index build_owner_mirror(const IT_ neighbor_rank,
          LAFEM::VectorMirror<DT_, IT_>& owner_mirror,
          const LAFEM::VectorMirror<DT_, IT_>& halo_mirror,
          const LAFEM::DenseVector<IT_, IT_>& dof_owners)
        {
          // get halo mirror indices
          const Index n = halo_mirror.num_indices();
          const IT_* halo_idx = halo_mirror.indices();
          const IT_* dof_own = dof_owners.elements();

          // count number of owner dofs
          Index num_owner = 0u;
          for(Index i = 0; i < n; ++i)
          {
            if(dof_own[halo_idx[i]] == neighbor_rank)
              ++num_owner;
          }

          // allocate mirror
          owner_mirror = LAFEM::VectorMirror<DT_, IT_>(dof_owners.size(), num_owner);

          if(num_owner <= Index(0))
            return Index(0);

          // store owner DOF indices
          IT_* own_idx = owner_mirror.indices();
          for(IT_ i = 0, k = 0; i < n; ++i)
          {
            if(dof_own[halo_idx[i]] == neighbor_rank)
            {
              own_idx[k++] = halo_idx[i];
            }
          }

          return num_owner;
        }

        template<typename DT_>
        static Index count_donee_dofs(
          const LAFEM::VectorMirror<DT_, IT_>& halo_mirror,
          const LAFEM::DenseVector<IT_, IT_>& own_dof_idx)
        {
          // get halo mirror indices
          const Index n = halo_mirror.num_indices();
          const IT_* halo_idx = halo_mirror.indices();
          const IT_* own_idx = own_dof_idx.elements();

          // count number of owner dofs
          Index num_donee = 0u;
          for(Index i = 0; i < n; ++i)
          {
            if(own_idx[halo_idx[i]] != ~IT_(0))
              ++num_donee;
          }

          return num_donee;
        }

        template<typename DT_>
        static Index build_donee_mirror(
          LAFEM::VectorMirror<DT_, IT_>& donee_mirror,
          const LAFEM::VectorMirror<DT_, IT_>& halo_mirror,
          const LAFEM::DenseVector<IT_, IT_>& own_dof_idx, Index offset)
        {
          // get halo mirror indices
          const Index n = halo_mirror.num_indices();
          const IT_* halo_idx = halo_mirror.indices();
          const IT_* own_idx = own_dof_idx.elements();
          IT_* donee_idx = donee_mirror.indices();

          // store donee DOF indices
          Index k = 0u;
          for(Index i = 0; i < n; ++i)
          {
            if(own_idx[halo_idx[i]] != ~IT_(0))
            {
              donee_idx[offset + k++] = own_idx[halo_idx[i]];
            }
          }

          return k;
        }

        template<typename DT_>
        static Index init_global_dof_idx(
          LAFEM::DenseVector<IT_, IT_>& global_dof_idx,
          const LAFEM::VectorMirror<DT_, IT_>& owned_mirror, Index offset)
        {
          IT_* g_dof_idx = global_dof_idx.elements();
          const IT_* own_idx = owned_mirror.indices();
          const Index n =  owned_mirror.num_indices();
          for(Index i = 0; i < n; ++i)
            g_dof_idx[own_idx[i]] = IT_(offset + i);
          return n;
        }

        template<typename DT_>
        static Index set_owner_dofs(
          LAFEM::DenseVector<IT_, IT_>& glob_dof_idx,
          const std::vector<IT_>& recv_buf,
          const LAFEM::VectorMirror<DT_, IT_>& mirror, Index offset)
        {
          const Index num_indices = mirror.num_indices();
          const IT_* mir_idx = mirror.indices();
          IT_* g_dof_idx = glob_dof_idx.elements();
          for(Index j(0); j < num_indices; ++j)
            g_dof_idx[mir_idx[j]] = recv_buf[offset + j];
          return num_indices;
        }

        template<typename DT_, typename DT2_>
        static Index gather(DT2_* buf,
          const LAFEM::DenseVector<DT_, IT_>& vector,
          const LAFEM::VectorMirror<DT_, IT_>& mirror)
        {
          XASSERT(mirror.size() == vector.size());
          const Index n = mirror.num_indices();
          const IT_* idx = mirror.indices();
          const DT_* val = vector.elements();
          for(Index i = 0u; i < n; ++i)
            buf[i] = DT2_(val[idx[i]]);
          return n;
        }

        template<typename DT_, typename DT2_>
        static Index scatter(const DT2_* buf,
          LAFEM::DenseVector<DT_, IT_>& vector,
          const LAFEM::VectorMirror<DT_, IT_>& mirror)
        {
          XASSERT(mirror.size() == vector.size());
          const Index n = mirror.num_indices();
          const IT_* idx = mirror.indices();
          DT_* val = vector.elements();
          for(Index i = 0u; i < n; ++i)
            val[idx[i]] = DT_(buf[i]);
          return n;
        }
      }; // class ADPAux<LAFEM::DenseVector<IT_, IT_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename IT_, int bs_>
      class ADPAux<LAFEM::DenseVectorBlocked<IT_, IT_, bs_>>
      {
      public:
        template<typename DT_>
        static Index get_local_dof_count(std::vector<Index>& v,
          const LAFEM::DenseVectorBlocked<IT_, IT_, bs_>& DOXY(global_dof_idx),
          const LAFEM::VectorMirror<DT_, IT_>& owned_mirror)
        {
          v.push_back(owned_mirror.num_indices() * Index(bs_));
          return v.back();
        }

        static String build_block_info(
          const LAFEM::DenseVectorBlocked<IT_, IT_, bs_>& DOXY(global_dof_idx), Index global_offset,
          std::deque<Index>& global_first, std::deque<Index>& global_num,
          std::deque<Index>& local_offset, std::deque<Index>& local_num)
        {
          String s = "<Blocked";
          s += " bs=\"" + stringify(bs_) + "\"";
          s += " gc=\"" + stringify(global_num.front()) + "\"";
          s += " gf=\"" + stringify(global_first.front()) + "\"";
          s += " go=\"" + stringify(global_offset + local_offset.front()) + "\"";
          s += " lc=\"" + stringify(local_num.front()) + "\"";
          s += " lo=\"" + stringify(local_offset.front()) + "\"";
          s += "/>";
          global_first.pop_front();
          global_num.pop_front();
          local_offset.pop_front();
          local_num.pop_front();
          return s;
        }

        template<typename DT_>
        static void alloc_idx_vector(LAFEM::DenseVectorBlocked<IT_, IT_, bs_>& idx_vec,
          const LAFEM::DenseVectorBlocked<DT_, IT_, bs_>& tmpl_vec, IT_ value)
        {
          idx_vec = LAFEM::DenseVectorBlocked<IT_, IT_, bs_>(tmpl_vec.size(), value);
        }

        template<typename DT_>
        static void update_owners(const IT_ neighbor_rank,
          LAFEM::DenseVectorBlocked<IT_, IT_, bs_>& dof_owners,
          const LAFEM::VectorMirror<DT_, IT_>& mirror)
        {
          const Index num_indices = mirror.num_indices();
          const IT_* mir_idx = mirror.indices();
          auto* dof_own = dof_owners.elements(); // dof_own is a Tiny::Vector<IT_,...>*
          for(Index j(0); j < num_indices; ++j)
          {
            auto& dof_owner = dof_own[mir_idx[j]];
            if(neighbor_rank < dof_owner[0])
            {
              for(int k = 0; k < bs_; ++k)
                dof_owner[k] = neighbor_rank;
            }
          }
        }

        static Index build_owned_dofs(const IT_ my_rank,
          LAFEM::DenseVectorBlocked<IT_, IT_, bs_>& owned_dofs,
          const LAFEM::DenseVectorBlocked<IT_, IT_, bs_>& dof_owners, Index offset)
        {
          const Index n = dof_owners.size();
          const auto* dof_own = dof_owners.elements();
          auto* own_idx = owned_dofs.elements();

          Index k = 0;
          for(Index i = 0; i < n; ++i)
          {
            if(dof_own[i][0] == my_rank)
            {
              for(int j = 0; j < bs_; ++j)
                own_idx[i][j] = IT_(offset + k++);
            }
          }

          return k;
        }

        template<typename DT_>
        static Index build_owned_mirror(const IT_ owner_rank,
          LAFEM::VectorMirror<DT_, IT_>& owner_mirror,
          const LAFEM::DenseVectorBlocked<IT_, IT_, bs_>& dof_owners)
        {
          // count the number of owned DOFs
          const Index n = dof_owners.size();
          const auto* dof_own = dof_owners.elements();
          Index num_owned = 0u;
          for(Index i = 0; i < n; ++i)
          {
            if(dof_own[i][0] == owner_rank)
              ++num_owned;
          }

          // allocate mirror
          owner_mirror = LAFEM::VectorMirror<DT_, IT_>(n, num_owned);
          IT_* mir_idx = owner_mirror.indices();

          // store owned DOF indices
          for(Index i = 0, k = 0; i < n; ++i)
          {
            if(dof_own[i][0] == owner_rank)
            {
              mir_idx[k++] = IT_(i);
            }
          }

          return num_owned * Index(bs_);
        }

        template<typename DT_>
        static Index build_owner_mirror(const IT_ neighbor_rank,
          LAFEM::VectorMirror<DT_, IT_>& owner_mirror,
          const LAFEM::VectorMirror<DT_, IT_>& halo_mirror,
          const LAFEM::DenseVectorBlocked<IT_, IT_, bs_>& dof_owners)
        {
          // get halo mirror indices
          const Index n = halo_mirror.num_indices();
          const IT_* halo_idx = halo_mirror.indices();
          const auto* dof_own = dof_owners.elements();

          // count number of owner dofs
          Index num_owner = 0u;
          for(Index i = 0; i < n; ++i)
          {
            if(dof_own[halo_idx[i]][0] == neighbor_rank)
              ++num_owner;
          }

          // allocate mirror
          owner_mirror = LAFEM::VectorMirror<DT_, IT_>(dof_owners.size(), num_owner);

          if(num_owner <= Index(0))
            return Index(0);

          // store owner DOF indices
          IT_* own_idx = owner_mirror.indices();
          for(IT_ i = 0, k = 0; i < n; ++i)
          {
            if(dof_own[halo_idx[i]][0] == neighbor_rank)
            {
              own_idx[k++] = halo_idx[i];
            }
          }

          return num_owner * Index(bs_);
        }

        template<typename DT_>
        static Index count_donee_dofs(
          const LAFEM::VectorMirror<DT_, IT_>& halo_mirror,
          const LAFEM::DenseVectorBlocked<IT_, IT_, bs_>& own_dof_idx)
        {
          // get halo mirror indices
          const Index n = halo_mirror.num_indices();
          const IT_* halo_idx = halo_mirror.indices();
          const auto* own_idx = own_dof_idx.elements();

          // count number of owner dofs
          Index num_donee = 0u;
          for(Index i = 0; i < n; ++i)
          {
            if(own_idx[halo_idx[i]][0] != ~IT_(0))
              num_donee += Index(bs_);
          }

          return num_donee;
        }

        template<typename DT_>
        static Index build_donee_mirror(
          LAFEM::VectorMirror<DT_, IT_>& donee_mirror,
          const LAFEM::VectorMirror<DT_, IT_>& halo_mirror,
          const LAFEM::DenseVectorBlocked<IT_, IT_, bs_>& own_dof_idx, Index offset)
        {
          // get halo mirror indices
          const Index n = halo_mirror.num_indices();
          const IT_* halo_idx = halo_mirror.indices();
          const auto* own_idx = own_dof_idx.elements();
          IT_* donee_idx = donee_mirror.indices();

          // store donee DOF indices
          Index k = 0u;
          for(Index i = 0; i < n; ++i)
          {
            if(own_idx[halo_idx[i]][0] != ~IT_(0))
            {
              for(int j = 0; j < bs_; ++j)
                donee_idx[offset + k++] = own_idx[halo_idx[i]][j];
            }
          }

          return k;
        }

        template<typename DT_>
        static Index init_global_dof_idx(
          LAFEM::DenseVectorBlocked<IT_, IT_, bs_>& global_dof_idx,
          const LAFEM::VectorMirror<DT_, IT_>& owned_mirror, Index offset)
        {
          auto* g_dof_idx = global_dof_idx.elements();
          const IT_* own_idx = owned_mirror.indices();
          const Index n =  owned_mirror.num_indices();
          for(Index i = 0; i < n; ++i)
          {
            for(int j = 0; j < bs_; ++j)
              g_dof_idx[own_idx[i]][j] = IT_(offset + i*Index(bs_) + Index(j));
          }
          return n * Index(bs_);
        }

        template<typename DT_>
        static Index set_owner_dofs(
          LAFEM::DenseVectorBlocked<IT_, IT_, bs_>& glob_dof_idx,
          const std::vector<IT_>& recv_buf,
          const LAFEM::VectorMirror<DT_, IT_>& mirror, Index offset)
        {
          const Index num_indices = mirror.num_indices();
          const IT_* mir_idx = mirror.indices();
          auto* g_dof_idx = glob_dof_idx.elements();
          for(Index i(0); i < num_indices; ++i)
          {
            for(int j = 0; j < bs_; ++j)
              g_dof_idx[mir_idx[i]][j] = recv_buf[offset + i*IT_(bs_) + IT_(j)];
          }
          return num_indices * Index(bs_);
        }

        template<typename DT_, typename DT2_>
        static Index gather(DT2_* buf,
          const LAFEM::DenseVectorBlocked<DT_, IT_, bs_>& vector,
          const LAFEM::VectorMirror<DT_, IT_>& mirror)
        {
          XASSERT(mirror.size() == vector.size());
          const Index n = mirror.num_indices();
          const IT_* idx = mirror.indices();
          const auto* val = vector.elements();
          for(Index i = 0u, k = 0u; i < n; ++i)
          {
            for(int j = 0; j < bs_; ++j, ++k)
              buf[k] = DT2_(val[idx[i]][j]);
          }
          return n * Index(bs_);
        }

        template<typename DT_, typename DT2_>
        static Index scatter(const DT2_* buf,
          LAFEM::DenseVectorBlocked<DT_, IT_, bs_>& vector,
          const LAFEM::VectorMirror<DT_, IT_>& mirror)
        {
          XASSERT(mirror.size() == vector.size());
          const Index n = mirror.num_indices();
          const IT_* idx = mirror.indices();
          auto* val = vector.elements();
          for(Index i = 0u, k = 0u; i < n; ++i)
          {
            for(int j = 0; j < bs_; ++j, ++k)
              val[idx[i]][j] = DT_(buf[k]);
          }
          return n * Index(bs_);
        }
      }; // class ADPAux<LAFEM::DenseVectorBlocked<IT_, IT_, bs_>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstIdxVec_, typename... RestIdxVec_>
      class ADPAux<LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>>
      {
      public:
        typedef typename FirstIdxVec_::IndexType IT_;
        typedef ADPAux<FirstIdxVec_> ADPAuxFirst;
        typedef ADPAux<LAFEM::TupleVector<RestIdxVec_...>> ADPAuxRest;

        template<typename FirstMirror_, typename... RestMirror_>
        static Index get_local_dof_count_raw(std::vector<Index>& v,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& owned_mirror)
        {
          Index n1 = ADPAuxFirst::get_local_dof_count(v, global_dof_idx.first(), owned_mirror.first());
          Index n2 = ADPAuxRest::get_local_dof_count_raw(v, global_dof_idx.rest(), owned_mirror.rest());
          return n1 + n2;
        }

        template<
          typename FirstMirror_, typename... RestMirror_>
        static Index get_local_dof_count(std::vector<Index>& v,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& owned_mirror)
        {
          v.push_back(get_local_dof_count_raw(v, global_dof_idx, owned_mirror));
          return v.back();
        }

        static String build_block_info_raw(
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx, Index global_offset,
          std::deque<Index>& global_first, std::deque<Index>& global_num,
          std::deque<Index>& local_offset, std::deque<Index>& local_num)
        {
          String s = ADPAuxFirst::build_block_info(global_dof_idx.first(), global_offset, global_first, global_num, local_offset, local_num);
          s += "\n";
          s += ADPAuxRest::build_block_info_raw(global_dof_idx.rest(), global_offset, global_first, global_num, local_offset, local_num);
          return s;
        }

        static String build_block_info(
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx, Index global_offset,
          std::deque<Index>& global_first, std::deque<Index>& global_num,
          std::deque<Index>& local_offset, std::deque<Index>& local_num)
        {
          // offsets need to be processed BEFORE the call to '_build_block_info_raw' because
          // the offsets of the tuple coincide with the offsets of the first tuple entry
          String slo = stringify(local_offset.front());
          String sgo = stringify(global_offset + local_offset.front());
          String sraw = build_block_info_raw(global_dof_idx, global_offset, global_first, global_num, local_offset, local_num);
          String s = "<Tuple";
          // counts and first need to be processed AFTER the call to '_build_block_info_raw' because
          // these correspond to the sums of the corresponding values of all the tuple entries
          s += " gc=\"" + stringify(global_num.front()) + "\"";
          s += " gf=\"" + stringify(global_first.front()) + "\"";
          s += " go=\"" + sgo + "\"";
          s += " lc=\"" + stringify(local_num.front()) + "\"";
          s += " lo=\"" + slo + "\">\n";
          s += sraw;
          s += "\n</Tuple>";
          global_first.pop_front();
          global_num.pop_front();
          local_offset.pop_front();
          local_num.pop_front();
          return s;
        }

        template<typename FirstTmplVector_, typename... RestTmplVector_>
        static void alloc_idx_vector(LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& idx_vec,
          const LAFEM::TupleVector<FirstTmplVector_, RestTmplVector_...>& tmpl_vec, IT_ value)
        {
          ADPAuxFirst::alloc_idx_vector(idx_vec.first(), tmpl_vec.first(), value);
          ADPAuxRest::alloc_idx_vector(idx_vec.rest(), tmpl_vec.rest(), value);
        }

        template<typename FirstMirror_, typename... RestMirror_>
        static void update_owners(const IT_ neighbor_rank,
          LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& dof_owners,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& mirror)
        {
          ADPAuxFirst::update_owners(neighbor_rank, dof_owners.first(), mirror.first());
          ADPAuxRest::update_owners(neighbor_rank, dof_owners.rest(), mirror.rest());
        }

        static Index build_owned_dofs(const IT_ my_rank,
          LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& owned_dofs,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& dof_owners, Index offset)
        {
          Index n1 = ADPAuxFirst::build_owned_dofs(my_rank, owned_dofs.first(), dof_owners.first(), offset);
          Index n2 = ADPAuxRest::build_owned_dofs(my_rank, owned_dofs.rest(), dof_owners.rest(), offset + n1);
          return n1 + n2;
        }

        template<typename FirstMirror_, typename... RestMirror_>
        static Index build_owned_mirror(const IT_ owner_rank,
          LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& owner_mirror,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& dof_owners)
        {
          Index n1 = ADPAuxFirst::build_owned_mirror(owner_rank, owner_mirror.first(), dof_owners.first());
          Index n2 = ADPAuxRest::build_owned_mirror(owner_rank, owner_mirror.rest(), dof_owners.rest());
          return n1 + n2;
        }

        template<typename FirstMirror_, typename... RestMirror_>
        static Index build_owner_mirror(const IT_ neighbor_rank,
          LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& owner_mirror,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& halo_mirror,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& dof_owners)
        {
          Index n1 = ADPAuxFirst::build_owner_mirror(neighbor_rank, owner_mirror.first(), halo_mirror.first(), dof_owners.first());
          Index n2 = ADPAuxRest::build_owner_mirror(neighbor_rank, owner_mirror.rest(), halo_mirror.rest(), dof_owners.rest());
          return n1 + n2;
        }

        template<typename FirstMirror_, typename... RestMirror_>
        static Index count_donee_dofs(
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& halo_mirror,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& own_dof_idx)
        {
          Index n1 = ADPAuxFirst::count_donee_dofs(halo_mirror.first(), own_dof_idx.first());
          Index n2 = ADPAuxRest::count_donee_dofs(halo_mirror.rest(), own_dof_idx.rest());
          return n1 + n2;
        }

        template<typename DT_, typename FirstMirror_, typename... RestMirror_>
        static Index build_donee_mirror(
          LAFEM::VectorMirror<DT_, IT_>& donee_mirror,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& halo_mirror,
          const LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& own_dof_idx, Index offset)
        {
          Index n1 = ADPAuxFirst::build_donee_mirror(donee_mirror, halo_mirror.first(), own_dof_idx.first(), offset);
          Index n2 = ADPAuxRest::build_donee_mirror(donee_mirror, halo_mirror.rest(), own_dof_idx.rest(), offset + n1);
          return n1 + n2;
        }

        template<typename FirstMirror_, typename... RestMirror_>
        static Index init_global_dof_idx(
          LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& global_dof_idx,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& owned_mirror, Index offset)
        {
          Index n1 = ADPAuxFirst::init_global_dof_idx(global_dof_idx.first(), owned_mirror.first(), offset);
          Index n2 = ADPAuxRest::init_global_dof_idx(global_dof_idx.rest(), owned_mirror.rest(), offset + n1);
          return n1 + n2;
        }

        template<typename FirstMirror_, typename... RestMirror_>
        static Index set_owner_dofs(
          LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>& glob_dof_idx,
          const std::vector<IT_>& recv_buf,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& mirror, Index offset)
        {
          Index n1 = ADPAuxFirst::set_owner_dofs(glob_dof_idx.first(), recv_buf, mirror.first(), offset);
          Index n2 = ADPAuxRest::set_owner_dofs(glob_dof_idx.rest(), recv_buf, mirror.rest(), offset + n1);
          return n1 + n2;
        }

        template<typename DT2_,
          typename FirstVector_, typename... RestVector_,
          typename FirstMirror_, typename... RestMirror_>
        static Index gather(DT2_* buf,
          const LAFEM::TupleVector<FirstVector_, RestVector_...>& vector,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& mirror)
        {
          Index n1 = ADPAuxFirst::gather(buf, vector.first(), mirror.first());
          Index n2 = ADPAuxRest::gather(&buf[n1], vector.rest(), mirror.rest());
          return n1 + n2;
        }

        template<typename DT2_,
          typename FirstVector_, typename... RestVector_,
          typename FirstMirror_, typename... RestMirror_>
        static Index scatter(const DT2_* buf,
          LAFEM::TupleVector<FirstVector_, RestVector_...>& vector,
          const LAFEM::TupleMirror<FirstMirror_, RestMirror_...>& mirror)
        {
          Index n1 = ADPAuxFirst::scatter(buf, vector.first(), mirror.first());
          Index n2 = ADPAuxRest::scatter(&buf[n1], vector.rest(), mirror.rest());
          return n1 + n2;
        }
      }; // class ADPAux<LAFEM::TupleVector<FirstIdxVec_, RestIdxVec_...>>

      /////////////////////////////////////////////////////////////////////////////////////////////////////////////////

      template<typename FirstIdxVec_>
      class ADPAux<LAFEM::TupleVector<FirstIdxVec_>>
      {
      public:
        typedef typename FirstIdxVec_::IndexType IT_;
        typedef ADPAux<FirstIdxVec_> ADPAuxFirst;

        template<typename FirstMirror_>
        static Index get_local_dof_count_raw(std::vector<Index>& v,
          const LAFEM::TupleVector<FirstIdxVec_>& global_dof_idx,
          const LAFEM::TupleMirror<FirstMirror_>& owned_mirror)
        {
          return ADPAuxFirst::get_local_dof_count(v, global_dof_idx.first(), owned_mirror.first());
        }

        template<typename FirstMirror_>
        static Index get_local_dof_count(std::vector<Index>& v,
          const LAFEM::TupleVector<FirstIdxVec_>& global_dof_idx,
          const LAFEM::TupleMirror<FirstMirror_>& owned_mirror)
        {
          v.push_back(get_local_dof_count_raw(v, global_dof_idx, owned_mirror));
          return v.back();
        }

        static String build_block_info_raw(
          const LAFEM::TupleVector<FirstIdxVec_>& global_dof_idx, Index global_offset,
          std::deque<Index>& global_first, std::deque<Index>& global_num,
          std::deque<Index>& local_offset, std::deque<Index>& local_num)
        {
          return ADPAuxFirst::build_block_info(global_dof_idx.first(), global_offset, global_first, global_num, local_offset, local_num);
        }

        static String build_block_info(
          const LAFEM::TupleVector<FirstIdxVec_>& global_dof_idx, Index global_offset,
          std::deque<Index>& global_first, std::deque<Index>& global_num,
          std::deque<Index>& local_offset, std::deque<Index>& local_num)
        {
          // offsets need to be processed BEFORE the call to '_build_block_info_raw' because
          // the offsets of the tuple coincide with the offsets of the first tuple entry
          String slo = stringify(local_offset.front());
          String sgo = stringify(global_offset + local_offset.front());
          String sraw = build_block_info_raw(global_dof_idx, global_offset, global_first, global_num, local_offset, local_num);
          String s = "<Tuple";
          // counts and first need to be processed AFTER the call to '_build_block_info_raw' because
          // these correspond to the sums of the corresponding values of all the tuple entries
          s += " gc=\"" + stringify(global_num.front()) + "\"";
          s += " gf=\"" + stringify(global_first.front()) + "\"";
          s += " go=\"" + sgo + "\"";
          s += " lc=\"" + stringify(local_num.front()) + "\"";
          s += " lo=\"" + slo + "\">\n";
          s += sraw;
          s += "\n</Tuple>";
          global_first.pop_front();
          global_num.pop_front();
          local_offset.pop_front();
          local_num.pop_front();
          return s;
        }

        template<typename FirstTmplVector_>
        static void alloc_idx_vector(LAFEM::TupleVector<FirstIdxVec_>& idx_vec,
          const LAFEM::TupleVector<FirstTmplVector_>& tmpl_vec, IT_ value)
        {
          ADPAuxFirst::alloc_idx_vector(idx_vec.first(), tmpl_vec.first(), value);
        }

        template<typename FirstMirror_>
        static void update_owners(const IT_ neighbor_rank,
          LAFEM::TupleVector<FirstIdxVec_>& dof_owners,
          const LAFEM::TupleMirror<FirstMirror_>& mirror)
        {
          ADPAuxFirst::update_owners(neighbor_rank, dof_owners.first(), mirror.first());
        }

        static Index build_owned_dofs(const IT_ my_rank,
          LAFEM::TupleVector<FirstIdxVec_>& owned_dofs,
          const LAFEM::TupleVector<FirstIdxVec_>& dof_owners, Index offset)
        {
          return ADPAuxFirst::build_owned_dofs(my_rank, owned_dofs.first(), dof_owners.first(), offset);
        }

        template<typename FirstMirror_>
        static Index build_owned_mirror(const IT_ owner_rank,
          LAFEM::TupleMirror<FirstMirror_>& owner_mirror,
          const LAFEM::TupleVector<FirstIdxVec_>& dof_owners)
        {
          return ADPAuxFirst::build_owned_mirror(owner_rank, owner_mirror.first(), dof_owners.first());
        }

        template<typename FirstMirror_>
        static Index build_owner_mirror(const IT_ neighbor_rank,
          LAFEM::TupleMirror<FirstMirror_>& owner_mirror,
          const LAFEM::TupleMirror<FirstMirror_>& halo_mirror,
          const LAFEM::TupleVector<FirstIdxVec_>& dof_owners)
        {
          return ADPAuxFirst::build_owner_mirror(neighbor_rank, owner_mirror.first(), halo_mirror.first(), dof_owners.first());
        }

        template<typename FirstMirror_>
        static Index count_donee_dofs(
          const LAFEM::TupleMirror<FirstMirror_>& halo_mirror,
          const LAFEM::TupleVector<FirstIdxVec_>& own_dof_idx)
        {
          return ADPAuxFirst::count_donee_dofs(halo_mirror.first(), own_dof_idx.first());
        }

        template<typename DT_, typename FirstMirror_>
        static Index build_donee_mirror(
          LAFEM::VectorMirror<DT_, IT_>& donee_mirror,
          const LAFEM::TupleMirror<FirstMirror_>& halo_mirror,
          const LAFEM::TupleVector<FirstIdxVec_>& own_dof_idx, Index offset)
        {
          return ADPAuxFirst::build_donee_mirror(donee_mirror, halo_mirror.first(), own_dof_idx.first(), offset);
        }

        template<typename FirstMirror_>
        static Index init_global_dof_idx(
          LAFEM::TupleVector<FirstIdxVec_>& global_dof_idx,
          const LAFEM::TupleMirror<FirstMirror_>& owned_mirror, Index offset)
        {
          return ADPAuxFirst::init_global_dof_idx(global_dof_idx.first(), owned_mirror.first(), offset);
        }

        template<typename FirstMirror_>
        static Index set_owner_dofs(
          LAFEM::TupleVector<FirstIdxVec_>& glob_dof_idx,
          const std::vector<IT_>& recv_buf,
          const LAFEM::TupleMirror<FirstMirror_>& mirror, Index offset)
        {
          return ADPAuxFirst::set_owner_dofs(glob_dof_idx.first(), recv_buf, mirror.first(), offset);
        }

        template<typename DT2_, typename FirstVector_, typename FirstMirror_>
        static Index gather(DT2_* buf,
          const LAFEM::TupleVector<FirstVector_>& vector,
          const LAFEM::TupleMirror<FirstMirror_>& mirror)
        {
          return ADPAuxFirst::gather(buf, vector.first(), mirror.first());
        }

        template<typename DT2_, typename FirstVector_, typename FirstMirror_>
        static Index scatter(const DT2_* buf,
          LAFEM::TupleVector<FirstVector_>& vector,
          const LAFEM::TupleMirror<FirstMirror_>& mirror)
        {
          return ADPAuxFirst::scatter(buf, vector.first(), mirror.first());
        }
      }; // class ADPAux<LAFEM::TupleVector<FirstIdxVec_>>
    } // namespace Intern
    /// \endcond
  } // namespace Global
} // namespace FEAT
