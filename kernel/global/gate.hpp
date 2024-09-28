// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GLOBAL_GATE_HPP
#define KERNEL_GLOBAL_GATE_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/global/synch_vec.hpp>
#include <kernel/global/synch_scal.hpp>

#include <vector>

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
     * This class provides the functionality for the synchronization of data across interfaces
     * between nearest neighbors (aka 'halos') in the overlapping domain decomposition approach
     * and is used by most other Global classes. Note that a gate is tied to one finite element
     * space (or a tuple thereof) and one vector type, i.e. if you need to synchronize different
     * vector types belonging to possibly different finite element spaces, then you also need one
     * gate object for each space-vector pair. In the case where you need two gate objects for the
     * same finite element space object, but for different vector classes, you can make use of the
     * convert functions to convert one gate into another one for the second vector type.
     *
     * \tparam LocalVector_
     * The type of the local vector container; may be any valid combination of LAFEM (meta-)vector types
     *
     * \tparam Mirror_
     * The type of the vector mirror; must be compatible to the local vector type
     *
     * To set up an object of this class, one has to perform three separate steps:
     * -# Set the Dist::Comm object for the gate by using the set_comm() function.
     * -# Add the mirror for each nearest neighbor by using the push() function.
     * -# Compile the gate by calling the compile() function and supplying it with a temporary local vector.
     *
     * \author Peter Zajac
     */
    template<typename LocalVector_, typename Mirror_>
    class Gate
    {
    public:
      /// the data type
      typedef typename LocalVector_::DataType DataType;
      /// the index type
      typedef typename LocalVector_::IndexType IndexType;
      /// the local vector type
      typedef LocalVector_ LocalVectorType;
      /// the mirror type
      typedef Mirror_ MirrorType;
      /// the internal buffer vector type
      typedef LAFEM::DenseVector<DataType, IndexType> BufferVectorType;

      typedef SynchScalarTicket<DataType> ScalarTicketType;
      typedef SynchVectorTicket<LocalVector_, Mirror_> VectorTicketType;

    public:
      /// our communicator
      const Dist::Comm* _comm;
      /// communication ranks
      std::vector<int> _ranks;
      /// vector mirrors
      std::vector<Mirror_> _mirrors;
      /// frequency vector
      LocalVector_ _freqs;

      /// Our 'base' class type
      template <typename LocalVector2_, typename Mirror2_>
      using GateType = Gate<LocalVector2_, Mirror2_>;

      /// this typedef lets you create a gate container with new Data and Index types
      template <typename DataType2_, typename IndexType2_>
      using GateTypeByDI = Gate<typename LocalVector_::template ContainerType<DataType2_, IndexType2_>, typename Mirror_::template MirrorType<DataType2_, IndexType2_> >;

    public:
      /// standard constructor
      explicit Gate() :
        _comm(nullptr)
      {
      }

      /**
       * \brief Constructor
       *
       * \param[in] comm
       * A \resident reference to the communicator to be used by the gate.
       */
      explicit Gate(const Dist::Comm& comm) :
        _comm(&comm)
      {
      }

      /// move constructor
      Gate(Gate&& other) :
        _comm(other._comm),
        _ranks(std::forward<std::vector<int>>(other._ranks)),
        _mirrors(std::forward<std::vector<Mirror_>>(other._mirrors)),
        _freqs(std::forward<LocalVector_>(other._freqs))
      {
      }

      /// move-assign operator
      Gate& operator=(Gate&& other)
      {
        if(this == &other)
        {
          return *this;
        }

        _comm = other._comm;
        _ranks = std::forward<std::vector<int>>(other._ranks);
        _mirrors = std::forward<std::vector<Mirror_>>(other._mirrors);
        _freqs = std::forward<LocalVector_>(other._freqs);

        return *this;
      }

      /// virtual destructor
      virtual ~Gate()
      {
      }

      /**
       * \brief Returns a const pointer to the underlying communicator
       *
       * \returns A const pointer to the underlying communicator
       */
      const Dist::Comm* get_comm() const
      {
        return this->_comm;
      }

      /**
       * \brief Sets the communicator for this gate
       *
       * \param[in] comm_
       * A \resident pointer to the communicator object to use.
       */
      void set_comm(const Dist::Comm* comm_)
      {
        this->_comm = comm_;
      }

      /**
       * \brief Returns a const reference to the neighbor ranks vector
       *
       * \returns A const reference to the neighbor ranks vector
       */
      const std::vector<int> get_ranks() const
      {
        return this->_ranks;
      }

      /**
       * \brief Returns a const reference to the neighbor mirrors vector
       *
       * \returns A const reference to the neighbor mirrors vector
       */
      const std::vector<Mirror_>& get_mirrors() const
      {
        return this->_mirrors;
      }

      /**
       * \brief Returns a const reference to the frequencies vector
       *
       * \returns A const reference to the frequencies vector
       */
      const LocalVector_& get_freqs() const
      {
        return this->_freqs;
      }

      /**
       * \brief Conversion function for same vector container type but with different MDI-Type
       *
       * \param[in] other
       * A \transient reference to the gate to convert from
       */
      template<typename LVT2_, typename MT2_>
      void convert(const Gate<LVT2_, MT2_>& other)
      {
        if((void*)this == (void*)&other)
          return;

        this->_comm = other._comm;
        this->_ranks = other._ranks;

        this->_mirrors.resize(other._mirrors.size());
        for(std::size_t i(0); i < other._mirrors.size(); ++i)
        {
          this->_mirrors.at(i).convert(other._mirrors.at(i));
        }

        this->_freqs.convert(other._freqs);
      }

      /**
       * \brief Conversion function for different vector container type
       *
       * This function (re)creates this gate from another gate with a different vector type, but the same mirrors.
       * This can be used to create a gate using DenseVectorBlocked from a gate using DenseVector or vice versa.
       *
       * \param[in] other
       * A (transient) reference to the gate using the other vector type to create this gate from.
       *
       * \param[in] vector
       * A temporary vector allocated to the correct size that is to be used for internal initialization.
       *
       * \param[in] mode
       * The clone-mode to be used for cloning the mirrors. Defaults to shallow clone.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       *
       */
      template<typename LVT2_>
      void convert(const Gate<LVT2_, Mirror_>& other, LocalVector_&& vector, LAFEM::CloneMode mode = LAFEM::CloneMode::Shallow)
      {
        if((void*)this == (void*)&other)
          return;

        this->_ranks.clear();
        this->_mirrors.clear();
        this->_mirrors.resize(other._mirrors.size());
        this->_freqs.clear(); // will be rebuild by compile function

        this->_comm = other._comm;
        this->_ranks = other._ranks;

        // shallow-clone mirrors
        for(std::size_t i(0); i < _mirrors.size(); ++i)
        {
          _mirrors.at(i).clone(other._mirrors.at(i), mode);
        }

        // compile this gate
        compile(std::forward<LocalVector_>(vector));
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
        temp += _ranks.size() * sizeof(int);

        return temp;
      }

      /**
       * \brief Adds a mirror for a neighbor process
       *
       * \param[in] rank
       * The rank of the neighbor with respect to our communicator
       *
       * \param[in] mirror
       * The mirror assembled for the halo of this neighbor process.
       */
      void push(int rank, Mirror_&& mirror)
      {
        XASSERT(this->_comm != nullptr);
        XASSERT(rank < this->_comm->size());

        // push rank and tags
        _ranks.push_back(rank);

        // push mirror
        _mirrors.push_back(std::move(mirror));
      }

      /**
       * \brief Compiles the gate to finish its setup
       *
       * \param[in] vector
       * A temporary vector allocated to the correct size which is used for initialization
       * of the internal frequencies vector. Its numerical contents are ignored.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      void compile(LocalVector_&& vector)
      {
        // initialize frequency vector
        _freqs = std::move(vector);
        _freqs.format(DataType(1));

        // loop over all mirrors
        for(std::size_t i(0); i < _mirrors.size(); ++i)
        {
          // sum up number of ranks per frequency entry, listed in different mirrors
          auto temp = _mirrors.at(i).create_buffer(_freqs);
          temp.format(DataType(1));

          // gather-axpy into frequency vector
          _mirrors.at(i).scatter_axpy(_freqs, temp);
        }

        // invert frequencies
        _freqs.component_invert(_freqs);
      }

      /**
       * \brief Returns the number of local DOFs
       *
       * \tparam perspective_
       * Specifies whether to compute the native or POD DOF count
       *
       * \returns The number of local DOFS; either in native or POD size
       */
      template<LAFEM::Perspective perspective_ = LAFEM::Perspective::native>
      Index get_num_local_dofs() const
      {
        return _freqs.template size<perspective_>();
      }

      /**
      * \brief Returns the number of global DOFs
      *
      * \tparam perspective_
      * Specifies whether to compute the native or POD DOF count
      *
      * \attention This function is collective, i.e. it must be called by all processes participating
      * in the gate's communicator, otherwise the application will deadlock.
      *
      * \returns The number of global DOFS; either in native or POD size
      */
      template<LAFEM::Perspective perspective_ = LAFEM::Perspective::native>
      Index get_num_global_dofs() const
      {
        XASSERT(this->_comm != nullptr);
        if(this->_comm->size() <= 1)
          return this->template get_num_local_dofs<perspective_>();

        // get my rank
        const int my_rank = this->_comm->rank();

        // get local number of DOFs
        const Index loc_dofs = this->template get_num_local_dofs<perspective_>();

        // create a local mask vector and format it to 1
        std::vector<int> mask(std::size_t(loc_dofs), 1);

        // set all DOFs, which are shared with a lower rank neighbor, to 0
        for(std::size_t i(0); i < _mirrors.size(); ++i)
        {
          if(this->_ranks.at(i) < my_rank)
            this->_mirrors.at(i).template mask_scatter<perspective_>(this->_freqs, mask, 0);
        }

        // count the number of DOFs that are still 1 and thus not owned by a lower rank neighbor
        Index owned_dofs(0u);
        for(const auto& k : mask)
          owned_dofs += Index(k);

        // now sum up the number of owned DOFs over all processes
        Index global_dofs(0u);
        this->_comm->allreduce(&owned_dofs, &global_dofs, std::size_t(1), Dist::op_sum);

        // done!
        return global_dofs;
      }

      /**
       * \brief Converts a type-1 vector into a type-0 vector.
       *
       * \param[inout] vector
       * On entry, the type-1 vector to be converted.
       * On exit, the converted type-0 vector.
       *
       * \note This function does not perform any synchronization.
       */
      void from_1_to_0(LocalVector_& vector) const
      {
        if(!_ranks.empty())
        {
          vector.component_product(vector, _freqs);
        }
      }

      /**
       * \brief Synchronizes a type-0 vector, resulting in a type-1 vector.
       *
       * This function performs a type-0 synchronization, i.e. sums up all local DOF contributions
       * for each DOF, by exchanging the DOF values of each shared DOF among all nearest neighbors.
       *
       * \param[inout] vector
       * On entry, the type-0 vector to be synchronized.\n
       * On exit, the synchronized type-1 vector.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      void sync_0(LocalVector_& vector) const
      {
        if(!_ranks.empty())
        {
          SynchVectorTicket<LocalVector_, Mirror_> ticket(vector, *_comm, _ranks, _mirrors);
          ticket.wait();
        }
      }

      /**
       * \brief Synchronizes a type-0 vector, resulting in a type-1 vector.
       *
       * This function performs a type-0 synchronization, i.e. sums up all local DOF contributions
       * for each DOF, by exchanging the DOF values of each shared DOF among all nearest neighbors.
       *
       * \param[inout] vector
       * On entry, the type-0 vector to be synchronized.\n
       * On exit, the synchronized type-1 vector.
       *
       * \returns A ticket that has to be waited upon to complete the sync operation.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      VectorTicketType sync_0_async(LocalVector_& vector) const
      {
        if(_ranks.empty())
          return SynchVectorTicket<LocalVector_, Mirror_>(); // empty ticket

        return SynchVectorTicket<LocalVector_, Mirror_>(vector, *_comm, _ranks, _mirrors);
      }

      /**
       * \brief Synchronizes a type-1 vector, resulting in a type-1 vector.
       *
       * This function performs a type-1 synchronization, i.e. averages all local DOF contributions
       * for each DOF, by exchanging the DOF values of each shared DOF among all nearest neighbors.
       *
       * \param[inout] vector
       * On entry, the type-1 vector to be synchronized.\n
       * On exit, the synchronized type-1 vector.
       *
       * \note
       * This function effectively applies the from_1_to_0() and sync_0()
       * functions onto the input vector.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      void sync_1(LocalVector_& vector) const
      {
        if(!_ranks.empty())
        {
          from_1_to_0(vector);
          SynchVectorTicket<LocalVector_, Mirror_> ticket(vector, *_comm, _ranks, _mirrors);
          ticket.wait();
        }
      }

      /**
       * \brief Synchronizes a type-1 vector, resulting in a type-1 vector.
       *
       * This function performs a type-1 synchronization, i.e. averages all local DOF contributions
       * for each DOF, by exchanging the DOF values of each shared DOF among all nearest neighbors.
       *
       * \param[inout] vector
       * On entry, the type-1 vector to be synchronized.\n
       * On exit, the synchronized type-1 vector.
       *
       * \returns A ticket that has to be waited upon to complete the sync operation.
       *
       * \note
       * This function effectively applies the from_1_to_0() and sync_0() functions onto the input vector.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      VectorTicketType sync_1_async(LocalVector_& vector) const
      {
        if(_ranks.empty())
          return SynchVectorTicket<LocalVector_, Mirror_>(); // empty ticket

        from_1_to_0(vector);
        return SynchVectorTicket<LocalVector_, Mirror_>(vector, *_comm, _ranks, _mirrors);
      }

      /**
       * \brief Computes a synchronized dot-product of two type-1 vectors.
       *
       * \param[in] x, y
       * The two type-1 vector whose dot-product is to be computed.
       *
       * \returns
       * The dot-product of \p x and \p y.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      DataType dot(const LocalVector_& x, const LocalVector_& y) const
      {
        // This is if there is only one process
        if(_comm == nullptr || _comm->size() == 1)
        {
          return x.dot(y);
        }
        // Even if there are no neighbors, we still need to sum up globally
        else if(_ranks.empty())
        {
          return sum(x.dot(y));
        }
        // If there are neighbors, we have to use the frequencies and sum up globally
        else
        {
          return sum(_freqs.triple_dot(x, y));
        }
      }

      /**
       * \brief Computes a synchronized dot-product of two type-1 vectors.
       *
       * \param[in] x, y
       * The two type-1 vector whose dot-product is to be computed.
       *
       * \returns A scalar ticket that has to be waited upon to complete the operation.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      ScalarTicketType dot_async(const LocalVector_& x, const LocalVector_& y, bool sqrt = false) const
      {
        return sum_async(_freqs.triple_dot(x, y), sqrt);
      }

      /**
       * \brief Computes a reduced sum over all processes.
       *
       * \param[in] x
       * The value that is to be summed over all processes.
       *
       * \returns
       * The reduced sum of all \p x.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      DataType sum(DataType x) const
      {
        return sum_async(x).wait();
      }

      /**
       * \brief Computes a reduced sum over all processes.
       *
       * \param[in] x
       * The value that is to be summed over all processes.
       *
       * \returns A scalar ticket that has to be waited upon to complete the operation.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      ScalarTicketType sum_async(DataType x, bool sqrt = false) const
      {
        return SynchScalarTicket<DataType>(x, *_comm, Dist::op_sum, sqrt);
      }

      /**
       * \brief Computes the minimum of a scalar variable over all processes.
       *
       * \param[in] x
       * What we want to compute the minimum over all processes of.
       *
       * \returns
       * The minimum of all \p x.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      DataType min(DataType x) const
      {
        return min_async(x).wait();
      }

      /**
       * \brief Computes the minimum of a scalar variable over all processes.
       *
       * \param[in] x
       * What we want to compute the minimum over all processes of.
       *
       * \returns A scalar ticket that has to be waited upon to complete the operation.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      ScalarTicketType min_async(DataType x) const
      {
        return SynchScalarTicket<DataType>(x, *_comm, Dist::op_min, false);
      }

      /**
       * \brief Computes the maximum of a scalar variable over all processes.
       *
       * \param[in] x
       * What we want to compute the maximum over all processes of.
       *
       * \returns
       * The maximum of all \p x.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      DataType max(DataType x) const
      {
        return max_async(x).wait();
      }

      /**
       * \brief Computes the maximum of a scalar variable over all processes.
       *
       * \param[in] x
       * What we want to compute the maximum over all processes of.
       *
       * \returns A scalar ticket that has to be waited upon to complete the operation.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      ScalarTicketType max_async(DataType x) const
      {
        return SynchScalarTicket<DataType>(x, *_comm, Dist::op_max, false);
      }

      /**
       * \brief Computes a reduced 2-norm over all processes.
       *
       * This function is equivalent to the call
       *    Math::sqrt(this->sum(x*x))
       *
       * \param[in] x
       * The value that is to be summarized over all processes.
       *
       * \returns
       * The reduced 2-norm of all \p x.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      DataType norm2(DataType x) const
      {
        return norm2_async(x).wait();
      }

      /**
       * \brief Computes a reduced 2-norm over all processes.
       *
       * This function is equivalent to the call
       *    Math::sqrt(this->sum(x*x))
       *
       * \param[in] x
       * The value that is to be summarized over all processes.
       *
       * \returns A scalar ticket that has to be waited upon to complete the operation.
       *
       * \attention This function is collective, i.e. it must be called by all processes participating
       * in the gate's communicator, otherwise the application will deadlock.
       */
      ScalarTicketType norm2_async(DataType x) const
      {
        return SynchScalarTicket<DataType>(x*x, *_comm, Dist::op_sum, true);
      }
    }; // class Gate<...>
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_GATE_HPP
