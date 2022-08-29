// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_GLOBAL_VECTOR_HPP
#define KERNEL_GLOBAL_VECTOR_HPP 1

#include <kernel/global/gate.hpp>
#include <kernel/util/math.hpp>
#include <kernel/lafem/container.hpp> // required for LAFEM::CloneMode

namespace FEAT
{
  namespace Global
  {
    /**
     * \brief Global vector wrapper class template
     *
     * This class implements a wrapper that contains a LAFEM vector as its core data object and
     * provides the necessary synchronization functions required in an MPI-parallel simulation
     * based on the overlapping domain decomposition approach. Effectively, this class only
     * couples a local LAFEM vector with its corresponding Global::Gate, which actually performs
     * all the dirty MPI work.
     *
     * This wrapper class implements all the 'black-box' functionality of a LAFEM vector and
     * includes the required synchronization for each operation. In the end, the only operations
     * which require direct synchronization/communication are the dot-product and all types of
     * norm functions. Most other functions (such as axpy, copy or scale) are parallel by nature
     * and do not require any type of synchronization, however, they should still be called on
     * all processes for the sake of consistency, unless a purely patch-local modification of
     * the vector is desired, of course.
     *
     * The class provides the two explicit synchronization functions sync_0() and synch_1(),
     * which synchronize the vector data between nearest neighbors in the overlapping domain
     * decomposition approach, thus ensuring that after all call to one of these functions,
     * each DOF, which is shared by more than one process, has the same value on each if these
     * processes, thus ensuring consistency of the DOF values across all processes.
     * The difference between these two function is as follows:
     * Each DOF, which is shared by more than one process, ...
     * - sync_0: ... is replaced by the \b sum of all its contributions from each process that
     *           shares this DOF.
     * - sync_1: ... is replaced by the \b average of all its contributions from each process that
     *           shares this DOF.
     *
     * The 'sync_0' function is used after each assembly of a dual (right-hand-side or defect)
     * vector to sum up all the linearform contributions for each DOF among all processes. This
     * function is also used to synchronize a dual vector after a matrix-vector multiplication.
     *
     * The 'sync_1' function is used whenever the local vectors are inconsistent over the processes,
     * e.g. after applying a patch-local preconditioner, by averaging each DOF over all processes
     * that contribute to that DOF.
     *
     * Finally, this class also offers a function named 'from_1_to_0', which effectively scales
     * each DOF by 1/K, where K is the number of processes that share the corresponding DOF. With
     * this function, a call to 'sync_1' is effectively equivalent to a call to 'from_1_to_0'
     * followed by a call to 'sync_0'.
     *
     * \tparam LocalVector_
     * The type of the local vector container; may be any valid combination of LAFEM (meta-)vector types
     *
     * \tparam Mirror_
     * The type of the vector mirror; must be compatible to the local vector type
     *
     * \author Peter Zajac
     */
    template<typename LocalVector_, typename Mirror_>
    class Vector
    {
    public:
      typedef Gate<LocalVector_, Mirror_> GateType;

      typedef typename LocalVector_::MemType MemType;
      typedef typename LocalVector_::DataType DataType;
      typedef typename LocalVector_::IndexType IndexType;
      typedef LocalVector_ LocalVectorType;

      /// Our 'base' class type
      template <typename LocalVector2_, typename Mirror2_ = Mirror_>
      using ContainerType = Vector<LocalVector2_, Mirror2_>;

      /// this typedef lets you create a vector container with new Memory, Datatype and Index types
      template <typename Mem2_, typename DataType2_, typename IndexType2_>
      using ContainerTypeByMDI = Vector<
        typename LocalVector_::template ContainerType<Mem2_, DataType2_, IndexType2_>,
        typename Mirror_::template MirrorType<Mem2_, DataType2_, IndexType2_> >;

    protected:
      /// a pointer to the gate responsible for synchronization
      const GateType* _gate;
      /// the internal local vector object
      LocalVector_ _vector;

    public:
      /// standard constructor
      Vector() :
        _gate(nullptr),
        _vector()
      {
      }

      /**
       * \brief Forwarding constructor
       *
       * \param[in] gate
       * A \resident pointer to the gate to be used for synchronization
       *
       * \param[in] args
       * The arguments that are to be passed to the local vector object constructor
       */
      template<typename... Args_>
      explicit Vector(const GateType* gate, Args_&&... args) :
        _gate(gate),
        _vector(std::forward<Args_>(args)...)
      {
      }

      /**
       * \brief Returns a reference to the internal local LAFEM vector object.
       *
       * \returns A reference to the internal local LAFEM vector object.
       */
      LocalVector_& local()
      {
        return _vector;
      }

      /**
       * \brief Returns a const reference to the internal local LAFEM vector object.
       *
       * \returns A const reference to the internal local LAFEM vector object.
       */
      const LocalVector_& local() const
      {
        return _vector;
      }

      template<typename OtherGlobalVector_>
      void convert(const GateType* gate, const OtherGlobalVector_ & other)
      {
        this->_gate = gate;
        this->_vector.convert(other.local());
      }

      /**
       * \brief Returns a const pointer to the internal gate of the vector.
       *
       * \returns A const pointer to the internal gate of the vector.
       */
      const GateType* get_gate() const
      {
        return _gate;
      }

      /**
       * \brief Returns a const pointer to the internal communicator of the gate of the vector.
       *
       * \returns a const pointer to the internal communicator of the gate of the vector.
       */
      const Dist::Comm* get_comm() const
      {
        return (_gate != nullptr ? _gate->get_comm() : nullptr);
      }

      /**
       * \brief Converts a type-1 vector into a type-0 vector
       *
       * This function scales each DOF (local vector value) x_i by 1/P_i, where P_i is the number
       * of processes that share the DOF x_i.
       *
       * \note This function is process-local and does not perform any communication, but it should
       * be not be called by individual processes only anyways.
       */
      void from_1_to_0()
      {
        if(_gate != nullptr)
        {
          _gate->from_1_to_0(this->_vector);
        }
      }

      /**
       * \brief Performs a type-0 synchronization of the vector, i.e. sums up all local DOF contributions
       *
       * This function performs a type-0 synchronization, i.e. sums up all local DOF contributions
       * for each DOF, by exchanging the DOF values of each shared DOF among all nearest neighbors.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      void sync_0()
      {
        if(_gate != nullptr)
          _gate->sync_0(_vector);
      }

      /**
       * \brief Performs a type-0 synchronization of the vector, i.e. sums up all local DOF contributions
       *
       * This function performs a type-0 synchronization, i.e. sums up all local DOF contributions
       * for each DOF, by exchanging the DOF values of each shared DOF among all nearest neighbors.
       *
       * \returns A SynchVectorTicket object that waits for the operation to complete.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      auto sync_0_async() -> decltype(_gate->sync_0_async(_vector))
      //decltype(_gate->sync_0(_vector)) sync_0_async()
      {
        return _gate->sync_0_async(_vector);
      }

      /**
       * \brief Performs a type-1 synchronization of the vector, i.e. averages all local DOF contributions
       *
       * This function performs a type-1 synchronization, i.e. averages all local DOF contributions
       * for each DOF, by exchanging the DOF values of each shared DOF among all nearest neighbors.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      void sync_1()
      {
        if(_gate != nullptr)
          _gate->sync_1(_vector);
      }

      /**
       * \brief Performs a type-1 synchronization of the vector, i.e. averages all local DOF contributions
       *
       * This function performs a type-1 synchronization, i.e. averages all local DOF contributions
       * for each DOF, by exchanging the DOF values of each shared DOF among all nearest neighbors.
       *
       * \returns A SynchVectorTicket object that waits for the operation to complete.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      auto sync_1_async() -> decltype(_gate->sync_1_async(_vector))
      //decltype(_gate->sync_1(_vector)) sync_1_async()
      {
        return _gate->sync_1_async(_vector);
      }

      /**
       * \brief Returns the total number of elements in this distributed vector
       *
       * \warning In parallel, this requires communication and is very expensive, so use sparingly!
       * \note This always returns the raw (or POD - Plain Old Data) size, as everything else is ambiguous.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       *
       * \returns The number of elements
       */
      Index size() const
      {
        // Compute total number of rows
        auto vec_l = clone();
        vec_l.format(DataType(1));

        return Index(vec_l.norm2sqr());
      }

      /**
       * \brief Creates and returns a clone of this global vector
       *
       * \param[in] mode
       * Specifies the clone mode for the internal vector object.
       *
       * \returns The created clone object
       */
      Vector clone(LAFEM::CloneMode mode = LAFEM::CloneMode::Weak) const
      {
        return Vector(_gate, _vector.clone(mode));
      }

      /**
       * \brief Creates this as a clone of another global vector
       *
       * \param[in] vector
       * A \transient reference to the vector that is to be cloned
       *
       * \param[in] mode
       * Specifies the clone mode for the internal vector object.
       */
      void clone(const Vector& other, LAFEM::CloneMode mode = LAFEM::CloneMode::Weak)
      {
        XASSERTM(&(other.local()) != &(this->local()), "Trying to self-clone a Global::Vector!");

        this->local() = other.local().clone(mode);
      }

      /**
       * \brief Clears the underlying vector.
       */
      void clear()
      {
        _vector.clear();
      }

      /**
       * \brief Reset all elements of the container to a given value or zero if missing.
       *
       * \param[in] alpha The value to be set (defaults to 0)
       */
      void format(DataType alpha = DataType(0))
      {
        _vector.format(alpha);
      }

      /**
       * \brief Reset all elements of the container to random values.
       *
       * \param[in] rng The random number generator.
       * \param[in] min Lower rng bound.
       * \param[in] max Upper rng bound.
       *
       * \note This function automatically synchronizes the resulting vector to ensure that the
       * shared DOFs are consistent over all processes.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      void format(Random & rng, DataType min, DataType max)
      {
        _vector.format(rng, min, max);
        sync_1();
      }

      /**
       * \brief Copies the contents of another vector into this vector
       *
       * \param[in] x
       * A \transient reference from which to copy the contents from.
       */
      void copy(const Vector& x)
      {
        // avoid self-copy
        if(this != &x)
        {
          _vector.copy(x.local());
        }
      }

      /**
       * \brief Performs an AXPY operation: this <- y + alpha*x
       *
       * \param[in] x, y
       * The \transient references to the two input vectors
       *
       * \param[in] alpha
       * The scaling factor for the input vector \p x
       */
      void axpy(const Vector& x, const Vector& y, const DataType alpha = DataType(1))
      {
        _vector.axpy(x.local(), y.local(), alpha);
      }

      /**
       * \brief Sets this to a scaled vector: this <- alpha*x
       *
       * \param[in] x
       * A \transient reference to the input vector that is to be scaled
       *
       * \param[in] alpha
       * The scaling factor for the input vector \p x
       */
      void scale(const Vector& x, const DataType alpha)
      {
        _vector.scale(x.local(), alpha);
      }

      /**
       * \brief Computes the dot-product of this vector and another vector
       *
       * \param[in] x
       * A \transient reference to the other vector for the dot-product
       *
       * \returns The dot-product of this and \p x
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      DataType dot(const Vector& x) const
      {
        if(_gate != nullptr)
          return _gate->dot(_vector, x.local());
        return _vector.dot(x.local());
      }

      /**
       * \brief Computes the dot-product of this vector and another vector
       *
       * \param[in] x
       * A \transient reference to the other vector for the dot-product
       *
       * \returns A scalar ticket that has to be waited upon to complete the operation.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      std::shared_ptr<SynchScalarTicket<DataType>> dot_async(const Vector& x) const
      {
        return _gate->dot_async(_vector, x.local());
      }

      /**
       * \brief Computes the squared Euclid norm of this vector
       *
       * \returns The squared Euclid norm of this vector
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      DataType norm2sqr() const
      {
        return dot(*this);
      }

      /**
       * \brief Computes the squared Euclid norm of this vector
       *
       * \returns A scalar ticket that has to be waited upon to complete the operation.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      std::shared_ptr<SynchScalarTicket<DataType>> norm2sqr_async() const
      {
        return dot_async(*this);
      }

      /**
       * \brief Computes the Euclid norm of this vector
       *
       * \returns The Euclid norm of this vector
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      DataType norm2() const
      {
        return Math::sqrt(norm2sqr());
      }

      /**
       * \brief Computes the Euclid norm of this vector
       *
       * \returns A scalar ticket that has to be waited upon to complete the operation.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      std::shared_ptr<SynchScalarTicket<DataType>> norm2_async() const
      {
        return _gate->dot_async(_vector, _vector, true);
      }

      /**
       * \brief Computes the component-wise inverse of a vector
       *
       * This function performs \f$ this_i \leftarrow \alpha / x_i \f$
       *
       * \param[in] x
       * A \transient reference to the vector whose components are to be inverted
       *
       * \param[in] alpha
       * The scaling factor for the component inversion
       */
      void component_invert(const Vector& x, const DataType alpha = DataType(1))
      {
        _vector.component_invert(x.local(), alpha);
      }

      /**
       * \brief Computes the component-wise product of two vector
       *
       * This function performs \f$ this_i \leftarrow x_i \cdot y_i \f$
       *
       * \param[in] x, y
       * The \transient references to the two vectors whose components are to be multiplied
       */
      void component_product(const Vector& x, const Vector& y)
      {
        _vector.component_product(x.local(), y.local());
      }

      /**
       * \brief Retrieve the absolute maximum value of this vector.
       *
       * \returns The largest absolute value of this vector.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      DataType max_abs_element() const
      {
        if (_gate != nullptr)
          return _gate->max_abs_element(_vector);
        return _vector.max_abs_element();
      }

      /**
       * \brief Retrieve the absolute maximum value of this vector.
       *
       * \returns A scalar ticket that has to be waited upon to complete the operation.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      std::shared_ptr<SynchScalarTicket<DataType>> max_abs_element_async() const
      {
        return _gate->max_abs_element_async(_vector);
      }

      /**
       * \brief Retrieve the absolute minimum value of this vector.
       *
       * \returns The smallest absolute value of this vector.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      DataType min_abs_element() const
      {
        if (_gate != nullptr)
          return _gate->min_abs_element(_vector);
        return _vector.min_abs_element();
      }

      /**
       * \brief Retrieve the absolute minimum value of this vector.
       *
       * \returns A scalar ticket that has to be waited upon to complete the operation.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      std::shared_ptr<SynchScalarTicket<DataType>> min_abs_element_async() const
      {
        return _gate->min_abs_element_async(_vector);
      }

      /**
       * \brief Retrieve the maximum value of this vector.
       *
       * \returns The largest value of this vector.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      DataType max_element() const
      {
        if (_gate != nullptr)
          return _gate->max_element(_vector);
        return _vector.max_element();
      }

      /**
       * \brief Retrieve the maximum value of this vector.
       *
       * \returns A scalar ticket that has to be waited upon to complete the operation.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      std::shared_ptr<SynchScalarTicket<DataType>> max_element_async() const
      {
        return _gate->max_element_async(_vector);
      }

      /**
       * \brief Retrieve the minimum value of this vector.
       *
       * \returns The smallest value of this vector.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      DataType min_element() const
      {
        if (_gate != nullptr)
          return _gate->min_element(_vector);
        return _vector.min_element();
      }

      /**
       * \brief Retrieve the minimum value of this vector.
       *
       * \returns A scalar ticket that has to be waited upon to complete the operation.
       *
       * \attention This function must be called by all processes participating in the gate's
       * communicator, otherwise the application will deadlock.
       */
      std::shared_ptr<SynchScalarTicket<DataType>> min_element_async() const
      {
        return _gate->min_element_async(_vector);
      }

      /// \copydoc FEAT::Control::Checkpointable::get_checkpoint_size()
      std::uint64_t get_checkpoint_size(LAFEM::SerialConfig& config)
      {
        return _vector.get_checkpoint_size(config);
      }

      /// \copydoc FEAT::Control::Checkpointable::restore_from_checkpoint_data(std::vector<char>&)
      void restore_from_checkpoint_data(std::vector<char>& data)
      {
        _vector.restore_from_checkpoint_data(data);
      }

      /// \copydoc FEAT::Control::Checkpointable::set_checkpoint_data(std::vector<char>&)
      std::uint64_t set_checkpoint_data(std::vector<char>& data, LAFEM::SerialConfig& config)
      {
        return _vector.set_checkpoint_data(data, config);
      }
    }; // class Vector<...>
  } // namespace Global
} // namespace FEAT

#endif // KERNEL_GLOBAL_VECTOR_HPP
