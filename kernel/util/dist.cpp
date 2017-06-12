// includes, FEAT
#include <kernel/util/dist.hpp>
#include <kernel/util/math.hpp> // for ilog10

// includes, system
#include <cstring> // for strcpy, memcpy
#include <cstdint>

namespace FEAT
{
  namespace Dist
  {
#if defined(FEAT_HAVE_MPI) || defined(DOXYGEN)

#ifdef FEAT_OVERRIDE_MPI_OPS
    // operation: x <- x + y
    struct OpSum
    {
      template<typename T_>
      static inline void eval(T_& x, T_& y)
      {
        x += y;
      }
    };

    /// operation: x <- max(x, y)
    struct OpMax
    {
      template<typename T_>
      static inline void eval(T_& x, T_& y)
      {
        if(x < y)
          x = y;
      }
    };

    /// operation: x <- min(x, y)
    struct OpMin
    {
      template<typename T_>
      static inline void eval(T_& x, T_& y)
      {
        if(y < x)
          x = y;
      }
    };

    // Helper function:
    // Interprets iv and iov as arrays of type T_ and
    // applies Op_ on each pair of array elements
    template<typename Op_, typename T_>
    inline void func_op_t(void* iv, void* iov, int* n)
    {
      T_* y = reinterpret_cast<T_*>(iv);
      T_* x = reinterpret_cast<T_*>(iov);
      for(int i(0); i < *n; ++i)
      {
        Op_::eval(x[i], y[i]);
      }
    }

    // Callback function for MPI operation overrides
    template<typename Op_>
    void op_callback(void* iv, void* iov, int* n, MPI_Datatype* dt)
    {
      if(*dt == MPI_SIGNED_CHAR)         func_op_t<Op_, signed char>       (iv, iov, n); else
      if(*dt == MPI_SHORT)               func_op_t<Op_, short>             (iv, iov, n); else
      if(*dt == MPI_INT)                 func_op_t<Op_, int>               (iv, iov, n); else
      if(*dt == MPI_LONG)                func_op_t<Op_, long>              (iv, iov, n); else
      if(*dt == MPI_LONG_LONG)           func_op_t<Op_, long long>         (iv, iov, n); else
      if(*dt == MPI_UNSIGNED_CHAR)       func_op_t<Op_, unsigned char>     (iv, iov, n); else
      if(*dt == MPI_UNSIGNED_SHORT)      func_op_t<Op_, unsigned short>    (iv, iov, n); else
      if(*dt == MPI_UNSIGNED)            func_op_t<Op_, unsigned int>      (iv, iov, n); else
      if(*dt == MPI_UNSIGNED_LONG)       func_op_t<Op_, unsigned long>     (iv, iov, n); else
      if(*dt == MPI_UNSIGNED_LONG_LONG)  func_op_t<Op_, unsigned long long>(iv, iov, n); else
      if(*dt == MPI_FLOAT)               func_op_t<Op_, float>             (iv, iov, n); else
      if(*dt == MPI_DOUBLE)              func_op_t<Op_, double>            (iv, iov, n); else
      if(*dt == MPI_LONG_DOUBLE)         func_op_t<Op_, long double>       (iv, iov, n); else
      if(*dt == MPI_INT8_T)              func_op_t<Op_, std::int8_t>       (iv, iov, n); else
      if(*dt == MPI_INT16_T)             func_op_t<Op_, std::int16_t>      (iv, iov, n); else
      if(*dt == MPI_INT32_T)             func_op_t<Op_, std::int32_t>      (iv, iov, n); else
      if(*dt == MPI_INT64_T)             func_op_t<Op_, std::int64_t>      (iv, iov, n); else
      if(*dt == MPI_UINT8_T)             func_op_t<Op_, std::uint8_t>      (iv, iov, n); else
      if(*dt == MPI_UINT16_T)            func_op_t<Op_, std::uint16_t>     (iv, iov, n); else
      if(*dt == MPI_UINT32_T)            func_op_t<Op_, std::uint32_t>     (iv, iov, n); else
      if(*dt == MPI_UINT64_T)            func_op_t<Op_, std::uint64_t>     (iv, iov, n); else
#ifdef FEAT_HAVE_QUADMATH
      if(*dt == dt__float128.dt)         func_op_t<Op_, __float128>        (iv, iov, n); else
#endif
      {
        // unsupported datatype
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
#endif // FEAT_OVERRIDE_MPI_OPS

    bool initialise(int& argc, char**& argv)
    {
      // initialise MPI runtime first
      if(MPI_Init(&argc, &argv) != MPI_SUCCESS)
        return false;

#ifdef FEAT_HAVE_QUADMATH
      // Create a custom MPI datatype for '__float128'
      Datatype& dt_f128 = const_cast<Datatype&>(Dist::dt__float128);
      MPI_Type_contiguous(int(sizeof(__float128)), MPI_BYTE, &dt_f128.dt);
      MPI_Type_commit(&dt_f128.dt);
#endif // FEAT_HAVE_QUADMATH

#ifdef FEAT_OVERRIDE_MPI_OPS
      // override MPI operations
      MPI_Op_create(op_callback<OpSum>, 1, &const_cast<Operation&>(Dist::op_sum).op);
      MPI_Op_create(op_callback<OpMax>, 1, &const_cast<Operation&>(Dist::op_max).op);
      MPI_Op_create(op_callback<OpMin>, 1, &const_cast<Operation&>(Dist::op_min).op);
#endif // FEAT_OVERRIDE_MPI_OPS

      return true;
    }

    void finalise()
    {

#ifdef FEAT_OVERRIDE_MPI_OPS
      Operation& my_op_sum = const_cast<Operation&>(Dist::op_sum);
      MPI_Op_free(&my_op_sum.op);
#endif // FEAT_OVERRIDE_MPI_OPS

#ifdef FEAT_HAVE_QUADMATH
      Datatype& dt_f128 = const_cast<Datatype&>(Dist::dt__float128);
      MPI_Type_free(&dt_f128.dt);
#endif // FEAT_HAVE_QUADMATH

      // finalise MPI
      MPI_Finalize();
    }

    // datatypes
    const Datatype dt_byte               (MPI_BYTE,               sizeof(char));
    const Datatype dt_char               (MPI_CHAR,               sizeof(char));
    const Datatype dt_wchar              (MPI_WCHAR,              sizeof(wchar_t));
    const Datatype dt_signed_char        (MPI_SIGNED_CHAR,        sizeof(signed char));
    const Datatype dt_signed_short       (MPI_SHORT,              sizeof(short));
    const Datatype dt_signed_int         (MPI_INT,                sizeof(int));
    const Datatype dt_signed_long        (MPI_LONG,               sizeof(long));
    const Datatype dt_signed_long_long   (MPI_LONG_LONG,          sizeof(long long));
    const Datatype dt_unsigned_char      (MPI_UNSIGNED_CHAR,      sizeof(unsigned char));
    const Datatype dt_unsigned_short     (MPI_UNSIGNED_SHORT,     sizeof(unsigned short));
    const Datatype dt_unsigned_int       (MPI_UNSIGNED,           sizeof(unsigned int));
    const Datatype dt_unsigned_long      (MPI_UNSIGNED_LONG,      sizeof(unsigned long));
    const Datatype dt_unsigned_long_long (MPI_UNSIGNED_LONG_LONG, sizeof(unsigned long long));
    const Datatype dt_float              (MPI_FLOAT,              sizeof(float));
    const Datatype dt_double             (MPI_DOUBLE,             sizeof(double));
    const Datatype dt_long_double        (MPI_LONG_DOUBLE,        sizeof(long double));
    const Datatype dt_signed_int8        (MPI_INT8_T,             sizeof(std::int8_t));
    const Datatype dt_signed_int16       (MPI_INT16_T,            sizeof(std::int16_t));
    const Datatype dt_signed_int32       (MPI_INT32_T,            sizeof(std::int32_t));
    const Datatype dt_signed_int64       (MPI_INT64_T,            sizeof(std::int64_t));
    const Datatype dt_unsigned_int8      (MPI_UINT8_T,            sizeof(std::uint8_t));
    const Datatype dt_unsigned_int16     (MPI_UINT16_T,           sizeof(std::uint16_t));
    const Datatype dt_unsigned_int32     (MPI_UINT32_T,           sizeof(std::uint32_t));
    const Datatype dt_unsigned_int64     (MPI_UINT64_T,           sizeof(std::uint64_t));
#ifdef FEAT_HAVE_QUADMATH
    // This needs to initialised by Dist::initialise() !
    const Datatype dt__float128          (0,                      sizeof(__float128));
#endif

    // operations
    const Operation op_sum(MPI_SUM);
    const Operation op_max(MPI_MAX);
    const Operation op_min(MPI_MIN);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* MPI Request wrapper implementation                                                        */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    Request::Request() :
      request(MPI_REQUEST_NULL)
    {
    }

    Request::Request(MPI_Request req_) :
      request(req_)
    {
    }

    Request::Request(Request&& other) :
      request(other.request)
    {
      other.request = MPI_REQUEST_NULL;
    }

    Request& Request::operator=(Request&& other)
    {
      if(this != &other)
      {
        XASSERT(is_null());
        this->request = other.request;
        other.request = MPI_REQUEST_NULL;
      }
      return *this;
    }

    bool Request::is_null() const
    {
      return request == MPI_REQUEST_NULL;
    }

    void Request::free()
    {
      if(request != MPI_REQUEST_NULL)
        MPI_Request_free(&request);
    }

    void Request::cancel()
    {
      if(request != MPI_REQUEST_NULL)
        MPI_Cancel(&request);
    }

    bool Request::test(Status& status)
    {
      int flag(0);
      MPI_Test(&request, &flag, status.mpi_status());
      return flag != 0;
    }

    bool Request::wait(Status& status)
    {
      if(request == MPI_REQUEST_NULL)
        return false;
      MPI_Wait(&request, status.mpi_status());
      return true;
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* MPI RequestVector implementation                                                          */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    bool RequestVector::test_all()
    {
      int flag = 0;
      MPI_Testall(_isize(), _reqs_array(), &flag, _stats_array());
      return flag != 0;
    }

    bool RequestVector::test_any(std::size_t& idx, Status& status)
    {
      int i = 0;
      int flag = 0;
      MPI_Testany(_isize(), _reqs_array(), &i, &flag, status.mpi_status());
      if(flag == 0)
        return false;
      idx = std::size_t(i);
      return true;
    }

    /*std::size_t RequestVector::test_some(std::size_t* indices, Status* statuses)
    {
      static_assert(sizeof(int) <= sizeof(std::size_t), "size_t must be at least size of int");
      XASSERT(indices != nullptr);
      XASSERT(statuses != nullptr);

      int outcount = -1;
      int* idx = reinterpret_cast<int*>(indices);
      MPI_Testsome(_isize(), _reqs_array(), &outcount, idx, (MPI_Status*)statuses);

      // unpack int to size_t ?
      if(sizeof(int) < sizeof(std::size_t))
      {
        for(std::size_t k = std::size_t(outcount); k > 0;)
        {
          --k;
          indices[k] = std::size_t(idx[k]);
        }
      }

      return std::size_t(outcount);
    }*/

    void RequestVector::wait_all()
    {
      MPI_Waitall(_isize(), _reqs_array(), _stats_array());
    }

    bool RequestVector::wait_any(std::size_t& idx, Status& status)
    {
      int i = -1;
      MPI_Waitany(_isize(), _reqs_array(), &i, status.mpi_status());
      if(i == MPI_UNDEFINED)
        return false;
      idx = std::size_t(i);
      return true;
    }

    /*std::size_t RequestVector::wait_some(std::size_t* indices, Status* statuses)
    {
      static_assert(sizeof(int) <= sizeof(std::size_t), "size_t must be at least size of int");
      XASSERT(indices != nullptr);
      XASSERT(statuses != nullptr);

      int outcount = -1;
      int* idx = reinterpret_cast<int*>(indices);
      MPI_Waitsome(_isize(), _reqs_array(), &outcount, idx, (MPI_Status*)statuses);

      // unpack int to size_t ?
      if(sizeof(int) < sizeof(std::size_t))
      {
        for(std::size_t k = std::size_t(outcount); k > 0;)
        {
          --k;
          indices[k] = std::size_t(idx[k]);
        }
      }

      // store statuses
      for(std::size_t k(0); k < std::size_t(outcount); ++k)
        _stats[indices[k]] = statuses[k];

      return std::size_t(outcount);
    }*/

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* MPI Comm wrapper implementation                                                           */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    Comm::Comm() :
      comm(MPI_COMM_NULL),
      _rank(0),
      _size(0)
    {
    }

    Comm::Comm(MPI_Comm comm_) :
      comm(comm_),
      _rank(0),
      _size(0)
    {
      if(comm_ != MPI_COMM_NULL)
      {
        MPI_Comm_rank(comm, &_rank);
        MPI_Comm_size(comm, &_size);
      }
    }

    Comm::Comm(Comm&& other) :
      comm(other.comm),
      _rank(other._rank),
      _size(other._size)
    {
      // do not free the world comm
      if(other.comm != MPI_COMM_WORLD)
      {
        other.comm = MPI_COMM_NULL;
        other._rank = other._size = 0;
      }
    }

    Comm& Comm::operator=(Comm&& other)
    {
      if(this != &other)
      {
        XASSERT(comm == MPI_COMM_NULL);
        comm = other.comm;
        _rank = other._rank;
        _size = other._size;
        // do not free the world comm
        if(other.comm != MPI_COMM_WORLD)
        {
          other.comm = MPI_COMM_NULL;
          other._rank = other._size = 0;
        }
      }
      return *this;
    }

    Comm::~Comm()
    {
      // do not free the world comm
      if((comm != MPI_COMM_WORLD) && (comm != MPI_COMM_NULL))
        MPI_Comm_free(&comm);
    }

    Comm Comm::world()
    {
      return Comm(MPI_COMM_WORLD);
    }

    bool Comm::is_world() const
    {
      return (comm == MPI_COMM_WORLD);
    }

    bool Comm::is_null() const
    {
      return (comm == MPI_COMM_NULL);
    }

    Comm Comm::comm_dup() const
    {
      // do not duplicate world and null comms
      if((comm == MPI_COMM_WORLD) || (comm == MPI_COMM_NULL))
        return Comm(comm);

      // create a real duplicate
      MPI_Comm newcomm = MPI_COMM_NULL;
      MPI_Comm_dup(comm, &newcomm);

      return Comm(newcomm);
    }

    Comm Comm::comm_create_range_incl(int count, int first, int stride) const
    {
      XASSERT(count > 0);
      XASSERT((first >= 0) && (first < _size));
      XASSERT(stride > 0);

      // get this comm's group
      MPI_Group group = MPI_GROUP_NULL;
      MPI_Comm_group(comm, &group);

      // create sub-group
      MPI_Group newgroup = MPI_GROUP_NULL;
      int ranges[3] = {first, first + (count-1)*stride, stride};
      MPI_Group_range_incl(group, 1, &ranges, &newgroup);

      // create new comm from group
      MPI_Comm newcomm = MPI_COMM_NULL;
      MPI_Comm_create(comm, newgroup, &newcomm);

      // free group handles
      MPI_Group_free(&newgroup);
      MPI_Group_free(&group);

      return Comm(newcomm);
    }

    Comm Comm::comm_create_incl(int n, const int* ranks) const
    {
      XASSERT(n > 0);
      XASSERT(ranks != nullptr);

      // get this comm's group
      MPI_Group group = MPI_GROUP_NULL;
      MPI_Comm_group(comm, &group);

      // create sub-group
      MPI_Group newgroup = MPI_GROUP_NULL;
      MPI_Group_incl(group, n, ranks, &newgroup);

      // create new comm from group
      MPI_Comm newcomm = MPI_COMM_NULL;
      MPI_Comm_create(comm, newgroup, &newcomm);

      // free group handles
      MPI_Group_free(&newgroup);
      MPI_Group_free(&group);

      return Comm(newcomm);
    }

    Comm Comm::comm_split(int color, int key) const
    {
      MPI_Comm newcomm = MPI_COMM_NULL;
      MPI_Comm_split(comm, color, key, &newcomm);
      return Comm(newcomm);
    }

    void Comm::barrier() const
    {
      MPI_Barrier(comm);
    }

    Request Comm::ibarrier() const
    {
      MPI_Request req(MPI_REQUEST_NULL);
      MPI_Ibarrier(comm, &req);
      return Request(req);
    }

    void Comm::bcast(void* buffer, std::size_t count, const Datatype& datatype, int root) const
    {
      MPI_Bcast(buffer, int(count), datatype.dt, root, comm);
    }

    Request Comm::ibcast(void* buffer, std::size_t count, const Datatype& datatype, int root) const
    {
      MPI_Request req(MPI_REQUEST_NULL);
      MPI_Ibcast(buffer, int(count), datatype.dt, root, comm, &req);
      return Request(req);
    }

    void Comm::gather(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype, int root) const
    {
      MPI_Gather(sendbuf, int(sendcount), sendtype.dt, recvbuf, int(recvcount), recvtype.dt, root, comm);
    }

    Request Comm::igather(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype, int root) const
    {
      MPI_Request req(MPI_REQUEST_NULL);
      MPI_Igather(sendbuf, int(sendcount), sendtype.dt, recvbuf, int(recvcount), recvtype.dt, root, comm, &req);
      return Request(req);
    }

    void Comm::scatter(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype, int root) const
    {
      MPI_Scatter(sendbuf, int(sendcount), sendtype.dt, recvbuf, int(recvcount), recvtype.dt, root, comm);
    }

    Request Comm::iscatter(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype, int root) const
    {
      MPI_Request req(MPI_REQUEST_NULL);
      MPI_Iscatter(sendbuf, int(sendcount), sendtype.dt, recvbuf, int(recvcount), recvtype.dt, root, comm, &req);
      return Request(req);
    }

    void Comm::allgather(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype) const
    {
      MPI_Allgather(sendbuf, int(sendcount), sendtype.dt, recvbuf, int(recvcount), recvtype.dt, comm);
    }

    Request Comm::iallgather(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype) const
    {
      MPI_Request req(MPI_REQUEST_NULL);
      MPI_Iallgather(sendbuf, int(sendcount), sendtype.dt, recvbuf, int(recvcount), recvtype.dt, comm, &req);
      return Request(req);
    }

    void Comm::allgatherv(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, const int* recvcounts, const int* displs, const Datatype& recvtype) const
    {
      MPI_Allgatherv(sendbuf, int(sendcount), sendtype.dt, recvbuf, recvcounts, displs, recvtype.dt, comm);
    }

    void Comm::alltoall(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype) const
    {
      MPI_Alltoall(sendbuf, int(sendcount), sendtype.dt, recvbuf, int(recvcount), recvtype.dt, comm);
    }

    Request Comm::ialltoall(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype) const
    {
      MPI_Request req(MPI_REQUEST_NULL);
#ifdef MSMPI_VER
      MPI_Alltoall(sendbuf, int(sendcount), sendtype.dt, recvbuf, int(recvcount), recvtype.dt, comm);
#else
      MPI_Ialltoall(sendbuf, int(sendcount), sendtype.dt, recvbuf, int(recvcount), recvtype.dt, comm, &req);
#endif
      return Request(req);
    }

    void Comm::alltoallv(const void* sendbuf, const int* sendcounts, const int* sdispls, const Datatype& sendtype, void* recvbuf, const int* recvcounts, const int* rdispls, const Datatype& recvtype) const
    {
      MPI_Alltoallv(sendbuf, sendcounts, sdispls, sendtype.dt, recvbuf, recvcounts, rdispls, recvtype.dt, comm);
    }

    void Comm::reduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op, int root) const
    {
      MPI_Reduce(sendbuf == recvbuf ? MPI_IN_PLACE : sendbuf, recvbuf, int(count), datatype.dt, op.op, root, comm);
    }

    Request Comm::ireduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op, int root) const
    {
      MPI_Request req(MPI_REQUEST_NULL);
      MPI_Ireduce(sendbuf == recvbuf ? MPI_IN_PLACE : sendbuf, recvbuf, int(count), datatype.dt, op.op, root, comm, &req);
      return Request(req);
    }

    void Comm::allreduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const
    {
      MPI_Allreduce(sendbuf == recvbuf ? MPI_IN_PLACE : sendbuf, recvbuf, int(count), datatype.dt, op.op, comm);
    }

    Request Comm::iallreduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const
    {
      MPI_Request req(MPI_REQUEST_NULL);
      MPI_Iallreduce(sendbuf == recvbuf ? MPI_IN_PLACE : sendbuf, recvbuf, int(count), datatype.dt, op.op, comm, &req);
      return Request(req);
    }

    void Comm::scan(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const
    {
      MPI_Scan(sendbuf == recvbuf ? MPI_IN_PLACE : sendbuf, recvbuf, int(count), datatype.dt, op.op, comm);
    }

    /*Request Comm::iscan(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const
    {
      MPI_Request req(MPI_REQUEST_NULL);
#ifdef MSMPI_VER
      MPI_Scan(sendbuf == recvbuf ? MPI_IN_PLACE : sendbuf, recvbuf, int(count), datatype.dt, op.op, comm);
#else
      MPI_Iscan(sendbuf == recvbuf ? MPI_IN_PLACE : sendbuf, recvbuf, int(count), datatype.dt, op.op, comm, &req);
#endif
      return Request(req);
    }*/

    void Comm::exscan(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const
    {
      MPI_Exscan(sendbuf == recvbuf ? MPI_IN_PLACE : sendbuf, recvbuf, int(count), datatype.dt, op.op, comm);
    }

    /*Request Comm::iexscan(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const
    {
      MPI_Request req(MPI_REQUEST_NULL);
#ifdef MSMPI_VER
      MPI_Exscan(sendbuf == recvbuf ? MPI_IN_PLACE : sendbuf, recvbuf, int(count), datatype.dt, op.op, comm);
#else
      MPI_Iexscan(sendbuf == recvbuf ? MPI_IN_PLACE : sendbuf, recvbuf, int(count), datatype.dt, op.op, comm, &req);
#endif
      return Request(req);
    }*/

    void Comm::send(const void* buffer, std::size_t count, const Datatype& datatype, int dest, int tag) const
    {
      MPI_Send(buffer, int(count), datatype.dt, dest, tag, comm);
    }

    Request Comm::isend(const void* buffer, std::size_t count, const Datatype& datatype, int dest, int tag) const
    {
      MPI_Request req(MPI_REQUEST_NULL);
      MPI_Isend(buffer, int(count), datatype.dt, dest, tag, comm, &req);
      return Request(req);
    }

    void Comm::recv(void* buffer, std::size_t count, const Datatype& datatype, int source, int tag, Status& status) const
    {
      MPI_Recv(buffer, int(count), datatype.dt, source, tag, comm, status.mpi_status());
    }

    Request Comm::irecv(void* buffer, std::size_t count, const Datatype& datatype, int source, int tag) const
    {
      MPI_Request req(MPI_REQUEST_NULL);
      MPI_Irecv(buffer, int(count), datatype.dt, source, tag, comm, &req);
      return Request(req);
    }

    void Comm::bcast_stringstream(std::stringstream& stream, int root) const
    {
      std::string str;
      std::size_t len(0);

      // Get size to broadcast
      if(rank() == root)
      {
        str = stream.str();
        len = str.length();
      }

      // broadcast string length
      bcast(&len, std::size_t(1), root);

      // empty string?
      if(len == std::size_t(0))
        return;

      // allocate buffer
      std::vector<char> buf(len + std::size_t(1));

      // fill buffer on root
      if(rank() == root)
        std::strcpy(buf.data(), str.c_str());

      // broadcast buffer
      bcast(buf.data(), buf.size(), root);

      // convert
      if(rank() != root)
        stream << std::string(buf.data(), len);
    }

    void Comm::bcast_binarystream(BinaryStream& stream, int root) const
    {
      // get the stream's internal data container
      std::vector<char>& data = stream.container();
      std::size_t len(0);

      // Get size to broadcast
      if(rank() == root)
        len = data.size();

      // broadcast vector length
      bcast(&len, std::size_t(1), root);

      // empty vector?
      if(len == std::size_t(0))
        return;

      // reallocate vector
      if(rank() != root)
        data.resize(len);

      // broadcast data
      bcast(data.data(), data.size(), root);
    }

    void Comm::print(std::ostream& os, const String& msg, int root) const
    {
      XASSERTM((0 <= root) && (root < _size), "invalid root rank argument");
      if(root == _rank)
        os << msg << std::endl;
    }

    void Comm::allprint(std::ostream& os, const String& msg, int root) const
    {
      XASSERTM((0 <= root) && (root < _size), "invalid root rank argument");

      if(_rank == root)
      {
        // determine the rank padding size first:
        const std::size_t n = std::size_t(_size);
        const std::size_t ndig = Math::ilog10(n);

        // gather message lengths
        std::size_t mylen = msg.size();
        std::vector<std::size_t> lengths(n);
        this->gather(&mylen, std::size_t(1), lengths.data(), std::size_t(1), root);

        // receive and print messages
        for(int i(0); i < _size; ++i)
        {
          // allocate message buffer
          std::size_t length = lengths.at(std::size_t(i));
          if(length == std::size_t(0))
            continue;

          // get the message to be printed
          String message;
          if(i == root)
            message = msg;
          else
          {
            // receive message
            std::vector<char> msgbuf(length);
            this->recv(msgbuf.data(), length, i);
            message.assign(msgbuf.data(), length);
          }

          // split message into single lines
          std::vector<String> lines;
          message.split_by_charset(lines, "\n");

          // set line prefix
          String prefix = stringify(i).pad_front(ndig);

          // print all lines prefixed
          for(auto it = lines.begin(); it != lines.end(); ++it)
            os << '[' << prefix << "] " << (*it) << std::endl;
        }
      }
      else // rank != root
      {
        // send the message length via gather
        std::size_t dummy(0), mylen = msg.size();
        this->gather(&mylen, std::size_t(1), &dummy, std::size_t(1), root);

        // send the message itself
        if(mylen > std::size_t(0))
          this->send(msg.data(), mylen, root);
      }
    }

    /* ######################################################################################### */
    /* ######################################################################################### */
    /* ######################################################################################### */
#else // non-MPI build
    /* ######################################################################################### */
    /* ######################################################################################### */
    /* ######################################################################################### */

    bool initialise(int& /*argc*/, char**& /*argv*/)
    {
      // nothing to do here
      return true;
    }

    void finalise()
    {
      // nothing to do here
    }

    // datatypes
    // Note: The numbers for the datatypes are arbitrary; it is just important that each datatype
    // has its own unique number.
    const Datatype dt_byte               ( 1, sizeof(char));
    const Datatype dt_char               ( 2, sizeof(char));
    const Datatype dt_wchar              ( 3, sizeof(wchar_t));
    const Datatype dt_signed_char        (11, sizeof(signed char));
    const Datatype dt_signed_short       (12, sizeof(short));
    const Datatype dt_signed_int         (13, sizeof(int));
    const Datatype dt_signed_long        (14, sizeof(long));
    const Datatype dt_signed_long_long   (15, sizeof(long long));
    const Datatype dt_unsigned_char      (21, sizeof(unsigned char));
    const Datatype dt_unsigned_short     (22, sizeof(unsigned short));
    const Datatype dt_unsigned_int       (23, sizeof(unsigned int));
    const Datatype dt_unsigned_long      (24, sizeof(unsigned long));
    const Datatype dt_unsigned_long_long (25, sizeof(unsigned long long));
    const Datatype dt_float              (31, sizeof(float));
    const Datatype dt_double             (32, sizeof(double));
    const Datatype dt_long_double        (33, sizeof(long double));
#ifdef FEAT_HAVE_QUADMATH
    const Datatype dt__float128          (34, sizeof(__float128));
#endif
    const Datatype dt_signed_int8        (41, sizeof(std::int8_t));
    const Datatype dt_signed_int16       (42, sizeof(std::int16_t));
    const Datatype dt_signed_int32       (43, sizeof(std::int32_t));
    const Datatype dt_signed_int64       (44, sizeof(std::int64_t));
    const Datatype dt_unsigned_int8      (51, sizeof(std::uint8_t));
    const Datatype dt_unsigned_int16     (52, sizeof(std::uint16_t));
    const Datatype dt_unsigned_int32     (53, sizeof(std::uint32_t));
    const Datatype dt_unsigned_int64     (54, sizeof(std::uint64_t));

    // operations
    const Operation op_sum(1);
    const Operation op_max(2);
    const Operation op_min(3);

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* Dummy Request wrapper implementation                                                      */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    Request::Request() :
      request(1)
    {
    }


    Request::Request(Request&& other) :
      request(other.request)
    {
      other.request = 0;
    }

    Request& Request::operator=(Request&& other)
    {
      if(this != &other)
      {
        this->request = other.request;
        other.request = 0;
      }
      return *this;
    }

    bool Request::is_null() const
    {
      return request == 0;
    }

    void Request::free()
    {
      request = 0;
    }

    void Request::cancel()
    {
    }

    bool Request::test(Status&)
    {
      request = 0;
      return true;
    }

    bool Request::wait(Status&)
    {
      bool ret = (request != 0);
      request = 0;
      return ret;
    }

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* MPI RequestVector implementation                                                          */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    bool RequestVector::test_all()
    {
      free();
      return true;
    }

    bool RequestVector::test_any(std::size_t& idx, Status&)
    {
      // try to find one non-null request
      for(std::size_t i(0); i < _reqs.size(); ++i)
      {
        if(!_reqs.at(i).is_null())
        {
          idx = i;
          _reqs.at(i).free();
          return true;
        }
      }
      // no active requests in vector
      return false;
    }

    /*std::size_t RequestVector::test_some(std::size_t* indices, Status* statuses)
    {
      XASSERT(indices != nullptr);
      XASSERT(statuses != nullptr);

      // return all non-null requests
      std::size_t k = 0;
      for(std::size_t i(0); i < _reqs.size(); ++i)
      {
        if(!_reqs.at(i).is_null())
        {
          indices[k++] = i;
          _reqs.at(i).free();
        }
      }
      return k;
    }*/

    void RequestVector::wait_all()
    {
      free();
    }

    bool RequestVector::wait_any(std::size_t& idx, Status& status)
    {
      return test_any(idx, status);
    }

    /*std::size_t RequestVector::wait_some(std::size_t* indices, Status* statuses)
    {
      return test_some(indices, statuses);
    }*/

    /* ***************************************************************************************** */
    /* ***************************************************************************************** */
    /* Dummy Comm wrapper implementation                                                         */
    /* ***************************************************************************************** */
    /* ***************************************************************************************** */

    Comm::Comm() :
      _rank(0),
      _size(0)
    {
    }

    Comm::Comm(int) :
      _rank(0),
      _size(1)
    {
    }

    Comm::Comm(Comm&& other) :
      _rank(other._rank),
      _size(other._size)
    {
    }

    Comm& Comm::operator=(Comm&& other)
    {
      if(this != &other)
      {
        _rank = other._rank;
        _size = other._size;
      }
      return *this;
    }

    Comm::~Comm()
    {
    }

    Comm Comm::world()
    {
      return Comm(1);
    }

    bool Comm::is_null() const
    {
      return _size == 0;
    }

    bool Comm::is_world() const
    {
      return _size == 1;
    }

    Comm Comm::comm_dup() const
    {
      return (_size == 0) ? Comm() : Comm(0);
    }

    Comm Comm::comm_create_range_incl(int, int, int) const
    {
      return Comm(1);
    }

    Comm Comm::comm_create_incl(int, const int*) const
    {
      return Comm(1);
    }

    Comm Comm::comm_split(int, int) const
    {
      return Comm(1);
    }

    void Comm::barrier() const
    {
      // nothing to do
    }

    Request Comm::ibarrier() const
    {
      // nothing to do
      return Request();
    }

    void Comm::bcast(void*, std::size_t, const Datatype&, int) const
    {
      // nothing to do
    }

    Request Comm::ibcast(void*, std::size_t, const Datatype&, int) const
    {
      // nothing to do
      return Request();
    }

    void Comm::gather(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype, int root) const
    {
      XASSERT(root == 0);
      alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype);
    }

    Request Comm::igather(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype, int root) const
    {
      XASSERT(root == 0);
      alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype);
      return Request();
    }

    void Comm::scatter(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype, int root) const
    {
      XASSERT(root == 0);
      alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype);
    }

    Request Comm::iscatter(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype, int root) const
    {
      XASSERT(root == 0);
      alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype);
      return Request();
    }

    void Comm::allgather(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype) const
    {
      alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype);
    }

    Request Comm::iallgather(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype) const
    {
      alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype);
      return Request();
    }

    void Comm::allgatherv(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, const int* recvcounts, const int* displs, const Datatype& recvtype) const
    {
      XASSERT(displs[0] == 0);
      alltoall(sendbuf, sendcount, sendtype, recvbuf, std::size_t(recvcounts[0]), recvtype);
    }

    void Comm::alltoall(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype) const
    {
      XASSERT(sendcount == recvcount);
      XASSERT(sendtype == recvtype);
      if(recvbuf != sendbuf)
        memcpy(recvbuf, sendbuf, sendcount * sendtype.size());
    }

    Request Comm::ialltoall(const void* sendbuf, std::size_t sendcount, const Datatype& sendtype, void* recvbuf, std::size_t recvcount, const Datatype& recvtype) const
    {
      alltoall(sendbuf, sendcount, sendtype, recvbuf, recvcount, recvtype);
      return Request();
    }

    void Comm::alltoallv(const void* sendbuf, const int* sendcounts, const int* sdispls, const Datatype& sendtype, void* recvbuf, const int* recvcounts, const int* rdispls, const Datatype& recvtype) const
    {
      XASSERT(sdispls[0] == 0);
      XASSERT(rdispls[0] == 0);
      alltoall(sendbuf, std::size_t(sendcounts[0]), sendtype, recvbuf, std::size_t(recvcounts[0]), recvtype);
    }

    void Comm::reduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op, int root) const
    {
      XASSERT(root == 0);
      allreduce(sendbuf, recvbuf, count, datatype, op);
    }

    Request Comm::ireduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op, int root) const
    {
      XASSERT(root == 0);
      allreduce(sendbuf, recvbuf, count, datatype, op);
      return Request();
    }

    void Comm::allreduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation&) const
    {
      if(sendbuf != recvbuf)
        memcpy(recvbuf, sendbuf, count * datatype.size());
    }

    Request Comm::iallreduce(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const
    {
      allreduce(sendbuf, recvbuf, count, datatype, op);
      return Request();
    }

    void Comm::scan(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const
    {
      allreduce(sendbuf, recvbuf, count, datatype, op);
    }

    /*Request Comm::iscan(const void* sendbuf, void* recvbuf, std::size_t count, const Datatype& datatype, const Operation& op) const
    {
      allreduce(sendbuf, recvbuf, count, datatype, op);
      return Request();
    }*/

    void Comm::exscan(const void*, void*, std::size_t, const Datatype&, const Operation&) const
    {
      // nothing to do
    }

    /*Request Comm::iexscan(const void*, void*, std::size_t, const Datatype&, const Operation&) const
    {
      // nothing to do
      return Request();
    }*/

    void Comm::send(const void*, std::size_t, const Datatype&, int, int) const
    {
      // nothing to do
    }

    Request Comm::isend(const void*, std::size_t, const Datatype&, int, int) const
    {
      // nothing to do
      return Request();
    }

    void Comm::recv(void*, std::size_t, const Datatype&, int, int, Status&) const
    {
      // nothing to do
    }

    Request Comm::irecv(void*, std::size_t, const Datatype&, int, int) const
    {
      // nothing to do
      return Request();
    }

    void Comm::bcast_stringstream(std::stringstream&, int) const
    {
      // nothing to do
    }

    void Comm::bcast_binarystream(BinaryStream&, int) const
    {
      // nothing to do
    }

    void Comm::print(std::ostream& os, const String& msg, int) const
    {
      os << msg << std::endl;
    }

    void Comm::allprint(std::ostream& os, const String& msg, int) const
    {
      os << msg << std::endl;
    }
#endif // FEAT_HAVE_MPI
  } // namespace Dist
} // namespace FEAT
