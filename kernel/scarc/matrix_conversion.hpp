#pragma once
#ifndef SCARC_GUARD_MATRIX_CONVERSION_HPP
#define SCARC_GUARD_MATRIX_CONVERSION_HPP 1

#include<kernel/base_header.hpp>
#include<kernel/scarc/scarc_error.hpp>
#include<kernel/scarc/scarc_data.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/vector_mirror.hpp>
#include <kernel/lafem/matrix_mirror.hpp>
#include <kernel/archs.hpp>
#include <kernel/assembly/mirror_assembler.hpp>

using namespace FEAST;
using namespace FEAST::Foundation;
using namespace FEAST::LAFEM;
using namespace FEAST::ScaRC;

namespace FEAST
{
  namespace ScaRC
  {
    ///type-0 to type-1 matrix conversion
    template<
             typename Mem_,
             typename DT_,
             typename IT_,
             template<typename, typename, typename> class MT_>
    struct MatrixConversion
    {
    };

    template<
             typename Mem_,
             typename DT_,
             typename IT_>
    struct MatrixConversion<Mem_, DT_, IT_, SparseMatrixCSR>
    {
      template<template<typename, typename> class ST_, typename VMT_>
      static SparseMatrixCSR<Mem_, DT_, IT_> value(const SparseMatrixCSR<Mem_, DT_, IT_>& origin, const ST_<VMT_, std::allocator<VMT_> >& vec_mirrors, const ST_<IT_, std::allocator<IT_> >& other_ranks)
      {
        SparseMatrixCSR<Mem_, DT_, IT_> result;
        result.clone(origin);
#ifndef SERIAL
        ST_<DenseVector<Mem_, DT_, IT_>, std::allocator<DenseVector<Mem_, DT_, IT_> > > val_sendbufs;
        ST_<DenseVector<Mem_, DT_, IT_>, std::allocator<DenseVector<Mem_, DT_, IT_> > > val_recvbufs;
        ST_<DenseVector<Mem_, IT_, IT_>, std::allocator<DenseVector<Mem_, IT_, IT_> > > colind_sendbufs;
        ST_<DenseVector<Mem_, IT_, IT_>, std::allocator<DenseVector<Mem_, IT_, IT_> > > colind_recvbufs;
        ST_<DenseVector<Mem_, IT_, IT_>, std::allocator<DenseVector<Mem_, IT_, IT_> > > rp_sendbufs;
        ST_<DenseVector<Mem_, IT_, IT_>, std::allocator<DenseVector<Mem_, IT_, IT_> > > rp_recvbufs;
        for(Index i(0) ; i < vec_mirrors.size() ; ++i)
        {
          MatrixMirror<Mem::Main, double> mat_mirror(vec_mirrors.at(i), vec_mirrors.at(i));
          SparseMatrixCSR<Mem::Main, double> buf_mat;
          Assembly::MirrorAssembler::assemble_buffer_matrix(buf_mat, mat_mirror, result);
          mat_mirror.gather_op(buf_mat, result);

          double* val(buf_mat.val());
          Index* row_ptr(buf_mat.row_ptr());
          Index* col_ind(buf_mat.col_ind());

          DenseVector<Mem::Main, double> val_sendbuf(buf_mat.used_elements());
          DenseVector<Mem::Main, double> val_recvbuf(buf_mat.used_elements());
          DenseVector<Mem::Main, Index> colind_sendbuf(buf_mat.used_elements());
          DenseVector<Mem::Main, Index> colind_recvbuf(buf_mat.used_elements());
          DenseVector<Mem::Main, Index> rp_sendbuf(buf_mat.rows() + 1);
          DenseVector<Mem::Main, Index> rp_recvbuf(buf_mat.rows() + 1);

          for(Index j(0) ; j < buf_mat.used_elements() ; ++j)
          {
            val_sendbuf(j, val[j]);
            colind_sendbuf(j, col_ind[j]);
          }
          for(Index j(0) ; j < buf_mat.rows() + 1 ; ++j)
          {
            rp_sendbuf(j, row_ptr[j]);
          }

          Status s1;
          Comm::send_recv(val_sendbuf.elements(),
              val_sendbuf.size(),
              other_ranks.at(i),
              val_recvbuf.elements(),
              val_recvbuf.size(),
              other_ranks.at(i),
              s1);
          Status s2;
          Comm::send_recv(colind_sendbuf.elements(),
              colind_sendbuf.size(),
              other_ranks.at(i),
              colind_recvbuf.elements(),
              colind_recvbuf.size(),
              other_ranks.at(i),
              s2);
          Status s3;
          Comm::send_recv(rp_sendbuf.elements(),
              rp_sendbuf.size(),
              other_ranks.at(i),
              rp_recvbuf.elements(),
              rp_recvbuf.size(),
              other_ranks.at(i),
              s3);

          val_sendbufs.push_back(std::move(val_sendbuf));
          val_recvbufs.push_back(std::move(val_recvbuf));
          colind_sendbufs.push_back(std::move(colind_sendbuf));
          colind_recvbufs.push_back(std::move(colind_recvbuf));
          rp_sendbufs.push_back(std::move(rp_sendbuf));
          rp_recvbufs.push_back(std::move(rp_recvbuf));
        }
        Comm::barrier();

        for(Index i(0) ; i < vec_mirrors.size() ; ++i)
        {
          MatrixMirror<Mem::Main, double> mat_mirror(vec_mirrors.at(i), vec_mirrors.at(i));
          SparseMatrixCSR<Mem::Main, double> buf_mat;
          Assembly::MirrorAssembler::assemble_buffer_matrix(buf_mat, mat_mirror, result);

          mat_mirror.gather_op(buf_mat, result);

          SparseMatrixCSR<Mem::Main, double> other_buf_mat(buf_mat.rows(),
              buf_mat.columns(),
              colind_recvbufs.at(i),
              val_recvbufs.at(i),
              rp_recvbufs.at(i));

          buf_mat.axpy<Algo::Generic>(buf_mat, other_buf_mat);
          mat_mirror.scatter_op(result, buf_mat);
        }
#endif
        return result;
      }
    };
  }
}

#endif
