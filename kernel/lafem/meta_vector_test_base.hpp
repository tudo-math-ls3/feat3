// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
// #include <kernel/archs.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/power_vector.hpp>
#include <kernel/lafem/tuple_vector.hpp>

namespace FEAT
{
  namespace LAFEM
  {
    /**
     * \brief Abstract base class for meta-vector tests.
     *
     * \author Peter Zajac
     */
    template<typename DataType_, typename IndexType_>
    class MetaVectorTestBase
      : public FEAT::TestSystem::UnitTest
    {
    public:
      typedef DataType_ DataType;
      typedef IndexType_ IndexType;

      typedef DenseVector<DataType, IndexType> ScalarVector;
      typedef PowerVector<ScalarVector, 2> PowerVector2;
      typedef TupleVector<PowerVector2, ScalarVector> MetaVector;

      explicit MetaVectorTestBase(const String& id_in, const String datatype_name = "none", const String index_name = "none", PreferredBackend preferred_backend = PreferredBackend::generic)
       : FEAT::TestSystem::UnitTest(id_in, datatype_name, index_name, preferred_backend)
      {
      }

      static DataType fx00(Index i)
      {
        return DataType(0.2) * DataType(i);
      }

      static DataType fx01(Index i)
      {
        return DataType(2) - DataType(i);
      }

      static DataType fx1(Index i)
      {
        return DataType(0.01) * Math::sqr(DataType(i+1));
      }

      static DataType fy00(Index i)
      {
        return -DataType(3) + DataType(2*i);
      }

      static DataType fy01(Index i)
      {
        return DataType(1) + Math::sqrt(DataType(i+1));
      }

      static DataType fy1(Index i)
      {
        return DataType(1) + DataType(7) / DataType(i+1);
      }

      static DataType fz00(Index i)
      {
        return Math::sqr(DataType(i+1));
      }

      static DataType fz01(Index i)
      {
        return -DataType(5) - Math::sqrt(DataType(i+1));
      }

      static DataType fz1(Index i)
      {
        return - DataType(7) / DataType(i+1);
      }

      // generate test-vector x
      static MetaVector gen_vector_x(Index n00, Index n01, Index n1)
      {
        DenseVector<DataType, IndexType> r00(n00), r01(n01), r1(n1);

        // fill vectors
        for(Index i(0); i < n00; ++i)
          r00(i, fx00(i));
        for(Index i(0); i < n01; ++i)
          r01(i, fx01(i));
        for(Index i(0); i < n1; ++i)
          r1(i, fx1(i));

        // construct data vector
        MetaVector x;
        x.template at<0>().template at<0>().convert(r00);
        x.template at<0>().template at<1>().convert(r01);
        x.template at<1>().convert(r1);
        return x;
      }

      // generate test-vector y
      static MetaVector gen_vector_y(Index n00, Index n01, Index n1)
      {
        DenseVector<DataType, IndexType> r00(n00), r01(n01), r1(n1);

        // fill vectors
        for(Index i(0); i < n00; ++i)
          r00(i, fy00(i));
        for(Index i(0); i < n01; ++i)
          r01(i, fy01(i));
        for(Index i(0); i < n1; ++i)
          r1(i, fy1(i));

        // construct data vector
        MetaVector y;
        y.template at<0>().template at<0>().convert(r00);
        y.template at<0>().template at<1>().convert(r01);
        y.template at<1>().convert(r1);
        return y;
      }

      // Generate test-vector z. It is crucial that z_i != 0
      static MetaVector gen_vector_z(Index n00, Index n01, Index n1)
      {
        DenseVector<DataType, IndexType> r00(n00), r01(n01), r1(n1);

        // fill vectors
        for(Index i(0); i < n00; ++i)
          r00(i, fz00(i));
        for(Index i(0); i < n01; ++i)
          r01(i, fz01(i));
        for(Index i(0); i < n1; ++i)
          r1(i, fz1(i));

        // construct data vector
        MetaVector z;
        z.template at<0>().template at<0>().convert(r00);
        z.template at<0>().template at<1>().convert(r01);
        z.template at<1>().convert(r1);
        return z;
      }

      // generate null vector
      static MetaVector gen_vector_null(Index n00, Index n01, Index n1)
      {
        MetaVector x;
        x.template at<0>().template at<0>() = ScalarVector(n00, DataType(0));
        x.template at<0>().template at<1>() = ScalarVector(n01, DataType(0));
        x.template at<1>() = ScalarVector(n1, DataType(0));
        return x;
      }
    }; // MetaVectorTestBase
  } // namespace LAFEM
} // namespace FEAT
