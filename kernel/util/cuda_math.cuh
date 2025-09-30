// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#ifndef __CUDACC__
static_assert(false, "Tyring to include a cuh header in a non nvcc compile unit!")
#endif

// This cuda specific header defines templated wrappers around useful mathematical operations.
// You should never include this into any .c, .cpp, .h or .hpp file.

// for numerical limits
#include <float.h>
#include <cooperative_groups.h>
#include <cooperative_groups/reduce.h>
#ifdef FEAT_HAVE_HALFMATH
#include "cuda_fp16.h"
#endif

namespace cg = cooperative_groups;

namespace FEAT
{
  /**
   * \brief FEAT CudaMath namespace
   *
   * This namespace encapsulates templated wrappers around useful cuda math functions.
   * This *does not* contain all cmath or cudamath functions available, yet.
   *
   * \author Maximilian Esser
   */
  namespace CudaMath
  {
    /** \brief Calculates the square root of a value.
     *
     * \tparam DT_ The used dataytpe.
     *
     * \param[in] val The positive value in. No runtime check if this is met.
     *
     * \return The square root of val.
     */
    template<typename DT_>
    __host__ __device__ DT_ __forceinline__ cuda_sqrt(DT_ val);

    // /** \copydoc cuda_sqrt() */
    template<> __host__ __device__ __forceinline__ double cuda_sqrt<double>(double val)
    {
      return sqrt(val);
    }

    // /** \copydoc cuda_sqrt() */
    template<> __host__ __device__ __forceinline__ float cuda_sqrt<float>(float val)
    {
      return sqrtf(val);
    }

    #ifdef FEAT_HAVE_HALFMATH
    // /** \copydoc cuda_sqrt() */
    template<> __host__ __device__ __forceinline__ Half cuda_sqrt<Half>(Half val)
    {
      #ifdef __CUDA_ARCH__
      return hsqrt(val);
      #else
      return __float2half(sqrtf(__half2float(val)));
      #endif
    }
    #endif

    /** \brief Calculates the inverse square root of a value.
     *
     * \tparam DT_ The used dataytpe.
     *
     * \param[in] val The positive value in. No runtime check if this is met.
     *
     * \return The positive value of 1/sqrt(val).
     */
    template<typename DT_>
    __host__ __device__ DT_ cuda_rsqrt(DT_ val);

    // /** \copydoc cuda_rsqrt() */
    template<> __host__ __device__ __forceinline__ double cuda_rsqrt<double>(double val)
    {
      return rsqrt(val);
    }

    // /** \copydoc cuda_rsqrt() */
    template<> __host__ __device__ __forceinline__ float cuda_rsqrt<float>(float val)
    {
      return rsqrtf(val);
    }

    #ifdef FEAT_HAVE_HALFMATH
    // /** \copydoc cuda_sqrt() */
    template<> __host__ __device__ __forceinline__ Half cuda_rsqrt<Half>(Half val)
    {
      #ifdef __CUDA_ARCH__
      return hrsqrt(val);
      #else
      return __float2half(rsqrtf(__half2float(val)));
      #endif
    }
    #endif

    /** \brief Calculates the square of a value.
     *
     * \tparam DT_ The used dataytpe.
     *
     * \param[in] val The value in.
     *
     * \return The value of val * val.
     */
    template<typename DT_>
    __host__ __device__ __forceinline__ DT_ cuda_sqr(DT_ val)
    {
      return val * val;
    }

    /** \brief Calculates the absolute value.
     *
     * \tparam DT_ The used dataytpe.
     *
     * \param[in] val The value in.
     *
     * \return The value of |val|.
     */
    template<typename DT_>
    __host__ __device__ __forceinline__ DT_ cuda_abs(DT_);

    // /** \copydoc cuda_abs() */
    template<>
    __host__ __device__ __forceinline__ double cuda_abs<double>(double val)
    {
      return fabs(val);
    }

    // /** \copydoc cuda_abs() */
    template<>
    __host__ __device__ __forceinline__ float cuda_abs<float>(float val)
    {
      return fabsf(val);
    }


    #ifdef FEAT_HAVE_HALFMATH
    template<>
    __host__ __device__ __forceinline__ Half cuda_abs<Half>(Half val)
    {
      return __habs(val);
    }
    #endif

    /** \brief Calculates the power operation of two values.
     *
     * \tparam DT_ The used dataytpe.
     *
     * \param[in] val The base value.
     *
     * \param[in] exp The exponent.
     *
     * \return The value of val^(exp).
     */
    template<typename DT_>
    __host__ __device__ DT_ cuda_pow(DT_ val, DT_ exp);

    // /** \copydoc cuda_pow() */
    template<> __host__ __device__ __forceinline__ double cuda_pow<double>(double val, double exp)
    {
      return pow(val, exp);
    }

    // /** \copydoc cuda_pow() */
    template<> __host__ __device__ __forceinline__ float cuda_pow<float>(float val, float exp)
    {
      return powf(val, exp);
    }

    #ifdef FEAT_HAVE_HALFMATH
    // /** \copydoc cuda_pow() */
    template<> __host__ __device__ __forceinline__ Half cuda_pow<Half>(Half val, Half exp)
    {
      return __float2half(powf(__half2float(val), __half2float(exp)));
    }
    #endif

    /**
     * \brief Checks if a value is normal.
     *
     * \tparam DT_ The datatype
     *
     * \param[in] val The value to be checked.
     *
     * \returns True if val is not NaN, inf, zero or subnormal.
     */
    template<typename DT_>
    __host__ __device__ __forceinline__ bool cuda_isnormal(DT_ val);

    template<> __host__ __device__ __forceinline__ bool cuda_isnormal<double>(double val)
    {
      return !((isfinite(val) == 0) || abs(val) < DBL_MIN);
    }

    template<> __host__ __device__ __forceinline__ bool cuda_isnormal<float>(float val)
    {
      return !((isfinite(val) == 0) || abs(val) < FLT_MIN);
    }

    #ifdef FEAT_HAVE_HALFMATH
    template<> __host__ __device__ __forceinline__ bool cuda_isnormal<Half>(Half val)
    {
      return ((__hisinf(val) || __hisnan(val)) && (__habs(val) >= CUDART_MIN_DENORM_FP16));
    }
    #endif


    template<typename DT_>
    __host__ __device__ __forceinline__ bool cuda_isnan(DT_ val)
    {
      return isnan(val);
    }

    #ifdef FEAT_HAVE_HALFMATH
    template<> __host__ __device__ __forceinline__ bool cuda_isnan<Half>(Half val)
    {
      return  __hisnan(val);
    }
    #endif


    /**
     * \brief Evaluates the sine function.
     *
     * \tparam DT_ The datatype.
     *
     * \param[in] val The value to be evaluated.
     *
     * \returns The value cos(val)
     */
    template<typename DT_>
    __host__ __device__ __forceinline__ DT_ cuda_sin(DT_ val);

    /**
     * \brief Evaluates the cosine function.
     *
     * \tparam DT_ The datatype.
     *
     * \param[in] val The value to be evaluated.
     *
     * \returns The value sin(val)
     */
    template<typename DT_>
    __host__ __device__ __forceinline__ DT_ cuda_cos(DT_ val);

    /**
     * \brief Evaluates the max function.
     *
     * \tparam DT_ The datatype.
     *
     * \param[in] a The first value.
     * \param[in] b The second value
     *
     * \returns The max of a and b.
     */
    template<typename DT_>
    __host__ __device__ __forceinline__ DT_ cuda_max(DT_ a, DT_ b);

    /**
     * \brief Evaluates the min function.
     *
     * \tparam DT_ The datatype.
     *
     * \param[in] a The first value.
     * \param[in] b The second value
     *
     * \returns The min of a and b.
     */
    template<typename DT_>
    __host__ __device__ __forceinline__ DT_ cuda_min(DT_ a, DT_ b);

    template<>
    __host__ __device__ __forceinline__ float cuda_sin(float val)
    {
      return sinf(val);
    }

    template<>
    __host__ __device__ __forceinline__ double cuda_sin(double val)
    {
      return sin(val);
    }

    #ifdef FEAT_HAVE_HALFMATH
    template<>
    __host__ __device__ __forceinline__ Half cuda_sin(Half val)
    {
      #ifdef __CUDA_ARCH__
      return hsin(val);
      #else
      return __float2half(sinf(__half2float(val)));
      #endif
    }
    #endif

    template<>
    __host__ __device__ __forceinline__ float cuda_cos(float val)
    {
      return cosf(val);
    }

    template<>
    __host__ __device__ __forceinline__ double cuda_cos(double val)
    {
      return cos(val);
    }

    #ifdef FEAT_HAVE_HALFMATH
    template<>
    __host__ __device__ __forceinline__ Half cuda_cos(Half val)
    {
      #ifdef __CUDA_ARCH__
      return hcos(val);
      #else
      return __float2half(cosf(__half2float(val)));
      #endif
    }
    #endif



    template<>
    __host__ __device__ __forceinline__ double cuda_max(double a, double b)
    {
      return fmax(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ float cuda_max(float a, float b)
    {
      return fmaxf(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ int cuda_max(int a, int b)
    {
      return max(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ unsigned int cuda_max(unsigned int a, unsigned int b)
    {
      return max(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ long cuda_max(long a, long b)
    {
      return max(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ unsigned long cuda_max(unsigned long a, unsigned long b)
    {
      return max(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ long long cuda_max(long long a, long long b)
    {
      return max(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ unsigned long long cuda_max(unsigned long long a, unsigned long long b)
    {
      return max(a, b);
    }
    #ifdef FEAT_HAVE_HALFMATH
    template<>
    __host__ __device__ __forceinline__ Half cuda_max(Half a, Half b)
    {
      return __hmax(a, b);
    }
    #endif

    template<>
    __host__ __device__ __forceinline__ double cuda_min(double a, double b)
    {
      return fmin(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ float cuda_min(float a, float b)
    {
      return fminf(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ int cuda_min(int a, int b)
    {
      return min(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ unsigned int cuda_min(unsigned int a, unsigned int b)
    {
      return min(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ long cuda_min(long a, long b)
    {
      return min(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ unsigned long cuda_min(unsigned long a, unsigned long b)
    {
      return min(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ long long cuda_min(long long a, long long b)
    {
      return min(a, b);
    }

    template<>
    __host__ __device__ __forceinline__ unsigned long long cuda_min(unsigned long long a, unsigned long long b)
    {
      return min(a, b);
    }

    #ifdef FEAT_HAVE_HALFMATH
    template<>
    __host__ __device__ __forceinline__ Half cuda_min(Half a, Half b)
    {
      return __hmin(a, b);
    }
    #endif

    /// TODO: There should be a std::limits implementation for cuda
    template<typename DT_>
    __host__ __device__ __forceinline__ DT_ cuda_get_eps();

    template<>
    __host__ __device__ __forceinline__ double cuda_get_eps<double>()
    {
      return DBL_EPSILON;
    }

    template<>
    __host__ __device__ __forceinline__ float cuda_get_eps<float>()
    {
      return FLT_EPSILON;
    }

    #ifdef FEAT_HAVE_HALFMATH
    template<>
    __host__ __device__ __forceinline__ Half cuda_get_eps<Half>()
    {
      return CUDART_ONE_FP16 - Half(0.99951171);
    }
    #endif

    template<typename DT_>
    __host__ __device__ __forceinline__ DT_ cuda_get_sqrt_eps();

    template<>
    __host__ __device__ __forceinline__ double cuda_get_sqrt_eps<double>()
    {
      return cuda_sqrt(DBL_EPSILON);
    }

    template<>
    __host__ __device__ __forceinline__ float cuda_get_sqrt_eps<float>()
    {
      return cuda_sqrt(FLT_EPSILON);
    }

    #ifdef FEAT_HAVE_HALFMATH
    template<>
    __host__ __device__ __forceinline__ Half cuda_get_sqrt_eps<Half>()
    {
      return cuda_sqrt(CUDART_ONE_FP16 - Half(0.99951171));
    }
    #endif

    // We have to provide our own implementation of atomicMax (and Min) for floating point types,
    // since cuda only provides these atomic operations for integer types.
    // The implementation are based on the IEEM standard and use the fact, that the exponent bits
    // have higher presedent in ordering than the mantisse. This of course depends on the implementation,
    // i.e. the endianness of floating point in relation to the integers, so prey to the nvidia gods, that
    // they do not decide to do something strange.
    // For this reason, we maybe should write a test function....
    // see https://stackoverflow.com/questions/17399119/how-do-i-use-atomicmax-on-floating-point-values-in-cuda
    // for a bit more discussion
    /**
     * \brief Atomic maximum evaluation.
     *
     * Writes the max of *address and value in *address in a threadsafe manner, returning the old value previously contained in *address.
     *
     * \tparam DT_ The datatype to be used.
     *
     * \param[in/out] address Device pointer to the value that is to be changend.
     * \param[in] value The value *address should be compared to.
     *
     * \return The old value address pointed to.
     */
    template<typename DT_>
    __device__ __forceinline__ DT_ cuda_atomic_max(DT_* address, DT_ value);

    template<>
    __device__ __forceinline__ int cuda_atomic_max(int* address, int value)
    {
      return atomicMax(address, value);
    }

    template<>
    __device__ __forceinline__ unsigned int cuda_atomic_max(unsigned int* address, unsigned int value)
    {
      return atomicMax(address, value);
    }

    // no atomic max provided for long ... if explicitly required, do it by cas loop
    // template<>
    // __device__ __forceinline__ long cuda_atomic_max(long* address, long value)
    // {
    //   return atomicMax(address, value);
    // }

    // template<>
    // __device__ __forceinline__ unsigned long cuda_atomic_max(unsigned long* address, unsigned long value)
    // {
    //   return atomicMax(address, value);
    // }

    template<>
    __device__ __forceinline__ long long cuda_atomic_max(long long* address, long long value)
    {
      return atomicMax(address, value);
    }

    template<>
    __device__ __forceinline__ unsigned long long cuda_atomic_max(unsigned long long* address, unsigned long long value)
    {
      return atomicMax(address, value);
    }

    /// This implemenation *should* be faster, but appears to introduce a floating point error... Dont ask me...
    // template<>
    // __device__ __forceinline__ float cuda_atomic_max(float* address, float value)
    // {
    //   //signbit to also check for negativ zero... this is not NaN save!
    //   //a bit of black magic by simply reinterpreting binary represantion... base idea:
    //   // a =  0   1 0    101                 b = 0  0 1    110
    //   //    sign  exp   mantisse
    //   // if sign(b) = 0  -> standard integer order (since first byte is also sign bit)
    //   // if sign(b) = 1  -> Smaller is better if we interpret sign bit as largest value
    //   //this uses the internal intrinsics for typcasting, which do not exists for other datatypes
    //   float old = !signbit(value) ? __int_as_float(atomicMax((int*)address, __float_as_int(value))) :
    //       __uint_as_float(atomicMin((unsigned int*)address, __float_as_uint(value)));

    //   return old;
    // }

    template<>
    __device__ __forceinline__ float cuda_atomic_max(float* address, float value)
    {
      static_assert(sizeof(float) == sizeof(int), "Short does not have 16 bits");

      //reinterpret our address, due to atomicCas, we have to use unsigned longs longs
      unsigned int* address_short = (unsigned int*)address;
      //get current value
      unsigned int old = *address_short;
      unsigned int assumed;
      //now loop until old and assumed value is the same, in this case, atomicCAS overwrites with the max value
      //neet thing: while atomicCAS executes, address can not we overwritten, so this "should" be threadsafe
      do
      {
        assumed = old;
        old = atomicCAS(address_short, assumed, __float_as_uint(fmax(value, __uint_as_float(assumed))));
      }while(assumed != old);

      return __uint_as_float(old);

    }

    // template<>
    // __device__ __forceinline__ double cuda_atomic_max(double* address, double value)
    // {
    //   static_assert(sizeof(double) == sizeof(long long), "Double does not have 64 bits");

    //   //we have to use c style casts... this "should" do the right thing....
    //   if(!signbit(value))
    //   {
    //     auto tmp = atomicMax((long long*)address, *(long long*)&value);
    //     return *(double*)&tmp;
    //   }
    //   auto tmp2 = atomicMax((unsigned long long*)address, *(unsigned long long*)&value);
    //   return *(double*)&tmp2;
    // }

    template<>
    __device__ __forceinline__ double cuda_atomic_max(double* address, double value)
    {
      static_assert(sizeof(double) == sizeof(long long), "Short does not have 16 bits");

      //reinterpret our address, due to atomicCas, we have to use unsigned longs longs
      unsigned long long* address_short = (unsigned long long*)address;
      //get current value
      unsigned long long old = *address_short;
      unsigned long long assumed;
      //now loop until old and assumed value is the same, in this case, atomicCAS overwrites with the max value
      //neet thing: while atomicCAS executes, address can not we overwritten, so this "should" be threadsafe
      do
      {
        assumed = old;
        old = atomicCAS(address_short, assumed, __double_as_longlong(fmax(value, __longlong_as_double(assumed))));
      }while(assumed != old);

      return __longlong_as_double(old);

    }

    #ifdef FEAT_HAVE_HALFMATH
    //no atomicMax for short, so we have to implement the atomic operation with a CAS loop
    template<>
    __device__ __forceinline__ Half cuda_atomic_max(Half* address, Half value)
    {
      static_assert(sizeof(Half) == sizeof(unsigned short), "Short does not have 16 bits");

      //reinterpret our address
      unsigned short* address_short = (unsigned short*)address;
      //get current value
      unsigned short old = *address_short;
      unsigned short assumed;
      //now loop until old and assumed value is the same, in this case, atomicCAS overwrites with the max value
      //neet thing: while atomicCAS executes, address can not we overwritten, so this "should" be threadsafe
      do
      {
        assumed = old;
        old = atomicCAS(address_short, assumed, __hmax(value, *(Half*)&assumed));
      }while(assumed != old);

      return old;

    }
    #endif

    template<typename DT_, typename IT_>
    __host__ __device__ inline DT_ cuda_invert_matrix(const IT_ n, const IT_ stride, DT_* __restrict__ a, IT_* __restrict__ p)
    {
      #ifndef __CUDA_ARCH__
      // make sure that the parameters are valid
      if((n <= IT_(0)) || (stride < n) || (a == nullptr) || (p == nullptr))
        return DT_(0);
      #endif

      #ifndef __CUDA_ARCH__
      // invert 1x1 explicitly
      if(n == IT_(1))
      {
        DT_ det = a[0];
        a[0] = DT_(1) / det;
        return det;
      }
      #endif

      // initialize identity permutation
      for(IT_ i = 0; i < n; ++i)
      {
        p[i] = i;
      }

      // initialize determinant to 1
      DT_ det = DT_(1);

      // primary column elimination loop
      for(IT_ k(0); k < n; ++k)
      {
        // step 1: find a pivot for the elimination of column k
        {
          // for this, we only check the rows p[j] with j >= k, as all
          // rows p[j] with j < k have already been eliminated and are
          // therefore not candidates for pivoting
          DT_ pivot = CudaMath::cuda_abs(a[p[k]*stride + p[k]]);
          IT_ i = k;

          // loop over all unprocessed rows
          for(IT_ j(k+1); j < n; ++j)
          {
            // get our matrix value and check whether it can be a pivot
            DT_ abs_val = CudaMath::cuda_abs(a[p[j]*stride + p[j]]);
            if(abs_val > pivot)
            {
              pivot = abs_val;
              i = j;
            }
          }

          // do we have to swap rows i and k?
          if(i > k)
          {
            // swap rows "virtually" by exchanging their permutation positions
            IT_ t = p[k];
            p[k] = p[i];
            p[i] = t;
          }
        }

        // compute pivot row offset
        const IT_ pk_off = p[k]*stride;

        // step 2: process pivot row
        {
          // update determinant by multiplying with the pivot element
          det *= a[pk_off + p[k]];

          // get our inverted pivot element
          const DT_ pivot = DT_(1) / a[pk_off + p[k]];

          // replace column entry by unit column entry
          a[pk_off + p[k]] = DT_(1);

          // divide the whole row by the inverted pivot
          for(IT_ j(0); j < n; ++j)
          {
            a[pk_off+j] *= pivot;
          }
        }

        // step 3: eliminate pivot column

        // loop over all rows of the matrix
        for(IT_ i(0); i < n; ++i)
        {
          // skip the pivot row
          if(i == p[k])
            continue;

          // compute row and pivot offsets
          const IT_ row_off = i*stride;

          // fetch elimination factor
          const DT_ factor =  a[row_off + p[k]];

          // replace by unit column entry
          a[row_off + p[k]] = DT_(0);

          // process the row
          for(IT_ j(0); j < n; ++j)
          {
            a[row_off + j] -= a[pk_off + j] * factor;
          }
        }
      }

      // return determinant
      return det;
    }


    // only callable by a coalesced, i.e. a sub warp group...
    template<typename DT_, typename IT_>
    __device__ __forceinline__ DT_ cuda_grouped_invert_matrix(const cg::coalesced_group& tg, const int n, const int stride, DT_* __restrict__ a, IT_* __restrict__ p)
    {
      const int idx = tg.thread_rank();
      const int g_size = tg.num_threads();

      // initialize identity permutation
      for(int i = idx; i < n; i += g_size)
      {
        p[i] = i;
      }

      // initialize determinant to 1
      DT_ det = DT_(1);

      // all threads work on one column, the could be more leeways to optimize by red-balck strategies
      // we have to see...
      for(int k = 0; k < n; ++k)
      {
        // step 1: use a sub-group to find the maximum diag pivot
        {
          for(int i = idx; (i < n); i += g_size)
          {
            auto active_g = cg::coalesced_threads();
            int k_prob = CudaMath::cuda_max(i,k);
            k_prob = reduce(active_g, k_prob,
              [&](const int& v1, const int& v2) -> int
              {return cuda_abs(a[p[v1]*stride + p[v1]]) >= cuda_abs(a[p[v2]*stride + p[v2]]) ?
                                      v1 : v2;});
            cg::invoke_one(active_g, [&](){
              if((k_prob > k) && (cuda_abs(a[p[k_prob]*stride + p[k_prob]]) > cuda_abs(a[p[k]*stride + p[k]])))
                cuda::std::swap(p[k], p[k_prob]);});
            active_g.sync();
          }
        }
        // wait for all threads to arrive
        tg.sync();
        // compute pivot row offset
        const int pk_off = p[k]*stride;

        // step 2: process pivot
        {
          // update determinant by multiplying with the pivot element
          det *= a[pk_off + p[k]];

          // get our inverted pivot element
          const DT_ pivot = DT_(1) / a[pk_off + p[k]];

          // replace column entry by unit column entry
          cg::invoke_one(tg, [&](){
                              a[pk_off + p[k]] = DT_(1);});
          tg.sync();
          // divide the whole row by the inverted pivot
          for(int i = idx; i < n; i += g_size)
          {
            a[pk_off+i] *= pivot;
          }
        }

        // step 3: eliminate pivot column

        // loop over all rows of the matrix
        for(int i = idx; i < n*n; i += g_size)
        {
          const int row = i/n;
          const int col = i%n;
          // skip the pivot row
          if(row == p[k])
            continue;

          // compute row and pivot offsets
          const int row_off = row*stride;

          // fetch elimination factor
          const DT_ factor =  a[row_off + p[k]];

          // combine update and replace of unit column entry in one step
          a[row_off + col] = (col == p[k] ? DT_(0) : a[row_off + col]) - a[pk_off + col] * factor;
        }

      }

      // return determinant
      return det;



    }

    template<typename DT_>
    __host__ __device__ constexpr inline std::size_t cuda_get_align();

    template<>
    __host__ __device__ constexpr inline std::size_t cuda_get_align<double>()
    {
      return std::size_t(8u);
    }

    template<>
    __host__ __device__ constexpr inline std::size_t cuda_get_align<float>()
    {
      return std::size_t(4u);
    }

    #ifdef FEAT_HAVE_HALFMATH
    template<>
    __host__ __device__ constexpr inline std::size_t cuda_get_align<Half>()
    {
      return std::size_t(2u);
    }
    #endif

    template<>
    __host__ __device__ constexpr inline std::size_t cuda_get_align<std::int64_t>()
    {
      return std::size_t(8u);
    }

    template<>
    __host__ __device__ constexpr inline std::size_t cuda_get_align<std::uint64_t>()
    {
      return std::size_t(8u);
    }

    template<>
    __host__ __device__ constexpr inline std::size_t cuda_get_align<std::int32_t>()
    {
      return std::size_t(4u);
    }

    template<>
    __host__ __device__ constexpr inline std::size_t cuda_get_align<std::uint32_t>()
    {
      return std::size_t(4u);
    }

    template<>
    __host__ __device__ constexpr inline std::size_t cuda_get_align<std::int16_t>()
    {
      return std::size_t(2u);
    }

    template<>
    __host__ __device__ constexpr inline std::size_t cuda_get_align<std::uint16_t>()
    {
      return std::size_t(2u);
    }

    /// this function returns the maximum alignment of the datatypes to allow for reinterpret casts
    template<typename DT_, typename IT_>
    __host__ __device__ constexpr inline std::size_t cuda_get_min_align()
    {
      return (cuda_get_align<DT_>() < cuda_get_align<IT_>()) ? cuda_get_align<DT_>() : cuda_get_align<IT_>();
    }
  }
}
