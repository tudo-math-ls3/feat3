
// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_SOLVER_VOXEL_AMAVANKA_HPP
#define KERNEL_SOLVER_VOXEL_AMAVANKA_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/backend.hpp>
#include <kernel/adjacency/coloring.hpp>
#include <kernel/lafem/null_matrix.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/saddle_point_matrix.hpp>
#include <kernel/lafem/tuple_matrix.hpp>
#include <kernel/global/matrix.hpp>
#include <kernel/solver/base.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/solver/amavanka.hpp>

#ifdef FEAT_HAVE_CUDA
#include <kernel/util/cuda_util.hpp>
#endif

namespace FEAT
{
  namespace Intern
  {
    /**
     * \brief Controls threading strategy of vanka assembly.
     *
     * \note For now, only controls the strategy on device side.
     */
    enum VankaAssemblyPolicy
    {
      /// Use one thread per macro
      oneThreadperBlock = 0,
      /// Use batched cublas kernels
      batchedAssembly = 1
    };

    /**
     * \brief Set the expected macro type.
     *
     * \attention Choosing uniform Macros if the macros are actually anisotropic
     *  leads to wrong results.
     *  While the other way round always works, uniform macros provide many
     *  performance benefits, so choose wisely.
     */
    enum VankaMacroPolicy
    {
      /// Macros can have different sizes
      anisotropicMacros = 0,
      /// All macros have the same size
      uniformMacros = 1
    };
  }
  namespace Solver
  {
    namespace Arch
    {
      /**
       * \brief Assembles vanka matrix on host.
       *
       * \tparam DT_ The datatype to be used.
       * \tparam IT_ The indextype to be used.
       * \tparam n_ The meta row/column size of the blocked matrix.
       *
       * This kernel assembles the vanka matrix on host side in a threadparallel fashion based
       * on a provided macro coloring. Paralleization done with OpenMP.
       *
       * \param[in] mat_wrap
       *  Wrapper object for the system matrix.
       *
       * \param[in/out] vanka_wrap
       *  Wrapper object for the preallocated vanka matrix.
       *
       * \param[in] macro_dofs
       *  Vector of Graphs mapping macros to dofs for each meta-block.
       *
       * \param[in] dof_macros
       *  Vector of Graphs mapping dofs to macros for each meta-block.
       *
       * \param[in/out] macro_mask
       *  Maskarray for singular macros. Ignored if not skip_singular.
       *
       * \param[in] coloring_data
       *  Handler for coloring data. Has to conform to the matrix objects.
       *
       * \param[in] stride
       *  The stride (i.e. actual data size) of the local matrix row.
       *
       * \param[in] omega
       *  Scaling factor for the vanka matrix.
       *
       * \param[in] eps
       *  Tolerance for singular matrix evaluation.
       *
       * \param[in] skip_singular
       *  Should singular local matrices be skipped and replaced by unit matrices?
       *
       * \author Maximilian Esser
       */
      template<typename DT_, typename IT_, int n_>
      void assemble_vanka_host(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
        Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& vanka_wrap, const std::vector<Adjacency::Graph>& macro_dofs,
        const std::vector<Adjacency::Graph>& dof_macros, std::vector<int>& macro_mask, const Adjacency::ColoringDataHandler& coloring_data,
        Index stride, DT_ omega, DT_ eps, bool skip_singular);

      /**
       * \brief Assembles vanka matrix on device.
       *
       * \tparam DT_ The datatype to be used.
       * \tparam IT_ The indextype to be used.
       * \tparam n_ The meta row/column size of the blocked matrix.
       *
       * This kernel assembles the vanka matrix on device side in a "one Macro per Thread" fashion based
       * on a provided macro coloring.
       *
       * \attention Only for testing purposes. Use the batched variant for good performance.
       *
       * \param[in] mat_wrap
       *  Wrapper object for the system matrix.
       *
       * \param[in/out] vanka_wrap
       *  Wrapper object for the preallocated vanka matrix.
       *
       * \param[in] d_macro_dofs
       *  Vector of device pointers of graphs mapping macros to dofs for each meta-block. See below
       *  for specification.
       *
       * \param[in] d_dof_macros
       *  Vector of device pointers of graphs mapping dofs to macros for each meta-block.
       *
       * \param[in/out] d_macro_mask
       *  Device ptr to maskarray for singular macros. Ignored if not skip_singular.
       *
       * \param[in] max_degree_dofs
       *  Vector of maximum degrees of the macro to dofs graph.
       *
       * \param[in] max_degree_macros
       *  Vector of maximum degrees of the dof to macros graph.
       *
       * \param[in] coloring_data
       *  Handler for coloring data. Has to conform to the matrix objects.
       *
       * \param[in] num_macros
       *  The number of macros used to assemble the local vanka blocks.
       *
       * \param[in] stride
       *  The stride (i.e. actual data size) of the local matrix rows.
       *
       * \param[in] omega
       *  Scaling factor for the vanka matrix.
       *
       * \param[in] eps
       *  Tolerance for singular matrix evaluation.
       *
       * \param[in] skip_singular
       *  Should singular local matrices be skipped and replaced by unit matrices?
       *
       * \param[in] uniform_macros
       *  Are the macros uniform in size? Based on this, compile-time optimizations are chosen.
       *
       *
       *  The specific format of d_macro_dofs depends on uniform_macros. If macros are uniform,
       *  for each pointer the very first entry holds the size of the uniform macro.
       *  Then there are num_macro times macro_size entries, each macro_sized subarray holding
       *  the indices of the mapped dofs:
       *  d_macro_dofs[i] -> [m_size, 0, 5, 29, 89, 1, 6, 21, 70, 4, .....]
       *                             <-----m------> <-----m-----> <--m-->   <- num_macros times
       *
       *  In the case of non-uniform macros and always for d_dof_macros, the i-th ptr maps to an array
       *  of size num_macros times max_degree_dofs[i]+1 entries, where each first entry holds
       *  the size of the local map, followed by the respective dof mapping, padded with ~Index(0)=:e
       *  to max_degree_dofs[i] entries:
       *  d_macros_dofs[i] -> [3, 0, 5, 14, e, 4, 1, 6, 19, 42, 2, 2, 6, e, e, 4,....]
       *                      n0  <--n0-->     n1 <----n1---->  n2 <n2->       n3         <- num_macros times
       *                          <---max_d-->    <---max_d-->    <---max_d-->
       *
       * \author Maximilian Esser
       */
      template<typename DT_, typename IT_, int n_>
      void assemble_vanka_device(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
        Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& vanka_wrap, const std::vector<Index*>& d_macro_dofs,
        const std::vector<Index*>& d_dof_macros, int* d_macro_mask, const std::vector<Index>& max_degree_dofs,
        const std::vector<Index>& max_degree_macros, const Adjacency::ColoringDataHandler& coloring_data,
        Index num_macros, Index stride, DT_ omega, DT_ eps, bool skip_singular, bool uniform_macros);

      /**
       * \brief Assembles vanka matrix on device.
       *
       * \tparam DT_ The datatype to be used.
       * \tparam IT_ The indextype to be used.
       * \tparam n_ The meta row/column size of the blocked matrix.
       *
       * This kernel assembles the vanka matrix on device side with the help of batchwise cublas kernels.
       * The macros-wise scattering is done based on a provided macro coloring.
       *
       * \param[in] mat_wrap
       *  Wrapper object for the system matrix.
       *
       * \param[in/out] vanka_wrap
       *  Wrapper object for the preallocated vanka matrix.
       *
       * \param[in] d_macro_dofs
       *  Vector of device pointers of graphs mapping macros to dofs for each meta-block. See below
       *  for specification.
       *
       * \param[in] d_dof_macros
       *  Vector of device pointers of graphs mapping dofs to macros for each meta-block.
       *
       * \param[in/out] d_macro_mask
       *  Device ptr to maskarray for singular macros. Ignored if not skip_singular.
       *
       * \param[in] max_degree_dofs
       *  Vector of maximum degrees of the macro to dofs graph.
       *
       * \param[in] max_degree_macros
       *  Vector of maximum degrees of the dof to macros graph.
       *
       * \param[in] coloring_data
       *  Handler for coloring data. Has to conform to the matrix objects.
       *
       * \param[in] num_macros
       *  The number of macros used to assemble the local vanka blocks.
       *
       * \param[in] stride
       *  The stride (i.e. actual data size) of the local matrix rows.
       *
       * \param[in] omega
       *  Scaling factor for the vanka matrix.
       *
       * \param[in] eps
       *  Tolerance for singular matrix evaluation.
       *
       * \param[in] skip_singular
       *  Should singular local matrices be skipped and replaced by unit matrices?
       *
       *  \attention For now, this kernel can only be used with uniform macros.
       *
       *  For d_macro_dofs the following has to hold:
       *  For each pointer the very first entry holds the size of the uniform macro.
       *  Then there are num_macro times macro_size entries, each macro_sized subarray holding
       *  the indices of the mapped dofs:
       *  d_macro_dofs[i] -> [m_size, 0, 5, 29, 89, 1, 6, 21, 70, 4, .....]
       *                             <-----m------> <-----m-----> <--m-->   <- num_macros times
       *
       *  For d_dof_macros, the i-th ptr maps to an array
       *  of size num_macros times max_degree_dofs[i]+1 entries, where each first entry holds
       *  the size of the local map, followed by the respective dof mapping, padded with ~Index(0)=:e
       *  to max_degree_dofs[i] entries:
       *  d_macros_dofs[i] -> [3, 0, 5, 14, e, 4, 1, 6, 19, 42, 2, 2, 6, e, e, 4,....]
       *                      n0  <--n0-->     n1 <----n1---->  n2 <n2->       n3         <- num_macros times
       *                          <---max_d-->    <---max_d-->    <---max_d-->
       *
       * \author Maximilian Esser
       */
      template<typename DT_, typename IT_, int n_>
      void assemble_vanka_device_batched(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
        Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& vanka_wrap, const std::vector<Index*>& d_macro_dofs,
        const std::vector<Index*>& d_dof_macros, int* d_macro_mask, const std::vector<Index>& max_degree_dofs,
        const std::vector<Index>& max_degree_macros, const Adjacency::ColoringDataHandler& coloring_data,
        Index num_macros, Index stride, Index actual_matrix_size, DT_ omega, bool skip_singular);
    }

  #ifndef __CUDACC__
    /**
     * \brief Additive Macro-wise Matrix-based Vanka preconditioner/smoother
     *
     * \tparam Matrix_ The MatrixType of the system matrix.
     * \tparam Filter_ The FilterType to be used.
     * \tparam pol_threading_ The threading policy to be used. Defaults to batched assembly.
     * \tparam macro_type_ The macro policy to be used. Defaults to uniform macros.
     *
     * This class implements an additive macro-wise Vanka smoother, which stores its
     * pre-computed operator as a sparse matrix, so that each application of the Vanka
     * smoother consists of only one sparse matrix-vector multiplication.
     *
     * This class is based on the standard Amavanka smoother and differs in one major
     * point:
     *  The backend of the numeric assembly of the underlying vanka matrix can be chosen
     *  at runtime by the PreferredBackend static variable and supports (for now) a,
     *  on OpenMP based, threadparallel generic assembly and a CUDA based device assembly.
     *
     * This class supports all matrix-types and macro distributions that the AmaVanka baseclass
     * does. For specific information of the AmaVanka solver refer to the baseclass documentation.
     *
     *
     * \attention
     * Due to the threadparallel nature of the assembly, this class requires a valid macro coloring
     * to be added either to the constructor or added by the fill_color() method before the
     * symbolic initialization is started.
     *
     *
     * \author Maximilian Esser
     */
    template<typename Matrix_,
    typename Filter_,
    FEAT::Intern::VankaAssemblyPolicy pol_threading_ = FEAT::Intern::VankaAssemblyPolicy::batchedAssembly,
    FEAT::Intern::VankaMacroPolicy macro_type_ = FEAT::Intern::VankaMacroPolicy::uniformMacros>
    class VoxelAmaVanka :
      public Solver::AmaVanka<Matrix_, Filter_>
    {
    public:
      /// our base-class
      typedef Solver::AmaVanka<Matrix_, Filter_> BaseClass;

      /// our data type
      typedef typename Matrix_::DataType DataType;
      /// our index type
      typedef typename Matrix_::IndexType IndexType;
      /// our vector type
      typedef typename Matrix_::VectorTypeL VectorType;

    protected:
      /// the type of our Vanka matrix
      typedef typename Intern::AmaVankaMatrixHelper<Matrix_>::VankaMatrix VankaMatrixType;
      /// our matrix data wrapper
      typedef Intern::CSRTupleMatrixWrapper<DataType, IndexType, Intern::AmaVankaMatrixHelper<Matrix_>::num_blocks> MatrixWrapper;
      /// our vanka data wrapper
      typedef Intern::CSRTupleMatrixWrapper<DataType, IndexType, Intern::AmaVankaMatrixHelper<VankaMatrixType>::num_blocks> VankaWrapper;
      /// coloring
      Adjacency::ColoringDataHandler _coloring_data;
      /// vector of graph arrays
      std::vector<Index*> _d_macro_dofs, _d_dof_macros;
      /// size data
      std::vector<Index> _max_degree_dofs, _max_degree_macros;
      /// array of macro mask
      int* _d_macro_mask;
      /// flag whether we should allocate additional device pointer
      bool _allocate_device = false;

      /// \brief Calculate the max degree of our graphs
      void _alloc_max_degrees()
      {
        _max_degree_dofs.resize(this->_macro_dofs.size());
        _max_degree_macros.resize(this->_macro_dofs.size());
        for(std::size_t i = 0; i < _max_degree_dofs.size(); ++i)
        {
          if constexpr(macro_type_ == FEAT::Intern::VankaMacroPolicy::uniformMacros)
            _max_degree_dofs[i] = this->_macro_dofs[i].degree(Index(0));
          else
            _max_degree_dofs[i] = this->_macro_dofs[i].degree();
          _max_degree_macros[i] = this->_dof_macros[i].degree();
        }
      }

      /// \brief Allocates device pointers, if required.
      void _alloc_device()
      {
        XASSERTM(_max_degree_dofs.size() == this->_macro_dofs.size(), "call _alloc_max_degrees beforehand");
        if(!_allocate_device)
          return;
      #ifdef FEAT_HAVE_CUDA
        _d_macro_dofs.resize(this->_macro_dofs.size());
        _d_dof_macros.resize(this->_dof_macros.size());
        for(int i = 0; i < int(this->_macro_dofs.size()); ++i)
        {
          Index malloc_size;
          if constexpr(macro_type_ == FEAT::Intern::VankaMacroPolicy::uniformMacros)
            malloc_size = this->_macro_dofs[i].get_num_nodes_domain() * _max_degree_dofs[i] * sizeof(Index);
          else
            malloc_size = this->_macro_dofs[i].get_num_nodes_domain() * (_max_degree_dofs[i]+1) * sizeof(Index);
          _d_macro_dofs[i] = (Index*)Util::cuda_malloc_managed(malloc_size);
          // prepare tmp array
          {
            Index* tmp_alias = _d_macro_dofs[i];
            const Index* dom_ptr = this->_macro_dofs[i].get_domain_ptr();
            const Index* img_ptr = this->_macro_dofs[i].get_image_idx();
            for(int k = 0; k < int(this->_macro_dofs[i].get_num_nodes_domain()); ++k)
            {
              if constexpr(macro_type_ == FEAT::Intern::VankaMacroPolicy::uniformMacros)
              {
                std::memcpy(tmp_alias + k*_max_degree_dofs[i], img_ptr + dom_ptr[k], _max_degree_dofs[i]*sizeof(Index));
              }
              else
              {
                const Index loc_size = dom_ptr[k+1] - dom_ptr[k];
                tmp_alias[k*(_max_degree_dofs[i]+1)] = loc_size;
                std::memcpy(tmp_alias + k*(_max_degree_dofs[i]+1) + 1, img_ptr + dom_ptr[k], loc_size*sizeof(Index));
                std::memset(tmp_alias + k*(_max_degree_dofs[i]+1) + loc_size + 1, ~int(0), (_max_degree_dofs[i] - loc_size)*sizeof(Index));
              }
            }
          }
          malloc_size = this->_dof_macros[i].get_num_nodes_domain()*(_max_degree_macros[i]+1)*sizeof(Index);
          _d_dof_macros[i] = (Index*)Util::cuda_malloc_managed(malloc_size);
          {
            Index* tmp_alias = _d_dof_macros[i];
            const Index* dom_ptr = this->_dof_macros[i].get_domain_ptr();
            const Index* img_ptr = this->_dof_macros[i].get_image_idx();
            for(int k = 0; k < int(this->_dof_macros[i].get_num_nodes_domain()); ++k)
            {
              const Index loc_size = dom_ptr[k+1] - dom_ptr[k];
              tmp_alias[k*(_max_degree_macros[i]+1)] = loc_size;
              std::memcpy(tmp_alias + k*(_max_degree_macros[i]+1) + 1, img_ptr + dom_ptr[k], loc_size*sizeof(Index));
              std::memset(tmp_alias + k*(_max_degree_macros[i]+1) + loc_size + 1, ~int(0), (_max_degree_macros[i] - loc_size)*sizeof(Index));
            }
          }
        }
        {
          if(this->_skip_singular)
          {
            _d_macro_mask = (int*)Util::cuda_malloc_managed(this->_macro_mask.size()*sizeof(int));
            Util::cuda_set_memory(_d_macro_mask, 0, this->_macro_mask.size());
          }
        }
      #endif
      }

      /// \brief Frees device pointers
      void _free_device()
      {
      #ifdef FEAT_HAVE_CUDA
        Util::cuda_free(_d_macro_mask);
        _d_macro_mask = nullptr;
        for(int i = 0; i < int(_d_macro_dofs.size()); ++i)
        {
          Util::cuda_free((void*)(_d_macro_dofs[i]));
          Util::cuda_free((void*)(_d_dof_macros[i]));
        }
        _d_dof_macros.clear();
        _d_macro_dofs.clear();
      #endif
      }

      /**
       * \brief Calls generic numeric kernel
       *
       * \param[in] mat_wrap
       *  Wrapper of the system matrix.
       *
       * \param[in/out] vanka_wrap
       *  Wrapper of the vanka matrix.
       *
       * \param[in] num_macros
       *  Number of macros. Ignored.
       *
       * \param[in] stride
       *  The local matrix stride.
       *
       * \param[in] eps
       *  Tolerance for singular matrix identifaication.
       */
      void _init_numeric_generic(const MatrixWrapper& mat_wrap, VankaWrapper& vanka_wrap, Index DOXY(num_macros), Index stride, DataType eps)
      {
        //call backend
        Arch::assemble_vanka_host(mat_wrap, vanka_wrap, this->_macro_dofs, this->_dof_macros, this->_macro_mask, _coloring_data,
                            stride, this->_omega, eps, this->_skip_singular);
      }

      #if defined(FEAT_HAVE_CUDA) || defined(DOXYGEN)
      /**
       * \brief Calls cuda numeric kernel
       *
       * \param[in] mat_wrap
       *  Wrapper of the system matrix.
       *
       * \param[in/out] vanka_wrap
       *  Wrapper of the vanka matrix.
       *
       * \param[in] num_macros
       *  Number of macros. Ignored.
       *
       * \param[in] stride
       *  The local matrix stride.
       *
       * \param[in] eps
       *  Tolerance for singular matrix identifaication.
       */
      void _init_numeric_cuda(const MatrixWrapper& mat_wrap, VankaWrapper& vanka_wrap, Index num_macros, Index stride, DataType eps)
      {
        XASSERTM(_allocate_device, "Allocate device disabled!");
        bool uniform_macros = (macro_type_ == FEAT::Intern::VankaMacroPolicy::uniformMacros);
        //call backend
        if constexpr(pol_threading_ == FEAT::Intern::VankaAssemblyPolicy::oneThreadperBlock)
        {
          Arch::assemble_vanka_device(mat_wrap, vanka_wrap, _d_macro_dofs, _d_dof_macros, _d_macro_mask, _max_degree_dofs, _max_degree_macros, _coloring_data, num_macros, stride, this->_omega, eps, this->_skip_singular, uniform_macros);
        }
        else if (pol_threading_ == FEAT::Intern::VankaAssemblyPolicy::batchedAssembly)
        {
          XASSERTM(uniform_macros, "Batched assembly only works with uniform macros!");
          // TODO: stride always actual local matrix size?
          Arch::assemble_vanka_device_batched(mat_wrap, vanka_wrap, _d_macro_dofs, _d_dof_macros, _d_macro_mask, _max_degree_dofs, _max_degree_macros, _coloring_data, num_macros, stride, stride, this->_omega, this->_skip_singular);
        }

      }
      #endif


    public:
      /**
       * \brief Constructor
       *
       * \tparam ColoringType_
       * The type of the coloring array, should either be a vector of ints or a Coloring object.
       *
       * \param[in] matrix
       * The saddle-point system matrix.
       *
       * \param[in] filter
       * The system filter.
       *
       * \param[in] coloring
       * The coloring of the macros to be added. This should be used, if automatic marco deduction is used.
       *
       * \param[in] omega
       * The damping parameter.
       *
       * \param[in] num_steps
       * The number of smoothing steps to be performed.
       */
      template<typename ColoringType_>
      explicit VoxelAmaVanka(const Matrix_& matrix, const Filter_& filter,
        const ColoringType_& coloring,
        const DataType omega = DataType(1), const Index num_steps = Index(1)) :
        BaseClass(matrix, filter, omega, num_steps),
        _coloring_data(),
        _d_macro_dofs(),
        _d_dof_macros(),
        _d_macro_mask(nullptr)
      {
        _coloring_data.fill_color(coloring);
        #ifdef FEAT_HAVE_CUDA
        _allocate_device = Util::cuda_get_device_count() > 0;
        #endif
      }

      // /**
      //  * \brief Sets whether device ptr should be allocated.
      //  *
      //  * \param[in] allocate Should device ptr be allocated?
      //  *
      //  * \warning Setting this to false while having PreferedBackend set to cuda during the assembly
      //  *          will lead to an error.
      //  */
      // void set_allocate_device(bool allocate)
      // {
      //   _allocate_device = allocate;
      // }

      /**
       * \brief Fills the coloring data.
       *
       * \tparam ColoringType_ Arraytype mapping index to a color
       *
       * \param[in] color The coloring data. Has to fit the macro dofs, i.e. call after pushing all macro dofs.
       * \param[in] hint Optionally give hint on the number of colors.
       *
       */
      template<typename ColoringType_>
      void fill_color(const ColoringType_& color, int hint = -1)
      {
        XASSERTM(color.size() == this->_macro_dofs.front().get_num_nodes_domain(), "Coloring does not fit macro dofs");
        XASSERTM(_coloring_data.initialized(), "Coloring data already initialized");
        _coloring_data.fill_color(color, hint);
      }

      /// \brief Returns the name of the solver.
      virtual String name() const override
      {
        return "VoxelAmaVanka";
      }

      /// \brief Initializes symbolic values and device pointers.
      virtual void init_symbolic() override
      {
        BaseClass::init_symbolic();
        this->watch_init_symbolic.start();
        _alloc_max_degrees();
        // _alloc_row_helper();
        _alloc_device();
        this->watch_init_symbolic.stop();
      }

      /// \brief Frees symbolic values and device pointers.
      virtual void done_symbolic() override
      {
        _free_device();
        // _accum_row_ctr.clear();
        // _accum_row_index.clear();
        _max_degree_dofs.clear();
        _max_degree_macros.clear();
        BaseClass::done_symbolic();
      }


      virtual void init_numeric() override
      {
        const DataType eps = Math::eps<DataType>();
        //call numeric init of BaseSolver, but not of Vanka
        this->watch_init_numeric.start();
        FEAT::Solver::SolverBase<typename Matrix_::VectorTypeL>::init_numeric();

        // get maximum macro size
        const Index num_macros = Index(this->_macro_dofs.front().get_num_nodes_domain());
        const Index stride = Intern::AmaVankaCore::calc_stride(this->_vanka, this->_macro_dofs);
        this->_vanka.format();
        //gather matrix wrappers
        auto matrix_wrapper = Solver::Intern::get_meta_matrix_wrapper(this->_matrix);
        // std::cout << matrix_wrapper.print();
        auto vanka_wrapper = Solver::Intern::get_meta_matrix_wrapper(this->_vanka);
        BACKEND_SKELETON_VOID(_init_numeric_cuda, _init_numeric_generic, _init_numeric_generic, matrix_wrapper, vanka_wrapper, num_macros, stride, eps)

        this->watch_init_numeric.stop();
      }


    }; // VoxelAmaVanka


    /**
     * \brief Creates a new VoxelAmaVanka smoother object
     *
     * \param[in] matrix
     * The system matrix.
     *
     * \param[in] filter
     * The system filter
     *
     * \param[in] coloring
     * The coloring of the voxel array.
     *
     * \param[in] omega
     * The damping parameter.
     *
     * \param[in] num_steps
     * The number of Vanka iterations to be performed.
     *
     * \returns
     * A shared pointer to a new AmaVanka object.
     */
    template<typename Matrix_, typename Filter_, typename ColoringType_,
    FEAT::Intern::VankaAssemblyPolicy pol_threading_ = FEAT::Intern::VankaAssemblyPolicy::batchedAssembly,
    FEAT::Intern::VankaMacroPolicy macro_type_ = FEAT::Intern::VankaMacroPolicy::uniformMacros>
    std::shared_ptr<VoxelAmaVanka<Matrix_, Filter_, pol_threading_, macro_type_>> new_voxel_amavanka(const Matrix_& matrix, const Filter_& filter, const ColoringType_& coloring,
                                                                        typename Matrix_::DataType omega = typename Matrix_::DataType(1), Index num_steps = Index(1))
    {
      return std::make_shared<VoxelAmaVanka<Matrix_, Filter_, pol_threading_, macro_type_>>(matrix, filter,coloring, omega, num_steps);
    }
  #endif //__CUDACC__

  }
}

#endif // KERNEL_SOLVER_VOXEL_AMAVANKA_HPP