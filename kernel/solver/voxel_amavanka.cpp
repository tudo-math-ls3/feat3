#include <kernel/base_header.hpp>
#include <kernel/adjacency/coloring.hpp>
#include <kernel/solver/amavanka_base.hpp>
#include <kernel/solver/voxel_amavanka.hpp>
#include <vector>

namespace FEAT
{
  namespace Solver
  {
    namespace Kernel
    {
      template<typename DT_, typename IT_, int n_, bool skip_singular_>
      void assemble_unscaled_vanka_host(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
        Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& vanka_wrap,
        const Adjacency::Graph** macro_dofs, int* macro_mask,
        const int* coloring_map, Index color_size,
        Index stride, DT_ eps)
      {
        typedef DT_ DataType;
        // typedef IT_ IndexType;
        if constexpr(skip_singular_)
        {
          #pragma omp parallel
          {
            // allocate arrays for local matrix
            std::unique_ptr<DataType[]> local(new DataType[stride*stride]);
            std::fill(local.get(), local.get() + stride*stride, DataType(0));
            std::unique_ptr<DataType[]> local_t(new DataType[stride*stride]);
            std::fill(local_t.get(), local_t.get() + stride*stride, DataType(0));
            std::unique_ptr<Index[]> pivot(new Index[stride]);
            std::fill(pivot.get(), pivot.get() + stride, Index(0));
            #pragma omp for
            for(int k = 0; k < int(color_size); ++k)
            {
              const int imacro = coloring_map[k];
              const std::pair<Index, Index> nrc = Intern::VoxelAmaVankaCore::gather(mat_wrap, local.get(), stride, Index(imacro), macro_dofs,
                                                                            Index(0), Index(0), Index(0), Index(0));
              // the approach used for checking the regularity of the local matrix is to check whether
              //
              //     || I - A*A^{-1} ||_F^2 < eps
              //
              // we could try to analyse the pivots returned by invert_matrix function instead, but
              // unfortunately this approach sometimes leads to false positives

              // make a backup if checking for singularity
              for(Index i(0); i < nrc.first; ++i)
                for(Index j(0); j < nrc.second; ++j)
                  local_t[i*stride+j] = local[i*stride+j];

              // invert matrix
              Math::invert_matrix(nrc.first, stride, local.get(), pivot.get());

              // compute (squared) Frobenius norm of (I - A*A^{-1})
              DataType norm = DataType(0);
              for(Index i(0); i < nrc.first; ++i)
              {
                for(Index j(0); j < nrc.first; ++j)
                {
                  DataType xij = DataType(i == j ? 1 : 0);
                  for(Index l(0); l < nrc.first; ++l)
                    xij -= local_t[i*stride+l] * local[l*stride+j]; // A_il * (A^{-1})_lj
                  norm += xij * xij;
                }
              }

              // is the matrix block singular?
              // Note: we check for !(norm < eps) instead of (norm >= eps),
              // because the latter one evaluates to false if norm is NaN,
              // which would result in a false negative
              const bool singular = !(norm < eps);

              // set macro regularity mask
              macro_mask[imacro] = (singular ? 0 : 1);

              // scatter local matrix
              if(!singular)
              {
                Intern::VoxelAmaVankaCore::scatter_add(vanka_wrap, local.get(), stride, Index(imacro), macro_dofs,
                  Index(0), Index(0), Index(0), Index(0));
              }

              // reformat local matrix TODO: necessary?
              for(Index i(0); i < nrc.first; ++i)
                for(Index j(0); j < nrc.second; ++j)
                  local[i*stride+j] = DataType(0);
            }
          }
        }
        else
        {
          #pragma omp parallel
          {
            // #pragma master
            // {
            //   std::cout << "Num threads " << omp_get_num_threads() << "\n";
            // }
            // allocate arrays for local matrix
            std::unique_ptr<DataType[]> local(new DataType[stride*stride]);
            std::fill(local.get(), local.get() + stride*stride, DataType(0));
            std::unique_ptr<Index[]> pivot(new Index[stride]);
            std::fill(pivot.get(), pivot.get() + stride, Index(0));
            #pragma omp for
            for(int k = 0; k < int(color_size); ++k)
            {
              const int imacro = coloring_map[k];
              const std::pair<Index, Index> nrc = Intern::VoxelAmaVankaCore::gather(mat_wrap, local.get(), stride, Index(imacro), macro_dofs,
                                                                            Index(0), Index(0), Index(0), Index(0));

              // invert matrix
              Math::invert_matrix(nrc.first, stride, local.get(), pivot.get());

              Intern::VoxelAmaVankaCore::scatter_add(vanka_wrap, local.get(), stride, Index(imacro), macro_dofs,
                Index(0), Index(0), Index(0), Index(0));

              // reformat local matrix TODO: necessary?
              for(Index i(0); i < nrc.first; ++i)
                for(Index j(0); j < nrc.second; ++j)
                  local[i*stride+j] = DataType(0);
            }
          }
        }
      }

      template<typename DT_, typename  IT_, int n_, bool skip_singular_>
      void scale_vanka_rows_host(Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& vanka_wrap, const DT_ omega,
                                  const Adjacency::Graph** graphs_dof_macros, const int* _m_mask)
      {
        for(int i = 0; i < n_; ++i)
        {
          const Index* row_dom_ptr = graphs_dof_macros[i]->get_domain_ptr();
          const Index* row_img_idx = graphs_dof_macros[i]->get_image_idx();
          const int hw = vanka_wrap.blocksizes[Index(1) + Math::min(Index(i),Index(1))];
          const int num_rows = int(vanka_wrap.tensor_counts[2*i]);
          for(int j = 0; j < n_; ++j)
          {
            DT_* vals = vanka_wrap.data_arrays[i*n_ + j];
            const IT_* row_ptr = vanka_wrap.row_arrays[i*n_+j];
            const IT_* col_idx = vanka_wrap.col_arrays[i*n_+j];
            const int hb = vanka_wrap.blocksizes[j +1];
            // now do the actual inner openmp loop over each column of the sub matrix
            #pragma omp parallel for
            for(int row = 0; row < num_rows; ++row)
            {
              // careful, num rows is counted against native elements, not raw elements
              Intern::VoxelAmaVankaCore::scale_row<DT_, IT_, skip_singular_>(vals, omega, row_ptr, col_idx, row_dom_ptr, row_img_idx, hw, hb, Index(row), _m_mask);
            }
          }
        }

      }
    }

    namespace Arch
    {
      template<typename DT_, typename IT_, int n_>
      void assemble_vanka_host(const Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& mat_wrap,
        Intern::CSRTupleMatrixWrapper<DT_, IT_, n_>& vanka_wrap, const std::vector<Adjacency::Graph>& macro_dofs,
        const std::vector<Adjacency::Graph>& dof_macros, std::vector<int>& macro_mask, const Adjacency::ColoringDataHandler& coloring_data,
        Index stride, DT_ omega, DT_ eps, bool skip_singular)
      {
        //further extract internal arrays
        const Adjacency::Graph* graphs_dof_macros[n_];
        const Adjacency::Graph* graphs_macro_dofs[n_];
        for(int i = 0; i < n_; ++i)
        {
          graphs_dof_macros[i] = &(dof_macros.at(std::size_t(i)));
          graphs_macro_dofs[i] = &(macro_dofs.at(std::size_t(i)));
        }
        int* _m_mask = macro_mask.data();
        if(!skip_singular)
          _m_mask = nullptr;
        auto& coloring_map = coloring_data.get_coloring_maps();
        for(Index k = 0; k < coloring_map.size(); ++k)
        {
          if(skip_singular)
            Solver::Kernel::template assemble_unscaled_vanka_host<DT_, IT_, n_, true>
                                            (mat_wrap, vanka_wrap, graphs_macro_dofs,
                                             _m_mask, coloring_map[k].data(), coloring_map[k].size(), stride, eps);
          else
            Solver::Kernel::template assemble_unscaled_vanka_host<DT_, IT_, n_, false>
                                            (mat_wrap, vanka_wrap, graphs_macro_dofs,
                                             _m_mask, coloring_map[k].data(), coloring_map[k].size(), stride, eps);
        }
        if(skip_singular)
          Solver::Kernel::template scale_vanka_rows_host<DT_, IT_, n_, true>(vanka_wrap, omega, graphs_dof_macros,_m_mask);
        else
          Solver::Kernel::template scale_vanka_rows_host<DT_, IT_, n_, false>(vanka_wrap, omega, graphs_dof_macros, _m_mask);

      }
    }

  }
}

//instantiate the templates
using namespace FEAT;
using namespace FEAT::Solver;


template void Arch::assemble_vanka_host(const Intern::CSRTupleMatrixWrapper<double, std::uint64_t, 2>&, Intern::CSRTupleMatrixWrapper<double, std::uint64_t, 2>&, const std::vector<Adjacency::Graph>&,
                                        const std::vector<Adjacency::Graph>&, std::vector<int>&, const Adjacency::ColoringDataHandler&, Index, double, double, bool);
template void Arch::assemble_vanka_host(const Intern::CSRTupleMatrixWrapper<double, std::uint32_t, 2>&, Intern::CSRTupleMatrixWrapper<double, std::uint32_t, 2>&, const std::vector<Adjacency::Graph>&,
                                        const std::vector<Adjacency::Graph>&, std::vector<int>&, const Adjacency::ColoringDataHandler&, Index, double, double, bool);
template void Arch::assemble_vanka_host(const Intern::CSRTupleMatrixWrapper<float, std::uint64_t, 2>&, Intern::CSRTupleMatrixWrapper<float, std::uint64_t, 2>&, const std::vector<Adjacency::Graph>&,
                                        const std::vector<Adjacency::Graph>&, std::vector<int>&, const Adjacency::ColoringDataHandler&, Index, float, float, bool);
template void Arch::assemble_vanka_host(const Intern::CSRTupleMatrixWrapper<float, std::uint32_t, 2>&, Intern::CSRTupleMatrixWrapper<float, std::uint32_t, 2>&, const std::vector<Adjacency::Graph>&,
                                        const std::vector<Adjacency::Graph>&, std::vector<int>&, const Adjacency::ColoringDataHandler&, Index, float, float, bool);
