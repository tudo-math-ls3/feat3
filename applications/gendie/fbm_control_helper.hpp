#pragma once

#include <kernel/runtime.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/string.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/assembly/slip_filter_assembler.hpp>
#include <kernel/assembly/error_computer.hpp>
#include <kernel/assembly/stokes_fbm_assembler.hpp>

#include <control/domain/parti_domain_control.hpp>
#include <control/stokes_blocked.hpp>

namespace FEAT::Control
{

  template<typename SpaceVelo_, typename SpacePres_>
  class StokesFBMDomainLevelBase :
    public Domain::StokesDomainLevel<typename SpaceVelo_::MeshType, typename SpaceVelo_::TrafoType, SpaceVelo_, SpacePres_>
  {
  public:
    typedef typename SpaceVelo_::MeshType MeshType;
    typedef typename SpaceVelo_::TrafoType TrafoType;
    typedef SpaceVelo_ SpaceVeloType;
    typedef SpacePres_ SpacePresType;
    typedef Control::Domain::StokesDomainLevel<MeshType, TrafoType, SpaceVeloType, SpacePresType> BaseClass;

    static constexpr bool fbm_support = true;

    /// the FBM assembler for this level
    std::unique_ptr<Assembly::StokesFBMAssembler<MeshType>> fbm_assembler;

    // inherit constructor
    using BaseClass::BaseClass;

    void create_fbm_assembler(const Dist::Comm& comm, const String& fbm_meshpart_names)
    {
      // get out mesh node
      const auto& mesh_node = *this->get_mesh_node();

      // allocate assembler
      fbm_assembler.reset(new Assembly::StokesFBMAssembler<MeshType>(*mesh_node.get_mesh()));

      // loop over all FBM meshparts
      std::deque<String> fbm_name_deque = fbm_meshpart_names.split_by_whitespaces();
      for(const String& name : fbm_name_deque)
      {
        // find meshpart node
        const auto* part_node = mesh_node.find_mesh_part_node(name);
        XASSERTM(part_node != nullptr, "FBM meshpart node '" + name + "'not found");

        // get the meshpart if it exists
        const auto* fbm_meshpart = part_node->get_mesh();
        if(fbm_meshpart)
          fbm_assembler->add_fbm_meshpart(*fbm_meshpart);
      }
      // synchronize over comm
      fbm_assembler->sync(mesh_node, comm);

      // compile the assembler
      fbm_assembler->compile();
    }
  }; // class DomainLevel

  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * \brief Navier-Stokes System Level base class
   */
  template<int dim_, typename LocalMatrixBlockA_, typename LocalMatrixBlockB_, typename LocalMatrixBlockD_, typename LocalScalarMatrix_>
  class StokesSystemLevelBase :
    public Control::StokesBlockedCombinedSystemLevel<dim_, typename LocalMatrixBlockA_::DataType, typename LocalMatrixBlockA_::IndexType, LocalMatrixBlockA_, LocalMatrixBlockB_, LocalMatrixBlockD_, LocalScalarMatrix_>
  {
  public:
    /// out base class
    typedef Control::StokesBlockedCombinedSystemLevel<dim_, typename LocalMatrixBlockA_::DataType, typename LocalMatrixBlockA_::IndexType, LocalMatrixBlockA_, LocalMatrixBlockB_, LocalMatrixBlockD_, LocalScalarMatrix_> BaseClass;

    /// the filtered local system matrix for Vanka
    typename BaseClass::LocalSystemMatrix local_matrix_sys;

    /// the velocity mass matrix
    typename BaseClass::GlobalMatrixBlockA velo_mass_matrix;

    typedef typename BaseClass::DataType DataType;
    typedef typename BaseClass::IndexType IndexType;
    typedef typename BaseClass::GlobalSystemVector GlobalSystemVector;
    typedef typename BaseClass::LocalVeloVector LocalVeloVector;
    typedef typename BaseClass::LocalPresVector LocalPresVector;
    static constexpr int dim = dim_;

    static constexpr bool fbm_support = false;

    template<typename Trafo_, typename SpaceVelo_>
    void assemble_velocity_laplace_matrix(Assembly::DomainAssembler<Trafo_>& dom_asm,
      const SpaceVelo_& space_velo, const DataType nu, bool defo, String cubature = "")
    {
      auto& loc_a = this->matrix_a.local();
      loc_a.format();
      Assembly::Common::LaplaceOperatorBlocked<dim> lapl_op;
      Assembly::Common::DuDvOperatorBlocked<dim> dudv_op;
      if(cubature.empty())
        cubature = "auto-degree:" + stringify(2*SpaceVelo_::local_degree+2);
      if(defo)
        Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_a, dudv_op, space_velo, cubature, nu);
      else
        Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_a, lapl_op, space_velo, cubature, nu);
    }

    void clear_velocity_laplace_matrix()
    {
      this->matrix_a.local() = BaseClass::LocalMatrixBlockA();
    }

    template<typename Trafo_, typename SpaceVelo_>
    void assemble_velocity_mass_matrix(Assembly::DomainAssembler<Trafo_>& dom_asm, const SpaceVelo_& space_velo, String cubature = "")
    {
      this->velo_mass_matrix = this->matrix_a.clone(LAFEM::CloneMode::Weak);
      auto& loc_m = this->velo_mass_matrix.local();
      loc_m.format();

      if(cubature.empty())
        cubature = "auto-degree:" + stringify(2*SpaceVelo_::local_degree+2);
      Assembly::Common::IdentityOperatorBlocked<dim> id_op;
      Assembly::assemble_bilinear_operator_matrix_1(dom_asm, loc_m, id_op, space_velo, cubature);
    }

    void compile_local_matrix()
    {
      // do we have to allocate the structures?
      if(this->local_matrix_sys.block_a().empty())
      {
        this->local_matrix_sys.block_a() = this->matrix_a.local().clone(LAFEM::CloneMode::Layout);
        this->local_matrix_sys.block_b() = this->matrix_b.local().clone(LAFEM::CloneMode::Layout);
        this->local_matrix_sys.block_d() = this->matrix_d.local().clone(LAFEM::CloneMode::Layout);
      }
      // copy local matrices (and sync 1 if necessary)
      this->matrix_a.convert_to_1(this->local_matrix_sys.block_a(), false);
      this->local_matrix_sys.block_b().copy(this->matrix_b.local());
      this->local_matrix_sys.block_d().copy(this->matrix_d.local());

      // apply velocity unit filters to A and B
      for(const auto& filter : this->get_local_velo_unit_filter_seq())
      {
        filter.second.filter_mat(this->local_matrix_sys.block_a());
        filter.second.filter_offdiag_row_mat(this->local_matrix_sys.block_b());
      }
    }
  }; // class SystemLevelBase

  /**
   * \brief Navier-Stokes System Level base class
   */
  template<int dim_, typename LocalMatrixBlockA_, typename LocalMatrixBlockB_, typename LocalMatrixBlockD_, typename LocalScalarMatrix_>
  class StokesFBMSystemLevelBase :
    public StokesSystemLevelBase<dim_, LocalMatrixBlockA_, LocalMatrixBlockB_, LocalMatrixBlockD_, LocalScalarMatrix_>
  {
  public:
    /// our base class
    typedef StokesSystemLevelBase<dim_, LocalMatrixBlockA_, LocalMatrixBlockB_, LocalMatrixBlockD_, LocalScalarMatrix_> BaseClass;

    /// the local interface filter
    typename BaseClass::LocalVeloUnitFilter filter_interface_fbm;

    typedef typename BaseClass::DataType DataType;
    typedef typename BaseClass::IndexType IndexType;
    typedef typename BaseClass::GlobalSystemVector GlobalSystemVector;
    typedef typename BaseClass::LocalVeloVector LocalVeloVector;
    typedef typename BaseClass::LocalPresVector LocalPresVector;
    static constexpr int dim = dim_;

    static constexpr bool fbm_support = true;

    /// FBM mask vectors of velocity and pressure
    std::vector<int> fbm_mask_velo, fbm_mask_pres;

    template<typename Mesh_, typename SpaceVelo_, typename SpacePres_>
    void assemble_fbm_filters(Assembly::StokesFBMAssembler<Mesh_>& fbm_asm, const SpaceVelo_& space_velo, const SpacePres_& space_pres, bool asm_mask, bool inner_fbm, bool no_scale)
    {
      // assemble velocity and pressure unit filters
      fbm_asm.assemble_inside_filter(this->get_local_velo_unit_filter_seq().find_or_add("fbm"), space_velo);
      fbm_asm.assemble_inside_filter(this->get_local_pres_unit_filter(), space_pres);

      // assemble interface filter for velocity
      if(inner_fbm)
        fbm_asm.assemble_interface_filter(filter_interface_fbm, space_velo, this->matrix_a, this->velo_mass_matrix, no_scale);
      else
        filter_interface_fbm = typename BaseClass::LocalVeloUnitFilter(space_velo.get_num_dofs());

      // assemble mask vectors on finest level
      if(asm_mask)
      {
        fbm_mask_velo.reserve(space_velo.get_num_dofs());
        for(int d(0); d <= dim; ++d)
        {
          for(auto k : fbm_asm.get_fbm_mask_vector(d))
            fbm_mask_velo.push_back(k);
        }
        fbm_mask_pres = fbm_asm.get_fbm_mask_vector(dim);
      }
    }

    template<typename SpacePres_>
    void assemble_pressure_mean_filter(const SpacePres_& space_pres, bool enable_fbm)
    {
      typename BaseClass::GlobalPresVector vec_glob_v(&this->gate_pres), vec_glob_w(&this->gate_pres);

      // fetch the local vectors
      LocalPresVector& vec_loc_v = vec_glob_v.local();
      LocalPresVector& vec_loc_w = vec_glob_w.local();

      // fetch the frequency vector of the pressure gate
      LocalPresVector& vec_loc_f = this->gate_pres._freqs;

      // assemble the mean filter
      Assembly::MeanFilterAssembler::assemble(vec_loc_v, vec_loc_w, space_pres, "gauss-legendre:2");

      // synchronize the vectors
      vec_glob_v.sync_1();
      vec_glob_w.sync_0();

      // apply pressure unit filter if FBM is enabled
      if(enable_fbm)
      {
        const auto& fil_p = this->get_local_pres_unit_filter();
        fil_p.filter_cor(vec_loc_v);
        fil_p.filter_def(vec_loc_w);
      }

      // create mean filter
      this->get_local_pres_mean_filter() = LocalPresMeanFilter(vec_loc_v.clone(), vec_loc_w.clone(), vec_loc_f.clone(), this->gate_pres.get_comm());
    }

    void apply_fbm_filter_to_rhs(GlobalSystemVector& vec_rhs) const
    {
      this->apply_fbm_filter_to_rhs(vec_rhs.local().first());
    }

    void apply_fbm_filter_to_rhs(LocalVeloVector& vec_rhs_v) const
    {
      this->filter_interface_fbm.filter_def(vec_rhs_v);
    }

    void apply_fbm_filter_to_def(GlobalSystemVector& vec_def_v, const GlobalSystemVector& vec_sol_v, const DataType factor) const
    {
      this->apply_fbm_filter_to_def(vec_def_v.local().first(), vec_sol_v.local().first(), factor);
    }

    void apply_fbm_filter_to_def(LocalVeloVector& vec_def_v, const LocalVeloVector& vec_sol_v, const DataType factor) const
    {
      if(this->filter_interface_fbm.used_elements() == Index(0))
        return;

      auto* vdef = vec_def_v.elements();
      const auto* vsol = vec_sol_v.elements();
      const IndexType* fidx = this->filter_interface_fbm.get_indices();
      const auto* fval = this->filter_interface_fbm.get_values();
      const IndexType* row_ptr = this->velo_mass_matrix.local().row_ptr();
      const IndexType* col_idx = this->velo_mass_matrix.local().col_ind();
      const auto* mval = this->velo_mass_matrix.local().val();

      Index n = this->filter_interface_fbm.used_elements();
      for(Index i(0); i < n; ++i)
      {
        IndexType row = fidx[i];
        vdef[row] = DataType(0);
        for(IndexType j(row_ptr[row]); j < row_ptr[row+1]; ++j)
        {
          vdef[row] += mval[j] * vsol[col_idx[j]];
        }
        vdef[row] *= fval[i][0] * factor;
      }
    }

  }; // class SystemLevelBase
}
