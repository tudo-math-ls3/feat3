// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once
#include <kernel/util/property_map.hpp>
#include <applications/gendie/materials.hpp>
#include <kernel/assembly/burgers_velo_material_assembly_job.hpp>
#include <kernel/assembly/domain_assembler_basic_jobs.hpp>
#include <kernel/assembly/common_operators.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/util/cuda_util.hpp>
#include <kernel/voxel_assembly/burgers_velo_material_assembler.hpp>
#include <kernel/cubature/dynamic_factory.hpp>

#include "format_helper.hpp"
#include "backend_plan.hpp"
#include "voxel_assembler_cache.hpp"
#include "template_helper.hpp"
#include "kernel/assembly/burgers_assembly_job.hpp"

#include <memory>
#include <type_traits>

namespace Gendie
{
  // CRTP interface for Flow System Assembler
  template<typename Derived_, typename Domain_>
  class SystemFlowAssemblerBaseCRTP
  {
  public:
    typedef Domain_ DomainType;
    typedef typename Domain_::WeightType SystemDataType;
    static constexpr int padlen = 30;
    static constexpr char pc = '.';
    String _cubature = "gauss-legendre:5";
    mutable FEAT::StopWatch watch_total, watch_restrict, watch_restrict_alloc;
    DomainType* domain;
    SystemDataType stiffness_factor = SystemDataType(1);
    SystemDataType theta = SystemDataType(0);
    FEAT::PreferredBackend asm_backend = FEAT::PreferredBackend::generic;
    /// per-level assembly backend selector (direct list and/or size heuristic)
    Intern::BackendSelector asm_backend_sel;
    /// use the voxel assembler for the velocity block (Carreau-family materials, Q2 hypercube spaces)
    bool voxel_asm = false;
    /// run the CUDA voxel assembly asynchronously (skip the end-of-kernel synchronization)
    bool async_asm = false;


  protected:
    /**
      * \brief Casts \c this to its true type
      *
      * \returns A Derived_ reference to \c this
      */
    Derived_& cast() {return static_cast<Derived_&>(*this);}

    /// \copydoc cast()
    const Derived_& cast() const {return static_cast<const Derived_&>(*this);}

    bool _parse([[maybe_unused]] const FEAT::PropertyMap* prop)
    {
      return true;
    }

    void _set_domain(DomainType& domain_)
    {
      domain = &domain_;
    }

    void _init()
    {

    }

    void _done()
    {

    }

    void _set_jacobian(double DOXY(weight))
    {
      XABORTM("Has to be implemented by derived class");
    }

    template<typename LevelType_, typename VecType_>
    void _assemble_matrices(std::deque<std::shared_ptr<LevelType_>>& DOXY(system_levels), const VecType_& DOXY(vec_sol)) const
    {
      XABORTM("Specialization has to be assembled by base class");
    }

    template<typename LevelType_>
    void _assemble_matrices_isotrop(std::deque<std::shared_ptr<LevelType_>>& DOXY(system_levels)) const
    {
      XABORTM("Specialization has to be assembled by base class");
    }

    void _set_theta(SystemDataType _theta)
    {
      this->theta = _theta;
    }

    template<typename LevelType_, typename VecType_>
    void _fill_convections(const std::deque<std::shared_ptr<LevelType_>>& system_levels, std::vector<typename LevelType_::GlobalVeloVector>& conv_vectors, const VecType_& vec_sol) const
    {
      this->watch_restrict_alloc.start();
      typedef typename LevelType_::GlobalVeloVector ConvVecType;
      constexpr bool same_vecs = std::is_same_v<VecType_, ConvVecType>;
      if constexpr (same_vecs)
      {
        conv_vectors.push_back(ConvVecType(system_levels.front()->gate_velo, vec_sol.local().template at<0>(), FEAT::LAFEM::CloneMode::Shallow));
      }
      else
      {
        conv_vectors.push_back(system_levels.front()->create_global_vector_velo());
        conv_vectors.back().local().convert(vec_sol.local().template at<0>());
      }
      for(std::size_t i = 1; i < system_levels.size(); ++i)
      {
        conv_vectors.push_back(system_levels.at(i)->create_global_vector_velo());
        conv_vectors.back().format();
      }
      this->watch_restrict_alloc.stop();
    }

    template<typename LevelType_>
    void _restrict_convections(const std::deque<std::shared_ptr<LevelType_>>& system_levels, std::vector<typename LevelType_::GlobalVeloVector>& conv_vectors) const
    {
      this->watch_restrict.start();
      for(std::size_t i = 0; i < system_levels.size(); ++i)
      {
        // do we need to participate in the communication?
        if((i+1) >= this->domain->size_virtual())
          break;
        // does this process have another system level?
        if((i+1) < system_levels.size())
        {
          const auto& vec_fine = conv_vectors.at(i);
          // create a coarse mesh velocity vector
          auto& vec_crs = conv_vectors.at(i+1);

          // truncate fine mesh velocity vector
          system_levels.at(i)->transfer_velo.trunc(vec_fine, vec_crs);
        }
        else
        {
          const auto& vec_fine = conv_vectors.at(i);
          // this process is a child, so send truncation to parent
          system_levels.at(i)->transfer_velo.trunc_send(vec_fine);
        }
      }
      this->watch_restrict.stop();
    }

    void _set_backend() const
    {
      FEAT::Backend::set_preferred_backend(this->asm_backend);
    }

    FEAT::String _format_timings() const
    {
      return FEAT::String{};
    }

    void _set_cubature(const String& cubature)
    {
      _cubature = cubature;
    }

    public:

    SystemFlowAssemblerBaseCRTP() = default;

    SystemFlowAssemblerBaseCRTP(DomainType& domain_, const FEAT::PropertyMap* prop)
     : domain(&domain_)
    {
      _parse(prop);
    }

    virtual ~SystemFlowAssemblerBaseCRTP() = default;

    SystemFlowAssemblerBaseCRTP(const SystemFlowAssemblerBaseCRTP&) = default;
    SystemFlowAssemblerBaseCRTP(SystemFlowAssemblerBaseCRTP&&) noexcept = default;
    SystemFlowAssemblerBaseCRTP& operator=(const SystemFlowAssemblerBaseCRTP&) = default;
    SystemFlowAssemblerBaseCRTP& operator=(SystemFlowAssemblerBaseCRTP&&) noexcept = default;

    void init()
    {
      watch_total.start();
      this->cast()._init();
      watch_total.stop();
    }

    void done()
    {
      watch_total.start();
      this->cast()._done();
      watch_total.stop();
    }

    /**
     * \brief Assembles the matrices defined by the system level
     *
     * \tparam VecType_ A global vectortype convertable to the inner system vectortype
     *
     * \param[in] vec_sol The local solution vector used as convection for assembly
     */
    template<typename LevelType_, typename VecType_>
    void assemble_matrices(std::deque<std::shared_ptr<LevelType_>>& system_levels, const VecType_& vec_sol) const
    {
      FEAT::Util::ScopedNvtxRange _nvtx("Asm SysMat");
      set_backend();
      watch_total.start();
      this->cast()._assemble_matrices(system_levels, vec_sol);
      watch_total.stop();
    }

    template<typename LevelType_>
    void assemble_matrices_isotrop(std::deque<std::shared_ptr<LevelType_>>& system_levels) const
    {
      FEAT::Util::ScopedNvtxRange _nvtx("Asm SysIso");
      set_backend();
      watch_total.start();
      this->cast()._assemble_matrices_isotrop(system_levels);
      watch_total.stop();
    }

    template<typename VecType_>
    void assemble_diagonal(VecType_& vec_diag, const VecType_& vec_sol, PreferredBackend backend = PreferredBackend::generic) const
    {
      FEAT::Util::ScopedNvtxRange _nvtx("Asm Diag");
      FEAT::Backend::set_preferred_backend(backend);
      watch_total.start();
      this->cast()._assemble_diagonal(vec_diag, vec_sol);
      watch_total.stop();
    }

    auto& get_domain()
    {
      return *domain;
    }

    /**
     * \brief Sets assembler such that the jacobian (if available) of the underlying system is assembled, weighted by a given parameter
     *
     * \param[in] weight Assembles the jacobian with a specific weight. How this parameter is interpreted exactly is specific to the actual system being assembled
     *                   but in general, weight == 0 means the standard system is assembled, weight > 0 some form ob (Pseudo)-Newton is assembled
     */
    void set_jacobian(double weight)
    {
      this->cast()._set_jacobian(weight);
    }

    SystemDataType calc_stiffness_factor() const
    {
      return this->cast()._calc_stiffness_factor();
    }

    void set_stiffness_factor(SystemDataType factor)
    {
      this->stiffness_factor = factor;
    }

    SystemDataType get_stiffness_factor() const
    {
      return this->stiffness_factor;
    }

    void set_theta(SystemDataType _theta)
    {
      return this->cast()._set_theta(_theta);
    }

    void set_backend() const
    {
      this->_set_backend();
    }

    /// sets the per-level assembly backend selector
    void set_asm_backend_selector(const Intern::BackendSelector& sel)
    {
      this->asm_backend_sel = sel;
    }

    /// enables/disables the voxel assembler for the velocity block
    void set_voxel_asm(bool use)
    {
      this->voxel_asm = use;
    }

    /// enables/disables asynchronous CUDA voxel assembly (no end-of-kernel synchronization)
    void set_async_asm(bool use)
    {
      this->async_asm = use;
    }

    FEAT::String format_timings() const
    {
      return this->cast()._format_timings();
    }

    double get_total_time_elapsed() const
    {
      return this->watch_total.elapsed();
    }

    void set_cubature(const String& cubature)
    {
      this->cast()._set_cubature(cubature);
    }

  }; //class SystemFlowAssemblerBaseCRTP

  /**
   * \brief Standard Burgers Carreau assembler
   *
   */
  template<typename Domain_>
  class CarreauSingleMatSystemFlowAssembler : public SystemFlowAssemblerBaseCRTP<CarreauSingleMatSystemFlowAssembler<Domain_>, Domain_>
  {
  public:
    typedef SystemFlowAssemblerBaseCRTP<CarreauSingleMatSystemFlowAssembler, Domain_> BaseClass;
    typedef typename BaseClass::DomainType DomainType;
    typedef typename BaseClass::SystemDataType SystemDataType;
    friend BaseClass;

    SystemDataType reg_eps = 1E-100;
    bool navier = false;

  protected:
    /// per-level cache of voxel assemblers, owned by the active solve (see _init/_done)
    mutable Intern::VoxelAssemblerCache _voxel_cache;

    bool _parse([[maybe_unused]] const FEAT::PropertyMap* prop)
    {

      return true;
    }

    void _init()
    {
      // resources are bound to the solve lifecycle; start from a clean slate
      _voxel_cache.clear();
    }

    void _done()
    {
      // free the voxel assembler resources when the solve tears down
      _voxel_cache.clear();
    }

    void _set_jacobian(double weight)
    {
      _weight_jacobian = weight;
    }

    SystemDataType _calc_stiffness_factor() const
    {
      const auto aT = this->_material->get_abiat_aT(this->_temperature);
      return this->_material->get_default_viscosity(_mesh_scale, aT);
    }

    template<typename Job_>
    void _set_job(Job_& job) const
    {
      typedef typename Job_::DataType AssemblyDataType;
      const auto aT = this->_material->get_abiat_aT(this->_temperature);
      job.deformation = true;
      job.nu = AssemblyDataType(this->_material->get_default_viscosity(_mesh_scale, aT));
      job.aT = AssemblyDataType(aT);
      job.sd_delta = AssemblyDataType(0);
      job.sd_nu = job.nu;
      job.beta = navier ? AssemblyDataType(_material->get_density_gram_per_unit(this->_mesh_scale)) : AssemblyDataType(0);
      job.reg_eps = AssemblyDataType(reg_eps);
      job.frechet_beta = AssemblyDataType(SystemDataType(job.beta) * this->_weight_jacobian);
      job.frechet_material = AssemblyDataType(SystemDataType(1) * this->_weight_jacobian);
      job.theta = AssemblyDataType(this->theta);
      job.alpha = AssemblyDataType(SystemDataType(1) / this->stiffness_factor);
    }

    /**
     * \brief Assembles the velocity A-block via the voxel assembler (CPU/CUDA by preferred backend).
     *
     * Only instantiated for Q2/Lagrange2 hypercube velocity spaces (see the \c if \c constexpr guard
     * in _assemble_matrices). The gendie Carreau-family material is mapped approximately onto the
     * voxel carreau functor mu = a_T * mu_0 * (1 + lambda*a_T*gamma_dot)^(-exp); this is exact only
     * for the plateau exponent d=1, unclamped shear rate. Parameters mirror \c _set_job.
     */
    template<typename MatrixA_, typename ConvVec_, typename Space_, typename Coloring_>
    void _assemble_matrix_a_voxel(std::size_t level, MatrixA_& matrix_a, const ConvVec_& convect, const Space_& space, const Coloring_& coloring) const
    {
      typedef typename MatrixA_::DataType DT;
      typedef typename MatrixA_::IndexType IT;
      // reuse the per-level voxel assembler (expensive to construct); only the scalar material
      // parameters below change between calls
      auto& voxel = _voxel_cache.template get<Space_, DT, IT>(level, space, coloring);
      // run the assembly asynchronously if requested (skips the device-wide sync/error-check at the end)
      voxel.runtime_params.synchronize = !this->async_asm;

      const auto aT = this->_material->get_abiat_aT(this->_temperature);
      const SystemDataType visco_scale = SystemDataType(1000) / this->_mesh_scale;
      const auto data = this->_material->get_viscosity_model()->get_data();

      voxel.deformation = true;
      voxel.nu = DT(this->_material->get_default_viscosity(this->_mesh_scale, aT));
      voxel.sd_nu = voxel.nu;
      voxel.sd_delta = DT(0);
      voxel.beta = navier ? DT(this->_material->get_density_gram_per_unit(this->_mesh_scale)) : DT(0);
      voxel.frechet_beta = DT(SystemDataType(voxel.beta) * this->_weight_jacobian);
      voxel.frechet_material = DT(SystemDataType(1) * this->_weight_jacobian);
      voxel.theta = DT(this->theta);
      voxel.reg_eps = DT(reg_eps);
      // approximate mapping of the Carreau-family model onto the voxel carreau functor
      voxel.material_type = FEAT::VoxelAssembly::MaterialType::carreau;
      voxel.mu_0 = DT(SystemDataType(data[0]) * visco_scale);
      voxel.lambda = DT(data[1]);
      voxel.exp = DT(SystemDataType(1) - SystemDataType(data[2]));
      voxel.a_T = DT(aT);

      const DT alpha = DT(SystemDataType(1) / this->stiffness_factor);
      voxel.assemble_matrix1(matrix_a, convect, space, FEAT::Cubature::DynamicFactory(this->_cubature), alpha);
    }

    /// always assume that matrix structure is preassambled
    template<typename LevelType_, typename VecType_>
    void _assemble_matrices(std::deque<std::shared_ptr<LevelType_>>& system_levels, const VecType_& vec_sol) const
    {
      typedef typename LevelType_::GlobalVeloVector ConvVecType;
      typedef typename LevelType_::DataType AssemblyDataType;

      // constexpr bool same_vecs = std::is_same_v<VecType_, ConvVecType>;
      // assemble convections
      std::vector<ConvVecType> conv_vecs;
      conv_vecs.reserve(system_levels.size());
      this->_fill_convections(system_levels, conv_vecs, vec_sol);
      this->_restrict_convections(system_levels, conv_vecs);

      this->watch_assembly.start();
      for(std::size_t i = 0; i < system_levels.size(); ++i)
      {
        if(i == 0) this->watch_assembly_fine.start();
        // setup
        auto& cur_sys = *system_levels[i];
        auto& cur_dom = *this->domain->at(i);
        ASSERTM(!cur_sys.matrix_a.local().empty(), "A Matrix not allocated");
        ASSERTM(!cur_sys.matrix_b.local().empty(), "B Matrix not allocated");
        ASSERTM(!cur_sys.matrix_d.local().empty(), "D Matrix not allocated");

        // resolve the per-level assembly backend (direct list or A-block size heuristic)
        FEAT::PreferredBackend lvl_backend = this->asm_backend_sel.backend_for(i, cur_sys.matrix_a.local().num_rows());

        bool did_voxel = false;
        // voxel assembler path: only for Q2/Lagrange2 hypercube velocity spaces on domains providing
        // an element coloring, only for Carreau-family materials, and only when this level actually
        // wants the CUDA backend (the voxel assembler exists for the CUDA path; CPU levels use the
        // regular domain-assembler job below)
        if constexpr(std::is_same_v<std::decay_t<decltype(cur_dom.space_velo)>, FEAT::VoxelAssembly::Q2StandardHyperCube<LevelType_::dim>> && Intern::supports_element_coloring<Domain_>)
        {
#ifdef FEAT_HAVE_CUDA
          if(this->voxel_asm && lvl_backend == FEAT::PreferredBackend::cuda && is_carreau_family(this->_material->get_current_visc_model_type()))
          {
            FEAT::Backend::set_preferred_backend(lvl_backend);
            cur_sys.matrix_a.local().format();
            this->_assemble_matrix_a_voxel(i, cur_sys.matrix_a.local(), conv_vecs[i].local(), cur_dom.space_velo, this->domain->at(i)->element_coloring);
            did_voxel = true;
          }
#endif // FEAT_HAVE_CUDA
        }

        if(!did_voxel)
        {
          if(lvl_backend == FEAT::PreferredBackend::cuda)
          {
            // the CPU domain-assembler Burgers job has no CUDA implementation; fall back to generic
            this->domain->comm().print("WARNING: assembly backend 'cuda' on level " + FEAT::stringify(i) + " needs voxel-assembly with a Carreau-family material and a Q2 space; falling back to 'generic'");
            lvl_backend = FEAT::PreferredBackend::generic;
          }
          FEAT::Backend::set_preferred_backend(lvl_backend);
          FEAT::Assembly::BurgersVeloMaterialBlockedMatrixAssemblyJob burgers_mat_job(
            cur_sys.matrix_a.local(), conv_vecs[i].local(), cur_dom.space_velo,
            this->_cubature, this->_material->get_visc_func(AssemblyDataType(_mesh_scale)), this->_material->get_visc_der_func(AssemblyDataType(_mesh_scale))
          );
          _set_job(burgers_mat_job);

          cur_sys.matrix_a.local().format();

          // actual assembly
          cur_dom.domain_asm.assemble(burgers_mat_job);
        }
        if(i == 0) this->watch_assembly_fine.stop();
      }
      this->watch_assembly.stop();
    }

    template<typename LevelType_>
    void _assemble_matrices_isotrop(std::deque<std::shared_ptr<LevelType_>>& system_levels) const
    {
      typedef typename LevelType_::DataType AssemblyDataType;
      this->watch_assembly.start();

      for(std::size_t i = 0; i < system_levels.size(); ++i)
      {
        if(i == 0) this->watch_assembly_fine.start();
        // setup
        auto& cur_sys = *system_levels[i];
        auto& cur_dom = *this->domain->at(i);
        ASSERTM(!cur_sys.matrix_a.local().empty(), "A Matrix not allocated");
        ASSERTM(!cur_sys.matrix_b.local().empty(), "B Matrix not allocated");
        ASSERTM(!cur_sys.matrix_d.local().empty(), "D Matrix not allocated");
        Assembly::Common::LaplaceOperatorBlocked<LevelType_::dim> laplace_op;

        cur_sys.matrix_a.local().format();

        // actual assembly
        Assembly::assemble_bilinear_operator_matrix_1(cur_dom.domain_asm, cur_sys.matrix_a.local(), laplace_op, cur_dom.space_velo, this->_cubature,
          AssemblyDataType(1));
        if(i == 0) this->watch_assembly_fine.stop();
      }
      this->watch_assembly.stop();
    }

    template<typename VecType_>
    void _assemble_diagonal(VecType_& vec_diag, const VecType_& vec_sol) const
    {
      typedef typename VecType_::DataType AssemblyDataType;
      this->watch_assembly.start();
      this->watch_assembly_fine.start();
      auto& cur_dom = *this->domain->at(0);
      FEAT::Assembly::BurgersVeloMaterialBlockedDiagonalAssemblyJob burgers_mat_job(
        vec_diag, vec_sol, cur_dom.space_velo,
        this->_cubature, this->_material->get_visc_func(AssemblyDataType(_mesh_scale)), this->_material->get_visc_der_func(AssemblyDataType(_mesh_scale))
      );
      _set_job(burgers_mat_job);

      // actual assembly
      cur_dom.domain_asm.assemble(burgers_mat_job);
      this->watch_assembly_fine.stop();
      this->watch_assembly.stop();
    }

    FEAT::String _format_timings() const
    {
      enum class _TimeID : std::size_t
      {
        total = 0,
        assembly,
        assembly_fine,
        restrict,
        restrict_alloc,
        num_entries
      };
      FEAT::String s;
      double timings_max[to_underlying(_TimeID::num_entries)], timings_min[to_underlying(_TimeID::num_entries)], timings[to_underlying(_TimeID::num_entries)];
      timings_min[to_underlying(_TimeID::total)] = timings_max[to_underlying(_TimeID::total)] = timings[to_underlying(_TimeID::total)] = this->watch_total.elapsed();
      timings_min[to_underlying(_TimeID::assembly)] = timings_max[to_underlying(_TimeID::assembly)] = timings[to_underlying(_TimeID::assembly)] = this->watch_assembly.elapsed();
      timings_min[to_underlying(_TimeID::assembly_fine)] = timings_max[to_underlying(_TimeID::assembly_fine)] = timings[to_underlying(_TimeID::assembly_fine)] = this->watch_assembly_fine.elapsed();
      timings_min[to_underlying(_TimeID::restrict)] = timings_max[to_underlying(_TimeID::restrict)] = timings[to_underlying(_TimeID::restrict)] = this->watch_restrict.elapsed();
      timings_min[to_underlying(_TimeID::restrict_alloc)] = timings_max[to_underlying(_TimeID::restrict_alloc)] = timings[to_underlying(_TimeID::restrict_alloc)] = this->watch_restrict_alloc.elapsed();
      // sync our timings
      this->domain->comm().allreduce(timings_max, timings_max, to_underlying(_TimeID::num_entries), FEAT::Dist::op_max);
      this->domain->comm().allreduce(timings_min, timings_min, to_underlying(_TimeID::num_entries), FEAT::Dist::op_min);

      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += FEAT::String("\n------------------------------------Timings Single Material System Assembler----------------------\n");
      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += format_subtime_mm("Allocate Restrictions", timings[to_underlying(_TimeID::restrict_alloc)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::restrict_alloc)], timings_max[to_underlying(_TimeID::restrict_alloc)], this->padlen);
      s += format_subtime_mm("Restrict Convection", timings[to_underlying(_TimeID::restrict)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::restrict)], timings_max[to_underlying(_TimeID::restrict)], this->padlen);
      s += format_subtime_mm("Assemble Systems", timings[to_underlying(_TimeID::assembly)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::assembly)], timings_max[to_underlying(_TimeID::assembly)], this->padlen);
      s += format_subtime_mm("Assemble Fine System", timings[to_underlying(_TimeID::assembly_fine)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::assembly_fine)], timings_max[to_underlying(_TimeID::assembly_fine)], this->padlen);
      s += format_subtime_mm("System Assembly Total Time", timings[to_underlying(_TimeID::total)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::total)], timings_max[to_underlying(_TimeID::total)], this->padlen);

      return s;
    }

    const Material<SystemDataType>* _material;
    mutable FEAT::StopWatch watch_assembly, watch_assembly_fine;
    SystemDataType _temperature;
    SystemDataType _mesh_scale;
    SystemDataType _weight_jacobian;


  public:
    CarreauSingleMatSystemFlowAssembler() = default;

    CarreauSingleMatSystemFlowAssembler(DomainType& domain_, const FEAT::PropertyMap* prop, const Material<SystemDataType>& material, SystemDataType temperature, SystemDataType mesh_scale)
      : BaseClass(domain_, prop),
      _material(&material),
      _temperature(temperature),
      _mesh_scale(mesh_scale),
      _weight_jacobian(SystemDataType(0))
    {
      _parse(prop);
    }

    virtual ~CarreauSingleMatSystemFlowAssembler() = default;

    CarreauSingleMatSystemFlowAssembler(const CarreauSingleMatSystemFlowAssembler&) = default;
    CarreauSingleMatSystemFlowAssembler(CarreauSingleMatSystemFlowAssembler&&) noexcept = default;
    CarreauSingleMatSystemFlowAssembler& operator=(const CarreauSingleMatSystemFlowAssembler&) = default;
    CarreauSingleMatSystemFlowAssembler& operator=(CarreauSingleMatSystemFlowAssembler&&) noexcept = default;

    void set_material(const Material<SystemDataType>& material, SystemDataType mesh_scale)
    {
      _material = &material;
      _mesh_scale = mesh_scale;
    }

  }; //class CarrauSystemFlowAssembler<...>

  /**
   * \brief Standard Burgers Carreau assembler
   *
   */
  template<typename Domain_>
  class BurgersSystemFlowAssembler : public SystemFlowAssemblerBaseCRTP<BurgersSystemFlowAssembler<Domain_>, Domain_>
  {
  public:
    typedef SystemFlowAssemblerBaseCRTP<BurgersSystemFlowAssembler, Domain_> BaseClass;
    typedef typename BaseClass::DomainType DomainType;
    typedef typename BaseClass::SystemDataType SystemDataType;
    friend BaseClass;

  protected:

    bool _parse([[maybe_unused]] const FEAT::PropertyMap* prop)
    {

      return true;
    }

    void _init()
    {
    }

    void _done()
    {
    }

    void _set_jacobian(double weight)
    {
      _weight_jacobian = weight;
    }

    SystemDataType _calc_stiffness_factor() const
    {
      return nu;
    }

    template<typename Job_>
    void _set_job(Job_& job) const
    {
      typedef typename Job_::DataType AssemblyDataType;
      job.deformation = deformation;
      job.nu = AssemblyDataType(nu);
      job.sd_delta = AssemblyDataType(0);
      job.sd_nu = AssemblyDataType(nu);
      job.beta = AssemblyDataType(beta);
      job.frechet_beta = AssemblyDataType(SystemDataType(job.beta) * this->_weight_jacobian);
      job.theta = AssemblyDataType(theta);
      job.alpha = AssemblyDataType(SystemDataType(1) / this->stiffness_factor);
    }

    /// always assume that matrix structure is preassambled
    template<typename LevelType_, typename VecType_>
    void _assemble_matrices(std::deque<std::shared_ptr<LevelType_>>& system_levels, const VecType_& vec_sol) const
    {
      typedef typename LevelType_::GlobalVeloVector ConvVecType;

      // constexpr bool same_vecs = std::is_same_v<VecType_, ConvVecType>;
      // assemble convections
      std::vector<ConvVecType> conv_vecs;
      conv_vecs.reserve(system_levels.size());
      this->_fill_convections(system_levels, conv_vecs, vec_sol);
      this->_restrict_convections(system_levels, conv_vecs);

      this->watch_assembly.start();
      for(std::size_t i = 0; i < system_levels.size(); ++i)
      {
        if(i == 0) this->watch_assembly_fine.start();
        // setup
        auto& cur_sys = *system_levels[i];
        auto& cur_dom = *this->domain->at(i);
        ASSERTM(!cur_sys.matrix_a.local().empty(), "A Matrix not allocated");
        ASSERTM(!cur_sys.matrix_b.local().empty(), "B Matrix not allocated");
        ASSERTM(!cur_sys.matrix_d.local().empty(), "D Matrix not allocated");
        FEAT::Assembly::BurgersBlockedMatrixAssemblyJob burgers_mat_job(
          cur_sys.matrix_a.local(), conv_vecs[i].local(), cur_dom.space_velo,
          this->_cubature);

        cur_sys.matrix_a.local().format();

        // actual assembly
        cur_dom.domain_asm.assemble(burgers_mat_job);
        if(i == 0) this->watch_assembly_fine.stop();
      }
      this->watch_assembly.stop();
    }

    template<typename LevelType_>
    void _assemble_matrices_isotrop(std::deque<std::shared_ptr<LevelType_>>& system_levels) const
    {
      typedef typename LevelType_::DataType AssemblyDataType;
      this->watch_assembly.start();

      for(std::size_t i = 0; i < system_levels.size(); ++i)
      {
        if(i == 0) this->watch_assembly_fine.start();
        // setup
        auto& cur_sys = *system_levels[i];
        auto& cur_dom = *this->domain->at(i);
        ASSERTM(!cur_sys.matrix_a.local().empty(), "A Matrix not allocated");
        ASSERTM(!cur_sys.matrix_b.local().empty(), "B Matrix not allocated");
        ASSERTM(!cur_sys.matrix_d.local().empty(), "D Matrix not allocated");
        Assembly::Common::LaplaceOperatorBlocked<LevelType_::dim> laplace_op;

        cur_sys.matrix_a.local().format();

        // actual assembly
        Assembly::assemble_bilinear_operator_matrix_1(cur_dom.domain_asm, cur_sys.matrix_a.local(), laplace_op, cur_dom.space_velo, this->_cubature,
          nu * AssemblyDataType(SystemDataType(1) / this->stiffness_factor));
        if(i == 0) this->watch_assembly_fine.stop();
      }
      this->watch_assembly.stop();
    }

    template<typename VecType_>
    void _assemble_diagonal(VecType_& vec_diag, const VecType_& vec_sol) const
    {
      this->watch_assembly.start();
      this->watch_assembly_fine.start();
      auto& cur_dom = *this->domain->at(0);
      FEAT::Assembly::BurgersBlockedDiagonalAssemblyJob burgers_mat_job(
        vec_diag, vec_sol, cur_dom.space_velo, this->_cubature);
      _set_job(burgers_mat_job);

      // actual assembly
      cur_dom.domain_asm.assemble(burgers_mat_job);
      this->watch_assembly_fine.stop();
      this->watch_assembly.stop();
    }

    FEAT::String _format_timings() const
    {
      enum class _TimeID : std::size_t
      {
        total = 0,
        assembly = 1,
        assembly_fine = 2,
        restrict = 3,
        restrict_alloc = 4,
        num_entries = 5
      };
      FEAT::String s;
      double timings_max[to_underlying(_TimeID::num_entries)], timings_min[to_underlying(_TimeID::num_entries)], timings[to_underlying(_TimeID::num_entries)];
      timings_min[to_underlying(_TimeID::total)] = timings_max[to_underlying(_TimeID::total)] = timings[to_underlying(_TimeID::total)] = this->watch_total.elapsed();
      timings_min[to_underlying(_TimeID::assembly)] = timings_max[to_underlying(_TimeID::assembly)] = timings[to_underlying(_TimeID::assembly)] = this->watch_assembly.elapsed();
      timings_min[to_underlying(_TimeID::assembly_fine)] = timings_max[to_underlying(_TimeID::assembly_fine)] = timings[to_underlying(_TimeID::assembly_fine)] = this->watch_assembly_fine.elapsed();
      timings_min[to_underlying(_TimeID::restrict)] = timings_max[to_underlying(_TimeID::restrict)] = timings[to_underlying(_TimeID::restrict)] = this->watch_restrict.elapsed();
      timings_min[to_underlying(_TimeID::restrict_alloc)] = timings_max[to_underlying(_TimeID::restrict_alloc)] = timings[to_underlying(_TimeID::restrict_alloc)] = this->watch_restrict_alloc.elapsed();
      // sync our timings
      this->domain->comm().allreduce(timings_max, timings_max, to_underlying(_TimeID::num_entries), FEAT::Dist::op_max);
      this->domain->comm().allreduce(timings_min, timings_min, to_underlying(_TimeID::num_entries), FEAT::Dist::op_min);

      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += FEAT::String("\n------------------------------------Timings Single Material System Assembler----------------------\n");
      s += FEAT::String("\n--------------------------------------------------------------------------------------------------\n");
      s += format_subtime_mm("Allocate Restrictions", timings[to_underlying(_TimeID::restrict_alloc)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::restrict_alloc)], timings_max[to_underlying(_TimeID::restrict_alloc)], this->padlen);
      s += format_subtime_mm("Restrict Convection", timings[to_underlying(_TimeID::restrict)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::restrict)], timings_max[to_underlying(_TimeID::restrict)], this->padlen);
      s += format_subtime_mm("Assemble Systems", timings[to_underlying(_TimeID::assembly)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::assembly)], timings_max[to_underlying(_TimeID::assembly)], this->padlen);
      s += format_subtime_mm("Assemble Fine System", timings[to_underlying(_TimeID::assembly_fine)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::assembly_fine)], timings_max[to_underlying(_TimeID::assembly_fine)], this->padlen);
      s += format_subtime_mm("System Assembly Total Time", timings[to_underlying(_TimeID::total)], timings[to_underlying(_TimeID::total)], timings_min[to_underlying(_TimeID::total)], timings_max[to_underlying(_TimeID::total)], this->padlen);

      return s;
    }

  public:
    mutable FEAT::StopWatch watch_assembly, watch_assembly_fine;
    SystemDataType nu, beta, theta, _weight_jacobian;
    bool deformation;


  public:
    BurgersSystemFlowAssembler() = default;

    BurgersSystemFlowAssembler(DomainType& domain_, const FEAT::PropertyMap* prop)
      : BaseClass(domain_, prop),
      nu(SystemDataType(1)),
      beta(SystemDataType(0)),
      theta(SystemDataType(0)),
      _weight_jacobian(SystemDataType(0)),
      deformation(false)
    {
      _parse(prop);
    }

    virtual ~BurgersSystemFlowAssembler() = default;

    BurgersSystemFlowAssembler(const BurgersSystemFlowAssembler&) = default;
    BurgersSystemFlowAssembler(BurgersSystemFlowAssembler&&) noexcept = default;
    BurgersSystemFlowAssembler& operator=(const BurgersSystemFlowAssembler&) = default;
    BurgersSystemFlowAssembler& operator=(BurgersSystemFlowAssembler&&) noexcept = default;

  }; //class BurgersSystemFlowAssembler<...>
}
