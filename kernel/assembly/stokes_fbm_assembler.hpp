// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

// includes, FEAT
#include <kernel/base_header.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/dist.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_bcsr.hpp>
#include <kernel/lafem/unit_filter.hpp>
#include <kernel/lafem/unit_filter_blocked.hpp>
#include <kernel/global/matrix.hpp>

// includes, system
#include <array>
#include <vector>
#include <memory>

namespace FEAT
{
  namespace Assembly
  {
    /**
     * \brief Stokes Fictitious Boundary Method assembler class
     *
     * This class can be used to assemble a triplet of unit filters for the velocity and pressure spaces of a Stokes system.
     * The three unit filters that result from this assembly are:
     * - An UnitFilterBlocked thats filter out all velocity DOFs which are wholly inside the FBM region by setting them to 0
     * - An UnitFilter that filters out all pressure DOFs which are wholly inside the FBM region by setting them to 0
     * - An UnitFilterBlocked for the velocity that represents the FBM interface DOFs which are handled by the div-free
     *   L2-projection, which is incorporated into the system matrix by replacing the corresponding A matrix rows by
     *   (properly scaled) mass matrix rows.
     *
     * To use this assembler for the proper handling of FBM boundary conditions, proceed as follows:
     * -# Assemble the matrix block A of the Stokes system (typically a Laplace matrix) as well as the velocity mass matrix
     * -# Create a StokesFBMAssembler object
     * -# Add all meshparts that represent the FBM region(s) to the assembler by calling add_fbm_meshpart()
     * -# In a parallel simulation, call the sync() function to synchronize the FBM data across all processes
     * -# Call the compile() function to finish the setup of the assembler
     * -# Use the assemble_inside_filter() function to assemble the unit-filters for the velocity and pressure spaces
     * -# Use the assemble_interface_filter() function to assemble a special unit filter that will represent the interface
     *    DOFs in the velocity space. Then, every time you (re)assemble the matrix block A of the system matrix, call the
     *    LAFEM::UnitFilterBlocked::filter_weak_matrix_rows() to replace the matrix A interface rows by the corresponding
     *    rows of the velocity mass matrix.
     * -# If you require a pressure mean filter, assemble the primal and dual vectors of the filter by using the
     *    Assembly::MeanFilterAssembler and filter them by applying the pressure unit filter's filter_cor function
     *    before pushing them into the mean filter.
     * -# If you require the characteristic function of the FBM region interface, e.g. for the volumetric computation of
     *    the drag and lift forces, then you can use the assemble_characteristic_vector() function for that.
     *
     * \author Peter Zajac
     */
    template<typename MeshType_>
    class StokesFBMAssembler
    {
    public:
      /// our mesh type
      typedef MeshType_ MeshType;
      /// our shape dimension
      static constexpr int shape_dim = MeshType::shape_dim;

    protected:
      /// a reference to our mesh object
      const MeshType& _mesh;

      /**
       * \brief Array of FBM mask vectors
       *
       * For each dimension, the mask vector contains an int value for each entity of that dimension, which can be equal
       * to one of the following four values:
       *
       * - 0: Indicates that the entity itself as well as all of its sub-dimensional entities are not contained in the
       *      FBM region, thus these are the entities that are <em>wholly outside</em> the FBM region.
       * - 1: Indicates that the entity itself is inside the FBM region, but at least one of its sub-dimensional faces
       *      is outside the FBM region. These are the entities that are <em>partially inside</em> the FBM region.
       * - 2: Indicates that the entity itself is outside the FBM region, but all of its sub-dimensional faces are
       *      inside in the FBM region. This is a borderline case which only happens if the mesh is too coarse to properly
       *      resolve a non-convex FBM region. These DOFs are treated as if they had the value 0.
       * - 3: Indicates that the entity itself as well as all of its sub-dimensional entities are inside the FBM region,
       *      thus these are the entities that are <em>wholly inside</em> the FBM region.
       */
      std::array<std::vector<int>, shape_dim+1> _fbm_masks;
      /// the meshpart representing the FBM inside region
      std::unique_ptr<Geometry::MeshPart<MeshType>> _fbm_meshpart_inside;
      /// the meshpart representing the FBM interface region
      std::unique_ptr<Geometry::MeshPart<MeshType>> _fbm_meshpart_interface;
      /// unit-filter assembler on inside mesh-part
      std::unique_ptr<Assembly::UnitFilterAssembler<MeshType_>> _unit_asm_inside;
      /// unit-filter assembler on interface mesh-part
      std::unique_ptr<Assembly::UnitFilterAssembler<MeshType_>> _unit_asm_interface;
      /// already compiled?
      bool _compiled;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] mesh_
       * A \resident reference to the underlying mesh that the FBM is to be assembled on
       */
      explicit StokesFBMAssembler(const MeshType& mesh_) :
        _mesh(mesh_),
        _fbm_masks(),
        _fbm_meshpart_inside(),
        _fbm_meshpart_interface(),
        _unit_asm_inside(),
        _unit_asm_interface(),
        _compiled(false)
      {
        // allocate FBM mask vectors
        for(int i(0); i <= shape_dim; ++i)
          _fbm_masks.at(std::size_t(i)).resize(_mesh.get_num_entities(i), 0);
      }

      // no move, no copy, no problems
      StokesFBMAssembler(const StokesFBMAssembler&) = delete;
      StokesFBMAssembler(StokesFBMAssembler&&) = delete;
      StokesFBMAssembler& operator=(const StokesFBMAssembler&) = delete;
      StokesFBMAssembler& operator=(StokesFBMAssembler&&) = delete;

      /// virtual destructor
      virtual ~StokesFBMAssembler()
      {
      }

      /**
       * \brief Clears the assembler
       */
      void clear()
      {
        for(int i(0); i <= shape_dim; ++i)
        {
          _fbm_masks.at(std::size_t(i)).clear();
          _fbm_masks.at(std::size_t(i)).resize(_mesh.get_num_entities(i), 0);
        }
        _fbm_meshpart_inside.reset();
        _fbm_meshpart_interface.reset();
        _unit_asm_inside.reset();
        _unit_asm_interface.reset();
        _compiled = false;
      }

      /**
       * \brief Returns a reference to an internal FBM mask vector
       *
       * \param[in] dim
       * The dimension of the mask vector that is to be returned.
       *
       * \returns
       * A const reference to the internal FBM mask vector
       */
      const std::vector<int>& get_fbm_mask_vector(int dim) const
      {
        XASSERT((dim >= 0) && (dim <= shape_dim));
        return _fbm_masks.at(std::size_t(dim));
      }

      /**
       * \brief Adds a FBM mesh-part to the assembler
       *
       * \param[in] fbm_part
       * A \transient reference to the FBM mesh-part to be added.
       */
      void add_fbm_meshpart(const Geometry::MeshPart<MeshType>& fbm_part)
      {
        XASSERTM(!_compiled, "assembler has already been compiled!");
        for(int i(0); i <= shape_dim; ++i)
          _add_fbm_target_set(_get_target_set(fbm_part.get_target_set_holder(), i), i);
      }

      /**
       * \brief Synchronizes the assembler over all processes in a parallel simulation
       *
       * \attention This is a collective function that must be called by all processes simultaneously!
       *
       * \param[in] mesh_node
       * A \transient reference to the mesh node containing the communication halos of this process's neighbors
       *
       * \param[in] comm
       * The distributed communicator that corresponds to the mesh node
       */
      void sync(const Geometry::RootMeshNode<MeshType>& mesh_node, const Dist::Comm& comm)
      {
        XASSERTM(!_compiled, "assembler has already been compiled!");

        // get halo map
        const auto& halo_map = mesh_node.get_halo_map();
        const std::size_t num_halos = halo_map.size();

        static constexpr std::size_t dim1(shape_dim + 1);

        // create buffers for all halos [halo][idx]
        std::vector<std::array<std::vector<int>, dim1>> send_bufs, recv_bufs;
        send_bufs.resize(num_halos);
        recv_bufs.resize(num_halos);

        // create request vectors
        Dist::RequestVector send_reqs, recv_reqs;
        send_reqs.reserve(num_halos*dim1);
        recv_reqs.reserve(num_halos*dim1);

        // loop over all halos
        auto it = halo_map.begin();
        for(Index h(0); h < num_halos; ++h, ++it)
        {
          // loop over all dimensions
          for(std::size_t d(0); d < dim1; ++d)
          {
            // get halo target set and its size
            const Geometry::TargetSet& halo_trg_set = _get_target_set(it->second->get_target_set_holder(), int(d));
            const Index n = halo_trg_set.get_num_entities();

            // nothing to do?
            if(n == Index(0))
              continue;

            // resize buffers
            send_bufs[h][d].resize(n);
            recv_bufs[h][d].resize(n);
            int* sbuf = send_bufs[h][d].data();

            // post receive request
            recv_reqs.push_back(comm.irecv(recv_bufs[h][d].data(), n, it->first));

            // gather send data via mirror
            for(Index i(0); i < n; ++i)
              sbuf[i] = _fbm_masks[d][halo_trg_set[i]];

            // send
            send_reqs.push_back(comm.isend(sbuf, n, it->first));
          }
        }

        // wait for all receives to finish
        recv_reqs.wait_all();

        // loop over all halos
        it = halo_map.begin();
        for(Index h(0); h < num_halos; ++h, ++it)
        {
          // loop over all dimensions
          for(std::size_t d(0); d < dim1; ++d)
          {
            // get halo target set and its size
            const Geometry::TargetSet& halo_trg_set = _get_target_set(it->second->get_target_set_holder(), int(d));
            const Index n = halo_trg_set.get_num_entities();

            // nothing to do?
            if(n == Index(0))
              continue;

            int* rbuf = recv_bufs[h][d].data();

            // scatter receive data via mirror
            for(Index i(0); i < n; ++i)
              _fbm_masks[d][halo_trg_set[i]] |= rbuf[i];
          }
        }

        // wait for all sends to finish
        send_reqs.wait_all();
      }

      /**
       * \brief Compiles the assembler
       */
      void compile()
      {
        XASSERTM(!_compiled, "assembler has already been compiled!");

        // process the mask
        _process_shapes(_fbm_masks, _mesh.get_index_set_holder());

        // count masked entities
        Index counts_i[shape_dim+1], counts_o[shape_dim+1];
        for(int d(0); d <= shape_dim; ++d)
        {
          Index c0(0u), c1(0u);
          for(int x : _fbm_masks.at(std::size_t(d)))
          {
            //c0 += Index( x       & 1); // count number of bit 0 set to 1
            //c1 += Index((x >> 1) & 1); // count number of bit 1 set to 1
            c0 += Index(x == 1);
            c1 += Index(x == 3);
          }
          counts_i[d] = c0;
          counts_o[d] = c1;
        }

        // allocate mesh-part
        this->_fbm_meshpart_interface.reset(new Geometry::MeshPart<MeshType>(counts_i, false));
        this->_fbm_meshpart_inside.reset(new Geometry::MeshPart<MeshType>(counts_o, false));

        // compute target sets
        for(int d(0); d <= shape_dim; ++d)
        {
          _build_tbm_target_set(_get_target_set(_fbm_meshpart_interface->get_target_set_holder(), d), d, 1);
          _build_tbm_target_set(_get_target_set(_fbm_meshpart_inside->get_target_set_holder(), d), d, 3);
        }

        // set up unit filter assemblers
        this->_unit_asm_interface.reset(new Assembly::UnitFilterAssembler<MeshType>());
        this->_unit_asm_inside.reset(new Assembly::UnitFilterAssembler<MeshType>());
        _unit_asm_interface->add_mesh_part(*this->_fbm_meshpart_interface);
        _unit_asm_inside->add_mesh_part(*this->_fbm_meshpart_inside);

        // done
        _compiled = true;
      }

      /// Returns a reference to the mesh-part representing the FBM interface
      const Geometry::MeshPart<MeshType>& get_meshpart_interface() const
      {
        XASSERTM(_compiled, "FBM assembler must be compiled first!");
        return *_fbm_meshpart_interface;
      }

      /// Returns a reference to the mesh-part representing the FBM inside region
      const Geometry::MeshPart<MeshType>& get_meshpart_inside() const
      {
        XASSERTM(_compiled, "FBM assembler must be compiled first!");
        return *_fbm_meshpart_inside;
      }

      /// Returns a reference to the internal unit-filter assembler on the FBM interface
      const Assembly::UnitFilterAssembler<MeshType>& get_unit_asm_interface() const
      {
        XASSERTM(_compiled, "FBM assembler must be compiled first!");
        return *_unit_asm_interface;
      }

      /// Returns a reference to the internal unit-filter assembler inside the FBM region
      const Assembly::UnitFilterAssembler<MeshType>& get_unit_asm_inside() const
      {
        XASSERTM(_compiled, "FBM assembler must be compiled first!");
        return *_unit_asm_inside;
      }

      /**
       * \brief Assembles a characteristic vector on the FBM region
       *
       * \param[out] vector
       * A \transient reference to the characteristic vector that is to be assembled
       *
       * \param[in] space
       * A \transient reference to the space that the characteristic vector is to be assembled for
       *
       * \param[in] val_outside
       * The vector value that is to be used outside the FBM interface/inside regions; defaults to 0
       *
       * \param[in] val_inside
       * The vector value that is to be used inside the FBM region; defaults to 1
       *
       * \param[in] val_intface
       * The vector value that is to be on the FBM interface; defaults to 1
       */
      template<typename DT_, typename IT_, typename Space_>
      void assemble_characteristic_vector(LAFEM::DenseVector<DT_, IT_>& vector, const Space_& space,
        const DT_ val_outside = DT_(0), const DT_ val_inside = DT_(1), const DT_ val_intface = DT_(1)) const
      {
        XASSERTM(_compiled, "FBM assembler must be compiled first!");
        if(vector.empty())
          vector = LAFEM::DenseVector<DT_, IT_>(space.get_num_dofs());

        LAFEM::UnitFilter<DT_, IT_> filter_i(space.get_num_dofs());
        LAFEM::UnitFilter<DT_, IT_> filter_o(space.get_num_dofs());
        _unit_asm_interface->assemble(filter_i, space);
        _unit_asm_inside->assemble(filter_o, space);

        filter_i.get_filter_vector().format(val_intface);
        filter_o.get_filter_vector().format(val_inside);

        vector.format(val_outside);
        filter_i.filter_sol(vector);
        filter_o.filter_sol(vector);
      }

      /**
       * \brief Assembles a UnitFilter inside the FBM region
       *
       * \param[inout] filter
       * A \transient reference to the filter to be assembled
       *
       * \param[in] space
       * A \transient reference to the space to assemble the filter for
       */
      template<typename DT_, typename IT_, typename Space_>
      void assemble_inside_filter(LAFEM::UnitFilter<DT_, IT_>& filter, const Space_& space) const
      {
        XASSERTM(_compiled, "FBM assembler must be compiled first!");
        this->_unit_asm_inside->assemble(filter, space);
      }

      /**
       * \brief Assembles a UnitFilterBlocked inside the FBM region
       *
       * \param[inout] filter
       * A \transient reference to the filter to be assembled
       *
       * \param[in] space
       * A \transient reference to the space to assemble the filter for
       */
      template<typename DT_, typename IT_, typename Space_, int dim_>
      void assemble_inside_filter(LAFEM::UnitFilterBlocked<DT_, IT_, dim_>& filter, const Space_& space) const
      {
        XASSERTM(_compiled, "FBM assembler must be compiled first!");
        this->_unit_asm_inside->assemble(filter, space);
      }

      /**
       * \brief Assembles a UnitFilter on the FBM interface (non-parallel version)
       *
       * \attention
       * This function can only be used in serial (non-parallel) simulations! If you intent to use FBM in a simulation
       * that uses Global::Matrix containers, then you need to call the overload that takes Global::Matrix parameters
       * for the A and M matrices, respectively, because of the synchronization that is necessary to compute the
       * correct filter entries, which itself requires access to the global matrix functionality!
       *
       * \param[inout] filter
       * A \transient reference to the filter to be assembled
       *
       * \param[in] space
       * A \transient reference to the space to assemble the filter for
       *
       * \param[in] matrix_a
       * A \transient reference to the system matrix block A that the FBM filter will be applied to later on
       *
       * \param[in] matrix_m
       * A \transient reference to the mass matrix that will be used to apply the FBM filter later on
       */
      template<typename DT_, typename IT_, typename Space_>
      void assemble_interface_filter(LAFEM::UnitFilter<DT_, IT_>& filter, const Space_& space,
        const LAFEM::SparseMatrixCSR<DT_, IT_>& matrix_a, const LAFEM::SparseMatrixCSR<DT_, IT_>& matrix_m) const
      {
        XASSERTM(_compiled, "FBM assembler must be compiled first!");
        this->_unit_asm_interface->assemble(filter, space);

        auto diag_a = matrix_a.create_vector_l();
        auto diag_m = matrix_m.create_vector_l();
        matrix_a.extract_diag(diag_a, true);
        matrix_m.extract_diag(diag_m, true);
        filter.get_filter_vector().format(diag_a.max_abs_element() / diag_m.max_abs_element());
      }

      /**
       * \brief Assembles a UnitFilterBlocked on the FBM interface (non-parallel version)
       *
       * \attention
       * This function can only be used in serial (non-parallel) simulations! If you intent to use FBM in a simulation
       * that uses Global::Matrix containers, then you need to call the overload that takes Global::Matrix parameters
       * for the A and M matrices, respectively, because of the synchronization that is necessary to compute the
       * correct filter entries, which itself requires access to the global matrix functionality!
       *
       * \param[inout] filter
       * A \transient reference to the filter to be assembled
       *
       * \param[in] space
       * A \transient reference to the space to assemble the filter for
       *
       * \param[in] matrix_a
       * A \transient reference to the system matrix block A that the FBM filter will be applied to later on
       *
       * \param[in] matrix_m
       * A \transient reference to the mass matrix that will be used to apply the FBM filter later on
       */
      template<typename DT_, typename IT_, typename Space_, int bh_, int bw_>
      void assemble_interface_filter(LAFEM::UnitFilterBlocked<DT_, IT_, bh_>& filter, const Space_& space,
        const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& matrix_a, const LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>& matrix_m) const
      {
        XASSERTM(_compiled, "FBM assembler must be compiled first!");
        this->_unit_asm_interface->assemble(filter, space);

        auto diag_a = matrix_a.create_vector_l();
        auto diag_m = matrix_m.create_vector_l();
        matrix_a.extract_diag(diag_a);
        matrix_m.extract_diag(diag_m);
        filter.get_filter_vector().format(diag_a.max_abs_element() / diag_m.max_abs_element());
      }

      /**
       * \brief Assembles a UnitFilter on the FBM interface (parallel version)
       *
       * \attention This is a collective function that must be called by all processes simultaneously!
       *
       * \param[inout] filter
       * A \transient reference to the filter to be assembled
       *
       * \param[in] space
       * A \transient reference to the space to assemble the filter for
       *
       * \param[in] matrix_a
       * A \transient reference to the global system matrix block A that the FBM filter will be applied to later on
       *
       * \param[in] matrix_m
       * A \transient reference to the global mass matrix that will be used to apply the FBM filter later on
       */
      template<typename DT_, typename IT_, typename Space_, typename Mirror_>
      void assemble_interface_filter(LAFEM::UnitFilter<DT_, IT_>& filter, const Space_& space,
        const Global::Matrix<LAFEM::SparseMatrixCSR<DT_, IT_>, Mirror_, Mirror_>& matrix_a,
        const Global::Matrix<LAFEM::SparseMatrixCSR<DT_, IT_>, Mirror_, Mirror_>& matrix_m) const
      {
        XASSERTM(_compiled, "FBM assembler must be compiled first!");
        this->_unit_asm_interface->assemble(filter, space);

        auto diag_a = matrix_a.create_vector_l();
        auto diag_m = matrix_m.create_vector_l();
        matrix_a.extract_diag(diag_a, true);
        matrix_m.extract_diag(diag_m, true);
        filter.get_filter_vector().format(diag_a.max_abs_element() / diag_m.max_abs_element());
      }

      /**
       * \brief Assembles a UnitFilterBlocked on the FBM interface (parallel version)
       *
       * \attention This is a collective function that must be called by all processes simultaneously!
       *
       * \param[inout] filter
       * A \transient reference to the filter to be assembled
       *
       * \param[in] space
       * A \transient reference to the space to assemble the filter for
       *
       * \param[in] matrix_a
       * A \transient reference to the global system matrix block A that the FBM filter will be applied to later on
       *
       * \param[in] matrix_m
       * A \transient reference to the global mass matrix that will be used to apply the FBM filter later on
       */
      template<typename DT_, typename IT_, typename Space_, typename Mirror_, int bh_, int bw_>
      void assemble_interface_filter(LAFEM::UnitFilterBlocked<DT_, IT_, bh_>& filter, const Space_& space,
        const Global::Matrix<LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>, Mirror_, Mirror_>& matrix_a,
        const Global::Matrix<LAFEM::SparseMatrixBCSR<DT_, IT_, bh_, bw_>, Mirror_, Mirror_>& matrix_m) const
      {
        XASSERTM(_compiled, "FBM assembler must be compiled first!");
        this->_unit_asm_interface->assemble(filter, space);

        auto diag_a = matrix_a.create_vector_l();
        auto diag_m = matrix_m.create_vector_l();
        matrix_a.extract_diag(diag_a, true);
        matrix_m.extract_diag(diag_m, true);
        filter.get_filter_vector().format(diag_a.max_abs_element() / diag_m.max_abs_element());
      }

    protected:
      /// auxiliary function: add a target set to the FBM mask vector
      void _add_fbm_target_set(const Geometry::TargetSet& target_set, int dim)
      {
        std::vector<int>& mask = _fbm_masks.at(std::size_t(dim));
        const Index n = target_set.get_num_entities();
        for(Index i(0); i < n; ++i)
          mask.at(target_set[i]) |= 3;
      }

      /// auxiliary function: build target set from mask vector
      void _build_tbm_target_set(Geometry::TargetSet& target_set, int dim, int mask_value)
      {
        const std::vector<int>& mask = _fbm_masks.at(std::size_t(dim));

        for(Index i(0), k(0); i < mask.size(); ++i)
        {
          if(mask[i] == mask_value)
          {
            target_set[k] = i;
            ++k;
          }
        }
      }

      /// auxiliary function: get a target set from a target set holder
      template<typename Shape_>
      static const Geometry::TargetSet& _get_target_set(const Geometry::TargetSetHolder<Shape_>& tsh, int dim)
      {
        XASSERT(Shape_::dimension >= dim);
        typedef typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType FacetType;
        if(Shape_::dimension == dim)
          return tsh.template get_target_set<Shape_::dimension>();
        else
          return _get_target_set(static_cast<const Geometry::TargetSetHolder<FacetType>&>(tsh), dim);
      }

      /// auxiliary function: get a target set from a target set holder (recursion end overload)
      static const Geometry::TargetSet& _get_target_set(const Geometry::TargetSetHolder<Shape::Vertex>& tsh, int dim)
      {
        XASSERT(dim == 0);
        return tsh.template get_target_set<0>();
      }

      /// auxiliary function: get a target set from a target set holder
      template<typename Shape_>
      static Geometry::TargetSet& _get_target_set(Geometry::TargetSetHolder<Shape_>& tsh, int dim)
      {
        XASSERT(Shape_::dimension >= dim);
        typedef typename Shape::FaceTraits<Shape_, Shape_::dimension - 1>::ShapeType FacetType;
        if(Shape_::dimension == dim)
          return tsh.template get_target_set<Shape_::dimension>();
        else
          return _get_target_set(static_cast<Geometry::TargetSetHolder<FacetType>&>(tsh), dim);
      }

      /// auxiliary function: get a target set from a target set holder (recursion end overload)
      static Geometry::TargetSet& _get_target_set(Geometry::TargetSetHolder<Shape::Vertex>& tsh, int dim)
      {
        XASSERT(dim == 0);
        return tsh.template get_target_set<0>();
      }

      /**
       * \brief Auxiliary function: update bit 1 of the mask vector for a shape based on a bit 0 of its face mask vector
       *
       * This function sets bit 1 of any shape entity to 0 if there is at least one face adjacent to that shape entity
       * whose bit 0 is equal to 0.
       *
       * \param[inout] shape_mask
       * The mask vector for the shape that whose bit 1 mask is to be updated
       *
       * \param[in] face_masj
       * The mask vector for the face that whose bit 0 mask is to be used for the update of the shape mask
       *
       * \param[in] faces_at_shape
       * The faces-at-shape index set of the underlying mesh which describes the faces that are adjacent to a shape
       */
      template<int m_>
      static void _process_mask(std::vector<int>& shape_mask, const std::vector<int>& face_mask, const Geometry::IndexSet<m_>& faces_at_shape)
      {
        XASSERT(Index(shape_mask.size()) == faces_at_shape.get_num_entities());
        XASSERT(Index(face_mask.size()) == faces_at_shape.get_index_bound());

        const Index n = faces_at_shape.get_num_entities();
        const int m = faces_at_shape.get_num_indices();
        for(Index i(0); i < n; ++i)
        {
          // bit 1 of shape mask is only 1 if bit 0 of all faces (of all dimensions) are 1
          int bit1 = (shape_mask[i] >> 1) & 1;
          for(int j(0); j < m; ++j)
            bit1 &= (face_mask[faces_at_shape(i,j)] & 1);
          shape_mask[i] = (shape_mask[i] & 1) | (bit1 << 1);
        }
      }

      /// auxiliary function: process masks for all face dimensions for a given shape
      template<typename Shape_, int face_dim_, std::size_t n_>
      static void _process_masks(std::array<std::vector<int>, n_>& masks, const Geometry::IndexSetWrapper<Shape_, face_dim_>& idx_wrapper)
      {
        _process_masks(masks, static_cast<const Geometry::IndexSetWrapper<Shape_, face_dim_-1>&>(idx_wrapper));
        _process_mask(masks.at(Shape_::dimension), masks.at(std::size_t(face_dim_)), idx_wrapper.template get_index_set<face_dim_>());
      }

      /// auxiliary function: process masks for all face dimensions for a given shape (recursion end)
      template<typename Shape_, std::size_t n_>
      static void _process_masks(std::array<std::vector<int>, n_>& masks, const Geometry::IndexSetWrapper<Shape_, 0>& idx_wrapper)
      {
        _process_mask(masks.at(Shape_::dimension), masks.front(), idx_wrapper.template get_index_set<0>());
      }

      /// auxiliary function: process masks shape dimensions
      template<typename Shape_, std::size_t n_>
      static void _process_shapes(std::array<std::vector<int>, n_>& masks, const Geometry::IndexSetHolder<Shape_>& idx_holder)
      {
        typedef typename Shape::FaceTraits<Shape_, Shape_::dimension-1>::ShapeType FacetType;
        _process_shapes(masks, static_cast<const Geometry::IndexSetHolder<FacetType>&>(idx_holder));
        _process_masks(masks, idx_holder.template get_index_set_wrapper<Shape_::dimension>());
      }

      /// auxiliary function: process masks shape dimensions (recursion end)
      template<std::size_t n_>
      static void _process_shapes(std::array<std::vector<int>, n_>&, const Geometry::IndexSetHolder<Shape::Vertex>&)
      {
        // end of recursion; nothing to do here
      }
    }; // class StokesFBMAssembler
  } // namespace Assembly
} // namespace FEAT
