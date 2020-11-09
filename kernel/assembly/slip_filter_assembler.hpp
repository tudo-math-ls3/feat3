// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#ifndef KERNEL_ASSEMBLY_SLIP_FILTER_ASSEMBLER_HPP
#define KERNEL_ASSEMBLY_SLIP_FILTER_ASSEMBLER_HPP 1

// includes, FEAT
#include <kernel/assembly/base.hpp>
#include <kernel/assembly/unit_filter_assembler.hpp>
#include <kernel/lafem/slip_filter.hpp>
#include <kernel/lafem/sparse_vector_blocked.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>
#include <kernel/space/lagrange1/element.hpp>
#include <kernel/space/lagrange2/element.hpp>

namespace FEAT
{
  namespace Assembly
  {

    /// \cond internal
    namespace Intern
    {
      template<typename Space_, int shape_dim_ = Space_::shape_dim>
      class Lagrange2InterpolatorWrapper;
    }
    /// \endcond

    /**
     * \brief Helper class that computes a weighted outer unit normal field
     *
     * \tparam Trafo_
     * Type of the transformation
     *
     * The transformation is a Finite Element function itself. Globally, it does not have enough regularity, but in
     * each cell (or rather, on each boundary facet), a local normal field can be computed. This can be interpolated
     * or projected onto a global space with more regularity so it can be used in computations.
     *
     * At the moment, only the standard P1/Q1 trafo is implemented, meaning the outer normal field is constant in
     * each boundary facet and discontinuous across facet boundaries. One method of mapping this function onto P1/Q1
     * on the boundary is to define the outer unit normal field pointwise in each mesh vertex by calculating the mean
     * of the normals of the surrounding facets, weighted by the facet's volume. The same can be done for P2/Q2
     * transformations, which are not implemented yet.
     *
     * In some circumstances, this leads to error estimates of optimal order, see \cite BD99
     *
     * This weighted mean outer unit normal is used in the SlipFilter.
     *
     * \see LAFEM::SlipFilter, SlipFilterAssembler, TraceAssembler
     *
     * \author Jordi Paul
     *
     */
    template<typename Trafo_>
    class OuterNormalComputer
    {
      public:
        /// Type for the transformation
        typedef Trafo_ TrafoType;
        /// The underlying mesh tupe
        typedef typename TrafoType::MeshType MeshType;
        /// Shape of the mesh's cells
        typedef typename MeshType::ShapeType ShapeType;
        /// Shape of the facets that make up the boundary we want to compute the normals for,
        // i.e. Simplex<2> faces for a Simplex<3> mesh
        typedef typename Shape::FaceTraits<ShapeType, ShapeType::dimension - 1>::ShapeType FacetType;

        /// The shape dimension of the boundary faces
        static constexpr int facet_dim = ShapeType::dimension-1;
        /// The number of coordinates for each vertex
        static constexpr int world_dim = MeshType::world_dim;

        /// The number of vertices that make up a Facet, i.e. 4 for a Hypercube<2>
        static constexpr int nvt_loc = TrafoType::num_coeff_facet;

        /// Floating point type
        typedef typename MeshType::CoordType CoordType;

        /**
         * \brief Computes facet orientations
         *
         * We need to figure out how a cell is oriented wrt. the reference cell. With the orientation code that we
         * get from the CongruencySampler, we know if the facets of the real cell are oriented like the corresponding
         * facets of the reference cell. Their orientation is known, so we can compute the orientation of the real
         * facets from all this.
         *
         * \param[in,out] orientation
         * Orientations in floating point format, -1.0 and 1.0 resp. for negatively or positively oriented facets
         *
         * \param[in] facets
         * Marker if a facet of the mesh is part of the boundary considered
         *
         * \param[in] mesh
         * The mesh of which we consider a part of the boundary
         *
         * orientation and facets are numbered like the facets in the original mesh (and not the corresponding
         * MeshPart for the boundary) because they need to be accessed from there.
         *
         */
        static void compute_orientations(CoordType* orientation, const Index* facets, const MeshType& mesh)
        {
          // Various index sets for the parent
          auto& idx_vert_at_facet(mesh.template get_index_set<facet_dim, 0>());
          auto& idx_vert_at_shape(mesh.template get_index_set<ShapeType::dimension, 0>());
          auto& idx_facet_at_shape(mesh.template get_index_set<ShapeType::dimension, facet_dim>());

          // This will contain the local vertex indices in the reference cell numbering
          typename MeshType::template IndexSet<facet_dim, 0>::Type::IndexTupleType reference_numbering;
          // This will contain the local vertex indices as they are numbered in the mesh
          typename MeshType::template IndexSet<facet_dim, 0>::Type::IndexTupleType my_numbering;

          // Go check all cells in the mesh
          for(Index k(0); k < mesh.get_num_entities(ShapeType::dimension); ++k)
          {
            // Check all facets
            for(int l(0); l < Shape::FaceTraits<ShapeType, facet_dim>::count; ++l)
            {
              // Index of the l-th facet in the parent mesh
              Index parent_l(idx_facet_at_shape[k][l]);

              // If the facet is in the meshpart ...
              if(facets[parent_l] != ~Index(0))
              {
                // ... check the numbering
                for(int i(0); i < nvt_loc; ++i)
                {
                  // The FaceIndexMapping gives us the reference cell numbering of the facet
                  int j(Geometry::Intern::FaceIndexMapping<ShapeType, facet_dim, 0>::map(l,i));
                  reference_numbering[i] = idx_vert_at_shape[k][j];

                  // This is how the facet is numbered in the mesh
                  my_numbering[i] = idx_vert_at_facet[parent_l][i];
                }

                // Get the orientation code. This describes if the facet given by my_numbering is a rotation of the
                // reference facet or a rotation of an inverted reference facet
                int orientation_code(
                  Geometry::Intern::CongruencySampler<FacetType>::compare(reference_numbering, my_numbering));

                // The orientation of the facet is the (possibly inverted) orientation of the appropriate facet in
                // the reference cell
                orientation[parent_l] = CoordType(
                  Geometry::Intern::CongruencySampler<FacetType>::orientation(orientation_code)
                  *Shape::ReferenceCell<ShapeType>::orientation(l));
              }
            } // facets
          } // cells
        } // compute_orientations

        /**
         * \brief Computes an outer unit normal field
         *
         * The outer unit normal is computed in all vertices of the mesh's facets that are identified by facets.
         * Because the real outer unit normal is discontinuous across facets, this is a mean where each contribution
         * is weighted by the corresponding facet's volume.
         *
         * \tparam VectorType_
         * Type for holding the value vector for a filter, so this will be a LAFEM::SparseVector.
         *
         * \tparam TrafoType_
         * Transformation, at the moment there is just the standard trafo.
         *
         * \param[in] nu
         * This will hold the weighted mean outer unit normal in the vertices
         *
         * \param[in] facets
         * All relevant facets are marked in this array
         *
         * \param[in] orientation
         * The orientation for all faces
         *
         * \param[in] trafo
         * The transformation
         *
         * orientation and facets are numbered like the facets in the original mesh (and not the corresponding
         * MeshPart for the boundary) because they need to be accessed from there.
         *
         */
        template<typename VectorType_, typename TrafoType_>
        static void compute_outer_unit_normal(
          VectorType_& nu,
          const Index* facets,
          const CoordType* orientation,
          const TrafoType_& trafo)
          {
            typedef typename VectorType_::DataType DataType;

            nu.format(DataType(0));
            // Vertex at facet index set from the parent
            auto& idx(trafo.get_mesh().template get_index_set<facet_dim,0>());
            // Temporary vector for holding one normal at a time
            Tiny::Vector<DataType, world_dim> tmp(DataType(0));

            // For every cell in the meshpart, compute the outer normal vectors in all local Lagrange points of the
            // element's transformation
            for(Index k(0); k < trafo.get_mesh().get_num_entities(facet_dim); ++k)
            {
              if(facets[k] != ~Index(0) )
              {
                CoordType vol(trafo.template compute_vol<FacetType>(k));

                // Compute all normals
                Tiny::Matrix<DataType, nvt_loc, world_dim> nu_loc(trafo.compute_oriented_normals(k));

                // Add the local contributions to the vector that is numbered according to the parent
                for(int l(0); l < nvt_loc; ++l)
                {
                  tmp = nu(idx[k][l]) + nu_loc[l]*DataType(orientation[k]*vol);
                  nu(idx[k][l],  tmp);
                }
              }
            }

            //// Normalize nu
            //for(Index i(0); i < nu.used_elements(); ++i)
            //{
            //  Index j(nu.indices()[i]);
            //  tmp = nu(j);
            //  tmp.normalize();

            //  nu(j,tmp);
            //}
          } // compute_outer_unit_normal

    }; // struct OuterNormalComputer

    /**
     * \brief Contains routines to assemble slip filters
     *
     * \tparam Mesh_
     * The mesh type.
     *
     * The assemble() routine is implemented for the homogenous slip condition for Lagrange 1/2 elements only. For
     * the many other limitations, see the SlipFilter class documentation.
     *
     * \see SlipFilter
     *
     */
    template<typename Mesh_>
    class SlipFilterAssembler
    {
      public:
        /// The mesh type
        typedef Mesh_ MeshType;
        /// Floating point type
        typedef typename Mesh_::CoordType CoordType;
        /// Type for the mapping of slip boundary entities to mesh entities
        typedef typename Geometry::TargetSetHolder<typename Mesh_::ShapeType> TargetSetHolderType;
        /// shape dimension
        static constexpr int facet_dim = Mesh_::shape_dim-1;
        /// The number of coordinates for each vertex
        static constexpr int world_dim = Mesh_::world_dim;

      private:
        /// Facets for the slip condition get marked by this
        Index* _facets;
        /// For every facet in the mesh, this is its orientation in the mesh
        CoordType* _orientation;
        /// Flag when to recompute _facets and _orientation (i.e. if MeshParts get added)
        bool recompute;
        /// Number of entities of various shape dimensions in the slip boundary represented by this filter
        Index _num_entities[Mesh_::shape_dim+1];
        /// Mapping of boundary entities to mesh entieties (i.e. vertices, edges, facets)
        TargetSetHolderType* _target_set_holder;

      public:
        /**
         * \brief Constructor from mesh
         *
         * \param[in] mesh
         * The mesh to assemble the slip filter for
         *
         * We need to know the number of vertices in the mesh for the final sizes of _facets and _orientation, or
         * some sort of size management has to be implemented
         *
         */
        explicit SlipFilterAssembler(const MeshType& mesh) :
          _facets(new Index[mesh.get_num_entities(facet_dim)]),
          _orientation(new CoordType[mesh.get_num_entities(facet_dim)]),
          recompute(false),
          _target_set_holder(nullptr)
          {
            for(Index i(0); i < mesh.get_num_entities(facet_dim); ++i)
            {
              _facets[i] = ~Index(0);
              _orientation[i] = CoordType(0);
            }

            for(int i(0); i < Mesh_::shape_dim+1; ++i)
              _num_entities[i] = Index(0);
          }

        /// Explicitly delete the copy constructor
        SlipFilterAssembler(const SlipFilterAssembler&) = delete;
        /// Explicitly delete the move constructor
        SlipFilterAssembler(SlipFilterAssembler&&) = delete;
        /// Explicitly delete the move assignment operator
        SlipFilterAssembler& operator=(SlipFilterAssembler&&) = delete;

        /**
         * \brief Virtual destructor
         */
        virtual ~SlipFilterAssembler()
        {
          if(_facets != nullptr)
            delete[] _facets;
          if(_orientation != nullptr)
            delete[] _orientation;
          if(_target_set_holder != nullptr)
            delete _target_set_holder;
        }

        /**
         * \brief Adds a MeshPart to the assembler
         *
         * This means that the boundary where the slip condition is enforced gets enlarged by the MeshPart.
         *
         * \param[in] meshpart
         * MeshPart identifying a part of the boundary where to enforce the slip condition
         *
         */
        void add_mesh_part(const Geometry::MeshPart<MeshType>& meshpart)
        {
          // Adding a MeshPart changes the structures
          recompute = true;
          const auto& facet_target(meshpart.template get_target_set<facet_dim>());

          for(Index l(0); l < facet_target.get_num_entities(); l++)
          {
            // Add the new facet if it is not already there
            if(_facets[facet_target[l]] == ~Index(0))
            {
              _facets[facet_target[l]] = _num_entities[facet_dim];
              _num_entities[facet_dim]++;
            }
          }
        }

        /**
         * \brief Recomputes the target mapping for geometric entities
         *
         * \param[in] mesh
         * The slip boundary this filter represents is part of this mesh.
         *
         */
        void recompute_target_set_holder(const Mesh_& mesh)
        {
          // We need a whole new TargetSetHolder
          if(_target_set_holder != nullptr)
            delete _target_set_holder;

          // We know how many facets we have from add_mesh_part
          _target_set_holder = new TargetSetHolderType(_num_entities);

          auto& facet_targets = _target_set_holder->template get_target_set<facet_dim>();

          // Copy the facet information from the marker array
          for(Index l(0); l < mesh.get_num_entities(facet_dim); ++l)
          {
            if(_facets[l] != ~Index(0))
              facet_targets[_facets[l]] = l;
          }

          // Now we deduct the vertices/edges/whatever at the slip boundary from the facet target mapping
          Geometry::Intern::template TargetSetComputer<facet_dim, 0>::top_to_bottom(
            *_target_set_holder, mesh.get_index_set_holder());

          // Update num_entities information from the possibly modified _target_set_holder
          for(int d(0); d <= facet_dim; ++d)
            _num_entities[d] = _target_set_holder->get_num_entities(d);
        }

        /**
         * \brief Assembles a homogeneous SlipFilter for a Lagrange 1 FE space
         *
         * \tparam Trafo_
         * The transformation
         *
         * \tparam MemType_
         * Memory type for the filter, i.e. Mem::Main or Mem::CUDA
         *
         * \tparam DataType_
         * Data type for the filter
         *
         * \tparam IndexType_
         * Index type for the filter
         *
         * \param[in,out] filter
         * The filter to be assembled
         *
         * \param[in,out] space
         * The filter gets applied of functions of this FE space
         */
        template<typename Trafo_, typename MemType_, typename DataType_, typename IndexType_>
        void assemble(LAFEM::SlipFilter<MemType_, DataType_, IndexType_, world_dim>& filter, const Space::Lagrange1::Element<Trafo_>& space)
        {
          // Allocate the filter if necessary
          if(filter.get_nu().size() == Index(0))
          {
            filter = LAFEM::SlipFilter<MemType_, DataType_, IndexType_, world_dim>
              (space.get_trafo().get_mesh().get_num_entities(0), space.get_num_dofs());
          }

          LAFEM::SlipFilter<Mem::Main, DataType_, IndexType_, world_dim> buffer;
          buffer.convert(filter);

          // Compute orientations if necessary
          if(recompute)
          {
            OuterNormalComputer<Trafo_>::compute_orientations(_orientation, _facets, space.get_trafo().get_mesh());
          }

          // Compute the weighted outer unit normal field
          OuterNormalComputer<Trafo_>::compute_outer_unit_normal(
            buffer.get_nu(), _facets, _orientation, space.get_trafo());

          // Generate the filter vector. For the Lagrange 1 with standard trafo case, this is just a clone operation
          buffer.get_filter_vector().clone(buffer.get_nu());

          // Upload assembled result
          filter.convert(buffer);
        }

        /**
         * \brief Assembles a homogeneous SlipFilter for a Lagrange 2 FE space
         *
         * \tparam Trafo_
         * The transformation
         *
         * \tparam MemType_
         * Memory type for the filter, i.e. Mem::Main or Mem::CUDA
         *
         * \tparam DataType_
         * Data type for the filter
         *
         * \tparam IndexType_
         * Index type for the filter
         *
         * \param[in,out] filter
         * The filter to be assembled
         *
         * \param[in,out] space
         * The filter gets applied of functions of this FE space
         *
         * This implements a primite P1/Q1 to P2/Q2 interpolation of the outer unit normal vector field. It would be
         * nicer to have a real inter FE space mapping, but that would require the application of node functionals
         * to FE functions. This is not implemented yet due to regularity issues (i.e. application of a Q1 node
         * functional requires point evaluation at mesh vertices, but Q1~ functions are not continuous so this is not
         * possible).
         *
         */
        template<typename Trafo_, typename MemType_, typename DataType_, typename IndexType_>
        void assemble(LAFEM::SlipFilter<MemType_, DataType_, IndexType_, world_dim>& filter, const Space::Lagrange2::Element<Trafo_>& space)
        {
          typedef Space::Lagrange2::Element<Trafo_> SpaceType;

          // The mesh the FE space lives on
          const auto& mesh = space.get_trafo().get_mesh();

          // Allocate the filter if necessary
          if(filter.get_nu().size() == Index(0))
          {
            filter = LAFEM::SlipFilter<MemType_, DataType_, IndexType_, world_dim>
              (mesh.get_num_entities(0), space.get_num_dofs());
          }

          // Always assemble the filter in Mem::Main and convert it to prevent excessive GPU up/downloading
          LAFEM::SlipFilter<Mem::Main, DataType_, IndexType_, world_dim> buffer;
          buffer.convert(filter);

          // Compute orientations if necessary
          if(recompute)
          {
            recompute_target_set_holder(space.get_trafo().get_mesh());
            OuterNormalComputer<Trafo_>::compute_orientations(_orientation, _facets, mesh);
          }

          // Compute the weighted outer unit normal field
          OuterNormalComputer<Trafo_>::compute_outer_unit_normal(
            buffer.get_nu(), _facets, _orientation, space.get_trafo());

          // We only have something to do if the filter is not empty after recompute_target_set_holder()
          if(buffer.get_nu().used_elements() > Index(0))
          {
            Intern::Lagrange2InterpolatorWrapper<SpaceType>::project(
              buffer.get_filter_vector(), buffer.get_nu(), space, _target_set_holder);
          }
          // Upload assembled result
          filter.convert(buffer);

        }

    }; // class SlipFilterAssembler

    /// \cond internal
    namespace Intern
    {
      /*
       * This is a hackish implementation of a Lagrange 1 to Lagrange 2 interpolator and mostly copied from
       * interpolator.hpp. The critical thing is that it requires the application of node functionals of one space
       * to FE functions from another space, so in general this will not be well-defined (like applying Q2 node
       * functionals (which are point evaluations) to Q1~ functions (which are not continuous). In the case of
       * conforming Lagrange elements, it is well-defined however, and implemented below.
       */
      template<typename Vector_, typename Space_, int shape_dim_>
      struct Lagrange2InterpolatorCore
      {
        public:
          template< typename TSH_>
          static void project(Vector_& lagrange_2_vector, const Vector_& lagrange_1_vector, const Space_& space, const TSH_* tsh)
          {
            typedef typename Vector_::DataType DataType;
            typedef typename Vector_::ValueType ValueType;

            // define dof assignment
            typedef typename Space_::template DofAssignment<shape_dim_, DataType>::Type DofAssignType;
            DofAssignType dof_assign(space);

            // skip empty dof assignments
            if(dof_assign.get_max_assigned_dofs() < 1)
              return;

            // Vertex at shape index set
            auto& idx(space.get_trafo().get_mesh().template get_index_set<shape_dim_, 0>());

            // loop over all entities
            if(tsh == nullptr)
            {
              const Index num_shapes = space.get_mesh().get_num_entities(shape_dim_);
              for(Index shape(0); shape < num_shapes; ++shape)
              {
                ValueType tmp(DataType(0));

                // evaluate the "node functional"
                for(int j(0); j < idx.get_num_indices(); ++j)
                  tmp += lagrange_1_vector(idx(shape, j));

                tmp *= DataType(1)/DataType(idx.get_num_indices());
                tmp.normalize();

                // prepare dof-assignment
                dof_assign.prepare(shape);

                // loop over all contributions
                Index dof_index(dof_assign.get_index(0));
                lagrange_2_vector(dof_index, tmp);

                // finish
                dof_assign.finish();
              }
            }
            else
            {
              auto& ts = tsh->template get_target_set<shape_dim_>();

              // loop over all entities
              const Index num_shapes = ts.get_num_entities();

              for(Index shape(0); shape < num_shapes; ++shape)
              {
                ValueType tmp(DataType(0));

                // evaluate the "node functional"
                for(int j(0); j < idx.get_num_indices(); ++j)
                  tmp += lagrange_1_vector(idx(ts[shape], j));

                tmp *= DataType(1)/DataType(idx.get_num_indices());
                tmp.normalize();

                // prepare dof-assignment
                dof_assign.prepare(ts[shape]);

                // loop over all contributions
                Index dof_index(dof_assign.get_index(0));
                lagrange_2_vector(dof_index, tmp);

                // finish
                dof_assign.finish();
              }
            }
          }
      };

      template<typename Vector_, typename Space_>
      struct Lagrange2InterpolatorCore<Vector_, Space_, 0>
      {
        public:
          template<typename TSH_>
          static void project(Vector_& lagrange_2_vector, const Vector_& lagrange_1_vector, const Space_& space, const TSH_* tsh)
          {
            typedef typename Vector_::DataType DataType;
            typedef typename Vector_::ValueType ValueType;

            // define dof assignment
            typedef typename Space_::template DofAssignment<0, DataType>::Type DofAssignType;
            DofAssignType dof_assign(space);

            if(tsh == nullptr)
            {
              // loop over all entities
              const Index num_verts = space.get_mesh().get_num_entities(0);
              for(Index vert(0); vert < num_verts; ++vert)
              {
                // prepare dof-assignment
                dof_assign.prepare(vert);

                // "loop over all contributions"
                Index dof_index(dof_assign.get_index(0));

                // evaluate the "node functional"
                ValueType tmp(lagrange_1_vector(vert));
                lagrange_2_vector(dof_index, tmp);

                // finish
                dof_assign.finish();
              }
            }
            else
            {
              auto& ts = tsh->template get_target_set<0>();
              // loop over all entities
              const Index num_verts = ts.get_num_entities();
              for(Index vert(0); vert < num_verts; ++vert)
              {
                // prepare dof-assignment
                dof_assign.prepare(ts[vert]);

                // "loop over all contributions"
                Index dof_index(dof_assign.get_index(0));

                // evaluate the "node functional"
                ValueType tmp(lagrange_1_vector(ts[vert]));
                lagrange_2_vector(dof_index, tmp);

                // finish
                dof_assign.finish();
              }
            }
          }
      };

      template<typename Space_, int shape_dim_>
      class Lagrange2InterpolatorWrapper
      {
        public:
          template<typename Vector_, typename TSH_>
          static void project(Vector_& lagrange_2_vector, const Vector_& lagrange_1_vector, const Space_& space, const TSH_* tsh)
          {
            // recurse down
            Lagrange2InterpolatorWrapper<Space_, shape_dim_ - 1>::project(lagrange_2_vector, lagrange_1_vector, space, tsh);

            // call interpolator core
            Lagrange2InterpolatorCore<Vector_, Space_, shape_dim_>::project(lagrange_2_vector, lagrange_1_vector, space, tsh);
          }
      };

      template<typename Space_>
      class Lagrange2InterpolatorWrapper<Space_, 0>
      {
        public:
          template<typename Vector_, typename TSH_>
          static void project(Vector_& lagrange_2_vector, const Vector_& lagrange_1_vector, const Space_& space, const TSH_* tsh)
          {
            // call interpolator core
            Lagrange2InterpolatorCore<Vector_, Space_, 0>::project(lagrange_2_vector, lagrange_1_vector, space, tsh);
          }
      };

    } // namespace Intern
    /// \endcond


  } //namespace Assembly
} // namespace FEAT

#endif // KERNEL_ASSEMBLY_SLIP_FILTER_ASSEMBLER_HPP
