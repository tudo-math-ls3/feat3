#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_SURFACE_MESH_HPP
#define KERNEL_GEOMETRY_ATLAS_SURFACE_MESH_HPP 1

#include <kernel/geometry/atlas/chart.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/index_calculator.hpp>

#include <kernel/util/math.hpp>

#include <deque>

namespace FEAT
{
  namespace Geometry
  {
    namespace Atlas
    {

      /**
       * Traits for the SurfaceMesh class
       */
      struct SurfaceMeshTraits
      {
        /// No explicit map is available in general
        static constexpr bool is_explicit = false;
        /// We support implicit projection
        static constexpr bool is_implicit = true;
        /// This is a world_dim dimensional object
        static constexpr int world_dim = 3;
        /// If there was a parametrisation, it would be the object's shape dim
        static constexpr int param_dim = 2;
      }; // struct SurfaceMeshTraits

      /**
       * \brief Boundary description by a surface mesh in 3d
       *
       * \tparam Mesh_
       * Type for the mesh this boundary description refers to
       *
       * \author Jordi Paul
       *
       */
      template<typename Mesh_>
      class SurfaceMesh :
        public ChartCRTP<SurfaceMesh<Mesh_>, Mesh_, SurfaceMeshTraits>
      {
      public:
        /// The CRTP base class
        typedef ChartCRTP<SurfaceMesh<Mesh_>, Mesh_, SurfaceMeshTraits> BaseClass;
        /// Floating point type for coordinates
        typedef typename BaseClass::CoordType CoordType;
        /// Vector type for world points
        typedef typename BaseClass::WorldPoint WorldPoint;
        /// The shape type is always Simplex<2>
        typedef Shape::Simplex<2> ShapeType;
        /// The type for the mesh defining the discrete surface
        typedef ConformalMesh<ShapeType, 3, 3, CoordType> SurfaceMeshType;
        /// Pointer to the surface mesh object
        SurfaceMeshType* _surface_mesh;

      private:
        /// Marker for edges
        // For computing point projections, we need to know if an edge was already traversed by adding the neighbour
        // across said edge to the search stack. If the surface is concave, it might happen that triangle A says
        // "go to my neighbour, B", and then B says "go to my neighbour, A". It's then clear that the wanted point
        // lies on the edge between A and B.
        bool* _edge_marker;
        /// Marker for triangles
        bool* _cell_marker;

      public:
        /// Explicitly delete empty default constructor
        SurfaceMesh() = delete;

        /**
         * \brief Constructor getting an object
         */
        explicit SurfaceMesh(SurfaceMeshType* surf_mesh) :
          _surface_mesh(surf_mesh),
          _edge_marker(new bool[_surface_mesh->get_num_entities(SurfaceMeshType::shape_dim-1)]),
          _cell_marker(new bool[_surface_mesh->get_num_entities(SurfaceMeshType::shape_dim)])
        {
        }

        /// Explicitly delete move constructor
        SurfaceMesh(SurfaceMeshType&&) = delete;

        /**
         * \brief Virtual destructor
         */
        virtual ~SurfaceMesh()
        {
          if(_surface_mesh != nullptr)
            delete _surface_mesh;

          delete[] _edge_marker;
          delete[] _cell_marker;
        }

        /// \copydoc ChartBase::bytes
        virtual std::size_t bytes() const override
        {
          if(_surface_mesh != nullptr)
            return _surface_mesh->bytes();
          else
            return std::size_t(0);
        }

        /** \copydoc ChartBase::get_type() */
        virtual String get_type() const override
        {
          return "SurfaceMesh";
        }

        /** \copydoc ChartBase::write */
        virtual void write(std::ostream& os, const String& sindent) const override
        {
          String indent(sindent);

          os << indent << "<SurfaceMesh";
          os << " verts=\" " << _surface_mesh->get_num_entities(0) << "\"";
          os << " trias=\" " << _surface_mesh->get_num_entities(ShapeType::dimension) << "\"";
          os << ">" << std::endl;

          // increase indent
          indent.resize(indent.size()+2, ' ');

          const auto& vtx = _surface_mesh->get_vertex_set();

          os << indent << "<Vertices>" << std::endl;
          indent.resize(indent.size()+2, ' ');

          for(Index i(0); i < vtx.get_num_vertices(); ++i)
          {
            const auto& v = vtx[i];
            os << indent << v[0];
            for(int j(1); j < SurfaceMeshType::world_dim; ++j)
              os << ' ' << v[j];
            os << std::endl;
          }

          indent.resize(indent.size()-2);
          os << indent << "</Vertices>" << std::endl;

          // Get vertex at cell index set
          const auto& idx = _surface_mesh->template get_index_set<ShapeType::dimension, 0>();

          os << indent << "<Triangles>" << std::endl;
          indent.resize(indent.size()+2, ' ');

          for(Index i(0); i < idx.get_num_entities(); ++i)
          {
            const auto& idx_loc = idx[i];
            os << indent << idx_loc[0];
            for(int j(1); j < idx.num_indices; ++j)
              os << ' ' << idx_loc[j];
            os << std::endl;
          }

          indent.resize(indent.size()-2);
          os << indent << "</Triangles>" << std::endl;

          indent.resize(indent.size()-2);
          os << indent << "</SurfaceMesh>";

        }

        /**
         * \brief Projects a single point to the surface given by the surface mesh
         *
         * \param[in,out] point
         * The point that gets projected.
         */
        void project_point(WorldPoint& point) const
        {
          CoordType signed_dist(0);
          WorldPoint grad_dist(CoordType(0));

          project_point(point, signed_dist, grad_dist);

        }

        /**
         * \brief Projects a single point to the surface given by the surface mesh
         *
         * \param[in,out] point
         * The point that gets projected.
         *
         * \param[in,out] signed_distance
         * The signed distance of the original to the projected point.
         *
         * \param[in,out] grad_dist
         * Gradient of the distance function.
         */
        void project_point(WorldPoint& point, CoordType& signed_distance, WorldPoint& grad_dist) const
        {
          // There is nothing to do if the surface mesh does not have any cells
          if(_surface_mesh->get_num_entities(SurfaceMeshType::shape_dim) == Index(0))
            return;

          // Vertex at cell information
          const auto& idx(_surface_mesh->template get_index_set<SurfaceMeshType::shape_dim, 0>());
          // The mesh's vertex set so we can get at the coordinates
          const auto& vtx(_surface_mesh->get_vertex_set());

          // When a point is found, this will hold the coefficients of its image under the projection wrt. the
          // standard P1/Q1 transformation
          Tiny::Vector<CoordType, SurfaceMeshType::shape_dim+1> coeffs(CoordType(0));
          // This will hold the coefficients of the transformation mapping
          Tiny::Vector<CoordType, SurfaceMeshType::shape_dim+1> bary(CoordType(0));

          // Since we want to project the point an arbitrary distance, we set this to somethin huge
          CoordType acceptable_distance(Math::huge<CoordType>());
          // Find the facet in the SurfaceMesh we project to and compute the coefficients
          Index best_facet =
            find_cell(coeffs, point, Index(0), acceptable_distance, _surface_mesh, true);

          // Clip the coefficients
          for(int j(0); j < SurfaceMeshType::shape_dim; ++j)
          {
            coeffs[j] = Math::max(coeffs[j], CoordType(0));
            coeffs[j] = Math::min(coeffs[j], CoordType(1));
          }
          // Now that we have the surface mesh cell with the lowest distance, we can compute the projected point
          WorldPoint projected_point(CoordType(0));

          // Evaluate the surface mesh trafo for the computed coefficients
          CoordType coeff0(CoordType(1));
          for(int j(0); j < SurfaceMeshType::shape_dim; ++j)
          {
            // Index of the local vertex i in the SurfaceMesh
            Index i(idx(best_facet, Index(j+1)));
            projected_point += coeffs[j]*vtx[i];
            coeff0 -= coeffs[j];
          }

          projected_point += coeff0*vtx[idx(best_facet, Index(0))];

          grad_dist = (projected_point - point);
          signed_distance = grad_dist.norm_euclid();

          // If the distance is too small, we set the gradient vector to zero
          if(signed_distance < Math::eps<CoordType>())
            grad_dist.format(CoordType(0));
          else
          {
            grad_dist.normalise();
            WorldPoint nu(get_normal_on_tria(best_facet, coeffs));
            signed_distance *= Math::signum(Tiny::dot(nu, grad_dist));
          }

          point = projected_point;

        }

        /**
         * \brief Computes the outer normal on a triangle for a single point
         *
         * \param[in] facet
         * Number of the triangle.
         *
         * \param[in] coeffs
         * Coefficients of the point where we want to compute the normal.
         *
         * \returns
         * The outer (unit) normal vector at the given point.
         */
        WorldPoint get_normal_on_tria(
          const Index facet, const Tiny::Vector<CoordType, SurfaceMeshType::shape_dim+1> DOXY(coeffs)) const
        {
          // Vertex at cell information
          const auto& idx(_surface_mesh->template get_index_set<SurfaceMeshType::shape_dim, 0>());
          // The mesh's vertex set so we can get at the coordinates
          const auto& vtx(_surface_mesh->get_vertex_set());

          Tiny::Matrix<CoordType, 2, SurfaceMeshType::world_dim> coords(CoordType(0));
          Tiny::Matrix<CoordType, SurfaceMeshType::world_dim, 2> coords_transpose(CoordType(0));

          Index i0(idx(facet, Index(0)));
          Index i1(idx(facet, Index(1)));
          Index i2(idx(facet, Index(2)));

          coords[0] = vtx[i1] - vtx[i0];
          coords[1] = vtx[i2] - vtx[i0];

          for(int i(0); i < coords.m; ++i)
          {
            for(int j(0); j < coords.n; ++j)
              coords_transpose(j,i) = coords(i,j);
          }

          WorldPoint nu(Tiny::orthogonal(coords_transpose));
          nu.normalise();

          return nu;
        }

        /**
         * \brief Orthogonally projects all vertices at facets of a MeshPart
         *
         * \param[in,out] mesh
         * The mesh the Meshpart refers to and whose vertices are to be projected
         *
         * \param[in] meshpart
         * The meshpart identifying the boundary of the mesh that is to be projected
         *
         * The MeshPart has to contain facets and a topology for them. To keep the search as local as possible,
         * the MeshPart is traversed facet-wise and after finishing the current facet, its neighbours are added to
         * the search stack.
         *
         */
        void project_meshpart(Mesh_& mesh, const MeshPart<Mesh_>& meshpart) const
        {
          // We need the topology of the meshpart because searching for points is cell-wise and we need to know
          // about the vertices in each cell
          if(meshpart.get_topology() == nullptr)
            throw InternalError(__func__,__FILE__,__LINE__,
            "The meshpart needs a topology for SurfaceMesh::project()!");

          Index num_facets(meshpart.get_num_entities(SurfaceMeshType::shape_dim));

          // There is nothing to do if the meshpart does not have any cells (which refer to facets of the original
          // mesh)
          if(num_facets == Index(0))
            return;

          // The number of vertices that need to be projected
          const Index num_verts(meshpart.get_num_entities(0));
          // Mapping of vertices from the meshpart to the real mesh
          const auto& ts_verts(meshpart.template get_target_set<0>());

          // Vertex at surface mesh cell information
          const auto& idx_sm(_surface_mesh->template get_index_set<SurfaceMeshType::shape_dim, 0>());
          // Vertex set of the surface mesh
          const auto& vtx_sm(_surface_mesh->get_vertex_set());

          // Vertex at facet information for the meshpart
          const auto& idx_mp(meshpart.template get_index_set<SurfaceMeshType::shape_dim, 0>());
          // Subfacet at facet information for the meshpart. i.e. if Mesh is a Hypercube<3> mesh, the boundary is a
          // Hypercube<2> mesh, the facets are quads and THEIR subfacets are edges, so this is edge@quad
          // Needed for computing neighbours.
          const auto& facet_idx_mp(meshpart.template get_index_set
            <SurfaceMeshType::shape_dim, SurfaceMeshType::shape_dim-1>());

          // Neighbour information for all facets in the MeshPart. Needed for adding facets to the search stack
          // so we can exploit local information
          typename MeshPart<Mesh_>::template IndexSet<SurfaceMeshType::shape_dim, SurfaceMeshType::shape_dim-1>::Type
            neigh_mp(meshpart.get_num_entities(SurfaceMeshType::shape_dim));

          Geometry::Intern::FacetNeighbours::compute(neigh_mp, facet_idx_mp);

          // The mesh's vertex set, since we modify the coordinates by projection
          auto& vtx(mesh.get_vertex_set());

          // When a point is found, this will hold the coefficients of its image under the projection wrt. the
          // standard P1/Q1 transformation
          Tiny::Vector<CoordType, SurfaceMeshType::shape_dim+1> coeffs(CoordType(0));

          // For every vertex, we need to know if it was already projected
          bool* vert_todo(new bool[num_verts]);
          for(Index i(0); i < num_verts; ++i)
            vert_todo[i] = true;

          // For every facet, we need to know if it was already added to the search stack
          bool* facet_todo(new bool[num_facets]);
          for(Index i(0); i < num_facets; ++i)
            facet_todo[i] = true;

          // For every facet we keep a guess, a cell number in the surface mesh where we start looking for vertices
          // from our facet.
          Index* guesses(new Index[num_facets]);
          for(Index i(0); i < num_facets; ++i)
            guesses[i] = Index(0);

          // The search stack for all facets. We add neighbours to it to keep the search as local as possible
          std::deque<Index> facet_stack;
          facet_stack.push_front(0);

          // Treat all facets
          while(!facet_stack.empty())
          {
            // This is the facet from the mesh part that's being worked on
            Index current_facet(facet_stack.front());
            // This will be the facet from the surface mesh were the vertex was found, set it to the best guess
            // initially
            Index facet_sm(guesses[current_facet]);

            facet_stack.pop_front();
            facet_todo[current_facet] = false;

            // We want to orthogonally project a point, so we need to figure out an acceptable distance for this
            // projection. Otherwise, it might happen that we project to the far side of the domain
            CoordType acceptable_distance(CoordType(1));

            // Compute acceptable distance. Because we do not have the transformation available, we cannot compute
            // i.e. the volume of the facet to deduce an acceptable distance from that. So instead, we pick a
            // random vertex of the facet, compute the distances to all other points, and take the geometric mean.
            for(int j(0); j < idx_mp.num_indices-1; ++j)
            {
              Index i0(idx_mp(current_facet, Index(0)));
              Index i_mp(idx_mp(current_facet, Index(j+1)));
              acceptable_distance *= ((vtx[i_mp] - vtx[i0]).norm_euclid());
            }
            acceptable_distance = Math::pow(
              acceptable_distance, CoordType(1)/CoordType(idx_mp.num_indices) - CoordType(1));

            // Find all vertices in the current facet
            for(int j(0); j < idx_mp.num_indices; ++j)
            {
              // The index of the local vertex j in the Meshpart
              Index i_mp(idx_mp(current_facet, Index(j)));

              if(vert_todo[i_mp])
              {
                vert_todo[i_mp] = false;
                const auto& x(vtx[ts_verts[i_mp]]);

                // Find the facet in the SurfaceMesh we project to and compute the coefficients. We break at the
                // first facet fulfilling the criteria
                facet_sm = find_cell(coeffs, x, guesses[current_facet], acceptable_distance,
                _surface_mesh, false);
                // Update the guess for this facet because there might be more points to project
                guesses[current_facet] = facet_sm;

                // This is the index of the local vertex j in the real mesh
                Index i_parent(ts_verts[i_mp]);
                // Clear the old vertex coordinates so we can just add onto that
                vtx[i_parent].format(CoordType(0));
                // Evaluate the surface mesh trafo for the computed coefficients
                CoordType coeff0(CoordType(1));
                for(int i(0); i < SurfaceMeshType::shape_dim; ++i)
                {
                  // Index of the local vertex i in the SurfaceMesh
                  Index i_sm(idx_sm(facet_sm, Index(i+1)));
                  vtx[i_parent] += coeffs[i]*vtx_sm[i_sm];
                  coeff0 -= coeffs[i];
                }
                vtx[i_parent] += coeff0*vtx_sm[idx_sm(facet_sm, Index(0))];

              } // vert_todo[i_mp]
            } // vertices in current facet

            // Now add all neighbours that are still on the todo list to the search stack
            for(int l(0); l < neigh_mp.num_indices; ++l)
            {
              Index k(neigh_mp(current_facet, Index(l)));
              if(k != ~Index(0) && facet_todo[k])
              {
                facet_todo[k] = false;
                facet_stack.push_back(k);
                // The last vertex was found in facet_sm, and since k is a neighbour of current, it cannot be too
                // far off
                guesses[k] = facet_sm;
              }
            } // adding neighbours

          } // search stack

          // Clean up
          delete[] facet_todo;
          delete[] vert_todo;
          delete[] guesses;
        }

        /// \copydoc ChartBase::dist()
        CoordType compute_dist(const WorldPoint& point) const
        {
          WorldPoint grad_dist(CoordType(0));
          return compute_dist(point, grad_dist);
        }

        /// \copydoc ChartBase::dist()
        CoordType compute_dist(const WorldPoint& point, WorldPoint& grad_dist) const
        {
          CoordType signed_distance(0);
          WorldPoint projected_point(point);
          project_point(projected_point, signed_distance, grad_dist);

          return Math::abs(signed_distance);
        }

        /// \copydoc ChartBase::signed_dist()
        CoordType compute_signed_dist(const WorldPoint& point) const
        {
          WorldPoint grad_dist(CoordType(0));
          return compute_signed_dist(point, grad_dist);
        }

        /// \copydoc ChartBase::signed_dist()
        CoordType compute_signed_dist(const WorldPoint& point, WorldPoint& grad_dist) const
        {
          CoordType signed_distance(0);
          WorldPoint projected_point(point);
          project_point(projected_point, signed_distance, grad_dist);

          return signed_distance;
        }

      private:
        /**
         * \brief Finds the cell in the surface mesh that contains a specific point
         *
         * \tparam world_dim
         * World dimension, meaning the number of coordinates.
         *
         * \tparam sc_
         * Stride for the coefficients vector.
         *
         * \tparam sx_
         * Stride for the vector containing the point.
         *
         * \param[out] coeffs
         * The coefficients that map x to its projection under the standard transformation.
         *
         * \param[in] x
         * The point we are looking for.
         *
         * \param[in] guess
         * The facet in the surface mesh we start looking for x in.
         *
         * \param[in] acceptable_distance
         * The maximum distance we want to move in orthogonal direction, so we do not project to the far side of
         * the domain by accident.
         *
         * \param[in] mesh
         * The (surface) mesh we are looking for x in.
         *
         * \param[in,out] traversed
         * Information about which facet was already traversed. No information is passed to or from this function,
         * this is just so it does not need to be reallocated for every vertex.
         *
         * For every subfacet, we need to know if it was already traversed by adding the neighbour across said
         * facet to the search stack. If the surface is concave, it might happen that facet A says "go to my
         * neighbour, B", and then B says "go to my neighbour, A". It's then clear that the wanted point lies on
         * the facet between A and B (c.f. determine_success()).
         *
         * \see determine_success
         * \returns The index of the cell in which x was found
         */
        template<int world_dim, int sc_, int sx_>
        Index find_cell(
          Tiny::Vector<CoordType, SurfaceMeshType::shape_dim+1, sc_>& coeffs,
          const Tiny::Vector<CoordType, world_dim, sx_>& x, Index guess, CoordType acceptable_distance,
          const SurfaceMeshType* mesh, bool find_best_approx) const
        {
          // For every subfacet, we need to know if it was already traversed by adding the neighbour across said
          // facet to the search stack. If the surface is concave, it might happen that facet A says "go to my
          // neighbour, B", and then B says "go to my neighbour, A". It's then clear that the wanted point lies on
          // the facet between A and B.

          // Clear traversed information
          for(Index facet(0); facet < mesh->get_num_entities(SurfaceMeshType::shape_dim-1); ++facet)
            _edge_marker[facet] = false;

          Index num_cells(mesh->get_num_entities(SurfaceMeshType::shape_dim));
          // Clear cell todo information
          for(Index cell(0); cell < num_cells; ++cell)
            _cell_marker[cell] = false;

          // Type of the vertex at cell index set
          typedef typename SurfaceMeshType::template IndexSet<SurfaceMeshType::shape_dim, 0>::Type VertAtCellIdxType;
          // Number of vertices per cell
          static constexpr int num_vert_loc = VertAtCellIdxType::num_indices;
          // Vertex at cell information
          const VertAtCellIdxType& idx(mesh->template get_index_set<SurfaceMeshType::shape_dim, 0>());
          // Facet at cell information
          const auto& facet_idx(mesh->template
              get_index_set<SurfaceMeshType::shape_dim, SurfaceMeshType::shape_dim-1>());
          // The mesh's vertex set so we can get at the coordinates
          const auto& vtx(mesh->get_vertex_set());

          const auto& neigh(mesh->get_neighbours());

          // This will contain all vertices making up a cell in the mesh. Since the mesh is of shape
          // Simplex<shape_dim>, there are shape_dim+1 vertices.
          Tiny::Matrix<CoordType, num_vert_loc, world_dim> coords(CoordType(0));

          // This is the search stack
          std::deque<Index> cell_stack;
          // We start with our guess
          cell_stack.push_front(guess);
          // This will be the cell we are currently working on
          Index current(~Index(0));

          bool success(false);

          Index best_approx_cell(~Index(0));

          CoordType best_approx_dist(Math::huge<CoordType>());

          Tiny::Vector<CoordType, SurfaceMeshType::shape_dim+1, sc_> current_coeffs(CoordType(0));

          // Search all cells
          while(!cell_stack.empty())
          {
            // Get new cell
            current = cell_stack.front();
            cell_stack.pop_front();
            _cell_marker[current] = true;

            // Collect all vertices that make up the cell in the coords matrix
            for(int j(0); j < SurfaceMeshType::shape_dim+1; ++j)
            {
              Index i(idx(current, Index(j)));
              coords[j] = vtx[i];
            }

            // Compute the coefficients of the inverse coordinate transformation. The last coefficient is the
            // signed distance of x to the facet made up by coords
            compute_inverse_coeffs(current_coeffs, x, coords);

            // From the coefficients, we can compute the barycentric coordinates (as the surface mesh is of Simplex
            // shape) so we get an idea about in which direction we need to search first
            Tiny::Vector<CoordType, SurfaceMeshType::shape_dim+1> bary(CoordType(1));
            CoordType ortho_dist = compute_bary_and_dist(bary, current_coeffs);

            // From all this, we can now determine if there is an orthogonal projection of x onto this facet. For
            // this, we need to know the traversed information as well and the neighbour information of the current
            // facet
            success = determine_success(ortho_dist, acceptable_distance, bary, facet_idx[current], _edge_marker);

            // Ok, maybe we are done
            if(success && !find_best_approx)
            {
              best_approx_cell = current;
              coeffs = current_coeffs;
              break;
            }
            // Check if we need to update the best approx data
            else if(success && find_best_approx)
            {
              if(ortho_dist < best_approx_dist)
              {
                best_approx_dist = ortho_dist;
                best_approx_cell = current;
                coeffs = current_coeffs;

              }
            }

            // If we are not done, add all neighbours that have not been searched to the search stack
            for(int l(0); l < SurfaceMeshType::shape_dim+1; ++l)
            {
              // The neighbour over the l-th facet
              Index n(neigh(current, Index(l)));
              // Check if the neighbour exists and is still on the to do list
              if(n != ~ Index(0) && !_cell_marker[n])
              {
                Index facet(facet_idx(current,Index(l)));
                // If the neighbour is still on the to do list and the facet it lies at is in the traversed list,
                // this means the barycentric coordinate for the corresponding direction was negative, meaning
                // that x lies in that direction. So we add that neighbour to the front.
                if(_edge_marker[facet])
                  cell_stack.push_front(n);
                else
                  cell_stack.push_back(n);

                // Since we added the neighbour to the stack, prevent it from being added again later.
                _cell_marker[n] = true;
              }
            } // adding neighbours

          } // search loop

          // If we get here without success being true, the point could not be found
          if(best_approx_cell == ~Index(0))
            throw InternalError(__func__,__FILE__,__LINE__,"Could not find point "+stringify(x)+" in SurfaceMesh");

          // Otherwise, it was found in current and coeffs is still set from that evaluation
          return best_approx_cell;
        } // find_cell()

        /**
         * \brief For given reference cell coordinates, this computes the barycentric coordinates and the distance
         *
         * \tparam DT_
         * The floating point type
         *
         * \tparam n_
         * World dimension of the points to compute the inverse mapping for
         *
         * \tparam sn_
         * Stride for the vector of coefficients
         *
         * \tparam sb_
         * Stride for the vector of barycentric coordinates
         *
         * \param[out] bary
         * The vector of barycentric coordinates.
         *
         * \param[in] coeffs
         * The coefficients for the mapping from the reference cell to the real cell
         *
         * \returns
         * The absolute value of the last coefficient, which is the absolute distance.
         *
         * Assume we have a non-degenerate \c Simplex<s> called \f$ S \f$ in \f$ \mathbb{R}^d \f$ defined by
         * vertices \f$ x^j \in \mathbb{R}^d, j = 0, \dots, s \f$. Assume for now that \f$ s = d \f$.
         *
         * Then
         * \f[
         *   \forall x \in \mathbb{R}^d: x = x^0 + \sum_{j=1}^s \lambda_j (x^j - x^0)
         * \f]
         * and \f$ \lambda_j, j=1,\dots,s \f$ with \f$ \lambda_0 := 1 - sum_{j=1}^s \f$ are called the
         * <b> barycentric coordinates of \f$ x \f$ wrt. S </b>. Note that the above equation can be rewritten as
         * \f[
         *   \forall x \in \mathbb{R}^d: x = \sum_{j=0}^s \lambda_j x^j.
         * \f]
         *
         * It is easy to see that
         * \f[
         *   x \in S \Leftrightarrow  \forall j = 0, \dots, d: \lambda_j \in [0, 1] .
         * \f]
         * If \f$ \exists j \in \{ 0, \dots, s \}: \lambda_j < 0 \f$, then \f$ x \notin S \f$ and \f$ x \f$ lies
         * on the far side of the plane defined by the facet opposite vertex \f$ j \f$. This makes the barycentric
         * coordinates very handy for finding out in which direction of a given simplex a point lies.
         *
         * This routine assumes that \f$ s = d-1 \f$, where \f$ d = 2, 3 \f$ and that the last coefficient is
         * actually the signed distance of the point \f$ x \f$ to which the coefficients belong to the hyperplane
         * defined by \f$ S \f$.
         *
         */
        template<typename DT_, int n_, int sn_, int sb_>
        static DT_ compute_bary_and_dist(Tiny::Vector<DT_, n_, sb_>& bary, const Tiny::Vector<DT_, n_, sn_>& coeffs)
        {
          bary(0) = CoordType(1);
          for(int j(0); j < n_-1; ++j)
          {
            bary(0) -= coeffs(j);
            bary(j+1) = coeffs(j);
          }

          return Math::abs(coeffs(n_-1));
        }

        /**
         * \brief Determines the success of the point search
         *
         * \tparam DT_
         * Floating point type for coordinates
         *
         * \tparam n_
         * Number of barycentric coordinates, meaning shape_dim+1
         *
         * \tparam sn_
         * Stride for the vector of barycentric coordinates
         *
         * \tparam IdxType_
         * Type of the IndexSet of appropriate size
         *
         * \param[in] ortho_dist
         * Orthogonal distance of the point to the appropriate plane
         *
         * \param[in] acceptable_distance
         * Tolerance for ortho_dist
         *
         * \param[in] bary
         * Barycentric coordinates of the point in the current cell
         *
         * \param[in] facet_idx
         * Facet at current cell information
         *
         * \param[in] traversed
         * Marker which facets have been traversed so far
         *
         * \returns
         * True if bary and traversed indicate that the point can be projected to the current cell.
         *
         * After the barycentric coordinates of a point wrt. a cell have been computed, we have to determine if the
         * point lies in that cell. This is definitely the case if
         * \f[
         *   \forall j=0, \dots, s: \lambda_j \in [0, 1] \wedge\mathrm{ortho_dist} = 0.
         * \f]
         *
         * If the first condition holds, the second is violated but \f$ \mathrm{ortho_dist} <
         * \mathrm{acceptable_distance} \f$, the search was successful as well.
         *
         * Because we are dealing with floats, we add some tolerance either way. The remaining problem is a concave
         * surface like this:
         *
         \verbatim

          x
          (exteriour)
          *--------------
          |            (1)
          | (0)       (interiour)
          |

         \endverbatim
         *
         * The barycentric coordinates of edge 0 will point into direction of edge 1 and vice versa. For this, we
         * keep the \c traversed information: If we first compute the barycentric coordinates wrt. edge (0), it will
         * mark the facet * (which is just a vertex in this case) as traversed. If we later computed the barycentric
         * coordinates of x wrt. edge (1), the coordinate belonging to vertex * will be negative, pointing us back
         * to edge (0). But * has been marked as traversed, so setting that negative barycentric coordinate to zero
         * means (if not other barycentric coordinate was lower than 0) we project to that facet, in this case: * .
         *
         */
        template<typename DT_, int n_, int sn_, typename IdxType_>
        static bool determine_success(
          DT_ ortho_dist,
          DT_ acceptable_distance,
          Tiny::Vector<DT_, n_, sn_>& bary,
          const IdxType_& facet_idx,
          bool* traversed)
        {
          static const DT_ tol = Math::sqrt(Math::eps<DT_>());
          bool success(true);

          // Determine success
          if(ortho_dist > acceptable_distance)
            success = false;

          int num_negative(0);

          for(int l(0); l < n_; ++l)
          {
            Index facet(facet_idx[Index(l)]);
            if(bary(l) < - tol)
            {
              ++num_negative;
              if(traversed[facet])
              {
                bary(l) = CoordType(0);
              }
              else
              {
                success = false;
                traversed[facet] = true;
              }
            }
          }

          if(num_negative > 1)
            success = false;

          return success;
        }

        /// \cond internal
        /**
         * \copydoc FEAT::Trafo::Standard::compute_inverse_coeffs()
         *
         * This Geometry. We no allowed dem tings from Trafo. So copy-paste.
         */
        template<typename DT_, int sb_, int sp_, int smx_, int snx_>
        void compute_inverse_coeffs(
          Tiny::Vector<DT_, 3, sb_>& coeffs,
          const Tiny::Vector<DT_, 3, sp_>& point,
          const Tiny::Matrix<DT_, 3, 3, smx_, snx_>& x) const
        {
          coeffs.format(DT_(0));

          // A will contain the transformation matrix
          Tiny::Matrix<DT_, 3, 3> A(DT_(0));
          // Fill the all rows in the first shape_dim = world_dim-1 columns
          for(int i(0); i < 3; ++i)
          {
            for(int j(0); j < 2; ++j)
              A(i,j) = x(j+1,i) - x(0,i);
          }

          // The last column is the additional direction for our augmented simplex and it is orthogonal to the rest
          Tiny::Vector<DT_, 3> ortho = Tiny::orthogonal(A.template size_cast<3, 2>());
          // Normalise this so the last coefficient is the signed distance
          ortho.normalise();

          // Set the last column in A
          for(int i(0); i < 3; ++i)
            A(i,2) = ortho(i);

          // This will be the inverse of the transformation matrix
          Tiny::Matrix<DT_, 3 , 3> Ainv(DT_(0));
          Ainv.set_inverse(A);

          // This is the solution of A u = point - x[0]
          coeffs = Ainv*(point - x[0]);
        }

        /**
         * \copydoc FEAT::Trafo::Standard::compute_inverse_coeffs()
         *
         * This Geometry. We no allowed dem tings from Trafo. So copy-paste.
         */
        template<typename DT_, int world_dim, int sc_, int sp_, int smx_, int snx_>
        void compute_inverse_coeffs(
          Tiny::Vector<DT_, 2, sc_>& coeffs,
          const Tiny::Vector<DT_, world_dim, sp_>& point,
          const Tiny::Matrix<DT_, 2, world_dim, smx_, snx_>& x) const
        {
          static_assert( (world_dim == 2 || world_dim == 3),
          "world dim has to be 2 or 3 for complementary barycentric coordinates");

          auto tmp = x[1]-x[0];
          DT_ sp(Tiny::dot(point - x[0],tmp));
          DT_ nsqr(Math::sqr(tmp.norm_euclid()));
          coeffs[0] = sp/nsqr;
          tmp = point - (x[0] + coeffs[0]*(x[1]-x[0]));
          coeffs[1] = tmp.norm_euclid();
        }
        /// \endcond
      }; // class SurfaceMesh

      /// \cond internal
      template<typename Mesh_>
      class SurfaceMeshVertsParser :
        public Xml::MarkupParser
      {
      public:
        typedef typename Mesh_::CoordType CoordType;
        typedef VertexSet<3, 3, CoordType> VertexSetType;

      protected:
        VertexSetType& _vertex_set;
        Index _read;

      public:
        explicit SurfaceMeshVertsParser(VertexSetType& vertex_set) :
          _vertex_set(vertex_set),
          _read(0)
        {
        }

        virtual bool attribs(std::map<String,bool>&) const override
        {
          return true;
        }

        virtual void create(int iline, const String& sline, const String&, const std::map<String, String>&, bool closed) override
        {
          if(closed)
            throw Xml::GrammarError(iline, sline, "Invalid closed markup");
        }

        virtual void close(int iline, const String& sline) override
        {
          // ensure that we have read all vertices
          if(_read < _vertex_set.get_num_vertices())
            throw Xml::GrammarError(iline, sline, "Invalid terminator; expected point");
        }

        virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String&) override
        {
          // no children allowed
          return nullptr;
        }

        virtual bool content(int iline, const String& sline) override
        {
          // make sure that we do not read more points than expected
          if(_read >= _vertex_set.get_num_vertices())
            throw Xml::ContentError(iline, sline, "Invalid content; expected terminator");

          // split line by whitespaces
          std::deque<String> scoords;
          sline.split_by_charset(scoords);

          // check size
          if(scoords.size() != std::size_t(3))
            throw Xml::ContentError(iline, sline, "Invalid number of coordinates");

          // get our current vertex
          auto& vtx = _vertex_set[_read];
          for(int i(0); i < 3; ++i)
          {
            if(!scoords.at(std::size_t(i)).parse(vtx[i]))
              throw Xml::ContentError(iline, sline, "Failed to parse vertex coordinate");
          }

          // okay, another point done
          ++_read;

          return true;
        }
      };

      template<typename Mesh_>
      class SurfaceMeshTriasParser :
        public Xml::MarkupParser
      {
      public:
        typedef IndexSet<3> IndexSetType;

      protected:
        IndexSetType& _index_set;
        Index _read;

      public:
        explicit SurfaceMeshTriasParser(IndexSetType& index_set) :
          _index_set(index_set),
          _read(0)
        {
        }

        virtual bool attribs(std::map<String,bool>&) const override
        {
          return true;
        }

        virtual void create(int iline, const String& sline, const String&, const std::map<String, String>&, bool closed) override
        {
          if(closed)
            throw Xml::GrammarError(iline, sline, "Invalid closed markup");
        }

        virtual void close(int iline, const String& sline) override
        {
          // ensure that we have read all triangles
          if(_read < _index_set.get_num_entities())
            throw Xml::GrammarError(iline, sline, "Invalid terminator; expected triangle indices");
        }

        virtual std::shared_ptr<MarkupParser> markup(int, const String&, const String&) override
        {
          // no children allowed
          return nullptr;
        }

        virtual bool content(int iline, const String& sline) override
        {
          // make sure that we do not read more points than expected
          if(_read >= _index_set.get_num_entities())
            throw Xml::ContentError(iline, sline, "Invalid content; exprected terminator");

          // split line by whitespaces
          std::deque<String> sidx;
          sline.split_by_charset(sidx);

          // check size
          if(sidx.size() != std::size_t(3))
            throw Xml::ContentError(iline, sline, "Invalid number of indices");

          // get our current vertex
          auto& idx = _index_set[_read];
          for(int i(0); i < 3; ++i)
          {
            if(!sidx.at(std::size_t(i)).parse(idx[i]))
              throw Xml::ContentError(iline, sline, "Failed to parse vertex index");
          }

          // okay, another point done
          ++_read;

          return true;
        }
      };

      template<typename Mesh_, bool enable_ = (Mesh_::shape_dim > 2)>
      class SurfaceMeshChartParser :
        public Xml::DummyParser
      {
      public:
        explicit SurfaceMeshChartParser(ChartBase<Mesh_>*&)
        {
          throw InternalError(__func__,__FILE__,__LINE__,"Thou shall not arrive here");
        }
      };

      template<typename Mesh_>
      class SurfaceMeshChartParser<Mesh_, true> :
        public Xml::MarkupParser
      {
      protected:
        typedef typename SurfaceMesh<Mesh_>::SurfaceMeshType SurfaceMeshType;
        ChartBase<Mesh_>*& _chart;
        SurfaceMeshType* _surfmesh;
        Index _nverts, _ntrias;
        bool _have_verts, _have_trias;

      public:
        explicit SurfaceMeshChartParser(ChartBase<Mesh_>*& chart) :
          _chart(chart),
          _surfmesh(nullptr),
          _nverts(0),
          _ntrias(0),
          _have_verts(false), _have_trias(false)
        {
        }

        virtual bool attribs(std::map<String,bool>& attrs) const override
        {
          attrs.emplace("verts", true);
          attrs.emplace("trias", true);
          return true;
        }

        virtual void create(
          int iline,
          const String& sline,
          const String&,
          const std::map<String, String>& attrs,
          bool closed) override
        {
          // make sure this one isn't closed
          if(closed)
            throw Xml::GrammarError(iline, sline, "Invalid closed SurfaceMesh markup");

          // try to parse the size
          if(!attrs.find("verts")->second.parse(_nverts))
            throw Xml::GrammarError(iline, sline, "Failed to parse number of vertices");
          if(!attrs.find("trias")->second.parse(_ntrias))
            throw Xml::GrammarError(iline, sline, "Failed to parse number of triangles");

          // create surface mesh
          Index num_entities[3] = {_nverts, 0, _ntrias};
          _surfmesh = new SurfaceMeshType(num_entities);
        }

        virtual void close(int iline, const String& sline) override
        {
          if(!_have_verts)
            throw Xml::GrammarError(iline, sline, "Vertices of SurfaceMesh are missing");
          if(!_have_trias)
            throw Xml::GrammarError(iline, sline, "Triangles of SurfaceMesh are missing");

          // deduct topology
          _surfmesh->deduct_topology_from_top();

          // okay
          _chart = new SurfaceMesh<Mesh_>(_surfmesh);
          _surfmesh = nullptr;
        }

        virtual bool content(int, const String&) override
        {
          return false;
        }

        virtual std::shared_ptr<Xml::MarkupParser> markup(int, const String&, const String& name) override
        {
          if(name == "Vertices")
          {
            _have_verts = true;
            return std::make_shared<SurfaceMeshVertsParser<Mesh_>>(_surfmesh->get_vertex_set());
          }
          if(name == "Triangles")
          {
            _have_trias = true;
            return std::make_shared<SurfaceMeshTriasParser<Mesh_>>(_surfmesh->template get_index_set<2,0>());
          }
          return nullptr;
        }
      };
      /// \endcond

    } // namespace Atlas
  } // namespace Geometry
} // namespace FEAT

#endif // KERNEL_GEOMETRY_ATLAS_SURFACE_MESH_HPP
