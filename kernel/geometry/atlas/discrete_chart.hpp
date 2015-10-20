#pragma once
#ifndef KERNEL_GEOMETRY_ATLAS_DISCRETE_CHART_HPP
#define KERNEL_GEOMETRY_ATLAS_DISCRETE_CHART_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/geometry/atlas/chart.hpp>
#include <kernel/geometry/conformal_mesh.hpp>
#include <kernel/geometry/mesh_streamer_factory.hpp>

#include <deque>

namespace FEAST
{
  namespace Geometry
  {
    namespace Atlas
    {
      /// \cond internal
      namespace Intern
      {
        /**
         * \brief Coefficient / facet numbering permutation
         *
         * For Simplex<d>, d != 1, facet i lies opposite of the local vertex i and nothing needs to be done
         */
        template<int shape_dim>
        inline int coeffperm(int i)
        {
          return i;
        }

        /**
         * \brief Coefficient / facet numbering permutation
         *
         * For Simplex<1>, facet i is the local vertex i, so the coefficients need to be permuted
         */
        template<>
        inline int coeffperm<1>(int i)
        {
          return 1-i;
        }
      } // namespace Intern
      /// \endcond

      /// Discrete chart traits
      template<typename MeshType>
      struct DiscreteChartTraits
      {
        /// No explicit map is available in general
        static constexpr bool is_explicit = false;
        /// We support implicit projection
        static constexpr bool is_implicit = true;
        /// This is a world_dim dimensional object
        static constexpr int world_dim = MeshType::world_dim;
        /// If there was a parametrisation, it would be the object's shape dim
        static constexpr int param_dim = MeshType::shape_dim;
      }; // struct DiscreteChartTraits

      /**
       * \brief Class template for discrete charts
       *
       * A DiscreteChart is just a (possibly very fine) mesh, that may even have a completely different shape type
       * than the mesh referring to it. It will in general never have an explicit form (a 1d DiscreteChart has it,
       * but then one could use a Polyline as well).
       *
       * \tparam Mesh_
       * Type of the (root) mesh this chart refers to.
       *
       * \tparam SurfaceMeshType
       * Type of the boundary mesh defining the DiscreteChart.
       *
       * Mesh_ could be ConformalMesh<Hypercube<3>> and SurfaceMeshType could be ConformalMesh<Simplex<2>>.
       *
       */
      template<typename Mesh_, typename SurfaceMesh_>
      class DiscreteChart:
        public ChartCRTP<DiscreteChart<Mesh_, SurfaceMesh_>, Mesh_, DiscreteChartTraits<SurfaceMesh_>>
      {
        public:
          /// The type for the mesh defining the discrete surface
          typedef SurfaceMesh_ SurfaceMeshType;
          /// The CRTP base class
          typedef ChartCRTP<DiscreteChart<Mesh_, SurfaceMesh_>, Mesh_, DiscreteChartTraits<SurfaceMesh_>> BaseClass;
          /// Floating point type for coordinates
          typedef typename BaseClass::CoordType CoordType;
          /// Vector type for world points
          typedef typename BaseClass::WorldPoint WorldPoint;

        public:
          /// Pointer to the surface mesh object
          SurfaceMesh_* _surface_mesh;

        private:
          /**
           * \brief Constructor
           * See the documentation of this class template for details.
           */
          explicit DiscreteChart(SurfaceMesh_* surface_mesh_) :
            _surface_mesh(surface_mesh_)
            {
              CONTEXT(name()+"::DiscreteChart(SurfaceMesh_*)");
            }

        public:
          /**
           * \brief Constructor
           * See the documentation of this class template for details.
           */
          explicit DiscreteChart() :
            _surface_mesh(nullptr)
            {
              CONTEXT(name()+"::DiscreteChart()");
            }

          /**
           * \brief Destructor
           */
          virtual ~DiscreteChart()
          {
            if(_surface_mesh != nullptr)
              delete _surface_mesh;
          }

          /** \copydoc ChartBase::get_type() */
          virtual String get_type() const override
          {
            CONTEXT(name()+"::get_type()");
            return "discrete";
          }

          /**
           * \brief Descriptive String
           *
           * \returns The class name as String
           */
          static String name()
          {
            return "DiscreteChart<"+Mesh_::name()+", "+SurfaceMesh_::name()+">";
          }

          /** \copydoc ChartBase::write_data_container() */
          virtual void write_data_container(MeshStreamer::ChartContainer& chart_container) const override
          {
            CONTEXT(name()+"::write_data_container()");
            Geometry::MeshWriter<SurfaceMeshType>::write_mesh(chart_container.mesh_data, *_surface_mesh);
          }

          /**
           * \brief Parses a MeshDataContainer
           *
           * \param[in] data
           * The MeshStreamer::MeshDataContainer holding the data
           *
           * \returns
           * A new object of type DiscreteChart containing the mesh data from the container
           */
          static DiscreteChart<Mesh_, SurfaceMesh_>* parse(FEAST::MeshStreamer::MeshDataContainer& data)
          {
            CONTEXT(name()+"::parse()");
            // Create a factory for our mesh
            MeshStreamerFactory<SurfaceMesh_> my_factory(&data);
            // Create our mesh
            SurfaceMeshType* surface_mesh(new SurfaceMeshType(my_factory));

            return new DiscreteChart(surface_mesh);
          }

          /**
           * \brief Projects a single point to the surface given by the discrete chart
           *
           * \warning Dumb and inefficient: This just moves the point to the nearest point of the surface mesh. This
           * means it is somewhat likely that different points get moved together and for every single point, and the
           * complexity is O(# surface_mesh_vertices) FOR EVERY SINGLE POINT PROJECTED.
           *
           * If possible, use project(Mesh_&, MeshPart&) version.
           *
           * \param[in,out] point
           * The point that gets projected.
           *
           */
          void project(WorldPoint& point) const
          {
            CONTEXT(name()+"::project(WorldPoint&)");
            Index min_index(0);
            CoordType min_dist(Math::Limits<CoordType>::max());

            typename SurfaceMesh_::VertexSetType& vtx(_surface_mesh->get_vertex_set());
            CoordType tol(Math::pow(Math::eps<CoordType>(), CoordType(0.5)));

            for(Index i(0); i < _surface_mesh->get_num_entities(0); ++i)
            {
              typename SurfaceMeshType::VertexSetType::VertexType tmp(point);
              tmp -= vtx[i];

              CoordType current_dist(tmp.norm_euclid());

              if(current_dist <= tol)
              {
                min_index = i;
                break;
              }

              if(current_dist < min_dist)
              {
                min_dist = current_dist;
                min_index = i;
              }

            }

            point = vtx[min_index];

          } // void project()

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
          void project(Mesh_& mesh, const MeshPart<Mesh_>& meshpart) const
          {
            CONTEXT(name()+"::project(Mesh_&, const MeshPart<Mesh>&)");

            // We need the topology of the meshpart because searching for points is cell-wise and we need to know
            // about the vertices in each cell
            if(meshpart.get_topology() == nullptr)
              throw InternalError("The meshpart needs a topology for DiscreteChart::project()!");

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
            const auto& idx_sm(_surface_mesh->template get_index_set<SurfaceMesh_::shape_dim, 0>());
            // Vertex set of the surface mesh
            const auto& vtx_sm(_surface_mesh->get_vertex_set());

            // Vertex at facet information for the meshpart
            const auto& idx_mp(meshpart.template get_index_set<SurfaceMeshType::shape_dim, 0>());
            // Subfacet at facet information for the meshpart. i.e. if Mesh is a Hypercube<3> mesh, the boundary is a
            // Hypercube<2> mesh, the facets are quads and THEIR subfacets are edges, so this is edge@quad
            // Needed for computing neighbours.
            const auto& facet_idx_mp(meshpart.template get_index_set
            <SurfaceMesh_::shape_dim, SurfaceMesh_::shape_dim-1>());

            // Neighbour information for all facets in the MeshPart. Needed for adding facets to the search stack
            // so we can exploit local information
            typename MeshPart<Mesh_>::template IndexSet<SurfaceMesh_::shape_dim, SurfaceMesh_::shape_dim-1>::Type
              neigh_mp(meshpart.get_num_entities(SurfaceMesh_::shape_dim));

            Geometry::Intern::FacetNeighbours::compute(neigh_mp, facet_idx_mp);

            // The mesh's vertex set, since we modify the coordinates by projection
            auto& vtx(mesh.get_vertex_set());

            // When a point is found, this will hold the coefficients of its image under the projection wrt. the
            // standard P1/Q1 transformation
            Tiny::Vector<CoordType, SurfaceMesh_::shape_dim+1> coeffs(CoordType(0));

            // For every vertex, we need to know if it was already projected
            bool* vert_todo(new bool[num_verts]);
            for(Index i(0); i < num_verts; ++i)
              vert_todo[i] = true;

            // For every facet, we need to know if it was already added to the search stack
            bool* facet_todo(new bool[num_facets]);
            for(Index i(0); i < num_facets; ++i)
              facet_todo[i] = true;

            // For every subfacet, we need to know if it was already traversed by adding the neighbour across said
            // facet to the search stack. If the surface is concave, it might happen that facet A says "go to my
            // neighbour, B", and then B says "go to my neighbour, A". It's then clear that the wanted point lies on
            // the facet between A and B.
            bool* traversed(new bool[_surface_mesh->get_num_entities(SurfaceMesh_::shape_dim-1)]);

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
              acceptable_distance = Math::pow(acceptable_distance, CoordType(1)/CoordType(idx_mp.num_indices));

              // Find all vertices in the current facet
              for(int j(0); j < idx_mp.num_indices; ++j)
              {
                // The index of the local vertex j in the Meshpart
                Index i_mp(idx_mp(current_facet, Index(j)));

                if(vert_todo[i_mp])
                {
                  vert_todo[i_mp] = false;
                  const auto& x(vtx[ts_verts[i_mp]]);

                  // Find the facet in the SurfaceMesh we project to and compute the coefficients
                  facet_sm = find_cell(coeffs, x, guesses[current_facet], acceptable_distance,
                  _surface_mesh, traversed);
                  // Update the guess for this facet because there might be more points to project
                  guesses[current_facet] = facet_sm;

                  // This is the index of the local vertex j in the real mesh
                  Index i_parent(ts_verts[i_mp]);
                  // Clear the old vertex coordinates so we can just add onto that
                  vtx[i_parent].format(CoordType(0));
                  // Evaluate the surface mesh trafo for the computed coefficients
                  CoordType coeff0(CoordType(1));
                  for(int i(0); i < SurfaceMesh_::shape_dim; ++i)
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
            delete[] traversed;
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
            Tiny::Vector<CoordType, SurfaceMesh_::shape_dim+1, sc_>& coeffs,
            const Tiny::Vector<CoordType, world_dim, sx_>& x,
            Index guess, CoordType acceptable_distance, const SurfaceMesh_* mesh, bool* traversed) const
          {
            CONTEXT(name()+"::find_cell()");

            // Clear traversed information
            for(Index facet(0); facet < mesh->get_num_entities(SurfaceMesh_::shape_dim-1); ++facet)
              traversed[facet] = false;

            Index num_cells(mesh->get_num_entities(SurfaceMeshType::shape_dim));
            // We need to recall which cells we already treated
            bool* cell_marker(new bool[num_cells]);
            for(Index cell(0); cell < num_cells; ++cell)
              cell_marker[cell] = false;

            // Type of the vertex at cell index set
            typedef typename SurfaceMesh_::template IndexSet<SurfaceMesh_::shape_dim, 0>::Type VertAtCellIdxType;
            // Number of vertices per cell
            static constexpr int num_vert_loc = VertAtCellIdxType::num_indices;
            // Vertex at cell information
            const VertAtCellIdxType& idx(mesh->template get_index_set<SurfaceMeshType::shape_dim, 0>());
            // Facet at cell information
            const auto& facet_idx(mesh->template
                get_index_set<SurfaceMeshType::shape_dim, SurfaceMesh_::shape_dim-1>());
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

            // Search all cells
            while(!cell_stack.empty())
            {
              // Get new cell
              current = cell_stack.front();
              cell_stack.pop_front();
              cell_marker[current] = true;

              // Collect all vertices that make up the cell in the coords matrix
              for(int j(0); j < SurfaceMesh_::shape_dim+1; ++j)
              {
                Index i(idx(current, Index(j)));
                coords[j] = vtx[i];
              }

              // Compute the coefficients of the inverse coordinate transformation. The last coefficient is the
              // signed distance of x to the facet made up by coords
              compute_inverse_coeffs(coeffs, x, coords);

              // From the coefficients, we can compute the barycentric coordinates (as the surface mesh is of Simplex
              // shape) so we get an idea about in which direction we need to search first
              Tiny::Vector<CoordType, SurfaceMesh_::shape_dim+1> bary(CoordType(1));
              CoordType ortho_dist = compute_bary_and_dist(bary, coeffs);

              // From all this, we can now determine if there is an orthogonal projection of x onto this facet. For
              // this, we need to know the traversed information as well and the neighbour information of the current
              // facet
              success = determine_success(ortho_dist, acceptable_distance, bary, facet_idx[current], traversed);

              // Ok, maybe we are done
              if(success)
              {
                break;
              }
              else
              {
                // If we are not done, add all neighbours that have not been searched to the search stack
                for(int l(0); l < SurfaceMesh_::shape_dim+1; ++l)
                {
                  // The neighbour over the l-th facet
                  Index n(neigh(current, Index(l)));
                  // Check if the neighbour exists and is still on the to do list
                  if(n != ~ Index(0) && !cell_marker[n])
                  {
                    Index facet(facet_idx(current,Index(l)));
                    // If the neighbour is still on the to do list and the facet it lies at is in the traversed list,
                    // this means the barycentric coordinate for the corresponding direction was negative, meaning
                    // that x lies in that direction. So we add that neighbour to the front.
                    if(traversed[facet])
                      cell_stack.push_front(n);
                    else
                      cell_stack.push_back(n);

                    // Since we added the neighbour to the stack, prevent it from being added again later.
                    cell_marker[n] = true;
                  }
                } // adding neighbours
              } // success check
            } // search loop

            // clean up
            delete[] cell_marker;

            // If we get here without success being true, the point could not be found
            if(!success)
            {
              String msg("Could not find point (");
              for(int i(0); i < world_dim; ++i)
                msg += (" "+stringify(x[i]));

              msg+=" ) in current mesh!";

              throw InternalError(msg);
            }

            // Otherwise, it was found in current and coeffs is still set from that evaluation
            return current;
          } // DiscreteMesh::find_cell()


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
            bary(Intern::coeffperm<n_-1>(0)) = CoordType(1);
            for(int j(0); j < n_-1; ++j)
            {
              bary(Intern::coeffperm<n_-1>(0)) -= coeffs(j);
              bary(Intern::coeffperm<n_-1>(j+1)) = coeffs(j);
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
          static bool determine_success(DT_ ortho_dist, DT_ acceptable_distance,
          Tiny::Vector<DT_, n_, sn_>& bary, const IdxType_& facet_idx, bool* traversed)
          {
            DT_ tol = Math::sqrt(Math::eps<DT_>());
            // Determine success
            bool success(true);
            if(ortho_dist > acceptable_distance)
              success = false;

            for(int l(0); l < n_; ++l)
            {
              Index facet(facet_idx[Index(l)]);
              if(bary(l) < - tol)
              {
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
            return success;
          }

          /// \cond internal
          /**
           * \copydoc FEAST::Trafo::Standard::compute_inverse_coeffs()
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
           * \copydoc FEAST::Trafo::Standard::compute_inverse_coeffs()
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

      }; // class DiscreteChart

    } // namespace Atlas

  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_ATLAS_DISCRETE_CHART_HPP
