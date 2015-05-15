#pragma once
#ifndef KERNEL_GEOMETRY_INDEX_CALCULATOR_HPP
#define KERNEL_GEOMETRY_INDEX_CALCULATOR_HPP 1

// includes, FEAST
#include <kernel/geometry/intern/index_representative.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>
#include <kernel/geometry/index_set.hpp>
#include <kernel/util/exception.hpp>

// includes, system
#include <set>
#include <vector>

namespace FEAST
{
  namespace Geometry
  {
    /**
     * \brief Stores the index representatives of an index set
     *
     * \todo detailed description
     *
     * \author Constantin Christof
     */
    template<typename Shape_>
    class IndexTree
    {
    public:
      /// number of indices per index vector
      static constexpr int num_indices = Shape::FaceTraits<Shape_, 0>::count;

      /// Index vector
      class IndexVector
      {
      public:
        Index idx[num_indices];

        IndexVector()
        {
        }

        IndexVector(const IndexVector& iv)
        {
          for(Index i(0); i < Index(num_indices); ++i)
            idx[i] = iv[i];
        }

        Index& operator[](Index i)
        {
          return idx[i];
        }

        const Index& operator[](Index i) const
        {
          return idx[i];
        }

        bool operator<(const IndexVector& other) const
        {
          // Lexicographical comparison ignoring the first entry
          for(Index i(1); i < Index(num_indices); ++i)
          {
            if (idx[i] < other[i])
            {
              return true;
            }
            else if (idx[i] > other[i])
            {
              return false;
            }
          }
          return false;
        }

      }; // class IndexTree::IndexVector

    private:
      /// Set of IndexVector representatives
      typedef std::set<IndexVector> RepSet;
      /// Vector of IV rep sets
      typedef std::vector<RepSet> RepSetVector;

      /// representative set vector
      RepSetVector _rep_set_vec;

    public:
      /**
       * \brief Constructor.
       *
       * \param[in] num_vertices
       * The total number of vertices in the mesh.
       */
      explicit IndexTree(Index num_vertices)
        : _rep_set_vec(num_vertices)
      {
        CONTEXT(name() + "::IndexTree()");
      }

      /// Destructor
      virtual ~IndexTree()
      {
        CONTEXT(name() + "::~IndexTree()");
      }

      /// returns number of indices of an index-representative
      Index get_num_indices() const
      {
        CONTEXT(name() + "::get_num_indices()");
        return Shape::FaceTraits<Shape_, 0>::count;
      }

      /// returns size of the i-th representative set
      Index get_set_size(Index i) const
      {
        CONTEXT(name() + "::get_set_size()");
        ASSERT_(i < _rep_set_vec.size());
        return Index(_rep_set_vec.at(i).size());
      }

      /// returns the value of the k-th component of the j-th index-representative in the i-th set
      Index get_index(Index i, Index j, Index k) const
      {
        CONTEXT(name() + "::get_index()");
        typename RepSet::const_iterator iter = _rep_set_vec[i].begin();
        advance(iter, j);
        return (*iter)[k];
      }

      /**
       * \brief Searches for an index vector within the tree.
       *
       * This function searches the index tree for an index vector's representative and, if found,
       * returns its id within the tree.
       *
       * \param[in] index_vector
       * The index vector whose representative is to be found.
       *
       * \returns
       * A bool-Index-pair indicating the result of the search.
       *   - If the first (bool) component is \c true then the second (Index) component contains the
       *     id of the index vector within the tree.
       *   - If the first (bool) component is \c false then the index vector's representative was
       *     not found within the index tree.
       */
      template<typename IndexVectorType_>
      std::pair<bool,Index> find(const IndexVectorType_& index_vector) const
      {
        CONTEXT(name() + "::find()");

        // calculate representative
        IndexVector representative;
        Intern::IndexRepresentative<Shape_>::compute(representative, index_vector);

        // get the corresponding representative set
        const RepSet& rep_set = _rep_set_vec[representative[0]];

        // try to find the representative
        typename RepSet::const_iterator iter = rep_set.find(representative);
        if(iter == rep_set.end())
          return std::make_pair(false, Index(0));
        else
          return std::make_pair(true, Index((*iter)[0]));
      }

      /**
       * \brief Inserts an index vector's representative into the index tree.
       *
       * \param[in] index_vector
       * The index vector whose representative is to be stored.
       *
       * \param[in] id
       * An id that is to be associated with the index vector.
       */
      template<typename IndexVectorType_>
      void insert(const IndexVectorType_& index_vector, Index id)
      {
        CONTEXT(name() + "::insert<...>()");

        // calculate representative
        IndexVector representative;
        Intern::IndexRepresentative<Shape_>::compute(representative, index_vector);

        // insert representative
        Index first_index = representative[0];
        ASSERT(first_index < _rep_set_vec.size(), "index out-of-range");

        representative[0] = id;
        _rep_set_vec[first_index].insert(representative);
      }

      /**
       * \brief Parses an index set into the tree.
       *
       * \param[in] index_set
       * The index set that is to be parsed.
       */
      template<typename IndexSet_>
      void parse(const IndexSet_& index_set)
      {
        CONTEXT(name() + "::parse()");

        static_assert(int(IndexSet_::num_indices) == int(num_indices), "index count mismatch");

        // fetch number of entities
        const Index num_entities = index_set.get_num_entities();

        // loop over all entities
        for(Index i(0); i < num_entities; ++i)
        {
          // insert the index vector
          insert(index_set[i], i);
        }
      }

      /**
       * \brief Enumerates the index vector representatives.
       *
       * This function loops over all index vector representatives in the index tree and assigns
       * an unique id to each representative. The ids are distributed in consecutive order beginning
       * from zero.
       *
       * \returns
       * The total number of representatives; coincides with the first unused id.
       */
      Index enumerate()
      {
        Index cur_id = 0;

        // loop over all index vector sets
        Index n = Index(_rep_set_vec.size());
        for(Index i(0); i < n; ++i)
        {
          typename RepSet::iterator it(_rep_set_vec[i].begin());
          typename RepSet::iterator jt(_rep_set_vec[i].end());
          for(; it != jt; ++it, ++cur_id)
          {
            // The RepSet is an std::set and its elements cannot be modified as this might screw up the ordering.
            // As we specifically abuse the 0th entry in each set element for saving an index and exclude that
            // from the comparison operator, it is ok to cast away the const qualifier and modify the 0th entry.
            const_cast<IndexVector&>(*it)[0] = cur_id;
          }
        }

        // return total number of index vectors
        return cur_id;
      }

      /// \brief Returns the class name
      static String name()
      {
        return "IndexTree<" + Shape_::name() + ">";
      }
    }; // class IndexTree


    /**
     * \brief Calculates the missing index sets if the vertex-at-shape
     * index sets are given.
     *
     * \author Constantin Christof
     */
    template<
      typename Shape_,
      int face_dim_>
    class IndexCalculator
    {
    public:
      /// Type of the subshape to calculate the missing information at
      typedef typename Shape::FaceTraits<Shape_, face_dim_>::ShapeType FaceType;
      /// Type for the IndexTree containing vertex\@subshape information
      typedef IndexTree<FaceType> IndexTreeType;

    public:
      /**
       * \brief Calculates an index set from vertex\@shape information
       *
       * Assume that for a Hypercube<3>, we have the vertex\@cell information. This routine can generate the
       * edge\@cell or face\@cell information from this.
       *
       * More generic: From vertex\@shape, generate subshape\@shape, where the type/dimension of subshape is
       * determined by face_dim_
       *
       * \tparam IndexSetIn_
       * Type for the vertex\@shape information IndexSet.
       *
       * \tparam IndexSetOut_
       * Type for the subshape\@shape information IndexSet.
       *
       * \param[in] index_tree
       * For every entity of subshape, this IndexTree holds the information which vertices this entity contains.
       *
       * \param[in] index_set_in
       * Provided vertex\@shape information.
       *
       * \param[out] index_set_out
       * subshape\@shape information.
       *
       */
      template<
        typename IndexSetIn_,
        typename IndexSetOut_>
      static bool compute(
        const IndexTreeType& index_tree,
        const IndexSetIn_& index_set_in,
        IndexSetOut_& index_set_out)
      {
        CONTEXT(name() + "::compute()");

        // index vector reference
        typedef typename IndexSetIn_::ConstIndexVectorReference ConstIndexVectorRefIn;
        typedef Intern::FaceIndexMapping<Shape_, face_dim_, 0> FimType;

        // fetch number of shapes
        const Index num_entities = index_set_in.get_num_entities();

        typename IndexTreeType::IndexVector current_face_indices;

        // loop over all shapes
        for(Index i(0); i < num_entities; ++i)
        {
          // get vertex-index-vector of shape i
          ConstIndexVectorRefIn current_cell_in = index_set_in[i];

          // loop over all cells of shape i
          for(int j(0); j < IndexSetOut_::num_indices; ++j)
          {
            // get the vertex-index-vector of cell j:
            // loop over all indices of the cell-vertex-vector
            for(int k(0); k < IndexTreeType::num_indices; ++k)
            {
              current_face_indices[Index(k)]  = current_cell_in[FimType::map(j, k)];
            }

            // try to find the index of the vector within the tree
            std::pair<bool, Index> bi = index_tree.find(current_face_indices);
            if(!bi.first)
              return false;
            index_set_out[i][j] = bi.second;
          }
        }

        // okay, all index vectors found
        return true;
      }

      /**
       * \brief For given vertex\@shape information, numbers subshapes and calculates vertex\@subshape
       */
      template<typename IndexSetIn_, typename IndexSetOut_>
      static void compute_vertex_subshape( const IndexSetIn_& index_set_in, IndexSetOut_& index_set_out)
      {
        CONTEXT(name() + "::compute()");

        // Type for shape to vertex@subshape mapping
        typedef Intern::FaceIndexMapping<Shape_, face_dim_, 0> FimType;
        // Number of shapes
        const Index num_shapes(index_set_in.get_num_entities());
        // Number of verticex
        const Index num_verts(index_set_in.get_index_bound());
        // Number of vertices per subshape
        const Index num_verts_subshape(Shape::FaceTraits<FaceType,0>::count);
        // Number of subshapes per shape, i.e. a Simplex<3> has four Simplex<2> as subshapes
        const Index num_subshapes_shape(Shape::FaceTraits<Shape_,face_dim_>::count);

        // For every vertex, the IndexTree saves which subshapes it is contained in
        IndexTreeType my_index_tree(num_verts);
        // This saves the global vertex numbers in the current subshape
        typename IndexTreeType::IndexVector current_face_indices;

        // Generate the IndexTree of subshapes
        for(Index k(0); k < num_shapes; ++k)
        {
          for(Index l(0); l < num_subshapes_shape; ++l)
          {
            // Get the ith index in the current subshape from the subshape to index mapping for shape k
            for(Index i(0); i < num_verts_subshape; ++i)
              current_face_indices[i] = index_set_in[k][FimType::map(int(l),int(i))];

            // Insert the current subshape into the IndexTree. We need to give an id as 2nd argument, this does not
            // get used in any way
            my_index_tree.insert(current_face_indices, num_verts+1);
          }
        }

        // Enumerate the subshape entities. my_index_tree[i] is the set of all subshapes having i as first vertex.
        // Each set k of those consists of IndexVectors v, where v[0] is the index of the subshape in the subshape
        // numbering and the subsequent entries are the vertex numbers.
        const Index num_subshapes(my_index_tree.enumerate());

        // The output IndexSet has num_subshapes entities and the maximum index is the number of vertices
        IndexSetOut_ my_index_set(num_subshapes, num_verts);
        index_set_out = std::move(my_index_set);

        for(Index i(0); i < num_verts; ++i)
        {
          // Size of the ith set in the IndexTree
          Index n = my_index_tree.get_set_size(i);

          // Iterate over the set
          for(Index j(0); j < n; ++j)
          {
            // This is the index of the subshape
            Index my_index = my_index_tree.get_index(i, j, 0);
            index_set_out[my_index][0] = i;

            // Iterate over the IndexVector that the set element j represents
            // Skip index 0 as this is contains the index of the subshape in the subshape numbering
            for(Index k(1); k < (IndexTreeType::num_indices); ++k)
              index_set_out[my_index][k] = my_index_tree.get_index(i, j, k);
          }

        }
        return;
      }

      /// \brief Returns the class name
      static String name()
      {
        return "IndexCalculator<" + Shape_::name() + "," + stringify(face_dim_) + ">";
      }
    }; // class IndexCalculator

    /// \cond internal
    namespace Intern
    {
      template<typename Shape_, int face_dim_, int cell_dim_ = Shape_::dimension>
      struct RisbHelper
      {
        typedef typename Shape::FaceTraits<Shape_, cell_dim_>::ShapeType CellType;
        typedef typename Shape::FaceTraits<Shape_, face_dim_>::ShapeType FaceType;

        static void compute(IndexSetHolder<Shape_>& ish, IndexTree<FaceType>& idx_tree)
        {
          // recurse down
          RisbHelper<Shape_, face_dim_, cell_dim_ - 1>::compute(ish, idx_tree);

          // create an index calculator
          IndexCalculator<CellType, face_dim_>::compute(
            idx_tree,
            ish.template get_index_set<cell_dim_, 0>(),
            ish.template get_index_set<cell_dim_, face_dim_>());
        }
      };

      template<typename Shape_, int face_dim_>
      struct RisbHelper<Shape_, face_dim_, face_dim_>
      {
        typedef typename Shape::FaceTraits<Shape_, face_dim_>::ShapeType FaceType;

        static void compute(IndexSetHolder<Shape_>&, IndexTree<FaceType>&)
        {
          // dummy
        }
      };

      /**
       * \brief Wrapper struct for generating subshape@shape information
       *
       * \tparam Shape_
       * Highest dimensional shape type
       *
       * \tparam face_dim_
       * Dimension of the subshape
       *
       * Using template recursion, this generates subshape@shape information for all subshapes of lower dimension,
       * i.e. for Simplex<3> it generates edge@face, edge@cell and face@cell
       */
      template<typename Shape_, int face_dim_ = Shape_::dimension - 1>
      struct RisbWrapper
      {
        typedef typename Shape::FaceTraits<Shape_, face_dim_>::ShapeType FaceType;
        static constexpr int num_verts = Shape::FaceTraits<FaceType, 0>::count;

        /**
         * \brief Generates subshape@shape information contained in an IndexSetHolder
         *
         * \param[in,out] ish
         * IndexsetHolder to be filled.
         *
         * The only information this really needs is vertex@shape. If vertex@subshape is missing for any subshape,
         * it is generated using IndexCalculator.
         *
         */
        static void wrap(IndexSetHolder<Shape_>& ish)
        {
          auto& vert_at_subshape_index_set(ish.template get_index_set<face_dim_,0>());

          // Check if vertex@subshape information is present (needed later on) and compute it if necessary
          if(vert_at_subshape_index_set.get_num_entities() == 0)
          {
            IndexCalculator<Shape_, face_dim_>::compute_vertex_subshape(
              ish.template get_index_set<Shape_::dimension, 0>(), vert_at_subshape_index_set);

            // Update dimensions of other IndexSets in the IndexSetHolder
            IndexSetHolderDimensionUpdater<Shape_::dimension, Shape_::dimension-face_dim_>::update(ish);
          }

          // recurse down
          RisbWrapper<Shape_, face_dim_ - 1>::wrap(ish);

          // get vertices-at-face index set
          IndexSet<num_verts>& vert_adj(vert_at_subshape_index_set);

          // build an index-tree from it
          IndexTree<FaceType> idx_tree(vert_adj.get_index_bound());
          idx_tree.parse(vert_adj);

          // Call the helper. This in turn calls IndexCalculater::compute, which needs vertex@shape and
          // vertex@subshape information for all subshapes.
          RisbHelper<Shape_, face_dim_>::compute(ish, idx_tree);
        }
      }; // RisbWrapper<Shape_, face_dim_>

      template<typename Shape_>
      struct RisbWrapper<Shape_, 0>
      {
        static void wrap(IndexSetHolder<Shape_>& /*ish*/)
        {
          // dummy
        }
      }; // RisbWrapper<Shape_, 0>
    } // namespace Intern
    /// \endcond

    /**
     * \brief Builder for redundant index sets
     *
     * This class builds all redundant index sets from the mandatory index sets.
     *
     * \author Peter Zajac
     */
    template<typename Shape_>
    class RedundantIndexSetBuilder
    {
    public:
      /// \brief Routine that does the actual work
      static void compute(IndexSetHolder<Shape_>& index_set_holder)
      {
        Intern::RisbWrapper<Shape_>::wrap(index_set_holder);
      }
    };
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_INDEX_CALCULATOR_HPP
