#pragma once
#ifndef KERNEL_GEOMETRY_INDEX_CALCULATOR_HPP
#define KERNEL_GEOMETRY_INDEX_CALCULATOR_HPP 1

// includes, FEAST
#include <kernel/geometry/congruency/index_representative.hpp>
#include <kernel/geometry/congruency/local_index_mapping.hpp>
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
      /// dummy enum
      enum
      {
        /// number of indices per index vector
        num_indices = Shape::FaceTraits<Shape_, 0>::count
      };

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
          for(int i(0); i < num_indices; ++i)
            idx[i] = iv[i];
        }

        Index& operator[](int i)
        {
          return idx[i];
        }

        const Index& operator[](int i) const
        {
          return idx[i];
        }

        bool operator<(const IndexVector& other) const
        {
          // Lexicographical comparison ignoring the first entry
          for(int i(1); i < num_indices; ++i)
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

      // representative set vector
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
      Index get_index(int i, int j, int k) const
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
        Congruency::IndexRepresentative<Shape_>::compute(representative, index_vector);

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
        Congruency::IndexRepresentative<Shape_>::compute(representative, index_vector);

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

        static_assert(IndexSet_::num_indices == num_indices, "index count mismatch");

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
       * This function loops over all index vector representatives in the index tree and assignes
       * an unique id to each representative. The id's are distributed in consecutive order beginning
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
            (*it)[0] = cur_id;
          }
        }

        // return total number of index vectors
        return cur_id;
      }

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
      int face_dim>
    class IndexCalculator
    {
    public:
      // cell-type (e.g. edge), Shape_ = shape-type (e.g quad)
      typedef typename Shape::FaceTraits<Shape_, face_dim>::ShapeType CellType;
      typedef IndexTree<CellType> IndexTreeType;

    public:
      /**
       * \brief Calculates an index set.
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
        typedef Congruency::LocalIndexMapping<Shape_, face_dim, 0> Lim;

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
              current_face_indices[k]  = current_cell_in[Lim::map(j, k)];
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

      static String name()
      {
        return "IndexCalculator<" + Shape_::name() + "," + stringify(face_dim_) + ">";
      }
    }; // class IndexCalculator
  } // namespace Geometry
} // namespace FEAST

#endif // KERNEL_GEOMETRY_INDEX_CALCULATOR_HPP
