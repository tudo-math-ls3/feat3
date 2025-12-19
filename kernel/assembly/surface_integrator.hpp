// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once

#include <kernel/base_header.hpp>
#include <kernel/cubature/dynamic_factory.hpp>
#include <kernel/cubature/rule.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/trafo/standard/inverse_mapping.hpp>
#include <kernel/util/tiny_algebra.hpp>

#include <vector>

namespace FEAT::Assembly
{

  template<typename DT_, int dim_>
  using BoundingBoxType = Tiny::Matrix<DT_, 2, dim_>;

  namespace Intern
  {
    template<typename DT_, int dim_, typename VertexSet_, typename IndexTuple_>
    static inline BoundingBoxType<DT_, dim_> get_boundingbox(const VertexSet_& vtx_set, const IndexTuple_& vidx, DT_ tol)
    {
      typedef DT_ DataType;
      BoundingBoxType<DT_, dim_> bbox;

      // build cell bounding box
      bbox[0] = bbox[1] = vtx_set[vidx[0]];
      for(int i(1); i < vidx.num_indices; ++i)
      {
        // get vertex and update bounding box
        const auto& vtx = vtx_set[vidx[i]];
        for(int j(0); j < dim_; ++j)
        {
          Math::minimax(DataType(vtx[j]), bbox[0][j], bbox[1][j]);
        }
      }

      // finally, incorporate tolerances
      for(int j(0); j < dim_; ++j)
      {
        // compute bbox dimension
        DataType bb_extra = tol * (bbox[1][j] - bbox[0][j]);
        bbox[0][j] -= bb_extra;
        bbox[1][j] += bb_extra;
      }

      return bbox;
    }

    template<typename DT_, int dim_>
    static inline bool check_point(const Tiny::Vector<DT_, dim_>& pt, const BoundingBoxType<DT_, dim_>& bb)
    {
      bool outside = false;
      for(int i = 0; i < dim_; ++i)
      {
        outside |= (pt[i] < bb[0][i]) || (pt[i] > bb[1][i]);
      }
      return !outside;
    }
  }

  /**
   * \brief Assemble class Discrete Integrale of an FE function over a d-1 manifold
   *
   * \tparam Trafo_ The trafo on which the FE spaces are defined on
   * \tparam SurfaceShape_ The shapetype of the surface elements
   *
   * \note For now, this only works for piecewise defined triangulations and we assume no cuvature
   *       on the single segments. Maybe we will add this in the future.
   *
   * Also any vertex mapping provided to this class is assumed to have the 0-th submanifold ordering in the sense of the FEAT3 mesh element mapping, i.e.
   * We always look onto the inner face and the ordering is always left to right:
   *
   *  Line2d:         0 ----- 1
   *
   *  Quad:           2 ----- 3
   *                  |       |
   *                  |       |
   *                  0 ----- 1
   *
   * Triang:          2
   *                  | -
   *                  |  |
   *                  |   -
   *                  0 --- 1
   *
   * \note It is assumed, that the normal orientation is the same for all elements.
   * We generally assume we view from the outside, so you have to set the _view_inside variable if this is not case.
   */
  template<typename Trafo_, typename SurfaceShape_>
  class SurfaceIntegrator
  {
  public:
    typedef Trafo_ TrafoType;
    typedef typename TrafoType::MeshType MeshType;
    typedef SurfaceShape_ SurfaceShapeType;
    static constexpr int domain_dim = Trafo_::ShapeType::dimension;
    static constexpr int sub_dim = SurfaceShapeType::dimension;
    typedef typename TrafoType::CoordType DataType;
    typedef Index IndexType;
    typedef Tiny::Vector<DataType, domain_dim> CoordVecType;
    typedef Tiny::Vector<DataType, domain_dim> NormalVecType;
    typedef Tiny::Vector<DataType, sub_dim> RefCoordVecType;
    typedef Tiny::Matrix<DataType, domain_dim, sub_dim> JacobianType;
    typedef Tiny::Matrix<DataType, sub_dim, domain_dim> JacobianInvType;
    typedef Tiny::Tensor3<DataType, domain_dim, sub_dim, sub_dim> HessianType;
    typedef Tiny::Matrix<DataType, domain_dim, sub_dim, sub_dim> HessianInvType;
    typedef Cubature::Rule<SurfaceShapeType, DataType, DataType, RefCoordVecType> CubatureRuleType;

    static constexpr int verts_per_face = Shape::FaceTraits<SurfaceShapeType, 0>::count;


  private:
    const TrafoType& _trafo;
    std::vector<Tiny::Vector<DataType, domain_dim>> _shape_vertices;
    Geometry::IndexSet<verts_per_face> _shape_indices;

    std::vector<Geometry::IndexTuple<verts_per_face>> _shape_index_tmp;

    std::vector<IndexType> _cell_helper;
    std::vector<IndexType> _cell_helper_offsets;
    std::vector<IndexType> _active_surfaces;
    std::vector<int> _cell_mask_vec;

    CubatureRuleType _cubature_rule;

    /// Normal defined from the inside or outside?
    bool _inside_view;

    /// tests whether two boundingboxes intersect in any way
    bool boundingbox_intersect_test(const BoundingBoxType<DataType, domain_dim>& bb_a, const BoundingBoxType<DataType, domain_dim>& bb_b)
    {
      bool no_intersect = false;
      for(int d = 0; d < domain_dim; ++d)
      {
        no_intersect |= (bb_a[0][d] > bb_b[1][d]) || (bb_b[0][d] > bb_a[1][d]);
      }

      return !no_intersect;
    }

  public:

    void set_normal_view(bool inside_view)
    {
      _inside_view = inside_view;
    }

    bool get_normal_view() const
    {
      return _inside_view;
    }

    void compile()
    {
      // move our helper vector
      _shape_indices.set_indices(std::move(_shape_index_tmp));
      // get mesh
      const MeshType& mesh = _trafo.get_mesh();
      const auto& vtx_set = mesh.get_vertex_set();
      // first of all, get bounding box of this mesh
      BoundingBoxType<DataType, domain_dim> mesh_bb;
      for(int k = 0; k < domain_dim; ++k)
      {
        mesh_bb[0][k] = Math::Limits<DataType>::max();
        mesh_bb[1][k] = Math::Limits<DataType>::min();
      }
      {
        for(const auto* verts = vtx_set.begin(); verts != vtx_set.end(); ++verts)
        {
          for(int k = 0; k < domain_dim; ++k)
            Math::minimax((*verts)[k], mesh_bb[0][k], mesh_bb[1][k]);
        }
      }

      // preprocess which surface elements should be checked at all
      _active_surfaces.clear();
      for(IndexType k = 0; k < _shape_indices.get_num_entities(); ++k)
      {
        const auto& loc_ind_tuple = _shape_indices[k];
        // domain_dim boundingbox of our shape
        auto shape_bb = Intern::get_boundingbox<DataType, domain_dim>(_shape_vertices, loc_ind_tuple, DataType(1E-2));
        // are we inside of our local domain boundingbox ?
        if(boundingbox_intersect_test(shape_bb, mesh_bb))
          _active_surfaces.push_back(k);
      }

      // now process which mesh cells are candidates for each surface cell
      _cell_helper.clear();
      _cell_helper_offsets.resize(_shape_indices.get_num_entities()+1);
      for(auto& c : _cell_helper_offsets)
      {
        c = IndexType(0);
      }


      // get mapping element to vertices
      const auto& ele_to_vert = mesh.template get_index_set<domain_dim, 0>();
      const Index num_elements = mesh.get_num_elements();
      constexpr Index batch_size = 1 << 19;
      std::vector<std::vector<IndexType>> gather_helper(std::min(batch_size, Index(_active_surfaces.size())));
      FEAT_PRAGMA_OMP(parallel for)
      for(std::size_t i = 0; i < gather_helper.size(); ++i)
      {
        // reserve 128 elements per entry, this should be enough without triggering a reallocation besides extrem cases
        gather_helper[i].reserve(128);
      }
      // now iterate over our surface elements in a batched manner
      for(Index num_batch = 0; num_batch*batch_size < _active_surfaces.size(); ++num_batch)
      {
        FEAT_PRAGMA_OMP(parallel for)
        for(Index t = 0; t < gather_helper.size(); ++t)
        {
          gather_helper.at(t).clear();
        }

        const Index offset = num_batch*batch_size;
        // now iterate through our mesh and test each cell for its boundingbox... TODO: we could speed this up by using macro elements
        FEAT_PRAGMA_OMP(parallel for)
        for(Index t = 0; t < num_elements; ++t)
        {
          // check if our current element is outside of our fluid domain, an thus should be ignored
          if(!_cell_mask_vec.empty() && _cell_mask_vec.at(t) == 3)
          {
            continue;
          }

          // create boundingbox
          auto cell_bb = Intern::get_boundingbox<DataType, domain_dim>(vtx_set, ele_to_vert[t], DataType(1E-4));

          // now test for each shape if the current cell interesects
          for(Index k = 0; k+offset < std::min(_active_surfaces.size(), offset+batch_size); ++k)
          {
            const auto& loc_ind_tuple = _shape_indices[_active_surfaces.at(k+offset)];
            // domain_dim boundingbox of our shape
            auto shape_bb = Intern::get_boundingbox<DataType, domain_dim>(_shape_vertices, loc_ind_tuple, DataType(1E-4));
            // are we inside of our local domain boundingbox ?
            if(boundingbox_intersect_test(shape_bb, cell_bb))
            {
              FEAT_PRAGMA_OMP(critical)
                gather_helper.at(k).push_back(t);
            }
          }
        }
        // now gather our information into our offset array
        for(std::size_t k = 0; k < gather_helper.size(); ++k)
        {
          const auto& v_info = gather_helper[k];
          _cell_helper_offsets.at(_active_surfaces.at(k+offset)+1) = v_info.size();
          for(const auto& val : v_info)
            _cell_helper.push_back(val);
        }
        // now gather the values
      }
      // now perform partial sum to get actual offsets
      std::partial_sum(_cell_helper_offsets.begin(), _cell_helper_offsets.end(), _cell_helper_offsets.begin());

      XASSERTM(_cell_helper_offsets.back() == _cell_helper.size(), "Offsets do not match size of cell helper array");

    }


  public:
    /**
     * \brief Initializes the integrator with a cell wise list of vertices
     *
     * \param[in] vertices Vector of tupled vertex data
     * \param[in] index The sparse index array into the vertex data, e.g. face 7 -> Verts (0, 7, 14, 36)
     */
    explicit SurfaceIntegrator(const Trafo_& trafo, const Cubature::DynamicFactory& cubature_factory_) :
     _trafo(trafo),
     _cubature_rule(Cubature::ctor_factory, cubature_factory_),
     _inside_view(false)
    {
    }

    /**
     * Careful, this overwrites any previous vertex info
     */
    template<typename Vtxs_>
    void set_vertices(const std::vector<Vtxs_>& vertices)
    {
      _shape_index_tmp.clear();
      _shape_vertices.resize(vertices.size());
      FEAT_PRAGMA_OMP(parallel for)
      for(std::size_t k = 0; k < std::size_t(vertices.size()); ++k)
      {
        _shape_vertices[k] = vertices[k];
      }
    }

    /**
     * \brief Simply pushes new vertices to the end of the current vertex array
     */
    template<typename Vtxs_>
    void add_face_vertices(const Vtxs_& face_vertices)
    {
      for(std::size_t i = 0; i < std::size_t(verts_per_face); ++i)
      {
        _shape_vertices.push_back(face_vertices[i]);
      }

      const IndexType num_indx = _shape_index_tmp.size();
      _shape_index_tmp.emplace_back();
      for(IndexType i = 0; i < IndexType(verts_per_face); ++i)
      {
        _shape_index_tmp.back()[int(i)] = num_indx*verts_per_face + IndexType(i);
      }
    }

    /**
     * Adds a new face, indexing into the current vertex set
     */
    template<typename VertIdx_>
    void add_face(const VertIdx_& vert_idx)
    {
      _shape_index_tmp.emplace_back();
      for(int i = 0; i < verts_per_face; ++i)
      {
        _shape_index_tmp.back()[i] = IndexType(vert_idx[std::size_t(i)]);
      }
    }


    template<typename Job_>
    void assemble(Job_& job)
    {
      if(this->_shape_indices.get_num_entities() == 0u)
      {
        return;
      }

      assemble_omp(job);
    }

    /**
     * \brief Sets the cell mask vector, which sorts out any cells completelly outside the fluid domain
     *
     * \param[in] cell_mask_vec Mask vector, for example generated by fbm_asm->get_fbm_mask_vector(dim)
     *
     * \note the mask vec has to match the size of all entities of the background mesh and holds the following encoding
     *  3 -> cell is complettly inside the fbm region, thus the cell is ignored
     *  any other -> cell has some subdimensional entitiy, which is outside the fbm region, i.e. part of the fluid domain
     *
     */
    template<typename MaskVec_>
    void set_mask_vector(const MaskVec_& cell_mask_vec)
    {
      _cell_mask_vec.resize(cell_mask_vec.size());
      FEAT_PRAGMA_OMP(parallel for)
      for(Index k = 0; k < cell_mask_vec.size(); ++k)
      {
        _cell_mask_vec[k] = int(cell_mask_vec[k]);
      }
    }

    template<typename Job_>
    void assemble_omp(Job_& job)
    {
      typedef typename Job_::Task TaskType;
      // typedef for our evalhelper
      typedef Trafo::Standard::EvalHelper<DataType, RefCoordVecType, CoordVecType, JacobianType, JacobianInvType, DataType, HessianType, HessianInvType, SurfaceShapeType, domain_dim> EvalHelper;

      FEAT_PRAGMA_OMP(parallel)
      {
        // create a task for this thread
        std::unique_ptr<TaskType> task(new TaskType(job));

        FEAT_PRAGMA_OMP(for)
        for(IndexType k = 0; k < _active_surfaces.size(); ++k)
        {
          // construct actual surface
          Tiny::Vector<CoordVecType, verts_per_face> surface_verts;
          for(int i = 0; i < verts_per_face; ++i)
          {
            surface_verts[i] = _shape_vertices[_shape_indices[_active_surfaces[k]][i]];
          }
          // we now have to construct the transformation from the surface to the reference element, here we use our trafo capabilities
          // we assume we have a id 0 orientation
          Tiny::Matrix<DataType, domain_dim, verts_per_face> coefficients;
          coefficients.format();
          EvalHelper::set_coefficients(coefficients, _shape_vertices, _shape_indices, _active_surfaces[k]);

          // our jacobian of the transformation
          JacobianType jac;
          jac.format();
          // since we have standard trafo, the jacobian is independent of our actual domain point
          // our determinant
          DataType jac_det = DataType(0);

          DataType face_volume = EvalHelper::volume(coefficients);

          // our cubature points
          std::vector<CoordVecType> points(std::size_t(_cubature_rule.get_num_points()));
          // our weights
          std::vector<DataType> weights(std::size_t(_cubature_rule.get_num_points()));
          // averaged normal (we assume we have planar normals, in theory this should also be a vector)
          NormalVecType normal(DataType(0));

          // now iterate through our cubature points
          for(int pt = 0; pt < _cubature_rule.get_num_points(); ++pt)
          {
            RefCoordVecType domain_point = _cubature_rule.get_point(pt);
            DataType weight = _cubature_rule.get_weight(pt);
            CoordVecType img_point;
            EvalHelper::map_point(img_point, domain_point, coefficients);
            EvalHelper::calc_jac_mat(jac, domain_point, coefficients);
            jac_det = jac.vol();
            // calculate normal <- normal should point outside
            // generally normal only makes sense for a dim-1 manifold
            if constexpr(sub_dim+1==domain_dim)
            {
              normal = Tiny::orthogonal(jac).normalize();
              if(_inside_view)
                normal = normal.negate();
            }
            else
            {
              normal.format();
            }
            DataType point_weight = jac_det * weight;
            points[std::size_t(pt)] = img_point;
            weights[std::size_t(pt)] = point_weight;
          }

          // prepare our job, we simply provide all the candidate cells, its the jobs job to handle the information
          task->prepare(_cell_helper, _cell_helper_offsets, points, weights, normal, face_volume, _active_surfaces[k]);

          // now, assemble our task
          task->assemble();

          // scatter
          FEAT_PRAGMA_OMP(critical)
          task->scatter();

        }
        // and combine
        FEAT_PRAGMA_OMP(critical)
        task->combine();
      }

    }

  }; // class SurfaceIntegrator
}
