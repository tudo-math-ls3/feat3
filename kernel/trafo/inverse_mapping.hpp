#pragma once
#ifndef KERNEL_TRAFO_INVERSE_MAPPING_HPP
#define KERNEL_TRAFO_INVERSE_MAPPING_HPP 1

#include <kernel/trafo/mapping_base.hpp>
#include <kernel/util/exception.hpp>

#include <vector>

namespace FEAT
{
  namespace Trafo
  {
    /// \cond internal
    namespace Intern
    {
      template<typename Shape_>
      struct InverseMappingHelper;
    } // namespace Intern
    /// \endcond

    /**
     * \brief Inverse Trafo Mapping Error class
     *
     * This exception is thrown by the InverseMapping class when the
     * point unmapping operation fails due to Newton iteration breakdown.
     *
     * \author Peter Zajac
     */
    class InverseMappingError :
      public FEAT::Exception
    {
    public:
      /**
       * \brief Constructor
       *
       * \param[in] point
       * The point that could not be unmapped.
       *
       * \param[in] cell
       * The index of the cell on which the unmapping failed.
       */
      template<typename DT_, int n_, int s_>
      explicit InverseMappingError(const Tiny::Vector<DT_, n_, s_>& point, Index cell) :
        Exception(String("Failed to unmap point ") + stringify(point) + " on cell " + stringify(cell))
      {
      }
    }; // class InverseMappingError

    /**
     * \brief Data structure for InverseMapping evaluations
     *
     * This class is used for storing the resulting data that arises
     * from an inverse mapping evaluatiopn.
     *
     * \tparam DataType_
     * The datatype that is used for coordinates.
     *
     * \tparam shape_dim_
     * The shape dimension of the underlying mesh.
     *
     * \tparam world_dim_
     * The world dimension of the underlying mesh.
     *
     * \author Peter Zajac
     */
    template<typename DataType_, int shape_dim_, int world_dim_ = shape_dim_>
    class InverseMappingData
    {
    public:
      /// the domain point type
      typedef Tiny::Vector<DataType_, shape_dim_> DomainPointType;
      /// the image point type
      typedef Tiny::Vector<DataType_, world_dim_> ImagePointType;

      /// the image point that was unmapped
      ImagePointType img_point;

      /// the indices of the cells that intersect with the image point
      std::vector<Index> cells;

      /// the domain points on each cell that map onto the image point
      std::vector<DomainPointType> dom_points;

      /// default constructor
      InverseMappingData()
      {
      }

      /// use default copy constructor
      InverseMappingData(const InverseMappingData& other) = default;

      /// move constructor
      InverseMappingData(InverseMappingData&& other) :
        img_point(other.img_point),
        cells(std::forward<std::vector<Index>>(other.cells)),
        dom_points(std::forward<std::vector<DomainPointType>>(other.dom_points))
      {
      }

      /// move-assign operator
      InverseMappingData& operator=(InverseMappingData&& other)
      {
        if(this != &other)
        {
          img_point = other.img_point;
          cells = std::forward<std::vector<Index>>(other.cells);
          dom_points = std::forward<std::vector<DomainPointType>>(other.dom_points);
        }
        return *this;
      }

      /**
       * \brief Checks whether no cells were found.
       */
      bool empty() const
      {
        return cells.empty();
      }

      /// \returns The numbers of cells that were found.
      std::size_t size() const
      {
        return cells.size();
      }
    }; // class InverseMappingData

    /**
     * \brief Inverse Trafo mapping class template
     *
     * This class template implements an "inverse trafo mapping", which
     * is capable of unmapping points, i.e. it can determine the cells
     * and reference coordinates that can be mapped onto a given input
     * point using the "normal" transformation.
     *
     * \tparam Trafo_
     * The type of the trafo mapping that is to be inverted.
     *
     * \tparam DataType_
     * The datatype that is used for internal calculations.
     *
     * \attention
     * <b>Note to implementors:</b> The implementation of this class works
     * only for the standard first-order transformations (both simplexes and
     * hypercube). If you want to add support for higher-order (iso-parametric)
     * transformations, you will need to find a way to compute the cell bounding
     * boxes for the candidate cell selection, as the algorithm implemented in
     * the #init_bounding_boxes() function works only for first-order trafos!
     *
     * \author Peter Zajac
     */
    template<typename Trafo_, typename DataType_>
    class InverseMapping
    {
    public:
      /// the underlying trafo type
      typedef Trafo_ TrafoType;
      /// the datatype to be used
      typedef DataType_ DataType;
      /// the shape type
      typedef typename TrafoType::ShapeType ShapeType;
      /// the shape dimension
      static constexpr int shape_dim = ShapeType::dimension;
      /// the world dimension
      static constexpr int world_dim = TrafoType::world_dim;

      /// the inverse mapping data type
      typedef InverseMappingData<DataType, shape_dim, world_dim> InvMapDataType;

      /// the image point type
      typedef typename InvMapDataType::ImagePointType ImagePointType;
      /// the domain point type
      typedef typename InvMapDataType::DomainPointType DomainPointType;

    protected:
      /// our bounding box type
      typedef Tiny::Matrix<DataType_, 2, world_dim> BBox;
      /// the trafo mapping
      const TrafoType& _trafo;
      /// the array of cell bounding boxes
      std::vector<BBox> _bboxes;
      /// tolerance for unmapped domain coordinates
      DataType _domain_tol;
      /// absolute tolerance for newton iteration
      DataType _newton_tol;
      /// maximum number of newton iterations for unmapping
      Index _newton_max_iter;

    public:
      /**
       * \brief Constructor
       *
       * \param[in] trafo
       * The trafo mapping that is to be inverted.
       *
       * \param[in] bbox_tol
       * The bounding box tolerance used for building the cell bounding boxes.
       *
       * \param[in] domain_tol
       * The coordinate tolerance for the unmapped domain points.
       *
       * This constructor automatically initialises the bounding boxes and
       * sets the default parameters for the Newton iteration to:
       * - tolerance: eps^(0.9)
       * - max. iterations: 10
       */
      explicit InverseMapping(const TrafoType& trafo,
        DataType bbox_tol = DataType(1E-2),
        DataType domain_tol = DataType(1E-4)
        ) :
        _trafo(trafo),
        _domain_tol(domain_tol),
        _newton_tol(Math::pow(Math::eps<DataType>(), DataType(0.9))),
        _newton_max_iter(10u)
      {
        init_bounding_boxes(bbox_tol);
      }

      /// virtual destructor
      virtual ~InverseMapping()
      {
      }

      /**
       * \brief Initialises the cell bounding boxes array
       *
       * \param[in] bbox_tol
       * The bounding box tolerance used for building the cell bounding boxes.
       *
       * \note
       * This function is automatically called by the constructor. You only need to
       * call this function explicitly if the underlying mesh or trafo has changed.
       */
      void init_bounding_boxes(DataType bbox_tol)
      {
        // get the mesh
        const auto& mesh = _trafo.get_mesh();
        // get the vertex set
        const auto& vtx_set = mesh.get_vertex_set();
        // get the vertices-at-cell index set
        const auto& vert_idx = mesh.template get_index_set<shape_dim, 0>();

        // resize bounding boxes array
        _bboxes.resize(std::size_t(mesh.get_num_elements()));

        // loop over all cells
        for(Index cell(0); cell < vert_idx.get_num_entities(); ++cell)
        {
          // get this cell's vertex indices
          const auto& vidx = vert_idx[cell];

          // get this cell's bounding box
          BBox& bbox = _bboxes.at(std::size_t(cell));

          // build cell bounding box
          bbox[0] = bbox[1] = vtx_set[vidx[0]];
          for(int i(1); i < vert_idx.get_num_indices(); ++i)
          {
            // get vertex and update bounding box
            const auto& vtx = vtx_set[vidx[i]];
            for(int j(0); j < world_dim; ++j)
            {
              Math::minimax(DataType(vtx[j]), bbox[0][j], bbox[1][j]);
            }
          }

          // finally, incorporate tolerances
          for(int j(0); j < world_dim; ++j)
          {
            // compute bbox dimension
            DataType bb_extra = bbox_tol * (bbox[1][j] - bbox[0][j]);
            bbox[0][j] -= bb_extra;
            bbox[1][j] += bb_extra;
          }
        }
      }

      /**
       * \brief Unmaps a given image point.
       *
       * This function "unmaps" a point given in real world coordinates, i.e. this
       * function determines the set of cells that intersect with the given point
       * and also computes the corresponding domain points for each of these cells.
       *
       * This function performs several steps:
       * - First, this function determines a set of candidate cells by using a
       *   simple bounding box test provided by the #find_candidate_cells() function.
       * - Then, for each of the candidate cells, this function tries to unmap
       *   the given point onto the cell to obtain the domain point coordinates
       *   by using a Newton iteration provided bythe #unmap_point_by_newton() function.
       * - Finally, this function checks whether the unmapped point is on the
       *   reference cell by using the #test_domain_point() function.
       *
       * \attention
       * This function may throw an InverseMappingError if the given input point
       * \p img_point could not be unmapped for a candidate cell, unless
       * \p ignore_failures was set to \c true. This usually indicates that
       * the mesh has quite ill-formed cell shapes, which cause the Newton
       * iteration to diverge -- especially, if the point to be tested is
       * not on the candidate cell, but outside of it.
       *
       * \param[in] img_point
       * The image point that is to be unmapped.
       *
       * \param[in] ignore_failures
       * Specifies whether to ignore cells on which the Newton iteration broke down.
       * If set to \c false, an InverseMappingError will be thrown if Newton fails
       * to converge, otherwise the corresponding cell will be skipped.
       *
       * \returns
       * An InverseMappingData object that contains the set of cells and
       * domain points, which are mapped onto the given \p img_point.
       */
      InvMapDataType unmap_point(const ImagePointType& img_point, bool ignore_failures = false) const
      {
        // store image point
        InvMapDataType inv_data;
        inv_data.img_point = img_point;

        // find candidate cells
        std::vector<Index> cells;
        if(!this->find_candidate_cells(cells, img_point))
          return inv_data;

        // loop over all candidate cells
        for(Index cell : cells)
        {
          // try to unmap the point by newton
          DomainPointType dom_point;
          if(!this->unmap_point_by_newton(dom_point, img_point, cell))
          {
            // point unmapping failed!
            if(ignore_failures)
              continue;

            // throw an inverse mapping error
            throw InverseMappingError(dom_point, cell);
          }

          // check if the domain point is on the reference element
          if(!this->test_domain_point(dom_point))
          {
            // nope ==> skip point
            continue;
          }

          // okay
          inv_data.cells.push_back(cell);
          inv_data.dom_points.push_back(dom_point);
        }

        // found at least one cell?
        return inv_data;
      }

      /**
       * \brief Determines a set of candidate cells by a bounding-box test.
       *
       * This function determines a set of "candidate" cells which may intersect
       * with a given image point by applying a simple bounding-box test.
       *
       * \note
       * This function may (and often will) give you "false positives", i.e.
       * the candidate cells that are selected by this function may not
       * intersect with the given input point.\n
       * However, this function \b never produces "false negatives", i.e.
       * if a cell actually intersects with the given input point, then
       * this particular cell will definitely be selected as a candidate.
       *
       * \param[out] cells
       * A vector that receives the indices of all found candidate cells.
       *
       * \param[in] img_point
       * The image point for which the candidates are to be found.
       *
       * \returns
       * \c true, if at least one candidate cell was found, otherwise \c false.
       */
      bool find_candidate_cells(std::vector<Index>& cells, const ImagePointType& img_point) const
      {
        // loop over all cells/bounding boxes
        for(std::size_t cell(0); cell < _bboxes.size(); ++cell)
        {
          const BBox& bbox = _bboxes.at(cell);

          // check whether the point is in the bounding box
          bool is_in = true;
          for(int j(0); j < shape_dim; ++j)
          {
            if((img_point[j] < bbox[0][j]) || (img_point[j] > bbox[1][j]))
            {
              is_in = false;
              break;
            }
          }

          // point inside bounding box?
          if(is_in)
            cells.push_back(Index(cell));
        }

        return !cells.empty();
      }

      /**
       * \brief Unmaps a given point on a particular cell by Newton iteration
       *
       * \param[out] dom_point
       * The domain point that is to be computed.
       *
       * \param[in] img_point
       * The image point that is to be unmapped.
       *
       * \param[in] cell
       * The index of the cell on which \p img_point is to be unmapped.
       *
       * \returns
       * \c true, if the Newton iteration successfully unmapped the input point,
       * or \c false, if the Newton iteration did not converge.
       */
      bool unmap_point_by_newton(DomainPointType& dom_point, const ImagePointType& img_point, const Index cell) const
      {
        // define an evaluator
        typedef typename TrafoType::template Evaluator<ShapeType, DataType>::Type TrafoEvaluator;
        TrafoEvaluator trafo_eval(_trafo);

        // for Newton iteration we need image points and inverse jacobians
        static constexpr TrafoTags trafo_tags = TrafoTags::img_point|TrafoTags::jac_inv;

        // define the trafo eval data
        typename TrafoEvaluator::template ConfigTraits<trafo_tags>::EvalDataType trafo_data;

        // initialise output domain point to cell centre
        for(int i(0); i < shape_dim; ++i)
          dom_point[i] = Shape::ReferenceCell<ShapeType>::template centre<DataType>(i);

        // compute squared tolerance
        const DataType tol_sqr = Math::sqr(_newton_tol);

        // no convergence yet
        bool converged = false;

        // prepare trafo evaluator
        trafo_eval.prepare(cell);

        // newton iteration
        for(Index iter(0); iter < _newton_max_iter; ++iter)
        {
          // evaluate trafo in current point
          trafo_eval(trafo_data, dom_point);

          // compute image point distance vector
          ImagePointType def = trafo_data.img_point - img_point;

          // check defect for tolerance
          if(def.norm_euclid_sqr() < tol_sqr)
          {
            converged = true;
            break;
          }

          // update reference point by newton
          dom_point -= trafo_data.jac_inv * def;
        }

        // finish trafo evaluator
        trafo_eval.finish();

        // Newton converged?
        return converged;
      }

      /**
       * \brief Tests whether a given domain point is on the reference element.
       *
       * \param[in] dom_point
       * The domain point that is to be tested.
       *
       * \returns
       * \c true, if \p dom_point is on the reference element, otherwise \c false.
       */
      bool test_domain_point(const DomainPointType& dom_point) const
      {
        return Intern::InverseMappingHelper<ShapeType>::is_on_ref(dom_point, _domain_tol);
      }
    }; // class InverseMapping

    /// \cond internal
    namespace Intern
    {
      template<int dim_>
      struct InverseMappingHelper<Shape::Simplex<dim_>>
      {
        // checks whether p is on the reference cell
        template<typename PT_, typename DT_>
        static bool is_on_ref(const PT_& p, const DT_ tol)
        {
          // reference simplex: all coords are in range [0,1]
          //                    and the sum of all coords is <= 1
          DT_ s = DT_(0);
          for(int i(0); i < dim_; ++i)
          {
            // coord outside range [-tol, 1+tol] ?
            if((p[i] < -tol) || (p[i] > DT_(1)+tol))
              return false;
            s += p[i];
          }
          // coord sum <= 1+tol ?
          return (s <= DT_(1) + tol);
        }
      };

      template<int dim_>
      struct InverseMappingHelper<Shape::Hypercube<dim_>>
      {
        // checks whether p is on the reference cell
        template<typename PT_, typename DT_>
        static bool is_on_ref(const PT_& p, const DT_ tol)
        {
          // reference hypercube: all coordinates are in range [-1,+1]
          for(int i(0); i < dim_; ++i)
          {
            // coord outside range [-1-tol, +1+tol] ?
            if((p[i] < -DT_(1)-tol) || (p[i] > DT_(1)+tol))
              return false;
          }
          return true;
        }
      };
    } // namespace Intern
    /// \endcond
  } // namespace Trafo
} // namespace FEAT

#endif // KERNEL_TRAFO_INVERSE_MAPPING_HPP
