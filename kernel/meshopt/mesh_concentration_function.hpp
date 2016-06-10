#pragma once
#ifndef KERNEL_MESHOPT_MESH_CONCENTRATION_FUNCTION_HPP
#define KERNEL_MESHOPT_MESH_CONCENTRATION_FUNCTION_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/meshopt/rumpf_trafo.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/comm_base.hpp>
#include <kernel/util/property_map.hpp>

#include <deque>

namespace FEAT
{
  namespace Meshopt
  {
    /**
     * \tparam Trafo_
     * Our transformation.
     *
     * \tparam RefCellTrafo_
     * Mesh optimisation reference cell transformation.
     *
     * \author Jordi Paul
     *
     **/
    template
    <
      typename Trafo_,
      typename RefCellTrafo_ = RumpfTrafo<Trafo_, typename Trafo_::CoordType>
    >
    class MeshConcentrationFunction
    {
      public:
        /// Type for the transformation
        typedef Trafo_ TrafoType;
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;

        /// Only Mem::Main is supported atm
        typedef Mem::Main MemType;
        /// The index type
        typedef Index IndexType;

        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;
        /// Vector type for element sizes etc.
        typedef LAFEM::DenseVector<MemType, CoordType, IndexType> ScalarVectorType;
        /// Vector type for element scales etc.
        typedef LAFEM::DenseVectorBlocked<MemType, CoordType, IndexType, MeshType::world_dim> VectorType;

        const Geometry::RootMeshNode<MeshType>* _mesh_node;
        ScalarVectorType _dist;
        VectorType _grad_dist;
        /// The sum of all entries in _conc
        CoordType _sum_conc;
        /// Vector for saving the mesh concentration, one entry per cell
        ScalarVectorType _conc;
        /// Gradient of the mesh concentration wrt. the world coordinates
        VectorType _grad_conc;
        VectorType _grad_sum_det;
        /// Gradient of the local optimal scales h wrt. the vertex coordinates.
        // Each entry in the DenseVectorBlocked represents one cell. Each cell's block contains
        // world_dim*(number of local vertices) entries. Has to be serialised like this because there is no
        // DenseVector that saves a Tiny::Matrix
        LAFEM::DenseVectorBlocked<MemType, CoordType, Index, MeshType::world_dim*Shape::FaceTraits<ShapeType,0>::count> _grad_h;
        typedef LAFEM::DenseVectorBlocked<MemType, CoordType, Index, MeshType::world_dim*Shape::FaceTraits<ShapeType,0>::count> GradHType;

        explicit MeshConcentrationFunction() :
          _mesh_node(nullptr),
          _dist(),
          _grad_dist(),
          _sum_conc(CoordType(0)),
          _conc(),
          _grad_conc(),
          _grad_sum_det(),
          _grad_h()
          {
          }
      protected:
        explicit MeshConcentrationFunction(const MeshConcentrationFunction& other) :
          _mesh_node(other._mesh_node),
          _sum_conc(other._sum_conc)
        {
          _dist.clone(other._dist, LAFEM::CloneMode::Deep);
          _grad_dist.clone(other._grad_dist, LAFEM::CloneMode::Deep);
          _grad_conc.clone(other._grad_conc, LAFEM::CloneMode::Deep);
          _grad_sum_det.clone(other._grad_sum_det, LAFEM::CloneMode::Deep);
          _grad_h.clone(other._grad_h, LAFEM::CloneMode::Deep);
        }
      public:

        /// \brief Virtual destructor
        virtual ~MeshConcentrationFunction()
        {
        }

        virtual void compute_dist() = 0;

        const ScalarVectorType& get_conc() const
        {
          return _conc;
        }

        void set_mesh_node(const Geometry::RootMeshNode<MeshType>* mesh_node_)
        {
          XASSERT(_mesh_node == nullptr);

          _mesh_node = mesh_node_;

          _sum_conc = CoordType(0);

          Index ndofs(_mesh_node->get_mesh()->get_num_entities(0));
          Index ncells(_mesh_node->get_mesh()->get_num_entities(ShapeType::dimension));

          _dist = ScalarVectorType(ndofs, CoordType(0));
          _grad_dist = VectorType(ndofs,CoordType(0));
          _conc = ScalarVectorType(ncells);
          _grad_conc = VectorType(ndofs, CoordType(0));
          _grad_h = GradHType(ncells,CoordType(0));
          _grad_sum_det = VectorType(ndofs, CoordType(0));
        }

        /**
         * \brief The class name
         *
         * \returns String with the class name
         */
        virtual String name() const
        {
          return "MeshConcentrationFunction<"+MeshType::name()+">";
        }

        /**
         * \brief Computes the mesh concentration function for each cell
         *
         * This means it computes
         * \f[
         *   \forall K \in \mathcal{T}_h: \underbrace{c(K)}_{ =: \mathtt{conc(cell)}} ~ \mathrm{and} ~ \underbrace{\sum_{K \in \mathcal{T}_h} c(K)}_{\mathtt{\_sum\_conc}}
         * \f]
         **/
        void compute_conc()
        {
          // Total number of cells in the mesh
          const Index ncells(_mesh_node->get_mesh()->get_num_entities(ShapeType::dimension));
          // Index set for local/global numbering
          const auto& idx = _mesh_node->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          _sum_conc = CoordType(0);
          for(Index cell(0); cell < ncells; ++cell)
          {
            // This will be the average of the levelset values at the vertices
            CoordType tmp(0);

            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              tmp += this->_dist(idx(cell,Index(j)));

            tmp = tmp/CoordType(Shape::FaceTraits<ShapeType,0>::count);

            this->_conc(cell, conc_val(tmp));
            _sum_conc+=this->_conc(cell);
          }

#ifdef FEAT_HAVE_MPI
          CoordType sum_conc_snd(_sum_conc);
          Util::Comm::allreduce(&sum_conc_snd, 1, &_sum_conc, MPI_SUM);
#endif
        }

        /**
         * \brief Computes the value of the concentration function
         *
         * \param[in] dist
         * The distance.
         *
         * \returns
         * The mesh concentration for the given distance.
         *
         */
        CoordType conc_val(CoordType dist) const
        {
          return dist; // Math::pow(this->_r_adapt_reg + Math::abs(dist),this->_r_adapt_pow);
        }

        /**
         * \brief Computes the derivative of the concentration function
         *
         * \param[in] dist
         * The distance.
         *
         * \returns
         * The derivative of mesh concentration wrt. the distance.
         *
         */
        CoordType conc_der(CoordType DOXY(dist)) const
        {
          return CoordType(0); //this->_r_adapt_pow * Math::pow(this->_r_adapt_reg + Math::abs(dist),this->_r_adapt_pow - CoordType(1))*Math::signum(dist);
        }

        /**
         * \brief Computes the local gradient of the concentration function wrt. the vertices
         *
         * \tparam Tgrad_
         * Type of the local gradient, i.e. Tiny::Matrix
         *
         * \tparam Tl_
         * Type for the vector of levelset values, i.e. Tiny::vector
         *
         * \tparam Tgradl_
         * Type for the local gradient of the levelset values wrt. the vertices, i.e. Tiny::Matrix
         *
         * \param[out] grad_loc_
         * The gradient of the concentration wrt. the vertices
         *
         * \param[in] dist_vals_
         * The levelset values at the vertices
         *
         * \param[in] dist_grad_vals_
         * The grandient of the levelset wrt. the vertices, evaluated at the vertices themselves
         *
         **/
        template<typename Tgrad_, typename Tl_, typename Tgradl_>
        void compute_grad_conc_local(Tgrad_& grad_loc_, const Tl_& dist_vals_, const Tgradl_& dist_grad_vals_) const
        {
          grad_loc_.format(CoordType(0));

          // This will be the average of the levelset values at the vertices
          CoordType val(0);
          for(int i(0); i < Shape::FaceTraits<ShapeType,0>::count; ++i)
            val += dist_vals_(i);

          val = val/CoordType(Shape::FaceTraits<ShapeType,0>::count);

          for(int i(0); i < Shape::FaceTraits<ShapeType,0>::count; ++i)
            grad_loc_[i] = conc_der(val)/CoordType(Shape::FaceTraits<ShapeType,0>::count) * dist_grad_vals_[i];
        }

        /**
         * \brief Computes the local gradient of the optimal scales
         *
         * The optimal scales \f$ h \f$ depend on the concentration \f$ c \f$ by the relation
         * \f[
         *   \forall K_k \in \mathcal{T}_h: h_k = \left( \frac{c (K_k)}{\sum_{l=1}^N c(K_l) \sum_{l=1}^N
         *   \mathrm{det} \nabla R_{T,l}(\phi)} \right)^{\frac{1}{d}}
         * \f]
         *
         * The concentration function in turn indirectly depends on the vertex locations
         * \f$ \{ x_i : i=1, \dots, n \} \f$ via
         * \f[ c(K_k) = (\alpha + \sum_{x_j \in K_k} \varphi(x_j))^\beta, \f]
         * so that we have to take this dependency into account for the full gradient. Define
         * \f[
         *   s_d := \sum_{l=1}^N \mathrm{det} \nabla R_{T,l}(\Phi), s_c :=  \frac{c (K_k)}{\sum_{l=1}^N c(K_l)}.
         * \f]
         * So for each \f$ x_j \f$ we arrive at
         * \f[
         *   \frac{\partial h(K_k)}{\partial x_j} = \frac{1}{d} \left( \frac{c(K_k)}{s_d s_c} \right)^
         *   {\frac{1}{d}-1} \left[ c(K_k) ( \frac{\partial s_d}{\partial x_j} s_c + s_d
         *   \frac{\partial s_c}{\partial x_j} ) + \frac{\partial c(K_k)}{\partial x_j} s_d s_c \right]
         *   (s_d s_c)^{-2}
         * \f]
         *
         */
        void compute_grad_h(const VectorType& coords)
        {
          XASSERT(_mesh_node != nullptr);

          _grad_h.format();

          CoordType sum_det = RefCellTrafo_::compute_sum_det(coords, *(_mesh_node->get_mesh()));
          RefCellTrafo_::compute_grad_sum_det(_grad_sum_det, coords, *(_mesh_node->get_mesh()));

          compute_grad_conc();

          // Index set for local/global numbering
          auto& idx = _mesh_node->get_mesh()->template get_index_set<ShapeType::dimension,0>();
          // This will hold the levelset values at the mesh vertices for one element
          Tiny::Vector<CoordType, Shape::FaceTraits<ShapeType,0>::count> dist_vals;
          // This will hold the levelset gradient values for one element for passing to other routines
          Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> dist_grad_vals;
          // This will hold the local gradient values for one element for passing to other routines
          Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc(0);
          // This will hold the computed local gradient values for one element for copy assigning to the blocked
          // datatype
          Tiny::Vector<CoordType, MeshType::world_dim*Shape::FaceTraits<ShapeType,0>::count> tmp(0);

          CoordType exponent = CoordType(1)/CoordType(MeshType::world_dim) - CoordType(1);

          for(Index cell(0); cell < _mesh_node->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            grad_loc.format();
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell, Index(j)));
              // Get levelset
              dist_vals(j) = this->_dist(i);
              // Get levelset gradient
              dist_grad_vals[j] = this->_grad_dist(i);
            }

            compute_grad_conc_local(grad_loc, dist_vals, dist_grad_vals);

            for(int d(0); d < MeshType::world_dim; ++d)
            {
              for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              {
                Index i(idx(cell, Index(j)));

                tmp(j*MeshType::world_dim +d) =
                  CoordType(1)/CoordType(MeshType::world_dim)*Math::pow(_conc(cell)/_sum_conc*sum_det,exponent)
                  *( _conc(cell)*(_grad_sum_det(i)(d)*_sum_conc + sum_det*_grad_conc(i)(d) )
                      + grad_loc(j,d) * sum_det *_sum_conc)
                  / Math::sqr(_sum_conc*sum_det);
              }
            }
            _grad_h(cell, tmp + _grad_h(cell));

          }

        } // compute_grad_h

        /**
         * \brief Computes the gradient of the sum of all mesh concentrations
         *
         * \f[
         *   \mathrm{grad\_conc}(k,i) = \frac{\partial}{\partial x_j} \sum_{l=1}^N c(K_l),
         * \f]
         * where \f$ i \f$ is the global index of the local vertex \f$ j \f$.
         *
         **/
        void compute_grad_conc()
        {
          // Clear the old gradient
          _grad_conc.format();

          // Index set for local/global numbering
          auto& idx = _mesh_node->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          Tiny::Vector <CoordType,MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc;
          // This will hold the levelset values at the mesh vertices for one element
          Tiny::Vector <CoordType, Shape::FaceTraits<ShapeType,0>::count> dist_loc;
          // This will hold the levelset gradient values for one element for passing to other routines
          Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_dist_loc;

          // Compute the functional value for each cell
          for(Index cell(0); cell < _mesh_node->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            // Collect levelset and levelset grad values
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              // Global vertex/dof index
              Index i(idx(cell, Index(j)));
              // Get levelset
              dist_loc(j) = this->_dist(i);
              // Get levelset gradient
              grad_dist_loc[j] = this->_grad_dist(i);
            }

            // Compute gradient of the concentration on this cell
            compute_grad_conc_local(grad_loc, dist_loc, grad_dist_loc);

            // Add local contributions to global gradient vector
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              Index i(idx(cell, Index(j)));
              this->_grad_conc(i, _grad_conc(i) + grad_loc[j]);
            }
          }

        } // compute_grad_conc()

    }; // class MeshConcentrationFunction

    /**
     * \tparam Trafo_
     * Our transformation.
     *
     * \tparam RefCellTrafo_
     * Mesh optimisation reference cell transformation.
     *
     * \author Jordi Paul
     *
     **/
    template
    <
      typename Trafo_,
      typename RefCellTrafo_ = RumpfTrafo<Trafo_, typename Trafo_::CoordType>
    >
    class ChartDistanceFunction : public MeshConcentrationFunction<Trafo_, RefCellTrafo_>
    {
      public:
        /// Type for the transformation
        typedef Trafo_ TrafoType;
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;

        /// Our base class
        typedef MeshConcentrationFunction<Trafo_, RefCellTrafo_> BaseClass;
        /// Only Mem::Main is supported atm
        typedef Mem::Main MemType;
        /// The index type
        typedef Index IndexType;

        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;
        /// Vector type for element sizes etc.
        typedef LAFEM::DenseVector<MemType, CoordType, IndexType> ScalarVectorType;
        /// Vector type for element scales etc.
        typedef LAFEM::DenseVectorBlocked<MemType, CoordType, IndexType, MeshType::world_dim> VectorType;

      protected:
        std::deque<String> _chart_list;

      public:

        ChartDistanceFunction() = delete;

        explicit ChartDistanceFunction(const std::deque<String>& chart_list_):
          _chart_list()
        {
          XASSERTM(chart_list_.size() > size_t(0), "Empty chart list.");

          for(const auto& it:chart_list_)
            _chart_list.push_back(it);
        }

        explicit ChartDistanceFunction(const ChartDistanceFunction& other) :
          BaseClass(other.BaseClass),
          _chart_list(other._chart_list)
          {
          }

        virtual ~ChartDistanceFunction()
        {
        }

        virtual void compute_dist() override
        {
          XASSERT(this->_mesh_node != nullptr);

          const auto& vtx = this->_mesh_node->get_mesh()->get_vertex_set();

          for(Index i(0); i < this->_mesh_node->get_mesh()->get_num_entities(0); ++i)
          {
            CoordType my_dist(0);
            for(const auto& it:_chart_list)
            {
              auto* chart = this->_mesh_node->get_atlas()->find_mesh_chart(it);
              if(chart == nullptr)
                throw InternalError(__func__,__FILE__,__LINE__,"Could not find chart "+it);

              my_dist += chart->dist(vtx[i]);

            }
            this->_dist(i, my_dist);
          }
        }
    };

    template<typename Trafo_, typename RefCellTrafo_>
    struct MeshConcentrationFunctionFactory
    {
      typedef typename Trafo_::MeshType MeshType;

      static std::shared_ptr<MeshConcentrationFunction<Trafo_, RefCellTrafo_>>
      create(const String& section_key, PropertyMap* config)
      {
        XASSERT(config != nullptr);

        std::shared_ptr<MeshConcentrationFunction<Trafo_, RefCellTrafo_>> result(nullptr);

        std::deque<String> chart_list;

        // Get the configuration section
        auto my_section = config->query_section(section_key);
        if(my_section == nullptr)
          throw InternalError(__func__,__FILE__,__LINE__,
          "Config is missing the referenced "+section_key+" section.");

        auto type_p = my_section->query("type");
        if(!type_p.second)
          throw InternalError(__func__,__FILE__,__LINE__,
          "Section "+section_key+" is missing the mandatory type key.");

        if(type_p.first== "ChartDistance")
        {
          auto chart_list_p = my_section->query("chart_list");
          if(!chart_list_p.second)
            throw InternalError(__func__,__FILE__,__LINE__,
            "Section "+section_key+" is missing the mandatory chart_list key.");

          chart_list_p.first.split_by_charset(chart_list, " ");

          auto real_result = std::make_shared<ChartDistanceFunction<Trafo_, RefCellTrafo_>>(chart_list);
          result = real_result;

        }

        return result;
      }
    };

  } // namespace Meshopt
} // namespace FEAST
#endif // KERNEL_MESHOPT_MESH_CONCENTRATION_FUNCTION_HPP
