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

    /// \cond internal

    // Forward declarations
    template<typename DT_>
    class ConcentrationFunctionDefault;

    template<typename DT_>
    class ConcentrationFunctionPowOfDist;
    /// \endcond


    /**
     * \brief Base class for mesh concentration functions
     *
     * \tparam Trafo_
     * Our transformation.
     *
     * \tparam RefCellTrafo_
     * Mesh optimisation reference cell transformation.
     *
     * This class implements the virtual interface all derived classes have to implement. A mesh concentration
     * function calculates some concentration for each cell and is then used by Meshopt::HyperelasticityFunctional
     * to compute the optimal scales h wrt. this concentration. As the concentration might depend on the DoF
     * (namely the vertex coordinates) of the mesh quality functional, this class also offers the computation of
     * the gradient of the optimal scales h wrt. these DoF.
     *
     * This class has no information about how the concentration is calculated. This is to be provided by derived
     * classes.
     *
     * \author Jordi Paul
     *
     */
    template
    <
      typename Trafo_,
      typename RefCellTrafo_ = RumpfTrafo<Trafo_, typename Trafo_::CoordType>
    >
    class MeshConcentrationFunctionBase
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
        /// Vector type for the gradient of h wrt. the DoF
        typedef LAFEM::DenseVectorBlocked<MemType, CoordType, IndexType,
        MeshType::world_dim*Shape::FaceTraits<ShapeType,0>::count> GradHType;
        /// Type of a mesh vertex
        typedef Tiny::Vector<CoordType, MeshType::world_dim> WorldPoint;

      public:
        /// \brief Virtual destructor
        virtual ~MeshConcentrationFunctionBase()
        {
        }

        /**
         * \brief Creates an empty clone of itself
         *
         * This function is handy for creating empty copies of objects of a derived class, which is used in
         * the HyperelasticityFunctional class.
         *
         * \returns A shared_ptr to a new object of this (or rather the derived) class.
         *
         */
        virtual std::shared_ptr<MeshConcentrationFunctionBase> create_empty_clone() const = 0;

        /**
         * \brief Computes the distance information
         */
        virtual void compute_dist() = 0;

        /**
         * \brief Computes the concentration for each cell according to the distance information.
         */
        virtual void compute_conc() = 0;

        /**
         * \brief Computes the local gradient of the optimal scales
         *
         * \param[in] coords
         * The current mesh vertex coordinates.
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
        virtual void compute_grad_h(const VectorType& DOXY(coords)) = 0;

        /**
         * \brief Returns a const reference to the concentration
         *
         * \returns A const reference to the concentration.
         */
        virtual const ScalarVectorType& get_conc() const = 0;

        /**
         * \brief Returns a const reference to the distance information
         *
         * \returns A const reference to the distance information in each vertex.
         */
        virtual const ScalarVectorType& get_dist() const = 0;

        /**
         * \brief Returns a const reference to the gradient of the distance function
         *
         * \returns A const reference to the gradient of the distance function in each vertex.
         */
        virtual const VectorType& get_grad_dist() const = 0;

        /**
         * \brief Returns a const reference to the concentration
         *
         * \returns A const reference to the concentration.
         */
        virtual const GradHType& get_grad_h() const = 0;

        /**
         * \brief Returns a descriptive String
         *
         * \returns The class name as String.
         */
        virtual String name() const = 0;

        /**
         * \brief Sets this object's mesh node
         *
         * \param[in] mesh_node_
         * The mesh node this function is to use.
         *
         * \note The mesh node cannot be set at construction time because this class is used from the
         * MeshConcentrationFuncionFactory, so this method is neccessary.
         */
        virtual void set_mesh_node(const Geometry::RootMeshNode<MeshType>* DOXY(mesh_node_)) = 0;

        /**
         * \brief Returns whether this function make use of the derivative of h wrt. the vertex coordinates.
         *
         * \returns True if the function makes use of the derivative of h wrt. the vertex coordinates.
         */
        virtual bool use_derivative() const = 0;

    };

    /**
     * \brief Class to compute a desired concentration for the mesh cell distribution
     *
     * \tparam ElementalFunction_
     * The scalar function that computes the concentration.
     *
     * \tparam Trafo_
     * Our transformation.
     *
     * \tparam RefCellTrafo_
     * Mesh optimisation reference cell transformation.
     *
     * \note To avoid virtual functions, this class explicitly knows which type of elemental function it has.
     *
     * The ElementalFunction_ is needed to compute grad_h etc., but we do not need to know from what quantity it maps
     * to the concentration. This is only needed in classes derived from this.
     */
    template
    <
      typename ElementalFunction_,
      typename Trafo_,
      typename RefCellTrafo_ = RumpfTrafo<Trafo_, typename Trafo_::CoordType>
    >
    class MeshConcentrationFunction : public MeshConcentrationFunctionBase<Trafo_, RefCellTrafo_>
    {
      public:
        /// The scalar function that computes the concentration from something
        typedef ElementalFunction_ ElementalFunction;
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
        /// Vector type for the gradient of h wrt. the DoF
        typedef LAFEM::DenseVectorBlocked
        <MemType, CoordType, Index, MeshType::world_dim*Shape::FaceTraits<ShapeType,0>::count> GradHType;
        /// Type for one mesh vertex
        typedef Tiny::Vector<CoordType, MeshType::world_dim> WorldPoint;

      protected:
        /// The mesh node this function works with
        const Geometry::RootMeshNode<MeshType>* _mesh_node;
        /// The scalar function mapping distance to concentration
        ElementalFunction _func;

        /// For all vertices, this holds their scalar "distance" to whatever
        ScalarVectorType _dist;
        /// For all vertices, this holds the gradient of the distance function
        VectorType _grad_dist;
        /// The sum of all entries in _conc
        CoordType _sum_conc;
        /// Vector for saving the mesh concentration, one entry per cell
        ScalarVectorType _conc;
        /// Gradient of the mesh concentration wrt. the world coordinates
        VectorType _grad_conc;
        /// Gradient of _sum_det wrt. the mesh vertices
        VectorType _grad_sum_det;
        /// Gradient of the local optimal scales h wrt. the vertex coordinates.
        // Each entry in the DenseVectorBlocked represents one cell. Each cell's block contains
        // world_dim*(number of local vertices) entries. Has to be serialised like this because there is no
        // DenseVector that saves a Tiny::Matrix
        GradHType _grad_h;

        /**
         * \brief Constructor setting the elemantal function
         *
         * \param[in] func_
         * The elemental function
         *
         * \note func_ gets copied because this constructor is called from the MeshConcentrationFunctionFactory
         */
        MeshConcentrationFunction(const ElementalFunction& func_) :
          _mesh_node(nullptr),
          _func(func_),
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
        ///**
        // *
        // */
        //explicit MeshConcentrationFunction(const MeshConcentrationFunction& other) :
        //  _mesh_node(other._mesh_node),
        //  _func(other._func),
        //  _sum_conc(other._sum_conc)
        //  {
        //    _dist.clone(other._dist, LAFEM::CloneMode::Deep);
        //    _grad_dist.clone(other._grad_dist, LAFEM::CloneMode::Deep);
        //    _grad_conc.clone(other._grad_conc, LAFEM::CloneMode::Deep);
        //    _grad_sum_det.clone(other._grad_sum_det, LAFEM::CloneMode::Deep);
        //    _grad_h.clone(other._grad_h, LAFEM::CloneMode::Deep);
        //  }

      public:
        /// \brief Virtual destructor
        virtual ~MeshConcentrationFunction()
        {
        }

        /// \copydoc BaseClass::get_conc()
        virtual const ScalarVectorType& get_conc() const override
        {
          return _conc;
        }

        /// \copydoc BaseClass::get_dist()
        virtual const ScalarVectorType& get_dist() const override
        {
          return _dist;
        }

        /// \copydoc BaseClass::get_grad_dist()
        virtual const VectorType& get_grad_dist() const override
        {
          return _grad_dist;
        }

        /// \copydoc BaseClass::get_grad_h()
        virtual const GradHType& get_grad_h() const override
        {
          return _grad_h;
        }

        /// \copydoc BaseClass::set_mesh_node()
        virtual void set_mesh_node(const Geometry::RootMeshNode<MeshType>* mesh_node_) override
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

        /// \copydoc BaseClass::name()
        virtual String name() const override
        {
          return "MeshConcentrationFunction<>";
        }

        /// \copydoc BaseClass::use_derivative()
        virtual bool use_derivative() const override
        {
          return _func.use_derivative;
        }

        /**
         * \brief Computes the mesh concentration function for each cell
         *
         * This means it computes
         * \f[
         *   \forall K \in \mathcal{T}_h: \underbrace{c(K)}_{ =: \mathtt{conc(cell)}} ~ \mathrm{and} ~ \underbrace{\sum_{K \in \mathcal{T}_h} c(K)}_{\mathtt{\_sum\_conc}}
         * \f]
         **/
        virtual void compute_conc() override
        {
          // Total number of cells in the mesh
          const Index ncells(_mesh_node->get_mesh()->get_num_entities(ShapeType::dimension));
          // Index set for local/global numbering
          const auto& idx = _mesh_node->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          this->_conc.format(CoordType(0));

          _sum_conc = CoordType(0);

          for(Index cell(0); cell < ncells; ++cell)
          {
            CoordType avg_dist(0);

            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              avg_dist += this->_dist(idx(cell,Index(j)));

            avg_dist = avg_dist/CoordType(Shape::FaceTraits<ShapeType,0>::count);

            this->_conc(cell, _func.conc_val(avg_dist));
            _sum_conc += this->_conc(cell);
          }

#ifdef FEAT_HAVE_MPI
          CoordType sum_conc_snd(_sum_conc);
          Util::Comm::allreduce(&sum_conc_snd, 1, &_sum_conc, MPI_SUM);
#endif
        }

        /**
         * \brief Computes the local gradient of the concentration function wrt. the vertices
         *
         * \tparam Tgrad_
         * Type of the local gradient, i.e. Tiny::Matrix
         *
         * \tparam Tl_
         * Type for the vector of distance values, i.e. Tiny::vector
         *
         * \tparam Tgradl_
         * Type for the local gradient of the distance wrt. the vertices, i.e. Tiny::Matrix
         *
         * \param[out] grad_loc_
         * The gradient of the concentration wrt. the vertices
         *
         * \param[in] grad_dist_loc_
         * The levelset values at the vertices
         *
         * \param[in] dist_loc_
         * The distances of the local vertices
         *
         **/
        template<typename Tgrad_, typename Tl_, typename Tgradl_>
        void compute_grad_conc_local(Tgrad_& grad_loc_, const Tl_& dist_loc_, const Tgradl_& grad_dist_loc_)
        {
          grad_loc_.format(CoordType(0));

          // This will be the average of the levelset values at the vertices
          CoordType val(0);
          for(int i(0); i < Shape::FaceTraits<ShapeType,0>::count; ++i)
            val += dist_loc_(i);

          val = val/CoordType(Shape::FaceTraits<ShapeType,0>::count);

          for(int i(0); i < Shape::FaceTraits<ShapeType,0>::count; ++i)
            grad_loc_[i] = _func.conc_der(val)/CoordType(Shape::FaceTraits<ShapeType,0>::count) * grad_dist_loc_[i];
        }

        /// \copydoc BaseClass::compute_grad_h()
        virtual void compute_grad_h(const VectorType& coords) override
        {
          XASSERT(_mesh_node != nullptr);

          _grad_h.format();

          if(_func.use_derivative)
          {

            CoordType sum_det = RefCellTrafo_::compute_sum_det(coords, *(_mesh_node->get_mesh()));
            RefCellTrafo_::compute_grad_sum_det(_grad_sum_det, coords, *(_mesh_node->get_mesh()));

            compute_grad_conc();

            // Index set for local/global numbering
            auto& idx = _mesh_node->get_mesh()->template get_index_set<ShapeType::dimension,0>();
            // This will hold the levelset values at the mesh vertices for one element
            FEAT::Tiny::Vector<CoordType, Shape::FaceTraits<ShapeType,0>::count> dist_loc;
            // This will hold the levelset gradient values for one element for passing to other routines
            FEAT::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_dist_loc;
            // This will hold the local gradient values for one element for passing to other routines
            FEAT::Tiny::Matrix<CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc(0);
            // This will hold the computed local gradient values for one element for copy assigning to the blocked
            // datatype
            FEAT::Tiny::Vector<CoordType, MeshType::world_dim*Shape::FaceTraits<ShapeType,0>::count> tmp(0);

            CoordType exponent = CoordType(1)/CoordType(MeshType::world_dim) - CoordType(1);

            for(Index cell(0); cell < _mesh_node->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
            {
              grad_loc.format();
              for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
              {
                Index i(idx(cell, Index(j)));
                // Get levelset
                dist_loc(j) = _dist(i);
                // Get levelset gradient
                grad_dist_loc[j] = _grad_dist(i);
              }

              compute_grad_conc_local(grad_loc, dist_loc, grad_dist_loc);

              for(int d(0); d < MeshType::world_dim; ++d)
              {
                for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
                {
                  Index i(idx(cell, Index(j)));

                  tmp(j*MeshType::world_dim +d) =
                    CoordType(1)/CoordType(MeshType::world_dim)*Math::pow(_conc(cell)/_sum_conc*sum_det,exponent)
                    *( _conc(cell)*(_grad_sum_det(i)(d)*_sum_conc + sum_det*_grad_conc(i)(d) )
                        + grad_loc(j,d) * sum_det *_sum_conc)
                    / Math::sqr(_sum_conc);
                }
              }
              _grad_h(cell, tmp + _grad_h(cell));

            }
          } // _func.use_derivative

        } // compute_grad_h

        /**
         * \brief Computes the gradient of the sum of all mesh concentrations
         *
         * \f[
         *   \mathrm{grad\_conc}(k,i) = \frac{\partial}{\partial x_j} \sum_{l=1}^N c(K_l),
         * \f]
         * where \f$ i \f$ is the global index of the local vertex \f$ j \f$.
         *
         */
        void compute_grad_conc()
        {
          // Clear the old gradient
          _grad_conc.format();

          // Index set for local/global numbering
          auto& idx = _mesh_node->get_mesh()->template get_index_set<ShapeType::dimension,0>();

          // This will hold the coordinates for one element for passing to other routines
          FEAT::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> x;
          // Local cell dimensions for passing to other routines
          FEAT::Tiny::Vector <CoordType,MeshType::world_dim> h;
          // This will hold the local gradient for one element for passing to other routines
          FEAT::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> grad_loc;
          // This will hold the levelset values at the mesh vertices for one element
          FEAT::Tiny::Vector <CoordType, Shape::FaceTraits<ShapeType,0>::count> grad_dist_loc;
          // This will hold the levelset gradient values for one element for passing to other routines
          FEAT::Tiny::Matrix <CoordType, Shape::FaceTraits<ShapeType,0>::count, MeshType::world_dim> dist_loc;

          // Compute the functional value for each cell
          for(Index cell(0); cell < _mesh_node->get_mesh()->get_num_entities(ShapeType::dimension); ++cell)
          {
            // Collect levelset and levelset grad values
            for(int j(0); j < Shape::FaceTraits<ShapeType,0>::count; ++j)
            {
              // Global vertex/dof index
              Index i(idx(cell, Index(j)));
              // Get levelset
              grad_dist_loc(j) = _dist(i);
              // Get levelset gradient
              dist_loc[j] = _grad_dist(i);
            }

            // Compute gradient of the concentration on this cell
            compute_grad_conc_local(grad_loc, grad_dist_loc, dist_loc);

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
     * \brief Class to compute a desired concentration for the mesh cell distribution base on distance to Charts
     *
     * \tparam ElementalFunction_
     * The scalar function that computes the concentration.
     *
     * \tparam Trafo_
     * Our transformation.
     *
     * \tparam RefCellTrafo_
     * Mesh optimisation reference cell transformation.
     *
     * \note To avoid virtual functions, this class explicitly knows which type of elemental function it has.
     *
     * The Charts are identified by a deque of Strings passed to the constructor and need to be present in the Atlas
     * of the object's MeshNode. At the moment, the sum of all distances is taken, which is to be revised in the
     * future.
     *
     */
    template
    <
      typename ElementalFunction_,
      typename Trafo_,
      typename RefCellTrafo_ = RumpfTrafo<Trafo_, typename Trafo_::CoordType>
    >
    class ChartDistanceFunction : public MeshConcentrationFunction<ElementalFunction_, Trafo_, RefCellTrafo_>
    {
      public:
        /// Type for the function mapping distance to concentration
        typedef ElementalFunction_ ElementalFunction;
        /// Type for the transformation
        typedef Trafo_ TrafoType;
        /// The mesh the transformation is defined on
        typedef typename TrafoType::MeshType MeshType;
        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;
        /// Type for one mesh vertex
        typedef Tiny::Vector<CoordType, MeshType::world_dim> WorldPoint;

        /// Our direct base class
        typedef MeshConcentrationFunction<ElementalFunction_, Trafo_, RefCellTrafo_> DirectBaseClass;
        /// Our base class
        typedef MeshConcentrationFunctionBase<Trafo_, RefCellTrafo_> BaseClass;
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
        /// List of charts to compute the distance to
        std::deque<String> _chart_list;

      public:
        /// Explicitly delete the empty default constructor
        ChartDistanceFunction() = delete;

        /**
         * \brief Constructor setting the ElementalFunction_ and the list of charts
         *
         * \param[in] func_
         * Elemental function.
         *
         * \param[in] chart_list_
         * The list of charts to compute the distance from.
         *
         */
        ChartDistanceFunction(const ElementalFunction& func_, const std::deque<String>& chart_list_):
          DirectBaseClass(func_),
          _chart_list()
          {
            XASSERTM(chart_list_.size() > size_t(0), "Empty chart list.");

            for(const auto& it:chart_list_)
              _chart_list.push_back(it);
          }

        //explicit ChartDistanceFunction(const ChartDistanceFunction& other) :
        //  BaseClass(other.BaseClass),
        //  _chart_list(other._chart_list)
        //  {
        //  }

        /// \copydoc BaseClass::~BaseClass
        virtual ~ChartDistanceFunction()
        {
        }

        /// \copydoc BaseClass::create_empty_clone()
        virtual std::shared_ptr<BaseClass> create_empty_clone() const override
        {
          std::shared_ptr<BaseClass> result(nullptr);
          result = std::make_shared<ChartDistanceFunction>(this->_func, _chart_list);
          return result;
        }

        /// \copydoc BaseClass::compute_dist()
        virtual void compute_dist() override
        {
          XASSERT(this->_mesh_node != nullptr);

          const auto& vtx = this->_mesh_node->get_mesh()->get_vertex_set();

          WorldPoint my_dist_vec(CoordType(0));
          WorldPoint tmp(CoordType(0));

          for(Index i(0); i < this->_mesh_node->get_mesh()->get_num_entities(0); ++i)
          {
            CoordType my_dist(0);
            my_dist_vec.format(CoordType(0));

            for(const auto& it:_chart_list)
            {
              auto* chart = this->_mesh_node->get_atlas()->find_mesh_chart(it);
              if(chart == nullptr)
                throw InternalError(__func__,__FILE__,__LINE__,"Could not find chart "+it);

              my_dist += Math::abs(chart->signed_dist(vtx[i], tmp));
              my_dist_vec += tmp;

            }
            this->_dist(i, my_dist);

            my_dist_vec.normalise();
            this->_grad_dist(i, my_dist_vec);
          }
        }
    }; // class ChartDistanceFunction

    /**
     * \brief Factory for MeshConcentrationFunctions
     *
     * \tparam Trafo_
     * Our transformation.
     *
     * \tparam RefCellTrafo_
     * Mesh optimisation reference cell transformation.
     *
     */
    template<typename Trafo_, typename RefCellTrafo_>
    struct MeshConcentrationFunctionFactory
    {
      /// The Meshtype
      typedef typename Trafo_::MeshType MeshType;
      /// Floating point precision of the mesh vertex coordinates
      typedef typename MeshType::CoordType CoordType;

      /**
       * \brief Creates a MeshConcentrationFunction according to a PropertyMap
       *
       * \param[in] section_key
       * Name of our section in config.
       *
       * \param[in] config
       * The PropertyMap holding our configuration
       *
       * \returns A shared_ptr to an object derived from MeshConcencentrationFunctionBase.
       */
      static std::shared_ptr<MeshConcentrationFunctionBase<Trafo_, RefCellTrafo_>>
        create(const String& section_key, PropertyMap* config)
        {
          XASSERT(config != nullptr);

          std::shared_ptr<MeshConcentrationFunctionBase<Trafo_, RefCellTrafo_>> result(nullptr);

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

            auto function_type_p = my_section->query("function_type");
            if(function_type_p.second)
            {
              // At the moment, there are only default and PowOfDist
              if(function_type_p.first == "PowOfDist")
              {
                CoordType minval(0);
                CoordType exponent(1);

                auto minval_p = my_section->query("minval");
                if(minval_p.second)
                  minval = std::stod(minval_p.first);

                auto exponent_p = my_section->query("exponent");
                if(exponent_p.second)
                  exponent = std::stod(exponent_p.first);

                typedef ConcentrationFunctionPowOfDist<CoordType> ElementalFunction;
                ElementalFunction my_func(minval, exponent, true);

                auto real_result = std::make_shared<ChartDistanceFunction<ElementalFunction, Trafo_, RefCellTrafo_>>
                  (my_func, chart_list);

                result = real_result;
              }
              else
                // Default case
              {
                typedef ConcentrationFunctionDefault<CoordType> ElementalFunction;
                ElementalFunction my_func;

                auto real_result = std::make_shared<ChartDistanceFunction<ElementalFunction, Trafo_, RefCellTrafo_>>
                  (my_func, chart_list);

                result = real_result;
              }
            }

          }

          XASSERTM(result != nullptr, "Unknown conc_function type!");

          return result;
        }
    };

    /**
     * \brief Default elemental distance concentration function
     *
     * \tparam DT_ Floating point precision
     *
     * This maps a signed distance to a concentration, leaving the derivative zero.
     */
    template<typename DT_>
    class ConcentrationFunctionDefault
    {
      public:
        /// Floating point precision
        typedef DT_ DataType;

        /// Does this utilise the derivative of the input wrt. some other DoF?
        const bool use_derivative;

        /**
         * \brief Default constructor
         *
         * This sets use_derivative to false.
         */
        explicit ConcentrationFunctionDefault():
          use_derivative(false)
        {
        }

        /**
         * \brief Copy constructor
         *
         * \param[in] other
         * Object to be copied.
         *
         */
        ConcentrationFunctionDefault(const ConcentrationFunctionDefault& other) :
          use_derivative(other.use_derivative)
        {
        }

        /**
         * \brief Empty destructor
         *
         * The shall be no classes derived from this, so this is non-virtual.
         *
         */
        ~ConcentrationFunctionDefault()
        {
        }

        /**
         * \brief Computes the concentration according to a distance.
         *
         * \param[in] dist
         * The input parameter
         *
         * \returns \f$ c(d) = |d|\f$
         */
        DataType conc_val(DT_ dist) const
        {
          return Math::abs(dist);
        }

        /**
         * \brief Computes the derivative of the concentration according to a distance.
         *
         * \param[in] dist
         * The input parameter
         *
         * \returns 0.
         */
        DataType conc_der(DT_ DOXY(dist)) const
        {
          return DataType(0);
        }
    };

    /**
     * \brief Default elemental distance concentration function
     *
     * \tparam DT_ Floating point precision
     *
     * This maps a signed distance to a concentration according to
     * \f[
     *    c(d) = (\alpha + |d|)^\beta, \alpha \geq 0.
     * \f]
     */
    template<typename DT_>
    class ConcentrationFunctionPowOfDist
    {
      public:
        /// Floating point precision
        typedef DT_ DataType;
        /// Does this utilise the derivative of the input wrt. some other DoF?
        const bool use_derivative;

      private:
        /// alpha
        const DataType _minval;
        /// beta
        const DataType _exponent;

      public:
        /**
         * \brief Constructor setting alpha, beta and use_derivative
         *
         * \param[in] minval_
         * alpha
         *
         * \param[in] exponent_
         * beta
         *
         * \param[in] use_derivative_
         * Whether to the derivative.
         *
         */
        explicit ConcentrationFunctionPowOfDist(DataType minval_, DataType exponent_, bool use_derivative_ = true):
          use_derivative(use_derivative_),
          _minval(minval_),
          _exponent(exponent_)
          {
          }

        //ConcentrationFunctionPowOfDist(const ConcentrationFunctionPowOfDist& other) :
        //  use_derivative(other.use_derivative),
        //  _minval(other._minval),
        //  _exponent(other._exponent)
        //  {
        //  }

        /**
         * \brief Empty destructor
         *
         * The shall be no classes derived from this, so this is non-virtual.
         *
         */
        ~ConcentrationFunctionPowOfDist()
        {
        }

        /**
         * \brief Computes the concentration according to a distance.
         *
         * \param[in] dist
         * The input parameter
         *
         * \returns \f$ c(d) = (\alpha |d|)^\beta \f$
         */
        DataType conc_val(DataType dist) const
        {
          return Math::pow(_minval + Math::abs(dist), _exponent);
        }

        /**
         * \brief Computes the concentration according to a distance.
         *
         * \param[in] dist
         * The input parameter
         *
         * \returns \f$ c'(d) = \beta (\alpha |d|)^{(\beta-1)} * sign(d) \f$
         */
        DataType conc_der(DataType dist) const
        {
          return _exponent*Math::pow(_minval + Math::abs(dist), _exponent - DataType(1))*Math::signum(dist);
        }
    };
  } // namespace Meshopt
} // namespace FEAST
#endif // KERNEL_MESHOPT_MESH_CONCENTRATION_FUNCTION_HPP
