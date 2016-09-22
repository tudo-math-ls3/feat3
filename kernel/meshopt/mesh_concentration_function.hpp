#pragma once
#ifndef KERNEL_MESHOPT_MESH_CONCENTRATION_FUNCTION_HPP
#define KERNEL_MESHOPT_MESH_CONCENTRATION_FUNCTION_HPP 1

#include <kernel/base_header.hpp>
#include <kernel/shape.hpp>

#include <kernel/analytic/common.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/intern/face_index_mapping.hpp>
#include <kernel/meshopt/rumpf_trafo.hpp>
#include <kernel/util/assertion.hpp>
#include <kernel/util/comm_base.hpp>
#include <kernel/util/property_map.hpp>
#include <kernel/util/mpi_cout.hpp>

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

    template<typename DT_, typename ShapeType_>
    struct AlignmentPenalty;

    template<typename DT_>
    struct AlignmentPenalty<DT_, Shape::Simplex<2>>
    {
      typedef DT_ DataType;
      typedef Shape::Simplex<2> ShapeType;

      template<typename Mesh_, typename Dist_, typename EdgeFreqs_>
      static DataType compute_constraint(const Mesh_& mesh, const Dist_& dist, const EdgeFreqs_& edge_freqs)
      {
        DataType constraint(0);
        const auto edge_idx = mesh.template get_index_set<1,0>();
        for(Index edge(0); edge < mesh.get_num_entities(1); ++edge)
        {
          constraint += edge_freqs(edge)*FEAT::Analytic::Common::template HeavisideRegStatic<DataType>::
            eval( - dist(edge_idx(edge,0)) * dist(edge_idx(edge,1)));
        }

#ifdef FEAT_HAVE_MPI
        Util::Comm::allreduce(&constraint, &constraint, 1, Util::CommOperationSum());
#endif
        // DEBUG
        return DataType(2)*constraint;
      }

      template<typename Mesh_, typename Dist_>
      static DataType compute_constraint(DataType* constraint_vec, const Mesh_& mesh, const Dist_& dist)
      {
        XASSERT(constraint_vec != nullptr);
        typedef Geometry::Intern::FaceIndexMapping<ShapeType, 1, 0> FimType;

        DataType constraint(0);
        const auto idx = mesh.template get_index_set<ShapeType::dimension,0>();
        for(Index cell(0); cell < mesh.get_num_entities(ShapeType::dimension); ++cell)
        {
          constraint_vec[cell] = DataType(0);
          for(int edge(0); edge < Shape::FaceTraits<ShapeType,1>::count; ++edge)
          {
            //These are the indices of the vertices on edge
            int i(FimType::map(edge,0));
            int j(FimType::map(edge,1));
            DataType my_constraint(FEAT::Analytic::Common::template HeavisideRegStatic<DataType>::
                eval(- DataType(1)*dist(idx(cell, Index(i))) * dist(idx(cell, Index(j)))) );

            constraint_vec[cell] += my_constraint;

            constraint += my_constraint;
          }

        }

#ifdef FEAT_HAVE_MPI
        Util::Comm::allreduce(&constraint, &constraint, 1, Util::CommOperationSum());
#endif
        return constraint;
      }

      template<typename Vector_, typename Mesh_, typename Dist_, typename GradDist_, typename EdgeFreqs_>
      static void add_constraint_grad(
        Vector_& grad, const DataType alignment_fval, const DataType fac, const Mesh_& mesh, const Dist_& dist, const GradDist_& grad_dist, const EdgeFreqs_& edge_freqs)
      {
        const auto edge_idx = mesh.template get_index_set<1,0>();
        /// Type of a mesh vertex
        typedef Tiny::Vector<DataType, Mesh_::world_dim> WorldPoint;

        WorldPoint grad_loc(DataType(0));

        for(Index edge(0); edge < mesh.get_num_entities(1); ++edge)
        {
          Index i(edge_idx(edge,0));
          Index j(edge_idx(edge,1));

          auto dist_prod =  dist(i) * dist(j);
          // Derivative of the heaviside function
          auto heaviside_der = FEAT::Analytic::Common::template HeavisideRegStatic<DataType>::der_x(-dist_prod);

          grad_loc =  (-fac * alignment_fval * heaviside_der * dist(j)) * grad_dist(i);
          // DEBUG
          grad(i, grad(i) + DataType(2)*edge_freqs(edge)*grad_loc);

          grad_loc =  (-fac * alignment_fval * heaviside_der * dist(i)) * grad_dist(j);
          // DEBUG
          grad(j, grad(j) + DataType(2)*edge_freqs(edge)*grad_loc);
        }
      }

    };

    template<typename DT_>
    struct AlignmentPenalty<DT_, Shape::Hypercube<2>>
    {

      typedef DT_ DataType;
      typedef Shape::Hypercube<2> ShapeType;

      template<typename Mesh_, typename Dist_, typename EdgeFreqs_>
      static DataType compute_constraint(const Mesh_& mesh, const Dist_& dist, const EdgeFreqs_& edge_freqs)
      {
        DataType constraint(0);
        const auto edge_idx = mesh.template get_index_set<1,0>();
        for(Index edge(0); edge < mesh.get_num_entities(1); ++edge)
        {
          constraint += edge_freqs(edge)*FEAT::Analytic::Common::template HeavisideRegStatic<DataType>::
            eval( - dist(edge_idx(edge,0)) * dist(edge_idx(edge,1)));
        }

        const auto cell_idx = mesh.template get_index_set<ShapeType::dimension,0>();
        for(Index cell(0); cell < mesh.get_num_entities(ShapeType::dimension); ++cell)
        {
          constraint += FEAT::Analytic::Common::template HeavisideRegStatic<DataType>::
            eval( - dist(cell_idx(cell,0)) * dist(cell_idx(cell,3)));
          constraint += FEAT::Analytic::Common::template HeavisideRegStatic<DataType>::
            eval( - dist(cell_idx(cell,1)) * dist(cell_idx(cell,2)));
        }

#ifdef FEAT_HAVE_MPI
        Util::Comm::allreduce(&constraint, &constraint, 1, Util::CommOperationSum());
#endif

        // DEBUG
        return DataType(2)*constraint;
      }

      template<typename Mesh_, typename Dist_>
      static DataType compute_constraint(DataType* constraint_vec, const Mesh_& mesh, const Dist_& dist)
      {
        XASSERT(constraint_vec != nullptr);
        typedef Geometry::Intern::FaceIndexMapping<ShapeType, 1, 0> FimType;

        DataType constraint(0);
        const auto idx = mesh.template get_index_set<ShapeType::dimension,0>();
        for(Index cell(0); cell < mesh.get_num_entities(ShapeType::dimension); ++cell)
        {
          constraint_vec[cell] = DataType(0);
          for(int edge(0); edge < Shape::FaceTraits<ShapeType,1>::count; ++edge)
          {
            //These are the indices of the vertices on edge
            int i(FimType::map(edge,0));
            int j(FimType::map(edge,1));
            DataType my_constraint(FEAT::Analytic::Common::template HeavisideRegStatic<DataType>::
                eval(- DataType(1)*dist(idx(cell, Index(i))) * dist(idx(cell, Index(j)))) );

            constraint_vec[cell] += my_constraint;

            constraint += my_constraint;
          }

          DataType my_constraint(FEAT::Analytic::Common::template HeavisideRegStatic<DataType>::
              eval( - dist(idx(cell,0)) * dist(idx(cell,3))));
          my_constraint = FEAT::Analytic::Common::template HeavisideRegStatic<DataType>::
            eval( - dist(idx(cell,1)) * dist(idx(cell,2)));

          constraint_vec[cell] += my_constraint;

          constraint += my_constraint;

        }

        return constraint;
      }

      template<typename Vector_, typename Mesh_, typename Dist_, typename GradDist_, typename EdgeFreqs_>
      static void add_constraint_grad(
        Vector_& grad, const DataType alignment_fval, const DataType fac, const Mesh_& mesh, const Dist_& dist,
        const GradDist_& grad_dist, const EdgeFreqs_& edge_freqs)
      {
        const auto edge_idx = mesh.template get_index_set<1,0>();
        /// Type of a mesh vertex
        typedef Tiny::Vector<DataType, Mesh_::world_dim> WorldPoint;

        WorldPoint grad_loc(DataType(0));

        for(Index edge(0); edge < mesh.get_num_entities(1); ++edge)
        {
          Index i(edge_idx(edge,0));
          Index j(edge_idx(edge,1));

          auto dist_prod =  dist(i) * dist(j);
          // Derivative of the heaviside function
          auto heaviside_der = FEAT::Analytic::Common::template HeavisideRegStatic<DataType>::der_x(-dist_prod);

          grad_loc =  (-fac * alignment_fval * heaviside_der * dist(j)) * grad_dist(i);
          // DEBUG
          grad(i, grad(i) + edge_freqs(edge)*DataType(2)*grad_loc);

          grad_loc =  (-fac * alignment_fval * heaviside_der * dist(i)) * grad_dist(j);
          // DEBUG
          grad(j, grad(j) + edge_freqs(edge)*DataType(2)*grad_loc);
        }

        const auto cell_idx = mesh.template get_index_set<ShapeType::dimension,0>();
        for(Index cell(0); cell < mesh.get_num_entities(ShapeType::dimension); ++cell)
        {
          Index i(cell_idx(cell,Index(0)));
          Index j(cell_idx(cell,Index(3)));

          auto dist_prod = dist(i) * dist(j);
          // Derivative of the heaviside function
          auto heaviside_der = FEAT::Analytic::Common::template HeavisideRegStatic<DataType>::der_x(-dist_prod);

          grad_loc =  (-fac * alignment_fval * heaviside_der * dist(j)) * grad_dist(i);
          // DEBUG
          grad(i, grad(i)+DataType(2)*grad_loc);

          grad_loc =  (-fac * alignment_fval * heaviside_der * dist(i)) * grad_dist(j);
          // DEBUG
          grad(j, grad(j)+DataType(2)*grad_loc);

          i = cell_idx(cell,Index(1));
          j = cell_idx(cell,Index(2));

          dist_prod = dist(i) * dist(j);
          // Derivative of the heaviside function
          heaviside_der = FEAT::Analytic::Common::template HeavisideRegStatic<DataType>::der_x(-dist_prod);

          grad_loc =  (-fac * alignment_fval * heaviside_der * dist(j)) * grad_dist(i);
          // DEBUG
          grad(i, grad(i)+DataType(2)*grad_loc);

          grad_loc =  (-fac * alignment_fval * heaviside_der * dist(i)) * grad_dist(j);
          // DEBUG
          grad(j, grad(j)+DataType(2)*grad_loc);
        }
      }
    };

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
        /// ShapeType of said mesh
        typedef typename MeshType::ShapeType ShapeType;
        /// The precision of the mesh coordinates
        typedef typename MeshType::CoordType CoordType;

        /// Only Mem::Main is supported atm
        typedef Mem::Main MemType;
        /// The index type
        typedef Index IndexType;

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
         * \brief Prints relevant information
         */
        virtual void print() const = 0;

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
         * \brief Computes the gradient of the mesh concentration for each cell according to the distance information.
         */
        virtual void compute_grad_conc() = 0;

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
         * \brief Returns a reference to gradient of the optimal scales wrt. the vertex coordinates
         *
         * \returns A reference to gradient of the optimal scales wrt. the vertex coordinates
         */
        virtual GradHType& get_grad_h() = 0;

        /**
         * \brief Returns a const reference to gradient of the optimal scales wrt. the vertex coordinates
         *
         * \returns A const reference to gradient of the optimal scales wrt. the vertex coordinates
         */
        virtual const GradHType& get_grad_h() const = 0;

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

        /**
         * \brief Computes the surface alignment constraint
         *
         * \returns The surface alignment constraint
         */
        virtual CoordType compute_constraint() const = 0;

        /**
         * \brief Computes the surface alignment constraint at every vertex
         *
         * \param[out] constraint_at_vtx
         * The vector of constraints at every vertex.
         *
         * The is for post processing and debugging purposes.
         *
         */
        virtual CoordType compute_constraint(CoordType* DOXY(constraint_at_vtx)) const = 0;

        /**
         * \brief Adds the scaled gradient of the constraint wrt. the vertex coordinates
         *
         * \param[in,out] grad
         * Vector to add the constraint gradient to.
         *
         * \param[in] constraint
         * The current constraint value.
         *
         * \param[in] fac
         * grad <- grad + fac * constraint_grad
         *
         */
        virtual void add_constraint_grad(VectorType& DOXY(grad), const CoordType DOXY(constraint), const CoordType DOXY(fac)) const = 0;

        /**
         * \brief Computes the gradient of the sum of all determinants of the trafo
         *
         * \param[in] coords
         * Set of coordinates to compute the gradient for
         *
         */
        virtual void compute_grad_sum_det(const VectorType& DOXY(coords)) = 0;

        /**
         * \brief Adds pointers to vectors that need synchronising (type-0 to type-1 vectors)
         *
         * \param[in,out] sync_vecs
         * deque holding pointers to all vector that need synchronising.
         *
         */
        virtual void add_sync_vecs(std::deque<VectorType*>& DOXY(sync_vecs)) = 0;

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
        /// Surface alignment penalty function
        typedef AlignmentPenalty<CoordType, ShapeType> PenaltyFunction;

      protected:
        /// The mesh node this function works with
        const Geometry::RootMeshNode<MeshType>* _mesh_node;
        /// The scalar function mapping distance to concentration
        ElementalFunction _func;

        /// For each edge, this contains 1/(# halos it is present in)
        ScalarVectorType _edge_freqs;
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
          _edge_freqs(),
          _dist(),
          _grad_dist(),
          _sum_conc(CoordType(0)),
          _conc(),
          _grad_conc(),
          _grad_sum_det(),
          _grad_h()
          {
          }

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

        /// \copydoc BaseClass::add_sync_vecs()
        virtual void add_sync_vecs(std::deque<VectorType*>& sync_vecs) override
        {
          if(_func.use_derivative)
          {
            sync_vecs.push_back(&_grad_sum_det);
            sync_vecs.push_back(&_grad_conc);
          }
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
        virtual GradHType& get_grad_h() override
        {
          return _grad_h;
        }

        /// \copydoc BaseClass::get_grad_h()
        virtual const GradHType& get_grad_h() const override
        {
          return _grad_h;
        }

        /// copydoc BaseClass::print()
        virtual void print() const override
        {
          _func.print();
        }

        /// \copydoc BaseClass::set_mesh_node()
        virtual void set_mesh_node(const Geometry::RootMeshNode<MeshType>* mesh_node_) override
        {
          XASSERT(_mesh_node == nullptr);

          _mesh_node = mesh_node_;

          _sum_conc = CoordType(0);

          Index ndofs(_mesh_node->get_mesh()->get_num_entities(0));
          Index nedges(_mesh_node->get_mesh()->get_num_entities(1));
          Index ncells(_mesh_node->get_mesh()->get_num_entities(ShapeType::dimension));

          _edge_freqs = ScalarVectorType(nedges, CoordType(1));
          _dist = ScalarVectorType(ndofs, CoordType(0));
          _grad_dist = VectorType(ndofs,CoordType(0));
          _conc = ScalarVectorType(ncells);
          _grad_conc = VectorType(ndofs, CoordType(0));
          _grad_h = GradHType(ncells,CoordType(0));
          _grad_sum_det = VectorType(ndofs, CoordType(0));

          for(const auto& it: _mesh_node->get_mesh_part_names())
          {
            if(it.starts_with("_halo"))
            {
              const auto& edge_ts = _mesh_node->find_mesh_part(it)->template get_target_set<1>();
              for(Index edge(0); edge < edge_ts.get_num_entities(); ++edge)
                _edge_freqs(edge_ts[edge], _edge_freqs(edge_ts[edge])+CoordType(1));
            }
          }
          _edge_freqs.component_invert(_edge_freqs, CoordType(1));
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
          Util::Comm::allreduce(&sum_conc_snd, &_sum_conc, 1, Util::CommOperationSum());
#endif
        }

        /// \copydoc BaseClass:compute_grad_sum_det()
        virtual void compute_grad_sum_det(const VectorType& coords) override
        {
          if(_func.use_derivative)
          {
            RefCellTrafo_::compute_grad_sum_det(_grad_sum_det, coords, *(_mesh_node->get_mesh()));
          }
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

            CoordType exponent(CoordType(1)/CoordType(MeshType::world_dim) - CoordType(1));

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
        virtual void compute_grad_conc() override
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

        /// \copydoc BaseClass::compute_constraint()
        virtual CoordType compute_constraint() const override
        {
            return PenaltyFunction::compute_constraint(*(_mesh_node->get_mesh()), this->_dist, this->_edge_freqs);
        }

        /// \copydoc BaseClass::compute_constraint()
        virtual CoordType compute_constraint(CoordType* constraint_vec) const override
        {
          XASSERT(constraint_vec != nullptr);
          return PenaltyFunction::compute_constraint(constraint_vec, *(_mesh_node->get_mesh()), this->_dist);
        }

        /// \copydoc BaseClass::add_constraint_grad()
        virtual void add_constraint_grad(
          VectorType& grad, const CoordType alignment_fval, const CoordType penalty_param) const override
        {
          PenaltyFunction::add_constraint_grad(
            grad, alignment_fval, penalty_param, *(_mesh_node->get_mesh()), _dist, _grad_dist, this->_edge_freqs);
        }

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
        /// How to mangle more than one distance (add, max, min)
        const String _operation;

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
        ChartDistanceFunction(const ElementalFunction& func_, const std::deque<String>& chart_list_, const String& operation_):
          DirectBaseClass(func_),
          _chart_list(),
          _operation(operation_)
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
          result = std::make_shared<ChartDistanceFunction>(this->_func, _chart_list, _operation);
          return result;
        }

        /// \copydoc BaseClass::compute_dist()
        virtual void compute_dist() override
        {
          if(_operation == "add")
            compute_dist_add();
          else if(_operation == "max")
            compute_dist_max();
          else if(_operation == "min")
            compute_dist_min();
          else
            throw InternalError(__func__,__FILE__,__LINE__,"Unknown operation "+_operation);
        }

        /// \copydoc BaseClass::compute_dist()
        void compute_dist_add()
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

            // Because we added distance function gradient vectors, we have to normalise again if possible
            CoordType my_norm(my_dist_vec.norm_euclid());
            if(my_norm > Math::eps<CoordType>())
              my_dist_vec *= (CoordType(1)/my_norm);

            this->_grad_dist(i, my_dist_vec);
          }
        }

        /// \copydoc BaseClass::compute_dist()
        void compute_dist_max()
        {
          XASSERT(this->_mesh_node != nullptr);

          const auto& vtx = this->_mesh_node->get_mesh()->get_vertex_set();

          WorldPoint my_dist_vec(CoordType(0));
          WorldPoint tmp(CoordType(0));

          for(Index i(0); i < this->_mesh_node->get_mesh()->get_num_entities(0); ++i)
          {
            CoordType my_dist(-Math::huge<CoordType>());
            my_dist_vec.format(CoordType(0));

            for(const auto& it:_chart_list)
            {
              auto* chart = this->_mesh_node->get_atlas()->find_mesh_chart(it);
              if(chart == nullptr)
                throw InternalError(__func__,__FILE__,__LINE__,"Could not find chart "+it);

              CoordType this_dist = chart->signed_dist(vtx[i], tmp);
              if(this_dist > my_dist)
              {
                my_dist = this_dist;
                my_dist_vec = tmp;
              }

            }

            this->_dist(i, my_dist);
            this->_grad_dist(i, my_dist_vec);
          }
        }
        /// \copydoc BaseClass::compute_dist()
        void compute_dist_min()
        {
          XASSERT(this->_mesh_node != nullptr);

          const auto& vtx = this->_mesh_node->get_mesh()->get_vertex_set();

          WorldPoint my_dist_vec(CoordType(0));
          WorldPoint tmp(CoordType(0));

          for(Index i(0); i < this->_mesh_node->get_mesh()->get_num_entities(0); ++i)
          {
            CoordType my_dist(Math::huge<CoordType>());
            my_dist_vec.format(CoordType(0));

            for(const auto& it:_chart_list)
            {
              auto* chart = this->_mesh_node->get_atlas()->find_mesh_chart(it);
              if(chart == nullptr)
                throw InternalError(__func__,__FILE__,__LINE__,"Could not find chart "+it);

              CoordType this_dist = chart->signed_dist(vtx[i], tmp);
              if(Math::abs(this_dist) < Math::abs(my_dist))
              {
                my_dist = this_dist;
                my_dist_vec = tmp;
              }

            }

            this->_dist(i, my_dist);
            this->_grad_dist(i, my_dist_vec);
          }
        }

        static String name()
        {
          return "ChartDistanceFunction<"+ElementalFunction::name()+">";
        }

        /// \copydoc BaseClass::print()
        virtual void print() const override
        {
          Util::mpi_cout(name()+" settings:\n");
          DirectBaseClass::print();
          Util::mpi_cout_pad_line("Operation:",_operation);
          for(const auto& it:_chart_list)
            Util::mpi_cout_pad_line("DistanceChart:",it);
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

            String operation("min");
            auto operation_p = my_section->query("operation");
            if(operation_p.second)
              operation = operation_p.first;

            chart_list_p.first.split_by_charset(chart_list, " ");

            auto function_type_p = my_section->query("function_type");
            if(function_type_p.second)
            {
              // At the moment, there are only default and PowOfDist
              if(function_type_p.first == "PowOfDist")
              {
                CoordType minval(0);
                CoordType exponent(1);
                bool use_derivative(true);

                auto minval_p = my_section->query("minval");
                if(minval_p.second)
                  minval = std::stod(minval_p.first);

                auto exponent_p = my_section->query("exponent");
                if(exponent_p.second)
                  exponent = std::stod(exponent_p.first);

                auto use_derivative_p = my_section->query("use_derivative");
                if(use_derivative_p.second)
                  use_derivative = (std::stoi(use_derivative_p.first) == 1);

                typedef ConcentrationFunctionPowOfDist<CoordType> ElementalFunction;
                ElementalFunction my_func(minval, exponent, use_derivative);

                auto real_result = std::make_shared<ChartDistanceFunction<ElementalFunction, Trafo_, RefCellTrafo_>>
                  (my_func, chart_list, operation);

                result = real_result;
              }
              // Default case
              else
              {
                typedef ConcentrationFunctionDefault<CoordType> ElementalFunction;
                ElementalFunction my_func;

                auto real_result = std::make_shared<ChartDistanceFunction<ElementalFunction, Trafo_, RefCellTrafo_>>
                  (my_func, chart_list, operation);

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
         */
        virtual ~ConcentrationFunctionDefault()
        {
        }

        /**
         * \brief Returns a descriptive String
         *
         * \returns The class name as String.
         */
        static String name()
        {
          return "ConcentrationFunctionDefault";
        }

        /**
         * \brief Prints relevant information
         */
        virtual void print() const
        {
          Util::mpi_cout(name()+" settings:\n");
          Util::mpi_cout_pad_line("Function:","c(d) = |d|");
          Util::mpi_cout_pad_line("use_derivative:",use_derivative);
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
         */
        virtual ~ConcentrationFunctionPowOfDist()
        {
        }

        /**
         * \brief Returns a descriptive String
         *
         * \returns The class name as String.
         */
        static String name()
        {
          return "ConcentrationFunctionPowOfDist";
        }

        /**
         * \brief Prints relevant information
         */
        virtual void print() const
        {
          Util::mpi_cout(name()+" settings:\n");
          Util::mpi_cout_pad_line("Function:","c(d) = (alpha + |d|)^beta");
          Util::mpi_cout_pad_line("alpha:",_minval);
          Util::mpi_cout_pad_line("beta:",_exponent);
          Util::mpi_cout_pad_line("use_derivative:",use_derivative);
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
