// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#pragma once
#include "kernel/assembly/asm_traits.hpp"
#include <kernel/assembly/surface_integrator.hpp>
#include <kernel/base_header.hpp>
#include <kernel/trafo/inverse_mapping.hpp>
#include <kernel/util/tiny_algebra.hpp>

#include <unordered_map>



namespace FEAT::Assembly
{

  template<typename AsmTraits_>
  class SurfaceIntegratorTaskBase
  {
  public:
    typedef typename AsmTraits_::TrafoType TrafoType;
    static constexpr int dim = TrafoType::ShapeType::dimension;
    typedef Index IndexType;
    typedef typename AsmTraits_::DataType DataType;
    typedef AsmTraits_ AsmTraits;
    typedef typename AsmTraits::TrafoEvaluator TrafoEval;
    typedef typename TrafoEval::DomainPointType DomainPointType;
    typedef typename TrafoEval::ImagePointType ImagePointType;

  protected:
    const TrafoType& _trafo;
    Trafo::InverseMapping<TrafoType, DataType> _inv_mapping;

    std::vector<IndexType> _cell_helper;
    std::vector<BoundingBoxType<DataType, dim>> _cell_bb;
    std::vector<ImagePointType> _domain_points;
    std::vector<DataType> _point_weights;
    std::vector<std::vector<IndexType>> _cell_to_domain_point;

    ImagePointType _cur_img_point;
    ImagePointType _cur_dom_point;
    ImagePointType _normal;
    IndexType _cur_surface_index;

    DataType _integration_weight;

  public:
    SurfaceIntegratorTaskBase(const TrafoType& trafo_)
     : _trafo(trafo_),
       _inv_mapping(_trafo, DataType(1E-2), DataType(1E-4), false),
      _cell_helper(),
      _cell_bb(),
      _domain_points(),
      _point_weights(),
      _cell_to_domain_point(),
      _cur_img_point(),
      _cur_dom_point(),
      _normal(),
      _cur_surface_index(~IndexType(0)),
      _integration_weight(DataType(0))
    {}

    template<typename DT_, typename ImgP_, typename IT_>
    void prepare(const std::vector<IT_>& cells, const std::vector<IT_>& cell_offsets, const std::vector<ImgP_>& points,
      const std::vector<DT_>& integration_weights, const ImgP_& normal, IT_ surface_ind)
    {
      _normal = ImagePointType::convert_new(normal);
      _cell_helper.resize(cell_offsets.at(surface_ind+1)-cell_offsets.at(surface_ind));
      std::copy(cells.begin()+cell_offsets.at(surface_ind), cells.begin()+cell_offsets.at(surface_ind+1), _cell_helper.begin());
      _cur_surface_index = surface_ind;
      _cell_bb.resize(_cell_helper.size());

      const auto& vtx_set = _trafo.get_mesh().get_vertex_set();
      const auto& vtx_ind = _trafo.get_mesh().template get_index_set<dim, 0>();
      for(IndexType k = 0; k < _cell_helper.size(); ++k)
      {
        _cell_bb[k] = Intern::get_boundingbox<DataType, dim>(vtx_set, vtx_ind[_cell_helper[k]], DataType(1E-4));
      }

      _cell_to_domain_point.resize(_cell_helper.size());
      std::for_each(_cell_to_domain_point.begin(), _cell_to_domain_point.end(), [](auto& ele){ele.clear();});
      _domain_points.clear();
      _point_weights.clear();

      std::unordered_map<IndexType, std::size_t> cell_to_index_map;
      for(std::size_t k = 0; k < _cell_helper.size(); ++k)
      {
        cell_to_index_map.insert({_cell_helper[k], k});
      }

      std::vector<IndexType> pt_helper;
      for(std::size_t pti = 0; pti < points.size(); ++pti)
      {
        pt_helper.clear();
        for(IndexType k = 0; k < _cell_helper.size(); ++k)
        {
          if(Intern::check_point(ImagePointType(points[pti]), _cell_bb[k]))
          {
            pt_helper.push_back(_cell_helper[k]);
          }
        }
        const auto& inv_mapping_data = _inv_mapping.unmap_point(ImagePointType(points[pti]), pt_helper, true);


        // extract the cell -> point information from the inv mapping
        const DataType point_weight = DataType(integration_weights[pti])/DataType(inv_mapping_data.size());
        for(std::size_t l = 0; l < inv_mapping_data.size(); ++l)
        {
          _domain_points.push_back(inv_mapping_data.dom_points[l]);
          _point_weights.push_back(point_weight);
          _cell_to_domain_point.at(cell_to_index_map[inv_mapping_data.cells[l]]).push_back(_point_weights.size()-1);
        }
      }
    }


    void prepare_point(const DomainPointType& dom_point, DataType point_weight)
    {
      _cur_dom_point = dom_point;
      _integration_weight = point_weight;
    }


    // to be implemented by child class
    void assemble()
    {
    }

    void scatter()
    {
    }

    void combine()
    {
    }


  };

  namespace Intern
  {
    template<typename AsmTraits_, typename DataType_, int domain_dim_, int image_dim_>
    struct LocalFEValueHolder
    {
      struct Empty {};
      std::conditional_t<*(AsmTraits_::space_config & SpaceTags::value), Tiny::Vector<DataType_, image_dim_>, Empty> value;
      std::conditional_t<*(AsmTraits_::space_config & SpaceTags::grad), Tiny::Matrix<DataType_, image_dim_, domain_dim_>, Empty> grad;
      std::conditional_t<*(AsmTraits_::space_config & SpaceTags::hess), Tiny::Tensor3<DataType_, image_dim_, domain_dim_, domain_dim_>, Empty> hess;

      static constexpr bool has_value = *(AsmTraits_::space_config & SpaceTags::value);
      static constexpr bool has_grad = *(AsmTraits_::space_config & SpaceTags::grad);
      static constexpr bool has_hess = *(AsmTraits_::space_config & SpaceTags::hess);
    };

    template<typename AsmTraits_, typename DataType_, int domain_dim_>
    struct LocalFEValueHolder<AsmTraits_, DataType_, domain_dim_, 1>
    {
      struct Empty {};
      std::conditional_t<*(AsmTraits_::space_config & SpaceTags::value), DataType_, Empty> value;
      std::conditional_t<*(AsmTraits_::space_config & SpaceTags::grad), Tiny::Vector<DataType_, domain_dim_>, Empty> grad;
      std::conditional_t<*(AsmTraits_::space_config & SpaceTags::hess), Tiny::Matrix<DataType_, domain_dim_, domain_dim_>, Empty> hess;

      static constexpr bool has_value = *(AsmTraits_::space_config & SpaceTags::value);
      static constexpr bool has_grad = *(AsmTraits_::space_config & SpaceTags::grad);
      static constexpr bool has_hess = *(AsmTraits_::space_config & SpaceTags::hess);
    };

    template<typename ValueType_>
    struct ValueTypeHelper
    {
      static constexpr int dim = 1;
    };

    template<typename DT_, int n_, int sn_>
    struct ValueTypeHelper<Tiny::Vector<DT_, n_, sn_>>
    {
      static constexpr int dim = n_;
    };
  }

  template<typename Derived_, typename AsmTraits_, typename Space_, typename FEVector_>
  class FEVectorIntegratorTaskCRTP :
    public SurfaceIntegratorTaskBase<AsmTraits_>
  {
  public:
    typedef SurfaceIntegratorTaskBase<AsmTraits_> BaseClass;
    typedef AsmTraits_ AsmTraits;
    typedef typename BaseClass::IndexType IndexType;
    typedef typename BaseClass::DomainPointType DomainPointType;
    static constexpr int dim = BaseClass::dim;
    typedef typename AsmTraits_::DataType DataType;
    typedef typename FEVector_::ValueType ValueType;
    static constexpr int image_dim = Intern::ValueTypeHelper<ValueType>::dim;

  protected:
    Derived_& cast()
    {
      return static_cast<Derived_&>(*this);
    }

    const Derived_& cast() const
    {
      return static_cast<const Derived_&>(*this);
    }

  public:
    /// the vector from which to gather the values
    const FEVector_& primal_vector;
    /// the test-/trial-space to be used
    const Space_& space;
    /// the trafo evaluator
    typename AsmTraits::TrafoEvaluator trafo_eval;
    /// the space evaluator
    typename AsmTraits::SpaceEvaluator space_eval;
    /// the space dof-mapping
    typename AsmTraits::DofMapping dof_mapping;
    /// the trafo evaluation data
    typename AsmTraits::TrafoEvalData trafo_data;
    /// the space evaluation data
    typename AsmTraits::SpaceEvalData space_data;
    /// the local vector to be gathered
    typename AsmTraits::template TLocalVector<ValueType> local_vector;
    /// the gather object
    typename FEVector_::GatherAxpy gather_axpy;
    /// local fe point/grad/hess values
    Intern::LocalFEValueHolder<AsmTraits, DataType, dim, image_dim> loc_value_holder;

  public:
    /**
      * \brief Constructor
      */
    explicit FEVectorIntegratorTaskCRTP(const Space_& space_, const FEVector_& primal_vec_) :
              BaseClass(space_.get_trafo()),
              primal_vector(primal_vec_),
              space(space_),
              trafo_eval(this->_trafo),
              space_eval(space),
              dof_mapping(space),
              gather_axpy(primal_vector)
    {
    }

  protected:
    void _prepare_cell(IndexType cell)
    {
      dof_mapping.prepare(cell);
      trafo_eval.prepare(cell);
      space_eval.prepare(trafo_eval);
      // gather our local vector
      local_vector.format();
      gather_axpy(local_vector, dof_mapping);
    }

    void _prepare_point(const DomainPointType& point)
    {
      trafo_eval(trafo_data, point);
      space_eval(space_data, trafo_data);
      const int num_loc_dofs = this->space_eval.get_num_local_dofs();
      if constexpr(loc_value_holder.has_value)
      {
        loc_value_holder.value = ValueType(0);
        for(int i = 0; i < num_loc_dofs; ++i)
        {
          Tiny::axpy(loc_value_holder.value, local_vector[i], space_data.phi[i].value);
        }
      }
      if constexpr(loc_value_holder.has_grad)
      {
        std::cout << "Has grad\n";
        loc_value_holder.grad.format();
        for(int i = 0; i < num_loc_dofs; ++i)
        {
          if constexpr(image_dim == 1)
            Tiny::axpy(loc_value_holder.grad, space_data.phi[i].grad, local_vector[i]);
          else
            loc_value_holder.grad.add_outer_product(local_vector[i], space_data.phi[i].grad);
        }
      }
      if constexpr(loc_value_holder.has_hess)
      {
        XABORTM("Hessian not implemented yet");
      }
    }

    void _integrate(DataType weight, IndexType point_idx)
    {
      XABORTM("Has to be implmented by derived class");
    }

    void _assemble_cell(const std::vector<IndexType>& domain_point_idx)
    {
      for(auto pti : domain_point_idx)
      {
        const auto dom_point = this->_domain_points[pti];
        this->cast()._prepare_point(dom_point);

        this->cast()._integrate(this->_point_weights[pti], pti);
      }
    }

    void _assemble()
    {
      for(IndexType k = 0; k < this->_cell_helper.size(); ++k)
      {
        const auto& domain_point_idx = this->_cell_to_domain_point.at(k);
        if(domain_point_idx.empty())
          continue;
        // prepare our cell
        this->cast()._prepare_cell(this->_cell_helper[k]);

        this->cast()._assemble_cell(domain_point_idx);
      }
    }


  public:
    void prepare_cell(IndexType cell)
    {
      this->cast()._prepare_cell(cell);
    }

    void prepare_point(const DomainPointType& point)
    {
      this->cast()._prepare_cell(point);
    }

    void assemble()
    {
      this->cast()._assemble();
    }

  };

  template<typename FaceVector_, typename Space_, typename FEVector_>
  class NormalValueSurfaceIntegratorJob
  {
  public:
    typedef Space_ SpaceType;
    typedef FEVector_ FEVector;
    typedef FaceVector_ FaceVector;
    /// the data-type of the vector
    typedef typename FEVector_::DataType DataType;
    /// the value-type of the vector
    typedef typename FEVector_::ValueType ValueType;

    static constexpr TrafoTags trafo_config_ = TrafoTags::img_point | TrafoTags::dom_point;
    static constexpr SpaceTags space_config_ = SpaceTags::value;

    /// our assembly traits
    typedef Assembly::AsmTraits1<
      DataType,
      Space_,
      trafo_config_,
      space_config_
    > AsmTraits;




    class Task :
      public FEVectorIntegratorTaskCRTP<Task, AsmTraits, Space_, FEVector_>
    {
    public:
      typedef FEVectorIntegratorTaskCRTP<Task, AsmTraits, Space_, FEVector_> BaseClass;
      typedef typename BaseClass::IndexType IndexType;
      typedef typename BaseClass::DomainPointType DomainPointType;
      static constexpr int dim = BaseClass::dim;
      /// the vector that is to be assembled
      FaceVector_& vector;
      DataType face_val;
      /// the scatter scaling factor
      DataType scatter_alpha;

    public:
      /**
       * \brief Constructor
       */
      explicit Task(NormalValueSurfaceIntegratorJob& job) :
               BaseClass(job.space, job.primal_vec),
               vector(job.face_vec),
               face_val(DataType(0)),
               scatter_alpha(job.scatter_alpha)
      {
      }

      void _assemble()
      {
        face_val = DataType(0);

        BaseClass::_assemble();

        this->vector(this->_cur_surface_index, face_val*scatter_alpha);
      }

      void _integrate(DataType weight, [[maybe_unused]] IndexType point_idx)
      {
        auto normal_val = Tiny::dot(this->loc_value_holder.value, this->_normal);
        face_val += normal_val * weight;
      }

    }; // class Task

    FaceVector_& face_vec;
    const FEVector_& primal_vec;
    const Space_& space;
    DataType scatter_alpha;

    explicit NormalValueSurfaceIntegratorJob(FaceVector_& face_vec_, const FEVector_& primal_vec_, const Space_& space_, DataType scatter_alpha_ = DataType(1)) :
      face_vec(face_vec_),
      primal_vec(primal_vec_),
      space(space_),
      scatter_alpha(scatter_alpha_)
    {
    }

  }; // class NormalValueSurfaceIntegrator

  template<typename FaceVector_, typename Space_, typename FEVector_>
  class NormalGradientSurfaceIntegratorJob
  {
  public:
    typedef Space_ SpaceType;
    typedef FEVector_ FEVector;
    typedef FaceVector_ FaceVector;
    /// the data-type of the vector
    typedef typename FEVector_::DataType DataType;
    /// the value-type of the vector
    typedef typename FEVector_::ValueType ValueType;

    static constexpr TrafoTags trafo_config_ = TrafoTags::img_point | TrafoTags::dom_point;
    static constexpr SpaceTags space_config_ = SpaceTags::grad;

    /// our assembly traits
    typedef Assembly::AsmTraits1<
      DataType,
      Space_,
      trafo_config_,
      space_config_
    > AsmTraits;




    class Task :
      public FEVectorIntegratorTaskCRTP<Task, AsmTraits, Space_, FEVector_>
    {
    public:
      typedef FEVectorIntegratorTaskCRTP<Task, AsmTraits, Space_, FEVector_> BaseClass;
      typedef typename BaseClass::IndexType IndexType;
      typedef typename BaseClass::DomainPointType DomainPointType;
      static constexpr int dim = BaseClass::dim;
      /// the vector that is to be assembled
      FaceVector_& vector;
      DataType face_val;
      /// the scatter scaling factor
      DataType scatter_alpha;

    public:
      /**
       * \brief Constructor
       */
      explicit Task(NormalGradientSurfaceIntegratorJob& job) :
               BaseClass(job.space, job.primal_vec),
               vector(job.face_vec),
               face_val(DataType(0)),
               scatter_alpha(job.scatter_alpha)
      {
      }

      void _assemble()
      {
        face_val = DataType(0);

        BaseClass::_assemble();

        this->vector(this->_cur_surface_index, face_val*scatter_alpha);
      }

      void _integrate(DataType weight, [[maybe_unused]] IndexType point_idx)
      {
        auto normal_grad = this->loc_value_holder.grad * this->_normal;
        face_val += Math::sqrt(Tiny::dot(normal_grad, normal_grad)) * weight;
      }

    }; // class Task

    FaceVector_& face_vec;
    const FEVector_& primal_vec;
    const Space_& space;
    DataType scatter_alpha;

    explicit NormalGradientSurfaceIntegratorJob(FaceVector_& face_vec_, const FEVector_& primal_vec_, const Space_& space_, DataType scatter_alpha_ = DataType(1)) :
      face_vec(face_vec_),
      primal_vec(primal_vec_),
      space(space_),
      scatter_alpha(scatter_alpha_)
    {
    }

  }; // class NormalGradientSurfaceIntegrator
}
