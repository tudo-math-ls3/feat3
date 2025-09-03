// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

// includes, FEAT
#include <kernel/geometry/hit_test_factory.hpp>
#include <kernel/analytic/function.hpp>
#include <kernel/util/exception.hpp>
#include <kernel/util/tiny_algebra.hpp>
#include <kernel/analytic/parsed_function.hpp>

// includes, system
#include<vector>

// The contents of this file require the 'fparser' third-party library.
#if defined(FEAT_HAVE_FPARSER) || defined(DOXYGEN)
#include <fparser.hh>

namespace FEAT
{
  namespace Geometry
  {
    /**
     * \brief Parsed-Hit-Test Factory class template
     *
     * This class template can be used to create a MeshPart for a particular mesh,
     * which consists of all entities that are inside the region characterized by a
     * function formula string.
     *
     * \tparam Mesh_
     * The type of the mesh for which the cell sub-set is to be computed.
     *
     * \author Gesa Pottbrock
     */
    template<typename Mesh_>
    class ParsedHitTestFactory :
      public Factory< MeshPart<Mesh_> >
    {
    public:
      /// The shape type of the mesh
      typedef typename Mesh_::ShapeType ShapeType;
      /// mesh part type
      typedef typename Mesh_::VertexType PointType;
      /// mesh part type
      typedef MeshPart<Mesh_> MeshType;
      /// target set holder type
      typedef typename MeshType::TargetSetHolderType TargetSetHolderType;

      /// mesh world dimension
      static constexpr int world_dim = Mesh_::world_dim;

    protected:
      class ParsedHitFunction
      {
      public:
        // mutable is mandatory! It cancels the const in ParsedHitTestFactory.
        mutable::FunctionParser _parser;

        ParsedHitFunction() :
          _parser()
        {
          _parser.AddConstant("pi", Math::pi<double>());
        }

        /**
         * \brief Creates a ParsedHitFunction
         *
         * This class creates a fparser Object from the given formula.
         * The ParsedHitFunction Object is a member of ParsedHitTestFactory.
         *
         * \param[in] formula
         * A string reference that represents the function formula.
         */
        void parse(const String& formula)
        {
          // add variables to our parser
          String vars("x");
          if (world_dim > 1) vars += ",y";
          if (world_dim > 2) vars += ",z";

          // try to parse the function
          const int ret = _parser.Parse(formula.c_str(), vars.c_str());
          if (ret >= 0)
          {
            String msg(_parser.ErrorMsg());
            msg.append("\n>>> '");
            msg.append(formula);
            msg.append("'");
            if (ret < int(formula.size()))
            {
              // ret contains the index of the first invalid character in the input string
              // append an additional line to mark the faulty character
              msg.append("\n>>>");
              msg.append(String(std::size_t(ret + 2), '-'));
              msg.append("^");
            }

            throw ParseError(msg);
          } //end : can't parse function

          // optimize the parsed function
          _parser.Optimize();
        }

        bool operator()(const PointType& point) const
        {
          //convert DataType to double
          const Tiny::Vector<double, world_dim> vars(point);

          // evaluate the parser
          const double val= _parser.Eval(vars.v);

          // check for errors
          switch (_parser.EvalError())
          {
          case 0: // no error
            break;

          case 1: // division by zero
            throw ParseError("Error in ParsedScalarFunction evaluation: division by zero");

          case 2: // sqrt error
            throw ParseError("Error in ParsedScalarFunction evaluation: sqrt of negative value");

          case 3: // log error
            throw ParseError("Error in ParsedScalarFunction evaluation:logarithm of a negative value");

          case 4: // trigonometric error
            throw ParseError("Error in ParsedScalarFunction evaluation: illegal input value");

          case 5: // recursion error
            throw ParseError("Error in ParsedScalarFunction evaluation: maximum recursion depth reached");

          default: // ???
            throw ParseError("Error in ParsedScalarFunction evaluation: unknown error");
          }

          // return true, if point is in the "inner" part of the mesh, otherwise return false
          return (val > 0.0);
        }
      };

      //Member Variables of ParsedHitTestFactory
      /// class for parsed_hit_fuction
      ParsedHitFunction _hit_func;
      /// reference to the input mesh
      const Mesh_& _mesh;
      /// internal data storing the indices
      std::vector<std::vector<Index>> _target_data;

    public:
      explicit ParsedHitTestFactory(const Mesh_& mesh) :
        _mesh(mesh),
        _target_data(std::size_t(_mesh.shape_dim + 1))
      {
      }

      /**
       * \brief Creates the ParsedHitTestFactory.
       *
       * The formula dictates the Meshpart.
       * Every point, for which the formula equals > 0 is in a "inner" point weighted with 1.
       * Otherwise the point will be weighted zero.
       *
       * \param[in] mesh
       * A \resident reference to the mesh for which the mesh part is to be computed.
       *
       * \param[in] formula for hit function
       * A string of the hit function.
       */
      explicit ParsedHitTestFactory(const Mesh_& mesh, const String& formula) :
        _hit_func(),
        _mesh(mesh),
        _target_data(std::size_t(_mesh.shape_dim + 1))
      {
        parse(formula);
      }

      /**
       * \brief Creates the ParsedHitTestFactory.
       *
       * The formula dictates the Meshpart.
       * Every point, for which the formula equals >= 0 is in a "inner" point weighted with 1.
       * Otherwise the point will be weighted zero.
       *
       * \param[in] mesh
       * A \resident reference to the mesh for which the mesh part is to be computed.
       *
       * \param[in] filter
       * A \resident reference to a filter meshpart. Only mesh entities of the filter
       * meshpart are tested. All other entities are considered as outside the region.
       *
       * \param[in] formula for hit function
       * A string of the hit function.
       */
      explicit ParsedHitTestFactory(const Mesh_& mesh, const MeshType& filter, const String& formula) :
        _hit_func(),
        _mesh(mesh),
        _target_data(std::size_t(_mesh.shape_dim + 1))
      {
        parse_filtered(formula, filter);
      }

      void parse(const String& formula)
      {
        _hit_func.parse(formula);
        Intern::hittest_compute_target_data<ParsedHitFunction, Mesh_>(_target_data, _mesh, _hit_func);
      }

      void parse_filtered(const String& formula, const MeshType& filter)
      {
        _hit_func.parse(formula);
        Intern::hittest_compute_filtered_target_data<ParsedHitFunction, Mesh_>(_target_data, _mesh, filter, _hit_func);
      }

      /// \copydoc Factory::get_num_entities()
      virtual Index get_num_entities(int dim) override
      {
        return Index(_target_data.at(std::size_t(dim)).size());
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        Intern::write_to_target_set<ShapeType>(target_set_holder, _target_data);
      }

      virtual void fill_attribute_sets(typename MeshType::AttributeSetContainer&) override
      {
        // do nothing as the object has no attribute sets
      }

      virtual void fill_index_sets(std::unique_ptr<typename MeshType::IndexSetHolderType>&) override
      {
        // do nothing as the object has no index sets
      }
    }; // class ParsedHitTestFactory

    /**
     * \brief Creates a new mesh-part from a formula hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] formula
     * A string of the hit function
     *
     * \returns
     * A mesh-part containing all entities for which the formula is >= 0
     */
    template<typename Mesh_>
    MeshPart<Mesh_> make_meshpart_by_formula_hit_test(const Mesh_& mesh, const String& formula)
    {
      ParsedHitTestFactory<Mesh_> factory(mesh, formula);
      return factory.make();
    }

    /**
     * \brief Creates a new mesh-part from a formula hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] filter
     * A \resident reference to a filter meshpart. Only mesh entities of the filter
     * meshpart are tested. All other entities are considered as outside the region.
     *
     * \param[in] formula
     * A string of the hit function
     *
     * \returns
     * A mesh-part containing all entities of filter for which the formula is >= 0
     */
    template<typename Mesh_>
    MeshPart<Mesh_> make_meshpart_by_filtered_formula_hit_test(const Mesh_& mesh, const MeshPart<Mesh_>& filter, const String& formula)
    {
      ParsedHitTestFactory<Mesh_> factory(mesh, filter, formula);
      return factory.make();
    }

    /**
     * \brief Creates a new mesh-part from a formula hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] formula
     * A string of the hit function
     *
     * \returns
     * A unique-pointer to mesh-part containing all entities for which the formula is >= 0
     */
    template<typename Mesh_>
    std::unique_ptr<MeshPart<Mesh_>> make_unique_meshpart_by_formula_hit_test(const Mesh_& mesh, const String& formula)
    {
      ParsedHitTestFactory<Mesh_> factory(mesh, formula);
      return factory.make_unique();
    }

    /**
     * \brief Creates a new mesh-part from a formula hit-test function
     *
     * \param[in] mesh
     * A \transient reference to the mesh for which a mesh-part is to be created.
     *
     * \param[in] filter
     * A \resident reference to a filter meshpart. Only mesh entities of the filter
     * meshpart are tested. All other entities are considered as outside the region.
     *
     * \param[in] formula
     * A string of the hit function
     *
     * \returns
     * A unique-pointer to a mesh-part containing all entities of filter for which the formula is >= 0
     */
    template<typename Mesh_>
    std::unique_ptr<MeshPart<Mesh_>> make_unique_meshpart_by_filtered_formula_hit_test(const Mesh_& mesh, const MeshPart<Mesh_>& filter, const String& formula)
    {
      ParsedHitTestFactory<Mesh_> factory(mesh, filter, formula);
      return factory.make_unique();
    }
  } // namespace Geometry
} //namespace FEAT
#endif // defined(FEAT_HAVE_FPARSER) || defined(DOXYGEN)
