// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2021 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.
#pragma once

#ifndef KERNEL_GEOMETRY_PARSED_HIT_TEST_FACTORY_HPP
#define KERNEL_GEOMETRY_PARSED_HIT_TEST_FACTORY_HPP 1

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
     * \int dim_
     * Dimension of the function formula.
     *
     * \author Gesa Pottbrock
     */
    template<typename Mesh_, int dim_>
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
          if (dim_ > 1) vars += ",y";
          if (dim_ > 2) vars += ",z";

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
          const Tiny::Vector<double, dim_> vars(point);

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
          //Returns true, if point is in the "inner" part of the Mesh
          //Otherwise returns false
          return (val>=0.0);
        }
      };
      //Member Variables of ParsedHitTestFactory
      ///class for parsed_hit_fuction
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
       * Every point, for which the formula equals >= 0 is in a "inner" point weighted with 1.
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

      void parse(const String& formula)
      {
        _hit_func.parse(formula);
        Intern::HitTestCompute<const ParsedHitFunction, Mesh_, ShapeType>::wrap(_target_data, _mesh, _hit_func);
      }

      /// \copydoc Factory::get_num_entities()
      virtual Index get_num_entities(int dim) override
      {
        return Index(_target_data.at(std::size_t(dim)).size());
      }

      virtual void fill_target_sets(TargetSetHolderType& target_set_holder) override
      {
        // call wrapper
        Intern::HitTestTargeter<ShapeType>::wrap(target_set_holder, _target_data);
        //includes functions apply, wrap
      }

      virtual void fill_attribute_sets(typename MeshType::AttributeSetContainer&) override
      {
        // do nothing as the object has no attribute sets
      }

      virtual void fill_index_sets(typename MeshType::IndexSetHolderType*&) override
      {
        // do nothing as the object has no index sets
      }
    }; // class ParsedHitTestFactory

  } // namespace Geometry
}//namespace FEAT
#endif // defined(FEAT_HAVE_FPARSER) || defined(DOXYGEN)
#endif //KERNEL_GEOMETRY_PARSED_HIT_TEST_FACTORY_HPP
