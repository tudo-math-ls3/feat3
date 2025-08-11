// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include "vtk_writer.hpp"

#include <kernel/assembly/fe_interpolator.hpp>
#include <kernel/assembly/discrete_projector.hpp>

namespace CCNDSimple
{
  VtkWriter::VtkWriter(const DomainControl& domain_,bool refined_) :
    domain(domain_),
    want_refined(refined_)
  {
    if(refined_)
    {
      Geometry::StandardRefinery<MeshType> refinery(domain.front()->get_mesh());
      refined_mesh = refinery.make_unique();
      exporter.reset(new Geometry::ExportVTK<MeshType>(*refined_mesh));
    }
    else
    {
      exporter.reset(new Geometry::ExportVTK<MeshType>(domain.front()->get_mesh()));
    }
  }

  VtkWriter::~VtkWriter()
  {
  }

  void VtkWriter::add_supported_args(SimpleArgParser& args)
  {
    args.support("vtk", "[<filename>]\n"
      "Specifies that a VTK output file is to be written.\n"
      "If no filename is given, a default one will be used.");
    args.support("vtk-step", "<stepping>\n"
      "Specifies the VTK output stepping for unsteady simulations.");
  }

  /// parse arguments from command line
  bool VtkWriter::parse_args(SimpleArgParser& args)
  {
    if(args.parse("vtk", name_prefix) >= 0)
      stepping = 1;
    args.parse("vtk-step", stepping);
    return true;
  }

  void VtkWriter::print_config()
  {
    print_pad(domain.comm(), "VTK Name Prefix", name_prefix);
    print_pad(domain.comm(), "VTK Stepping", stringify(stepping));
  }

  /// adds a vector reference to the checkpoint system
  bool VtkWriter::register_stokes_vector(const String& name, GlobalStokesVector& vector)
  {
    return stokes_vectors.emplace(name, &vector).second;
  }

  void VtkWriter::unregister_stokes_vectors()
  {
    stokes_vectors.clear();
  }

  bool VtkWriter::prepare_write()
  {
    // if stepping is 0, we'll never write anything
    if(stepping == Index(0))
      return false;

    // do we have to create an exporter?
    if(!exporter)
    {
      // refine mesh?
      if(want_refined)
      {
        Geometry::StandardRefinery<MeshType> refinery(domain.front()->get_mesh());
        refined_mesh = refinery.make_unique();
        exporter.reset(new Geometry::ExportVTK<MeshType>(*refined_mesh));
      }
      else
      {
        exporter.reset(new Geometry::ExportVTK<MeshType>(domain.front()->get_mesh()));
      }
    }

    // set VTK filename
    vtk_name = name_prefix;

    // yep, let's write some VTKs
    return true;
  }

  bool VtkWriter::prepare_write(Index step)
  {
    plot_line.clear();

    // VTK output, anyone?
    if(!prepare_write())
      return false;

    // no write this time?
    if((step % stepping) != Index(0))
      return false;

    // set VTK filename
    vtk_name = name_prefix + "." +  stringify(step).pad_front(5, '0');

    return true;
  }

  void VtkWriter::add_stokes_vector(const GlobalStokesVector& vector, const String& v_name, const String& p_name)
  {
    add_lagrange2_vector(vector.local().at<0>(), v_name);
    add_p1dc_vector(vector.local().at<1>(), p_name);
  }

  void VtkWriter::add_lagrange1_vector(const LocalScalarVectorType& vector, const String& name)
  {
    if(refined_mesh)
    {
      // interpolate form Lagrange-1 to Lagrange-2
      LocalScalarVectorType vector2(domain.front()->space_velo.get_num_dofs());
      Assembly::FEInterpolator<SpaceVeloType, SpaceTypeQ1>::interpolate(
          vector2, vector, domain.front()->space_velo, domain.front()->space_q1);
      exporter->add_vertex_scalar(name, vector2.elements());
    }
    else
    {
      exporter->add_vertex_scalar(name, vector.elements());
    }
  }

  void VtkWriter::add_lagrange2_vector(const LocalScalarVectorType& vector, const String& name)
  {
    // lagrange-2 dofs are lagrange-1 dofs on refined mesh, so we don't need to check whether we refined here
    exporter->add_vertex_scalar(name, vector.elements());
  }

  void VtkWriter::add_lagrange1_vector(const LocalFieldVectorType& vector, const String& name)
  {
    if(refined_mesh)
    {
      // interpolate form Lagrange-1 to Lagrange-2
      LocalFieldVectorType vector2(domain.front()->space_velo.get_num_dofs());
      Assembly::FEInterpolator<SpaceVeloType, SpaceTypeQ1>::interpolate(
        vector2, vector, domain.front()->space_velo, domain.front()->space_q1);
      //exporter->add_vertex_scalar(name, vector2.elements());
      exporter->add_vertex_vector(name, vector2);
    }
    else
    {
      exporter->add_vertex_vector(name, vector);
    }
  }

  void VtkWriter::add_lagrange2_vector(const LocalFieldVectorType& vector, const String& name)
  {
    // lagrange-2 dofs are lagrange-1 dofs on refined mesh, so we don't need to check whether we refined here
    exporter->add_vertex_vector(name, vector);
  }

  void VtkWriter::add_p0dc_vector(const LocalScalarVectorType& vector, const String& name)
  {
    if(refined_mesh)
    {
      const Index n = vector.size();
      const Index bs = Index(1 << dim); // = 2^dim

      // interpolate vector; each value is spread to a constant block of the same constant value
      LocalScalarVectorType vector2(n*bs);
      for(Index i(0), j(0); i < n; ++i)
      {
        for(Index k(0); k < bs; ++j, ++k)
          vector2(j, vector(i));
      }
      exporter->add_cell_scalar(name, vector2.elements());
    }
    else
    {
      exporter->add_cell_scalar(name, vector.elements());
    }
  }

  void VtkWriter::add_p1dc_vector(const LocalScalarVectorType& vector, const String& name)
  {
    if(refined_mesh)
    {
      const Index num_elems = domain.front()->get_mesh().get_num_elements();
      XASSERT(vector.size() == num_elems*(dim+1));

      Cubature::DynamicFactory cubature_factory("refine:midpoint");
      Cubature::Rule<ShapeType, DataType, DataType> cubature_rule(Cubature::ctor_factory, cubature_factory);
      const int num_pts = cubature_rule.get_num_points();

      LocalScalarVectorType vector2(num_elems*Index(num_pts));

      const SpacePresType& space = domain.front()->space_pres;

      // assembly traits
      typedef Assembly::AsmTraits1<DataType, SpacePresType, TrafoTags::none, SpaceTags::value> AsmTraits;

      // fetch the trafo
      const typename AsmTraits::TrafoType& trafo = space.get_trafo();

      // create a trafo evaluator
      typename AsmTraits::TrafoEvaluator trafo_eval(trafo);

      // create a space evaluator and evaluation data
      typename AsmTraits::SpaceEvaluator space_eval(space);

      // create a dof-mapping
      typename AsmTraits::DofMapping dof_mapping(space);

      // create trafo evaluation data
      typename AsmTraits::TrafoEvalData trafo_data;

      // create space evaluation data
      typename AsmTraits::SpaceEvalData space_data;

      // create local vector data
      typename AsmTraits::template TLocalVector<DataType> local_vector;

      // create matrix scatter-axpy
      typename LocalScalarVectorType::GatherAxpy gather_axpy(vector);

      // loop over all elements
      for(Index ielem(0), i(0); ielem < num_elems; ++ielem)
      {
        // format local vector
        local_vector.format();

        // initialize dof-mapping
        dof_mapping.prepare(ielem);

        // gather local vector data
        gather_axpy(local_vector, dof_mapping);

        // finish dof-mapping
        dof_mapping.finish();

        // prepare trafo evaluator
        trafo_eval.prepare(ielem);

        // prepare space evaluator
        space_eval.prepare(trafo_eval);

        // fetch number of local dofs
        const int num_loc_dofs = space_eval.get_num_local_dofs();

        for(int k(0); k < num_pts; ++k, ++i)
        {
          // compute trafo data
          trafo_eval(trafo_data, cubature_rule.get_point(k));

          // compute basis function data
          space_eval(space_data, trafo_data);

          // compute function value
          DataType value = DataType(0);
          for(int j(0); j < num_loc_dofs; ++j)
            value += local_vector[j] * space_data.phi[j].value;

          // save value
          vector2(i, value);
        }

        // finish evaluators
        space_eval.finish();
        trafo_eval.finish();
      }

      exporter->add_cell_scalar(name, vector2.elements());
    }
    else
    {
      exporter->add_cell_scalar(name, vector.elements());
    }
  }

  void VtkWriter::write()
  {
    plot_line += " ; > '" + vtk_name + ".pvtu'";
    exporter->write(vtk_name, domain.comm());
    exporter->clear();
  }

  bool VtkWriter::write_registered(const TimeStepping& time_stepping)
  {
    // write VTK file if desired
    if(prepare_write(time_stepping.time_step))
    {
      if(time_stepping.full_plot)
        domain.comm().print("\nWriting VTK output to '" + vtk_name + ".pvtu'");
      for(const auto& x : stokes_vectors)
        add_stokes_vector(*x.second, x.first + "_v", x.first + "_p");
      write();
      return true;
    }
    return false;
  }

  void VtkWriter::print_runtime(double total_time)
  {
    print_time(domain.comm(), "VTK Export Write Time", watch_vtk_write.elapsed(), total_time);
  }
} // namespace CCNDSimple
