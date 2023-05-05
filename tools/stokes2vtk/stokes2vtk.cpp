// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2023 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// This tool reads in either a single or an entire sequence of joined binary Q2/P1dc Stokes vectors, which have been
// written out by some application (typically one of the CCND applications), and generates a (sequence of) VTU file(s)
// that contain the velocity and pressure field data on a once refined mesh, which corresponds to Q1/Q0 data on that
// refined level. This way, only the binary Stokes vectors have to be saved and backup'ed instead of the actual
// VTU files, which can save a majority of the required disk space.
//
// \author Peter Zajac
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

#include <kernel/runtime.hpp>
#include <kernel/util/simple_arg_parser.hpp>
#include <kernel/util/stop_watch.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/geometry/export_vtk.hpp>
#include <kernel/geometry/mesh_atlas.hpp>
#include <kernel/geometry/mesh_node.hpp>
#include <kernel/geometry/mesh_part.hpp>
#include <kernel/geometry/mesh_file_reader.hpp>
#include <kernel/trafo/standard/mapping.hpp>
#include <kernel/space/lagrange2/element.hpp>
#include <kernel/space/discontinuous/element.hpp>
#include <kernel/assembly/asm_traits.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/dense_vector_blocked.hpp>
#include <kernel/lafem/dense_vector.hpp>
#include <kernel/lafem/tuple_vector.hpp>

namespace Stokes2Vtk
{
  using namespace FEAT;
  using namespace FEAT::Geometry;

  template<typename DT_, typename IT_, typename Trafo_>
  void asm_p1dc_to_4xp0dc(LAFEM::SparseMatrixCSR<DT_, IT_>& matrix,
    const Space::Discontinuous::Element<Trafo_, Space::Discontinuous::Variant::StdPolyP<1>>& space)
  {
    typedef DT_ DataType;

    // create refined midpoint cubature rule
    Cubature::DynamicFactory cubature_factory("refine:midpoint");
    Cubature::Rule<typename Trafo_::ShapeType, DataType, DataType> cubature_rule(Cubature::ctor_factory, cubature_factory);

    const Index row_bs = Index(cubature_rule.get_num_points());
    const Index col_bs = Index(space.num_local_dofs);
    const Index num_elems = space.get_mesh().get_num_elements();

    // assemble matrix structure if necessary
    if(matrix.used_elements() == 0u)
    {
      const Index num_rows = row_bs * num_elems;
      matrix = LAFEM::SparseMatrixCSR<DT_, IT_>(row_bs * num_elems, col_bs * num_elems, row_bs * col_bs * num_elems);
      IT_* row_ptr = matrix.row_ptr();
      IT_* col_idx = matrix.col_ind();

      for(Index irow(0), inze(0); irow < num_rows; ++irow)
      {
        row_ptr[irow] = irow * col_bs;
        for(Index j(0); j < col_bs; ++j, ++inze)
          col_idx[inze] = (irow / row_bs) * col_bs + j;
      }
      row_ptr[num_rows] = num_rows * col_bs;
    }

    // assembly traits
    typedef Assembly::AsmTraits1<DataType, Space::Discontinuous::Element<Trafo_, Space::Discontinuous::Variant::StdPolyP<1>>, TrafoTags::none, SpaceTags::value> AsmTraits;

    // create a trafo evaluator
    typename AsmTraits::TrafoEvaluator trafo_eval(space.get_trafo());

    // create a space evaluator and evaluation data
    typename AsmTraits::SpaceEvaluator space_eval(space);

    // create trafo evaluation data
    typename AsmTraits::TrafoEvalData trafo_data;

    // create space evaluation data
    typename AsmTraits::SpaceEvalData space_data;

    const IT_* row_ptr = matrix.row_ptr();
    DT_* val = matrix.val();

    // loop over all elements
    for(Index ielem(0); ielem < num_elems; ++ielem)
    {
      // prepare trafo evaluator
      trafo_eval.prepare(ielem);

      // prepare space evaluator
      space_eval.prepare(trafo_eval);

      // fetch number of local dofs
      XASSERT(space_eval.get_num_local_dofs() == int(col_bs));

      for(int k(0); k < int(row_bs); ++k)
      {
        // compute trafo data
        trafo_eval(trafo_data, cubature_rule.get_point(k));

        // compute basis function data
        space_eval(space_data, trafo_data);

        // fill matrix row
        IT_ inze = row_ptr[(ielem*row_bs) + Index(k)];
        for(int j(0); j < int(col_bs); ++j, ++inze)
          val[inze] = space_data.phi[j].value;
      }

      // finish evaluators
      space_eval.finish();
      trafo_eval.finish();
    }
  }

  template<typename Mesh_>
  int run_xml(SimpleArgParser& args, Geometry::MeshFileReader& mesh_reader)
  {
    static constexpr int dim = Mesh_::shape_dim;
    StopWatch watch_total, watch_mesh_in, watch_vec_in, watch_vtk_out, watch_mesh_format, watch_vec_format;
    watch_total.start();

    String stokes_name, velo_name/*, pres_name*/;
    int index_first(0), index_last(-1);

    // do we read in entire Stokes vectors?
    if(args.parse("stokes", stokes_name, index_first, index_last) < 0)
    {
      std::cout << "ERROR: You have to specify the input Stokes filename sequence via '--stokes <pattern> [<first-number> [<last-number>]]'" << std::endl;
      std::cout << "If <pattern> does not contain any asterisks, then it is interpreted as the name of a single file." << std::endl;
      std::cout << "If <pattern> contains a block of one or more asterisks ('*'), then this block servres as a" << std::endl;
      std::cout << "placeholder for the sequence index, e.g. the pattern 'foo.***.bin' will expand to 'foo.000.bin'," << std::endl;
      std::cout << "'foo.001.bin.', 'foo.002.bin', etc." << std::endl;
      std::cout << "If given, <first-index> specifies the first index for the pattern, otherwise the first index is 0." << std::endl;
      std::cout << "If given, <last-index> specifies the last index for the pattern, otherwise this tool will stop at" << std::endl;
      std::cout << "the first index for which no input file was found anymore." << std::endl << std::endl;
      return 1;
    }
    // do we read in only velocity vectors?
    if(stokes_name.empty())
    {
      if(args.parse("velo", velo_name, index_first, index_last) < 1)
      {
        std::cout << "ERROR: You have to specify the input velocity filename sequence via '--velo <pattern> [<first-number> [<last-number>]]'" << std::endl;
        std::cout << "If <pattern> does not contain any asterisks, then it is interpreted as the name of a single file." << std::endl;
        std::cout << "If <pattern> contains a block of one or more asterisks ('*'), then this block servres as a" << std::endl;
        std::cout << "placeholder for the sequence index, e.g. the pattern 'foo.***.bin' will expand to 'foo.000.bin'," << std::endl;
        std::cout << "'foo.001.bin.', 'foo.002.bin', etc." << std::endl;
        std::cout << "If given, <first-index> specifies the first index for the pattern, otherwise the first index is 0." << std::endl;
        std::cout << "If given, <last-index> specifies the last index for the pattern, otherwise this tool will stop at" << std::endl;
        std::cout << "the first index for which no input file was found anymore." << std::endl << std::endl;
        return 1;
      }
      /*if(args.parse("pres", pres_name) < 0)
      {
        std::cout << "ERROR: pres" << std::endl;;
        return 1;
      }*/
    }

    //const bool have_velo = (!stokes_name.empty() || !velo_name.empty());
    const bool have_pres = (!stokes_name.empty());

    // split input filename pattern into head and tail
    String name_head, name_tail;
    std::size_t name_n(0u);
    {
      String name_pattern = (!stokes_name.empty() ? stokes_name : velo_name);
      std::size_t i0 = name_pattern.find_first_of("*");
      std::size_t i1 = name_pattern.find_last_of("*");
      if(i0 != name_pattern.npos)
      {
        name_head = name_pattern.substr(0u, i0);
        name_tail = name_pattern.substr(i1+1u);
        name_n = i1 - i0 + 1u;
      }
      else // no pattern
      {
        name_head = name_pattern;
      }
    }

    // dump some info about our input
    if(!stokes_name.empty())
      std::cout << (!stokes_name.empty() ? "Stokes" : "Velocity");
    if(name_n > 0u)
    {
      std::cout << " filename sequence: " << name_head << String(name_n, '*') << name_tail << std::endl;
      std::cout << "First sequence index: " << index_first << std::endl;
      std::cout << "Last sequence index: " << index_last << std::endl;
    }
    else
      std::cout << " filename: " << name_head << std::endl;

    String vtk_name = name_head;
    if(name_n == 0u)
    {
      // try to split off extension
      std::size_t i0 = vtk_name.find_last_of('.');
      if(i0 != vtk_name.npos)
        vtk_name = vtk_name.substr(0u, i0);
    }
    const std::size_t pad_size = Math::max(std::size_t(4u), name_n);
    args.parse("vtk", vtk_name);
    if(name_n > 0u)
      std::cout << "VTU output sequence: " << vtk_name + "." + String(pad_size, '*') + ".vtu" << std::endl;
    else
      std::cout << "VTU output filename: " << vtk_name + ".vtu" << std::endl;

    bool want_diff = (args.check("diff") >= 0);
    std::cout << "Compute difference: " << (want_diff ? "yes" : "no") << std::endl;

    // parse levels
    Index level(0);
    args.parse("level", level);

    watch_mesh_in.start();

    // create an empty atlas and a root mesh node
    auto atlas = Geometry::MeshAtlas<Mesh_>::make_unique();
    auto node = Geometry::RootMeshNode<Mesh_>::make_unique(nullptr, atlas.get());

    // try to parse the mesh file
  #ifndef DEBUG
    try
  #endif
    {
      std::cout << "Parsing mesh files..." << std::endl;
      // Now parse the mesh file
      mesh_reader.parse(*node, *atlas, nullptr);
    }
  #ifndef DEBUG
    catch(std::exception& exc)
    {
      std::cerr << "ERROR: " << exc.what() << std::endl;
      return 1;
    }
    catch(...)
    {
      std::cerr << "ERROR: unknown exception" << std::endl;
      return 1;
    }
  #endif

    // adapt coarse mesh
    node->adapt();

    // refine up to level + 1
    std::cout << "Refining up to level " << level << " (+1)..." << std::endl;
    for(Index l(0); l < level; ++l)
    {
      node = node->refine_unique(AdaptMode::chart);
    }

    // get pressure DOF count
    const Index num_pres_dofs = node->get_mesh()->get_num_elements() * Index(dim+1);

    // do we need a pressure prolongation matrix?
    LAFEM::SparseMatrixCSR<double, Index> pres_pol;
    if(have_pres)
    {
      Trafo::Standard::Mapping<Mesh_> trafo(*node->get_mesh());
      Space::Discontinuous::ElementP1<decltype(trafo)> space(trafo);
      asm_p1dc_to_4xp0dc(pres_pol, space);
    }

    // refine once more
    node = node->refine_unique(AdaptMode::chart);

    // get our mesh
    Mesh_& mesh = *node->get_mesh();

    watch_mesh_in.stop();

    // get velocity DOF count
    const Index num_velo_dofs = mesh.get_num_vertices();

    // write our mesh into a string
    String vtu_head, vtu_tail;
    watch_mesh_format.start();
    {
      Geometry::ExportVTK<Mesh_> exporter(mesh);
      std::ostringstream oss;
      exporter.write_vtu(oss);
      String str_mesh = oss.str();

      // remove the tail of the file
      std::size_t off = str_mesh.rfind("</Piece>");
      XASSERT(off < str_mesh.size());

      // build head and tail
      vtu_head = str_mesh.substr(0u, off);
      vtu_tail = str_mesh.substr(off);
    }
    watch_mesh_format.stop();

    // free mesh node
    node.reset();

    std::uint64_t total_in(0u), total_out(0u), counter(0u);

    typedef LAFEM::DenseVectorBlocked<double, Index, dim> VeloVector;
    typedef LAFEM::DenseVector<double, Index> PresVector;
    typedef LAFEM::TupleVector<VeloVector, PresVector> StokesVector;

    StokesVector vec_prev;

    std::cout << std::endl;

    // loop over our input file sequence
    for(int iseq(index_first); (index_last < 0) || (iseq <= index_last); ++iseq)
    {
      // try to open input file
      String filename = name_head;
      if(name_n > 0u)
        filename += stringify(iseq).pad_front(name_n, '0') + name_tail;
      std::ifstream ifs(filename, std::ios_base::in|std::ios_base::binary);
      if(!ifs.is_open())
      {
        if((name_n > 0u) && (index_last < 0))
        {
          // not an error if index_last == -1
          std::cout << "\nVelocity file '" << filename << "' not found; so I assume that we're finished here" << std::endl;
          break;
        }
        else
        {
          std::cout << "\nERROR: Velocity file '" << filename << "' not found!" << std::endl;
          return 1;
        }
      }

      // get filesize
      ifs.seekg(0u, std::ios_base::end);
      total_in += std::size_t(ifs.tellg());
      ifs.seekg(0u, std::ios_base::beg);

      std::cout << "Processing '" << filename << "'...";

      // try to read vector
      watch_vec_in.start();
      StokesVector vector;
      try
      {
        if(!stokes_name.empty())
          vector.read_from(LAFEM::FileMode::fm_binary, ifs);
        else
          vector.first().read_from(LAFEM::FileMode::fm_binary, ifs);
        ifs.close();
      }
      catch(std::exception& exc)
      {
        std::cout << "\n\nERROR: " << exc.what() << std::endl;
        return 1;
      }
      watch_vec_in.stop();

      const Index vec_velo_size = vector.template at<0>().size();
      const Index vec_pres_size = vector.template at<1>().size(); // may be 0 if we process only velocity

      // check vector size
      if(vec_velo_size != num_velo_dofs)
      {
        std::cout << "\n\nERROR: invalid velocity vector size: found " << vec_velo_size << " DOFs but expected " << num_velo_dofs << std::endl;
        return 1;
      }
      if((vec_pres_size > 0u) && (vec_pres_size != num_pres_dofs))
      {
        std::cout << "\n\nERROR: invalid pressure vector size: found " << vec_pres_size << " DOFs but expected " << num_pres_dofs << std::endl;
        return 1;
      }

      if(vec_prev.first().empty() && want_diff)
        vec_prev = vector.clone(LAFEM::CloneMode::Shallow);

      // stringify vector
      watch_vec_format.start();
      String buffer;
      buffer.reserve(vec_velo_size * std::size_t(30*dim) + vec_pres_size * 30u);
      buffer += "<PointData>\n<DataArray type=\"Float64\" Name=\"velocity\" NumberOfComponents=\"" + stringify(dim) +"\" Format=\"ascii\">\n";
      const Tiny::Vector<double, dim>* velo_val = vector.template at<0>().elements();
      for(Index i(0); i < num_velo_dofs; ++i)
      {
        buffer += stringify(velo_val[i][0]);
        for(int j(1); j < dim; ++j)
        {
          buffer.push_back(' ');
          buffer += stringify(velo_val[i][j]);
        }
        buffer.push_back('\n');
      }

      // write diff if desired
      if(want_diff)
      {
        buffer += "</DataArray>\n<DataArray type=\"Float64\" Name=\"velocity_diff\" NumberOfComponents=\"" + stringify(dim) +"\" Format=\"ascii\">\n";
        const Tiny::Vector<double, dim>* velo_prev = vec_prev.template at<0>().elements();
        for(Index i(0); i < num_velo_dofs; ++i)
        {
          buffer += stringify(velo_val[i][0] - velo_prev[i][0]);
          for(int j(1); j < dim; ++j)
          {
            buffer.push_back(' ');
            buffer += stringify(velo_val[i][j] - velo_prev[i][j]);
          }
          buffer.push_back('\n');
        }
      }
      buffer += "</DataArray>\n</PointData>\n";

      if(have_pres)
      {
        // prolongate and write out pressure vector
        auto vec_p = pres_pol.create_vector_l();
        pres_pol.apply(vec_p, vector.template at<1>());
        double* pres_val = vec_p.elements();
        const Index np = vec_p.size();
        buffer += "<CellData>\n<DataArray type=\"Float64\" Name=\"pressure\" Format=\"ascii\">\n";
        for(Index i(0); i < np; ++i)
          (buffer += stringify(pres_val[i])).push_back('\n');
        buffer += "</DataArray>\n</CellData>\n";
      }
      watch_vec_format.stop();

      // write to file
      watch_vtk_out.start();
      String out_name = vtk_name;
      if(name_n > 0u)
        out_name += "." + stringify(iseq).pad_front(pad_size, '0');
      out_name += ".vtu";
      std::cout << " writing '" << out_name << "'...";
      std::ofstream ofs(out_name, std::ios_base::out);
      ofs << vtu_head;
      ofs << buffer;
      ofs << vtu_tail;
      ofs.close();
      std::cout << " done!" << std::endl;
      watch_vtk_out.stop();

      total_out += vtu_head.size() + vtu_tail.size() + buffer.size();
      ++counter;

      if(want_diff)
        vec_prev = std::move(vector);

      // no indexing?
      if(name_n == 0u)
        break;
    }

    // print out summary
    std::cout << "\nProcessed " << counter << " files, read in ";
    if(total_in < 1000000ull)
      std::cout << stringify_fp_fix(double(total_in) / 1024.0, 3) << " KiB";
    else if(total_in < 1000000000ull)
      std::cout << stringify_fp_fix(double(total_in) / (1024.0*1024.0), 3) << " MiB";
    else if(total_in < 1000000000000ull)
      std::cout << stringify_fp_fix(double(total_in) / (1024.0*1024.0*1024.0), 3) << " GiB";
    else
      std::cout << stringify_fp_fix(double(total_in) / (1024.0*1024.0*1024.0*1024.0), 3) << " TiB";
    std::cout << " and wrote out ";
    if(total_out < 1000000ull)
      std::cout << stringify_fp_fix(double(total_out) / 1024.0, 3) << " KiB";
    else if(total_out < 1000000000ull)
      std::cout << stringify_fp_fix(double(total_out) / (1024.0*1024.0), 3) << " MiB";
    else if(total_out < 1000000000000ull)
      std::cout << stringify_fp_fix(double(total_out) / (1024.0*1024.0*1024.0), 3) << " GiB";
    else
      std::cout << stringify_fp_fix(double(total_out) / (1024.0*1024.0*1024.0*1024.0), 3) << " TiB";
    std::cout << ", thus wasting " << stringify_fp_fix((double(total_out)/double(total_in) - 1.0)*100.0, 2);
    std::cout << "% additional disk space -- good job!" << std::endl << std::endl;

    watch_total.stop();
    double time_total = watch_total.elapsed();
    std::cout << "Total Runtime......: " << watch_total.elapsed_string().pad_front(8)  << " seconds" << std::endl;
    std::cout << "Mesh Refine Runtime: " << watch_mesh_in.elapsed_string().pad_front(8)     << " seconds ["
      << stringify_fp_fix(100.0*watch_mesh_in.elapsed() / time_total, 2, 6) << "%]" << std::endl;
    std::cout << "Mesh Format Runtime: " << watch_mesh_format.elapsed_string().pad_front(8) << " seconds ["
      << stringify_fp_fix(100.0*watch_mesh_format.elapsed() / time_total, 2, 6) << "%]" << std::endl;
    std::cout << "Vector Read Runtime: " << watch_vec_in.elapsed_string().pad_front(8)     << " seconds ["
      << stringify_fp_fix(100.0*watch_vec_in.elapsed() / time_total, 2, 6) << "%]" << std::endl;
    std::cout << "Vec Format Runtime.: " << watch_vec_format.elapsed_string().pad_front(8) << " seconds ["
      << stringify_fp_fix(100.0*watch_vec_format.elapsed() / time_total, 2, 6) << "%]" << std::endl;
    std::cout << "VTK Write Runtime..: " << watch_vtk_out.elapsed_string().pad_front(8)    << " seconds ["
      << stringify_fp_fix(100.0*watch_vtk_out.elapsed() / time_total, 2, 6) << "%]" << std::endl;

    return 0;
  }

  int run(int argc, char* argv[])
  {
    SimpleArgParser args(argc, argv);

    args.support("mesh");
    args.support("level");
    args.support("stokes");
    args.support("velo");
    //args.support("pres");
    args.support("vtk");
    args.support("diff");

    // need help?
    if((argc < 2) || (args.check("help") > -1))
    {
      std::cout << "\nUSAGE: stokes2vtu <options...>" << std::endl;
      std::cout << "\nMandatory Options:" << std::endl;
      std::cout <<   "------------------" << std::endl;
      std::cout << "--mesh <meshfiles...>" << std::endl;
      std::cout << "Specifies the mesh files to be read in." << std::endl << std::endl;
      std::cout << "--level <level>" << std::endl;
      std::cout << "Specifies the mesh refinement level of the velocity vectors." << std::endl << std::endl;
      std::cout << "--stokes <pattern> [<first-index> [<last-index>]]" << std::endl;
      std::cout << "Specifies the filename pattern and index bounds for the input Stokes vectors." << std::endl;
      std::cout << "If <pattern> contains a block of one or more asterisks ('*'), then this block servres as a" << std::endl;
      std::cout << "placeholder for the sequence index, e.g. the pattern 'foo.***.bin' will expand to 'foo.000.bin'," << std::endl;
      std::cout << "'foo.001.bin.', 'foo.002.bin', etc." << std::endl;
      std::cout << "If <pattern> does not contain any asterisks, then it is interpreted as the name of a single file." << std::endl;
      std::cout << "If given, <first-index> specifies the first index for the pattern, otherwise the first index is 0." << std::endl;
      std::cout << "If given, <last-index> specifies the last index for the pattern, otherwise this tool will stop at" << std::endl;
      std::cout << "the first index for which no input file was found anymore." << std::endl << std::endl;
      std::cout << "--velo <pattern> [<first-index> [<last-index>]]" << std::endl;
      std::cout << "Specifies the filename pattern and index bounds for the input Stokes vectors." << std::endl;
      std::cout << "This option is mutually exclusive with --stokes and allows to read in only velocity vectors rather" << std::endl;
      std::cout << "than Stokes vectors containing both velocity and pressure." << std::endl << std::endl;
      //           "123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-123456789-
      std::cout << "Optional Options:" << std::endl;
      std::cout << "-----------------" << std::endl;
      std::cout << "--vtk <name-prefix>" << std::endl;
      std::cout << "Specifies the prefix for the VTK files; if not given, then the prefix (that's everything in front" << std::endl;
      std::cout << "of the first asterisk) of the velocity input filename pattern is used." << std::endl << std::endl;
      std::cout << "--diff" << std::endl;
      std::cout << "If specified, then the difference of the current velocity field and the previous time step is also" << std::endl;
      std::cout << "written to the VTK file as a separate variable named 'velocity_diff'." << std::endl;
      return 0;
    }

    // check for unsupported options
    auto unsupported = args.query_unsupported();
    if( !unsupported.empty() )
    {
      // print all unsupported options to cerr
      for(auto it = unsupported.begin(); it != unsupported.end(); ++it)
        std::cerr << "ERROR: unsupported option '--" << (*it).second << "'" << std::endl;

      return 1;
    }

    int num_mesh_files = args.check("mesh");
    if(num_mesh_files < 1)
    {
      std::cerr << "ERROR: You have to specify at least one meshfile with --mesh <files...>" << std::endl;
      return 1;
    }

    // get our filename deque
    auto mpars = args.query("mesh");
    XASSERT(mpars != nullptr);

    // create an empty mesh file reader
    Geometry::MeshFileReader mesh_reader;
    mesh_reader.add_mesh_files(mpars->second);

    // read root markup
    mesh_reader.read_root_markup();

    // get mesh type
    const String mtype = mesh_reader.get_meshtype_string();

    std::cout << "Mesh Type: " << mtype << std::endl;

    // This is the list of all supported meshes that could appear in the mesh file
    typedef Geometry::ConformalMesh<Shape::Hypercube<2>, 2, Real> H2M2D;
    typedef Geometry::ConformalMesh<Shape::Hypercube<3>, 3, Real> H3M3D;

    if(mtype == "conformal:hypercube:2:2")
      return run_xml<H2M2D>(args, mesh_reader);
    if(mtype == "conformal:hypercube:3:3")
      return run_xml<H3M3D>(args, mesh_reader);

    std::cout << "ERROR: unsupported mesh type!" << std::endl;

    return 1;
  }
} // namespace Stokes2Vtk

int main(int argc, char* argv[])
{
  FEAT::Runtime::initialize(argc, argv);
  int ret = Stokes2Vtk::run(argc, argv);
  FEAT::Runtime::finalize();
  return ret;
}
