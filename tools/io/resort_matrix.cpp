#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/lafem/sparse_matrix_ell.hpp>
#include <kernel/adjacency/cuthill_mckee.hpp>
#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Adjacency;

int main(int argc, char ** argv)
{
    if (argc != 3 && argc != 2)
    {
      std::cout<<"Usage: 'resort-csr csr-file [csr-file-resorted]'"<<std::endl;
      exit(EXIT_FAILURE);
    }

    String input(argv[1]);
    if (input.size() < 5)
    {
      std::cout<<"Input Filetype not known: " << input << std::endl;
      exit(EXIT_FAILURE);
    }
    String input_extension(input.substr(input.size() - 4, 4));
    SparseMatrixCSR<Mem::Main, double, Index> orig;
    if (input_extension == ".mtx")
      orig.read_from_mtx(input);
    else if (input_extension == ".csr")
      orig.read_from_csr(input);
    else if (input_extension == ".ell")
    {
      SparseMatrixELL<Mem::Main, double, Index> temp(FileMode::fm_ell, input);
      orig.convert(temp);
    }
    else
    {
      std::cout<<"Input Filetype not known: " << input << std::endl;
      exit(EXIT_FAILURE);
    }

    Graph graph(Adjacency::rt_as_is, orig);
    Index best_radius;
    Index best_radius_index;
    orig.radius_row(best_radius, best_radius_index);
    std::cout<<"Initial: "<<best_radius<<std::endl;

    SparseMatrixCSR<Mem::Main, double, Index> best;
    best.clone(orig, CloneMode::Deep);

    auto csr = orig.clone(CloneMode::Deep);
    Permutation perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_default);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    Index test_radius;
    Index test_radius_index;
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"true root_minimum_degree sort_default: "<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_default);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"false root_minimum_degree sort_default: "<<best_radius<<std::endl;
    }

    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_maximum_degree, Adjacency::CuthillMcKee::sort_default);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"true root_maximum_degree sort_default: "<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_maximum_degree, Adjacency::CuthillMcKee::sort_default);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"false root_maximum_degree sort_default: "<<best_radius<<std::endl;
    }

    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_default, Adjacency::CuthillMcKee::sort_default);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"true root_default sort_default: "<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_default, Adjacency::CuthillMcKee::sort_default);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"false root_default sort_default: "<<best_radius<<std::endl;
    }


    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_asc);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"true root_minimum_degree sort_asc: "<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_asc);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"false root_minimum_degree sort_asc: "<<best_radius<<std::endl;
    }

    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_maximum_degree, Adjacency::CuthillMcKee::sort_asc);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"true root_maximum_degree sort_asc: "<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_maximum_degree, Adjacency::CuthillMcKee::sort_asc);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"false root_maximum_degree sort_asc: "<<best_radius<<std::endl;
    }

    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_default, Adjacency::CuthillMcKee::sort_asc);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"true root_default sort_asc: "<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_default, Adjacency::CuthillMcKee::sort_asc);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"false root_default sort_asc: "<<best_radius<<std::endl;
    }


    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_desc);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"true root_minimum_degree sort_desc: "<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_desc);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"false root_minimum_degree sort_desc: "<<best_radius<<std::endl;
    }

    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_maximum_degree, Adjacency::CuthillMcKee::sort_desc);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"true root_maximum_degree sort_desc: "<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_maximum_degree, Adjacency::CuthillMcKee::sort_desc);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"false root_maximum_degree sort_desc: "<<best_radius<<std::endl;
    }

    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_default, Adjacency::CuthillMcKee::sort_desc);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"true root_default sort_desc: "<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_default, Adjacency::CuthillMcKee::sort_desc);
    csr = orig.clone(CloneMode::Deep);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<"false root_default sort_desc: "<<best_radius<<std::endl;
    }

    if (argc != 2)
    {
      String output(argv[2]);
      if (output.size() < 5)
      {
        std::cout<<"Output Filetype not known: " << output << std::endl;
        exit(EXIT_FAILURE);
      }
      String output_extension(output.substr(output.size() - 4, 4));
      if (output_extension == ".mtx")
        best.write_out_mtx(output);
      else if (output_extension == ".csr")
        best.write_out_csr(output);
      else if (output_extension == ".ell")
      {
        SparseMatrixELL<Mem::Main, double, Index> temp;
        temp.convert(best);
        temp.write_out_ell(output);
      }
      else
      {
        std::cout<<"Output Filetype not known: " << output << std::endl;
        exit(EXIT_FAILURE);
      }
    }

}
