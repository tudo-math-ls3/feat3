// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/base_header.hpp>
#include <kernel/runtime.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/adjacency/cuthill_mckee.hpp>
#include <iostream>

using namespace FEAT;
using namespace FEAT::LAFEM;
using namespace FEAT::Adjacency;

int main(int argc, char ** argv)
{
  FEAT::Runtime::ScopeGuard runtime_scope_guard(argc, argv);
  if (argc != 3 && argc != 2)
  {
    std::cout<<"Usage: 'resort-csr csr-file [csr-file-resorted]'\n";
    FEAT::Runtime::abort();
  }

  String input(argv[1]);
  if (input.size() < 5)
  {
    std::cout<<"Input Filetype not known: " << input << "\n";
    FEAT::Runtime::abort();
  }
  String input_extension(input.substr(input.size() - 4, 4));
  SparseMatrixCSR<double, Index> orig;
  if (input_extension == ".mtx")
    orig.read_from(FileMode::fm_mtx, input);
  else if (input_extension == ".csr")
    orig.read_from(FileMode::fm_csr, input);
  else
  {
    std::cout<<"Input Filetype not known: " << input << "\n";
    FEAT::Runtime::abort();
  }

  Graph graph(Adjacency::RenderType::as_is, orig);
  Index best_radius;
  Index best_radius_index;
  orig.radius_row(best_radius, best_radius_index);
  std::cout<<"Initial: "<<best_radius<<"\n";

  SparseMatrixCSR<double, Index> best;
  best.clone(orig, CloneMode::Deep);

  Permutation perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::RootType::minimum_degree, Adjacency::CuthillMcKee::SortType::standard);
  auto csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  Index test_radius;
  Index test_radius_index;
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"true RootType::minimum_degree SortType::standard: "<<best_radius<<"\n";
  }
  perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::RootType::minimum_degree, Adjacency::CuthillMcKee::SortType::standard);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"false RootType::minimum_degree SortType::standard: "<<best_radius<<"\n";
  }

  perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::RootType::maximum_degree, Adjacency::CuthillMcKee::SortType::standard);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"true RootType::maximum_degree SortType::standard: "<<best_radius<<"\n";
  }
  perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::RootType::maximum_degree, Adjacency::CuthillMcKee::SortType::standard);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"false RootType::maximum_degree SortType::standard: "<<best_radius<<"\n";
  }

  perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::RootType::standard, Adjacency::CuthillMcKee::SortType::standard);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"true RootType::default SortType::standard: "<<best_radius<<"\n";
  }
  perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::RootType::standard, Adjacency::CuthillMcKee::SortType::standard);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"false RootType::default SortType::standard: "<<best_radius<<"\n";
  }


  perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::RootType::minimum_degree, Adjacency::CuthillMcKee::SortType::asc);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"true RootType::minimum_degree SortType::asc: "<<best_radius<<"\n";
  }
  perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::RootType::minimum_degree, Adjacency::CuthillMcKee::SortType::asc);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"false RootType::minimum_degree SortType::asc: "<<best_radius<<"\n";
  }

  perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::RootType::maximum_degree, Adjacency::CuthillMcKee::SortType::asc);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"true RootType::maximum_degree SortType::asc: "<<best_radius<<"\n";
  }
  perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::RootType::maximum_degree, Adjacency::CuthillMcKee::SortType::asc);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"false RootType::maximum_degree SortType::asc: "<<best_radius<<"\n";
  }

  perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::RootType::standard, Adjacency::CuthillMcKee::SortType::asc);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"true RootType::default SortType::asc: "<<best_radius<<"\n";
  }
  perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::RootType::standard, Adjacency::CuthillMcKee::SortType::asc);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"false RootType::standard SortType::asc: "<<best_radius<<"\n";
  }


  perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::RootType::minimum_degree, Adjacency::CuthillMcKee::SortType::desc);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"true RootType::minimum_degree SortType::desc: "<<best_radius<<"\n";
  }
  perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::RootType::minimum_degree, Adjacency::CuthillMcKee::SortType::desc);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"false RootType::minimum_degree SortType::desc: "<<best_radius<<"\n";
  }

  perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::RootType::maximum_degree, Adjacency::CuthillMcKee::SortType::desc);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"true RootType::maximum_degree SortType::desc: "<<best_radius<<"\n";
  }
  perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::RootType::maximum_degree, Adjacency::CuthillMcKee::SortType::desc);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"false RootType::maximum_degree SortType::desc: "<<best_radius<<"\n";
  }

  perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::RootType::standard, Adjacency::CuthillMcKee::SortType::desc);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"true RootType::standard SortType::desc: "<<best_radius<<"\n";
  }
  perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::RootType::standard, Adjacency::CuthillMcKee::SortType::desc);
  csr = orig.clone(CloneMode::Deep);
  csr.permute(perm, perm);
  csr.radius_row(test_radius, test_radius_index);
  if (test_radius < best_radius)
  {
    best_radius = test_radius;
    best.clone(csr, CloneMode::Deep);
    std::cout<<"false RootType::standard SortType::desc: "<<best_radius<<"\n";
  }

  if (argc != 2)
  {
    String output(argv[2]);
    if (output.size() < 5)
    {
      std::cout<<"Output Filetype not known: " << output << "\n";
      FEAT::Runtime::abort();
    }
    String output_extension(output.substr(output.size() - 4, 4));
    if (output_extension == ".mtx")
      best.write_out(FileMode::fm_mtx, output);
    else if (output_extension == ".csr")
      best.write_out(FileMode::fm_csr, output);
    else
    {
      std::cout<<"Output Filetype not known: " << output << "\n";
      FEAT::Runtime::abort();
    }
  }
}
