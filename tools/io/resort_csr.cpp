#include <kernel/base_header.hpp>
#include <kernel/lafem/sparse_matrix_csr.hpp>
#include <kernel/adjacency/cuthill_mckee.hpp>
#include <iostream>

using namespace FEAST;
using namespace FEAST::LAFEM;
using namespace FEAST::Adjacency;

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        std::cout<<"Usage 'resort-csr csr-file csr-file-resorted'"<<std::endl;
        exit(EXIT_FAILURE);
    }

    String input(argv[1]);
    String output(argv[2]);

    SparseMatrixCSR<Mem::Main, double, Index> csr(FileMode::fm_csr, input);
    Graph graph(Adjacency::rt_as_is, csr);
    Index best_radius;
    Index best_radius_index;
    csr.radius_row(best_radius, best_radius_index);
    std::cout<<best_radius<<std::endl;

    SparseMatrixCSR<Mem::Main, double, Index> best;
    best.clone(csr, CloneMode::Deep);

    Permutation perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_default);
    csr.permute(perm, perm);
    Index test_radius;
    Index test_radius_index;
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_default);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }

    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_maximum_degree, Adjacency::CuthillMcKee::sort_default);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_maximum_degree, Adjacency::CuthillMcKee::sort_default);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }

    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_default, Adjacency::CuthillMcKee::sort_default);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_default, Adjacency::CuthillMcKee::sort_default);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }


    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_asc);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_asc);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }

    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_maximum_degree, Adjacency::CuthillMcKee::sort_asc);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_maximum_degree, Adjacency::CuthillMcKee::sort_asc);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }

    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_default, Adjacency::CuthillMcKee::sort_asc);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_default, Adjacency::CuthillMcKee::sort_asc);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }


    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_desc);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_minimum_degree, Adjacency::CuthillMcKee::sort_desc);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }

    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_maximum_degree, Adjacency::CuthillMcKee::sort_desc);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_maximum_degree, Adjacency::CuthillMcKee::sort_desc);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }

    perm = CuthillMcKee::compute(graph, true, Adjacency::CuthillMcKee::root_default, Adjacency::CuthillMcKee::sort_desc);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }
    perm = CuthillMcKee::compute(graph, false, Adjacency::CuthillMcKee::root_default, Adjacency::CuthillMcKee::sort_desc);
    csr.permute(perm, perm);
    csr.radius_row(test_radius, test_radius_index);
    if (test_radius < best_radius)
    {
      best_radius = test_radius;
      best.clone(csr, CloneMode::Deep);
      std::cout<<best_radius<<std::endl;
    }

    best.write_out(FileMode::fm_csr, output);
}
