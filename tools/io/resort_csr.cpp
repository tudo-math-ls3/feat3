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

    SparseMatrixCSR<Mem::Main, double, Index> orig(FileMode::fm_csr, input);
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

    best.write_out(FileMode::fm_csr, output);
}
