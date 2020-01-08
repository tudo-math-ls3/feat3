// FEAT3: Finite Element Analysis Toolbox, Version 3
// Copyright (C) 2010 - 2020 by Stefan Turek & the FEAT group
// FEAT3 is released under the GNU General Public License version 3,
// see the file 'copyright.txt' in the top level directory for details.

#include <kernel/util/runtime.hpp>
#include <kernel/util/time_stamp.hpp>
#include <kernel/geometry/common_factories.hpp>
#include <kernel/adjacency/graph.hpp>

using namespace FEAT;

typedef Geometry::ConformalMesh<Shape::Quadrilateral> MeshType;

template<
  typename Adjactor1_,
  typename Adjactor2_>
Index _aux_inj(const Adjactor1_& adj1, const Adjactor2_& adj2, Index i, Index idx[])
{
  Index num_idx = 0;
  typename Adjactor1_::ImageIterator cur1(adj1.image_begin(i));
  typename Adjactor1_::ImageIterator end1(adj1.image_end(i));
  for(; cur1 != end1; ++cur1)
  {
    typename Adjactor2_::ImageIterator cur2(adj2.image_begin(*cur1));
    typename Adjactor2_::ImageIterator end2(adj2.image_end(*cur1));

    // loop over all image node indices
    for(; cur2 != end2; ++cur2)
    {
      Index jdx = *cur2;

      // check if we already have that image node index
      bool found = false;
      for(Index k(0); k < num_idx; ++k)
      {
        if(idx[k] == jdx)
        {
          found = true;
          break;
        }
      }
      // add the index if we don't have it already
      if(!found)
      {
        idx[num_idx] = jdx;
        ++num_idx;
      }
    }
  }
  return num_idx;
}

template<
  typename Adjactor1_,
  typename Adjactor2_>
Adjacency::Graph render_injectify_old(
  const Adjactor1_& adj1,
  const Adjactor2_& adj2)
{
  // validate adjactor dimensions
  XASSERTM(adj1.get_num_nodes_image() == adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

  // get counts
  Index _num_nodes_domain = adj1.get_num_nodes_domain();
  Index _num_nodes_image = adj2.get_num_nodes_image();
  Index _num_indices_image = 0;

  // allocate pointer array
  Index* _domain_ptr = new Index[_num_nodes_domain + 1];

  // allocate auxiliary array
  Index* aux = new Index[_num_nodes_image];

  // count number of adjacencies and build pointer array
  for(Index i(0); i < _num_nodes_domain; ++i)
  {
    _domain_ptr[i] = _num_indices_image;
    _num_indices_image += _aux_inj(adj1, adj2, i, aux);
  }
  _domain_ptr[_num_nodes_domain] = _num_indices_image;

  // delete auxiliary array
  delete [] aux;

  // allocate and build index array
  XASSERT(_num_indices_image > Index(0));
  /// \compilerhack gcc (up to version 9 at least) thinks, that _num_indices_image will be zero for sure at runtime.
#if defined(FEAT_COMPILER_GNU) && (FEAT_COMPILER_GNU < 100000)
  Index* _image_idx(nullptr);
  if (_num_indices_image > 0)
    _image_idx = new Index[_num_indices_image];
#else
  Index* _image_idx = new Index[_num_indices_image];
#endif

  for(Index i(0); i < _num_nodes_domain; ++i)
  {
    _aux_inj(adj1, adj2, i, &_image_idx[_domain_ptr[i]]);
  }

  return Adjacency::Graph(_num_nodes_domain, _num_nodes_image, _num_indices_image, _domain_ptr, _image_idx, false);
}

template<
  typename Adjactor1_,
  typename Adjactor2_>
Adjacency::Graph render_injectify_mask(
  const Adjactor1_& adj1,
  const Adjactor2_& adj2)
{
  // validate adjactor dimensions
  XASSERTM(adj1.get_num_nodes_image() == adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

  // get counts
  Index _num_nodes_domain = adj1.get_num_nodes_domain();
  Index _num_nodes_image = adj2.get_num_nodes_image();
  Index _num_indices_image = 0;

  // allocate pointer array
  Index* _domain_ptr = new Index[_num_nodes_domain + 1];

  // allocate auxiliary array
  std::vector<int> vidx_mask(_num_nodes_image, 0);
  int* idx_mask = vidx_mask.data();

  // count number of adjacencies and build pointer array
  for(Index i(0); i < _num_nodes_domain; ++i)
  {
    _domain_ptr[i] = _num_indices_image;
    for(auto it = adj1.image_begin(i); it != adj1.image_end(i); ++it)
    {
      for(auto jt = adj2.image_begin(*it); jt != adj2.image_end(*it); ++jt)
      {
        if(idx_mask[*jt] == 0)
        {
          ++_num_indices_image;
          idx_mask[*jt] = 1;
        }
      }
    }
    // reset mask
    for(auto it = adj1.image_begin(i); it != adj1.image_end(i); ++it)
      for(auto jt = adj2.image_begin(*it); jt != adj2.image_end(*it); ++jt)
        idx_mask[*jt] = 0;
  }
  _domain_ptr[_num_nodes_domain] = _num_indices_image;

  // allocate and build index array
  Index* _image_idx = new Index[_num_indices_image];
  Index cur = 0u;
  for(Index i(0); i < _num_nodes_domain; ++i)
  {
    for(auto it = adj1.image_begin(i); it != adj1.image_end(i); ++it)
    {
      for(auto jt = adj2.image_begin(*it); jt != adj2.image_end(*it); ++jt)
      {
        if(idx_mask[*jt] == 0)
        {
          _image_idx[cur] = *jt;
          ++cur;
          idx_mask[*jt] = 1;
        }
      }
    }
    // reset mask
    for(auto it = adj1.image_begin(i); it != adj1.image_end(i); ++it)
      for(auto jt = adj2.image_begin(*it); jt != adj2.image_end(*it); ++jt)
        idx_mask[*jt] = 0;
    XASSERT(cur == _domain_ptr[i+1]);
  }

  return Adjacency::Graph(_num_nodes_domain, _num_nodes_image, _num_indices_image, _domain_ptr, _image_idx, false);
}

/// renders injectified adjactor composition
template<
  typename Adjactor1_,
  typename Adjactor2_>
Adjacency::Graph render_injectify_new(
  const Adjactor1_& adj1,
  const Adjactor2_& adj2)
{
  // validate adjactor dimensions
  XASSERTM(adj1.get_num_nodes_image() == adj2.get_num_nodes_domain(), "Adjactor dimension mismatch!");

  // get counts
  Index _num_nodes_domain = adj1.get_num_nodes_domain();
  Index _num_nodes_image = adj2.get_num_nodes_image();
  Index _num_indices_image = 0;

  // allocate pointer array
  Index* _domain_ptr = new Index[_num_nodes_domain + 1];

  // allocate auxiliary index set
  std::set<Index> idx_set;

  // count number of adjacencies and build pointer array
  for(Index i(0); i < _num_nodes_domain; ++i)
  {
    _domain_ptr[i] = _num_indices_image;
    idx_set.clear();
    for(auto it = adj1.image_begin(i); it != adj1.image_end(i); ++it)
      idx_set.insert(adj2.image_begin(*it), adj2.image_end(*it));
    _num_indices_image += Index(idx_set.size());
  }
  _domain_ptr[_num_nodes_domain] = _num_indices_image;

  // allocate and build index array
  Index* _image_idx = new Index[_num_indices_image];
  for(Index i(0); i < _num_nodes_domain; ++i)
  {
    idx_set.clear();
    for(auto it = adj1.image_begin(i); it != adj1.image_end(i); ++it)
      idx_set.insert(adj2.image_begin(*it), adj2.image_end(*it));
    Index k = _domain_ptr[i];
    for(auto it = idx_set.begin(); it != idx_set.end(); ++it, ++k)
      _image_idx[k] = *it;
  }

  return Adjacency::Graph(_num_nodes_domain, _num_nodes_image, _num_indices_image, _domain_ptr, _image_idx, false);
}

void sort_new(Adjacency::Graph& graph)
{
  const Index* dom_ptr = graph.get_domain_ptr();
  Index* img_idx = graph.get_image_idx();

  for(Index i(0); i < graph.get_num_nodes_domain(); ++i)
  {
    std::sort(&img_idx[dom_ptr[i]], &img_idx[dom_ptr[i+1]]);
  }
}

void sort_old(Adjacency::Graph& graph)
{
  const Index* dom_ptr = graph.get_domain_ptr();
  Index* img_idx = graph.get_image_idx();

  // loop over all domain nodes
  for(Index i(0); i < graph.get_num_nodes_domain(); ++i)
  {
    // apply linear insertion sort onto the adjacency indices
    for(Index j(dom_ptr[i] + 1); j < dom_ptr[i+1]; ++j)
    {
      Index x = img_idx[j];
      Index k(j);
      for(; (k > dom_ptr[i]) && (x < img_idx[k-1]); --k)
      {
        img_idx[k] = img_idx[k-1];
      }
      img_idx[k] = x;
    }
  }
}


Adjacency::Graph build_prol(const MeshType& mesh_c)
{
  Index nc_v = mesh_c.get_num_entities(0);
  Index nc_e = mesh_c.get_num_entities(1);
  Index nc_q = mesh_c.get_num_entities(2);

  const auto& v_e = mesh_c.get_index_set<1,0>();
  const auto& v_q = mesh_c.get_index_set<2,0>();

  Index nf_v = nc_v + nc_e + nc_q;
  Index nnze = nc_v + 2*nc_e + 4*nc_q;

  Adjacency::Graph graph(nf_v, nc_v, nnze);

  Index* dom_ptr = graph.get_domain_ptr();
  Index* img_idx = graph.get_image_idx();

  for(Index i(0); i < nc_v; ++i)
  {
    dom_ptr[i] = i;
    img_idx[i] = i;
  }

  for(Index i(0); i < nc_e; ++i)
  {
    Index k = dom_ptr[nc_v+i] = nc_v + 2*i;
    img_idx[k+0] = v_e(i,0);
    img_idx[k+1] = v_e(i,1);
  }

  for(Index i(0); i < nc_q; ++i)
  {
    Index k = dom_ptr[nc_v+nc_e+i] = nc_v + 2*nc_e + 4*i;
    img_idx[k+0] = v_q(i,0);
    img_idx[k+1] = v_q(i,1);
    img_idx[k+2] = v_q(i,2);
    img_idx[k+3] = v_q(i,3);
  }
  dom_ptr[nf_v] = nnze;

  return graph;
}

Adjacency::Graph expand_prol(const Adjacency::Graph& prol, const Index n)
{
  XASSERT(prol.get_num_nodes_domain() <= n);

  // expand prolongation matrix
  const Index pnr = prol.get_num_nodes_domain();
  const Index pnc = prol.get_num_nodes_image();
  const Index pnz = prol.get_num_indices();
  const Index* dom_ptr_p = prol.get_domain_ptr();
  const Index* img_idx_p = prol.get_image_idx();

  // allocate expanded matrix
  Adjacency::Graph  ex_prol(n, n, pnz + n - pnc);
  Index* dom_ptr_q = ex_prol.get_domain_ptr();
  Index* img_idx_q = ex_prol.get_image_idx();

  // copy coarse identity
  for(Index i(0); i < pnc; ++i)
  {
    dom_ptr_q[i] = dom_ptr_p[i];
    for(Index j(dom_ptr_p[i]); j < dom_ptr_p[i+1]; ++j)
    {
      img_idx_q[j] = img_idx_p[j];
    }
  }
  dom_ptr_q[pnc] = dom_ptr_p[pnc];

  // copy interpolation + append identity
  for(Index i(pnc); i < pnr; ++i)
  {
    Index jq = dom_ptr_q[i];
    for(Index j(dom_ptr_p[i]); j < dom_ptr_p[i+1]; ++j, ++jq)
    {
      img_idx_q[jq] = img_idx_p[j];
    }
    img_idx_q[jq] = i;
    dom_ptr_q[i+1] = ++jq;
  }

  // append identity
  for(Index i(pnr); i < n; ++i)
  {
    Index jq = dom_ptr_q[i];
    img_idx_q[jq] = i;
    dom_ptr_q[i+1] = ++jq;
  }

  // export matrix
  return ex_prol;
}

void run()
{
  Index level_min = 5;
  Index level_max = 10;

  MeshType* mesh = nullptr;

  std::deque<Adjacency::Graph> prols;

  {
    Geometry::RefinedUnitCubeFactory<MeshType> factory(level_min);
    mesh = new MeshType(factory);
  }

  for(Index i(level_min); i < level_max; ++i)
  {
    MeshType* mesh_c = mesh;
    prols.push_front(build_prol(*mesh_c));
    {
      Geometry::StandardRefinery<MeshType> refinery(*mesh_c);
      mesh = new MeshType(refinery);
    }
    delete mesh_c;
  }

  const auto& vidx = mesh->get_index_set<2,0>();

  Adjacency::Graph neighbours;

  TimeStamp ts_0;
  Adjacency::Graph verts_at_elem(Adjacency::RenderType::as_is, vidx);
  TimeStamp ts_1;
  std::cout << "1 AS_IS        : " << ts_1.elapsed_string(ts_0).pad_front(7) << std::endl;

  Adjacency::Graph elems_at_vert(Adjacency::RenderType::transpose, vidx);
  TimeStamp ts_2;
  std::cout << "2 TRANSPOSE    : " << ts_2.elapsed_string(ts_1).pad_front(7) << std::endl;

  neighbours = render_injectify_new(elems_at_vert, verts_at_elem);
  TimeStamp ts_3;
  std::cout << "3 INJECT NEW   : " << ts_3.elapsed_string(ts_2).pad_front(7) << std::endl;

  neighbours = render_injectify_old(elems_at_vert, verts_at_elem);
  TimeStamp ts_4;
  std::cout << "4 INJECT OLD   : " << ts_4.elapsed_string(ts_3).pad_front(7) << std::endl;

  neighbours = render_injectify_mask(elems_at_vert, verts_at_elem);
  TimeStamp ts_5;
  std::cout << "5 INJECT MASK  : " << ts_5.elapsed_string(ts_4).pad_front(7) << std::endl;

  neighbours.sort_indices();
  TimeStamp ts_6;
  std::cout << "6 LIN INS SORT : " << ts_6.elapsed_string(ts_5).pad_front(7) << std::endl;

  Adjacency::Graph matrix = neighbours.clone();
  Adjacency::Graph tmp1, tmp2;

  Index lev(level_max);
  for(const auto& p : prols)
  {
    std::cout << String(40, '*') << std::endl;
    std::cout << "--> " << (--lev) << std::endl;

    // expand prol
    Adjacency::Graph prol = expand_prol(p, matrix.get_num_nodes_domain());
    Adjacency::Graph rest(Adjacency::RenderType::transpose, prol);

    ts_1.stamp();
    tmp1 = render_injectify_old(matrix, prol);
    tmp2 = render_injectify_old(rest, tmp1);
    ts_2.stamp();
    sort_old(tmp2);
    ts_3.stamp();
    std::cout << "1 INJECT OLD   : " << ts_2.elapsed_string(ts_1).pad_front(7) //<< std::endl;
      << " + " << ts_3.elapsed_string(ts_2).pad_front(7) //<< std::endl;
      << " = " << ts_3.elapsed_string(ts_1).pad_front(7) << std::endl;

    ts_2.stamp();
    tmp1 = render_injectify_mask(matrix, prol);
    tmp2 = render_injectify_mask(rest, tmp1);
    ts_3.stamp();
    sort_new(tmp2);
    ts_4.stamp();
    std::cout << "2 INJECT MASK  : " << ts_3.elapsed_string(ts_2).pad_front(7) //<< std::endl;
      << " + " << ts_4.elapsed_string(ts_3).pad_front(7) //<< std::endl;
      << " = " << ts_4.elapsed_string(ts_2).pad_front(7) << std::endl;

    tmp1 = render_injectify_new(matrix, prol);
    tmp2 = render_injectify_new(rest, tmp1);
    ts_4.stamp();
    std::cout << "3 INJECT NEW   : " << ts_4.elapsed_string(ts_3).pad_front(7) << std::endl;

    std::cout << "NNZE: " << matrix.get_num_indices() << "  >>  " << tmp2.get_num_indices() << std::endl;

    matrix = std::move(tmp2);
  }

  delete mesh;
}

int main(int argc, char** argv)
{
  Runtime::initialise(argc, argv);
  run();
  return Runtime::finalise();
}
