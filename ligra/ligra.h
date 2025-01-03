// This code is part of the project "Ligra: A Lightweight Graph Processing
// Framework for Shared Memory", presented at Principles and Practice of 
// Parallel Programming, 2013.
// Copyright (c) 2013 Julian Shun and Guy Blelloch
//
// Permission is hereby granted, free of charge, to any person obtaining a
// copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights (to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject to
// the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
// OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
// LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
// OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
// WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#ifndef LIGRA_H
#define LIGRA_H
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <cstring>
#include <string>
#include <algorithm>
#include <vector>
#include <cassert>
#include <random>    // For std::default_random_engine
#include <chrono>    // For std::chrono
#include "parallel.h"
#include "gettime.h"
#include "utils.h"
#include "vertex.h"
#include "compressedVertex.h"
#include "vertexSubset.h"
#include "graph.h"
#include "IO.h"
#include "parseCommandLine.h"
#include "gettime.h"
#include "index_map.h"
#include "edgeMap_utils.h"
#include "parseCommandLine.h"
using namespace std;

//*****START FRAMEWORK*****

typedef uint32_t flags;
const flags no_output = 1;
const flags pack_edges = 2;
const flags sparse_no_filter = 4;
const flags dense_forward = 8;
const flags dense_parallel = 16;
const flags remove_duplicates = 32;
const flags no_dense = 64;
inline bool should_output(const flags& fl) { return !(fl & no_output); }

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDense(graph<vertex> GA, VS& vertexSubset, F &f, const flags fl) {
  using D = tuple<bool, data>;
  long n = GA.n;
  vertex *G = GA.V;
  if (should_output(fl)) {
    D* next = newA(D, n);
    auto g = get_emdense_gen<data>(next);
    parallel_for (long v=0; v<n; v++) {
      std::get<0>(next[v]) = 0;
      if (f.cond(v)) {
        G[v].decodeInNghBreakEarly(v, vertexSubset, f, g, fl & dense_parallel);
      }
    }
    return vertexSubsetData<data>(n, next);
  } else {
    auto g = get_emdense_nooutput_gen<data>();
    parallel_for (long v=0; v<n; v++) {
      if (f.cond(v)) {
        G[v].decodeInNghBreakEarly(v, vertexSubset, f, g, fl & dense_parallel);
      }
    }
    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapDenseForward(graph<vertex> GA, VS& vertexSubset, F &f, const flags fl) {
  using D = tuple<bool, data>;
  long n = GA.n;
  vertex *G = GA.V;
  if (should_output(fl)) {
    D* next = newA(D, n);
    auto g = get_emdense_forward_gen<data>(next);
    parallel_for(long i=0;i<n;i++) { std::get<0>(next[i]) = 0; }
    parallel_for (long i=0; i<n; i++) {
      if (vertexSubset.isIn(i)) {
        G[i].decodeOutNgh(i, f, g);
      }
    }
    return vertexSubsetData<data>(n, next);
  } else {
    auto g = get_emdense_forward_nooutput_gen<data>();
    parallel_for (long i=0; i<n; i++) {
      if (vertexSubset.isIn(i)) {
        G[i].decodeOutNgh(i, f, g);
      }
    }
    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse(graph<vertex>& GA, vertex* frontierVertices, VS& indices,
        uintT* degrees, uintT m, F &f, const flags fl) {
  using S = tuple<uintE, data>;
  long n = indices.n;
  S* outEdges;
  long outEdgeCount = 0;

  if (should_output(fl)) {
    uintT* offsets = degrees;
    outEdgeCount = sequence::plusScan(offsets, offsets, m);
    outEdges = newA(S, outEdgeCount);
    auto g = get_emsparse_gen<data>(outEdges);
    parallel_for (size_t i = 0; i < m; i++) {
      uintT v = indices.vtx(i), o = offsets[i];
      vertex vert = frontierVertices[i];
      vert.decodeOutNghSparse(v, o, f, g);
    }
  } else {
    auto g = get_emsparse_nooutput_gen<data>();
    parallel_for (size_t i = 0; i < m; i++) {
      uintT v = indices.vtx(i);
      vertex vert = frontierVertices[i];
      vert.decodeOutNghSparse(v, 0, f, g);
    }
  }

  if (should_output(fl)) {
    S* nextIndices = newA(S, outEdgeCount);
    if (fl & remove_duplicates) {
      if (GA.flags == NULL) {
        GA.flags = newA(uintE, n);
        parallel_for(long i=0;i<n;i++) { GA.flags[i]=UINT_E_MAX; }
      }
      auto get_key = [&] (size_t i) -> uintE& { return std::get<0>(outEdges[i]); };
      remDuplicates(get_key, GA.flags, outEdgeCount, n);
    }
    auto p = [] (tuple<uintE, data>& v) { return std::get<0>(v) != UINT_E_MAX; };
    size_t nextM = pbbs::filterf(outEdges, nextIndices, outEdgeCount, p);
    free(outEdges);
    return vertexSubsetData<data>(n, nextM, nextIndices);
  } else {
    return vertexSubsetData<data>(n);
  }
}

template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapSparse_no_filter(graph<vertex>& GA,
    vertex* frontierVertices, VS& indices, uintT* offsets, uintT m, F& f,
    const flags fl) {
  using S = tuple<uintE, data>;
  long n = indices.n;
  long outEdgeCount = sequence::plusScan(offsets, offsets, m);
  S* outEdges = newA(S, outEdgeCount);

  auto g = get_emsparse_no_filter_gen<data>(outEdges);

  // binary-search into scan to map workers->chunks
  size_t b_size = 10000;
  size_t n_blocks = nblocks(outEdgeCount, b_size);

  uintE* cts = newA(uintE, n_blocks+1);
  size_t* block_offs = newA(size_t, n_blocks+1);

  auto offsets_m = make_in_imap<uintT>(m, [&] (size_t i) { return offsets[i]; });
  auto lt = [] (const uintT& l, const uintT& r) { return l < r; };
  parallel_for(size_t i=0; i<n_blocks; i++) {
    size_t s_val = i*b_size;
    block_offs[i] = pbbs::binary_search(offsets_m, s_val, lt);
  }
  block_offs[n_blocks] = m;
  parallel_for (size_t i=0; i<n_blocks; i++) {
    if ((i == n_blocks-1) || block_offs[i] != block_offs[i+1]) {
      // start and end are offsets in [m]
      size_t start = block_offs[i];
      size_t end = block_offs[i+1];
      uintT start_o = offsets[start];
      uintT k = start_o;
      for (size_t j=start; j<end; j++) {
        uintE v = indices.vtx(j);
        size_t num_in = frontierVertices[j].decodeOutNghSparseSeq(v, k, f, g);
        k += num_in;
      }
      cts[i] = (k - start_o);
    } else {
      cts[i] = 0;
    }
  }

  long outSize = sequence::plusScan(cts, cts, n_blocks);
  cts[n_blocks] = outSize;

  S* out = newA(S, outSize);

  parallel_for (size_t i=0; i<n_blocks; i++) {
    if ((i == n_blocks-1) || block_offs[i] != block_offs[i+1]) {
      size_t start = block_offs[i];
      size_t start_o = offsets[start];
      size_t out_off = cts[i];
      size_t block_size = cts[i+1] - out_off;
      for (size_t j=0; j<block_size; j++) {
        out[out_off + j] = outEdges[start_o + j];
      }
    }
  }
  free(outEdges); free(cts); free(block_offs);

  if (fl & remove_duplicates) {
    if (GA.flags == NULL) {
      GA.flags = newA(uintE, n);
      parallel_for(size_t i=0;i<n;i++) { GA.flags[i]=UINT_E_MAX; }
    }
    auto get_key = [&] (size_t i) -> uintE& { return std::get<0>(out[i]); };
    remDuplicates(get_key, GA.flags, outSize, n);
    S* nextIndices = newA(S, outSize);
    auto p = [] (tuple<uintE, data>& v) { return std::get<0>(v) != UINT_E_MAX; };
    size_t nextM = pbbs::filterf(out, nextIndices, outSize, p);
    free(out);
    return vertexSubsetData<data>(n, nextM, nextIndices);
  }
  return vertexSubsetData<data>(n, outSize, out);
}

// Decides on sparse or dense base on number of nonzeros in the active vertices.
template <class data, class vertex, class VS, class F>
vertexSubsetData<data> edgeMapData(graph<vertex>& GA, VS &vs, F f,
    intT threshold = -1, const flags& fl=0) {
  long numVertices = GA.n, numEdges = GA.m, m = vs.numNonzeros();
  if(threshold == -1) threshold = numEdges/20; //default threshold
  vertex *G = GA.V;
  if (numVertices != vs.numRows()) {
    cout << "edgeMap: Sizes Don't match" << endl;
    abort();
  }
  if (vs.size() == 0) return vertexSubsetData<data>(numVertices);
  uintT* degrees = NULL;
  vertex* frontierVertices = NULL;
  uintT outDegrees = 0;
  if(threshold > 0) { 
    vs.toSparse();
    degrees = newA(uintT, m);
    frontierVertices = newA(vertex,m);
    {parallel_for (size_t i=0; i < m; i++) {
	  uintE v_id = vs.vtx(i);
	  vertex v = G[v_id];
	  degrees[i] = v.getOutDegree();
	  frontierVertices[i] = v;
    }}
    outDegrees = sequence::plusReduce(degrees, m);
//    if (fl & no_dense) {
//      std::cout << "Out Degree " << outDegrees << std::endl;
//    }
    if (outDegrees == 0) return vertexSubsetData<data>(numVertices);
  }
  if (!(fl & no_dense) && (m + outDegrees > threshold)) {
    if(degrees) free(degrees);
    if(frontierVertices) free(frontierVertices);
    vs.toDense();
    return (fl & dense_forward) ?
      edgeMapDenseForward<data, vertex, VS, F>(GA, vs, f, fl) :
      edgeMapDense<data, vertex, VS, F>(GA, vs, f, fl);
  } else {
    auto vs_out =
      (should_output(fl) && fl & sparse_no_filter) ? // only call snof when we output
      edgeMapSparse_no_filter<data, vertex, VS, F>(GA, frontierVertices, vs, degrees, vs.numNonzeros(), f, fl) :
      edgeMapSparse<data, vertex, VS, F>(GA, frontierVertices, vs, degrees, vs.numNonzeros(), f, fl);
    free(degrees); free(frontierVertices);
    return vs_out;
  }
}

// Regular edgeMap, where no extra data is stored per vertex.
template <class vertex, class VS, class F>
vertexSubset edgeMap(graph<vertex> GA, VS& vs, F f, intT threshold = -1, const flags& fl=0) {
  return edgeMapData<pbbs::empty>(GA, vs, f, threshold, fl);
}

// Packs out the adjacency lists of all vertex in vs. A neighbor, ngh, is kept
// in the new adjacency list if p(ngh) is true.
// Weighted graphs are not yet supported, but this should be easy to do.
template <class vertex, class P>
vertexSubsetData<uintE> packEdges(graph<vertex>& GA, vertexSubset& vs, P& p, const flags& fl=0) {
  using S = tuple<uintE, uintE>;
  vs.toSparse();
  vertex* G = GA.V; long m = vs.numNonzeros(); long n = vs.numRows();
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  auto degrees = array_imap<uintT>(m);
  granular_for(i, 0, m, (m > 2000), {
    uintE v = vs.vtx(i);
    degrees[i] = G[v].getOutDegree();
  });
  long outEdgeCount = pbbs::scan_add(degrees, degrees);
  S* outV;
  if (should_output(fl)) {
    outV = newA(S, vs.size());
  }

  bool* bits = newA(bool, outEdgeCount);
  uintE* tmp1 = newA(uintE, outEdgeCount);
  uintE* tmp2 = newA(uintE, outEdgeCount);
  if (should_output(fl)) {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t offset = degrees[i];
      auto bitsOff = &(bits[offset]); auto tmp1Off = &(tmp1[offset]);
      auto tmp2Off = &(tmp2[offset]);
      size_t ct = G[v].packOutNgh(v, p, bitsOff, tmp1Off, tmp2Off);
      outV[i] = make_tuple(v, ct);
    }
  } else {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t offset = degrees[i];
      auto bitsOff = &(bits[offset]); auto tmp1Off = &(tmp1[offset]);
      auto tmp2Off = &(tmp2[offset]);
      size_t ct = G[v].packOutNgh(v, p, bitsOff, tmp1Off, tmp2Off);
    }
  }
  free(bits); free(tmp1); free(tmp2);
  if (should_output(fl)) {
    return vertexSubsetData<uintE>(n, m, outV);
  } else {
    return vertexSubsetData<uintE>(n);
  }
}

template <class vertex, class P>
vertexSubsetData<uintE> edgeMapFilter(graph<vertex>& GA, vertexSubset& vs, P& p, const flags& fl=0) {
  vs.toSparse();
  if (fl & pack_edges) {
    return packEdges<vertex, P>(GA, vs, p, fl);
  }
  vertex* G = GA.V; long m = vs.numNonzeros(); long n = vs.numRows();
  using S = tuple<uintE, uintE>;
  if (vs.size() == 0) {
    return vertexSubsetData<uintE>(n);
  }
  S* outV;
  if (should_output(fl)) {
    outV = newA(S, vs.size());
  }
  if (should_output(fl)) {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t ct = G[v].countOutNgh(v, p);
      outV[i] = make_tuple(v, ct);
    }
  } else {
    parallel_for (size_t i=0; i<m; i++) {
      uintE v = vs.vtx(i);
      size_t ct = G[v].countOutNgh(v, p);
    }
  }
  if (should_output(fl)) {
    return vertexSubsetData<uintE>(n, m, outV);
  } else {
    return vertexSubsetData<uintE>(n);
  }
}



//*****VERTEX FUNCTIONS*****

template <class F, class VS, typename std::enable_if<
  !std::is_same<VS, vertexSubset>::value, int>::type=0 >
void vertexMap(VS& V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if(V.dense()) {
    parallel_for(long i=0;i<n;i++) {
      if(V.isIn(i)) {
        f(i, V.ithData(i));
      }
    }
  } else {
    parallel_for(long i=0;i<m;i++) {
      f(V.vtx(i), V.vtxData(i));
    }
  }
}

template <class VS, class F, typename std::enable_if<
  std::is_same<VS, vertexSubset>::value, int>::type=0 >
void vertexMap(VS& V, F f) {
  size_t n = V.numRows(), m = V.numNonzeros();
  if(V.dense()) {
    parallel_for(long i=0;i<n;i++) {
      if(V.isIn(i)) {
        f(i);
      }
    }
  } else {
    parallel_for(long i=0;i<m;i++) {
      f(V.vtx(i));
    }
  }
}

//Note: this is the version of vertexMap in which only a subset of the
//input vertexSubset is returned
template <class F>
vertexSubset vertexFilter(vertexSubset V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  V.toDense();
  bool* d_out = newA(bool,n);
  {parallel_for(long i=0;i<n;i++) d_out[i] = 0;}
  {parallel_for(long i=0;i<n;i++)
      if(V.d[i]) d_out[i] = filter(i);}
  return vertexSubset(n,d_out);
}

template <class F>
vertexSubset vertexFilter2(vertexSubset V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  bool* bits = newA(bool, m);
  V.toSparse();
  {parallel_for(size_t i=0; i<m; i++) {
    uintE v = V.vtx(i);
    bits[i] = filter(v);
  }}
  auto v_imap = make_in_imap<uintE>(m, [&] (size_t i) { return V.vtx(i); });
  auto bits_m = make_in_imap<bool>(m, [&] (size_t i) { return bits[i]; });
  auto out = pbbs::pack(v_imap, bits_m);
  out.alloc = false;
  free(bits);
  return vertexSubset(n, out.size(), out.s);
}

template <class data, class F>
vertexSubset vertexFilter2(vertexSubsetData<data> V, F filter) {
  long n = V.numRows(), m = V.numNonzeros();
  if (m == 0) {
    return vertexSubset(n);
  }
  bool* bits = newA(bool, m);
  V.toSparse();
  parallel_for(size_t i=0; i<m; i++) {
    auto t = V.vtxAndData(i);
    bits[i] = filter(std::get<0>(t), std::get<1>(t));
  }
  auto v_imap = make_in_imap<uintE>(m, [&] (size_t i) { return V.vtx(i); });
  auto bits_m = make_in_imap<bool>(m, [&] (size_t i) { return bits[i]; });
  auto out = pbbs::pack(v_imap, bits_m);
  out.alloc = false;
  free(bits);
  return vertexSubset(n, out.size(), out.s);
}



//cond function that always returns true
inline bool cond_true (intT d) { return 1; }

template<class vertex>
void Compute(graph<vertex>&, commandLine);

template<class vertex>
void Compute(graph<vertex>&, graph<vertex> &, int, std::vector<int> &, commandLine);


//#define GRP

std::vector<int> id2part;
std::vector<int> idPartOffset;
std::vector<std::vector<int>> partNodes;

// std::vector<int> id2group;
// std::vector<long> part_left;
// std::vector<long> part_right;
// std::vector<long> vid;
// std::vector<std::vector<long>> start_l_l;
std::vector<long> start_l;
std::vector<long> pivot_l;

void load_parts(const std::string &filename, int NP) {
#ifdef PRINT_DEBUG_MESSAGES
  std::cout << filename << std::endl;
#endif
  std::ifstream infile;
  infile.open(filename);
  std::vector<std::pair<int, int>> id2part_t;
  // std::vector<int> id_o2n;
  // std::vector<int> id_n2o;

  partNodes.clear();
  for (int i = 0; i < NP; i++) {
    partNodes.emplace_back(std::vector<int>());
  }

  int id = 0;
  if (infile.is_open()) {
    int num;
    while (infile >> num) {
      id2part_t.emplace_back(id, num);
      id++;
    }
  } else {
    std::cout << "Cannot open partition file 1!!!" << std::endl;
    exit(-1);
  }
//  assert(id == vertices);
  std::sort(id2part_t.begin(),
            id2part_t.end(),
            [](const std::pair<int, int> &a, const std::pair<int, int> &b) -> bool {
              return (a.second < b.second) || (a.second == b.second && a.first < b.first);
            });
  int N = id2part_t.size();
  // vid.resize(N);
  // id_o2n.resize(N);
  // id_n2o.resize(N);
  id2part.resize(N);
  idPartOffset.resize(N);
  for (int i = 0; i < N; i++) {
    // vid[i] = i;
    auto p = id2part_t[i];
    // id_o2n[p.first] = i;
    // id_n2o[i] = p.first;
    // assert(i == p.first);
    
    id2part[p.first] = p.second;
    
    // id2part[p.first] = p.second;
    
    // print id2part
    // if(i < 100)
    // std::cout << i << " " << p.first << " " << p.second << std::endl;
    
    idPartOffset[p.first] = partNodes[p.second].size();
    partNodes[p.second].push_back(p.first);
  }
  // part_left.resize(NP);
  // part_right.resize(NP);
  // for (int p_i = 0; p_i < NP; p_i++) {
  //   part_left[p_i] = part_id[p_i].front();
  //   part_right[p_i] = part_id[p_i].back();
  // }
}

void read_starts(std::string filename, int queryNum) {
#ifdef PRINT_DEBUG_MESSAGES
  std::cout << filename << std::endl;
#endif
  std::ifstream infile;
  infile.open(filename);
  start_l.clear();
  int id = 0;
  if (infile.is_open()) {
    int num;
    while (infile >> num) {
      start_l.emplace_back(num);
      id++;
      if(id >= queryNum)
        break;
    }
  } else {
    std::cout << "Cannot open partition file 2!!!" << std::endl;
    exit(-1);
  }
}

// void read_pivots(std::string filename) {
// #ifdef PRINT_DEBUG_MESSAGES
//   std::cout << filename << std::endl;
// #endif
//   std::ifstream infile;
//   infile.open(filename);
//   pivot_l.clear();
//   int id = 0;
//   if (infile.is_open()) {
//     int num;
//     while (infile >> num) {
//       pivot_l.emplace_back(num);
//       id++;
//     }
//   } else {
//     std::cout << "Cannot open pivot file 2!!!" << std::endl;
//     exit(-1);
//   }
// }

// template <class vertex>
// void part_group(graph<vertex>& GB) {
//   int id = 0;
//   while (id < start_l_l.size()) {
//     auto t = start_l_l[id];
//     if (t.size() > 500) {
//       id++;
//       continue;
//     }
//     int pid = id2part[t[0]];
//     for (auto vi = part_left[pid]; vi < part_right[pid]; vi++) {
//       vertex &v = GB.V[vi];
//       auto degree = v.getOutDegree();
//       for (auto d = 0; d < degree; d++) {
//         uintE ngh = v.getOutNeighbor(d);
//         int ppid = id2part[ngh];
//         int next_part = id + 1;
//         while (next_part < start_l_l.size()) {
//           if (ppid == id2part[start_l_l[next_part][0]]) {
//             t.insert(t.begin(), start_l_l[next_part].begin(), start_l_l[next_part].end());
//             start_l_l.erase(start_l_l.begin() + next_part);
//           } else {
//             next_part++;
//           }
//         }
//       }
//     }
//     start_l_l[id] = t;
//     id++;
//   }
// }

// template <class vertex>
// void random_group(graph<vertex>& GB) {
//   int id = 0;
//   while (id < start_l_l.size()) {
//     auto t = start_l_l[id];
//     if (t.size() > 500) {
//       id++;
//       continue;
//     }
//     auto next_part = id + rand() % (start_l_l.size() - id + 1) + 1;
//     if (next_part < start_l_l.size()) {
//       t.insert(t.begin(), start_l_l[next_part].begin(), start_l_l[next_part].end());
//       start_l_l.erase(start_l_l.begin() + next_part);
//     }
//     start_l_l[id] = t;
//     id++;
//   }
// }

template <class vertex>
std::vector<int> select_pivots(graph<vertex> &GA, graph<vertex> &GB, int ps, int pn, int PartitionN)
{
  int pivotStrategy = ps;
  int pivotNum = pn;
  std::vector<int> pivots;
  
  int maxDegree = 0;
  int boundNum = 0;
  
  int *degree = new int[GA.n]();
  bool *isBound = new bool[GA.n]();
  std::vector<int> boundNodes;
  
  for (int i = 0; i < GA.n; i++)
  {
    vertex &va = GA.V[i];
    vertex &vb = GB.V[i];
    // degree[i] = va.getOutDegree() + vb.getOutDegree() + va.getInDegree() + vb.getInDegree();
    degree[i] = va.getInDegree() + vb.getInDegree();
    
    maxDegree = std::max(maxDegree, degree[i]);
    
    if (va.getOutDegree() > 0 || va.getInDegree() > 0)
    {
      isBound[i] = true;
      boundNodes.push_back(i);
      boundNum++;
    }
  }

  
  
  
  if (pivotStrategy == -1)
  {
    double avgDegree = (GA.m + GB.m) * 2.0 / (GA.n + GB.n);
    if (avgDegree <= 12.432) pivotStrategy = 6;
    else pivotStrategy = 5;
  }
  
  
  std::vector<int> ids;
  for (int i = 0; i < GA.n; ++i) ids.push_back(i);
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::shuffle(ids.begin(), ids.end(), std::default_random_engine(seed));
  
  int* selectedPerPart = new int[PartitionN]();
  
  if (pivotStrategy == 1)
  {
    pivots.resize(pivotNum);
    std::copy_n(ids.begin(), pivotNum, pivots.begin());
  }
  if (pivotStrategy == 2)
  {

    for (int i = 0; i < GA.n; i++)
    {
      int pid = id2part[ids[i]];
      if (selectedPerPart[pid] < pivotNum / PartitionN)
      {
        pivots.push_back(ids[i]);
        selectedPerPart[pid]++;
      }
    }
  }
  if (pivotStrategy == 3)
  {
    std::vector<int> boundNodesTemp = boundNodes;
    std::shuffle(boundNodesTemp.begin(), boundNodesTemp.end(), std::default_random_engine(seed));
    
    for (int i = 0; i < boundNodesTemp.size(); i++)
    {
      int pid = id2part[boundNodesTemp[i]];
      if (selectedPerPart[pid] < pivotNum / PartitionN)
      {
        pivots.push_back(boundNodesTemp[i]);
        selectedPerPart[pid]++;
      }
    }
    
  }
  
  std::vector<int> degreeOrder;
  for (int i = 0; i < GA.n; i++) degreeOrder.push_back(i);

  
  std::sort(degreeOrder.begin(), degreeOrder.end(), [&](const int &a, const int &b) -> bool {
    return degree[a] > degree[b];
  });
  
  if (pivotStrategy == 4)
  {
    
    for (int i = 0; i < GA.n; i++)
    {
      int pid = id2part[degreeOrder[i]];
      if (selectedPerPart[pid] < pivotNum / PartitionN)
      {
        pivots.push_back(degreeOrder[i]);
        selectedPerPart[pid]++;
      }
    }
    
  }
  if (pivotStrategy == 5)
  {
    if (pivotNum == -1)
    {
      if (maxDegree > 2414375)
      {
        pivotNum = GA.n * 0.002;
      }
      else
      {
        if (GA.n > 4207239)
        {
          pivotNum = GA.n * 0.007;
        }
        else
        {
          pivotNum = GA.n * 0.001;
        }
        
      }
    }
    pivots.resize(pivotNum);
    std::copy_n(degreeOrder.begin(), pivotNum, pivots.begin());
    
  }
  if (pivotStrategy == 6)
  {
    for (int i = 0; i < GA.n; i++)
    {
      if(isBound[degreeOrder[i]] == false)
        continue;
      
      if (pivotNum == -1)
      {
        pivots.push_back(degreeOrder[i]);
      }
      else
      {
        int pid = id2part[degreeOrder[i]];
        if (selectedPerPart[pid] < pivotNum / PartitionN)
        {
          pivots.push_back(degreeOrder[i]);
          selectedPerPart[pid]++;
        }
      }
    }
  }
  
  std::cout << "pivotStrategy = " << pivotStrategy << std::endl;
  std::cout << "pivotNum = " << pivots.size() << std::endl;
  
  delete[] degree;
  delete[] isBound;
  delete[] selectedPerPart;
  return pivots;
}



int parallel_main(int argc, char* argv[]) {
  commandLine P(argc,argv,"-p <#part> -ps <strategy> -pn <#pivot> -qn <#query> -qf <inSrc> <inFile_inter> <inFile_intra> <partFile>");

  char* iFile1 = P.getArgument(2); // get third last
  char* iFile2 = P.getArgument(1); // get second last
  char* iPart = P.getArgument(0); // get last
  int NP = P.getOptionIntValue("-p", -1);
  
  char *iSource = P.getOptionValue("-qf");
  int queryNum = P.getOptionIntValue("-qn", -1);
  
  int pivotStrategy = P.getOptionIntValue("-ps", -1);
  int pivotNum = P.getOptionIntValue("-pn", -1);
  

  load_parts(iPart, NP);

  bool symmetric = P.getOptionValue("-s");
  bool compressed = P.getOptionValue("-c");
  bool binary = P.getOptionValue("-b");
  bool mmap = P.getOptionValue("-m");

  graph<asymmetricVertex> Ga = readGraph<asymmetricVertex>(iFile1,compressed,symmetric,binary,mmap); //asymmetric graph
  graph<asymmetricVertex> Gb = readGraph<asymmetricVertex>(iFile2,compressed,symmetric,binary,mmap); //asymmetric graph

  std::vector<int> pivots = select_pivots(Ga, Gb, pivotStrategy, pivotNum, NP);

  read_starts(iSource, queryNum);

  Compute(Ga, Gb, NP, pivots, P);
  Ga.del();
  Gb.del();
}
#endif