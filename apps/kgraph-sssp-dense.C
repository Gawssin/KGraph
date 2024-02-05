#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#include "ligra.h"
#include <vector>
#include <queue>
#include "time.hpp"

int memoizationRange;

typedef std::pair<intE, uintE> WN;
typedef std::pair<int, std::pair<int, int>> WNH;

int schedule_partition(int *closest, int **closestPart, int NP, int k)
{
  int pi = -1;

  for (int i = 0; i < NP; i++)
  {
    closest[i] = INT_MAX;
  }

  for (int qi = 0; qi < k; ++qi)
  {
    for (int i = 0; i < NP; i++)
    {
      closest[i] = std::min(closest[i], closestPart[qi][i]);
    }
  }

  int minDis = INT_MAX;
  for (int i = 0; i < NP; i++)
  {
    if (closest[i] < minDis)
    {
      minDis = closest[i];
      pi = i;
    }
  }

  return pi;
}

template <class vertex>
void pivotSSP(graph<vertex> &GB, int pivotId, std::vector<WN> &pivotSP, int *pivotDist, std::priority_queue<WNH, std::vector<WNH>, std::greater<WNH>> &pq, int *vis)
{
  int partId = id2part[pivotId];
  int partN = partNodes[partId].size();
  int infDist = INT_MAX / 2;
  std::fill(pivotDist, pivotDist + partN, infDist);

  pivotDist[idPartOffset[pivotId]] = 0;

  pq.push(std::make_pair(0, std::make_pair(pivotId, 0)));

  vertex &vp = GB.V[pivotId];
  int outDegree = vp.getOutDegree();

  int tolNbr = outDegree * memoizationRange;
  pivotSP.reserve(tolNbr);

  std::fill(vis, vis + partN, -1);

  int nbrCnt = 0;

  while (!pq.empty())
  {
    intE curDist = pq.top().first;
    uintE vi = pq.top().second.first;
    intE level = pq.top().second.second;

    pq.pop();
    if (curDist > pivotDist[idPartOffset[vi]])
      continue;

    if (vi != pivotId)
    {
      if (vis[idPartOffset[vi]] != -1)
      {
        pivotSP[vis[idPartOffset[vi]]] = std::make_pair(curDist, vi);
      }
      else
      {
        vis[idPartOffset[vi]] = pivotSP.size();
        pivotSP.push_back(std::make_pair(curDist, vi));
        ++nbrCnt;
      }
    }

    if (nbrCnt >= tolNbr)
    {
      break;
    }

    vertex &v = GB.V[vi];
    auto degree = v.getOutDegree();
    for (int d = 0; d < degree; d++)
    {
      uintE ngh = v.getOutNeighbor(d);
      intE wgh = v.getOutWeight(d);

      if (pivotDist[idPartOffset[ngh]] > pivotDist[idPartOffset[vi]] + wgh)
      {
        pivotDist[idPartOffset[ngh]] = pivotDist[idPartOffset[vi]] + wgh;
        pq.push(std::make_pair(pivotDist[idPartOffset[ngh]], std::make_pair(ngh, level + 1)));
        if (level == 0)
        {
          vis[idPartOffset[ngh]] = pivotSP.size();
          pivotSP.push_back(std::make_pair(pivotDist[idPartOffset[ngh]], ngh));
          ++nbrCnt;
        }
      }
    }
  }

  while (!pq.empty())
    pq.pop();

  pivotSP.resize(pivotSP.size());
}

template <class vertex>
void Compute(graph<vertex> &GA, graph<vertex> &GB, int NP, std::vector<int> &pivot_l, commandLine P)
{
  std::cout << "start computing" << std::endl;

  memoizationRange = P.getOptionIntValue("-mr", 2);

  int pivotN = pivot_l.size();
  long n = GA.n;
  std::vector<long> s_l = start_l;
  long k = s_l.size();

  int *offsetGB = new int[NP + 1]();
  int maxPartN = 0;
  for (int i = 1; i < NP + 1; i++)
  {
    offsetGB[i] = offsetGB[i - 1] + partNodes[i - 1].size();
    maxPartN = std::max(maxPartN, int(partNodes[i - 1].size()));
  }

  int closest[NP];

  std::vector<intE *> ShortestPathLen_l;
  ShortestPathLen_l.reserve(k);

  bool **active = new bool *[k];

  int **closestPart = new int *[k]();

  parallel_for(int i = 0; i < k; i++)
  {
    closestPart[i] = new int[NP];
    fill(closestPart[i], closestPart[i] + NP, INT_MAX);
    active[i] = new bool[n]();
    active[i][s_l[i]] = true;

    ShortestPathLen_l[i] = newA(intE, n);
    fill(ShortestPathLen_l[i], ShortestPathLen_l[i] + n, INT_MAX / 2);

    ShortestPathLen_l[i][s_l[i]] = 0;
    int pid = id2part[s_l[i]];
    closestPart[i][pid] = 0;
  }

  int num_threads = omp_get_max_threads();
  std::cout << "num_threads: " << num_threads << std::endl;

  std::vector<WN> *pivotSP = new std::vector<WN>[pivotN];
  int *id2pivotIndex = new int[n]();
  std::fill(id2pivotIndex, id2pivotIndex + n, -1);

  int **pivotDist = new int *[num_threads];
  int **pivotVis = new int *[num_threads];
  for (int i = 0; i < num_threads; i++)
  {
    pivotDist[i] = new int[maxPartN];
    pivotVis[i] = new int[maxPartN];
  }

  std::priority_queue<WNH, std::vector<WNH>, std::greater<WNH>> *pivotPQ = new std::priority_queue<WNH, std::vector<WNH>, std::greater<WNH>>[num_threads];

  Timer t;
  t.tic();

  std::cout << "start memoizing pivot queries..." << std::endl;

  parallel_for(int i = 0; i < pivotN; i++)
  {
    int pivotId = pivot_l[i];
    id2pivotIndex[pivotId] = i;
    int partId = id2part[pivotId];
    int threadId = omp_get_thread_num();
    pivotSSP(GB, pivotId, pivotSP[i], pivotDist[threadId], pivotPQ[threadId], pivotVis[threadId]);
  }

  t.pause();
  t.print_ms("finish memoization");

  std::cout << "start evaluating user-input queries..." << std::endl;
  int pi = -1;
  int rounds = 0;

  while (pi = schedule_partition(closest, closestPart, NP, k), pi != -1)
  {
    rounds++;
    auto f = [&](int qi)
    {
      auto &ShortestPathLen = ShortestPathLen_l[qi];
      closestPart[qi][pi] = INT_MAX;

      for (int vi = offsetGB[pi]; vi < offsetGB[pi + 1]; ++vi)
      {
        if (active[qi][vi] == true)
        {
          active[qi][vi] = false;

          vertex &v = GA.V[vi];
          auto degree = v.getOutDegree();

          for (int d = 0; d < degree; ++d)
          {
            uintE ngh = v.getOutNeighbor(d);

            if (ShortestPathLen[ngh] > ShortestPathLen[vi] + v.getOutWeight(d))
            {
              ShortestPathLen[ngh] = ShortestPathLen[vi] + v.getOutWeight(d);
              active[qi][ngh] = true;
              int nghPart = id2part[ngh];

              closestPart[qi][nghPart] = std::min(closestPart[qi][nghPart], ShortestPathLen[ngh]);
            }
          }

          if (id2pivotIndex[vi] != -1)
          {
            int pivotIndex = id2pivotIndex[vi];

            for (int j = 0; j < pivotSP[pivotIndex].size(); ++j)
            {
              uintE ngh = pivotSP[pivotIndex][j].second;
              intE wgh = pivotSP[pivotIndex][j].first;
              if (ShortestPathLen[ngh] > ShortestPathLen[vi] + wgh)
              {
                ShortestPathLen[ngh] = ShortestPathLen[vi] + wgh;

                active[qi][ngh] = true;
                closestPart[qi][pi] = std::min(closestPart[qi][pi], ShortestPathLen[ngh]);
              }
            }
          }
          else
          {
            vertex &vb = GB.V[vi];
            auto degreeb = vb.getOutDegree();

            for (int d = 0; d < degreeb; ++d)
            {
              uintE ngb = vb.getOutNeighbor(d);
              if (ShortestPathLen[ngb] > ShortestPathLen[vi] + vb.getOutWeight(d))
              {
                ShortestPathLen[ngb] = ShortestPathLen[vi] + vb.getOutWeight(d);
                active[qi][ngb] = true;

                closestPart[qi][pi] = std::min(closestPart[qi][pi], ShortestPathLen[ngb]);
              }
            }
          }
        }
      }
    };

#pragma omp parallel for schedule(dynamic)
    for (int qi = 0; qi < k; ++qi)
    {
      if (closestPart[qi][pi] < INT_MAX)
      {
        f(qi);
      }
    }
  }

  std::cout << "finish evaluation" << std::endl;

  t.toc();
  t.print_ms("exuection time");
}

#pragma clang diagnostic pop