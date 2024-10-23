#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#include "ligra.h"
#include <vector>
#include <queue>
#include "time.hpp"

int maxLevel;
int maxPartN;

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
void pivotHopNeighbor(graph<vertex> &GB, int pivotId, std::vector<WN> &pivotSP, int *pivotDist, std::priority_queue<WNH, std::vector<WNH>, std::greater<WNH>> &pq, int *vis)
{
  int infDist = INT_MAX / 2;
  std::fill(pivotDist, pivotDist + maxPartN, infDist);

  pivotDist[idPartOffset[pivotId]] = 0;

  pivotSP.push_back(std::make_pair(0, pivotId));
  int curIndex = 0;
  while (curIndex < pivotSP.size())
  {
    int curLevel = pivotSP[curIndex].first;
    int curId = pivotSP[curIndex].second;

    if (curLevel + 1 > maxLevel)
    {
      break;
    }

    vertex &v = GB.V[curId];

    auto degree = v.getOutDegree();
    for (int d = 0; d < degree; d++)
    {
      uintE ngh = v.getOutNeighbor(d);

      if (pivotDist[idPartOffset[ngh]] > pivotDist[idPartOffset[curId]] + 1)
      {
        pivotDist[idPartOffset[ngh]] = pivotDist[idPartOffset[curId]] + 1;

        pivotSP.push_back(std::make_pair(pivotDist[idPartOffset[ngh]], ngh));
      }
    }
    ++curIndex;
  }
}

template <class vertex>
void Compute(graph<vertex> &GA, graph<vertex> &GB, int NP, std::vector<int> &pivot_l, commandLine P)
{

  std::cout << "start computing" << std::endl;

  maxLevel = P.getOptionIntValue("-ml", 2);

  int pivotN = pivot_l.size();
  long n = GA.n;
  std::vector<long> s_l = start_l;
  long k = s_l.size();

  int *offsetGB = new int[NP + 1]();
  maxPartN = 0;
  for (int i = 1; i < NP + 1; i++)
  {
    offsetGB[i] = offsetGB[i - 1] + partNodes[i - 1].size();
    maxPartN = std::max(maxPartN, int(partNodes[i - 1].size()));
  }

  int closest[NP];

  std::vector<intE *> BFSLevel_l;
  BFSLevel_l.reserve(k);

  bool **active = new bool *[k];

  int **closestPart = new int *[k]();

  parallel_for(int i = 0; i < k; i++)
  {
    closestPart[i] = new int[NP];
    fill(closestPart[i], closestPart[i] + NP, INT_MAX);
    active[i] = new bool[n]();
    active[i][s_l[i]] = true;

    BFSLevel_l[i] = newA(intE, n);
    fill(BFSLevel_l[i], BFSLevel_l[i] + n, INT_MAX);

    BFSLevel_l[i][s_l[i]] = 0;
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
    pivotHopNeighbor(GB, pivotId, pivotSP[i], pivotDist[threadId], pivotPQ[threadId], pivotVis[threadId]);
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
      auto &BFSLevel = BFSLevel_l[qi];
      closestPart[qi][pi] = INT_MAX;

      for (int vi : partNodes[pi])
      {
        if (active[qi][vi] == true)
        {
          active[qi][vi] = false;

          vertex &v = GA.V[vi];
          auto degree = v.getOutDegree();

          for (int d = 0; d < degree; ++d)
          {
            uintE ngh = v.getOutNeighbor(d);

            if (BFSLevel[ngh] > BFSLevel[vi] + 1)
            {
              BFSLevel[ngh] = BFSLevel[vi] + 1;
              active[qi][ngh] = true;
              int nghPart = id2part[ngh];

              closestPart[qi][nghPart] = std::min(closestPart[qi][nghPart], BFSLevel[ngh]);
            }
          }

          if (id2pivotIndex[vi] != -1)
          {

            int pivotIndex = id2pivotIndex[vi];
            for (int j = 0; j < pivotSP[pivotIndex].size(); ++j)
            {
              uintE ngh = pivotSP[pivotIndex][j].second;
              intE level = pivotSP[pivotIndex][j].first;
              if (BFSLevel[ngh] > BFSLevel[vi] + level)
              {
                BFSLevel[ngh] = BFSLevel[vi] + level;

                active[qi][ngh] = true;
                closestPart[qi][pi] = std::min(closestPart[qi][pi], BFSLevel[ngh]);
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
              if (BFSLevel[ngb] > BFSLevel[vi] + 1)
              {
                BFSLevel[ngb] = BFSLevel[vi] + 1;
                active[qi][ngb] = true;

                closestPart[qi][pi] = std::min(closestPart[qi][pi], BFSLevel[ngb]);
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
