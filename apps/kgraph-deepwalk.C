#pragma clang diagnostic push
#pragma ide diagnostic ignored "openmp-use-default-none"

#include "ligra.h"
#include <vector>
#include <queue>
#include <random>
#include "time.hpp"

typedef std::pair<intE, uintE> WN;
typedef std::pair<double, std::pair<intE, uintE>> PWN;

int maxLevel;

int schedule_partition(int *closest, int **closestPart, int NP, int k)
{
  int pi = -1;

  for (int i = 0; i < NP; i++)
  {
    closest[i] = INT_MAX;
  }

  for (int i = 0; i < NP; i++)
  {
    for (int qi = 0; qi < k; ++qi)
    {
      closest[i] = std::min(closest[i], closestPart[qi][i]);
    }
  }

  int minLevel = INT_MAX;
  for (int i = 0; i < NP; i++)
  {
    if (closest[i] < minLevel)
    {
      minLevel = closest[i];
      pi = i;
    }
  }

  return pi;
}

template <class vertex>
void memoizePivot(graph<vertex> &GA, graph<vertex> &GB, int curId, int curLevel, double curProb, std::vector<double> &probabilityPivot,
                  std::vector<int> &PivotMemOrigin, std::vector<std::vector<int>> &pivotPath, std::vector<int> &pivotPathIndex)
{

  intE d1 = GA.V[curId].getOutDegree();
  intE d2 = GB.V[curId].getOutDegree();

  pivotPathIndex.push_back(curId);

  if (curLevel == 0 || d1 + d2 == 0)
  {
    probabilityPivot.push_back(curProb);
    PivotMemOrigin.push_back(maxLevel - curLevel);

    pivotPath.push_back(pivotPathIndex);
    return;
  }

  double prob = curProb * 1.0 / (d1 + d2);

  for (int j = 0; j < d1; j++)
  {
    memoizePivot(GA, GB, GA.V[curId].getOutNeighbor(j), curLevel - 1, prob, probabilityPivot, PivotMemOrigin, pivotPath, pivotPathIndex);
    pivotPathIndex.pop_back();
  }

  for (int j = 0; j < d2; j++)
  {
    memoizePivot(GA, GB, GB.V[curId].getOutNeighbor(j), curLevel - 1, prob, probabilityPivot, PivotMemOrigin, pivotPath, pivotPathIndex);
    pivotPathIndex.pop_back();
  }
}

template <class vertex>
void Compute(graph<vertex> &GA, graph<vertex> &GB, int NP, std::vector<int> &pivot_l, commandLine P)
{

  int walk_length = P.getOptionIntValue("-wl", 1000);
  std::cout << "walk_length: " << walk_length << std::endl;
  maxLevel = P.getOptionIntValue("-ml", 2); // maxLevel = 2 for power-law graphs and 4 for road networks
  std::cout << "maxLevel: " << maxLevel << std::endl;
  long pivotN = pivot_l.size();
  std::random_device rd;
  std::default_random_engine generator(rd());

  long n = GA.n;
  std::vector<long> s_l = start_l;
  long k = s_l.size();
  int *closest = new int[NP]();

  int *id2pivotIndex = new int[n]();
  std::fill(id2pivotIndex, id2pivotIndex + n, -1);

  std::vector<int> *PivotMem = new std::vector<int>[pivotN];

  std::vector<double> *probabilityPivot = new std::vector<double>[pivotN]();

  Timer t;
  t.tic();

  int num_threads = omp_get_max_threads();
  std::cout << "num_threads: " << num_threads << std::endl;

  std::vector<std::vector<int>> *pivotPath = new std::vector<std::vector<int>>[pivotN]();

  std::vector<int> *pivotPathIndex = new std::vector<int>[pivotN]();

#pragma omp parallel for schedule(dynamic)
  for (int i = 0; i < pivotN; i++)
  {
    int pivotId = pivot_l[i];
    id2pivotIndex[pivotId] = i;

    memoizePivot(GA, GB, pivotId, maxLevel, 1.0, probabilityPivot[i], PivotMem[i], pivotPath[i], pivotPathIndex[i]);
    if (PivotMem[i].size() <= 1)
    {
      id2pivotIndex[pivotId] = -1;
    }
  }

  t.pause();
  t.print_ms("finish memoization");

  int **closestPart = new int *[k];
  closestPart[0] = new int[k * NP];
  fill(closestPart[0], closestPart[0] + k * NP, INT_MAX);

  int **walkPath = new int *[k];
  walkPath[0] = new int[k * (walk_length + maxLevel)];
  int *walkLevel = new int[k]();

  std::cout << "walk_length: " << walk_length << std::endl;

  parallel_for(int i = 0; i < k; i++)
  {
    closestPart[i] = closestPart[0] + i * NP;

    walkPath[i] = walkPath[0] + i * (long long)(walk_length + maxLevel);
    walkPath[i][walkLevel[i]++] = s_l[i];

    int pid = id2part[s_l[i]];

    closestPart[i][pid] = 1;
  }

  t.pause();
  t.print_ms("allocate buffers");

  std::cout << "start evaluating user-input queries..." << std::endl;
  int pi = -1;

  while (pi = schedule_partition(closest, closestPart, NP, k), pi != -1)
  {

    auto f = [&](int qi, int vi, int vlevel)
    {
      if (vlevel >= walk_length)
      {
        return;
      }

      int nextvindex, passLevel;
      int pid, nextId;

      if (id2pivotIndex[vi] != -1)
      {

        int pivotIndex = id2pivotIndex[vi];

        std::discrete_distribution<> distribution(probabilityPivot[pivotIndex].begin(), probabilityPivot[pivotIndex].end());
        nextvindex = distribution(generator);

        passLevel = PivotMem[pivotIndex][nextvindex];

        for (int i = 0; i < passLevel; i++)
        {
          walkPath[qi][walkLevel[qi]++] = pivotPath[pivotIndex][nextvindex][i + 1];
        }
      }
      else
      {

        int d1 = GA.V[vi].getOutDegree();
        int d2 = GB.V[vi].getOutDegree();
        int pid, nextId;
        if (d1 + d2 == 0)
        {

          std::uniform_int_distribution<> distribution(0, GA.n - 1);
          nextId = distribution(generator);
        }
        else
        {

          std::uniform_int_distribution<int> distribution(0, d1 + d2 - 1);

          int nextvindex = distribution(generator);

          if (nextvindex < d1)
          {
            nextId = GA.V[vi].getOutNeighbor(nextvindex);
          }
          else if (nextvindex < d1 + d2)
          {
            nextId = GB.V[vi].getOutNeighbor(nextvindex - d1);
          }
        }
        walkPath[qi][walkLevel[qi]++] = nextId;
      }
    };

#pragma omp parallel for schedule(dynamic)
    for (int qi = 0; qi < k; ++qi)
    {
      int curId = walkPath[qi][walkLevel[qi] - 1];
      int curPid = id2part[curId];
      int vlevel = walkLevel[qi];

      if (curPid == pi)
      {

        closestPart[qi][curPid] = INT_MAX;

        int cnt = walk_length;
        while (cnt--)
        {
          f(qi, curId, vlevel);

          curId = walkPath[qi][walkLevel[qi] - 1];
          vlevel = walkLevel[qi];
        }

        if (vlevel < walk_length)
        {
          curPid = id2part[curId];
          closestPart[qi][curPid] = vlevel;
        }
      }
    }
  }

  t.toc();
  t.print_ms("exuection time");
}

#pragma clang diagnostic pop