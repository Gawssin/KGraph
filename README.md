# KGraph: Concurrent Graph Query Processing with Memoization on Graph


Organization
--------

This repository contains the code for KGraph. It is implemented based on [**ForkGraph**](https://github.com/Xtra-Computing/ForkGraph).


Compilation
--------

**Compiler:**
* `g++ >= 7.5.0`


**Build system:**
* `CMake >= 3.12`


**To build:**
```sh
$ cd KGraph/
$ mkdir build
$ cd build && cmake ..
$ make -j
```


Input Graph Formats
-----------
The input graph comprises three files: `<inFile_inter>`, `<inFile_intra>` and `<partFile>`. The `<inFile_inter>` file contains the inter-partition edges, the `<inFile_intra>` file contains the intra-partition edges, and the `<partFile>` file contains the partition information. The number in the i-th line corresponds to the partition ID of the vertex with vertex ID i-1. 

`<inFile_inter>` and `<inFile_intra>` should be in the weighted adjacency graph format used by the [Problem Based Benchmark suite](http://www.cs.cmu.edu/~pbbs/benchmarks/graphIO.html) and [Ligra](https://github.com/jshun/ligra). We provide some example graphs in this format in the [shared OneDrive folder](https://1drv.ms/f/s!Avm-79-w1WHbgZMrEm3hMYBFnwY4vQ?e=JxtXju).


Usage
-----------
> **Note:**
> - For traversal-based applications, we implement two versions: `sparse` and `dense`. The `sparse` version is suitable for road networks and other graphs with low average degree, while the `dense` version is suitable for power-law graphs.
> - The choose of some parameters is crucial for the performance of KGraph and should be carefully tuned based on hardware resource (LLC size, number of processors and CPUs). We provide some recommended values for these parameters.

### Common parameters for all applications
```
$ ./kgraph-app -p <#part> -ps <strategy> -pn <#pivot> -qn <#query> -qf <inSrc> <inFile_inter> <inFile_intra> <partFile>
```
| Argument  | Description |
| :-----| :---- |
| **-p** | The number of partitions |
| **-ps<br>-pn** | The strategy for seleting pivots <br> The number of pivots |
| **-qn<br>-qf** | The number of queries to be evaluated <br> The file containing source vertices of queries |
| **\<inFile_inter\><br>\<inFile_intra\><br>\<partFile\>** | The input Graph |

### Specific parameters for each application

| Argument | Application&emsp;&emsp;&emsp;&emsp; | Description |
| :-----| :----- | :---- |
| **-h** | `kgraph-sssp-sparse`<br> `kgraph-sswp-sparse` <br> `kgraph-ssr-sparse` | The heuristic-Based Yielding used in [ForkGraph](https://github.com/Xtra-Computing/ForkGraph). We recomand setting it to 3200 for `Ca` and `Us`, and 204800 for `Eu`. |
| **-mr** | `kgraph-sssp-dense`<br> `kgraph-sswp-dense` <br> `kgraph-ssr-dense` | The extent of memoization for a pivot. The larger this value is, the more information the pivot will memoize, but it will also consume more memory. `[default: 2]` |
| **-wl** | `kgraph-deepwalk` <br>  `kgraph-ppr` | The walk length. `[default: 1000]` |
| **-ml** | `kgraph-deepwalk` <br>  `kgraph-ppr` | The length of a path memoized by a pivot that starts from itself. We recomand setting it to 4 for road networks and 2 for power-law graphs. |
