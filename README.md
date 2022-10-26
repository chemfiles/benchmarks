# Chemistry file readers benchmarks

This repository contains benchmarks comparing the performance of different
libraries when reading various file formats used in theoretical and
computational chemistry and biology. The benchmarks only test opening a file and
reading it in full, not running analysis on the resulting trajectory.

Libraries compared:

- [chemfiles](https://chemfiles.org)
- [MDAnalysis](https://mdanalysis.org)
- [OpenBabel](http://openbabel.org/)
- [ASE](https://wiki.fysik.dtu.dk/ase/)

## Running the benchmarks on your system

```bash
git clone https://github.com/chemfiles/benchmarks
cd benchmarks
./install-benchmarks.sh
./run-benchmarks.sh
```

## Caveats

Different libraries are able to read different formats, and can sometimes do
more work when opening and reading files.

We benchmark the code running from Python, even if some libraries have other
languages interfaces.

We use the version of each software as provided in conda, without special tuning
for a given system.

When in doubt, run benchmark using with your own code and use cases!

## Results

In most of the tested cases, chemfiles is faster than other softwares.
Exceptions to this seems to be the XTC and DCD readers in MDAnalysis, when
running on macOS.

Each library is tested with multiple format (text formats, direct reading of
compressed text, and binary formats). For each library, the time required to
open and read the full trajectory is reported; as well as the timing ratio with
respect to chemfiles (higher is slower).

### macOS 12.4, SSD APFS, M1 Max CPU

| file               | chemfiles (v0.10.3)   | ase (v3.22.1)         | MDAnalysis (v2.2.0)   | mdtraj (v1.9.7)        | openbabel (v3.1.0)    | parmed (v3.4.3)         | pytraj (v2.0.6.dev0)   |
|--------------------|-----------------------|-----------------------|-----------------------|------------------------|-----------------------|-------------------------|------------------------|
| water.xyz          | 4.87ms                | 46.00ms - **9.45x**   | 18.52ms - **3.80x**   |                        | 44.57ms - **9.16x**   |                         |                        |
| water.xyz.gz       | 12.75ms               | 505.79ms - **39.67x** | 25.41ms - **1.99x**   |                        | 411.59ms - **32.29x** |                         |                        |
| 1vln-triclinic.pdb | 11.47ms               | 66.54ms - **5.80x**   | 102.12ms - **8.90x**  | 334.51ms - **29.16x**  | 302.76ms - **26.39x** | 984.93ms - **85.86x**   | 20.38ms - **1.78x**    |
| lmsd.sdf           | 5.19ms                | 8.62ms - **1.66x**    |                       |                        | 41.62ms - **8.01x**   | 984.03ms - **189.44x**  | 7.04ms - **1.36x**     |
| vesicles.gro       | 12.91ms               | 181.34ms - **14.04x** | 311.04ms - **24.09x** | 215.69ms - **16.71x**  | 146.50ms - **11.35x** | 2707.65ms - **209.70x** |                        |
| molecules.mol2     | 5.18ms                |                       | 15.72ms - **3.03x**   |                        | 29.02ms - **5.60x**   | 325.65ms - **62.83x**   | 11.08ms - **2.14x**    |
| ubiquitin.xtc      | 79.41ms               |                       | 58.10ms - 0.73x       | 473.29ms - **5.96x**   |                       |                         | 69.61ms - 0.88x        |
| ubiquitin.trr      | 1.19ms                |                       | 24.43ms - **20.49x**  | 432.94ms - **363.21x** |                       |                         | 18.85ms - **15.81x**   |
| water.nc           | 0.92ms                | 36.44ms - **39.42x**  | 2.02ms - **2.18x**    | 4.21ms - **4.55x**     |                       |                         | 7.16ms - **7.75x**     |
| adk.dcd            | 7.18ms                |                       | 5.39ms - 0.75x        | 85.56ms - **11.92x**   |                       |                         | 12.74ms - **1.78x**    |

### Ubuntu Linux 20.04, SSD ext4, Intel Xeon 4214R CPU

| file               | chemfiles (v0.10.3)   | ase (v3.22.1)         | MDAnalysis (v2.2.0)   | mdtraj (v1.9.7)        | openbabel (v3.1.0)    | parmed (v3.4.3)         | pytraj (v2.0.6.dev0)   |
|--------------------|-----------------------|-----------------------|-----------------------|------------------------|-----------------------|-------------------------|------------------------|
| water.xyz          | 17.72ms               | 97.73ms - **5.51x**   | 121.17ms - **6.84x**  |                        | 143.75ms - **8.11x**  |                         |                        |
| water.xyz.gz       | 26.41ms               | 948.20ms - **35.90x** | 88.07ms - **3.33x**   |                        | 865.02ms - **32.75x** |                         |                        |
| 1vln-triclinic.pdb | 26.09ms               | 123.60ms - **4.74x**  | 315.89ms - **12.11x** | 620.56ms - **23.78x**  | 529.36ms - **20.29x** | 1860.58ms - **71.31x**  | 65.92ms - **2.53x**    |
| lmsd.sdf           | 15.92ms               | 20.88ms - **1.31x**   |                       |                        | 79.10ms - **4.97x**   | 1952.52ms - **122.66x** | 49.05ms - **3.08x**    |
| vesicles.gro       | 38.72ms               | 371.38ms - **9.59x**  | 568.72ms - **14.69x** | 393.73ms - **10.17x**  | 485.10ms - **12.53x** | 5229.44ms - **135.06x** |                        |
| molecules.mol2     | 14.66ms               |                       | 40.90ms - **2.79x**   |                        | 62.50ms - **4.26x**   | 602.32ms - **41.10x**   | 38.32ms - **2.61x**    |
| ubiquitin.xtc      | 169.95ms              |                       | 181.97ms - **1.07x**  | 895.10ms - **5.27x**   |                       |                         | 195.14ms - **1.15x**   |
| ubiquitin.trr      | 3.04ms                |                       | 123.98ms - **40.79x** | 793.42ms - **261.03x** |                       |                         | 81.19ms - **26.71x**   |
| water.nc           | 2.76ms                | 173.51ms - **62.83x** | 32.44ms - **11.75x**  | 30.65ms - **11.10x**   |                       |                         | 22.63ms - **8.19x**    |
| adk.dcd            | 15.87ms               |                       | 149.44ms - **9.42x**  | 149.02ms - **9.39x**   |                       |                         | 62.90ms - **3.96x**    |

## Contributions

Improvements to the benchmarks and more files to be tested are very welcome!
