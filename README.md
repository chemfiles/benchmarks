# Chemistry file readers benchmarks

This repository contains benchmarks comparing the performance of different
libraries when reading various file formats used in theoretical and
computational chemistry and biology. The benchmarks only test opening a file and
reading it in full, not runnning analysis on the resulting trajectory.

The code uses RAM disks (filesystems living in RAM) to try to abstract away as
much of the filesystem and SSD/HDD performance characteristics.

Libraries compared:

- [chemfiles](https://chemfiles.org)
- [MDAnalysis](https://mdanalysis.org)
- [OpenBabel](http://openbabel.org/)
- [ASE](https://wiki.fysik.dtu.dk/ase/)

## Caveats

Different libraries are able to read different formats, and can sometimes do
more work when opening and reading files.

We benchmark the code running from Python, even if some libraries have other
languages interfaces.

We use the version of each software as provided in PyPi or conda, without
special tuning for a given system.

When in doubt, run benchmark using with your own code and use cases!

## Results

In most of the tested cases, chemfiles is faster than other softwares.
Exceptions to this seems to be the SDF reader in ASE, and the XTC reader in
MDAnalysis.

### macOS 10.14, SSD APFS, Intel Core i7-4870HQ

| file               | chemfiles (v0.9.2)   | MDAnalysis (v0.20.1)   | openbabel (v3.0.0)    | ase (v3.19.0)         |
|--------------------|----------------------|------------------------|-----------------------|-----------------------|
| water.xyz          | 10.42ms              | 53.40ms - **5.13x**    | 81.04ms - **7.78x**   | 109.30ms - **10.49x** |
| water.xyz.gz       | 22.88ms              | 71.03ms - **3.10x**    | 479.78ms - **20.97x** | 808.83ms - **35.35x** |
| 1vln-triclinic.pdb | 26.08ms              | 237.38ms - **9.10x**   | 617.88ms - **23.69x** | 204.38ms - **7.84x**  |
| lmsd.sdf           | 13.58ms              |                        | 105.69ms - **7.79x**  | 7.13ms - 0.53x        |
| ubiquitin.xtc      | 267.82ms             | 115.60ms - 0.43x       |                       |                       |
| water.nc           | 3.92ms               | 5.25ms - **1.34x**     |                       | 161.53ms - **41.21x** |

### Ubuntu Linux 18.04, HDD ext4, Intel Core i5-4460S

| file               | chemfiles (v0.9.2)   | MDAnalysis (v0.20.1)   | openbabel (v3.0.0)    | ase (v3.19.0)          |
|--------------------|----------------------|------------------------|-----------------------|------------------------|
| water.xyz          | 22.75ms              | 44.11ms - **1.94x**    | 101.67ms - **4.47x**  | 98.48ms - **4.33x**    |
| water.xyz.gz       | 37.58ms              | 64.79ms - **1.72x**    | 538.10ms - **14.32x** | 1008.11ms - **26.82x** |
| 1vln-triclinic.pdb | 37.01ms              | 198.42ms - **5.36x**   | 626.34ms - **16.92x** | 166.15ms - **4.49x**   |
| lmsd.sdf           | 16.43ms              |                        | 86.51ms - **5.27x**   | 5.49ms - 0.33x         |
| ubiquitin.xtc      | 356.56ms             | 113.29ms - 0.32x       |                       |                        |
| water.nc           | 6.18ms               | 4.49ms - 0.73x         |                       | 155.26ms - **25.11x**  |

## Contributions

Improvements to the benchmarks and more files to be tested are very welcome!
