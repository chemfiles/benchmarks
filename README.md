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

## Running the benchmarks on your system

```bash
git clone https://github.com/chemfiles/benchmarks
cd benchmarks
pip install -r requirements.txt  # or use conda to install everything
python bench.py

# run without the ramdisk (Windows or Linux without sudo)
python bench.py --no-ramdisk
```

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

Each library is tested with multiple format (text formats, direct reading of
compressed text, and binary formats). For each library, the time required to
open and read the full trajectory is reported; as well as the timing ratio with
respect to chemfiles (higher is slower).

### macOS 10.14, SSD APFS, Intel Core i7-4870HQ

| file               | chemfiles (v0.10.1)   | MDAnalysis (v1.1.1)   | openbabel (v3.1.0)    | ase (v3.22.0)         |
|--------------------|-----------------------|-----------------------|-----------------------|-----------------------|
| water.xyz          | 10.92ms               | 50.54ms - **4.63x**   | 84.60ms - **7.75x**   | 124.99ms - **11.45x** |
| water.xyz.gz       | 28.14ms               | 61.55ms - **2.19x**   | 584.87ms - **20.78x** | 908.83ms - **32.30x** |
| 1vln-triclinic.pdb | 29.18ms               | 238.49ms - **8.17x**  | 519.53ms - **17.81x** | 159.65ms - **5.47x**  |
| lmsd.sdf           | 13.98ms               |                       | 99.53ms - **7.12x**   | 6.96ms - 0.50x        |
| vesicles.gro       | 31.08ms               | 523.99ms - **16.86x** | 339.95ms - **10.94x** | 456.95ms - **14.70x** |
| ubiquitin.xtc      | 206.21ms              | 107.77ms - 0.52x      |                       |                       |
| water.nc           | 3.78ms                | 5.11ms - **1.35x**    |                       | 126.16ms - **33.40x** |

### Ubuntu Linux 18.04, HDD ext4, Intel Core i5-4460S

| file               | chemfiles (v0.10.1)   | MDAnalysis (v1.1.1)   | openbabel (v3.1.0)    | ase (v3.22.0)         |
|--------------------|-----------------------|-----------------------|-----------------------|-----------------------|
| water.xyz          | 11.33ms               | 41.23ms - **3.64x**   | 93.63ms - **8.27x**   | 89.92ms - **7.94x**   |
| water.xyz.gz       | 24.40ms               | 56.87ms - **2.33x**   | 490.23ms - **20.09x** | 895.07ms - **36.69x** |
| 1vln-triclinic.pdb | 21.62ms               | 189.13ms - **8.75x**  | 559.01ms - **25.86x** | 127.91ms - **5.92x**  |
| lmsd.sdf           | 10.22ms               |                       | 79.04ms - **7.73x**   | 5.14ms - 0.50x        |
| vesicles.gro       | 25.92ms               | 437.74ms - **16.89x** | 507.86ms - **19.60x** | 364.06ms - **14.05x** |
| ubiquitin.xtc      | 203.99ms              | 119.48ms - 0.59x      |                       |                       |
| water.nc           | 3.49ms                | 3.65ms - **1.05x**    |                       | 132.95ms - **38.04x** |

## Contributions

Improvements to the benchmarks and more files to be tested are very welcome!
