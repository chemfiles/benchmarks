#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import json
import time
import sys
import warnings

from chemfiles import Trajectory, __version__ as CHEMFILES_VERSION
from MDAnalysis import Universe, __version__ as MDA_VERSION

warnings.filterwarnings("ignore")

MAX_REPETITIONS = 30
MAX_TIME_SECONDS = 10


def run_benchmark(path, function):
    # Warm up any file system caches
    for _ in range(3):
        function(path)

    # Repeat the benchmark multiple times
    repetitions = 0
    runtime = 0

    while (repetitions < MAX_REPETITIONS and runtime < MAX_TIME_SECONDS):
        repetitions += 1
        before = time.perf_counter()

        function(path)

        after = time.perf_counter()
        runtime += after - before

    runtime /= repetitions
    return {
        "path": path,
        "repetitions": repetitions,
        "time": runtime * 1e3
    }


def bench_chemfiles(path):
    file = Trajectory(path)
    for step in range(file.nsteps):
        file.read()


def bench_mdanalysis(path):
    universe = Universe(path)
    for ts in universe.trajectory:
        pass


def bench_chemfiles_open(path):
    file = Trajectory(path)


def bench_mdanalysis_open(path):
    universe = Universe(path)


def print_results(data):
    for code, res in data.items():
        print("{} ({})".format(code, res["version"]))
        for result in res["results"]:
            print("    {} {:.2f}ms".format(result["path"], result["time"]))
    print()


if __name__ == '__main__':
    FILES = [
        # "files/1aki.tng",
        # "files/1j8k.mmcif",
        "files/1vln-triclinic.pdb",
        "files/water.6.xyz.gz",
        "files/water.xyz",
    ]

    mda = []
    chemfiles = []
    for file in FILES:
        mda.append(run_benchmark(file, bench_mdanalysis_open))
        chemfiles.append(run_benchmark(file, bench_chemfiles_open))

    results = {
        "mda": {
            "version": MDA_VERSION,
            "results": mda,
        },
        "chemfiles": {
            "version": CHEMFILES_VERSION,
            "results": chemfiles,
        },
    }

    print("Opening files")
    print_results(results)

    mda = []
    chemfiles = []
    for file in FILES:
        mda.append(run_benchmark(file, bench_mdanalysis))
        chemfiles.append(run_benchmark(file, bench_chemfiles))

    results = {
        "mda": {
            "version": MDA_VERSION,
            "results": mda,
        },
        "chemfiles": {
            "version": CHEMFILES_VERSION,
            "results": chemfiles,
        },
    }

    print("Reading all steps")
    print_results(results)
