#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import math
import os
import shutil
import sys
import time
import warnings

import ase
import chemfiles
import MDAnalysis
import openbabel  # noqa
import pytraj
from ase import io  # noqa
from openbabel import openbabel as ob
from tabulate import tabulate

# Disable all kind of warnings
warnings.filterwarnings("ignore")
ob.obErrorLog.StopLogging()

MAX_REPETITIONS = 50
MAX_TIME_SECONDS = 10


def run_benchmark(package, path, function):
    print(f"    running {package}")
    # Warm up any file system caches
    runtime = 0
    for _ in range(3):
        before = time.perf_counter()
        try:
            function(path)
        except Exception as e:
            print(f"        can not read {os.path.basename(path)}: {e}")
            return {"path": path, "repetitions": 0, "time": float("nan")}
        after = time.perf_counter()
        runtime += after - before
        if runtime > MAX_TIME_SECONDS / 2:
            break

    # Repeat the benchmark multiple times
    repetitions = 0
    runtime = 0

    while repetitions < MAX_REPETITIONS and runtime < MAX_TIME_SECONDS:
        repetitions += 1
        before = time.perf_counter()

        function(path)

        after = time.perf_counter()
        runtime += after - before

    runtime /= repetitions
    return {
        "path": os.path.basename(path),
        "repetitions": repetitions,
        "time": runtime * 1e3,
    }


def bench_chemfiles(path):
    with chemfiles.Trajectory(path) as trajectory:
        for frame in trajectory:
            pass


def bench_mdanalysis(path):
    if path.endswith(".sdf"):
        raise Exception("format not supported")

    universe = MDAnalysis.Universe(path)
    for ts in universe.trajectory:
        pass


def bench_openbabel(path):
    if (
        path.endswith(".nc")
        or path.endswith(".xtc")
        or path.endswith(".dcd")
        or path.endswith(".trr")
    ):
        raise Exception("format not supported")

    mol = ob.OBMol()
    conversion = ob.OBConversion()
    conversion.ReadFile(mol, path)
    while conversion.Read(mol):
        pass


def bench_ase(path):
    if path.endswith(".dcd"):
        # ASE only supports DCD files created by CP2K
        raise Exception("format not supported")

    if path.endswith(".nc"):
        format = "netcdftrajectory"
    else:
        format = None

    for atoms in ase.io.iread(path, format=format):
        pass


def bench_pytraj(path):
    if path.endswith("adk.dcd"):
        topology = os.path.join(os.path.dirname(path), "adk.psf")
    elif path.endswith("water.nc"):
        topology = os.path.join(os.path.dirname(path), "water.pdb")
    elif path.endswith("ubiquitin.xtc") or path.endswith("ubiquitin.trr"):
        topology = os.path.join(os.path.dirname(path), "not-ubiquitin.pdb")

    elif path.endswith(".xyz") or path.endswith(".xyz.gz"):
        raise Exception("format not supported")
    elif path.endswith(".gro"):
        raise Exception("format not supported")
    else:
        topology = path

    trajectory = pytraj.TrajectoryIterator(path, topology)
    for step in trajectory:
        pass


def print_table(results):
    headers = [
        "{} (v{})".format(package, result["version"])
        for package, result in results.items()
    ]
    headers.insert(0, "file")

    data = []
    for i, file in enumerate(FILES):
        current = [os.path.basename(file)]
        chemfiles_time = 0
        for package, value in results.items():
            timing = value["results"][i]["time"]
            if math.isnan(timing):
                current.append("")
                continue

            if package == "chemfiles":
                chemfiles_time = timing
                current.append("{:.2f}ms".format(timing))
            else:
                ratio = timing / chemfiles_time
                if ratio > 1:
                    current.append("{:.2f}ms - **{:.2f}x**".format(timing, ratio))
                else:
                    current.append("{:.2f}ms - {:.2f}x".format(timing, ratio))

        data.append(current)

    print()
    print(tabulate(data, headers=headers, tablefmt="github"))


BENCHMARKS = {
    "chemfiles": bench_chemfiles,
    # "MDAnalysis": bench_mdanalysis,
    # "openbabel": bench_openbabel,
    # "ase": bench_ase,
    "pytraj": bench_pytraj,
}

FILES = [
    "files/water.xyz",
    "files/water.xyz.gz",
    "files/1vln-triclinic.pdb",
    "files/lmsd.sdf",
    "files/vesicles.gro",
    "files/ubiquitin.xtc",
    "files/ubiquitin.trr",
    "files/water.nc",
    "files/adk.dcd",
]


if __name__ == "__main__":
    results = {}
    for package in BENCHMARKS.keys():
        results[package] = {
            "version": eval(package).__version__,
            "results": [],
        }

    for file in FILES:
        print("reading {}".format(file))
        for package, function in BENCHMARKS.items():
            r = run_benchmark(package, file, function)
            results[package]["results"].append(r)

    print_table(results)
