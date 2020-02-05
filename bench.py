#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os
import gc
import sys
import math
import json
import time
import warnings
from shutil import copyfile

from ase import io
import ase
import chemfiles
import openbabel
import MDAnalysis

from ramdisk import ram_disk


warnings.filterwarnings("ignore")

MAX_REPETITIONS = 30
MAX_TIME_SECONDS = 10


def run_benchmark(path, function):
    # Warm up any file system caches
    runtime = 0
    for _ in range(3):
        before = time.perf_counter()
        try:
            function(path)
        except Exception as e:
            print(e)
            return {
                "path": path,
                "repetitions": 0,
                "time": float("nan")
            }
        after = time.perf_counter()
        runtime += after - before
        if runtime > MAX_TIME_SECONDS / 2:
            break

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
        "path": os.path.basename(path),
        "repetitions": repetitions,
        "time": runtime * 1e3,
    }


def bench_chemfiles(path):
    file = chemfiles.Trajectory(path)
    for step in range(file.nsteps):
        file.read()


def bench_mdanalysis(path):
    universe = MDAnalysis.Universe(path)
    for ts in universe.trajectory:
        pass


def bench_openbabel(path):
    mol = openbabel.OBMol()
    conversion = openbabel.OBConversion()
    conversion.ReadFile(mol, path)
    while conversion.Read(mol):
        pass


def bench_ase(path):
    for atoms in ase.io.iread(path):
        pass


def print_results(data):
    for package, value in data.items():
        print("{} ({})".format(package, value["version"]))
        for result in value["results"]:
            if math.isnan(result["time"]):
                continue
            print("    {:30} {:.2f}ms".format(result["path"], result["time"]))
    print()


BENCHMARKS = {
    "ase": bench_ase,
    "openbabel": bench_openbabel,
    "MDAnalysis": bench_mdanalysis,
    "chemfiles": bench_chemfiles,
}

FILES = [
    "files/1vln-triclinic.pdb",
    "files/water.xyz.gz",
    "files/water.xyz",
]


def main():
    results = {}
    for package in BENCHMARKS.keys():
        results[package] = {
            "version": eval(package).__version__,
            "results": [],
        }

    with ram_disk() as root:
        os.mkdir(os.path.join(root, "files"))
        for file in FILES:
            new_path = os.path.join(root, file)
            copyfile(src=file, dst=new_path)

            print("reading {}".format(file))
            for package, function in BENCHMARKS.items():
                r = run_benchmark(new_path, function)
                results[package]["results"].append(r)
        gc.collect()

    print_results(results)


if __name__ == '__main__':
    main()
