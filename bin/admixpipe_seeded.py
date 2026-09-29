#!/usr/bin/env python3
"""Run admixturePipeline.py with reproducible, distinct ADMIXTURE seeds.

AdmixPipe gives each ADMIXTURE run `-s np.random.randint(1000000)` from
numpy's unseeded global generator. This wrapper replaces that draw with the
next value of a permutation of 1..999999 shuffled by a generator seeded with
SEED, so every K and replicate gets its own seed and a run can be repeated
exactly.

Usage: admixpipe_seeded.py SEED /path/to/admixturePipeline.py [args...]
"""
import os
import runpy
import sys

import numpy as np

seed, script, *args = sys.argv[1:]
run_seeds = iter(np.random.default_rng(int(seed)).permutation(np.arange(1, 1_000_000)))
numpy_randint = np.random.randint


def randint(low, high=None, size=None, dtype=int):
    """Return the next run seed for AdmixPipe's randint(1000000) call."""
    if (low, high, size) == (1_000_000, None, None):
        return int(next(run_seeds))
    return numpy_randint(low, high, size, dtype)


np.random.randint = randint
sys.argv = [script, *args]
sys.path.insert(0, os.path.dirname(os.path.abspath(script)))
runpy.run_path(script, run_name="__main__")
