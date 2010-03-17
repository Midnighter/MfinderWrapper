#!/usr/bin/env python
# encoding: utf-8


"""
setup.py
 
Created by Moritz Beber 2010-02-09.
Copyright (c) 2010 Jacobs University of Bremen. All rights reserved.
"""


import os
from distutils.core import setup, Extension


mfinder_path = "mfinder"


setup(
    name = "mfinder_wrapper",
    version = "0.1",
    description = "Motif Analysis Wrapper with Custom Randomisation",
    author = "Moritz Beber",
    author_email = "moritz (dot) beber (at) googlemail (dot) com",
    url = "http://github.com/Midnighter/mfinder_wrapper",
    packages = [mfinder_path],
    ext_package = mfinder_path,
    ext_modules = [Extension("_mfinder",
        sources = [
            os.path.join(mfinder_path, "mfinder_wrap.c"),
            os.path.join(mfinder_path, "interface.c"),
            os.path.join(mfinder_path, "clustering.c"),
            os.path.join(mfinder_path, "grassberger.c"),
            os.path.join(mfinder_path, "hash.c"),
            os.path.join(mfinder_path, "list.c"),
            os.path.join(mfinder_path, "mat.c"),
            os.path.join(mfinder_path, "metropolis.c"),
            os.path.join(mfinder_path, "motif_ids.c"),
            os.path.join(mfinder_path, "output.c"),
            os.path.join(mfinder_path, "permutation.c"),
            os.path.join(mfinder_path, "prob.c"),
            os.path.join(mfinder_path, "random.c"),
            os.path.join(mfinder_path, "results.c"),
            os.path.join(mfinder_path, "role.c"),
            os.path.join(mfinder_path, "stubs.c"),
            os.path.join(mfinder_path, "switches.c")
            ],
        extra_compile_args = ["-O3"],
        include_dirs=[mfinder_path,],
        library_dirs=[mfinder_path,],
        libraries=['m'],
        define_macros=[("NDEBUG", None), ("UNIX", '1')],
        )],
        py_modules=["mfinder_wrapper", os.path.join(mfinder_path, "mfinder"),
            "interface_test", "randomise_flow", "randomise_metb"],
    )

