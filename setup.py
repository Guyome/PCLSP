#!/usr/bin/env python

from setuptools import setup
from os.path import join
from os import listdir

setup(name='PclspSolver',
    version='0.1',
    description='Heuristic solver for PCLSP',
    author='Guillaume Lanquepin',
    author_email='guillauem.lanquepin@himolde.net',
    package_dir = {'.': ''},
    py_modules=['solvers','tools','plots','pclsp_solver'],
    data_files=[('LUBP', ['LUBP/lupb.cpp']),
        ('LLBP',['LLBP/main.cpp', 'LLBP/LbIpopt.cpp', 'LLBP/LbIpopt.hpp'])],
    )

