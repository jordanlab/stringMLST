#!/usr/bin/env python
try:
  import os
  from setuptools import setup, find_packages
except ImportError:
  from distutils.core import setup
from os import path
here = path.abspath(path.dirname(__file__))
long_description="Fast k-mer based tool for alignment and assembly-free multi locus sequence typing (MLST) directly from genome sequencing reads."
with open(path.join(here, 'README.rst')) as f:
    long_description = f.read()

setup(
  name = 'stringMLST',
  scripts = ['stringMLST.py'], 
  version = '0.4',
  description = 'Fast k-mer based tool for alignment and assembly-free multi locus sequence typing (MLST) directly from genome sequencing reads.',
  long_description=long_description,
  author = 'Jordan Lab',
  author_email = 'pypi@atc.io',
  url = 'https://github.com/jordanlab/stringMLST',
  keywords = ['MLST', 'kmer', "NGS", "stringMSLT"], 
  classifiers = [
      'Programming Language :: Python :: 2.7',
      'Programming Language :: Python :: 3.5',
  ],
  install_requires=['lxml','pyfaidx'],
)
