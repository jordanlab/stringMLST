#!/usr/bin/env python
try:
  import os
  from setuptools import setup, find_packages
except ImportError:
  from distutils.core import setup

setup(
  name = 'stringMLST',
  scripts = ['stringMLST.py'], # this must be the same as the name above
  version = '0.3.5',
  description = 'Fast k-mer based tool for alignment and assembly-free multi locus sequence typing (MLST) directly from genome sequencing reads.',
  author = 'Jordan Lab',
  author_email = 'pypi@atc.io',
  url = 'https://github.com/anujg1991/stringMLST',
  keywords = ['testing', 'MLST', 'kmer', "NGS", "stringMSLT"], 
  classifiers = [
      'Programming Language :: Python :: 2.7',
      'Programming Language :: Python :: 3.5',
  ],
  install_requires=['lxml','pyfaidx'],
)