#!/usr/bin/env python
try:
  import os
  from setuptools import setup, find_packages
except ImportError:
  from distutils.core import setup
from os import path
here = path.abspath(path.dirname(__file__))
try:
  with open("README.md", "r") as fh:
      long_description = fh.read()
  except:
    pass
  else:
    long_description = 'Fast k-mer based tool for alignment and assembly-free multi locus sequence typing (MLST) directly from genome sequencing reads.'

setup(
  name = 'stringMLST',
  scripts = ['stringMLST.py'],
  version = '0.6',
  description = 'Fast k-mer based tool for alignment and assembly-free multi locus sequence typing (MLST) directly from genome sequencing reads.',
  long_description=long_description,
  long_description_content_type="text/markdown",
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
