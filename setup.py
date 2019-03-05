#!/usr/bin/env python
try:
  import os
  from setuptools import setup, find_packages
except ImportError:
  from distutils.core import setup
from os import path
here = path.abspath(path.dirname(__file__))
def readme(file):
  with open(path.join(here, 'README.md')) as fh:
      long_description_text = fh.read()
  return(long_description_text)

setup(
  name = 'stringMLST',
  scripts = ['stringMLST.py'],
  version = "0.6.1",
  description = 'Fast k-mer based tool for alignment and assembly-free multi locus sequence typing (MLST) directly from genome sequencing reads.',
  long_description=readme('README.md'),
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
