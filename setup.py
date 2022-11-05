#
# Copyright Tim Molteno 2017 tim@elec.ac.nz
#

from setuptools import setup, find_packages

import setuptools.command.test

with open('README.md') as f:
    readme = f.read()

setup(name='spotless',
      version='0.4.0',
      description='Grid-Free Deconvolution Directly From Visibilities',
      long_description=readme,
      long_description_content_type="text/markdown",
      url='http://github.com/tmolteno/TART',
      author='Tim Molteno',
      test_suite='nose.collector',
      tests_require=['nose'],
      author_email='tim@elec.ac.nz',
      license='GPLv3',
      install_requires=['numpy', 'matplotlib',
                        'healpy', 'astropy', 'tart', 'disko'],
      packages=['spotless'],
      scripts=['bin/spotless', 'bin/spotless_calibrate', 'bin/gridless'],
      classifiers=[
          "Development Status :: 4 - Beta",
          "Topic :: Scientific/Engineering",
          "Topic :: Communications :: Ham Radio",
          "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
          'Programming Language :: Python :: 3',
          'Programming Language :: Python :: 3.7',
          'Programming Language :: Python :: 3.8',
          'Programming Language :: Python :: 3.9',
          'Programming Language :: Python :: 3.10',
          'Programming Language :: Python :: 3 :: Only',
          "Intended Audience :: Science/Research"])
