"""Install RF-tools"""

import io
import os
import sys
from os import path

from setuptools import find_packages, setup
from setuptools.command.test import test as TestCommand

import rftools

root = path.abspath(path.dirname(__file__))

def read(*filenames, **kwargs):
    encoding = kwargs.get('encoding', 'utf-8')
    sep = kwargs.get('sep', '\n')
    buf = []
    for filename in filenames:
        with io.open(path.join(root, filename), encoding=encoding) as f:
            buf.append(f.read())
    return sep.join(buf)

long_description = read('README.md')

class PyTest(TestCommand):
    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        import pytest
        errcode = pytest.main(self.test_args)
        sys.exit(errcode)

setup(
    name = "rftools",
    version = rftools.__version__,
    author = "John Garrett",
    author_email = "garrettj403@gmail.com",
    description = "Tools to design RF components and networks.",
    license = "MIT",
    url = "https://github.com/garrettj403/RF-tools/",
    keywords = [
        "RF engineering",
        "Microwave engineering",
        "Electrical engineering",
        "RF components",
        "Touchstone",
        "HFSS",
    ],
    packages=find_packages(),
    install_requires=[
        'matplotlib',
        'numpy',
        'scipy',
    ],
    extras_require={'testing': ['pytest'],},
    tests_require=['pytest'],
    cmdclass={'test': PyTest},
    long_description=long_description,
    long_description_content_type='text/markdown',
    platforms='any',
    classifiers=[
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Topic :: Scientific/Engineering :: Physics",
    ],
    project_urls={
        'Changelog': 'https://github.com/garrettj403/RF-tools/blob/master/CHANGES.md',
        'Issue Tracker': 'https://github.com/garrettj403/RF-tools/issues',
    },
    scripts=[
        "bin/50ohm-line",
        "bin/waveguide",
        "bin/waveguide-att",
        "bin/cwaveguide",
        "bin/noisetemp",
        "bin/temp-rj",
        "bin/temp-cw",
    ],
)
