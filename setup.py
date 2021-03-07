# Copyright (c) 2020, Franck Lejzerowicz.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file LICENSE, distributed with this software.
# ----------------------------------------------------------------------------

import re, ast
from setuptools import find_packages, setup

classes = """
    License :: OSI Approved :: BSD License
    Topic :: Scientific/Engineering
    Topic :: Scientific/Engineering :: Bio-Informatics
    Programming Language :: Python :: 3.5
    Programming Language :: Python :: 3 :: Only
    Operating System :: Unix
    Operating System :: POSIX
    Operating System :: MacOS :: MacOS X
"""

classifiers = [s.strip() for s in classes.split('\n') if s]

description = (
    ""
)

with open("README.md") as f:
    long_description = f.read()

_version_re = re.compile(r"__version__\s+=\s+(.*)")

with open("Xrbfetch/__init__.py", "rb") as f:
    hit = _version_re.search(f.read().decode("utf-8")).group(1)
    version = str(ast.literal_eval(hit))

standalone = ['Xrbfetch=Xrbfetch.scripts._standalone_xrbfetch:standalone_xrbfetch']

setup(
    name="Xrbfetch",
    version=version,
    license="BSD",
    description=description,
    long_description=long_description,
    long_description_content_type="text/markdown",
    author="Franck Lejzerowicz",
    author_email="franck.lejzerowicz@gmail.com",
    maintainer="Franck Lejzerowicz",
    maintainer_email="franck.lejzerowicz@gmail.com",
    url="https://github.com/FranckLejzerowicz/Xrbfetch",
    packages=find_packages(),
    install_requires=[
        "click >= 6.7",
        'pandas >= 0.19.0',
        'numpy >= 1.12.1',
        'cython >= 0.29.15',
        'redbiom >= 0.3.5',
        'biom-format >= 2.1.5'
    ],
    classifiers=classifiers,
    entry_points={'console_scripts': standalone},
    package_data={'Xrbfetch': ['resources/newblooms.all.fasta']},
    python_requires='>=3.5',
)
