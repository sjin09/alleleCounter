# -*- coding: utf-8 -*-

# DO NOT EDIT THIS FILE!
# This file has been autogenerated by dephell <3
# https://github.com/dephell/dephell

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

import os.path

readme = ''
here = os.path.abspath(os.path.dirname(__file__))
readme_path = os.path.join(here, 'README.rst')
if os.path.exists(readme_path):
    with open(readme_path, 'rb') as stream:
        readme = stream.read().decode('utf8')

setup(
    long_description=readme,
    name='alleleCounter',
    version='0.0.1',
    description='alleleCounter',
    python_requires='==3.*,>=3.6.1',
    project_urls={
        "documentation": "https://alleleCounter.readthedocs.io",
        "homepage": "https://github.com/sjin09/alleleCounter",
        "repository": "https://github.com/sjin09/alleleCounter"
    },
    author='Sangjin Lee',
    author_email='sl17@sanger.ac.uk',
    license='MIT',
    classifiers=[
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],
    entry_points={
        "console_scripts": ["alleleCounter = alleleCounter.__main__:main"]
    },
    packages=['alleleCounter'],
    package_dir={"": "src"},
    package_data={"alleleCounter": ["*.typed"]},
    install_requires=[
        'argparse==1.*,>=1.4.0', 'biopython==1.*,>=1.79.0',
        'click==7.*,>=7.0.0', 'natsort==7.*,>=7.1.1', 'pysam==0.*,>=0.16.0'
    ],
    extras_require={
        "dev": [
            "black==20.*,>=20.8.0.b1", "coverage[toml]==5.*,>=5.3.0",
            "darglint==1.*,>=1.5.8", "flake8==3.*,>=3.8.4",
            "flake8-bandit==2.*,>=2.1.2", "flake8-bugbear==20.*,>=20.1.4",
            "flake8-docstrings==1.*,>=1.5.0",
            "flake8-rst-docstrings==0.*,>=0.0.14", "mypy==0.*,>=0.790.0",
            "pep8-naming==0.*,>=0.11.1", "pre-commit==2.*,>=2.8.2",
            "pre-commit-hooks==3.*,>=3.3.0", "pygments==2.*,>=2.7.2",
            "pytest==6.*,>=6.1.2", "reorder-python-imports==2.*,>=2.3.6",
            "safety==1.*,>=1.9.0", "sphinx==3.*,>=3.3.1",
            "sphinx-autobuild==2020.*,>=2020.9.1", "sphinx-click==2.*,>=2.5.0",
            "sphinx-rtd-theme==0.*,>=0.5.0", "typeguard==2.*,>=2.9.1",
            "xdoctest[colors]==0.*,>=0.15.0"
        ]
    },
)