#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    from setuptools import setup

except ImportError:
    from distutils.core import setup

setup(
    name='exocartographer',
    version='1.0',
    description='Map some exoplanets',
    author='Nick Cowan, Ben Farr',
    author_email='ncowan@amhurst.edu, farr@uchicago.edu',
    url='https://github.com/bfarr/exocartographer',
    include_package_data=True,
    packages=['exocartographer'],
    install_requires=['numpy', 'scipy'],
)