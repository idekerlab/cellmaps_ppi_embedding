#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""
import os
import re
from setuptools import setup, find_packages


with open(os.path.join('cellmaps_ppi_embedding', '__init__.py')) as ver_file:
    for line in ver_file:
        line = line.rstrip()
        if line.startswith('__version__'):
            version=re.sub("'", "", line[line.index("'"):])
        elif line.startswith('__description__'):
            desc = re.sub("'", "", line[line.index("'"):])
        elif line.startswith('__repo_url__'):
            repo_url = re.sub("'", "", line[line.index("'"):])
        elif line.startswith('__author__'):
            author = re.sub("'", "", line[line.index("'"):])
        elif line.startswith('__email__'):
            email = re.sub("'", "", line[line.index("'"):])

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = ['scipy>=1.10.1,<1.13',
                'cellmaps_utils>=0.4.0,<1.0.0',
                'node2vec>=0.4.6,<0.5.0',
                'networkx>=2.8,<2.9']

setup_requirements = [ ]

setup(
    author=author,
    author_email=email,
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    description=desc,
    install_requires=requirements,
    license="MIT license",
    long_description=readme + '\n\n' + history,
    long_description_content_type='text/x-rst',
    include_package_data=True,
    keywords='cellmaps_ppi_embedding',
    name='cellmaps_ppi_embedding',
    packages=find_packages(include=['cellmaps_ppi_embedding']),
    package_dir={'cellmaps_ppi_embedding': 'cellmaps_ppi_embedding'},
    package_data={'cellmaps_ppi_embedding': ['readme_outputs.txt']},
    scripts=['cellmaps_ppi_embedding/cellmaps_ppi_embeddingcmd.py'],
    setup_requires=setup_requirements,
    url=repo_url,
    version=version,
    zip_safe=False)
