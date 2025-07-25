==========================
Cell Maps PPI Embedder
==========================
The Cell Maps PPI Embedder is part of the Cell Mapping Toolkit

.. image:: https://img.shields.io/pypi/v/cellmaps_ppi_embedding.svg
        :target: https://pypi.python.org/pypi/cellmaps_ppi_embedding

.. image:: https://app.travis-ci.com/idekerlab/cellmaps_ppi_embedding.svg?branch=main
        :target: https://app.travis-ci.com/github/idekerlab/cellmaps_ppi_embedding

.. image:: https://readthedocs.org/projects/cellmaps-ppi-embedding/badge/?version=latest
        :target: https://cellmaps-ppi-embedding.readthedocs.io/en/latest/?badge=latest
        :alt: Documentation Status

.. image:: https://zenodo.org/badge/618128210.svg
        :target: https://zenodo.org/doi/10.5281/zenodo.10607492
        :alt: Zenodo DOI badge




Generate embeddings from Protein-Protein interaction networks using `node2vec <https://arxiv.org/abs/1607.00653>`__

* Free software: MIT license
* Documentation: https://cellmaps-ppi-embedding.readthedocs.io.



Dependencies
------------

* `cellmaps_utils <https://pypi.org/project/cellmaps-utils>`__
* `node2vec <https://pypi.org/project/node2vec>`__
* `networkx <https://pypi.org/project/networkx>`__

Compatibility
-------------

* Python 3.8+

Installation
------------

.. code-block::

   git clone https://github.com/idekerlab/cellmaps_ppi_embedding
   cd cellmaps_ppi_embedding
   pip install -r requirements_dev.txt
   make dist
   pip install dist/cellmaps_ppi_embedding*whl


Run **make** command with no arguments to see other build/deploy options including creation of Docker image

.. code-block::

   make

Output:

.. code-block::

   clean                remove all build, test, coverage and Python artifacts
   clean-build          remove build artifacts
   clean-pyc            remove Python file artifacts
   clean-test           remove test and coverage artifacts
   lint                 check style with flake8
   test                 run tests quickly with the default Python
   test-all             run tests on every Python version with tox
   coverage             check code coverage quickly with the default Python
   docs                 generate Sphinx HTML documentation, including API docs
   servedocs            compile the docs watching for changes
   testrelease          package and upload a TEST release
   release              package and upload a release
   dist                 builds source and wheel package
   install              install the package to the active Python's site-packages
   dockerbuild          build docker image and store in local repository
   dockerpush           push image to dockerhub


Before running tests, please install: ``pip install -r requirements_dev.txt``.

For developers
-------------------------------------------

To deploy development versions of this package
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Below are steps to make changes to this code base, deploy, and then run
against those changes.

#. Make changes

   Modify code in this repo as desired

#. Build and deploy

.. code-block::

    # From base directory of this repo cellmaps_ppi_embedding
    pip uninstall cellmaps_ppi_embedding -y ; make clean dist; pip install dist/cellmaps_ppi_embedding*whl



Needed files
------------

The output directory for the PPI download is required (see `Cell Maps PPI Downloader <https://github.com/idekerlab/cellmaps_ppidownloader/>`__).


Usage
-----

For information invoke :code:`cellmaps_ppi_embeddingcmd.py -h`

**Example usage**

.. code-block::

   cellmaps_ppi_embeddingcmd.py ./cellmaps_ppi_embedding_outdir --inputdir ./cellmaps_ppidownloader_outdir


Via Docker
~~~~~~~~~~~~~~~~~~~~~~

**Example usage**


.. code-block::

   Coming soon ... ...

Cite
-------

If you find this tool useful, please cite:

Lenkiewicz, J., Churas, C., Hu, M., Qian, G., Jain, M., Levinson, M. A., ... & Schaffer, L. V. (2025). Cell Mapping Toolkit: An end-to-end pipeline for mapping subcellular organization. Bioinformatics, 41(6), btaf205.

Credits
-------

This package was created with Cookiecutter_ and the `audreyr/cookiecutter-pypackage`_ project template.

.. _Cookiecutter: https://github.com/audreyr/cookiecutter
.. _`audreyr/cookiecutter-pypackage`: https://github.com/audreyr/cookiecutter-pypackage
.. _NDEx: http://www.ndexbio.org
