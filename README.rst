
.. image:: https://zenodo.org/badge/doi/10.5281/zenodo.57076.svg
   :target: https://zenodo.org/record/57076


.. contents:: Table of Contents

Introduction
============

CLASHChimeras is a Python_ package for analysing CLASH_ datasets. It takes
raw fastq files as input and provides comprehensive analysis of RNA
profiles and chimeric reads identification. The output is CSV_ and BED_ format
files for easy visualization in Genome Browsers.

Dependencies
============

CLASHChimeras requires certain software to be installed and setup before you
can use it completely. The software you need to explicitly install are the
following:

* Bowtie2_ - Fast and sensitive read alignment
* Tophat_ - A spliced read mapper for RNA-Seq

Usage
=====

The package can be used by three executable scripts

#. download-for-chimeras_
#. align-for-chimeras_
#. find-chimeras_

download-for-chimeras
---------------------

.. code-block:: bash

   $ download-for-chimeras -gor "H.sapiens" -mor hsa

Downloads required sequences and create bowtie2 indexes required for
alignment

align-for-chimeras
------------------

.. code-block:: bash

   $ align-for-chimeras -i input.fastq -si /path/to/database -ti
   /path/to/database
    -o output

Given a fastq file, this script executes bowtie2 and tophat aligners to generate
alignment files necessary for detecting chimeras in the reads


find-chimeras
-------------

.. code-block:: bash

   $ find-chimeras -s smallRNA.sam -t targetRNA.sam -o output

Given two SAM files, this script tries to find chimeras that
are observed between a smallRNA and a targetRNA

Check out the documentation_ for installation instructions and detailed usage
with example


.. _Python: https://www.python.org
.. _CLASH: http://www.nature.com/nprot/journal/v9/n3/abs/nprot.2014.043.html
.. _CSV: https://en.wikipedia.org/wiki/Tab-separated_values
.. _BED: http://www.genome.ucsc.edu/FAQ/FAQformat.html#format1
.. _pip: https://pypi.python.org/pypi/pip
.. _guide: https://docs.python.org/3.4/using/index.html
.. _Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
.. _Tophat: http://ccb.jhu.edu/software/tophat/index.shtml
.. _documentation: http://clashchimeras.readthedocs.org/en/latest/

