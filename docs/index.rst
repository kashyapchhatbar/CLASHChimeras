.. CLASHChimeras documentation master file, created by
   sphinx-quickstart on Sun Jun 28 11:30:45 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. contents:: Table of Contents

Introduction
============

CLASHChimeras is a Python_ package for analysing CLASH_ datasets. It takes
raw fastq files as input and provides comprehensive analysis of RNA
profiles and chimeric reads identification. The output is CSV_ and BED_ format
files for easy visualization in Genome Browsers.

Installation
============

You can install it using pip_ after you have setup Python version 3.4 or above.
Please use this guide_ for setting up Python_ if you have not done it already.
After setting up Python_ and pip_, you can run this on your shell

.. code-block:: bash

   $ pip3 install CLASHChimeras

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

.. argparse::
   :ref: clashchimeras.download.parseArguments
   :prog: download-for-chimeras

align-for-chimeras
------------------

.. argparse::
   :ref: clashchimeras.align.parseArguments
   :prog: align-for-chimeras

find-chimeras
-------------

.. argparse::
   :ref: clashchimeras.find.parseArguments
   :prog: find-chimeras

.. _Python: https://www.python.org
.. _CLASH: http://www.nature.com/nprot/journal/v9/n3/abs/nprot.2014.043.html
.. _CSV: https://en.wikipedia.org/wiki/Tab-separated_values
.. _BED: http://www.genome.ucsc.edu/FAQ/FAQformat.html#format1
.. _pip: https://pypi.python.org/pypi/pip
.. _guide: https://docs.python.org/3.4/using/index.html
.. _Bowtie2: http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml
.. _Tophat: http://ccb.jhu.edu/software/tophat/index.shtml