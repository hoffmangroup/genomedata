============
 Genomedata
============
by Michael Hoffman <michael.hoffman at utoronto dot ca>

.. image:: https://img.shields.io/pypi/v/genomedata.png
    :target: https://pypi.python.org/pypi/genomedata/
    :alt: Latest Version

Description
===========
Genomedata is both a format for efficient storage of multiple tracks of
numeric data anchored to a genome and a python interface to genomic datasets.
The file format allows fast random access to hundreds of gigabytes of data, 
while retaining a small disk space footprint. We have also developed utilities
to load data into this format.

Specifically, the genomedata package provides access to genome-scale data,
either using an HDF5_ container or a bigWig_ file.

.. _HDF5: http://www.hdfgroup.org/
.. _bigWig: http://www.genome.ucsc.edu/goldenPath/help/bigWig.html

Please see the following URL (and linked documentation) for information,
installation, and support: https://hoffmanlab.org/proj/genomedata/

Documentation
=============

Live documentation based on this repository can be found on `Read the Docs`_.

.. _Read the Docs: http://genomedata.readthedocs.io/en/latest/

License
========
Genomedata is free software: you can redistribute it and/or modify it under the terms of version 2 of the GNU General Public License as published by the Free Software Foundation.

Genomedata is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.