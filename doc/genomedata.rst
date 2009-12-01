========================
Genomedata documentation
========================
:Author: Michael M. Hoffman <mmh1 at washington dot edu>
:Organization: University of Washington
:Address: Department of Genome Sciences, PO Box 355065, Seattle, WA 98195-5065, United States of America
:Copyright: 2009 Michael M. Hoffman

For a broad overview, see the paper:

    Hoffman MM, Buske OJ, Noble WS, "The genomedata format for storing
    large-scale functional genomics data." In preparation.

Michael <mmh1 at washington dot edu> can send you a copy of the latest
manuscript. Please cite this paper if you use genomedata.

.. currentmodule:: genomedata

Installation
============

A simple, interactive script_ has been created to install genomedata
(and most dependencies) on any Linux platform. Installation is as
simple as downloading and running this script! For instance::

  wget http://noble.gs.washington.edu/proj/genomedata/install.py
  python install.py

.. _script: http://noble.gs.washington.edu/proj/genomedata/install.py

.. note:: 
  The following are prerequisites:

  - Python 2.X | X >= 5
  - Zlib

.. note:: Please send any install script bugs/issues/comments to:
          Orion Buske <stasis at uw dot edu>


.. _genomedata-overview:

Overview
========

Genomedata is a module to store and access large-scale functional
genomics data in a format which is both space-efficient and allows
efficient random-access.

Under the surface, genomedata is implemented as a collection of HDF5 files,
but genomedata provides a transparent interface to interact with your
underlying data without having to worry about the mess of repeatedly parsing
large data files or having to keep them in memory for random access.

The genomedata hierarchy:

  Each :class:`Genome` contains many :class:`Chromosomes <Chromosome>`
    Each :class:`Chromosome` contains many :class:`Supercontigs <Supercontig>`
      Each :class:`Supercontig` contains one ``continuous`` data set
        Each ``continuous`` data set is a numpy.array of floating
        point numbers with a column for each data track and a row
        for each base in the data set.

Why have :class:`Supercontigs <Supercontig>`?
  Genomic data seldom covers the entire genome but instead tends to be defined
  in large but scattered regions. In order to avoid storing the undefined
  data between the regions, chromosomes are divided into separate supercontigs
  when regions of defined data are far enough apart. They also serve as
  a convenient chunk since they can usually fit entirely in memory.


The workflow
============
A genomedata collection contains sequence and may also contain
numerical data associated with that sequence. You can easily load
sequence and numerical data into a genomedata collection with the
:ref:`genomedata-load` command (see command details additional details)::

    genomedata-load [-t trackname=signalfile]... [-s sequencefile]... GENOMEDATADIR


This command is a user-friendly shortcut to the typical workflow.
The underlying commands are still installed and may be used if more
fine-grained control is required. The commands and required ordering are:

1. :ref:`genomedata-load-seq`
#. :ref:`genomedata-open-data`
#. :ref:`genomedata-load-data`
#. :ref:`genomedata-close-data`

.. note:: A call to :program:`h5repack` after
          :ref:`genomedata-close-data` may be used to
          transparently compress the data.


Genomedata usage
================

Python interface
~~~~~~~~~~~~~~~~

The data in genomedata is accessed through the hierarchy described in
:ref:`genomedata-overview`. A full :ref:`Python API <python-api>` is 
also available. To appreciate the full benefit of genomedata,
it is most easily used as a contextmanager::

    from genomedata import Genome
    [...]
    genomedatadir = "/path/to/genomedata"
    with Genome(genomedatadir) as genome:
        [...]

.. note:: 
    If used as a context manager, chromosome access is memoized. 
    If not, chromosomes should be closed manually with 
    :meth:`Chromosome.close`.

Basic usage
-----------

Genomedata is designed to make it easy to get to the data you want.
Here are a few examples:

**Get arbitrary sequence** (10-bp sequence starting at chr2:1423):

>>> chromosome = genome["chr2"]
>>> seq = chromosome.seq[1423:1433]
>>> seq
array([116,  99,  99,  99,  99, 103, 103, 103, 103, 103], dtype=uint8)
>>> seq.tostring()
'tccccggggg'

**Get arbitrary data** (data from first 3 tracks for region chr8:999-1000):

>>> chromosome = genome["chr8"]
>>> chromosome[999:1001, 0:3]  # Note the half-open, zero-based indexing
array([[ NaN,  NaN,  NaN],
       [ 3. ,  5.5,  3.5], dtype=float32)

**Get data for a specific track** (specified data in first 5-bp of chr1):

>>> chromosome = genome["chr1"]
>>> data = chromosome[0:5, "sample_track"]
>>> data
array([ 47.,  NaN,  NaN,  NaN,  NaN], dtype=float32)

*Only specified data*:

>>> from numpy import isfinite
>>> data[isfinite(data)]
array([ 47.], dtype=float32)
    
.. note:: Specify a slice for the track to keep the data in column form:

          >>> col_index = chromosome.index_continuous("sample_track")
          >>> data = chromosome[0:5, col_index:col_index+1]


Command-line interface
~~~~~~~~~~~~~~~~~~~~~~

Genomedata collections can be created and loaded from the command line
with the :ref:`genomedata-load` command. 

.. _genomedata-load:

genomedata-load
---------------

Usage information follows, but in summary, this script takes as input:

- sequence files in |sequence file formats| format
- trackname, datafile pairs (specified as ``trackname=datafile``), where:
    * trackname is a ``string`` identifier (e.g. ``broad.h3k27me3``) 
    * datafile contains one column of data for this data track 
      in one of the following formats: |signal file formats|
- the name of the genomedata collection to create

For example, let's say you have sequence data for chrX (``chrX.fa``) and
chrY (``chrY.fa.gz``), as well as two signal tracks: high (``signal.high.wig``)
and low (``signal.low.bed.gz``). You could construct a genomedata collection
named ``mygenomedata`` in the current directory with the following command::

    genomedata-load -s chrX.fa -s chrY.fa.gz -t high=signal.high.wig -t low=signal.low.bed.gz mygenomedata

.. |signal file formats| replace:: |signal data formats|, or a gzip'd 
                         form of any of the preceding

.. |sequence file formats| replace:: FASTA_ (``.fa`` or ``.fa.gz``)

.. _FASTA: http://www.ncbi.nlm.nih.gov/blast/fasta.shtml

Command-line usage information::

 Usage: genomedata-load [OPTIONS] GENOMEDATADIR

 --track and --sequence may be repeated to specify multiple trackname=trackfile
 pairings and sequence files, respectively

 Options:
   --version             show program's version number and exit
   -h, --help            show this help message and exit
   -s SEQFILE, --sequence=SEQFILE
                         Add the sequence data in the specified file
   -t TRACK, --track=TRACK
                         Add data for the given track. TRACK should be
                         specified in the form: NAME=FILE, such as: -t
                         signal=signal.dat


Alternately, as described in :ref:`genomedata-overview`, the underlying
Python and C load scripts are also accessible for more finely-grained control.
This can be especially useful for parallelizing genomedata loading over a
cluster. 


.. _genomedata-load-seq:

genomedata-load-seq
-------------------

This command adds the provided sequence files to the specified genomedatadir,
creating it if it does not already exist. Sequence files should be in
|sequence file formats| format. Gaps of >= 100,000 base pairs 
(specified as :option:`gap-length`) in the reference sequence,
are used to divide the sequence into supercontigs.

::

 Usage: genomedata-load-seq [OPTION]... GENOMEDATADIR SEQFILE...
 
 Options:
   -g, --gap-length  XXX: Implement this.
   --version         show program's version number and exit
   -h, --help        show this help message and exit


.. _genomedata-open-data:

genomedata-open-data
--------------------

This command opens the specified tracknames in the genomedata object,
allowing data for those tracks to be added with :ref:`genomedata-load-data`.

::

 Usage: genomedata-open-data [OPTION]... GENOMEDATADIR TRACKNAME...
 
 Options:
   --version   show program's version number and exit
   -h, --help  show this help message and exit


.. _genomedata-load-data:

genomedata-load-data
--------------------

This command loads data from stdin into genomedata under the given trackname.
The input data must be in one of these supported datatypes: 
|signal data formats|.
A :option:`chunk-size` can be specified to control the size of hdf5 chunks
(the smallest data read size, like a page size). Larger values of 
:option:`chunk-size` can increase the level of compression, but they also
increase the minimum amount of data that must be read to access a single
value.

.. |signal data formats| replace:: WIG_, BED_, bedGraph_

.. _WIG: http://genome.ucsc.edu/FAQ/FAQformat#format6
.. _BED: http://genome.ucsc.edu/FAQ/FAQformat#format1
.. _bedGraph: http://genome.ucsc.edu/goldenPath/help/bedgraph.html

::

 Usage: genomedata-load-data [OPTION...] GENOMEDATADIR TRACKNAME
 Loads data into genomedata format
 Takes track data in on stdin
 
   -c, --chunk-size=NROWS     Chunk hdf5 data into blocks of NROWS. A higher
                              value increases compression but slows random
                              access. Must always be smaller than the max size
                              for a dataset. [default: 10000]
   -?, --help                 Give this help list
       --usage                Give a short usage message
   -V, --version              Print program version
 
 Mandatory or optional arguments to long options are also mandatory or optional
 for any corresponding short options.


.. _genomedata-close-data:

genomedata-close-data
---------------------

Closes the specified genomedata object.

::

 Usage: genomedata-close-data [OPTION]... GENOMEDATADIR
 
 Options:
   --version   show program's version number and exit
   -h, --help  show this help message and exit


.. _python-api:


Python API
~~~~~~~~~~

The genomedata package is designed to be used from a variety of scripting
languages, but currently only exports the following Python API.

.. module:: genomedata

.. autoclass:: Genome
   :members:
   :undoc-members:
   
   .. automethod:: __init__
   .. automethod:: __iter__
   .. automethod:: __getitem__

.. autoclass:: Chromosome
   :members:
   :undoc-members:
   
   .. automethod:: __iter__
   .. automethod:: __getitem__

.. autoclass:: Supercontig
   :members:
   :undoc-members:


Support
=======

There is a moderated genomedata-announce mailing list that you can 
subscribe to for information on new releases of genomedata:
    
https://mailman1.u.washington.edu/mailman/listinfo/genomedata-announce

If you want to report a bug or request a feature, please do so using the
Genomedata issue tracker:

http://code.google.com/p/genomedata/issues

For other support with Genomedata, or to provide feedback, please write
contact the authors directly. We are interested in all comments regarding the 
package and the ease of use of installation and documentation.

