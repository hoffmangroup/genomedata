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

The workflow
============
A genomedata collection contains sequence and may also contain
numerical data associated with that sequence. You can easily load
sequence and numerical data into a genomedata collection with the
:ref:`genomedata-load` command::

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

.. _genomedata-overview:

Overview
========

.. automodule:: genomedata


The genomedata hierarchy:

Each :class:`Genome` contains many :class:`Chromosomes <Chromosome>`
  Each :class:`Chromosome` contains many :class:`Supercontigs <Supercontig>`
    Each :class:`Supercontig` contains one ``continuous`` data set
      Each ``continuous`` data set is a numpy.array of floating
      point numbers with a column for each data track and a row
      for each base in the data set.

Genomic data seldom covers the entire genome but instead tends to be defined
in large but scattered regions. In order to avoid storing the undefined
data between the regions, chromosomes are divided into separate supercontigs
when regions of defined data are far enough apart.

Genomedata Usage
================

Command-line interface
~~~~~~~~~~~~~~~~~~~~~~

Genomedata collections can be created and loaded from the command line
with the :ref:`genomedata-load` command. 

.. _genomedata-load:

genomedata-load
---------------

Usage information follows, but in summary, this script accepts sequence files
(in ``.fa`` or ``.fa.gz`` format) and trackname, datafile pairs
(where trackname is a ``string`` and datafile is in ``.wig`` or ``wig.gz``
format). It outputs a genomedata collection as the specified directory.

::

 Usage: genomedata-load [OPTIONS] GENOMEDATADIR
 e.g. genomedata-load -t high=signal.high -t low=signal.low -s seq.X -s seq.Y outdir

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
creating it if it does not already exist. Gaps of >= 100,000
base pairs (specified as :option:`gap-length`) in the reference sequence,
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
A :option:`chunk-size` can be specified to control the size of hdf5 chunks
(the smallest data read size, like a page size). Larger values of 
:option:`chunk-size` can increase the level of compression, but they also
increase the minimum amount of data that must be read to access a single
value.

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




Python API
~~~~~~~~~~

The genomedata package is designed to be used from a variety of scripting
languages, but currently only exports the following Python API.


.. autoclass:: Genome
   :members:
   :undoc-members:
   
   .. automethod:: __init__
   .. automethod:: __iter__
   .. automethod:: __getitem__


.. autoclass:: Chromosome
   :members:
   :undoc-members:
   
   .. automethod:: __init__
   .. automethod:: __iter__
   .. automethod:: __getitem__


.. autoclass:: Supercontig
   :members:
   :undoc-members:
    
   .. automethod:: __init__

