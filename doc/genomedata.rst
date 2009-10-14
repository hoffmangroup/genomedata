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
:program:`genomedata-load` command::

    genomedata-load [-t trackname=signalfile]... [-s sequencefile]... GENOMEDATADIR


This command is a user-friendly shortcut to the typical workflow.
The underlying commands are still installed and may be used if more
fine-grained control is required. The commands and required ordering are:

  1. :program:`genomedata-load-seq`
  #. :program:`genomedata-open-data`
  #. :program:`genomedata-load-data`
  #. :program:`genomedata-close-data`

.. note:: A call to :program:`h5repack` after
          :program:`genomedata-close-data` may be used to
          transparently compress the data.

Overview
========

.. automodule:: genomedata


The important objects of the genomedata hierarchy are:

The Genome_ object
    Each ``genome`` contains many ``chromosomes``
The Chromosome_ object
    Each ``chromosome`` contains many ``supercontigs``
The Supercontig_ object
    Each ``supercontig`` contains one ``continuous`` data set
The continuous data set
    Each continuous data set is a numpy.array of floating point numbers
    with a column for each data track and a row for each base in the data set.


Genomedata Usage
================

Command-line interface
~~~~~~~~~~~~~~~~~~~~~~

Genomedata collections can be created and loaded from the command line
with the :program:`genomedata-load` command. 

Usage information follows, but in summary, it accepts sequence files
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


Python API
~~~~~~~~~~

The genomedata package is designed to be used from a variety of scripting
languages, but currently only exports the following Python API.

Genome
------

.. autoclass:: Genome
   :members:
   :undoc-members:
   
   .. automethod:: __init__
   .. automethod:: __iter__
   .. automethod:: __getitem__

Chromosome
----------

.. autoclass:: Chromosome
   :members:
   :undoc-members:
   
   .. automethod:: __init__
   .. automethod:: __iter__
   .. automethod:: __getitem__

Supercontig
-----------

.. autoclass:: Supercontig
   :members:
   :undoc-members:
    
   .. automethod:: __init__

