==================================
Genomedata |version| documentation
==================================
:Website: http://noble.gs.washington.edu/proj/genomedata
:Author: Michael M. Hoffman <mmh1 at washington dot edu>
:Organization: University of Washington
:Address: Department of Genome Sciences, PO Box 355065, Seattle, WA 98195-5065, United States of America
:Copyright: 2009 Michael M. Hoffman

For a broad overview, see the paper:

    Hoffman MM, Buske OJ, Noble WS, "The Genomedata format for storing
    large-scale functional genomics data." Submitted.

Michael <mmh1 at washington dot edu> can send you a copy of the latest
manuscript. Please cite this paper if you use Genomedata.

.. currentmodule:: genomedata

Installation
============

A simple, interactive script_ has been created to install Genomedata
(and most dependencies) on any \*nix platform. Installation is as
simple as downloading and running this script! For instance::

  wget http://noble.gs.washington.edu/proj/genomedata/install.py
  python install.py

.. _script: http://noble.gs.washington.edu/proj/genomedata/install.py

.. note::
  The following are prerequisites:

  - Linux/Unix
      This software has been tested on Linux and Mac OS X systems.
      We would love to add support for other systems in the future and
      will gladly accept any contributions toward this end.
  - Python 2.5 or 2.6
  - Zlib

.. note:: For questions, comments, or troubleshooting, please refer to
          the support_ section.

.. _genomedata-overview:

Overview
========

Genomedata provides a way to store and access large-scale functional
genomics data in a format which is both space-efficient and allows
efficient random-access. Genomedata archives are currently write-once,
although we are working to fix this.

Under the surface, Genomedata is implemented_ as one or more HDF5 files,
but Genomedata provides a transparent interface to interact with your
underlying data without having to worry about the mess of repeatedly parsing
large data files or having to keep them in memory for random access.

.. _implemented: :ref:`Implementation`

The Genomedata hierarchy:

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


Implementation
==============
Genomedata archives are implemented as one or more HDF5 files.
The :ref:`API <python-api>` handles both single-file
and directory archives transparently, but the implementation options
exist for several performance reasons.

Use a **directory** with few chromosomes/scaffolds:
  * Parallel load/access
  * Smaller file sizes

Use a **single file** with many chromosomes/scaffolds:
  * More efficient access with many chromosomes/scaffolds
  * Easier archive distribution

Implementing the archive as a directory makes it easier to
parallelize access to the data. In particular, it makes it easy to
create the archives in parallel with one chromosome on each
machine. It also reduces the likelihood of running into the
2 GB file limit applicable to older applications and older versions
of 32-bit UNIX. We are currently using an 81-track Genomedata
archive for our research which has a total size of 18 GB, but the
largest single file (chr1) is only 1.6 GB.

A directory-based Genomedata archive is not ideal for all circumstances,
however, such as when working with genomes with many chromosomes, contigs,
or scaffolds. In these situations, a single file implementation would be much
more efficient. Additionally, having the archive as a single file allows
the archive to be distributed much more easily (without tar/zip/etc).

.. note:: The default behavior is to implement the Genomedata archive as a
          directory if there are fewer than 100 sequences being loaded and as
          a single file otherwise.

.. versionadded:: 1.1
   Single-file-based Genomedata archives

Creation
========
A Genomedata archive contains sequence and may also contain
numerical data associated with that sequence. You can easily load
sequence and numerical data into a Genomedata archive with the
:ref:`genomedata-load` command (see command details additional details)::

    genomedata-load [-t trackname=signalfile]... [-s sequencefile]... GENOMEDATAFILE


This command is a user-friendly shortcut to the typical workflow.
The underlying commands are still installed and may be used if more
fine-grained control is required (for instance, parallel data loading).
The commands and required ordering are:

1. :ref:`genomedata-load-seq`
#. :ref:`genomedata-open-data`
#. :ref:`genomedata-load-data`
#. :ref:`genomedata-close-data`

Entire data tracks can be replaced with the following pipeline:

1. :ref:`genomedata-erase-data`
#. :ref:`genomedata-load-data`
#. :ref:`genomedata-close-data`

.. versionadded:: 1.1
   The ability to replace data tracks.

.. note:: A call to :program:`h5repack` after
          :ref:`genomedata-close-data` may be used to
          transparently compress the data.

.. _genomedata-load-example:

Example
~~~~~~~
The following is a brief example for creating a Genomedata archive from
sequence and signal files.

Given the following two sequence files:

1. chr1.fa::

     >chr1
     taaccctaaccctaaccctaaccctaaccctaaccctaaccctaacccta
     accctaaccctaaccctaaccctaaccct

#. chrY.fa.gz::

     >chrY
     ctaaccctaaccctaaccctaaccctaaccctaaccctCTGaaagtggac

and the following two signal files:

1. signal_low.wigFix::

     fixedStep chrom=chr1 start=5 step=1
     0.372
     -2.540
     0.371
     -2.611
     0.372
     -2.320

#. signal_high.bed.gz::

     chrY    0       12      4.67
     chrY    20      23      9.24
     chr1    1       3       2.71
     chr1    3       6       1.61
     chr1    6       24      3.14


A Genomedata archive (``genomedata.test``) could then be created with the
following command::

    genomedata-load -s chr1.fa -s chrY.fa.gz -t low=signal_low.wigFix -t high=signal_high.bed.gz genomedata.test

or the following pipeline::

   genomedata-load-seq genomedata.test chr1.fa chrY.fa.gz
   genomedata-open-data genomedata.test low high
   genomedata-load-data genomedata.test low < signal_low.wigFix
   zcat signal_high.bed.gz | genomedata-load-data genomedata.test high
   genomedata-close-data genomedata.test

.. note:: chr1.fa and chrY.fa.gz could als be combined into a single
          sequence file with two sequences

.. warning::
   It is important that the sequence names (`chrY`, `chr1`) in the signal files
   match the sequence identifiers in the sequence files exactly.

Genomedata usage
================

Python interface
~~~~~~~~~~~~~~~~
The data in Genomedata is accessed through the hierarchy described in
:ref:`genomedata-overview`. A full :ref:`Python API <python-api>` is
also available. To appreciate the full benefit of Genomedata,
it is most easily used as a contextmanager::

    from genomedata import Genome
    [...]
    gdfilename = "/path/to/genomedata/archive"
    with Genome(gdfilename) as genome:
        [...]

.. note::
    If Genome is used as a context manager, it will clean up any opened
    Chromosomes automatically.
    If not, the Genome object (and all opened chromosomes) should be
    closed manually with a call to :meth:`Genome.close`.

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

Genomedata archives can be created and loaded from the command line
with the :ref:`genomedata-load` command.

.. _genomedata-load:

genomedata-load
---------------

Usage information follows, but in summary, this script takes as input:

- sequence files in |sequence file formats| format, where the sequence
  identifiers are the names of the chromosomes/scaffolds to create.
- trackname, datafile pairs (specified as ``trackname=datafile``), where:
    * trackname is a ``string`` identifier (e.g. ``broad.h3k27me3``)
    * datafile contains signal data for this data track
      in one of the following formats: |signal file formats|
    * the chromosomes/scaffolds referred to in the datafile MUST be identical
      to those found in the sequence files
- the name of the Genomedata archive to create

See the :ref:`full example <genomedata-load-example>` for more details.

.. |signal file formats| replace:: |signal data formats|, or a gzip'd
                         form of any of the preceding

.. |sequence file formats| replace:: FASTA_ (``.fa`` or ``.fa.gz``)

.. _FASTA: http://www.ncbi.nlm.nih.gov/blast/fasta.shtml

Command-line usage information::

 Usage: genomedata-load [OPTIONS] GENOMEDATAFILE

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
This can be especially useful for parallelizing Genomedata loading over a
cluster.


.. _genomedata-load-seq:

genomedata-load-seq
-------------------

This command adds the provided sequence files to the specified Genomedata,
archive creating it if it does not already exist. Sequence files should be in
|sequence file formats| format. Gaps of >= 100,000 base pairs
(specified as :option:`gap-length`) in the reference sequence,
are used to divide the sequence into supercontigs. The FASTA sequence
identifier will be used as the name for the chromosomes/scaffolds
created within the Genomedata archive and must be consistent between these
sequence files and the data loaded later with :ref:`genomedata-load-data`.
See :ref:`this example <genomedata-load-example>` for details.

::

 Usage: genomedata-load-seq [OPTION]... GENOMEDATAFILE SEQFILE...

 Options:
   -g, --gap-length  XXX: Implement this.
   --version         show program's version number and exit
   -h, --help        show this help message and exit


.. _genomedata-open-data:

genomedata-open-data
--------------------

This command opens the specified tracknames in the Genomedata archive,
allowing data for those tracks to be added with :ref:`genomedata-load-data`.

::

 Usage: genomedata-open-data [OPTION]... GENOMEDATAFILE TRACKNAME...

 Options:
   --version   show program's version number and exit
   -h, --help  show this help message and exit


.. _genomedata-load-data:

genomedata-load-data
--------------------

This command loads data from stdin into Genomedata under the given trackname.
The input data must be in one of these supported datatypes:
|signal data formats|. The chromosome/scaffold references in these files must
match the sequence identifiers in the sequence files loaded with
:ref:`genomedata-load-seq`. See :ref:`this example <genomedata-load-example>`
for details. A :option:`chunk-size` can be specified to control the
size of hdf5 chunks (the smallest data read size, like a page size).
Larger values of :option:`chunk-size` can increase the level
of compression, but they also increase the minimum amount of data
that must be read to access a single value.

.. |signal data formats| replace:: WIG_, BED_, bedGraph_

.. _WIG: http://genome.ucsc.edu/FAQ/FAQformat#format6
.. _BED: http://genome.ucsc.edu/FAQ/FAQformat#format1
.. _bedGraph: http://genome.ucsc.edu/goldenPath/help/bedgraph.html

::

 Usage: genomedata-load-data [OPTION...] GENOMEDATAFILE TRACKNAME
 Loads data into Genomedata format
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

Closes the specified Genomedata arhive.

::

 Usage: genomedata-close-data [OPTION]... GENOMEDATAFILE

 Options:
   --version   show program's version number and exit
   -h, --help  show this help message and exit


.. _genomedata-erase-data:

genomedata-erase-data
---------------------

Erases all data associated with the specified tracks, allowing the data to
then be replaced. The pipeline for replacing a data track is:

1. :ref:`genomedata-erase-data`
#. :ref:`genomedata-load-data`
#. :ref:`genomedata-close-data`

::

 Usage: genomedata-erase-data [OPTION]... GENOMEDATAFILE TRACKNAME...

 Erase the specified tracks from the Genomedata archive in such a way that
 the track can be replaced (via genomedata-load-data).

 Options:
   --version      show program's version number and exit
   -h, --help     show this help message and exit
   -v, --verbose  Print status updates and diagnostic messages

.. _python-api:

Python API
~~~~~~~~~~

The Genomedata package is designed to be used from a variety of scripting
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


.. _support:

Support
=======

To stay informed of **new releases**, subscribe to the moderated 
``genomedata-announce`` mailing list (mail volume very low):

  https://mailman1.u.washington.edu/mailman/listinfo/genomedata-announce

For **discussion and questions** about the use of the Genomedata system,
there is a ``genomedata-users`` mailing list:
  
  https://mailman1.u.washington.edu/mailman/listinfo/genomedata-users

For issues related to the use of Genomedata on **Mac OS X**, 
please use the above mailing list or contact 
Jay Hesselberth <jay dot hesselberth at ucdenver dot edu>.

If you want to **report a bug or request a feature**, please do so using
our issue tracker:

  http://code.google.com/p/genomedata/issues

For other support with Genomedata, or to provide feedback, please write
contact the authors directly. We are interested in all comments regarding the
package and the ease of use of installation and documentation.

