==========================
 Genomedata documentation
==========================

For a broad overview, see Hoffman and Noble, "The genomedata format
for storing large-scale functional genomics data," in preparation.
Please cite that document if you use genomedata.

The workflow
============
A genomedata collection contains sequence and may also contain
numerical data associated with that sequence. You may load the
sequence and numerical data with these commands:

  1. ``genomedata-load-seq``
  2. ``genomedata-name-tracks``
  3. ``genomedata-load-data``
  4. ``genomedata-save-metadata``
