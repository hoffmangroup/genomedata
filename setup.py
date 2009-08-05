#!/usr/bin/env python

"""genomedata: tools for accessing large amounts of genomic data

LONG_DESCRIPTION
"""

__version__ = "0.1.4"

# Copyright 2008-2009 Michael M. Hoffman <mmh1@washington.edu>

from ez_setup import use_setuptools
use_setuptools()

from setuptools import find_packages, setup

doclines = __doc__.splitlines()
name, short_description = doclines[0].split(": ")
long_description = "\n".join(doclines[2:])

url = "http://noble.gs.washington.edu/~mmh1/software/%s/" % name.lower()
download_url = "%s%s-%s.tar.gz" % (url, name, __version__)

# XXX: remove these when the upstream packages are updated to fix these issues
dependency_links = ["http://pypi.python.org/packages/source/p/path.py/path-2.2.zip"]

classifiers = ["Natural Language :: English",
               "Programming Language :: Python"]

entry_points = """
[console_scripts]
genomedata-load-seq = genomedata.load_seq:main
genomedata-name-tracks = genomedata.name_tracks:main
genomedata-load-data = genomedata.load_data:main
genomedata-save-metadata = genomedata.save_metadata:main
genomedata-report = genomedata.report:main
"""

install_requires = ["numpy", "path", "tables>2.0.4"]

if __name__ == "__main__":
    setup(name=name,
          version=__version__,
          description=short_description,
          author="Michael Hoffman",
          author_email="mmh1@washington.edu",
          url=url,
          download_url=download_url,
          classifiers=classifiers,
          long_description=long_description,
          dependency_links=dependency_links,
          install_requires=install_requires,
          zip_safe=False,

          # XXX: this should be based off of __file__ instead
          packages=find_packages("."),
          include_package_data=True,
          entry_points=entry_points
          )
