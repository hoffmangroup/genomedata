#!/usr/bin/env python

"""genomedata: tools for accessing large amounts of genomic data

LONG_DESCRIPTION
"""

__version__ = "1.0.2"

# Copyright 2008-2009 Michael M. Hoffman <mmh1@washington.edu>

import os

from ez_setup import use_setuptools
use_setuptools()

from distutils.command.build_scripts import build_scripts
from platform import system, processor
from setuptools import find_packages, setup
from shutil import rmtree

doclines = __doc__.splitlines()
name, short_description = doclines[0].split(": ")
long_description = "\n".join(doclines[2:])

url = "http://noble.gs.washington.edu/proj/%s/" % name.lower()
download_url = "%s%s-%s.tar.gz" % (url, name, __version__)

# XXX: remove these when the upstream packages are updated to fix these issues
dependency_links = ["http://pypi.python.org/packages/source/p/path.py/path-2.2.zip"]

classifiers = ["Natural Language :: English",
               "Programming Language :: Python"]

entry_points = """
[console_scripts]
genomedata-load = genomedata.load_genomedata:main
genomedata-load-seq = genomedata._load_seq:main
genomedata-open-data = genomedata._open_data:main
genomedata-close-data = genomedata._close_data:main
genomedata-report = genomedata._report:main
"""
#genomedata-erase-track = genomedata._erase_track:main

install_requires = ["numpy", "path", "tables>2.0.4", "textinput"]

arch = "_".join([system(), processor()])

class DirSet(object):
    """Maintain a set of valid directories.

    add_dir: add the given directory to the set
    add_env: add the given ':'-separated environment variable to the set
    
    """
    def __init__(self):
        self._set = set()
    def add_dir(self, dir):
        if os.path.isdir(dir):
            self._set.add(dir)
    def add_env(self, env):
        if env in os.environ:
            for dir in os.environ[env].split(":"):
                self.add_dir(dir)
    def as_list(self):
        return list(self._set)

# Get compile flags/information from environment
library_dirs = DirSet()
include_dirs = DirSet()

library_dirs.add_env("LIBRARY_PATH")
library_dirs.add_env("LD_LIBRARY_PATH")
include_dirs.add_env("C_INCLUDE_PATH")
if "HDF5_DIR" in os.environ:
    hdf5_dir = os.environ["HDF5_DIR"]
    library_dirs.add_dir(os.path.join(hdf5_dir, "lib"))
    include_dirs.add_dir(os.path.join(hdf5_dir, "include"))

library_dirs = library_dirs.as_list()
include_dirs = include_dirs.as_list()
class BuildScriptWrapper(build_scripts):
    """Override the script-building machinery of distutils.

    Intercepts attempts to build scripts with the `scripts` argument to setup()
    If scripts is a list, this doesn't do anything
    If instead, it is a dict: BIN_NAME -> [SOURCE_FILES], and at least one
    of those source files is a .c file, then this script hops in.

    It extends the distutils compiler to compile the [SOURCE_FILES]
    and then links them into an executable called BIN_NAME, placing
    this executable into the appropriate directory to get added to the
    egg's bin directory.
    """
    def run(self):
        print "##################################################"
        from distutils.ccompiler import new_compiler
        from distutils.sysconfig import customize_compiler

        compiler = new_compiler()
        customize_compiler(compiler)

        # Customize compiler options
        compiler.add_library("hdf5")
        extra_postargs = ["-std=c99"]

        # Remove DNDEBUG flag from all compile statements
        bad_flag = "-DNDEBUG"
        for k, v, in compiler.__dict__.items():
            try:
                v.remove(bad_flag)
            except (AttributeError, TypeError, ValueError):
                pass

        # Compile and link any sources that are passed in
        output_dir = os.path.join(self.build_dir, arch)
        try:
            binaries = []
            for bin, srcs in self.scripts.iteritems():
                # Only compile srcs with c files
                try:
                    if not any([src.endswith(".c") for src in srcs]):
                        continue
                except AttributeError:
                    if not src.endswith(".c"):
                        continue

                objs = compiler.compile(srcs, output_dir=output_dir,
                                        include_dirs=include_dirs,
                                        extra_postargs=extra_postargs,
                                        debug=False)

                bin_path = os.path.join(output_dir, bin)

                compiler.link_executable(objs, bin_path,
                                         library_dirs=library_dirs)
                binaries.append(bin_path)

            # Replace dict script with actual before build_scripts.run() call
            self.scripts = binaries
        except AttributeError:  # Not a dict
            pass

        print "##################################################"

        build_scripts.run(self)  # Call actual script

        # If success, remove script build dir
        if os.path.isdir(output_dir):
            print "Removing script build dir: %s" % output_dir
            rmtree(output_dir)


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
          packages=find_packages(".", exclude="test"),
          entry_points=entry_points,
          scripts={"genomedata-load-data": ["src/genomedata_load_data.c"]},
          cmdclass={"build_scripts": BuildScriptWrapper}
          )
