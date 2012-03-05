#!/usr/bin/env python

"""genomedata: tools for accessing large amounts of genomic data

Genomedata is a format for efficient storage of multiple tracks of
numeric data anchored to a genome. The format allows fast random
access to hundreds of gigabytes of data, while retaining a small disk
space footprint. We have also developed utilities to load data into
this format.
"""

__version__ = "1.3.2"

# Copyright 2008-2011 Michael M. Hoffman <mmh1@washington.edu>

import os
import sys

from ez_setup import use_setuptools
use_setuptools()

from distutils.command.clean import clean
from distutils.command.build_scripts import build_scripts
from distutils.spawn import find_executable
from platform import system, processor
from setuptools import find_packages, setup
from shutil import rmtree
from subprocess import CalledProcessError, check_call

doclines = __doc__.splitlines()
name, short_description = doclines[0].split(": ")
long_description = "\n".join(doclines[2:])

url = "http://noble.gs.washington.edu/proj/%s/" % name.lower()
download_url = "%s%s-%s.tar.gz" % (url, name, __version__)

classifiers = ["Natural Language :: English",
               "Programming Language :: Python"]

entry_points = """
[console_scripts]
genomedata-info = genomedata._info:main
genomedata-load = genomedata.load_genomedata:main
genomedata-load-seq = genomedata._load_seq:main
genomedata-load-assembly = genomedata._load_seq:main
genomedata-open-data = genomedata._open_data:main
genomedata-close-data = genomedata._close_data:main
genomedata-report = genomedata._report:main
genomedata-erase-data = genomedata._erase_data:main
genomedata-test = test.run_tests:main
"""

install_requires = ["numpy", "forked-path", "tables>=2.2", "textinput"]

arch = "_".join([system(), processor()])

include_gnulib = (system() != "Linux")
GNULIB_BUILD_DIR = "src/build-deps"
GNULIB_LIB_DIR = "%s/gllib" % GNULIB_BUILD_DIR

class DirList(list):
    """Maintain a unique list of valid directories.

    This is a list, not a set to maintain order.

    add_dir: add the given directory to the set
    add_env: add the given ':'-separated environment variable to the set
    """
    def add_dir(self, dirname):
        if dirname not in self and os.path.isdir(dirname):
            self.append(dirname)

    def add_env(self, env):
        if env in os.environ:
            for dirname in os.environ[env].split(":"):
                self.add_dir(dirname)

    def append(self, item):
        if item not in self:
            list.append(self, item)

# Get compile flags/information from environment
library_dirnames = DirList()
include_dirnames = DirList()

if "HDF5_DIR" in os.environ:
    hdf5_dir = os.environ["HDF5_DIR"]
    library_dirnames.add_dir(os.path.join(hdf5_dir, "lib"))
    include_dirnames.add_dir(os.path.join(hdf5_dir, "include"))

if include_gnulib:
    # Gnulib for OS X dependencies
    library_dirnames.add_dir(GNULIB_LIB_DIR)
    include_dirnames.add_dir(GNULIB_LIB_DIR)

library_dirnames.add_env("LIBRARY_PATH")
library_dirnames.add_env("LD_LIBRARY_PATH")
include_dirnames.add_env("C_INCLUDE_PATH")

## fix types, since distutils does type-sniffing:
library_dirnames = list(library_dirnames)
include_dirnames = list(include_dirnames)

class InstallationError(Exception):
    pass

class CleanWrapper(clean):
    """Wraps `python setup.py clean` to also cleans Gnulib installation"""
    def run(self):
        clean.run(self)
        if include_gnulib:
            print >>sys.stderr, ">> Cleaning Gnulib build directory"
            try:
                check_call(["make", "clean"], cwd=GNULIB_BUILD_DIR)
            except CalledProcessError:
                print >>sys.stderr, ">> WARNING: Failed to clean Gnulib build!"

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

        # these two are necessary if a static HDF5 is installed:
        compiler.add_library("m")
        compiler.add_library("z")

        # is a static HDF5 is installed with --with-szip?
        h5dump_filename = find_executable("h5dump")
        try:
            # XXX: output should be redirected to /dev/null
            check_call("objdump -t %s | fgrep szip" % h5dump_filename, shell=True)
        except CalledProcessError:
            pass
        else:
            compiler.add_library("sz")

        if include_gnulib:
            compiler.add_library("gnu")

        extra_postargs = ["-std=c99", "-Wall"]

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
                                        include_dirs=include_dirnames,
                                        extra_postargs=extra_postargs,
                                        debug=False)

                bin_path = os.path.join(output_dir, bin)

                compiler.link_executable(objs, bin_path,
                                         library_dirs=library_dirnames)
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

def make_gnulib():
    print ">> Not a Linux system: configuring and making Gnulib libraries..."

    libfilename = "%s/libgnu.a" % GNULIB_LIB_DIR
    if os.path.isfile(libfilename):
        print ">> Found libgnu.a... skipping configure and make"
    else:
        commands = ["./configure", "make"]
        for command in commands:
            print ">> %s" % command
            try:
                check_call(command, cwd=GNULIB_BUILD_DIR)
            except CalledProcessError:
                raise InstallationError("Error compiling Gnulib")

    if not os.path.isfile(libfilename):
        raise InstallationError("Expected to find: %s" % libfilename)

    to_rm = ["%s/getopt.h" % GNULIB_LIB_DIR]
    for filename in to_rm:
        if os.path.isfile(filename):
            print ">> Removing: %s" % filename
            os.remove(filename)

    print ">> Gnulib libraries successfully created!"

if __name__ == "__main__":
    # Configure and make gnulib if not on Linux
    if include_gnulib:
        try:
            make_gnulib()
        except InstallationError, e:
            print >>sys.stderr, ">> ERROR: %s" % e
            sys.exit(1)

    setup(name=name,
          version=__version__,
          description=short_description,
          author="Michael Hoffman",
          author_email="mmh1@washington.edu",
          url=url,
          download_url=download_url,
          classifiers=classifiers,
          long_description=long_description,
          install_requires=install_requires,
          zip_safe=False,
          # XXX: this should be based off of __file__ instead
          packages=find_packages("."),  # including "test"
          include_package_data = True,
          entry_points=entry_points,
          scripts={"genomedata-load-data": ["src/genomedata_load_data.c"]},
          cmdclass={"build_scripts": BuildScriptWrapper,
                    "clean": CleanWrapper}
          )
