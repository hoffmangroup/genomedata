#!/usr/bin/env python
from __future__ import division, with_statement

import os
import re
import sys
import unittest

from pkg_resources import resource_listdir

def main(verbose=True):
    # Move to test directory to allow imports
    os.chdir(os.path.dirname(__file__))

    # Gather a list of unittest modules
    filenames = os.listdir(os.getcwd())
    regex = re.compile("^test_.*\.py$", re.IGNORECASE)
    module_filenames = filter(regex.search, filenames)
    make_module_name = lambda filename: filename[:-3]
    modulenames = map(make_module_name, module_filenames)
    if verbose:
        print "Found test modules: %r" % modulenames

    map(__import__, modulenames)  # Import them all
    modules = [sys.modules[modulename] for modulename in modulenames]

    # Run the test suite for each
    suite = unittest.TestSuite([module.suite() for module in modules])
    if verbose:
        verbosity = 2
    else:
        verbosity = 1

    unittest.TextTestRunner(verbosity=verbosity).run(suite)

if __name__ == "__main__":
    sys.exit(main())
