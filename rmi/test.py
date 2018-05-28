"""unittest entry point"""

import unittest
from . import test_xdist
from . import test_pp
from . import test_map
from . import test_em


def test_suite():
    """returns test suite with all IRP tests"""
    return unittest.TestSuite([test_xdist.test_suite(),
                               test_pp.test_suite(),
                               test_map.test_suite(),
                               test_em.test_suite()])
