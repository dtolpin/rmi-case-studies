"""Tests for MAP assignments
"""

import unittest
import rmi.map
from rmi.xdist import xgamma, xexpon
import scipy.stats


class TestMAP(unittest.TestCase):
    """Tests for MAP assignments
    """
    pass

def test_suite():
    """Returns testsuite for posterior probabilities.
    """
    return unittest.TestSuite([TestMAP(test)
                               for test
                               in []])
