"""Tests for extended distribution objects.
"""

import unittest
import rmi.xdist
import scipy.stats


class TestXDist(unittest.TestCase):
    """Tests for xdist objects.
    """
    def test_xexpon(self):
        pass

    def test_xgamma(self):
        F = rmi.xdist.xgamma(4, 2)
        self.assertAlmostEqual(F.if2(1.), 0.07825557100854585)
        self.assertAlmostEqual(F.ixf2(1.), 0.06696082973267764)


def test_suite():
    """Returns testsuite for posterior probabilities.
    """
    return unittest.TestSuite([TestXDist(test)
                               for test in ["test_xexpon",
                                            "test_xgamma"]])
