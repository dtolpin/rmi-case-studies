"""Tests for posterior probabilities
"""

import unittest
import scipy.stats
import rmi.pp
import rmi.xdist


class TestPP(unittest.TestCase):
    """Tests for posterior probabilities
    """
    def test_marginal(self):
        """Marginal probabilities.  """
        t_s = 0
        t_e = 10
        pi = 0.5
        F = rmi.xdist.xgamma(a=2, scale=1)
        G = rmi.xdist.xexpon(scale=1)

        P = rmi.pp.marginal(t_s, t_e, [], pi, F, G)
        self.assertEqual(P, [], "empty sequence")

        P = rmi.pp.marginal(t_s, t_e, [(5., 1.)], pi, F, G)
        self.assertEqual(len(P), 1, "singleton")

        P = rmi.pp.marginal(t_s, t_e, [(10. / 3., 1.), (20. / 3., 1.)],
                            pi, F, G)
        self.assertAlmostEqual(P[0], P[1], msg="equal probabilities")

        P = rmi.pp.marginal(t_s, t_e, [(3., 1.), (6., 2.)],
                            pi, F, G)
        self.assertEqual(len(P), 2, "general case")


def test_suite():
    """Returns testsuite for posterior probabilities.
    """
    return unittest.TestSuite([TestPP(test)
                               for test in ["test_marginal"]])
