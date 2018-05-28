"""Tests for EM algorithm to estimate process parameters.
"""

import unittest
import rmi.em
import scipy.stats


class TestEM(unittest.TestCase):
    """Tests for EM algorithm to estimate process parameters.
    """

    def test_short(self):
        """The sequence is too short.
        """
        with self.assertRaises(AssertionError):
            rmi.em.parameters(0, 10, [], 0.5)

        with self.assertRaises(AssertionError):
            rmi.em.parameters(0, 10, [(5, 1)], 0.5)

        with self.assertRaises(AssertionError):
            rmi.em.parameters(0, 15, [(5, 1), (10, 2)], 0.5)

    def test_none(self):
        """No intrusion.
        """
        F, G, labels = rmi.em.parameters(
            0, 10,
            [(1, 1), (3, 3), (6, 2), (9, 1)], 0.)
        self.assertEqual(labels, [0, 0, 0, 0])

    def test_all(self):
        """All events belong to the intrusion.
        """
        F, G, labels = rmi.em.parameters(
            0, 10,
            [(1, 1), (3, 3), (6, 2), (9, 1)], 1.)
        self.assertTrue(sum(labels) <= 2,
                        "at most a half of events can be marked")

    def test_some(self):
        """Some events belong to the intrusion.
        """
        F, G, labels = rmi.em.parameters(
            0, 10,
            [(2, 1), (3.5, 3), (4, 20), (9, 1)], 0.5)
        self.assertTrue(sum(labels) > 0,
                        "at least one positive label is anticipated")
        self.assertTrue(sum(labels) <= 2,
                        "at most a half of events can be marked")


def test_suite():
    """Returns testsuite for posterior probabilities.
    """
    return unittest.TestSuite([TestEM(test)
                               for test in ["test_short",
                                            "test_none",
                                            "test_all",
                                            "test_some"]])
