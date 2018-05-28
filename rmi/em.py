"""Expectaction-maximization algorithm for
intrusion detection.
"""

import math
import scipy.stats
import rmi.map
from rmi.xdist import xexpon, xgamma


def parameters(t_s, t_e, S, pe, Kmax=None, max_iter=0):
    """Estimates parameters of the process from a sequence
    with possible intrusion, using the EM algorithm.
    Assumes Gamma-distributed interarrival intervals
    and exponentially distributed marks. Arguments:
      t_s      --- start of interval
      t_e      --- end of interval
      S        --- sequence of events (t_i, y_i)
      pe       --- prior per-event intrusion probability
      Kmax     --- maximum number of intrusion events,
                   len(S)/2 by default.
      max_iter --- maximum number of iterations,
                   unbounded by default.
    Returns distributions of intervals (F) and marks (G), MAP labels.
    F.kwds, G.kwds hold distribution parameters.
    """
    N = len(S)
    assert N >= 3, "S must contain at least 3 events"

    if Kmax is None:
        Kmax = N // 2

    # Section 6.2, Algorithm 6

    F = None
    G = None
    labels = [0] * N
    iter = 1
    while True:
        # Extract normal subsequence according to current labels
        Sn = [s for l, s in zip(labels, S) if l == 0]

        # Fit parameters
        dts = [(t2 - t1) for (t1, _), (t2, _) in zip(Sn[:-1], Sn[1:])]
        tshape, _, tscale = scipy.stats.gamma.fit(dts, floc=0)
        _, yscale = scipy.stats.expon.fit([y for t, y in Sn], floc=0)

        # Create distribution objects
        F = xgamma(a=tshape, scale=tscale)
        G = xexpon(scale=yscale)

        # Compute new labels
        new_labels = rmi.map.labels(t_s, t_e, S, pe, F, G)

        # Stop when converged or stuck
        if (new_labels == labels or
                sum(new_labels) > Kmax or
                iter == max_iter):
            break

        labels = new_labels

    return F, G, labels
