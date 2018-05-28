"""Posterior probabilities
"""


def _ps(t_s, t_e, S, pe, F, G):
    """Marginal likelihood of S. Arguments:
      t_s --- start of interval
      t_e --- end of interval
      S   --- sequence of events (t_i, y_i)
      pe  --- prior per-event intrusion probability 
      F   --- interarrival distribution
      G   --- mark distribution
    Returns the marginal likelihood of S and a[], as a tuple.
    """
    N = len(S)

    # Section 5.3, Lemma 2

    a = [None] * N
    for k in range(N):
        tk, yk = S[k]
        a[k] = pe ** k * F.if2(tk - t_s) / F.meanpdf
        for j in range(k):
            tj = S[j][0]
            a[k] += a[j] * pe ** (k - j - 1) * F.pdf(tk - tj) / F.meanpdf
        a[k] *= (1 - pe) * G.pdf(yk) / G.meanpdf
    A = pe ** N * F.ixf2(t_e - t_s) / F.meanpdf
    for j in range(N):
        tj = S[j][0]
        A += a[j] * pe ** (N - j - 1) * F.if2(t_e - tj) / F.meanpdf

    return A, a


def _psnoi(t_s, t_e, S, pe, F, G):
    """Probability of S without intrusion. Arguments:
      t_s --- start of interval
      t_e --- end of interval
      S   --- sequence of events (t_i, y_i)
      pe  --- prior intrusion probability
      F   --- interarrival distribution
      G   --- mark distribution
    Returns the probability.
    """
    N = len(S)
    if N == 0:
        return F.ixf2(t_e - t_s) / F.mean() / F.meanpdf

    # Section 5.2, Lemma 2

    t = S[0][0]
    Psnoi = F.if2(t - t_s) / F.meanpdf
    for k in range(N - 1):
        t1, y1 = S[k]
        t2 = S[k + 1][0]
        Psnoi *= (1 - pe) * \
            F.pdf(t2 - t1) / F.meanpdf * \
            G.pdf(y1) / G.meanpdf
    t, y = S[-1]
    Psnoi *= (1 - pe) * \
             F.if2(t_e - t) / F.meanpdf * \
             G.pdf(y) / G.meanpdf

    return Psnoi


def intrusion(t_s, t_e, S, pe, F, G):
    """Computes the probability of an intrusion
    in the sequence, given the parameters. Arguments:
      t_s --- start of interval
      t_e --- end of interval
      S   --- sequence of events (t_i, y_i)
      pe  --- prior intrusion probability
      F   --- interarrival distribution
      G   --- mark distribution
    Returns the probability.
    """
    N = len(S)
    if N == 0:
        # An empty sequence cannot contain an intrusion
        return 0.

    # Compute probability of the sequence without intrusion
    Psnoi = _psnoi(t_s, t_e, S, pe, F, G)

    # Compute marginal likelihood of S
    Ps, _ = _ps(t_s, t_e, S, pe, F, G)

    # Section 5.4, Equation 17
    return 1. - Psnoi / Ps


def marginal(t_s, t_e, S, pe, F, G):
    """Computes the marginal probability of each
    event in the sequence to belong to the intrusion.
    Arguments:
      t_s --- start of interval
      t_e --- end of interval
      S   --- sequence of events (t_i, y_i)
      pe  --- prior intrusion probability
      F   --- interarrival distribution
      G   --- mark distribution
    Returns the list of probabilities.
    """
    N = len(S)
    if N == 0:
        return []

    # Section 5.4, Theorem 5

    A, af = _ps(t_s, t_e, S, pe, F, G)

    t_sb = -t_e
    t_eb = -t_s
    Sb = [(-t, y) for t, y in reversed(S)]

    _, ab = _ps(t_sb, t_eb, Sb, pe, F, G)
    ab.reverse()

    return [(1 - af[k] * ab[k] / ((1 - pe) * G.pdf(y) / G.meanpdf * A))
            for k, (_, y) in enumerate(S)]
