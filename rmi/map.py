"""Maximum a posteriori subsequence.
"""


def labels(t_s, t_e, S, pe, F, G):
    """Computes maximum aposteriori intrusion
    subsequence of the sequence given parameters.
    Returns the sequence as a list of pairs (t, y).
    Arguments:
      t_s --- start of interval
      t_e --- end of interval
      S   --- sequence of events (t_i, y_i)
      pe  --- prior per-event intrusion probability
      F   --- interarrival distribution
      G   --- mark distribution
    return []
    Returns the label vector: 1 - intrusion, 0 - normal behavior.
    """
    N = len(S)

    # Section 5.3

    # Allocate space for probabilities and back-links
    a = [None] * (N + 1)
    prev = [None] * (N + 1)   # back-linked list of normal events

    # Find the most likely predecessor of each event
    fnorm = F.if2(0.)
    gnorm = G.if2(0.)
    for k in range(N):
        tk, yk = S[k]
        a[k] = pe ** k * F.if2(tk - t_s) / F.meanpdf
        prev[k] = -1
        for j in range(k):
            tj = S[j][0]
            aj = a[j] * pe ** (k - j - 1) * F.pdf(tk - tj) / F.meanpdf
            if aj > a[k]:
                a[k] = aj
                prev[k] = j
        a[k] *= (1 - pe) * G.pdf(yk) / G.meanpdf

    # Find the most likely last event
    a[N] = pe ** N * F.ixf2(t_e - t_s) / F.meanpdf
    prev[N] = -1
    for j in range(N):
        tj = S[j][0]
        aj = a[j] * pe ** (N - j - 1) * F.if2(t_e - tj) / F.meanpdf
        if aj > a[N]:
            a[N] = aj
            prev[N] = j

    # Fill the list of labels by tracing the back links
    labels = [1] * N
    k = prev[N]
    while k != -1:
        labels[k] = 0
        k = prev[k]

    return labels
