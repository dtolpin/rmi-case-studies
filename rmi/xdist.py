import scipy.stats
from scipy.special import gamma, gammaincc

"""Extended distributions.

Distribution objects extended with methods required to compute P, Q, and R.
"""


def gamma_incomplete(a, scale):
    """Gamma incomplete as defined in Maxima.
    """
    return gamma(a) * gammaincc(a, scale)
    

class _Xdist(object):
    """Extended distribution wrapper.
    """
    def __init__(self, dist, if2, ixf2):
        self.dist = dist
        self.if2 = if2
        self.ixf2 = ixf2
        self.meanpdf = self.if2(0.)
        
    def __getattr__(self, name):
        return getattr(self.dist, name)
    
def xexpon(scale):
    """Extended exponential distribution.
    """
    return _Xdist(scipy.stats.expon(scale=scale),
                  if2=lambda x: scipy.exp(- x / scale) / (2. * scale),
                  ixf2=lambda x: scipy.exp(- x / scale) / (4. * scale))

def xgamma(a, scale):
    """Extended gamma distribution.
    """
    return _Xdist(scipy.stats.gamma(a=a, scale=scale),
                  if2=lambda x: gamma_incomplete(2 * a   - 1, 2 * x / scale) /
                                (2 ** (2 * a - 1) * gamma(a) * scale *
                                 gamma_incomplete(a, x / scale)),
                  ixf2=lambda x: (gamma_incomplete(2 * a, 2 * x / scale) -
                                  2 * x / scale * gamma_incomplete(2 * a - 1, 2 * x / scale)) /
                                  (2 ** (2 * a) * gamma(a) * scale *
                                   (gamma_incomplete(a + 1, x / scale) -
                                    x / scale * gamma_incomplete(a, x / scale))))
