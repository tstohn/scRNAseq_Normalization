"""NegativeBinomial module."""

import numpy as np
from scipy.stats import nbinom
from copulas.univariate.base import BoundedType, ParametricType, ScipyModel

class NegBinomUnivariate(ScipyModel):
    """Wrapper around scipy.stats.gamma.

    Documentation: https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gamma.html
    """

    PARAMETRIC = ParametricType.PARAMETRIC
    BOUNDED = BoundedType.SEMI_BOUNDED
    MODEL_CLASS = nbinom

    def _fit_constant(self, X):
        self._params = {
            'n': 5.0,
            'p': 0.5,
            'loc': np.unique(X)[0]
        }

    def _fit(self, X):
        n, p, loc = nbinom.fit(X)
        self._params = {
            'n': n,
            'p': p,
            'loc': loc
        }

    def _is_constant(self):
        return False

    def _extract_constant(self):
        return self._params['loc']