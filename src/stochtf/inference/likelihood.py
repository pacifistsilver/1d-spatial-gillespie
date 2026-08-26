"""Exact log-likelihood of observed counts under the promoter models.

The stationary distribution is available in closed form from
:mod:`stochtf.analytical.pgf`, so for iid counts the log-likelihood is

    log L(theta) = sum_i log P(y_i | theta)

evaluated deterministically over the whole distribution. There is no ABC
tolerance to tune and no simulation noise, so the sampler targets the true
posterior.
"""

import numpy as np

from stochtf.analytical import pgf

#: Probability floor. Keeps one unlucky observation from sending the whole
#: log-likelihood to -inf and stalling the sampler.
MIN_PROB = 1e-300


def prepare_counts(counts):
    """Validates observed counts and returns them as integers.

    Args:
        counts: Molecule numbers, one per cell. Flattened if not already 1-D.

    Returns:
        A 1-D int64 array of non-negative counts.

    Raises:
        ValueError: If counts are non-finite, non-integer, or negative.
    """
    y = np.asarray(counts)
    if y.ndim > 1:
        y = y.ravel()
    if not np.all(np.isfinite(y)):
        raise ValueError("counts contain non-finite values")
    y_int = np.rint(y).astype(np.int64)
    if np.any(np.abs(y - y_int) > 1e-8):
        raise ValueError("counts must be integers (molecule numbers)")
    if np.any(y_int < 0):
        raise ValueError("counts must be non-negative")
    return y_int


def log_likelihood(counts, a_s, b_s, a_n, b_n, k_y, gamma, model="dimer",
                   y_max=None):
    """Computes the exact log-likelihood of iid stationary counts.

    Args:
        counts: Molecule numbers, one per cell.
        a_s: SOX2 binding rate.
        b_s: SOX2 unbinding rate.
        a_n: NANOG binding rate.
        b_n: NANOG unbinding rate.
        k_y: Transcription rate in the active states.
        gamma: mRNA degradation rate.
        model: Promoter logic, selected via
          :data:`stochtf.analytical.pgf.MODEL_GATE`.
        y_max: Count grid bound. Chosen from the first two moments if None.

    Returns:
        The log-likelihood, or -inf for parameters outside the model's support.
        Out-of-support parameters return rather than raise so that a sampler can
        simply reject them.
    """
    y = prepare_counts(counts)

    if not all(np.isfinite([a_s, b_s, a_n, b_n, k_y, gamma])):
        return -np.inf
    if min(a_s, b_s, a_n, b_n, k_y) <= 0 or gamma <= 0:
        return -np.inf

    # The grid has to reach the largest observation, whatever the mean implies.
    if y_max is None:
        gate = pgf.MODEL_GATE[model]
        mean, second = pgf.factorial_moments(a_s, b_s, a_n, b_n, k_y, gamma,
                                             gate, order=2)
        var = second + mean - mean**2
        y_max = int(mean + 12.0 * np.sqrt(max(var, 1.0)) + 12)
    y_max = max(int(y.max()), y_max, 32)

    try:
        p = pgf.model_pmf(model, a_s, b_s, a_n, b_n, k_y, gamma, y_max=y_max)
    except (FloatingPointError, np.linalg.LinAlgError):
        return -np.inf

    if not np.all(np.isfinite(p)):
        return -np.inf

    return float(np.sum(np.log(np.maximum(p[y], MIN_PROB))))


def log_likelihood_factory(counts, model="dimer", k_y=None, gamma=None):
    """Builds a log-likelihood function bound to fixed data.

    The returned function takes its rates in prior-declaration order
    (alpha_s, alpha_n, beta_s, beta_n), which differs from the order
    :func:`log_likelihood` takes. Reordering happens here rather than at every
    call site.

    Args:
        counts: Molecule numbers, one per cell.
        model: Promoter logic to score against.
        k_y: Transcription rate to pin. Inferred if None.
        gamma: mRNA degradation rate to pin. Inferred if None.

    Returns:
        A callable taking four rates when k_y and gamma are both pinned, and six
        otherwise. Its arity is also exposed as the ``n_params`` attribute.
    """
    y = prepare_counts(counts)
    fixed = k_y is not None and gamma is not None

    if fixed:
        def f(a_s, a_n, b_s, b_n):
            return log_likelihood(y, a_s, b_s, a_n, b_n, k_y, gamma, model)
    else:
        def f(a_s, a_n, b_s, b_n, k_y_, gamma_):
            return log_likelihood(y, a_s, b_s, a_n, b_n, k_y_, gamma_, model)

    f.n_params = 4 if fixed else 6
    return f
